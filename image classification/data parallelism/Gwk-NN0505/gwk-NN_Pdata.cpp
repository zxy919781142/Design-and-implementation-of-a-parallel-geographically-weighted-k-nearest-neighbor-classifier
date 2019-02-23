#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "mpi.h"
using namespace std;
float Compute_Distance(float Raster_x,float Raster_y,float Sample_x,float Sample_y,float resulation)//计算两像元的空间距离
{
	return(sqrt((Raster_x-Sample_x)*(Raster_x-Sample_x)+(Raster_y-Sample_y)*(Raster_y-Sample_y))*resulation);
}
float Compute_SpecDist(float *RastScanline,long Raster_N,long Sample_N)//计算两像元的光谱距离
{
	return(fabs(RastScanline[Raster_N]-RastScanline[Sample_N]));
}

	/****************/
	/*计算X，Y坐标*/
	/****************/
void  Cmpute_XY(long Raster_N,long *Raster_X,long *Raster_Y,long X_size)
{
	*Raster_X =  Raster_N%X_size;
	*Raster_Y =  Raster_N/X_size;
}
	/****************/
	/*根据图像坐标计算存储位置Raster_N*/
	/****************/
long Compute_Raster_N(long Raster_X,long Raster_Y,long X_size)
{
	return(Raster_X+Raster_Y*X_size);
}

float Compute_Weight(long Raster_N,long *Sample_N,long Sample_n, float resulation,long X_size,float *RastScanline)//计算某个待分类像元归属权重
{
	float Raster_Weight = 0,Raster_Dis_Weight = 0;
	long Raster_x,Raster_y,Sample_x,Sample_y;

	Cmpute_XY(Raster_N,&Raster_x,&Raster_y,X_size);
	
	for(long i = 0;i < Sample_n;i++)
	{
		Cmpute_XY(Sample_N[i],&Sample_x,&Sample_y,X_size);
		Raster_Weight = Raster_Weight +(float)1/(Compute_Distance(Raster_x,Raster_y,Sample_x,Sample_y,resulation)*
			Compute_SpecDist(RastScanline,Raster_N,Sample_N[i]));
		Raster_Dis_Weight = Raster_Dis_Weight + 1/Compute_Distance(Raster_x,Raster_y,Sample_x,Sample_y,resulation);
	} 
		if((float)Raster_Weight/Raster_Dis_Weight<1){
		return ((float)Raster_Weight/Raster_Dis_Weight);}
	else return 1;
}

int fopen_txt_n(const char *fileName,long *raster_x,long *raster_y,long *Sample_N,long X_size)
{
	fstream out;
	char buffer[100];
	char buffer1[20];
	char buffer2[20];
	out.open(fileName,ios::in);
	int k = 0;
	while(!out.eof())
	{
		out.getline(buffer,100,'\n');
		int i = 0;
		while(buffer[i]!='	')
		{
			buffer1[i] = buffer[i];
			i++;
		}
		buffer1[i] = ' ';
		i++;
		int j = 0;
		while(buffer[i]!='\0')
		{
			buffer2[j] = buffer[i];
			i++;
			j++;
		}		
		buffer2[j] = ' ';


		long x = atof(buffer1);
		long y = atof(buffer2);
		raster_x[k] = x;
		raster_y[k] = y;
		Sample_N[k] =  Compute_Raster_N(raster_x[k],raster_y[k],X_size);
		k++;		
	}
	out.close();
	return (k);

}




int main(int nArgc, char *papszArgv[])
{
    GDALDataset  *poDataset;
	const char *pszFilename = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\nanjing2002.tif";
	GDALAllRegister();
	const char *pszFilename1 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\bareland.txt" ;
	const char *pszFilename2 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\fieldland.txt" ;
	const char *pszFilename3 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\grassland.txt" ;
	const char *pszFilename4 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\roadland.txt" ;
	const char *pszFilename5 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\waterland.txt" ;
	long raster_x[5][100],raster_y[5][100],Sample_N[5][100];

	int cp, np;
	MPI_Init(&nArgc, &papszArgv);//并行初始化
	MPI_Comm_rank(MPI_COMM_WORLD, &cp);//获取进程编号
	MPI_Comm_size(MPI_COMM_WORLD, &np);//获取进程数量
	double t1=0,t2=0,tz = 0;


	 /****************/
	/*打开栅格数据源*/
	/****************/
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    if( poDataset == NULL )
    {
        printf("open failed!\n");
    }
	/*********************/
	/*获取栅格数据源信息*/
	/********************/
	double        adfGeoTransform[6];

	const char *Driver_Name = poDataset->GetDriver()->GetDescription();//获取源数据driver name
	long RasterXsize =  poDataset->GetRasterXSize();//获取图像横向像素数
	long RasterYsize =  poDataset->GetRasterYSize();//获取图像纵向像素数
	long RasterCount =  poDataset->GetRasterCount();//获取图像层数

	OGRSpatialReference oSRS;
	if( poDataset->GetProjectionRef()  != NULL )//获取投影信息
	{
		oSRS = poDataset->GetProjectionRef();
	}
	float Origin_x,Origin_y;

	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )//获取图像起始点（左上角点）坐标，及分辨率
    {
        
        Origin_x = adfGeoTransform[0];
		Origin_y = adfGeoTransform[3] ;        
    }
	/*
	adfGeoTransform[0] 左上角X坐标
    adfGeoTransform[1] 东西方向像元分辨率
    adfGeoTransform[2] 旋转角度，如果是0则图像上为正北方向
    adfGeoTransform[3] 左上角Y坐标
    adfGeoTransform[4] 旋转角度，如果是0则图像上为正北方向
    adfGeoTransform[5] 南北方向像元分辨率
	*/
	/*********************/
	/*创建目标数据*/
	/********************/

	 GDALDataset *poDstDS;
	 char **papszOptions = NULL;
	 char *pszSRS_WKT = NULL;
	// GDALDriver *pDriver = GetGDALDriverManager()->GetDriverByName("GTiff"); //图像驱动
   //  papszOptions = CSLSetNameValue(papszOptions, "BIGTIFF", "IF_NEEDED"); //配置图像信息

	 const char *pszDstFilename1 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test1.tif";
	 const char *pszDstFilename2 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test2.tif";
	 const char *pszDstFilename3 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test3.tif";
	 const char *pszDstFilename4 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test4.tif";
	 const char *pszDstFilename5 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test5.tif";
	 const char *pszDstFilename6 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test6.tif";
	 const char *pszDstFilename7 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test7.tif";
	 const char *pszDstFilename8 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test8.tif";
	// GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(Driver_Name);//获取源数据driver
	 GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(Driver_Name);//获取源数据driver
	 GDALRasterBand  *poBand = poDataset->GetRasterBand(1);	
	 GDALRasterBand  *poDBand;
	 GDALDataType RasterDataType = poBand->GetRasterDataType();	 

	 switch(cp)
		{
			case 0:
			{
				poDstDS = poDriver->Create(pszDstFilename1, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 1:
			{
				poDstDS = poDriver->Create(pszDstFilename2, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 2:
			{
				poDstDS = poDriver->Create(pszDstFilename3, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 3:
			{
				poDstDS = poDriver->Create(pszDstFilename4, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 4:
			{
				poDstDS = poDriver->Create(pszDstFilename5, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}

			case 5:
			{
				poDstDS = poDriver->Create(pszDstFilename6, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 6:
			{
				poDstDS = poDriver->Create(pszDstFilename7, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}
			case 7:
			{
				poDstDS = poDriver->Create(pszDstFilename8, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
				break;
			}

			default:
			{
			}
		}
   //  poDstDS = poDriver->Create(pszDstFilename, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
                               
	 poDstDS->SetGeoTransform( adfGeoTransform );
	 oSRS.exportToWkt( &pszSRS_WKT );
     poDstDS->SetProjection( pszSRS_WKT );
     CPLFree( pszSRS_WKT );

	int bare_point_n = fopen_txt_n(pszFilename1,raster_x[0],raster_y[0],Sample_N[0],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	int field_point_n = fopen_txt_n(pszFilename2,raster_x[1],raster_y[1],Sample_N[1],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	int grass_point_n = fopen_txt_n(pszFilename3,raster_x[2],raster_y[2],Sample_N[2],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	int road_point_n = fopen_txt_n(pszFilename4,raster_x[3],raster_y[3],Sample_N[3],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	int water_point_n = fopen_txt_n(pszFilename5,raster_x[4],raster_y[4],Sample_N[4],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数

	

	/*********************/
	/*读源数据并写入目标数据*/
	/********************/
	int Bands = RasterCount;
	int p_ysize = RasterYsize/np;
	int y_size = p_ysize*cp;
	int p_lastYsize = RasterYsize -y_size;
	
	long raster_pixles_n = RasterXsize*RasterYsize;

		float *pafScanline;
		pafScanline = (float *) CPLMalloc(sizeof(float)*RasterXsize*RasterYsize);

	t1 = MPI_Wtime();
	for(; Bands>0; Bands--)
	{


		poBand = poDataset->GetRasterBand(Bands);
		
		if(cp!=np-1)
		{
			
			poBand->RasterIO( GF_Read, 0, y_size, RasterXsize, p_ysize, 
							  pafScanline, RasterXsize, p_ysize, GDT_Float32, 
							  0, 0 );
			raster_pixles_n =RasterXsize*p_ysize;
		}
		 if(cp==np-1)
		{
			
			
			poBand->RasterIO( GF_Read, 0, y_size, RasterXsize, p_lastYsize, 
							  pafScanline, RasterXsize, p_lastYsize, GDT_Float32, 
							  0, 0 );
			raster_pixles_n =RasterXsize*p_lastYsize;
		}

		

		/*基于像元的图像处理*/
		float pixles_weight[5];		
		
		long raster_n = 0;

		for(;raster_n < raster_pixles_n;raster_n++)
		{

		pixles_weight[0] = Compute_Weight(raster_n,Sample_N[0],bare_point_n, adfGeoTransform[1],RasterXsize,pafScanline);
		pixles_weight[1] = Compute_Weight(raster_n,Sample_N[1],field_point_n, adfGeoTransform[1],RasterXsize,pafScanline);
		pixles_weight[2] = Compute_Weight(raster_n,Sample_N[2],grass_point_n, adfGeoTransform[1],RasterXsize,pafScanline);
		pixles_weight[3] = Compute_Weight(raster_n,Sample_N[3],road_point_n, adfGeoTransform[1],RasterXsize,pafScanline);
		pixles_weight[4] = Compute_Weight(raster_n,Sample_N[4],water_point_n, adfGeoTransform[1],RasterXsize,pafScanline);
		// cout<<pixles_weight[0]<<" "<<pixles_weight[1]<<" "<<pixles_weight[2]<<" "<<pixles_weight[3]<<" "<<pixles_weight[4]<<" "<<endl;
		float max_weight=0;
		int sample_x = 10;
		for(int i = 0;i<5;i++){
		   if(pixles_weight[i]>max_weight)
			{
				max_weight = pixles_weight[i];
				sample_x = i;
			}
		   
		}
		//cout<<"class:"<<sample_x<<endl;
		switch(sample_x)
		{
			case 0:
			{
				pafScanline[raster_n] = 10;
				break;
			}
			case 1:
			{
				pafScanline[raster_n] = 50;
				break;
			}
			case 2:
			{
				pafScanline[raster_n] = 100;
				break;
			}
			case 3:
			{
				pafScanline[raster_n] = 200;
				break;
			}
			case 4:
			{
				pafScanline[raster_n] = 250;
				break;
			}
			default:
			{
			}
		}




		

		

		

		}
		if(cp!=np-1)
		{
			poDBand = poDstDS->GetRasterBand(Bands);
			poDBand->RasterIO( GF_Write, 0, y_size, RasterXsize, p_ysize, 
                      pafScanline, RasterXsize,p_ysize, GDT_Float32, 0, 0 ); 
			
		}
		if(cp==np-1)
		{	
		poDBand = poDstDS->GetRasterBand(Bands);	
		poDBand->RasterIO( GF_Write, 0, y_size, RasterXsize, p_lastYsize, pafScanline, RasterXsize, p_lastYsize, GDT_Float32, 0, 0 ); 
			
		}
		



	


	}
	t2 = MPI_Wtime();
	if(cp==0)
	printf("cosume time: %f",t2-t1);
	if( poDstDS != NULL )
	 GDALClose( (GDALDatasetH) poDstDS );
	MPI_Finalize();
}