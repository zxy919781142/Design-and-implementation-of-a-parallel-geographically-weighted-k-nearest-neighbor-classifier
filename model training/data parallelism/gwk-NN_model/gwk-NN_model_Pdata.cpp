#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "mpi.h"
#include "mpe.h"
#include "sstream"
using namespace std;
float Compute_Distance(float Raster_x,float Raster_y,float Sample_x,float Sample_y,float resulation)//计算两像元的空间距离
{
	return(sqrt((Raster_x-Sample_x)*(Raster_x-Sample_x)+(Raster_y-Sample_y)*(Raster_y-Sample_y))*resulation);
}
float Compute_SpecDist(float *RastScanline,long Raster_N,long Sample_N)//计算两像元的光谱距离
{
	return(fabs(RastScanline[Raster_N]-RastScanline[Sample_N]));
}

	
void  Cmpute_XY(long Raster_N,long *Raster_X,long *Raster_Y,long X_size)
{
	*Raster_X =  Raster_N%X_size;
	*Raster_Y =  Raster_N/X_size;
}
	
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
	return ((float)Raster_Weight/Raster_Dis_Weight);
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
float qiuxing(float *x, float xdata)
{
	
	if(xdata>1/x[0])
		return(1-x[1]);
	else
		return(1-x[1]*(3/2*x[0]*xdata-0.5*(x[0]*xdata)*(x[0]*xdata)*(x[0]*xdata)));
		
}

float gaosi(float *x, float xdata)
{
	return(1-x[0]*(1-exp(-3*(xdata/x[1])*(xdata/x[1]))));
}

float zhishu(float *x,float xdata)
{
	
	return(1-x[0]*(1-1/(exp(xdata/x[1]))));
	
}


float ydata_G(long h,long raster_X,long raster_Y,long *sample_X,long *sample_Y,long *zsample_X,long *zsample_Y,long n_sample,long n_zsample)//使用k近邻方法计算像元归属概率值
{
	long sampleN = 0,zsampleN=0;

	for(int i =0;i<n_sample;i++)
	if(abs(raster_X-sample_X[i])<h||abs(raster_Y-sample_Y[i])<h)
		sampleN++;
	for(int i =0;i<n_zsample;i++)
	if(abs(raster_X-zsample_X[i])<h||abs(raster_Y-zsample_Y[i])<h)
		zsampleN++;
	
	return((float)sampleN/zsampleN);
}
bool compute_ydata(float *ydata,long *xdata,long raster_Xsize,long raster_Ysize,long *raster_X,long *raster_Y,long *sample_x,long *sample_Y,long *zsample_X,long *zsample_Y,long n_sample,long n_zsample,int px_i)//计算某类样本归属概率值,获得xdata,ydata
{
	ofstream write; 
	if(px_i==0){
    write.open("C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\result\\bareland_1.csv");
	}
	if(px_i==1){
    write.open("C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\result\\fieldland_1.csv"); 
	}
	if(px_i==2){
    write.open("C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\result\\grassland_1.csv"); 
	}
	if(px_i==3){
    write.open("C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\result\\roadland_1.csv"); 
	}
	if(px_i==4){
    write.open("C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\result\\water_1.csv"); 
	}

	//*************************************************************
	long lh = 0;
	float data_g=0;
	if(raster_Xsize>raster_Ysize)
		lh = raster_Xsize;
	else
	lh = raster_Ysize;
	for(int i = 0;i<n_sample;i++)
	{
		float maxdata_g = 0;
		for(int h =1;h<lh;h = h+1)
		{
			data_g = ydata_G( h,raster_X[i],raster_Y[i],sample_x,sample_Y,zsample_X,zsample_Y,n_sample, n_zsample);
			if(data_g>maxdata_g)
			{
				maxdata_g = data_g;
				xdata[i] = h;
			}
		}
		ydata[i] = maxdata_g;
		write<<i<<","<<xdata[i]<<","<<ydata[i]<<endl;
	}

 
    write.close(); 




	return true;

}



void resnorm(long *xdata,float *ydata,long n_sample,int cp)	
{
	float x[] = {0.0001,0.0001};
	float qx1=0,qx2=0,gx1=0,gx2=0,zx1=0,zx2=0;
	float zresnorm1 = 0,zresnorm2 = 0,zresnorm3 = 0,zresnorm_min[] = {100000,100000,100000};

	for(;x[0]<1;x[0] = x[0]+0.001)
	{
		
		for(x[1]= 0.0001;x[1]<361;x[1] = x[1] + 0.5)
		{
			zresnorm1 = 0,zresnorm2 = 0,zresnorm3 = 0;
			for(int n = 0;n<n_sample;n++)
			{
			//	zresnorm1 =zresnorm1 + fabs(ydata[n] - qiuxing(x,xdata[n]));
				zresnorm2 =zresnorm2 + fabs(ydata[n] - gaosi(x,xdata[n]));
				zresnorm3 =zresnorm3 + fabs(ydata[n] - zhishu(x,xdata[n]));
			
			if(zresnorm_min[1]>zresnorm2)
			{
				zresnorm_min[1] = zresnorm2;
				gx1 = x[0];
				gx2 = x[1];
			}
			if(zresnorm_min[2]>zresnorm3)
			{
				zresnorm_min[2] = zresnorm3;
				zx1 = x[0];
				zx2 = x[1];
			}				
			
		}
	}

	cout<<"process: "<<cp<<" "<<"size: "<<n_sample<<" "<<"gsmin, zsmin, gc, gr, zc, zr: "<<" "<<zresnorm_min[1]<<" "<<zresnorm_min[2]<<" "<<gx1<<" "<<gx2<<" "<<zx1<<" "<<zx2<<endl;
	
	
}


int main(int nArgc, char *papszArgv[])
{
    GDALDataset  *poDataset;
	
	const char *pszFilename = "C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\nanjing2002.tif";


	GDALAllRegister();
	int N_sample = 5;
	
	int cp, np;
	MPI_Init(&nArgc, &papszArgv);//并行初始化
	MPI_Comm_rank(MPI_COMM_WORLD, &cp);//获取进程编号
	MPI_Comm_size(MPI_COMM_WORLD, &np);//获取进程数量
	



	
	const char *pszFilenameX[] = {"C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\bareland.txt","C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\fieldland.txt","C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\grassland.txt","C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\roadland.txt","C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\waterland.txt"};

	const char *zpszFilename = "C:\\Users\\zxy\\Desktop\\Gwk-NN_Program\\remote_data\\sample\\ztestsample.txt" ;
	long raster_x[5][300],raster_y[5][500],Sample_N[5][500],zraster_x[500],zraster_y[500],ZSample_N[500];



	
	
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    if( poDataset == NULL )
    {
        printf("open failed!\n");
		MPI_Finalize();
		return 0;
    }
	
	double        adfGeoTransform[6];

	const char *Driver_Name = poDataset->GetDriver()->GetDescription();
	long RasterXsize =  poDataset->GetRasterXSize();
	long RasterYsize =  poDataset->GetRasterYSize();
	long RasterCount =  poDataset->GetRasterCount();

	OGRSpatialReference oSRS;
	if( poDataset->GetProjectionRef()  != NULL )
	{
		oSRS = poDataset->GetProjectionRef();
	}
	float Origin_x,Origin_y;

	if( poDataset->GetGeoTransform( adfGeoTransform ) == CE_None )
    {
        
        Origin_x = adfGeoTransform[0];
		Origin_y = adfGeoTransform[3] ;        
    }
	
	int point_n = 0;
	int zsample_point_n = fopen_txt_n(zpszFilename,zraster_x,zraster_y,ZSample_N,RasterXsize);
	float ydata[500];long xdata[500];
	int Nxp= 0,p = 1,px_i = 0;



	double start_time=MPI_Wtime();

	
	if(np >= N_sample &&cp+1<=N_sample)
	{
		//MPE_Log_event( eventID_begin, 0, "start");
		point_n = fopen_txt_n(pszFilenameX[cp],raster_x[cp],raster_y[cp],Sample_N[cp],RasterXsize);
		if(point_n == 0)
		{
			printf("获取样本数据失败!\n");
			MPI_Finalize();
			return 0;
		}	

		compute_ydata(ydata,xdata,RasterXsize,RasterYsize,raster_x[cp],raster_y[cp],raster_x[cp],raster_y[cp],zraster_x,zraster_y,point_n,zsample_point_n,cp);//训练样本点随机概率计算
		if(!compute_ydata)
		{
			printf("计算样本点随机概率失败!\n");
			MPI_Finalize();
			return 0;
		}
		resnorm(xdata,ydata, point_n,cp);
	}
	else if(np < N_sample )
	{
   
	loop: 
		px_i = Nxp +cp;
		
		point_n = fopen_txt_n(pszFilenameX[px_i],raster_x[px_i],raster_y[px_i],Sample_N[px_i],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
		
		if(point_n == 0)
		{
			printf("获取样本数据失败!\n");
			MPI_Finalize();
			return 0;
		}
		

		
		compute_ydata(ydata,xdata,RasterXsize,RasterYsize,raster_x[px_i],raster_y[px_i],raster_x[px_i],raster_y[px_i],zraster_x,zraster_y,point_n,zsample_point_n,px_i);//训练样本点随机概率计算
		if(!compute_ydata)
		{
			printf("计算样本点随机概率失败!\n");
			MPI_Finalize();
			return 0;
		}

		resnorm(xdata,ydata, point_n,px_i);//模型拟合

		

		

		Nxp = N_sample-np*p;
		int Sxp = cp+np*p;
		
		if(Nxp <= np && Nxp>0)
		{
			if(Sxp<N_sample)
			{				
				Nxp = np*p;
				p++;
				goto loop;
			}
		}
			
		else if (Nxp >np)
		{
			Nxp = np*p;
			p++;
			goto loop;
		}		

	}
//	 MPE_Log_event( eventID_end, 0, "end");
	MPI_Barrier(MPI_COMM_WORLD);
	double end_time=MPI_Wtime();
	if(cp==0)
	printf("%d process consumes time:%f\n",np,end_time-start_time);


  
	MPI_Finalize();

	


}
