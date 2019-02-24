#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "windows.h"
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


	struct pp{
	 
		int a;
		float b;
	
	};

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



	
    poDataset = (GDALDataset *) GDALOpen( pszFilename, GA_ReadOnly );
    if( poDataset == NULL )
    {
        printf("open failed!\n");
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
	
    int cp, np;
	MPI_Init(&nArgc, &papszArgv);
	MPI_Comm_rank(MPI_COMM_WORLD, &cp);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	float t1 = 0,t2=0,t3=0;
	//float ts1=0,ts2 =0;
	float ts1=0;
	float ts2=0;
	float tsz=0;
	char* temp = (char *) malloc(2);

 	 int Bands = RasterCount;
	 float *pafScanline;
	 int point_n[5];

     GDALDataset *poDstDS;
	 GDALRasterBand  *poBand;
	 GDALRasterBand  *poDBand;
	 GDALDataType RasterDataType;
     poBand = poDataset->GetRasterBand(1);	
	 if(cp==0){
	 char **papszOptions = NULL;
	 char *pszSRS_WKT = NULL;
	 const char *pszDstFilename = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\nanjing2002_test.tif";
	 GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(Driver_Name);//获取源数据driver
		
	 GDALDataType RasterDataType = poBand->GetRasterDataType();	
     poDstDS = poDriver->Create( pszDstFilename, RasterXsize, RasterYsize, RasterCount, RasterDataType, 
                                papszOptions );
	 poDstDS->SetGeoTransform( adfGeoTransform );
	 oSRS.exportToWkt( &pszSRS_WKT );
     poDstDS->SetProjection( pszSRS_WKT );
     CPLFree( pszSRS_WKT );

	

	 }
	point_n[0] = fopen_txt_n(pszFilename1,raster_x[0],raster_y[0],Sample_N[0],RasterXsize);
	point_n[1] = fopen_txt_n(pszFilename2,raster_x[1],raster_y[1],Sample_N[1],RasterXsize);
	point_n[2] = fopen_txt_n(pszFilename3,raster_x[2],raster_y[2],Sample_N[2],RasterXsize);
	point_n[3] = fopen_txt_n(pszFilename4,raster_x[3],raster_y[3],Sample_N[3],RasterXsize);
	point_n[4] = fopen_txt_n(pszFilename5,raster_x[4],raster_y[4],Sample_N[4],RasterXsize);
	pafScanline = (float *) CPLMalloc(sizeof(float)*RasterXsize*RasterYsize);

	
	
	


    t1 = MPI_Wtime();


		
		
		long raster_pixles_n = RasterXsize*RasterYsize; 
		int raster_n =0;
		float pixles_weight1;


	for(; Bands>0; Bands--)
	{
		
		poBand = poDataset->GetRasterBand(Bands);

	
		
        poBand->RasterIO( GF_Read, 0, 0, RasterXsize, RasterYsize, pafScanline, RasterXsize, RasterYsize, GDT_Float32,  0, 0 );
	


	    raster_n =0;

		for(;raster_n <raster_pixles_n;raster_n++)
		{
			
		    //MPI_Status status;
			float pixles_weight[100];
			float pw[100];
			int Nxp= 0,p = 1,px_i = 0,N_sample=5, Sxp=0;
			p=1;
			
			if(np>=6 )
			{
				
			
				if(cp<=5 && cp>0 )
				{

               
				pixles_weight[cp-1] = Compute_Weight(raster_n,Sample_N[cp-1],point_n[cp-1], adfGeoTransform[1],RasterXsize,pafScanline);
				MPI_Send(&pixles_weight[cp-1],1,MPI_FLOAT,0,cp-1,MPI_COMM_WORLD);
			  
				}

				
			}			
			if (np>1 && np<=5)
			{
			    
			    if(cp>0 ){
					
				loop:
			 	   px_i = Nxp +cp-1;
				   if(px_i<=4){
				 
					   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);	
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,0,px_i,MPI_COMM_WORLD);
			  	
				   Nxp = 5-p*(np-1);
				            
				       if(Nxp>=np-1){
				        Nxp=(np-1)*p;
				         p++;
						 goto loop;
	                   }
				       else if(Nxp<np-1 && Nxp>0){
				   
				       Nxp=(np-1)*p;
				       p++;
					   goto loop;
				        }
					   
				   }
			   
				 }
		 
		
		}		 
				

	
			if (np==1)
			{

        pw[0] = Compute_Weight(raster_n,Sample_N[0],point_n[0], adfGeoTransform[1],RasterXsize,pafScanline);
		pw[1] = Compute_Weight(raster_n,Sample_N[1],point_n[1], adfGeoTransform[1],RasterXsize,pafScanline);
		pw[2] = Compute_Weight(raster_n,Sample_N[2],point_n[2], adfGeoTransform[1],RasterXsize,pafScanline);
		pw[3] = Compute_Weight(raster_n,Sample_N[3],point_n[3], adfGeoTransform[1],RasterXsize,pafScanline);
		pw[4] = Compute_Weight(raster_n,Sample_N[4],point_n[4], adfGeoTransform[1],RasterXsize,pafScanline);
		//cout<<pw[0]<<" "<<pw[1]<<" "<<pw[2]<<" "<<pw[3]<<" "<<pw[4]<<" "<<endl;
		 
			}
		
		MPI_Status status;


		MPI_Barrier(MPI_COMM_WORLD); 
		ts1=MPI_Wtime();
		if(cp==0 && np>1){
		for(int j=0;j<5  ;j++){

				
				
				MPI_Recv(&pixles_weight[j],1,MPI_FLOAT,MPI_ANY_SOURCE,j,MPI_COMM_WORLD,&status);
				



		}
		   ts2=MPI_Wtime();		
		   tsz=tsz+ts2-ts1;
		
		}

	


		
		
		if(cp==0 )

		{
		float max_weight=0;
		int sample_x = 10;
	


	
			for(int i = 0;i<5;i++){
			  // if(pw[i]>=0.9){pw[i]=1;}
			   if(pw[i]>max_weight)
				{
					max_weight =pw[i];
					sample_x = i;
				}
			}

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

	
	

		

	//	}

		}
	
		

			}


			if(cp==0){
			poDBand =poDstDS->GetRasterBand(Bands);
			poDBand->RasterIO( GF_Write, 0, 0, RasterXsize,  RasterYsize,  pafScanline, RasterXsize,  RasterYsize, GDT_Float32,  0, 0 );
			}  
				
}

			
	//	 printf("bands = %d\n",Bands);
    
	t2 = MPI_Wtime();

	if(cp==0){


	printf("consume time = %f",t2-t1);

	printf("send time = %f",tsz);
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize();

	 if( poDstDS != NULL )
	 GDALClose( (GDALDatasetH) poDstDS );
	return 0;
}
