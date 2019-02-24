// Gwk-NN_pdual.cpp : 定义控制台应用程序的入口点。
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include "cpl_string.h"
#include "ogr_spatialref.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "mpi.h"
#include "windows.h"
using namespace std;
float Compute_Distance(float Raster_x,float Raster_y,float Sample_x,float Sample_y,float resulation)
{
	return(sqrt((Raster_x-Sample_x)*(Raster_x-Sample_x)+(Raster_y-Sample_y)*(Raster_y-Sample_y))*resulation);
}
float Compute_SpecDist(float *RastScanline,long Raster_N,long Sample_N)
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

float Compute_Weight(long Raster_N,long *Sample_N,long Sample_n, float resulation,long X_size,float *RastScanline)
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
		int c;
	
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
	int point_n[5];

	int cp, np;
	MPI_Init(&nArgc, &papszArgv);
	MPI_Comm_rank(MPI_COMM_WORLD, &cp);
	MPI_Comm_size(MPI_COMM_WORLD, &np);



	double t1=0,t2=0,tz = 0;
	double ts1=0;
	double ts2=0;
	double tsz=0;
    double tts1=0;
	double tts2=0;
	double ttsz=0;


	
	int pz,ys;
	ys=np%4;
	if(ys==0){
	pz=int(np/4);
	}
	else
	{
	pz=int(np/4)+1;
	}

	int zh;
	
	/****************/
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
	

	 GDALDataset *poDstDS;
	 char **papszOptions = NULL;
	 char *pszSRS_WKT = NULL;
	

	 const char *pszDstFilename1 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test1.tif";
	 const char *pszDstFilename2 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test2.tif";
	 const char *pszDstFilename3 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test3.tif";
	 const char *pszDstFilename4 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test4.tif";
	 const char *pszDstFilename5 = "C:\\Users\\Administrator\\Desktop\\Gwk-NN_Program\\remote_data\\test5.tif";
	
	 GDALRasterBand  *poBand = poDataset->GetRasterBand(1);	
	 GDALRasterBand  *poDBand;
	 GDALDataType RasterDataType = poBand->GetRasterDataType();	 
	
	 switch(cp)
		{
			case 0:
		
			
			{
				poDstDS = poDriver->Create(pszDstFilename1, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
					 poDstDS->SetGeoTransform( adfGeoTransform );
	         oSRS.exportToWkt( &pszSRS_WKT );
             poDstDS->SetProjection( pszSRS_WKT );
             CPLFree( pszSRS_WKT );

				break;
			}

			case 4:
			
		
			{
				poDstDS = poDriver->Create(pszDstFilename2, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
					 poDstDS->SetGeoTransform( adfGeoTransform );
	 oSRS.exportToWkt( &pszSRS_WKT );
     poDstDS->SetProjection( pszSRS_WKT );
     CPLFree( pszSRS_WKT );
				break;
			}

			case 8:
			{
				poDstDS = poDriver->Create(pszDstFilename3, RasterXsize, RasterYsize, RasterCount, RasterDataType,  papszOptions );
		        poDstDS->SetGeoTransform( adfGeoTransform );
	           oSRS.exportToWkt( &pszSRS_WKT );
                poDstDS->SetProjection( pszSRS_WKT );
              CPLFree( pszSRS_WKT );
				break;
			}
	default:
			{
			}
		

		}

  

	point_n[0] = fopen_txt_n(pszFilename1,raster_x[0],raster_y[0],Sample_N[0],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	point_n[1] = fopen_txt_n(pszFilename2,raster_x[1],raster_y[1],Sample_N[1],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	point_n[2] = fopen_txt_n(pszFilename3,raster_x[2],raster_y[2],Sample_N[2],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	point_n[3] = fopen_txt_n(pszFilename4,raster_x[3],raster_y[3],Sample_N[3],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	point_n[4] = fopen_txt_n(pszFilename5,raster_x[4],raster_y[4],Sample_N[4],RasterXsize);//获取样本数据坐标点及样本数据坐标点个数
	
	

	int Bands = RasterCount;
	zh=int(cp/4);
	int lb;


	int p_ysize = RasterYsize/pz;
	int y_size = p_ysize*zh;
	int p_lastYsize = RasterYsize -y_size;
	
	long raster_pixles_n = RasterXsize*RasterYsize;

	
		float *pafScanline;
		pafScanline = (float *) CPLMalloc(sizeof(float)*RasterXsize*RasterYsize);


		float pixles_weight[5];	
		float pw[5];	
	



		MPI_Barrier(MPI_COMM_WORLD); 
		t1 = MPI_Wtime();


	
		if(np>4){

				for(; Bands>1; Bands--)
	{

			poBand = poDataset->GetRasterBand(Bands);
			int raster_n =0;

          if( int(cp/4)<=(pz-2))
		{

			poBand->RasterIO( GF_Read, 0, y_size, RasterXsize, p_ysize,  pafScanline, RasterXsize, p_ysize, GDT_Float32, 0, 0 );
		
			raster_pixles_n =RasterXsize*p_ysize;
			
		
			for(;raster_n < raster_pixles_n ;raster_n++)
			{
				
				int Nxp= 0,p = 1,px_i = 0, Sxp=0;

				loop:
			 	   px_i = Nxp +cp;
				   if(px_i<=4){
				   //pixles_weight[px_i].a=px_i;
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				// cout<<"send:"<<raster_n<<" "<<pixles_weight[px_i]<<endl;
				   if(cp!=0){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,0,px_i,MPI_COMM_WORLD);
				   }
			       //tag2=tag*(np-1);
				   Nxp = 5-p*4;
				            
				       if(Nxp>=4){
				        Nxp=4*p;
				         p++;
						 goto loop;
	                   }
				       else if(Nxp<4 && Nxp>0){
				   
				       Nxp=4*p;
				       p++;
					   goto loop;
				        }
					   
				   }	



			MPI_Status status;
			
			if (cp==0 ){
				

				ts1=MPI_Wtime();    
				for(int j=1;j<4;j++){
             
					MPI_Recv(&pw[j],1,MPI_FLOAT,j,j,MPI_COMM_WORLD,&status);
					
					
						}


            ts2=MPI_Wtime();		
            tsz=tsz+ts2-ts1;	
            float max_weight=0;
	    int sample_x = 10;
		 for(int i = 0;i<5;i++){
	              if(pw[i]>max_weight)
	                 {
			   max_weight = pw[i];
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
				

		                }
					

			}

		}

			if(cp%4==0){

			
			}  


   }


  else {
	          poBand->RasterIO( GF_Read, 0, y_size, RasterXsize, p_lastYsize, pafScanline, RasterXsize, p_lastYsize, GDT_Float32, 0, 0 );
			  raster_pixles_n =RasterXsize*p_lastYsize;
			  if(np%4==1 && int(cp/4) > (pz-2)){
			     for(;raster_n <raster_pixles_n ;raster_n++)
			        {
				for(int k=0;k<5;k++){
				pw[k] = Compute_Weight(raster_n,Sample_N[k],point_n[k], adfGeoTransform[1],RasterXsize,pafScanline);
				
				}


		               float max_weight=0;
		               int sample_x = 10;
			         for(int i = 0;i<5;i++){
			            if(pw[i]>max_weight)
				      {
					  max_weight = pw[i];
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

		                }
					
			
		}
	
		
            if(np%4==2 && int(cp/4) > (pz-2)){
           for(;raster_n <raster_pixles_n;raster_n++)
			{

				int Nxp= 0,p = 1,px_i = 0, Sxp=0;

				loop1:
			 	   px_i = Nxp +cp-4;
				   if(px_i<=4){
				   //pixles_weight[px_i].a=px_i;
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);	
				   if(cp!=4){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,4,1,MPI_COMM_WORLD);
				   }
			       //tag2=tag*(np-1);
				   Nxp = 5-p*(np-4);
				            
				       if(Nxp>=np-4){
				        Nxp=(np-4)*p;
				         p++;
						 goto loop1;
	                   }
				       else if(Nxp<np-4 && Nxp>0){
				   
				       Nxp=(np-4)*p;
				       p++;
					   goto loop1;
				        }
					   
				   }	
					

				MPI_Status status;
			
			        if (cp==4){
				
			            tts1=MPI_Wtime();	
			        

				    MPI_Recv(&pixles_weight[1],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
				    MPI_Recv(&pixles_weight[3],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
				    tts2=MPI_Wtime();
			            ttsz=ttsz+tts2-tts1;
				    float max_weight=0;
		                    int sample_x = 10;
				    for(int i = 0;i<5;i++){
			              if(pixles_weight[i]>max_weight)
				      {
					  max_weight = pixles_weight[i];
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
				

		                }
			
			}

		     }

		     if(cp==4){
			}  

			}
		
			if(np%4==3 && int(cp/4) > (pz-2)){
			
			
			for(;raster_n <raster_pixles_n;raster_n++)
			{
			  int Nxp= 0,p = 1,px_i = 0, Sxp=0;
			       

				loop2:
			 	   px_i = Nxp +cp-4;
				   if(px_i<=4){
				 
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				   
				   if(cp!=4){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,4,1,MPI_COMM_WORLD);
				   }
				   Nxp = 5-p*(np-4);
				            
				       if(Nxp>=np-4){
				        Nxp=(np-4)*p;
				         p++;
						 goto loop2;
	                   }
				       else if(Nxp<np-4 && Nxp>0){
				   
				       Nxp=(np-4)*p;
				       p++;
					   goto loop2;
				        }
					   
				   }	
			





			MPI_Status status;
			//MPI_Barrier(MPI_COMM_WORLD); 
			if (cp==np-3){
				
				   
				     tts1=MPI_Wtime();	
			        
					MPI_Recv(&pixles_weight[1],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
			        MPI_Recv(&pixles_weight[2],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
					 MPI_Recv(&pixles_weight[4],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
					  tts2=MPI_Wtime();	
					  ttsz=tts2-tts1+ttsz;
						float max_weight=0;
		               int sample_x = 10;
					   	for(int i = 0;i<5;i++){
			            if(pixles_weight[i]>max_weight)
				      {
					  max_weight = pixles_weight[i];
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
				

		                }
			
			}

			}

		     if(cp==np-3){
			//poDBand =poDstDS->GetRasterBand(Bands);
			//poDBand->RasterIO( GF_Write,0, y_size, RasterXsize, p_lastYsize, pafScanline, RasterXsize, p_lastYsize, GDT_Float32, 0, 0 );
			} 
		
			}

			    if( np==8 && int(cp/4)>(pz-2))
		{

			for(;raster_n < raster_pixles_n ;raster_n++)
			{
				
				int Nxp= 0,p = 1,px_i = 0, Sxp=0;

				loop3:
			 	   px_i = Nxp +cp-4;
				   if(px_i<=4){
				   //pixles_weight[px_i].a=px_i;
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				// cout<<"send:"<<raster_n<<" "<<pixles_weight[px_i]<<endl;
				   if(cp!=4){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,4,px_i,MPI_COMM_WORLD);
				   }
			       //tag2=tag*(np-1);
				   Nxp = 5-p*4;
				            
				       if(Nxp>=4){
				        Nxp=4*p;
				         p++;
						 goto loop3;
	                   }
				       else if(Nxp<4 && Nxp>0){
				   
				       Nxp=4*p;
				       p++;
					   goto loop3;
				        }
					   
				   }	



			MPI_Status status;
			
			if (cp==4 ){
				  //  cout<<raster_n<<endl;


				tts1=MPI_Wtime();    
				for(int j=1;j<4;j++){              
					MPI_Recv(&pw[j],1,MPI_FLOAT,4+j,j,MPI_COMM_WORLD,&status);
						    }


                       tts2=MPI_Wtime();		
		       ttsz=ttsz+tts2-tts1;
		   
		       float max_weight=0;
		       int sample_x = 10;
		       for(int i = 0;i<5;i++){
			            if(pw[i]>max_weight)
				      {
					  max_weight = pw[i];
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
				

		                }
					

			}

		}

					if(cp%4==0){

			
			}  







			}
			
		

			}
			}

		}

		if (np<=4){		

raster_pixles_n =RasterXsize*p_ysize;
		   
 for(Bands=3; Bands>1; Bands--)
	{
		
		poBand = poDataset->GetRasterBand(Bands);
         poBand->RasterIO( GF_Read, 0, y_size, RasterXsize, p_ysize,  pafScanline, RasterXsize, p_ysize, GDT_Float32, 0, 0 );
			
			int raster_n=0;
			
			if(np==1){
			  for(;raster_n <raster_pixles_n ;raster_n++)
			{
				
		
				for(int k=0;k<5;k++){
				pw[k] = Compute_Weight(raster_n,Sample_N[k],point_n[k], adfGeoTransform[1],RasterXsize,pafScanline);
				
				}


		               float max_weight=0;
		               int sample_x = 10;
					   	for(int i = 0;i<5;i++){
			            if(pw[i]>max_weight)
				      {
					  max_weight = pw[i];
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

		                }
					
			
		}
	
			
			}
	
            if(np==2){
			          for(;raster_n <raster_pixles_n;raster_n++)
			{

					int Nxp= 0,p = 1,px_i = 0, Sxp=0;

				loop41:
			 	   px_i = Nxp +cp;
				   if(px_i<=4){
				   
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				
				   if(cp!=0){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
				   }
			      
				   Nxp = 5-p*np;
				            
				       if(Nxp>=np){
				        Nxp=np*p;
				         p++;
						 goto loop41;
	                   }
				       else if(Nxp<np && Nxp>0){
				   
				       Nxp=np*p;
				       p++;
					   goto loop41;
				        }
					   
				   }	
					
			

				MPI_Status status;
			
			if (cp==0){
				
				    ts1=MPI_Wtime();	
			        

				MPI_Recv(&pixles_weight[1],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
				MPI_Recv(&pixles_weight[3],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
					
					
				
				     ts2=MPI_Wtime();
					tsz=tsz+ts2-ts1;
					


						float max_weight=0;
		               int sample_x = 10;
					   	for(int i = 0;i<5;i++){
			            if(pixles_weight[i]>max_weight)
				      {
					  max_weight = pixles_weight[i];
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
				

		                }
			
			}

		     }

		     if(cp==0){
	
			} 
			
			}

	
            if(np==3){
				for(;raster_n <raster_pixles_n;raster_n++)
			{
			  int Nxp= 0,p = 1,px_i = 0, Sxp=0;
			       

				loop42:
			 	   px_i = Nxp +cp;
				   if(px_i<=4){
				 
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				   
				   if(cp!=0){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,0,1,MPI_COMM_WORLD);
				   }
				   Nxp = 5-p*np;
				            
				       if(Nxp>=np){
				        Nxp=np*p;
				         p++;
						 goto loop42;
	                   }
				       else if(Nxp<np && Nxp>0){
				   
				       Nxp=np*p;
				       p++;
					   goto loop42;
				        }
					   
				   }	
			





			MPI_Status status;
			//MPI_Barrier(MPI_COMM_WORLD); 
			if (cp==0){
				
				   
				     ts1=MPI_Wtime();	
			        
					MPI_Recv(&pixles_weight[1],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
			        MPI_Recv(&pixles_weight[2],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
					 MPI_Recv(&pixles_weight[4],1,MPI_FLOAT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
					  ts2=MPI_Wtime();	
					  tsz=ts2-ts1+tsz;
						float max_weight=0;
		               int sample_x = 10;
					   	for(int i = 0;i<5;i++){
			            if(pixles_weight[i]>max_weight)
				      {
					  max_weight = pixles_weight[i];
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
				

		                }
			
			}

			}

		     if(cp==0){
			//poDBand =poDstDS->GetRasterBand(Bands);
			//poDBand->RasterIO( GF_Write,0, y_size, RasterXsize, p_lastYsize, pafScanline, RasterXsize, p_lastYsize, GDT_Float32, 0, 0 );
			} 
			
			}	

            if(np==4){

			
			for(;raster_n < raster_pixles_n ;raster_n++)
			{
				
				int Nxp= 0,p = 1,px_i = 0, Sxp=0;

				loop43:
			 	   px_i = Nxp +cp;
				   if(px_i<=4){
				   //pixles_weight[px_i].a=px_i;
				   pixles_weight[px_i]= Compute_Weight(raster_n,Sample_N[px_i],point_n[px_i], adfGeoTransform[1],RasterXsize,pafScanline);
				// cout<<"send:"<<raster_n<<" "<<pixles_weight[px_i]<<endl;
				   if(cp!=0){
				   MPI_Send(&pixles_weight[px_i],1,MPI_FLOAT,0,px_i,MPI_COMM_WORLD);
				   }
			       //tag2=tag*(np-1);
				   Nxp = 5-p*4;
				            
				       if(Nxp>=4){
				        Nxp=4*p;
				         p++;
						 goto loop43;
	                   }
				       else if(Nxp<4 && Nxp>0){
				   
				       Nxp=4*p;
				       p++;
					   goto loop43;
				        }
					   
				   }	



			MPI_Status status;
			
			if (cp==0 ){
				

				ts1=MPI_Wtime();    
				for(int j=1;j<4;j++){

                    
					MPI_Recv(&pw[j],1,MPI_FLOAT,j,j,MPI_COMM_WORLD,&status);
					
					
						}


            ts2=MPI_Wtime();		
		   tsz=tsz+ts2-ts1;
		   //cout<<raster_n<<"  "<<ts2-ts1<<"  "<<tsz<<endl;
	
						float max_weight=0;
		               int sample_x = 10;
					   	for(int i = 0;i<5;i++){
			            if(pw[i]>max_weight)
				      {
					  max_weight = pw[i];
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
				

		                }
					

			}

		}

					if(cp==0){

			
			}  







			}	









}

}
		MPI_Barrier(MPI_COMM_WORLD); 


	
t2 = MPI_Wtime();
	
	if(cp==0){
	printf("cosume time: %f",t2-t1);
	printf("send time = %f",tsz);
	}
	if(np==6 && cp==4){printf("send time(6) = %f",ttsz);}
	if(np==7 && cp==4){printf("send time(7) = %f",ttsz);}
	if(np==8 && cp==4){printf("send time(8) = %f",ttsz);}
	if(cp%4==0){
	if( poDstDS != NULL ){
	 GDALClose( (GDALDatasetH) poDstDS );
	}
	}
	MPI_Finalize();
}
