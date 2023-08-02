#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "reconstructImageOnGPU.h"
#include "mallocMD.h"


// kernel code

__constant__ double	devAnglesRad[MAX_Proj_Num];
__constant__ double	devCosTable[MAX_Proj_Num];
__constant__ double	devSinTable[MAX_Proj_Num];
__constant__ double	devXshift[MAX_Proj_Num];
__constant__ double	devYshift[MAX_Proj_Num];


__global__ void
backProjKernel( int z, int reconSize, int zstart, int reconStep, char* devMaskData, char* devMaskAngle, int xfactor,
			int projHeight, int projWidth, int sngrmHeight, float* devSngrmVoxData,
			double centerPos, double scale, double* devReconVoxData, double* devReconVoxData2, double DETECTOR_PITCH_MM_AT_ISO_CENT)
{
	int	x = blockIdx.x;
	int	y = blockIdx.y;
	int	th = blockIdx.z;

	int	horPos = z*reconSize*reconSize + y*reconSize + x;
	int	mas = y*reconSize + x;

	if( devMaskData[mas] != 0 )
	{
		double xx = (x - centerPos) * scale;
		double yy = -(y - centerPos) * scale;
		double zz = zstart +  z*reconStep;
		// Change for high resolution CBCT for middle area using M panel.
		//int maskangle = 1; 
		int maskangle = 0;
		if(devMaskAngle[th*reconSize*reconSize+y*reconSize+x] != 0) maskangle = 1;
		{
			double dbeta;
			double thdb0, thdb1, thdb2;
			
			if( th==0 )
			{
				thdb1 = devAnglesRad[th];
				thdb2 = devAnglesRad[th+1];
					
				dbeta = fabs(thdb2-thdb1)/2;	
			}
			else
			if( th == sngrmHeight-1 )
			{
				thdb0 = devAnglesRad[th-1];
				thdb1 = devAnglesRad[th];
					
				dbeta = fabs(thdb1-thdb0)/2;						
			}
			else
			{
			thdb0 = devAnglesRad[th-1];
			thdb2 = devAnglesRad[th+1];
			thdb1 = devAnglesRad[th+1]-2*PI;
			if(thdb2 < 0) thdb1 = devAnglesRad[th+1]+2*PI;
				
			dbeta = fabs(thdb2-thdb0)/2;
			if (dbeta > PI/2) dbeta = fabs(thdb1-thdb0)/2;
			}

			double	cos_thdb1 = devCosTable[th];
			double	sin_thdb1 = devSinTable[th];

			double X =  xx*cos_thdb1 + yy*sin_thdb1;
			double Y = -xx*sin_thdb1 + yy*cos_thdb1;
			double L2 =  1.0 - Y/DIST_BTWN_SRC_AND_ISOCENT_MM;
					
			double	U = X/L2/(DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM/ projWidth);
			double	intpWeight = fabs( fabs(U) - fabs(floor(U)) );

			double  UUU = (double)(projWidth/2 + U +devXshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM/ projWidth) - 0.5 );
			int	UU = (int)floor( UUU );

			double  ZETA = (double)(projHeight/2) + zz/L2/(DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM/ projWidth);	
			int	NZETA = (int)floor( ZETA );

			double	zintpWeight = fabs( fabs(ZETA) - fabs(floor(ZETA)) );

						
			if( (UU >= 0 && UU < projWidth))// &&  (NZETA >= 0 && NZETA < projHeight) )
			{
				double sum =
			                     devSngrmVoxData[th*projHeight*projWidth+NZETA*projWidth+UU]* (1.0-intpWeight) * (1.0-zintpWeight) +
					     devSngrmVoxData[th*projHeight*projWidth+NZETA*projWidth+UU+1] * intpWeight* (1.0-zintpWeight)  +
					     devSngrmVoxData[th*projHeight*projWidth+(NZETA+1)*projWidth+UU] * (1.0-intpWeight)*zintpWeight +
					     devSngrmVoxData[th*projHeight*projWidth+(NZETA+1)*projWidth+UU+1] * intpWeight*zintpWeight;

				devReconVoxData[horPos] += ( 1.0/L2/L2 * sum * dbeta * xfactor * maskangle);
				devReconVoxData2[horPos] = devReconVoxData[horPos];

			}

		}				
	}
	__syncthreads();
}

__global__ void
filterOptimizedMean( int reconSize, double* devReconVoxData, double*
		devReconVoxData2, char* devMaskData, int sngrmHeight)
{
	int	i = blockIdx.x;
	int	j = blockIdx.y;
	int	k = blockIdx.z;
	int     Nfilt = 3;
	int     filt = 50;
	int     count=0;
		  for(int L2 = 0; L2 < Nfilt; L2++)
		    {
		      for(int L3 = 0; L3 < Nfilt; L3++)
			if(devReconVoxData2[k*reconSize*reconSize+i*reconSize+j] < 
			   devReconVoxData2[k*reconSize*reconSize+(i-1+L2)*reconSize+j-1+L3] + filt  ||
			   devReconVoxData2[k*reconSize*reconSize+i*reconSize+j] > 
			   devReconVoxData2[k*reconSize*reconSize+(i-1+L2)*reconSize+j-1+L3] - filt) 
			  {
			    if(L2 * L3 != 1)
			    {
  			    devReconVoxData[k*reconSize*reconSize+i*reconSize+j] += devReconVoxData2[k*reconSize*reconSize+(i-1+L2)*reconSize+j-1+L3];
			    }
			    count++;
			  }
		    }
	devReconVoxData[k*reconSize*reconSize+i*reconSize+j]=devReconVoxData[k*reconSize*reconSize+i*reconSize+j]/(count-1.0);
	__syncthreads();
}

__global__ void
prepareConstMapGPU( int projHeight, int projWidth, double* devconstMap, double DETECTOR_PITCH_MM_AT_ISO_CENT)
{
	int	m = blockIdx.x;
	int	n = blockIdx.y;
	int     twod = m*projWidth + n;

	double theConst = 1.0/(8.0*PI*PI*DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM / projWidth);
	double gamma_x_cutoffFreq = GAMMA * CUTOFF_FREQ;
	double PI_x_cutoffFreq = PI * CUTOFF_FREQ;
	double one_minus_alpha = 1.0 - ALPHA;
	double omega = PI/GAMMA;
	double cos_omega = cos(omega);
	double sin_omega = sin(omega);

			if( m == n )
			{
				devconstMap[twod] = 0.0;
			}

			else
			{
				double MN = (double)(m-n);

				double theta = PI_x_cutoffFreq * MN;
				double delta = gamma_x_cutoffFreq * MN;
				double eta = 1.0/(1.0-1.0/delta);                  
				double xi = 1.0/(1.0+1.0/delta);
					
				double cos_theta = cos(theta);
				double sin_theta = sin(theta);

				devconstMap[twod] = theConst /MN/MN *
							( 2.0*ALPHA*theta*sin_theta
					 			+ 2.0*ALPHA*(cos_theta-1.0)
								+ xi*one_minus_alpha*theta*(sin_theta*cos_omega+cos_theta*sin_omega)
					 			+ eta*one_minus_alpha*theta*(sin_theta*cos_omega-cos_theta*sin_omega)
								+ xi*xi*one_minus_alpha*(cos_theta*cos_omega-sin_theta*sin_omega -1.0)
								+ eta*eta*one_minus_alpha*(cos_theta*cos_omega+sin_theta*sin_omega -1.0)
								);
			}
	__syncthreads();
}

__global__ void
makesinogram( int zstartpic, int zendpic, int zpitch, int projHeight, int projWidth, int sngrmHeight,
unsigned short* devprojVolume, float* devSngrmVoxData1, double DETECTOR_PITCH_MM_AT_ISO_CENT)
{
	int	th = blockIdx.x; //sngrmHeight
	int	n = blockIdx.y; //pitch
	int     m = blockIdx.z; //projWidth

	double	XcenterShift = -projWidth/2.0 - devXshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth);
	double	YcenterShift = -projHeight/2.0 + devYshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth);
	double  pos2 = fabs((YcenterShift + (n+zstartpic))) * (DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM /projWidth);
        double  pos1 = fabs((XcenterShift + m)) * (DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM /projWidth);
	double  Cm = DIST_BTWN_SRC_AND_ISOCENT_MM /sqrt(DIST_BTWN_SRC_AND_ISOCENT_MM * DIST_BTWN_SRC_AND_ISOCENT_MM + pos1 * pos1 + pos2 * pos2 );

	double  nw = (n+zstartpic + devYshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth));
	double	intpWeight = fabs( fabs(nw) - fabs(floor(nw)) );

	int     threed = th * projWidth * projHeight + (n+zstartpic) * projWidth + m;
	int     threed1 = th * projWidth * projHeight + (int)floor(nw) * projWidth + m;
	int     threed2 = th * projWidth * projHeight + (int)floor(nw) * projWidth + m-1;
	int     threed3 = th * projWidth * projHeight + (int)floor(nw+1) * projWidth + m;
	int     threed4 = th * projWidth * projHeight + (int)floor(nw+1) * projWidth + m-1;

	devSngrmVoxData1[threed] = (float)(devprojVolume[threed1] * (1-intpWeight)*0.5
				   	   +devprojVolume[threed2] * (1-intpWeight)*0.5
					   +devprojVolume[threed3] * intpWeight*0.5
                                           +devprojVolume[threed4] * intpWeight*0.5 );

	if(m == 0) devSngrmVoxData1[threed] = (float)(devprojVolume[threed1] * (1-intpWeight)
					             +devprojVolume[threed3] * intpWeight );

	devSngrmVoxData1[threed] = (float)(-log((MAX_PIXVALUE - devSngrmVoxData1[threed])/MAX_AIR_PIXDIFF) * Cm);
	if(devSngrmVoxData1[threed] < 0 ) devSngrmVoxData1[threed] = 0;

	__syncthreads();
}

__global__ void
makefiltersinogram( int zstartpic, int zendpic, int zpitch, int projHeight, int projWidth, int sngrmHeight,
float* devSngrmVoxData1, float* devSngrmVoxData, double* devconstMap, double coeff, double DETECTOR_PITCH_MM_AT_ISO_CENT)
{
	int	th = blockIdx.x; //sngrmHeight
	int	n = blockIdx.y; //pitch
	int     m = blockIdx.z; //projWidth
	int     threed = th * projWidth * projHeight + (n+zstartpic) * projWidth + m;

	double  XcenterShift = -projWidth/2.0 - devXshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth);
	double  YcenterShift = -projHeight/2.0 + devYshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth);
	double  pos2 = fabs((YcenterShift + (n+zstartpic))) * (DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM /projWidth);
        double  pos1 = fabs((XcenterShift + m)) * (DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM /projWidth);
	double  Cm = DIST_BTWN_SRC_AND_ISOCENT_MM /sqrt(DIST_BTWN_SRC_AND_ISOCENT_MM * DIST_BTWN_SRC_AND_ISOCENT_MM + pos1 * pos1 + pos2 * pos2 );

	devSngrmVoxData[threed] = coeff * devSngrmVoxData1[threed];
	double summ = 0.0;
	for(int k=0;k<projWidth;k++)
	{
		int	threed2 = th*projWidth*projHeight + (n+zstartpic)*projWidth + k;
		int	twodd = m * projWidth + k;
		summ += devconstMap[twodd] * devSngrmVoxData1[threed2];
	}
	devSngrmVoxData[threed] += float(summ);
	__syncthreads();
}

// hostCode

double*	devReconVoxData;
double*	devReconVoxData2;
char*	devMaskData;
char*	devMaskAngle;
float*	devSngrmVoxData;
float*	devSngrmVoxData1;
unsigned short* devprojVolume;
double* devconstMap;

void
initializeGPU(int reconSize, int zslices, int projHeight, int projWidth, int sngrmHeight, double* anglesRad, double* cosTable, double* sinTable, double* xshift, double* yshift)
{
	cudaSetDevice(1);

	cudaMalloc( (void**)&devReconVoxData, sizeof(double)*reconSize*reconSize*zslices );
	cudaMalloc( (void**)&devReconVoxData2, sizeof(double)*reconSize*reconSize*zslices );
	cudaMalloc( (void**)&devMaskData, sizeof(char)*reconSize*reconSize );
	cudaMalloc( (void**)&devMaskAngle, sizeof(char)*reconSize*reconSize*sngrmHeight );
	cudaMalloc( (void**)&devSngrmVoxData, sizeof(float)*projHeight*projWidth*sngrmHeight );
	cudaMalloc( (void**)&devSngrmVoxData1, sizeof(float)*projHeight*projWidth*sngrmHeight );
	cudaMalloc( (void**)&devprojVolume, sizeof(unsigned short)*projHeight*projWidth*sngrmHeight );
	cudaMalloc( (void**)&devconstMap, sizeof(double)*projHeight*projWidth );


	cudaMemcpyToSymbol( devAnglesRad, anglesRad, sizeof(double)*sngrmHeight );
	cudaMemcpyToSymbol( devCosTable, cosTable, sizeof(double)*sngrmHeight );
	cudaMemcpyToSymbol( devSinTable, sinTable, sizeof(double)*sngrmHeight );
	cudaMemcpyToSymbol( devXshift, xshift, sizeof(double)*sngrmHeight );
	cudaMemcpyToSymbol( devYshift, yshift, sizeof(double)*sngrmHeight );
}

void
terminateGPU()
{
	cudaFree(devReconVoxData);
	cudaFree(devReconVoxData2);
	cudaFree(devMaskData);
	cudaFree(devMaskAngle);
	cudaFree(devSngrmVoxData);
	cudaFree(devSngrmVoxData1);
	cudaFree(devprojVolume);
	cudaFree(devconstMap);

}


void
filterSinogramOnGPU( int zstart, int zslices, int projHeight, int projWidth, int sngrmHeight,
			unsigned short*** projVolume, float*** fltdVoxData,
			double* xshift, double* yshift, double DETECTOR_PITCH_MM_AT_ISO_CENT )
{
	if( projVolume == NULL )
	{
		fprintf(stderr,"ERROR: projVolume is NULL\n");
		return;
	}

	cudaMemcpy( devSngrmVoxData, (void*)&(fltdVoxData[0][0][0]), sizeof(float)*projHeight*projWidth*sngrmHeight, cudaMemcpyHostToDevice );
	cudaMemcpy( devSngrmVoxData1, (void*)&(fltdVoxData[0][0][0]), sizeof(float)*projHeight*projWidth*sngrmHeight, cudaMemcpyHostToDevice );
	cudaMemcpy( devprojVolume, (void*)&(projVolume[0][0][0]), sizeof(unsigned short)*projHeight*projWidth*sngrmHeight, cudaMemcpyHostToDevice );

	// ConstMap
	{
	   	dim3	DimGrid(projHeight,projWidth);
	   	dim3	DimBlk(1);
	   	prepareConstMapGPU<<<DimGrid, DimBlk>>>( projHeight, projWidth, devconstMap, DETECTOR_PITCH_MM_AT_ISO_CENT);
	   	cudaThreadSynchronize(); 
	}
	// sinogram
	{
		int zstartpic = (int)(zstart/(DETECTOR_PITCH_MM_AT_ISO_CENT*MAX_PIC_NUM/projWidth) + projHeight/2.0)-10;
		int zendpic = (int)((zstart+zslices)/(DETECTOR_PITCH_MM_AT_ISO_CENT*MAX_PIC_NUM/projWidth) + projHeight/2)+10;
		int zpitch = (int)(zendpic-zstartpic);
		dim3	DimGrid(sngrmHeight,zpitch,projWidth);
		//dim3	DimGrid(sngrmHeight,projHeight,projWidth);
		dim3	DimBlk(1);

		makesinogram<<<DimGrid, DimBlk>>>( 	zstartpic, zendpic, zpitch,
								projHeight,
								projWidth,
								sngrmHeight,
								devprojVolume,
								devSngrmVoxData1,
								DETECTOR_PITCH_MM_AT_ISO_CENT
								);
	}
	cudaMemcpy( (void*)&(fltdVoxData[0][0][0]), devSngrmVoxData1, sizeof(float)*projHeight*projWidth*sngrmHeight, cudaMemcpyDeviceToHost );
			/*char	filename[256];
			sprintf(filename,"fltsinogram%dx%d.float32.raw",projHeight,sngrmHeight);
			FILE*	fp2=fopen(filename,"wb");
			for(int th=0;th<sngrmHeight;th++)
			  for(int m=0;m<projWidth;m++)
			    {
			      float flt=fltdVoxData[th][projHeight/2][m];
			      fwrite( &flt, getByteSizeOfUnit(UNIT_FLOAT32),1,fp2);
			    }
			    fclose(fp2);*/
	// filtered sinogram
	{
		{
	  	//double	XcenterShift;// = -projWidth/2.0;// - xshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth) -1;
	  	//double	YcenterShift;// = -projWidth/2.0;// - yshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth) -1;
		int zstartpic = (int)(zstart/(DETECTOR_PITCH_MM_AT_ISO_CENT*MAX_PIC_NUM/projWidth) + projHeight/2.0)-10;
		int zendpic = (int)((zstart+zslices)/(DETECTOR_PITCH_MM_AT_ISO_CENT*MAX_PIC_NUM/projWidth) + projHeight/2)+10;
		int zpitch = (int)(zendpic-zstartpic);
		fprintf(stderr,"zstart = %d, zend = %d \n",zstartpic,zendpic);
		dim3	DimGrid(sngrmHeight,zpitch,projWidth);
		//dim3	DimGrid(sngrmHeight,projHeight,projWidth);
		dim3	DimBlk(1);
		double omega = PI/GAMMA;
	   	double coeff = CUTOFF_FREQ*CUTOFF_FREQ/4.0/
	  	       (DETECTOR_PITCH_MM_AT_ISO_CENT * MAX_PIC_NUM / projWidth) *(ALPHA*0.5+(1.0-ALPHA)*(1.0/omega*sin(omega)+1.0/omega/omega*(cos(omega)-1.0)));
		makefiltersinogram<<<DimGrid, DimBlk>>>( 	zstartpic, zendpic, zpitch,
								projHeight,
								projWidth,
								sngrmHeight,
								devSngrmVoxData1,
								devSngrmVoxData,
								devconstMap,
								coeff, DETECTOR_PITCH_MM_AT_ISO_CENT
								);
		}
	cudaThreadSynchronize();
	}
	cudaMemcpy( (void*)&(fltdVoxData[0][0][0]), devSngrmVoxData, sizeof(float)*projHeight*projWidth*sngrmHeight, cudaMemcpyDeviceToHost );
			/*char	filename2[256];
			sprintf(filename2,"fltsinogram%dx%d.float32_i.raw",projHeight,sngrmHeight);
			FILE*	fp3=fopen(filename2,"wb");
			for(int th=0;th<sngrmHeight;th++)
			  for(int m=0;m<projWidth;m++)
			    {
			      float flt=fltdVoxData[th][projHeight/2][m];
			      fwrite( &flt, getByteSizeOfUnit(UNIT_FLOAT32),1,fp3);
			    }
			    fclose(fp3);*/
	//fprintf(stderr,"fltdVoxData = %e \n",fltdVoxData[0][256][255]);
}


void
reconstructImageOnGPU( int reconSize, int zstart, int zslices, int reconStep, double*** reconVoxData, char** maskData,
		int xfactor, int projHeight, int projWidth, int sngrmHeight, float*** sngrmVoxData, char* maskAngle, double scale, double DETECTOR_PITCH_MM_AT_ISO_CENT )
{
	if( reconVoxData == NULL )
	{
		fprintf(stderr,"ERROR: reconVoxData is NULL\n");
		return;
	}

	cudaMemcpy( devReconVoxData, (void*)&(reconVoxData[0][0][0]), sizeof(double)*reconSize*reconSize*zslices, cudaMemcpyHostToDevice );
	cudaMemcpy( devReconVoxData2, (void*)&(reconVoxData[0][0][0]), sizeof(double)*reconSize*reconSize*zslices, cudaMemcpyHostToDevice );
	cudaMemcpy( devMaskData, (void*)&(maskData[0][0]), sizeof(char)*reconSize*reconSize, cudaMemcpyHostToDevice );
	cudaMemcpy( devMaskAngle, (void*)&(maskAngle[0]), sizeof(char)*reconSize*reconSize*sngrmHeight, cudaMemcpyHostToDevice );
	cudaMemcpy( devSngrmVoxData, (void*)&(sngrmVoxData[0][0][0]), sizeof(float)*projHeight*projWidth*sngrmHeight, cudaMemcpyHostToDevice );


	// back projection
	{
		dim3	DimGrid(reconSize,reconSize,sngrmHeight);
		dim3	DimBlk(1);
		double  centerPos = (double)reconSize/2.0 - 0.5;

		int z;
		for(z=0;z<zslices;z++)
		{
		backProjKernel<<<DimGrid, DimBlk>>>( z, reconSize, zstart, reconStep,
								devMaskData, devMaskAngle, xfactor,
								projHeight, projWidth, sngrmHeight, devSngrmVoxData,
								centerPos, scale, devReconVoxData, devReconVoxData2, DETECTOR_PITCH_MM_AT_ISO_CENT);
		}
		cudaThreadSynchronize();

	}

	// filter
	{
		dim3	DimGrid(reconSize, reconSize, zslices);
		dim3	DimBlk(1);
		filterOptimizedMean<<<DimGrid, DimBlk>>>( reconSize,
					  devReconVoxData, devReconVoxData2,
					  devMaskData, sngrmHeight);
		cudaThreadSynchronize();
	}
	cudaMemcpy( (void*)&(reconVoxData[0][0][0]), devReconVoxData, sizeof(double)*reconSize*reconSize*zslices, cudaMemcpyDeviceToHost );
	//fprintf(stderr,"reconVoxData = %lf %lf\n",reconVoxData[0][0][0],reconVoxData[0][135][135]);
}

























