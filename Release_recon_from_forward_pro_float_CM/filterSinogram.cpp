#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "filterSinogram.h"


double**
prepareConstMap(int projWidth, double DETECTOR_PITCH_CM_AT_ISO_CENT)
{
	double**	constMap = (double**)new2DArray( projWidth, projWidth, UNIT_FLOAT64 );

	if( constMap == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for prepareConstMap\n");
		return NULL;
	}

	double theConst = 1.0/(8.0*PI*PI*DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM / projWidth);
	double gamma_x_cutoffFreq = GAMMA * CUTOFF_FREQ;
	double PI_x_cutoffFreq = PI * CUTOFF_FREQ;
	double one_minus_alpha = 1.0 - ALPHA;
	double omega = PI/GAMMA;
	double cos_omega = cos(omega);
	double sin_omega = sin(omega);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int m=0;m<projWidth;m++)
	{
		for(int n=0;n<projWidth;n++)
		{
			if( m == n )
			{
				constMap[m][n] = 0.0;
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

				
				constMap[m][n] = theConst /MN/MN *
							( 2.0*ALPHA*theta*sin_theta
					 			+ 2.0*ALPHA*(cos_theta-1.0)
								+ xi*one_minus_alpha*theta*(sin_theta*cos_omega+cos_theta*sin_omega)
					 			+ eta*one_minus_alpha*theta*(sin_theta*cos_omega-cos_theta*sin_omega)
								+ xi*xi*one_minus_alpha*(cos_theta*cos_omega-sin_theta*sin_omega -1.0)
								+ eta*eta*one_minus_alpha*(cos_theta*cos_omega+sin_theta*sin_omega -1.0)
								);
			}
		}
	}

	return constMap;
}

void
filterSinogram( int projHeight, int projWidth, int sngrmHeight,
			float *** projVolume, float*** fltdVoxData,
			double** constMap, double* xshift, double* yshift, double DETECTOR_PITCH_CM_AT_ISO_CENT )
{
	if( fltdVoxData == NULL )
	{
		fprintf(stderr,"ERROR: fltdVoxData is NULL\n");
		return;
	}

	double omega = PI/GAMMA;
	double coeff = CUTOFF_FREQ*CUTOFF_FREQ/4.0/
	  (DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM / projWidth) *(ALPHA*0.5+(1.0-ALPHA)*(1.0/omega*sin(omega)+1.0/omega/omega*(cos(omega)-1.0)));
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int th=0;th<sngrmHeight;th++)
	{
	  double	XcenterShift = -projWidth/2.0 - xshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth) -1;
	  double	YcenterShift = -projWidth/2.0 - yshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM / projWidth) -1;
	  //fprintf(stderr,"X = %lf, Y = %lf \n", XcenterShift, YcenterShift);
	  for(int n=0;n<projHeight;n++)
	  //for(int n=256;n<257;n++)
	    {
	      double pos2 = fabs((YcenterShift + n)) * (DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM / projWidth);

		for(int m=0;m<projWidth;m++)
		{
		  double pos1 = fabs((XcenterShift + m)) * (DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM / projWidth);
		  double Cm = DIST_BTWN_SRC_AND_ISOCENT_CM /
					sqrt(	DIST_BTWN_SRC_AND_ISOCENT_CM*DIST_BTWN_SRC_AND_ISOCENT_CM +
						pos1 * pos1 + pos2 * pos2 );

		  
	  //########### if you use original  projection(bone has low value),and I0=1;###########
		  fltdVoxData[th][n][m] = -log( projVolume[th][n][m]) * Cm * coeff;//(poly energy)//MAX_AIR_PIXDIFF=1
		  
		  
		  for(int k = 0;k < projWidth;k++)
		  	{
		  	  fltdVoxData[th][n][m] += constMap[m][k] * (-log( projVolume[th][n][k])) * Cm;//(poly energy) //MAX_AIR_PIXDIFF=1

		  	}

       
		  
		}
		
		//fprintf(stderr,"th = %d, n = %d \n", th, n);
	    }
	}

	/*char	filename[256];
		sprintf(filename,"fltsinogram%dx%d.float32.raw",projWidth,sngrmHeight);
		FILE*	fp=fopen(filename,"wb");
		for(int th=0;th<sngrmHeight;th++)
		  for(int m=0;m<projWidth;m++)
		    {
		      float flt=fltdVoxData[th][256][m];
		fwrite( &flt, getByteSizeOfUnit(UNIT_FLOAT32),1,fp);
		    }
		    fclose(fp);*/
	
}


