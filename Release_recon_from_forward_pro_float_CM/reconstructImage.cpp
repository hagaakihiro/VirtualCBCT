#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "reconstructImage.h"


void
reconstructImage( int reconSize, int reconSlices, int zstart, double reconStep, double reconScale, double*** reconVoxData, char** maskData,
		  int projHeight, int projWidth, int sngrmHeight, float*** sngrmVoxData, int xfactor,
		  double* anglesRad, double* cosTable, double* sinTable,
		  double* xshift, double DETECTOR_PITCH_CM_AT_ISO_CENT )
{
  double scale = reconScale;
  double centerPos = (double)reconSize/2.0 - 0.5;
  
  if( reconVoxData == NULL )
    {
      fprintf(stderr,"ERROR: reconVoxData is NULL\n");
      return;
    }
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int z=0;z<reconSlices;z++)
    {
      double zz = zstart +  z*reconStep;
      //printf("zz = %lf %d\n", zz, xfactor);
      
      for(int y=0;y<reconSize;y++)
	{	
	  for(int x=0;x<reconSize;x++)
	    {	
	      if( maskData[y][x] != 0 )
		{
		  double xx = (x - centerPos) * scale;
		  double yy = - (y - centerPos) * scale;
		  
		  for(int th=0;th<sngrmHeight;th++)
		    {
		      double dbeta;
		      double thdb0, thdb1, thdb2;
		      
		      if( th==0 )
			{
			  thdb1 = anglesRad[th];
			  thdb2 = anglesRad[th+1];
			  
			  dbeta = fabs(thdb2-thdb1)/2;						
			}
		      else
			if( th == sngrmHeight-1 )
			  {
			    thdb0 = anglesRad[th-1];
			    thdb1 = anglesRad[th];
			    
			    dbeta = fabs(thdb1-thdb0)/2;						
			  }
			else
			  {
			    thdb0 = anglesRad[th-1];
			    thdb2 = anglesRad[th+1];
			    thdb1 = anglesRad[th+1]-2*PI;
			    if(thdb2 < 0) thdb1 = anglesRad[th+1]+2*PI;	
			    
			    dbeta = fabs(thdb2-thdb0)/2;
			    if (dbeta > PI/2) dbeta = fabs(thdb1-thdb0)/2;
			  }
		      
		      double	cos_thdb1 = cosTable[th];
		      double	sin_thdb1 = sinTable[th];
		      
		      double X =  xx*cos_thdb1 + yy*sin_thdb1;
		      double Y = -xx*sin_thdb1 + yy*cos_thdb1;
		      double L2 =  1.0 - Y/DIST_BTWN_SRC_AND_ISOCENT_CM;
		      
		      double U = X/L2/(DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM/ projWidth);	
		      double intpWeight = fabs( fabs(U) - fabs(floor(U)) );
		      
		      double  UUU = (double)(projWidth/2 + U + xshift[th]/(DETECTOR_PITCH_CM * MAX_PIC_NUM/ projWidth) - 0.5 );
		      int	UU = (int)floor( UUU );
		      
		      double  ZETA = (double)(projHeight/2) + zz/L2/(DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM/ projWidth);
		      //double  ZETA = (double)(projHeight/2 + zz/L2/(DETECTOR_PITCH_CM_AT_ISO_CENT * MAX_PIC_NUM/ projWidth));	//shimo
		      int	NZETA = (int)floor( ZETA );
		      double zintpWeight = fabs( fabs(ZETA) - fabs(floor(ZETA)) );
		      /*
		      if(z==reconSlices-1)
			{
			  printf("zz = %lf ZETA=%lf %d UUU=%lf %d\n", zz, ZETA,NZETA,UUU,UU);
			}
		      */
		      /*
		      if(z==0)
			{
			  printf("zz = %lf ZETA=%lf %d UUU=%lf %d\n", zz, ZETA,NZETA,UUU,UU);
			}
		      */
		      /*
		      //pre program
		      if( (UU >= 0 && UU < projWidth) &&  (NZETA >= 0 && NZETA < projHeight) )
			{
			  double sum =
			    sngrmVoxData[th][NZETA][UU] * (1.0-intpWeight) * (1.0-zintpWeight) +
			    sngrmVoxData[th][NZETA][UU+1] * intpWeight* (1.0-zintpWeight)  +
			    sngrmVoxData[th][NZETA+1][UU] * (1.0-intpWeight)*zintpWeight +
			    sngrmVoxData[th][NZETA+1][UU+1] * intpWeight*zintpWeight;
			  
			  reconVoxData[z][y][x] += ( 1.0/L2/L2 * sum * dbeta * xfactor);
			  //printf("sum = %lf %lf\n", sum, sngrmVoxData[th][NZETA][UU]);
			  //reconVoxData2[z][y][x] = reconVoxData[z][y][x];
			  
			}
		      */
		      //2022/2/16 shimomura program
		      if( (UU >= 0 && UU < projWidth-1) &&  (NZETA >= 0 && NZETA < projHeight-1) )
			{
			  double sum =
			    sngrmVoxData[th][NZETA][UU] * (1.0-intpWeight) * (1.0-zintpWeight) +
			    sngrmVoxData[th][NZETA][UU+1] * intpWeight* (1.0-zintpWeight)  +
			    sngrmVoxData[th][NZETA+1][UU] * (1.0-intpWeight)*zintpWeight +
			    sngrmVoxData[th][NZETA+1][UU+1] * intpWeight*zintpWeight;
			  
			  reconVoxData[z][y][x] += ( 1.0/L2/L2 * sum * dbeta * xfactor);
			  //printf("sum = %lf %lf\n", sum, sngrmVoxData[th][NZETA][UU]);
			  //reconVoxData2[z][y][x] = reconVoxData[z][y][x];
			  
			}
		      if( ( UU == projWidth-1) &&  (NZETA >= 0 && NZETA < projHeight-1) )
			{
			  //double sum =
			  //sngrmVoxData[th][NZETA][UU] * (1.0-intpWeight) * (1.0-zintpWeight) +
			  //sngrmVoxData[th][NZETA+1][UU] * (1.0-intpWeight)*zintpWeight;

			  double sum =
			    sngrmVoxData[th][NZETA][UU]  * (1.0-zintpWeight) +
			    sngrmVoxData[th][NZETA+1][UU] *zintpWeight;
  
			  
			  reconVoxData[z][y][x] += ( 1.0/L2/L2 * sum * dbeta * xfactor);
			  //printf("sum = %lf %lf\n", sum, sngrmVoxData[th][NZETA][UU]);
			  //reconVoxData2[z][y][x] = reconVoxData[z][y][x];
			  
			}
		      if( (UU >= 0 && UU < projWidth-1) &&  ( NZETA == projHeight-1) )
			{
			  //double sum =
			  //sngrmVoxData[th][NZETA][UU] * (1.0-intpWeight) * (1.0-zintpWeight) +
			  //sngrmVoxData[th][NZETA][UU+1] * intpWeight* (1.0-zintpWeight);//dame
			  double sum =
			  sngrmVoxData[th][NZETA][UU] * (1.0-intpWeight)  +
			  sngrmVoxData[th][NZETA][UU+1] * intpWeight;
			  
			  reconVoxData[z][y][x] += ( 1.0/L2/L2 * sum * dbeta * xfactor);
			  //printf("sum = %lf %lf\n", sum, sngrmVoxData[th][NZETA][UU]);
			  //reconVoxData2[z][y][x] = reconVoxData[z][y][x];
			  
			}

		      if( ( UU == projWidth-1) &&  ( NZETA == projHeight-1) )
			{
			  //double sum =
			  //sngrmVoxData[th][NZETA][UU]* (1.0-intpWeight) * (1.0-zintpWeight);
			
			  double sum =
			  sngrmVoxData[th][NZETA][UU];
			  
			  reconVoxData[z][y][x] += ( 1.0/L2/L2 * sum * dbeta * xfactor);
			  //printf("sum = %lf %lf\n", sum, sngrmVoxData[th][NZETA][UU]);
			  //reconVoxData2[z][y][x] = reconVoxData[z][y][x];
			  
			}
		      
		      
		      
		    }				
		}
	    }
	}
    }
}
