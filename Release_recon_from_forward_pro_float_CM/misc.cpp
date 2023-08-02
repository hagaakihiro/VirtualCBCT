#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mallocMD.h"
#include "physParams.h"


char**
prepareMaskData(int reconSize)
{
	char**	maskPixData = (char**)new2DArray( reconSize, reconSize, UNIT_UINT8 );

	if( maskPixData == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for renderMaskData\n");
		return NULL;
	}

	double maskRadius = (double)(reconSize-1)/2.0;
	double	sqMaskRad = maskRadius*maskRadius;

	for(int j=0;j<reconSize;j++)
	for(int i=0;i<reconSize;i++)
	{
		if( ((double)j-maskRadius)*((double)j-maskRadius) +
				((double)i-maskRadius)*((double)i-maskRadius) <= sqMaskRad )
		{
			maskPixData[j][i] = 1;
		}
		else
		{
			maskPixData[j][i] = 0;
		}
	}
	return maskPixData;
}

void
//clearReconstructedImage(int reconSize, int sngrmHeight, double** reconPixData)
clearReconstructedImage(int reconSize, int zslices, int sngrmHeight, double*** reconVoxData, short*** reconVoxDataShort)
{
	int length = reconSize*reconSize*zslices;
	double* pixPtr = (double*)get1DArrayOf3DArray((void***)reconVoxData,UNIT_FLOAT64);
	int length2 = reconSize*reconSize*zslices;
       	short* pixPtr2 = (short*)get1DArrayOf3DArray((void***)reconVoxDataShort,UNIT_SINT16);

	if( pixPtr == NULL )
	{
		fprintf(stderr,"ERROR: reconVoxDatas are NULL\n");
		return;
	}

	for(int i=0;i<length;i++) *pixPtr++ = 0.0;
       	for(int i=0;i<length2;i++) *pixPtr2++ = 0;

}

void
clearSinogramImage(int projHeight, int projWidth, int sngrmHeight, float*** fltdSinogram)
{
	int length = projHeight*projWidth*sngrmHeight;
	float* pixPtr = (float*)get1DArrayOf3DArray((void***)fltdSinogram,UNIT_FLOAT32);

	if( pixPtr == NULL )
	{
		fprintf(stderr,"ERROR: fltdSinogram is NULL\n");
		return;
	}

	for(int i=0;i<length;i++) *pixPtr++ = 0.0;

}


void
shiftAngles(int numAngles,double* angles,double dAngle)
{
	if( angles == NULL )
	{
		fprintf(stderr,"WARNING: angles is NULL\n");
		return;
	}
	for(int i=0;i<numAngles;i++)	*angles++ += dAngle;
}

void
prepareAngles(int numAngles,double* angles,int reconSize,char* reconMaskAngle, double scale, const char* SMLClassString)
{
	if( reconMaskAngle == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for reconMaskAngle\n");
		return;
	}
	double x,y;
	double th1, thstart, thend;
	int fac;
	if(strcmp(SMLClassString, "S") != 0) 
	  {
	    //	    fprintf(stderr,"%s\n",SMLClassString );
	    if(strcmp(SMLClassString, "M") != 0) fac = 0;
	    for(int k=0;k<numAngles;k++)
	      {
		for(int j=0;j<reconSize;j++)
		  {
		    y = - (-reconSize/2 + j + 0.5) * scale;
		    for(int i=0;i<reconSize;i++)
		      {
			x = (-reconSize/2 + i + 0.5) * scale;
			if(y > 0) 
			  {
			    th1 = -atan(x/y);
			    //fprintf(stderr,"%lf %lf %lf %lf\n",x,y, th1*180/3.14150265,angles[k]*180/3.14150265 );
			  }
			else if(x < 0 && y < 0) 
			  {
			    th1 = (-atan(x/y)+PI);
			    //fprintf(stderr,"%lf %lf %lf %lf\n",x,y, th1*180/3.14150265,angles[k]*180/3.14150265 );
			  }
			
			else if(x > 0 && y < 0) 
			  {
			    th1 = -(atan(x/y)+PI);
			    //fprintf(stderr,"%lf %lf %lf %lf\n",x,y, th1*180/3.14150265,angles[k]*180/3.14150265 );
			  }
			//if(x > 10 && x < 20 && y > -20 && y < -10) fprintf(stderr,"%lf %lf %lf %lf\n",x,y, th1*180/3.14150265,angles[k]*180/3.14150265 );
			thstart = th1 + PI - 10*PI/180;
			thend = th1 + 10*PI/180;
			if(thstart < -PI) thstart = thstart + 2*PI;
			else if(thstart > PI) thstart = thstart - 2*PI;
			if(thend > PI) thend = thend - 2*PI;
			else if(thend < -PI) thend = thend + 2*PI;
			//if(x > 10 && x < 20 && y > -20 && y < -10) fprintf(stderr,"%lf %lf %lf %lf %lf\n",x,y, thstart*180/PI, thend*180/PI ,angles[k]*180/PI );

			if(thstart > thend)
			  {
			    if(angles[k] > thstart && angles[k] < PI )
			      {
				reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 1;
			      }
			    else if(angles[k] > -PI && angles[k] < thend )
			      {
				reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 1;
			      }
			    else
			      {
				reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 0;
			      }
			  }
			else
			  {
			    if(angles[k] > thstart && angles[k] < thend )
			    //if(angles[k] > -PI/2 && angles[k] < PI/2 )
			      {
				reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 1;
			      }
			    else
			      {
				reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 0;
			      }
			  }
		      }
		  }
	      }
	    
	  }
	else
	  {
	    for(int k=0;k<numAngles;k++)
	      {
		for(int j=0;j<reconSize;j++)
		  {
		    for(int i=0;i<reconSize;i++)
		      {
			reconMaskAngle[k*reconSize*reconSize+j*reconSize+i] = 1;
		      }
		  }
	      }
	  }
}
int
prepareTrigonometricTable(int numAngles,double* anglesRad,double** cosTable,double** sinTable)
{
	*cosTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );
	*sinTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );

	if( *cosTable == NULL || *sinTable == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for renderMaskData\n");
		return 0;
	}
	
	for(int i=0;i<numAngles;i++)
	{
		(*cosTable)[i] = cos(anglesRad[i]);
		(*sinTable)[i] = sin(anglesRad[i]);
	}

	return numAngles;
}

static clock_t start, stop;


void
startTimer()
{
	start = clock();
}

void
stopTimer()
{
	stop = clock();
}

void
showLapTime(char* label)
{
	fprintf(stderr,"%s : %3.2lf sec.\n",label,((double)stop-(double)start)/(double)CLOCKS_PER_SEC);
}























