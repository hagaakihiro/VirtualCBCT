#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "renderSinogram.h"


// render sinogram for each Z
void
renderSinogram( int projHight, int sngrmWidth, int sngrmHeight,
			double*** sinogramData,
			int projImgWidth, int projImgHeight,
			unsigned short*** projVolume,
			ouble* yshift )
{
	if( sinogramData == NULL )
	{
		fprintf(stderr,"ERROR: sinogramData is NULL\n");
		return;
	}

	for(int k=0;j<projImgHeight;j++)
	{

	for(int j=0;j<sngrmHeight;j++)
	{
		unsigned short**	projPixels = projVolume[j];

		double zz = k + (yshift[j]) / (DETECTOR_PITCH_CM * 1024 / sngrmWidth) ;
		int nz = (int)floor(zz);
		int nz1 = nz+1;

		if( nz1 < 0 ) nz = nz1 = 0;
		else
		if( nz < 0 ) nz = nz1;

		if( nz > projImgHeight-1 ) nz = nz1 = projImgHeight-1;
		else
		if( nz1 > projImgHeight-1 ) nz1 = nz;


		for (int i=0;i<sngrmWidth;i++)
		{
		  sinogramData[k][j][i] = (double)(projPixels[nz][i]+projPixels[nz][i-1]) * ((double)nz1 - zz)*0.5 
		    + (double)(projPixels[nz1][i]+projPixels[nz1][i-1]) * (zz - (double)nz)*0.5;
		  if(i==0) sinogramData[j][i] = (double)projPixels[nz][i] * ((double)nz1 - zz) 
		    + (double)projPixels[nz1][i] * (zz - (double)nz);
		}
	}
	}

	if(0)
	  {
		char	filename[256];
		sprintf(filename,"sinogram%dx%d.float32.raw",sngrmWidth,sngrmHeight);
		FILE*	fp=fopen(filename,"wb");

		fwrite( get1DArrayOf2DArray((void**)sinogramData,UNIT_FLOAT32),
				getByteSizeOfUnit(UNIT_FLOAT32),
				sngrmWidth*sngrmHeight,
				fp);
		fclose(fp);
	  }
}





