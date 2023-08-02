#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mallocMD.h"
#include "misc.h"
#include "physParams.h"

int
prepareTrigonometricTable(int projImgWidth, int projImgHeight, int numAngles, double* anglesRad, double* xshift, double* yshift,
	double DETECTOR_PITCH_CM, double DETECTOR_PIXEL_NUM, double** cosTable, double** sinTable, double** XcenterShift, double** YcenterShift)
{
	*cosTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );
	*sinTable = (double*)new1DArray( numAngles, UNIT_FLOAT64 );
	*XcenterShift = (double*)new1DArray( numAngles, UNIT_FLOAT64 );
	*YcenterShift = (double*)new1DArray( numAngles, UNIT_FLOAT64 );

	if( *cosTable == NULL || *sinTable == NULL || *XcenterShift == NULL ||  *YcenterShift == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for renderMaskData\n");
		return 0;
	}
	
	for(int i=0;i<numAngles;i++)
	{
		(*cosTable)[i] = cos(anglesRad[i]);
		(*sinTable)[i] = sin(anglesRad[i]);

		//(*XcenterShift)[i] = -projImgWidth/2.0 - xshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgWidth));
		//(*XcenterShift)[i] = -projImgWidth/2.0 - xshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		(*XcenterShift)[i] = -projImgHeight/2.0 - xshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		(*YcenterShift)[i] = -projImgHeight/2.0 - yshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		//(*YcenterShift)[i] = -projImgHeight/2.0;// + yshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
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























