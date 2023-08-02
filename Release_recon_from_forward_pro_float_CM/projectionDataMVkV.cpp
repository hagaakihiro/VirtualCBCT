
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "projectionDataMVkV.h"


#define LINEBUF_SIZE	512

static char	lineBuf[LINEBUF_SIZE];

//////

static void
skipLines(FILE* fp,int nLines)
{
	for(int i=0;i<nLines;i++)	fgets(lineBuf,LINEBUF_SIZE,fp);
}

static void
loadInteger(FILE* fp,int* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%d",value);
}

static void
loadIntegerHex(FILE* fp,int* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%X",value);
}

static void
loadString(FILE* fp,char* string)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%s",string);
}

static void
loadDouble(FILE* fp,double* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf",value);
}

//////

int
loadData(char* infoFilename,PROJECTION_DATA* projData)
{
	FILE*	fp;

	if( (fp=fopen(infoFilename,"r")) == NULL )
	{
		fprintf(stderr,"info file: $s not found\n",infoFilename);
		return 0;
	}

	// skip comments (6lines)
	skipLines(fp,6);
	
	// load each element
	loadString( fp, projData->kVMVClassString );		// Energy label
	loadInteger( fp, &projData->totalProjNumber );		// total projections
	loadInteger( fp, &projData->projImgWidth );		// projection image width	
	loadInteger( fp, &projData->projImgHeight );		// projection image height
	loadInteger( fp, &projData->projImgFileOffset );	// projection image offset
	loadString( fp, projData->projImgFileDirname );		// projection image file directory name
	loadInteger( fp, &projData->projImgFileNumDigits );	// projection image file number digits
	if(strcmp(projData->kVMVClassString,"kV")==0)
	  {
	    loadIntegerHex( fp, &projData->projImgFile1stNumber );	// first projection image file number in hex
	    projData->DETECTOR_PITCH_CM_AT_ISO_CENT = DETECTOR_PITCH_CM_AT_ISO_CENT_kV;
	  }
	if(strcmp(projData->kVMVClassString,"MV")==0)
	  {
	    loadInteger( fp, &projData->projImgFile1stNumber );	// first projection image file number
	    projData->DETECTOR_PITCH_CM_AT_ISO_CENT = DETECTOR_PITCH_CM_AT_ISO_CENT_MV;
	  }
	loadString( fp, projData->projImgFileExt );		// projection image file extension
	loadString( fp, projData->angleDataFilename );		// angle data file name
	loadInteger( fp, &projData->reconSize );		// reconstruction size (square only)
	loadString( fp, projData->reconTgtClassString );	// reconstruction target class label
	loadInteger( fp, &projData->reconOffset );		// reconstruction start z position in cm
	loadDouble( fp, &projData->reconStep );	          	// reconstruction z step in cm
	loadInteger( fp, &projData->reconSlices );		// number of reconstructed slices
	loadDouble( fp, &projData->reconScale );		// Reconstruction size facor;
	loadInteger( fp, &projData->xfactor );		        // factor for CT value
	loadString( fp, projData->reconDataFilename );		// reconstruction file name

	fclose(fp);

	return loadProjectionData( projData );
}



PROJECTION_DATA*
newProjectionData()
{
	PROJECTION_DATA* ret = (PROJECTION_DATA*)malloc(sizeof(PROJECTION_DATA));

	if( ret == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for newProjectionData()\n");
		return NULL;
	}

	ret->anglesRad = NULL;
	ret->xshift = NULL;
	ret->yshift = NULL;

	ret->projVolume = NULL;

	return ret;
}

void
deleteProjectionData(PROJECTION_DATA* projData)
{
	if( projData == NULL )	return;

	delete1DArray( (void*)projData->anglesRad );
	delete1DArray( (void*)projData->xshift );
	delete1DArray( (void*)projData->yshift );

	delete3DArray( (void***)projData->projVolume,
				projData->projImgWidth, projData->projImgHeight, projData->usedProjNumber,
				UNIT_UINT16);

	free(projData);
}



// load angle, shift, pixels
int
loadProjectionData(PROJECTION_DATA* projData)
{
	FILE*	fp;
	int	lineCount=0, tgtCount=0, i=0;
	double	angle, x, y;
	char	classString[16];


	if( (fp=fopen(projData->angleDataFilename,"r")) == NULL )
	{
		fprintf(stderr,"angle data file: $s not found\n",projData->angleDataFilename);
		return 0;
	}

	// check target projection number
	double cwccw0,cwccw;
	double preangle[MAX_Proj_Num],prex[MAX_Proj_Num],prey[MAX_Proj_Num];
	int phasecharac[MAX_Proj_Num];
	while( fscanf(fp,"%lf %lf %lf %s", &angle,&x,&y,classString) != EOF ) 
	  {
	    i++;
	    if(i == 1) cwccw0 = angle;
	    if(i == 30) cwccw = (cwccw0 - angle);
	    // * means any
	    phasecharac[lineCount]=0;
	    if( !strcmp(projData->reconTgtClassString,"*") ||
		!strcmp(classString,projData->reconTgtClassString) )	
	      {
	       phasecharac[lineCount]=1;
	       tgtCount++;
	      }
	    preangle[lineCount]=angle;
	    prex[lineCount]=x;
	    prey[lineCount]=y;
	    //fprintf(stdout,"%lf %d\n",preangle[lineCount],phasecharac[lineCount]);
	    lineCount++;
	  }
	fprintf(stdout,"%d %d %lf %lf\n",tgtCount,lineCount, cwccw0, cwccw);
	// Phase0=1 for CCW while Phase0=-1 for CW
	double preangle1[MAX_Proj_Num],prex1[MAX_Proj_Num],prey1[MAX_Proj_Num];
	int phasecharac1[MAX_Proj_Num];
	int phase0 = 1, linemax=lineCount-1;
	i=0;
	if( cwccw < 0) // CW
	  {
	    phase0 = -1;
	    lineCount=0;
	    while ( lineCount <= linemax )
	      {
		preangle1[lineCount]=preangle[linemax-lineCount];
		prex1[lineCount]=prex[linemax-lineCount];
		prey1[lineCount]=prey[linemax-lineCount];
		phasecharac1[lineCount]=phasecharac[linemax-lineCount];
		//fprintf(stdout,"%lf %d %d\n",preangle1[lineCount],phasecharac1[lineCount],lineCount);
		lineCount++;
	      }
	  }
	else           // CCW
	  {
	    phase0 = 1;
	    lineCount=0;
	    while ( lineCount <= linemax )
	      {
		preangle1[lineCount]=preangle[lineCount];
		prex1[lineCount]=prex[lineCount];
		prey1[lineCount]=prey[lineCount];
		phasecharac1[lineCount]=phasecharac[lineCount];
		fprintf(stdout,"%lf %d %d\n",preangle1[lineCount],phasecharac1[lineCount],lineCount);
		lineCount++;
	      }
	  }

	//fprintf(stdout,"%lf %lf %d\n",cwccw,cwccw0,phase0);
	projData->usedProjNumber = projData->totalProjNumber;
	projData->totalProjNumber = tgtCount++;     //target count (mark)


	// check lines in angle data file
	if( lineCount != projData->totalProjNumber )
	{
		fprintf(stderr,"WARNING: total projection number mismatch (L=%d != %d)\n",lineCount,projData->totalProjNumber);
		//projData->usedProjNumber = projData->totalProjNumber;
	}
	else
	{
		fprintf(stderr,"%d of %d projections are used.\n",projData->usedProjNumber,projData->totalProjNumber);
		//projData->usedProjNumber = projData->totalProjNumber;
	}
	
	//projData->usedProjNumber = tgtCount;
	// prepare arrays
	projData->anglesRad = (double*)new1DArray(projData->usedProjNumber,UNIT_FLOAT64);
	projData->xshift = (double*)new1DArray(projData->usedProjNumber,UNIT_FLOAT64);
	projData->yshift = (double*)new1DArray(projData->usedProjNumber,UNIT_FLOAT64);
	

	if( projData->anglesRad == NULL ||
		projData->xshift == NULL || projData->yshift == NULL )
	{
		fprintf(stderr,"ERROR: not enough memory for loadAngleData()\n");
		return 0;
	}

	// rerurn to file head
	fseek(fp,0,SEEK_SET);


	// load angle, shift, and pixels
	lineCount = tgtCount = 0;
	double xaverage = 0;
	int indx0 = 0, indx00 = projData->usedProjNumber-1;
	//while( fscanf(fp,"%lf %lf %lf %s", &angle,&x,&y,classString) != EOF ) 
        //while ( tgtCount < projData->usedProjNumber )
	for( i=0; i <= linemax;i++ )
	{
	  if(phasecharac1[i]==1)
	    {
	      indx0 = tgtCount;
	      angle=preangle1[i];
	      x=prex1[i];
	      y=prey1[i];
	      // image file read
	      
	      
	      projData->anglesRad[indx0] = angle * PI / 180.0;
	      projData->xshift[indx0] = x;
	      projData->yshift[indx0] = y;
	      xaverage += x;
	      fprintf(stdout,"%lf %lf %lf %d %d \n",angle,x,y,indx0, i);
	      
	      tgtCount++;    //target projection number
	      lineCount++;// all projection number
	    }
	}
	fclose(fp);
	xaverage = xaverage/lineCount;
	if(xaverage < 6) sprintf(projData->SMLClassString, "%s",  "S");
	else if(xaverage > 15) sprintf(projData->SMLClassString, "%s",  "L");
	else sprintf(projData->SMLClassString, "%s",  "M");
	//	sprintf(projData->SMLClassString, "%s",  "S");
	fprintf(stderr,"Panel Offset %s\n",projData->SMLClassString);

	// reconstruction slice number check
	if( projData->reconStep > 0 &&
		projData->reconOffset + projData->reconStep * (projData->reconSlices-1) >= projData->projImgHeight )
	{
		fprintf(stderr,"WARNING: reconstruction slice number (%d) invalid\n",projData->reconSlices);

		projData->reconSlices = (projData->projImgHeight-projData->reconOffset)/projData->reconStep;

		fprintf(stderr,"         corrected to (%d)\n",projData->reconSlices);
	}
	else
	if( projData->reconStep < 0 &&
		projData->reconOffset + projData->reconStep * (projData->reconSlices-1) < 0 )
	{
		fprintf(stderr,"WARNING: reconstruction slice number (%d) invalid\n",projData->reconSlices);

		projData->reconSlices = (projData->reconOffset+1)/(-projData->reconStep);

		fprintf(stderr,"         corrected to (%d)\n",projData->reconSlices);
	}

	// Scatter correction
	//for(int th=0;th<lineCount;th++)
	//  {
	//  for(int ii=0;ii<projData->projImgHeight;ii++)
	//    {
	//    for(int i=0;i<projData->projImgWidth;i++)
	//      {
	//	//if(projData->projVolume[th][ii][i] > 25000) 
	//	  {
	//	    //projData->projVolume[th][ii][i] = (int)(projData->projVolume[th][ii][i]*0.33);
	//	  }
	//      }
	//    }
	//  }


	return tgtCount;
}


// load pixels
int
loadPixelData( char* filename, int offset, int length, unsigned short* pixels )
{
	FILE*	fp;
	int	readSize;
	//fprintf(stderr,"loading %s\n",filename);
	if( (fp=fopen(filename,"rb")) == NULL )
	{
		fprintf(stderr,"data file: %s not found\n",filename);
		return -1;
	}

	//fprintf(stderr,"loading %s\n",filename);

	fseek(fp,offset,SEEK_SET);
	readSize = fread( (void*)pixels, 1, length, fp );
	fclose(fp);

	if( readSize != length )
	{
		fprintf(stderr,"WARNING: pixel data size unmatch\n");
	}

	return  readSize;
}






