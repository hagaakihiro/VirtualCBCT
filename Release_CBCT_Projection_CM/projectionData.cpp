#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#include <math.h>
#endif

#include "physParams.h"
#include "mallocMD.h"
#include "projectionData.h"

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
loadDouble(FILE* fp,double* value)
{
	skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf",value);
}

static void
loadDouble3(FILE* fp, double* value1, double* value2, double* value3)
{
  skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%lf %lf %lf",value1,value2,value3);
}

static void
loadInteger10(FILE* fp, int* value1, int* value2, int* value3, int* value4, int* value5, int* value6, int* value7, int* value8, int* value9, int* value10)
{
  skipLines(fp,2);	fgets(lineBuf,LINEBUF_SIZE,fp);	sscanf(lineBuf,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",value1,value2,value3,value4,value5,value6,value7,value8,value9,value10);
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
  loadString( fp, projData->SpectrumString );		// Energy Spectrum
  loadInteger( fp, &projData->totalProjNumber );		// total projections
  loadInteger( fp, &projData->projImgWidth );		// projection image width	
  loadInteger( fp, &projData->projImgHeight );		// projection image height
  loadString( fp, projData->angleDataFilename );		// angle data file name
  loadString( fp, projData->matImgFileDirname );		// material image file directory name
  loadInteger( fp, &projData->zICposition );		// IC z position in z slices =207
  loadString( fp, projData->outputname );		// output projection name
  loadDouble3( fp, &projData->x_scale, &projData->y_scale, &projData->z_scale );// Pixel scale
  loadInteger( fp, &projData->reconSize );		// reconstruction size (square only)
  loadInteger( fp, &projData->reconSlices );		// reconstruction size (square only)
  loadInteger( fp, &projData->mat_num);	   	        // Number of Materials used
  loadString( fp, projData->mat_type );		        // human, gammex, catphan, etc
  loadDouble( fp, &projData->xfactor );		        // factor for density value
  loadDouble3( fp, &projData->x_reg, &projData->y_reg, &projData->z_reg );// Material image shift
  fclose(fp);

  projData->DETECTOR_PITCH_CM = DETECTOR_PITCH_CM_kV;
  projData->DETECTOR_PITCH_CM_AT_ISO_CENT = DETECTOR_PITCH_CM_AT_ISO_CENT_kV;
  projData->DIST_BTWN_SRC_AND_ISOCENT_CM = DIST_BTWN_SRC_AND_ISOCENT_CM_kV;
  projData->DETECTOR_PIXEL_NUM = DETECTOR_PIXEL_NUM_kV;
  projData->DIST_BTWN_DETECTOR_AND_ISOCENT_CM = DIST_BTWN_DETECTOR_AND_ISOCENT_CM_kV;
  
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
  
  ret->matVolume = NULL;
  
  return ret;
}

void
deleteProjectionData(PROJECTION_DATA* projData)
{
  if( projData == NULL )	return;

  delete1DArray( (void*)projData->anglesRad );
  delete1DArray( (void*)projData->xshift );
  delete1DArray( (void*)projData->yshift );
  delete1DArray( (void*)projData->matVolume );

  free(projData);
}



// load angle, shift, pixels
int
loadProjectionData(PROJECTION_DATA* projData)
{
  FILE*	fp;
  int	lineCount=0, tgtCount=0;
  double	angle, x, y;
  char	classString[16];
  char  matImgFilename0[256];


  if( (fp=fopen(projData->angleDataFilename,"r")) == NULL )
    {
      fprintf(stderr,"angle data file:%s not found\n",projData->angleDataFilename);
      return 0;
    }
  
  // prepare arrays
  projData->anglesRad = (double*)new1DArray(projData->totalProjNumber,UNIT_FLOAT64);
  projData->xshift = (double*)new1DArray(projData->totalProjNumber,UNIT_FLOAT64);
  projData->yshift = (double*)new1DArray(projData->totalProjNumber,UNIT_FLOAT64);
  // check target projection number
  int i = 0;
  while( fscanf(fp,"%lf %lf %lf %s", &angle,&x,&y,classString) != EOF ) 
    {
      projData->anglesRad[i]=angle * PI / 180.0;
      projData->xshift[i]=x;
      projData->yshift[i]=y;
      //fprintf(stdout,"%lf %d\n",projData->anglesRad[i], i);
      i++;
    }
 
  double z_length = (projData->DIST_BTWN_SRC_AND_ISOCENT_CM+projData->reconSize/2*projData->y_scale)
    /(projData->DIST_BTWN_DETECTOR_AND_ISOCENT_CM+projData->DIST_BTWN_SRC_AND_ISOCENT_CM)
    * projData->DETECTOR_PITCH_CM * projData->DETECTOR_PIXEL_NUM;//shimo 2022/1/20
  int iz_length = z_length/projData->z_scale;
  int linemin = projData->zICposition - (iz_length/2+1);//207-(40/2+1)=186
  int linemax = projData->zICposition + (iz_length/2+1);//207+(40/2+1)=228 -> 222 is max,so 222


  
  printf("iz_length %d linemin %d linemax %d \n", iz_length, linemin, linemax);//iz_length=40
  if(linemin < 0) linemin = 0;//186
  //if(linemax > projData->reconSlices) linemax = projData->reconSlices;//222
  if(linemax >= projData->reconSlices) linemax = projData->reconSlices;//222
  int data_size = projData->reconSize*projData->reconSize*projData->mat_num;
  unsigned short* pixels = (unsigned short*)new1DArray( data_size, UNIT_UINT16 );
  projData->matVolume = (unsigned short*)new1DArray( data_size*
						      linemax, UNIT_UINT16 );
  int t_count = 0;
  int k = 0;
  for( i = linemin; i < linemax;i++ )
    {
  
      
      //water cylinder
      sprintf(matImgFilename0, "%s_%03d.raw",projData->matImgFileDirname,i); //25cm 
       
      
      
      loadPixelData( matImgFilename0, 0, data_size * sizeof(unsigned short), pixels);
      printf("%s \n", matImgFilename0);


      
       for(int kk = 0; kk < data_size; kk++)
	{
	 
	  projData->matVolume[t_count*data_size+kk]= pixels[kk];// unsigned short* reconImageshort
	}


    
      t_count += 1;
    }
  fclose(fp);
  projData->usedZNumber = t_count;//int reconSlices & z_required_range
  projData->zICposition = projData->zICposition - linemin;  //207(IR.txtfile)-linemin
  // projection volume

  printf("projData->usedZNumber=%d \n",projData->usedZNumber);

    /*
     // Check the material image
  char	filename[256];
  sprintf(filename,"material_img.raw");
  FILE*	fp1=fopen(filename,"wb");
  
  fwrite( projData->matVolume,
	  getByteSizeOfUnit(UNIT_UINT16),
	  data_size*projData->usedZNumber,fp1);
  
  fclose(fp1);
	 */

  return t_count;
}


// load pixels
int
loadPixelData( char* filename, int offset, int length, unsigned short* pixels )
{
	FILE*	fp;
	int	readSize;
	if( (fp=fopen(filename,"rb")) == NULL )
	{
		fprintf(stderr,"data file: %s not found\n",filename);
		return -1;
	}

	fseek(fp,offset,SEEK_SET);
	readSize = fread( (void*)pixels, 1, length, fp );
	fclose(fp);

	if( readSize != length )
	{
		fprintf(stderr,"WARNING: pixel data size unmatch\n");
	}

	return  readSize;
}






