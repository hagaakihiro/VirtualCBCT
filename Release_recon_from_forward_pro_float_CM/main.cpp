#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>

#include "physParams.h"
#include "mallocMD.h"
#include "misc.h"
#include "projectionDataMVkV.h"
#include "filterSinogram.h"
#include "reconstructImage.h"

int
main(int argc,char* argv[])
{
  
	// show usage unless command is OK
	if( argc != 2 )
	{
		fprintf(stderr,"\n[usage]$ cbct <info filename>\n\n");
		exit(0);
	}

	// all in one as data structure
	PROJECTION_DATA* projData = newProjectionData();


	// load data required for reconstuction (angle, projection images, etc.)
	if( loadData(argv[1],projData) == 0 )
	{
		deleteProjectionData( projData );
		exit(0);
	}


	FILE*	fp;
	FILE*   fp1;
	FILE*   fp2;

	if( (fp=fopen(projData->reconDataFilename,"wb")) == NULL )
	{
		fprintf(stderr,"ERROR: output file NOT writable\n");
		exit(0);
	}
	else
	{
	        float***	fltdSinogram = (float***)new3DArray( projData->projImgWidth, projData->projImgHeight, projData->usedProjNumber, UNIT_FLOAT32 );
		double***	reconImage = (double***)new3DArray( projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_FLOAT64 );
		short*** 	reconImageshort = (short***)new3DArray( projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_SINT16 );
		//short*** 	reconImageshort = (short***)new3DArray( projData->reconSize, projData->reconSize, 240, UNIT_SINT16 );
		double**	constMap = prepareConstMap( projData->projImgWidth, projData->DETECTOR_PITCH_CM_AT_ISO_CENT );//filterSinogram.cpp

		char**		reconMaskImage = prepareMaskData(projData->reconSize);//misc.cpp
		double		*cosTable, *sinTable;
		char*	        reconMaskangle =  (char*)new1DArray(projData->reconSize*projData->reconSize*projData->usedProjNumber,UNIT_UINT8);

		prepareAngles( projData->usedProjNumber, projData->anglesRad, projData->reconSize, reconMaskangle, projData->reconScale, projData->SMLClassString);////misc.cpp
		prepareTrigonometricTable( projData->usedProjNumber, projData->anglesRad, &cosTable, &sinTable );////misc.cpp

		/////////////////////////
		////////////////////////
		/////////////////////
		{
			double zstart = projData->reconOffset;
			double zend = projData->reconOffset + projData->reconSlices * projData->reconStep;

			fprintf(stderr,"Reconstructing z = %lf to %lf, Step = %lf \n",zstart, zend, projData->reconStep);

			/*char	filename[256];
			sprintf(filename,"sinogram%dx%d.unsignedshort.raw",projData->projImgWidth,projData->usedProjNumber);
			FILE*	fp2=fopen(filename,"wb");
			for(int th=0;th<projData->usedProjNumber;th++)
			  for(int m=0;m<projData->projImgWidth;m++)
			    {
			      unsigned short flt=projData->projVolume[th][256][m];
			      fwrite( &flt, getByteSizeOfUnit(UNIT_UINT16),1,fp2);
			    }
			    fclose(fp2);*/

			startTimer();
			clearSinogramImage( projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber, fltdSinogram);

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			
			
			float*     PProjVolume = (float*)new1DArray( projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber, UNIT_FLOAT32 );
			///SL 64
			//if( (fp2=fopen("36ProjImageshort_SL64_64x64x64_int16_smooth.raw","rb")) == NULL )  //new ver


		 
			 
			
	       // 2022/1/3(shimo) test  generalized(autometed) "input projection file" //
		   char	projImgFilename[1024];
		  fprintf(stderr,"%s.%s\n",projData->projImgFileDirname,projData->projImgFileExt);
		       sprintf( projImgFilename,"%s.%s",
				projData->projImgFileDirname,
				 projData->projImgFileExt );
		     
	       if( (fp2=fopen(projImgFilename,"rb")) == NULL )

			  {
			    fprintf(stderr,"ERROR:projection file NOT writable\n");
			    exit(0);
			  }
			
			fread(PProjVolume,getByteSizeOfUnit(UNIT_FLOAT32), projData->projImgHeight*projData->projImgWidth*projData->usedProjNumber,fp2);
			
			fclose(fp2);

			float***     PPPProjVolume = (float***)new3DArray( projData->projImgWidth, projData->projImgHeight,
									projData->usedProjNumber, UNIT_FLOAT32 );

			for(int i=0;i<projData->usedProjNumber ;i++)
			  for(int ii =0;ii<projData->projImgHeight ;ii++)
			    for(int iii =0; iii<projData->projImgWidth;iii++)
			      {

int iiii = i*projData->projImgHeight*projData->projImgWidth +ii*projData->projImgWidth+iii;


 PPPProjVolume[i][ii][iii]=PProjVolume[iiii];


			      }

				//directory projection read

			fprintf(stderr,"kokodayo\n");
			
			filterSinogram( projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber,
					PPPProjVolume, fltdSinogram,
						constMap,
						projData->xshift,
						projData->yshift, projData->DETECTOR_PITCH_CM_AT_ISO_CENT );


			
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			// "his" format filterSinogram
			
			/*
			filterSinogram( projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber,
					projData->projVolume, fltdSinogram,
						constMap,
						projData->xshift,
						projData->yshift, projData->DETECTOR_PITCH_CM_AT_ISO_CENT );
			*/

			

			
			stopTimer();	showLapTime("filterSinogram");	startTimer();
			if( (fp1=fopen("fltsngrm.raw","wb")) == NULL )
			  {
			    fprintf(stderr,"ERROR: fltsngrm file NOT writable\n");
			    exit(0);
			  }
			fwrite( get1DArrayOf3DArray((void***)fltdSinogram,UNIT_FLOAT32),
				getByteSizeOfUnit(UNIT_FLOAT32),
				projData->projImgWidth*projData->projImgHeight*projData->usedProjNumber,
				fp1);
			fclose(fp1);
			// clear and reconstruct image
			clearReconstructedImage( projData->reconSize, projData->reconSlices, projData->usedProjNumber, reconImage, reconImageshort);

			reconstructImage( projData->reconSize, projData->reconSlices, projData->reconOffset, projData->reconStep, projData->reconScale, 
					  reconImage, reconMaskImage,
					  projData->projImgHeight, projData->projImgWidth, projData->usedProjNumber,
					  fltdSinogram, projData->xfactor,
					  projData->anglesRad, cosTable, sinTable,
					  projData->xshift, projData->DETECTOR_PITCH_CM_AT_ISO_CENT );
			stopTimer();	showLapTime("reconstructImage");

			
for(int k=0;k<projData->reconSlices;k++)
			  for(int j=0;j<projData->reconSize;j++)
			    for(int i=0;i<projData->reconSize;i++)
			      {
			    
				reconImageshort[k][j][i] = (short)reconImage[k][j][i];
			      }

	fwrite( get1DArrayOf3DArray((void***)reconImageshort,UNIT_SINT16),
					getByteSizeOfUnit(UNIT_SINT16),
					projData->reconSize*projData->reconSize*projData->reconSlices,
		fp);///**short image, not unsigned short image!!!!!

			
			/*     //0 paddig for shepp logan

			for(int k=0;k<projData->reconSlices;k++)
			  for(int j=0;j<projData->reconSize;j++)
			    for(int i=0;i<projData->reconSize;i++)
			      {
			    
				reconImageshort[k+15][j][i] = (short)reconImage[k][j][i];
			      }

			for(int k=0;k<15;k++)
			  for(int j=0;j<projData->reconSize;j++)
			    for(int i=0;i<projData->reconSize;i++)
			      {
			    
				reconImageshort[k][j][i] = 0;
			      }
			for(int k=225;k<240;k++)
			  for(int j=0;j<projData->reconSize;j++)
			    for(int i=0;i<projData->reconSize;i++)
			      {
			    
				reconImageshort[k][j][i] = 0;
			      }
	fwrite( get1DArrayOf3DArray((void***)reconImageshort,UNIT_SINT16),
					getByteSizeOfUnit(UNIT_SINT16),
					projData->reconSize*projData->reconSize*240,
					fp);

			*/

			fprintf(stderr,"done.\n");
		}

		/////////////////////////
		////////////////////////
		/////////////////////

		
		fclose(fp);

	/////////////////////////// 
	
		delete3DArray( (void***)fltdSinogram, projData->projImgWidth, projData->projImgHeight, projData->usedProjNumber, UNIT_FLOAT32 );
		delete3DArray( (void***)reconImage, projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_FLOAT64 );
		delete3DArray( (void***)reconImageshort, projData->reconSize, projData->reconSize, projData->reconSlices, UNIT_SINT16 );
		delete2DArray( (void**)constMap, projData->projImgWidth, projData->projImgWidth, UNIT_FLOAT64 );
		delete2DArray( (void**)reconMaskImage, projData->reconSize, projData->reconSize, UNIT_UINT8 );
		delete1DArray( (void*)reconMaskangle);
		delete1DArray( (void*)cosTable );
		delete1DArray( (void*)sinTable );
	}

	deleteProjectionData( projData );
}


