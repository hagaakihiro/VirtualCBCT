#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "physParams.h"
#include "mallocMD.h"
#include "misc.h"
#include "projectionData.h"
#include "reprojection.h"
#include "reprojectionOnGPU.h"
#include "lengthproj.h"
#include "IR_ImageCBCT.h"
#include "virtual_projection.h"
#define USE_GPU // Comment out when GPU is not used.

int
main(int argc,char* argv[])
{
  
    // show usage unless command is OK
    if( argc != 2 )
        {
            fprintf(stderr,"\n[usage]$ proj.exe <info filename>\n\n");
            exit(0);
        }

    // all in one as data structure
    PROJECTION_DATA* projData = newProjectionData();
  
  
    // load data required for reprojection (angle, spectrum, etc.)
    if( loadData(argv[1],projData) == 0 )
        {
            deleteProjectionData( projData );
            exit(0);
        }

  
    FILE*	fp;
    FILE*	fp1;
        
    FILE*	fp22;
    FILE*	fp3;
    FILE*   fp4;
    char filename[128];
  
    int startProjNumber, endProjNumber;
  
    //#################################
  
    char spectrum_data_name[128];
    sprintf(spectrum_data_name,"%s",projData->SpectrumString);
    printf("%s \n",spectrum_data_name);
  
    //#################################
    
    float* Attenuation = (float*)new1DArray( 1024*projData->mat_num, sizeof(float) );
    float* pdf_pe = (float*)new1DArray( 1024, sizeof(float) );
  
  
    int NE;
    int nthin=1;

  
    attenuation_coefficient( spectrum_data_name, projData->mat_num, Attenuation, pdf_pe, &NE, nthin, projData->mat_type);//prior_weight_production.cpp


  
    int z_required_range=projData->usedZNumber;
    printf("%d \n",z_required_range);
  
  
    short*    req_reconImageshort = (short*)new1DArray( projData->reconSize*projData->reconSize*z_required_range, UNIT_SINT16 );
    short*    reconImageshort = (short*)new1DArray( projData->reconSlices*projData->reconSize*projData->reconSize, UNIT_SINT16 );
    short*    reconImageshort_min = (short*)new1DArray( projData->reconSlices*projData->reconSize*projData->reconSize, UNIT_SINT16 );
    short*    FBP_reconImageshort = (short*)new1DArray( projData->reconSize*projData->reconSize*projData->FBP_slice, UNIT_SINT16 );
    float*    reconImagefloat = (float*)new1DArray( projData->reconSlices*projData->reconSize*projData->reconSize, UNIT_FLOAT32 );
    float*    req_npf_float = (float*)new1DArray( projData->reconSize*projData->reconSize*z_required_range, UNIT_FLOAT32 );
    float*    req_npf_float_re = (float*)new1DArray( projData->reconSize*projData->reconSize*z_required_range, UNIT_FLOAT32 );
    float*    npf_float = (float*)new1DArray( projData->reconSlices*projData->reconSize*projData->reconSize, UNIT_FLOAT32 );
    float*    npf_float_re = (float*)new1DArray( projData->reconSlices*projData->reconSize*projData->reconSize, UNIT_FLOAT32 );

    float* reprojection_float = (float*)new1DArray( projData->projImgWidth*projData->projImgHeight*projData->totalProjNumber, UNIT_FLOAT32 );
    float* reprojection_float_full = (float*)new1DArray( projData->projImgWidth*projData->projImgHeight*projData->totalProjNumber, UNIT_FLOAT32 );
    
    float* ProjVolume = (float*)new1DArray( projData->projImgWidth*projData->projImgHeight*projData->totalProjNumber, UNIT_FLOAT32 );//Haga shimo

 
    double	  *cosTable, *sinTable;
    double	  *XcenterShift, *YcenterShift;
    prepareTrigonometricTable( projData->projImgWidth, projData->projImgHeight, projData->totalProjNumber, projData->anglesRad, projData->xshift, projData->yshift,
			     projData->DETECTOR_PITCH_CM, projData->DETECTOR_PIXEL_NUM, &cosTable, &sinTable, &XcenterShift, &YcenterShift );//misc.cpp

    startProjNumber = 0;
    endProjNumber = projData->totalProjNumber;
    fprintf(stderr,"startProjNumber = %d,endProjNumber = %d \n",startProjNumber,endProjNumber);

    int ite = 0;
    double xalpha[1];

    fprintf(stderr,"matImage num = %d, center = %d \n",projData->usedZNumber, projData->zICposition);
    projData->reconScale = projData->x_scale;
    projData->reconStep = projData->z_scale;
    projData->z_reg = (projData->usedZNumber/2 - projData->zICposition)*projData->z_scale;//projData->zICposition = projData->zICposition - linemin;
  
    fprintf(stderr,"projData->usedZNumber/2 = %d, projData->zICposition = %d \n",projData->usedZNumber/2, projData->zICposition);
    fprintf(stderr,"xreg = %lf, zreg = %lf \n",projData->x_reg, projData->z_reg);

 
  
#ifdef USE_GPU
    initializeGPU( projData->reconSize, projData->reconScale, z_required_range, projData->reconStep,
		 projData->projImgHeight, projData->projImgWidth, projData->totalProjNumber,
		 cosTable, sinTable, XcenterShift, YcenterShift,
		 projData->x_reg, projData->y_reg, projData->z_reg,
		 projData->DETECTOR_PITCH_CM, projData->DETECTOR_PITCH_CM_AT_ISO_CENT,
		 projData->DIST_BTWN_SRC_AND_ISOCENT_CM, projData->DETECTOR_PIXEL_NUM,
		 projData->mat_num);
  
    reprojectionOnGPU( ite, projData->reconSize, z_required_range, projData->matVolume, reprojection_float, req_npf_float_re, ProjVolume,
		projData->reconScale, projData->reconStep,
		projData->projImgWidth, projData->projImgHeight, startProjNumber, endProjNumber,
		projData->anglesRad, cosTable, sinTable, XcenterShift, YcenterShift, projData->xfactor, xalpha,
		projData->DETECTOR_PITCH_CM, projData->DETECTOR_PITCH_CM_AT_ISO_CENT,
		projData->DIST_BTWN_SRC_AND_ISOCENT_CM, projData->DETECTOR_PIXEL_NUM,
		     projData->x_reg, projData->y_reg, projData->z_reg, Attenuation, pdf_pe, projData->mat_num, NE, projData->outputname);

#else
    for(int ke = 0; ke < NE; ke++)//energy
        {
            printf("%d ke/ %d NE \n",ke,NE);
            reprojection( ite, projData->reconSize, z_required_range, projData->matVolume, reprojection_float, req_npf_float_re, ProjVolume,
                         projData->reconScale, projData->reconStep,
                         projData->projImgWidth, projData->projImgHeight, startProjNumber, endProjNumber,
                         projData->anglesRad, cosTable, sinTable, XcenterShift, YcenterShift, projData->xfactor, xalpha,
                         projData->DETECTOR_PITCH_CM, projData->DETECTOR_PITCH_CM_AT_ISO_CENT,
                         projData->DIST_BTWN_SRC_AND_ISOCENT_CM, projData->DETECTOR_PIXEL_NUM,
                         projData->x_reg, projData->y_reg, projData->z_reg, Attenuation, pdf_pe, projData->mat_num, NE, ke );

            printf("%f %f %f %f \n",reprojection_float[0],reprojection_float[262144],reprojection_float[1],reprojection_float[130000]);
      
            for(int i = 0; i < projData->projImgWidth*projData->projImgHeight*projData->totalProjNumber; i++) reprojection_float_full[i] += reprojection_float[i];
        }
    ///#######CPU
    sprintf(filename,"ProjImage_float.raw");
    fp3=fopen(filename,"wb");
    fwrite( reprojection_float_full,
        getByteSizeOfUnit(UNIT_FLOAT32),
        projData->projImgWidth*projData->projImgHeight*projData->totalProjNumber,fp3);//Reprojection.raw
    fclose(fp3);

#endif
  

  //########

 

  
  /////////////////////////// 
	
    delete1DArray( (void*)reconImageshort);
    delete1DArray( (void*)reconImageshort_min);
    delete1DArray( (void*)reconImagefloat);
    delete1DArray( (void*)req_reconImageshort);
    delete1DArray( (void*)FBP_reconImageshort);
    delete1DArray( (void*)reprojection_float);
    delete1DArray( (void*)reprojection_float_full);
    //delete1DArray( (void*)npf_double);
    //delete1DArray( (void*)npf_double_re);
    delete1DArray( (void*)npf_float);
    delete1DArray( (void*)npf_float_re);

    deleteProjectionData( projData );
}


