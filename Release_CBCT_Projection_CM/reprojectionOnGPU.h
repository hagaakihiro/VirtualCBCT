
void initializeGPU(int reconSize, double reconScale, int reconSlices, double reconStep,
	      int projImgHeight, int projImgWidth, int totalProjNumber,
	      double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift,
	      double x_reg, double y_reg, double z_reg,
	      double DETECTOR_PITCH_CM, double DETECTOR_PITCH_CM_AT_ISO_CENT,
	      double DIST_BTWN_SRC_AND_ISOCENT_CM, double DETECTOR_PIXEL_NUM,
		   int num_material);

void	terminateGPU();

void	
reprojectionOnGPU( int ite, int reconSize, int reconSlices, unsigned short* reconImageshort, float* reprojection_float, float* npf_double,float* projVolume,//float* projVolume,
		      double reconScale, double reconStep,
		      int projImgWidth, int projImgHeight, int startProjNumber, int endProjNumber,
		      double* anglesRad, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor, double* xalpha,
		      double DETECTOR_PITCH_CM, double DETECTOR_PITCH_CM_AT_ISO_CENT,
		      double DIST_BTWN_SRC_AND_ISOCENT_CM, double DETECTOR_PIXEL_NUM,
		   double x_reg, double y_reg, double z_reg, float* Attenuation, float* pdf_pe, int num_material, int NE,char* outputname);



