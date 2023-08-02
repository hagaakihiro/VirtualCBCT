
void	
reprojection( int ite, int reconSize, int reconSlices, unsigned short* reconImageshort, float* reprojection_float, float* npf_double, float* projVolume,
		      double reconScale, double reconStep,
		      int projImgWidth, int projImgHeight, int startProjNumber, int endProjNumber,
		      double* anglesRad, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor, double* xalpha,
		      double DETECTOR_PITCH_MM, double DETECTOR_PITCH_MM_AT_ISO_CENT,
		      double DIST_BTWN_SRC_AND_ISOCENT_MM, double DETECTOR_PIXEL_NUM,
		      double x_reg, double y_reg, double z_reg, float* Attenuation, float* pdf_pe, int num_material, int NE, int ke );






