

void	initializeGPU(int reconSize, int zslices, int projHeight, int projWidth, int sngrmHeight, double* anglesRad, 
                      double* cosTable, double* sinTable, double* xshift, double* yshift );
void	terminateGPU();

void	reconstructImageOnGPU( int reconSize, int zstart, int zslices, int reconStep, double*** reconVoxData, char** maskData, int xfactor,
			       int projHeight, int projWidth, int sngrmHeight, float*** sngrmVoxData, char* maskAngle, double reconScale, double DETECTOR_PITCH_MM_AT_ISO_CENT  );

void    filterSinogramOnGPU( int zstart, int zslices, int projHeight, int projWidth, int sngrmHeight,
			unsigned short*** projVolume, float*** fltdVoxData,
			     double* xshift, double* yshift, double DETECTOR_PITCH_MM_AT_ISO_CENT );
