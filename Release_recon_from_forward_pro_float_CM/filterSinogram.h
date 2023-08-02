
double**	prepareConstMap(int projWidth, double DETECTOR_PITCH_CM_AT_ISO_CENT);

void	filterSinogram( int projHeight,	int projWidth, int sngrmHeight,
				float*** projVolume, float*** fltdVoxData,
				double** constMap, double* xshift, double* yshift, double DETECTOR_PITCH_CM_AT_ISO_CENT );
