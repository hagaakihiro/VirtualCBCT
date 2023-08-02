int	prepareTrigonometricTable(int projImgWidth, int projImgHeight, int numAngles, double* anglesRad, double* xshift, double* yshift,
				  double DETECTOR_PITCH_CM, double DETECTOR_PIXEL_NUM, double** cosTable, double** sinTable, double** XcenterShift, double** YcenterShift);

void	startTimer();
void	stopTimer();
void	showLapTime(char* label);



