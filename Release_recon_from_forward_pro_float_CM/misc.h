
char**	prepareMaskData(int reconSize);

void	clearReconstructedImage(int reconSize, int zslices, int sngrmHeight, double*** reconVoxData, short*** reconVoxDataShort);
void    clearSinogramImage(int projHeight, int projWidth, int sngrmHeight, float*** fltdSinogram);

void	shiftAngles(int numAngles,double* angles,double dAngle);
void	prepareAngles(int numAngles,double* angles,int reconSize,char* reconMaskAngle, double reconScale, const char* SMLClassString);
int	prepareTrigonometricTable(int numAngles,double* anglesRad,double** cosTable,double** sinTable);

void	startTimer();
void	stopTimer();
void	showLapTime(char* label);



