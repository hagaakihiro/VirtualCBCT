
typedef struct projectionData
{
	int	totalProjNumber;
        char    GantRotationDir[2];
	int	usedProjNumber;

	int	projImgWidth;
	int	projImgHeight;

	char	projImgFileDirname[256];

	int	projImgFileOffset;

	int	projImgFileNumDigits;
	int	projImgFile1stNumber;
	char	projImgFileExt[8];

	char	angleDataFilename[256];

	int	reconSize;
	char	reconTgtClassString[16];
	int	reconOffset;
	double	reconStep;
	int	reconSlices;
	int	xfactor;
	char	kVMVClassString[16];
	char	SMLClassString[8];
	double	reconScale;

	char	reconDataFilename[256];

	double*	anglesRad;
	double*	xshift;
	double*	yshift;
        int*    inteangle;

        double DETECTOR_PITCH_CM_AT_ISO_CENT;

	unsigned short***	projVolume;

} PROJECTION_DATA;


int	loadData(char* infoFilename,PROJECTION_DATA* projData);

PROJECTION_DATA*	newProjectionData();
void	deleteProjectionData(PROJECTION_DATA* projData);

int	loadProjectionData(PROJECTION_DATA* projData);
int	loadPixelData( char* filename, int offset, int length, unsigned short* pixels );

