
typedef struct projectionData
{
  int     totalProjNumber;
  
  int     usedZNumber;

  int     projImgWidth;
  int     projImgHeight;

  char	matImgFileDirname[256];
  char	outputname[256];
  int     zICposition;

  int     projImgFileNumDigits;
  int     projImgFile1stNumber;
  char	projImgFileExt[8];
  
  char	angleDataFilename[256];
  
  int     reconSize;
  char	reconTgtClassString[16];
  int     reconOffset;
  double  reconStep;
  int     reconSlices;
  int     FBP_slice;
  double  reconScale;
  double  xfactor;
  int     mat_num;
  
  char	reconDataFilename[256];
  
  double*	anglesRad;
  double*	xshift;
  double*	yshift;
  int*    inteangle;
  
  // Material image scales
  double	x_scale;
  double	y_scale;
  double        z_scale;
  // Material image scales
  double	x_reg;
  double	y_reg;
  double  z_reg;
  
  
  int     CM_NUM_Phase;
  double  angle_interval;

  int     startNum[5];
  int     endNum[5];

  double  DETECTOR_PITCH_CM;
  double  DETECTOR_PITCH_CM_AT_ISO_CENT;
  double  DIST_BTWN_SRC_AND_ISOCENT_CM;
  double  DETECTOR_PIXEL_NUM;
  double  DIST_BTWN_DETECTOR_AND_ISOCENT_CM;

  unsigned short*	matVolume;
  char    SpectrumString[256];
  char    mat_type[16];
} PROJECTION_DATA;


int	loadData(char* infoFilename,PROJECTION_DATA* projData);

PROJECTION_DATA*	newProjectionData();
void	deleteProjectionData(PROJECTION_DATA* projData);

int	loadProjectionData(PROJECTION_DATA* projData);
int	loadPixelData( char* filename, int offset, int length, unsigned short* pixels );

