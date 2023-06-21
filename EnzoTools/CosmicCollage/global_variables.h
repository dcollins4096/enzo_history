#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif


// from UserDefines.C
EXTERN int debug, verbosedebug, read_in_random_seed, read_in_num_projections, read_in_random_shifts;
EXTERN int NumberOfProjections, OutputProjectionGridSize, ProjectionType, 
  InputDatasetPrecision, NumberOfClusters;
EXTERN unsigned int random_seed;
EXTERN double AngularProjSize, ComovingBoxSize, Hubble0;

EXTERN char *OutputProjectionFileName, *InputDataFileName, *InputProjectionFileDatasetName,
  *RebinJustClustersFileNameList, *RebinNoClustersFileNameList, *ClusterPositionFileNameList,
  *NoClustersProjectionFileName, *JustClustersProjectionFileName,
  *OutputProjectionFileDatasetName, *InputProjectionFileNameList, *LogFileName, *RebinFileNameList;

// other global variables
EXTERN float *MasterProjectionArray;
EXTERN float *RebinArray;
EXTERN float *InputProjectionArrayFloat;
EXTERN double *InputProjectionArrayDouble;

EXTERN double *InputProjectionArrayNoClustersDouble;
EXTERN double *InputProjectionArrayJustClustersDouble;
EXTERN float *InputProjectionArrayNoClustersFloat;
EXTERN float *InputProjectionArrayJustClustersFloat;
EXTERN float *RebinArrayNoClusters;
EXTERN float *RebinArrayJustClusters;
EXTERN float *MasterProjectionArrayNoClusters;
EXTERN float *MasterProjectionArrayJustClusters;

EXTERN double *xshifts, *yshifts;
EXTERN int *random_axes;

EXTERN double *haloxpos, *haloypos, *halozpos, *haloRvir, *haloYdec[HALO_YDEC_BINS], *haloredshift;

// MAX_PROJECTION_FILES is set in standard_includes.h
EXTERN char *InputProjectionFileNames[3][MAX_PROJECTION_FILES];  
EXTERN char *ClusterPositionFileNames[3][MAX_PROJECTION_FILES];  
EXTERN char *RebinnedProjectionFileNames[MAX_PROJECTION_FILES];  
EXTERN char *RebinnedProjectionNoClustersFileNames[MAX_PROJECTION_FILES];  
EXTERN char *RebinnedProjectionJustClustersFileNames[MAX_PROJECTION_FILES];  

EXTERN double gridratio[MAX_PROJECTION_FILES], projredshift[MAX_PROJECTION_FILES];

EXTERN int InputProjectionGridSize[MAX_PROJECTION_FILES];


