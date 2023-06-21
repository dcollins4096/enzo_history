#ifndef ANALYSIS_BASE_CLASS_H
#define ANALYSIS_BASE_CLASS_H

/***********************************************************************
/  
/ ANALYSIS SUPERCLASS
/
************************************************************************/

class AnalysisBaseClass{
  
 public:
  AnalysisBaseClass( TopGridData *metadata,
		HierarchyEntry *topgrid,
		int level = -1,
		FLOAT *Left = NULL,
		FLOAT *Right = NULL);

  ~AnalysisBaseClass();

  void PrintParameters();

  void PrintGridInfo( grid *Grid);

  void SetAnalysisRegion( FLOAT Left[MAX_DIMENSION],
			  FLOAT Right[MAX_DIMENSION] );

  int OpenData();
  int CloseData();

  int NumberOfGridsInVolume( );
  void InitializeFastFind();
  int IndexGrids();
  void IndexGridsByProcessor();
  void BuildMyGridArray();

  HierarchyEntry *FindGrid( HierarchyEntry *Grid, float *point);
  HierarchyEntry *FastFindGrid(float *point);

  AnalysisBaseClass *Next;

#ifdef UNIT_TESTING
 public:
#else 
 protected:
#endif

  const float G;
  const float MPC_CM;
  const float MSOLAR_G;

  TopGridData *MetaData;
  HierarchyEntry *TopGrid;
  int MaximumAnalysisLevel;
  FLOAT AnalysisRegionLeftEdge[MAX_DIMENSION];
  FLOAT AnalysisRegionRightEdge[MAX_DIMENSION];
  FLOAT DomainDelta[MAX_DIMENSION];


  void PrintTrueFalse(FILE *buffer, bool testvalue,
		      const char *fmt="%s\n");


  void OpenCloseDatasets( HierarchyEntry *Grid,
			  int current_level,
			  int Open,
			  int *NumberOfGrids );

  void LoopNumberOfGrids( HierarchyEntry *Grid,
				     int current_level);
  void LoopIndexGrids( HierarchyEntry *Grid,
		       int current_level,
		       int *current_index );
  void LoopIndexGridsByProcessor(HierarchyEntry *Grid,
				 int current_level,
				 int *current_indexes);
  void LoopCountGridsPerProcessor(HierarchyEntry *Grid,
				  int current_level);
  void LoopBuildMyGridArray(HierarchyEntry *Grid,
			    int current_level,
			    int *current_index);
  float GetParticleMass( HierarchyEntry *Grid );
  void FlagGridCells(HierarchyEntry *Grid);

  float CurrentRedshift;

  //hid_t       HDF5FloatType;

  void HDF5CreateFile( char *name );
  hid_t HDF5OpenFile( char *name );
  void HDF5CloseFile( hid_t file_id );
  hid_t HDF5CreateGroup( hid_t file_id, char *name );
  hid_t HDF5OpenGroup( hid_t file_id, char *name );
  void HDF5CloseGroup( hid_t group_id );
  void HDF5MakeDataset( hid_t loc_id, 
			// char *group_name,
			char *dset_name,
			int rank, hsize_t dims[], 
			float *data,
			char *units = NULL);
  //<dcc> the int version.
  void HDF5MakeDataset( hid_t group_id, 
			char *dset_name,
			int rank, hsize_t dims[], 
			int *data,
			char *units = NULL);
  //</dcc>
  herr_t      h5_status;
  void HDF5ReadDataset( hid_t group_id,
			char *dset_name,
			float *data );
 
  void HDF5WriteStringAttr(hid_t dset_id, char *Alabel, char *String);

  HierarchyEntry *ContainingGrid( HierarchyEntry *Grid, float *point);

  HierarchyEntry **FastFindGridArray;
  int *FastFindDeltaN;
  float *FastFindCellWidth;

  HierarchyEntry **MyGridArray;

  int NumberOfGrids;
  Eint32 *GridsPerProcessor;
  void CountGridsPerProcessor();

 public:
  inline float ReturnCurrentRedshift(void){ return CurrentRedshift;};
};

#endif
