#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"

#include "AnalysisBaseClass.h"

void my_exit(int exit_status);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

void FlagLevel( HierarchyEntry *Grid, HierarchyEntry ***GridArray, int *dx,
		float cell_width );


AnalysisBaseClass::AnalysisBaseClass( TopGridData *metadata,
			    HierarchyEntry *topgrid,
			    int level,
			    FLOAT *Left,
			    FLOAT *Right):
#ifdef USE_BAD_VALUES
  G(6.67e-8),
  MPC_CM(3.0856e+24),
  MSOLAR_G(1.989e+33)
#else
  G(6.6742e-8),
  MPC_CM(3.0856776e+24),
  MSOLAR_G(1.989e+33)
#endif
{

  int i;
  MetaData = metadata;
  TopGrid = topgrid;

  Next = NULL;
  
  /* Create some calculated data to save flops. */
  for( i=0; i < MAX_DIMENSION; i++)
    DomainDelta[i] = DomainRightEdge[i] - DomainLeftEdge[i];  

  if( Left ){
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionLeftEdge[i] = Left[i];
  }else{
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionLeftEdge[i] = DomainLeftEdge[i];
  }

  if( Right ){
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionRightEdge[i] = Right[i];
  }else{
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionRightEdge[i] = DomainRightEdge[i];
  }

  if( level < 0 ){
    MaximumAnalysisLevel = MaximumRefinementLevel;
  }else{
    MaximumAnalysisLevel = level;
  }

  CurrentRedshift = 0.0;
  if (ComovingCoordinates) {
    FLOAT a =1.0, dadt, t = topgrid->GridData->ReturnTime();
    
    if (CosmologyComputeExpansionFactor(t, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in ComputeExpansionFactor.\n");
      my_exit(EXIT_FAILURE);
    }

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    /* Thought this was a good idea, but it isn't 
       if(CurrentRedshift < 0)
       CurrentRedshift = 0.0;
    */
  }
  
  FastFindGridArray = NULL;
  FastFindDeltaN = NULL;
  FastFindCellWidth = NULL;

  MyGridArray = NULL;

  GridsPerProcessor = NULL;

}

AnalysisBaseClass::~AnalysisBaseClass(){

  delete[] FastFindGridArray;
  delete[] FastFindDeltaN;
  delete[] FastFindCellWidth;
  delete[] GridsPerProcessor;
  delete[] MyGridArray;
}

void AnalysisBaseClass::HDF5CreateFile( char *name ){

  hid_t       file_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;

  //Don't create existing files.
  FILE * test = fopen(name,"r");
  if( test ){
    fclose( test );
    return;
  }

  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if( file_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateFile, cannot create file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Fclose(file_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateFile, cannot close file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

}

hid_t AnalysisBaseClass::HDF5OpenFile( char *name ){

  hid_t       file_id;
  herr_t      h5_error = -1;

  file_id = H5Fopen(name, H5F_ACC_RDWR, H5P_DEFAULT);

  if( file_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5OpenFile, cannot open file %s\n", name);
    my_exit(EXIT_FAILURE);
  }

  return file_id;
}

void AnalysisBaseClass::HDF5CloseFile( hid_t file_id ){

  herr_t      h5_status;
  herr_t      h5_error = -1;

  h5_status = H5Fclose(file_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CloseFile, cannot close file\n" );
    my_exit(EXIT_FAILURE);
  }

}

hid_t AnalysisBaseClass::HDF5CreateGroup( hid_t loc_id, char *name ){

  hid_t group_id;
  herr_t h5_status;
  herr_t      h5_error = -1;

  group_id = H5Gcreate(loc_id, name, 0);

  if( group_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CreateGroup, cannot create group %s\n", name);
    my_exit(EXIT_FAILURE);
  }
  if( debug )
    fprintf(stderr,"HDF5CreateGroup: created %s\n",name);
  return group_id;
}


hid_t AnalysisBaseClass::HDF5OpenGroup( hid_t loc_id, char *name ){

  hid_t group_id;
  herr_t h5_status;
  herr_t      h5_error = -1;

  group_id = H5Gopen(loc_id, name);

  if( group_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5GetGroup, cannot open group %s\n", name);
    my_exit(EXIT_FAILURE);
  }
  return group_id;
}

void AnalysisBaseClass::HDF5CloseGroup( hid_t group_id ){

  herr_t      h5_status;
  herr_t      h5_error = -1;

  h5_status = H5Gclose(group_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5CloseGroup, cannot close group\n" );
    my_exit(EXIT_FAILURE);
  }
}

void AnalysisBaseClass::HDF5ReadDataset( hid_t group_id,
				     char *dset_name,
				     float *data ){

  herr_t h5_status;
  herr_t      h5_error = -1;

  hid_t dset_id =  H5Dopen(group_id, dset_name);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot open read dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot read dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  if( debug )
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, successfully wrote %s.\n",dset_name);
  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5ReadDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
}

void AnalysisBaseClass::PrintParameters(){

  //printf("Enzo Analysis class parameters:\n");
  printf("MaximumAnalysisLevel                 = %"ISYM"\n", MaximumAnalysisLevel);
  printf("CurrentRedshift                      = %"FSYM"\n", CurrentRedshift); 
}

void AnalysisBaseClass::PrintGridInfo(grid *Grid){

  int i;
  int grid_rank;
  int grid_dims[MAX_DIMENSION];
  FLOAT grid_left[MAX_DIMENSION];
  FLOAT grid_right[MAX_DIMENSION];

  Grid->ReturnGridInfo(&grid_rank, grid_dims, grid_left, grid_right);

  printf("Grid information:\n");
  printf("  Rank      = %"ISYM"\n", grid_rank);
  printf("  Dims      = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %4"ISYM, grid_dims[i]);
  printf("\n" );

  printf("  LeftEdge  = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %8.6"FSYM, grid_left[i]);
  printf("\n" );

  printf("  RightEdge = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %8.6"FSYM, grid_right[i]);
  printf("\n" );
}

void AnalysisBaseClass::PrintTrueFalse(FILE *buffer, bool testvalue,
				  const char *fmt){

  if(testvalue){
    fprintf(buffer, fmt, "TRUE");
  }else{
    fprintf(buffer, fmt, "FALSE");
  }
}
 
void AnalysisBaseClass::HDF5MakeDataset( hid_t group_id, 
				    //char *group_name,
				    char *dset_name,
				    int rank, hsize_t dims[], 
				    float *data,
				    char *units ){
  herr_t      h5_status;
  herr_t      h5_error = -1;

  // WARN WARN WARN
  // harcoding 32-bit IO
  hid_t float_type_id = HDF5_R4;
  hid_t file_type_id = HDF5_FILE_R4;

  int i, size = 1;
  for(i = 0; i < rank; i++)
    size *= int(dims[i]);

  Eflt32 *buffer = new Eflt32[size];

  for(i = 0; i < size; i++)  
    buffer[i] = Eflt32(data[i]);

  hid_t file_dsp_id = H5Screate_simple(rank, dims, NULL);

  if( file_dsp_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot open dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  hid_t dset_id =  H5Dcreate(group_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot create dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  
  /* Write the dataset. */
  h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot write dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  if(units)
    HDF5WriteStringAttr(dset_id, "Units", units);
  
  h5_status = H5Sclose(file_dsp_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  delete [] buffer;
}


//<dcc> made one for ints.
 
void AnalysisBaseClass::HDF5MakeDataset( hid_t group_id, 
					 //char *group_name,
					 char *dset_name,
					 int rank, hsize_t dims[], 
					 int *data,
					 char *units ){
  herr_t      h5_status;
  herr_t      h5_error = -1;

  // WARN WARN WARN
  // harcoding 32-bit IO
  hid_t float_type_id = HDF5_I4;
  hid_t file_type_id = HDF5_FILE_I4;

  int i, size = 1;
  for(i = 0; i < rank; i++)
    size *= int(dims[i]);

  Eint32 *buffer = new Eint32[size];

  for(i = 0; i < size; i++)  
    buffer[i] = Eint32(data[i]);

  hid_t file_dsp_id = H5Screate_simple(rank, dims, NULL);

  if( file_dsp_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot open dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  hid_t dset_id =  H5Dcreate(group_id, dset_name, file_type_id, file_dsp_id, H5P_DEFAULT);

  if( dset_id == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot create dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }
  
  /* Write the dataset. */
  h5_status = H5Dwrite(dset_id, float_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot write dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  if(units)
    HDF5WriteStringAttr(dset_id, "Units", units);
  
  h5_status = H5Sclose(file_dsp_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataspace\n");
    my_exit(EXIT_FAILURE);
  }

  h5_status = H5Dclose(dset_id);

  if( h5_status == h5_error ){
    fprintf(stderr, "AnalysisBaseClass::HDF5MakeDataset, cannot close dataset %s\n", dset_name);
    my_exit(EXIT_FAILURE);
  }

  delete [] buffer;
}

//</dcc>

void AnalysisBaseClass::HDF5WriteStringAttr(hid_t dset_id, char *Alabel, char *String){

  hid_t       attr_id;
  hid_t       attr_dsp_id;
  hid_t       attr_type_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;

  const char  *NoString = "none";

  attr_dsp_id = H5Screate(H5S_SCALAR);

  if( attr_dsp_id == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to open dataspace\n"); 
    my_exit(EXIT_FAILURE);
  }

  attr_type_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(attr_type_id, 80);

  attr_id = H5Acreate(dset_id, Alabel, attr_type_id,  attr_dsp_id, H5P_DEFAULT);

  if(  attr_id == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to create attribute %s\n", Alabel ); 
    my_exit(EXIT_FAILURE);
  }

  if( strlen(String) > 0 ){
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) String);
    if(  h5_status == h5_error ){ 
      fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to write string %s\n", String); 
      my_exit(EXIT_FAILURE); 
    }
  }else{
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) NoString);
    if(  h5_status == h5_error ){ 
      fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to write string %s\n", NoString); 
      my_exit(EXIT_FAILURE); 
    }
  }

  h5_status = H5Aclose(attr_id);

  if(  h5_status == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to close attribute\n"); 
    my_exit(EXIT_FAILURE); 
  }

  h5_status = H5Sclose(attr_dsp_id);

  if(  h5_status == h5_error ){ 
    fprintf(stderr, "AnalysisBaseClass::HDF5WriteStringAttr, unable to close dataspace\n"); 
    my_exit(EXIT_FAILURE); 
  }
}

void AnalysisBaseClass::LoopNumberOfGrids( HierarchyEntry *Grid,
					   int current_level){
  
  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){
      if(current_level < MaximumAnalysisLevel){
	LoopNumberOfGrids(Grid->NextGridNextLevel, current_level + 1);
      }    
      NumberOfGrids++;
    }
    
    Grid = Grid->NextGridThisLevel;
  }
}

int AnalysisBaseClass::NumberOfGridsInVolume(  ){
  NumberOfGrids = 0;
  LoopNumberOfGrids( TopGrid, 0 );
  return NumberOfGrids;
}

void AnalysisBaseClass::LoopCountGridsPerProcessor( HierarchyEntry *Grid,
						    int current_level){
  
  int proc_index;

  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){      
      if(current_level < MaximumAnalysisLevel){
	LoopCountGridsPerProcessor(Grid->NextGridNextLevel, current_level + 1);
      }
      proc_index = Grid->GridData->ReturnProcessorNumber();
      GridsPerProcessor[proc_index]++;      
    }

    Grid = Grid->NextGridThisLevel;
  }
}

void AnalysisBaseClass::CountGridsPerProcessor(){
  if(!GridsPerProcessor)
    GridsPerProcessor = new Eint32[NumberOfProcessors];
  
  for(int i = 0; i < NumberOfProcessors; i++)
    GridsPerProcessor[i] = 0;

  LoopCountGridsPerProcessor( TopGrid, 0 );
}

void AnalysisBaseClass::LoopBuildMyGridArray(HierarchyEntry *Grid,
					     int current_level,
					     int *current_index){

  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){      
      if(current_level < MaximumAnalysisLevel){
	LoopBuildMyGridArray(Grid->NextGridNextLevel, current_level + 1, 
			     current_index );
      }    
      if(Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
	MyGridArray[(*current_index)++] = Grid;	
      }
    }    
    Grid = Grid->NextGridThisLevel;
  }
}


void AnalysisBaseClass::BuildMyGridArray(  ){
  int index = 0;

  if(GridsPerProcessor == NULL)
    CountGridsPerProcessor();

  delete[] MyGridArray;
  MyGridArray = new HierarchyEntry*[GridsPerProcessor[MyProcessorNumber]];
  
  LoopBuildMyGridArray( TopGrid, 0, &index);

}


void AnalysisBaseClass::FlagGridCells(HierarchyEntry *Grid){
  
  Grid->GridData->ClearFlaggingField();

  HierarchyEntry *SubGrid = Grid->NextGridNextLevel;
  
  while(SubGrid != NULL){
    if(Grid->GridData->FlagRefinedCells(SubGrid->GridData)==FAIL)
      my_exit( EXIT_FAILURE );
    SubGrid = SubGrid->NextGridThisLevel;
  }
}

void AnalysisBaseClass::SetAnalysisRegion( FLOAT Left[MAX_DIMENSION],
			  FLOAT Right[MAX_DIMENSION] ){
  
  for( int i = 0; i < MAX_DIMENSION; i++ ){
    AnalysisRegionLeftEdge[i] = Left[i];
    AnalysisRegionRightEdge[i] = Right[i];
  }
}

void AnalysisBaseClass::InitializeFastFind(){
  int i, j, num_cells = 1;

  HierarchyEntry *this_grid, *next_level;
  
  FastFindDeltaN = new int[MAX_DIMENSION];
  FastFindCellWidth = new float[MAX_DIMENSION];
  
  // hard coded to 16 for unigrid case
  
  for( i=0; i < MAX_DIMENSION; i++){
    if(MaximumAnalysisLevel == 0){
      num_cells*= 16;
      FastFindCellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i])/16.0;      
      FastFindDeltaN[i] = 16;
    }else{
      num_cells*= MetaData->TopGridDims[i];
      FastFindCellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i])/((float)(MetaData->TopGridDims[i]));
      FastFindDeltaN[i] = MetaData->TopGridDims[i];
    }
  }

  FastFindGridArray = new HierarchyEntry*[num_cells];

  this_grid = TopGrid;

  // hard coded one level, since that's all we can resolve with topgriddims

  while( this_grid ){
    this_grid->GridData->FlagGridArray( &FastFindGridArray, FastFindDeltaN, FastFindCellWidth, this_grid );
    if(MaximumAnalysisLevel > 0){
      next_level = this_grid->NextGridNextLevel;
      while( next_level ){
	next_level->GridData->FlagGridArray( &FastFindGridArray, FastFindDeltaN, FastFindCellWidth, next_level );
	next_level = next_level->NextGridThisLevel;
      }
    }
    this_grid = this_grid->NextGridThisLevel;
  }

  return;
}


HierarchyEntry *AnalysisBaseClass::ContainingGrid( HierarchyEntry *Grid, float *point){

  while( Grid != NULL ){
    if( Grid->GridData->PointInGrid( point ))
      break;
    Grid = Grid->NextGridThisLevel;
  }

  return Grid;
}

HierarchyEntry *AnalysisBaseClass::FindGrid( HierarchyEntry *Grid, float *point){

  HierarchyEntry *NextGrid;
  Grid = ContainingGrid( Grid, point );

  if( Grid )
    while( Grid->NextGridNextLevel ){
      NextGrid = ContainingGrid( Grid->NextGridNextLevel, point );
      if( NextGrid )
	Grid = NextGrid;
      else
	break;
    }
  else
    my_exit(EXIT_FAILURE);
  
  return Grid;
}


HierarchyEntry *AnalysisBaseClass::FastFindGrid(float *point){

  int i, i_l[MAX_DIMENSION];
  HierarchyEntry *NextGrid;
  HierarchyEntry *Grid;

  if(!FastFindGridArray)
    InitializeFastFind();

  for(i=0; i < MAX_DIMENSION; i++){
    i_l[i] = int((point[i] - DomainLeftEdge[i])/((float)(FastFindCellWidth[i])));
  }

  Grid = FastFindGridArray[i_l[0] + i_l[1]*FastFindDeltaN[0] +
				   i_l[2]*FastFindDeltaN[0]*FastFindDeltaN[1]];

  if(MaximumAnalysisLevel > 1){
    
    if( Grid ){
      while( Grid->NextGridNextLevel ){
	NextGrid = ContainingGrid( Grid->NextGridNextLevel, point );
	if( NextGrid ){
	  Grid = NextGrid;
	}else{
	  break;
	}
      }
    }else{
      my_exit(EXIT_FAILURE);
    }
  }

  return Grid;
}
