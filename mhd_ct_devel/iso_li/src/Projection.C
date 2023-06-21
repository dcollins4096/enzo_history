#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <hdf5.h>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_PNG  
#include "pngwriter.h"
#endif
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
#ifdef USE_PNG  
#include "ColorTables.h"
#endif
#include "ealInt.h"
#include "ealFloat.h"
#include "AnalysisBaseClass.h"
#include "Projection.h"

void my_exit(int status);

int FindField(int field, int farray[], int numfields);
field_type_int get_field_id(char *field_name);
char *get_field_name(field_type_int field);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);

Projection::Projection( TopGridData *mdata, HierarchyEntry *topgrid, 
			int num_fields, field_type_int *fields, int axis, 
			int max_level , char * Title, FLOAT *left, FLOAT *right):
  
  AnalysisBaseClass( mdata, topgrid, max_level, left, right){

  int i, max_dim[2];
  FLOAT proj_frac[2];

  if( Title != NULL ){
    ProjectionTitle = new char[ strlen(Title) + 1];
    strcpy( ProjectionTitle, Title);
  }else{
    ProjectionTitle = new char[30];
    sprintf(ProjectionTitle,"%s","NoGroupName");
  }


  // projection axis and permutations
  // x -> (y, z)
  // y -> (x, z)
  // z -> (x, y)

  ProjectionAxis = axis;

  switch(ProjectionAxis){
  case 0:
    Dim1 = 1;
    Dim2 = 2;
    break;
  case 1:
    Dim1 = 0;
    Dim2 = 2;
    break;
  case 2:
    Dim1 = 0;
    Dim2 = 1;
    break;
  }

  proj_frac[0] = (AnalysisRegionRightEdge[Dim1] - AnalysisRegionLeftEdge[Dim1])/DomainDelta[Dim1];
  proj_frac[1] = (AnalysisRegionRightEdge[Dim2] - AnalysisRegionLeftEdge[Dim2])/DomainDelta[Dim2];

  //fprintf(stderr,"PJ: proj_frac = %f %f\n",proj_frac[0],proj_frac[1]);


  max_dim[0] = MetaData->TopGridDims[Dim1];
  max_dim[1] = MetaData->TopGridDims[Dim2];

  for(i = 0; i < MaximumAnalysisLevel; i++){
    max_dim[0] *= RefineBy;
    max_dim[1] *= RefineBy;
  }


  ProjectionDims[0] = max( (int)(proj_frac[0]*max_dim[0]), 1);
  ProjectionDims[1] = max( (int)(proj_frac[1]*max_dim[1]), 1);

  //fprintf(stderr,"max_dim[0] = %"ISYM", max_dim[1] = %"ISYM"\n", max_dim[0], max_dim[1]);
  //fprintf(stderr,"proj_frac[0] = %"FSYM", proj_frac[1] = %"FSYM"\n", proj_frac[0], proj_frac[1]);
  //fprintf(stderr,"PD: %d %d\n",ProjectionDims[0], ProjectionDims[1]);

  
  ProjectionDX[0] = (AnalysisRegionRightEdge[Dim1] - AnalysisRegionLeftEdge[Dim1])/ProjectionDims[0];
  ProjectionDX[1] = (AnalysisRegionRightEdge[Dim2] - AnalysisRegionLeftEdge[Dim2])/ProjectionDims[1];

  BufferSize = (ProjectionDims[0])*(ProjectionDims[1]);

  LevelField = new ealInt(BufferSize);

  NumberOfFields = num_fields;

  if(NumberOfFields){
    FieldType = new int[NumberOfFields];
    ProjectionField = new ealFloat*[NumberOfFields];
    
    for(i = 0; i < NumberOfFields; i++){

      FieldType[i] = fields[i];
      ProjectionField[i] = new ealFloat(BufferSize);
    }
  }else{
    FieldType = NULL;
    ProjectionField = NULL;
  }

}

Projection::~Projection(void){

  int i;
  delete ProjectionTitle;
  delete LevelField;
  if(FieldType)
    delete [] FieldType;

  if(ProjectionField){
    for(i = 0; i < NumberOfFields; i++){
      delete ProjectionField[i];
    }
    
    delete[] ProjectionField;
  }
}

int proj_counter = 0;
void Projection::Project(){
  if( debug ){
    fprintf(stderr, "  proj:");
    proj_counter = 0;
  }
  LoopProject(TopGrid, 0);
  int field_index;
  for(field_index = 0; field_index < NumberOfFields; field_index++){
    if( debug )
      fprintf(stderr, "  reduce field %d\n",field_index);
    ProjectionField[field_index]->ReduceSum();
  }
  LevelField->ReduceMax();
}

void Projection::LoopProject( HierarchyEntry *Grid,
			      int current_level){
  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){

      if(Grid->NextGridNextLevel && current_level <= MaximumAnalysisLevel)
	LoopProject(Grid->NextGridNextLevel, current_level + 1);

      if(Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
	ProjectGridParticles(Grid->GridData);

	if(current_level < MaximumAnalysisLevel){
	  FlagGridCells( Grid );
	}else{
	  Grid->GridData->ClearFlaggingField();
	}
#ifdef NOT_YET
	if(LoadGridDataAtStart == FALSE){
	  if(Grid->GridData->ReadGrid(NULL, -1, FALSE, TRUE) == FAIL){
	    fprintf(stderr, "Failed to read in grid data, exiting. Grid information follows.\n");
	    PrintGridInfo(Grid->GridData);
	    my_exit(EXIT_FAILURE);
	  }
	}
#endif //NOT_YET
	if( debug ){
	  fprintf(stderr, "%d\n ",proj_counter++);
	  fflush(stderr);
	}
	ProjectGrid(Grid->GridData, current_level);

	Grid->GridData->DeleteFlaggingField();
#ifdef NOT_YET
	if(LoadGridDataAtStart == FALSE){
	  Grid->GridData->DeleteAllFields();
	}
#endif //NOT_YET
      }
    }
    
    Grid = Grid->NextGridThisLevel;
  }
}

void Projection::ProjectGrid(grid *Grid, int CurrentLevel){

  int field_index, i, j, cell_index, proj_index;

  FLOAT x_cell[MAX_DIMENSION];
  FLOAT proj_delta[MAX_DIMENSION];
  int i_proj[MAX_DIMENSION];
  int i_grid[MAX_DIMENSION];
  int i_start[MAX_DIMENSION], i_stop[MAX_DIMENSION];

  FLOAT left[MAX_DIMENSION];
  FLOAT right[MAX_DIMENSION];

  EnzoArrayFloat *baryon_field;

  EnzoArrayInt *flagging_field = Grid->CreateFieldArrayInt(gFlaggingField);

  proj_delta[Dim1] = ProjectionDX[0];
  proj_delta[Dim2] = ProjectionDX[1];
  proj_delta[ProjectionAxis] = flagging_field->CellSize[ProjectionAxis];
  
  for( i = 0; i < MAX_DIMENSION; i++ ){
    left[i] = max( Grid->GetGridLeftEdge(i), AnalysisRegionLeftEdge[i] );
    right[i] = min( Grid->GetGridRightEdge(i), AnalysisRegionRightEdge[i] );
    i_start[i] = (int)((left[i] - AnalysisRegionLeftEdge[i])/proj_delta[i]);
    i_stop[i] = max( (int)((right[i] - AnalysisRegionLeftEdge[i])/proj_delta[i]), i_start[i]+1);
  }
  for(field_index = 0; field_index < NumberOfFields; field_index++){
    baryon_field = Grid->CreateFieldArrayFloat(FieldType[field_index]);
    if(baryon_field){
      for( i_proj[Dim2] = i_start[Dim2]; i_proj[Dim2] < i_stop[Dim2]; i_proj[Dim2]++){  
	for( i_proj[Dim1] = i_start[Dim1]; i_proj[Dim1] < i_stop[Dim1]; i_proj[Dim1]++){       

	  proj_index = i_proj[Dim2]*(ProjectionDims[0]) + i_proj[Dim1];
      
	  for (i_proj[ProjectionAxis] = i_start[ProjectionAxis]; 
	       i_proj[ProjectionAxis] < i_stop[ProjectionAxis]; 
	       i_proj[ProjectionAxis]++){
	    for( i = 0; i < MAX_DIMENSION; i++ ){
	      x_cell[i] = ( ((float) i_proj[i]) + 0.5 ) * proj_delta[i] + AnalysisRegionLeftEdge[i];
	      i_grid[i] = (int) ((x_cell[i] - Grid->GetGridLeftEdge(i))/(baryon_field->CellSize[i])) + baryon_field->StartIndex[i];
	    }
	    cell_index = i_grid[0] + i_grid[1]*(baryon_field->Dimension[0]) +
	      i_grid[2]*(baryon_field->Dimension[0])*(baryon_field->Dimension[1]);

	    if(!(flagging_field->Array[cell_index])){
	      ProjectionField[field_index]->Array[proj_index] += baryon_field->Array[cell_index]*(proj_delta[ProjectionAxis]);
	  
	      LevelField->Array[proj_index] = max(LevelField->Array[proj_index], CurrentLevel);
 
	    }
	  }
	}
      }
      delete baryon_field;
      baryon_field = NULL;
    }
  }

  delete flagging_field;
}

// Future-proofing, implement when needed
void Projection::ProjectGridParticles(grid *Grid){

}

void Projection::Write(char *outputname){

  if(MyProcessorNumber == ROOT_PROCESSOR){

    int i;

    char basename[MAX_LINE_LENGTH];
    if(outputname == NULL){
      strcpy( basename, MetaData->AnalysisName);
      strcat( basename, ".proj" );
    }else{
      strcpy( basename, outputname);
    }

    char *fieldname;

    FILE *fptr = fopen(basename, "a");

    fprintf(fptr, "ProjectionTitle        = %s\n",ProjectionTitle);
    fprintf(fptr, "ProjectionAxis         = %"ISYM"\n", ProjectionAxis);

    fprintf(fptr, "Dimension              = ");
    WriteListOfInts(fptr, 2, ProjectionDims);  
  
    fprintf(fptr, "LeftEdge               = ");
    WriteListOfFloats(fptr, 3, AnalysisRegionLeftEdge);

    fprintf(fptr, "RightEdge              = ");
    WriteListOfFloats(fptr, 3, AnalysisRegionRightEdge);
    
    fprintf(fptr, "Time                   = %"GOUTSYM"\n", 
	    TopGrid->GridData->ReturnTime());

    strcat( basename, ".h5" );

    HDF5CreateFile(basename);

    hid_t hdf_file = HDF5OpenFile(basename);

    hid_t hdf_group = HDF5CreateGroup(hdf_file,ProjectionTitle);

    hsize_t dims[2];
    dims[0] = ProjectionDims[0];
    dims[1] = ProjectionDims[1];

    for(i = 0; i < NumberOfFields; i++){
      fieldname = get_field_name(FieldType[i]);
      if(fieldname){
	fprintf(fptr, "Field[%"ISYM"]               = %s\n", i, fieldname);
	HDF5MakeDataset( hdf_group, fieldname, 2, dims, ProjectionField[i]->Array);
	delete[] fieldname;
      }
    }

    HDF5MakeDataset( hdf_group, "MaxLevel",2,dims,LevelField->Array);
    HDF5CloseGroup(hdf_group);
    HDF5CloseFile(hdf_file);
    int file_status = fclose(fptr);
    assert(file_status == 0);
  }
}

EnzoArrayFloat *Projection::CreateFieldArrayFloat(char *field_name){

  field_type_int field_id = get_field_id(field_name);

  return this->CreateFieldArrayFloat(field_id);
}

EnzoArrayFloat *Projection::CreateFieldArrayFloat(field_type_int field){

  EnzoArrayFloat *array = NULL;
  int sindex[] = {0, 0};
  int eindex[2];
  eindex[0] = ProjectionDims[0] -1;
  eindex[1] = ProjectionDims[1] -1;

  field_type_int field_index = FindField(field, 
				     this->FieldType, 
				     this->NumberOfFields);  

  if(field_index != -1){
    array = new EnzoArrayFloat(2, ProjectionDims,
			       sindex, eindex,
			       ProjectionDX);    
    array->Array = ProjectionField[field_index]->Array;
  }

  return array;
}


void Projection::Plot(char *outputname, int color_table){
  
    
    Eint32 field_index, index;
    Eint32 i, j, r, g, b, color_index;
    float field_min, field_max, field_delta, field_val;

    char basename[MAX_LINE_LENGTH];
    char imagename[MAX_LINE_LENGTH];
    char *fieldname;

#ifdef USE_PNG
  if(MyProcessorNumber == ROOT_PROCESSOR){

    if(outputname == NULL){
      strcpy( basename, MetaData->AnalysisName);
    }else{
      strcpy( basename, outputname);
    }

    for(field_index = 0; field_index < NumberOfFields; field_index++){
      fieldname = get_field_name(FieldType[field_index]);
      if(fieldname){
	sprintf(imagename, "%s.%s.png", basename, fieldname);
	pngwriter png_image(ProjectionDims[0], ProjectionDims[1], 1.0, imagename);

	field_min = ProjectionField[field_index]->Min();
	field_max = ProjectionField[field_index]->Max();
	
	field_delta = max(field_max - field_min, tiny_number);
	
	for (j = 0; j < ProjectionDims[1]; j++) {
	  index = j*ProjectionDims[0];
	  
	  for (i = 0; i < ProjectionDims[0]; i++, index++) {
	    field_val = ProjectionField[field_index]->Array[index];
	    
	    color_index = nint(255.0*((field_val - field_min)/field_delta));
	    
	    r = 256*RedArrays[color_table][color_index]+color_index;
	    g = 256*BlueArrays[color_table][color_index]+color_index;
	    b = 256*GreenArrays[color_table][color_index]+color_index;
	    
	    png_image.plot(i, j, r, g, b);	     
	  }
	}

	png_image.close();
	
	delete[] fieldname;
      }
    }
  }
#endif
}
