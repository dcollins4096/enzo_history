/***********************************************************************
/
/  GRID CLASS (READ GRID)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness, July 2002
/  modified2:  Alexei Kritsuk, Jan 2004   a trick for RandomForcing //AK
/  modified3:  Robert Harkness, Jan 2007 for HDF5 memory buffering
/  modified4:  Robert Harkness, April 2008
/  modified5:  Matthew Turk, September 2009 for refactoring and removing IO_TYPE
/
/  PURPOSE:
/
************************************************************************/
 
//  Input a grid from file pointer fpt
 
#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "h5utilities.h"
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
void my_exit(int status);
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif
 
// HDF5 function prototypes
 

 
// function prototypes
 
int ReadListOfFloats(FILE *fptr, int N, FLOAT floats[]);
int ReadListOfInts(FILE *fptr, int N, int nums[]);
 
static int GridReadDataGridCounter = 0;
 
 
int grid::Group_ReadGrid(FILE *fptr, int GridID, HDF5_hid_t file_id, 
			 int ReadText, int ReadData, int ReadEverything)
{
 
  int i, j, k, field, size, active_size, dim;
  char name[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
  char logname[MAX_LINE_LENGTH];
  char procfilename[MAX_LINE_LENGTH];
 
  char id[MAX_GROUP_TAG_SIZE];
  char pid[MAX_TASK_TAG_SIZE];
  char gpid[MAX_TASK_TAG_SIZE];
 
  int ActiveDim[MAX_DIMENSION];
 
  FILE *log_fptr;
 
  hid_t       group_id, dset_id, old_fields;
  hid_t       file_dsp_id;
  hid_t       num_type;
 
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[MAX_DIMENSION];
 
  herr_t      h5_status;
  herr_t      h5_error = -1;
 
  int         num_size;
 
  char *ParticlePositionLabel[] =
    {"particle_position_x", "particle_position_y", "particle_position_z"};
  char *ParticleVelocityLabel[] =
    {"particle_velocity_x", "particle_velocity_y", "particle_velocity_z"};
  char *ParticleAttributeLabel[] = {"creation_time", "dynamical_time",
                                    "metallicity_fraction", "alpha_fraction"};
 
  if(ReadText){

    /* Read general grid class data */

    /* make sure quantities defined at least for 3d */
 
    for (int dim = GridRank; dim < 3; dim++) {
      GridDimension[dim] = 1;
      GridStartIndex[dim] = 0;
      GridEndIndex[dim] = 0;
    }
    if (fscanf(fptr, "GridRank = %"ISYM"\n", &GridRank) != 1) {
            ENZO_FAIL("Error reading GridRank.");
    }
 
    if (fscanf(fptr, "GridDimension = ") != 0) {
            ENZO_FAIL("Error reading GridDimension(0).");
    }
 
    if (ReadListOfInts(fptr, GridRank, GridDimension) == FAIL) {
            ENZO_FAIL("Error reading GridDimension(1).");
    }
 
    fscanf(fptr, "GridStartIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridStartIndex) == FAIL) {
            ENZO_FAIL("Error reading GridStartIndex.");
    }
 
    fscanf(fptr, "GridEndIndex = ");
 
    if (ReadListOfInts(fptr, GridRank, GridEndIndex) == FAIL) {
            ENZO_FAIL("Error reading GridEndIndex.");
    }
 
    fscanf(fptr, "GridLeftEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridLeftEdge) == FAIL) {
            ENZO_FAIL("Error reading GridLeftEdge.");
    }
 
    fscanf(fptr, "GridRightEdge = ");
 
    if (ReadListOfFloats(fptr, GridRank, GridRightEdge) == FAIL) {
            ENZO_FAIL("Error reading GridRightEdge.");
    }
 
    if (fscanf(fptr, "Time = %"PSYM"\n", &Time) != 1) {
            ENZO_FAIL("Error reading Time.");
    }
 
    if (fscanf(fptr, "SubgridsAreStatic = %"ISYM"\n", &SubgridsAreStatic) != 1) {
            ENZO_FAIL("Error reading SubgridsAreStatic.");
    }
 
    /* Read baryon field quantities. */
 
    if (fscanf(fptr, "NumberOfBaryonFields = %"ISYM"\n",
	       &NumberOfBaryonFields) != 1) {
            ENZO_FAIL("Error reading NumberOfBaryonFields.");
    }
    if (NumberOfBaryonFields > 0) {
 
      fscanf(fptr, "FieldType = ");
 
      if (ReadListOfInts(fptr, NumberOfBaryonFields, FieldType) == FAIL) {
		ENZO_FAIL("Error reading FieldType.");
      }
 
      fgetpos(fptr, &BaryonFileNamePosition); //AK
 
      if (fscanf(fptr, "BaryonFileName = %s\n", procfilename) != 1) {
		ENZO_FAIL("Error reading BaryonFileName.");
      }
 
      fscanf(fptr, "CourantSafetyNumber    = %"FSYM"\n", &CourantSafetyNumber);
      fscanf(fptr, "PPMFlatteningParameter = %"ISYM"\n", &PPMFlatteningParameter);
      fscanf(fptr, "PPMDiffusionParameter  = %"ISYM"\n", &PPMDiffusionParameter);
      fscanf(fptr, "PPMSteepeningParameter = %"ISYM"\n", &PPMSteepeningParameter);
    }

    /* 3) Read particle info */
 
    if (fscanf(fptr, "NumberOfParticles = %"ISYM"\n", &NumberOfParticles) != 1) {
            ENZO_FAIL("error reading NumberOfParticles.");
    }
 
    if (NumberOfParticles > 0) {
 
      /* Read particle file name. */
    
      if (fscanf(fptr, "ParticleFileName = %s\n", procfilename) != 1) {
		ENZO_FAIL("Error reading ParticleFileName.");
      }
    }
 
    /* 4) Read gravity info */
 
    if (SelfGravity)
      if (fscanf(fptr, "GravityBoundaryType = %"ISYM"\n",&GravityBoundaryType) != 1) {
		ENZO_FAIL("Error reading GravityBoundaryType.");
      }

    // If HierarchyFile has different Ghostzones (which should be a parameter not a macro ...)
    // (useful in a restart with different hydro/mhd solvers) 
    int ghosts =DEFAULT_GHOST_ZONES;
    if (GridStartIndex[0] != ghosts)  {
	if (GridID < 2)
     fprintf(stderr,"Grid_Group_ReadGrid: Adjusting Ghostzones which in the hierarchy file did not match the selected HydroMethod.\n");
      for (int dim=0; dim < GridRank; dim++) {
	GridDimension[dim]  = GridEndIndex[dim]-GridStartIndex[dim]+1+2*ghosts;
	GridStartIndex[dim] = ghosts;
	GridEndIndex[dim]   = GridStartIndex[dim]+GridDimension[dim]-1-2*ghosts;
	 if (GridID < 2) fprintf(stderr, "dim: GridStart,GridEnd,GridDim:  %i: %i %i %i\n",
				  dim, GridStartIndex[dim], GridEndIndex[dim], GridDimension[dim]);
      }
    }

  } // (if (ReadText) )

  snprintf(name, MAX_LINE_LENGTH-1, "/Grid%"GROUP_TAG_FORMAT""ISYM, GridID);

  if (NumberOfBaryonFields > 0 && ReadData &&
      (MyProcessorNumber == ProcessorNumber)) {

#ifndef SINGLE_HDF5_OPEN_ON_INPUT
    file_id = H5Fopen(procfilename,  H5F_ACC_RDONLY, H5P_DEFAULT);
    if( file_id == h5_error )ENZO_VFAIL("Error opening %s", procfilename)
#endif
 
    group_id = H5Gopen(file_id, name);
    if( group_id == h5_error )ENZO_VFAIL("Error opening group %s", name)
 
    /* fill in ActiveDim for dims up to 3d */
 
    for (int dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;
 
    /* check dimensions of HDF file against this grid
       (note: we don't bother to check the coordinate arrays)  */
 
    size = 1;
    active_size = 1;
 
    for (int dim = 0; dim < GridRank; dim++) {
      size *= GridDimension[dim];
      active_size *= ActiveDim[dim];
    }
 
    //  CAUTION - are the coordinates reversed?
 
    for (int dim = 0; dim < GridRank; dim++) {
      OutDims[GridRank-dim-1] = ActiveDim[dim];
    }
 
    /* allocate temporary space */
 
    float *temp = new float[active_size];
 
    if(ReadEverything == TRUE) {
      old_fields = H5Gopen(group_id, "OldFields");
    }
 
    /* loop over fields, reading each one */

    for (field = 0; field < NumberOfBaryonFields; field++) {

      BaryonField[field] = new float[size];
      for (i = 0; i < size; i++)
        BaryonField[field][i] = 0;

      this->read_dataset(GridRank, OutDims, DataLabel[field],
          group_id, HDF5_REAL, (VOIDP) temp,
          TRUE, BaryonField[field], ActiveDim);

      if(ReadEverything == TRUE) {
        OldBaryonField[field] = new float[size];
        for (i = 0; i < size; i++)
          OldBaryonField[field][i] = 0;

        this->read_dataset(GridRank, OutDims, DataLabel[field],
            old_fields, HDF5_REAL, (VOIDP) temp,
            TRUE, OldBaryonField[field], ActiveDim);

      }
 
    } // end: loop over fields
 

    if (HydroMethod == MHD_RK) { // This is the MHD with Dedner divergence cleaning that needs an extra field
      // 

   
      int activesize = 1;
      for (int dim = 0; dim < GridRank; dim++)
	activesize *= (GridDimension[dim]-2*DEFAULT_GHOST_ZONES);
      
      if (divB == NULL) 
	divB = new float[activesize];
      
      /* if we restart from a different solvers output without a Phi Field create here and set to zero */
    int PhiNum; 
    if ((PhiNum = FindField(PhiField, FieldType, NumberOfBaryonFields)) < 0) {
      fprintf(stderr, "Starting with Dedner MHD method with no Phi field. \n");
      fprintf(stderr, "Adding it in Grid_ReadGrid.C \n");
      char *PhiName = "Phi";
      PhiNum = NumberOfBaryonFields;
      int PhiToAdd = PhiField;
      this->AddFields(&PhiToAdd, 1);
      DataLabel[PhiNum] = PhiName;
    } else { 
      if (0) 
	for (int n = 0; n < size; n++)
	  BaryonField[PhiNum][n] = 0.;
    }

      for (int dim = 0; dim < 3; dim++)
	if (gradPhi[dim] == NULL)
	  gradPhi[dim] = new float[activesize];
      
      for (int dim = GridRank; dim < 3; dim++)
	for (int n = 0; n < activesize; n++)
	  gradPhi[dim][n] = 0.0;
      
    } /* if HydroMethod == MHD */


    delete [] temp;
 
  }  // end:   if (NumberOfBaryonFields > 0 && ReadData &&
  //      (MyProcessorNumber == ProcessorNumber)) {

  /* Compute Flux quantities */

  this->PrepareGridDerivedQuantities();
 
 
  if (NumberOfParticles > 0 && ReadData &&
      (MyProcessorNumber == ProcessorNumber)) {
  
 
    /* Open file if not already done (note: particle name must = grid name). */
 
    if (NumberOfBaryonFields == 0) {
 
#ifndef SINGLE_HDF5_OPEN_ON_INPUT 
      file_id = H5Fopen(procfilename, H5F_ACC_RDONLY, H5P_DEFAULT);
      if( file_id == h5_error )ENZO_VFAIL("Error opening file %s", name)
#endif
 
      group_id = H5Gopen(file_id, name);
      if( group_id == h5_error )ENZO_VFAIL("Error opening group %s", name)
 
    } // end: if (NumberOfBaryonFields == 0)
 
    /* Allocate room for particles. */
 
    this->AllocateNewParticles(NumberOfParticles);
 
    TempIntArray[0] = NumberOfParticles;
 
    FLOAT *temp = new FLOAT[NumberOfParticles];
 
    /* Read ParticlePosition (use temporary buffer). */
 
    for (int dim = 0; dim < GridRank; dim++) {
      this->read_dataset(1, TempIntArray, ParticlePositionLabel[dim],
            group_id, HDF5_FILE_PREC, (VOIDP) ParticlePosition[dim], FALSE);
    }
 
    /* Read ParticleVelocity. */
 
    for (int dim = 0; dim < GridRank; dim++) {
      this->read_dataset(1, TempIntArray, ParticleVelocityLabel[dim],
            group_id, HDF5_REAL, (VOIDP) ParticleVelocity[dim], FALSE);
    }
 
    this->read_dataset(1, TempIntArray, "particle_mass",
          group_id, HDF5_REAL, (VOIDP) ParticleMass, FALSE);

    /* Read ParticleNumber into temporary buffer and Copy to ParticleNumber. */
 
    this->read_dataset(1, TempIntArray, "particle_index",
          group_id, HDF5_INT, (VOIDP) ParticleNumber, FALSE);

    // Read ParticleType if present
 
    if (ParticleTypeInFile == TRUE) {
 
      /* Read ParticleType into temporary buffer and Copy to ParticleType. */
      this->read_dataset(1, TempIntArray, "particle_type",
            group_id, HDF5_INT, (VOIDP) ParticleType, FALSE);
 
      int abs_type;
      for (i = 0; i < NumberOfParticles; i++) {
	abs_type = ABS(ParticleType[i]);
        if (abs_type < PARTICLE_TYPE_GAS ||
            abs_type > NUM_PARTICLE_TYPES-1) {
          ENZO_VFAIL("file: %s: particle %"ISYM" has unknown type %"ISYM"\n", name, i, ParticleType[i])
        }
      }
 
    } else {
 
      /* Otherwise create the type. */
 
      for (i = 0; i < NumberOfParticles; i++)
        ParticleType[i] = ReturnParticleType(i);
 
    }
 
 
    /* Read ParticleAttributes. */
    if (AddParticleAttributes) {
      for (j = 0; j < NumberOfParticleAttributes; j++) {
	ParticleAttribute[j] = new float[NumberOfParticles];
	for (i=0; i < NumberOfParticles; i++)
	  ParticleAttribute[j][i] = 0;
      }
    } else {
    for (j = 0; j < NumberOfParticleAttributes; j++) {
 
      this->read_dataset(1, TempIntArray, ParticleAttributeLabel[j],
            group_id, HDF5_REAL, (VOIDP) ParticleAttribute[j], FALSE);

    }
    } // ENDELSE add particle attributes
 
    delete [] temp;
 

  } // end: if (NumberOfParticles > 0) && ReadData && (MyProcessorNumber == ProcessorNumber)
 
  /* Close file. */
 
  if ( (MyProcessorNumber == ProcessorNumber) &&
       (NumberOfParticles > 0 || NumberOfBaryonFields > 0)
       && ReadData ){
 
    if (ReadEverything == TRUE) {
      hid_t acc_node;
      H5E_BEGIN_TRY{
        acc_node = H5Gopen(group_id, "Acceleration");
      }H5E_END_TRY
      /* We just check for existence, because for SOME REASON grids don't
         know their own level. */
      if(acc_node != h5_error){
        char acc_name[255];
        size = 1;
        for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];
        float *temp = new float[size];
        for (dim = 0; dim < GridRank; dim++) {
          if(this->AccelerationField[dim] != NULL) {
            delete this->AccelerationField[dim];
          }
          snprintf(acc_name, 254, "AccelerationField%d", dim);
          this->read_dataset(GridRank, OutDims, acc_name,
              acc_node, HDF5_REAL, (VOIDP) temp,
              TRUE, AccelerationField[dim], ActiveDim);
        }
        delete temp;
        H5Gclose(acc_node);
      }

      this->ReadAllFluxes(group_id);
    }
    h5_status = H5Gclose(group_id);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#ifndef SINGLE_HDF5_OPEN_ON_INPUT 

    h5_status = H5Fclose(file_id);
    if( h5_status == h5_error ){ENZO_FAIL("Error in IO");}

#endif
  }
 
  return SUCCESS;
 
}

int grid::read_dataset(int ndims, hsize_t *dims, char *name, hid_t group,
                  hid_t data_type, void *read_to, int copy_back_active,
                  float *copy_to, int *active_dims)
{
  hid_t file_dsp_id;
  hid_t dset_id;
  hid_t h5_status;
  herr_t      h5_error = -1;
  int i, j, k, dim;
  /* get data into temporary array */

  file_dsp_id = H5Screate_simple((Eint32) ndims, dims, NULL);
  if( file_dsp_id == h5_error ){ENZO_FAIL("Error creating file dataspace");}

  dset_id =  H5Dopen(group, name);
  if( dset_id == h5_error )ENZO_VFAIL("Error opening %s", name)

  h5_status = H5Dread(dset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) read_to);
  if( dset_id == h5_error )ENZO_VFAIL("Error reading %s", name)

  h5_status = H5Sclose(file_dsp_id);
  if( dset_id == h5_error )ENZO_VFAIL("Error closing dataspace %s", name)

  h5_status = H5Dclose(dset_id);
  if( dset_id == h5_error )ENZO_VFAIL("Error closing %s", name)

  if(copy_back_active == TRUE) {
    /* copy active region into whole grid */

    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++)
          copy_to[i + j*GridDimension[0] +
            k*GridDimension[0]*GridDimension[1]] =
	      ((float *)read_to)[(i-GridStartIndex[0])                             +
	                         (j-GridStartIndex[1])*active_dims[0]              +
	                         (k-GridStartIndex[2])*active_dims[0]*active_dims[1] ];

  }
  return SUCCESS;
}

int grid::ReadAllFluxes(hid_t grid_node)
{
  /* We get the attribute describing to us the number of subgrids. */

  int i;
  hid_t flux_group, subgrid_group;
  hid_t h5_error = -1;
  char name[255];

  readAttribute(grid_node, HDF5_INT, "NumberOfSubgrids", 
            (void *) &this->NumberOfSubgrids, 1);

  /* Now for every subgrid, we read a flux group, and all of its associated
     baryon fields. */

  fprintf(stderr, "Received NumberOfSubgrids = %d\n", this->NumberOfSubgrids);

  this->SubgridFluxStorage = new fluxes*[this->NumberOfSubgrids];

  flux_group = H5Gopen(grid_node, "Fluxes");
  if(flux_group == h5_error) ENZO_FAIL("Can't open Fluxes group");

  for(i = 0; i < this->NumberOfSubgrids; i++) {
    snprintf(name, 254, "Subgrid%08d", i);
    subgrid_group = H5Gopen(flux_group, name);
    if(subgrid_group == h5_error)ENZO_VFAIL("IO Problem opening %s", name)

      this->SubgridFluxStorage[i] = new fluxes;
    this->ReadFluxGroup(subgrid_group, this->SubgridFluxStorage[i]);
    H5Gclose(subgrid_group);
  }
  subgrid_group = H5Gopen(flux_group, "BoundaryFluxes");
  this->BoundaryFluxes = new fluxes;
  this->ReadFluxGroup(subgrid_group, this->BoundaryFluxes);

  H5Gclose(subgrid_group);
  H5Gclose(flux_group);

}

int grid::ReadFluxGroup(hid_t flux_group, fluxes *fluxgroup)
{
  hid_t h5_error = -1;
  hid_t axis_group = h5_error;
  hid_t left_group, right_group;
  int i, j, field, dim;
  hsize_t size;

  char name[255];

  for (dim = 0; dim < GridRank; dim++) {
    /* compute size (in floats) of flux storage */

    snprintf(name, 254, "Axis%d", dim);
    axis_group = H5Gopen(flux_group, name);
    if(axis_group == h5_error)ENZO_VFAIL("Can't open %s", name)

    size = 1;

    left_group = H5Gopen(axis_group, "Left");
    if(left_group == h5_error){ENZO_FAIL("IO Problem with Left");}

    right_group = H5Gopen(axis_group, "Right");
    if(right_group == h5_error){ENZO_FAIL("IO Problem with Right");}

    readAttribute(left_group, HDF5_INT, "StartIndex",
        fluxgroup->LeftFluxStartGlobalIndex[dim], TRUE);
    readAttribute(left_group, HDF5_INT, "EndIndex",
        fluxgroup->LeftFluxEndGlobalIndex[dim], TRUE);

    readAttribute(right_group, HDF5_INT, "StartIndex",
        fluxgroup->RightFluxStartGlobalIndex[dim], TRUE);
    readAttribute(right_group, HDF5_INT, "EndIndex",
        fluxgroup->RightFluxEndGlobalIndex[dim], TRUE);

    for (j = 0; j < GridRank; j++) {
      size *= fluxgroup->LeftFluxEndGlobalIndex[dim][j] -
        fluxgroup->LeftFluxStartGlobalIndex[dim][j] + 1;
    }

    for (field = 0; field < NumberOfBaryonFields; field++) {
      /* For now our use case ensures these will always exist forever
         and if they don't, we need a hard failure. */
      /* Note also that if you pass a pre-initialized fluxgroup, this will leak
         memory. */
      fluxgroup->LeftFluxes[field][dim]  = new float[size];
      fluxgroup->RightFluxes[field][dim]  = new float[size];

      this->read_dataset(1, &size, DataLabel[field], left_group,
          HDF5_REAL, (void *) fluxgroup->LeftFluxes[field][dim],
          FALSE);
    
      this->read_dataset(1, &size, DataLabel[field], right_group,
          HDF5_REAL, (void *) fluxgroup->RightFluxes[field][dim],
          FALSE);

    }
	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	     field++) {
          fluxgroup->LeftFluxes[field][dim] = NULL;
          fluxgroup->RightFluxes[field][dim] = NULL;
	}

    H5Gclose(left_group);
    H5Gclose(right_group);
    H5Gclose(axis_group);
  }
}
