/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (INITIALIZE BOUNDARY FACE TO A CONSTANT VALUE)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <hdf5.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>


#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

// Set external boundaries to a face-constant value.
//   (Note: we do not need BoundaryValue if the boundary type is inflow)

#ifdef OOC_BOUNDARY

int WRITE_BT(boundary_type *bt_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);
int WRITE_BV(float         *bv_buffer, int field, int dim, int face, int slabsize, int BoundaryDimension[], int BoundaryRank, int Nfields);

#endif

int ExternalBoundary::InitializeExternalBoundaryFace(int dim, 
					     boundary_type LeftBoundaryType,
					     boundary_type RightBoundaryType,
					     float LeftBoundaryValue[],
					     float RightBoundaryValue[])
{
  int field, index;
  
#ifdef OOC_BOUNDARY
  hid_t       file_id, dset_id, attr_id;
  hid_t       file_dsp_id, mem_dsp_id, attr_dsp_id;
  hid_t       file_type_id, mem_type_id;
  hid_t       int_file_type_id, int_mem_type_id;
  
  hsize_t     mem_stride, mem_count, mem_block;
  hsize_t     file_stride[4], file_count[4], file_block[4];
  
  hssize_t    mem_offset;
  hssize_t    file_offset[4];
  
  herr_t      h5_status;
  herr_t      h5_error = -1;
  
  int slabsize;
  float *bv_buffer;
  boundary_type *bt_buffer;
#endif
  
  /* Error check */
  
  if (dim > BoundaryRank) {
    fprintf(stderr, "Dimension %d > BoundaryRank %d.\n", dim, BoundaryRank);
    return FAIL;
  }
  
  /* compute size of entire mesh */
  
  int size = 1;
  for (int i = 0; i < BoundaryRank; i++)
    size = size*BoundaryDimension[i];
  
  /* set BoundaryType faces to a constant */
  
  if (BoundaryDimension[dim] != 1){
    
#ifdef OOC_BOUNDARY
    
    slabsize = size/BoundaryDimension[dim];
    bt_buffer  = new boundary_type[size/BoundaryDimension[dim]]; 
    
    for (field = 0; field < NumberOfBaryonFields; field++) {
      
      for (index = 0; index < size/BoundaryDimension[dim]; index++) {
        bt_buffer[index] = LeftBoundaryType;
      }
      WRITE_BT(bt_buffer, field, dim, 0, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      
      for (index = 0; index < size/BoundaryDimension[dim]; index++) {
        bt_buffer[index] = RightBoundaryType;
      }
      WRITE_BT(bt_buffer, field, dim, 1, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
      
    } // end of loop over fields
    
    ExternalBoundaryTypeIO = TRUE;
    
    delete [] bt_buffer;
    
#else

    for (field = 0; field < NumberOfBaryonFields; field++) {
      BoundaryType[field][dim][0] = 
	new boundary_type[size/BoundaryDimension[dim]];
      BoundaryType[field][dim][1] = 
	new boundary_type[size/BoundaryDimension[dim]];
      for (index = 0; index < size/BoundaryDimension[dim]; index++) {
	BoundaryType[field][dim][0][index] = LeftBoundaryType;
	BoundaryType[field][dim][1][index] = RightBoundaryType;
      }
    }
    //ExternalBoundaryTypeIO = FALSE; //This line is in the 'other' copy, but I don't understand why.
#endif
  }//BoundaryDim !=1
  /* Set Magnetic Boundary conditions */
  if(MHD_Used)
    for(field=0;field<3;field++){
      MagneticBoundaryType[field][dim][0]=LeftBoundaryType;
      MagneticBoundaryType[field][dim][1]=RightBoundaryType;
    }
  
  
  /* Set MHD Inflow BC's here.  To be done later. */
  
  /* If required, set BoundaryType faces to a constant (usually inflow) */
  
  if (BoundaryDimension[dim] != 1){
    
#ifdef OOC_BOUNDARY
    
    slabsize = size/BoundaryDimension[dim];
    bv_buffer = new float[size/BoundaryDimension[dim]];
    
    for (field = 0; field < NumberOfBaryonFields; field++) {
      
      if (LeftBoundaryType == inflow) {
        for (index = 0; index < size/BoundaryDimension[dim]; index++) {
          bv_buffer[index] = LeftBoundaryValue[field];
        }
        WRITE_BV(bv_buffer, field, dim, 0, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
        ExternalBoundaryValueIO = TRUE;
      }
      
      if (RightBoundaryType == inflow) {
        for (index = 0; index < size/BoundaryDimension[dim]; index++) {
          bv_buffer[index] = RightBoundaryValue[field];
        }
        WRITE_BV(bv_buffer, field, dim, 1, slabsize, BoundaryDimension, BoundaryRank, NumberOfBaryonFields);
        ExternalBoundaryValueIO = TRUE;
      }
      
    } // end of loop over fields
    
    delete [] bv_buffer;
    
#else
    
    for (field = 0; field < NumberOfBaryonFields; field++) {
      
      if (LeftBoundaryType == inflow) {
	BoundaryValue[field][dim][0] = new float[size/BoundaryDimension[dim]];
	for (index = 0; index < size/BoundaryDimension[dim]; index++)
	  BoundaryValue[field][dim][0][index] = LeftBoundaryValue[field];
      }
      
      if (RightBoundaryType == inflow) {
	BoundaryValue[field][dim][1] = new float[size/BoundaryDimension[dim]];
	for (index = 0; index < size/BoundaryDimension[dim]; index++)
	  BoundaryValue[field][dim][1][index] = RightBoundaryValue[field];
      }
      
    } // end of loop over fields
    
    //ExternalBoundaryValueIO = FALSE; //See comment above.
    
#endif
    
  } // end of for loop over fields
  
  return SUCCESS;
  
}
