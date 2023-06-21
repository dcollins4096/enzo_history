#ifdef MHDF

//
// Allocate, Deallocte Magnetic and Electric Quantities.
// This is probably a temporary routine: this functionally belongs
// in Grid_AllocateGrids.C
//

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


void grid::MHD_Allocate(){
  int field, dim;
  
  for(field=0; field<3; field++){
    MagneticSize[field] = 1;
    ElectricSize[field] = 1;
    
    for(dim=0; dim<3; dim++){
      MagneticDims[field][dim] = GridDimension[dim];
      ElectricDims[field][dim] = GridDimension[dim] +1;
      
      MHDStartIndex[field][dim] = GridStartIndex[dim];
      MHDEndIndex[field][dim] = GridEndIndex[dim];
      
      MHDeStartIndex[field][dim] = GridStartIndex[dim];
      MHDeEndIndex[field][dim] = GridEndIndex[dim]+1;
      
      if( field == dim )
	{
	  MagneticDims[field][dim]++;
	  ElectricDims[field][dim]--;
	  MHDEndIndex[field][dim]++;
	  MHDeEndIndex[field][dim]--;
	  
	}
      
      MagneticSize[field] *= MagneticDims[field][dim];
      ElectricSize[field] *= ElectricDims[field][dim];
    }
  }

  int i;
  for( field=0;field<NFaces;field++){
    if( MagneticField[field] != NULL ){
      fprintf(stderr,"Grid_MHD_Allocate: Non Null Magnetic Field\n");
    }else{
      MagneticField[field] = new float[ MagneticSize[field] ];
      //fprintf(stderr," Make B %d %p\n", field,  MagneticField[field]);
      //for( i=0;i<MagneticSize[field]; i++)
      //MagneticField[field][i] = 0;
    }
  }
  for( field=0; field<NEdges;field++){
    if( ElectricField[field] != NULL ){
      fprintf(stderr,"Grid_MHD_Allocate: Non Null Electric Field\n");      
    }else{
      ElectricField[field] = new float[ ElectricSize[field] ];
      //fprintf(stderr," Make E %d %p\n", field,  ElectricField[field]);
      for( i=0; i<ElectricSize[field]; i++)
	ElectricField[field][i] = 0;
    }
  }
  //<dbg>
  /*
    int index;
    
    for(field=0; field<3; field++)
    for( int k=0;k<MagneticDims[field][2]; k++)
    for( int j=0;j<MagneticDims[field][1]; j++)
    for( int i=0;i<MagneticDims[field][0]; i++){
    index = i + MagneticDims[field][0] * (j + MagneticDims[field][1] * k );
    MagneticField[ field ][index] = 0.0;
    }
  */
  //</dbg>
  fprintf(stderr,"Allocate some stuff %d %d\n", MagneticSize[0], ElectricSize[0]);
}
void grid::MHD_Deallocate(){

  int field;
  fprintf(stderr,"Delete some stuff\n");
  for( field=0;field<NFaces;field++)
    if( MagneticField[field] != NULL ){
      delete [] MagneticField[field];
      MagneticField[field]=NULL;
      //fprintf(stderr,"Null B\n");
    }
  for( field=0; field<NEdges;field++)
    if( ElectricField[field] != NULL ){
      delete [] ElectricField[field];
      ElectricField[field] = NULL;
      //fprintf(stderr,"Null E\n");
    }
}

#endif //MHDF
