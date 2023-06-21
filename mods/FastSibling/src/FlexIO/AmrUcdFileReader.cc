#include <stdio.h>
#include <stdlib.h>
#include "AmrUcdFileReader.hh"

void AmrUcdFileReader::getUcd(FlexArray<AmrNode*> &nodes, 
			      FlexArray<int> &cells){
  genUCD.purge(); // Auto-Purges AmrNode's
  for(int i=0;i<activeGrids.getSize();i++){
    int lev = grids[activeGrids[i]].level;
    if(levelmask[lev]){
      printf("Add Grid %u to the hierarchy [%u]\n",activeGrids[i],i);
      // genUCD.addGrid(grid->level, grid->origin, grid->dx, grid->rank,
      // grid->dims, (float *)(grid->data));
      genUCD.addGrid(grids[activeGrids[i]]);
    }
    else {
      //printf("Grid is active, but level mask is OFF for grid [%u] level %u\n",
      //	     activeGrids[i],lev);
    }
  }
  printf("****buildNodeHierarchy\n");
  genUCD.buildNodeHierarchy();
  printf("****buildUCD\n");
  genUCD.buildUCD(nodes,cells);
  printf("nnodes = %u ncells = %u\n",nodes.getSize(),cells.getSize());
  puts("*** done ucd build");
}
