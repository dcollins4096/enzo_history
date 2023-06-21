/***********************************************************************
/
/  GENERATE GRID TO MPI TASK MAP
/
/  written by: Robert Harkness
/  date:       January, 2006
/
/  PURPOSE:
/
************************************************************************/

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <assert.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#define MAX_TASKS 2048
#define MAX_NODES 256
#define TASKS_PER_NODE 8

int AssignGridToTaskMap(Eint64 GridIndex[], Eint64 Memory[], int Ntask)
{

  int i, j, k, l, m, n;
  int err;
  int node;
  int gtemp;
  int mnode;
  double maxmem;

  double grid_size[MAX_TASKS];
  int grid[MAX_TASKS];
  int task_node[MAX_TASKS];
  int task[MAX_TASKS];

  double freemem[MAX_NODES];
  int freecpu[MAX_NODES];
  int actual_node[MAX_NODES];
  int node_input[MAX_NODES];
  int tasks_per_node[MAX_NODES];
  double nodemem[MAX_NODES];

  char node_name[80];


// Initialize MPI

  for ( i = 0; i < MAX_TASKS; i++ ) {
    grid_size[i] = 0.0;
    task_node[i] = -1;
  }

  for ( i = 0; i < MAX_NODES; i++ ) {
    actual_node[i] = 0;
    node_input[i] = 0;
    freemem[i] = 0.0;
    freecpu[i] = 0;
    nodemem[i] = 0.0;
  }

  MPI_Arg nt;
  MPI_Arg id;
  MPI_Arg nn;
  MPI_Arg lname;
  MPI_Arg node_number;
  MPI_Datatype DataTypeInt = (sizeof(Eint64) == 4) ? MPI_INT : MPI_LONG_LONG_INT;

  err = MPI_Comm_size(MPI_COMM_WORLD, &nt);
    assert( err == 0 );
  err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    assert( err == 0 );
  err = MPI_Get_processor_name(node_name, &lname);
    assert( err == 0 );

  sscanf(node_name, "ds%3d", &node_number);

  // TG  sscanf(node_name, "tg-c%3d", &node_number);

  fprintf(stderr, "Proc %d of %d is on node %s [%d]\n", id, nt, node_name, node_number);

  MPI_Barrier(MPI_COMM_WORLD);

//  Get a list of the actual node names to define properties
//  TASKS_PER_NODE is the number running in this job and NOT
//  necessaarily the same as in the target job

  nn = (nt-1)/TASKS_PER_NODE+1;

  if ( id == 0 )
    fprintf(stderr, "Number of nodes %d\n", nn);

  if ( id % TASKS_PER_NODE == 0 )
    node_input[id/TASKS_PER_NODE] = node_number;

  MPI_Allreduce(node_input, actual_node, nn, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);

// Define node properties
// As long as the total number of freecpu[] is > nt this should work

  for ( i = 0; i < nn; i++) {
    if ( actual_node[i] < 300 ) {
// TG    if ( actual_node[i] < 130 ) {
      freemem[i] = 12.5;
      tasks_per_node[i] = 8;
      freecpu[i] = tasks_per_node[i];
    } else {
      freemem[i] = 27.0;
      tasks_per_node[i] = 8;
      freecpu[i] = tasks_per_node[i];
    }
  }

  if ( id == 0 ) 
    for ( i = 0; i < nn; i++) {
      fprintf(stderr, "Logical Node %3d  Actual Node %3"ISYM"  Memory %6.2f\n", i, actual_node[i], freemem[i]);
    }


// Generate mapping

  for ( i = 0; i < nt; i++ ) {
//    fscanf(fp, "Grid %4d  PN %4d  Memory %16d\n", &grid[i], &old_pn[i], &gtemp);
    grid[i] = GridIndex[i];
    gtemp = Memory[i];
    grid_size[i] = ((double) gtemp) / 1073741824.0;  // in GBytes
  }

  if ( id == 0 ) {
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "%"ISYM"  %"ISYM"  %6.2f\n", i, grid[i], grid_size[i]);
    }
  }

// Sort grid[] and grid_size[] to monotonically decreasing grid_size

  int last1, save;
  double dsave;

  for ( m = 1; m < nt; m++ ) {
    last1 = nt-m+1;
    for ( l = 1; l < last1; l++) {
      if( grid_size[l] > grid_size[l-1] ) {
        dsave = grid_size[l-1];
        grid_size[l-1] = grid_size[l];
        grid_size[l] = dsave;
        save = grid[l-1];
        grid[l-1] = grid[l];
        grid[l] = save;
      }
    }
  }

  if ( id == 0 ) {
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "%"ISYM"  %"ISYM"  %6.2f\n", i, grid[i], grid_size[i]);
    }
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

/* Here we have on every task
   grid[i]               the grid number
   grid_size[i]          the size of the grid in memory
                         strictly monotonic decreasing

   nt              number of MPI tasks
   nn              number of physical nodes
   freemem[n]      free memory on node n
   freecpu[n]      free cpus on node n
*/

// assign grid[i] to task[j] on node[m]

  for ( i = 0; i < nt; i++ ) {

  maxmem = 0.0;
  mnode = -1;
  for ( n = 0; n < nn; n++ ) {
    if ( freecpu[n] > 0 ) {
      if ( freemem[n] > maxmem ) {
        maxmem = freemem[n];
        mnode = n;
      }
    }
  }

  freecpu[mnode] = freecpu[mnode] - 1;
  freemem[mnode] = freemem[mnode] - grid_size[i];
  task_node[i] = mnode;
  if ( freemem[mnode] < 0.0 )
    fprintf(stderr, "memory < 0 task %"ISYM" mnode %"ISYM" freemem %6.2f\n", i, mnode, freemem[mnode]);

  }

  MPI_Barrier(MPI_COMM_WORLD);

  for ( i = 0; i < nt; i++ ) {
     task[i] = i;
    if ( id == 0 ) fprintf(stderr, "Task %"ISYM"  Grid %"ISYM"  Node %"ISYM"  Size %6.2f\n", i, grid[i], task_node[i], grid_size[i]);
  }
  if ( id == 0 ) fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");

  MPI_Barrier(MPI_COMM_WORLD);

  for ( m = 1; m < nt; m++ ) {
    last1 = nt-m+1;
    for ( l = 1; l < last1; l++) {
      if( grid[l] < grid[l-1] ) {
        save = grid[l-1];
        grid[l-1] = grid[l];
        grid[l] = save;
        save = task_node[l-1];
        task_node[l-1] = task_node[l];
        task_node[l] = save;
        save = task[l-1];
        task[l-1] = task[l];
        task[l] = save;
        dsave = grid_size[l-1];
        grid_size[l-1] = grid_size[l];
        grid_size[l] = dsave;
      }
    }
  }

  for ( i = 0; i < nt; i++ ) {
    TaskMap[i] = task[i];
    TaskMemory[i] = 0;
  }

  if ( id == 0 ) {
    for ( i = 0; i < nt; i++ ) {
      fprintf(stderr, "Grid %"ISYM"  Task %"ISYM"  Node %"ISYM"  Size %6.2f\n", grid[i], task[i], task_node[i], grid_size[i]);
    }
    fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for ( i = 0; i < nt; i++ ) {
    nodemem[task_node[i]] = nodemem[task_node[i]] + grid_size[i];
  }

  if ( id == 0 )
  for ( n = 0; n < nn; n++) {
    fprintf(stderr, "Node %"ISYM"  Memory %8.2f GBytes\n", n, nodemem[n]);
  }
  if ( id == 0 ) fprintf(stderr, "+++++++++++++++++++++++++++++++++++\n");

  MPI_Barrier(MPI_COMM_WORLD);

  return SUCCESS;
}

