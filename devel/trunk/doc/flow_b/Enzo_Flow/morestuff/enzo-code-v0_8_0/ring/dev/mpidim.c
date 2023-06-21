#include <stdio.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
int me, us;
int rank;
int mpi_layout[3] = {0,0,0};

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &me);
MPI_Comm_size(MPI_COMM_WORLD, &us);

printf("Processor %d of %d\n",me,us);
rank = 3;
scanf("%d",&us);
printf("PEs %d\n",us);
MPI_Dims_create(us, rank, mpi_layout);
printf("MPI_layout %d %d %d\n",mpi_layout[0],mpi_layout[1],mpi_layout[2]);

MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
}
