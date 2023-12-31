/***********************************************************************
/
/  UPDATE HYDRODYNAMICS EQUATIONS USING GPU
/
/  written by: Peng Wang, KIPAC/Stanford
/  date:       January, 2009
/  modified1: Tom Abel July 2009
/
/  description: 
/      This is the heart of the Hydro solver where fluxes at cell
/    interfaces are calculated using primitive variables at cell centers. 
/    The reconstruction is done using piecewise linear method (PLM) and 
/    the Riemann solver is Local Lax-Friedrichs (LLF) solver.
/      The GPU implementation uses one thread per cell. The key goal
/    is to allow arbitrary grid dimension in the parallization model, 
/    which is crucial for adaptive mesh calculations where arbitrary grid 
/    size can emerge. This is achieved by the use of method of line (MOL) spatial
/    discretization, in which flux calculation is made purely local (only need 4 
/    continuous cells for a flux is needed). This makes the parallelization
/    model easily extend to 3D.
/
************************************************************************/

#include "../macros_and_parameters.h"
#include "../typedefs.h"
#include "../global_data.h"
//#include <stdio.h>
#include <cutil.h>

// hack for making things compile
#define CUDA_BLOCK_SIZE 64
#define CUDA_GRID_SIZE 640
#define Gamma 1.001
#define Theta_Limiter 1.5
#define EOSType 3
#define EOSSoundSpeed 1.
#define NEQ_HYDRO 5
#define iD 0
#define iS1 1
#define iS2 2
#define iS3 3
#define iEtot 4


// forward declarations
__global__ void ComputeInternalEnergy_kernel(float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, int size );
__global__ void HydroSweepX_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
				         float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
					 int size);
__global__ void HydroSweepY_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
			   	         float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				         float *FluxTau, int size, int dim0, int dim1, int dim2);
__global__ void HydroSweepZ_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
					 float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				         float *FluxTau, int size, int dim0, int dim1, int dim2);
__global__ void HydroComputedUx_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau,
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size);
__global__ void HydroComputedUy_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size, int dim0, int dim1, int dim2);
__global__ void HydroComputedUz_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size, int dim0, int dim1, int dim2);
__global__ void HydroUpdatePrim_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau,
					     float dt, int size);
__device__ void LLF_PLM_CUDA3(float *Prim, float *Flux, const int &tx);
__device__ void plm_point(const float &vm1, const float &v, const float &vp1, float &vl_plm);
__device__ void EOS(float &p, float &rho, float &e, float &cs, const int &eostype, const int &mode);
__device__ float minmod(const float &a, const float &b, const float &c);
__device__ float Max(const float &a, const float &b, const float &c);
__device__ float Min(const float &a, const float &b, const float &c);


/******************************************************************************
/
/   Function HydroTimeUpdate_CUDA: update primitive variables by one time step
/     using method of lines by calling CUDA kernels
/
/   Input : 
/     Prim[NEQ_HYDRO][GridDimension^3] : primitive variables at cell center
/     GridDimension[3] : grid dimension
/     GridStartIndex[3] : the starting index of active data
/     GridRank : dimension of the problem
/     dtdx : dt/dx
/   Output: 
/     Prim[NEQ_HYDRO][GridDimension^3] get update by one time step
/
*******************************************************************************/

int HydroTimeUpdate_CUDA(float **Prim, int GridDimension[], 
			 int GridStartIndex[], int GridEndIndex[], int GridRank,
		          float dtdx, float dt)
{

  // compute grid size
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);

  // allocate device memory for primitive variables, with an additional field for internal energy
  float *PrimDevice = 0;
  int totalsize = sizeof(float)*(NEQ_HYDRO+1)*size;
  if (cudaMalloc((void**)&PrimDevice, totalsize) != cudaSuccess) {
    printf("cudaMalloc for PrimDevice with size %d failed.\n", totalsize);
    return FAIL;
  }

  // copy primitives to device memory
  cudaMemcpy(PrimDevice,        Prim[0], sizeof(float)*size, cudaMemcpyHostToDevice); // density
  cudaMemcpy(PrimDevice+size,   Prim[1], sizeof(float)*size, cudaMemcpyHostToDevice); // vx
  cudaMemcpy(PrimDevice+2*size, Prim[2], sizeof(float)*size, cudaMemcpyHostToDevice); // vy
  cudaMemcpy(PrimDevice+3*size, Prim[3], sizeof(float)*size, cudaMemcpyHostToDevice); // vz
  cudaMemcpy(PrimDevice+4*size, Prim[4], sizeof(float)*size, cudaMemcpyHostToDevice); // energy

  // set pointers of primitives on device
  float *Rho_Device  = PrimDevice, 
        *Vx_Device   = PrimDevice + size, 
	*Vy_Device   = PrimDevice + 2*size, 
	*Vz_Device   = PrimDevice + 3*size, 
	*Etot_Device = PrimDevice + 4*size, 
	*Eint_Device = PrimDevice + 5*size;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  cudaEventRecord(start, 0);

  // allocate device memory for flux
  float *FluxDevice = 0;
  if (cudaMalloc((void**)&FluxDevice, sizeof(float)*NEQ_HYDRO*size) != cudaSuccess) {
    cudaFree(PrimDevice);
    printf("cuda Malloc failed for FluxDevice with size %d\n", sizeof(float)*NEQ_HYDRO*size);
    return FAIL;
  }

  float *FluxD_Device   = FluxDevice, 
        *FluxS1_Device  = FluxDevice + size, 
	*FluxS2_Device  = FluxDevice + 2*size, 
	*FluxS3_Device  = FluxDevice + 3*size, 
	*FluxTau_Device = FluxDevice + 4*size;

  // allocate device memory for dU
  float *dUDevice = 0;
  if (cudaMalloc((void**)&dUDevice, sizeof(float)*NEQ_HYDRO*size) != cudaSuccess) {
    cudaFree(PrimDevice);
    cudaFree(FluxDevice);
    printf("cuda Malloc failed for dUDevice with size %d\n", sizeof(float)*NEQ_HYDRO*size);
    return FAIL;
  }

  float *dUD_Device   = dUDevice, 
        *dUS1_Device  = dUDevice + size, 
	*dUS2_Device  = dUDevice + 2*size, 
	*dUS3_Device  = dUDevice + 3*size, 
	*dUTau_Device = dUDevice + 4*size;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&elapsedTime, start, stop);

  //  PerformanceTimers[35] += elapsedTime/1e3;

  // compute gridSize
  dim3 dimBlock(CUDA_BLOCK_SIZE);
  dim3 dimGrid;
  if (size <= CUDA_BLOCK_SIZE*CUDA_GRID_SIZE) {
    dimGrid.x = size/CUDA_BLOCK_SIZE+1;
    dimGrid.y = 1;
  } else {
    dimGrid.x = CUDA_GRID_SIZE;
    dimGrid.y = size/(CUDA_BLOCK_SIZE*CUDA_GRID_SIZE) + 1;
  }

  // call the kernel functions to do the computation
  cudaEventRecord(start, 0);

  // computer internel energy
  ComputeInternalEnergy_kernel<<<dimGrid, dimBlock>>>(Vx_Device, Vy_Device, Vz_Device, Etot_Device, Eint_Device, size);
  cudaThreadSynchronize();

  // compute flux in x direction
  HydroSweepX_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, Eint_Device, 
					         FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device, 
					         size);
  // compute dU for x direction
  HydroComputedUx_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
					             dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device, dUTau_Device, 
					             dtdx, size);

  if (GridRank > 1) {
    HydroSweepY_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, Eint_Device, 
		 			           FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device, 
					           size, GridDimension[0], GridDimension[1], GridDimension[2]);

    HydroComputedUy_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
  						       dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device, dUTau_Device,
						       dtdx, size, GridDimension[0], GridDimension[1], GridDimension[2]);
  }

  if (GridRank > 2) {
    HydroSweepZ_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, Eint_Device,
		 			           FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device, 
					           size, GridDimension[0], GridDimension[1], GridDimension[2]);

    HydroComputedUz_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
  						       dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device, dUTau_Device,
					               dtdx, size, GridDimension[0], GridDimension[1], GridDimension[2]);
  }

  // update prim
  HydroUpdatePrim_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device,
   					             dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device, dUTau_Device,
						     dt, size);

  cudaEventRecord(stop, 0);
  cudaThreadSynchronize();
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);


  // copy prim back to cpu
  cudaEventRecord(start, 0);
  
  CUDA_SAFE_CALL(cudaMemcpy(Prim[0], PrimDevice,        sizeof(float)*size, cudaMemcpyDeviceToHost)); // density
  CUDA_SAFE_CALL(cudaMemcpy(Prim[1], PrimDevice+size,   sizeof(float)*size, cudaMemcpyDeviceToHost)); // vx
  CUDA_SAFE_CALL(cudaMemcpy(Prim[2], PrimDevice+2*size, sizeof(float)*size, cudaMemcpyDeviceToHost)); // vy
  CUDA_SAFE_CALL(cudaMemcpy(Prim[3], PrimDevice+3*size, sizeof(float)*size, cudaMemcpyDeviceToHost)); // vz
  CUDA_SAFE_CALL(cudaMemcpy(Prim[4], PrimDevice+4*size, sizeof(float)*size, cudaMemcpyDeviceToHost)); // energy
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //  PerformanceTimers[37] += elapsedTime/1e3;    

  cudaEventRecord(start, 0);

  CUDA_SAFE_CALL(cudaFree(PrimDevice));
  CUDA_SAFE_CALL(cudaFree(FluxDevice));
  CUDA_SAFE_CALL(cudaFree(dUDevice));

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //  PerformanceTimers[38] += elapsedTime/1e3;

  cudaEventDestroy( start ); 
  cudaEventDestroy( stop ); 

  return SUCCESS;

}

__global__ void ComputeInternalEnergy_kernel(float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid;
  igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igrid >= size)
    return;

  // internal energy = total energy - kinetic energy
  Eneint[igrid] = Etot[igrid] - 0.5*(Vx[igrid]*Vx[igrid] + Vy[igrid]*Vy[igrid] + Vz[igrid]*Vz[igrid]);

}
  
// kernel: compute flux in the x direction
__global__ void HydroSweepX_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
				         float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				         float *FluxTau, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_HYDRO*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_HYDRO*CUDA_BLOCK_SIZE]; // output flux

  if (igrid >= 2 && igrid <= size - 2) { // only do flux computation for active cells

    // load data from device to shared.
    int idx_prim1 = (tx+2)*NEQ_HYDRO;
    PrimLine[idx_prim1++] = Rho [igrid];
    PrimLine[idx_prim1++] = Eneint[igrid];
    PrimLine[idx_prim1++] = Vx  [igrid];
    PrimLine[idx_prim1++] = Vy  [igrid];
    PrimLine[idx_prim1  ] = Vz  [igrid];

    // if the first, load in two more cells for boundary condition
    if (tx == 0 || igrid == 2) {
      for (int i = -2; i <=-1; i++) {
        const int idx_prim  = igrid + i;
        int idx_prim1 = (i+tx+2)*NEQ_HYDRO;
        PrimLine[idx_prim1++] = Rho [idx_prim];
        PrimLine[idx_prim1++] = Eneint[idx_prim];
        PrimLine[idx_prim1++] = Vx  [idx_prim];
        PrimLine[idx_prim1++] = Vy  [idx_prim];
        PrimLine[idx_prim1  ] = Vz  [idx_prim];
      }
    }

    // if the last, load in one more cell for boundary condition
    if (tx == CUDA_BLOCK_SIZE - 1 || igrid == size - 2) {
      const int idx_prim  = igrid + 1;
      int idx_prim1 = (tx+3)*NEQ_HYDRO;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vx  [idx_prim];
      PrimLine[idx_prim1++] = Vy  [idx_prim];
      PrimLine[idx_prim1  ] = Vz  [idx_prim];
    }

  }

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igrid >= 2 && igrid <= size - 2) {
    // the main computation: calculating the flux at tx
    LLF_PLM_CUDA3(PrimLine, FluxLine, tx);

    // copy 1D Flux back to Flux
    int idx_prim1 = tx*NEQ_HYDRO;
    FluxD  [igrid] = FluxLine[idx_prim1++];
    FluxS1 [igrid] = FluxLine[idx_prim1++];
    FluxS2 [igrid] = FluxLine[idx_prim1++];
    FluxS3 [igrid] = FluxLine[idx_prim1++];
    FluxTau[igrid] = FluxLine[idx_prim1  ];
  }
}

// kernel: compute flux in the y direction
__global__ void HydroSweepY_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
			   	         float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				         float *FluxTau, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridy = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  // convert to array index
  int k = igridy/(dim0*dim1);
  int i = (igridy - k*dim0*dim1)/dim1;
  int j = igridy - k*dim0*dim1 - i*dim1;
  int igrid = i + (j + k*dim1) * dim0;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_HYDRO*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_HYDRO*CUDA_BLOCK_SIZE]; // output flux

  if (igridy >= 2 && igridy <= size - 2) {

    // load data from device to shared.
    int idx_prim1 = (tx+2)*NEQ_HYDRO;
    PrimLine[idx_prim1++] = Rho [igrid];
    PrimLine[idx_prim1++] = Eneint[igrid];
    PrimLine[idx_prim1++] = Vy  [igrid]; // vx = vy
    PrimLine[idx_prim1++] = Vz  [igrid]; // vy = vz
    PrimLine[idx_prim1  ] = Vx  [igrid]; // vz = vx

    // if the first, load in two more cells for boundary condition
    if (tx == 0 || igridy == 2) {
      for (int di = -2; di <=-1; di++) {
        int igridypi = igridy + di;
        int k = igridypi/(dim0*dim1);
        int i = (igridypi - k*dim0*dim1)/dim1;
        int j = igridypi - k*dim0*dim1 - i*dim1;
        int idx_prim = i + (j + k*dim1) * dim0;

        int idx_prim1 = (di+tx+2)*NEQ_HYDRO;
        PrimLine[idx_prim1++] = Rho [idx_prim];
        PrimLine[idx_prim1++] = Eneint[idx_prim];
        PrimLine[idx_prim1++] = Vy  [idx_prim];
        PrimLine[idx_prim1++] = Vz  [idx_prim];
        PrimLine[idx_prim1  ] = Vx  [idx_prim];
      }
    }

    // if the last, load in one more cell for boundary condition
    if (tx == CUDA_BLOCK_SIZE - 1 || igridy == size - 2) {
      int igridyp1 = igridy + 1;
      int k = igridyp1/(dim0*dim1);
      int i = (igridyp1 - k*dim0*dim1)/dim1;
      int j = igridyp1 - k*dim0*dim1 - i*dim1;
      int idx_prim = i + (j + k*dim1) * dim0;

      int idx_prim1 = (tx+3)*NEQ_HYDRO;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vy  [idx_prim];
      PrimLine[idx_prim1++] = Vz  [idx_prim];
      PrimLine[idx_prim1  ] = Vx  [idx_prim];
    }
  }

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igridy >= 2 && igridy <= size - 2) {
    // the main computation: calculating the flux at tx
    LLF_PLM_CUDA3(PrimLine, FluxLine, tx);

    // copy 1D Flux back to Flux
    int idx_prim1 = tx*NEQ_HYDRO;
    FluxD  [igrid] = FluxLine[idx_prim1++];
    FluxS2 [igrid] = FluxLine[idx_prim1++]; // vy = vx
    FluxS3 [igrid] = FluxLine[idx_prim1++]; // vz = vy
    FluxS1 [igrid] = FluxLine[idx_prim1++]; // vx = vz
    FluxTau[igrid] = FluxLine[idx_prim1  ];
  }

}

// kernel: compute flux in z direction
__global__ void HydroSweepZ_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot, float *Eneint, 
					 float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				         float *FluxTau, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridz = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  // convert to array index
  int j = igridz / (dim0*dim2);
  int i = (igridz - j*dim0*dim2) / dim2;
  int k = igridz - j*dim0*dim2 - i*dim2;
  int igrid = i + (j + k*dim1) * dim0;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_HYDRO*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_HYDRO*CUDA_BLOCK_SIZE]; // output flux

  if (igridz >= 2 && igridz <= size - 2) {

    // load data from device to shared.
    int idx_prim1 = (tx+2)*NEQ_HYDRO;
    PrimLine[idx_prim1++] = Rho [igrid];
    PrimLine[idx_prim1++] = Eneint[igrid];
    PrimLine[idx_prim1++] = Vz  [igrid]; // vx = vz
    PrimLine[idx_prim1++] = Vx  [igrid]; // vy = vx
    PrimLine[idx_prim1  ] = Vy  [igrid]; // vz = vy
  
    // if the first, load in two more cells for boundary condition
    if (tx == 0 || igridz == 2) {
      for (int di = -2; di <=-1; di++) {

        int igridzpi = igridz + di;
        int j = igridzpi / (dim0*dim2);
        int i = (igridzpi - j*dim0*dim2) / dim2;
        int k = igridzpi - j*dim0*dim2 - i*dim2;
        int idx_prim = i + (j + k*dim1) * dim0;

        int idx_prim1 = (di+tx+2)*NEQ_HYDRO;
        PrimLine[idx_prim1++] = Rho [idx_prim];
        PrimLine[idx_prim1++] = Eneint[idx_prim];
        PrimLine[idx_prim1++] = Vz  [idx_prim];
        PrimLine[idx_prim1++] = Vx  [idx_prim];
        PrimLine[idx_prim1  ] = Vy  [idx_prim];
      }
    }

    // if the last, load in one more cell for boundary condition
    if (tx == CUDA_BLOCK_SIZE - 1 || igridz == size - 2) {
      int igridzp1 = igridz + 1;
      int j = igridzp1 / (dim0*dim2);
      int i = (igridzp1 - j*dim0*dim2) / dim2;
      int k = igridzp1 - j*dim0*dim2 - i*dim2;
      int idx_prim = i + (j + k*dim1) * dim0;

      int idx_prim1 = (tx+3)*NEQ_HYDRO;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vz  [idx_prim];
      PrimLine[idx_prim1++] = Vx  [idx_prim];
      PrimLine[idx_prim1  ] = Vy  [idx_prim];
    }
  } 

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igridz >= 2 && igridz <= size - 2) {
    // the main computation: calculating the flux at tx
    LLF_PLM_CUDA3(PrimLine, FluxLine, tx);

    // copy 1D Flux back to Flux
    int idx_prim1 = tx*NEQ_HYDRO;
    FluxD  [igrid] = FluxLine[idx_prim1++];
    FluxS3 [igrid] = FluxLine[idx_prim1++]; // vz = vx
    FluxS1 [igrid] = FluxLine[idx_prim1++]; // vx = vy
    FluxS2 [igrid] = FluxLine[idx_prim1++]; // vy = vz
    FluxTau[igrid] = FluxLine[idx_prim1  ];
  }
}

// compute flux difference in x direction
__global__ void HydroComputedUx_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau,
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igrid < 2 || igrid > size - 3)
    return;

  int igridp1 = igrid + 1;
  dUD  [igrid] = (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] = (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] = (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] = (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] = (FluxTau[igrid] - FluxTau[igridp1])*dtdx;

}

// compute flux difference in y direction
__global__ void HydroComputedUy_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridy = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igridy < 2 || igridy > size - 3)
    return;

  int k = igridy/(dim0*dim1);
  int i = (igridy - k*dim0*dim1)/dim1;
  int j = igridy - k*dim0*dim1 - i*dim1;
  int igrid = i + (j + k*dim1) * dim0;

  int igridyp1 = igridy + 1;
  k = igridyp1/(dim0*dim1);
  i = (igridyp1 - k*dim0*dim1)/dim1;
  j = igridyp1 - k*dim0*dim1 - i*dim1;
  int igridp1 = i + (j + k*dim1) * dim0;


  dUD  [igrid] += (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] += (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] += (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] += (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] += (FluxTau[igrid] - FluxTau[igridp1])*dtdx;

}

// compute flux difference in z direction
__global__ void HydroComputedUz_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3, float *FluxTau, 
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau, 
					     float dtdx, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridz = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igridz < 2 || igridz > size - 3)
    return;

  int j = igridz / (dim0*dim2);
  int i = (igridz - j*dim0*dim2) / dim2;
  int k = igridz - j*dim0*dim2 - i*dim2;
  int igrid = i + (j + k*dim1) * dim0;

  int igridzp1 = igridz + 1;
  j = igridzp1 / (dim0*dim2);
  i = (igridzp1 - j*dim0*dim2) / dim2;
  k = igridzp1 - j*dim0*dim2 - i*dim2;
  int igridp1 = i + (j + k*dim1) * dim0;

  dUD  [igrid] += (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] += (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] += (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] += (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] += (FluxTau[igrid] - FluxTau[igridp1])*dtdx;

}

__global__ void HydroUpdatePrim_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
					     float *dUD, float *dUS1, float *dUS2, float *dUS3, float *dUTau,
					     float dt, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igrid < 2 || igrid > size - 3)
    return;

  float D, S1, S2, S3, Tau;
  D   = Rho[igrid];
  S1  = D*Vx[igrid];
  S2  = D*Vy[igrid];
  S3  = D*Vz[igrid];
  Tau = D*Etot[igrid];

  D   += dUD[igrid];
  S1  += dUS1[igrid];
  S2  += dUS2[igrid];
  S3  += dUS3[igrid];
  Tau += dUTau[igrid];

  Rho[igrid] = D;
  Vx[igrid] = S1/D;
  Vy[igrid] = S2/D;
  Vz[igrid] = S3/D;
  Etot[igrid] = Tau/D;
  
}

// the main computation routine: compute Flux at the left cell interface given Prim at the center
__device__ void LLF_PLM_CUDA3(float *Prim, float *Flux, const int &tx)
{
  // those crazily many fields are defined in stead of using arrays
  // to avoid letting the compiler putting arrays in local memory
  float Priml_rho, Priml_eint, Priml_vx, Priml_vy, Priml_vz;
  float Primr_rho, Primr_eint, Primr_vx, Primr_vy, Primr_vz;
  float Ul_D, Ul_S1, Ul_S2, Ul_S3, Ul_Etot;
  float Ur_D, Ur_S1, Ur_S2, Ur_S3, Ur_Etot;
  float Fl_D, Fl_S1, Fl_S2, Fl_S3, Fl_Etot;
  float Fr_D, Fr_S1, Fr_S2, Fr_S3, Fr_Etot;
  float cs, v2, p, lm_l, lp_l, lm_r, lp_r;

  // 1. do PLM reconstruction for all primitive fields

  plm_point(Prim[tx*NEQ_HYDRO], Prim[(tx+1)*NEQ_HYDRO], Prim[(tx+2)*NEQ_HYDRO], Priml_rho);
  plm_point(Prim[(tx+3)*NEQ_HYDRO], Prim[(tx+2)*NEQ_HYDRO], Prim[(tx+1)*NEQ_HYDRO], Primr_rho);

  plm_point(Prim[tx*NEQ_HYDRO+1], Prim[(tx+1)*NEQ_HYDRO+1], Prim[(tx+2)*NEQ_HYDRO+1], Priml_eint);
  plm_point(Prim[(tx+3)*NEQ_HYDRO+1], Prim[(tx+2)*NEQ_HYDRO+1], Prim[(tx+1)*NEQ_HYDRO+1], Primr_eint);

  plm_point(Prim[tx*NEQ_HYDRO+2], Prim[(tx+1)*NEQ_HYDRO+2], Prim[(tx+2)*NEQ_HYDRO+2], Priml_vx);
  plm_point(Prim[(tx+3)*NEQ_HYDRO+2], Prim[(tx+2)*NEQ_HYDRO+2], Prim[(tx+1)*NEQ_HYDRO+2], Primr_vx);

  plm_point(Prim[tx*NEQ_HYDRO+3], Prim[(tx+1)*NEQ_HYDRO+3], Prim[(tx+2)*NEQ_HYDRO+3], Priml_vy);
  plm_point(Prim[(tx+3)*NEQ_HYDRO+3], Prim[(tx+2)*NEQ_HYDRO+3], Prim[(tx+1)*NEQ_HYDRO+3], Primr_vy);

  plm_point(Prim[tx*NEQ_HYDRO+4], Prim[(tx+1)*NEQ_HYDRO+4], Prim[(tx+2)*NEQ_HYDRO+4], Priml_vz);
  plm_point(Prim[(tx+3)*NEQ_HYDRO+4], Prim[(tx+2)*NEQ_HYDRO+4], Prim[(tx+1)*NEQ_HYDRO+4], Primr_vz);

  // 2. use LLF Riemann solver to compute flux  

  // 2.1, compute Fl and Ul
  v2 = pow(Priml_vx,2) + pow(Priml_vy,2) + pow(Priml_vz,2);

  EOS(p, Priml_rho, Priml_eint, cs, EOSType, 2);    

  Ul_D    = Priml_rho;
  Ul_S1   = Priml_rho * Priml_vx;
  Ul_S2   = Priml_rho * Priml_vy;
  Ul_S3   = Priml_rho * Priml_vz;
  Ul_Etot = Priml_rho * (Priml_eint + 0.5*v2);

  Fl_D    = Priml_rho * Priml_vx;
  Fl_S1   = Ul_S1 * Priml_vx + p;
  Fl_S2   = Ul_S2 * Priml_vx;
  Fl_S3   = Ul_S3 * Priml_vx;
  Fl_Etot = Priml_rho*(0.5*v2 + Priml_eint + p/Priml_rho)*Priml_vx;

  // largest and smallest eigenvectors, reuse the space of cs

  lp_l = Priml_vx + cs;
  lm_l = Priml_vx - cs;

  // 2.2 Fr and Ur
  v2 = pow(Primr_vx,2) + pow(Primr_vy,2) + pow(Primr_vz,2);

  EOS(p, Primr_rho, Primr_eint, cs, EOSType, 2);

  Ur_D    = Primr_rho;
  Ur_S1   = Primr_rho * Primr_vx;
  Ur_S2   = Primr_rho * Primr_vy;
  Ur_S3   = Primr_rho * Primr_vz;
  Ur_Etot = Primr_rho * (Primr_eint + 0.5*v2);

  Fr_D    = Primr_rho * Primr_vx;
  Fr_S1   = Ur_S1 * Primr_vx + p;
  Fr_S2   = Ur_S2 * Primr_vx;
  Fr_S3   = Ur_S3 * Primr_vx;
  Fr_Etot = Primr_rho*(0.5*v2 + Primr_eint + p/Primr_rho)*Primr_vx;

  // largest and smallest eigenvectors, reuse the space of cs
  lp_r = Primr_vx + cs;
  lm_r = Primr_vx - cs;

  // 2.3. compute the maximum local wave speed
  lp_l = Max(0, lp_l, lp_r);
  lm_l = Max(0, -lm_l, -lm_r);
  lp_l = (lp_l > lm_l) ? lp_l : lm_l;

  // 2.4. compute the flux
  Flux[tx*NEQ_HYDRO+iD   ] = 0.5*(Fl_D + Fr_D - lp_l*(Ur_D - Ul_D));
  Flux[tx*NEQ_HYDRO+iS1  ] = 0.5*(Fl_S1 + Fr_S1 - lp_l*(Ur_S1 - Ul_S1));
  Flux[tx*NEQ_HYDRO+iS2  ] = 0.5*(Fl_S2 + Fr_S2 - lp_l*(Ur_S2 - Ul_S2));
  Flux[tx*NEQ_HYDRO+iS3  ] = 0.5*(Fl_S3 + Fr_S3 - lp_l*(Ur_S3 - Ul_S3));
  Flux[tx*NEQ_HYDRO+iEtot] = 0.5*(Fl_Etot + Fr_Etot - lp_l*(Ur_Etot - Ul_Etot));

}

// minmod limiter function
__device__ float minmod(const float &a, const float &b, const float &c)
{
  return 0.25*(sign(a)+sign(b))*ABS((sign(a)+sign(c)))*Min(ABS(a), ABS(b), ABS(c));
}

// do PLM reconstruction for a point
__device__ void plm_point(const float &vm1, const float &v, const float &vp1, float &vl_plm)
{

  float dv_l = (v-vm1) * Theta_Limiter;
  float dv_r = (vp1-v) * Theta_Limiter;
  float dv_m = 0.5*(vp1-vm1);
  
  float dv = minmod(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

// return the maximum of a, b, c
__device__ float Max(const float &a, const float &b, const float &c)  
{
  if (a > b) {
    if (a > c)
      return a;
    else 
      return c;
  } else {
    if (b > c)
      return b;
    else
      return c;
  }
}

// return the minimum of a, b, c
__device__ float Min(const float &a, const float &b, const float &c)
{
  if (a<b) {
    if (c<a)
      return c;
    else 
      return a;
  } else {
    if (c<b)
      return c;
    else 
      return b;
  }
}

/***********************************************************
/
/   Function EOS: compute equation of state
/
/   Input: 
/     eostype: 
/       0: ideal gas
/     mode:  
/       1: given p and rho, calculate others.
/       2: given rho and e, calculate others.
/
************************************************************/

__device__ void EOS(float &p, float &rho, float &e, float &cs, const int &eostype, const int &mode)
{

  if (eostype == 0) {
    
    if (mode == 1) {
      e = p / rho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e;
    }

    cs = sqrt(Gamma*p/rho);

  }

  if (eostype == 3) { // straight isothermal
    cs = EOSSoundSpeed;
    p = rho*cs*cs;
    e = p / ((Gamma-1.0)*rho);
  }


}
