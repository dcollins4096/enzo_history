/*****************************************************************************
 *                                                                           *
 * Copyright 2005 Daniel R. Reynolds
 * Copyright 2005 Laboratory for Computational Astrophysics                  *
 * Copyright 2005 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  This routine fills in ghost cells for EnzoVectors grids based on 
/  neighbor information.
/
/  Note: we perform all x0 communication first, then all x1 
/  communication, then all x2 communication.  This ensures that values 
/  at edges and corners of the box are communicated appropriately.
/
/  written by: Daniel R. Reynolds
/  date:       May, 2006
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "EnzoVector.h"

//  Vector Boundary Communication Routine
//  (may be used for parallelism, or even for single-proc. periodic BCs)
int EnzoVector::exchange()
{
#ifdef USE_MPI

  // some local variables
  int i, j, k, idx;
  MPI_Arg myrank;
  MPI_Arg one=1;
  MPI_Datatype FDataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Datatype IDataType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
  int x0len = Nx0 + Ng0l + Ng0r;
  int x1len = Nx1 + Ng1l + Ng1r;
  int x2len = Nx2 + Ng2l + Ng2r;

  // Get MPI processor rank
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

//   fprintf(stdout,"p%i: entering EnzoVector::exchange routine\n",myrank);
//   fprintf(stdout,"p%i:   lengths = (%"ISYM",%"ISYM",%"ISYM")\n",
// 	  myrank,x0len,x1len,x2len);

  // set exchange tags
  int msg_xch_x0l = 10000;
  int msg_xch_x0r = 10001;
  int msg_xch_x1l = 10002;
  int msg_xch_x1r = 10003;
  int msg_xch_x2l = 10004;
  int msg_xch_x2r = 10005;

  // allocate request IDs
  MPI_Request id_recv_x0l, id_send_x0l;
  MPI_Request id_recv_x0r, id_send_x0r;
  MPI_Request id_recv_x1l, id_send_x1l;
  MPI_Request id_recv_x1r, id_send_x1r;
  MPI_Request id_recv_x2l, id_send_x2l;
  MPI_Request id_recv_x2r, id_send_x2r;
  
  // allocate MPI status object
  MPI_Status status;

//   fprintf(stdout,"  p%i: determining send buffers, Nbors = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM"), MPI_PROC_NULL = %"ISYM"\n",myrank,Nbors[0][0],Nbors[0][1],Nbors[1][0],Nbors[1][1],Nbors[2][0],Nbors[2][1],int(MPI_PROC_NULL));

  /////////////////////////////////////////////////////////
  // Determine send buffer sizes (initialize to Dirichlet value)
  int NborGhosts[3][2];
  NborGhosts[0][0] = 1;  NborGhosts[0][1] = 1;
  NborGhosts[1][0] = 1;  NborGhosts[1][1] = 1;
  NborGhosts[2][0] = 1;  NborGhosts[2][1] = 1;


  //   open receive value for x0L proc ghosts
  if (Nbors[0][0] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[0][0]), one, IDataType, MPI_Arg(Nbors[0][0]), 
		  msg_xch_x0l, MPI_COMM_WORLD, &id_recv_x0l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0L receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x0R proc ghosts
  if (Nbors[0][1] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng0r), one, IDataType, MPI_Arg(Nbors[0][1]), 
		  msg_xch_x0l, MPI_COMM_WORLD, &id_send_x0l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0R send error\n",myrank);
      return FAIL;
    }
  //   wait for x0L proc ghosts
  if (Nbors[0][0] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x0l, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0L wait error\n",myrank);
      return FAIL;
    }

  //   open receive value for x0R proc ghosts
  if (Nbors[0][1] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[0][1]), one, IDataType, MPI_Arg(Nbors[0][1]), 
		  msg_xch_x0r, MPI_COMM_WORLD, &id_recv_x0r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0R receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x0L proc ghosts
  if (Nbors[0][0] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng0l), one, IDataType, MPI_Arg(Nbors[0][0]), 
		  msg_xch_x0r, MPI_COMM_WORLD, &id_send_x0r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0L send error\n",myrank);
      return FAIL;
    }
  //   wait for x0R proc ghosts
  if (Nbors[0][1] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x0r, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x0R wait error\n",myrank);
      return FAIL;
    }

  //   open receive value for x1L proc ghosts
  if (Nbors[1][0] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[1][0]), one, IDataType, MPI_Arg(Nbors[1][0]), 
		  msg_xch_x1l, MPI_COMM_WORLD, &id_recv_x1l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1L receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x1R proc ghosts
  if (Nbors[1][1] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng1r), one, IDataType, MPI_Arg(Nbors[1][1]), 
		  msg_xch_x1l, MPI_COMM_WORLD, &id_send_x1l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1R send error\n",myrank);
      return FAIL;
    }
  //   wait for x1L proc ghosts
  if (Nbors[1][0] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x1l, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1L wait error\n",myrank);
      return FAIL;
    }

  //   open receive value for x1R proc ghosts
  if (Nbors[1][1] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[1][1]), one, IDataType, MPI_Arg(Nbors[1][1]), 
		  msg_xch_x1r, MPI_COMM_WORLD, &id_recv_x1r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1R receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x1L proc ghosts
  if (Nbors[1][0] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng1l), one, IDataType, MPI_Arg(Nbors[1][0]), 
		  msg_xch_x1r, MPI_COMM_WORLD, &id_send_x1r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1L send error\n",myrank);
      return FAIL;
    }
  //   wait for x1R proc ghosts
  if (Nbors[1][1] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x1r, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x1R wait error\n",myrank);
      return FAIL;
    }

  //   open receive value for x2L proc ghosts
  if (Nbors[2][0] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[2][0]), one, IDataType, MPI_Arg(Nbors[2][0]), 
		  msg_xch_x2l, MPI_COMM_WORLD, &id_recv_x2l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2L receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x2R proc ghosts
  if (Nbors[2][1] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng2r), one, IDataType, MPI_Arg(Nbors[2][1]), 
		  msg_xch_x2l, MPI_COMM_WORLD, &id_send_x2l) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2R send error\n",myrank);
      return FAIL;
    }
  //   wait for x2L proc ghosts
  if (Nbors[2][0] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x2l, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2L wait error\n",myrank);
      return FAIL;
    }

  //   open receive value for x2R proc ghosts
  if (Nbors[2][1] != MPI_PROC_NULL)
    if (MPI_Irecv(&(NborGhosts[2][1]), one, IDataType, MPI_Arg(Nbors[2][1]), 
		  msg_xch_x2r, MPI_COMM_WORLD, &id_recv_x2r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2R receive error\n",myrank);
      return FAIL;
    }
  //   fill and send value for x2L proc ghosts
  if (Nbors[2][0] != MPI_PROC_NULL)
    if (MPI_Isend(&(Ng2l), one, IDataType, MPI_Arg(Nbors[2][0]), 
		  msg_xch_x2r, MPI_COMM_WORLD, &id_send_x2r) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2L send error\n",myrank);
      return FAIL;
    }
  //   wait for x2R proc ghosts
  if (Nbors[2][1] != MPI_PROC_NULL)
    if (MPI_Wait(&id_recv_x2r, &status) != 0) {
      fprintf(stderr,"GravityBdryExchange p%i: x2R wait error\n",myrank);
      return FAIL;
    }

//   fprintf(stdout,"  p%"ISYM": my ghosts = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,Ng0l,Ng0r,Ng1l,Ng1r,Ng2l,Ng2r);
//   fprintf(stdout,"  p%"ISYM": neighbor ghosts = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,NborGhosts[0][0],NborGhosts[0][1],NborGhosts[1][0],NborGhosts[1][1],NborGhosts[2][0],NborGhosts[2][1]);

  /////////////////////////////////////////////////////////
  // set buffer sizes, allocate send/receive buffers
  MPI_Arg x0Lsendsize = x1len*x2len*NborGhosts[0][0];
  MPI_Arg x0Rsendsize = x1len*x2len*NborGhosts[0][1];
  MPI_Arg x0Lrecvsize = x1len*x2len*Ng0l;
  MPI_Arg x0Rrecvsize = x1len*x2len*Ng0r;
  
  MPI_Arg x1Lsendsize = x0len*x2len*NborGhosts[1][0];
  MPI_Arg x1Rsendsize = x0len*x2len*NborGhosts[1][1];
  MPI_Arg x1Lrecvsize = x0len*x2len*Ng1l;
  MPI_Arg x1Rrecvsize = x0len*x2len*Ng1r;
  
  MPI_Arg x2Lsendsize = x0len*x1len*NborGhosts[2][0];
  MPI_Arg x2Rsendsize = x0len*x1len*NborGhosts[2][1];
  MPI_Arg x2Lrecvsize = x0len*x1len*Ng2l;
  MPI_Arg x2Rrecvsize = x0len*x1len*Ng2r;

  int SBuffSize = (x0Lsendsize > x0Rsendsize) ? x0Lsendsize : x0Rsendsize;
  SBuffSize = (x1Lsendsize > SBuffSize) ? x1Lsendsize : SBuffSize;
  SBuffSize = (x1Rsendsize > SBuffSize) ? x1Rsendsize : SBuffSize;
  SBuffSize = (x2Lsendsize > SBuffSize) ? x2Lsendsize : SBuffSize;
  SBuffSize = (x2Rsendsize > SBuffSize) ? x2Rsendsize : SBuffSize;

  int RBuffSize = (x0Lrecvsize > x0Rrecvsize) ? x0Lrecvsize : x0Rrecvsize;
  RBuffSize = (x1Lrecvsize > RBuffSize) ? x1Lrecvsize : RBuffSize;
  RBuffSize = (x1Rrecvsize > RBuffSize) ? x1Rrecvsize : RBuffSize;
  RBuffSize = (x2Lrecvsize > RBuffSize) ? x2Lrecvsize : RBuffSize;
  RBuffSize = (x2Rrecvsize > RBuffSize) ? x2Rrecvsize : RBuffSize;
  float *SendBuf = new float[SBuffSize];
  for (int i=0; i<SBuffSize; i++)  SendBuf[i] = 0.0;
  float *RecvBuf = new float[RBuffSize];
  for (int i=0; i<RBuffSize; i++)  RecvBuf[i] = -1.0;

//   fprintf(stdout,"  p%"ISYM": send buffer sizes = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,x0Lsendsize,x0Rsendsize,x1Lsendsize,x1Rsendsize,x2Lsendsize,x2Rsendsize);
//   fprintf(stdout,"  p%"ISYM": recv buffer sizes = (%"ISYM":%"ISYM",%"ISYM":%"ISYM",%"ISYM":%"ISYM")\n",MyProcessorNumber,x0Lrecvsize,x0Rrecvsize,x1Lrecvsize,x1Rrecvsize,x2Lrecvsize,x2Rrecvsize);


  /////////////////////////////////////////////////////////
  // iterate over all of the species, communicating data to neighbors
  // (do one species at a time)
  float *mydata;
  for (int idat=0; idat<Nspecies; idat++) {

//     fprintf(stdout,"  p%i: updating species %"ISYM"\n",
// 	    myrank,idat);
  
    // extract this variable from the vector array
    mydata = data[idat];
    
//     fprintf(stdout,"  p%i: updating x0L bdry\n",myrank);
  
    /////////////////////////
    // Update x0L boundaries
    {
      // open receive buffer for x0L boundary
      if (Nbors[0][0] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x0Lrecvsize, FDataType, MPI_Arg(Nbors[0][0]), 
		      msg_xch_x0l, MPI_COMM_WORLD, &id_recv_x0l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0L receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x0R boundary
      if (Nbors[0][1] != MPI_PROC_NULL) {
	idx = 0;
	for (k=0; k<x2len; k++)
	  for (j=0; j<x1len; j++)
	    for (i=x0len-Ng0r-NborGhosts[0][1]; i<x0len-Ng0r; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];
      
	if (MPI_Isend(SendBuf, x0Rsendsize, FDataType, MPI_Arg(Nbors[0][1]), 
		      msg_xch_x0l, MPI_COMM_WORLD, &id_send_x0l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0R send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x0L data, and update ghost cells
      if (Nbors[0][0] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x0l, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0L wait error\n",myrank);
	  return FAIL;
	}
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=0; j<x1len; j++)
	    for (i=0; i<Ng0l; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    

//     fprintf(stdout,"  p%i: updating x0R bdry\n",myrank);  

    /////////////////////////
    // Update x0R boundaries
    {
      // open receive buffer for x0R boundary
      if (Nbors[0][1] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x0Rrecvsize, FDataType, MPI_Arg(Nbors[0][1]), 
		      msg_xch_x0r, MPI_COMM_WORLD, &id_recv_x0r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0R receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x0L boundary
      if (Nbors[0][0] != MPI_PROC_NULL) {
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=0; j<x1len; j++)
	    for (i=Ng0l; i<Ng0l+NborGhosts[0][0]; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];
      
	if (MPI_Isend(SendBuf, x0Lsendsize, FDataType, MPI_Arg(Nbors[0][0]), 
		      msg_xch_x0r, MPI_COMM_WORLD, &id_send_x0r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0L send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x0R data, and update ghost cells
      if (Nbors[0][1] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x0r, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x0R wait error\n",myrank);
	  return FAIL;
	}
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=0; j<x1len; j++)
	    for (i=x0len-Ng0r; i<x0len; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }


//     fprintf(stdout,"  p%i: updating x1L bdry\n",myrank);  

    /////////////////////////
    // Update x1L boundaries
    {
      // open receive buffer for x1L boundary

      if (Nbors[1][0] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x1Lrecvsize, FDataType, MPI_Arg(Nbors[1][0]), 
		      msg_xch_x1l, MPI_COMM_WORLD, &id_recv_x1l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1L receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x1R boundary
      if (Nbors[1][1] != MPI_PROC_NULL) {
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=x1len-Ng1r-NborGhosts[1][1]; j<x1len-Ng1r; j++)
	    for (i=0; i<x0len; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];
      
	if (MPI_Isend(SendBuf, x1Rsendsize, FDataType, MPI_Arg(Nbors[1][1]), 
		      msg_xch_x1l, MPI_COMM_WORLD, &id_send_x1l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1R send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x1L data, and update ghost cells
      if (Nbors[1][0] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x1l, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1L wait error\n",myrank);
	  return FAIL;
	}
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=0; j<Ng1l; j++)
	    for (i=0; i<x0len; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }


//     fprintf(stdout,"  p%i: updating x1R bdry\n",myrank);  

    /////////////////////////
    // Update x1R boundaries
    {
      // open receive buffer for x1R boundary
      if (Nbors[1][1] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x1Rrecvsize, FDataType, MPI_Arg(Nbors[1][1]), 
		      msg_xch_x1r, MPI_COMM_WORLD, &id_recv_x1r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1R receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x1L boundary
      if (Nbors[1][0] != MPI_PROC_NULL) {
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=Ng1l; j<Ng1l+NborGhosts[1][0]; j++)
	    for (i=0; i<x0len; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];
      
	if (MPI_Isend(SendBuf, x1Lsendsize, FDataType, MPI_Arg(Nbors[1][0]), 
		      msg_xch_x1r, MPI_COMM_WORLD, &id_send_x1r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1L send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x1R data, and update ghost cells
      if (Nbors[1][1] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x1r, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x1R wait error\n",myrank);
	  return FAIL;
	}
	idx=0;
	for (k=0; k<x2len; k++)
	  for (j=x1len-Ng1r; j<x1len; j++)
	    for (i=0; i<x0len; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }


//     fprintf(stdout,"  p%i: updating x2L bdry\n",myrank);  

    /////////////////////////
    // Update x2L boundaries
    {
      // open receive buffer for x2L boundary
      if (Nbors[2][0] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x2Lrecvsize, FDataType, MPI_Arg(Nbors[2][0]), 
		      msg_xch_x2l, MPI_COMM_WORLD, &id_recv_x2l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2L receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x2R boundary
      if (Nbors[2][1] != MPI_PROC_NULL) {
	idx=0;
	for (k=x2len-Ng2r-NborGhosts[2][1]; k<x2len-Ng2r; k++)
	  for (j=0; j<x1len; j++)
	    for (i=0; i<x0len; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];
      
	if (MPI_Isend(SendBuf, x2Rsendsize, FDataType, MPI_Arg(Nbors[2][1]), 
		      msg_xch_x2l, MPI_COMM_WORLD, &id_send_x2l) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2R send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x2L data, and update ghost cells
      if (Nbors[2][0] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x2l, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2L wait error\n",myrank);
	  return FAIL;
	}
	idx=0;
	for (k=0; k<Ng2l; k++)
	  for (j=0; j<x1len; j++)
	    for (i=0; i<x0len; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }


//     fprintf(stdout,"  p%i: updating x2R bdry\n",myrank);  

    /////////////////////////
    // Update x2R boundaries
    {
      // open receive buffer for x2R boundary
      if (Nbors[2][1] != MPI_PROC_NULL) {
	if (MPI_Irecv(RecvBuf, x2Rrecvsize, FDataType, MPI_Arg(Nbors[2][1]), 
		      msg_xch_x2r, MPI_COMM_WORLD, &id_recv_x2r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2R receive error\n",myrank);
	  return FAIL;
	}
      }
    
      // fill and send buffer for x2L boundary
      if (Nbors[2][0] != MPI_PROC_NULL) {
	idx=0;
	for (k=Ng2l; k<Ng2l+NborGhosts[2][0]; k++)
	  for (j=0; j<x1len; j++)
	    for (i=0; i<x0len; i++)
	      SendBuf[idx++] = mydata[(k*x1len + j)*x0len + i];

	if (MPI_Isend(SendBuf, x2Lsendsize, FDataType, MPI_Arg(Nbors[2][0]), 
		      msg_xch_x2r, MPI_COMM_WORLD, &id_send_x2r) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2L send error\n",myrank);
	  return FAIL;
	}
      }
    
      // wait for x2R data, and update ghost cells
      if (Nbors[2][1] != MPI_PROC_NULL) {
	if (MPI_Wait(&id_recv_x2r, &status) != 0) {
	  fprintf(stderr,"GravityBdryExchange p%i: x2R wait error\n",myrank);
	  return FAIL;
	}
	
	idx=0;
	for (k=x2len-Ng2r; k<x2len; k++)
	  for (j=0; j<x1len; j++)
	    for (i=0; i<x0len; i++)
	      mydata[(k*x1len + j)*x0len + i] = RecvBuf[idx++];
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

  // Delete file exchange buffers
  delete[] SendBuf;
  delete[] RecvBuf;

//   fprintf(stdout,"  p%i: returning from EnzoVector::exchange\n",myrank);  


  return SUCCESS;


#else
  // If MPI not used, just copy one boundary to ghosts of other boundary

  // first check that right and left procs are the same 
  // (should be either MPI_PROC_NULL or 0, but definitely the same)
  if (Nbors[0][0] != Nbors[0][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x0 neighbors");
    return FAIL;
  }
  if (Nbors[1][0] != Nbors[1][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x1 neighbors");
    return FAIL;
  }
  if (Nbors[2][0] != Nbors[2][1]) {
    fprintf(stderr,"EnzoVector::exchange error in x2 neighbors");
    return FAIL;
  }


  // iterate over all of the species, filling ghosts with opposite face data
  float *mydata;
  for (idat=0; idat<Nspecies; idat++) {

    // extract this variable from the vector array
    mydata = data[idat];
    
    // Update x0L boundaries
    if (Nbors[0][0] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<Ng0l; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[(k*x1len + j)*x0len + x0len-Ng0r-Ng0l + i];
    
    // Update x0R boundaries
    if (Nbors[0][1] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<Ng0r; i++)
	    mydata[(k*x1len + j)*x0len + x0len-Ng0r + i] = 
	      mydata[(k*x1len + j)*x0len + Ng0l + i];
    
    // Update x1L boundaries
    if (Nbors[1][0] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<Ng1l; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[(k*x1len + x1len-Ng0r-Ng0r + j)*x0len + i];
    
    // Update x1R boundaries
    if (Nbors[1][1] != MPI_PROC_NULL) 
      for (k=0; k<x2len; k++)
	for (j=0; j<Ng1r; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + x1len-Ng1r + j)*x0len + i] = 
	      mydata[(k*x1len + Ng1l + j)*x0len + i];
    
    // Update x2L boundaries
    if (Nbors[2][0] != MPI_PROC_NULL) 
      for (k=0; k<Ng2l; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[(k*x1len + j)*x0len + i] = 
	      mydata[((k + x2len-Ng2r-Ng2l)*x1len + j)*x0len + i];

    // Update x2R boundaries
    if (Nbors[2][1] != MPI_PROC_NULL) 
      for (k=0; k<Ng2r; k++)
	for (j=0; j<x1len; j++)
	  for (i=0; i<x0len; i++)
	    mydata[((k+x2len-Ng2r)*x1len + j)*x0len + i] = 
	      mydata[((k+Ng2l)*x1len + j)*x0len + i];
    
  }  // end for idat

  return SUCCESS;

#endif  // end if USE_MPI
}

/******************************************************************/
