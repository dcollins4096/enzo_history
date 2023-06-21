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
//======================================================================
//
// File:        performance.C
//
// Description: Performance-related  declarations
//
//----------------------------------------------------------------------
//
// Namespaces:  jb
//
//----------------------------------------------------------------------
//
// James Bordner
// UCSD
//
//======================================================================

namespace jb {
  int mpi_send_start,mpi_send_stop;
  int mpi_recv_start,mpi_recv_stop;
  int mpi_sendrecv_start,mpi_sendrecv_stop;
  int mpi_barrier_start,mpi_barrier_stop;
  int mpi_gather_start,mpi_gather_stop;
  int mpi_reduce_start,mpi_reduce_stop;
};
