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
#ifndef MESSAGE_H
#define MESSAGE_H
void c_error (char *, int );
void c_warning (char *, int );

#define ERROR_MESSAGE c_error(__FILE__,__LINE__);
#define WARNING_MESSAGE c_warning(__FILE__,__LINE__);
#define DEBUG_MESSAGE c_warning(__FILE__,__LINE__);

      
#endif

