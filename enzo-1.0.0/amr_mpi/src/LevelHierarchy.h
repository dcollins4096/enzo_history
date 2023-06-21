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
/  LEVEL HIERARCHY STRUCTURE AND ROUTINES
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

struct LevelHierarchyEntry
{
  LevelHierarchyEntry *NextGridThisLevel;  /* next entry on this level */
  grid                *GridData;           /* pointer to this entry's grid */
  HierarchyEntry      *GridHierarchyEntry; /* pointer into hierarchy */
};
