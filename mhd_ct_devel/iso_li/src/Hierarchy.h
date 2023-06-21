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
/  HIERARCHY STRUCTURE
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

struct HierarchyEntry
{
  HierarchyEntry *NextGridThisLevel; /* pointer to the next grid on level */
  HierarchyEntry *NextGridNextLevel; /* pointer to first child of this grid */
  HierarchyEntry *ParentGrid;        /* pointer to this grid's parent */
  grid           *GridData;          /* pointer to this grid's data */
};
