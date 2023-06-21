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
/  SORT A LIST OF INTS AND DRAG A NUMBER OF FIELDS WITH IT
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/

#include <stdio.h>
#include "macros_and_parameters.h"

#define SWAP(A, n1, n2, tmp) {tmp = A[n1]; A[n1] = A[n2]; A[n2] = tmp;}

void QuickSortAndDrag(int List[], int left, int right, 
		      int NumberToDrag1, float *DragList1[], 
		      int NumberToDrag2, FLOAT *DragList2[])
{

  int i, n, last, leftright2, temp1;
  float temp2;
  FLOAT temp3;

  if (left >= right)
    return;

  leftright2 = (left + right)/2;

  SWAP(List, left, leftright2, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, leftright2, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, leftright2, temp3)

  last = left;

  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp1)
      for (n = 0; n < NumberToDrag1; n++)
	SWAP(DragList1[n], last, i, temp2)
      for (n = 0; n < NumberToDrag2; n++)
	SWAP(DragList2[n], last, i, temp3)
    }

  SWAP(List, left, last, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, last, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, last, temp3)

  QuickSortAndDrag(List, left  , last-1, 
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2);
  QuickSortAndDrag(List, last+1, right ,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2);
  
}

/* Float version of the above w/o FLOAT (should just have one version). */

void QuickSortAndDragFloat(float List[], int left, int right, int NumberToDrag,
                           float *DragList[])
{

  int i, n, last, leftright2;
  float temp;

  if (left >= right)
    return;

  leftright2 = (left + right)/2;

  SWAP(List, left, leftright2, temp)
  for (n = 0; n < NumberToDrag; n++)
    SWAP(DragList[n], left, leftright2, temp)

  last = left;

  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp)
      for (n = 0; n < NumberToDrag; n++)
	SWAP(DragList[n], last, i, temp)
    }

  SWAP(List, left, last, temp)
  for (n = 0; n < NumberToDrag; n++)
    SWAP(DragList[n], left, last, temp)

  QuickSortAndDragFloat(List, left  , last-1, NumberToDrag, DragList);
  QuickSortAndDragFloat(List, last+1, right , NumberToDrag, DragList);

}
