/****************************************************************************/
/* This file contains declarations used by the FastSiblingLocator routines. */
/****************************************************************************/

/* This is the declaration for a link in the link list. */

struct ChainingMeshLink {
  grid *GridData;
  ChainingMeshLink *NextLink;
};

/* This is the chaining mesh structure declarataion. */

struct ChainingMeshStructure {
  int Rank;
  int Dimension[MAX_DIMENSION];
  FLOAT LeftEdge[MAX_DIMENSION];
  FLOAT RightEdge[MAX_DIMENSION];
  FLOAT CellSize[MAX_DIMENSION];
  boundary_type BoundaryType[MAX_DIMENSION];
  ChainingMeshLink **HeadOfChain;
  ChainingMeshLink *OutsideMesh;
};

/* This is the declaration for the list of sibling grids. */

struct SiblingGridList {
  int NumberOfSiblings;
  grid **GridList;
};
