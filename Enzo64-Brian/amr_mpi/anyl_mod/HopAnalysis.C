/***********************************************************************
/
/  PERFORMS A "HOP" CLUSTERING ANALYSIS (astro-ph/9712200)
/
/  written by: Greg Bryan
/  date:       February, 1999
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
//
 
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
 
#include "hdf4.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#undef DEFINE_STORAGE
extern "C" {
#include "kd.h"
}
 
/* function prototypes */
 
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		    ExternalBoundary *Exterior);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int  DepositParticleMassField(HierarchyEntry *Grid, float RequestTime = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData,
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid,
				      TopGridData *MetaData,
				      LevelHierarchyEntry *LevelArray[],
				      int level);
int CommunicationInitialize(int *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);
extern "C" void hop_main(KD kd);
extern "C" void regroup_main(float dens_outer);
 
main(int argc, char *argv[])
{
  CommunicationInitialize(&argc, &argv);
 
  /* Main declarations */
 
  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
 
  /* Initialize */
 
  float HopDensityThreshold    = 160;   // reproduces FOF b=0.2
  int   UseBaryons             = 0;
  debug                        = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];
  FLOAT RegionLeft[MAX_DIMENSION], RegionRight[MAX_DIMENSION];
  int level, dim, i, j, part;
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    RegionLeft[dim] = RegionRight[dim] = FLOAT_UNDEFINED;
 
  /* --------------------------------------------------------------- */
  /* Interpret command-line arguments. */
 
  char c;
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {
 
	/* get beginning of region selection (float coordinate). */
 
      case 'b':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionLeft[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Begin coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* debug */
 
      case 'd':
	debug = TRUE;
	break;
 
	/* get finish of region selection (float coordinate). */
 
      case 'f':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"PSYM, &RegionRight[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Finish coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* gas and dm */
 
      case 'g':
	UseBaryons = 1;
	break;
 
	/* get density threshold. */
 
      case 't':
	if (--argc > 0) {
	  sscanf((*++argv), "%"PSYM, &HopDensityThreshold);
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;
 
	/* Unknown */
 
      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	my_exit(EXIT_FAILURE);
	
      } // end of switch(c)
 
  /* Error check for number of parameters, and set parameter file. */
 
  if (argc != 1) {
    fprintf(stderr, "enzohop [-b #] [-f #] [-t #] [-g] [-d] amr_file\n");
    fprintf(stderr, "  -b)egin region\n");
    fprintf(stderr, "  -f)inish region\n");
    fprintf(stderr, "  -t)hreshold for hop (default 160)\n");
    fprintf(stderr, "  -g)as particles also used (normally just dm)\n");
    fprintf(stderr, "  -d)ebug\n");
    my_exit(FAIL);
  }
  ParameterFile = argv[0];
 
  /* Read the saved file. */
 
  SetDefaultGlobalValues(MetaData);
  if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior) == FAIL) {
    fprintf(stderr, "Error in ParameterFile %s.\n", ParameterFile);
    my_exit(EXIT_FAILURE);
  }
  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels
 
  /* Set the Cell width of the root grid. */
 
  float BaseRadius = (DomainRightEdge[0] - DomainLeftEdge[0])/
      float(MetaData.TopGridDims[0]);
 
  /* If undefined, set parameters. */
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (RegionLeft[dim] == FLOAT_UNDEFINED)
      RegionLeft[dim] = DomainLeftEdge[dim];
    if (RegionRight[dim] == FLOAT_UNDEFINED)
      RegionRight[dim] = DomainRightEdge[dim];
  }
 
  if (debug)
    printf("HopAnalysis region: Left = %"GSYM" %"GSYM" %"GSYM"   Right = %"GSYM" %"GSYM" %"GSYM"\n",
	   RegionLeft[0], RegionLeft[1], RegionLeft[2],
	   RegionRight[0], RegionRight[1], RegionRight[2]);
 
  /* Initialize Particle List info */
 
  ListOfParticles *ListOfParticlesHead[2];
  int TotalNumberOfParticles[2];
  for (i = 0; i < 2; i++) {
    ListOfParticlesHead[i] = NULL;
    TotalNumberOfParticles[i] = 0;
  }
 
  /* --------------------------------------------------------------- */
  /* Loop over all the levels, and collect particles */
 
  printf("Collecting particles...\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
 
    /* If SelfGravity, set all the particle mass fields. */
 
    LevelHierarchyEntry *Temp = LevelArray[level];
#ifdef UNUSED
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }
#endif /* UNUSED */
 
    /* Loop over all the grids. */
 
    Temp = LevelArray[level];
    while (Temp != NULL) {
 
      /* Allocate a new ListOfParticles for this grid. */
 
      for (i = 0; i < 2; i++) {
	ListOfParticles *Temp = ListOfParticlesHead[i];
	ListOfParticlesHead[i] = new ListOfParticles;
	ListOfParticlesHead[i]->NextList = Temp;
      }
 
      /* Set particle density. */
 
#ifdef UNUSED
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData,
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
#endif /* UNUSED */
 
      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */
 
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
 
      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */
 
      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
        					 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
 
      /* Generate particle list for this grid. */
 
      Temp->GridData->OutputAsParticleData(RegionLeft, RegionRight,
					   ListOfParticlesHead, BaseRadius);
      if (UseBaryons)
	TotalNumberOfParticles[0] += ListOfParticlesHead[0]->NumberOfParticles;
      else {
	if (ListOfParticlesHead[0]->NumberOfParticles > 0) {
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    delete [] ListOfParticlesHead[0]->ParticlePosition[dim];
	    delete [] ListOfParticlesHead[0]->ParticleVelocity[dim];
	  }
	  delete [] ListOfParticlesHead[0]->ParticleRadius;
	  for (j = 0; j < ListOfParticlesHead[0]->NumberOfValues; j++)
	    delete [] ListOfParticlesHead[0]->ParticleValue[j];
	}
      }
      TotalNumberOfParticles[1] += ListOfParticlesHead[1]->NumberOfParticles;
 
      /* Delete grid. */
 
      delete Temp->GridData;
 
      /* Next grid on this level. */
 
      Temp = Temp->NextGridThisLevel;
 
    } // end loop over grids
 
  } // end loop over levels
 
  /* --------------------------------------------------------------- */
  /* Convert particles into kd structure. */
 
  if (debug)
    printf("TotalNumberOfParticles = %"ISYM" %"ISYM"\n", TotalNumberOfParticles[0],
	   TotalNumberOfParticles[1]);
 
  /* Initialize the kd structure. */
 
  if (sizeof(float) != 4) {
    fprintf(stderr, "error: hop is hard-coded for 4-byte floats.\n");
    my_exit(EXIT_FAILURE);
  }
 
  KD kd;
  int nBucket = 16, kdcount = 0;
  kdInit(&kd, nBucket);
  kd->nActive = TotalNumberOfParticles[0]+TotalNumberOfParticles[1];
  kd->p = new PARTICLE[kd->nActive];
  if (kd->p == NULL) {
    fprintf(stderr, "failed allocating particles.\n");
    my_exit(EXIT_FAILURE);
  }
 
  float DensityFactor = 1.0/(2.78e11*OmegaMatterNow*POW(HubbleConstantNow, 2) *
			     POW(ComovingBoxSize/HubbleConstantNow, 3));
  ListOfParticles FullList[2];
 
  for (i = 0; i < 2; i++) {
 
    if (TotalNumberOfParticles[i] > 0) {
 
      /* Allocate a field for these particles. */
 
      FullList[i].NumberOfValues = ListOfParticlesHead[i]->NumberOfValues;
      if (debug)
	printf("NumberOfValues[%"ISYM"] = %"ISYM"\n", i, FullList[i].NumberOfValues);
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	FullList[i].ParticlePosition[dim] = new
	  float[TotalNumberOfParticles[i]];
	FullList[i].ParticleVelocity[dim] = new
	  float[TotalNumberOfParticles[i]];
      }
      FullList[i].ParticleRadius = new float[TotalNumberOfParticles[i]];
      for (j = 0; j < FullList[i].NumberOfValues; j++)
	FullList[i].ParticleValue[j] = new float[TotalNumberOfParticles[i]];
      if (i == 1)
	FullList[i].ParticleIndex = new int[TotalNumberOfParticles[i]];
 
      /* Copy grid lists into the full list. */
 
      ListOfParticles *Temp = ListOfParticlesHead[i];
      int count = 0;
      while (Temp != NULL) {
//     printf("i=%"ISYM" part=%"ISYM" count=%"ISYM"\n", i, Temp->NumberOfParticles, count);
	for (dim = 0; dim < MetaData.TopGridRank; dim++)
	  for (part = 0; part < Temp->NumberOfParticles; part++) {
	    FullList[i].ParticlePosition[dim][count+part] =
	      Temp->ParticlePosition[dim][part];
	    FullList[i].ParticleVelocity[dim][count+part] =
	      Temp->ParticleVelocity[dim][part];
	  }
	for (part = 0; part < Temp->NumberOfParticles; part++)
	  FullList[i].ParticleRadius[count+part] = Temp->ParticleRadius[part];
	for (j = 0; j < FullList[i].NumberOfValues; j++)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList[i].ParticleValue[j][count+part] =
	      Temp->ParticleValue[j][part];
	if (i == 1)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList[i].ParticleIndex[count+part] = Temp->ParticleIndex[part];
 
	/* Copy positions into kd structure. */
 
	int izero = SelfGravity - 1;  // always zero
	for (part = 0; part < Temp->NumberOfParticles; part++) {
	  kd->p[kdcount].r[0] = Temp->ParticlePosition[0][part];
	  kd->p[kdcount].r[1] = Temp->ParticlePosition[1][part];
	  kd->p[kdcount].r[2] = Temp->ParticlePosition[2][part];
	  for (dim = 0; dim < MetaData.TopGridRank; dim++)
	    if (kd->p[kdcount].r[dim] !=
		kd->p[kdcount].r[dim+izero]) {
	      printf("warn: %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"GSYM"\n", kdcount, part, i, dim,
		     kd->p[kdcount].r[dim]);
	      kd->p[kdcount].r[dim] = 0;
	    }
	  kd->p[kdcount].fMass = Temp->ParticleValue[0][part]*DensityFactor;
	  kdcount++;
	}
 
	/* Delete this grid's particles. */
 
	if (Temp->NumberOfParticles > 0) {
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    delete [] Temp->ParticlePosition[dim];
	    delete [] Temp->ParticleVelocity[dim];
	  }
	  delete [] Temp->ParticleRadius;
	  if (i == 1)
	    delete [] Temp->ParticleIndex;
	  for (j = 0; j < FullList[i].NumberOfValues; j++)
	    delete [] Temp->ParticleValue[j];
	}
 
	count += Temp->NumberOfParticles;
	Temp = Temp->NextList;
      }
 
    } // end: if (TotalNumberOfParticles[i] > 0)
 
  } // end: loop over dark matter/baryon particle lists
 
  /* --------------------------------------------------------------- */
  /* Call hop. */
 
  fprintf(stderr, "Calling hop...\n");
  hop_main(kd);
 
  fprintf(stderr, "Calling regroup...\n");
  regroup_main(HopDensityThreshold);
 
 
  /* --------------------------------------------------------------- */
  /* Read the group membership and compute group properties. */
 
  FILE *fptr;
  if ((fptr = fopen("zregroup.tag", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output zregroup.hop\n");
    my_exit(EXIT_FAILURE);
  }
 
  int nActive, nGroups;
  fread(&nActive, 4, 1, fptr);
  fread(&nGroups, 4, 1, fptr);
  printf("nActive = %"ISYM"(=%"ISYM")   nGroups = %"ISYM"\n", nActive, kd->nActive, nGroups);
 
  /* Allocate space and read group memberships for the particles. */
 
  int *GroupID = new int[nActive];
  if (fread(GroupID, 4, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file zregroup.hop\n");
    my_exit(EXIT_FAILURE);
  }
  fclose(fptr);
 
  /* Allocate space and read group memberships for the particles. */
 
  float *Density = new float[nActive];
  if ((fptr = fopen("output_hop.den", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output output_hop.den\n");
    my_exit(EXIT_FAILURE);
  }
  fread(&nActive, 4, 1, fptr);
  if (fread(Density, 4, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file output_hop.den\n");
    my_exit(EXIT_FAILURE);
  }
  fclose(fptr);
 
  /* Allocate and initialize group information. */
 
  int const NumberOfGroupProperties = 10;
  float *GroupProperties[NumberOfGroupProperties];
  for (i = 0; i < NumberOfGroupProperties; i++) {
    GroupProperties[i] = new float[nGroups];
    for (j = 0; j < nGroups; j++)
      GroupProperties[i][j] = 0;
  }
 
  /* Loop over particles, adding to group properties. */
 
  int i1, i2;
  float Luminosity;
  for (i = 0; i < nActive; i++)
    if ((j = GroupID[i]) >= 0) {
 
      /* Convert GroupID into twopart index for FullList. */
 
      i1 = (i < TotalNumberOfParticles[0]) ? 0 : 1;
      i2 = i - i1*TotalNumberOfParticles[0];
 
      /* Total mass. */
 
      GroupProperties[0][j] += FullList[i1].ParticleValue[0][i2];
 
      /* Free-free emission. */
 
      if (i1 == 0)
	GroupProperties[1][j] += (Luminosity =
				  FullList[i1].ParticleValue[1][i2] *
				  FullList[i1].ParticleValue[0][i2] *
				  sqrt(FullList[i1].ParticleValue[2][i2]));
 
      /* Luminosity-weighted Temperature. */
 
      if (i1 == 0)
	GroupProperties[2][j] += Luminosity*FullList[i1].ParticleValue[2][i2];
 
      /* Number of points. */
 
      GroupProperties[3][j] += 1.0;
 
      /* Max density (and position). */
 
      if (Density[i] > GroupProperties[4][j]) {
	GroupProperties[4][j] = Density[i];
	GroupProperties[5][j] = FullList[i1].ParticlePosition[0][i2];
	GroupProperties[6][j] = FullList[i1].ParticlePosition[1][i2];
	GroupProperties[7][j] = FullList[i1].ParticlePosition[2][i2];
      }
 
    }
 
  /* Normalize temperature. */
 
  for (j = 0; j < nGroups; j++)
    if (GroupProperties[1][j] > 0)
      GroupProperties[2][j] /= GroupProperties[1][j];
 
  /* Output group properties. */
 
  if ((fptr = fopen("HopAnalysis.out", "w")) == NULL) {
    fprintf(stderr, "Error opening regroup output HopAnalysis.out\n");
    my_exit(EXIT_FAILURE);
  }
 
  fprintf(fptr, "#Group     Mass      # part     max dens     x        y       z\n");
  for (j = 0; j < nGroups; j++) {
    fprintf(fptr, "%"ISYM"      ", j);
    for (i = 0; i < 8; i++)
      if (i < 1 || i > 2)
	fprintf(fptr, " %.9"GSYM" ", GroupProperties[i][j]);
    fprintf(fptr, "\n");
  }
 
  fclose(fptr);
 
  my_exit(EXIT_SUCCESS);
}
 
void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}
