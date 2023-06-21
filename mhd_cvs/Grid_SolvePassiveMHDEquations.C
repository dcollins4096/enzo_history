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
/  GRID CLASS (Evolve Magentic Fields Passively)
/
/  written by: Hao Xu
/  date:       May, 2006
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

#include "performance.h"
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "pout.h"
/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);
extern "C" void FORTRAN_NAME(ppm_de)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual,
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
			  int colindex[]);
extern "C" void FORTRAN_NAME(ppm_lr)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual, 
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
			  int colindex[]);
extern "C" void FORTRAN_NAME(zeus_main)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge, float *C1, float *C2,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual, int *igamfield, 
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
                            int colindex[], int *bottom,
			  float *minsupecoef);

extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
                                        float *ex, float *ey, float *ez,
                                        float *dx, float *dy, float *dz,
                                        int *idim, int *jdim, int *kdim,
                                        int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                        float *dt, MHD_Centering *method);


extern "C" void FORTRAN_NAME(create_e)(float *fx1, float *fy1, float *fz1,
                                       float *fx2, float *fy2, float *fz2,
                                       float *ex, float *ey, float *ez,
                                       int *idim, int *jdim, int *kdim,   
                                       int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                       float * dt, int * ProjectE);

extern "C" void FORTRAN_NAME(divb_rj)(float *fx1, float *fy1, float *fz1, 
                                      float *fx2, float *fy2, float *fz2,
                                      float *bxb, float *byb, float *bzb,
                                      float *bxc, float *byc, float *bzc,
                                      float *dx, float *dy, float *dz, float *dt,
                                      int *nx, int *ny, int *nz);

extern "C" void FORTRAN_NAME (passivemhd_create_e)(float *density, float *newvx, float *newvy, float *newvz,
         float *oldvx, float *oldvy, float *oldvz,
         float *pre,
         float *bxf, float *byf, float *bzf, float *oldbxf, float *oldbyf, float *oldbzf,
         float *ex, float *ey, float *ez,
         float * dx, float *dy, float *dz, int *idim, int *jdim, int *kdim, float *dtIn, float *olddt, int *ProjectE,
         int *i1, int *i2, int *j1,int *j2,int *k1, int *k2, FLOAT *a, int *biermann,
         float *c, float *mh, float *e, float *chi
      );

extern "C" void FORTRAN_NAME (passivemhd_computer_b)(float *density, float *newvx, float *newvy, float *newvz,
         float *oldvx, float *oldvy, float *oldvz,
         float *pre,
         float *bx, float *by, float *bz,
         float * dx, float *dy, float *dz, int *idim, int *jdim, int *kdim, float *dtIn, int *ProjectE,
         int *i1, int *i2, int *j1,int *j2,int *k1, int *k2, FLOAT *a, int *biermann
      );

#ifdef BIERMANN
int MHDCosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                       float *TemperatureUnits, float *TimeUnits,
                       float *VelocityUnits, FLOAT Time,
                       float *BFieldUnits);
int  CosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                       float *TemperatureUnits, float *TimeUnits,
                       float *VelocityUnits, FLOAT Time);
#endif


#ifdef PASSIVE
int grid::SolvePassiveMHDEquations(int CycleNumber, int NumberOfSubgrids, 
			      fluxes *SubgridFluxes[], int level, int grid)
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;


  //dcc
  fprintf(stderr,"=== GPFS SolvePassiveEquations CycleNumber = %d Level = %d proc= %d dtFixed = %15.12f ===\n",
  CycleNumber, level, MyProcessorNumber, dtFixed);
  Pout("shd: level, grid", level, grid);
  PoutF("shd: left", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
  PoutF("shd: rigt", GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
  if( level == 1 && grid == 1 ) {
    fprintf(stderr,"gonna stop.\n");
  }
  // /dcc

  this->DebugCheck("SolveHydroEquations");

  if (NumberOfBaryonFields > 0) {

    /* initialize */

    int dim, i, idim, j, jdim, field, size, subgrid, n;
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, coloff[MAX_COLOR];
    long_int GridGlobalStart[MAX_DIMENSION];
    FLOAT a = 1, dadt;

    int k, face, coord, axis;
    int  Dim[3]={1,1,1};
    int Findex, Cindex;
    int * FluxDims[3][3];
#define fluxe(dim,coord,face,end,sub) FluxExtents[dim+3*(coord+3*(face+2*(end+2*sub)))]
int * FluxExtents = new int[3*3*2*2*NumberOfSubgrids];

    float *colourpt = NULL;

//#ifdef BIERMANN
    // Biermann battery constants in cgs units, convert to enzo units later
    float speedoflight = 3.0e10;
    float hydrogenmass = 1.6733e-24;
    float electroncharge = 4.803e-10;
    float chi=1.0;
//#endif //BIERMANN    


if(Passive_Used ==1){
    for(dim=0;dim<3;dim++)
    for(coord=0;coord<3;coord++)
      for(face=0;face<2;face++)
	for(n=0;n<2;n++){
	  FluxDims[dim][coord] = new int[NumberOfSubgrids];
	  for(subgrid=0;subgrid<NumberOfSubgrids;subgrid++){
	    fluxe(dim,coord,face,n,subgrid)=0;
	    FluxDims[dim][coord][subgrid]=0;
	  }
	}
}//Passive==1

    /* Compute size (in floats) of the current grid. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
    
    /* Find fields: density, total energy, velocity1-3. */

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
      return FAIL;
    }

    /* If multi-species being used, then treat them as colour variables
       (note: the solver has been modified to treat these as density vars). */

    int NumberOfColours = 0, ColourNum;

    if (MultiSpecies > 0) {

      NumberOfColours = 6 + 3*(MultiSpecies-1);

      if ((ColourNum = 
	   FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
	fprintf(stderr, "Could not find ElectronDensity.\n");
	return FAIL;
      }

      /* Set Offsets from the first field (here assumed to be the electron
	 density) to the other multi-species fields. */

      colourpt = BaryonField[ColourNum];

      for (i = 0; i < NumberOfColours; i++)
	coloff[i] = BaryonField[ColourNum+i] - colourpt;

    }
 
    /* Add metallicity as a colour variable. */

    int MetalNum;

    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
      if (colourpt == NULL) {
	colourpt = BaryonField[MetalNum];
	ColourNum = MetalNum;
      }

      coloff[NumberOfColours++] = BaryonField[MetalNum  ] - colourpt;
      coloff[NumberOfColours++] = BaryonField[MetalNum+1] - colourpt;
      coloff[NumberOfColours++] = BaryonField[MetalNum+2] - colourpt;

    }

    /* Get easy to handle pointers for each variable. */

    float *density     = BaryonField[DensNum];
    float *totalenergy = BaryonField[TENum];
    float *gasenergy   = BaryonField[GENum];

    /* Velocity1 must exist, but if 2 & 3 aren't present, then create blank
       buffers for them (since the solver needs to advect something). */

    float *velocity1, *velocity2, *velocity3;
    velocity1 = BaryonField[Vel1Num];

    if (GridRank > 1)
      velocity2 = BaryonField[Vel2Num];    
    else {
      velocity2 = new float[size];
      for (i = 0; i < size; i++)
        velocity2[i] = 0.0;
    }

    if (GridRank > 2)
      velocity3 = BaryonField[Vel3Num];    
    else {
      velocity3 = new float[size];
      for (i = 0; i < size; i++)
        velocity3[i] = 0.0;
    }


   //
   // Allocate Electric Field.
   //   
       
   for(field=0;field<3;field++){   
       
     if(ElectricField[field] != NULL ) {
       delete [] ElectricField[field];
     }
     
     ElectricField[field] = new float[ElectricSize[field]];
     for(i=0;i<ElectricSize[field]; i++) ElectricField[field][i] =
0.0+level*(1+CycleNumber);
     
       
     if(MagneticField[field]==NULL){
       fprintf(stderr, "========== Solve MHD create Magnetic Field\nShit,that's not good.\n");
       return FAIL;
     }  
  
   }


    /* Allocate space for magnetic fluxes */

   float *MagneticFlux[3][2];
//each magnetic field has components from the flux of two 'other' axis
  
  for( field=0;field<3;field++)
    for(axis=0;axis<2;axis++)
      MagneticFlux[field][axis] = new float[MagneticSize[field]];
       
       
  for(field=0;field<3;field++)
    for( i=0;i<MagneticSize[field]; i++){
      MagneticFlux[field][0][i] = 0.0;
      MagneticFlux[field][1][i] = 0.0;
    }



    /* Allocate temporary space for Zeus_Hydro. */

    float *zeus_temp;
    if (HydroMethod == Zeus_Hydro)
      zeus_temp = new float[size];
    int LowestLevel = (level > MaximumRefinementLevel-1) ? TRUE : FALSE;

    /* Determine if Gamma should be a scalar or a field. */

    int UseGammaField = FALSE;
    float *GammaField;
    if (HydroMethod == Zeus_Hydro && MultiSpecies > 1) {
      UseGammaField = TRUE;
      GammaField = new float[size];
      if (this->ComputeGammaField(GammaField) == FAIL) {
	fprintf(stderr, "Error in grid->ComputeGammaField.\n");
	return FAIL;
      }
    } else {
      GammaField = new float[1];
      GammaField[0] = Gamma;
    }

    /* Set minimum support. */

    float MinimumSupportEnergyCoefficient = 0;
    if (UseMinimumPressureSupport == TRUE && level > MaximumRefinementLevel-1)
      if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
	fprintf(stderr, "Error in grid->SetMinimumSupport,\n");
	return FAIL;
      }

    /* allocate space for fluxes */

    for (i = 0; i < NumberOfSubgrids; i++) {
      for (dim = 0; dim < GridRank; dim++)  {

	/* compute size (in floats) of flux storage */

        size = 1;
        for (j = 0; j < GridRank; j++)
          size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
                  SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

	/* set unused dims (for the solver, which is hardwired for 3d). */

        for (j = GridRank; j < 3; j++) {
          SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
        }

	/* Allocate space (if necessary). */

        for (field = 0; field < NumberOfBaryonFields; field++) {
	  if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new float[size];
	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->RightFluxes[field][dim] = new float[size];
	  for (n = 0; n < size; n++) {
	    SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
	    SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
	  }
        }

	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	     field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

      }  // next dimension

      /* make things pretty */

      for (dim = GridRank; dim < 3; dim++)
        for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

    } // end of loop over subgrids

    /* compute global start index for left edge of entire grid 
       (including boundary zones) */

    for (dim = 0; dim < GridRank; dim++) {
      GridGlobalStart[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
	GridStartIndex[dim];
    }

    /* fix grid quantities so they are defined to at least 3 dims */
   
    for (i = GridRank; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridVelocity[i]    = 0.0;
      GridGlobalStart[i] = 0;
    }

    /* allocate temporary space for solver (enough to fit 31 of the largest
       possible 2d slices plus 4*NumberOfColours). */

    int tempsize = max(max(GridDimension[0]*GridDimension[1], 
                           GridDimension[1]*GridDimension[2]),
		           GridDimension[2]*GridDimension[0]  );
    float *temp = new float[tempsize*(31+NumberOfColours*4)];

    /* create and fill in arrays which are easiler for the solver to
       understand. */

    int *leftface  = new int[NumberOfSubgrids*3*(18+2*NumberOfColours)];
    int *rightface = leftface + NumberOfSubgrids*3*1;
    int *istart    = leftface + NumberOfSubgrids*3*2;
    int *jstart    = leftface + NumberOfSubgrids*3*3;
    int *iend      = leftface + NumberOfSubgrids*3*4;
    int *jend      = leftface + NumberOfSubgrids*3*5;
    int *dindex    = leftface + NumberOfSubgrids*3*6;
    int *Eindex    = leftface + NumberOfSubgrids*3*8;
    int *uindex    = leftface + NumberOfSubgrids*3*10;
    int *vindex    = leftface + NumberOfSubgrids*3*12;
    int *windex    = leftface + NumberOfSubgrids*3*14;
    int *geindex   = leftface + NumberOfSubgrids*3*16;
    int *colindex  = leftface + NumberOfSubgrids*3*18;
    float *standard = NULL;
    if (NumberOfSubgrids > 0) standard = SubgridFluxes[0]->LeftFluxes[0][0];

    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
      for (dim = 0; dim < GridRank; dim++) {

        /* Set i,j dimensions of 2d flux slice (this works even if we 
           are in 1 or 2d) the correspond to the dimensions of the global
           indicies.  I.e. for dim = 0, the plane is dims 1,2
                           for dim = 1, the plane is dims 0,2
                           for dim = 2, the plane is dims 0,1 */

	idim = (dim == 0) ? 1 : 0;
	jdim = (dim == 2) ? 1 : 2;

        /* Set the index (along the dimension perpindicular to the flux
           plane) of the left and right flux planes.  The index is zero
           based from the left side of the entire grid. */

	leftface[subgrid*3+dim] = 
	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][dim] - 
	    GridGlobalStart[dim];
	rightface[subgrid*3+dim] = 
	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][dim] -
	    GridGlobalStart[dim];   // (+1 done by fortran code)

        /* set the start and end indicies (zero based on entire grid)
           of the 2d flux plane. */

	istart[subgrid*3+dim] =
	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][idim] -
	    GridGlobalStart[idim];
	jstart[subgrid*3+dim] =
	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][jdim] -
	    GridGlobalStart[jdim];
	iend[subgrid*3+dim] =
	  SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][idim] -
	    GridGlobalStart[idim];
	jend[subgrid*3+dim] =
	  SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][jdim] -
	    GridGlobalStart[jdim];

        /* Compute offset from the standard pointer to the start of 
           each set of flux data.  This is done to compensate for
           fortran's inability to handle arrays of pointers or structs.
	   NOTE: This pointer arithmatic is illegal; some other way should
	   be found to do it (like write higher level ppm stuff in c++). */

	dindex[subgrid*6+dim*2] = 
	  SubgridFluxes[subgrid]->LeftFluxes[DensNum][dim] - standard;
	dindex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[DensNum][dim] - standard;
	Eindex[subgrid*6+dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[TENum][dim] - standard;
	Eindex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[TENum][dim] - standard;
	uindex[subgrid*6+dim*2] =
	  SubgridFluxes[subgrid]->LeftFluxes[Vel1Num][dim] - standard;
	uindex[subgrid*6+dim*2+1] =
	  SubgridFluxes[subgrid]->RightFluxes[Vel1Num][dim] - standard;

	if (GridRank > 1) {
          vindex[subgrid*6+dim*2] =
	    SubgridFluxes[subgrid]->LeftFluxes[Vel2Num][dim] - standard;
          vindex[subgrid*6+dim*2+1] =
	    SubgridFluxes[subgrid]->RightFluxes[Vel2Num][dim] - standard;
        }
	if (GridRank > 2) {
          windex[subgrid*6+dim*2] =
	    SubgridFluxes[subgrid]->LeftFluxes[Vel3Num][dim] - standard;
          windex[subgrid*6+dim*2+1] =
	    SubgridFluxes[subgrid]->RightFluxes[Vel3Num][dim] - standard;
        }

        if (DualEnergyFormalism) {
          geindex[subgrid*6+dim*2] =
            SubgridFluxes[subgrid]->LeftFluxes[GENum][dim] - standard;
          geindex[subgrid*6+dim*2+1] =
            SubgridFluxes[subgrid]->RightFluxes[GENum][dim] - standard;
        }

	for (i = 0; i < NumberOfColours; i++) {
	  colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2] =
	    SubgridFluxes[subgrid]->LeftFluxes[ColourNum+i][dim] - standard;
	  colindex[i*NumberOfSubgrids*6 + subgrid*6 + dim*2 + 1] =
	    SubgridFluxes[subgrid]->RightFluxes[ColourNum+i][dim] - standard;
	}

      }

    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */

    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	  == FAIL) {
	fprintf(stderr, "Error in CsomologyComputeExpansionFactors.\n");
	return FAIL;
      }

/* Create a cell width array to pass (and convert to absolute coords). */

     float *CellWidthTemp[MAX_DIMENSION];
     for (dim = 0; dim < MAX_DIMENSION; dim++) {
       CellWidthTemp[dim] = new float[GridDimension[dim]];
       for (i = 0; i < GridDimension[dim]; i++)
       if (dim < GridRank)
         CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
       else
         CellWidthTemp[dim][i] = 1.0;
     }



if(Passive_Used==1){
#ifdef OLD_CENTER
   
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField, first call \n");
    return FAIL;
  }
#endif //OLD_CENTER
}//Passive_Used==1


    /* Prepare Gravity. */

    int GravityOn = 0, FloatSize = sizeof(float);
    if (SelfGravity || UniformGravity || PointSourceGravity)
      GravityOn = 1;

    /* call a Fortran routine to actually compute the hydro equations
       on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */

#ifdef BIERMANN    
   int Biermann = 1;
   
    /* Compute Units. */
                                                                                                                                                             
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
        VelocityUnits = 1, BFieldUnits = 1;
                                                                                                                                                             
  if(ComovingCoordinates){
    if(Passive_Used){
        if (MHDCosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                          &TimeUnits, &VelocityUnits, Time,&BFieldUnits) == FAIL) {
      fprintf(stderr, "Error in MHD CosmologyGetUnits.\n");
      return FAIL;
    }
   }else{
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                          &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
    }// MHD_Used
    }//ComovingCoordinates

   /* Transform speed of light, hydrogen mass and electron change in units of the ENZO */

    electroncharge *= TimeUnits*BFieldUnits/(speedoflight*DensityUnits*pow(LengthUnits,3));
    speedoflight /= VelocityUnits;
    hydrogenmass /= DensityUnits*pow(LengthUnits,3);
   

#else
   int Biermann = 0;
#endif 
if(0){        
      int UseDT = 1;
//      MHD_ComputeCurrent();
      FORTRAN_NAME(passivemhd_create_e)(density,OldBaryonField[Vel1Num],OldBaryonField[Vel2Num], OldBaryonField[Vel3Num],
          OldBaryonField[Vel1Num],OldBaryonField[Vel2Num], OldBaryonField[Vel3Num],
          BaryonField[GENum],
          MagneticField[0], MagneticField[1], MagneticField[2],
          OldMagneticField[0], OldMagneticField[1], OldMagneticField[2],
          ElectricField[0], ElectricField[1], ElectricField[2],
          CellWidth[0], CellWidth[1], CellWidth[2],
          GridDimension, GridDimension +1, GridDimension +2,
          &dtFixed, &dtFixed, &UseDT,
          GridStartIndex, GridEndIndex,
          GridStartIndex+1, GridEndIndex+1,
          GridStartIndex+2, GridEndIndex+2,
          &a,&Biermann,
          &speedoflight,&hydrogenmass, &electroncharge, &chi
      );
  
    // Update Magnetic Fields
    float dtUsed = 1.0;
  
    FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
                                ElectricField[0], ElectricField[1], ElectricField[2],
                                CellWidth[0], CellWidth[1], CellWidth[2],
                                GridDimension, GridDimension +1, GridDimension +2,
                                GridStartIndex, GridEndIndex, 
                                GridStartIndex+1, GridEndIndex+1,
                                GridStartIndex+2, GridEndIndex+2,
                                &dtUsed, &MHD_CenteringMethod);
          
    for(i=0;i<MagneticSize[0];i++) 
       MagneticField[0][i]=0.5*(OldMagneticField[0][i]+MagneticField[0][i]);
    for(i=0;i<MagneticSize[1];i++)                                          
       MagneticField[1][i]=0.5*(OldMagneticField[1][i]+MagneticField[1][i]);
    for(i=0;i<MagneticSize[2];i++)                                          
       MagneticField[2][i]=0.5*(OldMagneticField[2][i]+MagneticField[2][i]);

            
#ifdef OLD_CENTER
              
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField, first call \n");
    return FAIL;
  }
#endif //OLD_CENTER
          
                 
          
}



    
    if (HydroMethod == PPM_DirectEuler) {
      FORTRAN_NAME(ppm_de)(
			density, totalenergy, velocity1, velocity2, velocity3,
                          gasenergy,
			&GravityOn, AccelerationField[0],
                           AccelerationField[1],
                           AccelerationField[2],
			&Gamma, &dtFixed, &CycleNumber,
                          CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			&GridRank, &GridDimension[0], &GridDimension[1],
                           &GridDimension[2], GridStartIndex, GridEndIndex,
			GridVelocity, &PPMFlatteningParameter,
                           &PressureFree,
			&PPMDiffusionParameter, &PPMSteepeningParameter,
                           &DualEnergyFormalism, &DualEnergyFormalismEta1,
			   &DualEnergyFormalismEta2,
			&NumberOfSubgrids, leftface, rightface,
			istart, iend, jstart, jend,
			standard, dindex, Eindex, uindex, vindex, windex,
			  geindex, temp,
                        &NumberOfColours, colourpt, coloff, colindex);
    }

    if (HydroMethod == PPM_LagrangeRemap) {
      FORTRAN_NAME(ppm_lr)(
			density, totalenergy, velocity1, velocity2, velocity3,
                          gasenergy,
			&GravityOn, AccelerationField[0],
                           AccelerationField[1],
                           AccelerationField[2],
			&Gamma, &dtFixed, &CycleNumber,
                          CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			&GridRank, &GridDimension[0], &GridDimension[1],
                           &GridDimension[2], GridStartIndex, GridEndIndex,
			GridVelocity, &PPMFlatteningParameter,
                           &PressureFree,
			&PPMDiffusionParameter, &PPMSteepeningParameter,
                           &DualEnergyFormalism, &DualEnergyFormalismEta1,
                           &DualEnergyFormalismEta2,
			&NumberOfSubgrids, leftface, rightface,
			istart, iend, jstart, jend,
			standard, dindex, Eindex, uindex, vindex, windex,
			  geindex, temp,
                        &NumberOfColours, colourpt, coloff, colindex);
    }

    if (HydroMethod == Zeus_Hydro) {
      FORTRAN_NAME(zeus_main)(
			density, totalenergy, velocity1, velocity2, velocity3,
                          zeus_temp, &ZEUSLinearArtificialViscosity,
			  &ZEUSQuadraticArtificialViscosity,
			&GravityOn, AccelerationField[0],
                           AccelerationField[1],
                           AccelerationField[2],
			GammaField, &dtFixed, &CycleNumber,
                          CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			&GridRank, &GridDimension[0], &GridDimension[1],
                           &GridDimension[2], GridStartIndex, GridEndIndex,
			GridVelocity, &PPMFlatteningParameter,
                           &PressureFree,
			&PPMDiffusionParameter, &PPMSteepeningParameter,
                           &DualEnergyFormalism, &UseGammaField,
                           &DualEnergyFormalismEta1, &DualEnergyFormalismEta2,
			&NumberOfSubgrids, leftface, rightface,
			istart, iend, jstart, jend,
			standard, dindex, Eindex, uindex, vindex, windex,
			  geindex, temp,
                        &NumberOfColours, colourpt, coloff, colindex,
			&LowestLevel, &MinimumSupportEnergyCoefficient);
    }

/* Put the update of magnetic fields here */



if (Passive_Used==1)  {

    // Create Electric fields from velocity and magnetic fields
    int UseDT = 1;
    float Olddt = Time - OldTime;
    
if(0){
 FORTRAN_NAME(passivemhd_create_e)(density,BaryonField[Vel1Num],BaryonField[Vel2Num], BaryonField[Vel3Num],
          OldBaryonField[Vel1Num],OldBaryonField[Vel2Num], OldBaryonField[Vel3Num],
          BaryonField[GENum],
          MagneticField[0], MagneticField[1], MagneticField[2],
          OldMagneticField[0], OldMagneticField[1], OldMagneticField[2],
          ElectricField[0], ElectricField[1], ElectricField[2],
          CellWidth[0], CellWidth[1], CellWidth[2],
          GridDimension, GridDimension +1, GridDimension +2,
          &dtFixed, &dtFixed, &UseDT,
          GridStartIndex, GridEndIndex,
          GridStartIndex+1, GridEndIndex+1,
          GridStartIndex+2, GridEndIndex+2,
          &a,&Biermann,
          &speedoflight,&hydrogenmass, &electroncharge, &chi
      );

    // Update Magnetic Fields
    float dtUsed = 1.0;

    FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
                                ElectricField[0], ElectricField[1], ElectricField[2],
                                CellWidth[0], CellWidth[1], CellWidth[2],
                                GridDimension, GridDimension +1, GridDimension +2,
                                GridStartIndex, GridEndIndex,
                                GridStartIndex+1, GridEndIndex+1,
                                GridStartIndex+2, GridEndIndex+2,
                                &dtUsed, &MHD_CenteringMethod);

    for(i=0;i<MagneticSize[0];i++)
       MagneticField[0][i]=0.5*(OldMagneticField[0][i]+MagneticField[0][i]);
    for(i=0;i<MagneticSize[1];i++)
       MagneticField[1][i]=0.5*(OldMagneticField[1][i]+MagneticField[1][i]);
    for(i=0;i<MagneticSize[2];i++)
       MagneticField[2][i]=0.5*(OldMagneticField[2][i]+MagneticField[2][i]);

}
    
if(1){
      FORTRAN_NAME(passivemhd_create_e)(BaryonField[DensNum],OldBaryonField[Vel1Num],OldBaryonField[Vel2Num],OldBaryonField[Vel3Num],
          OldBaryonField[Vel1Num],OldBaryonField[Vel2Num], OldBaryonField[Vel3Num],
          BaryonField[GENum],
          MagneticField[0], MagneticField[1], MagneticField[2],
          OldMagneticField[0], OldMagneticField[1], OldMagneticField[2],
          ElectricField[0], ElectricField[1], ElectricField[2],
          CellWidth[0], CellWidth[1], CellWidth[2],
          GridDimension, GridDimension +1, GridDimension +2,
          &dtFixed, &Olddt, &UseDT,
          GridStartIndex, GridEndIndex,
          GridStartIndex+1, GridEndIndex+1,
          GridStartIndex+2, GridEndIndex+2,
          &a,&Biermann,
          &speedoflight,&hydrogenmass, &electroncharge, &chi
      );

    // Update Magnetic Fields
    float dtUsed = 1.0;

    FORTRAN_NAME(curl_of_e)(OldMagneticField[0], OldMagneticField[1], OldMagneticField[2],
                                ElectricField[0], ElectricField[1], ElectricField[2],
                                CellWidth[0], CellWidth[1], CellWidth[2],
                                GridDimension, GridDimension +1, GridDimension +2,
                                GridStartIndex, GridEndIndex,
                                GridStartIndex+1, GridEndIndex+1,
                                GridStartIndex+2, GridEndIndex+2, 
                                &dtUsed, &MHD_CenteringMethod);

   for(i=0;i<MagneticSize[0];i++)
       MagneticField[0][i]=OldMagneticField[0][i];
    for(i=0;i<MagneticSize[1];i++)
       MagneticField[1][i]=OldMagneticField[1][i];
    for(i=0;i<MagneticSize[2];i++)
       MagneticField[2][i]=OldMagneticField[2][i];

/*
FORTRAN_NAME(passivemhd_computer_b)(OldBaryonField[DensNum],OldBaryonField[Vel1Num],OldBaryonField[Vel2Num],OldBaryonField[Vel3Num],
          OldBaryonField[Vel1Num],OldBaryonField[Vel2Num], OldBaryonField[Vel3Num],
          OldBaryonField[GENum],
          CenteredB[0], CenteredB[1], CenteredB[2],            
          CellWidth[0], CellWidth[1], CellWidth[2],            
          GridDimension, GridDimension +1, GridDimension +2,
          &dtFixed, &UseDT,
          GridStartIndex, GridEndIndex,
          GridStartIndex+1, GridEndIndex+1,
          GridStartIndex+2, GridEndIndex+2,
          &a,&Biermann
      );
*/
#ifdef OLD_CENTER
   
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField, first call \n");
    return FAIL;
  }
#endif //OLD_CENTER
   
}

/* fill Subgrid->ElectricField for MHD boundary correction */
  //define fluxe(dim,coord,face,end,sub) FluxExtents[dim+3*(coord+3*(face+2*(end+2*sub)))]
  
  
  for(subgrid=0;subgrid<NumberOfSubgrids; subgrid++)
    for(dim=0;dim<3;dim++)
      for( field=0;field<3;field++){
	
	SubgridFluxes[subgrid]->RightElectric[field][dim]=NULL;
	SubgridFluxes[subgrid]->LeftElectric[field][dim]=NULL;

     size=1; 
	if( field != dim){
	  for( coord=0;coord<3;coord++){
	    Dim[coord]=FluxDims[dim][coord][subgrid] ;
	    if(coord != field && coord != dim )
	      Dim[coord]++;
	    size *= Dim[coord];
	  }
	  //fprintf(stderr, "sgfeea: level %d sub %d dim %d field %d size %d\n",
	  //level, subgrid, dim, field, size);
	  SubgridFluxes[subgrid]->LeftElectric[field][dim]=new float[size];
	  SubgridFluxes[subgrid]->RightElectric[field][dim]=new float[size];
	  
	  for(i=0;i<size;i++){
	    SubgridFluxes[subgrid]->LeftElectric[field][dim][i]=0.0;
	    SubgridFluxes[subgrid]->RightElectric[field][dim][i] = 0.0;
	  }

       //for your reference:
	//fluxe(dim, coord, left or right, start or end, subgrid)
	  /*
	fprintf(stderr, "sgfee2: level %d fluxes(%d, 2, 0, 0, %d ) %d %d\n",
		level, dim, subgrid, fluxe(dim,2,0,0,subgrid), fluxe(dim,2,0,1,subgrid));
	fprintf(stderr, "sgfee2: level %d fluxes(%d, 1, 0, 0, %d ) %d %d\n",
		level, dim, subgrid, fluxe(dim,1,0,0,subgrid), fluxe(dim,1,0,1,subgrid));
	fprintf(stderr, "sgfee2: level %d fluxes(%d, 0, 0, 0, %d ) %d %d\n",
		level, dim, subgrid, fluxe(dim,0,0,0,subgrid), fluxe(dim,0,0,1,subgrid));
	  */
	for(k=0;k<Dim[2];k++)
	  for(j=0;j<Dim[1];j++)
	    for(i=0;i<Dim[0];i++){
	      Findex=i+Dim[0]*(j+Dim[1]*k);
	      Cindex=(i+fluxe(dim,0,0,0,subgrid)-1)
		+ElectricDims[field][0]*(j+fluxe(dim,1,0,0,subgrid)-1
		 +ElectricDims[field][1]*(k+fluxe(dim,2,0,0,subgrid)-1));
	      SubgridFluxes[subgrid]->LeftElectric[field][dim][Findex]= 
		ElectricField[field][Cindex];//*dtFixed; taken care of in the definition of ElectricField
	      Cindex=(i+fluxe(dim,0,1,0,subgrid)-1 )
		+ElectricDims[field][0]*(j+fluxe(dim,1,1,0,subgrid)-1
		 +ElectricDims[field][1]*(k+fluxe(dim,2,1,0,subgrid)-1));
	      SubgridFluxes[subgrid]->RightElectric[field][dim][Findex]=
		ElectricField[field][Cindex];// *dtFixed; taken care of in the definition of E.
            }
        }//field !=dim
       }          



}//end Passive_used



/* fill Subgrid->ElectricField for MHD boundary correction */
  //define fluxe(dim,coord,face,end,sub) FluxExtents[dim+3*(coord+3*(face+2*(end+2*sub)))]

   JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDEFlux0");
  for(subgrid=0;subgrid<NumberOfSubgrids; subgrid++)
    for(dim=0;dim<3;dim++)
      for( field=0;field<3;field++){
                                       
        SubgridFluxes[subgrid]->RightElectric[field][dim]=NULL;
        SubgridFluxes[subgrid]->LeftElectric[field][dim]=NULL;


        size=1;
        if( field != dim){
  
          for( coord=0;coord<3;coord++){
            Dim[coord]=FluxDims[dim][coord][subgrid] ;
            if(coord != field && coord != dim )
              Dim[coord]++;
            size *= Dim[coord];
          }                            
          //fprintf(stderr, "sgfeea: level %d sub %d dim %d field %d size %d\n",
          //level, subgrid, dim, field, size);
          SubgridFluxes[subgrid]->LeftElectric[field][dim]=new float[size];
          SubgridFluxes[subgrid]->RightElectric[field][dim]=new float[size];
          
          for(i=0;i<size;i++){
            SubgridFluxes[subgrid]->LeftElectric[field][dim][i]=0.0;
            SubgridFluxes[subgrid]->RightElectric[field][dim][i] = 0.0;
          }
          //Global memory leak here, too.


      for(k=0;k<Dim[2];k++)          
          for(j=0;j<Dim[1];j++)
            for(i=0;i<Dim[0];i++){
              Findex=i+Dim[0]*(j+Dim[1]*k);
          
              Cindex=(i+fluxe(dim,0,0,0,subgrid)-1)
                +ElectricDims[field][0]*(j+fluxe(dim,1,0,0,subgrid)-1
                 +ElectricDims[field][1]*(k+fluxe(dim,2,0,0,subgrid)-1));
            
              SubgridFluxes[subgrid]->LeftElectric[field][dim][Findex]=
                ElectricField[field][Cindex];//*dtFixed; taken care of in the definition of ElectricField
          
              Cindex=(i+fluxe(dim,0,1,0,subgrid)-1 )
                +ElectricDims[field][0]*(j+fluxe(dim,1,1,0,subgrid)-1
                 +ElectricDims[field][1]*(k+fluxe(dim,2,1,0,subgrid)-1));
          
              SubgridFluxes[subgrid]->RightElectric[field][dim][Findex]=
                ElectricField[field][Cindex];// *dtFixed; taken care of in the definition of E.

          }
        }//field==dim
          
      }//E field

      JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeMagFluxDelete");
  for(field=0;field<3;field++)
    for(axis=0;axis<2;axis++){
      delete [] MagneticFlux[field][axis];
    }
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDAfterLocalFluxDelete");   





    /* deallocate temporary space for solver */

    delete [] temp;
    if (GridRank < 3) delete [] velocity3;
    if (GridRank < 2) delete [] velocity2;
    if (HydroMethod == Zeus_Hydro)
      delete [] zeus_temp;

    delete [] leftface;
    delete [] GammaField;

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] CellWidthTemp[dim];

      /* deallocate temporary space for magnetic fluxes */
    for(field=0;field<3;field++)
    for(axis=0;axis<2;axis++){
      delete [] MagneticFlux[field][axis];
    }    

  }  // end: if (NumberOfBaryonFields > 0)



  this->DebugCheck("SolveHydroEquations (after)");
  if( level == 1 && grid == 1 ) {
    fprintf(stderr,"gonna stop.\n");
  }
  return SUCCESS;

}
#endif /* PASSIVE*/
