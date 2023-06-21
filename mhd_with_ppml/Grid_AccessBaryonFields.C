/***********************************************************************
/
/  GRID CLASS (ACCESS THE BARYON FIELDS)
/
/  written by: Daniel R. Reynolds
/  date:       October, 2006
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
// Density field 
float* grid::AccessDensity() {
  int DensNum = 0;
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessDensity cannot find density field\n");
    return NULL;
  }
  return BaryonField[DensNum];
}

// Total Energy field 
float* grid::AccessTotalEnergy() {
  int EnergyNum = 0;
  if ((EnergyNum = FindField(TotalEnergy,FieldType,NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessTotalEnergy cannot find total energy field\n");
    return NULL;
  }
  return BaryonField[EnergyNum];
}

// Gas Energy field 
float* grid::AccessGasEnergy() {
  int GENum = 0;
  if (((GENum = FindField(InternalEnergy,FieldType,NumberOfBaryonFields))<0) 
      || (DualEnergyFormalism == FALSE)) {
    fprintf(stderr,"Grid:AccessGasEnergy cannot find internal energy field\n");
    return NULL;
  }
  return BaryonField[GENum];
}

// Velocity1 field 
float* grid::AccessVelocity1() {
  int Vel1Num = 0;
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessVelocity1 cannot find velocity field\n");
    return NULL;
  }
  return BaryonField[Vel1Num];
}

// Velocity2 field 
float* grid::AccessVelocity2() {
  int Vel2Num = 0;
  if ((Vel2Num = FindField(Velocity2, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessVelocity2 cannot find velocity field\n");
    return NULL;
  }
  return BaryonField[Vel2Num];
}

// Velocity3 field 
float* grid::AccessVelocity3() {
  int Vel3Num = 0;
  if ((Vel3Num = FindField(Velocity3, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessVelocity3 cannot find velocity field\n");
    return NULL;
  }
  return BaryonField[Vel3Num];
}

// Electron Density field 
float* grid::AccessElectronDensity() {
  int ENum = 0;
  if ((ENum = FindField(ElectronDensity,FieldType,NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessElectronDensity cannot find electron field\n");
    return NULL;
  }
  return BaryonField[ENum];
}

// Hydrogen-I Density field 
float* grid::AccessHIDensity() {
  int HINum = 0;
  if ((HINum = FindField(HIDensity, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHIDensity cannot find HI field\n");
    return NULL;
  }
  return BaryonField[HINum];
}

// Helium-I Density field 
float* grid::AccessHeIDensity() {
  int HeINum = 0;
  if ((HeINum = FindField(HeIDensity, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIDensity cannot find HeI field\n");
    return NULL;
  }
  return BaryonField[HeINum];
}

// Helium-II Density field 
float* grid::AccessHeIIDensity() {
  int HeIINum = 0;
  if ((HeIINum = FindField(HeIIDensity, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find HeII field\n");
    return NULL;
  }
  return BaryonField[HeIINum];
}

// Radiation Energy (Grey, or 0th bin)
float* grid::AccessRadiationFrequency0()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq0, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq0 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

// Radiation Energy (1st bin)
float* grid::AccessRadiationFrequency1()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq1, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq1 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

// Radiation Energy (2nd bin)
float* grid::AccessRadiationFrequency2()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq2, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq2 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

// Radiation Energy (3rd bin)
float* grid::AccessRadiationFrequency3()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq3, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq3 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

// Radiation Energy (4th bin)
float* grid::AccessRadiationFrequency4()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq4, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq4 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

// Radiation Energy (5th bin)
float* grid::AccessRadiationFrequency5()
 {
  int RadNum = 0;
  if ((RadNum = FindField(RadiationFreq5, FieldType, NumberOfBaryonFields))<0) {
    fprintf(stderr,"Grid:AccessHeIIDensity cannot find RadiationFreq5 field\n");
    return NULL;
  }
  return BaryonField[RadNum];
}

