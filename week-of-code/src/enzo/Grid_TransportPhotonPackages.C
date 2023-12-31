#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (TRANSPORT PHOTON PACKAGES)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This is the heart of the radiative transfer algorithm.
/    On each Grid we initialize photo and heating rates and then call
/    WalkPhotonPackage so all photon packages are transported along their
/    own directions and the photo-ionization and heating rates on 
/    on the grid are updated on the fly. 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

#ifdef CONFIG_BFLOAT_4
#define ROUNDOFF 1e-6
#endif
#ifdef CONFIG_BFLOAT_8
#define ROUNDOFF 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define ROUNDOFF 1e-16
#endif

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
PhotonPackageEntry *PopPhoton(PhotonPackageEntry * &Node);
PhotonPackageEntry *DeletePhotonPackage(PhotonPackageEntry *PP);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::TransportPhotonPackages(int level, ListOfPhotonsToMove **PhotonsToMove, 
				  int GridNum, grid **Grids0, int nGrids0, 
				  grid *ParentGrid, grid *CurrentGrid)
{

  int i,j,k, dim, index, count;
  grid *MoveToGrid;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || MultiSpecies < 1) 
    return SUCCESS;

  if (RadiativeTransfer < 1) 
    return SUCCESS;

  if (RadiativeTransfer > 0 && GridRank < 3) {
    fprintf(stderr, "Grid_TransportPhotonPackage: failed\n");
    fprintf(stderr, "Grid_TransportPhotonPackage: "
	    "Transfer in less than 3D is not implemented.\n");
    ENZO_FAIL("");
  }

  if (PhotonPackages->NextPackage == NULL)
    return SUCCESS;

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stdout, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /* Find radiative transfer fields. */

  int kphHINum, gammaHINum, kphHeINum, gammaHeINum, kphHeIINum, gammaHeIINum,
    kdissH2INum;
  if (IdentifyRadiativeTransferFields(kphHINum, gammaHINum, kphHeINum, 
				      gammaHeINum, kphHeIINum, gammaHeIINum, 
				      kdissH2INum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifyRadiativeTransferFields.\n");
    ENZO_FAIL("");
  }

  int RPresNum1, RPresNum2, RPresNum3;
  if (RadiationPressure) {
    if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3)
	== FAIL) {
      fprintf(stdout, "Error in IdentifyRadiationPressureFields.\n");
      ENZO_FAIL("");
    }
  }

  /* Get units. */

  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    fprintf(stdout, "Error in GetUnits.\n");
    ENZO_FAIL("");
  }

  if (DEBUG) fprintf(stdout,"TransportPhotonPackage: initialize fields.\n");
  if (DEBUG) fprintf(stdout,"TransportPhotonPackage: %"ISYM" %"ISYM" .\n",
		     GridStartIndex[0], GridEndIndex[0]);

  PhotonPackageEntry *PP, *FPP, *SavedPP, *PausedPP;
  PP = PhotonPackages;

  if (DEBUG) {
    count = 0;
    while ((PP->NextPackage) != NULL) { 
      count++;
      PP=PP->NextPackage;
    }
    fprintf(stdout, "TransportPhotonPackage: done initializing.\n");
    fprintf(stdout, "counted %"ISYM" packages\n", count);
  }

  /* If requested, make vertex centered field (only when it doesn't
     exist ... see inside routine). */

  if (RadiativeTransferInterpolateField)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (FieldsToInterpolate[i] == TRUE)
	if (this->ComputeVertexCenteredField(i) == FAIL) {
	  fprintf(stderr, "Error in grid->ComputeVertexCenteredField "
		  "(field %"ISYM").\n", i);
	  ENZO_FAIL("");
	}

  count = 0;
  PP = PhotonPackages->NextPackage;
  FPP = this->FinishedPhotonPackages;
  PausedPP = this->PausedPhotonPackages;
  
  int dcount = 0;
  int tcount = 0;
  int pcount = 0;
  int trcount = 0;
  int AdvancePhotonPointer;
  int DeleteMe, DeltaLevel, PauseMe;

  FLOAT EndTime = PhotonTime+dtPhoton-ROUNDOFF;

  while (PP != NULL) {

    DeleteMe = FALSE;
    PauseMe = FALSE;
    MoveToGrid = NULL;
    AdvancePhotonPointer = TRUE;

    if ((PP->CurrentTime) < EndTime) {
      WalkPhotonPackage(&PP,
			&MoveToGrid, ParentGrid, CurrentGrid, Grids0, nGrids0,
			DensNum, HINum, HeINum, HeIINum, H2INum,
			kphHINum, gammaHINum, kphHeINum, gammaHeINum,
			kphHeIINum, gammaHeIINum, kdissH2INum, RPresNum1,
			RPresNum2, RPresNum3, DeleteMe, PauseMe, DeltaLevel, 
			DensityUnits, TemperatureUnits, VelocityUnits, 
			LengthUnits, TimeUnits);
      tcount++;
    } else {

      /* If all work is finished, store in FinishedPhotonPackages and
	 don't check for work until next timestep */

      SavedPP = PopPhoton(PP);
      PP = PP->NextPackage;
      InsertPhotonAfter(FPP, SavedPP);
      AdvancePhotonPointer = FALSE;

    }

    if (DEBUG > 1) 
      fprintf(stdout, "photon #%"ISYM" %x %x %x\n",
	      tcount,  PP,  PhotonPackages, 
	       MoveToGrid); 

    if (PauseMe == TRUE) {
      if (DEBUG > 1) fprintf(stdout, "paused photon %x\n", PP);
      SavedPP = PopPhoton(PP);
//      printf("paused photon %x lvl %"ISYM" ipix %"ISYM" CSRC %x leafID %"ISYM"\n",
//	     SavedPP, SavedPP->level, SavedPP->ipix, SavedPP->CurrentSource,
//	     SavedPP->CurrentSource->LeafID);
      PP = PP->NextPackage;
      InsertPhotonAfter(PausedPP, SavedPP);
      AdvancePhotonPointer = FALSE;
      MoveToGrid = NULL;
      pcount++;
    }

    if (DeleteMe == TRUE) {   
      if (DEBUG > 1) fprintf(stdout, "delete photon %x\n", PP);
      dcount++;
      PP = DeletePhotonPackage(PP);
      MoveToGrid = NULL;
    } 

    if (MoveToGrid != NULL) {
      if (DEBUG) {
	fprintf(stdout, "moving photon from %x to %x\n", 
		 CurrentGrid,  MoveToGrid);
	fprintf(stdout, "moving photon %x %x %x %x\n", 
		 PP,  PP->PreviousPackage, 
		 PP->NextPackage,  PhotonPackages);
      }
      ListOfPhotonsToMove *NewEntry = new ListOfPhotonsToMove;
      NewEntry->NextPackageToMove = (*PhotonsToMove)->NextPackageToMove;
      (*PhotonsToMove)->NextPackageToMove = NewEntry;
      NewEntry->PhotonPackage = PP;
      NewEntry->FromGrid = CurrentGrid;
      NewEntry->ToGrid   = MoveToGrid;
      NewEntry->ToGridNum= MoveToGrid->GetGridID();
      NewEntry->ToLevel  = level + DeltaLevel;
      NewEntry->ToProcessor = MoveToGrid->ReturnProcessorNumber();

      if (NewEntry->ToProcessor >= NumberOfProcessors)
	printf("TransportPH(P%"ISYM" :: G%"ISYM"): WARNING BAD TO_PROC -- P%"ISYM"->P%"ISYM".\n",
	       MyProcessorNumber, GridNum, ProcessorNumber, 
	       NewEntry->ToProcessor);

      if (PP->PreviousPackage != NULL) 
	PP->PreviousPackage->NextPackage = PP->NextPackage;
      if (PP->NextPackage != NULL) 
	PP->NextPackage->PreviousPackage = PP->PreviousPackage;
      trcount++;
    } // ENDIF MoveToGrid

    if (AdvancePhotonPointer == TRUE)
      PP = PP->NextPackage;

    // Merge "paused" photons only when all photons have been transported
    if (PP == NULL && PausedPP->NextPackage != NULL) {
      if (this->MergePausedPhotonPackages() == FAIL) {
	fprintf(stderr, "Error in grid::MergePausedPhotonPackages.\n");
	ENZO_FAIL("");
      }
      // Reset temp pointers
      PP = PhotonPackages->NextPackage;
      //FPP = this->FinishedPhotonPackages;
      PausedPP = this->PausedPhotonPackages;
      this->PausedPhotonPackages->NextPackage = NULL;
      this->PausedPhotonPackages->PreviousPackage = NULL;
    }

  } // ENDWHILE photons

  if (DEBUG) {
    fprintf(stdout, "grid::TransportPhotonPackage: "
	    "transported %"ISYM" deleted %"ISYM" paused %"ISYM"\n",
	    tcount, dcount, pcount);
    printf("L%d/G%d (%x): tr %d, del %d, move %d (NumberOfPhotons = %d/%d/%d)\n",
	   level, this->ID, this, tcount, dcount, trcount, NumberOfPhotonPackages,
	   this->ReturnRealPhotonCount(), NumberOfPhotonPackages-dcount);
  }
  NumberOfPhotonPackages -= dcount;

  /* For safety, clean up paused photon list */

#ifdef UNUSED
  if (PausedPhotonPackages->NextPackage != NULL) {
    PausedPP = PausedPhotonPackages->NextPackage;
    while (PausedPP != NULL) {
      PausedPP = DeletePhotonPackage(PausedPP);
      PausedPP = PausedPP->NextPackage;
    }
    PausedPhotonPackages->NextPackage = NULL;
    PausedPhotonPackages->PreviousPackage = NULL;
  }
#endif

#ifdef UNUSED
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    if (HasRadiation == TRUE) break;
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      if (HasRadiation == TRUE) break;
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	if (BaryonField[kphHINum][index] > 0) {
	  HasRadiation = TRUE;
	  break;
	}
      } // ENDFOR i
    }  // ENDFOR j
  } // ENDFOR k
#endif /* UNUSED */

  // Debug xyz-axis for a unigrid 64^3 with a source in the corner.
#define NO_DEBUG_AXES
#ifdef DEBUG_AXES
  printf("PHDebug(x): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][14914], BaryonField[kphHINum][14915], 
	 BaryonField[kphHINum][14916], 
	 BaryonField[kphHeIINum][14914], BaryonField[kphHeIINum][14915], 
	 BaryonField[kphHeIINum][14916], 
	 BaryonField[HINum][14914], BaryonField[HINum][14915], 
	 BaryonField[HINum][14916]);
  printf("PHDebug(y): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][14983], BaryonField[kphHINum][15053], 
	 BaryonField[kphHINum][15123], 
	 BaryonField[kphHeIINum][14983], BaryonField[kphHeIINum][15053], 
	 BaryonField[kphHeIINum][15123], 
	 BaryonField[HINum][14983], BaryonField[HINum][15053], 
	 BaryonField[HINum][15123]);
  printf("PHDebug(z): kph= %"GSYM" %"GSYM" %"GSYM", Nph = %"GSYM" %"GSYM" %"GSYM"\n, HI = %"GSYM" %"GSYM" %"GSYM"\n",
	 BaryonField[kphHINum][19813], BaryonField[kphHINum][24713], 
	 BaryonField[kphHINum][29613], 
	 BaryonField[kphHeIINum][19813], BaryonField[kphHeIINum][24713], 
	 BaryonField[kphHeIINum][29613], 
	 BaryonField[HINum][19813], BaryonField[HINum][24713], 
	 BaryonField[HINum][29613]);
#endif /* DEBUG_AXES */
	 
  return SUCCESS;
}
