/***********************************************************************
/
/  INITIALIZE THE METAL-COOLING RATES
/
/  written by: John Wise
/  date:       May, 2008
/  modified1:
/
/  PURPOSE:
/    For runs with chemical feedback, initialize the CoolData.metals 
/    rate table.  The table must be tabulated against temperature (rows) 
/    and electron fraction (columns).  Number of temperature bins must be 
/    equal to CoolData.NumberOfTemperatureBins.
/
/  RETURNS: ENZO_SUCCESS or FAIL
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

int ReadMetalCoolingRates(float TemperatureUnits, float LengthUnits, 
			  float aUnits, float DensityUnits, float TimeUnits, 
			  float aye)
{

  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int i, NumberOfTemperatureBins;
  float TemperatureRange[2];

  if ((fptr = fopen(MetalCoolingTable, "r")) == NULL) {
    fprintf(stderr, "Error opening metal cooling table %s\n", MetalCoolingTable);
    ENZO_FAIL("");
  }

  // The second and third lines have the number of bins and temp/x_e
  // ranges, respectively.

  fgets(line, MAX_LINE_LENGTH, fptr);
  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"ISYM" %"ISYM, &NumberOfTemperatureBins, 
	      &CoolData.NumberOfElectronFracBins)) != 2) {
    fprintf(stderr, "Error reading number of bins (line 2)\n");
    ENZO_FAIL("");
  }

  if (NumberOfTemperatureBins != CoolData.NumberOfTemperatureBins) {
    fprintf(stderr, 
	    "Number of temperature bins (=%"ISYM") in metal cooling table MUST equal\n"
	    "NumberOfTemperatureBins in other rate tables (=%"ISYM")\n",
	    NumberOfTemperatureBins, CoolData.NumberOfTemperatureBins);
    ENZO_FAIL("");
  }

  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"FSYM" %"FSYM" %"FSYM" %"FSYM, 
  //if ((sscanf(line, "# %e %e %e %e", 
	      TemperatureRange, TemperatureRange+1,
	      &CoolData.ElectronFracStart, &CoolData.ElectronFracEnd)) != 4) {
    fprintf(stderr, "Error reading number of ranges (line 3)\n");
    ENZO_FAIL("");
  }

  if (TemperatureRange[0] != CoolData.TemperatureStart ||
      TemperatureRange[1] != CoolData.TemperatureEnd) {
    fprintf(stderr, 
	    "Temperature range [%"GSYM", %"GSYM"] in metal cooling table MUST equal the\n"
	    "temperature range [%"GSYM", %"GSYM"] in the other rate tables.\n",
	    TemperatureRange[0], TemperatureRange[1], CoolData.TemperatureStart,
	    CoolData.TemperatureEnd);
    ENZO_FAIL("");
  }

  int prev_pos;
  int ixe, itemp, index, nbins;

  nbins = CoolData.NumberOfTemperatureBins * CoolData.NumberOfElectronFracBins;
  CoolData.metals = new float[nbins];

  // Pass comments
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (line[0] != '#') {
      fseek(fptr, prev_pos, SEEK_SET);
      break;
    } // ENDIF not comment
    prev_pos = ftell(fptr);
  } // ENDWHILE file

  // Read the rate table
  index = 0;
  for (itemp = 0; itemp < CoolData.NumberOfTemperatureBins; itemp++)
    for (ixe = 0; ixe < CoolData.NumberOfElectronFracBins; ixe++) {
      index = ixe*CoolData.NumberOfTemperatureBins + itemp;
      if ((fscanf(fptr, "%"FSYM, &CoolData.metals[index])) == EOF) {
	fprintf(stderr, "EOF reached at itemp = %"ISYM", ixe = %"ISYM"\n", 
		itemp, ixe);
	ENZO_FAIL("");
      }
    }

  fclose(fptr);

  /* Convert to code units (a la calc_rates.src) */

  const double mh = 1.673e-24;
  double tbase1, xbase1, dbase1, coolunit;

  tbase1 = TimeUnits;
  xbase1 = LengthUnits / (aye * aUnits);
  dbase1 = DensityUnits * pow(aye*aUnits, 3.0);
  coolunit = (pow(aUnits, 5) * xbase1*xbase1 * mh*mh) / 
    (pow(tbase1, 3) * dbase1);

  for (i = 0; i < nbins; i++)
    CoolData.metals[i] /= coolunit;

  return SUCCESS;

}
