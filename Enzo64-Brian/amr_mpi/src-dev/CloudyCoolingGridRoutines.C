/********************************************************************************
 *                       CloudyCoolingGridRoutines.C                            *
 *                      Britton Smith - November, 2005                          *
 *                                                                              *
 * Read in heating, cooling, and mean molecular weight values from files        *
 *                                                                              *
 * For CloudyCoolingGridRank = 0, values are only a function of temperature.    *
 * For CloudyCoolingGridRank = 1, values are function of density and temp.      *
 * For CloudyCoolingGridRank = 2, values are function of density, metallicity,  *
 *                                and temperature.                              *
 * Heating, cooling, and mmw are stored in 1D arrays.                           *
 * Cooling and heating are stored as log10 of the values.                       *
 * Linear interpolation is used to calculate non-gridpoint values.              *
 *                                                                              *
 * Routines in this file:                                                       *
 *          InitializeCloudyCooling - initialize data arrays and call functions *
 *                                    to read in data.                          *
 *          coolingGridReadGrid1D - read in data for interpolation over         *
 *                                  density and temperature.                    *
 *          coolingGridReadGrid2D - read in data for interpolation over         *
 *                                  density, metalliticy, and temperature.      *
 *          coolingGridReadCoolingMapFirstMap - read first cooling map in a set *
 *                                              and fill temperature, heating,  *
 *                                              cooling, and mmw arrays.        *
 *          coolingGridReadCoolingMap - read cooling map from a set and fill    *
 *                                      heating, cooling, and mmw array. Don't  *
 *                                      fill temperature array (already filled).*
 *          coolingGridGetParameterValuesFromLine - read in a list of values    *
 *                                                  from a line and place into  *
 *                                                  an array.                   *
 *          coolingGridInterpolate0D - interpolate over temperature and return  *
 *                                     value from thrown data field.  Returns   *
 *                                     either double or float depending on type *
 *                                     of data field.                           *
 *          coolingGridInterpolate1D - interpolate over density and temperature *
 *                                     and return value from thrown data field. *
 *                                     Returns double or float depending on     *
 *                                     type of data field.                      *
 *          coolingGridInterpolate2D - interpolate over density, metallicity,   *
 *                                     and temperature and return value from    *
 *                                     thrown data field.  Returns double or    *
 *                                     float depending on type of data field.   *
 ********************************************************************************/

#include <math.h>
#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
using namespace std;

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

/**************************** Functions Prototypes ******************************/

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

// Take file containing list of files produced by Cloudy and some other information.
// Read all files on list and fill heating, cooling, and mmw arrays.
// For 1D parameter space
int coolingGridReadGrid1D(char* coolingGridRunFile);

// For 2D parameter space
int coolingGridReadGrid2D(char* coolingGridRunFile);

// Read cooling map file and fill data arrays
int coolingGridReadCoolingMap(char *coolingMapFilename,float *heating,
			      float *cooling,float *mmw,int *offset);

// Read first cooling map file.
// Get same stuff as above, but also fill temperature array.
int coolingGridReadCoolingMapFirstMap(char *coolingMapFilename,float *temperature,
				      float *heating,float *cooling,float *mmw,int *offset);


// Take a string containing list of parameter values and put those values into float array.
int coolingGridGetParameterValuesFromLine(string lineString,char separator,int maxValues,
					  float *valuesPtr,int *numberOfValues);


// Initialize Cloudy Cooling
int InitializeCloudyCooling(FLOAT Time)
{

  CoolData.HydrogenFractionByMass = 0.76;

  FLOAT a = 1, dadt;

  /* If using cosmology, compute the expansion factor and get units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, MassUnits = 1, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }

  if (ComovingCoordinates) {

    if (CosmologyComputeExpansionFactor(Time, &a, &dadt) 
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
      return FAIL;
    }

    aUnits = 1.0/(1.0 + InitialRedshift);

  }

  /* Get conversion units (ripped from InitializeEquilibriumCoolData.C) */

  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(a*aUnits);
  double dbase1 = DensityUnits * POW(a*aUnits, 3);
  double mh = 1.67e-24;
  double CoolUnit = (POW(aUnits,5) * POW(xbase1,2) * POW(mh,2)) /
                    (POW(tbase1,3) * dbase1);

  int q;

  if (CloudyCoolingData.CloudyCoolingGridRank > MAX_CLOUDY_COOLING_DIMENSION) {
    fprintf(stderr,"Error in InitializeCloudyCooling: CloudyCoolingGridRank can not be greater than MAX_CLOUDY_COOLING_DIMENSION.\n");
    return FAIL;
  }

  // Interpolate only over temperature.

  else if (CloudyCoolingData.CloudyCoolingGridRank == 0) {

    // Set maxes for cooling data
    CloudyCoolingData.coolingGridMaxTemperatureValues = 122;

    // Initialize storage arrays
    CloudyCoolingData.coolingGridTemperature = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridHeating = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridCooling = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridMeanMolecularWeight = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];

    int coolingGridDataIndex = 0;
    if (coolingGridReadCoolingMapFirstMap(CloudyCoolingData.CloudyCoolingGridRunFile,
					  CloudyCoolingData.coolingGridTemperature,
					  CloudyCoolingData.coolingGridHeating,
					  CloudyCoolingData.coolingGridCooling,
					  CloudyCoolingData.coolingGridMeanMolecularWeight,
					  &coolingGridDataIndex) == FAIL) {
      fprintf(stderr, "Error in coolingGridReadGrid.\n");
      return FAIL;
    }
    CloudyCoolingData.coolingGridNumberTemperatureValues = coolingGridDataIndex;

    /* Convert cooling and heating values into code units. */

    for (q = 0;q < CloudyCoolingData.coolingGridNumberTemperatureValues;q++) {
      CloudyCoolingData.coolingGridHeating[q] -= log10(CoolUnit);
      CloudyCoolingData.coolingGridCooling[q] -= log10(CoolUnit);
    }

  }

  // Interpolate over density and temperature.

  else if (CloudyCoolingData.CloudyCoolingGridRank == 1) {

    // Set maxes for cooling data
    CloudyCoolingData.coolingGridMaxParameterValuesDimension1 = 41;
    CloudyCoolingData.coolingGridMaxTemperatureValues = 122;

    // Initialize storage arrays
    CloudyCoolingData.coolingGridParameterValuesDimension1 = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1];
    CloudyCoolingData.coolingGridTemperature = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridHeating = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		 CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridCooling = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		 CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridMeanMolecularWeight = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		CloudyCoolingData.coolingGridMaxTemperatureValues];

    // Read in cooling grid data
    if (coolingGridReadGrid1D(CloudyCoolingData.CloudyCoolingGridRunFile) == FAIL) {
      fprintf(stderr,"Error in InitializeCloudyCooling.\n");
      return FAIL;
    }

    /* Convert cooling and heating values into code units. */

    for (q = 0;q < (CloudyCoolingData.coolingGridNumberParameterValuesDim1*
		    CloudyCoolingData.coolingGridNumberTemperatureValues);q++) {
      CloudyCoolingData.coolingGridHeating[q] -= log10(CoolUnit);
      CloudyCoolingData.coolingGridCooling[q] -= log10(CoolUnit);
    }

  }

  // Interpolate over density, metallicity, and temperature.

  else if (CloudyCoolingData.CloudyCoolingGridRank == 2) {

    // Set maxes for cooling data
    CloudyCoolingData.coolingGridMaxParameterValuesDimension1 = 41;
    CloudyCoolingData.coolingGridMaxParameterValuesDimension2 = 41;
    CloudyCoolingData.coolingGridMaxTemperatureValues = 122;

    // Initialize storage arrays
    CloudyCoolingData.coolingGridParameterValuesDimension1 = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1];
    CloudyCoolingData.coolingGridParameterValuesDimension2 = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension2];
    CloudyCoolingData.coolingGridTemperature = 
      new float[CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridHeating = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		 CloudyCoolingData.coolingGridMaxParameterValuesDimension2*
		 CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridCooling = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		 CloudyCoolingData.coolingGridMaxParameterValuesDimension2*
		 CloudyCoolingData.coolingGridMaxTemperatureValues];
    CloudyCoolingData.coolingGridMeanMolecularWeight = 
      new float[CloudyCoolingData.coolingGridMaxParameterValuesDimension1*
		CloudyCoolingData.coolingGridMaxParameterValuesDimension2*
		CloudyCoolingData.coolingGridMaxTemperatureValues];

    // Read in cooling grid data
    if (coolingGridReadGrid2D(CloudyCoolingData.CloudyCoolingGridRunFile) == FAIL) {
      fprintf(stderr,"Error in InitializeCloudyCooling.\n");
      return FAIL;
    }

    /* Convert cooling and heating values into code units. */

    for (q = 0;q < (CloudyCoolingData.coolingGridNumberParameterValuesDim1*
		    CloudyCoolingData.coolingGridNumberParameterValuesDim2*
		    CloudyCoolingData.coolingGridNumberTemperatureValues);q++) {
      CloudyCoolingData.coolingGridHeating[q] -= log10(CoolUnit);
      CloudyCoolingData.coolingGridCooling[q] -= log10(CoolUnit);
    }

  }

  else {
    fprintf(stderr,"Error in InitializeCloudyCooling: CloudyCoolingGridRank must be 0, 1, or 2.\n");
    return FAIL;
  }

  return SUCCESS;
}

// Read cooling map file and fill data arrays for 1D + temperature interpolation
int coolingGridReadGrid1D(char* coolingGridRunFile)
{
  FILE *runptr;
  char line[MAX_LINE_LENGTH];
  char coolingGridFilePrefix[MAX_LINE_LENGTH];
  char coolingGridDataDir[MAX_LINE_LENGTH];
  char coolingMapFilename[MAX_LINE_LENGTH];
  string lineString,dataDirString;
  float parameter1;
  int fileIndex;
  int coolingGridDataIndex;

  // Initialize data array offset variable
  coolingGridDataIndex = 0;

  // Get data dir from path to run file
  dataDirString = coolingGridRunFile;
  dataDirString = dataDirString.substr(0,dataDirString.rfind("/")+1);
  memset(coolingGridDataDir,'\0',MAX_LINE_LENGTH);
  dataDirString.copy(coolingGridDataDir,dataDirString.length());

  // Open run file to read in file numbers and parameters
  if ((runptr = fopen(coolingGridRunFile, "r")) == NULL) {
    fprintf (stderr,"Error opening cooling grid run file %s.\n", coolingGridRunFile);  
    return FAIL;
  }

  while (fgets(line, MAX_LINE_LENGTH, runptr) != NULL) {
    lineString = line;

    // Get output file prefix.
    sscanf(line, "# outputFilePrefix = %s", coolingGridFilePrefix);

    // Get parameter values for cooling maps.
    if (lineString.compare("# Loop commands and values:\n") == 0) {
      // Get parameter values for first dimension.
      if (fgets(line, MAX_LINE_LENGTH, runptr) == NULL) {
	fprintf (stderr,"Error: cooling grid file is corrupt.\n");  
	return FAIL;
      }
      lineString = line;
      lineString = lineString.substr(lineString.find(":")+2,
				     lineString.length()-lineString.find(":")-3);
      if (coolingGridGetParameterValuesFromLine(lineString,' ',
						CloudyCoolingData.coolingGridMaxParameterValuesDimension1,
						CloudyCoolingData.coolingGridParameterValuesDimension1,
						&CloudyCoolingData.coolingGridNumberParameterValuesDim1) == FAIL) {
	fprintf (stderr,"Error in coolingGridReadGrid.\n");
	return FAIL;
      }

    }

    // Get cooling map file indices and parameter values
    if (sscanf(line,"%"ISYM"\t%"FSYM"\n",&fileIndex,&parameter1) == 2) {
      sprintf(coolingMapFilename,"%s%s_run%"ISYM".dat",coolingGridDataDir,coolingGridFilePrefix,fileIndex);

      if (coolingGridDataIndex) {
	if (coolingGridReadCoolingMap(coolingMapFilename,
				      CloudyCoolingData.coolingGridHeating,
				      CloudyCoolingData.coolingGridCooling,
				      CloudyCoolingData.coolingGridMeanMolecularWeight,
				      &coolingGridDataIndex) == FAIL) {
	  fprintf(stderr, "Error in coolingGridReadGrid.\n");
	  return FAIL;
	}
      }
      else {
	if (coolingGridReadCoolingMapFirstMap(coolingMapFilename,
					      CloudyCoolingData.coolingGridTemperature,
					      CloudyCoolingData.coolingGridHeating,
					      CloudyCoolingData.coolingGridCooling,
					      CloudyCoolingData.coolingGridMeanMolecularWeight,
					      &coolingGridDataIndex) == FAIL) {
	  fprintf(stderr, "Error in coolingGridReadGrid.\n");
	  return FAIL;
	}
	CloudyCoolingData.coolingGridNumberTemperatureValues = coolingGridDataIndex;
      }
    }
  }

  fclose(runptr);

  return SUCCESS;
}

// Read cooling map file and fill data arrays for 2D + temperature interpolation
int coolingGridReadGrid2D(char* coolingGridRunFile)
{
  FILE *runptr;
  char line[MAX_LINE_LENGTH];
  char coolingGridFilePrefix[MAX_LINE_LENGTH];
  char coolingGridDataDir[MAX_LINE_LENGTH];
  char coolingMapFilename[MAX_LINE_LENGTH];
  string lineString,dataDirString;
  float parameter1,parameter2;
  int fileIndex;
  int coolingGridDataIndex;

  // Initialize data array offset variable
  coolingGridDataIndex = 0;

  // Get data dir from path to run file
  dataDirString = coolingGridRunFile;
  dataDirString = dataDirString.substr(0,dataDirString.rfind("/")+1);
  memset(coolingGridDataDir,'\0',MAX_LINE_LENGTH);
  dataDirString.copy(coolingGridDataDir,dataDirString.length());

  // Open run file to read in file numbers and parameters
  if ((runptr = fopen(coolingGridRunFile, "r")) == NULL) {
    fprintf (stderr,"Error opening cooling grid run file %s.\n", coolingGridRunFile);  
    return FAIL;
  }

  while (fgets(line, MAX_LINE_LENGTH, runptr) != NULL) {
    lineString = line;

    // Get output file prefix.
    sscanf(line, "# outputFilePrefix = %s", coolingGridFilePrefix);

    // Get parameter values for cooling maps.
    if (lineString.compare("# Loop commands and values:\n") == 0) {
      // Get parameter values for first dimension.
      if (fgets(line, MAX_LINE_LENGTH, runptr) == NULL) {
	fprintf (stderr,"Error: cooling grid file is corrupt.\n");  
	return FAIL;
      }
      lineString = line;
      lineString = lineString.substr(lineString.find(":")+2,
				     lineString.length()-lineString.find(":")-3);
      if (coolingGridGetParameterValuesFromLine(lineString,' ',
						CloudyCoolingData.coolingGridMaxParameterValuesDimension1,
						CloudyCoolingData.coolingGridParameterValuesDimension1,
						&CloudyCoolingData.coolingGridNumberParameterValuesDim1) == FAIL) {
	fprintf (stderr,"Error in coolingGridReadGrid.\n");
	return FAIL;
      }
      
      // Get parameter values for second dimension.
      if (fgets(line, MAX_LINE_LENGTH, runptr) == NULL) {
	fprintf (stderr,"Error in coolingGridReadGrid: cooling grid file is corrupt.\n");  
	return FAIL;
      }
      lineString = line;
      lineString = lineString.substr(lineString.find(":")+2,
				     lineString.length()-lineString.find(":")-3);
      if (coolingGridGetParameterValuesFromLine(lineString,' ',
						CloudyCoolingData.coolingGridMaxParameterValuesDimension1,
						CloudyCoolingData.coolingGridParameterValuesDimension2,
						&CloudyCoolingData.coolingGridNumberParameterValuesDim2) == FAIL) {
	fprintf (stderr,"Error in coolingGridReadGrid.\n");
	return FAIL;
      }
    }

    // Get cooling map file indices and parameter values
    if (sscanf(line,"%"ISYM"\t%"FSYM"\t%"FSYM"\n",&fileIndex,&parameter1,&parameter2) == 3) {
      sprintf(coolingMapFilename,"%s%s_run%"ISYM".dat",coolingGridDataDir,coolingGridFilePrefix,fileIndex);

      if (coolingGridDataIndex) {
	if (coolingGridReadCoolingMap(coolingMapFilename,
				      CloudyCoolingData.coolingGridHeating,
				      CloudyCoolingData.coolingGridCooling,
				      CloudyCoolingData.coolingGridMeanMolecularWeight,
				      &coolingGridDataIndex) == FAIL) {
	  fprintf(stderr, "Error in coolingGridReadGrid.\n");
	  return FAIL;
	}
      }
      else {
	if (coolingGridReadCoolingMapFirstMap(coolingMapFilename,
					      CloudyCoolingData.coolingGridTemperature,
					      CloudyCoolingData.coolingGridHeating,
					      CloudyCoolingData.coolingGridCooling,
					      CloudyCoolingData.coolingGridMeanMolecularWeight,
					      &coolingGridDataIndex) == FAIL) {
	  fprintf(stderr, "Error in coolingGridReadGrid.\n");
	  return FAIL;
	}
	CloudyCoolingData.coolingGridNumberTemperatureValues = coolingGridDataIndex;
      }
    }
  }

  fclose(runptr);

  return SUCCESS;
}

// Read cooling map file and fill data arrays
int coolingGridReadCoolingMap(char *coolingMapFilename,float *heating,float *cooling,float *mmw,int *offset)
{
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  float floatValue1;
  float floatValue2,floatValue3;
  float floatValue4;
  int numberTemperatureValues = 0;

  // Open cooling map file to read in heating, cooling, mmw
  if ((fptr = fopen(coolingMapFilename, "r")) == NULL) {
    fprintf (stderr,"Error opening cooling map file %s in coolingGridReadCoolingMap.\n", coolingMapFilename);  
    return FAIL;
  }

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (sscanf(line,"%"ESYM"\t%"ESYM"\t%"ESYM"\t%"FSYM"\n",&floatValue1,&floatValue2,&floatValue3,&floatValue4) == 4) {
      heating[*offset] = (floatValue2 > 0) ? log10(floatValue2) : -99.0;
      cooling[*offset] = (floatValue3 > 0) ? log10(floatValue3) : -99.0;
      mmw[*offset] = floatValue4;
      (*offset)++;
      numberTemperatureValues++;
    }

  }

  fclose (fptr);

  if (numberTemperatureValues != CloudyCoolingData.coolingGridNumberTemperatureValues) {
    fprintf (stderr,"Error in coolingGridReadCoolingMap: cooling maps contain different numbers of points.\n"); 
    return FAIL;
  }

  return SUCCESS;
}

// Read first cooling map file.
// Get same stuff as above, but also fill temperature array.
int coolingGridReadCoolingMapFirstMap(char *coolingMapFilename,float *temperature,float *heating,
				      float *cooling,float *mmw,int *offset)
{
  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  float floatValue1;
  float floatValue2,floatValue3;
  float floatValue4;

  // Open cooling map file to read in heating, cooling, mmw
  if ((fptr = fopen(coolingMapFilename, "r")) == NULL) {
    fprintf (stderr,"Error opening cooling map file %s in coolingGridReadCoolingMapFirstMap.\n", coolingMapFilename);  
    return FAIL;
  }

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (*offset >= CloudyCoolingData.coolingGridMaxTemperatureValues) {
      fprintf (stderr,"Error in coolingGridReadCoolingMapFirstMap: number of cooling map temperature points exceeds maximum value.\n");
	return FAIL;
    }
    if (sscanf(line,"%"ESYM"\t%"ESYM"\t%"ESYM"\t%"FSYM"\n",&floatValue1,&floatValue2,&floatValue3,&floatValue4) == 4) {
      temperature[*offset] = floatValue1;
      heating[*offset] = (floatValue2 > 0) ? log10(floatValue2) : -99.0;
      cooling[*offset] = (floatValue3 > 0) ? log10(floatValue3) : -99.0;
      mmw[*offset] = floatValue4;
      (*offset)++;
    }
  }

  fclose (fptr);
  return SUCCESS;
}

// Take a string containing list of parameter values and put those values into float array.
int coolingGridGetParameterValuesFromLine(string lineString,char separator,int maxValues,
					  float *valuesPtr,int *numberOfValues)
{
  string value;
  int nextSpace;
  int index = 0;
  char charValue[MAX_LINE_LENGTH];
  while (lineString.length() > 0) {
    if (index >= maxValues) {
      fprintf (stderr,"Error in coolingGridGetParameterValuesFromLine: number of parameters exceeds maximum allowed.\n");
      return FAIL;
    }
    nextSpace = lineString.find(separator);
    if (nextSpace == string::npos) {
      value = lineString;
      lineString = "";
    }
    else {
      value = lineString.substr(0,nextSpace);
      lineString = lineString.substr(nextSpace+1);
    }
    memset(charValue,'\0',MAX_LINE_LENGTH);
    value.copy(charValue,value.length());
    valuesPtr[index] = atof(charValue);
    index++;
  }
  *numberOfValues = index;
  return SUCCESS;
}

// Take values for first and second parameters (probably density and metallicity) and temperature.
double coolingGridInterpolate2D_double(float parameter1,float parameter2,float temperature,double *dataField)
{
  int parameter1Index, parameter2Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q,w;
  double slope,temperatureValue[2],p2Value[2],p1Value;


  // find parameter1 index
  if (parameter1 <= CloudyCoolingData.coolingGridParameterValuesDimension1[0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1]) {
    parameter1Index = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter2 index
  if (parameter2 <= CloudyCoolingData.coolingGridParameterValuesDimension2[0]) {
    parameter2Index = 0;
  }
  else if (parameter2 >= CloudyCoolingData.coolingGridParameterValuesDimension2[CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 1]) {
    parameter2Index = CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 2;
  }
  else {
    parameter2Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 1;
    while (highpoint - parameter2Index > 1) {
      midpoint = (highpoint + parameter2Index) >> 1;
      if (parameter2 >= CloudyCoolingData.coolingGridParameterValuesDimension2[midpoint]) parameter2Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


  // interpolate over parameter 1
  for (q = 0;q < 2;q++) {
    // inerpolate over parameter 2
    for (w = 0;w < 2;w++) {
      // interpolate over temperature
      interpolateIndex = ((q+parameter1Index)*CloudyCoolingData.coolingGridNumberParameterValuesDim2 + (w+parameter2Index))
	* CloudyCoolingData.coolingGridNumberTemperatureValues + temperatureIndex;

      slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
	(CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);

      temperatureValue[w] = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
	dataField[interpolateIndex];
    }
    slope = (temperatureValue[1]-temperatureValue[0]) /
      (CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index+1]-CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index]);
    p2Value[q] = (parameter2-CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index]) * slope + temperatureValue[0];
  }
  slope = (p2Value[1]-p2Value[0]) /
    (CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index+1]-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]);
  p1Value = (parameter1-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]) * slope + p2Value[0];

  return p1Value;
}

// Float returning version of 2D interpolation routine.
// Do interpolation from specified data field.
float coolingGridInterpolate2D_float(float parameter1,float parameter2,float temperature,float *dataField)
{
  int parameter1Index, parameter2Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q,w;
  float slope,temperatureValue[2],p2Value[2],p1Value;


  // find parameter1 index
  if (parameter1 <= CloudyCoolingData.coolingGridParameterValuesDimension1[0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1]) {
    parameter1Index = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find parameter2 index
  if (parameter2 <= CloudyCoolingData.coolingGridParameterValuesDimension2[0]) {
    parameter2Index = 0;
  }
  else if (parameter2 >= CloudyCoolingData.coolingGridParameterValuesDimension2[CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 1]) {
    parameter2Index = CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 2;
  }
  else {
    parameter2Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim2 - 1;
    while (highpoint - parameter2Index > 1) {
      midpoint = (highpoint + parameter2Index) >> 1;
      if (parameter2 >= CloudyCoolingData.coolingGridParameterValuesDimension2[midpoint]) parameter2Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


  // interpolate over parameter 1
  for (q = 0;q < 2;q++) {
    // inerpolate over parameter 2
    for (w = 0;w < 2;w++) {
      // interpolate over temperature
      interpolateIndex = ((q+parameter1Index)*CloudyCoolingData.coolingGridNumberParameterValuesDim2 + (w+parameter2Index))
	* CloudyCoolingData.coolingGridNumberTemperatureValues + temperatureIndex;

      slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
	(CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);

      temperatureValue[w] = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
	dataField[interpolateIndex];
    }
    slope = (temperatureValue[1]-temperatureValue[0]) /
      (CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index+1]-CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index]);
    p2Value[q] = (parameter2-CloudyCoolingData.coolingGridParameterValuesDimension2[parameter2Index]) * slope + temperatureValue[0];
  }
  slope = (p2Value[1]-p2Value[0]) /
    (CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index+1]-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]);
  p1Value = (parameter1-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]) * slope + p2Value[0];

  return p1Value;
}


// Take values for first parameter (probably densityy) and temperature.
double coolingGridInterpolate1D_double(float parameter1,float temperature,double *dataField)
{
  int parameter1Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q;
  double slope,temperatureValue[2],p1Value;


  // find parameter1 index
  if (parameter1 <= CloudyCoolingData.coolingGridParameterValuesDimension1[0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1]) {
    parameter1Index = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


  // interpolate over parameter 1
  for (q = 0;q < 2;q++) {
    // interpolate over temperature
    interpolateIndex = (q+parameter1Index) * CloudyCoolingData.coolingGridNumberTemperatureValues + temperatureIndex;
    
    slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
      (CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);
    
    temperatureValue[q] = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
      dataField[interpolateIndex];
  }
  slope = (temperatureValue[1]-temperatureValue[0]) /
    (CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index+1]-
     CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]);

  p1Value = (parameter1-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]) * slope + temperatureValue[0];
  
  return p1Value;
}

// Float returning version of 1D interpolation routine.
// Take values for first parameter (probably densityy) and temperature.
float coolingGridInterpolate1D_float(float parameter1,float temperature,float *dataField)
{
  int parameter1Index, temperatureIndex;
  int midpoint,highpoint;
  int interpolateIndex,q;
  float slope,temperatureValue[2],p1Value;


  // find parameter1 index
  if (parameter1 <= CloudyCoolingData.coolingGridParameterValuesDimension1[0]) {
    parameter1Index = 0;
  }
  else if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1]) {
    parameter1Index = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 2;
  }
  else {
    parameter1Index = 0;
    highpoint = CloudyCoolingData.coolingGridNumberParameterValuesDim1 - 1;
    while (highpoint - parameter1Index > 1) {
      midpoint = (highpoint + parameter1Index) >> 1;
      if (parameter1 >= CloudyCoolingData.coolingGridParameterValuesDimension1[midpoint]) parameter1Index = midpoint;
      else highpoint = midpoint;
    }
  }

  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


  // interpolate over parameter 1
  for (q = 0;q < 2;q++) {
    // interpolate over temperature
    interpolateIndex = (q+parameter1Index) * CloudyCoolingData.coolingGridNumberTemperatureValues + temperatureIndex;
    
    slope = (dataField[interpolateIndex+1] - dataField[interpolateIndex]) /
      (CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);
    
    temperatureValue[q] = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
      dataField[interpolateIndex];
  }
  slope = (temperatureValue[1]-temperatureValue[0]) /
    (CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index+1]-
     CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]);

  p1Value = (parameter1-CloudyCoolingData.coolingGridParameterValuesDimension1[parameter1Index]) * slope + temperatureValue[0];
  
  return p1Value;
}


// Take value for temperature and interpolate.
double coolingGridInterpolate0D_double(float temperature,double *dataField)
{
  int temperatureIndex;
  int midpoint,highpoint;
  double slope,temperatureValue;


  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


    // interpolate over temperature
    
    slope = (dataField[temperatureIndex+1] - dataField[temperatureIndex]) /
      (CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);
    
    temperatureValue = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
      dataField[temperatureIndex];

    return temperatureValue;
}


// Float returning version of 0D interpolation function.
// Take value for temperature and interpolate.
float coolingGridInterpolate0D_float(float temperature,float *dataField)
{
  int temperatureIndex;
  int midpoint,highpoint;
  float slope,temperatureValue;


  // find temperature index
  if (temperature <= CloudyCoolingData.coolingGridTemperature[0]) {
    temperatureIndex = 0;
  }
  else if (temperature >= CloudyCoolingData.coolingGridTemperature[CloudyCoolingData.coolingGridNumberTemperatureValues - 1]) {
    temperatureIndex = CloudyCoolingData.coolingGridNumberTemperatureValues - 2;
  }
  else {
    temperatureIndex = 0;
    highpoint = CloudyCoolingData.coolingGridNumberTemperatureValues - 1;
    while (highpoint - temperatureIndex > 1) {
      midpoint = (highpoint + temperatureIndex) >> 1;
      if (temperature >= CloudyCoolingData.coolingGridTemperature[midpoint]) temperatureIndex = midpoint;
      else highpoint = midpoint;
    }
  }


  // interpolate over temperature
    
  slope = (dataField[temperatureIndex+1] - dataField[temperatureIndex]) /
    (CloudyCoolingData.coolingGridTemperature[temperatureIndex+1] - CloudyCoolingData.coolingGridTemperature[temperatureIndex]);

  temperatureValue = (temperature - CloudyCoolingData.coolingGridTemperature[temperatureIndex]) * slope +
    dataField[temperatureIndex];

    return temperatureValue;
}
