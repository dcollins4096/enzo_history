/***********************************************************************
/
/  Cloudy Cooling Data
/
***********************************************************************/

#define MAX_CLOUDY_COOLING_DIMENSION 2

struct CloudyCoolingDataType
{
  // Number of dimension to interpolate over
  int CloudyCoolingGridRank;

  // Flag to control whether or not to include heating 
  // from Cloudy.
  int IncludeCloudyHeating;

  // Maximum number of values for first parameter space
  int coolingGridMaxParameterValuesDimension1;
  // Maximum number of values for second parameter space
  int coolingGridMaxParameterValuesDimension2;
  // Maximum number of temperature points in cooling maps
  int coolingGridMaxTemperatureValues;

  // Cooling grid run file
  char *CloudyCoolingGridRunFile;

  // Values for first parameter space (density)
  float *coolingGridParameterValuesDimension1;
  // Values for second parameter space (metallicity)
  float *coolingGridParameterValuesDimension2;
  // Temperature values
  float *coolingGridTemperature;

  // Number of values for first parameter space
  int coolingGridNumberParameterValuesDim1;
  // Number of values for second parameter space
  int coolingGridNumberParameterValuesDim2;
  // Number of temperature points in cooling maps
  int coolingGridNumberTemperatureValues;

  // Heating values
  float *coolingGridHeating;
  // Cooling values
  float *coolingGridCooling;
  // Array holding mean molecular weight values
  float *coolingGridMeanMolecularWeight;

};
