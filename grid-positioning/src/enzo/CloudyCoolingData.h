/***********************************************************************
/
/  Cloudy Cooling Data
/
***********************************************************************/

#define CLOUDY_COOLING_MAX_DIMENSION 5

struct CloudyCoolingDataType
{

  // Use a CMB temperature floor.
  int CMBTemperatureFloor;

  // Redshift independent temperature floor.
  float ConstantTemperatureFloor;

  // Flag to control whether or not to include heating from Cloudy.
  int IncludeCloudyHeating;

  // Flag to control whether or not to include mean molecular weight from Cloudy.
  int IncludeCloudyMMW;

  // To convert from mass fraction to metallicity.
  /*
    x = SUM { A_i * m_i}, for i = 3 to N.
    A_i = solar number abundance with respect H.
    m_i = atomic weight.
    N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), x = 0.018477.
   */
  float CloudyMetallicityNormalization;

  // Factor to account for extra electrons from metals.
  /* 
     f = SUM { A_i * i }, for i = 3 to N.
     N = Atomic number of heaviest element in cooling model.
     For solar abundance patters and N = 30 (Zn), f = 9.153959e-3.
   */
  float CloudyElectronFractionFactor;

  // Cooling grid file.
  char *CloudyCoolingGridFile;

  // Rank of Cloudy dataset.
  int CloudyCoolingGridRank;

  // Dimension of Cloudy dataset.
  int *CloudyCoolingGridDimension;

  // Dataset parameter values.
  float **CloudyCoolingGridParameters;

  // Heating values
  float *CloudyHeating;

  // Cooling values
  float *CloudyCooling;

  // Array holding mean molecular weight values
  float *CloudyMeanMolecularWeight;

};
