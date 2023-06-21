/* Region structure used for Parallel FFT */

struct region {
  int StartIndex[MAX_DIMENSION];
  int RegionDim[MAX_DIMENSION];
  int Processor;
  float *Data;
};
