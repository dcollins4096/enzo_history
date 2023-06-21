/* Region structure used for Parallel FFT */

struct region {
  int StartIndex[MAX_DIMENSION];
  int RegionDim[MAX_DIMENSION];
  int Processor;
  float *Data;
};

/*
struct commSndRcv {
 struct region *Sends, *Receives;
 int sends, receives, SendSize, ReceiveSize;
} *cSndRcv;
*/

struct commSndRcv {
 struct region *Sends, *Receives;
 int sends, receives, SendSize, ReceiveSize;
};

/* int first_pass = 0; */
