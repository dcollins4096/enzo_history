#ifndef ARRAY_H_
#define ARRAY_H_

template<typename T>
class EnzoArray
{
public:

  EnzoArray(int rank, int *dims, int *start, int *end,
	    float *cell_size=NULL){
    
    int i;
    this->Rank = rank;
    
    for(i = 0; i < MAX_DIMENSION; i++){
      this->Dimension[i] = 0;
      this->StartIndex[i] = 0;
      this->EndIndex[i] = 0;
      this->CellSize[i] = 0.0;
    }
    
    for(i = 0; i < rank; i++){
      this->Dimension[i] = dims[i];
      this->StartIndex[i] = start[i];
      this->EndIndex[i] = end[i];
      if(cell_size)
	this->CellSize[i] = cell_size[i];
    }
    
    Array = NULL;
    
    for(i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++){
      Vector[i] = NULL;
    }
  };

  // No real destructor, because it doesn't allocate data.
  ~EnzoArray(){};

  int Rank;                        // number of dimensions
  int Dimension[MAX_DIMENSION];    // total dimensions of all grids
  int StartIndex[MAX_DIMENSION];   // starting index of the active region
                                   //   (zero based)
  int EndIndex[MAX_DIMENSION];     // stoping index of the active region
                                   //   (zero based)
  float CellSize[MAX_DIMENSION];
  
  T *Array;

  // used for velocities and positions
  T *Vector[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];

};

#define EnzoArrayFLOAT EnzoArray<FLOAT>
#define EnzoArrayFloat EnzoArray<float>
#define EnzoArrayInt   EnzoArray<int>

#endif
