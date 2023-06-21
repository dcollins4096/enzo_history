#include "Bounds.hh"

void Bounds::setFromExtents(double *p1,double *p2){
  for(int i=0;i<3;i++){
    if(p1[i]<p2[i]){ 
      min[i]=p1[i];
      max[i]=p2[i];
    }
    else {
      min[i]=p2[i];
      max[i]=p1[i];
    }
  }
}
void Bounds::setFromOriginDx(double *origin,double *dx,int *dims){
  for(int i=0;i<3;i++){
    min[i]=origin[i];
    if(dims[i]>1)
      max[i]=origin[i]+dx[i]*(double)(dims[i]-1);
    else max[i]=origin[i]; // failsafe for cheating on dims
  }
}
int Bounds::contains(double point[3]) {
  for(int i=0;i<3;i++){
    if(min[i]>point[i]) return 0;
    if(max[i]<point[i]) return 0;
  }
  return 1;
}
