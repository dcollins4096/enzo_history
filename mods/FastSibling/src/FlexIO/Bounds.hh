#ifndef __BOUNDS_HH_
#define __BOUNDS_HH_

struct Bounds {
  double min[3],max[3];
  void setFromExtents(double *p1,double *p2);
  void setFromOriginDx(double *origin,double *dx,int *dims);
  int contains(double point[3]);
  // uses restricted assumption that box is entirely inside or 
  // entirely outside its parent
  inline int contains(Bounds &b){
    if(contains(b.min) && contains(b.max)) return 1;
    else return 0;
  }
};

#endif

