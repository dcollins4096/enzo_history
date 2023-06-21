#ifndef __AMRUCDGRIDHIERARCHY_HH_
#define __AMRUCDGRIDHIERARCHY_HH_

#include "FlexArrayTmpl.H"
#include "OrderedList.H"
#include "AmrNode.hh"
#include "Bounds.hh"
#include "IO.hh"
#include "AmrGrid.h"
#include "Vec3f.hh"

class AmrUcdGridHierarchy {
  struct AmrUcdGrid {
    int level;
    int rank;
    int dims[3];
    double origin[3],delta[3];
    IObase::DataType datatype;
    union {
      void *vdata;
      int *idata;
      float *fdata;
      double *ddata;
      short *sdata;
      char *cdata;
    };
    float *data; // local data
    // Attributes
    long nelements;
    int gridID,complete;
    // DATATYPE data; // assume f77 order (should be void*)
    int dataveclen;
    Bounds bounds;
    FlexArray<AmrNode> points; // assume f77 order
    AmrUcdGrid *parent;
    FlexArray<AmrUcdGrid*> children;
    // Constructors
    AmrUcdGrid(AmrGrid &g):complete(0),parent(0){
      copyFromAmrGrid(g);
      regenPoints();
    }
    ~AmrUcdGrid(){
      if(fdata!=data) delete data;
      data=0;
    }
  private:
    // utilities
    void copyFromAmrGrid(AmrGrid &g);
    void regenPoints();
  };
  struct LevelInfo {
    double delta[3]; // find differences in delta
    int childstride[3]; // stride through children
    void setDelta(double *dx){
      for(int i=0;i<3;i++){
	// childstride[i]=1; // why?
	delta[i]=dx[i]; // CRASH !!!
      }
      printf("setDelta=%lf\n",delta[0]);
    }
    void setStride(double *childdx){
      for(int i=0;i<3;i++)
	childstride[i] = (int)(delta[i]/childdx[i] + 0.5);
      
      printf("Stride=%u from dx=%lf\n",childstride[0],childdx[0]);
    }
  };
public:
  typedef OrderedList<AmrUcdGrid*> GridList;
  FlexArray<GridList> grids;
  
  FlexArray<LevelInfo> levelinfo;
// public:
  ~AmrUcdGridHierarchy();
  void addGrid(AmrGrid &grid); // copy from
  //*** Now do parenting of the nodes (can be combined with the addGrid
#if 0 // disabled
  struct CellInfo {
    struct Face {
      Vec3f normal;
      // int neighbor[3]; return as -1,0,1 values
      Vec3f minext[3];
      Vec3f maxext[3];
    };
    Face faces[6];
    float Vec3f[3],Vec3f[3];
  };
  static inline AmrUcdGrid *Grid2Cell(int coord[3],AmrUcdGrid *gridinfo ,CellInfo *cellinfo){
    // automatically hoists grid reference if intersection with denser region
    // automatically deprecates grid reference if intersection with boundary
    // returns null if encounter with outer boundary
    // Ordering of faces are {+x -x +y -y +z -z}

    // actually totally wrong... need to compute cell min max!!!!!!
    Vec3f min(gridinfo->bounds.min),max(gridinfo->bounds.max);

    cellinfo->face[0].minext.set(min[0],min[1],min[2]); // m000
    cellinfo->face[0].maxext.set(min[0],max[1],max[2]); // M011
    cellinfo->face[1].minext.set(max[0],min[1],min[2]); // m100
    cellinfo->face[1].maxext.set(max[0],max[1],max[2]); // M111
    cellinfo->face[2].minext.set(min[0],min[1],min[2]); // m000
    cellinfo->face[2].maxext.set(max[0],min[1],max[2]); // M101
    cellinfo->face[3].minext.set(min[0],max[1],min[2]); // m010
    cellinfo->face[3].maxext.set(max[0],max[1],max[2]); // M111
    cellinfo->face[4].minext.set(min[0],min[1],min[2]); // m000
    cellinfo->face[4].maxext.set(max[0],max[1],min[2]); // M110
    cellinfo->face[5].minext.set(min[0],min[1],max[2]); // m001
    cellinfo->face[5].maxext.set(max[0],max[1],max[2]); // M111
    
    Vec3f p1(min[0],min[1],min[2]),// (0,0,0),
      p2(max[0],min[1],min[2]),// (1,0,0),
      p3(max[0],max[1],min[2]),// (1,1,0),
      p4(min[0],max[1],min[2]),// (0,1,0),
      p5(min[0],min[1],max[2]),// (0,0,1),
      p6(max[0],min[1],max[2]),// (1,0,1),
      p7(max[0],max[1],max[2]),// (1,1,1),
      p8(min[0],max[1],max[2]); // (0,1,1)
    cellinfo->face[0].normal.norm(p1,p4,p8);
    cellinfo->face[1].normal.norm(p3,p2,p6);
    cellinfo->face[2].normal.norm(p2,p1,p5);
    cellinfo->face[3].normal.norm(p4,p3,p7);
    cellinfo->face[4].normal.norm(p1,p4,p3);
    cellinfo->face[5].normal.norm(p8,p7,p6);
  }
#endif // disabled
  void purge();
  void buildNodeHierarchy();	
  void print();
  void buildUCD(FlexArray<AmrNode*> &nodes,FlexArray<int> &cells);
};


#endif

