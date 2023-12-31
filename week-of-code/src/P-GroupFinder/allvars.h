#include "nrsrc/nrutil.h"
#include "ngbtree.h"


#define  KERNEL_TABLE 10000
#define  PI               3.1415927
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS  9.10953e-28
#define  THOMPSON    6.65245e-25

#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

#define GADGET 0
#define ENZO 1

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define nint(A) (int) ((A) + 0.5*sign(A))

extern int     ThisTask, NTask;

extern double  Time;
extern double  BoxSize;
extern double  leftEdge[3], rightEdge[3];

extern double  SearchRadius;

extern int     NumPart;   /* total particle number */
extern int     *NpartInGrids;

extern int     *Nslab, *Nshadow;
extern int     Nlocal, *Noffset;
extern int     Nlocal_in_file;
extern int     *NtoLeft, *NtoRight;
extern int     *Nslab_local, *NtoLeft_local, *NtoRight_local;

extern int     Ncontrib, *ContribID, *ContribHead;

extern int     NLinkAccross;

extern double  GridExtension, GridCorner[3];
extern int     ***GridFirst, ***GridLast, ***GridFlag;
extern int     *GridNext;
extern int     *Head,*Next;   /* for link-lists of groups */
extern int     *Tail,*Len;
extern int     Ngroups, NgroupsAll, *NgroupsList;


extern struct  gr_data
{
  int Len, Tag;
} *GroupDat, *GroupDatAll;

extern int    Nx,Ny,Nz;

/* Quantities for all particles */

extern struct particle_data 
{
  double  	Pos[3];
  float  	Vel[3];
  int    	Type;
  int    	ID;
  int    	MinID;
  int    	GrLen;
  float  	Mass;
  float  	Mfs, Mclouds, Sfr;
  float  	Energy;
  float  	Rho;
  int    	PartID;
  int           slab;
} *P, *Pbuf_local;

struct id_data 
{
  int    ID;
  int    index;
};
 

struct idmin_data 
{
  int    minID;
  int    index;
  int    len;
};
 

extern double  UnitLength_in_cm,
               UnitMass_in_g,
               UnitVelocity_in_cm_per_s,
               UnitTime_in_s,
               UnitTime_in_Megayears,
               UnitDensity_in_cgs,
               UnitPressure_in_cgs,
               UnitCoolingRate_in_cgs,
               UnitEnergy_in_cgs,
               G,
               Hubble,
               EnzoMassUnit,
               EnzoVelocityUnit,
               EnzoTimeUnit;

/* ----------------------------------------------------------
 *
 *  Variables for subfind
 *
 */

extern int    DesDensityNgb;
extern int    DesLinkNgb;
extern float  Theta;
 
extern double Time;
extern double Omega;
extern double OmegaLambda;
extern double G;
extern double H0;
extern float  Epsilon;


extern struct grouptree_data
{
  int lefthead, leftlen;
  int righthead, rightlen;
  int left, right;
} *GroupTree;


extern int   NumInGroup;
extern int   NSubGroups, NSubGroupsAll;

extern int    AnzNodes;
extern int    MaxNodes;


extern int   *Id;
extern int   *SubGroupLen, *SubGroupTag;
extern int   *GroupTag, *NextFinal;
extern int   *Index;
extern int   *NewHead;
extern int   *Node;
extern int   *SortId;
extern float *Density;
extern float *Potential;
extern float *Energy;


/* tabulated smoothing kernel */

extern double  Kernel[KERNEL_TABLE+1],
               KernelDer[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];
