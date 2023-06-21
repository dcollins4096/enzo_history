#ifndef __AMRCELL_HH_
#define __AMRCELL_HH_
typedef float Coord[3];
typedef char ICoord[3];
enum BrickIndex {Min=0,Max,Mid};
#define B_MIN 0
#define B_MAX 1
#define B_MID 2

#define NEGX 0
#define POSX 1
#define MIDX 2
#define NEGY 0
#define POSY 1
#define MIDY 2
#define NEGZ 0
#define POSZ 1
#define MIDZ 2
#if 0
struct AmrCell {
  typedef AmrCell *AmrCellp;
  struct AmrCellp neighbors[6];
  typedef float Coord[3];
  typedef int ICoord[3];
  Coord vertices[2][2][2]; // vertices the form the 8 corners of the brick 
  ICoord normverts[2][2][2];
};

struct AmrCellOps {
  // returns the next ray to intersect
  // or null if we are at the end of the line
  AmrCell *rayIntersect(AmrCell *cell, // current cell
			double entrypoint[3], // exitpoint for prev cell
			double &length, // length of intersection
			double &val, // data value contained by node
			/* later val can take into account interpolation
			/ of neighboring values.  This would combine 
			with the length argument if trilinear */
			double rayroot[3], double rayunitvec[3]);
  //double aperture); // for cone-based tracing.. 
  //must have beam aperture.  We can use the algorithm to trace down the
  // centerline, but must use beam aperture for interpolation.
  // computes which face is intersected by the ray
  // returns integer facenumber of exit + exitpoint
  int faceIntersect(AmrCell *cell,
		    double exitpoint[3], // exitpoint from cell
		    double rayroot[3],double rayunitvec[3]);
};
#endif
//CAPI:struct from Performer
struct Vec3f
{
  //   PFSTRUCT_DECLARE

public:
    float vec[3];
  
public:
  // constructors and destructors
  //CAPI:private
  Vec3f(float _x, float _y, float _z) { set(_x, _y, _z); }
  Vec3f() {};
  
public:
  // sets and gets
  //CAPI:arrayclass
  //CAPI:verb SetVec3
  void set(float _x, float _y, float _z) {
    vec[0] = _x;
    vec[1] = _y;
    vec[2] = _z; 
  }
  void set(float *_x){
    vec[0]=_x[0];
    vec[1]=_x[1];
    vec[2]=_x[2];
  }

public:
    // other functions
    //CAPI:verb
    void copy(const Vec3f&  _v) { *this = _v; }
    int equal(const Vec3f&  _v) const { 
	return (vec[0] == _v[0] && 
		vec[1] == _v[1] &&
		vec[2] == _v[2]);
    }
    int almostEqual(const Vec3f& _v, float _tol) const;

    void negate(const Vec3f& _v) { 
	vec[0] = -_v[0];
	vec[1] = -_v[1]; 
	vec[2] = -_v[2]; 
    }

    float dot(const Vec3f&  _v) const {
	return (vec[0] * _v[0] + 
		vec[1] * _v[1] +
		vec[2] * _v[2]);
    }

    void add(const Vec3f& _v1, const Vec3f& _v2) { 
	vec[0] = _v1[0] + _v2[0]; 
	vec[1] = _v1[1] + _v2[1]; 
	vec[2] = _v1[2] + _v2[2]; 
    }

    void sub(const Vec3f& _v1, const Vec3f& _v2) { 
	vec[0] = _v1[0] - _v2[0]; 
	vec[1] = _v1[1] - _v2[1]; 
	vec[2] = _v1[2] - _v2[2]; 
    }

    void scale(float _s, const Vec3f& _v) { 
	vec[0] = _s * _v[0]; 
	vec[1] = _s * _v[1]; 
	vec[2] = _s * _v[2]; 
    }

    void addScaled(const Vec3f& _v1, float _s, const Vec3f& _v2) { 
	vec[0] = _v1[0] + _s * _v2[0]; 
	vec[1] = _v1[1] + _s * _v2[1]; 
	vec[2] = _v1[2] + _s * _v2[2]; 
    }

    void combine(float _a, const Vec3f& _v1, float _b, const Vec3f& _v2) { 
	vec[0] = _a * _v1[0] + _b * _v2[0]; 
	vec[1] = _a * _v1[1] + _b * _v2[1]; 
	vec[2] = _a * _v1[2] + _b * _v2[2]; 
    }

    //CAPI:verb SqrDistancePt3
    float sqrDistance(const Vec3f& _v) const { 
	return (PF_SQUARE(vec[0] - _v[0]) +
		PF_SQUARE(vec[1] - _v[1]) +
		PF_SQUARE(vec[2] - _v[2]));
    }

    float normalize();
    float length() const;
    //CAPI:verb DistancePt3
    float distance(const Vec3f& _v) const;
    void  cross(const Vec3f&  _v1, const Vec3f&  _v2);

    //CAPI:verb XformVec3
    void xformVec(const Vec3f& _v, const pfMatrix& _m);

    //CAPI:verb XformPt3
    void xformPt(const Vec3f& _v, const pfMatrix& _m);

    //CAPI:verb FullXformPt3
    void fullXformPt(const Vec3f& _v, const pfMatrix& _m);

public:
    // Operators
    float&  operator [](int i) { return vec[i]; }

    const float&  operator [](int i) const { return vec[i]; }

    int operator ==(const Vec3f& _v) const {
        return vec[0] == _v[0] && vec[1] == _v[1] && vec[2] == _v[2];
    }
    int operator !=(const Vec3f& _v) const {
        return !(*this == _v);
    }

public:
    // Vec3f operators (N.B. return by value can be slow)

    Vec3f operator -() const {
        return Vec3f(-vec[0], -vec[1], -vec[2]);
    }

    Vec3f operator +(const Vec3f& _v) const {
        return Vec3f(vec[0]+_v[0], vec[1]+_v[1], vec[2]+_v[2]);
    }

    Vec3f operator -(const Vec3f& _v) const {
        return Vec3f(vec[0]-_v[0], vec[1]-_v[1], vec[2]-_v[2]);
    }

    friend inline Vec3f operator *(float _s, const Vec3f&);
    friend inline Vec3f operator *(const Vec3f& _v, float _s);
    friend inline Vec3f operator /(const Vec3f& _v, float _s);
    friend inline Vec3f operator *(const Vec3f& _v, const pfMatrix& _m);
    
public:
    // Assignment Operators
    Vec3f&  operator =(const Vec3f& _v) {
        vec[0] = _v[0]; 
	vec[1] = _v[1];
	vec[2] = _v[2]; 
	return *this;
    }

    Vec3f& operator *=(float _s) {
        vec[0] *= _s; 
	vec[1] *= _s; 
	vec[2] *= _s; 
	return *this;
    }
    
    Vec3f& operator /=(float _s) {
	_s = 1.0/_s; 
	return *this *= _s;
    }
    
    Vec3f& operator +=(const Vec3f& _v) {
        vec[0] += _v[0]; 
	vec[1] += _v[1]; 
	vec[2] += _v[2];
	return *this;
    }

    Vec3f& operator -=(const Vec3f& _v) {
        vec[0] -= _v[0]; 
	vec[1] -= _v[1]; 
	vec[2] -= _v[2];
	return *this;
    }
};


inline Vec3f operator *(float _s, const Vec3f& _v) {
    return Vec3f(_v[0]*_s, _v[1]*_s, _v[2]*_s);
}

inline Vec3f operator *(const Vec3f& _v, float _s) {
    return Vec3f(_v[0]*_s, _v[1]*_s, _v[2]*_s);
}

inline Vec3f operator /(const Vec3f& _v, float _s) {
    _s = 1.0f/_s;
    return Vec3f(_v[0]*_s, _v[1]*_s, _v[2]*_s);
}

inline Vec3f operator *(const Vec3f& _v, const pfMatrix&  _m) {
    // transform as point (w=1), assuming affine transformation
    // i.e. does not use slower dst.xformFullPt().
    Vec3f dst; dst.xformPt(_v, _m); return dst;
}

struct GridInfo {
  float minext[3],maxext[3];
  int dims[3];
  float *data;
  GridInfo *children,*parent;
  int nchildren;
};

void FindChildren(float minext[3],float maxext[3],FlexArray<AmrGrid> &allgrids,
		  FlexArray<AmrGrid> outgrids){
  outgrids.setSize(0);
  for(int i=0;i<allgrids.getSize();i++){
    if(allgrids[i].minext[0] < minext[0] ||
       allgrids[i].minext[0] < minext[0] ||
       allgrids[i].minext[0] < minext[0] ||
       
  }
}

GridInfo *BuildGridInfoHierarchy(FlexArray<AmrGrid> &grids){
  GridInfo *rootgrid = new GridInfo;
  
}

struct CellInfo {
  struct Face {
    Vec3f normal;
    int neighbor[3];
    Vec3f minext[3];
    Vec3f maxext[3];
  };
  Face faces[6];
  float Vec3f[3],Vec3f[3];
};

int Grid2Cell(int coord[3],GridInfo *gridinfo,CellInfo *cellinfo){
  
}

#endif
