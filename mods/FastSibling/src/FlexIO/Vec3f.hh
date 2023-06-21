#ifndef __VEC3F_HH_
#define __VEC3F_HH_

#include <stdlib.h>
#include <math.h>

    //CAPI:verb SqrDistancePt3
#define PF_SQUARE(_x) ((_x)*(_x))
#define PF_MIN2(a,b) ((a) < (b) ? (a) : (b))
#define PF_MAX2(a,b) ((a) > (b) ? (a) : (b))
#define PF_MIN3(a,b,c) ((a) < (b) ? PF_MIN2(a,c) : PF_MIN2(b,c))
#define PF_MAX3(a,b,c) ((a) > (b) ? PF_MAX2(a,c) : PF_MAX2(b,c))
#define PF_MIN4(a,b,c,d) ((a) < (b) ? PF_MIN3(a,c,d) : PF_MIN3(b,c,d))
#define PF_MAX4(a,b,c,d) ((a) > (b) ? PF_MAX3(a,c,d) : PF_MAX3(b,c,d))
#define PF_CLAMP(_x, _lo, _hi) \
        (((_x) < (_lo)) ? (_lo) : (_x) > (_hi) ? (_hi) : (_x))

/*
 * PF_ABS is faster than calling fabsf
 * PF_ABSLT etc are faster than (PF_ABS(x1) < x2)
 */
#define PF_ABS(_x1) ((_x1 < 0) ? -(_x1) : (_x1))
#define PF_ABSLT(_x1,_x2) ((_x1) < (_x2) && -(_x1) < (_x2))
#define PF_ABSGT(_x1,_x2) ((_x1) > (_x2) || -(_x1) > (_x2))
#define PF_ABSLE(_x1,_x2) ((_x1) <= (_x2) && -(_x1) <= (_x2))
#define PF_ABSGE(_x1,_x2) ((_x1) >= (_x2) || -(_x1) >= (_x2))
/*
 * Speed oriented macros
 */
#define PF_PI       3.14159265358979323846f /* F for SP float */
#define PF_PI_D     3.14159265358979323846  /* slower DP for more precision */

#define PF_DEG2RAD(x)   ((x)*PF_PI  /180.0f)
#define PF_DEG2RAD_D(x) ((x)*PF_PI_D/180.0)

#define PF_RAD2DEG(x)   ((x)*180.0f/PF_PI)
#define PF_RAD2DEG_D(x) ((x)*180.0 /PF_PI_D)

#define PF_HUGEVAL  3.40282347e+37f


/* macro for fast square roots */
/* thresholds chosen so it's no worse than pfSqrt() */
#define PF_SQRT1(_x) \
        (((_x) > 0.9996f && (_x) < 1.001f) ? \
          0.5f + 0.5f*(_x) : \
          pfSqrt(_x))

#define PF_1OVERSQRT1(_x) \
        (((_x) > 0.9996f && (_x) < 1.001f) ? \
          1.5f - 0.5f*(_x) : \
          1.0f/pfSqrt(_x))

struct Vec3f {
  //   PFSTRUCT_DECLARE

public:
    float vec[3];
  
public:
  // constructors and destructors
  //CAPI:private
  Vec3f(float _x, float _y, float _z) { set(_x, _y, _z); }
  Vec3f(float *_p) { set(_p); }
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
  void set(float *_p){
    vec[0]=_p[0];
    vec[1]=_p[1];
    vec[2]=_p[2];
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


    float sqrDistance(const Vec3f& _v) const { 
	return (PF_SQUARE(vec[0] - _v[0]) +
		PF_SQUARE(vec[1] - _v[1]) +
		PF_SQUARE(vec[2] - _v[2]));
    }
  void norm(const Vec3f& _v1, const Vec3f& _v2, const Vec3f& _v3){
    Vec3f s1,s2;
    s1 = _v1-_v2;
    s2 = _v3-_v2;
    cross(s1,s2);
    normalize();
  }
  float normalize(){
    float len=length();
    if(len==0) this->set(len,len,len);
    vec[0]/=len; vec[1]/=len; vec[2]/=len;
    return len;
  }
  float length() const{
    return sqrt(PF_SQUARE(vec[0]) +
		PF_SQUARE(vec[1]) +
		PF_SQUARE(vec[2]));
  }
    //CAPI:verb DistancePt3
  float distance(const Vec3f& _v) const{
    return sqrt(this->sqrDistance(_v));
  }
  void  cross(const Vec3f&  v1, const Vec3f&  v2){
    float temp[3];

    temp[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
    temp[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
    temp[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);
    set(temp);
  }
  void  cross(const Vec3f&  v2){
    float temp[3];

    temp[0] = (vec[1] * v2[2]) - (vec[2] * v2[1]);
    temp[1] = (vec[2] * v2[0]) - (vec[0] * v2[2]);
    temp[2] = (vec[0] * v2[1]) - (vec[1] * v2[0]);
    set(temp);
  }
  /*CAPI:verb XformVec3
     void xformVec(const Vec3f& _v, const pfMatrix& _m);

    //CAPI:verb XformPt3
     void xformPt(const Vec3f& _v, const pfMatrix& _m);

    //CAPI:verb FullXformPt3
    void fullXformPt(const Vec3f& _v, const pfMatrix& _m);*/

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
  // friend inline Vec3f operator *(const Vec3f& _v, const pfMatrix& _m);
    
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

//inline Vec3f operator *(const Vec3f& _v, const pfMatrix&  _m) {
    // transform as point (w=1), assuming affine transformation
    // i.e. does not use slower dst.xformFullPt().
//    Vec3f dst; dst.xformPt(_v, _m); return dst;
//}

#endif

