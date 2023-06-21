/* Macro definitions for HDF4 compatibility */

#if defined(SP2)
typedef void         *VOIDP;
typedef float        float32;
typedef double       float64;
#endif

#if defined(SUN)
typedef void         *VOIDP;
typedef int          int32;
typedef float        float32;
typedef double       float64;
#endif

#if defined(IRIX) || defined(IRIS4) || defined(COMPAQ) || defined(IA64)
typedef void         *VOIDP;
typedef int          int32;
typedef float        float32;
typedef double       float64;
#endif
