
#ifdef SMALL_INTS
#define ISYM "d"
#endif

#ifdef LARGE_INTS
#define ISYM "lld"
#endif

#ifdef r4
#define FSYM "f"
#endif

#ifdef r8
#define FSYM "lf"
#endif

#ifdef p8
#define PSYM "lf"
#define GSYM "g"
#endif

#ifdef p16
#define PSYM "Lf"
#define GSYM "g"
#endif
