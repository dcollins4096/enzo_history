/* Define units in CGS */

#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

/* Units of length is centimetres */

EXTERN float GlobalLengthUnits;

/* Units of mass in grams */

EXTERN float GlobalMassUnits;

/* Units of density in g/cm^3 */

EXTERN float GlobalDensityUnits;

/* Units of time in seconds */

EXTERN float GlobalTimeUnits;

