#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

EXTERN int trace;
EXTERN int trace_level;
EXTERN FILE *trace_fptr;


void flow1( const char *name );
void flow2( const char *name );
void print_trace( const char *io, const char *name );

