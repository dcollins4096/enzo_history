/*
   This is just a C wrapper for init_sprng() called from Fortran code.
*/
 
#include <stdio.h>
#include <stdlib.h>

#define SIMPLE_SPRNG    /* simple interface   */
#include "sprng.h"      /* SPRNG header file  */
 
#include "macros_and_parameters.h"
 
/*
  amd    trailing underscore_            nec    trailing underscore_
  dec    trailing underscore_            sgi    trailing underscore_
  intel  trailing underscore_            sun    trailing underscore_
  linux  trailing underscore_
 
  hp     NO trailing underscore          cray   NO trailing underscore
  ibm    NO trailing underscore
 
  Cray uses UPPERCASE names.
*/
 
#if defined(IRIS4) || defined(SUN)   || defined(COMPAQ) || \
    defined(IA64) || defined(LINUX) || defined(NEC) || defined(CRAYX1)
extern "C" void enzo_seed_(int *seed)
#endif
 
#if defined(SP2) || defined(HP)
extern "C" void enzo_seed(int *seed)
#endif
 
{
  Eint32 seed_value;
  Eint32 gtype = 0;
  seed_value = (*seed);
  printf("Seed value %ld\n", seed_value);
  init_sprng(gtype, seed_value, SPRNG_DEFAULT);
  print_sprng();
}
