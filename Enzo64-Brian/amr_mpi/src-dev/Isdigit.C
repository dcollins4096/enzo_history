#include<stdio.h>
#include<ctype.h>

#include "macros_and_parameters.h"

Eint32 hide_isdigit(Eint32 c)
{
  int i;

//  i = (Eint32) isdigit((unsigned char) c);

  i = 0;
  fprintf(stderr, "What is this shit? %c\n", c);

  if ( c >= '0' && c <= '9' ) i = 1;

  return(i);
}
