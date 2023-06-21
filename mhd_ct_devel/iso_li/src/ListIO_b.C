
//
// List IO that partakes of macros and parameters overloading.
//

#include <stdio.h>
#include <string.h>
#include <cstdlib>

#include "macros_and_parameters.h"

void ReadTokensFromString(char * string, char ** array, char * delim,int size_of_array){
  char * tmp;
  int counter = 0;
  int lenToCopy;
  tmp = strtok(string, delim);
  while( tmp != NULL && counter < size_of_array){
    lenToCopy = strlen(tmp);
    //Ensure the last character isn't a newline.
    if( strstr( tmp, "\n" ) != NULL )
      lenToCopy -= 1;
    array[counter] = new char[lenToCopy+2];
    strncpy(array[counter],tmp,lenToCopy+2);

    //This puts the end of the string where it needs to be.  Odd things
    //were happening on AIX...
    sprintf(array[counter]+lenToCopy,"\0");
    counter++;
    tmp = strtok(NULL,delim);
  }
}

void ReadIntsFromString(char * string, int * array, int size_of_array){
  char * pointer = string;
  for( int counter = 0;  counter < size_of_array  ; counter++)
    array[counter] =  strtol(pointer, &pointer,0);
}
void ReadFloatFromString(char * string, float * array, int size_of_array){
  char * pointer = string;
  for( int counter = 0;  counter < size_of_array  ; counter++)
    array[counter] =  strtod(pointer, &pointer);
}
