/***********************************************************************
/
/  EXCEPTION CLASS
/
/  written by: Matthew Turk
/  date:       June, 2009
/  modified:   
/
/  PURPOSE:
/
************************************************************************/
#ifndef __ENZO_EXCEPTIONS__
#define __ENZO_EXCEPTIONS__

#ifndef OSX10_4
#include <execinfo.h>
#else
#include "execinfo_local.h"
#endif
#include <signal.h>
#include <stdio.h>

#include <exception>

// Example from 
//
// http://www.ibm.com/developerworks/linux/library/l-cppexcep.html
//

// This must be included BEFORE macros_and_parameters.h
// so we use int here

extern char current_error[255];

 class EnzoFatalException
 {
 public:
     EnzoFatalException(const char *error_msg,
                        const char *filename = NULL,
                        int line_number = 0)
     {
         void * array[25];
         int nSize = backtrace(array, 25);
         char ** symbols = backtrace_symbols(array, nSize);
         fprintf(stderr, "Caught fatal exception:\n\n");
         fprintf(stderr, "   '%s'\n", error_msg);
         if(filename != NULL)
            fprintf(stderr, "at %s:%d\n\n", filename, line_number);
         fprintf(stderr, "Backtrace:\n\n");

         for (int i = 0; i < nSize; i++)
         {
             fprintf(stderr, "BT symbol: %s\n", symbols[i]);
         }

         delete [] symbols;
     }
    void WriteDebuggingOutput()
    {

    }
 };

#endif
