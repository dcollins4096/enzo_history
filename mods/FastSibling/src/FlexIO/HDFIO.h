#ifndef __HDFIO_H_
#define __HDFIO_H_

#include "IO.h"
#include "Arch.h"

IOFile HDFIOopen PROTO((char *filename,char *accessname));
IOFile HDFIOopenRead PROTO((char *filename));
IOFile HDFIOopenWrite PROTO((char *filename));
IOFile HDFIOopenAppend PROTO((char *filename));

#endif
