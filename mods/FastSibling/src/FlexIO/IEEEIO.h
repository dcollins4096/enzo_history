#ifndef __IEEEIO_H_
#define __IEEEIO_H_

#include "IO.h"
#include "Arch.h"

IOFile IEEEopen PROTO((char *filename,char *accessname));
IOFile IEEEopenRead PROTO((char *filename));
IOFile IEEEopenWrite PROTO((char *filename));
IOFile IEEEopenAppend PROTO((char *filename));
void IEEEbufferOn PROTO((IOFile fileID,int bufsize));
void IEEEbufferOff PROTO((IOFile fileID));
#endif
