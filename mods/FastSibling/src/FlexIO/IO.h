/* IO.h C header */
#ifndef __IO_H_
#define __IO_H_
#ifndef UNITYPE
#define UNITYPE
#define BYTE 0
#define INT8 0
#define INT16 1
#define INT32 2
#define INT64 3
#define FLOAT32 4
#define FLOAT64 5
#define UCHAR 6
#define UINT8 6
#define UINT16 7
#define UINT32 8
#define UINT64 9
/* special string types */
#define CHAR 10
#define CHAR8 10
#define STRING 10
#define UNICODE 11
#define CHAR16 11
#endif
/* error */
#define TYPEERROR -1

#define IOREAD 0
#define IOWRITE 1
#define IOAPPEND 2
#define IOCREATE 1 /* Same as IOWRITE */

#include "Arch.h"

int IOclose PROTO((IOFile deviceID));
int IOisValid PROTO((IOFile deviceID));
int IOsizeOf PROTO((int datatype));
int IOnElements PROTO((int rank,int *dims));
int IOnBytes PROTO((int datatype,int rank,int *dims));
int IOwrite PROTO((IOFile deviceID,int typeID,int rank,int *dims,void *data));
int IOreadInfo PROTO((IOFile deviceID,int *typeID,int *rank,int *dims,int maxdims));
int IOread PROTO((IOFile deviceID,void *data));
int IOseek PROTO((IOFile deviceID,int dataset_index));
int IOnDatasets PROTO((IOFile deviceID));
int IOwriteAnnotation PROTO((IOFile deviceID,char *annotation));
int IOreadAnnotation PROTO((IOFile deviceID,int index,char *annotation,int maxsize));
int IOreadAnnotationInfo PROTO((IOFile deviceID,int index,int *size));
int IOnAnnotations PROTO((IOFile deviceID));

int IOwriteAttribute PROTO((IOFile deviceID,char *name,int type,Long length,void *data));
int IOreadIndexedAttributeInfo PROTO((IOFile deviceID,int number,
				      char *name,int *type,Long *nelem,int maxnamelen));
int IOreadAttributeInfo PROTO((IOFile deviceID,char *name,int *type,Long *nelem));
int IOreadAttribute PROTO((IOFile deviceID,int number,void *data));
int IOnAttributes PROTO((IOFile deviceID));
int IOreserveChunk PROTO((IOFile deviceID,int typeID,int rank,int *dims));
int IOwriteChunk PROTO((IOFile deviceID,int *chunkdims,int *chunkorigin,void *data));
int IOreadChunk PROTO((IOFile deviceID,int *chunkdims,int *chunkorigin,void *data));
int IOwriteStream PROTO((IOFile deviceID,void *data,int length));
int IOreadStream PROTO((IOFile deviceID,void *data,int length));
int IOpause PROTO((IOFile deviceID));
int IOresume PROTO((IOFile deviceID));
#endif
