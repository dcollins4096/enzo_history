/***********************************************************************
/
/  WRITE HDF5 STRING ATTRIBUTE
/
/  written by: Robert Harkness
/  date:       July 2002
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAILS HARD
/
************************************************************************/

#ifdef USE_HDF5
#include <hdf5.h>
#include <stdio.h>
#include <assert.h>

#include "macros_and_parameters.h"

// HDF5 function prototypes

#include "extern_hdf5.h"




int HDF5_WriteStringAttr(hid_t dset_id, char *Alabel, char *String, FILE *log_fptr)
{

  hid_t       attr_id;
  hid_t       attr_dsp_id;
  hid_t       attr_type_id;
  herr_t      h5_status;
  herr_t      h5_error = -1;

  const char  *NoString = "none";

#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  if (io_log) fprintf(log_fptr, "Enter WSA\n");

//  if (io_log) fprintf(log_fptr, "  Alabel: %d %s\n", strlen(Alabel), Alabel);
//  if (io_log) fprintf(log_fptr, "  String: %d %s\n", strlen(String), String);

  attr_dsp_id = H5Screate(H5S_SCALAR);
    if (io_log) fprintf(log_fptr, "  H5Screate attr_dsp_id: %d\n", attr_dsp_id);
    assert( attr_dsp_id != h5_error );

  attr_type_id = H5Tcopy(H5T_C_S1);
                 H5Tset_size(attr_type_id, 80);

  attr_id = H5Acreate(dset_id, Alabel, attr_type_id,  attr_dsp_id, H5P_DEFAULT);
    if (io_log) fprintf(log_fptr, "  H5Acreate attr_id: %d\n", attr_id);
    assert( attr_id != h5_error );

  if( strlen(String) > 0 )
  {
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) String);
  }
  else
  {
    h5_status = H5Awrite(attr_id, attr_type_id, (void *) NoString);
  }

    if (io_log) fprintf(log_fptr, "  H5Awrite: %d\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Aclose(attr_id);
    if (io_log) fprintf(log_fptr, "  H5Aclose: %d\n", h5_status);
    assert( h5_status != h5_error );

  h5_status = H5Sclose(attr_dsp_id);
    if (io_log) fprintf(log_fptr, "  H5Sclose: %d\n", h5_status);
    assert( h5_status != h5_error );

  return SUCCESS;

}
#endif /* USE_HDF5 */
