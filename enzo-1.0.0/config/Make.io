#======================================================================
#
# Make.io for HDF5
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-05-13  Created
#
#======================================================================

TEXT_IO = $(TEXT_IO_HDF5)

DEFINES_IO = -DUSE_HDF5
OBJS_IO = $(OBJS_HDF5)
FLAGS_IO = $(FLAGS_HDF5_ARCH)
LIBS_IO = $(LIBS_HDF5_ARCH)
