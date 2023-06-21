#======================================================================
#
# Make.fft for old F77 FFT
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-12-15  Created
#
#======================================================================

TEXT_FFT = $(TEXT_FFT_F77)

DEFINES_FFT = -DFFT_F77
OBJS_FFT = fourn.o power_of_2.o


