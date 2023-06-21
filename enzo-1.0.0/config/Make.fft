#======================================================================
#
# Make.fft for new F90 FFT
#
# James Bordner (jbordner@cosmos.ucsd.edu)
#
# 2003-12-15  Created
#
#======================================================================

TEXT_FFT = $(TEXT_FFT_F90)

DEFINES_FFT = -DFFT_F90
OBJS_FFT = fft_f90.o fft90.o fft_rotate.o

#  fft90.o contains singleton module used by fft_f90.o--requires this dependency:

fft_f90.o: fft90.o

