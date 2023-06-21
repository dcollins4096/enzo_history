#
# macros for c3
#
CPP           = /lib/cpp
DEFINES       = -DCONVEX
OPTIONS       = -DVECLIB

## Fortran compiler
FC            = fc
FC_OPTIMIZE   = -O2 -g -fpp

## C Compiler
CC            = cc
CC_OPTIMIZE   = -O1 -g

## C++ Compiler
CCxx	      = CC
CCxx_OPTIMIZE = -O1 -g +w

## Linker/loader
LDR           = fc
LDR_OPTIMIZE  = -g -O2
LIB           = -L/usr/ncsa/lib -L/usr/convex/c++.3.0.1/lib \
	        -lc -lC -ltask -lm -ldf -lveclib
INCLUDE       = -I/usr/ncsa/include 

OPTS          = -c -DCONVEX ${OPTIONS} ${INCLUDE}
DEBUG         = 

