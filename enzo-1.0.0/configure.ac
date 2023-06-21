#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT(ENZO, 0.2.0, jbordner@cosmos.ucsd.edu)
AC_CONFIG_SRCDIR([uniglue/extern_hdf5.h])

AC_CANONICAL_TARGET
AC_PROG_CC(cc gcc)
# Checks for programs.
# AC_PROG_CXX
# AC_PROG_CC

# Checks for libraries.

# Checks for header files.
# AC_HEADER_STDC
# AC_CHECK_HEADERS([limits.h malloc.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
# AC_C_CONST
# AC_TYPE_SIZE_T
# AC_HEADER_STDBOOL

# Checks for library functions.
# AC_FUNC_MALLOC
# AC_CHECK_FUNCS([sqrt strchr strrchr strstr])

AC_CHECK_LIB(ncurses,initscr,[AC_SUBST(defines_ncurses,[-DUSE_NCURSES])
                              AC_SUBST(libs_ncurses,[-lncurses])])
AC_CHECK_LIB(curses,initscr,[AC_SUBST(defines_curses,[-DUSE_CURSES])
                             AC_SUBST(libs_curses,[-lcurses])])

AC_CONFIG_FILES([Makefile
                 config/Makefile
                 config/MAKE.CONFIG
                 config/MAKE.GLOBAL
                 amr_mpi/Makefile
                 amr_mpi/anyl/Makefile
                 amr_mpi/src/Makefile
                 hdf5_tools/Makefile
                 inits/Makefile
                 ring/Makefile
                 uniglue/Makefile
                 uniglue/dev/Makefile])

AC_CONFIG_COMMANDS(Make.arch,     [touch config/Make.arch])
AC_CONFIG_COMMANDS(Make.mach,     [touch config/Make.mach])
AC_CONFIG_COMMANDS(Make.mach.machine,[touch config/Make.mach.`bin/gethost`])
AC_CONFIG_COMMANDS(Make.mpi,      [touch config/Make.mpi])
AC_CONFIG_COMMANDS(Make.perf,     [touch config/Make.perf])
AC_CONFIG_COMMANDS(Make.io,       [touch config/Make.io])
AC_CONFIG_COMMANDS(Make.precision,[touch config/Make.precision])
AC_CONFIG_COMMANDS(Make.particles,[touch config/Make.particles])
AC_CONFIG_COMMANDS(Make.obj,      [touch config/Make.obj])
AC_CONFIG_COMMANDS(Make.opt,      [touch config/Make.opt])
AC_CONFIG_COMMANDS(Make.endian,   [touch config/Make.endian])
AC_CONFIG_COMMANDS(Make.fft,      [touch config/Make.fft])
AC_CONFIG_COMMANDS(Make.jbopt,    [touch config/Make.jbopt])
AC_CONFIG_COMMANDS(Make.test,     [touch config/Make.test])
AC_CONFIG_COMMANDS(default,       [gmake default])

AC_OUTPUT


