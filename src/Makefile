##intel
#F90=mpiifort
#GFORTRAN
F90=mpif90
LD=$(F90)
# ifdef HPCTOOLKIT
#    $(info HPCToolkit is loaded - building version for profiling)
#    DBGFLAGS=-g -debug inline-debug-info
# else
#    $(info HPCToolkit is not loaded - building release version)
#    DBGFLAGS=
# endif



OFLAGS= 


#gfortran
#debug run
# FFLAGS=-ffree-line-length-none -fimplicit-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero -fbacktrace -fdefault-double-8 -fbackslash  -Og -fcray-pointer -Wno-lto-type-mismatch #-fdump-core -fdump-parse-tree
FFLAGS=-ffree-line-length-none -fimplicit-none -fdefault-real-8 -fopenmp -finit-local-zero -fdefault-double-8  -Og -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch -fbacktrace -fbounds-check
#production run -lstdc++
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -finit-local-zero -fbackslash -fimplicit-none -fopenmp -O2  -flto -march=native -ffree-line-length-none -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch
# FFLAGS=-ffree-line-length-none -fimplicit-none -fdefault-real-8 -fopenmp -finit-local-zero -fdefault-double-8  -O2 -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch

LIBS = libparmetis.dylib libmetis.dylib libtecio.dylib

include Makefile_common


 
