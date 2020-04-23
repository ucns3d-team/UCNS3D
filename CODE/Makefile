##intel fortran
F90=mpiifort

## gfortran 
#F90=mpifc

LD=$(F90)
ifdef HPCTOOLKIT
   $(info HPCToolkit is loaded - building version for profiling)
   DBGFLAGS=-g -debug inline-debug-info
else
   $(info HPCToolkit is not loaded - building release version)
   DBGFLAGS=
endif

ifdef MKLROOT


###!GFORTRAN WITH MKL
#BLASFLAGS=-I${F95ROOT}/include/intel64/lp64 -m64 -I${MKLROOT}/include
#BLASLIBS=${F95ROOT}/lib/intel64/libmkl_blas95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

##INTEL FORTRAN WITH MKL
BLASFLAGS=-I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
BLASLIBS= ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl

else
   $(error MKLROOT not set. Typing 'module load mkl' may help)
endif

OFLAGS= 
##INTEL COMPILER INSTRUCTIONS
#FFLAGS=-i4 -r8 -O0  -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static $(BLASFLAGS) -check bounds -fp-model strict -debug inline-debug-info -zero 
FFLAGS=-i4 -r8 $(DBGFLAGS) -O2 -ipo -xHost -fp-model fast=2 -zero -qopenmp  -qopenmp-link=static $(BLASFLAGS)   #NEW WAY 30% FASTER

##GFORTRAN COMPILE INSTRUCTIONS
#FFLAGS=-ffree-line-length-none -fdefault-real-8 -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch -DGNU -static $(BLASFLAGS)

LIBS =$(BLASLIBS) -lstdc++ 

include Makefile.common


 
