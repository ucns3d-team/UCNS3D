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

#intel MKL

#ifdef MKLROOT

#BLASFLAGS=  -I${MKLROOT}/include
#BLASLIBS=  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
#else
#   $(error MKLROOT not set. Typing 'module load mkl' may help)
#endif

#open blas

 BLASFLAGS=  -I /usr/local/Cellar/openblas/0.3.10_1/include	###Ensure that you include your BLAS library path here
 BLASLIBS=  /usr/local/Cellar/openblas/0.3.10_1/lib/libopenblas.a   ###Ensure that you include your BLAS library path here

OFLAGS= 

#intel compilers
# debug run
#FFLAGS=-i4 -r8 -O0  -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static $(BLASFLAGS) -check bounds -fp-model strict -debug inline-debug-info -zero 
# production run
#FFLAGS=-i4 -r8 -O3 -ipo -xHost -fp-model precise -zero -qopenmp  -qopenmp-link=static $(BLASFLAGS)  

##
#gfortran
#debug run
#FFLAGS=-ffree-line-length-none -fimplicit-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O1 -fcray-pointer -Wno-lto-type-mismatch $(BLASFLAGS)
#production run
FFLAGS=-fdefault-real-8 -fdefault-double-8 -finit-local-zero -flto -fbackslash -fimplicit-none -fopenmp -O3 -march=native -ffree-line-length-none -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch $(BLASFLAGS)


LIBS =$(BLASLIBS) -lstdc++

include Makefile.common


 
