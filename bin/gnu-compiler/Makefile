## Makefile to be build locally using GNU compilers.

##intel
# F90=mpiifort
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





#gfortran
#debug run
FFLAGS=-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch 
#production run
# FFLAGS=-fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -O3 -march=native -Wno-lto-type-mismatch   

#aocl
#FFLAGS= -O3 -Mfreeform -Mallocatable=03 -Mdaz -fdefault-real-8 -fopenmp -mavx2 -march=znver2 
#FFLAGS= -O2 -mp -fdefault-real-8 -flto -fopenmp -mavx2 -march=znver2 

#gnu hybrid static dynamic
LIBS = -Wl,-Bstatic libparmetis.a libmetis.a libtecio.a -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi 

#intel static linking
# LIBS = libparmetis.a libmetis.a libtecio.a -lstdc++ 


include Makefile_common


 
