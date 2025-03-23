##intel
#F90=mpiifort
#GFORTRAN
F90=ftn
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
#FFLAGS=-g -fbacktrace -ffree-line-length-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -fcheck=all -Wno-lto-type-mismatch -fallow-argument-mismatch


#this is the most stable for reproducability (without the combination of (-fno-lto -fno-tree-sra) no reproducible results beyond 00)
#FFLAGS= -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none  -O2 -march=native -fcray-pointer -fno-lto -fno-tree-sra -Wno-lto-type-mismatch -fallow-argument-mismatch 

FFLAGS= -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch -fallow-argument-mismatch

#aocl
#FFLAGS= -O3 -Mfreeform -Mallocatable=03 -Mdaz -fdefault-real-8 -fopenmp -mavx2 -march=znver2
#FFLAGS= -O2 -mp -fdefault-real-8 -flto -fopenmp -mavx2 -march=znver2

#gnu hybrid static dynamic
LIBS = -Wl,-Bstatic libparmetis_gnu_mpi.a libmetis.a libtecio.a  -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi
###-Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc 

#intel static linking
#LIBS = libparmetis.a libmetis.a libtecio.a  -lstdc++


include Makefile_common


 
