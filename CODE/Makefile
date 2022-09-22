##intel
F90=mpiifort
#GFORTRAN
#F90=mpif90
LD=$(F90)
# ifdef HPCTOOLKIT
#    $(info HPCToolkit is loaded - building version for profiling)
#    DBGFLAGS=-g -debug inline-debug-info
# else
#    $(info HPCToolkit is not loaded - building release version)
#    DBGFLAGS=
# endif


#intel MKL
BLASFLAGS=  -I${MKLROOT}/include
BLASLIBS=  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

#BLASFLAGS= -I  ${BLIS_ROOT}/include     #Ensure that you include your BLAS library path here
#BLASLIBS=  ${BLIS_ROOT}/lib/libblis.a   #Ensure that you include your BLAS library path here


#amdblis aocl
#BLASFLAGS=  -I/opt/hlrs/non-spack/libraries/aocl-gcc/2.1.0/amd-blis/include/
#BLASLIBS= /opt/hlrs/non-spack/libraries/aocl-gcc/2.1.0/amd-blis/lib/libblis.a

OFLAGS= 
#intel compilers
# debug run
FFLAGS=-i4 -r8 -O0  -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static $(BLASFLAGS) -check bounds -fp-model strict -debug inline-debug-info -zero 
# production run
#FFLAGS=-i4 -r8 -ipo -march=core-avx2 -mtune=core-avx2 -O3 -fp-model precise -zero -qopenmp  -qopenmp-link=static $(BLASFLAGS)  #amd zen
FFLAGS=-i4 -r8 -ipo -xHost -O2 -fp-model precise -zero -qopenmp  -qopenmp-link=static $(BLASFLAGS)


#gfortran
#debug run
#FFLAGS=-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch $(BLASFLAGS)
#production run
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -O3 -march=native -Wno-lto-type-mismatch   $(BLASFLAGS) 

#aocl
#FFLAGS= -O3 -Mfreeform -Mallocatable=03 -Mdaz -fdefault-real-8 -fopenmp -mavx2 -march=znver2 $(BLASFLAGS)
#FFLAGS= -O2 -mp -fdefault-real-8 -flto -fopenmp -mavx2 -march=znver2 $(BLASFLAGS)

#gnu hybrid static dynamic
#LIBS = -Wl,-Bstatic libparmetis.a libmetis.a libtecio.a $(BLASLIBS) -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi 

#intel static linking
LIBS = libparmetis.a libmetis.a libtecio.a $(BLASLIBS) -lstdc++ 


include Makefile.common


 
