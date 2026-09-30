##intel
#F90=mpiifort
#GFORTRAN
#F90=mpifort
#Archer2 fortran compiler wrapper
F90=amdflang
LD=$(F90)
# ifdef HPCTOOLKIT
#    $(info HPCToolkit is loaded - building version for profiling)
#    DBGFLAGS=-g -debug inline-debug-info
# else
#    $(info HPCToolkit is not loaded - building release version)
#    DBGFLAGS=
# endif
#
EXTRALIBSPATH=/scratch/project_465003108/HACKATHON/UCNS3D/libst

#EXTRALIBSPATH=/pfs/lustrep1/scratch/project_465001996/cce_build/local/lib

#EXTRALIBSPATH=/work/e01/e01/takist/DG_RUN_TGV/GPUX/finalito/libs/
#EXTRALIB2=/work/e01/e01/takist/DG_RUN_TGV/GPUX/TEC/360ex_2025r2/bin/

#intel MKL
#BLASFLAGS=  -I${MKLROOT}/include
#BLASLIBS=  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

#BLASFLAGS= -I  ${BLIS_ROOT}/include     #Ensure that you include your BLAS library path here
#BLASLIBS=  ${BLIS_ROOT}/lib/libblis.a   #Ensure that you include your BLAS library path here


#amdblis aocl

#cray blas
#We can rely on compiler wrappers on Archer2 to link to the correct blas libraries
#BLASFLAGS=  -I    /opt/cray/pe/libsci/22.12.1.1/CRAY/9.0/x86_64/include/
#BLASLIBS=/opt/cray/pe/libsci/22.12.1.1/CRAY/9.0/x86_64/lib/libsci_cray.a
#BLASLIBS=/opt/cray/pe/libsci/22.12.1.1/CRAY/9.0/x86_64/lib/libsci_cray_mp.a

OFLAGS=
#intel compilers
# debug run
#FFLAGS=-i4 -r8 -O1  -xHost -ipo -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static $(BLASFLAGS) -check bounds -fp-model strict -debug inline-debug-info -zero
# production run
#FFLAGS=-i4 -r8 -O2 -xHost -ipo -fp-model precise -zero -qopenmp  -qopenmp-link=static $(BLASFLAGS)


#gfortran
#debug run
#FFLAGS=-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch $(BLASFLAGS)
#production run
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch -fallow-argument-mismatch $(BLASFLAGS)

#crayftn
#debug run
#FFLAGS=-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch -fallow-argument-mismatch $(BLASFLAGS)
#production run
#FFLAGS=-eZ -s real64 -fbackslash -fopenmp -e 0 -e I -O2 -Dxpu 
#FFLAGS=  \
  -Mpreprocess \
  -fdefault-real-8 \
  -Mnobackslash \
  -fimplicit-none \
  -fopenmp --offload-arch=gfx90a \
  -O1 -fno-lto \
  -Dgpk


 # $mpich$(which amdflang) -O3 -g \
            -fopenmp -fopenmp-version=51 --offload-arch=gfx90a \
            ${t}.f90 -o ${t}_amd.exe \
            -lomp -lmpifort \
            -I$MPICH_DIR/include \
            $PE_MPICH_GTL_DIR_amd_gfx90a $PE_MPICH_GTL_LIBS_amd_gfx90a

FFLAGS=-cpp -fdefault-real-8 -fbackslash \
         -fopenmp --offload-arch=gfx90a\
         -O3 -Dxpu \
	 -lomp -lmpifort \
            -I${MPICH_DIR}/include \
            ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a}

#FFLAGS=-Mpreprocess -fdefault-real-8 -Mnobackslash -Mrecursive -fopenmp -fno-lto -O1 -Dgpu

#FFLAGS=-cpp -fdefault-real-8 -fopenmp  -Ofast  -Dxpu 
        


# Optional compile-time upper bounds for fixed GPU work arrays.
# By default these are derived from UCNS3D.DAT before compilation.
# Override manually with, e.g. make GPU_MAX_FLAGS="-DGPU_MAX_IORDER=5 -DGPU_MAX_SPECIES=5".
PYTHON ?= python3
GPU_MAX_CONFIG ?= UCNS3D.DAT
GPU_MAX_FLAGS_AUTO := $(shell $(PYTHON) gpu_max_flags.py $(GPU_MAX_CONFIG))
GPU_MAX_FLAGS ?= $(GPU_MAX_FLAGS_AUTO)
GPU_MAX_DEPS := gpu_max_flags.py
ifneq ($(wildcard $(GPU_MAX_CONFIG)),)
GPU_MAX_DEPS += $(GPU_MAX_CONFIG)
endif
FFLAGS += $(GPU_MAX_FLAGS)

#aocl
#FFLAGS= -O3 -Mfreeform -Mallocatable=03 -Mdaz -fdefault-real-8 -fopenmp -mavx2 -march=znver2 $(BLASFLAGS)
#FFLAGS= -O2 -mp -fdefault-real-8 -flto -fopenmp -mavx2 -march=znver2 $(BLASFLAGS)

#gnu hybrid static dynamic
#LIBS = -Wl,-Bstatic libparmetis_gnu_mpi.a libmetis.a libtecio.a $(BLASLIBS) -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi
###-Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc

#intel static linking
#LIBS = libGKlib.a libparmetis.a libmetis.a libtecio.a $(BLASLIBS) -lstdc++
#LIBS = -Wl,-Bstatic ./libtecio.a ./libparmetis.a ./libmetis.a -Wl,-Bdynamic \
        -lpthread -lm -ldl 

#LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic ./libtecio.a -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lmetis -lparmetis
LINKFLAGS=-fopenmp


#LIBS = -Wl, -L$(EXTRALIB2) -L$(EXTRALIBSPATH) -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lmpi -ltecio -lparmetis -lmetis
#LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic libtecio.a libparmetis.a libmetis.a -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl 

LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic libtecio.a -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lmpi -lmetis -lparmetis -L/opt/cray/pe/cce/15.0.0/cce/x86_64/lib/


include Makefile_common

