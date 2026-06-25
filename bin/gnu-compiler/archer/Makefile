# UCNS3D GNU/Archer CPU build template
#
# Copy this Makefile and Makefile_common into src, or run from src with:
#   make -f ../bin/gnu-compiler/archer/Makefile
#
# This template defaults to the GNU profile through the Cray ftn wrapper.

COMPILER ?= gnu

PYTHON ?= python3
UNAME_S := $(shell uname -s)

UCNS3D_LIB_ROOT_CANDIDATES ?= ../../bin/lib ../bin/lib
UCNS3D_LIB_ROOT ?= $(firstword $(wildcard $(UCNS3D_LIB_ROOT_CANDIDATES)))
ifeq ($(strip $(UCNS3D_LIB_ROOT)),)
UCNS3D_LIB_ROOT := ../../bin/lib
endif

TECPLOT_LIB ?= $(UCNS3D_LIB_ROOT)/tecplot/libtecio.a
PARMETIS_LIB ?= $(UCNS3D_LIB_ROOT)/parmetis/libparmetis.a
METIS_LIB ?= $(UCNS3D_LIB_ROOT)/metis/libmetis.a

GNU_LIB_ROOT ?= /root/local/lib
GNU_TECPLOT_LIB ?= $(firstword $(wildcard libtecio.a) $(TECPLOT_LIB))
GNU_PARMETIS_LIB ?= $(firstword $(wildcard $(GNU_LIB_ROOT)/libparmetis.a) $(PARMETIS_LIB))
GNU_METIS_LIB ?= $(firstword $(wildcard $(GNU_LIB_ROOT)/libmetis.a) $(METIS_LIB))
GNU_GKLIB ?= $(firstword $(wildcard $(GNU_LIB_ROOT)/libGKlib.a))

MAC_LIB_ROOT ?= .
MAC_TECPLOT_LIB ?= $(firstword $(wildcard $(MAC_LIB_ROOT)/libtecio.dylib) $(wildcard $(MAC_LIB_ROOT)/libtecio.a))
MAC_PARMETIS_LIB ?= $(firstword $(wildcard $(MAC_LIB_ROOT)/libparmetis.dylib) $(wildcard $(MAC_LIB_ROOT)/libparmetis.a))
MAC_METIS_LIB ?= $(firstword $(wildcard $(MAC_LIB_ROOT)/libmetis.dylib) $(wildcard $(MAC_LIB_ROOT)/libmetis.a))

EXTRALIBSPATH = /pfs/lustrep1/scratch/project_x/cce_build/local/lib


# -----------------------------------------------------------------------------
# Compiler profiles
# -----------------------------------------------------------------------------

ifeq ($(COMPILER),cray)
F90 = ftn
LD  = $(F90)

# Cray Fortran compiler: production run
FFLAGS = -eZ -s real64 -fbackslash -fopenmp -e 0 -e I -O2 -Dxpu

LINKFLAGS = -fopenmp

LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic libtecio.a -Wl,-Bdynamic \
       -lstdc++ -lpthread -lm -ldl -lmpi -lmetis -lparmetis \
       -L/opt/cray/pe/cce/15.0.0/cce/x86_64/lib/

else ifeq ($(COMPILER),intel)
F90 = mpiifx
LD  = $(F90)

# Intel Fortran compiler: production run
FFLAGS = -i4 -r8 -ipo -O2 -fpp -march=core-avx2 -fp-model strict \
         -zero -qopenmp -qopenmp-link=static

LINKFLAGS = -qopenmp

LIBS = -Wl,-Bstatic $(TECPLOT_LIB) $(PARMETIS_LIB) $(METIS_LIB) \
       -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl

else ifeq ($(COMPILER),gnu)
F90 ?= ftn
LD  = $(F90)

# GNU Fortran: CPU/OpenMP build, no XPU/OpenMP-target acceleration.
ifeq ($(UNAME_S),Darwin)
GNU_ARCH_FLAGS ?= -mcpu=native
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp \
         -ffree-line-length-none -finit-local-zero -fimplicit-none \
         -O2 $(GNU_ARCH_FLAGS) -fcray-pointer -fno-lto -fno-tree-sra \
         -Wno-lto-type-mismatch -fallow-argument-mismatch

MAC_MISSING_LIBS :=
ifeq ($(strip $(MAC_TECPLOT_LIB)),)
MAC_MISSING_LIBS += libtecio
endif
ifeq ($(strip $(MAC_PARMETIS_LIB)),)
MAC_MISSING_LIBS += libparmetis
endif
ifeq ($(strip $(MAC_METIS_LIB)),)
MAC_MISSING_LIBS += libmetis
endif
ifneq ($(strip $(MAC_MISSING_LIBS)),)
$(error macOS GNU build needs $(MAC_MISSING_LIBS) in MAC_LIB_ROOT='$(MAC_LIB_ROOT)' or explicit MAC_*_LIB paths)
endif

LINKFLAGS = -fopenmp

LIBS = $(MAC_PARMETIS_LIB) $(MAC_METIS_LIB) $(MAC_TECPLOT_LIB)
else
GNU_ARCH_FLAGS ?= -march=native
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp \
         -ffree-line-length-none -finit-local-zero -fimplicit-none \
         -O2 $(GNU_ARCH_FLAGS) -fcray-pointer -fno-lto -fno-tree-sra \
         -Wno-lto-type-mismatch -fallow-argument-mismatch -fPIE

LINKFLAGS = -fopenmp

LIBS = -Wl,-Bstatic $(GNU_TECPLOT_LIB) $(GNU_PARMETIS_LIB) \
       $(GNU_METIS_LIB) $(GNU_GKLIB) \
       -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi
endif

else
$(error Unknown COMPILER='$(COMPILER)'; expected 'cray', 'intel', or 'gnu')
endif


# -----------------------------------------------------------------------------
# Optional compiler recipes
# -----------------------------------------------------------------------------

# Intel Fortran
# FFLAGS = -i4 -r8 -O0 -g -fpp -debug all -traceback -qopenmp \
#          -qopenmp-link=static -fp-model strict -debug inline-debug-info -zero
# FFLAGS = -i4 -r8 -ipo -fpp -xHost -O2 -fp-model precise -zero \
#          -qopenmp -qopenmp-link=static

# GNU Fortran
# FFLAGS = -ffree-line-length-none -fdefault-real-8 -fopenmp -Wall \
#          -fbounds-check -finit-local-zero -fdump-parse-tree -fdump-core \
#          -fbacktrace -fdefault-double-8 -fbackslash -O0 -fcray-pointer \
#          -Wno-lto-type-mismatch -fallow-argument-mismatch
# FFLAGS = -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp \
#          -ffree-line-length-none -finit-local-zero -fimplicit-none -flto \
#          -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch \
#          -fallow-argument-mismatch

# AMD Fortran
# FFLAGS = -Mpreprocess -fdefault-real-8 -Mnobackslash -Mrecursive -fopenmp \
#          -fno-lto -O2 -Dxpu


# -----------------------------------------------------------------------------
# GPU compile-time bounds
# -----------------------------------------------------------------------------

# By default the helper script looks in the current build directory for
# UCNS3D.DAT and companion files such as REALGAS.DAT. If they are absent,
# the helper emits conservative default bounds.
# Override manually with, for example:
#   make GPU_MAX_FLAGS="-DGPU_MAX_IORDER=5 -DGPU_MAX_SPECIES=5"
# Or point the helper at an input deck outside the build directory:
#   make GPU_MAX_CONFIG=/path/to/UCNS3D.DAT
GPU_MAX_CONFIG ?=
GPU_MAX_SCRIPT ?= gpu_max_flags.py
GPU_MAX_FLAGS_AUTO :=
ifneq ($(wildcard $(GPU_MAX_SCRIPT)),)
GPU_MAX_FLAGS_AUTO := $(shell $(PYTHON) $(GPU_MAX_SCRIPT) $(GPU_MAX_CONFIG))
endif
GPU_MAX_FLAGS ?= $(GPU_MAX_FLAGS_AUTO)

GPU_MAX_DEPS :=
ifneq ($(wildcard $(GPU_MAX_SCRIPT)),)
GPU_MAX_DEPS += $(GPU_MAX_SCRIPT)
endif
ifneq ($(strip $(GPU_MAX_CONFIG)),)
ifneq ($(wildcard $(GPU_MAX_CONFIG)),)
GPU_MAX_DEPS += $(GPU_MAX_CONFIG)
endif
else
ifneq ($(wildcard UCNS3D.DAT),)
GPU_MAX_DEPS += UCNS3D.DAT
endif
ifneq ($(wildcard REALGAS.DAT),)
GPU_MAX_DEPS += REALGAS.DAT
endif
ifneq ($(wildcard MULTISPECIES.DAT),)
GPU_MAX_DEPS += MULTISPECIES.DAT
endif
ifneq ($(wildcard MULTISPECIES_DIFF.DAT),)
GPU_MAX_DEPS += MULTISPECIES_DIFF.DAT
endif
endif

FFLAGS += $(GPU_MAX_FLAGS)


# Optional link recipes

# GNU hybrid static/dynamic
# LIBS = -Wl,-Bstatic libparmetis_gnu_mpi.a libmetis.a libtecio.a \
#        -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi

# Intel static linking
# LIBS = libGKlib.a libparmetis.a libmetis.a libtecio.a -lstdc++

# Local static libraries with dynamic system libraries
# LIBS = -Wl,-Bstatic ./libtecio.a ./libparmetis.a ./libmetis.a \
#        -Wl,-Bdynamic -lpthread -lm -ldl

# External METIS/PARMETIS path
# LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic ./libtecio.a -Wl,-Bdynamic \
#        -lstdc++ -lpthread -lm -ldl -lmetis -lparmetis

# Dynamic system/MPI libraries
# LIBS = -L$(EXTRALIBSPATH) -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl \
#        -lmpi -ltecio -lparmetis -lmetis

# Static local libraries with dynamic system libraries
# LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic libtecio.a libparmetis.a libmetis.a \
#        -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl


include Makefile_common
