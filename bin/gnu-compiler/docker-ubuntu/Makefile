# UCNS3D: option 6 / 64-bit ParMETIS only.
# Select COMPILER=cray|amd|intel|gnu; clean before changing compilers or flags.
# All partition libraries must be built together with idx_t=64, real_t=32,
# using the solver's MPI environment. Default solver/MPI integers stay 32-bit.
# METIS_LINK=static|shared selects partition-library linkage only.
# METIS_PREFIX=/path/to/install selects its lib directory; alternatively set
# METIS_LIB_DIR directly. See README_PARMETIS64.md for compiler examples.

COMPILER ?= cray
PYTHON ?= python3
# An empty goal list means the default build; "clean all" still needs bounds.
UCNS3D_BUILD_GOALS := $(filter-out clean,$(if $(strip $(MAKECMDGOALS)),$(MAKECMDGOALS),all))
UNAME_S := $(shell uname -s)

# Standard static-library layout, relative to the source/build directory.
UCNS3D_LIB_ROOT ?= ../bin/lib
TECPLOT_LIB ?= $(UCNS3D_LIB_ROOT)/tecplot/libtecio.a
GNU_LIB_ROOT ?= .
GNU_TECPLOT_LIB ?= $(TECPLOT_LIB)
MAC_LIB_ROOT ?= .
MAC_TECPLOT_LIB ?= $(if $(filter static,$(METIS_LINK)),$(TECPLOT_LIB),$(firstword $(wildcard $(MAC_LIB_ROOT)/libtecio.dylib) $(TECPLOT_LIB)))
EXTRALIBSPATH = /pfs/lustrep1/scratch/project_x/cce_build/local/lib

# Kept as a compatibility check for older build commands. The Fortran kind
# is fixed to c_int64_t, even when FFLAGS is overridden or this flag omitted.
METIS_INDEX_BITS ?= 64
ifneq ($(METIS_INDEX_BITS),64)
$(error Only METIS_INDEX_BITS=64 is supported; use matched 64-bit METIS/ParMETIS libraries)
endif

# Preserve static Intel/GNU-Linux and shared Cray/AMD/GNU-macOS defaults.
ifneq ($(filter $(COMPILER),cray amd),)
METIS_LINK ?= shared
else ifeq ($(COMPILER):$(UNAME_S),gnu:Darwin)
METIS_LINK ?= shared
else
METIS_LINK ?= static
endif
ifeq ($(METIS_LINK),static)
METIS_LIB_SUFFIX := a
else ifeq ($(METIS_LINK),shared)
ifeq ($(UNAME_S),Darwin)
METIS_LIB_SUFFIX := dylib
else
METIS_LIB_SUFFIX := so
endif
else
$(error METIS_LINK must be static or shared)
endif

# Explicit installation/directory overrides use a flat library directory.
# Otherwise static libraries use the standard ../bin/lib subdirectories,
# while shared libraries retain the compiler-profile locations.
METIS_PREFIX ?=
ifneq ($(strip $(METIS_PREFIX)),)
METIS_LIB_DIR ?= $(METIS_PREFIX)/lib
else ifeq ($(METIS_LINK),shared)
ifneq ($(filter $(COMPILER),cray amd),)
METIS_LIB_DIR ?= $(EXTRALIBSPATH)
else ifeq ($(COMPILER):$(UNAME_S),gnu:Darwin)
METIS_LIB_DIR ?= $(MAC_LIB_ROOT)
else ifeq ($(COMPILER),gnu)
METIS_LIB_DIR ?= $(GNU_LIB_ROOT)
else
METIS_LIB_DIR ?= .
endif
endif

# Explicit filenames prevent a shared build from silently selecting archives,
# or a static build from accidentally finding system 32-bit shared libraries.
PARMETIS_LIB ?= $(if $(strip $(METIS_LIB_DIR)),$(METIS_LIB_DIR),$(UCNS3D_LIB_ROOT)/parmetis)/libparmetis.$(METIS_LIB_SUFFIX)
METIS_LIB ?= $(if $(strip $(METIS_LIB_DIR)),$(METIS_LIB_DIR),$(UCNS3D_LIB_ROOT)/metis)/libmetis.$(METIS_LIB_SUFFIX)
GKLIB ?= $(if $(strip $(METIS_LIB_DIR)),$(METIS_LIB_DIR),$(UCNS3D_LIB_ROOT)/parmetis)/libGKlib.$(METIS_LIB_SUFFIX)
# Retain existing per-platform path overrides.
GNU_PARMETIS_LIB ?= $(PARMETIS_LIB)
GNU_METIS_LIB ?= $(METIS_LIB)
GNU_GKLIB ?= $(GKLIB)
MAC_PARMETIS_LIB ?= $(PARMETIS_LIB)
MAC_METIS_LIB ?= $(METIS_LIB)
MAC_GKLIB ?= $(GKLIB)
ifeq ($(COMPILER):$(UNAME_S),gnu:Darwin)
PARTITION_LIBRARY_FILES = $(MAC_PARMETIS_LIB) $(MAC_METIS_LIB) $(MAC_GKLIB)
else ifeq ($(COMPILER),gnu)
PARTITION_LIBRARY_FILES = $(GNU_PARMETIS_LIB) $(GNU_METIS_LIB) $(GNU_GKLIB)
else
PARTITION_LIBRARY_FILES = $(PARMETIS_LIB) $(METIS_LIB) $(GKLIB)
endif
# Dependent libraries follow their users: ParMETIS -> METIS -> GKlib.
PARTITION_LIBS = $(PARTITION_LIBRARY_FILES)
ifeq ($(METIS_LINK),shared)
PARTITION_RPATH ?= $(foreach dir,$(sort $(dir $(abspath $(PARTITION_LIBRARY_FILES)))),-Wl,-rpath,$(patsubst %/,%,$(dir)))
PARTITION_LIBS += $(PARTITION_RPATH)
endif

# -----------------------------------------------------------------------------
# Compiler profiles
# -----------------------------------------------------------------------------

ifeq ($(COMPILER),cray)
F90 = ftn
LD  = $(F90)

# Cray Fortran compiler: production run
FFLAGS = -eZ -s real64 -fbackslash -fopenmp -e 0 -e I -O2 -Dxpu

LINKFLAGS = -fopenmp

LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic $(TECPLOT_LIB) -Wl,-Bdynamic \
       $(PARTITION_LIBS) -lstdc++ -lpthread -lm -ldl -lmpi \
       -L/opt/cray/pe/cce/15.0.0/cce/x86_64/lib/

else ifeq ($(COMPILER),amd)
F90 = amdflang
LD  = $(F90)

# AMD Flang: LUMI gfx90a XPU build, matching the AMD Makefile.
# Load the AMD/Cray MPICH environment to provide MPICH_DIR and the GTL flags.
EXTRALIBSPATH = /scratch/project_465003108/HACKATHON/UCNS3D/libst

FFLAGS = -cpp -fdefault-real-8 -fbackslash \
         -fopenmp --offload-arch=gfx90a \
         -O3 -fopenmp-target-fast -Dxpu \
         -I$(MPICH_DIR)/include

LINKFLAGS = -fopenmp

LIBS = -L$(EXTRALIBSPATH) -Wl,-Bstatic $(TECPLOT_LIB) -Wl,-Bdynamic \
       $(PARTITION_LIBS) -lstdc++ -lpthread -lm -ldl -lomp -lmpifort -lmpi \
       $(PE_MPICH_GTL_DIR_amd_gfx90a) $(PE_MPICH_GTL_LIBS_amd_gfx90a)

else ifeq ($(COMPILER),intel)
F90 = mpiifx
LD  = $(F90)

# Intel Fortran compiler: production run
 FFLAGS = -i4 -r8 -ipo -O2 -fpp -march=core-avx2 -fp-model strict \
          -zero -qopenmp -qopenmp-link=static
#FFLAGS =         -i4 -r8 -O0 -g -fpp -debug all -traceback -qopenmp \
 #         -qopenmp-link=static -fp-model strict -debug inline-debug-info -zero

LINKFLAGS = -qopenmp

LIBS = -Wl,-Bstatic $(TECPLOT_LIB) -Wl,-Bdynamic \
       $(PARTITION_LIBS) -lstdc++ -lpthread -lm -ldl

else ifeq ($(COMPILER),gnu)
F90 ?= mpif90
LD  = $(F90)

# GNU Fortran: CPU/OpenMP build, no XPU/OpenMP-target acceleration.
ifeq ($(UNAME_S),Darwin)
GNU_ARCH_FLAGS ?= -mcpu=native
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp \
         -ffree-line-length-none -finit-local-zero -fimplicit-none \
         -O2 $(GNU_ARCH_FLAGS) -fcray-pointer -fno-lto -fno-tree-sra \
         -Wno-lto-type-mismatch -fallow-argument-mismatch

LINKFLAGS = -fopenmp

LIBS = $(MAC_TECPLOT_LIB) $(PARTITION_LIBS) -lstdc++ -lpthread -lm
else
GNU_ARCH_FLAGS ?= -march=native
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp \
         -ffree-line-length-none -finit-local-zero -fimplicit-none \
         -O2 $(GNU_ARCH_FLAGS) -fcray-pointer -fno-lto -fno-tree-sra \
         -Wno-lto-type-mismatch -fallow-argument-mismatch -fPIE

LINKFLAGS = -fopenmp

# Partition files select their own linkage; MPI/system libraries stay dynamic.
LIBS = -Wl,-Bstatic $(GNU_TECPLOT_LIB) -Wl,-Bdynamic \
       $(PARTITION_LIBS) -lstdc++ -lpthread -lm -ldl -lc -lmpi
endif

else
$(error Unknown COMPILER='$(COMPILER)'; expected 'cray', 'amd', 'intel', or 'gnu')
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
# Cleaning needs neither Python nor case data. Explicit flags bypass generation.
ifneq ($(strip $(UCNS3D_BUILD_GOALS)),)
ifeq ($(origin GPU_MAX_FLAGS),undefined)
ifeq ($(wildcard $(GPU_MAX_SCRIPT)),)
$(error Missing GPU_MAX_SCRIPT='$(GPU_MAX_SCRIPT)'; restore the helper or provide GPU_MAX_FLAGS explicitly)
endif
# Use an exit marker rather than .SHELLSTATUS so older GNU make also detects failure.
GPU_MAX_FLAGS_AUTO := $(shell $(PYTHON) "$(GPU_MAX_SCRIPT)" $(if $(strip $(GPU_MAX_CONFIG)),"$(GPU_MAX_CONFIG)") || printf '%s\n' __UCNS3D_BOUNDS_FAILED__)
ifneq ($(filter __UCNS3D_BOUNDS_FAILED__,$(GPU_MAX_FLAGS_AUTO)),)
$(error GPU work-array limit generation failed; use PYTHON=/path/to/python3 (Python 3.6+) and check GPU_MAX_CONFIG)
endif
ifeq ($(strip $(GPU_MAX_FLAGS_AUTO)),)
$(error GPU_MAX_SCRIPT produced no work-array limits; check PYTHON and GPU_MAX_CONFIG)
endif
endif
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


include Makefile_common
