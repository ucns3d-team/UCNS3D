#!/bin/sh
#
# UCNS3D : An Open Source, High-Order, Finite-Volume CFD Solver
#
# Main contact:  Panagiotis (Takis) Tsoutsanis (Cranfield University)
#
# Licensed under the GNU General Public License, Version 3.0 (the "License");
# you may not use this file except in compliance with the License.
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "as is" basis,
# without warranties or conditions of any kind, either express or implied.
# see the License for the specific language governing permissions and
# limitations under the License.
#
# --------------------------------------------------------------------
#
# Makefile for building UCNS3D PROTOTYPE on HPC systems 
#
# Load correct Programming Environment for build (see below)
#


#> ===== Building on Cirrus

#> Executable and source files

EXECUTABLE=test
SOURCES=prototype_cpu.f90

#> Programming Environment

# --- Nvidia
# module load gcc/10.2.0
# module avail nvidia/nvhpc-nompi/24.5
# module openmpi/4.1.6-cuda-12.4-nvfortran
# module load intel-20.4/cmkl 

# --- GNU
# module load gcc/10.2.0
# module load intel-20.4/cmkl

#> Compiler flags and libs

# --- Nvidia
# Note: -mp=multicore for multicore CPU offload, -mp-gpu for GPU offload

#BLASFLAGS= -I${MKLROOT}/include
#FFLAGS= -mp -mp=multicore -Minfo -Mnovect -r8 -Mbackslash -march=native $(BLASFLAGS) # $(MPIFLAGS)
#FFLAGS= -Mcray=pointer $(BLASFLAGS)

# --- GNU
# Macro definition picks "include mpif.h" instead of "use mpi" when compiling under GNU on Cirrus (see Cirrus user guide) 
CPPFLAGS= -cpp -D__CIRRUS_GNU__

# Note: switching on LTO can slow down compilation
FFLAGS= -fopenmp -fimplicit-none -march=native -g -finit-local-zero -ffree-line-length-none -fbackslash -fdefault-real-8 -fdefault-double-8 -Wno-lto-type-mismatch -Wno-argument-mismatch -flto -fno-tree-sra
#fno-lto

#-fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch -Wno-argument-mismatch
#-fno-lto -fno-tree-sra

BLASLIBS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group
LIBS= -Wl,-Bdynamic -lstdc++ -lpthread -lm -ldl -lc -lmpi $(BLASLIBS) 

#> Compiler wrapper
F90=mpif90   # -v # verbose output


all: $(EXE)

$(EXE):
	$(F90) $(CPPFLAGS) $(FFLAGS) $(SRC) $(LIBS) -o $(EXE)

clean:
	rm $(EXE) *.o *.mod
