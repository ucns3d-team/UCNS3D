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
SOURCES=prototype_gpu_matmu_no_dtypes.f90

#> Programming Environment

# --- Cray/Archer2
# module load PrgEnv-cray
# module load rocm
# module load craype-accel-amd-gfx90a
# module load craype-x86-milan

#FFLAGS= -s real64 -fbackslash -fopenmp -e 0 -e I -g #-O2 -Ofast

# --- Nvidia
#nvfortran
FFLAGS=  -gpu=mem:managed,maxregcount:64 -Minfo -mp=gpu -r8 -Mbackslash -march=native -Mcray=pointer -Ofast

# --- AMD
# AMD flang ? (GNU?)
#FFLAGS= -Mfree -mp -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -fno-lto -Wno-lto-type-mismatch 

# --- GNU
#FFLAGS= -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch -fallow-argument-mismatch

#> Compiler wrapper

# Check compiler wrapper with: $(F90) -craype-verbose <source file>  | tr " " "\n"
#F90=ftn  

# --- Nvidia/Cirrus
F90=mpifort

all: $(EXECUTABLE)

$(EXECUTABLE):
	$(F90) $(FFLAGS) $(SOURCES) -o $(EXECUTABLE) 


clean:
	rm $(EXECUTABLE) *.o *.mod
