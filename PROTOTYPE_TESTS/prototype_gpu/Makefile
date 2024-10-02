
#GFORTRAN
#F90=mpifort
#Archer2 fortran compiler wrapper
F90=ftn

#Cray compiler

prototype.x:
	$(F90) -s real64 -fbackslash -fopenmp -e 0 -e I -Ofast prototype_gpu1.f90 -o prototype.x


#GNU compiler
#$(F90) -fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -finit-local-zero -fimplicit-none -flto -fcray-pointer -Ofast -march=native -Wno-lto-type-mismatch -fallow-argument-mismatch prototype_gpu1.f90 -o prototype.x
