UCNS3D CODE Instructions


1. Requirements

a) Linux x86-64 (Tested on Redhat, Ubuntu, Centos, Suse)
b) Intel Parallel Studio Version 17 or newer
c) Intel MKL library
d) Tecplot or Paraview for visualisation



2. Compiling

a) Open a terminal window in the CODE directory
b) Type "make -f Makefile clean all" for a clean installation or "make -f Makefile" for recompiling the changed files and their dependencies
c) the name of the executable is ucns3d_p


3. Grids





4. Running

For running ucns3d you will need the following files in a directory of your choice:
4.1 a grid file generated with any software packages exported in Ansys fluent format (ASCII *.msh extension), given the name grid.msh (more on this at section 4 of this file)
4.2 the UCNS3D.DAT parameter file responsible for all the settings of the code (more on README_UCNS3D.txt file)
4.3 the executable ucns3d_p
4.4 For interactively running the code specify the number of threads to be used by typing in the terminal window:
    "export OMP_NUM_THREADS=N", N being the number of threads to be used (use 1 for MPI only mode)
    in the same terminal window run the code by typing:
    "mpirun -np M ./ucns3d_p",  M being the number of MPI processes (at least 2 are required)
4.5 For running at different HPC systems sample scripts are provided




