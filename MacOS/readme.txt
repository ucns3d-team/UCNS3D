. Use the following dynamic libraries:
libmetis.dylib  (Compiled from source)
libparmetis.dylib (Compiled from source)
libtecio.dylib (Copied from Tecplot executable)
The above are also provided in this folder but probably won't work, you will have to compile on your system.

--------------------------------------------------
Locate the required libraries:

1)	Download and install the latest free trial of tecpot, then the dynamic library (libtecio.dylib), is available in 360’s folder, Tecplot 360 EX 2021 r1 /Contents/Frameworks/    (you can run ucns3d without tecplot, and you can use a vtk output but this is required for compiling only)

2)	Download the latest version of parmetis and in the CMakeLists.txt of the parmetis directory, in the section starting with: if(SHARED) you should add:
set(METIS_LIBRARY_TYPE SHARED) 

a.	Then open a terminal window in the parmetis directory and execute the following commands
i.	make config shared=1
ii.	make
b.	Then you can use the generated libmetis.dylib and libparmetis.dylib from your directory parmetis/build/Darwin…./libmetis/ and parmetis/build/Darwin…./libparmetis/ respectively

3)	Install OpenBlas with Homebrew and modify the Makefile to point to the correct location of the libraries.

--------------------------------------------------------


2. Use the Makefile.common and Makefile for MacOS from this folder and copy to the CODE directory

3. Open a terminal window and compile as:
make -f Makefile clean all (this is for new make and clean)
or
make -f Makefile


4. Type the following in a terminal window prior to running the application

install_name_tool -change @rpath/libtecio.dylib /Users/Username/code_directory/libtecio.dylib /Users/Username/executable_directory/executable_name



5. Type the following in a terminal window (values for A greater than or equal to 1, and B greater than 1)
export OMP_NUM_THREADS=A
mpirun -np B ./ucns3d_p

6. Enjoy!
