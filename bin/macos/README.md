# Setup on MacOS

Use the following dynamic libraries:
libmetis.dylib  (Compiled from source)
libparmetis.dylib (Compiled from source)
libtecio.dylib (Copied from Tecplot executable)
The above are also provided in this folder but probably won't work, you will have to compile on your system.

--------------------------------------------------
Locate the required libraries:

1. Download and install the latest free trial of tecpot, then the dynamic library (libtecio.dylib), is available in 360’s folder, Tecplot 360 EX 2021 r1 /Contents/Frameworks/    (you can run ucns3d without tecplot, and you can use a vtk output but this is required for compiling only)

2. Download the latest version of parmetis and in the CMakeLists.txt of the parmetis directory, in the section starting with: if(SHARED) you should add:
set(METIS_LIBRARY_TYPE SHARED) 

    2.1. Then open a terminal window in the parmetis directory and execute the following commands

```
make config shared=1
```
```
make
```

2.2. Then you can use the generated libmetis.dylib and libparmetis.dylib from your directory parmetis/build/Darwin…./libmetis/ and parmetis/build/Darwin…./libparmetis/ respectively

3. Use the Makefile and Makefile_common for MacOS from this folder and copy them to the src directory.
   The build does not require BLAS, LAPACK, OpenBLAS, or MKL.

4. Open a terminal window in the src directory and compile as:

```
make -f Makefile clean all
```

```
make -f Makefile
```

If the libraries are not in the src directory, point the Makefile to them:

```
make -f Makefile MAC_LIB_ROOT=/path/to/macos/libs
```

5. Type the following in a terminal window prior to running the application if libtecio.dylib is not found automatically:

```
install_name_tool -change @rpath/libtecio.dylib /Users/Username/code_directory/libtecio.dylib /Users/Username/executable_directory/executable_name
```


6. Type the following in a terminal window (values for A greater than or equal to 1, and B greater than 1)
```
export OMP_NUM_THREADS=A
```
```
mpirun -np B ./ucns3d_p
```
