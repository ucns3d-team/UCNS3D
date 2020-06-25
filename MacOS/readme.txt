. Use the following dynamic libraries:
libmetis.dylib  (Compiled from source)
libparmetis.dylib (Compiled from source)
libtecio.dylib (Copied from Tecplot executable)
The above are also provided in this folder but probably won't work, you will have to compile on your system.

2. Use the Makefile.common and Makefile for MacOS from this folder and copy to the CODE directory

3. Open a terminal window and compile as:
make -f Makefile clean all (this is for new make and clean)
or
make -f Makefile


4. Type the following in a terminal window prior to running the application

install_name_tool -change @rpath/libtecio.dylib /Users/Username/code_directory/libtecio.dylib /Users/Username/executable_directory/executable_name


5. Type the following in a terminal window (values for A and B is down to you but B>1)
export OMP_NUM_THREADS=A
mpirun -np B ./ucns3d_p

6. Enjoy!
