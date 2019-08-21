Readme file for the UCNS3D.DAT parameter file


The user should select the respective options for the configuration of
the code that is desired such as the type of equations to be solved, 
numerical methods, output frequencies etc.


Notes:

1) the use must ensure to rite the coordinates for all the probe positions that are desired.
so if 4 probe positions are required the code will expect to read 4 lines with the x,y,z coordinates of the probes at lines 68,69,70, and 71 of the dat file

2) the output rates of the output files are all given in wallclock time (seconds), therefore controlling how many files will be written in a given time

3)for problems with inflows/outflows or pressure farfields the free stream conditions need to be specified, along with initial conditions 4 profile

4) the total simulation time is carrying the units of the flow variables and grid dimensions as well.

5) the user can develop various code configurations by modifying the READ_UCNS3D subroutine at the io.f90 module

6) the user can create initial solution profiles and assign them a number at the profile.f90 module, where some examples of well established problems are included such as 2D explosion, Taylor Green Vortex etc.

7) a series of test problems is provided for testing and getting familiar with ucns3d



