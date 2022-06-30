
# TESTS: Test problems description


Overview
===========================

The test problems include 2D, and 3D using Linear Advection, Euler, Navier-Stokes, RANS equations and Diffuse Interface Model with a
five-equation paradigm.

• All the directores include the grid files (MESHES), the pointwise files for generating the meshes,
the parameter files UCNS3D.DAT (DATS), and sample run directories with some representative outputs.


• It has to be noted that for each test problem a different initial condition is used, and the user should
start using the configurations of the parameter file UCNS3D.DAT provided prior to changing any of the variables
to ensure that the code runs as intended.



TESTS 1
===========================

References for the first set of tests [TESTS1](https://doi.org/10.5281/zenodo.3375432).


2D Double Mach Reflection
----------------------------
2D/EULER/DMR/   test 3.5 from the reference below

[dmr](https://doi.org/10.1016/j.jcp.2019.07.039)


2D Explosion
----------------------------
2D/EULER/EXPLOSION/  test 5.1.1 from the reference below

[explosion](http://www.global-sci.com/intro/article_detail/cicp/7586.html)


2D Kelvin Helmholtz Instability
--------------------------------
2D/EULER/KHI/ test 3.7 from the reference below

[Kelvin-Helmholtz](https://doi.org/10.1016/j.jcp.2019.07.039)


2D Shu-Osher vortex
--------------------------------
2D/EULER/SHU_VORTEX/ test 3.1 from the reference below

[Shu-Osher-vortex](https://doi.org/10.1016/j.jcp.2018.02.009)




2D Linear Advection
--------------------------------
2D/LINEAR_ADVECTION/SMOOTH_PROFILE/ test 3.2 from the reference below

[Linear Advection](https://doi.org/10.1016/j.jcp.2019.07.039)




2D Solid Body Rotation
--------------------------------

2D/LINEAR_ADVECTION/SOLID_BODY_ROTATION/ test 3.3 from the reference below

[Solid Body Rotation](https://doi.org/10.1016/j.jcp.2019.07.039)


2D RANS of NACA0012
--------------------------------

2D/RANS/NACA0012/   test 4.1 from  the reference below

[NACA0012](https://doi.org/10.1016/j.compfluid.2017.01.002)


3D iLES of Taylor-Green Vortex
--------------------------------

3D/ILES/TGV/  test 3.3 from   the reference below

[iLES TGV](https://doi.org/10.1016/j.amc.2018.04.076)


3D RANS of ONERA M6 Wing
--------------------------------

3D/RANS/ONERA_M6/  test 4.4 from  the reference below

[ONERA_M6](https://doi.org/10.1016/j.compfluid.2017.01.002)




TESTS 2
===========================

References for the second set of tests [TESTS2](https://doi.org/10.5281/zenodo.6538622).



2D Shock-Interaction with Water Droplet
--------------------------------

CAVITY/ test 4.5 from the reference below

[droplet-shock](https://doi.org/10.1007/s10915-021-01673-y)


2D Shock-Interaction with Helium Bubble
--------------------------------

HELIUM/ test 4.6 from from the reference below

[helium-shock](https://doi.org/10.1007/s10915-021-01673-y)


2D Shock-induced collabpse of array of bubbles
--------------------------------

BUBBLES/ 2D CASE I  from the reference below

[BUBBLES_ARRAY](https://doi:10.1017/jfm.2020.535)

2D NOH Infinite Strong Shock
--------------------------------

NOH/ 2  test problem from the reference below

[noh](https://doi.org/10.1016/0021-9991(87)90074-X)






