# UCNS3D Tutorial: PSP (Pressure Sensitive Paint) Rotor Simulation

This tutorial will guide you through the process of setting up and running a Multiple Reference Frame (MRF) simulation for the PSP rotor using the UCNS3D solver.

## 1. Case Overview

The Pressure Sensitive Paint (PSP) rotor is a model scale set of rotor blades designed jointly by the U.S. Army Aviation Development Directorate (ADD) and NASA, based on the High Lift Rotor (HLR) blades.

Key characteristics:
- Four-bladed rotor
- Linear twist: -14 degrees (starting at Blade Station 16.75)
- Main chord length: 5.45 inches
- 30-degree tip sweep
- Tip chord length: 3.27 inches
- Government RC-series airfoils
- Dynamic pressure sensors at r/R = 93% and r/R = 99%

<p align="center">
<img width="912" height="811" src="results/psp_experimental.png">
</p>
 This case is part of the [AIAA Rotorcraft Hover Prediction Workshop (HPW)](https://www.aiaa-hpw.org/)

## 2. File Setup

### 2.1 Mesh Files

You can download the mesh from [AIAA Rotorcraft Hover Prediction Workshop (HPW)](https://www.aiaa-hpw.org/file-share) and create your own compuational domain.

Alternatively you can download [full helicopter pointwise](https://dspace.lib.cranfield.ac.uk/handle/1826/22474 ) or the [periodic mesh.](https://www.dropbox.com/scl/fi/sndjeqp3xdq6y0a827euy/grid.msh?rlkey=fy2fqykmasf31dmccc3vu6v8p&st=b4c07jgq&dl=0)


Required files in working directory:
- `grid.msh`: Main grid file

<p align="center">
<img width="912" height="811" src="results/psp_mesh.png">
</p>

### 2.2 ROTFRAME.dat File

Create a file named `ROTFRAME.dat` with the following content for periodic simulation:

```
!ROTATIONAL REFERENCE FRAME MODE ( 0 DEACTIVE,1 ACTIVE, 2-MRF)
2
!SRF ORIGIN POINT
0.0 0.0 0.0
!SRF VELOCITY
0.0 0.0 0.0
!Periodic: 1: ACTIVE  0: DEACTIVE |\| Angle: (deg) |\| Reference Velocity: (WxR)
1                   90.0            218.072
============= Multiple Reference Frame ==================
1                  !Number of rotors
0.0 0.0 -0.5      !Point 1 (MRF- Bottom cylinder face center)
0.0 0.0 0.5       !Point 2 (MRF- Top cylinder face center)
2.0 120.427       !MRF Radius - MRF Rotational velocity (rad/s)
```

### 2.3 Planform Definition

Create a file to define the blade sections according to Table 1:

```
r/R     Twist(deg)  Chord(in)  Sweep(Deg)  Airfoil
0.12    8.2         5.45       0           -
0.17    8.2         5.45       0           RC(4)-12
0.25    7.01        5.45       0           RC(4)-12
0.65    1.4         5.45       0           RC(4)-12
0.70    0.7         5.45       0           RC(4)-10
0.80    -0.7        5.45       0           RC(4)-10
0.85    -1.4        5.45       0           RC(6)-8
0.95    -2.8        5.45       0           RC(6)-8
1.00    -3.5        3.27       30          RC(6)-8
```

### 2.4 UCNS3D Input File

Create your main UCNS3D input file  with the following settings. Adjust the settings as needed for different simulation runs. For this coarse mesh we will deactivate the Turbulence as the mesh is coarse and with no Boundary layer.

```
====================================================================================================================================================================================================|
----------------------------------------------------------------------------------------UCNS3D PARAMETERS-------------------------------------------------------------------------------------------|
====================================================================================================================================================================================================|
|2:2D, 3:3D ||\|| I(STATISTICS 1=ENABLED, 0=DISABLED)  ||\||    CODE CONFIGURATION (0=DEFAULT)
3				 0   					0
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
EQUATIONS: |1: Navier-Stokes |2: Euler | 3: Linear-sinewave | 4: Linear-step  ||\||  Initial conditions profile (4 DEFAULT)
		1	   								4
|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
TURBULENCE MODEL ACTIVATION:1:Active 0:Deactive   ||COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED  ||SCALAR TRANSPORT COMPUTATION (number of passive scalars [only output for the first one])
0							0							0
```



## 3. Running the Simulation

1. Ensure all required files are in your working directory.
2. Run the UCNS3D solver using the command:
   ```
   mpirun -np <number_of_processes> ./ucns3d_p
   ```
3. Monitor the simulation progress through the output files (Force.dat, Moment.dat. and residual.dat)

## 4. Post-processing

The post-processing data for the full helicopter out ground effect can be downloaded at [https://dspace.lib.cranfield.ac.uk/handle/1826/22472](https://dspace.lib.cranfield.ac.uk/handle/1826/22472)

<p align="center">
<img width="912" height="811" src="results/psp_IGE_OGE.jpg">
</p>

### 4.1 Pressure Distribution Analysis

Focus on pressure measurements at key stations:
- r/R = 93% (pressure sensor location))

<p align="center">
<img width="912" height="811" src="results/psp_medium_Cp_rr_093
.pdf">
</p>
#- r/R = 99% (pressure sensor location)
- r/R = 93% (pressure sensor location))

<p align="center">
<img width="912" height="811" src="results/psp_medium_Cp_rr_099
.pdf">
</p>



### 4.2 Performance Parameters
Calculate and compare:
- Thrust coefficient
- Torque coefficient
- Figure of Merit
- Blade loading distribution

<p align="center">
<img width="912" height="811" src="results/FoMxCt.pdf">
</p>

## 5. Validation Data

Compare results with experimental data from:
1. NASA 14x22 Rotor Test Cell data
2. PSP validation test data
3. Available pressure sensor measurements at r/R = 93% and 99%


## 6. References

1. PSP Rotor Workshop documentation
2. Original NASA/ADD design documentation
3. Experimental validation datasets

By following this tutorial, you should be able to set up and run a simulation of the PSP rotor using UCNS3D with the Multiple Reference Frame approach and periodic boundary conditions.
