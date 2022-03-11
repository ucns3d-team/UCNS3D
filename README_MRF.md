---
title: Multiple Reference Frame on UCNS3D
permalink: tutorials-MRF-UCNS3D.html
keywords: Helicopter, Hovering, rotor, Wind turbine, 
summary: Tutorial for MRF application using UCNS3D
---

This document describes how to run a Multiple Reference Frame simulation on UCNS3D solver. The files for this tutorial are located in the Cranfield University repository Cord  (https://cord.cranfield.ac.uk/)


# Setup
 
The file ROTFRAME.dat contains the inputs used by the solver on the simulation. It includes the flag for the rotational reference frame, axis origin, rotational velocity, rotational peridiocity flag, periodic angle, reference velocity (tip blade velocity), number of rotors and points and radius  of the MRF domain. 


Detailed Parameter Values Settings
===========================

ROTATIONAL REFERENCE FRAME MODE: customisable reference frame settings
-----------------------------------------------------
 ```
        POSSIBLE VALUES= 0 --> DEACTIVE
                         1 --> Single Rotating Frame
                         2 --> Multiple Rotating Frame
```
# Single Rotating Frame
SRF ORIGIN POINT: coordinates of the rotation center for SRF case
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  xyz coordinate e.g.(0.0 0.0 0.0)
```
SRF VELOCITY: rotatinal velocity vector
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  rotational u v w velocity:  e.g: Caradonna (0.0 0.0 265.9)
```
# Rotational Peridiocity
Periodic: Rotational peridiocity (check the rotational_periodic_readme.md for mesh instructions)
-----------------------------------------------------
 ```
        POSSIBLE VALUES= 0 --> Deactivated
                         1 --> Active
```
Angle: rotatinal peridiocity angle
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  Any value   e.g: PSP 4blade(90.0)
```
Reference Velocity: Reference velocity for turbulence calculation (blade tip velocity (WxR))
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  Any value
```
# Multiple Reference Frame

<p align="center">
<img width="1000" height="500" src="docs/MRF_schematic.png">
</p>

Number of rotors: Total number of rotational subdomais
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  Any integer value
```
For each rotor, we have the inputs point1, point2,radius and rotational velocity.  For two or more rotors, you must ensure that you provide the same inputs for each rotor in sequence.

Point1: coordinates of the center of the bottom cylinder face of the rotational domain
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  xyz coordinate e.g.(0.0 0.0 0.0)
```
Point2: coordinates of the center of the Top cylinder face of the rotational domain
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  xyz coordinate e.g.(0.0 0.0 0.0)
```
Radius: Radii of the rotor in metres
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  any value
```
Rotational velocity: Rotational velocity (rad/s)
-----------------------------------------------------
 ```
        POSSIBLE VALUES=  any value
```


The test cases consists of four helicopter cases: 

* Caradonna and Tung (out of ground effect): The mesh files and results of the test case from section 4 used on the publication ["Hovering rotor solutions by high-order methods on unstructured grids"](https://doi.org/10.1016/j.ast.2019.105648)  can be download at the [cranfield reposity cord](https://figshare.com/s/616cc1b31e4cda09ac2c).




* Paint Sensitive Pressure rotor (out of ground effect): The periodic mesh files and results of the test case from section 5 used on the publication ["Hovering rotor solutions by high-order methods on unstructured grids"](https://doi.org/10.1016/j.ast.2019.105648) can be download at the [cranfield reposity cord](https://figshare.com/s/616cc1b31e4cda09ac2c).

* `Black Hawk` UH-60a (out of ground effect): The periodic mesh files and results of the test case from section 6 used on the publication ["Hovering rotor solutions by high-order methods on unstructured grids"](https://doi.org/10.1016/j.ast.2019.105648)  can be download at the [cranfield reposity cord](https://figshare.com/s/616cc1b31e4cda09ac2c).

* Caradonna and Tung (with and without ground effect): The mesh files and results of the test case used on the publication ["Simple multiple reference frame for high-order solution of hovering rotors with and without ground effect"](https://doi.org/10.1016/j.ast.2021.106518)  can be download at the [cranfield reposity cord](https://doi.org/10.17862/cranfield.rd.c.5284412.v1).

* S-76 (with and without ground effect): The mesh files and results of the test case used on the publication ["Simple multiple reference frame for high-order solution of hovering rotors with and without ground effect"](https://doi.org/10.1016/j.ast.2021.106518)   can be download at the [cranfield reposity cord](https://doi.org/10.17862/cranfield.rd.c.5284412.v1).

* Paint Sensitive Pressure full-helicopter (with and without ground effect): The periodic mesh files and results of the test case from section 5 used on the publication ["Numerical Investigation of full helicopter with and without the ground effect"](https://doi.org/10.1016/j.ast.2022.107401)and the  Data underlying this study can be accessed through the [Cranfield University repository](https://doi.org/10.17862/cranfield.rd.
19146182.v1). For this study we used the surface unstructured Blade grid downloaded at (https://www.aiaa-hpw.org/psp-rotor) and the robin fuselage geometry  downloaded at (https://github.com/Applied-Scientific-Research/robin-surface-mesh) whith changes accordingly to the mod7 geometry as exposed at (https://ntrs.nasa.gov/citations/20100033120).


## Running and post processing the Simulation

The detail for running the simulation are provided on instructions given at (https://github.com/ucns3d-team/UCNS3D)




