# `UCNS3D.DAT` — input parameter reference

`UCNS3D.DAT` is the main control file of UCNS3D. It selects the equations, the numerical methods, the flow conditions and the output settings of a run.

This page lists every parameter in the order it appears in the file, with all the values it accepts and a short explanation of each.

## Contents

- [How the file is laid out](#how-the-file-is-laid-out)
- [Quick reference: which line holds what](#quick-reference-which-line-holds-what)
- [1. General setup](#1-general-setup)
- [2. Flow conditions](#2-flow-conditions)
- [3. Spatial discretisation](#3-spatial-discretisation)
- [4. Time integration](#4-time-integration)
- [5. Boundaries, mesh scaling and gradients](#5-boundaries-mesh-scaling-and-gradients)
- [6. Run length and output](#6-run-length-and-output)
- [7. Restart information](#7-restart-information)
- [8. Probes](#8-probes)
- [Optional companion files](#optional-companion-files)
- [Notes](#notes)

## How the file is laid out

The file is read **by line position**, not by keyword.

- Values are read from fixed line numbers. The lines in between are headers and separators that the code skips.
- **Never add or delete a line.** You can change the text of a header line freely, but the number of lines before each value line must stay the same.
- Values on one line are separated by spaces.
- Probe coordinates, if any, start on line 68, one probe per line.

## Quick reference: which line holds what

| Line | Values, in order |
|---:|---|
| 5 | `dimensiona` `statistics` `code_profile` |
| 8 | `governingequations` `initcond` |
| 11 | `turbulence` `icoupleturb` `passivescalar` |
| 18 | `rres` `ufreestream` `vvel` `wvel` `pres` \[`Mach_in`\] |
| 21 | `aoa` `vectorx` `vectory` `vectorz` |
| 24 | `gamma` `prandtl` `reynolds` `charlength` |
| 30 | `spatiladiscret` `iriemann` `spatialorder` `limiter` `poly` |
| 33 | `wenocnschar` `ees` `wenoz` `wenocentralweight` |
| 36 | `temporder` `cfl` `timestep` `upperlimit` `reslimit` |
| 39 | `iboundary` `boundtype` `scaler` |
| 42 | `greengo` `lmach` |
| 48 | `out_time` `ntmax` `wallc` |
| 51 | `tecplot` `ievery` `ievery2` `ieveryav` `stencil_io` |
| 54 | `averaging` |
| 57 | `outsurf` `iforce` `surfshear` |
| 60 | `ires_turb` `ires_unsteady` `lamps` `prev_turbmodel` |
| 66 | `nprobes` |
| 68+ | probe coordinates, one probe per line |

---

## 1. General setup

### `dimensiona` — dimensions of the problem (line 5)

| Value | Meaning |
|---|---|
| `2` | 2D problem |
| `3` | 3D problem |

### `statistics` — parallel scalability statistics (line 5)

| Value | Meaning |
|---|---|
| `0` | Statistics collection disabled (default) |
| `1` | Statistics collection enabled. This is expensive; use it only for performance studies |

### `code_profile` — code configuration (line 5)

A code profile is a named set of internal settings (stencil sizes, implicit solver, robustness options and so on). It lets you switch the whole configuration with one number.

| Value | Meaning |
|---|---|
| `0` | Default settings |
| `1` | Jacobi implicit solver, set up for DDES |
| `11` | Matrix-free LU-SGS implicit solver, set up for DDES |
| `9` | Jacobi implicit solver, robust settings for RANS |
| `91` | Matrix-free LU-SGS implicit solver, robust settings for RANS |
| `777` | Robust settings (as `9`) |
| `888` | Default settings with first-order reconstruction in the cells next to boundaries |
| `100`, `101` | Hybrid DG–FV method |
| `102` | Pure DG method |

Other profiles can be added or edited in `parameters.f90`.

### `governingequations` — equations to solve (line 8)

| Value | Meaning |
|---|---|
| `1` | Navier–Stokes |
| `2` | Euler |
| `3` | Linear advection equation |
| `4` | Linear advection equation, second test variant (used for gradient-approximation tests) |
| `-1` | Multicomponent (multi-phase) Euler equations |

### `initcond` — initial condition profile (line 8)

| Value | Meaning |
|---|---|
| `4` | Free-stream values everywhere (default). Use this for problems with inflow, outflow or pressure far-field boundaries |
| other | A numbered test-case profile, for example the 2D explosion or the Taylor–Green vortex. Profiles are defined in `profile.f90`, where you can add your own |

### `turbulence` — turbulence model (line 11)

| Value | Meaning |
|---|---|
| `0` | Deactivated (laminar or inviscid) |
| `1` | Active. The model itself (Spalart–Allmaras or k-ω SST) is chosen by the code profile |

### `icoupleturb` — coupling of the turbulence model (line 11)

| Value | Meaning |
|---|---|
| `0` | Decoupled: the turbulence equations are solved separately from the mean flow (default) |
| `1` | Coupled: the turbulence equations are solved together with the mean flow |

### `passivescalar` — number of passive scalars (line 11)

| Value | Meaning |
|---|---|
| `0` | No passive scalars |
| `1, 2, …, N` | Number of passive scalars to transport. Only the first one is written to the output files |

---

## 2. Flow conditions

### Free-stream state (line 18)

| Parameter | Meaning | Values |
|---|---|---|
| `rres` | Free-stream density | Any positive value |
| `ufreestream` | Free-stream velocity, x component | Any value |
| `vvel` | Free-stream velocity, y component | Any value |
| `wvel` | Free-stream velocity, z component | Any value |
| `pres` | Free-stream pressure | Any positive value, or `-1` (see below) |
| `Mach_in` | Prescribed inlet Mach number | Optional sixth value; leave it out if not needed |

Setting `pres = -1` sets the pressure to `rres / gamma`. The speed of sound is then 1, so `ufreestream` is the Mach number.

### `aoa` — angle of attack (line 21)

Any value, in degrees.

### `vectorx`, `vectory`, `vectorz` — plane of the angle of attack (line 21)

These three values select the plane in which the angle of attack is measured.

| Values | Plane |
|---|---|
| `1 1 0` | XY |
| `1 0 1` | XZ |
| `0 1 1` | YZ |

### Fluid properties and reference values (line 24)

| Parameter | Meaning | Values |
|---|---|---|
| `gamma` | Ratio of specific heats, used for single-component fluids | Any value (1.4 for air) |
| `prandtl` | Laminar Prandtl number | Any value (0.72 for air) |
| `reynolds` | Reynolds number, `Re = rres × ufreestream × charlength / viscosity` | Any positive value |
| `charlength` | Characteristic length | Any positive value |

The Reynolds number is used only to set the free-stream viscosity. If you know the Reynolds number per metre, multiply it by `charlength`.

---

## 3. Spatial discretisation

### `spatiladiscret` — type of scheme (line 30)

| Value | Meaning |
|---|---|
| `1` | Central (linear) scheme without limiter. For smooth flows only |
| `2` | MUSCL with a slope limiter (default) |
| `3` | WENO variants |

### `iriemann` — approximate Riemann solver (line 30)

| Value | Meaning |
|---|---|
| `1` | HLLC. Resolves contact and shear waves; the usual choice for viscous and low-speed flows |
| `2` | Rusanov (local Lax–Friedrichs). The most dissipative and the most robust |
| `3` | Roe |
| `4` | Rotated Roe–HLL hybrid, designed to avoid the carbuncle at strong shocks |
| `5` | T-Roe (3D only) |
| `8` | Hybrid HLL–HLLC: HLL in and next to shocks, HLLC elsewhere. For strong shocks combined with boundary layers |
| `9` | HLL. Robust at strong shocks, but smears contact and shear layers |
| `12` | Roe with half dissipation (3D, ideal gas only) |

Roe-type solvers (`3`, `4`, `5`, `12`) are for ideal-gas cases; use HLLC, Rusanov, HLL or the hybrid with real-gas or multicomponent cases.

### `spatialorder` — order of spatial accuracy (line 30)

Any integer from `1` to `7`.

### `limiter` — slope limiter for MUSCL (line 30)

The limiter matters even with central or WENO schemes. Some cells may not have enough directional stencils, or may be switched back to MUSCL by the MOOD technique, and those cells use this limiter.

| Value | Meaning |
|---|---|
| `1` | Minmod, equivalent to Barth–Jespersen |
| `2` | MOG limiter |
| `3` | MOGE limiter |
| `4` | MOGV limiter |
| `5` | Van Albada |
| `6` | Van Leer |
| `7` | Venkatakrishnan |
| `8` | Venkatakrishnan with smooth-region blending |
| `9` | MOG variant with extended-stencil blending |
| `10` | Restrictive Barth–Jespersen |
| `11` | Vector Barth–Jespersen |
| `21`–`31` | The same limiters (value minus 20), evaluated at the cell vertices instead of the face quadrature points. For example `27` is Venkatakrishnan at the vertices |

### `poly` — basis-function polynomial (line 30)

| Value | Meaning |
|---|---|
| `1` | Generic polynomials: x, y, z, x², y², z², xy, yz, xz, … (default) |
| `2` | Legendre polynomials, shifted to the interval 0 to 1 |
| `4` | DG / modal basis |

### `wenocnschar` — reconstruction variables (line 33)

| Value | Meaning |
|---|---|
| `1` | Conserved variables (default) |
| `2` | Characteristic variables. Works only with WENO-type schemes |
| `3` | Primitive variables. Suitable for multicomponent flows |

### `ees` — directional stencil algorithm (line 33)

| Value | Meaning |
|---|---|
| `0` | Default stencils |
| `1` | Restrictive stencils |
| `2` | Symmetrical stencils |
| `5` | Compact stencils. **Required** to activate the compact WENO (CWENO) and CWENOZ schemes |

### `wenoz` — weight normalisation (line 33)

Applies to WENO schemes only.

| Value | Meaning |
|---|---|
| `0` | Default weights. With `ees = 5` this gives the CWENO variant |
| `1` | Z-type weights. With `ees = 5` this gives the CWENOZ variant |

### `wenocentralweight` — linear weight of the central stencil (line 33)

Any positive value. Higher values favour the central stencil, which suits smooth problems.

| Scheme | Suggested range |
|---|---|
| CWENO | 10³ to 10⁶ (10⁵ works across many problems) |
| CWENOZ | 2 to 100 |
| WENO | 100 to 10⁵ |

---

## 4. Time integration

### `temporder` — time discretisation method (line 36)

| Value | Meaning | CFL limit |
|---|---|---|
| `1` | Forward Euler | below 1.0 |
| `2` | 2nd-order Runge–Kutta (SSP) | below 1.0 |
| `3` | 3rd-order Runge–Kutta (SSP) | below 1.0 |
| `4` | 4th-order Runge–Kutta (SSP) | below 1.5 |
| `5` | Forward Euler with local time stepping, for steady-state problems | below 1.0 |
| `10` | Implicit backward Euler, for steady-state problems | no formal limit |
| `11` | Implicit dual time stepping, 2nd order, for unsteady problems | no formal limit |
| `12` | Explicit dual time stepping, 2nd order, for unsteady problems | below 1.0 |

### `cfl` — CFL number (line 36)

Any value allowed by the chosen time method.

For dual time stepping (`11`, `12`) the CFL number applies only to the pseudo-steady problem solved at each physical time step. For the implicit options (`10`, `11`) a large value speeds up convergence; reduce it if the run diverges.

### `timestep` — physical time-step size (line 36)

Any positive value. Used only by the dual time stepping options (`11`, `12`).

### `upperlimit` — maximum inner iterations (line 36)

The maximum number of pseudo-time iterations per physical time step in dual time stepping. If the inner problem has not converged by then, the code moves on to the next time step.

For option `12`, a value of 20 is normally enough for a three-orders-of-magnitude drop in the inner residual.

### `reslimit` — residual convergence criterion (line 36)

The normalised residual below which the run (or the inner problem of dual time stepping) is considered converged.

| Time method | Suggested value |
|---|---|
| `11`, `12` | about 0.001 is enough for a wide range of problems |
| `5`, `10` | a value near 0.00001 may be needed |

---

## 5. Boundaries, mesh scaling and gradients

### `iboundary` — periodic boundaries (line 39)

| Value | Meaning |
|---|---|
| `0` | No periodic boundaries |
| `1` | The domain has periodic boundaries |

### `boundtype` — far-field type (line 39)

| Value | Meaning |
|---|---|
| `0` | Supersonic |
| `1` | Subsonic. By default the far-field state is determined automatically inside the code |

### `scaler` — mesh scaling (line 39)

Any positive value. The grid coordinates are **divided** by this number; use `1` if the mesh is already in the units you want.

### `greengo` — gradient approximation for the diffusive fluxes (line 42)

| Value | Meaning |
|---|---|
| `0` | Least squares (default). Poor-quality cells automatically use Green–Gauss instead |
| `1` | Green–Gauss in all cells |

This setting affects only the gradients used for the viscous (diffusion) fluxes.

### `lmach` — low-Mach-number correction (line 42)

| Value | Meaning |
|---|---|
| `0` | No correction |
| `1` | Correction active. It mainly improves low-order MUSCL and WENO schemes up to 3rd order; artefacts may appear with higher-order methods |

---

## 6. Run length and output

### Stopping criteria (line 48)

The run stops when the first of these limits is reached.

| Parameter | Meaning | Values |
|---|---|---|
| `out_time` | Simulation time at which to finish. It carries the units of the flow variables and the grid | Any positive value |
| `ntmax` | Maximum number of iterations | Any positive integer |
| `wallc` | Wall-clock time limit in seconds. A checkpoint file and an output file are written when it is reached | Any positive value |

On a restart the iteration counter continues from the restart file, so `ntmax` must be larger than the iteration stored in that file.

### `tecplot` — output file format (line 51)

| Value | Meaning |
|---|---|
| `0` | Tecplot ASCII |
| `1` | Tecplot binary, one file for the whole domain |
| `2` | VTK binary, one file for the whole domain |
| `3` | VTK binary, partitioned output |
| `4` | Tecplot binary, partitioned output |
| `5` | VTU binary, one file for the whole domain. Use only for meshes up to about 50 million cells |
| `6` | pVTU binary, partitioned parallel output |

### Output frequencies (line 51)

| Parameter | Meaning |
|---|---|
| `ievery` | How often to write a solution file, in wall-clock seconds. For unsteady code profiles it is the interval in simulation time instead |
| `ievery2` | How often to write a restart (checkpoint) file, in wall-clock seconds |
| `ieveryav` | How often to write an averaged solution file, in wall-clock seconds |

Because these are wall-clock intervals, they control how many files are written in a given run time, whatever the speed of the machine.

### `stencil_io` — write stencils (line 51)

| Value | Meaning |
|---|---|
| `0` | Disabled |
| `1` | Write the reconstruction stencils at each probe location |

### `averaging` — time averaging (line 54)

| Value | Meaning |
|---|---|
| `0` | Deactivated (default) |
| `1` | Compute the mean and RMS of the solution. Use only for unsteady ILES, DES, DDES or URANS simulations |

### Surface output and forces (line 57)

| Parameter | Value | Meaning |
|---|---|---|
| `outsurf` | `0` / `1` | Write the surface solution on wall boundaries: off / on |
| `iforce` | `0` / `1` | Compute forces (CL, CD): off / on |
| `surfshear` | `0` / `1` | Include shear stresses in the surface output: off / on |

---

## 7. Restart information

These four values on line 60 describe the **restart file being read**, not the new run. They tell the code what the file contains.

| Parameter | Value | Meaning |
|---|---|---|
| `ires_turb` | `0` / `1` | Previous run was without / with a turbulence model |
| `ires_unsteady` | `0` / `1` | Previous run was steady / unsteady |
| `lamps` | `0` / `1` | No passive scalar / a passive scalar is present in the restart file |
| `prev_turbmodel` | `0` / `1` / `2` | Previous turbulence model: none / Spalart–Allmaras / k-ω |

A run restarts automatically when a restart file is present in the run folder.

---

## 8. Probes

### `nprobes` — number of probes (line 66)

Any non-negative integer. Probes record the solution at fixed points during the run.

If `nprobes` is greater than 0, give the coordinates of every probe starting on line 68, one probe per line (`x y` in 2D, `x y z` in 3D). For example, four probes need four lines: 68, 69, 70 and 71.

---

## Optional companion files

These files are read from the run folder only if they are present.

| File | Purpose |
|---|---|
| `THERMAL.DAT` | Makes no-slip walls isothermal. It holds a flag and the wall temperature in kelvin. Without it, walls are adiabatic |
| `AXISSYMMETRIC.DAT` | Switches a 2D case to axisymmetric and gives the position of the axis |
| `RESTART.dat` | Solution to restart from |

## Notes

1. **Mesh formats.** The solver reads meshes in [ANSYS Fluent](https://www.ansys.com/en-gb/products/fluids/ansys-fluent) and [STAR-CD](https://www.plm.automation.siemens.com/global/en/products/simcenter/STAR-CCM.html) formats. The STAR-CD format is deprecated by Siemens and has been integrated into STAR-CCM+.
2. **Probes.** Make sure the coordinates of every probe are written; the code expects exactly `nprobes` coordinate lines.
3. **Output rates.** All output intervals are in wall-clock seconds, except where stated above.
4. **Inflow, outflow and far-field problems.** Specify the free-stream conditions and use initial condition profile `4`.
5. **Units.** The total simulation time has the units implied by the flow variables and the grid dimensions.
6. **Custom configurations.** New code profiles can be created by editing the parameter-reading routine (`parameters.f90`).
7. **Custom initial conditions.** New initial solution profiles can be added in `profile.f90` and given a number. Examples of well-known problems, such as the 2D explosion and the Taylor–Green vortex, are included.
8. **Test problems.** A set of test cases is provided for getting familiar with UCNS3D.
