DETAILED PARAMETERS VALUES

DIMENSIONA = DIMENSIONS OF PROBLEM
           
        POSSIBLE VALUES= 2 --> 2D problem
                         3 --> 3D problem

STATISTICS = PARALLEL SCALABILITY STATISTICS
           
        POSSIBLE VALUES= 0 --> disabled statistics collection
                         1 --> enabled statistics collection (expensive)
                         
 CODE_PROFILE = customisable CODE PROFILE selected
           
        POSSIBLE VALUES= 0  --> default values
                         1  --> jacobi for DDES
                         11 --> matrix free lu-sgs for DDES
                         9  --> jacobi for RANS
                         91 --> matrix free lu-sgs for RANS

                         
governingequations = TYPE OF EQUATIONS TO BE SOLVED
           
        POSSIBLE VALUES= 1 --> Navier-Stokes
                         2 --> Euler
                         3 --> Linear Advection equation
                         4 --> Gradient approximation sample equation
                        -1 --> Multicomponent Euler equations
                        
turbulence=  TURBULENCE MODEL
        POSSIBLE VALUES= 0 --> Deactivated
                         1 --> Active
                         
                         
icoupleturb= COUPLING OF TURBULENCE MODEL
        POSSIBLE VALUES= 0 --> DECOUPLED (DEFAULT)
                         1 --> COUPLED
                         
PASSIVESCALAR= NUMBER OF PASSIVE SCALARS
        POSSIBLE VALUES= 0 --> DEACTIVATED
                         1,2,3,..,N --> As many required, but only the first one is written in output file

RRES=           DENSITY VALUE AT FREE STREAM
        POSSIBLE VALUES=  --> Any positive value
                         
                         
                         
ufreestream=    U FREESTREAM VALUE 
        POSSIBLE VALUES= --> Any value
                         

VVEL=           V FREESTREAM VALUE 
        POSSIBLE VALUES= --> Any value

WVEL=           W FREESTREAM VALUE 
        POSSIBLE VALUES= --> Any value

PRES=           PRESSURE FREESTREAM VALUE 
        POSSIBLE VALUES= --> Any positive value
                         -1 --> It will set pressure at P=RRES/GAMMA, 
                         resulting in SPEED OF SOUND=1, AND Ufreestream=Mach number
                         

AOA=  ANGLE OF ATTACK
        POSSIBLE VALUES=  --> ANY VALUE
                        
                         

vectorx,vectory,vectorz=  SELECT WITH RESPECT TO WHICH AXIS THE ANGLE OF ATTACK IS DEFINED (XY,XZ AND SET ACCORDINGLY THE VALUES)
        
GAMMA=   USED FOR SINGLE COMPONENT FLUIDS
        POSSIBLE VALUES=  --> ANY VALUE
                         
                         

PRANDTL= PRANDTL CONSTANT
        POSSIBLE VALUES= --> ANY VALUE

Reynolds= REYNOLDS NUMBER
        POSSIBLE VALUES= The value is defined as (Re=(RRES*UFREESTREAM*CHARLENGTH)/(VISC)),
        and it is used to determing the freestream value of viscosity

 CharLength= CHARACTERISTIC LENGTH
        POSSIBLE VALUES= --> ANY VALUE
                         

spatiladiscret= DEFINE THE TYPE OF SCHEME THAT WILL BE USED
        POSSIBLE VALUES= 1 --> CENTRAL SCHEME NO LIMITER
                         2 --> MUSCL (DEFAULT)
                         3 --> WENO VARIANTS
                         
                         
iRiemann=       APPROXIMATE RIEMANN SOLVERS
        POSSIBLE VALUES= 1 --> HLLC
                         2 --> RUSANOV(LLF)
                         3 --> ROE
                         4 --> HYBRID ROE-HLL (CARBUNCLE FREE)

spatialorder=  ORDER OF SPATIAL DISCRETISATION
        POSSIBLE VALUES= 1,2,3,..,7 --> SPATIAL ORDER OF ACCURACY
                         

LIMITER=    TYPE OF LIMITER FOR MUSCL-SCHEMES (USEFUL EVEN WHEN USING CENTRAL OR WENO LIMITERS, SINCE
            SOME CELLS MIGHT NOT HAVE SUFFICIENT NUMBER OF DIRECTIONAL STENCILS, OR DUE TO MOOD TECHNIQUE 
            THEY MIGHT REVERT TO MUSCL METHOD)
        POSSIBLE VALUES= 1 --> MINMOD (BARTH AND JESPERSEN EQUIVALENT LIMITER)
                         2 --> MOG    (MOG LIMITER)
                         3 --> MOGE
                         4 --> MOGV
                         5 --> VAN ALBADA
                         6 --> VAN LEER 
                         7 --> VENKATAKRISHNAN

POLY=       BASIS FUNCTION POLYNOMIAL
        POSSIBLE VALUES= 1 --> GENERIC (DEFAULT x+y+z+x^2+y^2+z^2+xy+zy+xz)
                         2 --> LEGENDRE (SHIFTED FROM O TO 1)
                         
                         

wenocnschar= RECONSTRUCTION VARIABLES
        POSSIBLE VALUES= 1 --> CONSERVED (DEFAULT)
                         2 --> PRIMITIVE (SUITABLE FOR MULTICOMPONENT FLOWS)
                         3 --> CHARACTERISTICS (WORKS ONLY FOR WENO TYPE OF SCHEMES)
                         

EES=    DIRECTIONAL STENCILS ALGORITHMS
        POSSIBLE VALUES= 0 --> DEFAULT
                         1 --> RESTRICTIVE
                         2 --> SYMMETRICAL ONES
                         5 --> COMPACT WENO/WENOZ (YOU MUST USE THIS SETTING FOR 
                         ACTIVATING COMPACT WENO/WENOZ SCHEMES)

wenoz=  WEIGHTS NORMALISATION(APPLICABLE TO WENO METHOD ONLY)
        POSSIBLE VALUES= 0 --> DEFAULT (WHEN EES=5 IT ACTIVATES THE CWENO VARIANT)
                         1 --> CWENOZ WHEN EES=5
                         
                         
wenocentralweight=  LINEAR WEIGHT FOR CENTRAL STENCIL
        POSSIBLE VALUES=  --> ANY VALUE 
                            (USE 10^3-10^6 FOR CWENO (higher values more suitable for smooth problems, 10^5 works across many problems))
                            (USE 2-100 FOR CWENOZ)
                            (USE 100-10^5 FOR WENO)
                         

temporder=  TEMPORAL DISCRETISATION METHOD
        POSSIBLE VALUES= 1 --> FORWARD EULER (CFL LIMIT <1.0)
                         2 --> 2ND-ORDER RUNGE-KUTTA (SSP) (CFL LIMIT <1.0)
                         3 --> 3RD-ORDER RUNGE-KUTTA (SSP) (CFL LIMIT <1.0)
                         4 --> 4TH-ORDER RUNGE-KUTTA (SSP) (CFL LIMIT <1.5)
                         5 --> FORWARD EULER WITH LOCAL TIME STEPPING FOR STEADY STATE PROBLEMS  (CFL LIMIT <1.0)
                         10 --> IMPLICIT BDF-EULER FOR STEADY STATE PROBLEMS (NO CFL LIMIT)
                         11 --> IMPLICIT DUAL TIME STEPPING SECOND ORDER FOR UNSTEADY PROBLEMS  (NO CFL LIMIT)
                         12 --> EXPLICIT DUAL TIME STEPPING SECOND ORDER FOR UNSTEADY PROBLEMS  (CFL LIMIT <1.0)

CFL=            CFL NUMBER
        POSSIBLE VALUES= --> ANY VALUE (ACCORDING TO THE TEMPORAL DISCRETISATION METHOD)
                        FOR DUAL TIME STEPPING PROBLEMS THE CFL NUMBER CORRESPONDS ONLY TO THE
                        CFL NUMBER USED FOR THE PSEUDO STEADY-STATE PROBLEM AT EACH NEWTON ITERATION
                        (HENCE A LARGE VALUE SHOULD BE ASSIGNED FOR OPTION 10,11 TO ACCELERATE CONVERGENCE)
                         
                         
                         
timestep=  EXPLICIT DEFINITION OF TIME STEP SIZE  
        POSSIBLE VALUES= --> ANY VALUE
                         USED ONLY BY OPTION (11,12) DUAL TIME STEPPING FOR ADVANCING THE SOLUTION.
                         
                         
upperlimit= UPPER LIMIT OF ITERATIONS FOR PSEUDO-STEADY STATE PART OF DUAL-TIME STEPPING
        POSSIBLE VALUES= --> ANY VALUE (FOR OPTION 12 A VALUE OF 20 IS MORE THAN ENOUGH FOR 
                            CONVERGENCE TO THREE ORDERS OF MAGNITUDE REDUCTION IN RESIDUAL)
                            IF THIS NUMBER OF ITERATIONS IS NOT SUFFICIENT THE DUAL TIME WILL PROCEED TO THE
                            NEXT STEP
        
                          
                         
reSLIMIT=       NORMALISED RESIDUAL CONVERGENCE CRITERION FOR STEADY STATE PROBLEMS (OR PSEUDO-STEADY STATE COMPONENT OF DTS)
        POSSIBLE VALUES= --> ANY VALUE (FOR OPTION 11, 12 A VALUE OF 0.001) IS MORE THAN ENOUGH FOR 
                            A WIDE RANGE OF PROBLEMS
                           FOR OPTION (5, 10) A VALUE CLOSE TO 0.00001 MIGHT BE REQUIRED

iboundary= PRESENCE OF PERIDIC BOUNDARY IN THE DOMAIN
        POSSIBLE VALUES= 1 --> PERIODIC
                         0 --> NON PERIODIC 


boundtype=
        POSSIBLE VALUES= 0 --> SUPERSONIC
                         1 --> SUBSONIC (BY DEFAULT FARFIELD IS DETERMINED AUTOMATICALLY WITHIN THE CODE)

SCALER= SCALE THE MESH 
        POSSIBLE VALUES= --> ANY VALUE (DIVIDES THE GRID COORDINATES BY THE SCALER VALUE)
                         

GREENGO=    GRADIENTS APPROXIMATION
        POSSIBLE VALUES= 0 --> LSQ (DEFAULT, EVERYTHING IS COMPUTED USING LSQ EXCEPT BAD QUALITY CELLS THAT USE GREEN GAUSS
                                    ONLY FOR THE APPROXIMATION OF THE GRADIENTS FOR THE DIFFUSION FLUXES)
                         1 --> GREEN GAUSS (GREEN GAUSS ONLY FOR THE APPROXIMATION OF THE GRADIENTS FOR THE DIFFUSION FLUXES)


LMACH=      LOW MACH NUMBER CORRECTION
        POSSIBLE VALUES= 0 --> NO CORRECTION
                         1 --> LMACH CORRECTION (IMPROVES MAINLY THE LOW-ORDER MUSCL AND WENO SCHEMES UP TO 3RD-ORDER, ARTIFACTS
                         MAY APPEAR WHEN ENGAGED WITH HIGHER-ORDER METHODS)

OUT_TIME=  TIME TO FINISH THE SIMULATION
        POSSIBLE VALUES= --> ANY VALUE
                         

                         
NTMAX=    MAXIMUM NUMBER OF ITERATIONS TO FINISH THE SIMULATION
        POSSIBLE VALUES= --> ANY VALUE

WALLC=   WALLCLOCK TIME LIMIT (A CHECKPOINT FILE AND OUTPUT FILE WILL BE WRITTEN WHEN THIS TIME IS MET)
        POSSIBLE VALUES= --> ANY VALUE
                         

TECPLOT=   OUTPUT FILE FORMAT 
        POSSIBLE VALUES= 1 --> TECPLOT BINARY (ONE FILE FOR THE ENTIRE DOMAIN)
                         2 --> VTK BINARY (ONE FILE FOR THE ENTIRE DOMAIN)
                         3 --> VTK BINARY PARTITIONED OUTPUT  
                         4 --> TECPLOT BINARY PARTITIONED OUTPUT


IEVERY=  HOW OFTEN (WALLCLOCK TIME IN SECONDS) TO WRITE AN OUTPUT FILE
         POSSIBLE VALUES= --> ANY VALUE


IEVERY2=HOW OFTEN (WALLCLOCK TIME IN SECONDS) TO WRITE A RESTART/CHECKPOINT FILE
         POSSIBLE VALUES= --> ANY VALUE


IEVERYAV=HOW OFTEN (WALLCLOCK TIME IN SECONDS) TO WRITE AN AVERAGED OUTPUT FILE
         POSSIBLE VALUES= --> ANY VALUE


STENCIL_IO= ENABLE WRITING OF THE OUTPUT FILES FOR THE STENCILS FOR EACH OF THE PROBE LOCATIONS
        POSSIBLE VALUES= 1 --> ENABLED
                         0 --> DISABLED

Averaging= ENABLE AVERAGING WITHIN THE CODE
        POSSIBLE VALUES= 1 --> ACTIVATED (SHOULD ONLY BE USED FOR UNSTEADY ILES, DDES,DES,URANS SIMULATIONS)
                         0 --> DEACTIVATED (DEFAULT)



OUTSURF= ENABLE WRITING SURFACE OUTPUT SOLUTION FOR WALL BOUNDARIES
        POSSIBLE VALUES= 1 --> ACTIVE
                         0 --> DEACTIVATED


IFORCE= COMPUTE FORCES (CL,CD)
        POSSIBLE VALUES= 1 --> ACTIVE
                         0 --> DEACTIVATED
                         
surfshear=  ENABLE WRITING SHEAR STRESSES ON SURFACE OUTPUT SOLUTION FOR WALL BOUNDARIES
        POSSIBLE VALUES= 1 --> ACTIVE
                         0 --> DEACTIVATED


IRES_TURB= PREVIOUS SIMULATION TYPE (RESTART)
        POSSIBLE VALUES= 0 --> WITHOUT TURBULENCE MODEL
                         1 --> WITH TURBULENCE MODEL


IRES_UNSTEADY=  PREVIOUS SIMULATION TYPE (RESTART)
        POSSIBLE VALUES= 0 --> STEADY
                         1 --> UNSTEADY


LAMPS=  PREVIOUS PASSIVE SCALAR PRESENT IN RESTART FILE 
        POSSIBLE VALUES= 0 --> NO PREVIOUS PASSIVE SCALAR
                         1 --> PREVIOUS PASSIVE SCALAR PRESENT IN RESTART FILE


Prev_turbmodel=     PREVIOUS TURBULENCE MODEL USED
        POSSIBLE VALUES= 0 --> NO PREVIOUS TURBULENCE MODEL
                         1 --> SPALART-ALLMARAS 
                         2 --> K-OMEGA

NPROBES=  NUMBER OF PROBES IN THE DOMAIN
        POSSIBLE VALUES= ANY VALUE --> ENSURE THAT YOU PROVIDE THEIR COORDINATES BELOW
                         
                         
                         
