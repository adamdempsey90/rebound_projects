c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c      USER_PARAMETERS.INC
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c Author: Matthew John Payne
c
c Define various additional quantities / parameters / switches etc
c
c These are used to define and control the various user-defined forces 
c 
c ------------------------------------------------------------------------------
      integer NPARAMS
c 
      integer USEROUTPUT
c 
      integer DRIFTSWITCH
      real*8  DRIFTTOL
c 
      integer MIGRATIONSWITCH
      real*8  MIGRATIONSTART, MIGRATIONEND, MIGRATIONRATIO
      character*20 MIGRATIONMECHANISM
c 
      integer TIMEDEPSWITCH,GASDRAGSWITCH
      real*8  GASDRAGSTART,GASDRAGEND,GASDRAGCUTOFF
      real*8  TAU_DISC,DISSIPATION_OFFSET
      real*8  ALPHA_DISC,EPSILON_DISC,GAMMA_DISC,BETA_DISC
      real*8  H_DISC,H2_DISC,Z00_DISC,Z002_DISC
      real*8  RHO_GCM3,RHO_CGS,RHO_DISC
      real*8  ONE28_OVER_9PI,C_S0,SIXSQRTPIOVER8,MU_OVER_SIGMA
      integer NODE_PRECESS_MODEL,PERIASTRON_PRECESS_MODEL
      integer GAS_ECCENTRICITY_MODEL
      real*8  A_BIN,E_BIN,MAX_FACTOR,A_DISK_MAX
      real*8  I_BIN_GAS,MIN_INC,COS_IBG,SIN_IBG
      real*8  MASS_PRIM,MASS_SEC,MASS_RATIO
      real*8  LARWOOD1,LARWOOD2,LARWOOD3,LARWOOD,LAR_SCALE
      real*8  PRECESS_TIMESCALE
      real*8  FORCED_ECC_SCALING,DISK_FORCED_ECC,DISK_FORCED_ECC2
      real*8  dN_INITIAL,dW_INITIAL,dN_ARBITRARY,dW_ARBITRARY
c 
      integer MASSGROWTHSWITCH
      real*8  MASSGROWTHSTART, MASSGROWTHEND
      character*20 MASSGROWTHMODEL
c 
      integer HYBRIDTOBSSWITCH
      character*20 HYBRIDTOBSCRITERION
      real*8  HYBRIDTOBSDIST
      integer BSTOHYBRIDSWITCH
      character*20 BSTOHYBRIDCRITERION
      real*8  BSTOHYBRIDTIME
c
      integer RELATIVITYSWITCH
      real*8  c2
c
c
      integer RECORDCLOSEENCOUNTERS
      integer CENTRALSTARCOLLISIONS
      integer BIGBIGCOLLISIONSONLY
c
      integer TIDALSWITCH
      integer TIDALBINARYSWITCH
      real*8  TIDALCUT
      real*8  TIDALSTARLOVE
      real*8  TIDALSTARQ
      real*8  TIDALSTARMOI
      real*8  TIDALSTARSPINX,TIDALSTARSPINY,TIDALSTARSPINZ
      real*8  TIDALSTARSPINMAG
c
      integer DISKPOTENTIALSWITCH
c
      integer T4NAGASAWATIDALSWITCH
      integer T4NAGASAWABINARYSWITCH
c
      integer DISKGRAVITYSWITCH
c
c
c
c
c

c------------------------------------------------------------------------------
c How many additional parameters can be associated with each individual body
      parameter (NPARAMS = 18)
c
c------------------------------------------------------------------------------
c Do you wish to switch on the limited output to user.out? (1 = ON, 0 = OFF)
      parameter (USEROUTPUT = 1)
c
c ------------------------------------------------------------------------------
c Set the parameters for the drift-checking
c DRIFTSWITCH = 1/0 ==>> Turns it ON/OFF
      parameter (DRIFTSWITCH = 0)
      parameter (DRIFTTOL = 1.001)
c
c ------------------------------------------------------------------------------
c Set migration switches/parameters
c MIGRATIONSWITCH ( 1 = ON,   0 = OFF)
c MIGRATIONSTART & MIGRATIONEND in years
c Note that the individual migration rates are specified in big.in, etc
c Chose different migration mechanisms... 
c ...'LPSIMPLE_EXP'
      parameter (MIGRATIONSWITCH    = 1)
      parameter (MIGRATIONSTART     = 0)
      parameter (MIGRATIONEND       = 1458000)
      parameter (MIGRATIONMECHANISM = 'LP03')
      parameter (MIGRATIONRATIO = 0)
c
c ------------------------------------------------------------------------------
c Set Gas Drag Switches/parameters
c 
c Timing___________
      parameter (TIMEDEPSWITCH         = 0)                       ! ( 0 = NO GAS DISS,   1 = GAS CAN DISS, SO HAVE TIME DEP)
      parameter (GASDRAGSWITCH         = 0)                       ! ( 0 = OFF, 1 = Simple (circular) gas drag, 2 = Eccentric Gas Drag Model)
      parameter (GASDRAGSTART          = 0)			  ! When to start gas drag during the simulation
      parameter (GASDRAGEND            = 1000000000)		  ! When to stop gas drag during the simulation
      parameter (GASDRAGCUTOFF         = 1e-07)			  ! Don't simulate gas drag for very tenuous disks
      parameter (DISSIPATION_OFFSET    = 0)			  ! Can be +ve or -ve
      				       	 			  ! +ve implies "time disc is ON FULL" (i.e. does not dissipate)
								  ! -ve implies "disc has started to dissipate" (i.e. systems have evolved)
      parameter (TAU_DISC              = 1e6)			  ! Timescale for exponential disc dissipation (years)
c Disk_____________
      parameter (ALPHA_DISC            = 1.5) 			  ! Disc surface density, Sigma ~ r^-ALPHA = r^-(EPSILON-GAMMA)
      parameter (GAMMA_DISC            = 1.25)			  ! Z00 propto r^GAMMA_DISC, where Z0 is the exponential decay const in exp{z^2 / Z0^2}
      parameter (EPSILON_DISC          = ALPHA_DISC + GAMMA_DISC) ! Disc midplane density,      rho_0 ~ r^-EPSILON_DISC
      parameter (BETA_DISC             = 0.5)			  ! Disc temperature structure, T     ~ r^-BETA_DISC
      parameter (H_DISC                = 0.05)		  	  ! Disc vertical scale height, H = 0.05
      parameter (H2_DISC               = H_DISC*H_DISC)		  !                             H^2
      parameter (Z00_DISC              = sqrt(2.0)*H_DISC) 	  ! Alternative Disc vertical scale height, Z00 = sqrt{2}*H ~ 0.07
      parameter (Z002_DISC             = Z00_DISC*Z00_DISC) 	  !                                         Z00^2
      parameter (RHO_CGS               = 0.49806095E+03)	  ! RHO_CGS = AU*AU*AU*K2 / MSUN (=0.49806095E+03)
      parameter (RHO_GCM3              = 1.4e-09)		  ! RHO_GCM3 = Midplane density of disc in g cm^-3
      parameter (RHO_DISC              = RHO_GCM3 * RHO_CGS) 	  ! RHO_DISC = Midplane density of disc in correct units for Mercury
      				       	 			  ! 1.4*10^-9 g cm^-3 * rhocgs  = 6.9728533e-07
								  ! N.B. 1.4*10^-9 g cm^-3 should correspond to Sigma~1,800 g cm^-2 for h/a=0.05
c Constants________
      parameter (C_S0		         = 1.0e5)		  	  ! C_S0 = ~1km	= 1e5 cm	(see Lyra '09, Rafokov '04)
      parameter (MU_OVER_SIGMA	       = 3.87e-09)		  ! MU_OVER_SIGMA = mean molecular mass / cross section = 2.3 * mH / 10^-15 cm^2 
      				       	 			  ! (see Hayashi '81, Rafikov 04)
      parameter (ONE28_OVER_9PI        = 4.52707394)		  ! ONE28_OVER_9PI	   = 128 / (9*PI)
      parameter (SIXSQRTPIOVER8        = 3.75994241)		  ! 6sqrt{pi/8}
c Binary___________      
      parameter (NODE_PRECESS_MODEL		= 0)		  ! 0=arbitrary; 1=Larwood    
      parameter (PERIASTRON_PRECESS_MODEL    	= 1)		  ! 0=arbitrary; 1=minmod; 2=superbee
      parameter (GAS_ECCENTRICITY_MODEL    	= -2)		  ! <0= simplest model(-1=minmod; -2=superbee): 0, ed=constant
      						  		  ! N.B. minmod and superbee are tested only in Gamma Cephei by Paardekooper et al (2008)
      parameter (A_BIN    	       = 0.5)		  	  ! Semi-Major Axis of the Binary 
      parameter (MAX_FACTOR    	       = 0.2)		  	  ! Max Semi-Major Axis of the Disk in a Binary System as a fraction of the Binary Semi-Maj Axis
      parameter (A_DISK_MAX    	       = MAX_FACTOR*A_BIN) 	  ! Max Semi-Major Axis of the Disk in a Binary System
      parameter (E_BIN    	       = 0.01)		  	  ! Eccentricity of the Binary
      parameter (I_BIN_GAS    	       = 0*0.0174532925)	  ! Inclination of the Gas Disk relative to the Binary Orbit (Radians [0.0174532925Rad = 1deg])
      				       	 			  ! There is an implicit assumption that the binary must be in the i=0 plane (to get other bodies correct)
      parameter (MIN_INC	       = 1e-05)		  	  ! Minimum inclination of I_BIN_GAS necessary to bother with calculations
      parameter (COS_IBG    	       = dcos(I_BIN_GAS))	  !           cosine of I_BIN_GAS
      parameter (SIN_IBG    	       = dsin(I_BIN_GAS))	  !           sine of   I_BIN_GAS
      parameter (MASS_PRIM    	       = 1)		  	  ! Mass of the Primary Star		
      parameter (MASS_SEC    	       = 9.09090909090909e-05)		  	  ! Mass of the Secondary Star	
      parameter (MASS_RATIO    	       = MASS_SEC/MASS_PRIM)	  ! Mass Ratio: Secondary / Primary
      parameter (LAR_SCALE 	       = 1)			  ! Allow a scaling of the analytic Larwood presession rate.
      parameter (LARWOOD1    	       = LAR_SCALE*(15.0/32.0)*MASS_RATIO)  ! Precession rate (Larwood et al. 1996)
      parameter (LARWOOD2    	       = LARWOOD1*MAX_FACTOR**3)  !                  ||
      parameter (LARWOOD3 	       = LARWOOD2*COS_IBG)  	  !                  ||
      parameter (LARWOOD 	       = LARWOOD3*sqrt(MASS_PRIM))!		     || 
      parameter (PRECESS_TIMESCALE     = 2000)	   		  ! Precession timescale (years) for the disk in the XXX model from Paardekooper)
      parameter (FORCED_ECC_SCALING    = 1)	   		  ! GAS_ECCENTRICITY_MODEL = +1, need to scale eccentricity to simple forced eccentricity from binary
      parameter (DISK_FORCED_ECC=FORCED_ECC_SCALING*1.25*A_BIN/E_BIN) 
      parameter (DISK_FORCED_ECC2      = DISK_FORCED_ECC**2)
      parameter (dN_INITIAL            = 0)	  		  ! Initial precession for the node in NODE_PRECESS_MODEL = 1
      parameter (dN_ARBITRARY          = 0)		  	  ! Allow arbitrary scaling of the precession for the periastron in NODE_PRECESS_MODEL = 0 & 1
      parameter (dW_INITIAL            = 0)		  	  ! Initial precession for the periastron in  PERIASTRON_PRECESS_MODEL = 0
      parameter (dW_ARBITRARY          = 0)		  	  ! Allow arbitrary scaling of the precession for the periastron in PERIASTRON_PRECESS_MODEL = 0 & 2
c
c
c ------------------------------------------------------------------------------
c Set mass growth switches/parameters
c MASSGROWTHSWITCH ( 1 = ON,   0 = OFF)
c MASSGROWTHSTART & MASSGROWTHEND in years
c Note that the individual migration rates are specified in big.in, etc
c Specify growth curve: 'EXPONENTIAL' or 'LINEAR'
      parameter (MASSGROWTHSWITCH    = 0)
      parameter (MASSGROWTHSTART     = 1)
      parameter (MASSGROWTHEND       = 5000000)
      parameter (MASSGROWTHMODEL     = 'EXPONENTIAL')
c
c ------------------------------------------------------------------------------
c Set integration method switches/parameters
      parameter (HYBRIDTOBSSWITCH    = 0)
      parameter (HYBRIDTOBSCRITERION = 'DISTANCE')
      parameter (HYBRIDTOBSDIST      = 0.1)
c
      parameter (BSTOHYBRIDSWITCH    = 0)
      parameter (BSTOHYBRIDCRITERION = 'TIME')
      parameter (BSTOHYBRIDTIME      = 1000000)
c
c ------------------------------------------------------------------------------
c Set Relativity Switch
      parameter (RELATIVITYSWITCH    = 0)
      parameter (c2                  = 2.98185398e4)
c
c ------------------------------------------------------------------------------
c Set Close Encounter Switch
      parameter (RECORDCLOSEENCOUNTERS    = 0)
      parameter (CENTRALSTARCOLLISIONS    = 0)
      parameter (BIGBIGCOLLISIONSONLY     = 0)
c
c
c ------------------------------------------------------------------------------
c Set Tidal Switches/parameters
      parameter (TIDALSWITCH = 0)
      parameter (TIDALBINARYSWITCH = 0)
      parameter (TIDALCUT = 0.5)
      parameter (TIDALSTARLOVE = 0.028)
      parameter (TIDALSTARQ = 1e6)
      parameter (TIDALSTARMOI = 0.076)
      parameter (TIDALSTARSPINX = 0)
      parameter (TIDALSTARSPINY = 0)
      parameter (TIDALSTARSPINZ = 1)
      parameter (TIDALSTARSPINMAG = 2e-6)
c
c ------------------------------------------------------------------------------
c Set Nagasawa Tidal Switches/parameters
      parameter (T4NAGASAWATIDALSWITCH  = 0)
      parameter (T4NAGASAWABINARYSWITCH = 0)
c
c 
c 
c
c ------------------------------------------------------------------------------
c Set Parameters related to the non-central disk potential
      parameter (DISKGRAVITYSWITCH = 0)
c
