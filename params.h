// v.1.2.1

//---- exactly one method should be defined------
//#define GILLESPIE
//#define FASTER_KMC
#define NORMAL  // standard simulation method described in the paper (non-KMC)
//-----------------------------------------------

//#define PUSHING // if defined, cells can push away other cells as they grow

//#define CONST_BIRTH_RATE // if defined, birth rate does not depend on the no. of empty neighbours as long as there is at least one


//#define MAKE_TREATMENT_N // if defined, simulate treatment after reaching given size
//#define MAKE_TREATMENT_T // if defined, simulate treatment after reaching given time

const float gama=1e-2, gama_res=5e-8 ; // these are rates per daughter cell. Rates per diploid exome will be 2x higher (these values are given in the paper)

//#define MIGRATION_MATRIX

// used when MIGRATION_MATRIX is not defined
const float migr=0 ;
// used only when MIGRATION_MATRIX is defined
//const float migr[2][2]={{0,0} // before treatment: WT/resistant
//                        ,{0,1e-5}}; // after treatment: WT/resistant

//#define DEATH_ON_SURFACE ; // if defined then cells die on surface only upon treatment, if not then cells die also in the volume
//#define CORE_IS_DEAD  // when set, core cells are set to dead 
#define SHOW_ONLY_DRIVERS
//#define NO_MECHANICS  // if defined, no mechanics is simulated (this speeds up everything but looses info on spatial distribution)
const float timescale=1./log(2.) ; // calculates the timescale factor from the division time [days] (first number, here 1.)
const float death0=0.95, growth0=1.0 ;   // before treatment

// if death on surface:
//float death1=/*0.1*/0.99, growth1=0.0 ;   // after treatment
// if death in volume
const float death1=1., growth1=0.5 ; // after treatment

const float driver_adv=0.05 ; 
const float driver_prob=2e-5 ; // driver probability per haploid genome (should be 2e-5)
const float driver_balance=1 ; // 1==drivers decrease death, 0==drivers increase growth, intermediate values are also possible
const float driver_migr_adv=0.0,  max_migr=1e-5 ; // maximum migration prob. is taken into account only when driver_migr_adv>0
const int driver_mode = 0 ; // 0== drivers affect bd only, 1==drivers affect simultaneously bd and migr, 2==drivers affect bd xor migr with equal prob.

const float cutoff=0.1 ;
#ifdef __MAIN
int max_size=int(1e6) ; // this is ignored when MAKE_TREATMENT_T is defined
float time_to_treat=10*(365./12) ; // this is ignored when MAKE_TREATMENT_N is defined

#endif

const int _resol=1 ; // spatial resolution of sampling [cells]
const int _bins=10000 ; // max number of bins