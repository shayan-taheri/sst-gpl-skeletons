#ifndef _BENCH_GTC_OPT
#define _BENCH_GTC_OPT

/* measure the time in kernels */
#define FINE_TIMER           0

/* put a barrier in timer */
#define BARRIER_ON           1

/* Set to 1 when ntoroidal = mzetamax in experiments */
#define ASSUME_MZETA_EQUALS1 1

/* Use SSE intrinsics for some steps in charge/push */
/* Set this to 1 only when mzeta=1 and on x86 platforms */
#define SIMD_CODE            0

/* Avoid lookups to the pgyro/gyro arrays in address calculation step.
 * This reduces the memory footprint while increasing the floating 
 * point computation. Should be a good strategy on GPUs. */
#define GYRO_LOCAL_COMPUTE   0

/* This setting avoids loads from auxiliary arrays in push,
   and stores in charge */
#define ONTHEFLY_PUSHAUX     0

/* Useful when npartdom > 1 */
#define UNIFORM_PARTICLE_LOADING 0

/* Charge density grid isn't updated if this is set to 0.
   Used to quantify performance impact of irregular accesses 
   to the charge density grid. */
#define DO_CHARGE_UPDATES    1

#define IDEAL_ALIGNMENT      16

/* Do not modify below this point */
/*****************************************************/

/* Currently unsupported options, carried 
   over from pthreads code  */

#ifndef SINGLE_PRECISION
#define SINGLE_PRECISION     0
#else
#define SINGLE_PRECISION     1
#endif

#ifndef MIXED_PRECISION  
#define MIXED_PRECISION 0
#else
#define MIXED_PRECISION 1
#endif

#define RADIAL_FINEPART      0

#define PRINT_CHECKSUM       0

/* Permute radial coordinate of particles. 
 * This is set to 0 by default.*/
#define PERMUTE_ZION_R       0

#define SQRT_PRECOMPUTED     0

#define LOCKING              1

#if defined(__x86_64__)
#if (SINGLE_PRECISION | MIXED_PRECISION)
#define SIMD_INCR_DOUBLEPAIR 0
#else
#define SIMD_INCR_DOUBLEPAIR 1
#endif
#else
#define SIMD_INCR_DOUBLEPAIR 0
#endif

#define NUM_TRIALS 10
#define RNG_SEED 232323

#endif
