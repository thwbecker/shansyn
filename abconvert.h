/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: abconvert.h,v 1.16 2009/04/16 00:55:53 becker Exp becker $ */
/*
  header file for abconvert.c. conversion of
  spherical harmonics expansion

  $Id: abconvert.h,v 1.16 2009/04/16 00:55:53 becker Exp becker $

*/

#include "shansyn.h"
#ifndef BE_VERBOSE
 #define VERBOSE 0
#else
 #define VERBOSE 1
#endif

/* conversion factors from physical to geodetic 
   normalization */

#define GEODETIC_FACTOR(l, m) ((TWO_SQRT_PI)*pow(-1.0,(COMP_PRECISION)(-m)))
#define RICK_PRE_F (4.0*SQUARE(PI))
#define RICK_FACTOR(l, m) (RICK_PRE_F/pow(-1.0,(COMP_PRECISION)(m)))
#define RICK_SCALAR_FACTOR(l, m) (pow(-1.0,(COMP_PRECISION)(-m))*2.0*SQRT_TWO_PI)

/* 
   output switches 
*/
#define ABPHYS_OUT 0
#define ABGEOD_OUT 1
#define POWER_OUT 2
#define TAPER_OUT 3
#define LMAB_OUT 4
#define CORRL_OUT 5
#define CORRT_OUT 6
#define VECABAB_OUT 7
#define VECAABBR_OUT 8
#define LMAB_GEODETIC_OUT 9
#define TRMS_OUT 10
#define TPOWER_OUT 11
#define VECABAB_POL_OUT 12
#define VECABAB_TOR_OUT 13
#define VECABAB_NEW_OUT 14
#define GSH_OUT 15
#define MEAN_OUT 16
#define CORRTH_OUT 17
#define ABPHYS_NONZERO_OUT 18
#define CCL_COUPLING 19
#define ABPHYS_OUT_LONG 20
#define VECABAB_OUT_LONG 21
#define SPEAR_CORRL_OUT 22
#define SPEAR_CORRT_OUT 23
#define SPEAR_CORRTH_OUT 24
#define ABPHYS_NONZERO_OUT_ONE_COLUMN 25
#define POWER_OUT_NN 26
#define ADMITTANCE_OUT 27
#define ADMITTANCE_OUT_NN 28

/* 
   taper switches 
*/
#define NO_TAPER 0
#define COSSQR_TAPER 1
#define ONEML_TAPER 2
#define ONEMLSQR_TAPER 3
#define SINOVRL_TAPER 4
#define SINOVRLM_TAPER 5
#define NNR_TAPER 6
#define ZERO_TAPER 7
#define L0_TAPER 8
#define FROM_FILE_TAPER 9
#define PHI_ROTATE 10
#define COSP4_TAPER 11
#define NR_TAPER 12
#define FROM_SH_FILE_TAPER 13
#define PASS_POL_TAPER 14
#define PASS_TOR_TAPER 15
#define SET_L_UNITY 16
#define LAPLACIAN 17
// file to read filter from 
#define FILTER_FILE "filter.dat"
#define FILTER_SH_FILE "filter.ab"
#define ONLY_L_ZERO 18
#define EXP_TAPER 19
#define HIGHP_TAPER 20
#define ONLY_M0_TERMS 21
/* 
   input switches 
*/
#define AB_INPUT 0
#define AB_INPUT_HC 14
#define LAB_GEOD_INPUT 1
#define AB_GEOD_INPUT 2
#define AB_RICK_INPUT 3
#define ABAB_INPUT 4
#define ABAB_GEOD_INPUT 5
#define AABBR_INPUT 6
#define LMAB_GEOD_INPUT 7
#define INTERPOLATE 8
#define INTERPOLATE_ABAB 9
#define LMAB_FNORM_INPUT 10
#define AB_MASTERS_INPUT 11
#define GSH_INPUT 12
#define AB_NONZERO_INPUT 13
#define INTERPOLATE_GSH 15


/* defaults */
#define IN_FORMAT_DEFAULT AB_INPUT
#define OUT_FORMAT_DEFAULT ABPHYS_OUT
#define TAPERING_DEFAULT NO_TAPER 
#define AMPLITUDE_DEFAULT 1.0
#define LC_DEFAULT 0.0

