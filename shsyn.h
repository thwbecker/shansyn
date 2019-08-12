/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: shsyn.h,v 1.5 2009/04/16 00:55:53 becker Exp becker $ */
#include "shansyn.h"

#ifndef BE_VERBOSE
 #define BE_VERBOSE FALSE
#endif
#ifndef WPF
 #define WPF TRUE
#endif 

#define MY_FLT_MAX FLT_MAX
#define MY_FLT_MIN FLT_MIN

#define EPS_SMALL 1.0e-16
#define FOURIER_FACTOR HALF_SQRT_TWO 


/* names of binary files for velocity output */

#define FIRST_VEL_OUT  "vec_p.grd" /* is in x or phi direction */
#define SECOND_VEL_OUT "vec_t.grd" /* is in -y or theta direction  */

/* io matters */
#define IN_FILE stdin
#define OUT_FILE stdout

/* TP and UP coefficients assignment rules */
#define TP(l, m, j) (*(r+ (j)* (rslmsize) +  ((l)*(l+1)) + 2*(m) + 0 ))
#define UP(l, m, j) (*(r+ (j)* (rslmsize) +  ((l)*(l+1)) + 2*(m) + 1 ))

#define AP(l, m)    ( (*(ap+POSLM(l, m))) )
#define BP(l, m)    ( (*(bp+POSLM(l, m))) )
#define AT(l, m)    ( (*(ap+lmsize+POSLM(l, m))) )
#define BT(l, m)    ( (*(bp+lmsize+POSLM(l, m))) )

/* sum modes */
#define SUM_OVER_L_AND_M 0
#define FFT_SUM 1
/* output modes */
#define ASCII_STDOUT 0
#define ONE_GRD 1
#define BINARY_STDOUT 2
#define TWO_GRDS 3
#define GRADIENT 4
/* default methods */
#define DEF_SUM_MODE SUM_OVER_L_AND_M
#define DEF_OUT_MODE ASCII_STDOUT
/* default values for expansion  */
#define XMINS 0.0
#define XMAXS 360.0
#define YMINS -90.0
#define YMAXS 90.0
#define DEF_INC 1.0
/* 
   memory limit for the P array 
   after which we will evaluate 
   Legendre functions each time....
   in MB
*/
#define P_MEM_LIMIT 4200.0




