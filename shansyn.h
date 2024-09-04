/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: spherical_harmonics_functions.h,v 1.17 2009/04/16 00:55:53 becker Exp becker $ */
/*

  function declarations for spherical harmonics package

*/
// structure for a model consisting of spherical harmonic
// expansions at different depths between dmin and dmax
// can be of type DISCRETE, CHEBYSHEV or SPLINES
//
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <zlib.h>


#ifdef USE_GMT4
#include "gmt.h"
#else
#include "gmt_dev.h"
//#include "gmt_private.h"
#endif

#include "precision.h"

struct mod{
  int n,lmax,lms;
  int vector_field,izero;
  int radial_type;			/* radial representation */
  int is_gsh;				/* generalized spherical
					   harmonics 
					   0: regular scalar
					   1: GSH format scalar
					   3: two phi (2 on input)
					   5: four phi (4 in input)
					*/
  COMP_PRECISION **a,**b,*d;
  COMP_PRECISION **amp,**amt,**bmp,**bmt;
  COMP_PRECISION dmin,dmax;
};
/*
  Chebyshev interpolation related, model output types
*/
#define DISCRETE 0
#define CHEBYSHEV 1
#define SPLINES 2
/* model output types */
#define HC_TYPE 3
//
// modes of the linear regression
//
#define REGRESS_XDETERMINED 0
#define REGRESS_ITERATIVE 1
//
// real earth constants
//
#define RADIUS_EARTH 6371.e3
#define REARTH_KM (RADIUS_EARTH/1e3)
#define RCMB_KM 2871.

#include "determine_coeff.h"
// macros
#define MEMERROR {fprintf(stderr,"memory allocation error\n");exit(-1);}

#ifndef SQUARE	
#define SQUARE(x) ((x)*(x))
#endif

#define SIGNIFICANTLY_DIFFERENT(x, y) ((fabs((x)-(y))>1.0e-04)?(1):(0))

#ifdef BOOLEAN
// the upper case BOOLEAN gets defined by old GMT as int, so 
typedef int boolean;
#else
typedef unsigned short boolean;
#define BOOLEAN unsigned short
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef MAX
#define MAX(x, y) (((x)<(y))?(y):(x))
#endif 
#ifndef MIN
#define MIN(x, y) (((x)<(y))?(x):(y))
#endif 


#ifndef SIGN
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif
#ifndef FMAX
#define FMAX MAX
#endif

#define ONE_MEGABYTE 1048576.0

#ifndef POSLM
 #define POSLM(l, m) ( (((l) * (1+(l))) / 2) + (m))
#endif
#ifndef LPOSLM
 #define LPOSLM(l, m, j) (((j)*(lmsize)+POSLM(l,m)))
#endif
#ifndef P
 #define P(l, m, j)        (p[LPOSLM(l, m, j)])
#endif
#ifndef DPTHETA
 #define DPTHETA(l, m, j)  (dptheta[LPOSLM(l, m, j)])
#endif
#ifndef DP2THETA
#define DP2THETA(l, m, j) (dp2theta[LPOSLM(l, m, j)])
#endif

#ifndef PI
#define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#endif
#ifndef TWOPI
#define TWOPI 6.2831853071795864769252867665590057683943387987502116419498891846
#endif
#ifndef TWO_PI
#define TWO_PI TWOPI
#endif
#ifndef PIOVERONEEIGHTY
#define PIOVERONEEIGHTY 0.0174532925199432957692369076848861271344287188854172545609719144
#endif
#ifndef PIHALF
#define PIHALF 1.5707963267948966192313216916397514420985846996875529104874722961
#endif
#ifndef SQRT_TWO
#define SQRT_TWO 1.4142135623730950488016887242096980785696718753769480731766797379
#endif
#ifndef  TWO_PISQR
#define TWO_PISQR 19.739208802178717237668981999752302270627398814481581252826698752
#endif
#ifndef ONEEIGHTYOVERPI 
#define ONEEIGHTYOVERPI  57.295779513082320876798154814105170332405472466564321549160243861
#endif
#ifndef HALF_SQRT_TWO 
#define HALF_SQRT_TWO    0.7071067811865475244008443621048490392848359376884740365883398689
#endif
#ifndef TWO_SQRT_PI
#define TWO_SQRT_PI 3.5449077018110320545963349666822903655950989122447742564276155797
#endif
#ifndef ONE_OVER_SQRT_FOUR_PI
#define ONE_OVER_SQRT_FOUR_PI 0.2820947917738781434740397257803862929220253146644994284220428608
#endif
#ifndef SQRT_TWO_PI
#define SQRT_TWO_PI 2.5066282746310005024157652848110452530069867406099383166299235763
#endif
#ifndef SQRT_THREE_FOURTHS
#define SQRT_THREE_FOURTHS 0.8660254037844386467637231707529361834714026269051903140279034897
#endif
#ifndef SQRT_THREE
#define SQRT_THREE 1.7320508075688772935274463415058723669428052538103806280558069794
#endif

#define THETA2LATITUDE(x) ( (90.0 - (x)*ONEEIGHTYOVERPI) )
#define LATITUDE2THETA(x) ( (90.0 - (x))*PIOVERONEEIGHTY )
#define DEG2RAD(x) ( (x)*PIOVERONEEIGHTY )
#define PHI2LONGITUDE(x) ( (x) * ONEEIGHTYOVERPI )
#define LONGITUDE2PHI(x) ( (x) * PIOVERONEEIGHTY )
#define RAD2DEG(x) ( (x) * ONEEIGHTYOVERPI)

#include "auto_proto.h"


