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
#include <math.h>
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
#define REARTH 6371.
#define RCMB 2871.

#include "determine_coeff.h"
// macros
#define MEMERROR {fprintf(stderr,"memory allocation error\n");exit(-1);}



