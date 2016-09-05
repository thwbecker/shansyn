/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: shana.h,v 1.16 2009/04/16 00:55:53 becker Exp becker $ */
//
// header file for shana.c, expansion routine
// 
// $Id: shana.h,v 1.16 2009/04/16 00:55:53 becker Exp becker $
//
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include "shansyn.h"

#ifndef BE_VERBOSE
 #define BE_VERBOSE 0
#else
 #define BE_BERBOSE 1
#endif
/* reduces lmax automatically for trapezoidal integration
   to suit lmax+1 <= nlat */
//#define REDUCE_LMAX

/* #define CHECK_FLOATING_EXCEPTION */
#define EPS_SMALL EPS_COMP_PREC
/* filenames */
#define FST_VEC_FILE_START "vec_p" /* name of first vector component
				      file */
#define SCD_VEC_FILE "vec_t.grd"

/* limit for P array size 
   if P array gets larger, the routine switches to single 
   Legendre function evaluation
*/
#define MEMLIMIT 1000

#ifndef GMT_PRECISION
 #define GMT_PRECISION DATA_PRECISION
#endif

#define STRING_LENGTH 100

/* integration modes */
#define TRAPEZOIDAL 0
#define GAUSSIAN 2
#define LEAST_SQUARES 3

/* operational modes */
#define BLOCK 0
#define CONTOUR 1
#define ASCII_BLOCK 2
#define ASCII_BLOCK_HEADER 3
#define GMT_BLOCK 4

#define POINTS_FOR_A_MATRIX 5
#define POINTS_FOR_A_MATRIX_ASCII 7
#define POINTS_FOR_AVEC_MATRIX 50
#define POINTS_FOR_AVEC_MATRIX_ASCII 70

#define POINTS_FOR_AB 6




// default settings
#define DEF_INT_MODE TRAPEZOIDAL
#define DEF_LMAX 31
#define DEF_OPMODE 2
#define DEF_DAMPING 0.0

#ifdef USE_LAPACK
/* LAPACK */
extern void sgels_( char *, int *, int *,int *, float *, int *, float *, int *, float *, int *, int *);
void solver_ab_lls(float *, int , int ,float *,char *);
#endif
