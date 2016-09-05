#ifndef SHANSYN_READ_PRECISION
/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: precision.h,v 1.9 2009/04/16 00:55:53 becker Exp becker $ */

// this shouldn't be needed since we are including limits.h
// but somehow we needed to fix that for LINUX
#ifndef FLT_MAX
 #define FLT_MAX 3.40282347E+38F
#endif
#ifndef FLT_MIN
 #define FLT_MIN 1.17549435E-38F
#endif


/*

  data precision for input is always single precision
  since we want to interact with GMT

*/
#define DATA_PRECISION float
#define GMT_PRECISION float
#define THREE_GMTDATA_FSCAN_FORMAT "%f %f %f"

/* 
   precision for internal calculations 

   it turns out that FFT and sum methods return 
   EXACTLY the same results for single and 
   COMP_PRECISION precision settings lm 10 mode deviations 
   from the analytical function are 2.8e-07 and 5.6e-16 
   per node respectively

*/
/*
  
  default is double precision calculation, if SINGLE_PREC
  is not defined

  run 

  cproto -f2 *.c | grep -v main > auto_proto.h

  if precision is changed


*/
#ifndef SINGLE_PREC

 #define COMPILE_PREC "double"
 #define COMP_PRECISION double
 #define TWO_TWO_DATA_FSCAN_FORMAT "%*g %*g %lf %lf"
 #define FOUR_DATA_FSCAN_FORMAT "%lf %lf %lf %lf"
 #define TWO_DATA_FSCAN_FORMAT "%lf %lf"
 #define DATA_FSCAN_FORMAT "%lf"
 #define ASCII_DATA_FORMAT "%22.15e " 
 #define EPS_COMP_PREC 5.0e-15
 #define EPS_DATA_PREC  5.0e-7
 #define COEFF_DATA_FORMAT "%22.15e %22.15e"
#else

 #define COMPILE_PREC "float"
 #define COMP_PRECISION float
 #define TWO_TWO_DATA_FSCAN_FORMAT "%*g %*g %f %f"
 #define FOUR_DATA_FSCAN_FORMAT "%f %f %f %f"
 #define TWO_DATA_FSCAN_FORMAT "%f %f"
 #define DATA_FSCAN_FORMAT "%f" 
 #define ASCII_DATA_FORMAT "%15.7e "
 #define EPS_COMP_PREC 5.0e-7
 #define EPS_DATA_PREC  5.0e-7
 #define COEFF_DATA_FORMAT "%15.7e %15.7e"
#endif

#define SHANSYN_READ_PRECISION
#endif
