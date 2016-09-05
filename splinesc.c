/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: splinesc.c,v 1.3 2001/03/11 18:20:04 becker Exp becker $ */
#include "shansyn.h"

/*
  
  spline interpolation routine, calls FORTRAN code 
  due to Ritsema as distributed with the s20rts model

  spline knots are equally spaced from -1 .. 1 in m steps 
  when first called

  x1, x2:     left and right boundaries of integral, don't change
              once initialized
  c[0...m-1]: array with weights for splines 0 ... m-1
  m:          number of splines
  x:          location for interpolation

  $Id: splinesc.c,v 1.3 2001/03/11 18:20:04 becker Exp becker $
  
*/

#define MAX_SPLINE_N 50
COMP_PRECISION spline_base(COMP_PRECISION x1,COMP_PRECISION x2,
			   COMP_PRECISION *c,int m,
			   COMP_PRECISION x)
{
  static int init=0;
  int i;
  // static arrays for spline setup
  static COMP_PRECISION 
    qq0[MAX_SPLINE_N*MAX_SPLINE_N],
    qq[3*MAX_SPLINE_N*MAX_SPLINE_N],
    spknt[MAX_SPLINE_N],
    qqwk[3*MAX_SPLINE_N],dx;
  COMP_PRECISION xp,val;
  // check range 
  if ((x-x1)*(x-x2) > 0.0){
    fprintf(stderr,"spline_base: x (%g) not in range (%g - %g)\n",
	    x,x1,x2);
    exit(-1);
  }
  if(!init){
    if(m > MAX_SPLINE_N){
      fprintf(stderr,"spline_base: works only for order <= %i (m: %i)\n",
	      MAX_SPLINE_N,m);
      exit(-1);
    }
    splhsetup(spknt,qq0,qq,qqwk,&m);
    dx=(x2 - x1)/2.0;
    init=1;
  }
  // local coordinate from -1 to 1
  xp=-1.0+(x-x1)/dx;
  val=0.0;
  for(i=0;i<m;i++)// add up contribution from all base functions
    val += c[i]*splh(&m,spknt,qq0,qq,&i,&xp);
  return(val);

}
#undef MAX_SPLINE_N
// norm damping for splines
COMP_PRECISION spline_norm_damping(int i,int j,int n)
{
  if(i!=j)
    return 0.0;
  else{
    // all damped the same way
    return 1.0;
  }
}


