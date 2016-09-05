/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: gauleg.c,v 1.6 2001/03/08 03:04:49 becker Exp $ */
#include <math.h>

/* 
   find weights and integration points for 
   Gauss-Legendre integration of order n
   from Numerical Recipes.

   C-type x[0...n-1] arrays have to call this 
   function like ...,x-1,...

   $Id: gauleg.c,v 1.6 2001/03/08 03:04:49 becker Exp $
   
*/
#include "precision.h"
#include "spherical_harmonics_functions.h"
// relative precision for weight locations
#ifdef DOUBLE_PREC
#define EPS_GAUL 3.0e-15
#else
#define EPS_GAUL 3.0e-11
#endif

void gauleg_orig(COMP_PRECISION x1,COMP_PRECISION x2,
		 COMP_PRECISION *x,COMP_PRECISION *w,
		 int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  
  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++) {
    z=cos(3.1415926535897932384626434*(i-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS_GAUL);
    x[i]=xm-xl*z;
    x[n+1-i]=xm+xl*z;
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i]=w[i];
  }
}
#undef EPS_GAUL

