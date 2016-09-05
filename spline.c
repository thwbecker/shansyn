/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: spline.c,v 1.5 2001/03/08 03:04:49 becker Exp $ */
#include "shana.h"
/*
  
  spline interpolation routine from Numerical Recipes

  used for theta-interpolation in shana.c


  as opposed to spline base functions used by 
  determine_coefficients, these are in splinesc.c


  $Id: spline.c,v 1.5 2001/03/08 03:04:49 becker Exp $


 */


#define NR_END 1
#define FREE_ARG char*

void spline(COMP_PRECISION *x,COMP_PRECISION *y,
	    int n,COMP_PRECISION yp1 ,COMP_PRECISION ypn,
	    COMP_PRECISION *y2)
{
  int i,k;
  COMP_PRECISION p,qn,sig,un,*u;
  
  u=vector(1,n-1);
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30){
    qn=un=0.0;
  }else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free_vector(u,1,n-1);
} 

void splint(COMP_PRECISION *xa,COMP_PRECISION *ya,
	    COMP_PRECISION *y2a,int n,COMP_PRECISION x,
	    DATA_PRECISION *y)
{
	int klo,khi,k;
	COMP_PRECISION h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
	  fprintf(stderr,"Bad xa input to routine splint\n");
	  exit(-1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

COMP_PRECISION *vector(long  nl,long nh)
/* allocate a COMP_PRECISION vector with subscript range v[nl..nh] */
{
	COMP_PRECISION *v;

	v=(COMP_PRECISION *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
	if (!v) {
	  fprintf(stderr,"memerror\n");
	  exit(-1);
	}
	return v-nl+NR_END;
}


void free_vector(COMP_PRECISION *v,long nl,long nh)
/* free a COMP_PRECISION vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
