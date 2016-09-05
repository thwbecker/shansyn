/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: plgndr.c,v 1.16 2009/04/16 00:55:29 becker Exp becker $ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "shansyn.h"
/* 
   this source file holds all Legendre funciton routines,
   some are very specialized and operate on blocks of 
   data or include Gauss integration weights

   part of shana and shsyn

   (C) Thorsten Becker, becker@eps.harvard.edu 


   Calculate associated Legendre polynomials of max order 
   lmax for vector x. Recursion formula modified from, e.g., 
   Numerical Recipes. Normalization follows, 
   e.g., Dahlen and Tromp, for fully normalized, physics 
   convention (integral over unit sphere = 1).
   
   X_l^m = ( (2l+1)/4/pi (l-m)!/(l+m)! )^0.5 P'_l^m 
   P'_l^m = (-1)^m * (1-z^2)^{m/2} d^m/dz^m P_l(z)
   P_l(z)= 1 / 2^l / l! d^l/dz^l (z^2-1)^l 

   Routine calculates and stores all functions in P, 
   intended to be called once.
*/

/* 1/(2*sqrt(pi)) */
#define NORM_FAC_FOR_ALL_M 0.28209479177387814347403972578039

/* 
   Legendre functions 
*/

/*
  
  routine to create a whole array P(l,m,j) that holds
  the Legendre function of order l,m for the j-th point
  in latitudinal (y) direction

  y is given as cos(theta)

*/

void  plgndr(COMP_PRECISION *y, int nlat, 
	     COMP_PRECISION *p, int lmax)
{
  COMP_PRECISION som,tmp,fact1,fact2;
  int i,l,j,m,intfact1,intfact2,lmsize,twol,lpm,lmm;
  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*
		(((COMP_PRECISION)lmax)+2)/2.0);
  
  for(j=0;j<nlat;j++){  /* loop over all y */

    for(m=0;m<=lmax;m++){ /* loop over all m */

      /* compute P m m */

      P(m, m, j) = NORM_FAC_FOR_ALL_M;

      if (m > 0) {
	/* |sin(theta)| */
       	som=sqrt((1.0- y[j])*(1.0+ y[j]));  
       	for (i=1,intfact1=1,intfact2=2,tmp=1.0; 
       	     i<=m; 
       	     i++,intfact1 += 2,intfact2 += 2){ 
       	  P(m,m,j) *= -som;   
       	  tmp  *= (((COMP_PRECISION)intfact1)/
		   ((COMP_PRECISION)intfact2)); 
       	} 
       	P(m,m,j) *= sqrt(tmp); 
       	P(m,m,j) *= sqrt((COMP_PRECISION)(1+2*m)); 
      }  

      if(m==lmax)
	continue;

      /* compute p m+1 m */
      /* the factor is Norm_{m+1}^m/Norm_m^m (2m+1) */
	 
      P(m+1,m,j) = y[j] * sqrt((COMP_PRECISION)(3+2*m)) * P(m,m,j);
    

      /* compute the rest, i.e. P l m up to P lmax m  from 
	 stable recurrence formula P_l^m= fact1 * x * P_{l-1}^m - fact2 * P_{l-2}^m */
      
      for (l=m+2,twol=2*m+4;
	   l<=lmax;
	   l++,twol += 2) {
	//
	lpm=l+m;
	lmm=l-m;
	/* this is Norm_l^m / Norm_{l-1}^m * (l+m-1)/(l-m), ergo
	   Sqrt((-1 + 4*Power(L,2))/(Power(L,2) - Power(M,2))) */
	fact1  = (COMP_PRECISION)((twol+1)*(twol-1));
	fact1 /= (COMP_PRECISION)(lpm*lmm);
	fact1  = sqrt(fact1);
	
	/* this is Norm_l^m / Norm_{l-2}^m * (l+m-1)/(l-m), ergo
	   Sqrt(((1 + 2*L)*(-1 + L - M)*(-1 + L + M))/
	   ((-3 + 2*L)*(L - M)*(L + M))) */

	fact2  = ((COMP_PRECISION)(twol+1)*(COMP_PRECISION)(lmm-1)*(COMP_PRECISION)(lpm-1));
	fact2 /= ((COMP_PRECISION)(twol-3)*(COMP_PRECISION)(lmm)*(COMP_PRECISION)(lpm));
	fact2  = sqrt(fact2);

	P(l,m,j)  = (y[j]* fact1 * P(l-1,m,j));
	P(l,m,j) -= (      fact2 * P(l-2,m,j));
      }
    }
  }

}
/*

  same as above but we calculate the Legendre functions
  at the actual Gauss integration abcissae values
  also include the weighting factor for Gaussian quadrature
  
  nr_gauss_pts should be set to lmax+1

  the abscissae values are output to absc(0...nr_gauss_pts-1)

*/
void  plgndr_g(COMP_PRECISION *p, int lmax,
	       int nr_gauss_pts,
	       COMP_PRECISION *absc)
{
  COMP_PRECISION *w;
  int i,l,m,lmsize;
  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*
		(((COMP_PRECISION)lmax)+2)/2.0);
  if((w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nr_gauss_pts))
     ==NULL){
    fprintf(stderr,"memerror in plgndr_g for w\n");
    exit(-1);
  }
  // obtain abscissae values
  // call numerical recipes style, with shifted index
  gauleg_orig(-1.0,1.0,absc-1,w-1,nr_gauss_pts);
  // get Legendre functions at the Gauss points
  plgndr(absc, nr_gauss_pts,p,lmax);
  // scale with weight
  for(i=0;i<nr_gauss_pts;i++)
    for(l=0;l<=lmax;l++)
      for(m=0;m<=l;m++)
	P(l,m,i) *= w[i];
  free(w);
}

/* 
   this version takes a twodimensional array of 
   coordinates and uses the second
   for the y WHILE THE COSINE HAS TO BE TAKEN HERE 
   (special case of plgndr)

   y is given as a 2D array with phi, theta
*/

void  plgndr2(DATA_PRECISION *y, int nlat, COMP_PRECISION *p, int lmax)
{ 
  COMP_PRECISION som,tmp,fact1,fact2,sinus;
  int i,l,j,m,intfact1,intfact2,lmsize,twol,lpm,lmm;
  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*(((COMP_PRECISION)lmax)+2)/2.0);
  
  for(j=0;j<nlat;j++){  /* loop over all y */

    for(m=0;m<=lmax;m++){ /* loop over all m */

      /* compute P m m */      

      P(m, m, j)=NORM_FAC_FOR_ALL_M;

      if (m > 0) {  
	/* |sin(theta)| */
	sinus=sin(*(y+j*2+1));
       	som=(sinus >= 0)?(sinus):(-sinus);
       	for (i=1,intfact1=1,intfact2=2,tmp=1.0; 
       	     i<=m; 
       	     i++,intfact1 += 2,intfact2 += 2){ 
       	  P(m,m,j) *= -som;   
       	  tmp  *= (((COMP_PRECISION)intfact1)/
		   ((COMP_PRECISION)intfact2)); 
       	} 
       	P(m,m,j) *= sqrt(tmp); 
       	P(m,m,j) *= sqrt((COMP_PRECISION)(1+2*m)); 
      }  
      if(m==lmax)
	continue;
      /* compute p m+1 m */
      P(m+1,m,j) = cos(*(y+j*2+1)) * sqrt((COMP_PRECISION)(3+2*m)) * P(m,m,j);
      /* compute the rest, i.e. P l m up to P lmax m  from 
	 stable recurrence formula P_l^m= fact1 * x * P_{l-1}^m - fact2 * P_{l-2}^m */
      for (l=m+2,twol=2*m+4;
	   l<=lmax;
	   l++,twol += 2) {
	lpm=l+m;
	lmm=l-m;
	/* this is Norm_l^m / Norm_{l-1}^m * (l+m-1)/(l-m), ergo
	   Sqrt((-1 + 4*Power(L,2))/(Power(L,2) - Power(M,2))) */
	fact1  = (COMP_PRECISION)((twol+1)*(twol-1));
	fact1 /= (COMP_PRECISION)(lpm*lmm);
	fact1  = sqrt(fact1);
	
	/* this is Norm_l^m / Norm_{l-2}^m * (l+m-1)/(l-m), ergo
	   Sqrt(((1 + 2*L)*(-1 + L - M)*(-1 + L + M))/
	   ((-3 + 2*L)*(L - M)*(L + M))) */

	fact2  = ((COMP_PRECISION)(twol+1)*(COMP_PRECISION)(lmm-1)*(COMP_PRECISION)(lpm-1));
	fact2 /= ((COMP_PRECISION)(twol-3)*(COMP_PRECISION)(lmm  )*(COMP_PRECISION)(lpm  ));
	fact2  = sqrt(fact2);

	P(l,m,j)  = (cos(y[j*2+1])* fact1 * P(l-1,m,j));
	P(l,m,j) -= (               fact2 * P(l-2,m,j));
      }
    }
  }

}
/*

  routine calculates single value of legendre function
  of order l,m at point x

*/

COMP_PRECISION slgndr(int l,int m, COMP_PRECISION x)
{
  COMP_PRECISION pll,pmm,pmmp1,som,tmp,fact1,fact2;
  int i,ll,intfact1,intfact2,twol,lpm,lmm;
#ifdef DEBUG
  if (m < 0 || m > l || fabs(x) > 1.0){
    fprintf(stderr,"slgndr: bad arguments");
    exit(-1);
  }
#endif
  pll = 0.;			/* for compiler */
  pmm=NORM_FAC_FOR_ALL_M;

  if (m > 0) {
    som=sqrt((1.0-x)*(1.0+x));
    for (i=1,intfact1=1,intfact2=2,tmp=1.0; 
	 i<=m; 
	 i++,intfact1 += 2,intfact2 += 2){ 
      pmm *= -som;   
      tmp *= (((COMP_PRECISION)intfact1)/
	      ((COMP_PRECISION)intfact2)); 
    }
    pmm *= sqrt(tmp);
    pmm *= sqrt((COMP_PRECISION)(1+2*m)); 
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*sqrt((COMP_PRECISION)(3+2*m))*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2,twol=2*m+4;
	   ll<=l;
	   ll++,twol += 2) {
	lpm=ll+m;
	lmm=ll-m;
	fact1  = (COMP_PRECISION)((twol+1)*(twol-1));
	fact1 /= (COMP_PRECISION)(lpm*lmm);
	fact1  = sqrt(fact1);
	fact2  = ((COMP_PRECISION)(twol+1)*(COMP_PRECISION)(lmm-1)*(COMP_PRECISION)(lpm-1));
	fact2 /= ((COMP_PRECISION)(twol-3)*(COMP_PRECISION)(lmm)*(COMP_PRECISION)(lpm));
	fact2  = sqrt(fact2);

	pll  = (x * fact1 * pmmp1);
	pll -= (    fact2 * pmm  );

	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}

/* 
   calculate the derivative factors for the Legendre 
   functions 
   
   sort of, since we are really calculating

   d X_lm / d theta where 

   X_lm = NormFac(l,m) P_lm(cos(theta))

   as in B.58 of Dahlen and Tromp

   therefore, P_lm=X_lm/NormFac(l,m) and dP_lm(\mu)/d \mu = 
   - 1/(NormFac(l,m) *sin(theta)) d X_lm(theta)/dtheta
   
   if calc_second is set, calculate second derivatives

   y is given as cos(theta)

*/
/*

  this is not well tested

  - we want a more elegant way to deal with y= +/-1 arguments
  
  - stability of this routine has to be established

*/


void  pdtheta_lgndr(COMP_PRECISION *y, int nlat, COMP_PRECISION *p, 
		    COMP_PRECISION *dptheta, int lmax, 
		    COMP_PRECISION *dp2theta,int calc_second)
{
  COMP_PRECISION fact1,fact2,*trigf;
  int l,j,m,lmsize,twol,lmm,lpm;
  lmsize= (int)((((COMP_PRECISION)lmax)+1.0)*(((COMP_PRECISION)lmax)+2)/2.0);
  trigf=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nlat);
  if(!trigf)
    MEMERROR;
  for(j=0;j<nlat;j++){
    trigf[j]=1.0-SQUARE(y[j]);
    if(trigf[j] <= EPS_COMP_PREC)
      trigf[j]=EPS_COMP_PREC;
    trigf[j]= 1.0/sqrt(trigf[j]);
    //fprintf(stderr,"%22.16e %22.16e \n",y[j],trigf[j]);
  }
  if(lmax < 2){
    for(j=0;j<nlat;j++){
      DPTHETA(0,0,j) = 0.0;
      if(lmax == 1){
	DPTHETA(1,0,j) =  -trigf[j] * ( SQRT_THREE * P(0,0,j) - y[j] * P(1,0,j));
	DPTHETA(1,1,j) =   trigf[j] * y[j] * P(1,1,j);
      }
      if(calc_second){/* for formula, see below */
	DP2THETA(0,0,j) = 0.0;
	if(lmax == 1){
	  DP2THETA(1,0,j)=   trigf[j]  * y[j] * DPTHETA(1,0,j) - P(1,0,j);
	  DP2THETA(1,0,j)=   trigf[j]  * y[j] * DPTHETA(1,0,j) - P(1,0,j);
	}
      }
    }
  }else{
    for(j=0;j<nlat;j++){  /* loop over all y */
      DPTHETA(0,0,j) = 0.0;
      DPTHETA(1,0,j) =  -trigf[j] * ( SQRT_THREE * P(0,0,j) - y[j] * P(1,0,j));
      DPTHETA(1,1,j) =   trigf[j] * y[j] * P(1,1,j);
      for(l=2,twol=4;
	  l<=lmax;
	  l++,twol += 2){ /* loop over all l for m=0 and m=l */
	/* 
	   for l 0 use  
	   dX_l0/dtheta= -l/sqrt(1-y^2) {sqrt((2l+1)(2l-1)) X_l-1,0 - y X_l,0}
	*/
	DPTHETA(l,0,j)  = sqrt(((COMP_PRECISION)(twol+1))/((COMP_PRECISION)(twol-1))) * P(l-1,0,j);
	DPTHETA(l,0,j) -= (y[j] * P(l,0,j));
	DPTHETA(l,0,j) *= -trigf[j] * ((COMP_PRECISION)l);
	/*

	  d X_ll /d theta = l cot(theta) X_ll 
	                  = l X_ll y/sqrt(1-y^2)
			  from B.123
	*/
	DPTHETA(l,l,j) =  trigf[j] * y[j] * ((COMP_PRECISION)l) * P(l,l,j);
	/* 
	   l m  from B.120, which says
	   d X_lm/d theta = 1/2 sqrt((l-m)(l+m+1)) X_l,m+1 - 1/2 sqrt ((l+m)(l-m+1)) X_l, m-1
	*/
	for(m=1;m<l;m++){ /* loop over all m from 1 to lmax-1 */
	  lmm=l-m;
	  lpm=l+m;
	  fact1=sqrt((COMP_PRECISION)(lmm*(lpm+1)));
	  fact2=sqrt((COMP_PRECISION)(lpm*(lmm+1)));
	  DPTHETA(l,m,j)  = 0.5 * ( fact1 * P(l,m+1,j) - fact2 * P(l,m-1,j));
	}
#ifdef DEBUG
	/* check floating exceptions */
	for(m=0;m<=l;m++)
	  if(!finite(DPTHETA(l,m,j))){
	    fprintf(stderr,"pdtheta_lgndr: j: %i y: %g l: %i m: %i DP: %g Pmp1: %g Pmm1: %g \n",
		    j,y[j],l,m,DPTHETA(l,m,j),P(l,m+1,j),P(l,m-1,j));
	    exit(-1);
	  }
#endif
      }
    }
  }
  free(trigf);
}
