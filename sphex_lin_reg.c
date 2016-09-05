/*

  calculates linear regression between spherical harmonic coefficients

  calls numerical recipes routines for data with erorrs only in sig_y
  (simple linear regression), or errors in both sig_x and sig_y
  (iterative least squares)

  $Id: sphex_lin_reg.c,v 1.4 2001/07/26 19:56:07 becker Exp becker $

*/
#include "shansyn.h"


void sphex_lin_reg(COMP_PRECISION *slope,COMP_PRECISION *sigma_slope,
		   COMP_PRECISION *a, COMP_PRECISION *b,COMP_PRECISION *c, 
		   COMP_PRECISION *d,int l,int regmode,int wmode,
		   COMP_PRECISION *sigma)
{
  COMP_PRECISION *x,*y,*sigx,*sigy,chi2,q;
  int m,lmin,lmax,n,lmsize,i,not_weighted=0;
  if(l<0){// sum over all l if l given negative
    lmin=  0;
    lmax= -l;
  }else{
    lmin = l;
    lmax = l;
  }
  // for all routines, we assume that we do not have the 
  // sigma values for the data
  //
  // first, create x and y value vectors, sigma will be unity
  //
  lmsize= (int)((((COMP_PRECISION)(lmax-lmin))+1.0)*
		  (((COMP_PRECISION)(lmax-lmin))+2.0));
  x=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lmsize);
  y=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lmsize);
  sigx=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lmsize);
  sigy=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lmsize);
  if(!x || !y || !sigx|| !sigy){
    fprintf(stderr,"memerror\n");exit(-1);
  }
  // assign x and y arrays
  for(n=-1,l=lmin;l <= lmax;l++){
    for(m=0;m<=l;m++){
      n++;
      x[n]= *(a+POSLM(l, m));
      y[n] = *(c+POSLM(l, m));
      if(m != 0){
	n++;
	x[n]= *(b+POSLM(l, m));
	y[n]= *(d+POSLM(l, m));
      }
    }
  }
  n++;
  // select a weighing
  switch(wmode){
  case 0:{// no weights
    for(i=0;i<n;i++)
      sigx[i]=sigy[i]=1.0;
    not_weighted=1;
    break;
  }
  case 1:{// some weights
    for(i=0;i<n;i++){
      sigx[i]=sigma[0];
      sigy[i]=sigma[1];
    }
    not_weighted=0;
    break;
  }
  case 2:{
    for(i=-1,l=lmin;l <= lmax;l++){
      for(m=0;m<=l;m++){
	i++;
	sigx[i]=((COMP_PRECISION)l+1)*sigma[0];
	sigy[i]=((COMP_PRECISION)l+1)*sigma[1];
	if(m != 0){
	  i++;
	  sigx[i]=((COMP_PRECISION)l+1)*sigma[0];
	  sigy[i]=((COMP_PRECISION)l+1)*sigma[1];
	}
      }
    }
    not_weighted=0;
    break;
  }
  default:{
    fprintf(stderr,"sphex_lin_reg: weighting mode %i is undefined\n",
	    wmode);
    exit(-1);
  }}
  //
  // call numrec style
  //
  switch(regmode){
  case REGRESS_XDETERMINED:{
    linreg_fit((x-1),(y-1),n,(sigx-1),not_weighted,
	       slope,(slope+1),sigma_slope,(sigma_slope+1),
	       &chi2,&q);
    break;
  }
  case REGRESS_ITERATIVE:{
    linreg_fitexy((x-1),(y-1),n,(sigx-1),(sigy-1),
		  slope,(slope+1),sigma_slope,(sigma_slope+1),
		  &chi2,&q);
    break;
  }
  default:{
    fprintf(stderr,"linear regression mode %i is undefined\n",
	    regmode);
    exit(-1);
    break;
  }}
  free(x);free(y);free(sigx);free(sigy);
}
