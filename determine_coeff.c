/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: determine_coeff.c,v 1.7 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"
/*

  determine coefficients for base function base_func
  with damping weights dfunc(i,n) where i the local mode
  (e.g, all unity, or 0 for i==0 (constant))

  c = coefficients [0...m-1]

  WARNING:
  
  y(x) is function given, y[0...ndata-1...ndata+2m-1] 

  (last 2m parameters set to zero for damping)

  x1, x2: limits of interval to be interpolated
  lamda: norm damping parameter, give in fractions of RMS
  omega: roughness damping
  size: output model size, length of c
  vr: variance reduction 1-misfit^2/length(y)^2


 */
#define FUNC(x1,x2,x3,x4,x5) ((*func)(x1,x2,x3,x4,x5))
#define DFUNC(i,j,k) ((*dfunc)(i,j,k))

void determine_coeff(COMP_PRECISION *c,int m,
		     COMP_PRECISION *x,COMP_PRECISION *y, 
		     int ndata,
		     COMP_PRECISION x1, COMP_PRECISION x2,
		     COMP_PRECISION lambda,
		     COMP_PRECISION omega,
		     COMP_PRECISION *modelsize, 
		     COMP_PRECISION *datasize,
		     COMP_PRECISION *chi2,
		     COMP_PRECISION (*func)(COMP_PRECISION,
					    COMP_PRECISION,
					    COMP_PRECISION *,
					    int,
					    COMP_PRECISION),
		     COMP_PRECISION (*dfunc)(int,int,int))
						 
{
  COMP_PRECISION *a,tmp,rms;
  int i,j,k,ndatam,ndatamm;
  static int warned=0;
  if((ndata < m)&&(!warned)){
    fprintf(stderr,"determine_coeff: nr of data (%i) < order (%i)\n",
	    ndata,m);
    warned=1;
  } 
  if(ndata<1){
    fprintf(stderr,"determine_coeff: need at least one data point, n: %i\n",
	    ndata);
    exit(-1);
  }
  ndatam=ndata+m;// for norm
  ndatamm=ndatam+m;// and roughness
  a=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m*ndatamm);
  if(!a)MEMERROR;
  // assemble matrix for least squares fit
  *datasize = 0.0;
  for(i=0;i<ndata;i++){
    for(j=0;j<m;j++){
      for(k=0;k<m;k++)
	c[k]=(k == j)?(1.0):(0.0);
      *(a+j*ndatamm+i) = FUNC(x1,x2,c,m,x[i]);
    }
    *datasize += y[i]*y[i];
  }
  rms= sqrt(*datasize/(COMP_PRECISION)ndata);
  // scale lambda and omega scaling factors to RMS of signal
  lambda *= rms;
  omega  *= rms;
  // norm damping
  for(k=0,i=ndata;i<ndatam;i++,k++)
    for(j=0;j<m;j++){
      *(a+j*ndatamm+i) = DFUNC(k,j,m)*lambda;
    }
  // roughness damping
  for(k=0,i=ndatam;i<ndatamm;i++,k++)
    for(j=0;j<m;j++){
      if(j == k)
	*(a+j*ndatamm+i) = omega*CHEBEV_DER_INT_FUNC(x1,x2,j);
      else
	*(a+j*ndatamm+i) = 0.0;
    }
  svd_solver(a,c,y,ndatamm,m);
  free(a);
  // solution vector length
  for(*modelsize=0.0,i=0;i<m;i++)
    *modelsize += c[i]*c[i];
  // chi2 
  for(*chi2=0.0,i=0;i<ndata;i++){
    tmp=FUNC(x1,x2,c,m,x[i]) - y[i];
    *chi2 += tmp*tmp;
  }
}

/*
  
  this is a driver for numerical recipes SVD routines
  used for least squares

 */
void svd_solver(COMP_PRECISION *a, 
		COMP_PRECISION *x,
		COMP_PRECISION *y,
		int ndata, int m)
{
  COMP_PRECISION *v,*w,wmax,wminlim,*work;
  int i;
  work=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  v=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m*m);
  w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  if((!v) || (!w) || (!work))
    MEMERROR;
  svd_driver(a,ndata,m,v,w);
  for(wmax=0.0,i=0;i<m;i++){
    if(w[i] > wmax)
      wmax=w[i];
  }
  // cutoff value for singular values
  wminlim=wmax * 1.0e-8;
  for(i=0;i<m;i++){
    if(w[i] < wminlim)
      w[i]=0.0;
  }
  svbksb(a,w,v,&ndata,&m,&ndata,&m,y,x,work);
  free(v);free(w);free(work);
}

void svd_driver(COMP_PRECISION *a,int m,int n,
		COMP_PRECISION *v,
		COMP_PRECISION *w)
{
  COMP_PRECISION *wrkarr;
  wrkarr=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n);
  if(!wrkarr)
    MEMERROR;
  svdcmp(a,&m,&n,&m,&n,w,v,wrkarr);
  free(wrkarr);
}

