/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: plotlegendre.c,v 1.7 2001/05/19 15:51:57 becker Exp becker $ */
#include "shansyn.h"
/*
  
  driver routine for various Legendre function implementations
  main purpose is visualization and benchmarking

  
*/

void gauleg_orig(COMP_PRECISION ,COMP_PRECISION ,
		COMP_PRECISION *,COMP_PRECISION *,
		int );
#define PLG_OUTPUT_FORMAT "%22.15e %22.15e\n"

int main(int argc, char *argv[] )
{
  COMP_PRECISION *y, *p,fac,*abscissas,*weights,z,*dptheta,dummy,dz;
  int lmax,l,m,i,j,nlat,lmsize,form;
  switch(argc){
  case 5:{
    sscanf(argv[1],"%i",&l);
    sscanf(argv[2],"%i",&m);
    sscanf(argv[3],"%i",&nlat); 
    sscanf(argv[4],"%i",&form);
    break;
  }case 4:{
    sscanf(argv[1],"%i",&l);
    sscanf(argv[2],"%i",&m);
    sscanf(argv[3],"%i",&nlat);
    form=0;
    break;
  }
  case 3:{
    sscanf(argv[1],"%i",&l);
    sscanf(argv[2],"%i",&m);
    nlat=101;form=0;
    break;
  }
  default:{
    fprintf(stderr,"%s l m [n, 100] [form, 0]\n",argv[0]);
    fprintf(stderr,"\t prints associated Legendre function X_l^m (cos(theta))\n");
    fprintf(stderr,"\t at n evenly spaced points from -1 <= cos(theta) <= 1\n");
    fprintf(stderr,"\t output is in the format\n\n");
    fprintf(stderr,"\t cos(theta) X_l^m\n\n");
    fprintf(stderr,"\t form: 0 for physical normalization, X_l^m as in Dahlen and Tromp\n");
    fprintf(stderr,"\t       1 for geodetical normalization\n");
    fprintf(stderr,"\t       2 for Gauss quadrature points with weights\n");
    fprintf(stderr,"\t       3 for physical normalization, single points\n");
    fprintf(stderr,"\t       4 dP_l^m/d (ArcCos(x)), first derivatives\n");
    exit(-1);
  }}
  if(m > l || m < 0 || l < 0){
    fprintf(stderr,"%s: l=%i and m=%i is nonsense\n",argv[0],l,m);
    exit(-1);
  }
  if(nlat<2){
    nlat=2;
    fprintf(stderr,"%s: increasing nlat to %i\n",argv[0],nlat);
  }
  lmax=l;
  lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
  if(form == 2){
    fprintf(stderr,"%s: Assoc. Legendre Pol_{%i}^{%i}, Gauss quad. points\n",argv[0],l,m);
    if(((abscissas=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(lmax+1)))==NULL)||
       ((weights=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(lmax+1)))==NULL)){
      fprintf(stderr,"%s: memerror for Gauss\n",argv[0]);exit(-1);}
    gauleg_orig(1.0,-1.0,abscissas-1,weights-1,lmax+1);
    for(i=0;i<lmax+1;i++)
      fprintf(stdout,PLG_OUTPUT_FORMAT,*(abscissas+i),*(weights+i));
  }else{
    if(form != 3){
      y=(COMP_PRECISION *)calloc(sizeof(COMP_PRECISION),nlat);
      /* y coordinate */
      if((y=(COMP_PRECISION *)calloc(nlat,sizeof(COMP_PRECISION)))==NULL)
	{fprintf(stderr,"memerror\n");exit(-1);};
      j=nlat-1;
      for(i=0;i<nlat;i++)
	*(y+i)=-1.0+(COMP_PRECISION)i/(COMP_PRECISION)j*2.0;
      
      if((p=(COMP_PRECISION *)calloc(lmsize*nlat,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: switching to single function mode!\n",argv[0]);
	fprintf(stderr,"%s: Assoc. Legendre Pol_{%i}^{%i}, physical convention\n",argv[0],l,m);
	if(form == 4){
	  fprintf(stderr,"%s: will only calculate derivatives for Legendre table\n",argv[0]);
	  exit(-1);
	}
	dz=2.0/(COMP_PRECISION)(nlat-1);
	for(z= -1.0;z <= 1.0+1.0e-14;z += dz) 
	  fprintf(stdout,PLG_OUTPUT_FORMAT,z, slgndr(l,m,z)); 
	exit(-1);
      }else{
	plgndr(y,nlat,p,lmax);
	if(form == 4){
	  fprintf(stderr,"%s: calculating derivatives, dP_lm/d theta, for l=%i m=%i\n",argv[0],l,m);
	  if((dptheta=(COMP_PRECISION *)calloc(lmsize*nlat,sizeof(COMP_PRECISION)))==NULL){
	    fprintf(stderr,"%s: not enough memory for second table for dP/dtheta\n",argv[0]);
	    exit(-1);
	  }
	  pdtheta_lgndr(y,nlat,p,dptheta,lmax,&dummy,0);
	  for(j=0;j<nlat;j++)
	    fprintf(stdout,PLG_OUTPUT_FORMAT,*(y+j),DPTHETA(l, m, j));
	  exit(-1);
	}
      }
      if(form == 1){
	fac  = pow(-1.0,(COMP_PRECISION)m);
	fac *= TWO_SQRT_PI;
	if(m != 0)
	  fac *= sqrt(2.0);
	fprintf(stderr,"%s: Assoc. Legendre Pol_{%i}^{%i}, geodetic convention\n",argv[0],l,m);
      }else{
	fprintf(stderr,"%s: Assoc. Legendre Pol_{%i}^{%i}, physical convention\n",argv[0],l,m);
	fac=1.0;
      }
      for(j=0;j<nlat;j++)
	fprintf(stdout,PLG_OUTPUT_FORMAT,*(y+j),fac*P(l, m, j));
    }else{
      fprintf(stderr,"%s: Assoc. Legendre Pol_{%i}^{%i}, physical convention\n",argv[0],l,m);
      dz=2.0/(COMP_PRECISION)(nlat-1);
      for(z= -1.0;z <= 1.0+1.0e-14;z += dz)
	fprintf(stdout,PLG_OUTPUT_FORMAT,z, slgndr(l,m,z));
    }
  } 
  return 0;
}








