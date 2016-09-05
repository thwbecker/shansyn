/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: powercorr.c,v 1.10 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

/* 

   model[2]: input models, possibly with several layers

   l < 0 : sum from lmin ... -lmax
   l >= 0 : compute corr per degree
   
   lmin L minimum correlation
   cmode : 1 linear 2 pearson

   layer[2] selects which layer of the respective models to use

*/
COMP_PRECISION calc_correlation_model(struct mod *model, 
				      int l,int lmin,int cmode,
				      int weighted_for_degree,
				      int *layer)
{
  COMP_PRECISION corr;
  int i;
#ifdef DEBUG
  for(i=0;i<2;i++){
    if((layer[i] < 0) || (layer[i] >= model[i].n)){
      fprintf(stderr,"calc_correlation_model: ERROR: model %i only has %i layers, selected %i\n",
	      i+1,model[i].n,layer[i]);
      exit(-1);
    }
  }
  if(model[0].vector_field != model[1].vector_field){
    fprintf(stderr,"calc_correlatio_model: ERROR: vector field type mismath, %i vs. %i\n",
	    model[0].vector_field, model[1].vector_field);
    exit(-1);
  }
  if(model[0].is_gsh != model[1].is_gsh){
    fprintf(stderr,"calc_correlatio_model: ERROR: GSH type mismath, %i vs. %i\n",
	    model[0].is_gsh, model[1].is_gsh);
    exit(-1);
  }
  
#endif  
  switch(model[0].vector_field){
  case 0:
    corr = correlation(model[0].a[layer[0]],model[0].b[layer[0]],
		       model[1].a[layer[1]],model[1].b[layer[1]],
		       l,weighted_for_degree,lmin,cmode);
    break;
  case 1:
    corr = correlation_pt(model[0].amp[layer[0]],model[0].bmp[layer[0]],
			  model[0].amt[layer[0]],model[0].bmt[layer[0]],
			  model[1].amp[layer[1]],model[1].bmp[layer[1]],
			  model[1].amt[layer[1]],model[1].bmt[layer[1]],
			  l,weighted_for_degree,lmin,cmode);
    break;
  case 2:
    corr = correlation_gsh(model[0].amp[layer[0]],model[0].amt[layer[0]],
			   model[0].bmp[layer[0]],model[0].bmt[layer[0]],
			   model[1].amp[layer[1]],model[1].amt[layer[1]],
			   model[1].bmp[layer[1]],model[1].bmt[layer[1]],
			   l,weighted_for_degree,
			   (model[0].is_gsh-1),lmin,cmode);
    break;
  default:
    fprintf(stderr,"vector field type %i not implemented for out mode %i\n",
	    model[0].vector_field,cmode);
    exit(-1);
    break;
  }
  return corr;
}


COMP_PRECISION calc_total_power_model(struct mod *model,int lmax,
				      int layer)
{
  COMP_PRECISION power;
  switch(model->vector_field){	
  case 0:
    /* scalar */
    power = calc_total_power(model->a[layer],model->b[layer],lmax);
    break;
  case 1:
    /* vector */
    power = calc_total_power(model->amp[layer],model->bmp[layer],lmax) + 
      calc_total_power(model->amt[layer],model->bmt[layer],lmax);
    break;
  case 2:
    /* GSH */
    power = calc_total_power_gsh(model->amp[layer],model->amt[layer],
				 model->bmp[layer],model->bmt[layer],lmax);
    break;
  default:
    fprintf(stderr,"vector field type %i not implemented\n",
	    model->vector_field);
    exit(-1);
    break;
  }
  return power; 
}
COMP_PRECISION calc_rms_model(struct mod *model,int lmax,
			      int layer)
{
  COMP_PRECISION rms;
  switch(model->vector_field){	/* scalar */
  case 0:
    rms = calc_rms(model->a[layer],model->b[layer],lmax);
    break;
  case 1:
    rms = calc_rms(model->amp[layer],model->bmp[layer],lmax) + 
      calc_rms(model->amt[layer],model->bmt[layer],lmax);
    break;
  case 2:			/* GSH */
    rms = calc_rms_gsh(model->amp[layer],model->amt[layer],
		       model->bmp[layer],model->bmt[layer],lmax);
    break;
  default:
    fprintf(stderr,"vector field type %i not implemented \n",
	    model->vector_field);
    exit(-1);
    break;
  }
  return rms; 
}

/* power per degree */
COMP_PRECISION degree_power_model(struct mod *model,int l,int layer)
{
  COMP_PRECISION power;
  switch(model->vector_field){	/* scalar */
  case 0:
    power = degree_power(model->a[layer],model->b[layer],l);
    break;
  case 1:
     /* vector field using P + T, should
	really do differently
     */
      power = degree_power(model->amp[layer],model->bmp[layer],l) + 
	degree_power(model->amt[layer],model->bmt[layer],l);
      break;
  case 2:
    /* GSH */
    power = degree_power_gsh(model->amp[layer],model->amt[layer],
			     model->bmp[layer],model->bmt[layer],l);
    break;
  default:
    fprintf(stderr,"vector field type %i not implemented for out mode\n",
	    model->vector_field);
    exit(-1);
    break;
  }
  return power; 
}


/*

calculates power per degree and unit area,
normalized by the 2l+1 entries

(note: this is in units of value^2, i.e. take sqrt later

*/



COMP_PRECISION degree_power(COMP_PRECISION *a,
			    COMP_PRECISION *b, int l)
{
  COMP_PRECISION tmp;
  int m,nfac,os;
  // m=0
  tmp = SQUARE(a[POSLM(l, 0)]);
  // 1<=m<=l
  for(m=1;m<=l;m++){
    os = POSLM(l, m);
    tmp += SQUARE(a[os]);
    tmp += SQUARE(b[os]);
  }
  nfac = 2*l + 1;
  tmp /= (COMP_PRECISION)nfac;
  return tmp;
}

/*
  
  calculates power per degree and unit area for real/imaginary coeffs
  normalized by 4l+2 coefficients
  
*/
COMP_PRECISION degree_power_gsh(COMP_PRECISION *ar,
				COMP_PRECISION *ai,
				COMP_PRECISION *br,
				COMP_PRECISION *bi,
				int l)
{
  COMP_PRECISION tmp;
  int m,os,nfac;
  // m=0
  os = POSLM(l,0);
  tmp = SQUARE(ar[os]) + SQUARE(ai[os]);
  // 1<=m<=l
  for(m=1;m<=l;m++){
    os = POSLM(l, m);
    tmp += SQUARE(ar[os]) + SQUARE(ai[os]);
    tmp += SQUARE(br[os]) + SQUARE(bi[os]);
  }
  nfac = 4*l+2;
  tmp /= (COMP_PRECISION)nfac;
  return tmp;
}

/*

  calculate linear correlation coefficient between 
  two GSH 2phi or 4phi model expansions

  this sums up the dot product of all l >= 2 terms for 2phi (ialpha == 2), 
  and l >= 4 for 4phi (ialpha == 4)

*/
COMP_PRECISION correlation_gsh(COMP_PRECISION *ar, COMP_PRECISION *ai,
			       COMP_PRECISION *br, COMP_PRECISION *bi,
			       COMP_PRECISION *cr, COMP_PRECISION *ci, 
			       COMP_PRECISION *dr, COMP_PRECISION *di, 
			       int l, int weighted_for_degree,
			       int ialpha, int lmin,
			       int cmode)
{
  COMP_PRECISION w,*x=NULL,*y=NULL,prob,corr;
  int m,lmax,os,n=0;
  if(l < 0){
    //
    // sum over all nonzero ls, lmax = -l, if l given negative
    //
    switch(ialpha){
    case 2:
      lmin = MAX(2,lmin);
      break;
    case 4:
      lmin = MAX(4,lmin);
      break;
    default:
      fprintf(stderr,"correlation_gsh: ialpha %i undefined\n",ialpha);
      exit(-1);
    }
    lmax= -l;
  }else{
    lmin = l;
    lmax = l;
  }
  /* assemble parameters for correlation */
  for(l=lmin;l <= lmax;l++){
    if(weighted_for_degree)	/* weighted by the number of entries per degree (no good) */
      w = 1.0/(4.0*(COMP_PRECISION)l+2.0);
    else
      w = 1.0;
    for(m=0;m <= l;m++){
      //fprintf(stderr,"%11g %11g\t%11g %11g\n",ar[os],ai[os],cr[os],ci[os]);
      os = POSLM(l, m);
      add_to_xy(&x,&y,&n,ar[os] * w,cr[os] * w);
      add_to_xy(&x,&y,&n,ai[os] * w,ci[os] * w);
      if(m != 0){
	//fprintf(stderr,"%11g %11g\t%11g %11g\n",br[os],bi[os],dr[os],di[os]);
	add_to_xy(&x,&y,&n,br[os] * w,dr[os] * w);
	add_to_xy(&x,&y,&n,bi[os] * w,di[os] * w);
      }
    }
  }
  /* actually calculate correlation */
  corr = corr_sub(x,y,n,&prob,cmode);
  free(x);free(y);
  return(corr);
}
void add_to_xy(COMP_PRECISION **x, COMP_PRECISION **y, 
	       int *n, COMP_PRECISION xval, COMP_PRECISION yval)
{
  int n1;
  n1 = (*n) + 1;
  if(((*x = (COMP_PRECISION *)realloc(*x,sizeof(COMP_PRECISION)*n1))==NULL)||
     ((*y = (COMP_PRECISION *)realloc(*y,sizeof(COMP_PRECISION)*n1))==NULL)){
    fprintf(stderr,"add_to_xy: memory error\n");
    exit(-1);
  }
  *(*x + (*n)) = xval;
  *(*y + (*n)) = yval;
  *n = n1;
}
/*

calculate linear correlation coefficient between two model
expansions, weighted by each degree if weighted_by_degree is set
(this is no good)

*/
COMP_PRECISION correlation(COMP_PRECISION *a, COMP_PRECISION *b,
			   COMP_PRECISION *c, COMP_PRECISION *d, 
			   int l, int weighted_for_degree,
			   int lmin,int cmode)
{
  COMP_PRECISION w,corr,prob,*x=NULL,*y=NULL;
  int m,lmax,os,n=0;
  if(l < 0){// sum over all l if l given negative
    lmin=  MAX(1,lmin);
    lmax= -l;
  }else{
    lmin = l;
    lmax = l;
  }
  for(l=lmin;l <= lmax;l++){
    if(weighted_for_degree)	/* weighted by the number of entries per degree */
      w = 1.0/(2.0*(COMP_PRECISION)l+1.0);
    else
      w = 1.0;			/* better */
    for(m=0;m<=l;m++){
      os = POSLM(l, m);
      add_to_xy(&x,&y,&n,a[os] * w,c[os] * w);
      if(m != 0)
	add_to_xy(&x,&y,&n,b[os] * w,d[os] * w);
    }
  }
  corr = corr_sub(x,y,n,&prob,cmode);
  free(x);free(y);
  return(corr);
}

/* correlation for a vector field */
COMP_PRECISION correlation_pt(COMP_PRECISION *ap, COMP_PRECISION *bp, /* first file poloidal */
			      COMP_PRECISION *at, COMP_PRECISION *bt, /* first file toroidal */
			      COMP_PRECISION *cp, COMP_PRECISION *dp, /* second file pol  */
			      COMP_PRECISION *ct, COMP_PRECISION *dt, /* second file tor */
			      int l, int weighted_for_degree, int lmin,
			      int cmode)
{
  COMP_PRECISION w,corr,prob,*x=NULL,*y=NULL;
  int m,lmax,os,n=0;
  if(l < 0){// sum over all l if l given negative
    lmin=  MAX(1,lmin);
    lmax= -l;
  }else{
    lmin = l;
    lmax = l;
  }
  for(l=lmin;l <= lmax;l++){
    if(weighted_for_degree)	/* weighted by the number of entries per degree (no good) */
      w = 1.0/(4.0*(COMP_PRECISION)l+2.0);
    else
      w = 1.0;	
    for(m=0;m <= l;m++){
      os = POSLM(l, m);
      add_to_xy(&x,&y,&n,ap[os] * w,cp[os] * w);
      add_to_xy(&x,&y,&n,at[os] * w,ct[os] * w);
      if(m != 0){
	add_to_xy(&x,&y,&n,bp[os] * w,dp[os] * w);
	add_to_xy(&x,&y,&n,bt[os] * w,dt[os] * w);
      }
    }
  }
  corr = corr_sub(x,y,n,&prob,cmode);
  free(x);free(y);
  return(corr);
}

COMP_PRECISION ccl_correlation(COMP_PRECISION *a, COMP_PRECISION *b,
			       int l,int lmax, int *dof,int cmode)
{
  COMP_PRECISION corr,prob,*x=NULL,*y=NULL;
  int m,os,os2,n=0;
  if(l >= lmax){
    fprintf(stderr,"ccl_corr: can only compute r_l,l+1 up to %i-1=%i\n",
	    lmax,lmax-1);
    exit(-1);
  }
  for(m=0;m <= l;m++){
    os =  POSLM(l, m);		/* l,  m */
    os2 = POSLM(l+1, m);	/* l+1,m */
    add_to_xy(&x,&y,&n,a[os],a[os2]);
    if(m != 0)
      add_to_xy(&x,&y,&n,b[os],b[os2]);
  }
  *dof = n;
  corr = corr_sub(x,y,n,&prob,cmode);
  free(x);free(y);
  return(corr);
}


/*

  calculates correlation based on 1/(2l+1) coefficients
  (not good)

 */
COMP_PRECISION weighted_correlation(COMP_PRECISION *a, COMP_PRECISION *b,
				    COMP_PRECISION *c, COMP_PRECISION *d,
				    int l,int lmin,
				    int cmode)
{
  return correlation(a,b,c,d,l,1,lmin,cmode);
}


/*
  calculates RMS normalized by surface area of sphere, without l = 0 term
*/
COMP_PRECISION calc_rms(COMP_PRECISION *a, COMP_PRECISION *b,
			int lmax)
{
  int l;
  COMP_PRECISION rms;
  rms=0.0;
  for(l=1;l<=lmax;l++)
    rms += ((COMP_PRECISION)(2*l+1))*degree_power(a,b,l);
  return(sqrt(rms)/TWO_SQRT_PI);
}
/*
  same for gsh
*/
COMP_PRECISION calc_rms_gsh(COMP_PRECISION *ar, 
			    COMP_PRECISION *ai,
			    COMP_PRECISION *br, 
			    COMP_PRECISION *bi,int lmax)
{
  int l;
  COMP_PRECISION rms;
  rms=0.0;
  for(l=1;l<=lmax;l++){		/* this is OK, we're summing over zero
				   terms */
    rms += ((COMP_PRECISION)(4*l+2))*
      degree_power_gsh(ar,ai,br,bi,l);
  }
  return(sqrt(rms)/TWO_SQRT_PI);
}
/*
  total power normalized by surface area of sphere, including l =0 term
*/

COMP_PRECISION calc_total_power(COMP_PRECISION *a, 
				COMP_PRECISION *b,int lmax)
{
  int l;
  COMP_PRECISION rms;
  rms=0.0;
  for(l=0;l<=lmax;l++)
    rms += ((COMP_PRECISION)(2*l+1))*degree_power(a,b,l);
  return(sqrt(rms)/TWO_SQRT_PI);
}

/*
  total power normalized by surface area of sphere
*/

COMP_PRECISION calc_total_power_gsh(COMP_PRECISION *ar, 
				    COMP_PRECISION *ai,
				    COMP_PRECISION *br, 
				    COMP_PRECISION *bi,int lmax)
{
  int l;
  COMP_PRECISION rms;
  rms=0.0;
  for(l=0;l<=lmax;l++)
    rms += ((COMP_PRECISION)(4*l+2))*
      degree_power_gsh(ar,ai,br,bi,l);
  return(sqrt(rms)/TWO_SQRT_PI);
}

