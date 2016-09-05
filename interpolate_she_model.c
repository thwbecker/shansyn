/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: interpolate_she_model.c,v 1.11 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

/*

  interpolate a SHE at depth z

*/

int interpolate_she_model(struct mod *out_model,
			  struct mod *model, 
			  COMP_PRECISION z, 
			  int lmax,
			  int allow_lin_model_extrapolation) /*  */
{
  int j,k,l,m,i,interpolate;
  COMP_PRECISION fac,fac2,*ca,*cb;
  
  switch(model->radial_type){
  case DISCRETE:
    // discrete layers, interpolate layers somehow
    
    /* 


       linear interpolation 

    */
    // check for levels
    for(i=1;i<model->n;i++){
      if(model->d[i] >= model->d[i-1]){
	fprintf(stderr,"interpolate_she_model: depth levels in discrete model not ordered in descending order\n");
	exit(-1);
      }
    }
    
    if((z > model->d[0]) || 
       (z < model->d[model->n-1])){ /* outside bounds */
      if(allow_lin_model_extrapolation){
	interpolate = TRUE;
	fprintf(stderr,"interpolate_she_model: WARNING: extrapolating at z: %g\n",z);
      }else
	interpolate = FALSE;
    }else{
      interpolate = TRUE;
    }
    if(interpolate){		/* regular operation */
      i = model->n-1;
      while((i>0)&&(model->d[i]<z))
	i--;
      if(i == (model->n - 1))
	i=model->n-2;
      j=i+1;
      fac=(z - model->d[j])/(model->d[i]-model->d[j]);
      fac2=1.0-fac;
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  k=POSLM(l,m);
	  if(out_model->vector_field){
	    *(out_model->amp[0]+k) =
	      fac  * *(model->amp[i]+k) + fac2 * *(model->amp[j]+k);
	    *(out_model->bmp[0]+k) =
	      fac  * *(model->bmp[i]+k) + fac2 * *(model->bmp[j]+k);
	    *(out_model->amt[0]+k) =
	      fac  * *(model->amt[i]+k) + fac2 * *(model->amt[j]+k);
	    *(out_model->bmt[0]+k) =
	      fac  * *(model->bmt[i]+k) + fac2 * *(model->bmt[j]+k);
	  }else{
	    *(out_model->a[0]+k) =
	      fac  * *(model->a[i]+k) + fac2 * *(model->a[j]+k);
	    *(out_model->b[0]+k) =
	      fac  * *(model->b[i]+k) + fac2 * *(model->b[j]+k);
	  }
	}
    }else{
      //fprintf(stderr,"out of bounds %g %g %g\n",z,model->d[0],model->d[model->n-1]);
      /* return NaN */
      for(l=0;l <= lmax;l++)
	for(m=0;m <= l;m++){
	  k = POSLM(l,m);
	  if(out_model->vector_field){
	    *(out_model->amp[0]+k) = my_make_nan();
	    *(out_model->bmp[0]+k) = my_make_nan();
	    *(out_model->amt[0]+k) = my_make_nan();
	    *(out_model->bmt[0]+k) = my_make_nan();
	  }else{
	    *(out_model->a[0]+k) = my_make_nan();
	    *(out_model->b[0]+k) = my_make_nan();
	  }
	}
      return 1;
    }

    break;
  case CHEBYSHEV:
    if(out_model->vector_field){
      fprintf(stderr,"interpolate: vector field Chebychev not implemented yet\n");
      exit(-1);
    }
    if((ca=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*
				   model->n))==NULL)MEMERROR;
    if((cb=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*
				   model->n))==NULL)MEMERROR;
    for(l=0;l <= lmax;l++){
      for(m=0;m <= l;m++){
	k = POSLM(l,m);
	for(i=0;i<model->n;i++){
	  ca[i]= *(model->a[i]+k);
	  cb[i]= *(model->b[i]+k);
	}
	*(out_model->a[0]+k) = 
	  CHEBEV_FUNC(model->dmin,model->dmax,ca,model->n,z);
	*(out_model->b[0]+k) = 
	  CHEBEV_FUNC(model->dmin,model->dmax,cb,model->n,z);
      }
    }
    free(ca);free(cb);
    break;
  case SPLINES:
    if(out_model->vector_field){
      fprintf(stderr,"interpolate: vector field Spliines not implemented yet\n");
      exit(-1);
    }
    if((ca=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*
				    model->n))==NULL)MEMERROR;
    if((cb=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*
				    model->n))==NULL)MEMERROR;
    for(l=0;l <= lmax;l++)
      for(m=0;m <= l;m++){
	k = POSLM(l,m);
	for(i=0;i<model->n;i++){
	  ca[i]= *(model->a[i]+k);
	  cb[i]= *(model->b[i]+k);
	}
	*(out_model->a[0]+k) = 
	  spline_base(model->dmin,model->dmax,ca,model->n,z);
	*(out_model->b[0]+k) = 
	  spline_base(model->dmin,model->dmax,cb,model->n,z);
      }
    free(ca);free(cb);
    break;
  default:
    fprintf(stderr,"interpolate_she_model: can not deal with model type %i\n",
	    model->radial_type);
    exit(-1);
    break;
  }
  return 0;
} 

/*

  obtain nmodel mean expansions by averaging 
  nmodel spherical harmonic models 
  from z1 to z2 in steps steps
  

  as opposed to interpolate_she_model, arrays are allocated here



*/

void mean_model(struct mod *out_model,struct mod *model, 
		COMP_PRECISION z1, COMP_PRECISION z2,
		int steps,int use_r2_weights)
{
  COMP_PRECISION *a,*b,dz,z,r,w,ws=0.;
  int i,l,m,lmsize,j=0,iread;
  struct mod int_model[1];
  /* make room for temp model */
  copy_model_par(model,int_model);
  int_model->n = 1;
  allocate_model_coefficients(int_model);

  if(z2 < z1){
    fprintf(stderr,"mean_expansions: z2 should be bigger than z1, %g %g\n",
	    z1,z2);
    exit(-1);
  }
  // init mean expansions with zeroes
  for(l=0;l <= int_model->lmax;l++){
    for(m=0;m <= l;m++){
      if(int_model->vector_field){
	*(out_model->amp[0]+POSLM(l,m)) = *(out_model->bmp[0]+POSLM(l,m)) = 0.0;
	*(out_model->amt[0]+POSLM(l,m)) = *(out_model->bmt[0]+POSLM(l,m)) = 0.0;
      }else{
	*(out_model->a[0]+POSLM(l,m)) = *(out_model->b[0]+POSLM(l,m)) = 0.0;
      }
    }
  }
  //
  // delta z
  //
  dz=(z2 - z1)/(COMP_PRECISION)(steps-1);
  if(dz > 0){
    ws = 0.0;// sum of weights    
    // average from z1 to z2
    for(j=0,z=z1;z <= z2+1e-5;z += dz,j++){
      interpolate_she_model(int_model,model,z,int_model->lmax,TRUE); /* allow extrapolation */
      if(use_r2_weights == 0){// simple mean
	for(l=0;l <= int_model->lmax;l++)// add to mean
	  for(m=0;m <= l;m++){
	    if(int_model->vector_field){
	      out_model->amp[0][POSLM(l,m)] += int_model->amp[0][POSLM(l,m)];
	      out_model->bmp[0][POSLM(l,m)] += int_model->bmp[0][POSLM(l,m)];
	      out_model->amt[0][POSLM(l,m)] += int_model->amt[0][POSLM(l,m)];
	      out_model->bmt[0][POSLM(l,m)] += int_model->bmt[0][POSLM(l,m)];
	    }else{
	      out_model->a[0][POSLM(l,m)] += int_model->a[0][POSLM(l,m)];
	      out_model->b[0][POSLM(l,m)] += int_model->b[0][POSLM(l,m)];
	    }
	  }
	ws += 1.0;
      }else{// weight by radius^2, assuming constant Delta z
	r=(REARTH - z)/REARTH;// radius
	if((r > 1) || (r < 0)){
	  fprintf(stderr,"mean_models: error, r (%g) out of range, given depth %g\n",r,z);
	  exit(-1);
	}
	w = SQUARE(r);// weight
	for(l=0;l <= int_model->lmax;l++)// add to mean
	  for(m=0;m <= l;m++){
	    if(int_model->vector_field){
	      out_model->amp[0][POSLM(l,m)] += int_model->amp[0][POSLM(l,m)] * w;
	      out_model->bmp[0][POSLM(l,m)] += int_model->bmp[0][POSLM(l,m)] * w;
	      out_model->amt[0][POSLM(l,m)] += int_model->amt[0][POSLM(l,m)] * w;
	      out_model->bmt[0][POSLM(l,m)] += int_model->bmt[0][POSLM(l,m)] * w;
	    }else{
	      out_model->a[0][POSLM(l,m)] += int_model->a[0][POSLM(l,m)] * w;
	      out_model->b[0][POSLM(l,m)] += int_model->b[0][POSLM(l,m)] * w;
	    }
	  }
	ws += w;
      }
    }
  }else{
    ws=0.0;
  }
  if(ws == 0.0){
    fprintf(stderr,"mean_models: error, no averaging, j: %i ws: %g\n",j,ws);
    exit(-1);
  }
  for(l=0;l <= int_model->lmax;l++){
    for(m=0;m <= l;m++){
      if(int_model->vector_field){
	out_model->amp[0][POSLM(l,m)] /= ws;
	out_model->bmp[0][POSLM(l,m)] /= ws;
	out_model->amt[0][POSLM(l,m)] /= ws;
	out_model->bmt[0][POSLM(l,m)] /= ws;
      }else{
	out_model->a[0][POSLM(l,m)] /= ws;
	out_model->b[0][POSLM(l,m)] /= ws;
      }
    }
  }
  deallocate_model_coefficients(int_model);
}

COMP_PRECISION my_make_nan(void)
{
  static int been_here = 0;
  static COMP_PRECISION val;
  if(!been_here){
    //GMT_make_fnan(val);
    val = sqrt(-1.0);
    been_here = 1;
  }
  return val;
  
}
