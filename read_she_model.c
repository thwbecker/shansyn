/* 

   part of the shansyn spherical harmonics package, see COPYRIGHT for license 

   $Id: read_she_model.c,v 1.10 2009/04/16 00:55:29 becker Exp becker $ 

*/
#include "shansyn.h"


/*
  read in spherical harmonic expansion model with depth layers
  from file

  the first integer, the layer number n, determines the radial type of
  model. 

  0<n<1000: discrete layer, will need to read depths
  n < 0   : Chebyshev expansion of order -n
  1000<=n : splines interpolation of order n-1000
  
  fname: filename
  lmax: 
         if the input value lmax is set to <=0, no change in l_max
         if > 0, the model resolution will be changed to l_max,
         adding zeroes if necessary

  nexp:  number of expansions given in a column (1,2,3..; normally 1)

  expect_gsh: 1: expect generalized spherical harmonics 0: regular format

*/
void read_she_model(char *fname, /* filename */
		    struct mod *model,
		    int lmax,int nexp,
		    int expect_gsh)
{
  FILE *in;
  int j,l,m,i,iread;
  model->dmin =  FLT_MAX;
  model->dmax = -FLT_MAX;
  in=fopen(fname,"r");
  if(!in){
    fprintf(stderr,"read_she_model: error opening spherical harmonics model file \"%s\"\n",
	    fname);
    exit(-1);
  }
  if(fscanf(in,"%i",&model[0].n)!=1){ /* read number of layers/radial type */
    fprintf(stderr,"read_she_model: read error model file %s, could not read nl\n",fname);
    exit(-1);
  }
  if(model[0].n < 0){
    model[0].n *= -1;
    model[0].radial_type=CHEBYSHEV;
    fprintf(stderr,"WARNING: assuming Chebychev radial parameterization\n");
  }else{
    if(model[0].n > 1000){
      fprintf(stderr,"WARNING: assuming splines, if you have more than 1000 linear layers, modified code!\n");
      model[0].radial_type = SPLINES;
      model[0].n -= 1000.0;
    }else
      model[0].radial_type = DISCRETE;
  }
  if(model[0].n < 1){
    fprintf(stderr,"read_she_model: read error model file %s, number of layers is %i\n",
	    fname,model[0].n);
    exit(-1);
  }
  model[0].lmax = 0;
  for(i=1;i < nexp;i++)
    copy_model_par(model,(model+i));

  for(i=0;i < nexp;i++)
    if(!(model[i].d = (COMP_PRECISION *)
	 malloc(sizeof(COMP_PRECISION)*model[i].n)))
    MEMERROR;
  for(j=0;j < model[0].n;j++){
    /* 
       j -- depth / layer loop
    */
    if(model[0].radial_type == DISCRETE){// for discrete model, need to read in depth 
      if(fscanf(in,DATA_FSCAN_FORMAT,&model[0].d[j])!=1){
	fprintf(stderr,"read_she_model: read error model file %s, layer %i, could not read d, nl=%i, lmax=%i\n",
		fname,j,model[0].n,model[0].lmax);exit(-1);}
      if(model[0].d[j] < 0){
	fprintf(stderr,"read_she_model: depth argument has to be positive, %g\n",
		model[0].d[j]);
	exit(-1);
      }
      if(j && (model[0].d[j] > model[0].d[j-1])){
	fprintf(stderr,
		"read_she_model: model %s\nread_she_model: discrete layers have to sorted be from bottom to top,\nread_she_model: largest z (positive) first, %i: %g %i: %g\n",
		fname,j,model[0].d[j],j-1,model[0].d[j-1]);
	exit(-1);
      }
      for(i=1;i<nexp;i++)
	model[i].d[j] = model[0].d[j];
    }else{// Chebyshev or splines model, read dmin and dmax only once
      if(j == 0){
	if(fscanf(in,TWO_DATA_FSCAN_FORMAT,&model[0].dmin,
		  &model[0].dmax)!=2){
	  fprintf(stderr,"read_she_model: read error dmin / dmax Chebyshev model\n");
	}
	for(i=1;i<nexp;i++){
	  model[i].dmin = model[0].dmin;
	  model[i].dmax = model[0].dmax;
	}
      }
    }
    if(model[0].radial_type == DISCRETE){
      if(model[0].d[j] > model[0].dmax)
	model[0].dmax = model[0].d[j];
      if(model[0].d[j] < model[0].dmin)
	model[0].dmin = model[0].d[j];
      for(i=1;i<nexp;i++){
	model[i].dmin = model[0].dmin;
	model[i].dmax = model[0].dmax;
      }
    }
    //fprintf(stderr,"%i/%i at %g, exp  %i, expect gsh %i\n",j+1,model[0].n,model[0].d[j],nexp,expect_gsh);
    if(fscanf(in,"%i",&model[0].lmax) != 1){
      fprintf(stderr,"read_she_model: read error file %s, layer %i could not read lmax\n",
	      fname,j);exit(-1);}
    for(i=1;i < nexp;i++){
      model[i].lmax = model[0].lmax;
    }
    /* 
       determine GSH type 
    */
    if(expect_gsh){
      if(fscanf(in,"%i",&model[0].is_gsh) != 1){
	fprintf(stderr,"read_she_model: read error file %s, layer %i could not read GSH type\n",
		fname,j);exit(-1);}
      if((model[0].is_gsh != 0)&&(model[0].is_gsh != 2)&&(model[0].is_gsh != 4)){
	fprintf(stderr,"read_she_model: GSH type on input has to be 0, 2, or 4, but is %i\n",model[0].is_gsh);
	exit(-1);
      }
      model[0].is_gsh += 1;	/* internal storage is 0, 1, 3, 5 */
      for(i=1;i < nexp;i++)
	model[i].is_gsh = model[0].is_gsh;
    }else{
      for(i=0;i < nexp;i++)
	model[i].is_gsh = 0;
    }
    if(j==0){
      for(i=0;i < nexp;i++){
	// first layer, allocate pointers for coefficients
	// and determine offesets and the like
	allocate_model_coefficients((model+i));
	//fprintf(stderr,"allocating %i is_gsh: %i lms %i\n",i,model[i].is_gsh,model[i].lms);
      }
    }
    if(!model[0].is_gsh){
      /* regular */
      for(l=0;l <= model[0].lmax;l++){
	for(m=0;m <= l;m++){
	  for(i=0;i < nexp;i++){
	    if((fscanf(in,TWO_DATA_FSCAN_FORMAT,
		       (model[i].a[j]+POSLM(l, m)),
		       (model[i].b[j]+POSLM(l, m)))) != 2){
	      fprintf(stderr,"read_she_model: read error file %s, layer %i, l %i, m %i, model %i out of %i\n",
		      fname,j,l,m,i,nexp);exit(-1);
	    }
	  }
	}
      }
    }else{
      /* GSH */
      iread = 0;
      //fprintf(stderr,"GSH vec %i izero %i is_gsh %i lmax %i lms %i\n",model[0].vector_field, model[0].izero, model[0].is_gsh,model[0].lmax, model[0].lms);
      
      for(l=model[0].izero;l <= model[0].lmax;l++){
	for(m=0;m <= l;m++){
	  for(i=0;i < nexp;i++){
	    read_gsh_coeff(in,
			   (model[i].is_gsh)-1, l, m, 
			   model[i].a[j],   model[i].b[j], 
			   model[i].amp[j], model[i].amt[j], 
			   model[i].bmp[j], model[i].bmt[j],
			   &iread);
	  }
	}
      }
      for(i=0;i < nexp;i++){
	if(iread != model[i].lms * nexp){
	  fprintf(stderr,"gsh read error, expected %i, read %i non-zero coefficients\n",
		  model[i].lms,iread);
	  exit(-1);
	}
      }

    }
  }
  fclose(in);
  for(i=0;i < nexp;i++){
    // modify l_max
    if(lmax > 0){
      fprintf(stderr,"modifying lmax to %i\n",lmax);
      select_lms(model[i].is_gsh,lmax,
		 &(model[i].lms),&(model[i].vector_field),&(model[i].izero));
      /* shrink / expand */
      for(j=0;j < model[i].n;j++){
	model[i].a[j]=(COMP_PRECISION *)
	  realloc(model[i].a[j],
		  sizeof(COMP_PRECISION)*model[i].lms*model[i].n);
	model[i].b[j]=(COMP_PRECISION *)
	  realloc(model[i].b[j],
		  sizeof(COMP_PRECISION)*model[i].lms*model[i].n);
	if(!model[i].a[j] || !model[i].b[j])MEMERROR;
	/* fill in with zeroes */
	if(lmax > model[i].lmax){
	  if(!model[i].is_gsh){
	    for(l=model[i].lmax+1;l<=lmax;l++)
	      for(m=0;m<=l;m++)
		*(model[i].a[j]+POSLM(l, m))= 
		  *(model[i].b[j]+POSLM(l, m))=0.0;
	  }else{
	    fprintf(stderr,"ERROR: GSH lmax modification not implemented yet\n");
	    exit(-1);
	  }
	}
	if((j==0) && (i==0))
	  fprintf(stderr,"read_she_model: changed l_max from %i to %i\n",model[i].lmax,lmax);
      }
      model[i].lmax=lmax;
    }
    if(i==0)
      switch(model[i].radial_type){
      case DISCRETE:{
	fprintf(stderr,"read_she_model: discrete layer model %s, %i layers, L= %i, dmin/max: %g %g GSH: %i vector: %i\n",
		fname,model[i].n,model[i].lmax,model[i].dmin,model[i].dmax,model[i].is_gsh,model[i].vector_field);
	break;
      }
      case CHEBYSHEV:{
	fprintf(stderr,"read_she_model: Chebyshev model %s, order %i, L= %i, dmin/max: %g %g GSH: %i vector: %i\n",
	      fname,model[i].n,model[i].lmax,model[i].dmin,model[i].dmax,model[i].is_gsh,model[i].vector_field);
	break;
      }
      case SPLINES:{
	fprintf(stderr,"read_she_model: spline model %s, order %i, L= %i, dmin/max: %g %g GSH: %i vector: %i\n",
		fname,model[i].n,model[i].lmax,model[i].dmin,model[i].dmax,model[i].is_gsh,model[i].vector_field);
	break;
      }
      default:{
	fprintf(stderr,"read_she_model: ended up with wrong model code: %i\n",
		model[i].radial_type);
	exit(-1);
      }}
  }
}

/* 
   depending on the type of model, allocate scalar or vector field
   coefficients

   for this, lmax, is_gsh, n have to be defined
 */
void allocate_model_coefficients(struct mod *model)
{
  int j;
  /* detect how much room is needed */
  select_lms(model->is_gsh,model->lmax,
	     &(model->lms),&(model->vector_field),
	     &(model->izero));
  if(model->lms == 0){
    fprintf(stderr,"allocate_model_coefficients: lms is %i\n",
	    model->lms);
    exit(-1);
  }
  if(model->vector_field){
    /* vector */
    model->amp = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    model->bmp = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    model->amt = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    model->bmt = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    if(!model->amp || !model->bmp || 
       !model->amt || !model->bmt)MEMERROR;

    /* need to have empty ones */
    model->a = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
    model->b = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
 
  }else{
    /* scalar */
    model->a = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    model->b = (COMP_PRECISION **)
      malloc(sizeof(COMP_PRECISION *)*model->n);
    if(!model->a || !model->b)MEMERROR;
    
    /* empty */
    model->amp = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
    model->bmp = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
    model->amt = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
    model->bmt = (COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *));
 
  }

  /* layers */
  for(j=0;j < model->n;j++){
    if(model->vector_field){
      //fprintf(stderr,"allocate %i %i %i\n",j,model->n,model->lms);
      model->amp[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      model->bmp[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      if(!model->amp[j] || !model->bmp[j])MEMERROR;
      
      model->amt[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      model->bmt[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      if(!model->amt[j] || !model->bmt[j])MEMERROR;
    }else{
      model->a[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      model->b[j]=(COMP_PRECISION *)
	malloc(sizeof(COMP_PRECISION)*model->lms*model->n);
      if(!model->a[j] || !model->b[j])MEMERROR;
    }
  }
}
void deallocate_model_coefficients(struct mod *model)
{
  int j;
  /* layers */
  for(j=0;j < model->n;j++){
    if(model->vector_field){
      free(model->amp[j]);
      free(model->bmp[j]);
      free(model->amt[j]);
      free(model->bmt[j]);
    }else{
      free(model->a[j]);
      free(model->b[j]);
    }
  }
  if(model->vector_field){
    /* vector */
    free(model->amp);
    free(model->bmp);
    free(model->amt);
    free(model->bmt);
  }else{
    free(model->a);
    free(model->b);
  }
}


/* copy from a to b */
void copy_model_par(struct mod *a,struct mod *b)
{
  b->n = a->n;
  b->lmax = a->lmax;
  b->lms = a->lms;
  b->vector_field = a->vector_field;
  b->radial_type = a->radial_type;
  b->is_gsh = a->is_gsh;
  b->izero = a->izero;
  b->dmin = a->dmin;
  b->dmax = a->dmax;
  
}
