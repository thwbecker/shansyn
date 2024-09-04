/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: abconvert.c,v 1.32 2009/04/16 00:55:29 becker Exp becker $ */
/* 
   program does tapering and format conversion for
   spherical harmonic coefficients 

   handles physical convention 
   and geodetic normalization

*/
#include "abconvert.h"

void write_nonzero_coefficients(COMP_PRECISION **,COMP_PRECISION **,int,COMP_PRECISION,int ,FILE *, int ,int );


int main(int argc, char *argv[] )
{
  int lmax,lmax1,l,m,tapering,out_format,lmsize,in_format,j,lstart,
    md,llmax,layers=0,i,verbose=VERBOSE,highp = -1,rc,
    vector_field=0,lmin = 0,nset,shps,
    os1,ialpha=-1,icoeff=0,izero=0,iread;
  int cmode = 0;
  BOOLEAN normalize_power_by_ncoeff = TRUE, /* default is power per
					       degree and unit area */
    interpolate_mode = FALSE,
    use_correlation=FALSE,
    use_admittance=FALSE,
    use_second = FALSE;
  long seed;
  COMP_PRECISION *a=NULL,*b=NULL,tmp,amplitudefactor,
    *filter=NULL,rf1,rf2,dp,atmp,btmp,*amt=NULL,*amp=NULL,
    *bmp=NULL,*bmt=NULL,*c=NULL,*d=NULL,*cmt=NULL,*cmp=NULL,
    *dmp=NULL,*dmt=NULL,lc,weight1,weight2,*af=NULL,*bf=NULL,
    fac1,fac2,fac3,tfac[2];
  FILE *in=NULL;
  /* help page?> */
  if(argc > 1 && strcmp(argv[1],"-h")==0 ){
    phelp(argv[0]);exit(-1);
  }
  /* 
     default values 
  */
  in_format=IN_FORMAT_DEFAULT;
  out_format=OUT_FORMAT_DEFAULT;
  tapering=TAPERING_DEFAULT;
  amplitudefactor=AMPLITUDE_DEFAULT;
  lc=LC_DEFAULT;
  md = -1;
  llmax=99999;
  /* deal with input values */
  if(argc > 1)
    sscanf(argv[1],"%i",&out_format);
  if(argc>2)
    sscanf(argv[2],"%i",&tapering);
  if(argc>3)
    sscanf(argv[3],DATA_FSCAN_FORMAT,&amplitudefactor);
  if(argc>4)
    sscanf(argv[4],"%i",&in_format);
  if(argc>5){
    sscanf(argv[5],"%i",&llmax);
    if(llmax == 0){		/* use all coefficients */
      llmax = 99999;
    }else if(llmax < 0){	/* high pass */
      if(tapering != 0)  
	fprintf(stderr,"%s: overriding taper %i with high pass\n",
		argv[0],tapering);
      tapering = HIGHP_TAPER;
      highp = - llmax;
      llmax = 99999;
    }
  }
  if(argc > 6){
    if(tapering == HIGHP_TAPER){
      fprintf(stderr,"%s: WARNING: overriding lc setting\n",argv[0]);
    }else{
      sscanf(argv[6],DATA_FSCAN_FORMAT,&lc);
      if((lc >= 1.0)&&(tapering != EXP_TAPER)){
	fprintf(stderr,"%s: WARNING: lc >= 1: switching tapering off. is that what you wanted?\n",
		argv[0]);
	lc=0.0;
	tapering=NO_TAPER;
      }
    }
  }
  if(argc>7){
    sscanf(argv[7],"%i",&md);		/* debugging m level */
  }
  if(argc>8){
    phelp(argv[0]);exit(-1);
  }
  /* 

     done input values

  */
  if((in_format == INTERPOLATE) || (in_format == INTERPOLATE_ABAB) || (in_format == INTERPOLATE_GSH))
    interpolate_mode = TRUE;
  else
    interpolate_mode = FALSE;
  
 
  /* 
     decide if we need correlations, and which ones, or if we compute admittance 

  */
  switch(out_format){
  case CORRL_OUT:
  case CORRTH_OUT:
  case CORRT_OUT:		/* linear correlation */
    use_correlation = TRUE;
    cmode = 1;
    break;
  case SPEAR_CORRL_OUT:
  case SPEAR_CORRTH_OUT:
  case SPEAR_CORRT_OUT:		/* pearson's correlation */
    use_correlation = TRUE;
    cmode = 2;
    break;
  case ADMITTANCE_OUT:
  case ADMITTANCE_OUT_NN:
    use_admittance = TRUE;
    break;
  default:
    break;
  }

  /*
    
    DATA INPUT

    
    read in set(s) of coefficients and convert to internal 
    convention 

    start of input loop 
    
  */
  while(fscanf(stdin,"%i",&lmax) == 1){
    /* 

    
       all formats need to start with an integer specifying the degree of the expansion


    */
    lmax1 = lmax+1;
    switch(in_format){// start input format case structure
    case INTERPOLATE:
    case AB_INPUT:
    case AB_INPUT_HC:
    case AB_NONZERO_INPUT:
      /* 
	 regular SH format 
      */
      if(verbose)
	fprintf(stderr,"%s: reading A B format in physical convention%s lmax: %i\n",
		argv[0],(in_format == AB_INPUT_HC)?(", long hc format,"):(","),lmax);
      if(in_format == AB_INPUT_HC){
	rc=fscanf(stdin,"%*i %*f %i %i %*i",&nset,&shps);
	if((nset != 1)||(shps!=1)){
	  fprintf(stderr,"%s: error with long (hc) format: nset: %i shps: %i\n",
		  argv[0],nset,shps);
	  exit(-1);
	}
      }
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",
		argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);
      }
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if(in_format == AB_NONZERO_INPUT){
	    /* only non-zero terms are listed */
	    if(m==0){		/* only read A */
	      if((fscanf(stdin,DATA_FSCAN_FORMAT,(a+os1)))!=1){
		fprintf(stderr,"%s: read error, l=%i m=%i\n\n",
			argv[0],l,m);exit(-1);}
	      b[os1] = 0.0;
	    }else{
	      if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+os1),(b+os1)))!=2){
		fprintf(stderr,"%s: read error, l=%i m=%i\n\n",
			argv[0],l,m);exit(-1);}
	    }
	  }else{
	    /* regular input */
	    if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+os1),(b+os1)))!=2){
	      fprintf(stderr,"%s: read error, l=%i m=%i\n\n",
		      argv[0],l,m);exit(-1);}
	  }
	}
      break;
    case LAB_GEOD_INPUT:// layer lm A B (different from lm A B geodetic!!
      /* plates/rick_expansions/

	 geodetic normalization l m format 
      
      */

      if(verbose)
	fprintf(stderr,"%s: reading l m A B LAYER FORMAT in GEODETIC convention\n",argv[0]);
      
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",
		argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      if(layers!=1)
	if(verbose)
	  fprintf(stderr,"%s: data file has %i layers, reading only layer 1\n",
		  argv[0],layers);
      for(rc=m=0;m<layers;m++) {
	rc+=fscanf(stdin,"%*g");
      }
      rc=fscanf(stdin,"%*i %*g");
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1=POSLM(l, m);
	  if((fscanf(stdin,TWO_TWO_DATA_FSCAN_FORMAT,
		     (a+os1),(b+os1)))!=2){
	    fprintf(stderr,"%s: read error, l=%i m=%i\n\n",
		    argv[0],l,m);exit(-1);}
	  *(a+os1) *= (GEODETIC_FACTOR(l,m));
	  *(b+os1) *= (GEODETIC_FACTOR(l,m));
	}
      break;
    case AB_GEOD_INPUT:{
      /* 
	 
	 geodetic convention A B format
      
      */
      if(verbose)
	fprintf(stderr,"%s: reading A B format in geodetic convention lmax=%i\n",
		argv[0],lmax);
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+os1),(b+os1)))!=2){
	    fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);exit(-1);}
	  *(a+os1) *= (GEODETIC_FACTOR(l,m));
	  *(b+os1) *= (GEODETIC_FACTOR(l,m));
	}
      break;
    }
    case AB_RICK_INPUT:{
      /* 

	 "fully normalized", A B format

      */
      if(verbose)
	fprintf(stderr,"%s: reading A B format in Rick's convention, lmax: %i\n",argv[0],lmax);
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+os1),(b+os1)))!=2){
	    fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);exit(-1);}
	  *(a+os1) *= (RICK_SCALAR_FACTOR(l,m));
	  *(b+os1) *= (RICK_SCALAR_FACTOR(l,m));
	}
      break;
    }
    case AB_MASTERS_INPUT:{
      /* 

	 Masters/Edmonds nomalization AB format

      */
      if(verbose)
	fprintf(stderr,"%s: reading A B format in G. Master's (Edmonds, 1960) convention \n",argv[0]);
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+os1),(b+os1)))!=2){
	    fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);exit(-1);}
	  *(a+os1) *=  SQRT_TWO; 
	  *(b+os1) *= -SQRT_TWO;
	}
      break;
    }
    case INTERPOLATE_ABAB:// velocity ABAB format input
    case ABAB_INPUT:{
      /* 

	 poloidal/toroidal expansion in Dahlen & Tromp format

      */
      if(verbose)
	fprintf(stderr,
		"%s: reading poltor vector AB AB format in physical convention, lmax: %i\n",
		argv[0],lmax);
      /* allocate */
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((amt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      if((amp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      /* input */
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){ 
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,FOUR_DATA_FSCAN_FORMAT,
		     (amp+os1),(bmp+os1),(amt+os1),(bmt+os1)))!=4){
	    fprintf(stderr,"%s: read error, A l=%i m=%i\n\n",argv[0],l,m);
	    exit(-1);
	  }
	}
      vector_field = 1;		/* pol tor */
      break;
    }
    case ABAB_GEOD_INPUT:{
      /* 

	 poloidal/toroidal input in geodetic normalization

      */
      if(verbose)
	fprintf(stderr,"%s: reading poltor vector AB AB format in geodetic convention\n",argv[0]);
      
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((amt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      if((amp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){	  
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,FOUR_DATA_FSCAN_FORMAT,
		     (amp+os1),(bmp+os1),(amt+os1),
		     (bmt+os1)))!=4){
	    fprintf(stderr,"%s: read error, A l=%i m=%i\n\n",argv[0],l,m);exit(-1);}
	  *(amp+os1) *= GEODETIC_FACTOR(l, m);
	  *(bmp+os1) *= GEODETIC_FACTOR(l, m);
	  *(amt+os1) *= GEODETIC_FACTOR(l, m);
	  *(bmt+os1) *= GEODETIC_FACTOR(l, m);
	}
      vector_field = 1;
      break;
    }
    case AABBR_INPUT:{
      /* 


	 poloidal/toroidal input in rick's "fully normalized" convention

      */
      if(verbose)
	fprintf(stderr,"%s: reading poltor vector AA\\nBB format in Rick's convention\n",argv[0]);
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((amt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      if((amp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if(fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(amp+os1),(amt+os1))!=2){
	    if(verbose)fprintf(stderr,"%s: read error, A l=%i m=%i\n\n",
			       argv[0],l,m);exit(-1);}
	  if(fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(bmp+os1),(bmt+os1))!=2){
	    if(verbose)fprintf(stderr,"%s: read error, B l=%i m=%i\n\n",
			       argv[0],l,m);exit(-1);}
	  *(amp+os1) *= RICK_FACTOR(l, m);
	  *(bmp+os1) *= RICK_FACTOR(l, m);
	  *(amt+os1) *= RICK_FACTOR(l, m);
	  *(bmt+os1) *= RICK_FACTOR(l, m);
	}
      vector_field =  1;
      break;
    }
    case LMAB_GEOD_INPUT:{
      /* 

	 geodetic format, lm A B format

      */
      if(verbose)
	fprintf(stderr,"%s: reading l m A B format in geodetic convention \n",argv[0]);
      
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,TWO_TWO_DATA_FSCAN_FORMAT,
		     (a+os1),(b+os1)))!=2){
	    if(verbose)
	      fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);
	    exit(-1);
	  }
	  *(a+os1) *= GEODETIC_FACTOR(l, m);
	  *(b+os1) *= GEODETIC_FACTOR(l, m);
	}
      break;
    }
    case LMAB_FNORM_INPUT:{
      /* 


	 poloidal/toroidal Edmonds convention

      */
      if(verbose)
	fprintf(stderr,"%s: reading l m A B in fully normalized (Edmonds, 1960) format\n",argv[0]);
      
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      if((a=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (b=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(stdin,TWO_TWO_DATA_FSCAN_FORMAT,
		     (a+os1),(b+os1)))!=2){
	    if(verbose)fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);
	    exit(-1);
	  }
	}
      // go from "fully normalized" (Edmonds or old Harvard)  to Rick's convention
      // (this from advect code....)
      fac1 = 1.0/sqrt(4.0*PI);
      fac2 = fac1/sqrt(2.0);
      for(l=0;l<=lmax;l++){
	a[POSLM(l,0)] *= fac1;
	if(l>0){
	  for(m=1;m<=l;m++){
	    os1 = POSLM(l,m);
	    fac3 = fac2 * pow(-1.0,(COMP_PRECISION)(m));
	    a[os1] *= fac3;
	    b[os1] *= fac3;
	  }
	}
      }
      // go from Rick's convention to our physical norm
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  a[os1] *= RICK_SCALAR_FACTOR(l,m);
	  b[os1] *= RICK_SCALAR_FACTOR(l,m);
	}
      break;
    }
    case INTERPOLATE_GSH:
    case GSH_INPUT:{
      /* 
	 
	 generalized spherical harmonics input, those can be scalar,
	 2phi, or 4phi
	 

      */
      rc=fscanf(stdin,"%i",&ialpha); /* read type flag
				     0: scalar, 2: twophi, 4: four phi
				  */
      fprintf(stderr,"%s: reading lmax = %i GSH, BW normalization, %i type ... ",
	      argv[0],lmax,ialpha);
      /* make room in poloidal/toroidal storage */
      lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
      /* 
	 
	 compute non-zero coefficient number 
	 
      */
      select_lms((ialpha+1),lmax,&icoeff,&vector_field,&izero);
      if(lmax < izero){
	fprintf(stderr,"\n%s: gsh ialpha %i requires lmax > %i, lmax = %i\n",
		argv[0],ialpha,izero,lmax);
	exit(-1);
      }
      /* make room for expansion */
      if(!vector_field){
	/* scalar */
	if((a=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (b=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  if(verbose)fprintf(stderr,"\n%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	amt = bmt = amp = bmp = NULL;
      }else{
	/* real and imag parts */
	if((amp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (bmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  if(verbose)fprintf(stderr,"\n%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	if((amt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (bmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  if(verbose)fprintf(stderr,"\n%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	a = NULL; b = NULL;
      }
      /* we don't need to zero out the first coefficients, as those are 
	 set to zero by calloc */
      /* read in non-zero coefficients */
      iread = 0;
      for(l=izero;l <= lmax;l++)
	for(m=0;m <= l;m++)	/* read GSH format and rescale */
	  read_gsh_coeff(stdin,ialpha, l, m, a, b, 
			 amp, amt, bmp, bmt,&iread);
      if(iread != icoeff){
	fprintf(stderr,"\n%s: gsh read error, expected %i, read %i non-zero coefficients\n",
		argv[0],icoeff,iread);
	exit(-1);
      }
      fprintf(stderr,"done\n");
      break;
    }
    default:{
      fprintf(stderr,"%s: input format %i is undefined\n",
	      argv[0],in_format);
      phelp(argv[0]);exit(-1);
      break;
    }
    }
    /* 

       done with first file/expansion input. now check for second expansion
    

    */
    /* 

       check if we want correlations or interpolations 
       
       in this case, will need another file
    

    */
    if(use_correlation || interpolate_mode || use_admittance){
      /* 

	 input or output format needs a second file

      */
      use_second = TRUE;
      if(verbose){
	if(interpolate_mode){
	  switch(in_format){
	  case INTERPOLATE:
	    fprintf(stderr,"%s: second AB scalar file for interpolation, lmax=%i\n",argv[0],lmax);
	    break;
	  case INTERPOLATE_ABAB:
	    fprintf(stderr,"%s: second ABAB vector field file for interpolation, lmax=%i\n",argv[0],lmax);
	    break;
	  case INTERPOLATE_GSH:
	    fprintf(stderr,"%s: second GSH file for interpolation, lmax=%i\n",argv[0],lmax);
	    break;
	  default:
	    fprintf(stderr,"%s: internal error, interpolate mode not defined\n",argv[0]);
	    exit(-1);
	  }
	}else{
	  if(use_admittance)
	    fprintf(stderr,"%s: second file for admittance, lmax=%i\n",
		    argv[0],lmax);

	  else
	    fprintf(stderr,"%s: second file for correlation coefficient, lmax=%i\n",
		    argv[0],lmax);
	}
      }
      if(fscanf(stdin,"%i",&i) != 1){ /* second lmax */
	fprintf(stderr,"%s: read error second file lmax for mode %i \n",
		argv[0],out_format);
	exit(-1);
      }
      if(in_format == AB_INPUT_HC){
	rc = fscanf(stdin,"%*i %*f %i %i %*i",&nset,&shps);
	if((nset != 1)||(shps!=1)){
	  fprintf(stderr,"%s: error with long (hc) format: nset: %i shps: %i\n",
		  argv[0],nset,shps);
	  exit(-1);
	}
      }
      if(i != lmax){
	if(verbose)fprintf(stderr,"%s: second lmax=%i not equal to first file lmax=%i\n\n",argv[0],i,lmax);
	exit(-1);
      }
      if((in_format == GSH_INPUT)||(in_format == INTERPOLATE_GSH)){ 
	/* gsh format */
	rc = fscanf(stdin,"%i",&j);
	if(j != ialpha){
	  if(verbose)fprintf(stderr,"%s: gsh second ialpha=%i (0/2/4 types) not equal to first file ialphs=%i\n\n",
			     argv[0],j,ialpha);
	  exit(-1);
	}
      }
      /* allocate */
      if(!vector_field){
	//
	// read in second scalar expansion
	//
	/* this is kinda bad and should be standardized with first
	   expansion read */
	if((c=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (d=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	cmt = dmt = cmp = dmp = NULL;
	switch(in_format){
	case INTERPOLATE:
	case AB_INPUT:		/* scalar, physical convention */
	case AB_INPUT_HC:		/* scalar, physical convention */
	  for(l=0;l<=lmax;l++)
	    for(m=0;m<=l;m++){
	      os1 = POSLM(l, m);
	      if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(c+os1),(d+os1)))!=2){
		if(verbose)fprintf(stderr,"%s: second file read error, l=%i m=%i\n\n",argv[0],l,m);
		exit(-1);
	      }
	    }
	  break;
	case GSH_INPUT:		/* scalar, gsh */
	case INTERPOLATE_GSH:		/* scalar, gsh */
	  iread = 0;
	  for(l=izero;l<=lmax;l++)
	    for(m=0;m<=l;m++)
	      read_gsh_coeff(stdin,
			     ialpha, l, m, 
			     c, d, cmp, cmt, dmp, dmt, &iread);
	  if(iread != icoeff){
	    fprintf(stderr,"\n%s: gsh read error, expected %i, read %i non-zero coefficients\n",
		    argv[0],icoeff,iread);
	    exit(-1);
	  }
	  fprintf(stderr,"%s: second GSH file done\n",argv[0]);
	  break;
	default:
	  fprintf(stderr,"%s: expecting second scalar for outmode %i but input type %i undefined in second read\n",
		  argv[0],out_format,in_format);
	  exit(-1);
	  break;
	}
      }else{
	/* 
	   vector field

	   poloidal/toroidal coefficients or GSH 2phi/4phi terms
	

	*/
	if((cmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (dmt=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  if(verbose)fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	if((cmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	   (dmp=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	  fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",argv[0],lmax,lmsize); 
	  phelp(argv[0]);exit(-1);}
	c = d = NULL;
	switch(in_format){
	case INTERPOLATE_ABAB:
	case ABAB_INPUT:		/* vector, physical convection*/
	  for(l=0;l<=lmax;l++)
	    for(m=0;m<=l;m++){
	      os1 = POSLM(l, m);
	      if((fscanf(stdin,FOUR_DATA_FSCAN_FORMAT,
			 (cmp+os1),(dmp+os1),
			 (cmt+os1),(dmt+os1)))!=4){
		fprintf(stderr,"%s: read error, second ABAB file l=%i m=%i\n\n",argv[0],l,m);
		exit(-1);
	      }
	    }
	  break;
	case INTERPOLATE_GSH:
	case GSH_INPUT:
	  /* read in non-zero coefficients */
	  iread = 0;
	  for(l=izero;l<=lmax;l++)
	    for(m=0;m<=l;m++)
	      read_gsh_coeff(stdin,ialpha,l,m,
			     c,d,cmp,cmt,dmp,dmt,&iread);
	  if(iread != icoeff){
	    fprintf(stderr,"\n%s: gsh second file read error, expected %i, read %i non-zero coefficients\n",
		    argv[0],icoeff,iread);
	    exit(-1);
	  }
	  break; 
	default:
	  fprintf(stderr,"%s: expecting second AB AB velocity/GSH file for outmode %i but input type %i undefined in second read\n",
		  argv[0],out_format,in_format);
	  exit(-1);
	  break;
	}
      }	/* end velocity type */
    }// end second expansion file input loop
    //
    // 
    //
    /*
      open filter file, if required

    */
    if(tapering == FROM_FILE_TAPER){
      lc=0.0;
      fprintf(stderr,"%s: resetting l_c to zero for filtering with %i weights read from %s\n",
	      argv[0],lmax1,FILTER_FILE);
      in=fopen(FILTER_FILE,"r");
      if(!in){
	fprintf(stderr,"%s: could not read %s for filter weights\n",
		argv[0],FILTER_FILE);
	fprintf(stderr,"%s: taper option %i was set\n",
		argv[0],tapering);
	exit(-1);
      }
      //
      // read in filter weights
      //
      filter=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(lmax1));
      if(!filter)
	MEMERROR;
      for(i=0;i <= lmax;i++)
	if(fscanf(in,DATA_FSCAN_FORMAT,(filter+i))!=1){
	  fprintf(stderr,"%s: read error, weight l=%i of %s, expecting %i\n",
		  argv[0],i,FILTER_FILE,lmax);
	  exit(-1);
	}
      fclose(in);
      if(verbose)
	fprintf(stderr,"%s: succesfully read weights\n",argv[0]);
    }else if(tapering == FROM_SH_FILE_TAPER){
      /* 
	 spherical harmonics from file scaling 
      */
      
      lc=0.0;
      fprintf(stderr,"%s: resetting l_c to zero for filtering with SH file %s, lmax = %i\n",
	      argv[0],FILTER_SH_FILE,lmax);
      in=fopen(FILTER_SH_FILE,"r");
      if(!in){
	fprintf(stderr,"%s: could not read %s\n",
		argv[0],FILTER_SH_FILE);
	fprintf(stderr,"%s: taper option %i was set\n",
		argv[0],tapering);
	exit(-1);
      }
      /* read in a spherical harmonics file for tapering of l,m
	 coefficients (this might screw up normalizations) */
      if((af=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
	 (bf=(COMP_PRECISION *)
	  calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
	fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",
		argv[0],lmax,lmsize); 
	phelp(argv[0]);exit(-1);
      }
      if(fscanf(in,"%i",&l)!=1){ /* lmax */
	fprintf(stderr,"%s: read error SH weights file\n",argv[0]);exit(-1);}
      if(l < lmax){
	fprintf(stderr,"%s: need lmax %i but weights SH file is only lmax %i\n",
		argv[0],lmax,l);
	exit(-1);
      }
      /* read in weights */
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  if((fscanf(in,TWO_DATA_FSCAN_FORMAT,
		     (af+os1),(bf+os1)))!=2){
	    fprintf(stderr,"%s: read error SH weights, l=%i m=%i\n\n",
		    argv[0],l,m);exit(-1);}
	}
      fclose(in);
      if(verbose)
	fprintf(stderr,"%s: succesfully read SH file for lmax = %i weights\n",
		argv[0],lmax);
    }
    /* 
       

       APPLY TAPER
     

    */
    if(verbose)
      switch(tapering){
      case NO_TAPER:{/* no taper */
	break;
      }
      case COSSQR_TAPER:{
	fprintf(stderr,"%s: applying cos^2 taper  from l' > %g\n",argv[0],lc);
	break;
      }
      case EXP_TAPER:{
	fprintf(stderr,"%s: applying exp taper with width %g\n",argv[0],lc);
	if(lc == 1)
	  fprintf(stderr,"%s: WARNING: width %g, is that what you wanted?\n",argv[0],lc);
	if(lc == 0){
	  fprintf(stderr,"%s: ERROR: set lc to something but zero for exp taper\n",argv[0]);
	  exit(-1);
	}
	break;
      }
      case COSP4_TAPER:{
	fprintf(stderr,"%s: applying cos^4 taper  from l' > %g\n",argv[0],lc);
	break;
      }
      case ONEML_TAPER:{
	fprintf(stderr,"%s: applying linear taper from l' > %g\n",argv[0],lc);
	break;
      }
      case ONEMLSQR_TAPER:{
	fprintf(stderr,"%s: applying square taper from l' > %g\n",argv[0],lc);
	break;
      }
      case SINOVRL_TAPER:{
	fprintf(stderr,"%s: applying Lanczos taper  from l' > %g\n",argv[0],lc);
	break;
      }
      case NNR_TAPER:{
	fprintf(stderr,"%s: *removing* net rotation part of pol/tor AB set\n",argv[0]);
	break;
      }
      case PASS_POL_TAPER:{
	fprintf(stderr,"%s: removing all toroidal terms of an ABAB set\n",argv[0]);
	break;
      }
      case PASS_TOR_TAPER:{
	fprintf(stderr,"%s: removing all poloidal terms of an ABAB set\n",argv[0]);
	break;
      }
      case NR_TAPER:{
	fprintf(stderr,"%s: passing *only net rotation* part of pol/tor AB set\n",argv[0]);
	break;
      }
      case ZERO_TAPER:{
	fprintf(stderr,"%s: setting all coefficients >  l' (%g = %i) to zero\n",argv[0],
		lc,(int)(lc*(float)lmax1));
	break;
      }
      case HIGHP_TAPER:{
	lc = MAX(0,(float)highp/(float)lmax1);
	fprintf(stderr,"%s: setting all coefficients <  %i (l' = %g) to zero\n",argv[0],
		(int)(lc*(float)lmax1),lc);
	break;
      }
      case L0_TAPER:{
	fprintf(stderr,"%s: removing l=0 term (mean=0)\n",argv[0]);
	break;
      }
      case ONLY_L_ZERO:{
	fprintf(stderr,"%s: removing all but l=0 term (RMS=0)\n",argv[0]);
	break;
      }
      case ONLY_M0_TERMS:{
	fprintf(stderr,"%s: only keeping m=0 terms\n",argv[0]);
	break;
      }
      case FROM_FILE_TAPER:{
	fprintf(stderr,"%s: applying filter from weights file\n",argv[0]);
	break;
      }
      case FROM_SH_FILE_TAPER:{
	fprintf(stderr,"%s: applying filter from SH file\n",argv[0]);
	break;
      }
      case LAPLACIAN:{
	fprintf(stderr,"%s: compute Laplacian\n",argv[0]);
	break;
      }
      case DIVERGENCE:{
	fprintf(stderr,"%s: compute divergence, switching to pol out\n",argv[0]);
	if(!vector_field){
	  fprintf(stderr,"%s: ERROR: need vector field input\n",argv[0]);
	  exit(-1);
	}
	out_format = VECABAB_POL_OUT;
	break;
      }
      case VORTICITY:{
       fprintf(stderr,"%s: compute |vorticity|, switching to tor out\n",argv[0]);
	if(!vector_field){
	  fprintf(stderr,"%s: ERROR: need vector field input\n",argv[0]);
	  exit(-1);
	}
	out_format = VECABAB_TOR_OUT;
	break;
      }
      case PHI_ROTATE:{
	fprintf(stderr,"%s: rotating field by %g degrees eastward\n",
		argv[0],amplitudefactor);
	break;
      }
      case SET_L_UNITY:{	/* debugging mode */
	if(md < 0){
	  seed = (long) md;
	  fprintf(stderr,"%s: WARNING: setting l <= %i terms to random (seed: %i), rest to zero\n",
		  argv[0],llmax,(int)seed);
	}else{
	  seed = -2;
	  fprintf(stderr,"%s: WARNING: setting l = %i m = %i terms to random (seed: %i), rest to zero\n",
		  argv[0],llmax,md,(int)seed);
	}
	break;
      }
      default:{
	fprintf(stderr,"%s: tapering mode %i is not defined\n",argv[0],tapering);
	phelp(argv[0]);
	exit(-1);
	break;
      }}				/* end tapering switch */
    if(interpolate_mode){
      /*
	
	input format needs interpolation between files
	in this case, amp factor is the weight

      */
      weight1 = amplitudefactor;
      weight2 = 1.0 - amplitudefactor;
      if(verbose)
	fprintf(stderr,"%s: first/second set are weighted by %g/%g\n",
		argv[0],weight1,weight2);
      if(weight1 > 1 || weight1 < 0){
	fprintf(stderr,"%s: first weight not between 0 and 1\n",
		argv[0]);
      }
      switch(in_format){
	
      case INTERPOLATE:// AB FORMAT
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    *(a+os1) *= weight1;
	    *(a+os1) += *(c+os1) * weight2;
	    *(b+os1) *= weight1;
	    *(b+os1) += *(d+os1) * weight2;
	  }
	break;
      case INTERPOLATE_GSH:
      case INTERPOLATE_ABAB:
	lstart = (in_format == INTERPOLATE_ABAB)?(0):(izero);
	for(l=lstart;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    *(amp+os1) *= weight1;
	    *(amp+os1) += *(cmp+os1) * weight2;
	    *(bmp+os1) *= weight1;
	    *(bmp+os1) += *(dmp+os1) * weight2;
	    *(amt+os1) *= weight1;
	    *(amt+os1) += *(cmt+os1) * weight2;
	    *(bmt+os1) *= weight1;
	    *(bmt+os1) += *(dmt+os1) * weight2;
	  }
	break;
      default:
	fprintf(stderr,"%s: internal error\n",argv[0]);
	break;
      }
      /* reset amplitude factor */
      amplitudefactor = 1.0;
    }else if(tapering == PHI_ROTATE){
      /* 
	 
	 rotation mode

      */
      if(in_format != AB_INPUT){
	fprintf(stderr,"%s: rotate taper only works with AB input\n",argv[0]);
	exit(-1);
      }
      /*
	
	if tapering is rotation, then amplitude factor is 
	the phi shift in degrees
	apply rotation here

      */
      dp=amplitudefactor*PIOVERONEEIGHTY;
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  tmp = (COMP_PRECISION)m*dp;
	  rf1=cos(tmp);
	  rf2=sin(tmp);
	  atmp=rf1* *(a+os1) - rf2 * *(b+os1);
	  btmp=rf2* *(a+os1) + rf1 * *(b+os1);
	  *(a+os1)=atmp;
	  *(b+os1)=btmp;
	}
      /* reset */
      amplitudefactor=1.0;
    }else{
      /* amplification/scaling?  */
      if(verbose && (amplitudefactor != 1.0))
	fprintf(stderr,"%s: multiplying A and B by %g\n",argv[0],amplitudefactor);
      if(amplitudefactor == 0.0)
	fprintf(stderr,"%s: WARNING: scaling by zero!\n",argv[0]);
    }
    /*
      
      now apply taper
      
    */     
    if(tapering == SET_L_UNITY){
      ran1(&seed);
      /* 
	 debugging filter, set one set of l to random values, the
	 other to zero
      */
      for(l=0;l<=lmax;l++)
	for(m=0;m <= l;m++){
	  os1 = POSLM(l, m);
	  if(l <= llmax){
	    if((md < 0) || (m == md))
	      fac1 = 1.0;
	    else
	      fac1 = 0.0;
	  }else{
	    fac1 = 0.0;
	  }
	  if(vector_field){
	    *(amp+os1) = fac1 * (-.5+ran1(&seed));	/* should be zero anyway */
	    *(amt+os1) = fac1 * (-.5+ran1(&seed));
	    if(m != 0){
	      *(bmp+os1)=fac1* (-.5+ran1(&seed));
	      *(bmt+os1)=fac1* (-.5+ran1(&seed));
	    }else{
	      *(bmp+os1)=0;
	      *(bmt+os1)=0;
	    }
	  }else{
	    a[os1] =   fac1* (-.5+ran1(&seed));	/* should be zero anyway */
	    if(m!=0)
	      b[os1] = fac1* (-.5+ran1(&seed));
	    else
	      b[os1] = 0;
	  }
	}
      /* reset  */
      llmax = lmax;
      /* end debugging mode */
    }else{
      switch(vector_field){
      case 0:
	// apply taper for AB scalar
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    taperf(tfac,l,lmax,lc,tapering,m,filter,af,bf,
		   amplitudefactor);
	    *(a+os1) *= tfac[0]; *(b+os1) *= tfac[1];
	    if((!interpolate_mode) && use_second){	/* also scale second expansion */
	      *(c+os1) *= tfac[0]; *(d+os1) *= tfac[1];
	    }
	  }
	break;
      case 1: 			/* poloidal/toroidal coefficients */
      case 2:			/* or Real / Imaginary  from GSH */
	switch(tapering){
	case NR_TAPER:
	  /* 
	     remove all but NR components , ie. pass net rotation 
	     only 
	  */
	  *(amp+POSLM(0, 0))=0.0;	/* should be zero anyway */
	  *(amt+POSLM(0, 0))=0.0;
	  /* poloidal l=1 terms */
	  *(amp+POSLM(1, 0))=0.0;
	  *(amp+POSLM(1, 1))=0.0;
	  *(bmp+POSLM(1, 1))=0.0;
	  for(l=2;l<=lmax;l++)	/* all l>= 2 terms are zero */
	    for(m=0;m<=l;m++){
	      os1 = POSLM(l, m);
	      *(amp+os1) = 0.0;
	      *(bmp+os1) = 0.0;
	      *(amt+os1) = 0.0;
	      *(bmt+os1) = 0.0;
	    }
	  if((!interpolate_mode) && (use_second)){
	    *(cmp+POSLM(0, 0))=0.0;	/* should be zero anyway */
	    *(cmt+POSLM(0, 0))=0.0;
	    /* poloidal l=1 terms */
	    *(cmp+POSLM(1, 0))=0.0;
	    *(cmp+POSLM(1, 1))=0.0;
	    *(dmp+POSLM(1, 1))=0.0;
	    for(l=2;l<=lmax;l++)	/* l>=2 terms, reset zero */
	      for(m=0;m<=l;m++){
		os1 = POSLM(l, m);
		*(cmp+os1) = 0.0;
		*(dmp+os1) = 0.0;
		*(cmt+os1) = 0.0;
		*(dmt+os1) = 0.0;
	      }
	  }
	  break;
	case NNR_TAPER:/* 
			  NNR and other tapers are NOT
			  EXCLUSIVE since highpass might 
			  still apply! 
		       */
	
	  /* remove the net rotation component */
	  *(amt+POSLM(1, 0))=0.0;
	  *(amt+POSLM(1, 1))=0.0;
	  *(bmt+POSLM(1, 1))=0.0;
	  if((!interpolate_mode) && (use_second)){
	    *(cmt+POSLM(1, 0))=0.0;
	    *(cmt+POSLM(1, 1))=0.0;
	    *(dmt+POSLM(1, 1))=0.0;
	  }
	  break;
	case PASS_POL_TAPER:
	  if(vector_field == 2){
	    fprintf(stderr,"%s: pass pol doesn't make sense for GSH\n",
		    argv[0]);exit(-1);}
	  /* 
	     remove all toroidal components
	  */
	  for(l=0;l<=lmax;l++)	/* all l>= 2 terms are zero */
	    for(m=0;m <= l;m++){
	      os1 = POSLM(l, m);*(amt+os1) = 0.0;
	      *(bmt+os1) = 0.0;
	    }
	  if((!interpolate_mode) && (use_second)){
	    for(l=0;l<=lmax;l++)	/* l>=2 terms, reset zero */
	      for(m=0;m<=l;m++){
		os1 = POSLM(l, m);*(cmt+os1) = 0.0;
		*(dmt+os1) = 0.0;
	      }
	  }
	  break;
	case PASS_TOR_TAPER:
	  if(vector_field == 2){
	    fprintf(stderr,"%s: pass tor doesn't make sense for GSH\n",
		    argv[0]);exit(-1);}
	  /* 
	     remove all poloidal components
	  */
	  for(l=0;l<=lmax;l++)	/* all l>= 2 terms are zero */
	    for(m=0;m <= l;m++){
	      os1 = POSLM(l, m);*(amp+os1) = 0.0;
	      *(bmp+os1) = 0.0;
	    }
	  if((!interpolate_mode) && (use_second)){
	    for(l=0;l<=lmax;l++)	/* l>=2 terms, reset zero */
	      for(m=0;m<=l;m++){
		os1 = POSLM(l, m);*(cmp+os1) = 0.0;
		*(dmp+os1) = 0.0;
	      }
	  }
	  break;
	default:
	  break;
	}	/* end tapering switch */
      
	/* 
	   scale the rest 
	*/
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    taperf(tfac,l,lmax,lc,tapering,m,filter,af,bf,amplitudefactor);
	    *(amp+os1) *= tfac[0];
	    *(bmp+os1) *= tfac[1];
	    *(amt+os1) *= tfac[0];
	    *(bmt+os1) *= tfac[1];
	    if((!interpolate_mode) && (use_second)){
	      *(cmp+os1) *= tfac[0];
	      *(dmp+os1) *= tfac[1];
	      *(cmt+os1) *= tfac[0];
	      *(dmt+os1) *= tfac[1];
	    }
	  }
	break;			/* end vector field != 0 */
      default:
	fprintf(stderr,"%s: taper for vector_field %i not implemented\n",
		argv[0],vector_field);
	exit(-1);
	break;
      }
    }
    /* 
       

       OUTPUT 


       change l_max

    */
    if(llmax != 99999){
      if(llmax > lmax){
	if(vector_field){
	  fprintf(stderr,"%s: expanding coefficients to new l_max only implemented for single set\n",
		  argv[0]);
	  exit(-1);
	}
	if(verbose)
	  fprintf(stderr,"%s: new lmax value (%i) with original lmax (%i), fill up with zeroes\n",
		  argv[0],llmax,lmax);
	lmsize=(int)((((float)llmax)+1.0)*(((float)llmax)+2)/2.0);
	if((a=(COMP_PRECISION *)realloc(a,lmsize*sizeof(COMP_PRECISION)))==NULL ||
	   (b=(COMP_PRECISION *)realloc(b,lmsize*sizeof(COMP_PRECISION)))==NULL){
	  fprintf(stderr,
		  "%s: memerror while resizing lmax=%i lmsize=%i\n",
		  argv[0],llmax,lmsize); 
	  exit(-1);
	}
	for(l=lmax+1;l<=llmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l,m);
	    a[os1]=b[os1]=0.0;
	  }
	lmax=llmax;
      }else if(llmax<0){
	fprintf(stderr,"%s: new lmax is negative (%i), that's nonsense\n",argv[0],lmax);
	phelp(argv[0]);exit(-1);
      }else if(lmax!=llmax){
	if(verbose)fprintf(stderr,"%s: limiting output from %i to lmax=%i\n",argv[0],lmax,llmax);    
	lmax=llmax;
      }
    } /* end lmax != 99999 */
    
    /* 
       

       output
    
    
    */
    switch(out_format){// begin output format case structure
    case ABPHYS_OUT_LONG:
    case ABPHYS_OUT:{
      if(vector_field){
	/*  */
	if(verbose)
	  fprintf(stderr,"%s: physical convention %s AB format output, lmax=%i\n",
		  argv[0],(out_format==ABPHYS_OUT_LONG)?("long(hc)"):("regular (short)"),lmax);
	switch(tapering){
	case PASS_POL_TAPER:	/* output of poloidal only */
	  write_coefficients(&amp,&bmp,lmax,1.0,
			     (out_format==ABPHYS_OUT_LONG)?(2):(1),
			     stdout,1);
	  break;
	case PASS_TOR_TAPER:	/* output of toroidal */
	  write_coefficients(&amt,&bmt,lmax,1.0,
			     (out_format==ABPHYS_OUT_LONG)?(2):(1),
			     stdout,1);
	  break;
	default:
	  fprintf(stderr,"%s: error physical/DT convention AB format output, but vector_field %i input\n",
		  argv[0],vector_field);
	  exit(-1);
	}
      }else{
	write_coefficients(&a,&b,lmax,1.0,
			   (out_format==ABPHYS_OUT_LONG)?(2):(1),
			   stdout,1);
      }
      break;
    }
    case ABPHYS_NONZERO_OUT_ONE_COLUMN:
    case ABPHYS_NONZERO_OUT:{
      if(vector_field){
	fprintf(stderr,"%s: error physical/DT convention AB format output, but vector_field %i input\n",
		argv[0],vector_field);
	exit(-1);
      }
      if(verbose)
	fprintf(stderr,"%s: physical convention AB-nonzero format output, lmax=%i\n",argv[0],lmax);
      write_nonzero_coefficients(&a,&b,lmax,1.0,TRUE,stdout,1,(out_format==ABPHYS_NONZERO_OUT_ONE_COLUMN));
      break;
    }
    case ABGEOD_OUT:{
      if(verbose)fprintf(stderr,"%s: geodetic convention A'B' format output, lmax=%i\n",argv[0],lmax);
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  fprintf(stdout,"%21.14e\n",
		  *(a+os1)/GEODETIC_FACTOR(l,m));
	}
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  fprintf(stdout,"%21.14e\n",
		  *(b+os1)/GEODETIC_FACTOR(l,m));
	}
      break;
    }
    case POWER_OUT:
    case POWER_OUT_NN:{
      if(out_format == POWER_OUT_NN){
	normalize_power_by_ncoeff=FALSE;
	if(verbose)
	  fprintf(stderr,"%s: not using number of coefficient normalization for power (i.e. not per area)\n",
		  argv[0]);
      }else{
	normalize_power_by_ncoeff=TRUE;
      }
      if(!vector_field){
	if(verbose){
	  if(out_format == POWER_OUT_NN)
	    fprintf(stderr,"%s: output: l power_per_degree, lmax=%i\n",
		    argv[0],lmax);
	  else
	    fprintf(stderr,"%s: output: l power_per_degree_per_unit_area, lmax=%i\n",
		    argv[0],lmax);
	}
	    
	for(l=0;l<=lmax;l++)
	  fprintf(stdout,"%i %21.14e\n",l,degree_power(a,b,l,normalize_power_by_ncoeff));
      }else{
	if(vector_field == 1){
	  if(verbose)
	    fprintf(stderr,"%s: output: l pwr/degree/unit_area|_pol  pwr/degree/unit_area|_tor, lmax=%i\n",
		    argv[0],lmax);
	  for(l=0;l<=lmax;l++){
	    fprintf(stdout,"%i %21.14e %21.14e\n",
		    l,degree_power(amp,bmp,l,normalize_power_by_ncoeff),
		    degree_power(amt,bmt,l,normalize_power_by_ncoeff));
	  }
	}else{			/* power for GSH  */

	  if(verbose)
	    fprintf(stderr,"%s: output: l GSH pwr/degree/unit_area , lmax=%i\n",
		    argv[0],lmax);
	  for(l=0;l<=lmax;l++){
	    fprintf(stdout,"%i %21.14e\n",
		    l,degree_power_gsh(amp,amt,bmp,bmt,l,normalize_power_by_ncoeff));
	  }
	}
      }
      break;
    }
    case MEAN_OUT:{
      if(vector_field){
	fprintf(stderr,"%s: mean output undefined for vector fields\n",argv[0]);
	exit(-1);
      }else{
	if(verbose)
	  fprintf(stderr,"%s: output of mean (scaled l = 0 term)\n",argv[0]);
	fprintf(stdout,"%21.14e\n",a[POSLM(0,0)]/TWO_SQRT_PI);
      }
      break;
    }
    case TRMS_OUT:{ 
      if(!vector_field){
	if(verbose)
	  fprintf(stderr,"%s: output: RMS of expansion (no l = 0 terms!), lmax=%i\n",
		  argv[0],lmax);
	fprintf(stdout,"%21.14e\n",calc_rms(a,b,lmax));
      }else{
	if(vector_field == 1){
	  if(verbose)
	    fprintf(stderr,"%s: output: RMS/unit_area|_pol  RMS/unit_area|_tor, lmax=%i\n",
		    argv[0],lmax);
	  fprintf(stdout,"%21.14e %21.14e\n",
		  calc_rms(amp,bmp,lmax),calc_rms(amt,bmt,lmax));
	}else{
	  if(verbose)
	    fprintf(stderr,"%s: output: GSH RMS/unit_area|, lmax=%i\n",
		    argv[0],lmax);
	  fprintf(stdout,"%21.14e\n",
		  calc_rms_gsh(amp,amt,bmp,bmt,lmax));
	}
      }
      break;
    }
    case TPOWER_OUT:{ 
      if(!vector_field){
	if(verbose)
	  fprintf(stderr,"%s: output: magnitude (sqrt(total power)) of expansion, lmax=%i\n",
		  argv[0],lmax);
	fprintf(stdout,"%21.14e\n",calc_total_power(a,b,lmax));
      }else{
	if(vector_field == 1){
	  if(verbose)
	    fprintf(stderr,"%s: output: total power/unit_area|_pol  total power/unit_area|_tor, lmax=%i\n",
		    argv[0],lmax);
	  fprintf(stdout,"%21.14e %21.14e\n",
		  calc_total_power(amp,bmp,lmax),calc_total_power(amt,bmt,lmax));
	}else{
	  if(verbose)
	    fprintf(stderr,"%s: output: total power/unit_area, lmax=%i\n",
		    argv[0],lmax);
	  fprintf(stdout,"%21.14e\n",
		  calc_total_power_gsh(amp,amt,bmp,bmt,lmax));
	}
      
      }
      break;
    }
    case TAPER_OUT:{
      if(verbose)
	fprintf(stderr,"%s: shape of the taper function\n",argv[0]);
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  taperf(tfac,l,lmax,lc,tapering,m,filter,af,bf,
		 amplitudefactor);	  
	  fprintf(stdout,"%i %i %g %g %g\n",
		  l,m,(COMP_PRECISION)l/(COMP_PRECISION)lmax,
		  tfac[0],tfac[1]);
	}
      break;
    }
    case VECABAB_POL_OUT:{
      if(!vector_field){
	fprintf(stderr,"%s: did not read pol/tor coefficients, output %i does not work\n",
		argv[0],VECABAB_POL_OUT);
	exit(-1);
      }else{
	if(vector_field == 2){
	  fprintf(stderr,"%s: GSH format not implemented for out mode %i\n",
		  argv[0],out_format);
	  exit(-1);
	}
	fprintf(stderr,"%s: output of poloidal part of expansion\n",argv[0]);
	fprintf(stdout,"%i\n",lmax);
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    fprintf(stdout,"%21.14e %21.14e\n",
		    *(amp+os1),*(bmp+os1));
	  }
      }
      break;
    }
    case VECABAB_TOR_OUT:{
      if(!vector_field){
	fprintf(stderr,"%s: did not read pol/tor coefficients, output %i does not work\n",
		argv[0],VECABAB_POL_OUT);
	exit(-1);
      }else{
	if(vector_field == 2){
	  fprintf(stderr,"%s: GSH format not implemented for out mode %i\n",
		  argv[0],out_format);
	  exit(-1);
	}
	fprintf(stderr,"%s: output of toroidal part of expansion\n",argv[0]);
	fprintf(stdout,"%i\n",lmax);
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    fprintf(stdout,"%21.14e %21.14e\n",
		    *(amt+os1),*(bmt+os1));
	  }
      }
      break;
    }
    case LMAB_OUT:{
      if(!vector_field){
	if(verbose)
	  fprintf(stderr,"%s: physical convention l m A B format output, lmax=%i\n",argv[0],lmax);
 	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    fprintf(stdout,"%i %i %21.14e %21.14e\n",l,m,
		    *(a+os1),*(b+os1));
	  }
      }else{
	if(vector_field == 1){
	  if(verbose)
	    fprintf(stderr,"%s: physical convention l m Ap Bp At Bt format output, lmax=%i\n",argv[0],lmax);
 	  fprintf(stdout,"%i\n",lmax);
	  for(l=0;l<=lmax;l++)
	    for(m=0;m<=l;m++){
	      os1 = POSLM(l, m);
	      fprintf(stdout,"%i %i %21.14e %21.14e %21.14e %21.14e\n",
		      l,m,
		      *(amp+os1),*(bmp+os1),
		      *(amt+os1),*(bmt+os1));
	    }
	}else{
	  if(verbose)
	    fprintf(stderr,"%s: physical convention l m Ar Ai Br Bi format output, lmax=%i\n",argv[0],lmax);
 	  fprintf(stdout,"%i %i\n",lmax,ialpha);
	  for(l=0;l<=lmax;l++)
	    for(m=0;m<=l;m++){
	      os1 = POSLM(l, m);
	      fprintf(stdout,"%i %i %21.14e %21.14e %21.14e %21.14e\n",
		      l,m,
		      *(amp+os1),*(amt+os1),
		      *(bmp+os1),*(bmt+os1));
	    }

	}
      }
      break;
    } 
    case LMAB_GEODETIC_OUT:{
      if(verbose)
	fprintf(stderr,"%s: geodetic convention l m A B format output, lmax=%i\n",argv[0],lmax);
      if(!vector_field){
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    fprintf(stdout,"%i %i %21.14e %21.14e\n",l,m,
		    *(a+os1)/GEODETIC_FACTOR(l,m), 
		    *(b+os1)/GEODETIC_FACTOR(l,m));
	  }
      }else{
	if(vector_field == 2){
	  fprintf(stderr,"%s: GSH format not implemented for out mode %i\n",
		  argv[0],out_format);
	  exit(-1);
	}
	fprintf(stdout,"%i\n",lmax);
	for(l=0;l<=lmax;l++)
	  for(m=0;m<=l;m++){
	    os1 = POSLM(l, m);
	    fprintf(stdout,"%i %i %21.14e %21.14e %21.14e %21.14e\n",
		    l,m,
		    *(amp+os1)/GEODETIC_FACTOR(l,m),*(bmp+os1)/GEODETIC_FACTOR(l,m),
		    *(amt+os1)/GEODETIC_FACTOR(l,m),*(bmt+os1)/GEODETIC_FACTOR(l,m));
	  }
      }
      break;
    }
    case ADMITTANCE_OUT_NN:
    case ADMITTANCE_OUT:{
      if(vector_field){
	fprintf(stderr,"%s: admittance not implemented for vector field\n",argv[0]);
	exit(-1);
      }
      if(verbose)fprintf(stderr,"%s: computing admittance power(S1,S2)/power(S2), lmax=%i\n",
			 argv[0],lmax);
      if(out_format == ADMITTANCE_OUT_NN){
	normalize_power_by_ncoeff = FALSE;
	if(verbose)
	  fprintf(stderr,"%s: not normalizing by number of coefficients\n",argv[0]);
      }else{
	normalize_power_by_ncoeff = TRUE;
      }
      for(l=1;l <= lmax;l++){
	tmp = admittance(a,b,c,d,l,normalize_power_by_ncoeff);
	if(finite(tmp))
	  fprintf(stdout,"%5i %20.10lf\n",l,tmp);
	else
	  fprintf(stdout,"%5i         nan\n",l);
      }
      break;
    }
    case SPEAR_CORRL_OUT:
    case CORRL_OUT:{ 
      /* 

	 calculate correlation coefficient as a function of l between two
	 files
      
      */
      if(verbose)fprintf(stderr,"%s: %s correlation coefficients of each degree, lmax=%i\n",
			 argv[0],(cmode==1)?("linear"):("rank"),lmax);
      for(l=1;l <= lmax;l++){
	switch(vector_field){
	case 0:			/* scalar */
	  tmp = correlation(a,b,c,d,l,0,lmin,cmode);
	  break;
	case 1:			/* poloidal/toroidal */
	  tmp = correlation_pt(amp,bmp,amt,bmt,
			       cmp,dmp,cmt,dmt,l,0,lmin,cmode);
	  break;
	case 2:			/* GSH 2phi or 4phi  */
	  tmp = correlation_gsh(amp,amt,bmp,bmt,
				cmp,cmt,dmp,dmt,l,0,ialpha,lmin,cmode);
	  break;
	default:
	  fprintf(stderr,"%s: vector field type %i  not implemented for out mode %i\n",
		  argv[0],vector_field,out_format);
	  exit(-1);
	  break;
	}
	if(finite(tmp))
	  fprintf(stdout,"%5i %20.10lf\n",l,tmp);
	else
	  fprintf(stdout,"%5i         nan\n",l);
      }
      break;
    }
      /* hmmm, why?  */
    case CCL_COUPLING:{
      if(verbose){fprintf(stderr,"%s: cross l coupling\n",argv[0]);}
      switch(vector_field){
      case 0:
	for(l=1;l<lmax;l++){
	  fac1 = ccl_correlation(a,b,l,lmax,&i,1);
	  fprintf(stdout,"%6.1f %20.10lf %i\n",(float)l+0.5,fac1,i); /* l+0.5 r_l,l+1 n(r_l,l+1) */
	}
	break;
      default:
	fprintf(stderr,"%s: vector field type %i  not implemented for out mode %i\n",
		argv[0],vector_field,out_format);
	exit(-1);
	break;
      }

      break;
    }
    case SPEAR_CORRT_OUT:
    case CORRT_OUT:{ /* calculate total correlation coefficient
			between two files */
      if(verbose)fprintf(stderr,"%s: total %s correlation coefficient, summed over all l, lmax=%i\n",
			 argv[0],(cmode==1)?("linear"):("rank"),lmax);
      switch(vector_field){
      case 0:			/* scalar */
	tmp = correlation(a,b,c,d,-lmax,0,lmin,cmode);
	break;
      case 1:			/* poloidal/toroidal */
	tmp = correlation_pt(amp,bmp,amt,bmt,cmp,dmp,cmt,dmt,-lmax,0,lmin,cmode);
	break;
      case 2:			/* GSH 2phi/4phi */
	tmp = correlation_gsh(amp,amt,bmp,bmt,cmp,cmt,dmp,dmt,-lmax,0,ialpha,lmin,cmode);
	break;
      default:
	fprintf(stderr,"%s: vector field type %i  not implemented for out mode %i\n",
		argv[0],vector_field,out_format);
	exit(-1);
	break;
      }
      fprintf(stdout,"%20.10lf\n",tmp);
      break;
    }
    case SPEAR_CORRTH_OUT:
    case CORRTH_OUT:{ /* calculate total correlation coefficient
			 between two files, limited to degrees higher
			 than mc*/
      lmin = md;
      if(verbose)fprintf(stderr,"%s: total %s correlation coefficient, summed from %i to %i\n",
			 argv[0],(cmode==1)?("linear"):("rank"),lmin,lmax);
      switch(vector_field){
      case 0:
	fprintf(stdout,"%20.10lf\n",correlation(a,b,c,d,-lmax,0,lmin,cmode));
	break;
      case 2:
	fprintf(stdout,"%20.10lf\n",correlation_gsh(amp,amt,bmp,bmt,cmp,cmt,dmp,dmt,
						    -lmax,0,ialpha,lmin,cmode));
	break;
      default:
	fprintf(stderr,"%s: vector field type %i  not implemented for out mode %i\n",
		argv[0],vector_field,out_format);
	exit(-1);
	break;
      }
      break;
    }
    case VECABAB_OUT_LONG:
    case VECABAB_OUT:{
      if(!vector_field){
	fprintf(stderr,"%s: velocity output requested, yet no velocity input read\n",argv[0]);
	exit(-1);
      }
      if(vector_field == 1){
	if(verbose)
	  fprintf(stderr,"%s: poltor vector output in ABAB format (physical, %s), lmax=%i\n",argv[0],
		  (out_format==VECABAB_OUT_LONG)?("long(hc)"):("regular (short)"),lmax);
	write_vector_coefficients(&amp,&amt,&bmp,&bmt,lmax,1.0,(out_format==VECABAB_OUT_LONG)?(2):(1),stdout,1);
      }else{			/* GSH */
	if(verbose)
	  fprintf(stderr,"%s: GSH output in Ar Ai Br Bi format, lmax=%i\n",argv[0],lmax);
	fprintf(stdout,"%i %i\n",lmax,ialpha);
	/* note this has mixed arguments on purpose */
	write_vector_coefficients(&amp,&bmp,&amt,&bmt,lmax,1.0,0,stdout,1);
      }
      break;
    }
    case VECABAB_NEW_OUT:{
      if(!vector_field){
	fprintf(stderr,"%s: velocity output requested, yet no velocity input read\n",argv[0]);
	exit(-1);
      }
      if(vector_field == 2){
	fprintf(stderr,"%s: GSH format not implemented for out mode %i\n",
		argv[0],out_format);
	exit(-1);
      }
      if(verbose)
	fprintf(stderr,"%s: poltor vector output in new ABAB format (physical), lmax=%i\n",argv[0],lmax);
      write_vector_coefficients(&amp,&amt,&bmp,&bmt,lmax,1.0,1,stdout,1);
      break;
    }
    case VECAABBR_OUT:{
      if(!vector_field){
	fprintf(stderr,"%s: velocity output requested, yet no velocity input read\n",argv[0]);
	exit(-1);
      }
      if(vector_field == 2){
	fprintf(stderr,"%s: GSH format not implemented for out mode %i\n",
		argv[0],out_format);
	exit(-1);
      }
      if(verbose)fprintf(stderr,"%s: poltor vector output in ABAB format (Rick's convention), lmax=%i\n",argv[0],lmax);
      fprintf(stdout,"%i\n",lmax);
      for(l=0;l<=lmax;l++)
	for(m=0;m<=l;m++){
	  os1 = POSLM(l, m);
	  fprintf(stdout,"%21.14e %21.14e\n%21.14e %21.14e\n",
		  *(amp+os1)/RICK_FACTOR(l,m),
		  *(amt+os1)/RICK_FACTOR(l,m),
		  *(bmp+os1)/RICK_FACTOR(l,m),
		  *(bmt+os1)/RICK_FACTOR(l,m));
	}
      break;
    }
    case GSH_OUT:{
      if(ialpha < 0){
	if(!vector_field){
	  fprintf(stderr,"%s: assuming we have scalar fields\n",argv[0]);
	  ialpha = 0;
	}else{
	  fprintf(stderr,"%s: GSH output selected, but no GSH read\n",
		  argv[0]);
	  exit(-1);
	}
      }
      if(verbose)
	fprintf(stderr,"%s: GSH type %i BW convention output, lmax=%i\n",
		argv[0],ialpha,lmax);
      write_gsh_coeff_set(ialpha, lmax, a, b, amp,amt, bmp, bmt,
			  1.0,stdout);
      break;
    }
    default:{
      phelp(argv[0]);
      exit(-1);
      break;
    }}// end switch
  }// end scan of multiple input files loop
  return 0;
}
//
// this function returns tapering factors for A and B
// as a function of l and m
//
// normally, the two tapering factors will be identially, and
// only a function of l, not m. the exception is when we're scaling
// with another SH file
//
// output: fac[2] for a and b
//
void taperf(COMP_PRECISION *fac,
	    int l, int lmax, 
	    COMP_PRECISION limit, int mode, int m,
	    COMP_PRECISION *filter,
	    COMP_PRECISION *af, COMP_PRECISION *bf,
	    COMP_PRECISION amplitudefactor)
{
  /* should not taper in m since they form a full set */
  COMP_PRECISION lp,tmp,exp_scale=10000;
  /*  */
  static COMP_PRECISION r2fac = 1/(REARTH_KM*REARTH_KM);
  static COMP_PRECISION r2fac_vel = 0.01/(REARTH_KM)*1e6;
  
  /* fractional l */
  lp=(COMP_PRECISION)l/(COMP_PRECISION)lmax;
  if((lp > 1.0)||(lp < 0)){
    fprintf(stderr,"taperf: error l (%i) > lmax (%i) or smaller than zero\n",
	    l,lmax);
    exit(-1);
  }
  if(mode == EXP_TAPER){	/* different interpretation of limit
				   factor */
    exp_scale = limit;
    limit = 0.0;
  }

  if(lp < limit){
    //
    // are we to the left of the tapering limit?
    //
    switch(mode){
    case L0_TAPER: // only l=0 taper does something then
      if(l == 0)
	fac[0] = fac[1] = 0.0;
      else
	fac[0] = fac[1] = 1.0;
      break;
    case ONLY_L_ZERO:
      if(l != 0)
	fac[0] = fac[1] = 0.0;
      else
	fac[0] = fac[1] = 1.0;
      break;
    case ONLY_M0_TERMS:
      if(m == 0){
	fac[0] = fac[1] = 1.0;
      }else{
	fac[0] = fac[1] = 0.0;
      }
      break;
    case ZERO_TAPER:
      fac[0] = fac[1] = 1.0;
      break;
    case HIGHP_TAPER:
      fac[0] = fac[1] = 0.0;
      break;
    default:
      fac[0] = fac[1] = 1.0;
      break;
    }
  }else{
    // rescale lp 
    if(limit == 1.0){
      fprintf(stderr,"taperf: error: limit = 1.0\n");
      exit(-1);
    }
    lp = (lp-limit)/(1.0-limit);
    switch(mode){
    case PHI_ROTATE:
    case NNR_TAPER:
    case NR_TAPER:
    case PASS_TOR_TAPER:
    case PASS_POL_TAPER:
    case NO_TAPER:{
      fac[0] = fac[1] = 1.0;
      break;
    }
    case COSSQR_TAPER:{
      tmp=cos(lp*PIHALF);
      fac[0] = fac[1] = SQUARE(tmp);
      break;
    }
    case COSP4_TAPER:{
      tmp=cos(lp*PIHALF);
      tmp= SQUARE(tmp);
      fac[0] = fac[1] = SQUARE(tmp);
      break;
    }
    case  ONEML_TAPER:{
      fac[0] = fac[1] =  1.0-lp;
      break;
    }
    case ONEMLSQR_TAPER:{
      fac[0] = fac[1] = 1.0-SQUARE(lp);
      break;
    }
    case SINOVRL_TAPER:{
      tmp=(lp != 0.0)?(sin(PI*lp)/(PI*lp)):(1.0);
      fac[0] = fac[1] = tmp;
      break;
    }
    case EXP_TAPER:{
      fac[0] = fac[1] = exp(-pow((COMP_PRECISION)l/exp_scale,2));
      break;
    }
    case ZERO_TAPER:{
      fac[0] = fac[1] = 0.0;
      break;
    }
    case HIGHP_TAPER:{
      fac[0] = fac[1] = 1.0;
      break;
    }
    case L0_TAPER:{
      if(l==0){
	fac[0] = fac[1] = 0.0;
      }else{
	fac[0] = fac[1] = 1.0;
      }
      break;
    }
    case ONLY_M0_TERMS:
      if(m == 0){
	fac[0] = fac[1] = 1.0;
      }else{
	fac[0] = fac[1] = 0.0;
      }
      break;
 
    case ONLY_L_ZERO:{
      if(l != 0){
	fac[0] = fac[1] = 0.0;
      }else{
	fac[0] = fac[1] = 1.0;
      }
      break;
    }
    case FROM_FILE_TAPER:{
      fac[0] = fac[1] = filter[l];
      break;
    }  
    case LAPLACIAN:{
      fac[0] = fac[1] = -((COMP_PRECISION)l*((COMP_PRECISION)l+1.0))*r2fac;
      break;
    }  
    case VORTICITY:
    case DIVERGENCE:{		/* divergence of pol/tor field */
      fac[0] = fac[1] = -sqrt(((COMP_PRECISION)l*((COMP_PRECISION)l+1.0)))*r2fac_vel;
      break;
    }  
    case FROM_SH_FILE_TAPER:{
      fac[0] = af[POSLM(l,m)];
      fac[1] = bf[POSLM(l,m)];
      break;
    }
    default:{
      fprintf(stderr,"don't know taper number %i\n",
	      mode);
      exit(-1);
    }
    }
  }
  fac[0] *= amplitudefactor;
  fac[1] *= amplitudefactor;
  return;
}




void phelp(char *name)
{
  fprintf(stderr,"%s  [output fmt, %i] [tapering, %i] [ampl, %g] [input fmt, %i] [lmax, lmax] [lc, %g]\n",
	  name,OUT_FORMAT_DEFAULT,TAPERING_DEFAULT,AMPLITUDE_DEFAULT,IN_FORMAT_DEFAULT,LC_DEFAULT);
  fprintf(stderr,"\t  Reads spherical harmonic coefficents A B  from stdin\n");
  fprintf(stderr,"\t  and converts the coefficients into different format.\n");
  fprintf(stderr,"\t  Works for various combinations of normalizations and\n");
  fprintf(stderr,"\t  for scalar and vector fields.\n");
  fprintf(stderr,"\t  Our internal convention is ""theoretical physics"", or ""physical norm"" as\n");
  fprintf(stderr,"\t  in Dahlen and Tromp (1998, B.8) (DT)\n");
  fprintf(stderr,"\t  output fmt:    Output format to stdout\n");
  fprintf(stderr,"\t                 %i: AB format, physical norm. (Dahlen & Tromp conv.)\n",
	  ABPHYS_OUT);
  fprintf(stderr,"\t                 %i: AB format, physical norm, long header format as in hc\n",
	  ABPHYS_OUT_LONG);
  
  fprintf(stderr,"\t                 %i: AB format, geodetic norm. (NASA conv.)\n\n",ABGEOD_OUT);
  fprintf(stderr,"\t                 %i: power per unit area and degree, as in DT (delta = flat)\n",POWER_OUT);
  fprintf(stderr,"\t                 %i: power per degree, -not- normalizing by number of coeff.\n",POWER_OUT_NN);
  fprintf(stderr,"\t                 %i: RMS of expansion (no mean, l=0, terms) \n",TRMS_OUT);
  fprintf(stderr,"\t                 %i: mean of expansion (scaled l=0 term) \n",MEAN_OUT);
  fprintf(stderr,"\t                 %i: total power or magnitude of expansion (include l=0 term)\n\n",
	  TPOWER_OUT);

  fprintf(stderr,"\t                 %i: taper as a function of l (for debugging)\n",TAPER_OUT);
  fprintf(stderr,"\t                 %i: l m AB format, physical/DT norm\n",LMAB_OUT);
  fprintf(stderr,"\t                     or l m Ap Bp At Bt for vector\n");
  fprintf(stderr,"\t                     or l m Ar Ai Br Bi for GSH\n\n");

  fprintf(stderr,"\t                 %i: correlation coefficient per degree\n",CORRL_OUT);
  fprintf(stderr,"\t                     using two input files in AB or GSH format (linear, Pearson)\n");
  fprintf(stderr,"\t                 %i: total correlation coefficient, \n",CORRT_OUT);
  fprintf(stderr,"\t                     using two input files in AB or GSH format (linear, Pearson)\n");
  fprintf(stderr,"\t                 %i: total correlation coefficient, \n",CORRTH_OUT);
  fprintf(stderr,"\t                     but restricted from md to lmax (linear, Pearson)\n\n");

  fprintf(stderr,"\t                 %i: correlation coefficient per degree\n",SPEAR_CORRL_OUT);
  fprintf(stderr,"\t                     using two input files in AB or GSH format (ranked, Spear)\n");
  fprintf(stderr,"\t                 %i: total correlation coefficient, \n",SPEAR_CORRT_OUT);
  fprintf(stderr,"\t                     using two input files in AB or GSH format (ranked, Spear)\n");
  fprintf(stderr,"\t                 %i: total correlation coefficient, \n",SPEAR_CORRTH_OUT);
  fprintf(stderr,"\t                     but restricted from md to lmax (ranked, Spear)\n\n");

  fprintf(stderr,"\t                 %i: admittance - power(C1,C2)/power(C2,C2)\n",ADMITTANCE_OUT);
  fprintf(stderr,"\t                 %i: admittance - not normalizing power by numbers of coefficients\n\n",ADMITTANCE_OUT_NN);
  
  fprintf(stderr,"\t                 %i: cross degree correlation within one expansion (linear, Pearson)\n\n",CCL_COUPLING);

  fprintf(stderr,"\t                 %i: vector field, A_p B_p A_t B_t format (phys./DT norm)\n",VECABAB_OUT);
  fprintf(stderr,"\t                     or GSH 2phi/4phi in  A_r A_i B_r B_i format \n");
  fprintf(stderr,"\t                 %i: vector field, phys./DT norm, long (hc) header\n",VECABAB_OUT_LONG);
  fprintf(stderr,"\t                 %i: vector field, A_p A_t\\nB_p B_t format (Rick's norm)\n",
	  VECAABBR_OUT);
  fprintf(stderr,"\t                 %i: l m AB format (geodetic norm)\n",LMAB_GEODETIC_OUT);
  fprintf(stderr,"\t                 %i: vector field, physical norm, poloidal part as scalar only (DT)\n",
	  VECABAB_POL_OUT);
  fprintf(stderr,"\t                 %i: vector field, physical norm, toroidal part as scalar only (DT)\n",
	  VECABAB_TOR_OUT);
  fprintf(stderr,"\t                 %i: vector field, A_p B_p A_t B_t format (phys./DT norm), new format for HC\n",
	  VECABAB_NEW_OUT);
  fprintf(stderr,"\t                 %i: GSH format, BW normalization\n",GSH_OUT);
  fprintf(stderr,"\t                 %i: AB format, physical norm. (Dahlen & Tromp conv.) skip B term for m=0\n",ABPHYS_NONZERO_OUT);
  fprintf(stderr,"\t                 %i: \t same, use one column instead of two entries per line\n\n",ABPHYS_NONZERO_OUT_ONE_COLUMN);
  fprintf(stderr,"\t  tapering:      %i:  no tapering, l'=l/l_{max}\n",NO_TAPER);
  fprintf(stderr,"\t                 %i:  cos(-pi/2*l')**2  tapering for l' >= lc \n",COSSQR_TAPER);
  fprintf(stderr,"\t                 %i:  cos(-pi/2*l')**4  tapering for l' >= lc \n",COSP4_TAPER);
  fprintf(stderr,"\t                 %i:  1-l'              tapering for l' >= lc \n",ONEML_TAPER);
  fprintf(stderr,"\t                 %i:  1-(l')**2         tapering for l' >= lc \n",ONEMLSQR_TAPER);
  fprintf(stderr,"\t                 %i:  sin(pi*l')/(pi*l')tapering for l' >= lc \n",SINOVRL_TAPER);
  fprintf(stderr,"\t                 %i:  exp(-(l/lc)**2) tapering, set lc to something but unity!\n\n",EXP_TAPER);
  fprintf(stderr,"\t                 %i:  no net rotation for pol/tor AB\n",NNR_TAPER);
  fprintf(stderr,"\t                 %i:  only net rotation for pol/tor AB\n",NR_TAPER);
  fprintf(stderr,"\t                 %i:  pass only the poloidal field for ABAB\n",PASS_POL_TAPER);
  fprintf(stderr,"\t                 %i:  pass only the toroidal field for ABAB\n\n",PASS_TOR_TAPER);

  fprintf(stderr,"\t                 %i:  all coefficients l'>=lc set to zero (see -lmax for highpass)\n",ZERO_TAPER);
  fprintf(stderr,"\t                 %i:  l=0 term is set to zero (mean=0)\n",L0_TAPER);
  fprintf(stderr,"\t                 %i:  only l=0 term is passed  (RMS=0)\n",ONLY_L_ZERO);
  fprintf(stderr,"\t                 %i:  only m=0 terms are passed\n\n",ONLY_M0_TERMS);
  
  fprintf(stderr,"\t                 %i:  compute Laplacian (multiply with - l(l+1))\n\n",LAPLACIAN);
  fprintf(stderr,"\t                 %i:  compute  divergence         of PT field (multiply pol with - (l(l+1))/Re^2)\n",DIVERGENCE);
  fprintf(stderr,"\t                 %i:  compute  vertical vorticity of PT field (multiply tor with - (l(l+1))/Re^2\n\n",VORTICITY);

  fprintf(stderr,"\t                 %i:  read in w_0 ... w_{l_{max}} weights from file \"%s\"\n",
	  FROM_FILE_TAPER,FILTER_FILE);
  fprintf(stderr,"\t                      and scale all coefficients at l with the weight.\n");
  fprintf(stderr,"\t                      In this case, lc is set to zero by default.\n");
  fprintf(stderr,"\t                 %i:  read in phys norm AB SH file which will multiply input, from file \"%s\"\n",
	  FROM_SH_FILE_TAPER,FILTER_SH_FILE);
  fprintf(stderr,"\t                      This will multiply the expansion with those coefficients.\n");
  fprintf(stderr,"\t                      In this case, lc is set to zero by default.\n");

  fprintf(stderr,"\t                 %i:  rotate field by amount dp, where dp is given as the multiplication factor\n",
	  PHI_ROTATE);
  fprintf(stderr,"\t                      factor, in degrees. E.g. amp=10 means rotate by 10 degrees around\n");
  fprintf(stderr,"\t                      the rotation axis in easterly direction\n\n");
  
  fprintf(stderr,"\t                 %i:  debugging mode, sets all coefficients for l <= lmax to random (md is seed) values, the rest to zero\n",
	  SET_L_UNITY);
  fprintf(stderr,"\t                      if md >= 0, will only set l = llmax and m = md coefficients to random values\n\n");
  fprintf(stderr,"\t  ampl:          multiplication (scaling) factor, scales SH for most filters\n");
  fprintf(stderr,"\t                     does have a different meaning for rotations, and interpolation\n\n");
  fprintf(stderr,"\t  input fmt:     %i: AB format       physical/DT norm.\n",AB_INPUT);
  fprintf(stderr,"\t                 %i: AB format, physical, long (hc) format\n",AB_INPUT_HC);
  fprintf(stderr,"\t                 %i: layer AB format geodetic norm.\n",LAB_GEOD_INPUT);
  fprintf(stderr,"\t                 %i: AB format       geodetic norm.\n",AB_GEOD_INPUT);
  fprintf(stderr,"\t                 %i: AB format       Rick ""fully norm""\n",AB_RICK_INPUT);
  fprintf(stderr,"\t                 %i: vector field, A_p B_p A_t B_t format (physical/DT conv.)\n",
	  ABAB_INPUT);
  fprintf(stderr,"\t                 %i: vector field, A_p B_p A_t B_t format (geodetic conv.)\n",
	  ABAB_GEOD_INPUT);
  fprintf(stderr,"\t                 %i: vector field, A_p A_t \\n B_p B_t format (Rick's conv.)\n",
	  AABBR_INPUT);
  fprintf(stderr,"\t                 %i: l m A B format  geodetic (Rick)\n",LMAB_GEOD_INPUT);
  fprintf(stderr,"\t                 %i: interpolate two sets of coefficients\n",INTERPOLATE);
  fprintf(stderr,"\t                     using two input files in AB/DT format\n");
  fprintf(stderr,"\t                     in this case, the ampl factor gives the weight of the\n");
  fprintf(stderr,"\t                     first file, the second is weighted with 1-ampl\n");
  fprintf(stderr,"\t                 %i: interpolate two sets of vector coefficients (DT)\n",INTERPOLATE_ABAB);
  fprintf(stderr,"\t                     using two input files in A_p B_p A_t B_t  format as in %i\n",
	  ABAB_INPUT);
  fprintf(stderr,"\t                     ampl factor meaning same as for input mode %i\n",INTERPOLATE);
  fprintf(stderr,"\t                 %i: interpolate two sets of GSH vector coefficients \n",INTERPOLATE_GSH);
  fprintf(stderr,"\t                     using two input files in GSH format as in %i\n",GSH_INPUT);
  fprintf(stderr,"\t                     ampl factor meaning same as for input mode %i\n",INTERPOLATE);
  fprintf(stderr,"\t                 %i: l m A B ""fully normalized"" (Edmonds, 1960) aka old Harvard\n",LMAB_FNORM_INPUT);
  fprintf(stderr,"\t                 %i: AB format  G. Master's ""fully norm"" (Edmonds, 1960\n",AB_MASTERS_INPUT);
  fprintf(stderr,"\t                 %i: generalized spherical harmonics input (0, 2, or 4 type, BW convention)\n",GSH_INPUT);
  fprintf(stderr,"\t                 %i: AB format physical/DT norm, B not listed for m =0\n\n",AB_NONZERO_INPUT);

  fprintf(stderr,"\t  lmax:          maximum degree of output, by default = input lmax\n");
  fprintf(stderr,"\t                 (if zero, use input lmax)\n");
  fprintf(stderr,"\t                 (if negative, set all coefficients SMALLER than lmax zero(highpass))\n\n");

  fprintf(stderr,"\t  lc:            start of tapering (l/lmax has to be >=lc), lc is by default = %g,\n",LC_DEFAULT);
  fprintf(stderr,"\t                 for filter type tapering, it is always set to zero.\n");            
  fprintf(stderr,"\t                 for exp type tapering, lc is sort of the half width.\n\n");            
  fprintf(stderr,"\t  md:            set only l = lmax, m = md to random values, for taper mode %i (default: <0, all m's, used as seed)\n",
	  SET_L_UNITY);
  fprintf(stderr,"\t                 also used for limited corrleation mode %i\n",CORRTH_OUT);

  fprintf(stderr,"\n\n");
}

