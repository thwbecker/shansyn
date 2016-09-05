/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: mod_modelbase.c,v 1.8 2009/04/16 00:55:29 becker Exp $ */
#include "shansyn.h"

/*
  convert a model given by spherical harmonics coefficients
  at different, discrete depth layers to a chebyshev polynomial model 

*/

#define PWR_OUT_FILE "pwr.dat"
void fit_base_functions(COMP_PRECISION,COMP_PRECISION,
			COMP_PRECISION *,COMP_PRECISION *,
			struct mod *, 
			COMP_PRECISION (*)(COMP_PRECISION ,
					   COMP_PRECISION ,
					   COMP_PRECISION *,
					   int ,COMP_PRECISION),
			COMP_PRECISION (*)(int,int,int));

int main(int argc, char **argv)
{
 
  int i,lms,order,mode,iter,maxiter,hit;
  COMP_PRECISION lambda,omega,vr,totmsize,
    targetvr;
  int expect_gsh = 0;
  struct mod model[2];
  FILE *out;
  // default values
  lambda=0;
  omega=0.0;
  order=20;
  mode=SPLINES;
  
  maxiter=15;
  targetvr=0.95;

  switch(argc){
  case 2:{
    break;
  }
  case 3:{
    sscanf(argv[2],"%i",&mode);
    break;
  }
  case 4:{
    sscanf(argv[2],"%i",&mode);
    sscanf(argv[3],"%i",&order);
    break;
  }
  case 5:{
    sscanf(argv[2],"%i",&mode);
    sscanf(argv[3],"%i",&order);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&lambda);
    break;
  }
  case 6:{
    sscanf(argv[2],"%i",&mode);
    sscanf(argv[3],"%i",&order);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&lambda);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&omega);
    break;
  }
  case 7:{
    sscanf(argv[2],"%i",&mode);
    sscanf(argv[3],"%i",&order);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&lambda);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&omega);
    sscanf(argv[6],"%i",&maxiter);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 [mode, %i] [order, %i] [lambda, %g] [omega, %g] [maxiter, %i\n",
	    argv[0],mode,order,lambda,omega,maxiter);
    fprintf(stderr,"\treads a discrete layer spherical harmonics model\n");
    fprintf(stderr,"\tand converts it into a different base function model of order \"order\"\n\n");
    fprintf(stderr,"\tmode=%i: Chebyshev polynomials (normalized)\n",CHEBYSHEV);
    fprintf(stderr,"\tmode=%i: splines with equally spaced knots\n\n",
	    SPLINES);
    fprintf(stderr,"\tlambda and omega are the norm/roughness damping factor as a fraction of the RMS signal\n");
    fprintf(stderr,"\twrites to \"%s\" the power of each base function sorted by order\n",
	    PWR_OUT_FILE);
    fprintf(stderr,"\tmaxiter: if set to > 1, will try to reduce the damping maiter times until vr: %g is met\n",
	    targetvr);
    exit(-1);
    break;
  }}
  // read in model
  read_she_model(argv[1],model,-1,1,expect_gsh);
  // allocate memory for new model
  model[1].radial_type=mode;
  model[1].n=order;
  lms= (int)((((COMP_PRECISION)model[0].lmax)+1.0)
	     *(((COMP_PRECISION)model[0].lmax)+2)/2.0);
  model[1].a=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*model[1].n);
  model[1].b=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*model[1].n);
  if(!model[1].a || !model[1].b)MEMERROR;
  for(i=0;i<model[1].n;i++){
    model[1].a[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lms);
    model[1].b[i]=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*lms);
    if(!model[1].a[i] || !model[1].b[i])MEMERROR;
  }
  fprintf(stderr,"%s: converting to type %s model of order %i, damping: %g/%g times RMS\n",
	  argv[0],(mode==CHEBYSHEV?"Chebyshev":"splines"),
	  model[1].n,lambda,omega);
  model[1].lmax=model[0].lmax;
  // mantle models in km
  model[1].dmin=0;
  model[1].dmax=2871.0;
  fprintf(stderr,"%s: adjusting boundaries to %g and %g\n",
	  argv[0],model[1].dmin,model[1].dmax);
  iter=hit=0;
  do{
    iter++;
    if(iter > 1){
      if(lambda<0)
	lambda=0;
      if(omega<0)
	omega=0;
      fprintf(stderr,"%s: iteration %2i: lambda: %10g omega: %10g old vr: %10g\n",
	      argv[0],iter,lambda,omega,vr);
      if((omega <= 0)&&(lambda <= 0))
	hit++;
    }
    switch(mode){
    case CHEBYSHEV:{
      fit_base_functions(lambda,omega,&totmsize,&vr,model, 
			 CHEBEV_FUNC,cheb_norm_damping);
      // decrease damping for next possible iteration
      lambda -= 0.15;
      omega  -= 0.15;
      break;
    }
    case SPLINES:{
      if(omega != 0.0){
	fprintf(stderr,"%s: setting roughness damping to zero for splines\n",
		argv[0]);
	omega=0.0;
      }
      fit_base_functions(lambda,omega,&totmsize,&vr,model, 
			 spline_base,spline_norm_damping);
      // for next iteration
      lambda -= 0.15;
      break;
    }
    default:{
      fprintf(stderr,"%s: alternative base function mode %i undefined\n",
	      argv[0],mode);
      exit(-1);
    }}
  }while((vr<targetvr) && (iter < maxiter) && (!hit));
  if((iter==maxiter)&&(vr<targetvr)){
    fprintf(stderr,"%s: max number of iterations (%i) reached, vr target (%g) missed\n",
	    argv[0],maxiter,targetvr);
  }
  out=myopen(PWR_OUT_FILE,"w");
  fprintf(out,"# tot_s^2: %g vr: %g\n",totmsize,vr);
  for(i=0;i<model[1].n;i++)
    fprintf(out,"%i %g\n",
	    i,calc_rms(model[1].a[i],model[1].b[i],
		       model[1].lmax));
  fclose(out);
  fprintf(stderr,"%s: written fit values and coeffifient \"power\" to \"%s\"\n",
	  argv[0],PWR_OUT_FILE);
  fprintf(stderr,"%s: tot model size^2: %g model vr: %g\n",
	  argv[0],totmsize,vr);
  write_model((model+1),1.0,1,stdout);
  return 0;
}



void fit_base_functions(COMP_PRECISION lambda,
			COMP_PRECISION omega,
			COMP_PRECISION *totmodelsize,
			COMP_PRECISION *vr,
			struct mod *model, 
			COMP_PRECISION (*f1)(COMP_PRECISION ,
					     COMP_PRECISION ,
					     COMP_PRECISION *,
					     int ,COMP_PRECISION),
			COMP_PRECISION (*f2)(int,int,int))
{
  int i,j,l,m,nm,n,order;
  COMP_PRECISION *y,*c,*x,modelsize,chi2,datasize,totdatasize;
  n=model[0].n;
  order=model[1].n;
  nm=n+2*order;/* length interpolation vectors plus 2 * 
		  number of parameters for  damping */
  //interpolate at n resorted depth layers
  x= (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION) * n);
  if(!x)
    MEMERROR;
  // depths
  for(j=0,i=n-1;i>=0;i--,j++)
    x[j] = model[0].d[i];
  y=(COMP_PRECISION *)calloc(nm,sizeof(COMP_PRECISION));
  c=(COMP_PRECISION *)calloc(order,sizeof(COMP_PRECISION));
  if(!y || !c)
    MEMERROR;
  *totmodelsize=0.0;
  *vr          =0.0;
  totdatasize  =0.0;
  for(l=0;l<=model[0].lmax;l++)
    for(m=0;m<=l;m++){
      // A coefficients
      // first 0 .. model[0].n-1 values are coefficients at depth
      for(j=0,i=n-1;i>=0;i--,j++)
	y[j]= *(model[0].a[i]+POSLM(l,m));
      // interpolate
      determine_coeff(c,order,x,y,n,model[1].dmin,model[1].dmax,
		      lambda,omega,&modelsize,&datasize,&chi2,
		      f1,f2);
      for(i=0;i < order;i++)
	*(model[1].a[i]+POSLM(l,m)) = c[i];
      if(finite(chi2))
	*vr += chi2;
      *totmodelsize += modelsize;
      totdatasize   += datasize;
      if(m!=0){// B coefficients
	for(j=0,i=n-1;i>=0;i--,j++)
	  y[j]= *(model[0].b[i]+POSLM(l,m));
	determine_coeff(c,order,x,y,n,model[1].dmin,model[1].dmax,
			lambda,omega,&modelsize,&datasize,&chi2,
			f1,f2);
	if(finite(chi2))
	  *vr += chi2;
	*totmodelsize += modelsize;
	totdatasize += datasize;
	for(i=0;i < order;i++){
	  *(model[1].b[i]+POSLM(l,m)) = c[i];
	}
      }else
	for(i=0;i < order;i++)
	  *(model[1].b[i]+POSLM(l,m)) = 0.0;
    }
  *vr = 1.0 - *vr/totdatasize;
  free(c);free(y);free(x);
}
