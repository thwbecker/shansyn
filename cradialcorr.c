/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: cradialcorr.c,v 1.6 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

//
// calculates the radial correlation function
//
// there is also a newer program called calc_radial_corr
//
// $Id: cradialcorr.c,v 1.6 2009/04/16 00:55:29 becker Exp becker $
//

int main(int argc, char **argv)
{
 
  int mode,i,j,cmode;
  COMP_PRECISION zmin,zmax,dz,z1,z2,tmp[3],rdist=0,dmode,z,r,dr,a;
  struct mod model[1],out_model[2];
  int layers[2]={0,0};
  int lmin = 0;
#ifndef SHANA_EXPECT_GSH
  int expect_gsh = 0 ;
#else
  int expect_gsh = 1;
#endif

  zmin = 0;
  zmax = 2850;
  dz = 50.0;
  dmode = 0;			/* will be mode */
  cmode = 1;
  switch(argc){
  case 2:{
    break;
  }
  case 3:{
    sscanf(argv[2],"%lf",&dmode);
    break;
  }
  case 6:{
    sscanf(argv[2],"%lf",&dmode);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    break;
  }
   case 7:{
    sscanf(argv[2],"%lf",&dmode);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],"%i",&expect_gsh);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 [mode, %i] [zmin(%g) zmax(%g) dz(%g) expect_gsh(%i)]\n",
	    argv[0],(int)dmode,zmin,zmax,dz,expect_gsh);
    fprintf(stderr,"calculates the radial correlation function of model file1\n");
    fprintf(stderr,"output is:\nz_1 z_2 r_8 r_20 r_total\n");
    fprintf(stderr,"mode 0: uses only the input model depths\n");
    fprintf(stderr,"mode 1: uses z values from zmin to zmax in dz steps\n");
    fprintf(stderr,"mode<0: find distance diagonal to correlation of -mode\n"); /* this is actually mode 2 internally */
    exit(-1);
    break;
  }}
  /* deal with mode */
  if(dmode < 0){
    rdist = -dmode;
    fprintf(stderr,"%s: finding distance from diagonal to correlation %g\n",argv[0],rdist);
    mode = 2;
  }else{
    mode=(int)dmode;
  }
  read_she_model(argv[1],model,-1,1,expect_gsh);
  for(i=0;i < 2;i++){
    copy_model_par(model,(out_model+i));
    if(mode != 0)
      out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
    if(mode == 0)		/* make copies of the original model,
				   this should work as long as we
				   don't destroy the original (?!) */
      memcpy((out_model+i), model, sizeof(struct mod ));
  }


  switch(mode){
  case 0:{			/* compute matrix at model depth levels */
    if(model->radial_type != DISCRETE){
      fprintf(stderr,"%s: mode 0 only works for discrete input level models\n",
	      argv[0]);
      exit(-1);
    }
    for(i=0;i< model->n;i++)
      for(j=i;j < model->n;j++){
	if(i==j){
	  fprintf(stdout,"%g %g %g %g %g\n",
		  model->d[i],model->d[j],1.0,1.0,1.0);
	}else{
	  layers[0] = i;layers[1] = j;
	  
	  tmp[0] = calc_correlation_model(out_model,
					  -((model->lmax>=8)?(8):(model->lmax)),lmin,cmode,FALSE,layers);
	  tmp[1] = calc_correlation_model(out_model,
					  -((model->lmax>=20)?(20):(model->lmax)),lmin,cmode,FALSE,layers);
	  tmp[2] = calc_correlation_model(out_model,-model->lmax,lmin,cmode,FALSE,layers);
	  
	  fprintf(stdout,"%g %g %g %g %g\n",
		  model->d[i],model->d[j],tmp[0],tmp[1],tmp[2]);
	  fprintf(stdout,"%g %g %g %g %g\n",
		  model->d[j],model->d[i],tmp[0],tmp[1],tmp[2]);
	}
      }
    break;
  }
  case 1:{
    layers[0] = 0;layers[1] = 0;
    /* 
       compute matrix at interpolated values 
    */
    for(z1=zmin;z1 <= zmax+EPS_COMP_PREC;z1+=dz)
      for(z2=z1;z2 <= zmax+EPS_COMP_PREC;z2+=dz)
	if(z1==z2){
	  fprintf(stdout,"%g %g %g %g %g\n",z1,z2,1.0,1.0,1.0);
	}else{
	  interpolate_she_model(out_model,
				model,z1,model->lmax,FALSE);
	  interpolate_she_model((out_model+1),
				model,z2,model->lmax,FALSE);
	  tmp[0]=calc_correlation_model(out_model,
					-((model->lmax>=8)?(8):(model->lmax)),lmin,cmode,FALSE,layers);
	  
	  tmp[1]=calc_correlation_model(out_model,
					-((model->lmax>=8)?(20):(model->lmax)),lmin,cmode,FALSE,layers);
	  tmp[2]=calc_correlation_model(out_model,-model->lmax,lmin,cmode,FALSE,layers);
	  fprintf(stdout,"%g %g %g %g %g\n",z1,z2,tmp[0],tmp[1],tmp[2]);
	  fprintf(stdout,"%g %g %g %g %g\n",z2,z1,tmp[0],tmp[1],tmp[2]);
	}
    break;
  } 
  case 2:{
    /* find distance to a certain correlation */
    layers[0] = 0;layers[1] = 0;

    dr=0.1;
    for(z=zmin;z <= zmax+EPS_COMP_PREC;z+=dz){
      r = 0;tmp[0]=1;z1=zmin;z2=zmin;
      while((tmp[0] > rdist)&&(z1>=zmin)&&(z1<=zmax)&& 
	    (z2>=zmin)&&(z2<=zmax)){
	r += dr;
	a = r * HALF_SQRT_TWO;
	z1 = z - a;
	z2 = z + a;
	interpolate_she_model(out_model,
			      model,z1,model->lmax,FALSE);
	interpolate_she_model((out_model+1),
			      model,z2,model->lmax,FALSE);
	tmp[0] = calc_correlation_model(out_model,-model->lmax,lmin,cmode,FALSE,layers);

      }
      if(tmp[0] <= rdist)
	fprintf(stdout,"%11g %11g\n",z,r);
      else
	fprintf(stdout,"%11g NaN\n",z);
    }
    break;
  }
  default:{
    fprintf(stderr,"%s: mode %i is undefined\n",argv[0],mode);
    exit(-1);
  }}
  return 0;
}


