/*
  part of the shansyn spherical harmonics package, see COPYRIGHT for license 
  $Id: cmodelpower.c,v 1.13 2009/04/16 00:55:29 becker Exp becker $ 

*/
#include "shansyn.h"
/*

  calculates power per degree at different depths for 
  a spherical harmonic expansion model

*/

int main(int argc, char **argv)
{
 
  int l,lmaxlim,allow_extra;
  COMP_PRECISION dz,zmin,zmax,z,*rs,tmpd,
    weight,ws;
  struct mod model[1],out_model[1];

#ifndef SHANA_EXPECT_GSH
  int expect_gsh = 0 ;
#else
  int expect_gsh = 1;
#endif


  // boundaries
  zmin=0;
  zmax=2850;
  dz=50.0;
  lmaxlim=-1;
  allow_extra = 0;
  switch(argc){
  case 2:{
    break;
  }
  case 5:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    break;
  }
  case 6:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&lmaxlim);
    break;
  }
  case 7:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&lmaxlim);
    sscanf(argv[6],"%i",&allow_extra);
    break;
  }
  case 8:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&lmaxlim);
    sscanf(argv[6],"%i",&allow_extra);
    sscanf(argv[6],"%i",&expect_gsh);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 [zmin(%g) zmax(%g) dz(%g) lmaxlim(%i) allow_extra(%i) expect_gsh(%i)]\n",
	    argv[0],zmin,zmax,dz,lmaxlim,allow_extra,expect_gsh);
    fprintf(stderr,"calculates power per degree and unit area at depth intervals dz from zmin to zmax\n");
    fprintf(stderr,"based on a spherical harmonic model in file1\n");
    fprintf(stderr,"output is:\n\n depth total_power total_rms pwr_{l=0}^2 pwr_{l=1}^2 pwr_{l=2}^2 ... pwr_{l=lmax}^2\n\n");
    fprintf(stderr,"the last line has the r^2 averaged total values, indicated by depth=total\n");
    fprintf(stderr,"will limit total_rms of spherical harmonics to l_{max}=lmaxlim. if lmaxlim==-1, will use lmax\n");
    fprintf(stderr,"note that total_power and rms are in data units while pwr_l^2 is data^2\n");
    fprintf(stderr,"rms is total_power without the l=0 term\n");
    fprintf(stderr,"both RMS and total power are normalized by the surface area of the sphere\n");
    fprintf(stderr,"allow_extra: if set, will allow extrapolation (%i)\n",allow_extra);
    fprintf(stderr,"expect_gsh: expect a GSH format model (%i)\n",expect_gsh);
    exit(-1);
    break;
  }}
  if(dz<=0){ 
    dz=50;
    fprintf(stderr,"%s: adjusting dz to %g\n",argv[0],dz);
  }
  // read in model
  read_she_model(argv[1],model,-1,1,expect_gsh);
  copy_model_par(model,out_model);
  out_model->n = 1;
  allocate_model_coefficients(out_model);

  if(!(rs=(COMP_PRECISION *)
       calloc(model->lmax+3,sizeof(COMP_PRECISION))))MEMERROR;

  // loop through all depths
  if(zmin < 0 || zmax  < 0 || zmax > 6371){
    fprintf(stderr,"%s: negative depths not good, use positive < 6371, assuming km\n",
	    argv[0]);
    exit(-1);
  }
  if(lmaxlim == -1){// no limit in l_max on RMS
    lmaxlim = model->lmax;
  }else{// limit RMS 
    if(lmaxlim < 0){
      fprintf(stderr,"%s: RMS limiting l_max needs to be >0\n",
	      argv[0]);
      exit(-1);
    }else{
      fprintf(stderr,"%s: limiting RMS calculation to l_max: %i\n",
	      argv[0],lmaxlim);
    }
  }
  if(dz == 0)
    dz=1.0;
  ws=0.0;
  for(z=zmin;z <= zmax+EPS_COMP_PREC;z+=dz){
    // weigh by radius^2
    weight= (6371.0-z)/6371.0;
    weight *= weight;
    // get model layer
    interpolate_she_model(out_model,model,z,
			  model->lmax,allow_extra);

    
    // total power 
    tmpd = calc_total_power_model(out_model,lmaxlim,0);
    fprintf(stdout,"%12f %12.5e ",z,tmpd);
    rs[0] += tmpd*weight;
    // RMS
    tmpd=calc_rms_model(out_model,lmaxlim,0);
    fprintf(stdout,"%12.5e ",tmpd);
    rs[1] += tmpd*weight;
    
    // \sigma^2_l
    for(l=0;l <= model->lmax;l++){
      tmpd = degree_power_model(out_model,l,0);
      fprintf(stdout,"%12.5e ",tmpd);
      rs[l+2] += tmpd*weight;
    }
    ws += weight;
    fprintf(stdout,"\n");
  }
  // print out average numbers
  printf("total         %12.5e %12.5e ",rs[0]/ws,rs[1]/ws);
  for(l=0;l <= model->lmax;l++)
    printf("%12.5e ",rs[l+2]/ws);
  printf("\n");
  return 0;
}


