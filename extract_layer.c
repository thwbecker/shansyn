/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */

#include "shansyn.h"

int main(int argc, char **argv)
{
  COMP_PRECISION z,**a,**b,dz,zmin,zmax,tmp;
  int lmax,steps,use_r2_weighting,nexp=1,extrapolate=1,i;
  struct mod *model,*out_model;
#ifndef SHANA_EXPECT_GSH
  int expect_gsh = 0 ;
#else
  int expect_gsh = 1;
#endif
  lmax=-1;
  steps=30;
  dz=0.0;
  use_r2_weighting=0;


  switch(argc){
  case 3:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    break;
  }
  case 4:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    break;
  }
  case 5:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    break;
  }
  case 6:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&steps);
    break;
  }
  case 7:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&steps);
    sscanf(argv[6],"%i",&use_r2_weighting);
    break;
  }
  case 8:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&steps);
    sscanf(argv[6],"%i",&use_r2_weighting);
    sscanf(argv[7],"%i",&nexp);
    break;
  }
  case 9:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&steps);
    sscanf(argv[6],"%i",&use_r2_weighting);
    sscanf(argv[7],"%i",&nexp);
    sscanf(argv[8],"%i",&extrapolate);
    break;
  }
  case 10:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&z);
    sscanf(argv[3],"%i",&lmax);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[5],"%i",&steps);
    sscanf(argv[6],"%i",&use_r2_weighting);
    sscanf(argv[7],"%i",&nexp);
    sscanf(argv[8],"%i",&extrapolate);
    sscanf(argv[9],"%i",&expect_gsh);
    break;
  }
   default:{
    fprintf(stderr,"%s model_file depth [lmax, %i] [dz, %g] [steps, %i] [use_r2_w, %i] [nexp, %i] [extrapolate, %i] [gsh, %i]\n",
	    argv[0],lmax,dz,steps,use_r2_weighting,nexp,extrapolate,expect_gsh);
    fprintf(stderr,"\textracts spherical harmonic expansion from model model_file\n\tat depth depth\n");
    fprintf(stderr,"\tif lmax is set >0, will limit output to lmax, else uses original lmax\n");
    fprintf(stderr,"\tif dz is set to != 0, will average from z-dz/2 to z+dz/2 in %i steps\n",steps);
    fprintf(stderr,"\tsteps can be changed from %i as argument five\n",steps);
    fprintf(stderr,"\tif use_r2_w is set to unity, will weight the mean by radius^2\n");
    fprintf(stderr,"\tnexp: number of SHE in a row, normally unity\n");
    fprintf(stderr,"\textrapolate: if set to zero, will not extrapolate\n");
    fprintf(stderr,"\tgsh: if set to unity, will expect GSH expansion format\n");
    exit(-1);
  }}
  model = (struct mod *)malloc(sizeof(struct mod)*nexp);
  out_model = (struct mod *)malloc(sizeof(struct mod)*nexp);

  read_she_model(argv[1],model,-1,nexp,expect_gsh);
  if(lmax > 0){
    if(lmax > model[0].lmax){
      lmax=model[0].lmax;
      fprintf(stderr,"%s: limited by original l_max: %i\n",
	      argv[0],lmax);
    }
  }else
    lmax=model[0].lmax;
  fprintf(stderr,"%s: interpolating %s at z=%g lmax=%i, %i expansions, GSH mode: %i\n",
	  argv[0],argv[1],z,lmax,nexp,expect_gsh);
  
  /* copy over the model parameters  */
  
  for(i=0;i < nexp;i++){
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }
  
  if(dz == 0){
    for(i=0;i < nexp;i++){
      interpolate_she_model((out_model+i),
			    (model+i),z,lmax,extrapolate); /* allowing extrapolation */
    }
  }else{
    dz /= 2.0;
    zmin = z-dz;
    zmax = z+dz;
    if(zmin>zmax){
      tmp=zmin;zmin=zmax;zmax=tmp;dz=-dz;
    }
    fprintf(stderr,"%s: averaging with dz: %g (from %g to %g), in %i steps, rw2: %i\n",
	    argv[0],dz*2,zmin,zmax,steps,use_r2_weighting);
    fprintf(stderr,"%s: that is %g%% of the mantle volume\n",argv[0],
	    (pow(REARTH-zmin,3.0)-pow(REARTH-zmax,3.0))/
	    (pow(REARTH,3.0)-pow(RCMB,3.0))*100.0);
    for(i=0;i < nexp;i++){
      mean_model((out_model+i),(model+i),zmin, zmax,
		 steps,use_r2_weighting);
    }
  }
  write_model_layer(out_model,0,1.0,TRUE,nexp,stdout);
  deallocate_model_coefficients(model);
  deallocate_model_coefficients(out_model);
  free(model);
  free(out_model);
  return 0;			/*  */
}
