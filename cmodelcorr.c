/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */

#include "shansyn.h"

//
// calculates correlation coefficients for each l
// for two different spherical harmonic expansion
// model files at various depths
//
// NOTE: there is also cmodelmeancorr for boxcar mean correlations
//

int main(int argc, char **argv)
{
  int i,lmax,l,llim,extrapolate=0,cmode;
  int layers[2] = {0,0};
  COMP_PRECISION dz,zmin,zmax,z,rms[2],
    shift,stretch;
  struct mod model[2],out_model[2];

#ifndef SHANA_EXPECT_GSH
  int expect_gsh = 0 ;
#else
  int expect_gsh = 1;
#endif

  // boundaries
  zmin=50.0;
  zmax=2850;
  dz=50.0;
  // shift of model 2 with depth
  shift = 0.0;
  // stretch of depth layers of model 2 
  stretch=1.0;
  // llim, if set to -1, will use all coefficients, ie. l_max
  llim=-1;
  /* default type of correlation */
  cmode = 1;			/* 1: linear 2: parametric */

  switch(argc){
  case 3:
    break;
  case 6:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    break;
  case 7:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    break;
  case 8:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&stretch);
    break;
  case 9:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&stretch);
    sscanf(argv[8],"%i",&llim);
    break;
  case 10:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&stretch);
    sscanf(argv[8],"%i",&llim);
    sscanf(argv[9],"%i",&extrapolate);
    break;
  case 11:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&stretch);
    sscanf(argv[8],"%i",&llim);
    sscanf(argv[9],"%i",&extrapolate);
    sscanf(argv[10],"%i",&cmode);
    break;
  case 12:
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&shift);
    sscanf(argv[7],DATA_FSCAN_FORMAT,&stretch);
    sscanf(argv[8],"%i",&llim);
    sscanf(argv[9],"%i",&extrapolate);
    sscanf(argv[10],"%i",&cmode);
    sscanf(argv[11],"%i",&expect_gsh);
    break;
  default:
    fprintf(stderr,"%s file1 file2 [zmin(%g) zmax(%g) dz(%g)] [shift, %g] [stretch, %g] [L, %i] [extrapolate, %i] [cmode, %i] [expect_gsh, %i]\n",
	    argv[0],zmin,zmax,dz,shift,stretch,llim,extrapolate,cmode,expect_gsh);
    fprintf(stderr,"calculates correlation at depth intervals dz from zmin to zmax\n");
    fprintf(stderr,"based on two spherical harmonic model sets in file1 and file2\n");
    fprintf(stderr,"output is:\n\n depth RMS_mod_1 RMS_mod_2 r_{tot/L} r_{l_max=8} r_{l_max=20} r_tot^w r_{l=1} r_{l=2} ... r_{l=lmax} \n\n");
    fprintf(stderr,"if shift is set to non-zero, will shift all z's for the second model by shift in depth\n");
    fprintf(stderr,"if stretch is set to non unity, will stretch depth levels of model 2 by strech\n");
    fprintf(stderr,"if L is set to a positive number, will limit r_{tot} to r_L, else r_{tot} will be up to the l_max of the models\n");
    fprintf(stderr,"if extrapolate = 1, will extrapolate models, else not\n");
    fprintf(stderr,"cmode: 1=linear, Pearson correlation 2=parametric, Spearman rank\n");
    exit(-1);
    break;
  }
  if(dz <= 0){ 
    dz=50;
    fprintf(stderr,"%s: adjusting dz to %g\n",argv[0],dz);
  }
  // read in models
  for(i=0;i < 2;i++){
    read_she_model(argv[1+i],(model+i),-1,1,expect_gsh);
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }
  // pick smaller lmax 
  lmax=(model[0].lmax > model[1].lmax)?(model[1].lmax):(model[0].lmax);
  // loop through all depths
  if(shift != 0.0){
    fprintf(stderr,"%s: WARNING: shift is non-zero, but %g\n",
	    argv[0],shift);
    
  }
  // stretch depth layers
  if(stretch != 1.0){
    fprintf(stderr,"%s: WARNING: stretch is non-unity, but %g\n",
	    argv[0],stretch);
    if(stretch < 1.0){
      zmax = model[1].d[0] * stretch;
      fprintf(stderr,"%s: limiting zmax to %g\n",
	      argv[0],zmax);
    }
    for(i=0;i<model[1].n;i++)
      model[1].d[i] *= stretch;
  }
  if(llim >= 1)
    fprintf(stderr,"%s: using l_max=%i for r_{tot}, not %i\n",
	    argv[0],llim,lmax);
  for(z=zmin;z <= zmax+EPS_COMP_PREC;z+=dz){
    /* no extrapolation */
    // model A
    interpolate_she_model(out_model,model,z,
			  lmax,extrapolate);
    rms[0]=calc_rms_model(out_model,lmax,0);
    // model B
    interpolate_she_model((out_model+1),(model+1),z+shift,
			  lmax,extrapolate);
    rms[1]=calc_rms_model((out_model+1),lmax,0);
    /* 
       1      2    3       4   5   6     7 
       depth rms1 rms2  r_lim r_8 r_20 r_total 
    */
    fprintf(stdout,"%12f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
	    z,rms[0],rms[1],
	    calc_correlation_model(out_model,-((llim >=   1)?(llim):(lmax)),0,cmode,FALSE,layers),
	    calc_correlation_model(out_model,-((lmax >=   8)?(8):(lmax)),   0,cmode,FALSE,layers),
	    calc_correlation_model(out_model,-((lmax >=  20)?(20):(lmax)),  0,cmode,FALSE,layers),
	    calc_correlation_model(out_model,-lmax,                         0,cmode,TRUE,layers));
    for(l=1;l<=lmax;l++)	/*  8    9    10   11
				   r(1) r(2) r(3) r(4) 
				*/
      fprintf(stdout,"%12.5e ",calc_correlation_model(out_model,l,0,cmode,FALSE,layers));
    fprintf(stdout,"\n");
  }
  return 0;
}



