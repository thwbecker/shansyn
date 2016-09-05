/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: cmodellinreg.c,v 1.9 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

//
// calculates correlation best fit slopes 
// for two different spherical harmonic expansion
// model files at various depths
//

int main(int argc, char **argv)
{
 
  int i,j,lmax,mode,wmode,zcnt,expect_gsh;
  COMP_PRECISION dz,zmin,zmax,z,slope8[2],rms[2],
    slope[2],slope20[2],sigslope8[2],sigma0[2]={0.1,0.1},
    sigslope[2],sigslope20[2],sigma[2];
  struct mod model[2];
  struct mod out_model[2];

  expect_gsh = 0;
  
  // boundaries
  zmin=50.0;
  zmax=2850;
  dz=50.0;
  mode= REGRESS_XDETERMINED;
  wmode=0;
  switch(argc){
  case 3:{
    break;
  }
  case 6:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    break;
  }
  case 7:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],"%i",&mode);
    break;
  }
  case 9:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],"%i",&mode);
    sscanf(argv[7],DATA_FSCAN_FORMAT,sigma0);
    sscanf(argv[8],DATA_FSCAN_FORMAT,(sigma0+1));
    break;
  }
  default:{
    fprintf(stderr,"%s file1 file2 [zmin(%g) zmax(%g) dz(%g) mode(%i)] [sigma0x(%g) sigma0y(%g)]\n",
	    argv[0],zmin,zmax,dz,mode,sigma0[0],sigma0[1]);
    fprintf(stderr,"calculates  linear best fit abscissae values (a) and slopes (b) \nat depth intervals dz from zmin to zmax\n");
    fprintf(stderr,"(y=a+bx) based on two spherical harmonic model sets in file1 and file2\n");
    fprintf(stderr,"output is:\n\n depth a_tot a_8 a_20 b_tot b_8 b_20 sig_a_tot sig_a_8 sig_a_20 sig_b_tot sig_b_8 sig_b_20\n\n");
    fprintf(stderr,
	    "mode %i: assume that x(file1) is without errors\n",
	    REGRESS_XDETERMINED);
    fprintf(stderr,
	    "mode %i: assume that x and y have errors\n",
	    REGRESS_ITERATIVE);
    fprintf(stderr,"sigma0x/y: set to zero for no sigma0 x/y (will be scaled for normal lin reg, assumed to be unity for iterative\n");
    fprintf(stderr,"sigma0x/y: if set to != zero, assumed that input data has constant sigma0 x/y\n");
    fprintf(stderr,"all sigma0 values are given as fractions of the RMS of the `data' at that depth\n");
    exit(-1);
    break;
  }}
  // read in models
  for(i=0;i < 2;i++){
    read_she_model(argv[1+i],(model+i),-1,1,expect_gsh);
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }
  // pick smaller lmax 
  lmax=(model[0].lmax > model[1].lmax)?(model[1].lmax):(model[0].lmax);

  fprintf(stderr,"%s: lmax models: %i/%i, using min lmax=%i\n",
	  argv[0],model[0].lmax,model[1].lmax,lmax);
  if(mode == REGRESS_XDETERMINED)
    fprintf(stderr,"%s: assuming first model without errors\n",argv[0]);
  else
    fprintf(stderr,"%s: assuming both models with errors\n",argv[0]);
  if(sigma0[0] == 0.0)
    wmode=0;
  else
    wmode=1;
  fprintf(stderr,"%s: weight mode: %i, sigma01: %g, sigma02: %g\n",
	  argv[0],wmode,sigma0[0],sigma0[1]);

  if(dz==0)
    dz=1;

  // loop through all depths
  zcnt=(int)((zmax-zmin)/dz)+1;
  if(zcnt<=0){
    fprintf(stderr,"%s: bounds error, zmin: %g dz: %g zmax: %g\n",
	    argv[0],zmin,dz,zmax);
    exit(-1);
  }
  for(z=zmin,i=0;i<zcnt;i++,z+=dz){
    for(j=0;j<2;j++){
      interpolate_she_model((out_model+j),(model+j),
			    z,lmax,FALSE);
      if(out_model[j].is_gsh){
	fprintf(stderr,"%s: ERROR: GSH not implemented yet\n",argv[0]);
	exit(-1);
      }
      rms[j] = calc_rms(*out_model[j].a,*out_model[j].b,lmax);
      // scale sigma by RMS
      sigma[j]=sigma0[j]*rms[j];
    }
    sphex_lin_reg(slope,sigslope,
		  *out_model[0].a,*out_model[0].b,
		  *out_model[1].a,*out_model[1].b,
		  -lmax,mode,wmode,sigma);
    sphex_lin_reg(slope8,sigslope8,
		  *out_model[0].a,*out_model[0].b,
		  *out_model[1].a,*out_model[1].b,
		  -((lmax>=8)?(8):(lmax)),mode,wmode,sigma);
    sphex_lin_reg(slope20,sigslope20,
		  *out_model[0].a,*out_model[0].b,
		  *out_model[1].a,*out_model[1].b,
		  -((lmax>=20)?(20):(lmax)),mode,wmode,sigma);
    fprintf(stdout,"%12f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    z,
	    slope[0],slope8[0],slope20[0],
	    slope[1],slope8[1],slope20[1],
	    sigslope[0],sigslope8[0],sigslope20[0],
	    sigslope[1],sigslope8[1],sigslope20[1]);
  }
  return 0;
}



