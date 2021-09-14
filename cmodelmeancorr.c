/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: cmodelmeancorr.c,v 1.6 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

//
// calculates correlation coefficients for each l
// for two different spherical harmonic expansion
// model files at various depths, the fields are however 
// first averaged over a certain depth interval
//

int main(int argc, char **argv)
{
 
  int lmax,l,steps,cmode,i;
  COMP_PRECISION dz,zmin,zmax,z,z1,z2,
    rms[2],deltaz;
  struct mod model[2],out_model[2];
  
  int expect_gsh = 0;

  zmin=0;
  zmax=2900;
  deltaz=500;
  dz=50;
  cmode=1;
  switch(argc){
  case 3:{
    break;
  }
  case 4:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    break;
  }
  case 5:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
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
    sscanf(argv[6],DATA_FSCAN_FORMAT,&deltaz);
    break;
  }
  case 8:{
    sscanf(argv[3],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&dz);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&deltaz);
    sscanf(argv[7],"%i",&cmode);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 file2 [zmin, %g] [zmax, %g] [dz, %g] [delta z, %g] [cmode, %i]\n",
	    argv[0],zmin,zmax,dz,deltaz,cmode);
    fprintf(stderr,"calculates correlation between fields\n");
    fprintf(stderr,"based on two spherical harmonic model sets in file1 and file2\n");
    fprintf(stderr,"output is:\n\n depth rms_mod_1 rms_mod_2 r_tot r_{l_max=8} r_{l_max=20} r_tot^w r_{l=1} r_{l=2} ... r_{l=lmax} \n\n");
    fprintf(stderr,"the depth will move from zmin (%g) to zmax (%g) in dz (%g) steps,\n",
	    zmin,zmax,dz);
    fprintf(stderr,"while fields will be average from z-deltaz/2 to z+deltaz/2, deltaz/2: %g\n",
	    deltaz/2.0);
    fprintf(stderr,"cmode: 1=linear, Pearson correlation 2=parametric, Spearman rank\n");
    exit(-1);
    break;
  }}
  fprintf(stderr,"%s: %g - %g - %g, using delta z: %g\n",
	  argv[0],zmin,dz,zmax,deltaz);
  if(deltaz<0){
    fprintf(stderr,"%s: delta z should be >= 0, %g\n",
	    argv[0],deltaz);
    exit(-1);
  }
  // read in models
  for(i=0;i<2;i++){
    read_she_model(argv[1+i],(model+i),-1,1,expect_gsh);
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }

  // pick smaller lmax 
  lmax=(model[0].lmax > model[1].lmax)?(model[1].lmax):(model[0].lmax);
  steps = (int)deltaz/20;
  deltaz /= 2.0;
  // depth loop
  for(z=zmin;z<=zmax+1.0e-5;z+=dz){
    z1=z-deltaz;
    if(z1<zmin)
      z1=zmin;
    z2=z+deltaz;
    if(z2>zmax)
      z2=zmax;
    for(i=0;i<2;i++){
      mean_model((out_model+i),(model+i),z1, z2,steps,TRUE);
      rms[i] = calc_rms(out_model[i].a[0],out_model[i].b[0],lmax);
    }
    fprintf(stdout,"%12f %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e ",
	    z1+deltaz/2.,rms[0],rms[1],
	    correlation(out_model[0].a[0],out_model[0].b[0],
			out_model[1].a[0],out_model[1].b[0],-lmax,0,1,cmode),
	    correlation(out_model[0].a[0],out_model[0].b[0],
			out_model[1].a[0],out_model[1].b[0],-((lmax>=8)?(8):(lmax)),0,1,cmode),
	    correlation(out_model[0].a[0],out_model[0].b[0],
			out_model[1].a[0],out_model[1].b[0],-((lmax>=20)?(20):(lmax)),0,1,cmode),
	    weighted_correlation(out_model[0].a[0],out_model[0].b[0],
				 out_model[1].a[0],out_model[1].b[0],-lmax,1,cmode));
    for(l=1;l<=lmax;l++)
      fprintf(stdout,"%12.5e ",
	      correlation(out_model[0].a[0],out_model[0].b[0],
			  out_model[1].a[0],out_model[1].b[0],l,0,1,cmode));
    fprintf(stdout,"\n");
  }
  return 0;
}
