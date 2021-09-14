/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: model2scatter.c,v 1.3 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

//
// reads two spherical harmonics models  and writes all of their
// coefficients as x y scatter data to stdout (linear regression, e.g.)
//
int main(int argc, char **argv)
{
 
  int i,lmax,l,m,weighted=0;
  COMP_PRECISION dz,zmin,zmax,z,r,w;
  struct mod model[2];
  struct mod out_model[2];
  int expect_gsh =0;

  // boundaries
  zmin=50.0;
  zmax=2850;
  dz=50.0;

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
    sscanf(argv[6],"%i",&weighted);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 file2 [zmin(%g) zmax(%g) dz(%g)] [weighted, 0]\n",
	    argv[0],zmin,zmax,dz);
    fprintf(stderr,"writes all non zero coefficients as x y data for both models to stdout\n");
    fprintf(stderr,"(can be used for linear regression, say). if weighted is set, will weigh by r^2\n");
    exit(-1);
    break;
  }}
  fprintf(stderr,"%s: weighted: %i\n",argv[0],weighted);
  // read in models
  for(i=0;i<2;i++){
    read_she_model(argv[1+i],(model+i),-1,1,expect_gsh);
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }
  // pick smaller lmax 
  lmax=(model[0].lmax > model[1].lmax)?
    (model[1].lmax):(model[0].lmax);
  fprintf(stderr,"%s: lmax models: %i/%i, using min lmax=%i\n",
	  argv[0],model[0].lmax,model[1].lmax,lmax);
  for(z=zmin;z<=zmax+EPS_COMP_PREC;z+=dz){
    r=6371.0-z;

    /* don't allow extrapolation */
    interpolate_she_model(out_model,model,z,lmax,FALSE);
    interpolate_she_model((out_model+1),(model+1),z,lmax,FALSE);
    /*  */
    w=(weighted)?(r*r/40589641.):(1.0);
    for(l=0;l<=lmax;l++)
      for(m=0;m<=l;m++){
	fprintf(stdout,"%12.5e %12.5e\n",
		w*out_model[0].a[0][POSLM(l,m)],
		w*out_model[1].a[0][POSLM(l,m)]);
	if(m!=0)
	  fprintf(stdout,"%12.5e %12.5e\n",
		  w*out_model[0].b[0][POSLM(l,m)],
		  w*out_model[1].b[0][POSLM(l,m)]);
      }
  }
  return 0;
}


