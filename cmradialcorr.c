/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: cmradialcorr.c,v 1.1 2011/04/01 20:43:39 becker Exp becker $ */
#include "shansyn.h"


//
// calculates the cross-model radial correlation function
//
//
// $Id: cmradialcorr.c,v 1.1 2011/04/01 20:43:39 becker Exp becker $
//

int main(int argc, char **argv)
{
 
  int mode,i,cmode,lmax_min;
  COMP_PRECISION zmin,zmax,dz,z1,z2,tmp[3],rdist,dmode;
  struct mod model[2],out_model[2];
  int expect_gsh = 0;
  zmin=0;
  zmax=2850;
  dz=50.0;
  dmode = 1;			/* will be mode */
  cmode = 1;
  switch(argc){
  case 3:{
    break;
  }
  case 4:{
    sscanf(argv[3],"%lf",&dmode);
    break;
  }
  case 5:{
    sscanf(argv[3],"%lf",&dmode);
    sscanf(argv[4],DATA_FSCAN_FORMAT,&zmin);
    sscanf(argv[5],DATA_FSCAN_FORMAT,&zmax);
    sscanf(argv[6],DATA_FSCAN_FORMAT,&dz);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 file2 [mode, %i] [zmin(%g) zmax(%g) dz(%g)]\n",
	    argv[0],(int)dmode,zmin,zmax,dz);
    fprintf(stderr,"calculates the cross-model radial correlation function of models file1 and file2\n");
    fprintf(stderr,"output is:\nz_1 z_2 r_8 r_20 r_total\n");
    fprintf(stderr,"mode 1: uses z values from zmin to zmax in dz steps\n");
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
  for(i=0;i<2;i++){
    read_she_model(argv[i+1],(model+i),    -1,1,expect_gsh);
    copy_model_par((model+i),(out_model+i));
    out_model[i].n = 1;
    allocate_model_coefficients((out_model+i));
  }
  lmax_min = (model[0].lmax < model[1].lmax)?(model[0].lmax):(model[1].lmax);
  fprintf(stderr,"%s: using %s and %s, lmax %i mode %i\n",argv[0],argv[1],argv[2],lmax_min,mode);


  switch(mode){
  case 1:{
    for(z1=zmin;z1<=zmax+EPS_COMP_PREC;z1+=dz)
      for(z2=z1;z2<=zmax+EPS_COMP_PREC;z2+=dz){
	interpolate_she_model(out_model,    model,    z1,model[0].lmax,FALSE);
	interpolate_she_model((out_model+1),(model+1),z2,model[1].lmax,FALSE);

	tmp[0]=correlation(out_model[0].a[0],out_model[0].b[0],out_model[1].a[0],out_model[1].b[0],-((lmax_min >=  8)?(8): (lmax_min)),0,1,cmode);
	tmp[1]=correlation(out_model[0].a[0],out_model[0].b[0],out_model[1].a[0],out_model[1].b[0],-((lmax_min >= 20)?(20):(lmax_min)),0,1,cmode);
	tmp[2]=correlation(out_model[0].a[0],out_model[0].b[0],out_model[1].a[0],out_model[1].b[0],-lmax_min,0,1,cmode);

	if(z1==z2)
	  fprintf(stdout,"%g %g %g %g %g\n",z1,z2,tmp[0],tmp[1],tmp[2]);
	else{
	  fprintf(stdout,"%g %g %g %g %g\n",z1,z2,tmp[0],tmp[1],tmp[2]);
	  fprintf(stdout,"%g %g %g %g %g\n",z2,z1,tmp[0],tmp[1],tmp[2]);
	}
      }
    break;
  } 
   default:{
    fprintf(stderr,"%s: mode %i is undefined\n",argv[0],mode);
    exit(-1);
  }}
  return 0;
}


