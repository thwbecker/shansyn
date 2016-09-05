/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: modmodellmax.c,v 1.5 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"
//
// add additional zero value spherical harmonic coefficients or reduce nominal
// resolution
//
//
int main(int argc, char **argv)
{
 
  int newlmax = 31,nexp = 1;
  int expect_gsh = 0;
  struct mod model[3];
  switch(argc){
  case 2:break;
  case 3:{
    sscanf(argv[2],"%i",&newlmax);
    break;
  }
  case 4:{
    sscanf(argv[2],"%i",&newlmax);
    sscanf(argv[3],"%i",&nexp);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 [new_l_max, %i] [nexp, %i]\n",
	    argv[0],newlmax,nexp);
    fprintf(stderr,"increases or decreases l_max of model\n");
    fprintf(stderr,"if decreased, simple cutoff\n");
    fprintf(stderr,"if increased, adds zeros to model\n");
    exit(-1);
    break;
  }}
  if(nexp > 3){
    fprintf(stderr,"%s: nexp too high at %i, max is 3\n",
	    argv[0],nexp);
    exit(-1);
  }
  if(newlmax<=0){
    fprintf(stderr,"%s: new l_max of %i is not good\n",
	    argv[0],newlmax);
    exit(-1);
  }
  // read in model
  read_she_model(argv[1],model,newlmax,nexp,expect_gsh);
  write_model(model,1.0,nexp,stdout);
  return 0;
}

