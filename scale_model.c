/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: scale_model.c,v 1.4 2002/01/07 22:09:26 becker Exp becker $ */
#include "shansyn.h"

//
// scales model with factor 
//

int main(int argc, char **argv)
{
  COMP_PRECISION factor=1.0;
  int nexp=1;
  int expect_gsh = 0;
  struct mod model[3];
  switch(argc){
  case 2:{
    break;
  }
  case 3:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&factor);
    break;
  }
  case 4:{
    sscanf(argv[2],DATA_FSCAN_FORMAT,&factor);
    sscanf(argv[3],"%i",&nexp);
    break;
  }
  default:{
    fprintf(stderr,"%s file1 [factor (%g)] [nexp, %i]\n",
	    argv[0],factor,nexp);
    fprintf(stderr,"scales spherical harmonic model by factor\n");
    fprintf(stderr,"nexp is the number of expansions in a row\n");
    exit(-1);
    break;
  }}
  if(nexp > 3){
    fprintf(stderr,"%s: error, too many expansions: %i (<3)\n",
	    argv[0],nexp);
    exit(-1);
  }
  read_she_model(argv[1],model,-1,nexp,expect_gsh);
  if(factor != 1.0){
    fprintf(stderr,"%s: scaling model with factor %g\n",argv[0],factor);
  }else{
    fprintf(stderr,"%s: no scaling, factor: %g\n",argv[0],factor);
  }
  write_model(model,factor,nexp,stdout);
  return 0;
}


