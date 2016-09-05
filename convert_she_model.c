/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* 

converts old model format to new spherical harmonics format


*/
/* $Id: convert_she_model.c,v 1.3 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"

int main(int argc, char **argv)
{
  struct mod model[1];
  COMP_PRECISION fac=1.0;
  int expect_gsh = 0;
  switch(argc){
  case 2:
    break;
  case 3:
    sscanf(argv[2],DATA_FSCAN_FORMAT,&fac);
    break;
  default:
    fprintf(stderr,"%s model_file [fac, %g]\n",argv[0],fac);
    fprintf(stderr,"\tconverts spherical harmonic expansion from model model_file\n");
    fprintf(stderr,"\tto the new format as used by the hc Hager & O'Connell code\n");
    exit(-1);
    break;
  }
  read_she_model(argv[1],model,-1,1,expect_gsh);

  model[0].radial_type = HC_TYPE;
  write_model(model,fac,1,stdout);
  return 0;
}

