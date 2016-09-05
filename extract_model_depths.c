/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: extract_model_depths.c,v 1.5 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"
//
// extract depth layers of a spherical harmonic model
//
int main(int argc, char **argv)
{
 
  int i,nexp=1;
#ifndef SHANA_EXPECT_GSH
  int expect_gsh = 0 ;
#else
  int expect_gsh = 1;
#endif

  struct mod model[3];
  switch(argc){
  case 2:
    break;
  case 3:
    sscanf(argv[2],"%i",&nexp);
    break;
  case 4:
    sscanf(argv[2],"%i",&expect_gsh);
    break;
  default:{
    fprintf(stderr,"%s file1 [nexp, %i] [expect_gsh, %i]\n",argv[0],nexp,expect_gsh);
    fprintf(stderr,"writes depth layers of model to stdout\n");
    exit(-1);
    break;
  }}
  if(nexp > 3){
    fprintf(stderr,"%s: error, too many expansions: %i (<3)\n",
	    argv[0],nexp);
    exit(-1);
  }
  read_she_model(argv[1],model,-1,nexp,expect_gsh);
  for(i=0;i < model[0].n;i++)
    printf("%3i %11g\n",i+1,model[0].d[i]);
  return 0;
}

