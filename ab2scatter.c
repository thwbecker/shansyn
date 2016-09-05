/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: ab2scatter.c,v 1.2 2005/10/12 19:16:19 becker Exp becker $ */
#include "shansyn.h"


int main(int argc, char **argv)
{
  COMP_PRECISION *a,*b;
  int lmsize,l,m,lmax;
  if(argc!=1){
    fprintf(stderr,"%s: reads in spherical harmonic expansion and writes AB coefficients (no B_l0 terms) \n\tto stdout\n",
	    argv[0]);
    exit(-1);
  }  
  fscanf(stdin,"%i",&lmax);
  lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
  if((a=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL ||
     (b=(COMP_PRECISION *)calloc(lmsize,sizeof(COMP_PRECISION)))==NULL){
    fprintf(stderr,"%s: memerror, lmax=%i lmsize=%i\n",
	    argv[0],lmax,lmsize); 
    exit(-1);
  }
  for(l=0;l<=lmax;l++)
    for(m=0;m<=l;m++)
      if((fscanf(stdin,TWO_DATA_FSCAN_FORMAT,(a+POSLM(l, m)),
		 (b+POSLM(l, m))))!=2){
	fprintf(stderr,
		"%s: read error, l=%i m=%i\n\n",
		argv[0],l,m);exit(-1);}
  for(l=0;l<=lmax;l++){
    printf("%g\n",a[POSLM(l,0)]);
    for(m=1;m<=l;m++){
      printf("%g\n",a[POSLM(l,m)]);
      printf("%g\n",b[POSLM(l,m)]);
    }
  }
  return 0;
}
