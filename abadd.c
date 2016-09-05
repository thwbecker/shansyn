/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: abadd.c,v 1.5 2009/04/16 00:55:29 becker Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "precision.h"

#define POSLM(l, m) ( (((l) * (1+(l))) / 2) + (m))

int main(int argc, char *argv[] )
{
  int lmax,l,m,lmsize,lmax2,n;
  double *a,*b,atmp,btmp;
 

  if(argc > 1 && strcmp(argv[1],"-h")==0 ){
    fprintf(stderr,"%s\n",argv[0]);
    fprintf(stderr,"\t reads sequentially spherical coefficents A B  from stdin\n");
    fprintf(stderr,"\t (produced by shana) and adds them up to stdout\n");
    fprintf(stderr,"\n\n");
    exit(-1);
  }
  fscanf(stdin,"%i",&lmax);

  lmsize= (int)((((float)lmax)+1.0)*(((float)lmax)+2)/2.0);
  a=(double *)calloc(lmsize,sizeof(double));
  b=(double *)calloc(lmsize,sizeof(double));

  for(l=0;l<=lmax;l++)
    for(m=0;m<=l;m++)
      if((fscanf(stdin,"%lf %lf",(a+POSLM(l, m)),(b+POSLM(l, m))))!=2){
	fprintf(stderr,"%s: read error, l=%i m=%i\n\n",argv[0],l,m);exit(-1);}
  
  n=1;
  while(fscanf(stdin,"%i",&lmax2) == 1){
    if(lmax2 != lmax){
      fprintf(stderr,"%s: need identical lmax (%i) in sequential AB files\n",argv[0],lmax);
      exit(-1);
    }
    for(l=0;l<=lmax;l++)
      for(m=0;m<=l;m++){
	if((fscanf(stdin,"%lf %lf",&atmp,&btmp))!=2){
	  fprintf(stderr,"%s: read error, l=%i m=%i AB file %i\n",argv[0],l,m,n+1);exit(-1);}
	else{
	  *(a+POSLM(l, m)) += atmp;
	  *(b+POSLM(l, m)) += btmp;
	}
      }
    n++;
  }
  if(n > 1)
    fprintf(stderr,"%s: added %i files, lmax=%i\n",argv[0],n,lmax);
  else
    fprintf(stderr,"%s: one file streamed through, lmax=%i\n",argv[0],lmax);
  
  fprintf(stdout,"%i\n",lmax);
  for(l=0;l<=lmax;l++)
    for(m=0;m<=l;m++){
      fprintf(stdout,COEFF_DATA_FORMAT,
	      *(a+POSLM(l, m)),*(b+POSLM(l, m)));
      fprintf(stdout,"\n");
    }

  return 0;
}
