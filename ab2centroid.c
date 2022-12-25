/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: ab2centroid.c,v 1.3 2005/04/05 19:06:14 becker Exp becker $ */
#include "shansyn.h"


int main(int argc, char **argv)
{
  COMP_PRECISION phi,theta,*a,*b,area;
  int lmsize,l,m,lmax,rc;
  if(argc!=1){
    fprintf(stderr,"%s: reads in spherical harmonic expansion of 1/0 field and returns area and centroid coordinates\n\tto stdout\n",
	    argv[0]);
    fprintf(stderr,"\tinput is in physical normalization, output is area (4pi is max) lon lat in degrees\n");
    exit(-1);
  }  
  rc = fscanf(stdin,"%i",&lmax);
  if(!rc){
    fprintf(stderr,"%s: read error header\n",argv[0]);
    exit(-1);
  }
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
  // area 
  area=2.0*a[POSLM(0,0)]*sqrt(PI);
  // center of mass
  theta=((4.*a[POSLM(0,0)] - sqrt(3.0)*a[POSLM(1,0)])*pow(PI,1.5))/4.;
  theta/=area;
  phi=  ((4.*a[POSLM(0,0)] + sqrt(3.0)*b[POSLM(1,1)])*pow(PI,1.5))/2.;
  phi/=area;

  printf("%g %g %g\n",area,PHI2LONGITUDE(phi),THETA2LATITUDE(theta));
  return 0;
}
