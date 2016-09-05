/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: ones.c,v 1.2 2001/03/08 03:04:49 becker Exp $ */
#include <stdio.h>
float checker(int ,int ,int ,int ,int);

void main(int argc, char *argv[] )
{
  int i,j,n,m,che;
  double val;
  switch(argc){
  case 3:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&m);
    val=1;che=0;
    break;
  }
  case 4:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&m);
    sscanf(argv[3],"%lf",&val);
    che=0;
    break;
  }
  case 5:{
    sscanf(argv[1],"%i",&n);
    sscanf(argv[2],"%i",&m);
    sscanf(argv[3],"%lf",&val);
    sscanf(argv[4],"%i",&che);
    break;
  }
  
  default:{
    fprintf(stderr,"%s n m [val, 1.0] [checker, 0]\n",argv[0]);
    fprintf(stderr,"\tprints out ascii array of dimensions n times n\n");
    fprintf(stderr,"\tand value val. If checker is set, prints checkerboard.\n");
    exit(-1);
  }}
  for(i=0;i<n;i++){
    for(j=0;j<m;j++)
      fprintf(stdout,"%15.10e ",val*checker(i,j,n,m,che));
    fprintf(stdout,"\n");
  }
}
float checker(int i ,int j ,int n ,int m,int che)
{
  float tmp1,tmp2;
  if(che == 0)return 1.0;

  tmp1=4.*(float)i/(float)n+0.75;
  tmp2=4.*(float)j/(float)m+0.75;

  if((int)(tmp1-0.5) >= (int)tmp1){
    if(tmp2 == (int)tmp2)
      return 1.0;
    else
      return 0.0;
  }else{
    if((int)(tmp2-0.5) >= (int)tmp2)
      return 0.0;
    else
      return 1.0;
  }
}
