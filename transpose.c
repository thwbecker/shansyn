/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: transpose.c,v 1.2 2001/03/08 03:04:49 becker Exp $ */
#include <stdio.h>
#ifndef FORMAT
 #define FORMAT "%20.10e\n"
#endif
void main(void)
{
  double t;
  while(fscanf(stdin,"%lf",&t)==1)
    fprintf(stdout,FORMAT,t);
}
