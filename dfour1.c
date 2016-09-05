/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: dfour1.c,v 1.4 2001/05/19 15:51:57 becker Exp $ */
/* fast fourier transform from numerical recipes 

*/
#include <math.h>
#include <stdio.h>
#define SWAP(a,b) tempr=(double)(a);(a)=(b);(b)=(COMP_PRECISION)tempr
#include "precision.h"
#include "spherical_harmonics_functions.h"
void dfour1(COMP_PRECISION *data,unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*((double)data[j])-wi*((double)data[j+1]);
	tempi=wr*((double)data[j+1])+wi*((double)data[j]);
	data[j]=data[i]-(COMP_PRECISION)tempr;
	data[j+1]=data[i+1]-(COMP_PRECISION)tempi;
	data[i] += (COMP_PRECISION)tempr;
	data[i+1] += (COMP_PRECISION)tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP
