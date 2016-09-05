/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: four1.c,v 1.2 2001/03/08 03:04:49 becker Exp $ */
#include <math.h>
#ifndef COMP_PRECISION
 #define COMP_PRECISION double
#endif


/* fast fourier transform from numerical recipes 


   changed index notation to from 0 to n-1

*/


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(COMP_PRECISION *data,unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  COMP_PRECISION wtemp,wr,wpr,wpi,wi,theta;
  COMP_PRECISION tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j-1],data[i-1]);
      SWAP(data[j],data[i]);
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
	tempr=wr*data[j-1]-wi*data[j];
	tempi=wr*data[j]+wi*data[j-1];
	data[j-1]=data[i-1]-tempr;
	data[j]=data[i]-tempi;
	data[i-1] += tempr;
	data[i] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
#undef SWAP

