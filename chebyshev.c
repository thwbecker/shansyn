/* part of the shansyn spherical harmonics package, see COPYRIGHT for license */
/* $Id: chebyshev.c,v 1.7 2001/03/15 22:50:39 becker Exp becker $ */
#include "shansyn.h"

//
// Chebyshev base functions
// $Id: chebyshev.c,v 1.7 2001/03/15 22:50:39 becker Exp becker $ 
//


/*
  evaluate Chebyshev polynomials at a <= x <= b 
  given coefficients c[0]...c[m-1]

*/
COMP_PRECISION chebev(COMP_PRECISION a,COMP_PRECISION b,
		      COMP_PRECISION *c,int m,COMP_PRECISION x)
{
  COMP_PRECISION d=0.0,dd=0.0,sv,y,y2;
  int j;
  if ((x-a)*(x-b) > 0.0){
    fprintf(stderr,"chebev: x (%g) not in range (%g - %g)\n",
	    x,a,b);
    exit(-1);
  }
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>=1;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}
// normalized
COMP_PRECISION norm_chebev(COMP_PRECISION a,COMP_PRECISION b,
			   COMP_PRECISION *c,int m,
			   COMP_PRECISION x)
{
  COMP_PRECISION d=0.0,dd=0.0,sv,y,y2,tmp;
  int j;
  if ((x-a)*(x-b) > 0.0){
    fprintf(stderr,"chebev: x (%g) not in range (%g - %g)\n",
	    x,a,b);
    exit(-1);
  }
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>=1;j--) {
    sv=d;
    d=y2*d-dd+c[j]*chebnorm(j);
    dd=sv;
  }
  tmp=y*d-dd+0.5*c[0]*chebnorm(0);

  return tmp;
}

/*
  
  normalization factor sqrt((4k^2-1)/(4k^2-2))

*/

COMP_PRECISION chebnorm(int n)
{
  static COMP_PRECISION factor[51]={0.70710678118654752, 1.2247448713915890, 1.0350983390135313, 
				    1.0145993123917847, 1.0080322575483706, 1.0050890913907349, 
				    1.0035149493261806, 1.0025740068320432, 1.0019665702377579, 
				    1.0015515913132542, 1.0012554932753529, 1.0010368069140517, 
				    1.0008707010791882, 1.0007415648034324, 1.0006391819124997, 
				    1.0005566379501475, 1.0004891171728023, 1.0004331817400483, 
				    1.0003863241403532, 1.0003466805443029, 1.0003128421787780, 
				    1.0002837281941049, 1.0002584981302063, 1.0002364904845643, 
				    1.0002171788493409, 1.0002001401000727, 1.0001850309942748, 
				    1.0001715707312960, 1.0001595277987336, 1.0001487099420817, 
				    1.0001389564378277, 1.0001301320846162, 1.0001221224893116, 
				    1.0001148303385539, 1.0001081724271648, 1.0001020772727612, 
				    1.0000964831880291, 1.0000913367129822, 1.0000865913323700, 
				    1.0000822064204603, 1.0000781463682668, 1.0000743798580373, 
				    1.0000708792572744, 1.0000676201102981, 1.0000645807098093, 
				    1.0000617417343842, 1.0000590859405520, 1.0000565979002632, 
				    1.0000542637762561, 1.0000520711291977, 1.0000500087515628};
  COMP_PRECISION tmp;
  if(n<50)
    return factor[n];
  else{
    tmp=2.0*(COMP_PRECISION)n;
    tmp *= tmp;
    return sqrt((tmp-1.0)/(tmp-2.0));
  }
}

//
// norm damping function 
// effect of base function j on base function i out of n
// normally, \delta_ij 
COMP_PRECISION cheb_norm_damping(int i,int j,int n)
{
  // best results for simple norm damping
  // for all modes but the zero (constant) polynomial
  if(i!=j)
    return 0.0;
  else{
    return (i==0)?(0.0):(1.0);// don't damp constant term
    //return((COMP_PRECISION)(i)/(COMP_PRECISION)n);
  }
}

COMP_PRECISION chebev_der_int(COMP_PRECISION x1,
			      COMP_PRECISION x2,int i)
{
  static COMP_PRECISION 
    f[101]={0., 1.41421, 3.26599, 5.25357, 7.3238, 9.45331, 11.6289, 
	    13.8421, 16.087, 18.3592, 20.6555, 22.9733, 25.3103, 27.665, 
	    30.0357, 32.4213, 34.8207, 37.233, 39.6573, 42.0929, 44.5392, 
	    46.9956, 49.4616, 51.9366, 54.4204, 56.9125, 59.4125, 61.9202, 
	    64.4351, 66.9571, 69.4858, 72.0211, 74.5627, 77.1103, 79.6639, 
	    82.2231, 84.7879, 87.3581, 89.9334, 92.5139, 95.0992, 97.6894, 
	    100.284, 102.884, 105.487, 108.096, 110.708, 113.324, 115.945, 
	    118.569, 121.198, 123.83, 126.466, 129.105, 131.748, 134.394, 
	    137.044, 139.697, 142.354, 145.014, 147.676, 150.342, 153.011, 
	    155.683, 158.358, 161.036, 163.717, 166.4, 169.087, 171.776, 
	    174.467, 177.162, 179.859, 182.558, 185.26, 187.965, 190.671, 
	    193.381, 196.092, 198.807, 201.523, 204.242, 206.962, 209.686, 
	    212.411, 215.138, 217.868, 220.6, 223.334, 226.07, 228.808, 
	    231.548, 234.29, 237.034, 239.779, 242.527, 245.277, 248.029, 
	    250.782, 253.537, 256.294};
  if(i>100)
    return (COMP_PRECISION)i * 2.7525 -18.9621;
  else
    return f[i];
}
// for normalized
COMP_PRECISION norm_chebev_der_int(COMP_PRECISION x1,COMP_PRECISION x2,int i)
{
  static COMP_PRECISION f[51]={0., 1.73205, 3.38062, 5.33027, 7.38263, 9.50142, 11.6698, 13.8777, 16.1186, 
			       18.3877, 20.6815, 22.9971, 25.3324, 27.6855, 30.0549, 32.4394, 34.8378,
			       37.2491, 39.6726, 42.1075, 44.5531, 47.0089, 49.4743, 51.9489, 54.4322,
			       56.9239, 59.4235, 61.9308, 64.4454, 66.967, 69.4955, 72.0305, 74.5718, 		
			       77.1192, 79.6725, 82.2315, 84.7961, 87.3661, 89.9412, 92.5215, 95.1067, 
			       97.6967, 100.291, 102.891, 105.494, 108.102, 110.714, 113.331, 115.951, 
			       118.576, 121.204};
  if(i>50)
    return (COMP_PRECISION)i * 2.6225 -9.9248;
  else
    return f[i];
}
