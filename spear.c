/* 
   part of the shansyn spherical harmonics package, see COPYRIGHT for license 
*/
/* $Id: spear.c,v 1.6 2009/04/16 00:55:29 becker Exp becker $ */
#include "shansyn.h"




COMP_PRECISION corr_sub(COMP_PRECISION *x, COMP_PRECISION *y, int n, 
			COMP_PRECISION *prob, int mode)
{
  COMP_PRECISION zd,probd,probrs,corr,d;
  static int remove_avg = 0;	/* typically, this is not removed */
  switch(mode){
  case 1:			/* pearsons (regular) correlation */
    pearsn((x-1),(y-1),(unsigned long)n,&corr,&probd,&zd,remove_avg);
    *prob = 1 - probd;
    break;
  case 2:			/* spearman rank */
    spear((x-1),(y-1),(unsigned long)n,&d,&zd, &probd, &corr,&probrs);
    *prob = 1 - probrs;		/* probability of rejection of null hypothesis, higher is better */
    break;
  default:
    fprintf(stderr,"corr_sub: mode %i undefined\n",mode);
    exit(-1);
  }
  
  
  return(corr);
}


#define TINY 1.0e-20

void pearsn(COMP_PRECISION x[],COMP_PRECISION y[], unsigned long n,
	    COMP_PRECISION *r,COMP_PRECISION *prob, COMP_PRECISION *z,
	    int remove_mean)
{
  unsigned long j;
  COMP_PRECISION yt,xt,t,df;
  COMP_PRECISION syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;
  if(remove_mean){
    for (j=1;j<=n;j++) {
      //fprintf(stderr,"%12g %12g\n",x[j],y[j]);
      ax += x[j];
      ay += y[j];
    }
    ax /= n;
    ay /= n;
  }else{
    ax = 0.0;
    ay = 0.0;
  }
  for (j=1;j<=n;j++) {
    xt=x[j]-ax;
    yt=y[j]-ay;
    //fprintf(stderr,"%12g %12g\n",xt,yt);
    sxx += xt*xt;
    syy += yt*yt;
    sxy += xt*yt;
  }
  *r=sxy/sqrt(sxx*syy);
  *z=0.5*log((1.0+(*r)+TINY)/(1.0-(*r)+TINY));
  df=n-2;
  t=(*r)*sqrt(df/((1.0-(*r)+TINY)*(1.0+(*r)+TINY)));
  *prob=betai(0.5*df,0.5,df/(df+t*t));
}
#undef TINY


void spear(COMP_PRECISION data1[], COMP_PRECISION data2[], unsigned long n, 
	   COMP_PRECISION *d, COMP_PRECISION *zd, 
	   COMP_PRECISION *probd, COMP_PRECISION *rs, COMP_PRECISION *probrs)
{
	unsigned long j;
	COMP_PRECISION vard,t,sg,sf,fac,en3n,en,df,aved,*wksp1,*wksp2;

	wksp1=vector(1,n);
	wksp2=vector(1,n);
	for (j=1;j<=n;j++) {
	  //fprintf(stderr,"%20g %20g\n",data1[j],data2[j]);
		wksp1[j]=data1[j];
		wksp2[j]=data2[j];
	}
	sort2(n,wksp1,wksp2);
	crank(n,wksp1,&sf);
	sort2(n,wksp2,wksp1);
	crank(n,wksp2,&sg);
	*d=0.0;
	for (j=1;j<=n;j++)
		*d += SQUARE(wksp1[j]-wksp2[j]);
	en=n;
	en3n=en*en*en-en;
	aved=en3n/6.0-(sf+sg)/12.0;
	fac=(1.0-sf/en3n)*(1.0-sg/en3n);
	vard=((en-1.0)*en*en*SQUARE(en+1.0)/36.0)*fac;
	*zd=(*d-aved)/sqrt(vard);
	*probd=erfcc(fabs(*zd)/1.4142136);
	*rs=(1.0-(6.0/en3n)*(*d+(sf+sg)/12.0))/sqrt(fac);
	fac=(*rs+1.0)*(1.0-(*rs));
	if (fac > 0.0) {
		t=(*rs)*sqrt((en-2.0)/fac);
		df=en-2.0;
		*probrs=betai(0.5*df,0.5,df/(df+t*t));
	} else
		*probrs=0.0;
	free_vector(wksp2,1,n);
	free_vector(wksp1,1,n);
}



void 
crank (unsigned long n, COMP_PRECISION w[], COMP_PRECISION *s)
{
	unsigned long j=1,ji,jt;
	COMP_PRECISION t,rank;

	*s=0.0;
	while (j < n) {
		if (w[j+1] != w[j]) {
			w[j]=j;
			++j;
		} else {
			for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
			rank=0.5*(j+jt-1);
			for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
			t=jt-j;
			*s += t*t*t-t;
			j=jt;
		}
	}
	if (j == n) w[n]=n;
}



#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void 
sort2 (unsigned long n, COMP_PRECISION arr[], COMP_PRECISION brr[])
{
	unsigned long i,ir=n,j,k,l=1;
	int *istack,jstack=0;
	COMP_PRECISION a,b,temp;

	istack=ivector(1,NSTACK);
	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
				arr[i+1]=a;
				brr[i+1]=b;
			}
			if (!jstack) {
				free_ivector(istack,1,NSTACK);
				return;
			}
			ir=istack[jstack];
			l=istack[jstack-1];
			jstack -= 2;
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			SWAP(brr[k],brr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
				SWAP(brr[l+1],brr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
				SWAP(brr[l],brr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
				SWAP(brr[l+1],brr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			b=brr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
				SWAP(brr[i],brr[j])
			}
			arr[l]=arr[j];
			arr[j]=a;
			brr[l]=brr[j];
			brr[j]=b;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
}
#undef M
#undef NSTACK
#undef SWAP

/* CAUTION: This is the traditional K&R C (only) version of the Numerical
   Recipes utility file nrutil.c.  Do not confuse this file with the
   same-named file nrutil.c that is supplied in the same subdirectory or
   archive as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only traditional K&R.           */





COMP_PRECISION 
erfcc (double x)
{
	COMP_PRECISION t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}

COMP_PRECISION 
betai (double a, double b, double x)
{
	void nrerror();
	COMP_PRECISION bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}



#define MAXIT 1000
#ifdef SINGLE_PREC
#define EPS 3.0e-7
#else
#define EPS 5.0e-15
#endif
#define FPMIN 1.0e-30

COMP_PRECISION 
betacf (double a, double b, double x)
{
	int m,m2;
	COMP_PRECISION aa,c,d,del,h,qab,qam,qap;
	if(!finite(x))
	  return x;
	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) {
	  fprintf(stderr,"a %g b %g x %g\n",a,b,x);
	  nrerror("a or b too big, or MAXIT too small in betacf");
	}
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN
