#include "shansyn.h"

/*

  adapted from Numerical Recipes in C, p.686

  iterative linear regression routine


  $Id: fitxy.c,v 1.5 2009/04/16 00:55:29 becker Exp becker $
*/


void linreg_fit(COMP_PRECISION *x,COMP_PRECISION *y,int ndata,
		COMP_PRECISION *sig,int mwt,
		COMP_PRECISION *a,COMP_PRECISION *b,
		COMP_PRECISION *siga,COMP_PRECISION *sigb,
		COMP_PRECISION *chi2,COMP_PRECISION *q)
{
  int i;
  COMP_PRECISION  wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat,
    tmpd;
  
  *b=0.0;
  if (mwt) {
    ss=0.0;
    for (i=1;i<=ndata;i++) {
      wt=1.0/SQUARE(sig[i]);
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
  } else {
    for (i=1;i<=ndata;i++) {
      sx += x[i];
      sy += y[i];
    }
    ss=ndata;
  }
  sxoss=sx/ss;
  if (mwt) {
    for (i=1;i<=ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      *b += t*y[i]/sig[i];
    }
  } else {
    for (i=1;i<=ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      *b += t*y[i];
    }
  }
  *b /= st2;
  *a=(sy-sx*(*b))/ss;
  *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
  *sigb=sqrt(1.0/st2);
  *chi2=0.0;
  if (mwt == 0) {
    for (i=1;i<=ndata;i++){
      tmpd = y[i]-(*a)-(*b)*x[i];
      *chi2 += SQUARE(tmpd);
    }
    *q=1.0;
    sigdat=sqrt((*chi2)/(ndata-2));
    *siga *= sigdat;
    *sigb *= sigdat;
  } else {
    for (i=1;i<=ndata;i++){
      tmpd = (y[i]-(*a)-(*b)*x[i])/sig[i];
      *chi2 += SQUARE(tmpd);
    }
    *q=gammq(0.5*(ndata-2),0.5*(*chi2));
  }
}



/* 

   iterative routine

*/



/*

  these are global variables for the chi2 type 
  routines

*/
int nn;
COMP_PRECISION *xx,*yy,*sx,*sy,*ww,aa,offs;


#define POTN 1.571000
#define BIG 1.0e30
#ifndef PI
#define PI 3.14159265358979324
#endif
#define ACC 1.0e-4


void linreg_fitexy(COMP_PRECISION *x,COMP_PRECISION *y,
		   int ndat,
		   COMP_PRECISION *sigx,COMP_PRECISION *sigy,
		   COMP_PRECISION *a,COMP_PRECISION *b,
		   COMP_PRECISION *siga,COMP_PRECISION *sigb,
		   COMP_PRECISION *chi2,COMP_PRECISION *q)
{
  int j;
  COMP_PRECISION swap,amx,amn,varx,vary,ang[7],ch[7],
    scale,bmn,bmx,d1,d2,r2,tmpd,
    dum1,dum2,dum3,dum4,dum5;

  xx=vector(1,ndat);
  yy=vector(1,ndat);
  sx=vector(1,ndat);
  sy=vector(1,ndat);
  ww=vector(1,ndat);
  avevar(x,ndat,&dum1,&varx);
  avevar(y,ndat,&dum1,&vary);
  scale=sqrt(varx/vary);
  nn=ndat;
  for (j=1;j<=ndat;j++) {
    xx[j]=x[j];
    yy[j]=y[j]*scale;
    sx[j]=sigx[j];
    sy[j]=sigy[j]*scale;
    ww[j]=sqrt(SQUARE(sx[j])+SQUARE(sy[j]));
  }
  fit(xx,yy,nn,ww,1,&dum1,b,&dum2,&dum3,&dum4,&dum5);
  offs=ang[1]=0.0;
  ang[2]=atan(*b);
  ang[4]=0.0;
  ang[5]=ang[2];
  ang[6]=POTN;
  for (j=4;j<=6;j++) ch[j]=chixy(ang[j]);
  mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],
	 (COMP_PRECISION (*)())chixy);
  *chi2=brent(ang[1],ang[2],ang[3],(COMP_PRECISION (*)())chixy,ACC,b);
  *chi2=chixy(*b);
  *a=aa;
  *q=gammq(0.5*(nn-2),*chi2*0.5);
  for (r2=0.0,j=1;j<=nn;j++) r2 += ww[j];
  r2=1.0/r2;
  bmx=BIG;
  bmn=BIG;
  offs=(*chi2)+1.0;
  for (j=1;j<=6;j++) {
    if (ch[j] > offs) {
      d1=fabs(ang[j]-(*b));
      while (d1 >= PI) d1 -= PI;
      d2=PI-d1;
      if (ang[j] < *b) {
	swap=d1;
	d1=d2;
	d2=swap;
      }
      if (d1 < bmx) bmx=d1;
      if (d2 < bmn) bmn=d2;
    }
  }
  if (bmx < BIG) {
    bmx=zbrent((COMP_PRECISION (*)())chixy,*b,*b+bmx,ACC)-(*b);
    amx=aa-(*a);
    bmn=zbrent((COMP_PRECISION (*)())chixy,*b,*b-bmn,ACC)-(*b);
    amn=aa-(*a);
    tmpd = cos(*b);
    *sigb=sqrt(0.5*(bmx*bmx+bmn*bmn))/(scale*SQUARE(tmpd));
    *siga=sqrt(0.5*(amx*amx+amn*amn)+r2)/scale;
  } else (*sigb)=(*siga)=BIG;
  *a /= scale;
  *b=tan(*b)/scale;
  free_vector(ww,1,ndat);
  free_vector(sy,1,ndat);
  free_vector(sx,1,ndat);
  free_vector(yy,1,ndat);
  free_vector(xx,1,ndat);
}
#undef POTN
#undef BIG
#undef PI
#undef ACC
#define BIG 1.0e30
COMP_PRECISION chixy(COMP_PRECISION bang)
{
  int j;
  COMP_PRECISION ans,avex=0.0,avey=0.0,sumw=0.0,b,tmpd;
  
  b=tan(bang);
  for (j=1;j<=nn;j++) {
    tmpd = b*sx[j];
    ww[j] = SQUARE(tmpd)+SQUARE(sy[j]);
    sumw += (ww[j] = (ww[j] == 0.0 ? BIG : 1.0/ww[j]));
    avex += ww[j]*xx[j];
    avey += ww[j]*yy[j];
  }
  if (sumw == 0.0) sumw = BIG;
  avex /= sumw;
  avey /= sumw;
  aa=avey-b*avex;
  for (ans = -offs,j=1;j<=nn;j++){
    tmpd = yy[j]-aa-b*xx[j];
    ans += ww[j]*SQUARE(tmpd);
  }
  return ans;
}
#undef BIG


#define ITMAX 1000
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

COMP_PRECISION brent(COMP_PRECISION ax,COMP_PRECISION bx,
		     COMP_PRECISION cx,COMP_PRECISION (*f)(),
		     COMP_PRECISION tol,COMP_PRECISION *xmin)
{
  int iter;
  COMP_PRECISION a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  COMP_PRECISION e=0.0;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
	  if (u < x) a=u; else b=u;
	  if (fu <= fw || w == x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	  } else if (fu <= fv || v == x || v == w) {
	    v=u;
	    fv=fu;
	  }
	}
  }
  nrerror("Too many iterations in brent");
  *xmin=x;
  return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT



COMP_PRECISION gammq(COMP_PRECISION a,COMP_PRECISION x)
{
  COMP_PRECISION gamser,gammcf,gln;
  
  if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return gammcf;
  }
}
#define ITMAX 1000
#define EPS EPS_COMP_PREC
#define FPMIN 1.0e-30

void gcf(COMP_PRECISION *gammcf,COMP_PRECISION a,
	 COMP_PRECISION x,COMP_PRECISION *gln)
{
  int i;
  COMP_PRECISION an,b,c,d,del,h;
  
  *gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN


#define ITMAX 1000
#define EPS EPS_COMP_PREC

void gser(COMP_PRECISION *gamser,COMP_PRECISION a,
	  COMP_PRECISION x,COMP_PRECISION *gln)
{
  int n;
  COMP_PRECISION sum,del,ap;
  
  *gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) nrerror("x less than 0 in routine gser");
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}
#undef ITMAX
#undef EPS


void avevar(COMP_PRECISION *data,unsigned long n,
	    COMP_PRECISION *ave,COMP_PRECISION *var)
{
  unsigned long j;
  COMP_PRECISION s,ep;
  
  for (*ave=0.0,j=1;j<=n;j++) *ave += data[j];
  *ave /= n;
  *var=ep=0.0;
  for (j=1;j<=n;j++) {
    s=data[j]-(*ave);
    ep += s;
    *var += s*s;
  }
  *var=(*var-ep*ep/n)/(n-1);
}


void fit(COMP_PRECISION *x,COMP_PRECISION *y,
	 int ndata,COMP_PRECISION *sig,int mwt,
	 COMP_PRECISION *a,COMP_PRECISION *b,
	 COMP_PRECISION *siga,COMP_PRECISION *sigb,
	 COMP_PRECISION *chi2,COMP_PRECISION *q)
{
  int i;
  COMP_PRECISION wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat,tmpd;
	
	*b=0.0;
	if (mwt) {
		ss=0.0;
		for (i=1;i<=ndata;i++) {
			wt=1.0/SQUARE(sig[i]);
			ss += wt;
			sx += x[i]*wt;
			sy += y[i]*wt;
		}
	} else {
		for (i=1;i<=ndata;i++) {
			sx += x[i];
			sy += y[i];
		}
		ss=ndata;
	}
	sxoss=sx/ss;
	if (mwt) {
		for (i=1;i<=ndata;i++) {
			t=(x[i]-sxoss)/sig[i];
			st2 += t*t;
			*b += t*y[i]/sig[i];
		}
	} else {
		for (i=1;i<=ndata;i++) {
			t=x[i]-sxoss;
			st2 += t*t;
			*b += t*y[i];
		}
	}
	*b /= st2;
	*a=(sy-sx*(*b))/ss;
	*siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
	*sigb=sqrt(1.0/st2);
	*chi2=0.0;
	if (mwt == 0) {
	  for (i=1;i<=ndata;i++){
	    tmpd = y[i]-(*a)-(*b)*x[i];
	      *chi2 += SQUARE(tmpd);
	  }
		*q=1.0;
		sigdat=sqrt((*chi2)/(ndata-2));
		*siga *= sigdat;
		*sigb *= sigdat;
	} else {
	  for (i=1;i<=ndata;i++){
	    tmpd= (y[i]-(*a)-(*b)*x[i])/sig[i];
	    *chi2 += SQUARE(tmpd);
	  }
		*q=gammq(0.5*(ndata-2),0.5*(*chi2));
	}
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(COMP_PRECISION *ax,COMP_PRECISION *bx,
	    COMP_PRECISION *cx,COMP_PRECISION *fa,
	    COMP_PRECISION *fb,COMP_PRECISION *fc,
	    COMP_PRECISION (*func)())
{
	COMP_PRECISION ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#define ITMAX 1000
#define EPS EPS_COMP_PREC

COMP_PRECISION zbrent(COMP_PRECISION (*func)(),COMP_PRECISION x1,
		      COMP_PRECISION x2,COMP_PRECISION tol)
{
	int iter;
	COMP_PRECISION a,b,c,d,e=0.,min1,min2;
	COMP_PRECISION fa,fb,fc,p,q,r,s,tol1,xm;
	d =0.;
	a = x1;
	b = c = x2;
	fa=(*func)(a);
	fb=(*func)(b);
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
#undef ITMAX
#undef EPS

