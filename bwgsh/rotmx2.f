c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c     $Id: rotmx2.f,v 1.3 2006/08/03 00:20:24 becker Exp $

c---------------------------------------------------------------------------
cprog rotmx2
cxref
#include "prec.h"

      subroutine rotmx2(nmax,l,theta,d,id1,id2)
      implicit double precision (a-h,o-z)
      double precision d,theta,xtmp
      dimension d(id1,id2)
c     data big,small,dlbig,dlsml/1.d35,1.d-35,35.d0,-35.d0/
      data big,small,dlbig,dlsml/1.d25,1.d-25,25.d0,-25.d0/
      data pi/3.14159265358979d0/
      dfloat(n)=n
      th=theta
      if((th.gt.pi).or.(th.lt.0.d0)) stop 'illegal arg in rotmx2'
      if(l.ne.0) goto 350
      d(1+nmax,l+1)=1.d0
      return
350   isup=1
      if(th.le.pi/2.d0) goto 310
      th=pi-th
      isup=-1
310   nm=2*l+1
      nmp1=nm+1
      lp1=l+1
      lm1=l-1
      lp2=l+2
      nrow=2*nmax+1
      nmaxp1=nmax+1
      lmn=l-nmax
      if(th.ne.0.d0) goto 320
      do 330 im1ct=1,nrow
      im1=im1ct+lmn
      do 330 im2=lp1,nm
      d(im1ct,im2)=0.d0
      if(im1.eq.im2) d(im1ct,im2)=1.d0
330   continue
      goto 400
320   continue
c
c     zero l.h.s. of matrix
c
      do 340 im1=1,nrow
      do 340 im2=1,lp1
340   d(im1,im2)=0.d0
c
c        set up parameters
c
      xtmp = 0.5d0*th
      shth=dsin(xtmp)
      chth=dcos(xtmp)
      xtmp = 2.d0*chth
      sth=xtmp*shth
      cth=xtmp*chth-1.d0

      dlogf=dlog10(chth/shth)
      dlogs=dlog10(shth)
c
c       iterate from last column using 1. as starting value
c
      do 10 im1ct=1,nrow
      im1=im1ct+lmn
      m1=im1-lp1
      rm1=m1
      nm2=min0(im1-1,nm-im1)
      d(im1ct,nm)=1.d0
      if(nm2.eq.0) goto 10
      do 20 nit=1,nm2
      m2=l-nit
      im2=m2+lp1
      if(m2.ne.lm1) goto 70
      t1=0.d0
      goto 30
70    t1=-dsqrt(dfloat((im2+1)*(l-m2-1)))*d(im1ct,im2+2)
30    d(im1ct,im2)=t1-(2.d0/sth)*(cth*dfloat(m2+1)-rm1)
     1    *d(im1ct,im2+1)
      d(im1ct,im2)=d(im1ct,im2)/dsqrt(dfloat(im2*(l-m2)))
      temp=d(im1ct,im2)
      rmod=dabs(temp)
      if(rmod.lt.big) goto 20
      if(nit.eq.nm2) goto 20
      d(im1ct,nit+1)=dlbig
      d(im1ct,im2)=d(im1ct,im2)/big
      d(im1ct,im2+1)=d(im1ct,im2+1)/big
20    continue
10    continue
c
c        set up normalization for rightmost column
c
      t1=dfloat(2*l)*dlogs
      if(lmn.eq.0) goto 720
      do 710 i=1,lmn
      m1=i-l
      t1=dlogf+0.5d0*dlog10(dfloat(lp1-m1)/dfloat(l+m1))+t1
710   continue
720   d(1,1)=t1
      if(nrow.eq.1) goto 730
      do 110 im1ct=2,nrow
      m1=im1ct-nmaxp1
110   d(im1ct,1)=dlogf+0.5d0*dlog10(dfloat(l-m1+1)/dfloat(l+m1))
     1     +d(im1ct-1,1)
730   sgn=-1.d0
      if((lmn/2)*2.ne.lmn) sgn=1.d0
c
c       renormalize rows
c
      do 120 im1ct=1,nrow
      im1=im1ct+lmn
      sgn=-sgn
      csum=d(im1ct,1)
      mult=1
520   if(dabs(csum).lt.dlbig) goto 510
      mult=mult*2
      csum=0.5*csum
      goto 520
510   fac=10.d0**csum
      sfac=small/fac
      nm2=min0(im1-1,nm-im1)
      nm2p1=nm2+1
      do 130 im2=1,nm2p1
      if((d(im1ct,im2+1).eq.0.d0).or.(im2.ge.nm2)) goto 250
      csum=csum*dfloat(mult)+d(im1ct,im2+1)
      mult=1
220   if(dabs(csum).lt.dlbig) goto 210
      mult=mult*2
      csum=0.5d0*csum
      goto 220
210   fac=10.d0**csum
      sfac=small/fac
250   in2=nmp1-im2
      do 270 i=1,mult
      temp=d(im1ct,in2)
      rmod=dabs(temp)
      if(rmod.gt.sfac) goto 260
      d(im1ct,in2)=0.d0
      goto 130
260   d(im1ct,in2)=d(im1ct,in2)*fac
270   continue
      d(im1ct,in2)=sgn*d(im1ct,in2)
130   continue
120   continue
c
c       fill rest of matrix
c
400   if(isup.gt.0) goto 410
      sgn=-1.d0
      if((lmn/2)*2.ne.lmn) sgn=1.d0
      do 420 im1ct=1,nrow
      sgn=-sgn
      im1=im1ct+lmn
      nm2=min0(im1,nmp1-im1)
      do 420 in2=1,nm2
      im2=nmp1-in2
420   d(im1ct,in2)=sgn*d(im1ct,im2)
      do 430 im1ct=1,nrow
      im1=im1ct+lmn
      in1=nmp1-im1
      in1ct=in1-lmn
      sgn=-1.d0
      nm2=min0(im1,in1)
      do 440 nit=1,nm2
      sgn=-sgn
      im2=1+nm2-nit
      in2=nmp1-im2
      im2ct=im2-lmn
      in2ct=in2-lmn
      d(in1ct,in2)=sgn*d(im1ct,im2)
      if(in2ct.gt.nrow) goto 440
      d(im2ct,im1)=d(in1ct,in2)
      d(in2ct,in1)=d(im1ct,im2)
440   continue
430   continue
      return
410   do 450 im1ct=1,nrow
      im1=im1ct+lmn
      in1=nmp1-im1
      in1ct=in1-lmn
      sgn=-1.d0
      nm2=min0(im1,in1)
      do 460 nit=1,nm2
      sgn=-sgn
      im2=nm-nm2+nit
      in2=nmp1-im2
      im2ct=im2-lmn
      in2ct=in2-lmn
      d(in1ct,in2)=sgn*d(im1ct,im2)
      if(im2ct.gt.nrow) goto 460
      d(im2ct,im1)=d(in1ct,in2)
      d(in2ct,in1)=d(im1ct,im2)
460   continue
450   continue
      return
      end
