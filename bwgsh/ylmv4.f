c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c       $Id: ylmv4.f,v 1.3 2006/08/03 23:34:53 becker Exp $
c---------------------------------------------------------------------------
c
c
c
c        subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
c
c  evaluates the vectors ypp,ytp, of length     lenv= (lmax-1)*(2*lmax+6)
c  and vectors           ypppp,ytppp, of length lenv4=(lmax-3)*(2*lmax+10).
c
c  ypp,ytp represent contributions to the symmetric trace-free
c  second rank tensor c.
c
c  ypppp,ytppp represent contributions to the completely symmetric trace-free
c  fourth rank tensor e.
c
c    In spherical coordinates ('t'=theta=colatitiude, 'p'=phi=longitude)
c
c    c  =  ( c_tt  c_tp )
c          ( c_pt  c_pp )
c
c    and  c_pp   = -c_tt   = ypp . coefs
c         c_tp    = c_pt   = ytp . coefs
c
c     Similarly for e:
c     e_pppp =  e_tttt = -e_pptt   = ypppp . coefs
c     e_tppp = -e_tttp             = ytppp . coefs
c       
c 
c     Scalar spherical harmonics are also calculated and
c     placed in y(1) -- y( (lmax+1)**2 )
c
c------------------------------------------------------------------
c I believe scalar harmonics here are orthogonal but NOT
c orthonormal: a factor sqrt(2) is missing from harmonics with
c nonzero m. see Dahlen and Tromp eq. B.72
c------------------------------------------------------------------
c
c     The companion routine ylmavv4() (q.v.) calculates the
c     contribution to the average value of  c_ll (= l.c.l) and e_llll
c     along the minor arc and the complete great circle,
c     where 'l' represents the tangent to the path:
c         
c           l_t = -cos(az)      l_p = sin(az)
c 
c     and az is the local azimuth of the path.
c
c     Thus:
c            c_ll   =  - c_pp * cos(2*az)  - c_tp * sin(2*az)
c     and  
c            e_llll =  e_pppp * cos(4*az) + e_tppp * sin(4*az)
c----------------------------------------------------------------
c     'd' is d.p. workspace. (notice that it's bigger than in ylmv)
c----------------------------------------------------------------
#include "prec.h"

      subroutine ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,d)
      dimension y((lmax+1)**2)
     1         ,ypp((lmax-1)*(2*lmax+6))
     1         ,ytp((lmax-1)*(2*lmax+6))
     1         ,ypppp((lmax-3)*(2*lmax+10))
     1         ,ytppp((lmax-3)*(2*lmax+10))
      double precision theta,d(9*(2*lmax+1))
      complex cfac,dfac
      data radian/57.295779513082321d0/,rr4pi/0.282094791773878d0/
      theta=(90.d0-xlat)/radian
      dfac=cexp(cmplx(0.d0,xlon/radian))
      k=0
      ka=0
      ka4=0
      do l=0,lmax               ! l loop
         i2lp1 = 2*l+1
         if(l.lt.2) then
c     l < 2
            call rotmx2(0,l,theta,d,1,i2lp1)
            ind=l
            cfac=rr4pi*sqrt(dfloat(i2lp1))
            do m=0,l            ! m loop
               k=k+1
               ind=ind+1
               y(k) = d(ind) * real(cfac) ! scalar a
               if(m.ne.0) then ! m != 0 
                  k=k+1
                  y(k)=d(ind)* aimag(cfac) ! scalar b
               endif
               cfac=cfac*dfac
            enddo
         else if(l.lt.4) then
            call rotmx2(2,l,theta,d,5,i2lp1)
            ind=5*l+3
            indp=ind+2
            indm=indp
            cfac=rr4pi * sqrt(dfloat(i2lp1))
            do m=0,l            ! m loop
               k=k+1
               y(k)=d(ind)*real(cfac) ! scalar a 
               ! 2phi
               ka=ka+1
               ypp(ka) = -d(indp)*real(cfac)  ! real a
               ytp(ka) = -d(indp)*aimag(cfac) ! imag a
               ka1=ka+1
               ypp(ka1)= -ytp(ka)
               ytp(ka1)=  ypp(ka)
               ka=ka+1
               if(m.ne.0) then ! m != 0 
                  k=k+1
                  y(k)=d(ind)*aimag(cfac)
                  ka=ka+1
                  ypp(ka)= -d(indm)*real(cfac)   ! real b 
                  ytp(ka)= +d(indm)*aimag(cfac)  ! imag b
                  ka1=ka+1
                  ypp(ka1)=-ytp(ka)
                  ytp(ka1)= ypp(ka)
                  ka=ka+1
               endif
               ind=ind+5
               indp=indp+5
               indm=indm-5
               cfac=cfac*dfac
          enddo
       else
          call rotmx2(4,l,theta,d,9,i2lp1)
          ind=9*l+5
          indp=ind+2
          indm=indp
          indp4=ind+4
          indm4=indp4
          cfac=rr4pi*sqrt(dfloat(i2lp1))
          do m=0,l
             k=k+1
             y(k)=d(ind)*real(cfac)
             ! 2phi
             ka=ka+1
             ypp(ka)= -d(indp)*real(cfac)
             ytp(ka)= -d(indp)*aimag(cfac)
             ka1=ka+1
             ypp(ka1)= -ytp(ka)
             ytp(ka1)=  ypp(ka)
             ka=ka+1
             ! 4phi
             ka4=ka4+1
             ypppp(ka4)= -d(indp4)*real(cfac)
             ytppp(ka4)= -d(indp4)*aimag(cfac)
             ka41=ka4+1
             ypppp(ka41)= -ytppp(ka4)
             ytppp(ka41)=  ypppp(ka4)
             ka4=ka4+1
             if(m.ne.0) then
                ! scalar
                k=k+1
                y(k)=d(ind)*aimag(cfac)
                ! 2phi
                ka=ka+1
                ypp(ka)= -d(indm)*real(cfac)
                ytp(ka)= +d(indm)*aimag(cfac)
                ka1=ka+1
                ypp(ka1)= -ytp(ka)
                ytp(ka1)=  ypp(ka)
                ka=ka+1
                ! 4phi
                ka4=ka4+1
                ypppp(ka4)= -d(indm4)*real(cfac)
                ytppp(ka4)= +d(indm4)*aimag(cfac)
                ka41=ka4+1
                ypppp(ka41)= -ytppp(ka4)
                ytppp(ka41)=  ypppp(ka4)
                ka4=ka4+1
             endif
             ind=ind+9
             indp=indp+9
             indm=indm-9
             indp4=indp4+9
             indm4=indm4-9
             cfac=cfac*dfac
          enddo
       endif
      enddo
      return
      end
