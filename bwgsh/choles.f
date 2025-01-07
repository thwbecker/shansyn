c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c     $Id: choles.f,v 1.2 2006/08/02 21:09:49 becker Exp becker $
c
c---------------------------------------------------------------------------
c------------------------single-processor Cholesky factorization

#include "prec.h"      
      subroutine choles(a,g,b,y,x,n,nono)

      implicit CPREC (a-h, o-z)
      CPREC a(1),g(1)
      CPREC b(1),y(1),x(1)
c     
c        a= row-wise p.d. symm. system  n*(n+1)/2
c        g= cholesky storage
c        b= r.h.s. vector               n
c        y= temp. vector
c        x= answer vector
c        n= system dimension
c        nono .gt. 0 is the level at which p.d. failed
c
c        (a,g) and (b,y,x) may be equivalenced.
c
c----------------------------------------------------------
c-----first compute cholesky decomposition

      nono=0
      
      if(a(1).le.0.) then
         nono=1
         return
      endif
      
      g(1)=sqrt(a(1))
      y(1)=b(1)/g(1)
      
      do 400 i=2,n
         
         kz=(i*(i-1))/2
         g(kz+1)=a(kz+1)/g(1)
         sg=g(kz+1)**2
         y(i)=b(i)-g(kz+1)*y(1)
         
         if(i.gt.2) then
            
            jmax=i-1
            
            do 200 j=2,jmax
               
               gkz=a(kz+j)
               kj=(j*(j-1))/2
               kmax=j-1
               
               do 100 k=1,kmax
                  gkz=gkz-g(kz+k)*g(kj+k)
 100           continue

               g(kz+j)=gkz/g(kj+j)
               y(i)=y(i)-g(kz+j)*y(j)
               sg=sg+g(kz+j)**2
               
 200        continue

         endif
         gkz=a(kz+i)-sg
         
         if(gkz.le.0.) then
            nono=i
            return
         endif
         
         g(kz+i)=sqrt(gkz)
         y(i)=y(i)/g(kz+i)
         
 400  continue
      
      kz=(n*(n-1))/2
      x(n)=y(n)/g(kz+n)
      if(n.le.1) return
      
c-----
c     compute solution for particular rhs
      
      do 600 k=2,n
      
         i=n+1-k
         x(i)=y(i)
         jmin=i+1
         
         do 500 j=jmin,n
            kj=(j*(j-1))/2
            x(i)=x(i)-g(kj+i)*x(j)
 500     continue

         kz=(i*(i+1))/2
         x(i)=x(i)/g(kz)
         
 600  continue
      
      return
      end

c---------------------------------------------------------------------------
      subroutine atamask(a,v,aout,vout,ilata,nunk,maskata,maskv,kvout)
c---  removes from ata matrix (a) to be inverted columns/rows that are
c---  systematically 0. corresponding atd (v) entries are also removed.
      CPREC a(ilata),aout(ilata),v(nunk),vout(nunk)
      integer maskv(nunk),maskata(ilata)
      if(ilata.ne.nunk*(nunk+1)/2)stop "wrong dimensions"
      
      kvout=0
      do i=1,ilata
         maskata(i)=1
      enddo
      do i=1,nunk
         index=i*(i+1)/2
         if(a(index).eq.0.)then
            print*,i," diagonal entry is 0"
            maskv(i)=0
            do j=1,i
               index=i*(i-1)/2+j
               
c     TEST
c     print*,i,j,index,a(index)
               
               if(a(index).ne.0.)stop "this does not make sense"
               maskata(index)=0
            enddo
            do j=i+1,nunk
               index=j*(j-1)/2+i
               
c     TEST
c     print*,i,j,index,a(index)
               
               if(a(index).ne.0.)stop "this does not make sense"
               maskata(index)=0
            enddo
            
c     pause
            
         else
            maskv(i)=1
            kvout=kvout+1
            vout(kvout)=v(i)
         endif
      enddo
      kaout=0
      do i=1,ilata
         if(maskata(i).ne.0)then
            kaout=kaout+1
            aout(kaout)=a(i)
            
c     TEST
c     if(a(i).eq.0.)stop
            
         endif
      enddo
      return
      end
c---------------------------------------------------------------------------
      subroutine contribution_ata(rowa,dat,ncoef,nunk,ilata,ata,atd)
c---- given a row of A (rowa) and the corresponding datum (dat), add their
c---- contribution to the matrix A^tA and the vector A^td
c---
      integer nunk,ilata,ncoef
      
      CPREC rowa(nunk),ata(ilata),atd(nunk),dat

      k1=0
      do i1=1,ncoef
         atd(i1)=atd(i1)+rowa(i1)*dat
         do j1=1,i1
            k1=k1+1 
            ata(k1)=ata(k1)+rowa(i1)*rowa(j1)
         enddo
      enddo
      return
      end
c---------------------------------------------------------------------------
