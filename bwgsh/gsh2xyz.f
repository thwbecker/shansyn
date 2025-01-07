c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c       $Id: gsh2xyz.f,v 1.4 2006/08/25 20:56:47 becker Exp $
c
c expand a spherical harmonic expansion to xyz
c

	parameter(lmx=55)
	character namein*80,nameout*80,chlmax*3,chialpha*1
	dimension y((lmx+1)**2)
     &         ,ypp((lmx-1)*(2*lmx+6))
     &         ,ytp((lmx-1)*(2*lmx+6))
     &         ,ypppp((lmx-3)*(2*lmx+10))
     &         ,ytppp((lmx-3)*(2*lmx+10))
	double precision wkspc(9*(2*lmx+1))
	dimension cmod((LMX+1)**2)
c
c	input parameters
c
	write(0,*)"what input model?"
	read*,namein
	write(0,*)"what increment for grid?"
	read*,xincr
c--------------------following line useful e.g. if model given as
c--------------------relative perturbation and we need percent
	write(0,*)"scale model by? (1 if no scale)"
	read*,scale

	open(1,file=namein,status='old')
	read(1,*)lmax,ialpha
	if((ialpha.ne.0).and.(ialpha.ne.2).and.(ialpha.ne.4))then
	   print *,'ialpha read: ',ialpha
	   stop "non-existing ialpha term"
	endif

	write(0,*)"maximum l=",lmax
	write(0,*)"interrupt expansion at l=? (-1 for default)"
	read*,lmaxuser
	if(lmaxuser.eq.-1)then
	   lmaxuser = lmax
	endif
!
!       compute number of coefficients
!
	if(ialpha.eq.0)then
	   ncoef=(lmax+1)**2
	elseif(ialpha.le.2)then
	   ncoef=(lmax-1)*(2*lmax+6)
	else
	   ncoef=(lmax-3)*(2*lmax+10)
	endif
c
c       read in coefficients
c
	read(1,*)(cmod(i),i=1,ncoef)
	close(1)
	do k=1,80
	   if(namein(k:k).eq.' ')goto2
	enddo
 2	continue
c
c       limit lmax?
c
	lmax=lmaxuser		! user can filter-lapo 2004-04-03
	write(chlmax,'(i3.3)')lmax
	write(chialpha,"(i1.1)")ialpha
c	nameout=namein(1:k-1)//".eps"//chialpha//"."//chlmax//'.xyz'
	nameout=namein(1:k-1)//'.xyz'
	open(2,file=nameout)

	igrid=0
	if(lmax.gt.lmx)stop "l too big"
	xtmp = xincr/2.d0
	do xlat=-90.d0+xtmp,90.d0-xtmp,xincr ! lon loop
	   do xlon=0.d0,360.d0-xincr,xincr ! lat loop
	      igrid=igrid+1
c	      if(mod(igrid,1000).eq.0)write(0,*)igrid," grid points done"
c       compute spherical harmonics at lon,lat 
	      call ylmv4(xlat,xlon,lmx,y,ypp,ytp,ypppp,ytppp,wkspc)
c       !
c       sum up contributions
c       
	      z=0.d0
	      z2=0.0d0
	      if(ialpha.eq.0)then
c-------------------------------scalar SH
		 do k=1,ncoef
		    z=z+y(k)*cmod(k)
		 enddo
		 write(2,*)xlon,xlat,z*scale
c-------------------------------generalized SH for second order tensor
	      elseif(ialpha.eq.2)then
		 do k=1,ncoef
		    z=  z-ypp(k)*cmod(k)
		    z2=z2-ytp(k)*cmod(k)
		 enddo
		 write(2,*)xlon,xlat,z*scale,z2*scale
	      elseif(ialpha.eq.4)then
c-------------------------------generalized SH for fourth order tensor
		 do k=1,ncoef
		    z=z  +ypppp(k)*cmod(k)
		    z2=z2+ytppp(k)*cmod(k)
		 enddo
		 write(2,*)xlon,xlat,z*scale,z2*scale
	      else
		 stop "error in variable ialpha"
	      endif

	   enddo
	enddo

	close(2)
	end
