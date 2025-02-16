c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c     $Id: xyz2gsh.f,v 1.5 2006/08/25 20:56:38 becker Exp becker $
c
c---find expansion of xyz function over generalized (N=2 or N=4) or scalar spherical
c---harmonics (N=0) by solution of least squares problem. needed to express
c---maps of azimuthal anisotropy in terms of generalized s.h. coefficients.
c 
c       CHANGED INPUT FORMAT
c
c       reads in lon lat dat    for scalar
c
c       and
c
c       lon lat e1 e2              for 2phi and 4phi 
c
c
c
#include "prec.h"
	parameter(lmx=50)
c	parameter(lmx=63)
	
	parameter(nunk=(lmx-1)*(2*lmx+6))
	parameter(ilata=(nunk*(nunk+1))/2)
c nx and ny are actuallly not used, idsize is the important part 
c	parameter(nx=182,ny=91,idsize=nx*ny) !increase idsize if spacing < 2 deg
	parameter(idsize=30000,idsize2=idsize*2)
c	parameter(nx=364,ny=182,idsize=nx*ny) !increase idsize if spacing < 1 deg (THIS LINE WILL MAKE HUGE BINARIES)

	character namein*300,nameout*300,namein2*300
	dimension y((lmx+1)**2)
     &         ,ypp((lmx-1)*(2*lmx+6))
     &         ,ytp((lmx-1)*(2*lmx+6))
     &         ,ypppp((lmx-3)*(2*lmx+10))
     &         ,ytppp((lmx-3)*(2*lmx+10))
	double precision wkspc(9*(2*lmx+1))
	dimension dat(idsize2),rowa(nunk),ata(ilata),atd(nunk),x(nunk)
     &         ,gstore(ilata),ytemp(nunk),a(idsize2,nunk)
	CPREC ataout(ilata),atdout(nunk),xmasked(nunk)
	integer maskv(nunk),maskata(ilata)
	integer*4 ldep_damp,use_mc
	double precision dampcoef(nunk),damp_scale,damp_fac

	write(0,*),'lmax?'
	read *,lmax
	if(lmax.gt.lmx)then
	   write(0,*),'lmax out of bounds',lmax,lmx
	   stop
	endif

	write(0,*),"term 0, 2, or 4?"
	read*,ialpha
	if(ialpha.eq.0)then
c       isotropic
	   write(0,*),"input model?"
	   read*,namein
	   open(1,file=namein,status='old')
	   ncoef=(lmax+1)**2
	elseif(ialpha.eq.2.or.ialpha.eq.4)then
c       2phi/4phi
	   if(ialpha.eq.2)then ! 2phi
	      ncoef=(lmax-1)*(2*lmax+6)
	      if(lmax.lt.2)then
		 stop 'lmax needs to > 1 for 2phi'
	      endif
	   else ! 4phi
	      ncoef=(lmax-3)*(2*lmax+10)
	      if(lmax.lt.4)then
		 stop 'lmax needs to > 3 for 4phi'
	      endif
	   endif
	      
	   write(0,*),"a1 a2 or a3 a4 input model?"
	   read*,namein
	   open(1,file=namein,status='old')
	else
	   stop "non-existing term"
	endif

	write(0,*),"output file?"
	read*,nameout

	write(0,*),"damping factor?"
	read*,damp
	
	write(0,*),"degree dependent damping (0/1)?"
	read*,ldep_damp

	do i=1,ilata
	   ata(i)=0.d0
	enddo

c----------------read in all the gridpoints and increment matrices
	i=1
	
1	continue

	if(i.gt.idsize)then
	   write(0,*),'too many data points'
	   write(0,*),i,idsize
	   stop
	endif
	i1=i+1
	
	if(ialpha.eq.0)then
c-------------------------------scalar SH
	   read(1,*,end=20)xlon,xlat,dat(i)
	   call ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,wkspc)
	   do k=1,ncoef
	      rowa(k)=y(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)
	elseif(ialpha.eq.2)then
c-------------------------------generalized SH for second order tensor
	   read(1,*,end=20)xlon,xlat,dat(i),dat(i1)
	   dat(i) = -dat(i)
	   dat(i1)= -dat(i1)

	   call ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,wkspc)

	   do k=1,ncoef
	      rowa(k)=ypp(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)

	   do k=1,ncoef
	      rowa(k)=ytp(k)
	      a(i1,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i1),ncoef,nunk,ilata,ata,atd)
	   i=i+1
	   
	elseif(ialpha.eq.4)then
c-------------------------------generalized SH for fourth order tensor

	   read(1,*,end=20)xlon,xlat,dat(i),dat(i1)
	   
	   call ylmv4(xlat,xlon,lmax,y,ypp,ytp,ypppp,ytppp,wkspc)
	   do k=1,ncoef
	      rowa(k)=ypppp(k)
	      a(i,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i),ncoef,nunk,ilata,ata,atd)

	   do k=1,ncoef
	      rowa(k)=ytppp(k)
	      a(i1,k)=rowa(k)
	   enddo
	   call contribution_ata(rowa,dat(i1),ncoef,nunk,ilata,ata,atd)
	   i=i+1
	else
	   stop "error in variable ialpha"
	endif
	
c	if(mod(i,1000).eq.0)write(0,*),i," gridpoints read"
c--go back to read the next datum
	i=i+1

c	print *,i,idsize
	goto 1

 20	continue 

	ndat=i-1
	
	if(ialpha.eq.0)then
	   write(0,*),ndat,' gridpoints read, scalar, ndat ',ndat
	else
	   write(0,*),ndat/2,' gridpoints read, 2phi or 4phi, ndat ',ndat
	endif

	close(1)
	
c---------------------------------------regularization
	k=0
	il0=ialpha
c----------------------MUST BE CHANGED FOR GENERALIZED SH PARAMETERIZATION
	damp_scale = float(lmax)/2.0d0*(float(lmax)/2.d0+1.0d0)
	do il=il0,lmax
	   if(ialpha.eq.0)then
C       scalara
	      immax=2*il+1
	   else
C       2psi/4psi
	      immax=2*(2*il+1)
	   endif
	   damp_fac = float(il*(il+1))/damp_scale + 0.1d0
	   ! print *,il,damp_fac
	   do im=1,immax
	      k=k+1
	      if(k.gt.nunk)then
		 print *,'error, k ',k,' nunk ',nunk
		 stop
	      endif
	      dampcoef(k)= damp_fac 
	   enddo
	enddo

! assign damping to the diagonal elements of the ATA matrix?
	k=0
	icount=0
	do ix=1,ncoef
	   do iy=1,ix
	      k=k+1
	      if(iy.eq.ix)then
	         icount=icount+1
		 if(ldep_damp.ne.0)then
c       wave length dependent damping
		    ata(k)=ata(k)+damp*dampcoef(icount)
		 else
c------------NON-L-DEPENDENT DAMPING 
		    ata(k)=ata(k)+damp
		 endif
	      endif
	   enddo
	enddo

cc--mask coefficients associated with basis functions that cannot be constrained
c SHOULD NOT BE NEEDED NOW
c	write(0,*),"mask ata matrix and atd vector..."
c	call atamask(ata,atd,ataout,atdout,ilata,nunk,maskata,maskv,nonzero)

c--find least squares solution x via single-processor cholesky factorization of ATA
c	write(0,*),"cholesky factorization..."

	use_mc=0
	if(use_mc.eq.0)then
c	print *,'using own choles'
c	solve (A^TA+d I) x = A^T d
c       A-packed storage, storage, right hand side, temp vector, solution vector dimenions, error code
	   print *,'using home made Cholesky solve'
	   call choles(ata,gstore,atd,ytemp,x,ncoef,nono)
c	old,not sure
c	call choles(ataout,gstore,atdout,ytemp,xmasked,nonzero,nono)

	
	   if(nono.ne.0)then
	      write(0,*),"cholesky factorization failed ",nono
c	do j=0,nonzero-1
c	istart=j*(j+1)/2+1
c	j1=j+1
c	istop=j1*(j1+1)/2
c       c	write(*,*)j+1,(ata(i),i=istart,istop)
c       c	pause
c	enddo
	      stop
	   endif
	else
	   print *,'using LAPACK Cholesky'
	   call LAPACK_CFACTOR('U',ncoef,ata,nono) ! cholesky factorization
	   if(nono.ne.0)then
	      print *,'LAPACK factorization failed',nono
	      stop
	   endif
	   call LAPACK_CSOLVER('U',ncoef,1,ata,atd,ncoef,nono)
	   if(nono.ne.0)then
	      print *,'LAPACK solve failed',nono
	      stop
	   endif
	   do icl=1,ncoef
	      x(icl)=atd(icl)
	   enddo
	endif
	


c--unmask solution coefficients that were eliminated by atamask
c	write(0,*),"unmask solution..."
c	k=0
c	do i=1,nunk
c	   if(maskv(i).eq.1)then
c	      k=k+1
c	      x(i)=xmasked(k)
c	   elseif(maskv(i).eq.0)then
c	      x(i)=0.
c	   else
c	      write(0,*),"illegal maskv value for coeff. ",i
c	      stop
c	   endif
c	enddo

c--see how well the least squares solution fits the gridpoints
c	write(0,*),"computing variance reduction..."
	truerms=0.d0
	realchisq=0.d0
	denom=0.d0
	do irow=1,ndat
	   error=0.d0
	   syndat=0.d0
	   do icl=1,ncoef
	      syndat=syndat+a(irow,icl)*x(icl)
	   enddo
	   error=(syndat-dat(irow))**2
	   realchisq=error+realchisq
	   truerms=truerms+abs(syndat-dat(irow))
	   denom=denom+(dat(irow)*dat(irow))
	enddo
	truerms=truerms/float(ndat)
	varred=realchisq/denom
	varred=1.-varred
	write(0,*),'NUMBER OF DATA:',ndat
	write(0,*),'CHI-SQUARE:',realchisq
	write(0,*),'RMS MISFIT OBTAINED:',truerms
	write(0,*),'VARIANCE REDUCTION:',varred
c
c       compute approximate norm
c
	xnorm = 0.d0
	do i=1,ncoef
	   xnorm = xnorm + x(i)**2
	enddo
	xnorm = sqrt(xnorm)
	write(6,'(i1,1x,i3,1x,e12.4,1x,i6,1x,2(e20.7,1x),f10.7,1x,e20.7,1x,i1)')
     &        ialpha,lmax,damp,ndat,realchisq,truerms,varred,
     &        xnorm,ldep_damp

c
c       output
c
	open(2,file=nameout)

	write(2,*)lmax,ialpha
	write(2,'(5(e20.7))')(x(i),i=1,ncoef)
	close(2)
	
	end

