c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c       $Id: xyz2aniso.f,v 1.4 2006/08/03 01:10:42 becker Exp $
c
c-----find fast 2psi or 4psi direction, and azimuthal anisotropy variation 
c-----amplitude, over a grid, from gridded alpha_i files
#include "prec.h"
	character*80 namein1,namein2
	character chipsi*1
	pi=2.d0*asin(1.d0)
	print*,"2psi or 4psi term? (2/4)"
	read*,ipsi
	if(ipsi.ne.2.and.ipsi.ne.4)stop "non-existing term"
	write(chipsi,"(i1.1)")ipsi
	print*,"input file? (a1 and a2 or a3 and a4)"
	read*,namein1
c	print*,"second input file? (a2 or a4)"
c	read*,namein2
	print*,"what output? (1=azi, am; 2=vector comps)"
	read*,ivec
	do k=1,80
	   if(namein1(k:k).eq." ")goto1
	enddo
 1	open(1,file=namein1,status="old")
c	open(2,file=namein2,status="old")
	if(ivec.eq.1)then
	   open(91,file=namein1(1:k-1)//"."//chipsi//"pazi")
	   open(92,file=namein1(1:k-1)//"."//chipsi//"pam")
	else
	   open(91,file=namein1(1:k-1)//"."//chipsi//".x")
	   open(92,file=namein1(1:k-1)//"."//chipsi//".y")
	   if(ipsi.eq.4)then
	      open(93,file=namein1(1:k-1)//"."//chipsi//".x2")
	      open(94,file=namein1(1:k-1)//"."//chipsi//".y2")
	   endif
	endif
2	read(1,*,end=3)xlon,xlat,v1,v2
c	read(2,*)xlon1,xlat1,v2
c	if(xlon.ne.xlon1.or.xlat.ne.xlat1)stop "grids are not compatible"
	am=-999999.
cTEST
	do i=1,ipsi
	   psi=(atan(v2/v1)+PCAST(i-1)*pi)/PCAST(ipsi)
	   c=v1*cos(PCAST(ipsi)*psi)+v2*sin(PCAST(ipsi)*psi)
	   if(c.gt.am)then
	      psi_fast=psi
	      am=c
	   endif
	enddo
	if(psi_fast.gt.pi/2.)psi_fast=psi_fast-pi
	if(psi_fast.lt.-pi/2.)psi_fast=psi_fast+pi
	
	if(ivec.eq.1)then
	   write(91,*)xlon,xlat,psi_fast
	   write(92,*)xlon,xlat,am
	else
	   write(91,*)xlon,xlat,am*sin(psi_fast)
	   write(92,*)xlon,xlat,am*cos(psi_fast)
	   if(ipsi.eq.4)then
	      write(93,*)xlon,xlat,-am*cos(psi_fast)
	      write(94,*)xlon,xlat,am*sin(psi_fast)
	   endif
	endif
	goto2
3	continue
	close(1)
	close(2)
	close(91)
	close(92)
	end
