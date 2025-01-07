c
c     code by boschi & woodhouse (GJI, 2006). distribution from SPICE website/Lapo Boschi as of 08/2006
c
c     minor changes by twb@usc.edu
c
c     $Id: aniso2xyz.f,v 1.4 2006/08/03 01:10:42 becker Exp $
c
c--find eps1 and eps2 or eps3 and eps4 from fast azim. and ampl. of anisotropy
      character fileamp*80,fileang*80,stri1*4,stri2*4
      
      
      print*,"2psi or 4psi (2/4)"
      read*,ipsi
      if(ipsi.ne.2.and.ipsi.ne.4)stop "non existing term"
      print*,"xyz file with max amplitude"
      read*,fileamp
      print*,"xyz file with fast azimuth"
      read*,fileang
      open(1,file=fileamp,status="old")
      open(2,file=fileang,status="old")
      do k=1,80
         if(fileamp(k:k).eq." ")goto10
      enddo
 10   if(ipsi.eq.2)then 
         stri1="EPS1"
         stri2="EPS2"
      else
         stri1="EPS3"
         stri2="EPS4"
      endif
      open(91,file=fileamp(1:k-1)//"."//stri1)
      open(92,file=fileamp(1:k-1)//"."//stri2)
 1    read(1,*,end=2)x,y,amp
      read(2,*,end=3)x1,y1,ang
      if(x.ne.x1.or.y.ne.y1)stop "files are not compatible"
      if(ipsi.eq.2)then
         xtmp = 2.d0*ang
         a2_a1=tan(xtmp)
         a1=amp/(cos(xtmp)+(a2_a1*sin(xtmp)))
      else
         xtmp = 4.d0*ang
         a2_a1=tan(xtmp)
         a1=amp/(cos(xtmp)+(a2_a1*sin(xtmp)))
      endif
      a2=a1*a2_a1
      write(91,*)x,y,a1
      write(92,*)x,y,a2
      goto1
 3    print*,"end of file 2 reached"
 2    print*,"end of file 1 reached (OK)"
      close(1)
      close(2)
      close(91)
      close(92)
      end
