       subroutine gmesh(nyg,nxg,yg,xg,ixg,iyg)
       implicit double precision(a-h,o-z)
       parameter(ig=1999)
       dimension xg(ixg),yg(iyg)
       dimension xmi(21),ymi(21),xma(21),yma(21)
c       dimension xy(501)	   !格子線の補正で用いる
       integer*2 jg(2,ig)
       common /dgrd/dx,dy,in0,in1,jn0,jn1
       common /region/nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &     neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
       common /sufcood/xb(21,1001),yb(21,1001)
       common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
       common /gmovedata/ inmin,jnmin
c
       write(*,*)"in gmesh nr=",nr     !
       open(1,file="col.inf")
       rewind(1)
       read(1,*)ny0,ny1,eta,rr,dwwy
       close(1)
       write(*,*)"ny0,ny1,eta,rr,dwwy"
       write(*,*)ny0,ny1,eta,rr,dwwy
c
       iyy=0
       write(*,*)"in gmesh iyy=",iyy
       xma(1)=-9999d0  
       xma(2)=-9999d0
       xmi(1)=9999d0
       xmi(2)=9999d0
       xmin=9999d0
       ymin=9999d0
       xmax=-9999d0
       ymax=-9999d0
       do 20 i=1,nr
         xmi(i)=9999d0
         ymi(i)=9999d0
         xma(i)=-9999d0
         yma(i)=-9999d0
         do 20 j=1,nb(i)
           xmi(i)=dmin1(xmi(i),yb(i,j))
           ymi(i)=dmin1(ymi(i),xb(i,j))
           xma(i)=dmax1(xma(i),yb(i,j))
           yma(i)=dmax1(yma(i),xb(i,j))
 20    continue
       do 30 i=1,nr
       xmin=dmin1(xmi(i),xmin)
       xmax=dmax1(xma(i),xmax)
       ymin=dmin1(ymi(i),ymin)
       ymax=dmax1(yma(i),ymax)
 30    continue
       write(*,*)"yma(2)=",yma(2),"yma(1)=",yma(1)
       write(*,*)"ymi(2)=",ymi(2),"ymi(1)=",ymi(1)
       write(*,*)"xma(2)=",xma(2),"xma(1)=",xma(1)
       write(*,*)"xmi(2)=",xmi(2),"xmi(1)=",xmi(1)
       write(*,*)"ymin=",ymin,"xmin=",xmin
       write(*,*)"ymax=",ymax,"xmax=",xmax
c
       ny=int(yma(2)/(rr*2d0)*dble(ny1-ny0)+5d-1)+ny0
       ww=dmin1(0.8d0+(xma(2)-1d0)*2d-1,1.5d0)
       nx=min0(40,int(xma(2)/yma(2)*dble(ny)+5d-1))
       dx=xma(2)/dble(nx)*ww
       dy=yma(2)/dble(ny)
       write(*,*) 'dx,dy',dx,dy
c       dxg=dy*0.15d0
c       dxh=dy*0.08d0
c       dyh=dx/2d0
       dxg=dy*0.15d0*2d0
       dxh=dy*0.08d0*2d0
       dyh=dx/2d0*2d0
c--------------------------mesh for xaxis-----------------
c      write(*,*) 'mesh start yzmin0=',yzmin0
c -----
       x=xmin
       wx=0.011d0
       xg(1)=x-wx
       ii=2
       iorg=ii 
       xg(ii)=x
       xg(ii+1)=-5d-2
       xg(ii+2)=0d0
       x=0d0
       do 100 i=ii+3,ixg-1
c         wx=dxg*1.1d0      !下半分
         wx=dxg*0.65d0      !下半分
c
         if(i.le.(ii+2)+2)wx=wx/2d0    !最初の2本
c
c         if(x.ge.(xma(2)-xmi(2))*5d-1+xmi(2))wx=dxh/1.5d0 
c         if(x.ge.(xma(2)-xmi(2))*5d-1+xmi(2))wx=dxg*0.35d0  !上半分
         if(x.ge.(xma(2)-xmi(2))*5d-1+xmi(2))wx=dxg*0.45d0    !上半分
         x=x+wx
         xg(i)=x
         if(x.gt.xma(2)-1.5d0*wx) then
           nxs=i
           goto 200
         endif
 100  continue
c
200      write(*,*)"xg(",i,")",xg(i)
      xg(i+1)=xma(2)
      xg(i+2)=xma(2)+wx
      xg(i+3)=xma(2)+wx*2d0
      xg(i+4)=xma(2)+wx*3d0
      xg(i+5)=xma(2)+wx*4d0
      write(*,*)"xg",xg(i+1),xg(i+2)
      i00=i+5
      x=xg(i00)
      wx=0.11d0
      do i=i00+1,ixg
	  x=x+wx
	  xg(i)=x
         write(*,*)"xg(",i,")",xg(i)
         if(x.gt.xmax) then
           nxs=i
           if(dabs((xmax-xg(i-1))/(xg(i)-xg(i-1))).le.0.5d0) then
             xg(i-1)=xg(i-2)+wx/2d0
             xg(i)=xmax
             i01=i
             goto 250
            else
             xg(i)=xmax
             i01=i
             goto 250
           endif
         endif
      enddo
250   xg(i01+1)=xg(i01)+wx
c  -----
      nx=i01+1
c       if(xma(2).gt.5d-1*1.2d0) then
c         n0=in1
c         nm=(n0+nxs)/2
c         x0=xg(n0)
c         x1=xg(nxs)
cc         xm=xg(nm)
c         xm=(x0+x1)/2d0
c         k0=nm-n0+1
c         k1=nxs-nm+1
c         if(k0.gt.2 .and. k1.gt.2) then
c           etta=dmax1(eta,.5d0)
c           call fcd(k0,x0,(xm-x0)*(1d0-etta)+x0,xm,xy)
c         do 500 i=2,k0-1
c 500     xg(i+n0-1)=xy(i)
c         call fcd(k1,xm,(x1-xm)*etta+xm,x1,xy)
c         do 600 i=1,k1-1
c 600     xg(i+nm-1)=xy(i)
c         x0=xg(nxs)
c         x1=xg(nx)
c         k=nx-nxs+1
c         call fcd(k,x0,(x1-x0)*(1-etta)+x0,x1,xy)
c         do 700 i=1,k-1
c 700     xg(i+nxs-1)=xy(i)
c         endif
c       endif
c
       wy=9d-2
       yg(1)=ymin-wy       !j region is between 1 and j00
       yg(2)=ymin
       yg(3)=-5d-2
       yg(4)=0d0
       jorg=2
       y=yg(4)
       if(iyy.eq.1)then
         wyy=(yma(2)-ymi(2))/4d0/8d0/2d0*8d-1
c         wyy=(yma(2)-ymi(2))/4d0/8d0*2d0
         wy=wyy*2d0
       else
         wyy=(yma(2)-ymi(2))/8d0/dble(dwwy)/2d0*8d-1
c         wyy=(yma(2)-ymi(2))/8d0/dble(dwwy)*2d0
         wy=wyy/2d0
       endif
       write(*,*)"wy y yma(2)",wy,y,yma(2)
       do 800 j=5,iyg     !j region is between 1 and j01
       wytt=wy
c
       if(j.le.4+2)wytt=wytt/2d0   !最初の2本
c
       y=y+wytt
       yg(j)=y
c       if(yg(j).gt.yma(2)*0.75d0) goto 1000
        if(yg(j).gt.yma(2)*0.6d0) goto 1000
 800   continue
1000   j01=j
       write(*,*)"j01",j01
       write(*,*)"until yma(2)*0.6d0"
       write(*,*)"yg",yg(j01)
c
       if(iyy.eq.1)then
         wy=wyy
       else
         wy=wyy/2d0/2d0
       endif
       do 900 j=j01+1,iyg            !j region is between j04+1 and j05
       y=y+wy
       yg(j)=y
       if (y.ge.yma(2)+1d-12) then
         if (dabs((yma(2)-yg(j-1))/(yg(j)-yg(j-1))).le.0.5d0) then
           yg(j-1)=yg(j-2)+wy/2d0
           yg(j)=yma(2)
           j02=j
           goto 1050
         else
           yg(j)=yma(2)
           j02=j
           goto 1050
         endif
       endif
 900   continue
 1050  y=yma(2)
c
       write(*,*)"j02",j02
       write(*,*)yma(2)+1.5d0,yma(1)-1.5d0*wy
       wy=9d-2
       do 1100 j=j02+1,iyg          !j region is between j05+1 and j06
       y=y+wy
       yg(j)=y
       if (y.gt.dmin1(yma(2)+1.5d0,yma(1)-1.5d0*wy)) goto 1200
c       if(y.gt.(yma(1)-2d-1)) goto 1200
 1100  continue
 1200  j03=j
       yg(j03+1)=yma(1)
       yg(j03+2)=yma(1)+wy
       ny=j03+2
c
c        do 1300 ii=1,2
c        do 1400 i=2,nx
c 1400   if (xg(i).ge.xma(ii)) goto 1410
c 1410   if (dabs((xma(ii)-xg(i-1))/(xg(i)-xg(i-1))).le.0.5d0) then
c          xg(i-1)=xma(ii)
c          write(*,*)"xg(",i-1,")",xg(i-1)
c          goto 1420
c        else
c          xg(i)=xma(ii)
c          write(*,*)"xg(",i,")",xg(i)
c          goto 1420
c        endif
c 1420   continue
c 1300   continue
c-------------------------------------------------------------
       nxg=nx
       nyg=ny
       write(*,*)"nxg,nyg",nxg,nyg
c
       open(1,file="xgyg1.res")
       rewind(1)
       write(1,*)"about xg nxg=",nxg
       do i=1,nxg
       write(1,*)"xg(",i,")=",xg(i)
       enddo
       write(1,*)"about yg nyg=",nyg
       do i=1,nyg
       write(1,*)"yg(",i,")=",yg(i)
       enddo
       close(1)
c
       do 1500 i=1,ig
       do 1500 k=1,2
1500   jg(k,i)=0
       jg(1,iorg)=1
       jg(2,jorg)=1
c
       write(*,*)"before insert/nxg,nyg",nxg,nyg
       do 1600 i=1,2
       call insertgi(nxg,xg,xmi(i),1,jg,ixg,ig)
       call insertgi(nyg,yg,ymi(i),2,jg,iyg,ig)
       call insertgi(nxg,xg,xma(i),1,jg,ixg,ig)
       call insertgi(nyg,yg,yma(i),2,jg,iyg,ig)
1600    continue
       write(*,*)"after insert/nxg,nyg",nxg,nyg
c
       open(1,file="xgyg2.res")
       rewind(1)
       write(1,*)"about xg nxg=",nxg
       do i=1,nxg
       write(1,*)"xg(",i,")=",xg(i)
       enddo
       write(1,*)"about yg nyg=",nyg
       do i=1,nyg
       write(1,*)"yg(",i,")=",yg(i)
       enddo
       close(1)
c
       return
       end
c

