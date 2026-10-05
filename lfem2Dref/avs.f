        subroutine avs3dm(nelt,npt,ne,idx,xx,yy,zz,uu,vv,ww,pp,
     &                  ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        integer*2 nes(21,1001,1001)    !
        dimension ne(ine,3),idx(0:ine),xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
c        dimension nes(21,1001,1001)
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),     !%(6.25)
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
        write(*,*)"in avs3dm"    !
c        write(lpp,'(f8.1)') ti*1d3
        write(lpp,'(f8.1)') ti*1d4
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-3d-'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
        do 300 i=1,npt
        zz(i)=0d0
 300    ww(i)=0d0
c surface
c        pi=3.1415926/2.
        pi=3.1415926535897932384d0/2d0   !
        ist=11
c droplet surface
        nps=0
        nels=0

        do 320 i=1,nr
          if(ndx(i).eq.-1) goto 320
          nels=nels+(ist-1)*nb(i)
         do 330 j=1,ist
c            s=pi/float(ist-1)*float(j-1)
            s=pi/dble(ist-1)*dble(j-1)   !
            do 330 k=1,nb(i)
              nps=nps+1
              n=ne(nelb(i,k),neb(i,k))
              nes(i,j,k)=nps+npt
              xx(nps+npt)=dcos(s)*xx(n)
              yy(nps+npt)=yy(n)
              zz(nps+npt)=-dsin(s)*xx(n)
              uu(nps+npt)=dcos(s)*uu(n)
              vv(nps+npt)=vv(n)
              ww(nps+npt)=-dsin(s)*uu(n)
 330        pp(nps+npt)=1d0
 320    continue
        if(npt+nps.gt.ind .or. nelt+nels.gt.ine) then
        write(*,*) 'avs/error ind <',npt+nps,' or ine <',nelt+nels
        stop
        endif
c
        write(*,*)"logging started in avs3dm"   !
        open(3,file=vfile)
c        rewind(3)
        rewind(3)
c        write(3,*) '#',1         !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt+nps,nelt+nels
        write(3,*)  npt+nps,nelt+nels,4,0,0     !
        do 510 i=1,npt+nps
 510    write(3,'(i4,3e12.4)') i,-sngl(xx(i)),sngl(zz(i)),sngl(yy(i))
        do 540 i=1,nelt
 540    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=1,3)
        n=nelt
        do 560 i=1,nr
        if(ndx(i).eq.-1) goto 560
        do 570 j=1,ist-1
        do 570 k=1,nb(i)
        n=n+1
 570    write(3,'(2i4,a5,4i5)') n,ndx(i)+10-idmin,'quad'
     &          ,nes(i,j,k),nes(i,j,mcy(k+1,nb(i)))
     &          ,nes(i,j+1,mcy(k+1,nb(i))),nes(i,j+1,k)
 560    continue
c        write(3,*) 4,0      !vanished 
        write(3,*) 2,3,1
        write(3,*) 'velocity,'
        write(3,*) 'pressure,'
        do 580 i=1,npt+nps                                  !%hara(6.25)
 580    write(3,'(i4,10e12.4)')
     &          i,-sngl(uu(i)),sngl(ww(i)),sngl(vv(i)),sngl(pp(i))
        close(3)
        write(*,*)"loggind endedn in avs3dm"   !
c
        return
        end
c
        subroutine avs2dm(nelt,npt,ne,idx,xx,yy,zz,uu,vv,ww,pp,
     &                  ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),idx(0:ine),xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
        write(*,*)"in avs2dm"    !
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-2d-'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
c       do 300 i=1,npt
c       zz(i)=0d0
c 300   ww(i)=0d0
c
        write(*,*)"logging started in avs2dm"  !
        open(3,file=vfile)
        write(*,*)"file opened in avs2dm"   !
        rewind(3)
c        write(3,*) '#',1            !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt,nelt,
        write(3,*)  npt,nelt,4,0,0    !
        do 510 i=1,npt
 510    write(3,'(i4,3e12.4)') i,-sngl(xx(i)),0e0,sngl(yy(i))

        do 540 i=1,nelt
 540    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=1,3)
c       write(3,*) 4,0				!vanished 
        write(3,*) 2,3,1
        write(3,*) 'velocity,'
        write(3,*) 'pressure,'
        do 580 i=1,npt
 580    write(3,'(i4,10e12.4)')
     &          i,-sngl(uu(i)),0e0,sngl(vv(i)),sngl(pp(i))

        close(3)
        write(*,*)"logging ended in avs2dm"   !
c
        return
        end
c
        subroutine avs3dmmas(nelt,np2,ne,idx,xx,yy,zz,tt,
     &									ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        integer*2 nes(21,1001,1001)    !
        dimension ne(ine,6),idx(0:ine),xx(ind),yy(ind),zz(ind)
        dimension tt(ind),ip(ind)	        !
	  dimension xtmp(ind),ytmp(ind),ztmp(ind),ttmp(ind)
c        dimension nes(21,1001,1001)
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),     !%(6.25)
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                   ncomw(21,1001),ncopw(21,1001)

c
        write(*,*)"in avs3dmmas"    !
c        write(lpp,'(f8.1)') ti*1d3
        write(lpp,'(f8.1)') ti*1d4
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-3d-mas'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
c surface
        pi=3.1415926/2.	     !保存するデータは単精度
c        pi=3.14159265358979d0/2d0     !32384d0/2d0  !倍精度の仮数部は約17桁
        ist=11
c get nodes
        ntmp=np2
	  do i=1,ntmp
	    xtmp(i)=xx(i)
	    ytmp(i)=yy(i)
	    ztmp(i)=zz(i)
	    ttmp(i)=tt(i)
	  enddo
        np2=0
        do 1100 i=1,nelt
	  do 1000 j=1,3
	    do k=1,np2
	    if(ip(k).eq.ne(i,j))goto 1000
	    enddo
	    np2=np2+1
	    ip(np2)=ne(i,j)
1000    continue
1100    continue
c
c droplet surface
        write(*,*)"got nodes"
        nps=0
        nels=0
        do 320 i=1,nr
          if(ndx(i).eq.-1) goto 320
          nels=nels+(ist-1)*nb(i)
         do 330 j=1,ist
c            s=pi/float(ist-1)*float(j-1)
            s=pi/dble(ist-1)*dble(j-1)   !
            do 330 k=1,nbw(i),2
              nps=nps+1
              n=ne(nelbw(i,k),nebw(i,k))
              nes(i,j,k)=nps+np2
              xx(nps+np2)=dcos(s)*xx(n)
              yy(nps+np2)=yy(n)
              zz(nps+np2)=-dsin(s)*xx(n)
              tt(nps+np2)=tt(n)
 330          ip(nps+np2)=nps+np2
 320    continue
        if(np2+nps.gt.ind .or. nelt+nels.gt.ine) then
        write(*,*) 'avs/error ind <',np2+nps,' or ine <',nelt+nels
        stop
        endif
c
        write(*,*)"loggin started in avs3dmmas"
	  write(*,*)"np2,nps,nelt,nels",np2,nps,nelt,nels
        open(3,file=vfile)
        rewind(3)
c        write(3,*) '#',1         !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt+nps,nelt+nels
        write(3,*)  np2+nps,nelt+nels,1,0,0     !
        do 510 i=1,np2+nps
510     write(3,'(i4,3e12.4)') ip(i),-sngl(xx(ip(i))),sngl(zz(ip(i))),
     &                                sngl(yy(ip(i)))
        do 540 i=1,nelt
 540    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=1,3)
        n=nelt
        do 560 i=1,nr
        if(ndx(i).eq.-1) goto 560
        do 570 j=1,ist-1
        do 570 k=1,nbw(i),2
        n=n+1
 570    write(3,'(2i4,a5,4i5)') n,ndx(i)+10-idmin,'quad'
     &          ,nes(i,j,k),nes(i,j,mcy(k+1,nb(i)))
     &          ,nes(i,j+1,mcy(k+1,nb(i))),nes(i,j+1,k)
 560    continue
c        write(3,*) 4,0      !vanished 
        write(3,*) 1,1
        write(3,*) 'concentratin,'
        do 580 i=1,np2+nps                                  !%hara(6.25)
 580    write(3,'(i4,10e12.4)') ip(i),sngl(tt(ip(i)))
        close(3)
c
	  np2=ntmp
	  do i=1,np2
	    xx(i)=xtmp(i)
	    yy(i)=ytmp(i)
	    zz(i)=ztmp(i)
	    tt(i)=ttmp(i)
	  enddo
	  write(*,*)"loggin ended in avs3dmmas"
c
        return
        end
c
        subroutine avs2dmmas(nelt,np2,ne,idx,xx,yy,zz,tt,
     &                  ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,6),idx(0:ine),xx(ind),yy(ind),zz(ind)
        dimension tt(ind),ip(ind)    !
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
        write(*,*)"in avs2dmmas"    !
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-2d-mas'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
c       do 300 i=1,npt
c       zz(i)=0d0
c 300   ww(i)=0d0
c get nodes
        ntmp=np2
        np2=0
        do 1100 i=1,nelt
	  do 1000 j=1,3
	    do k=1,np2
	    if(ip(k).eq.ne(i,j))goto 1000
	    enddo
	    np2=np2+1
	    ip(np2)=ne(i,j)
1000    continue
1100    continue
        open(1,file="avs.res")
	  rewind(1)
	  write(1,*)"np2=",np2
	  do i=1,np2				 
	  write(1,*)i,",",ip(i),",",sngl(xx(ip(i))),",",sngl(yy(ip(i)))
     &                                           ,",",sngl(tt(ip(i)))
	  enddo
	  close(1)
c
        open(3,file=vfile)
        rewind(3)
c        write(3,*) '#',1            !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt,nelt,
        write(3,*)  np2,nelt,1,0,0    !
        do 510 i=1,np2
 510    write(3,'(i4,3e12.4)')ip(i),-sngl(xx(ip(i))),0e0,sngl(yy(ip(i)))

        do 540 i=1,nelt
 540    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=1,3)
c       write(3,*) 4,0				!vanished 
        write(3,*) 1,1
        write(3,*) 'concentration,'
        do 580 i=1,np2
 580    write(3,'(i4,10e12.4)')ip(i),sngl(tt(ip(i)))
c
        close(3)
	  np2=ntmp
c
        return
        end
c
        subroutine save2el(nel,np2,ne2,idx,xx2,yy2,tt,ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne2(ine,6),idx(0:ine),xx2(ind),yy2(ind),tt(ind)
        character lpp*8,rfile*30,vfile*30
	  common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
	  common /dprop/dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
c
        write(*,*)"in save2el"
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-2el'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        write(*,*)"save started"    !
	  open(3,file=vfile)
	  rewind(3)
	  write(3,*)ti,np2,nel
        do 510 i=1,np2
 510    write(3,'(i4,3e12.4)')i,sngl(xx2(i)),sngl(yy2(i)),sngl(tt(i))
        do 540 i=1,nel
 540    write(3,'(2i4,6i6)') i,idx(i),(ne2(i,j),j=1,6)
        close(3)
c
        lr=2     !自由表面の平均濃度を求める
	  ccc=0d0
	  nb=0
	  do 200 i=1,nbw(lr)
	  if(ncomw(lr,i).ne.0)goto 200
        ccc=ccc+tt(ne2(nelbw(lr,i),nebw(lr,i)))
	  nb=nb+1
	  if(ncomw(lr,i+i).ne.0)then   !最後の節点
          ccc=ccc+tt(ne2(nelbw(lr,i+1),nebw(lr,i+1)))
          nb=nb+1
	  endif
200     continue
        ccc=ccc/dble(nb)
        open(3,file='suf-cav.res',access='append')
	  write(3,*)sngl(ti),nb,sngl(ccc)
	  close(3)
c
        return
        end
c
        subroutine avs2dmmasr(nelt,np2,ne,idx,xx,yy,zz,tt,
     &                  ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,6),idx(0:ine),xx(ind),yy(ind),zz(ind)
        dimension tt(ind),ip(ind)    !
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
        write(*,*)"in avs2dmmasr"    !
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-2d-masr'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
c       do 300 i=1,npt
c       zz(i)=0d0
c 300   ww(i)=0d0
c get nodes
        ntmp=np2
        np2=0
        do 1100 i=1,nelt
	  do 1000 j=4,6
	    do k=1,np2
	    if(ip(k).eq.ne(i,j))goto 1000
	    enddo
	    np2=np2+1
	    ip(np2)=ne(i,j)
1000    continue
1100    continue
        open(1,file="avs.res")
	  rewind(1)
	  write(1,*)"np2=",np2
	  do i=1,np2				 
	  write(1,*)i,",",ip(i),",",sngl(xx(ip(i))),",",sngl(yy(ip(i)))
     &                                           ,",",sngl(tt(ip(i)))
	  enddo
	  close(1)
c
        open(3,file=vfile)
        rewind(3)
c        write(3,*) '#',1            !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt,nelt,
        write(3,*)  np2,nelt,1,0,0    !
        do 510 i=1,np2
 510    write(3,'(i4,3e12.4)')ip(i),-sngl(xx(ip(i))),0e0,sngl(yy(ip(i)))
        do 540 i=1,nelt
 540    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=4,6)
c       write(3,*) 4,0				!vanished 
        write(3,*) 1,1
        write(3,*) 'concentration,'
        do 580 i=1,np2
 580    write(3,'(i4,10e12.4)')ip(i),sngl(tt(ip(i)))
c
        close(3)
	  np2=ntmp
c
        return
        end