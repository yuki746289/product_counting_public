c        include "\usr\flib\addchr.for"
        include 'addchr.for'	       !ファイル名に拡張子を加える
c        include 'matrix1.f'
c        include 'flow.f'
c        include 'grid.f'
c        include 'gmesh-dropletSpread2.f' 
c        include 'gmesh-dropletSpread3.f' 
c        include 'gmesh-dropletVapor.f' 
c        include 'thinmesh3.f'
c        include 'renum.f'
c        include 'coalesce.f'
c        include 'plotmesh.f'
c        include 'avs.f'
c        include 'renum2.f'
c        include 'choleski.f'
c        include 'mass.f'
c        include 'therm.f'	         !物性値の計算
c        include 'temp.f'
        implicit double precision (a-h,o-z)
        parameter(ind=28001,ine=18001,ibw=501,ibt=251,lbw=27,ica=3)
c        parameter(ind=8001,ine=5001,ibw=501,ibt=251,lbw=27,ica=3)
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine),irx(0:ine)
        dimension nen(ine,3),xn(ind),yn(ind),idn(0:ine),irn(0:ine)
        dimension new(ine,3),xw(ind),yw(ind),idw(0:ine),irw(0:ine)
        dimension jnt(ica,5,ind),jnw(ind),ar0(ine)
        dimension mcr(ica),ncr(ica,5),kcr(ica,5,21)
        dimension npc(ica,5),nbwc(ica,5),nelrs(21),nelre(21)
        dimension kcrg(21),kbd(1001),bdv(1001)
        dimension sa(ind,ibw),sf(ind),dsf(ind)
        dimension x0(ind),y0(ind)
        dimension uu(ind),vv(ind),pp(ind),uvp(ind)
        dimension u0(ind),v0(ind),p0(ind),uvp0(ind)
        character*30 rfile,rfile1,rfile2,rfile3,tmfile
        character*4 fd1,fd2,fd3
c*******物質移動現象***************************************************
	  dimension tt(ind),tt0(ind),tt22(ind),ttk(ind)     !
	  dimension ne2(ine,6),xx2(ind),yy2(ind)
	  dimension ne22(ine,6),xx22(ind),yy22(ind)
	  dimension nek(ine,6),xxk(ind),yyk(ind)
	  dimension kuv(2001,10),uvb(2001,10),nuv(10)
	  dimension sa2(ind,ibw),sf2(ind)
	  dimension jww(ind),jwwb(ind),ip(ind),jp(ind)
	  dimension jbw(ind),px2(ind),py2(ind)  !boundfsで法線ベクトルの方向余弦を格納する時に用いる
        dimension js2(3)
c	  dimension tmpx(ind),tmpy(ind)
        character lpp*8      !tmp fileを各時間ごとに違う名前で保存する為に用いる        
c*********************************************************************
c*******熱移動,1次要素の変数******************************************
        dimension net(ine,3),xxt(ind),yyt(ind)    !温度補間に用いる
        dimension tp(ind),tp0(ind),tpn(ind)       !tpn(ind):温度補間に用いる
        dimension sat(ind,ibw),sft(ind)
        dimension jwt(ind),jwtt(ind),jwtr(ind),jwtrr(ind),jwtrrr(ind)
        character ltp*8
c*********************************************************************
        common /dplt/ lvp
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /regcon/ ncon(21),kcon(21,21)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
        common /dgrd/ dx,dy,in0,in1,jn0,jn1
	  common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)   !added sentence
c-------物質移動,2次要素のcommon文----------------------
        common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
	  common /mass/ dmas(9)
c-------------------------------------------------------
c        data msw/4/,ratio/2/,zkasu/1.0/
        data msw/4/,ratio/2/,zkasu/1.0/
        data fd1,fd2,fd3/'.dat','.res','.tmp'/
        data rfile/'impac.dat'/
c-------変数初期化--------------------------------------
        data js2/2,3,1/
        ti=0d0
        idx(0)=0
        irx(0)=0
c-------物質移動----------------------
        area0=0d0
	  area1=0d0
	  ttw=1d0
	  tt00=0d0       !初期濃度:CA0 [-]
	  irec=1
	  cdv=0d0
	  cdvp=0d0        
	  pi=3.145926535897d0
c-------変数初期化(熱移動)-----------------
        tpini=tini(1)/tplate(1)        !初期液滴内温度:Θ=Tini/Tw
c---------------------------------- MS/CQ Fortran only
        open(1,file='input.txt')
        read(1,'(a)') rfile
        close(1)
        open(5,file=rfile)
c----------------------------------
        write(*,*) '0=plot, 1=no-plot'
        read(5,*) lvp
        write(*,*) 'avs (0, 2, 3 or 5)'
        read(5,*) iavs
        write(*,*) 'start at t=0 or not ? (1/0)'
        read(5,*) iyy
        write(*,*) 'dt,wdt1(dat),wdt2(res),wdt3(tmp),timax'
        read(5,*) dt,wdt1,wdt2,wdt3,timax
        write(*,*) 'dx,dy,in0,in1,jn0,jn1'
        read(5,*) dx,dy,in0,in1,jn0,jn1
        write(*,*) 'results file'
        read(5,'(a)') rfile
c
        call data(nel,np,ne,xx,yy,idx,ind,ine,ibw,dltm)
c-------------物質移動,2次要素------------------------
        id=1
        pe=(1d0/visc(id))*dmas(id)
c-----------------------------------------------------
        if(iyy.eq.1) then
          close(5)
          call sufele(nel,ne,xx,yy,idx,xx,yy,ine,ind)
          write(*,*)"sufele ended"     !
          call genirx(nel,irx,idx,ine)
          write(*,*)"genirx ended"     !
          icoal=1
          wtime1=0d0
          wtime2=0d0
          wtime3=0d0
        else
          write(*,*) 'temporary file'
          read(5,'(a)') tmfile
          close(5)
          call addchr(tmfile,fd3,rfile3)
          open(4,file=rfile3,status='unknown',access='sequential',
     $            form='unformatted')
          read(4) ti,hmax,icoal,wtime1,wtime2,wtime3
          read(4) np,(xx(i),yy(i),i=1,np)
          read(4) np,(uu(i),vv(i),pp(i),i=1,np)
          read(4) nel,((ne(i,j),j=1,3),idx(i),irx(i),ar0(i),i=1,nel)
          read(4) nr,(nb(i),ndx(i),(nelb(i,j),neb(i,j),j=-1,nb(i)+2),
     &          (ncom(i,j),ncop(i,j),j=1,nb(i)+1),i=1,nr)
          read(4) (ncon(i),(kcon(i,j),j=1,nr),i=1,nr)
          read(4) nset,(kset(i),i=1,nr)
          read(4) (nelrs(i),nelre(i),i=1,nr)
          read(4) msw,(mcr(i),(ncr(i,j),npc(i,j),nbwc(i,j),j=1,mcr(i))
     &          ,i=1,msw)
          read(4) (((kcr(i,j,k),k=1,ncr(i,j)),j=1,mcr(i)),i=1,msw)
          read(4) (((jnt(i,j,k),k=1,np),j=1,mcr(i)),i=1,msw)
          read(4) ncp,(kcp(i),icp(i),jcp(i),st(i),ancl(i),
     &                  icl(i),icla(i),staa(i),strr(i),i=1,ncp)
c-----------------物質移動,2次要素のデータ-------------------------
	    read(4) np2,(xx2(i),yy2(i),i=1,np2)  
	    read(4) np2,(tt(i),i=1,np2)
	    read(4) nel,((ne2(i,j),j=1,6),i=1,nel)
          read(4) nr,(nbw(i),(nelbw(i,j),nebw(i,j),j=-1,nbw(i)+2),
     &           (ncomw(i,j),ncopw(i,j),j=1,nbw(i)+2),i=1,nr)
	    read(4) tt00,ttw,csta,time1												   
	    read(4) cdv,area0,aaa,ccc1,nbw2
c------------------------------------------------------------------
c-----------------熱移動-------------------------------------------
          read(4) np,(tp(i),i=1,np)
c------------------------------------------------------------------
c-----------------計算の時間間隔dt---------------------------------
          read(4) imesh,kmesh,dt
c------------------------------------------------------------------
          close(4)
c
c--------蒸発を考慮に入れた時間をti=0d0とする(初回のみ呼ぶ)---------------------------
c--------sub initialic(速度が0になる),sub cloalesce(icoal=0になる)にコメントをつける--
c         ti=0d0
c          wtime1=0d0
c          wtime2=0d0
c          wtime3=0d0
c          icoal=1	     !mesh切り直しの為
c	    iyy=1        !初期濃度分布，吸収量を計算する為
c	    time1=0d0    !まだ濃度計算は行っていないのでtime1(前回のど計算を行った時間)を0d0とする
c
c----------------------------------------------------------------------
c
c          call pltnt(nel,ne,xx,yy,idx,ind,ine)
          call coalesce(nel,ne,xx,yy,irx,idx,icoal,ind,ine)
        endif
c
          wtime1=wtime1+wdt1
          wtime2=wtime2+wdt2
          wtime3=wtime3+wdt3
c
cc        write(*,*)"start pltnt"      !
cc        call pltnt(nel,ne,xx,yy,idx,ind,ine)
cc        write(*,*)"ended pltnt"      !
c        write(*,*) ti,nel,np,hmax,icoal
c       pause
c
c       icoal=1
c
c-----------------------mesh再発生の処理-----------------------------
c
        if (hmax.gt.ratio .or. icoal.eq.1) then
c         call sufele(nel,ne,xx,yy,idx,ine,ind)
c-----メッシュ再発生------------------------------
          write(*,*)"grd started in if"   !
          call grd(nel,np,ne,xx,yy,idx,irx,jnt,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          new,xw,yw,idw,irw,jnw,
     &          nelrs,nelre,npc,nbwc,mcr,ncr,kcr,
     &          ar0,ind,ine,ibw,ica)
c-----メッシュ再発生終了----------------------------
c-----熱移動(iyy=0の場合における温度補間に用いる)---
          nelt=nel
          npt=np
          do i=1,nelt
          do j=1,3
            net(i,j)=ne(i,j)
            xxt(net(i,j))=xx(net(i,j))
            yyt(net(i,j))=yy(net(i,j))
          enddo
          enddo
c-------熱移動終了----------------------------------
c-------速度補間------------------------------------
          write(*,*)"grdfcn started in if"      !
          call grdfcn(nel,np,ne,xx,yy,idx,irx,jnw,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          uu,pp,vv,u0,v0,p0,ind,ine)
c-------速度補間終了--------------------------------
c
c---------mesh切り直し後の時間間隔を決定--------------------------------
          imesh=1
          kmesh=1
c-----------------------------------------------------------------------
c
c---------物質移動，2次要素--------------------------------
          write(*,*)"generate2el started in if"       !
c	    call generate2el(np,ne,nel,xx,yy,ine,ind,
c     $                 np22,ne22,xx22,yy22,nbw22)
          call genet2el(np,nel,ne,xx,yy,
     &                  np22,ne22,xx22,yy22,nbw22)
	    call sufelew
	    write(*,*)"np22=",np22,"nbw(nr)=",nbw(1),nbw(2)
c---------物質移動,2次要素の処理終了-----------------------
c
c---------iyy=1(ti=0の場合)---------------------------------------------
	    if(iyy.eq.1)then      !iyy=1
c-----------------------------------------------------------
c
c---------物質移動,2次要素--------------------------------------
c
	      write(*,*)"calc mass in iyy=1"
	      np2=np22
	      do i=1,nel
	      do j=1,6
	        ne2(i,j)=ne22(i,j)
	      enddo
	      enddo
	      do i=1,np2
	        xx2(i)=xx22(i)
	        yy2(i)=yy22(i)
	        tt(i)=tt00
	      enddo
c
c
          open(1,file="ne2.res")
	    rewind(1)
	    write(1,*)"np2=",np2
	    do i=1,nel
	    do j=1,6
	    write(1,*)i,",",j,",",ne2(i,j)
	    enddo
	    enddo
	    close(1)
c
c
	      write(*,*)"sub bunpu"
c	      if(dltm.ne.0d0) call bunpu(np2,tt,xx2,yy2,dltm,ind)
            call bunpu2(np2,tt,ind)     !液滴蒸発
	      write(*,*)"sub btw"
c	      call btw(ne2,nuv,kuv,uvb,tt,ine,ind)
	      write(*,*)"sub calc"
	      call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,csta,ind,ine)
	      write(*,*)"csta=",sngl(3d0/8d0/pi*csta)
	      call sufar(ne,xx,yy,aaa,ind,ine)
	      area0=aaa
            open(33,file="callog.res",access='append')
	      write(33,*)"ti=",sngl(ti),"   iyy=1"
	      close(33)
	      call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,tt,   !
     &                  ti,rfile,ine,ind)
	      call avs2dmmasr(nel,np2,ne2,idx,xx2,yy2,xn,tt,   !
     &                  ti,rfile,ine,ind)
	      call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                                   ti,rfile,ine,ind)
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !蒸発開始時からの物質の増加量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),   !aaa(2πRH):液滴表面積, vvv(πR^2H):液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)          !ccc:増加量, ccct:物質量
	  write(4,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	 !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)    !ccc/vvv:平均濃度の初期濃度との差，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c---------------------------------------------------------------
c---------物質移動,2次要素の処理終了----------------------------
c
c--------------熱移動-------------------------------------------
c
        do i=1,nel
        do j=1,3
          tp(ne(i,j))=0d0
          if(idx(i).eq.1)tp(ne(i,j))=tpini              !index=1の領域の節点に初期温度を代入
          if(idx(i).eq.-1 .and. xx(ne(i,j)).gt.-1d-15   !平板に温度代入(R≧0,Z≦0)
     &                    .and. yy(ne(i,j)).lt.1d-15) 
     &                   tp(ne(i,j))=tplate(1)/tplate(1)
        enddo
        enddo
        call avs2dmtmp(nel,np,ne,idx,xx,yy,tp,tt,
     &                  0d0,rfile,ine,ind)
c
c--------------熱移動終了の処理終了-----------------------------
c
c-------------iyy=0(ti≠0の場合)---------------------------------------------
	    else				  !iyy=0
c---------------------------------------------------------------
c
c
c--------------熱移動(温度補間---------------------------
c
          write(*,*)"grdfct started"
          call grdfct(nelt,npt,net,xxt,yyt,tp,
     &                    neln,npn,nen,xn,yn,tpn,irx,ind,ine)
          do i=1,nel
          do j=1,3
            if(idx(i).eq.-1 .and. xx(ne(i,j)).gt.-1d-15      !平板に温度代入(R≧0,Z≦0)
     &                      .and. yy(ne(i,j)).lt.1d-15) 
     &                     tp(ne(i,j))=tplate(1)/tplate(1)
          enddo
          enddo
c
c--------------熱移動終了の処理終了----------------------
c
c---------物質移動，2次要素--------------------------------
	      write(*,*)"calc mass in iyy=0"
	      neal=nel
	      call grdfc5(neal,np2,nbw2,ne2,xx2,yy2,tt,
     &                  nel,np22,nbw22,ne22,xx22,yy22,tt22,irx,ind,ine)
	      nel=neal
            irec=0
	      call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
c           -------recovery of concentration-------
            if(ccc2.lt.ccc1)then
		  dtt1=dt
		  ccc3=ccc2
1400        continue
	      call getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &                         ip,jp,ine,ind)
	      call renum2(2,npk,nek,nelk,jww,nbw2)
	      do i=1,nelk
	      do j=1,6
	      xxk(jww(nek(i,j)))=xx2(ip(nek(i,j)))
	      yyk(jww(nek(i,j)))=yy2(ip(nek(i,j)))
	      ttk(jww(nek(i,j)))=tt(ip(nek(i,j)))
	      jwwb(jww(nek(i,j)))=ip(nek(i,j))
	      nek(i,j)=jww(nek(i,j))
	      enddo
	      enddo
c
	      do 1450 j=nbw2*2
		  do 1450 i=1,npk
1450        sa2(i,j)=0d0
            do 1500 i=1,npk
1500        sf2(i)=0d0
	    call noudot(nelk,nbw2,nek,xxk,yyk,tp,
     &                 ne2,jww,jp,sa2,sf2,dtt1,ttk,tt0,ind,ine,ibw)	    	   
c	      write(*,*)"sub boundfsw"
c            call boundfsw(ne2,uu,vv,tt,px2,py2,ind,ine)    !自由表面に境界濃度を与える
c            call boundfsw2(nel,ne,xx,yy,ne2,                 !time1:前回濃度計算を行った時間
c     &                        tt,irx,time1,ti,ind,ine)
c	      call btw2(ne2,nuv,kuv,uvb,ttk,jww,jp,ine,ind)
c            write(*,*)"sub btwvap"
c            call btwvap(ne2,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)     !境界条件[K]{C}={f}の[K]を求める
c	      call bound(npk,nbw2,sa2,sf2,nuv,kuv,uvb,ind,ibw)        ![K],{f}をマトリックスに組み込む
	      call choleski2(npk,nbw2,sa2,sf2,ind,ibw)
c            call gauss(npk,nbw2,sa2,sf2,ind,ibw)
	      do 1550 i=npk
1550	      ttk(i)=sf2(i)
            do i=1,npk
	      tt(jwwb(i))=ttk(i)
	      enddo
	      call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
	      if(ccc2.ge.ccc1 .or. ccc2.lt.ccc3)goto 1700
	      dtt1=(ccc1-ccc3)/(ccc2-ccc3)*dtt1/2d0
	      if(dtt1.lt.dt)dtt1=dt
	      write(*,*)"ccc1,ccc2,ccc3=",sngl(ccc1),sngl(ccc2),sngl(ccc3)
	      write(*,*)"dt,dtt1=",sngl(dt),sngl(dtt1)
	      write(*,*)"-----------------"
	      ccc3=ccc2
	      goto 1400
	      endif
1700        continue	    
	      irec=1
	      time1=ti
	      call btw(ne2,nuv,kuv,uvb,tt,ine,ind)
	      inin=1
	    endif
c---------物質移動,2次要素の処理終了---------------------
cc                call pltnt(nel,ne,xx,yy,idx,ind,ine)
cc        write(*,*)"ended pltnt and if sentence"  !
c       pause
        endif
c
c---------mesh再発生の処理終了-----------------------------
c
        if(iyy.eq.1) call initialc(nel,ne,xx,yy,idx,uu,vv,pp,ind,ine)
c               call pltv(np,ne,xx,yy,uu,vv,ine,ind)
c               pause
c
        write(*,*) ' ------- grid data -------'
        write(*,*) 'nset=',nset,'   msw=',msw
        write(*,*) 'element=',nel,'  nodes=',np
c
        write(*,*)"ti=",ti,"nb(1)=",nb(1),"nb(2)=",nb(2)       !
        open(1,file="initial.res")         !
        rewind(1)
        write(1,*)"ti=0d0"
        do i=1,nr
        write(1,*)"nr=",i
        do j=1,nb(i)
        write(1,*)"xx(",ne(nelb(i,j),neb(i,j)),")=",
     &                                     xx(ne(nelb(i,j),neb(i,j))),
     &   "yy(",ne(nelb(i,j),neb(i,j)),")=",yy(ne(nelb(i,j),neb(i,j))),
c     &   "vv(",ne(nelb(i,j),neb(i,j)),")=",vv(ne(nelb(i,j),neb(i,j))),
     &   "ncom(",i,",",j,")=",ncom(i,j)
        enddo
        enddo
        close(1)
c        pause	 
c
        call addchr(rfile,fd1,rfile1)
        call addchr(rfile,fd2,rfile2)
        call addchr(rfile,fd3,rfile3)
c......
        do 60 i=1,np
        u0(i)=uu(i)
        v0(i)=vv(i)
        p0(i)=pp(i)
        x0(i)=xx(i)
        y0(i)=yy(i)
c-------熱移動------------------------------
        tp0(i)=tp(i)
c-------熱移動終了--------------------------
 60     continue
c-------物質移動,2次要素------------------
        do 1800 i=1,np2
1800    tt0(i)=tt(i)
c-----------------------------------------
	  open(1,file="vv.res")
	  rewind(1)
	  do i=1,np
	  write(1,*)sngl(uu(i)),sngl(vv(i)),sngl(pp(i))
	  enddo
	  close(1)
c
c......
!          open(1,file=rfile1,status='unknown',access='sequential',
!     $            form='unformatted')
!          open(2,file=rfile2)
!          rewind(1)
!          rewind(2)
!            write(1) sngl(dt),intv,sngl(timax)
!            write(1) sngl(re),sngl(we),sngl(fr),sngl(sta),
!     &                 sngl(str),sngl(sv2),sngl(rr)
!          if(iyy.eq.1) then
!            write(1) sngl(ti),npd,neld,np,nel
!            write(1) sngl(2*rr),sngl(rr),sngl(xx(ksf(nsf0))),
!     &               sngl(xx(ksf(nsfs))),sngl(vv(kct(nct)))
!            do 80 i=1,nel
! 80         write(1) (int2(ne(i,j)),j=1,3)
!            do 82 i=1,npd
! 82         write(1) sngl(uu(i)),sngl(vv(i)),sngl(pp(i))
!            do 83 i=1,np
! 83         write(1) sngl(xx(i)),sngl(yy(i)),sngl(tt(i))
!            write(1) nsf,nsf0,nsfs,(int2(ksf(i)),i=1,nsf)
!            write(1) nct,(int2(kct(i)),i=1,nct)
c            write(1) nrd,(sngl(sdi(i)),sngl(rdi(i)),i=1,nrd)
c            write(1) nqc,hav,tav,qtt,qht,(xc(i),hc(i),qc(i),i=1,nqc)
!          endif
!          close(1)
!          close(2)
c.....
       do 200 it=1,100000000
c          if(ti.gt.5d-1)dt=1d-4     !
c------mesh 切り直し後の時間間隔-----------------------------
          if(imesh.eq.1)then
            if(kmesh.eq.1)dt=dt/20d0
            if(kmesh.eq.2001)dt=dt*20d0
            if(kmesh.eq.2001)imesh=0
            kmesh=kmesh+1
          endif
c--------------------------------------------------------------
          ti=ti+dt
c
c         do 220 i=1,np
c         xx(i)=x0(i)
c 220     yy(i)=y0(i)
c
          call fst(ne,xx,yy,uu,vv,ine,ind)
c*****************************************************************************************
c
c         各オブジェクトkcを構成する子オブジェクトmcr(kc)ごとに流動，圧力場を計算する
c
c*****************************************************************************************
          do 230 m=1,mcr(1)     !mcr(kc):オブジェクトkcを構成する子オブジェクトの数
            npf=npc(1,m)        !npc(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)を構成する節点の数
            nbwf=nbwc(1,m)      !nbwc(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)のバンド幅
            ncrg=ncr(1,m)       !ncr(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)を構成する領域の数
            nelf=0              !nelf:オブジェクトkcの子オブジェクトmcr(kc)を構成する要素の数
            do 240 i=1,ind
 240        jnw(i)=jnt(1,m,i)   !jnt(kc,mcr(kc),np):オブジェクトkcを構成する最適化した節点の番号
            do 250 n=1,ncrg
              kcrg(n)=kcr(1,m,n)    !kcr(kc,mcr(kc),ncr(kc,mcr(kc))
                                    !オブジェクトkcを構成する領域ncr(kc,mcr(kc)))の番号
              do 260 i=nelrs(kcrg(n)),nelre(kcrg(n))    !nelrs(kr):領域krを構成する要素の最初の番号
                                                        !nelrs(kr):領域krを構成する要素の最後の番号
              nelf=nelf+1
              idw(nelf)=idx(i)
              irw(nelf)=irx(i)
              do 260 j=1,3
                new(nelf,j)=jnw(ne(i,j))               !jnw(古)=新
                xw(new(nelf,j))=xx(ne(i,j))
                yw(new(nelf,j))=yy(ne(i,j))
                uvp0(new(nelf,j)*3-2)=uu(ne(i,j))
                uvp0(new(nelf,j)*3-1)=vv(ne(i,j))
                uvp0(new(nelf,j)*3)=pp(ne(i,j))
                jbw(jnw(ne(i,j)))=ne(i,j)              !jbw(新)=古
 260          continue
 250        continue
c             write(*,*) npf,nbwf,ncrg,nelf
c               call pltv(np,ne,xx,yy,uu,vv,ine,ind)
c               call pltnt(nelf,new,xw,yw,idw,ind,ine)
c               pause
c
c           配列xn(ind),yn(ind)に節点座標を代入する
            do 270 i=1,npf  !オブジェクトkcの子オブジェクトmcr(kc)を構成する節点の数
             xn(i)=xw(i)
             yn(i)=yw(i)
c           原点(0d0,0d0)上に存在する節点の半径方向の速度を零にする
c            if(xn(i).eq.0d0 .and. yn(i).eq.0d0) then
c                nbd=1
c                kbd(nbd)=i*3-1
c                bdv(nbd)=0d0
c            endif
 270        continue
c           配列uvp(ind*3)に速度，圧力を代入する
            do 272 i=1,npf*3
 272        uvp(i)=uvp0(i)
c           配列uvp0(ind*3)に境界条件を代入する
c            call boundf1(nelf,ne,irw,jnw,ncrg,kcrg,
c     &                          nbd,kbd,bdv,uvp0,ind,ine)
            call boundf2(nelf,ne,irw,jnw,ncrg,kcrg,    !境界速度を与える条件追加
     &                          nbd,kbd,bdv,uvp0,ind,ine)
c
            nbd1=0                   !
            do i=1,nb(2)
            if(ncom(2,i).eq.1)then
              nbd1=nbd1+1
              if(ncom(2,i+1).eq.0)nbd1=nbd1+1
            endif
            enddo
c            write(*,*)"nbd=",nbd,"nbd1=",nbd1
c            write(*,*) 'nbd=',nbd
c
            write(*,*)"----------------------"
            do 280 lk=1,26      !xn,yn,xw,yw.new,uvpは節点番号を並べ替えている
c             マトリクスsa(3*ind,ibw),sf(3*ind)を初期化する
              call setmat(npf*3,nbwf*6,sa,sf,dsf,ind,ibw)
c             流動マトリクスを発生させる
c             前回計算した圧力，速度uvp0と移動後の節点xn(npf),yn(npf)を用いて計算する
              call flow(nelf,nbwf*3,sa,sf,new,xn,yn,uvp0,idw,dt,
     &                  ind,ine,ibw)
cc             call bounduvp(nbd,kbd,bdv,uvp,ind)
c             自由表面上に境界条件を格納した配列dsf(ind)を求める
              call boundfs(nelf,ne,xn,yn,tp0,uvp,jnw,tt0,
     &                                 ncrg,kcrg,dsf,ind,ine)
c              call boundfs(nelf,ne,xn,yn,uvp,jnw,ncrg,kcrg,dsf,
c     &                                jbw,px2,py2,ind,ine)
c             境界条件をマトリクスsa(3*ind,ibw)および列ベクトルsf(3*ind)に代入する
c             ↑自由表面および速度零の条件
              call boundasym(npf*3,nbwf*3,sa,sf,dsf,nbd,kbd,bdv,ind,ibw)
c             マトリクスsa(3*ind,ibw)および列ベクトルsa(3*ind)を解く
              call gauss(npf*3,nbwf*3,sa,sf,ind,ibw)
              err=0d0
              do 290 i=1,npf
                 errt=dabs(sf(3*i-2)-uvp(3*i-2))+
     &             dabs(sf(3*i-1)-uvp(3*i-1))+dabs(sf(3*i)-uvp(3*i))
                 if(errt.gt.err)then
                 err=errt
                 kkk=i
c                 write(*,*)kkk,err
c                 write(*,*)dabs(sf(3*i-2)-uvp(3*i-2)),
c     &                     dabs(sf(3*i-1)-uvp(3*i-1)),
c     &                     dabs(sf(3*i-0)-uvp(3*i-0))
c                 write(*,*)"-----------------------------"
c                 write(*,*)sngl(sf(3*i-2)),sngl(uvp(3*i-2))
c                 write(*,*)sngl(sf(3*i-1)),sngl(uvp(3*i-1))
c                 write(*,*)sngl(sf(3*i-0)),sngl(uvp(3*i-0))
c                 write(*,*)sngl(xn(i)),sngl(yn(i))
c                 pause
                endif
                zkar=zkasu
                if(err.lt.1d-6)zkar=zkasu*4.5d-1                    !移動後の節点位置において速度を計算し，
                uvp(3*i-2)=uvp(3*i-2)+(sf(3*i-2)-uvp(3*i-2))*zkar   !速度から求めた新たな移動後の節点位置
                uvp(3*i-1)=uvp(3*i-1)+(sf(3*i-1)-uvp(3*i-1))*zkar   !との差が小さくなった場合，流動ループ
                uvp(3*i)=uvp(3*i)+(sf(3*i)-uvp(3*i))*zkar           !を抜ける.
                xn(i)=xw(i)+(uvp0(3*i-2)+uvp(3*i-2))*dt/2d0
                yn(i)=yw(i)+(uvp0(3*i-1)+uvp(3*i-1))*dt/2d0
 290          continue
c             表示
c               call pltnt(nelf,new,xn,yn,idw,ind,ine)
               write(*,'(a3,i3,a3,i5,a5,e20.10)') 
     &                          'lk',lk,"kp",kkk,"err",sngl(err)
c           収束判定
            if(err.lt.1d-7) then
              open(10,file="err.res",access='append')
              write(10,*)"ti",sngl(ti),"  lk",lk,"  err",sngl(err)
              close(10)
              goto 300
	      endif
280         continue
            open(10,file="err.res",access='append')
            write(10,*)"ti",sngl(ti),"  lk",lk,"  err",sngl(err),
     &                 "   No convergion"
            close(10)
                write(*,*) ' !!!No convergion in lk loop/err=',sngl(err)
 300        continue
c
c           速度，圧力を格納した配列uvp(3*ind)に第1種の境界条件を代入する
c
            call bounduvp(nbd,kbd,bdv,uvp,ind) 
c
c           オブジェクトkcの子オブジェクトmcr(kc)を構成する節点の座標，速度，圧力を更新する
c
            do 310 n=1,ncrg
            do 310 i=nelrs(kcrg(n)),nelre(kcrg(n))
            do 310 j=1,3
            xx(ne(i,j))=xn(jnw(ne(i,j)))
            yy(ne(i,j))=yn(jnw(ne(i,j)))
            uu(ne(i,j))=uvp(jnw(ne(i,j))*3-2)
            vv(ne(i,j))=uvp(jnw(ne(i,j))*3-1)
310         pp(ne(i,j))=uvp(jnw(ne(i,j))*3)
230      continue			   !流動計算終了
c
c--------物質移動解析--------------------------------------------
c        2次要素を構成する節点の座標の更新
         do i=1,nel
	   do j=1,6
	     if(j.le.3)then
	       xx2(ne2(i,j))=xx(ne(i,j))
	       yy2(ne2(i,j))=yy(ne(i,j))
	     else
	       xx2(ne2(i,j))=(xx(ne(i,j-3))+xx(ne(i,js2(j-3))))/2d0
	       yy2(ne2(i,j))=(yy(ne(i,j-3))+yy(ne(i,js2(j-3))))/2d0
	     endif
	   enddo
	   enddo
c-----------------------------------------------------------------
c         オブジェクトkcにおけるx,y座標の最大値(xfmax,yfmax)
c         およびx=0d0におけるy座標の最大値yfasyを計算する
          xfmax=0d0
          yfmax=0d0
          yfasy=0d0
          do 330 m=1,mcr(1)
            do 330 n=1,ncr(1,m)
              k=kcr(1,m,n)
              do 330 j=1,nb(k)
              xfmax=dmax1(xx(ne(nelb(k,j),neb(k,j))),xfmax)
              yfmax=dmax1(yy(ne(nelb(k,j),neb(k,j))),yfmax)
              if(xx(ne(nelb(k,j),neb(k,j))).lt.1d-12)
     &           yfasy=dmax1(yy(ne(nelb(k,j),neb(k,j))),yfasy)
 330      continue
          xcl=xx(kcp(1))
          ucl=uu(kcp(1))
c*****************************************************************
c
c        熱移動解析
c
c*****************************************************************
         do 4000 m=1,mcr(1)    !mcr(kc):オブジェクトkcを構成する子オブジェクトの数
           npf=npc(1,m)        !npc(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)を構成する節点の数
           nbwf=nbwc(1,m)      !nbwc(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)のバンド幅
           ncrg=ncr(1,m)       !ncr(kc,mcr(kc)):オブジェクトkcの子オブジェクトmcr(kc)を構成する領域の数
           nelf=0              !nelf:オブジェクトkcの子オブジェクトmcr(kc)を構成する要素の数
           do 4100 i=1,ind
 4100      jnw(i)=jnt(1,m,i)   !jnt(kc,mcr(kc),np):オブジェクトkcを構成する最適化した節点の番号
           do 4200 n=1,ncrg
             kcrg(n)=kcr(1,m,n)    !kcr(kc,mcr(kc),ncr(kc,mcr(kc))
                                   !オブジェクトkcを構成する領域ncr(kc,mcr(kc)))の番号
             do 4300 i=nelrs(kcrg(n)),nelre(kcrg(n))  !nelrs(kr):領域krを構成する要素の最初の番号
                                                       !nelrs(kr):領域krを構成する要素の最後の番号
             nelf=nelf+1
             idw(nelf)=idx(i)
             irw(nelf)=irx(i)
             do 4300 j=1,3
               new(nelf,j)=jnw(ne(i,j))               !jnw(古)=新
               xw(new(nelf,j))=xx(ne(i,j))
               yw(new(nelf,j))=yy(ne(i,j))
               uvp0(new(nelf,j)*3-2)=uu(ne(i,j))
               uvp0(new(nelf,j)*3-1)=vv(ne(i,j))
               uvp0(new(nelf,j)*3)=pp(ne(i,j))
               jbw(jnw(ne(i,j)))=ne(i,j)              !jbw(新)=古
 4300        continue
 4200      continue
c
           nelt=0
           npt=0
           do i=1,npf
           jwt(i)=0
           enddo
c          npt,nelt,net(nelt,3)取得
           do i=1,nelf
             if(idw(i).eq.1)then      !****インデクス番号1の領域について温度場を計算する****
               nelt=nelt+1
               do j=1,3
               if(jwt(new(i,j)).eq.0)then   !まだ考慮されていない節点
                 npt=npt+1
                 net(nelt,j)=npt
                 jwt(new(i,j))=npt           !jwt(古)=新
                 jwtt(npt)=new(i,j)          !jwtt(新)=古
               else
                 net(nelt,j)=jwt(new(i,j))   !既に考慮されている節点
               endif
               enddo
             endif
           enddo
c          npt,nelt,net(nelt,3)の最適化
           call renum2(1,npt,net,nelt,jwtr,nbt)
           do i=1,npt   !i:新
           jwtrrr(jwtt(i))=jwtr(i)       !jwtrrr(古)=新_新
           jwtrr(jwtr(i))=jwtt(i)        !jwtrr(新_新)=古
           enddo
c
c          write(*,*)"npt,nelt,nbt",npt,nelt,nbt
c
           do i=1,npt
             do j=1,nbt*2
             sat(i,j)=0d0
             enddo
             sft(i)=0d0
           enddo
           write(*,*)"calc temperature"
           call temp(nel,nbt,ne,xx,yy,sat,sft,dt,tp0       !マトリックスを組み立てる
     &                ,pp,p0,tt0,idx,jwtrrr,ibw,ind,ine)
            call boundtmp(nel,nbt,ne,xx,yy,sat,sft,        !自由表面に第2種,蒸発潜熱の境界条件を与える
     &                      tp0,tt0,idx,jwtrrr,ibw,ind,ine)
            call boundtmp2(nel,npt,nbt,ne,sat,sft,idx
     &                      ,jwtrrr,ibw,ind,ine)            !平板との境界面に第1種の境界条件を与える
c            call choleski2(npt,nbt,sat,sft,ind,ibw)        !マトリックスを解く
            call gauss(npt,nbt,sat,sft,ind,ibw)
            do i=1,npt   !i:新_新
            tp(jwtrr(i))=sft(i)
            enddo
            call boundtmp3(nel,ne,tp,idx,ind,ine)          !平板との境界面に第1種の境界条件を与える
4000     continue
c*****************************************************************
c
c        物質移動解析
c
c*****************************************************************
       if(irec.eq.1)goto 2100      !irec=1の場合，濃度場の計算を行う
         time2=ti
         write(*,*)"recovery of concentration ti",sngl(ti)
         write(*,*)"ccc1,ccc2=",sngl(ccc1),sngl(ccc2)
         write(*,*)"-------------------------------"
         call getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &                 ip,jp,ine,ind)
         call renum2(2,npk,nek,nelk,jww,nbw2)
         do i=1,nelk
         do j=1,6
         xxk(jww(nek(i,j)))=xx2(ip(nek(i,j)))         !ip(新)=古
         yyk(jww(nek(i,j)))=yy2(ip(nek(i,j)))         !jp(古)=新
         ttk(jww(nek(i,j)))=tt0(ip(nek(i,j)))         !jww(新)=新_新
         jwwb(jww(nek(i,j)))=ip(nek(i,j))             !jwwb(新_新)=古
         nek(i,j)=jww(nek(i,j))                       !jww(jp(古))=新_新
         enddo
         enddo
c
	   do 1900 j=1,nbw2*2
	   do 1900 i=1,npk
1900     sa2(i,j)=0d0
         do 1950 i=1,npk
1950     sf2(i)=0d0
         dtt=time2-time1
		 if(dtt.lt.1d-10)then
		 write(*,*)"dtt=",dtt
		 pause
		 endif
	    call noudot(nelk,nbw2,nek,xxk,yyk,tp,
     &                 ne2,jww,jp,sa2,sf2,dtt,ttk,tt0,ind,ine,ibw)
c            call boundfsw(ne2,uu,vv,tt,px2,py2,ind,ine)    !自由表面に境界濃度を与える
c         call boundfsw2(nel,ne,xx,yy,ne2,                  !time1:前回濃度計算を行った時間
c     &                        tt,irx,time1,time2,ind,ine)  !time2:今回濃度計算を行う時間
c	    call btw2(ne2,nuv,kuv,uvb,ttk,jww,jp,ine,ind)
c         call btwvap(ne2,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)
c	   call bound(npk,nbw2,sa2,sf2,nuv,kuv,uvb,ind,ibw)
	   call choleski2(npk,nbw2,sa2,sf2,ind,ibw)
c       call gauss(npk,nbw2,sa2,sf2,ind,ibw)
	   do 2000 i=1,npk
2000   ttk(i)=sf2(i)
       do i=1,npk
	   tt(jwwb(i))=ttk(i)
	   enddo
	   call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
	   if(ccc2.ge.ccc1)irec=1
c
        open(33,file="callog.res",access='append')
        write(33,*)"dt=",sngl(time2-time1)
        close(33)
c
	   time1=time2
c
	   open(33,file="callog.res",access='append')
	   write(33,*)"ti=",sngl(ti),"   in time loop dt=",time2-time1
	   close(33)
	   call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,tt,   !
     &                  ti,rfile,ine,ind)
	   call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                            ti,rfile,ine,ind)
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !初期濃度からの増分から求めた物質量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	   !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)       !ccc:増加量, ccct:吸収量
	  write(4,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	   !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)      !ccc/vvv:増加した濃度の平均値，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c---------------------------------------------------------------

c
2100   continue
c---------------------------------------------------------
c
c       call pltv(np,ne,xx,yy,uu,vv,ine,ind) 
cc        write(*,*)"pltnt started in main"     !
cc        call pltnt(nel,ne,xx,yy,idx,ind,ine) 
cc        write(*,*)"pltnt ended in main"       !
c
c**********************************************************************
c
c     dataの保存
c
c**********************************************************************
        if(ti+1d-6*wdt1.ge.wtime1) then
	    time2=ti
	    if(time2.lt.time1+1d-10)goto 3100    !既に濃度場が計算されている場合，濃度場を計算しない
	    write(*,*)"save avs data"
          call getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &                  ip,jp,ine,ind)
	    call renum2(2,npk,nek,nelk,jww,nbw2)
	    do i=1,nelk
	    do j=1,6
	    xxk(jww(nek(i,j)))=xx2(ip(nek(i,j)))
	    yyk(jww(nek(i,j)))=yy2(ip(nek(i,j)))
	    ttk(jww(nek(i,j)))=tt0(ip(nek(i,j)))
	    jwwb(jww(nek(i,j)))=ip(nek(i,j))
	    nek(i,j)=jww(nek(i,j))
	    enddo
	    enddo
c
	    do 2900 j=1,nbw2*2
	    do 2900 i=1,npk
2900      sa2(i,j)=0d0
          do 2950 i=1,npk
2950      sf2(i)=0d0
          dtt=time2-time1
		if(dtt.lt.1d-10)then
		write(*,*)"dtt=",dtt
		pause
		endif
	    call noudot(nelk,nbw2,nek,xxk,yyk,tp,      !jww(jp(古))=新
     &		          ne2,jww,jp,sa2,sf2,dtt,ttk,tt0,ind,ine,ibw)	    	   
c	      write(*,*)"sub boundfsw"
c            call boundfsw(ne2,uu,vv,tt,px2,py2,ind,ine)    !自由表面に境界濃度を与える
c         call boundfsw2(nel,ne,xx,yy,ne2,                  !time1:前回濃度計算を行った時間
c     &                        tt,irx,time1,time2,ind,ine)  !time2:今回濃度計算を行う時間
c	     call btw2(ne2,nuv,kuv,uvb,ttk,jww,jp,ine,ind)
c          write(*,*)"sub btwvap"
c          call btwvap(ne2,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)
c	    call bound(npk,nbw2,sa2,sf2,nuv,kuv,uvb,ind,ibw)
	    call choleski2(npk,nbw2,sa2,sf2,ind,ibw)
c          call gauss(npk,nbw2,sa2,sf2,ind,ibw)
	    do 3000 i=1,npk
3000      ttk(i)=sf2(i)
          do i=1,npk
	    tt(jwwb(i))=ttk(i)
	    enddo
	    call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
	    if(ccc2.ge.ccc1)irec=1
c
          open(33,file="callog.res",access='append')
          write(33,*)"dt=",sngl(time2-time1)
          close(33)
c
	    time1=time2
3100	    continue
          open(1,file="noudo.res")
	    rewind(1)
		write(1,*)"np2=",np2
		do i=1,np2
		write(1,*)xx2(i),",",yy2(i),",",tt(i)
		enddo
		write(1,*)"nuv(2)=",nuv(2)
		do i=1,nuv(2)
		write(1,*)i,",",kuv(i,2),",",sngl(xxk(kuv(i,2))),",",
     &                 sngl(yyk(kuv(i,2))),",",sngl(uvb(i,2))
		enddo
		close(1) 
c
          wtime1=ti+wdt1
          if(iavs.eq.2 .or. iavs.eq.5)then 
          write(*,*)"avs2dm started in main"  !
c---------速度保存------------------------------------------
          call avs2dm(nel,np,ne,idx,xx,yy,xn,uu,vv,yn,pp,
     &                  ti,rfile,ine,ind)
c---------速度保存終了--------------------------------------
c---------温度保存------------------------------------------
          call avs2dmtmp(nel,np,ne,idx,xx,yy,tp,tt,
     &                  ti,rfile,ine,ind)
c
	    write(ltp,'(f8.1)') ti*1d4
          do i=1,8
            if(ltp(i:i).eq.' ') ltp(i:i)='0'
          enddo
          open(1,file='tempr-'//ltp//'.res')
          rewind(1)
          do i=1,np
          write(1,*)sngl(xx(i)),",",sngl(yy(i)),",",sngl(tp(i))
          enddo
          close(1)
c---------温度保存終了--------------------------------------
c---------濃度保存保存--------------------------------------
	    open(33,file="callog.res",access='append')
	    write(33,*)"ti=",sngl(ti),"   in logging"
	    close(33)
	    call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,tt,ti,rfile,ine,ind)
	    call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                            ti,rfile,ine,ind)
c---------濃度保存終了--------------------------------------
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !初期濃度からの増分から求めた物質量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	       !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)            !ccc:増加量, ccct:吸収量
	  write(4,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	       !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)          !ccc/vvv:増加した濃度の平均値，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c---------------------------------------------------------------
           endif
c
          if(iavs.eq.3 .or. iavs.eq.5)then
          write(*,*)"avs3dm started in main"  !
          call avs3dm(nel,np,ne,idx,xx,yy,xn,uu,vv,yn,pp,
     &                  ti,rfile,ine,ind)
	    call avs3dmmas(nel,np2,ne2,idx,xx2,yy2,xn,tt,
     &				  ti,rfile,ine,ind)
          write(*,*)"avs3dm ended in main"    !
          endif
c-----------------中心軸における界面位置を求める-------------------
          ymax=-9999			 
          do i=1,nb(2)
	    if(dabs(xx(ne(nelb(2,i),neb(2,i)))).lt.1d-12 .and. 
     &                             idx(nelb(2,i)).eq.1)
     &               ymax=dmax1(ymax,yy(ne(nelb(2,i),neb(2,i))))
	     enddo
	     open(1,file='t-z.res',access='append')
	     write(1,*)sngl(ti),ymax
	     close(1)
c-------------------------------------------------------------------
c
        endif
c------------------------------------------------------------------
c       save res data
c------------------------------------------------------------------
        if(ti+1d-6*wdt2.ge.wtime2) then       !save res data
          write(*,*)"save res"
          wtime2=ti+wdt2
          open(2,file=rfile2,access='append')
          write(2,'(20E12.4)') ti,xfmax,xcl,yfmax,yfasy,ucl,
     &                          real(icl(1)),ancl(1),st(1)
          close(2)
        endif
c------------------------------------------------------------------
c       save tmp data
c------------------------------------------------------------------
        if(ti+1d-6*wdt3.ge.wtime3) then       !save tmp data
          write(*,*)"save tmp"
          wtime3=ti+wdt3
c
          write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
          do i=1,8
            if(lpp(i:i).eq.' ') lpp(i:i)='0'
	    enddo
          call addchr(rfile,lpp//fd3,rfile3)
c
          open(4,file=rfile3,status='unknown',access='sequential',
     $            form='unformatted')
          rewind(4)
          write(4) ti,hmax,icoal,wtime1-wdt1,wtime2-wdt2,wtime3-wdt3
          write(4) np,(xx(i),yy(i),i=1,np)
          write(4) np,(uu(i),vv(i),pp(i),i=1,np)
          write(4) nel,((ne(i,j),j=1,3),idx(i),irx(i),ar0(i),i=1,nel)
          write(4) nr,(nb(i),ndx(i),(nelb(i,j),neb(i,j),j=-1,nb(i)+2),
     &          (ncom(i,j),ncop(i,j),j=1,nb(i)+1),i=1,nr)
          write(4) (ncon(i),(kcon(i,j),j=1,nr),i=1,nr)
          write(4) nset,(kset(i),i=1,nr)
          write(4) (nelrs(i),nelre(i),i=1,nr)
          write(4) msw,(mcr(i),(ncr(i,j),npc(i,j),nbwc(i,j),j=1,mcr(i))
     &          ,i=1,msw)
          write(4) (((kcr(i,j,k),k=1,ncr(i,j)),j=1,mcr(i)),i=1,msw)
          write(4) (((jnt(i,j,k),k=1,np),j=1,mcr(i)),i=1,msw)
          write(4) ncp,(kcp(i),icp(i),jcp(i),st(i),ancl(i),
     &                  icl(i),icla(i),staa(i),strr(i),i=1,ncp)
c-------------------物質移動,2次要素のデータ------------------------
	    write(4) np2,(xx2(i),yy2(i),i=1,np2)
	    write(4) np2,(tt(i),i=1,np2)
	    write(4) nel,((ne2(i,j),j=1,6),i=1,nel)
          write(4) nr,(nbw(i),(nelbw(i,j),nebw(i,j),j=-1,nbw(i)+2),
     &             (ncomw(i,j),ncopw(i,j),j=1,nbw(i)+2),i=1,nr)
	    write(4) tt00,ttw,csta,time1
	    write(4) cdv,area0,aaa,ccc1,nbw2
c-------------------------------------------------------------------
c-----------------熱移動-------------------------------------------
          write(4) np,(tp(i),i=1,np)
c------------------------------------------------------------------
c-----------------計算の時間間隔dt---------------------------------
          write(4) imesh,kmesh,dt
c------------------------------------------------------------------
          close(4)
        endif
c
        write(*,*)"coalesce"
        call coalesce(nel,ne,xx,yy,irx,idx,icoal,ind,ine) 
        if(icoal.eq.0)
     &          call checkelemnt(nel,ne,xx,yy,idx,ar0,hmax,ind,ine)
c
c---------------------mesh再発生----------------------------------
c
        if (hmax.gt.ratio .or. icoal.eq.1) then
		write(*,*)"regenerate mesh"
c
c-------物質移動,2次要素の計算(dttt,ccc1,tt(np2)----------------------------
          time2=ti
	    if(inin.eq.0)then
	    ccc=csta
	    inin=1
	    else
	    ccc=ccc1
	    endif
c
c
c--------------------------------------
          open(4,file=rfile3,status='unknown',access='sequential',
     $            form='unformatted')
          rewind(4)
          write(4) ti,hmax,icoal,wtime1-wdt1,wtime2-wdt2,wtime3-wdt3
          write(4) np,(xx(i),yy(i),i=1,np)
          write(4) np,(uu(i),vv(i),pp(i),i=1,np)
          write(4) nel,((ne(i,j),j=1,3),idx(i),irx(i),ar0(i),i=1,nel)
          write(4) nr,(nb(i),ndx(i),(nelb(i,j),neb(i,j),j=-1,nb(i)+2),
     &          (ncom(i,j),ncop(i,j),j=1,nb(i)+1),i=1,nr)
          write(4) (ncon(i),(kcon(i,j),j=1,nr),i=1,nr)
          write(4) nset,(kset(i),i=1,nr)
          write(4) (nelrs(i),nelre(i),i=1,nr)
          write(4) msw,(mcr(i),(ncr(i,j),npc(i,j),nbwc(i,j),j=1,mcr(i))
     &          ,i=1,msw)
          write(4) (((kcr(i,j,k),k=1,ncr(i,j)),j=1,mcr(i)),i=1,msw)
          write(4) (((jnt(i,j,k),k=1,np),j=1,mcr(i)),i=1,msw)
          write(4) ncp,(kcp(i),icp(i),jcp(i),st(i),ancl(i),
     &                  icl(i),icla(i),staa(i),strr(i),i=1,ncp)
c-------------------物質移動,2次要素のデータ------------------------
	    write(4) np2,(xx2(i),yy2(i),i=1,np2)
	    write(4) np2,(tt(i),i=1,np2)
	    write(4) nel,((ne2(i,j),j=1,6),i=1,nel)
          write(4) nr,(nbw(i),(nelbw(i,j),nebw(i,j),j=-1,nbw(i)+2),
     &             (ncomw(i,j),ncopw(i,j),j=1,nbw(i)+2),i=1,nr)
	    write(4) tt00,ttw,csta,time1
	    write(4) cdv,area0,aaa,ccc1,nbw2
c-------------------------------------------------------------------
c-----------------熱移動-------------------------------------------
          write(4) np,(tp(i),i=1,np)
c------------------------------------------------------------------
c-----------------計算の時間間隔dt---------------------------------
          write(4) imesh,kmesh,dt
c------------------------------------------------------------------
          close(4)
c
c
	    call getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &               ip,jp,ine,ind)
	    call renum2(2,npk,nek,nelk,jww,nbw2)
	    do i=1,nelk
	      do j=1,6
	      xxk(jww(nek(i,j)))=xx2(ip(nek(i,j)))
	      yyk(jww(nek(i,j)))=yy2(ip(nek(i,j)))
	      ttk(jww(nek(i,j)))=tt0(ip(nek(i,j)))
	      jwwb(jww(nek(i,j)))=ip(nek(i,j))
	      nek(i,j)=jww(nek(i,j))
	      enddo
	    enddo
c
	    do 2200 j=1,nbw2*2
		do 2200 i=1,npk
2200      sa2(i,j)=0d0
          do 2250 i=1,npk
2250      sf2(i)=0d0
          dttt=time2-time1
c
c
c-------------------------------------------------------------------
c----------------既にこの時間の計算が終了している場合--------------
c
		if(dttt.lt.1d-10)then
          write(*,*)"already calcrated 
     &                          at this time in regenerating mesh"
	    call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc1,ind,ine)
	    open(33,file="callog.res",access='append')
          write(33,*)"ti=",sngl(ti),
     &                 "   finished counting at this time before grd"
          close(33)
          write(*,*)"dttt=",dttt
	    goto 9000		  !切り直し前の濃度計算を行わない
		endif
c
c-------------------------------------------------------------------
c-------------------------------------------------------------------
c
c
		write(*,*)"sub noudot"
	    call noudot(nelk,nbw2,nek,xxk,yyk,tp,
     &                 ne2,jww,jp,sa2,sf2,dttt,ttk,tt0,ind,ine,ibw)	    	   
c	    write(*,*)"sub btw"
c          call boundfsw(ne2,uu,vv,tt,px2,py2,ind,ine)    !自由表面に境界濃度を与える
c         call boundfsw2(nel,ne,xx,yy,ne2,                  !time1:前回濃度計算を行った時間
c     &                        tt,irx,time1,time2,ind,ine)  !time2:前回濃度計算を行った時間
c	     call btw2(ne2,nuv,kuv,uvb,ttk,jww,jp,ine,ind)
c          write(*,*)"sub btwvap"
c          call btwvap(ne2,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)
c		write(*,*)"sub bound"
c	    call bound(npk,nbw2,sa2,sf2,nuv,kuv,uvb,ind,ibw)
		write(*,*)"sub choleski2 np2=",np2,"nbw2=",nbw2
	    call choleski2(npk,nbw2,sa2,sf2,ind,ibw)
c		 call gauss(npk,nbw2,sa2,sf2,ind,ibw)
		write(*,*)"input data into tt(np2)"
	    do 2300 i=1,npk
2300      ttk(i)=sf2(i)
          do i=1,npk
	    tt(jwwb(i))=ttk(i)
	    enddo
          write(*,*)"calc ccc1 in regenerating mesh"
	    call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc1,ind,ine)
c
	    open(33,file="callog.res",access='append')
          write(33,*)"ti=",sngl(ti),
     &                 "   in regenerating mesh before grd"
          close(33)
          call avs2dm(nel,np,ne,idx,xx,yy,xn,uu,vv,yn,pp,
     &                  ti,rfile,ine,ind)
          call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,tt,ti,rfile,ine,ind)   !
         call avs2dmmasr(nel,np2,ne2,idx,xx2,yy2,xn,tt,ti,rfile,ine,ind)   !
  	    call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                            ti,rfile,ine,ind)
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !初期濃度からの増分から求めた物質量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	       !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)            !ccc:増加量, ccct:吸収量
	  write(4,'(20E15.6)')sngl(ti),sngl(dpaa),sngl(dpvv),	       !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)          !ccc/vvv:増加した濃度の平均値，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c
c---------------------------------------------------------------
c---------------------------------------------------------------
c
c-------熱移動(データ保存)-------------------------
c
          write(*,*)"avs2dmtmp started"
          call avs2dmtmp(nel,np,ne,idx,xx,yy,tp,tt,
     &                  ti,rfile,ine,ind)
c
c-------熱移動終了-----------------------------------------------
c
9000    continue
c
          neal=nel
          write(*,*)"grd started in the end of main"    !
          call grd(nel,np,ne,xx,yy,idx,irx,jnt,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          new,xw,yw,idw,irw,jnw,
     &          nelrs,nelre,npc,nbwc,mcr,ncr,kcr,
     &          ar0,ind,ine,ibw,ica)
c-----熱移動(iyy=0の場合における温度補間に用いる)----------------
          nelt=nel
          npt=np
          do i=1,nelt
          do j=1,3
            net(i,j)=ne(i,j)
            xxt(net(i,j))=xx(net(i,j))
            yyt(net(i,j))=yy(net(i,j))
          enddo
          enddo
c-------熱移動終了-----------------------------------------------
c-------速度補間-------------------------------------------------
          write(*,*)"grdfcn started in the end of main"  !
          call grdfcn(nel,np,ne,xx,yy,idx,irx,jnw,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          uu,pp,vv,u0,v0,p0,ind,ine)
c-------速度補間終了---------------------------------------------
c
c---------mesh切り直し後の時間間隔を決定--------------------------------
          imesh=1
          kmesh=1
c-----------------------------------------------------------------------
c
c-------熱移動(温度補間)-----------------------------------------
c
          write(*,*)"grdfnt started inthe end of main"
          call grdfct(nelt,npt,net,xxt,yyt,tp,
     &                    neln,npn,nen,xn,yn,tpn,irx,ind,ine)
          do i=1,nel
          do j=1,3
            if(idx(i).eq.-1 .and. xx(ne(i,j)).gt.-1d-15      !平板に温度代入(R≧0,Z≦0)
     &                      .and. yy(ne(i,j)).lt.1d-15) 
     &                     tp(ne(i,j))=tplate(1)/tplate(1)
          enddo
          enddo
c
c-------熱移動終了-----------------------------------------------
c
c---------物質移動,2次要素の計算(ccc2,tt(np2))----------------------------
c
c          call generate2el(np,ne,nel,xx,yy,ine,ind,
c     $                 np22,ne22,xx22,yy22,nbw22)
          call genet2el(np,nel,ne,xx,yy,
     &                np22,ne22,xx22,yy22,nbw22)
	    call sufelew
	    call grdfc5(neal,np2,nbw2,ne2,xx2,yy2,tt,
     &                nel,np22,nbw22,ne22,xx22,yy22,tt22,irx,ind,ine)
	    nel=neal
c	    time1=time2
          irec=1
	    write(*,*)"calc ccc2 in regenerating mesh"
	    call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
		if(ccc2.lt.ccc1)irec=0
c
          open(1,file="ne22.res")
	    rewind(1)
	    write(1,*)"nel=",nel
	    do i=1,nel
	    do j=1,6
	    write(1,*)i,",",ne2(i,j),",",xx2(ne2(i,j)),",",yy2(ne2(i,j))
	    enddo
	    enddo
	    write(1,*)"-----"
	    write(1,*)"np2=",np2
	    do i=1,np2
	    write(1,*)i,",",xx2(i),",",yy2(i)
	    enddo
	    close(1)
c
          open(33,file="callog.res",access='append')
          write(33,*)"dt=",sngl(time2-time1)
          close(33)
c
	    open(33,file="callog.res",access='append')
          write(33,*)"ti=",sngl(ti),
     &                 "   in regenerating mesh before recov"
          close(33)
          call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,
     &                   tt,ti+4d1,rfile,ine,ind)   !
         call avs2dmmasr(nel,np2,ne2,idx,xx2,yy2,xn,
     &                   tt,ti+4d1,rfile,ine,ind)   !
  	    call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                            ti+4d1,rfile,ine,ind)
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !初期濃度からの増分から求めた物質量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti+4d1),sngl(dpaa),sngl(dpvv),    !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)            !ccc:増加量, ccct:吸収量
	  write(4,'(20E15.6)')sngl(ti+4d1),sngl(dpaa),sngl(dpvv),	   !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)          !ccc/vvv:増加した濃度の平均値，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c
c---------recovery of concentration(dtt1,ccc2,tt(np2))--------------------------
          if(ccc2.lt.ccc1)then
	      write(*,*)"recovery of concentartion"
		  dtt1=dt
		  ccc3=ccc2
2550        continue
	      call getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &                       ip,jp,ine,ind)
	      call renum2(2,npk,nek,nelk,jww,nbw2)
            do i=1,nelk
	        do j=1,6
	          xxk(jww(nek(i,j)))=xx2(ip(nek(i,j)))  !jp(古)=新
	          yyk(jww(nek(i,j)))=yy2(ip(nek(i,j)))  !ip(新)=古
	          ttk(jww(nek(i,j)))=tt(ip(nek(i,j)))   !jww(新)=新_新
	          jwwb(jww(nek(i,j)))=ip(nek(i,j))      !jwwb(新_新)=古
	          nek(i,j)=jww(nek(i,j))
	        enddo
	      enddo
c
            do 2600 j=1,nbw2*2
  		  do 2600 i=1,npk
2600        sa2(i,j)=0d0
            do 2650 i=1,npk
2650        sf2(i)=0d0
            write(*,*)"sub noudot"
 	       call noudot(nelk,nbw2,nek,xxk,yyk,tp,
     &                 ne2,jww,jp,sa2,sf2,dtt1,ttk,tt0,ind,ine,ibw)	    	   
c	       write(*,*)"sub btw"
c            call boundfsw(ne2,uu,vv,tt,px2,py2,ind,ine)    !自由表面に境界濃度を与える
c             call boundfsw2(nel,ne,xx,yy,ne2,              !time1:前回濃度計算を行った時間
c     &                        tt,irx,time1,time2,ind,ine)  !time2:前回濃度計算を行った時間
c	        call btw2(ne2,nuv,kuv,uvb,ttk,jww,jp,ine,ind)
c             write(*,*)"sub btwvap"
c             call btwvap(ne2,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)
c	       write(*,*)"sub bound"
c	       call bound(npk,nbw2,sa2,sf2,nuv,kuv,uvb,ind,ibw)
             call choleski2(npk,nbw2,sa2,sf2,ind,ibw)
c            call gauss(npk,nbw2,sa2,sf2,ind,ibw)
	        do 2700 i=1,npk
2700        ttk(i)=sf2(i)
            do i=1,npk
	        tt(jwwb(i))=ttk(i)
	        enddo
            write(*,*)"calc recov ccc2"
	        call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,ccc2,ind,ine)
  	  	    if(ccc2.ge.ccc1 .or. ccc2.lt.ccc3)goto 2750		  !ccc1:mesh切り直し前の吸収量
  		    dtt1=(ccc1-ccc3)/(ccc2-ccc3)*dtt1/2d0
		    if(dtt1.lt.dt)dtt1=dt
c
            open(33,file="callog.res",access='append')
            write(33,*)"dt=",sngl(dtt1),
     &                      "ccc1-3",sngl(ccc1),sngl(ccc2),sngl(ccc3)
            close(33)
c
		    write(*,*)"ccc1,ccc2,ccc3=",
     &                              sngl(ccc1),sngl(ccc2),sngl(ccc3)
	        write(*,*)"dt,dtt1=",sngl(dt),sngl(dtt1)
	        write(*,*)"-----------------"
	        ccc3=ccc2
	        goto 2550
	      endif
2750      continue
	    time1=time2	    
	      irec=1
c----------物質移動,2次要素の計算終了-----------------------------------
          open(33,file="callog.res",access='append')
          write(33,*)"ti=",sngl(ti),"   in regenerating mesh"
	    close(33)
          call avs2dmmas(nel,np2,ne2,idx,xx2,yy2,xn,
     &                            tt,ti+5d1,rfile,ine,ind)   !
          call avs2dmmasr(nel,np2,ne2,idx,xx2,yy2,xn,
     &                   tt,ti+5d1,rfile,ine,ind)   !
 	    call save2el(nel,np2,ne2,idx,xx2,yy2,tt,
     &                            ti+5d1,rfile,ine,ind)
c---------------吸収量の保存-----------------------------
        write(*,*)"save amount"    !
        call sufar(ne,xx,yy,dpaa,ind,ine)
	  call vdrop(nel,ne,xx,yy,idx,dpvv,ind,ine)
	  call calw2(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcc,ind,ine)  !初期濃度からの増分から求めた物質量
	  call calw3(nel,ne2,xx2,yy2,tt,nelrs,nelre,dpcct,ind,ine) !液滴の全物質量
c	  call riron(ti,pe,dltm,flth,amth)
	  open(3,file='ampart.res',access='append')
	  open(4,file='amtotal.res',access='append')
	  write(3,'(20E15.6)')sngl(ti+5d1),sngl(dpaa),sngl(dpvv),    !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcc),sngl(dpcc/dpvv)            !ccc:増加量, ccct:吸収量
	  write(4,'(20E15.6)')sngl(ti+5d1),sngl(dpaa),sngl(dpvv),    !aaa:液滴表面積, vvv:液滴体積
     &                     sngl(dpcct),sngl(dpcct/dpvv)          !ccc/vvv:増加した濃度の平均値，ccct/vvv:平均濃度
	  close(4)
	  close(3)
c---------------------------------------------------------------
c
c               call pltv(np,ne,xx,yy,uu,vv,ine,ind)
c               call pltnt(nel,ne,xx,yy,idx,ind,ine)
c               pause
          write(*,*)"grdfcn ended in the end of main"    !
 9999   continue
c
        endif
c
c--------------------mesh再発生の処理終了----------------------------
c
        write(*,*) '*ti,hmax,xcl',sngl(ti),sngl(hmax),sngl(xcl)
        write(*,*) '  icl,st,ag,ucl='
     &  ,icl(1),sngl(st(1)),sngl(ancl(1)),sngl(ucl)
c       write(*,*) 'icl,icp,jcp=',icl(1),icp(1),jcp(1)
c
        do 610 i=1,np
        u0(i)=uu(i)
        v0(i)=vv(i)
        p0(i)=pp(i)
        x0(i)=xx(i)
 610    y0(i)=yy(i)
c-------物質移動--------------------------------
        do 2800 i=1,np2
2800    tt0(i)=tt(i)
c-------物質移動終了----------------------------
c-------熱移動----------------------------------
        do i=1,np
        tp0(i)=tp(i)
        enddo
c-------熱移動終了------------------------------
c
        open(1,file="uvp.res")
	  rewind(1)
	  do i=1,np
	  write(1,*)sngl(uu(i)),sngl(vv(i)),sngl(pp(i))
	  enddo
	  close(1)
c
        if(ti.ge.timax) goto 1100
 200    continue
c
 1100 call clos
      write(*,*) 'ti .ge. timax'
      stop
      end
c
        subroutine data(nel,np,ne,xx,yy,idx,ind,ine,ibw,dltm)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine)
        dimension zz(6)
        common /doper/ vi,di,sta,str
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        common /dbound/ kfbnd(-9:9,-9:9)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /dconhy/ stsa(21),zksa(21),zmsa(21),       !%yama(6.24)
     &          stsr(21),zksr(21),zmsr(21)
        common /thininf/ dthin(-9:9,-9:9),athin(-9:9,-9:9),
     &                  mthin(-9:9,-9:9)
	  common /mass/ dmas(9)
        character*50 gfile,gfile1,gfile2,pfile
        character*2 fd1,fd2,f(3)*1,ff*1
        data fd1,fd2/'.1','.2'/, f/'u','v','p'/
c
        write(*,*) 'Characteristic velosity and length (vi,ri)'
c        read(5,*) vi,di
        read(5,*) vi,di,dltm     !
        write(*,*) 'grid file'
        read(5,'(a30)') gfile
c       write(*,*) 'propety data'
c       read(5,'(a30)') pfile
c
        call addchr(gfile,fd1,gfile1)
        call addchr(gfile,fd2,gfile2)
c
        open(1,file=gfile1)
        rewind(1)
        read(1,*) np
        do 10 i=1,np
 10     read(1,*) j,k,xx(k),yy(k)
        close(1)
c
        open(1,file=gfile2)
        rewind(1)
        read(1,*) nel,np,nbwt
        do 20 i=1,nel
 20     read(1,*) k,(ne(k,j),j=1,3),zz,idx(k)
        close(1)
c
        write(*,*) '             ------- grid data -------'
        write(*,*) 'element=',nel,'  nodes=',np,'  band width=',nbwt
        if(np*3.gt.ind .or. nel.gt.ine .or. nbwt*6.gt.ibw) then
                write(*,*) 'ine,ind,ibw=',ine,ind,ibw
                stop
        endif
        call arrangegrid(nel,ne,idx,ind,ine)
        write(*,*)"ne(",nel,",1)=",ne(nel,1),"ne(",nel,",2)=",ne(nel,2),
     &        "ne(",nel,",3)=",ne(nel,3)
        write(*,*)"after arranging grid nel=",nel
c
        write(*,*)"input start"   !
        msw=0
        read(5,*) n
c        do 100 i=1,n
c 100    read(5,*) k,dens(k),visc(k),gy(k),svc(k)
        do 100 i=1,n									 !
 100    read(5,*) k,dens(k),visc(k),gy(k),svc(k),dmas(k)
        write(*,*)"finish 1st n=",n
        read(5,*) n
        do 110 i=1,n
        read(5,*) j,k,sten(j,k)
 110    sten(k,j)=sten(j,k)
        write(*,*)"end 2st n=",n
        read(5,*) n
        do 120 i=1,n
        read(5,*) k,mswi(k),mfin(k),etas(k)
 120    msw=max(mswi(k),msw)
        write(*,*)"end 3st n=",n
        read(5,*) n
        do 130 i=1,n
 130    read(5,*) j,k,kfbnd(j,k)
        write(*,*)"end 4st n=",n
        read(5,*) n
        do 140 i=1,n
        read(5,*) j,k,l
        icoml(j,k,l)=i
        read(5,*) stsa(i),zksa(i),zmsa(i),stsr(i),zksr(i),zmsr(i),v0
        zksa(i)=zksa(i)/v0
        zksr(i)=zksr(i)/v0
 140    continue
        write(*,*)"end 5st n=",n
        read(5,*) n
        if(n.eq.0)go to 161
        do 160 i=1,n
 160    read(5,*) j,k,mthin(j,k),dthin(j,k),athin(j,k)
 161    write(*,*)"input fnished"     !
c        close(5)
      return
      end
c
        subroutine initialc(nelt,ne,xx,yy,idx,uu,vv,pp,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),idx(0:ine),uu(ind),vv(ind),pp(ind)
        dimension xx(ind),yy(ind)
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        common /doper/ vi,di,sta,str
        do 20 i=1,nelt
        do 20 j=1,3
          k=ne(i,j)
          uu(k)=0d0
          vv(k)=0d0
          pp(k)=0d0
 20     continue
        do 40 i=1,nelt
          if(idx(i).eq.1) then
            do 50 j=1,3
            k=ne(i,j)
c            vv(k)=-vi     !
            pp(k)=4d0*sten(0,1)/di
 50         continue
          elseif(idx(i).eq.2) then
            do 60 j=1,3
            k=ne(i,j)
            vv(k)=0d0
            pp(k)=dens(2)*gy(2)*yy(k)
 60         continue
          endif
 40     continue

        return
        end
c
        function fconhy(im,icl,ucl)
        implicit double precision (a-h,o-z)
        common /dconhy/ stsa(21),zksa(21),zmsa(21),       !%yama(6.24)
     &          stsr(21),zksr(21),zmsr(21)
        if(icl.gt.0) then
          if(zksa(im).lt.1d10) then
            u=dmax1(0d0,ucl)
            fconhy=(u/zksa(im))**(1d0/zmsa(im))+stsa(im)
          else
            fconhy=stsa(im)
          endif     
        elseif(icl.lt.0) then
          if(zksa(im).lt.1d10) then
            u=dmin1(0d0,ucl)
            fconhy=(u/zksr(im))**(1d0/zmsr(im))+stsr(im)
          else
            fconhy=stsr(im)
          endif
        endif
        fconhy=dmax1(0.1d0,fconhy)
        fconhy=dmin1(179.9d0,fconhy)
        return
        end
c
        subroutine arrangegrid(nel,ne,idx,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,3),nec(ine,3),idx(0:ine)
        dimension js2(3)
        data js2/2,3,1/
c
        do 100 i=1,nel
        do 100 j=1,3
 100    nec(i,j)=0
c
        do 120 i=1,nel
          do 130 j=1,nel
            if(j.eq.i)goto 130
            do 140 k1=1,3
              if(nec(j,k1).ne.0)goto 140
                k2=js2(k1)
                do 150 jc=1,3
 150              if((ne(i,jc).eq.ne(j,k2)).and.
     &             (ne(i,js2(jc)).eq.ne(j,k1)))goto 190
 140        continue
            goto 130
 190        nec(i,jc)=j
            nec(j,k1)=i
            if(nec(i,1).ne.0.and.nec(i,2).ne.0.and.ne(i,3).ne.0)
     &          goto 120
 130      continue
 120    continue
c
        do 10 i=1,nel
        if(idx(i).ne.1)goto 10
        do 20 j=1,3
          if(idx(nec(i,j)).eq.0)then
          do 30 k=1,nel
          if(idx(k).ne.1)goto 30
          do 40 l=1,3
            if(ne(i,j).eq.ne(k,js2(l)).and.
     &         idx(nec(k,js2(l))).eq.-1)goto 50
 40       continue
 30       continue
          endif
 20     continue
 10     continue
        write(*,*)"didn't find nodes between droplets"
        pause
	  return
 50     nel=nel+1
        ne(nel,1)=ne(i,j)
        ne(nel,2)=ne(k,l)
        ne(nel,3)=ne(i,js2(j))
        idx(nel)=1
c
        return
        end
