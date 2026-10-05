       include "head.for"
       dimension ne(ine,3),net(ine,3),new(ine,3)
       dimension idx(0:ine),idt(0:ine),idw(0:ine)
       dimension xx(ind),yy(ind)
       dimension xt(ind),yt(ind)
       dimension xw(ind),yw(ind)
       dimension xn(ind),yn(ind)
c
c      物性値
       dimension vre(-5:ica),vma(-5:ica),vfr(-5:ica)
       dimension vwe(-5:ica,-5:ica)    !※※※※※※※
c
c       境界
        dimension nc(ica),nelc(ica,-1:ibd)
        dimension ndx(0:ica)
        dimension neled(0:ica)
        dimension mr(ica)
        dimension nb(ica)
        dimension nelb(ica,-1:ibd)
        dimension neb(ica,-1:ibd)
        dimension ncom(ica,-1:ibd)
c
        dimension mb(ica)
        dimension melb(ica,-1:ibd)
        dimension meb(ica,-1:ibd)
        dimension mcom(ica,-1:ibd)
c
c       バンド幅
        dimension ip(ind),jp(ind)
        dimension iq(ine),jq(ine)
        dimension jww(ind),jvv(ind)
c
c       物理量
        dimension uu(ind),vv(ind),pp(ind)
        dimension ut(ind),vt(ind),pt(ind)
        dimension uw(ind),vw(ind),pw(ind)
c
c       マトリックス
        dimension uvp(ind*3),uvp0(ind*3)
        dimension mbd(ibd),vbd(ibd)        
        dimension sa(ind,ibw),sf(ind),dsf(ind)
c
        character*40 tfile
        character*12 lpp
c
c       収束係数
        data zkasu/1.0/
c
c     ファイルを読み込む
       call data(iyy,dt,wdt1,wdt2,tmax,vi,ri,vre,vma,vfr,vwe
     &           ,np,nele,ne,idx,xx,yy)
c
c      初期化
       if(iyy.eq.0)then
         ti=0d0
         tsav=wdt1     !保存時間
         ttmp=wdt2     !保存時間
         do kp=1,np
           uu(kp)=0d0
           vv(kp)=0d0
           pp(kp)=0d0
         enddo
       else
         open(1,file="tmp.inp")
         rewind(1)
         read(1,*)ti,tsav,ttmp
         read(1,*)np,nele
         read(1,*)(ne(i,1),ne(i,2),ne(i,3),i=1,nele)
         read(1,*)(idx(i),i=0,nele)
         read(1,*)(xx(i),yy(i),i=1,np)
         read(1,*)(uu(i),vv(i),pp(i),i=1,np)
         close(1)
       endif
c
c      界面を探す
       call sufele(nele,ne,idx,nr,nc,nelc
     &                    ,ndx,neled,na,mr,nb,nelb,neb,ncom)
c
c        call pltne(np,nele,ne,xx,yy)
c        call pltnb(na,nb,nelb,neb,ncom)
c-----------------------------------------------------------------------------
c
c      計算ループ
c
c-----------------------------------------------------------------------------
100    continue
c
       write(*,*)"------------------------------"
       write(*,*)"ti:",sngl(ti)
c
c      領域ループ
       do kr=1,nr
c
c        領域krの要素
c        nt,nelt
c        net(nelt,3)
c        xt(nt),yt(nt)
c        ut(nt),vt(nt),pt(nt)
         nelt=0
         do kc=1,nc(kr)
           nelt=nelt+1
           net(nelt,1)=ne(nelc(kr,kc),1)
           net(nelt,2)=ne(nelc(kr,kc),2)
           net(nelt,3)=ne(nelc(kr,kc),3)
           idt(nelt)=idx(nelc(kr,kc))
           iq(nelt)=nelc(kr,kc)     !iq(新)=古
           jq(nelc(kr,kc))=nelt     !jq(古)=新
c           write(*,*)"kc:",kc,nc(kr),nelc(kr,kc)
c           write(*,*)"kc:",kc,nc(kr)
c     &                      ,net(nelt,1),net(nelt,2),net(nelt,3)
c           pause
         enddo
c
         nt=0
         do 200 i=1,np
           ich=0
           do kelt=1,nelt
             if(net(kelt,1).eq.i)ich=1
             if(net(kelt,2).eq.i)ich=1
             if(net(kelt,3).eq.i)ich=1
           enddo
           if(ich.eq.0)goto 200
           nt=nt+1
           xt(nt)=xx(i)
           yt(nt)=yy(i)
           ut(nt)=uu(i)
           vt(nt)=vv(i)
           pt(nt)=pp(i)
           ip(nt)=i     !ip(新)=古
           jp(i)=nt     !ip(古)=新
200      continue
c
         do kelt=1,nelt
           net(kelt,1)=jp(net(kelt,1))
           net(kelt,2)=jp(net(kelt,2))
           net(kelt,3)=jp(net(kelt,3))
c           write(*,*)"kelt:",kelt,nelt
c     &                      ,net(kelt,1),net(kelt,2),net(kelt,3)
c           pause
         enddo
c
c       バンド幅
c       nw,nelw
c       new(nelw,3)
c       xw(nw),yw(nw)
c       uw(nw),vw(nw),pw(nw)
c        call renum2(index,np,ne,nel,jww,ine,ind,nband)
c        do kt=1,nt
c          jvv(jww(kt))=kt    !jvv(新_新)=新
c        enddo
c
        do kt=1,nt
          jww(kt)=nt-kt+1    !jww(新)=新_新
          jvv(jww(kt))=kt    !jvv(新_新)=新
        enddo
c
        nbm=-99999
        do kele=1,nelt
        do i=1,3
        do j=1,3
          kbm=iabs(net(kele,i)-net(kele,j))
          if(kbm.gt.nbm)nbm=kbm     !nbm:現在考慮しているtreeのバンド幅
        enddo
        enddo
        enddo
        nbw=nbm+1
        write(*,*)'band',nbw
c
        nw=nt
        do kw=1,nw  !新_新
          xw(kw)=xt(jvv(kw))
          yw(kw)=yt(jvv(kw))
          uw(kw)=ut(jvv(kw))
          vw(kw)=vt(jvv(kw))
          pw(kw)=pt(jvv(kw))
        enddo
c
        nelw=nelt
        do kelw=1,nelw
          new(kelw,1)=jww(net(kelw,1))
          new(kelw,2)=jww(net(kelw,2))
          new(kelw,3)=jww(net(kelw,3))
          idw(kelw)=idt(kelw)
        enddo
c
c       境界の節点
c       la,mb(la),melb(la,mb(la)),meb(la,mb(la))
        la=0
        do 300 ka=1,na
          if(kr.ne.mr(ka))goto 300
          la=la+1
          mb(la)=nb(ka)
          do kb=-1,nb(ka)+2
            melb(la,kb)=jq(nelb(ka,kb))
            meb(la,kb)=neb(ka,kb)
            mcom(la,kb)=ncom(ka,kb)
          enddo
300     continue
c
c        call pltne(nw,nelw,new,xw,yw)
c        call pltnb(la,mb,melb,meb,mcom)
c
c       マトリックス
c       uvp0(nw*3)
c       xn(nw),yn(nw)
        do kw=1,nw
          uvp0(kw*3-2)=uw(kw)
          uvp0(kw*3-1)=vw(kw)
          uvp0(kw*3-0)=pw(kw)
          uvp(kw*3-2)=uw(kw)
          uvp(kw*3-1)=vw(kw)
          uvp(kw*3-0)=pw(kw)
          xn(kw)=xw(kw)
          yn(kw)=yw(kw)
        enddo
c
c       第1種の境界条件
c       nbd,mbd(nbd),vbd(nbd)
        call boundf1(new,la,mb,melb,meb,mcom,nbd,mbd,vbd)
        do i=1,nbd
          uvp0(mbd(i))=vbd(i)
        enddo
c
c       流動解析
c       nelw,new(nelw,3),idw(nelw)
c       nw,xw(nw),yw(nw)
c       uw(nw),vw(nw),pw(nw)
c       la,mb(la),melb(la,mb(la)),meb(la,mb(la))
        do 400 lk=1,26
c         マトリックスの初期化
          call setmat(nw,nbw,sa,sf,dsf)
c         マトリックスの作成
          call flow(nelw,nbw*3,new,idw,xn,yn,uvp0,dt,vre,vma,vfr,sa,sf)
          
c          nbwm=nbw*3-1           !境界上の節点速度0の条件
c          do ii=1,nw*3
c          do j=max(1,ii-nbwm),min(ii+nbwm,nw*3)
c            write(*,*)"sa:",ii,j,sa(ii,nbw+j-ii)
c            pause
c          enddo
c          enddo
c          do i=1,nw
c            write(*,*)"sf:",i,sngl(sf(i*3-2))
c     &                        ,sngl(sf(i*3-1))
c     &                        ,sngl(sf(i*3-0))
c            pause
c          enddo
c
c         第2種境界条件(自由表面)
c         dsf(nw*3)
          call boundfs(kr,ndx,new,xn,yn,vwe
     &                      ,la,mb,melb,meb,mcom,dsf)
c         境界条件代入
          call boundasym(nw*3,nbw*3,nbd,mbd,vbd,dsf,sa,sf)
          
c         ガウスの消去法
          call gauss(nw*3,nbw*3,sa,sf)
c
c         収束誤差の計算
          err=0d0
          do 290 i=1,nw
            errt=dabs(sf(3*i-2)-uvp(3*i-2))
     &          +dabs(sf(3*i-1)-uvp(3*i-1))
     &          +dabs(sf(3*i-0)-uvp(3*i-0))
            if(errt.gt.err)then
              err=errt
              kkk=i
c              write(*,*)kkk,err
c              write(*,*)dabs(sf(3*i-2)-uvp(3*i-2)),
c     &                   dabs(sf(3*i-1)-uvp(3*i-1)),
c     &                   dabs(sf(3*i-0)-uvp(3*i-0))
c              write(*,*)"-----------------------------"
c              write(*,*)sngl(sf(3*i-2)),sngl(uvp(3*i-2))
c              write(*,*)sngl(sf(3*i-1)),sngl(uvp(3*i-1))
c              write(*,*)sngl(sf(3*i-0)),sngl(uvp(3*i-0))
c              write(*,*)sngl(xn(i)),sngl(yn(i))
c              pause
            endif
            zkar=zkasu
            if(err.lt.1d-6)zkar=zkasu*4.5d-1                    !移動後の節点位置において速度を計算し，
            uvp(3*i-2)=uvp(3*i-2)+(sf(3*i-2)-uvp(3*i-2))*zkar   !速度から求めた新たな移動後の節点位置
            uvp(3*i-1)=uvp(3*i-1)+(sf(3*i-1)-uvp(3*i-1))*zkar   !との差が小さくなった場合，流動ループ
            uvp(3*i-0)=uvp(3*i-0)+(sf(3*i-0)-uvp(3*i-0))*zkar   !を抜ける.
            xn(i)=xw(i)+(uvp0(3*i-2)+uvp(3*i-2))*dt/2d0
            yn(i)=yw(i)+(uvp0(3*i-1)+uvp(3*i-1))*dt/2d0
 290      continue
c
          write(*,'(a3,i3,a3,i5,a5,e20.10)') 
     &                     'lk',lk,"kp",kkk,"err",sngl(err)
c         収束判定
          if(err.lt.1d-7) then
            open(10,file="err.res",access='append')
            write(10,*)"ti",sngl(ti),"  lk",lk,"  err",sngl(err)
            close(10)
            goto 500
          endif
400    continue    !収束ループ
c
       open(10,file="err.res",access='append')
       write(10,*)"ti",sngl(ti),"  lk",lk,"  err",sngl(err),
     &               "   No convergion"
       close(10)
       write(*,*) ' !!!No convergion in lk loop/err=',sngl(err)
500    continue
c
c      値コピー
       do kw=1,nw   !新_新
         xx(ip(jvv(kw)))=xn(kw)  !ip(新)=古,jvv(新_新)=新
         yy(ip(jvv(kw)))=yn(kw)
         uu(ip(jvv(kw)))=uvp(3*kw-2)
         vv(ip(jvv(kw)))=uvp(3*kw-1)
         pp(ip(jvv(kw)))=uvp(3*kw-0)
       enddo
c
       call pltne(np,nele,ne,xx,yy)
       call pltnb(na,nb,nelb,neb,ncom)
c-------------------------------------------------------------------------------
c
c      要素の再発生
c
c-------------------------------------------------------------------------------
       if(icheck(nele,ne,xx,yy).eq.1)then
         write(*,*)"mesh is regenerated."
c        要素の再発生
         call grid(np,nele,ne,idx,xx,yy,uu,vv,pp,nr
     &                 ,nc,nelc,ndx,mr,na,nb,neb,nelb)
c        界面を探す
         call sufele(nele,ne,idx,nr,nc,nelc
     &                      ,ndx,neled,na,mr,nb,nelb,neb,ncom)
c
         open(10,file="err.res",access='append')
         write(10,*)"ti",sngl(ti),"  remesh"
         close(10)
c         pause
       endif
c-------------------------------------------------------------------------------
c
c      データ保存
c
c-------------------------------------------------------------------------------
c      tmp fileの保存
       if(ti.gt.(ttmp-dt*1d-2) .and. ti.lt.(ttmp+dt*1d-2))then
         write(*,*)'**** save tmp file data ****'
c         write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
         write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
         do i=1,12
           if(lpp(i:i).eq.' ') lpp(i:i)='0'
	   enddo
         call addchr("tmp",lpp//'.inp',tfile)
         open(1,file=tfile)
         rewind(1)
         write(1,*)ti,tsav,ttmp
         write(1,*)np,nele
         write(1,*)(ne(i,1),ne(i,2),ne(i,3),i=1,nele)
         write(1,*)(idx(i),i=0,nele)
         write(1,*)(xx(i),yy(i),i=1,np)
         write(1,*)(uu(i),vv(i),pp(i),i=1,np)
         close(1)
         ttmp=ttmp+wdt2
       endif
c-------------------------------------------------------------------------------
c
c      タイムステップ
c
c-------------------------------------------------------------------------------
       enddo   !タイムステップ
c
       ti=ti+dt
       if(ti.lt.tmax)goto 100
c
       stop
       end
