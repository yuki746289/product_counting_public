       include 'addchr.for'
       include "header.h"
       dimension ne(ine,4),new(ine,4),nen(ine,4)
       dimension xx(ind),yy(ind),zz(ind)
       dimension xn(ind),yn(ind),zn(ind)
       dimension idx(0:ine),idn(0:ine)
       dimension ip(ind),jp(ind),ie(ine),je(ine)
       dimension jww(ind),jmm(ind)
       dimension uu(ind),vv(ind),ww(ind),pp(ind)
       dimension u0(ind),v0(ind),w0(ind),p0(ind)
       dimension un0(ind),vn0(ind),wn0(ind),pn0(ind)
       dimension uvwp(4*ind),uvwp0(4*ind)
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension mn(ind)  !mn(nn):境界条件ne(nelb(nb(1)),neb(nb(1))を代入する節点番号
       character*40 rfile,tfile
       character*12 lpp
c------記録用--------------------------------------------------------
       dimension pnx(ind),pny(ind),pnz(ind)
       dimension vnx(ind),vny(ind),vnz(ind),hh(ind)
c------物質移動--------------------------------------------------------
       dimension nec(ine,4)        !新
       dimension nea(ine,4),xa(ind),ya(ind),za(ind)    !新_新
       dimension cc(ind),c0(ind),ca0(ind),ta0(ind),idc(0:ine)
       dimension ipc(ind),jpc(ind),iec(ine),jec(ine)
       dimension jwc(ind),jmc(ind)
       dimension mc(ind)
c------熱移動------------------------------------------
       dimension tn0(ind)      !流動領域
       dimension net(ine,4)    !新
       dimension neb(ine,4),xb(ind),yb(ind),zb(ind)    !新_新
       dimension tt(ind),t0(ind),cb0(ind),tb0(ind),idt(0:ine)
       dimension ipt(ind),jpt(ind),iet(ine),jet(ine)
       dimension jwt(ind),jmt(ind)
       dimension mt(ind)
c------メッシュ切り直し------------------------------------------
       dimension un(ind),vn(ind),wn(ind)
       dimension cn(ind),tn(ind),pn(ind)
       dimension ms(ine,3)
c------------------------------------------------------
       common /time/dt,timax,tsv,ttmp
       common /radian/phi2(ind),psi2(ind),eta2(ind)    !曲率の計算
       data zcft/1d0/   !緩和係数
       data emax/1d-7/  !収束判定:時間間隔よりも小さくする
c----------------------------------------------------------------------
c
c      初期化およびdata読み込み
c
c----------------------------------------------------------------------
       open(10,file='err.res')  !err.resファイルの作成
       rewind(10)
c
       do 10 i=0,ine   !インデックス番号の初期化
10     idx(i)=0     !idx(0:nele)
       call data(np,nele,ne,xx,yy,zz,nbw,iyy,ivs,idx,rfile)
c------座標系の選択----------------------------------------------------
       icoord=1      !0:オイラー座標,1:ラグラジアン座標
c------形状関数の計算--------------------------------------------------
       imt=1        !0:形状関数の計算を行わない,1:形状関数の計算を行う
c      オイラー座標:形状関数を1回計算する(対流項除く),ラグラジアン座標:形状関数を毎時間計算する
c------記録時間--------------------------------------------------------
       tsav=tsv     !dataを保存する最初の時間
       ttp=ttmp     !tmp fileを保存する最初の時間
c------流動解析--------------------------------------------------------
       icalf=1      !0:流動解析を行わない,1:流動解析を行う
       icdf=1       !0:流動領域に節点を割り振り直さない,1:流動領域に節点を割り振り直す
       call bdmodify(1,nn,ne,mn,ine,ind)     !流動解析に関する境界条件の修正
c------物質移動--------------------------------------------------------
       icalc=0      !0:物質移動解析を行わない,1:物質移動解析を行う
       icnt=1       !0:物質移動領域に節点を割り振り直さない,1:物質移動領域に節点を割り振り直す
       call bdmodify(2,nc,ne,mc,ine,ind)     !物質移動解析に関する境界条件の修正
c------熱移動----------------------------------------------------------
       icalt=0      !0:熱移動解析を行わない,1:熱移動解析を行う
       iheat=1      !0:熱移動領域に節点を割り振り直さない,1:熱移動領域に節点を割り振り直す
       call bdmodify(3,nt,ne,mt,ine,ind)     !熱移動解析に関する境界条件の修正
c----------------------------------------------------------------------
c
c      最初から計算
c
c----------------------------------------------------------------------
       if(iyy.eq.0)then
c         write(*,*)'iyy:0'
         ti=0d0
         call initf(np,nele,ne,u0,v0,w0,p0,uu,vv,ww,pp,idx,ind,ine)
         call initc(np,nele,ne,c0,cc,idx,ind,ine)
         call initt(np,nele,ne,t0,tt,idx,ind,ine)
c----------------------------------------------------------------------
c
c      途中から計算
c
c----------------------------------------------------------------------
       elseif(iyy.eq.1)then
c         write(*,*)'iyy:1'
         call addchr(rfile,'.tmp',tfile)
c         open(1,file='calc.tmp',status='unknown',access='sequential',
c     $            form='unformatted')
         open(1,file=tfile)
         rewind(1)
         read(1,*)ti,tsav,ttp
         read(1,*)np,nele,nbw
         read(1,*)(ne(i,1),ne(i,2),ne(i,3),ne(i,4),i=1,nele)
         read(1,*)(idx(i),i=0,nele)
         read(1,*)(xx(i),yy(i),zz(i),i=1,np)
c-------流動----------
         read(1,*)(u0(i),v0(i),w0(i),p0(i),i=1,np)
c------物質移動-------
         read(1,*)(c0(i),i=1,np)
c-------熱移動--------
         read(1,*)(t0(i),i=1,np)
c---------------------
         close(1)
         tsav=tsav+tsv
         ttp=ttp+ttmp
       endif
c---------------------------------------------------------------------
c
c      表面を探す
c
c---------------------------------------------------------------------
        call findms(nele,ne,nels,ms)
c---------------------------------------------------------------------
c
c      節点位置の修正
c
c---------------------------------------------------------------------
        icount=0
3700    continue
        if(elshape2(np,nele,ne,xx,yy,zz,ind,ine).gt.7d-1)then
          write(*,*)"start remesh:",np,nele
          call remesh2(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
          write(*,*)"end remesh:",np,nele
          dr=elshape2(np,nele,ne,xx,yy,zz,ind,ine)
          if(dr.gt.7d-1)then
            icount=icount+1
            write(*,*)"icount:",icount
            if(icount.eq.10)stop
            goto 3700
          endif
          write(10,*)'remesh:',icount
c          pause
        endif
c---------------------------------------------------------------------
c
c      曲率計算で使用する角度の初期化
c
c---------------------------------------------------------------------
       do i=1,ind
         phi2(i)=-999d0
         psi2(i)=-999d0
         eta2(i)=-999d0
       enddo
c*******************************************************************
c
c      時間ループの開始
c
c*******************************************************************
       write(*,*)'started time loop'
       do 2004 it=1,100000
       ti=ti+dt
       write(*,*)'*************************************************'
       write(*,*)'ti',ti
c------形状関数の計算-----------------------------------------------
       if(icoord.eq.0 .and. it.ge.2)imt=0  !オイラー座標において，形状関数を1回計算する
c-------------------------------------------------------------------
c*******************************************************************
c
c      流動解析
c
c*******************************************************************
c----------------------------------------------------------------------
       if(icalf.eq.0)goto 40000    !0:流動解析を行わない,1:流動解析を行う
c----------------------------------------------------------------------
       write(*,*)'---- calc flow ----'
       if(icdf.eq.0)goto 10000     !既に流動領域に節点番号が割り振り直されている場合
c----------------------------------------------------------------------
c
c      流動領域の節点に新しい番号を割り振る
c
c----------------------------------------------------------------------
c       write(*,*)"np"
c       call pltrset()
c       call pltsetnp(np,xx,yy,zz,ind)
c       pause
c
       nelw=0
       npw=0
       do 1000 i=1,np
         ip(i)=0
         jp(i)=0
1000   continue
c
       do 1100 kele=1,nele  !古
         if(idx(kele).eq.1 .or. idx(kele).eq.2 .or.  !index番号が流動領域の場合
     &      idx(kele).eq.3 .or. idx(kele).eq.4)then
           nelw=nelw+1
           do 1200 i=1,4
             if(ip(ne(kele,i)).eq.0)then
               npw=npw+1
               new(nelw,i)=npw
               ip(ne(kele,i))=npw     !ip(古)=新
               jp(npw)=ne(kele,i)     !jp(新)=古
             elseif(ip(ne(kele,i)).ne.0)then
               new(nelw,i)=ip(ne(kele,i))
             endif
1200       continue
           ie(kele)=nelw          !ie(古)=新
           je(nelw)=kele          !je(新)=古
         endif
1100   continue
c----------------------------------------------------------------------
c
c      バンド幅最適化
c
c----------------------------------------------------------------------
       call renum3(npw,nelw,new,jww,jmm,nbw)   !jww(新)=新_新
c***************************************************
c        do i=1,npw
c          jww(i)=i
c          jmm(i)=i
c        enddo
c
c        バンド幅を求める
c         write(*,*)'calc band'
c         nbm=-99999
c         do kele=1,nelw
c           do i=1,4
c             do j=1,4
c               kbm=iabs(jww(new(kele,i))-jww(new(kele,j)))
c               if(kbm.gt.nbm)nbm=kbm     !nbm:現在考慮しているtreeのバンド幅
c             enddo
c           enddo
c         enddo
c         nbw=nbm+1
c         write(*,*)'band',nbw
c***************************************************
       idn(0)=0
       do 1300 kele=1,nelw   !新_新(要素は新=新_新)
         idn(kele)=idx(je(kele))    !je(新_新)=古
         do 1300 i=1,4              !idn(0:nelw)
         nen(kele,i)=jww(new(kele,i))
1300   continue
       do 1400 kp=1,npw     !新_新
         xn(kp)=xx(jp(jmm(kp)))     !jp(新)=古
         yn(kp)=yy(jp(jmm(kp)))     !jmm(新_新)=新
         zn(kp)=zz(jp(jmm(kp)))
c         pn0(kp)=p0(jp(jmm(kp)))    !jp(jmm(新_新))=古
c         un0(kp)=u0(jp(jmm(kp)))
c         vn0(kp)=v0(jp(jmm(kp)))
c         wn0(kp)=w0(jp(jmm(kp)))
         tn0(kp)=t0(jp(jmm(kp)))    !物性値を計算するために用いる
c         write(*,*)kp,jmm(kp),jp(jmm(kp))
c         pause
1400  continue
      npn=npw          !npn:流動領域の節点数
      neln=nelw        !neln:流動領域の要素数
c      icdf=0           !流動領域に節点番号を割り振り直さない
      write(*,*)"npn,nelw",npn,nelw
10000 continue
c
c      write(*,*)"npn"
c      call pltrset()
c      call pltsetnp(npn,xn,yn,zn,ind)
c      call pltsetne(neln,nen,2d0,2,ine)
c      pause
c----------------------------------------------------------------------
c
c     マトリクスに初期値を代入する
c
c----------------------------------------------------------------------
      do 1500 kp=1,npn
        uvwp0(4*kp-3)=u0(jp(jmm(kp)))       !uvwp0(4*新_新-i)
        uvwp0(4*kp-2)=v0(jp(jmm(kp)))
        uvwp0(4*kp-1)=w0(jp(jmm(kp)))
        uvwp0(4*kp-0)=p0(jp(jmm(kp)))
c        sf(kp)=p0(jp(jmm(kp)))
1500  continue
c
      do 1600 kp=1,npn
        uvwp(4*kp-3)=uvwp0(4*kp-3)
        uvwp(4*kp-2)=uvwp0(4*kp-2)
        uvwp(4*kp-1)=uvwp0(4*kp-1)
        uvwp(4*kp-0)=uvwp0(4*kp-0)
1600  continue
c----------------------------------------------------------------------
c
c     繰り返し計算
c
c----------------------------------------------------------------------
      do 1700 lk=1,26
c
c       N-S式のマトリクスを解く  sa(4*ind,ibw),sf(4*ind)
        call setmat(sa,sf,npn,4*nbw)                                   !マトリクスを初期化する
        call flow(neln,4*nbw,sa,sf,nen,xn,yn,zn                        !速度u,v,pのマトリクスsa,列ベクトルsfを発生させる
     &             ,uvwp0,tn0,idn,ind,ine,ibw)                          !nbw=4*nbw
        call boundf1(npn,ne,idn,sa,sf,4*nbw,ip,jww,ind,ine,ibw)        !マトリクスsaにu,v,w,pの境界条件を代入
c        call boundfs(npn,np,nele,ne,xx,yy,zz,idn,sa,sf,4*nbw,ip,jww)   !マトリクスsaに自由表面の境界条件を代入
        call boundfs(npn,neln,nen,xn,yn,zn,idn,sa,sf,4*nbw,ip,jww
     &                                              ,hh,vnx,vny,vnz)    !マトリクスsaに自由表面の境界条件を代入
        call gauss(4*npn,4*nbw,sa,sf,4*ind,ibw)                        !ガウスの消去法でマトリクスを解く
        call boundf1r(ne,sf,ip,jww,ind,ine)                            !列ベクトルsfにu,v,pの境界条件を代入
c
c       速度u,v,w,pの残差を求める
        err=-999d0
        umx=-999d0                                           !列ベクトルに用いる速度:uvwp0(4*ind)
        vmx=-999d0                                           !対流項に用いる速度:uvwp(4*ind)
        wmx=-999d0
        pmx=-999d0
        do i=1,npn
          ee=dabs(uvwp(4*i-3)-sf(4*i-3))+dabs(uvwp(4*i-2)-sf(4*i-2))
     &      +dabs(uvwp(4*i-1)-sf(4*i-1))+dabs(uvwp(4*i-0)-sf(4*i-0))
          if(err.lt.ee)err=ee
          umx=dmax1(umx,dabs(sf(4*i-3)))
          vmx=dmax1(vmx,dabs(sf(4*i-2)))
          wmx=dmax1(wmx,dabs(sf(4*i-1)))
          pmx=dmax1(pmx,dabs(sf(4*i-0)))
          uvwp(4*i-3)=uvwp(4*i-3)+(sf(4*i-3)-uvwp(4*i-3))*zcft
          uvwp(4*i-2)=uvwp(4*i-2)+(sf(4*i-2)-uvwp(4*i-2))*zcft   !速度から求めた新たな移動後の節点位置
          uvwp(4*i-1)=uvwp(4*i-1)+(sf(4*i-1)-uvwp(4*i-1))*zcft   !との差が小さくなった場合，流動ループ
          uvwp(4*i-0)=uvwp(4*i-0)+(sf(4*i-0)-uvwp(4*i-0))*zcft           !を抜ける.
          xn(i)=xx(jp(jmm(i)))+(uvwp0(4*i-3)+uvwp(4*i-3))*dt/2d0
          yn(i)=yy(jp(jmm(i)))+(uvwp0(4*i-2)+uvwp(4*i-2))*dt/2d0
          zn(i)=zz(jp(jmm(i)))+(uvwp0(4*i-1)+uvwp(4*i-1))*dt/2d0
c          write(*,*)"i:",i,"/",npn
c          write(*,*)"uvwp:",sngl(uvwp(4*i-3)),sngl(uvwp(4*i-2))
c     &                      ,sngl(uvwp(4*i-1)),sngl(uvwp(4*i-0))
c          write(*,*)"xx:",sngl(xn(i)),sngl(yn(i)),sngl(zn(i))
c          pause
        enddo
c
c       表示
        write(*,*)"----------------------------------------------------"
        write(*,*)'lk,err',lk,err
        write(*,*)'umx,vmx,wmx,pmx',sngl(umx),sngl(vmx),
     &                               sngl(wmx),sngl(pmx)
c
        if(err.lt.emax)goto 2000   !終了判定
1700  continue
      write(*,*)'***No convergion lk=',lk,'err=',sngl(err),'***'
2000  continue     !繰り返し計算終了
c
c     残差の記録
      if(err.ge.emax)
     &     write(10,*)'ti',sngl(ti),'lk',lk,'err',err,'No convergion'
      if(err.lt.emax)write(10,*)'ti',sngl(ti),'lk',lk,'err',err
c----------------------------------------------------------------------
c
c     求めた値を配列に代入する
c
c----------------------------------------------------------------------
      do 2100 kp=1,npn     !新_新
        xx(jp(jmm(kp)))=xn(kp)
        yy(jp(jmm(kp)))=yn(kp)
        zz(jp(jmm(kp)))=zn(kp)
        uu(jp(jmm(kp)))=uvwp(4*kp-3)     !jp(新)=古
        vv(jp(jmm(kp)))=uvwp(4*kp-2)     !jmm(新_新)=新
        ww(jp(jmm(kp)))=uvwp(4*kp-1)
        pp(jp(jmm(kp)))=uvwp(4*kp-0)     !jp(jmm(新_新))=古
c        pp(jp(jmm(kp)))=sf(kp)
        pnx(jp(jmm(kp)))=vnx(kp)*hh(kp)     !表面張力のベクトル
        pny(jp(jmm(kp)))=vny(kp)*hh(kp)     !表面張力のベクトル
        pnz(jp(jmm(kp)))=vnz(kp)*hh(kp)     !表面張力のベクトル
2100  continue
40000 continue
c*******************************************************************
c
c      物質移動解析
c
c*******************************************************************
c----------------------------------------------------------------------
       if(icalc.eq.0)goto 30000    !0:物質移動解析を行わない,1:物質移動解析を行う
c----------------------------------------------------------------------
       write(*,*)'---- calc mass ----'
       if(icnt.eq.0)goto 20000     !既に物質移動領域に節点番号が割り振り直されている場合
c----------------------------------------------------------------------
c
c      物質移動領域の節点に新しい番号を割り振る
c
c----------------------------------------------------------------------
       nelc=0
       npc=0
       do 2400 i=1,np
         ipc(i)=0
         jpc(i)=0
2400   continue
c
       do 2500 kele=1,nele  !古
         if(idx(kele).eq.2 .or. idx(kele).eq.4 .or.  !index番号が物質移動領域の場合
     &      idx(kele).eq.5 .or. idx(kele).eq.6)then
           nelc=nelc+1
           do 2600 i=1,4
             if(ipc(ne(kele,i)).eq.0)then
               npc=npc+1
               nec(nelc,i)=npc
               ipc(ne(kele,i))=npc     !ipc(古)=新
               jpc(npc)=ne(kele,i)     !jpc(新)=古
             elseif(ipc(ne(kele,i)).ne.0)then
               nec(nelc,i)=ipc(ne(kele,i))
             endif
2600       continue
           iec(kele)=nelc          !iec(古)=新
           jec(nelc)=kele          !jec(新)=古
         endif
2500   continue
c----------------------------------------------------------------------
c
c      バンド幅最適化
c
c----------------------------------------------------------------------
       call renum3(npc,nelc,nec,jwc,jmc,nbc)   !jwc(新)=新_新
c                                               !jmc(新_新)=新
       idc(0)=0
       do 2700 kele=1,nelc   !新_新(要素は新=新_新)
         idc(kele)=idx(jec(kele))   !jec(新_新)=古
         do 2700 i=1,4              !idc(0:nelw)
         nea(kele,i)=jwc(nec(kele,i))
2700   continue
       do 2800 kp=1,npc     !新_新
         xa(kp)=xx(jpc(jmc(kp)))     !jpc(新)=古
         ya(kp)=yy(jpc(jmc(kp)))     !jmc(新_新)=新
         za(kp)=zz(jpc(jmc(kp)))
c         ca0(kp)=c0(jpc(jmc(kp)))     !jpc(jmc(新_新))=古
c         ta0(kp)=t0(jpc(jmc(kp)))     !物性値を計算するために用いる
c         write(*,*)kp,jmc(kp),jpc(jmc(kp))
c         pause
2800  continue
      npa=npc          !npa:物質移動領域の節点数
      nela=nelc        !nela:物質移動領域の要素数
      icnt=0           !物質移動領域に節点番号を割り振り直さない
20000 continue
c----------------------------------------------------------------------
      do kp=1,npa     !新_新
        ca0(kp)=c0(jpc(jmc(kp)))     !jpc(jmc(新_新))=古
        ta0(kp)=t0(jpc(jmc(kp)))     !物性値を計算するために用いる
      enddo
c----------------------------------------------------------------------
c
c      濃度マトリクスを発生させ，濃度を計算する
c
c----------------------------------------------------------------------
       call setmat(sa,sf,npa,nbc)                              !マトリクスを初期化する
       call mtmass(nela,nbc,sa,sf,nea,xa,ya,za                 !nea(ine,4):物質移動領域の要素のみ格納
     &               ,uu,vv,ww,ca0,ta0,jmc,jpc,idc,ind,ine,ibw) !マトリクスsa, 列ベクトルsfを発生させる
       call boundm(npa,ne,idc,sa,sf,nbc,ipc,jwc,nc,mc)         !マトリクスsaに境界条件を代入する
       call gauss(npa,nbc,sa,sf,4*ind,ibw)                     !ガウスの消去法でマトリクスsaを解く
       call boundmr(ne,sf,ipc,jwc,nc,mc,ind,ine)               !列ベクトルsfに境界条件を代入する
c
      do 2900 kp=1,npa     !新_新
        cc(jpc(jmc(kp)))=sf(kp)           !jpc(jmc(新_新))=古
2900  continue
c
30000 continue
c*******************************************************************
c
c      熱移動解析
c
c*******************************************************************
c----------------------------------------------------------------------
       if(icalt.eq.0)goto 50000    !0:熱移動解析を行わない,1:熱移動解析を行う
c----------------------------------------------------------------------
       write(*,*)'---- calc heat ----'
       if(iheat.eq.0)goto 60000     !既に熱移動領域に節点番号が割り振り直されている場合
c----------------------------------------------------------------------
c
c      熱移動領域の節点に新しい番号を割り振る
c
c----------------------------------------------------------------------
       nelt=0
       npt=0
       do 3000 i=1,np
         ipt(i)=0
         jpt(i)=0
3000   continue
c
       do 3100 kele=1,nele  !古
         if(idx(kele).eq.3 .or. idx(kele).eq.4 .or.  !index番号が熱移動領域の場合
     &      idx(kele).eq.5 .or. idx(kele).eq.7)then
           nelt=nelt+1
           do 3200 i=1,4
             if(ipt(ne(kele,i)).eq.0)then
               npt=npt+1
               net(nelt,i)=npt
               ipt(ne(kele,i))=npt     !ipt(古)=新
               jpt(npt)=ne(kele,i)     !jpt(新)=古
             elseif(ipt(ne(kele,i)).ne.0)then
               net(nelt,i)=ipt(ne(kele,i))
             endif
3200       continue
           iet(kele)=nelt          !iet(古)=新
           jet(nelc)=kele          !jet(新)=古
         endif
3100   continue
c----------------------------------------------------------------------
c
c      バンド幅最適化
c
c----------------------------------------------------------------------
       call renum3(npt,nelt,net,jwt,jmt,nbt)   !jwt(新)=新_新
c                                               !jmt(新_新)=新
       idt(0)=0
       do 3300 kele=1,nelt   !新_新(要素は新=新_新)
         idt(kele)=idx(jet(kele))   !jet(新_新)=古
         do 3300 i=1,4              !idt(0:nelw)
         neb(kele,i)=jwt(net(kele,i))
3300   continue
       do 3400 kp=1,npt     !新_新
         xb(kp)=xx(jpt(jmt(kp)))     !jpt(新)=古
         yb(kp)=yy(jpt(jmt(kp)))     !jmt(新_新)=新
         zb(kp)=zz(jpt(jmt(kp)))
c         cb0(kp)=c0(jpt(jmt(kp)))     !物性値を計算するために用いる
c         tb0(kp)=t0(jpt(jmt(kp)))     !jpt(jmt(新_新))=古
c         write(*,*)kp,jmt(kp),jpt(jmt(kp))
c         pause
3400  continue
      npb=npt          !npb:熱移動領域の節点数
      nelb=nelt        !nelb:熱移動領域の要素数
      iheat=0          !熱移動領域に節点番号を割り振り直さない
60000 continue
c----------------------------------------------------------------------
      do kp=1,npb     !新_新
        cb0(kp)=c0(jpt(jmt(kp)))     !jpt(jmt(新_新))=古
        tb0(kp)=t0(jpt(jmt(kp)))     !物性値を計算するために用いる
      enddo
c----------------------------------------------------------------------
c
c      温度マトリクスを発生させ，温度を計算する
c
c----------------------------------------------------------------------
       call setmat(sa,sf,npb,nbt)                              !マトリクスを初期化する
       call mtheat(nelb,nbt,sa,sf,neb,xb,yb,zb                 !nea(ine,4):物質移動領域の要素のみ格納
     &               ,uu,vv,ww,cb0,tb0,jmt,jpt,idt,ind,ine,ibw) !マトリクスsa, 列ベクトルsfを発生させる
       call boundt(npb,ne,idt,sa,sf,nbt,ipt,jwt,nt,mt)         !マトリクスsaに境界条件を代入する
       call gauss(npb,nbt,sa,sf,4*ind,ibw)                     !ガウスの消去法でマトリクスsaを解く
       call boundtr(ne,sf,ipt,jwt,nt,mt,ind,ine)               !列ベクトルsfに境界条件を代入する
c
      do 3500 kp=1,npb     !新_新
        tt(jpt(jmt(kp)))=sf(kp)           !jpt(jmt(新_新))=古
3500  continue
c
50000 continue
c**********************************************************************
c
c     メッシュの切り直し
c
c**********************************************************************
c        メッシュ切り直し
c        dr=sfshape(nels,ms,xx,yy,zz)
c        dummy=elshape2(np,nele,ne,xx,yy,zz,ind,ine)
c        if(elshape2(np,nele,ne,xx,yy,zz,ind,ine).gt.7d-1)then
c          write(*,*)"start remesh:",np,nele
c          call copydata(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt,idx
c     &                ,npn,nen,neln,xn,yn,zn,un,vn,wn,pn,cn,tn,idn) !変数のコピー
c          call remesh(np,ne,nele,xx,yy,zz)
c          write(*,*)"end remesh:",np,nele
c          call interpolation(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt,idx
c     &                    ,npn,nen,neln,xn,yn,zn,un,vn,wn,pn,cn,tn,idn)
c          call findms(nele,ne,nels,ms)
c          dummy=elshape(np,nele,ne,xx,yy,zz,ind,ine)
c          pause
c        endif
c
c       節点移動
        icount=0
3600    continue
        write(*,*)"vlmin:",calvlmin(np,ne,nele,xx,yy,zz)
        if(elshape2(np,nele,ne,xx,yy,zz,ind,ine).gt.7d-1)then
          vlmin=calvlmin(np,ne,nele,xx,yy,zz)   !最小体積を計算する
          write(*,*)"start remesh:",np,nele,sngl(vlmin)
          call remesh2(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
          vlmin=calvlmin(np,ne,nele,xx,yy,zz)   !最小体積を計算する
          write(*,*)"end remesh:",np,nele,sngl(vlmin)
          dr=elshape2(np,nele,ne,xx,yy,zz,ind,ine)
          if(dr.gt.7d-1)then
            icount=icount+1
            write(*,*)"icount:",icount
            if(icount.eq.10)stop
            goto 3600
          endif
          write(10,*)'remesh:',icount
c
c         曲率計算で使用する角度の初期化
          do i=1,ind
            phi2(i)=-999d0
            psi2(i)=-999d0
            eta2(i)=-999d0
          enddo
c          pause
        endif
c**********************************************************************
c
c     dataの保存
c
c**********************************************************************
c     tmp fileの保存
      if(ti.gt.(ttp-dt*1d-2) .and. ti.lt.(ttp+dt*1d-2))then
        write(*,*)'**** save tmp file data ****'
c        write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
        write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
        do i=1,12
          if(lpp(i:i).eq.' ') lpp(i:i)='0'
	  enddo
        call addchr(rfile,lpp//'.tmp',tfile)
c        open(1,file=tfile,status='unknown',access='sequential',
c     $           form='unformatted')
        open(1,file=tfile)
        rewind(1)
        write(1,*)ti,tsav,ttp
        write(1,*)np,nele,nbw
        write(1,*)(ne(i,1),ne(i,2),ne(i,3),ne(i,4),i=1,nele)
        write(1,*)(idx(i),i=0,nele)
        write(1,*)(xx(i),yy(i),zz(i),i=1,np)
        write(1,*)(uu(i),vv(i),ww(i),pp(i),i=1,np)
        write(1,*)(cc(i),i=1,np)   !物質移動
        write(1,*)(tt(i),i=1,np)   !熱移動
        close(1)
        ttp=ttp+ttmp
      endif
c     avs fileの保存
      if(ti.gt.(tsav-dt*1d-2) .and. ti.lt.(tsav+dt*1d-2))then
        write(*,*)'**** save avs file data ****'
        call svdata(np,nele,ne,xx,yy,zz,uu,vv,ww,pp,pnx,pny,pnz,idx,ti,
     &                                           rfile,ind,ine)
        if(icalc.eq.1)
     &    call svdatam(np,nele,ne,xx,yy,zz,cc,idx,ti,rfile,ind,ine)    !物質移動
        if(icalt.eq.1)
     &    call svdatat(np,nele,ne,xx,yy,zz,tt,idx,ti,rfile,ind,ine)    !熱移動
        tsav=tsav+tsv       !tsv:dataを保存する時間間隔
      endif
c----------------------------------------------------------------------
c
c     速度，圧力を更新する
c
c----------------------------------------------------------------------
      do 2300 i=1,np
        u0(i)=uu(i)
        v0(i)=vv(i)
        w0(i)=ww(i)
        p0(i)=pp(i)
c------物質移動-------
        c0(i)=cc(i)
c---------------------
c------熱移動-------
        t0(i)=tt(i)
c---------------------
2300  continue
c--------------------------------------------------------------------
c
c      時間ループ終了
c
c---------------------------------------------------------------------
2004  continue
      write(*,*)'calc end'
      close(10)    !err.resファイルを閉じる
      stop
      end
c --
c -