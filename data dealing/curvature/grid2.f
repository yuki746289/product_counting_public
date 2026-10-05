
c
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c       表面の節点を作る
c       表面:ns,xs(ns),ys(ns),zs(ns)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c

c---------------------------------------------------------------------------
c
c       球
c
c---------------------------------------------------------------------------
        subroutine gsuf1(ns,xs,ys,zs,dl0,ind)
        implicit double precision(a-h,o-z)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        dl0=2.2d0   !領域の幅
        r=1d0   !直径2d0
        ns=0    !節点の数
        nn=0    !節点にならなかった数
        do 10 while(1.eq.1)
c
c         乱数
          phi=drand(0)*2d0*PI   !0d0 <= drand(0) <= 1d0
          psi=drand(0)*2d0*PI
          x0=r*dcos(phi)*dcos(psi)
          y0=r*dsin(phi)
          z0=r*dcos(phi)*dsin(psi)
c
c         調べる
          do ks=1,ns
            dl=dsqrt((xs(ks)-x0)**2+(ys(ks)-y0)**2+(zs(ks)-z0)**2)
            if(dl.lt.1d-1)then
              nn=nn+1   !節点にならなかった数
              if(nn.eq.19000)return
              goto 10
            endif
          enddo
c
c         登録
          nn=0  !節点にならなかった数
          ns=ns+1
          xs(ns)=x0
          ys(ns)=y0
          zs(ns)=z0
10      continue
c
        end
c---------------------------------------------------------------------------
c
c       トーラス
c
c---------------------------------------------------------------------------
        subroutine gsuf2(ns,xs,ys,zs,dl0,ind)
        implicit double precision(a-h,o-z)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        dl0=2.2d0   !領域の幅
        r1=0.75d0   !トーラス中心円の半径
        r2=0.25d0   !トーラス断面の半径
        ns=0    !節点の数
        nn=0    !節点にならなかった数
        do 10 while(1.eq.1)
c
c         乱数
          phi=drand(0)*2d0*PI   !0d0 <= drand(0) <= 1d0 →  0<=φ<2π
          psi=drand(0)*2d0*PI   !0d0 <= drand(0) <= 1d0 →  0<=ψ<2π
          x0=r2*dcos(phi)*dcos(psi)+r1*dcos(psi)
          y0=r2*dsin(phi)
          z0=r2*dcos(phi)*dsin(psi)+r1*dsin(psi)
c
c         調べる
          do ks=1,ns
            dl=dsqrt((xs(ks)-x0)**2+(ys(ks)-y0)**2+(zs(ks)-z0)**2)
            if(dl.lt.1d-1)then
              nn=nn+1   !節点にならなかった数
              if(nn.eq.19000)return
              goto 10
            endif
          enddo
c
c         登録
          nn=0  !節点にならなかった数
          ns=ns+1
          xs(ns)=x0
          ys(ns)=y0
          zs(ns)=z0
10      continue
c
        end

c
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c       格子を発生させる(8分木)
c       格子:nt,xt(nt),yt(nt),zt(nt),dl(nt)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c




c---------------------------------------------------------------------------
c
c       8分木:格子を発生させる
c
c---------------------------------------------------------------------------
        subroutine octtree(ns,xs,ys,zs,nt,xt,yt,zt,dl0,ind)
        implicit double precision(a-h,o-z)
        dimension mt(ind)
        dimension xs(ind),ys(ind),zs(ind)          !表面の節点
        dimension xt(ind),yt(ind),zt(ind),dl(ind)  !格子
        dimension x0(ind),y0(ind),z0(ind),d0(ind)  !格子:分割分
        dimension xu(ind),yu(ind),zu(ind),du(ind)  !コピー用
c
        nt=1        !格子の数
        xt(1)=0d0   !格子の座標
        yt(1)=0d0
        zt(1)=0d0
        dl(1)=dl0   !格子の幅
c
300     continue
c
        do 10 kt=1,nt
10      mt(kt)=0
        do 100 kt=1,nt
          xb1=xt(kt)-dl(kt)/2d0     !xb1 <= x <= xb2
          xb2=xt(kt)+dl(kt)/2d0
          yb1=yt(kt)-dl(kt)/2d0     !yb1 <= y <= yb2
          yb2=yt(kt)+dl(kt)/2d0
          zb1=zt(kt)-dl(kt)/2d0     !zb1 <= z <= zb2
          zb2=zt(kt)+dl(kt)/2d0
          do ks=1,ns
            if(xs(ks).gt.xb1-1d-15 .and. xs(ks).lt.xb2+1d-15 .and.
     &         ys(ks).gt.yb1-1d-15 .and. ys(ks).lt.yb2+1d-15 .and.
     &         zs(ks).gt.zb1-1d-15 .and. zs(ks).lt.zb2+1d-15)then
              mt(kt)=mt(kt)+1  !格子内の表面節点の数
            endif
          enddo
100     continue
        do kt=1,nt
          if(mt(kt).ge.2)goto 30
        enddo
        goto 500 !ループ抜ける
c
c       格子を分割する
30      continue
        n0=0
        do 200 kt=1,nt
          if(mt(kt).lt.2)goto 200
          call addoct(xt(kt),yt(kt),zt(kt),dl(kt),n0,x0,y0,z0,d0,ind)
200     continue
c       格子の登録:分割していないもの
        nu=0
        do 400 kt=1,nt
          if(mt(kt).ge.2)goto 400  !分割したものを除く
          nu=nu+1
          xu(nu)=xt(kt)
          yu(nu)=yt(kt)
          zu(nu)=zt(kt)
          du(nu)=dl(kt)
400     continue
c       格子の登録:分割したもの
        do k0=1,n0
          nu=nu+1
          xu(nu)=x0(k0)
          yu(nu)=y0(k0)
          zu(nu)=z0(k0)
          du(nu)=d0(k0)
        enddo
c       移し変え
        nt=0
        do ku=1,nu
          nt=nt+1
          xt(nt)=xu(ku)
          yt(nt)=yu(ku)
          zt(nt)=zu(ku)
          dl(nt)=du(ku)
        enddo
c
        write(*,*)"nt",nt
        goto 300
c
c       表面の節点がある格子を除く
500     continue
        nu=0
        do 600 kt=1,nt
          if(mt(kt).eq.1)goto 600
          nu=nu+1   !表面の節点を含まない格子の数
          xu(nu)=xt(kt)
          yu(nu)=yt(kt)
          zu(nu)=zt(kt)
          du(nu)=dl(kt)
600     continue
        nt=0
        do ku=1,nu
          nt=nt+1
          xt(nt)=xu(ku)
          yt(nt)=yu(ku)
          zt(nt)=zu(ku)
          dl(nt)=du(ku)
        enddo
c
        end
c---------------------------------------------------------------------------
c
c       8分木:格子を分割する
c
c---------------------------------------------------------------------------
        subroutine addoct(xt,yt,zt,dl,n0,x0,y0,z0,d0,ind)
        implicit double precision(a-h,o-z)
        dimension x0(ind),y0(ind),z0(ind),d0(ind)
c
        n0=n0+1
        x0(n0)=xt-dl/4d0  !p1(-dx,-dy,-dz)
        y0(n0)=yt-dl/4d0  !p1(-dx,-dy,-dz)
        z0(n0)=zt-dl/4d0  !p1(-dx,-dy,-dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt+dl/4d0  !p2( dx,-dy,-dz)
        y0(n0)=yt-dl/4d0  !p2( dx,-dy,-dz)
        z0(n0)=zt-dl/4d0  !p2( dx,-dy,-dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt+dl/4d0  !p3( dx,-dy, dz)
        y0(n0)=yt-dl/4d0  !p3( dx,-dy, dz)
        z0(n0)=zt+dl/4d0  !p3( dx,-dy, dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt-dl/4d0  !p4(-dx,-dy, dz)
        y0(n0)=yt-dl/4d0  !p4(-dx,-dy, dz)
        z0(n0)=zt+dl/4d0  !p4(-dx,-dy, dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt-dl/4d0  !p5(-dx, dy,-dz)
        y0(n0)=yt+dl/4d0  !p5(-dx, dy,-dz)
        z0(n0)=zt-dl/4d0  !p5(-dx, dy,-dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt+dl/4d0  !p6( dx, dy,-dz)
        y0(n0)=yt+dl/4d0  !p6( dx, dy,-dz)
        z0(n0)=zt-dl/4d0  !p6( dx, dy,-dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt+dl/4d0  !p7( dx, dy, dz)
        y0(n0)=yt+dl/4d0  !p7( dx, dy, dz)
        z0(n0)=zt+dl/4d0  !p7( dx, dy, dz)
        d0(n0)=dl/2d0

        n0=n0+1
        x0(n0)=xt-dl/4d0  !p8(-dx, dy, dz)
        y0(n0)=yt+dl/4d0  !p8(-dx, dy, dz)
        z0(n0)=zt+dl/4d0  !p8(-dx, dy, dz)
        d0(n0)=dl/2d0
c
        return
        end
c---------------------------------------------------------------------------
c
c       並べ替え : 原点に近い節点 → その節点にに近い節点 → その節点にに近い節点
c
c---------------------------------------------------------------------------
        subroutine sort(nt,xt,yt,zt,mmp,ind)
        implicit double precision(a-h,o-z)
        dimension xt(ind),yt(ind),zt(ind),mt(ind)
        dimension mmp(ind)
c
        x0=0d0
        y0=0d0
        z0=0d0
        do 100 kt=1,nt
100     mt(kt)=0
c
        do mp=1,nt      !順番
          dmin=999d0
          do 200 kt=1,nt    !節点
            if(mt(kt).eq.1)goto 200
            dd=dsqrt((xt(kt)-x0)**2+(yt(kt)-y0)**2+(zt(kt)-z0)**2)
            if(dd.lt.dmin)kmin=kt
            if(dd.lt.dmin)dmin=dd
200       continue
          mt(kmin)=1
          mmp(mp)=kmin
c          x0=xt(kmin)
c          y0=yt(kmin)
c          z0=zt(kmin)
        enddo
c
        return
        end
c---------------------------------------------------------------------------
c
c       並べ替え : それまでの節点の重心に近い節点
c
c---------------------------------------------------------------------------
        subroutine sort2(nt,xt,yt,zt,mmp,ind)
        implicit double precision(a-h,o-z)
        dimension xt(ind),yt(ind),zt(ind),mt(ind)
        dimension mmp(ind)
c
c       初期化
        do 100 kt=1,nt
100     mt(kt)=0
c
c       最初の節点(原点に近い節点)
        dmin=999d0
        do kt=1,nt
          dd=dsqrt((xt(kt))**2+(yt(kt))**2+(zt(kt))**2)
          if(dd.lt.dmin)nmin=kt
          if(dd.lt.dmin)dmin=dd
        enddo
        nn=1    !登録した節点の数
        mmp(nn)=nmin
        mt(nmin)=1  !登録した節点
c
        do mp=2,nt      !順番
c         登録した節点の重心
          xg=0
          yg=0
          zg=0
          do kn=1,nn
            xg=xg+xt(mmp(kn))/dble(nn)
            yg=yg+yt(mmp(kn))/dble(nn)
            zg=zg+zt(mmp(kn))/dble(nn)
          enddo
c         重心に最も近い節点を探す
          dmin=999d0
          do 200 kt=1,nt    !節点
            if(mt(kt).eq.1)goto 200
            dd=dsqrt((xt(kt)-xg)**2+(yt(kt)-yg)**2+(zt(kt)-zg)**2)
            if(dd.lt.dmin)nmin=kt
            if(dd.lt.dmin)dmin=dd
200       continue
          nn=nn+1
          mmp(nn)=nmin
          mt(nmin)=1
        enddo
c
        write(*,*)"sort2,nn,nt",nn,nt
        pause
c
        return
        end
c
c
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c       面を作る
c       表面:nels,ms(nels,3)
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
c



c---------------------------------------------------------------------------
c
c       表面を作る
c
c---------------------------------------------------------------------------
c        subroutine initsuf(ns,xs,ys,zs,nels,ms)
        subroutine initsuf(ns,xs,ys,zs,nels,ms,mgr)
        include "header.h"
        dimension xs(ind),ys(ind),zs(ind)
        dimension ms(ine,3),ms0(ine,3) !ms0(nel0,3)    :今回登録した要素
        dimension le(ine,2),le0(ine,2) !le(ie,2)       :辺
        dimension mn(ind)
        dimension mgr(ind)     !節点候補のグループ 1-n, 0は境界
        data PI/3.1415926535897932384626433832795/
c
c       内側の基準点
        x0=0d0
c        x0=0.75d0
        y0=0d0
        z0=0d0
c
c       最初の面
        call stsuf(x0,y0,z0,ns,xs,ys,zs,n1,n2,n3,ind)
        write(*,*)"n1,n2,n3",n1,n2,n3

        nels=1
        ms(1,1)=n1
        ms(1,2)=n2
        ms(1,3)=n3

        ie=3
        le(1,1)=n2
        le(1,2)=n1
        le(2,1)=n3
        le(2,2)=n2
        le(3,1)=n1
        le(3,2)=n3
c
c       面の法線は内向き

100     continue

        write(*,*)"================================================="

        ie0=0   !今回登録した辺の数
        do 300 ke=1,ie
          k1=le(ke,1)
          k2=le(ke,2)


          write(*,*)"---------------"
          write(*,*)"ke",ke,"/",ie
          write(*,*)"k1,k2",k1,k2
          write(*,*)"mgr",mgr(k1),mgr(k2)
          write(*,*)"ie0",ie0
          write(*,*)"nels",nels
c          do kels=1,nels
c            write(*,*)"ms",kels,":",(ms(kels,j),j=1,3)
c          enddo


c         既に登録した要素
          do kels=1,nels
            if(k1.eq.ms(kels,1) .and. k2.eq.ms(kels,2))goto 300
            if(k1.eq.ms(kels,2) .and. k2.eq.ms(kels,3))goto 300
            if(k1.eq.ms(kels,3) .and. k2.eq.ms(kels,1))goto 300
          enddo
c         要素内部の要素を探す
          call ichkeg2(ns,nels,ms,mn,ind,ine)  !mn(ns) 0:外側の節点,1:要素内部の節点
c         節点を探す
          smax=-999d0
          do 200 ks=1,ns
            if(mn(ks).eq.1)goto 200                                !要素内部の節点
            call edge(k1,k2,ks,xs,ys,zs,sht,ind)                   !角度の計算
            if(sht.le.smax+1d-10)goto 200
            if(isame(k1,k2,ks,nels,ms,ine).eq.1)goto 200           !既にある面は飛ばす
            if(iacute(k1,k2,ks,nels,ms,ns,xs,ys,zs).eq.1)goto 200  !鋭角(acute angle)になる面は飛ばす
            if(inele(k1,k2,ks,nels,ms).eq.1)goto 200               !1辺の要素の数は最大2個
            if((mgr(k1).ne.0 .or.  mgr(k2).ne.0) .and.
     &         (mgr(k1).ne.0 .and. mgr(k1).ne.mgr(ks) .and.
     &          mgr(ks).ne.0))goto 200                             !同じ節点グループではない
            if((mgr(k1).ne.0 .or.  mgr(k2).ne.0) .and.
     &         (mgr(k2).ne.0 .and. mgr(k2).ne.mgr(ks) .and. 
     &          mgr(ks).ne.0))goto 200                             !同じ節点グループではない
            if(icube(k1,k2,ks,xs,ys,zs,smax).eq.1)goto 200              !立方体を点対称に作成する条件
            write(*,*)"ks",ks,dble(sht/PI*180d0),"/",sngl(smax/PI*180d0)
            kmax=ks
            smax=sht
c            call pltsetline(k1,k2,ks,5d0,4)
c            pause
200       continue
          if(smax.eq.-999d0)write(*,*)"節点が見つかりません"
          if(smax.eq.-999d0)stop
          write(*,*)"kmax,ns",kmax,ns
          write(*,*)"mgr:",mgr(kmax)
c         要素を登録
          nels=nels+1
          ms(nels,1)=k1
          ms(nels,2)=k2
          ms(nels,3)=kmax
c         次の辺を登録
          ie0=ie0+1
          le0(ie0,1)=kmax
          le0(ie0,2)=k2
          ie0=ie0+1
          le0(ie0,1)=k1
          le0(ie0,2)=kmax
          call pltsetme(nels,ms,3d0,2,ind,ine)
c          pause
300     continue
        if(ie0.eq.0)return     !辺がない
c
        ie=ie0
        do 500 ke=1,ie0
        do 500 j=1,2
500     le(ke,j)=le0(ke,j)
c
        do 1 kels=1,nels
1       write(*,*)"kels",kels,":",ms(kels,1),ms(kels,2),ms(kels,3)

        goto 100
c
        end
c---------------------------------------------------------------------------
c
c       角度を計算する
c
c--------------------------------------------------------------------------
        subroutine edge(k1,k2,ks,xs,ys,zs,sht,ind)
        implicit double precision(a-h,o-z)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        v1x=xs(k1)-xs(ks)
        v1y=ys(k1)-ys(ks)
        v1z=zs(k1)-zs(ks)

        v2x=xs(k2)-xs(ks)
        v2y=ys(k2)-ys(ks)
        v2z=zs(k2)-zs(ks)

        sht=finpr(v1x,v1y,v1z,v2x,v2y,v2z)
c        sht=PI-sht

        return
        end
c---------------------------------------------------------------------------
c
c       最初の面
c       p0(x0,y0,z0)要素内部の節点    面は内向き
c---------------------------------------------------------------------------
        subroutine stsuf(x0,y0,z0,ns,xs,ys,zs,n1,n2,n3,ind)
        implicit double precision(a-h,o-z)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        dmin=999d0
        do ks=1,ns
          dd=dsqrt((xs(ks)-x0)**2+(ys(ks)-y0)**2+(zs(ks)-z0)**2)
          if(dd.lt.dmin)kmin=ks
          if(dd.lt.dmin)dmin=dd
        enddo
        n1=kmin
c
        dmin=999d0
        do 10 ks=1,ns
          if(ks.eq.n1)goto 10
          dd=dsqrt((xs(ks)-xs(n1))**2+(ys(ks)-ys(n1))**2
     &                               +(zs(ks)-zs(n1))**2)
          if(dd.lt.dmin)kmin=ks
          if(dd.lt.dmin)dmin=dd
10      continue
        n2=kmin
c
        dmin=999d0
        do 20 ks=1,ns
          if(ks.eq.n1 .or. ks.eq.n2)goto 20
          if(artri(xs(n1),ys(n1),zs(n1),xs(n2),ys(n2),zs(n2)
     &                      ,xs(ks),ys(ks),zs(ks)).eq.0d0)goto 20
          dd=dsqrt((xs(ks)-xs(n2))**2+(ys(ks)-ys(n2))**2
     &                               +(zs(ks)-zs(n2))**2)
          if(dd.lt.dmin)kmin=ks
          if(dd.lt.dmin)dmin=dd
20      continue
        n3=kmin
c
c       面の法線との角度
        call chkne(x0,y0,z0,xs(n1),ys(n1),zs(n1),xs(n2),ys(n2),zs(n2)
     &                      ,xs(n3),ys(n3),zs(n3),shita)
        if( shita .gt. (89d0/180d0*PI) )then   !節点入れ換え
          nn=n2
          n2=n3
          n3=nn
          call chkne(x0,y0,z0,xs(n1),ys(n1),zs(n1),xs(n2),ys(n2),zs(n2)
     &                        ,xs(n3),ys(n3),zs(n3),shita)
          write(*,*)"stsuf,shita",shita/PI*180d0
          pause
          return
        endif
        write(*,*)"stsuf,shita",shita/PI*180d0
        pause
c
        return
        end
c---------------------------------------------------------------------------
c
c       要素の中の節点は選ばない
c       作った面の内部の節点
c
c---------------------------------------------------------------------------
        subroutine ichkeg2(ns,nels,ms,mn,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3)
        dimension nec(ine,3),mn(ind)
        dimension js(4)
        data js/1,2,3,1/
c
c       初期化
        do 100 i=1,ns
100     mn(i)=0         !0:節点候補,1:作った面内部の要素
c
c       面の節点を探す
        do 200 kels=1,nels
        do 200 j=1,3
          mn(ms(kels,j))=1      !1:作った面内部の要素(後で縁の要素除く)
200     continue
c
c       nec(nelt,3)を探す(隣接する要素)
        do 300 kels=1,nels
        do 300 j=1,3
300     nec(kels,j)=0
c
        do 400 i=1,nels-1
        do 400 j=1,3    !辺の数
          i1=ms(i,js(j))
          i2=ms(i,js(j+1))
          do 500 k=i+1,nels
          do 500 l=1,3  !辺の数
            j1=ms(k,js(l))
            j2=ms(k,js(l+1))
            if(i1.eq.j2 .and. i2.eq.j1)then    !辺を共有している
              nec(i,j)=k
              nec(k,l)=i
            endif
500       continue
400     continue
c
c       内部の節点を探す
c       辺1の節点番号:1-2
c       辺2の節点番号:2-3
c       辺3の節点番号:3-1
        do 600 ks=1,ns
600     mn(ks)=0
        do 700 kels=1,nels
        do 700 j=1,3    !辺の番号
          if(nec(kels,j).eq.0)then !外側の要素
            mn(ms(kels,js(j)))=0
            mn(ms(kels,js(j+1)))=0
          endif
700     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       既にある面は飛ばす
c
c---------------------------------------------------------------------------
        function isame(k1,k2,ks,nels,ms,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3)
c
        do k=1,nels
          if((k1.eq.ms(k,1) .or. k1.eq.ms(k,2) .or. k1.eq.ms(k,3)) .and.
     &       (k2.eq.ms(k,1) .or. k2.eq.ms(k,2) .or. k2.eq.ms(k,3)) .and.
     &       (ks.eq.ms(k,1) .or. ks.eq.ms(k,2) .or. ks.eq.ms(k,3)) )then
            isame=1
            return
          endif
        enddo
        isame=0
        return
        end
c---------------------------------------------------------------------------
c
c       鋭角(acute angle)になる面は飛ばす
c         面が折れ曲がる場合
c         辺k1-k2を共有する面との角度.k1,k2の中心を基準とした角度
c         面が合わさっていても横向きの角度で0°にならない
c---------------------------------------------------------------------------
        function iacute(k1,k2,ks,nels,ms,ns,xs,ys,zs)
        include "header.h"
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        do kels=1,nels
          if(ms(kels,1).eq.k2 .and. ms(kels,2).eq.k1)then
            kt=ms(kels,3)
            goto 10
          endif
          if(ms(kels,2).eq.k2 .and. ms(kels,3).eq.k1)then
            kt=ms(kels,1)
            goto 10
          endif
          if(ms(kels,3).eq.k2 .and. ms(kels,1).eq.k1)then
            kt=ms(kels,2)
            goto 10
          endif
        enddo
c
10      continue
        x1=xs(k1)
        y1=ys(k1)
        z1=zs(k1)
        x2=xs(k2)
        y2=ys(k2)
        z2=zs(k2)
        xt=xs(kt)
        yt=ys(kt)
        zt=zs(kt)
        xr=xs(ks)
        yr=ys(ks)
        zr=zs(ks)
        xg=(xs(k1)+xs(k2))/2d0
        yg=(ys(k1)+ys(k2))/2d0
        zg=(zs(k1)+zs(k2))/2d0
c
        ck1=-(xt-xg)*(x2-x1)-(yt-yg)*(y2-y1)-(zt-zg)*(z2-z1)
        ck2=(x2-x1)**2+(y2-y1)**2+(z2-z1)**2
        ck=ck1/ck2
        cl1=-(xr-xg)*(x2-x1)-(yr-yg)*(y2-y1)-(zr-zg)*(z2-z1)
        cl2=(x2-x1)**2+(y2-y1)**2+(z2-z1)**2
        cl=cl1/cl2
c
c
        v1x=xt+ck*(x2-x1)-xg
        v1y=yt+ck*(y2-y1)-yg
        v1z=zt+ck*(z2-z1)-zg
        v2x=xr+cl*(x2-x1)-xg
        v2y=yr+cl*(y2-y1)-yg
        v2z=zr+cl*(z2-z1)-zg
c
        sht=finpr(v1x,v1y,v1z,v2x,v2y,v2z)        
c        write(*,*)"sht:",sngl(sht/PI*180d0)
c        write(*,*)"ck,cl",sngl(ck),sngl(cl)
c        write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
c        write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
c        write(*,*)"xr:",sngl(xr),sngl(yr),sngl(zr)
c        write(*,*)"xt:",sngl(xt),sngl(yt),sngl(zt)
c        write(*,*)"v1:",sngl(v1x),sngl(v1y),sngl(v1z)
c        write(*,*)"v2:",sngl(v2x),sngl(v2y),sngl(v2z)
c
c       面の角度が<10°の時は飛ばす
        iacute=0
        if(sht.lt.10d0/180d0*PI)iacute=1
c
        return
        end
c---------------------------------------------------------------------------
c
c       1辺の要素の数を探す.最大2個.
c
c---------------------------------------------------------------------------
        function inele(k1,k2,ks,nels,ms)
        include "header.h"
        dimension ms(ine,3)
c
        ifg=0
c       違う向きでは2つまで
c       辺k1-ksが存在する要素
        n=0
        do kels=1,nels
          if(ms(kels,1).eq.k1 .and. ms(kels,2).eq.ks)n=n+1
          if(ms(kels,2).eq.k1 .and. ms(kels,3).eq.ks)n=n+1
          if(ms(kels,3).eq.k1 .and. ms(kels,1).eq.ks)n=n+1
          if(ms(kels,1).eq.ks .and. ms(kels,2).eq.k1)n=n+1
          if(ms(kels,2).eq.ks .and. ms(kels,3).eq.k1)n=n+1
          if(ms(kels,3).eq.ks .and. ms(kels,1).eq.k1)n=n+1
        enddo
        if(n.ge.2)ifg=1
c       辺k2-ksが存在する要素
        n=0
        do kels=1,nels
          if(ms(kels,1).eq.k2 .and. ms(kels,2).eq.ks)n=n+1
          if(ms(kels,2).eq.k2 .and. ms(kels,3).eq.ks)n=n+1
          if(ms(kels,3).eq.k2 .and. ms(kels,1).eq.ks)n=n+1
          if(ms(kels,1).eq.ks .and. ms(kels,2).eq.k2)n=n+1
          if(ms(kels,2).eq.ks .and. ms(kels,3).eq.k2)n=n+1
          if(ms(kels,3).eq.ks .and. ms(kels,1).eq.k2)n=n+1
        enddo
        if(n.ge.2)ifg=1
c
        inele=ifg
        return
        end
c---------------------------------------------------------------------------
c
c       立方体を点対称に作成する条件
c
c---------------------------------------------------------------------------
        function icube(k1,k2,ks,xs,ys,zs,smax)
        include "header.h"
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.1415926535897932384626433832795/
c
        if((xs(k1).eq.0d0 .and. xs(k2).eq.0d0 .and. xs(ks).eq.0d0).or.
     &     (xs(k1).eq.2d0 .and. xs(k2).eq.2d0 .and. xs(ks).eq.2d0))then
          xc1=zs(k1)
          xc2=zs(k2)
          xc3=zs(ks)
          yc1=ys(k1)
          yc2=ys(k2)
          yc3=ys(ks)
      elseif((ys(k1).eq.0d0 .and. ys(k2).eq.0d0 .and. ys(ks).eq.0d0).or.
     &     (ys(k1).eq.2d0 .and. ys(k2).eq.2d0 .and. ys(ks).eq.2d0))then
          xc1=xs(k1)
          xc2=xs(k2)
          xc3=xs(ks)
          yc1=zs(k1)
          yc2=zs(k2)
          yc3=zs(ks)
      elseif((zs(k1).eq.0d0 .and. zs(k2).eq.0d0 .and. zs(ks).eq.0d0).or.
     &     (zs(k1).eq.2d0 .and. zs(k2).eq.2d0 .and. zs(ks).eq.2d0))then
          xc1=xs(k1)
          xc2=xs(k2)
          xc3=xs(ks)
          yc1=ys(k1)
          yc2=ys(k2)
          yc3=ys(ks)
        else
c          write(*,*)"in icube no surface was found."
c          write(*,*)"x1:",dble(xs(k1)),dble(ys(k1)),dble(zs(k1))
c          write(*,*)"x2:",dble(xs(k2)),dble(ys(k2)),dble(zs(k2))
c          write(*,*)"xs:",dble(xs(ks)),dble(ys(ks)),dble(zs(ks))
          icube=1
          return
        endif
c        call edge(k1,k2,ks,xs,ys,zs,sht,ind)                   !角度の計算
c        write(*,*)"sht:",sngl(sht/PI*180d0),"/",sngl(smax/PI*180d0)
c        write(*,*)"xc1:",sngl(xc1),sngl(yc1)
c        write(*,*)"xc2:",sngl(xc2),sngl(yc2)
c        write(*,*)"xc3:",sngl(xc3),sngl(yc3)
c        write(*,*)"x1:",dble(xs(k1)),dble(ys(k1)),dble(zs(k1))
c        write(*,*)"x2:",dble(xs(k2)),dble(ys(k2)),dble(zs(k2))
c        write(*,*)"xs:",dble(xs(ks)),dble(ys(ks)),dble(zs(ks))
c
c         右
          ifg=1
          x1=1d0
          y1=1d0
          x2=2d0
          y2=0d0
          x3=2d0
          y3=2d0
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc1,yc1,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc1,yc1,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc1,yc1,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn1:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 100
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc2,yc2,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc2,yc2,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc2,yc2,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn1:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 100
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc3,yc3,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc3,yc3,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc3,yc3,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn1:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 100
          goto 400
c
c         上
100       continue
          ifg=2
          x1=1d0
          y1=1d0
          x2=2d0
          y2=2d0
          x3=0d0
          y3=2d0
          ar=artri(x21,y21,0d0,x22,y22,0d0,x23,y23,0d0)
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc1,yc1,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc1,yc1,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc1,yc1,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn2:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 200
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc2,yc2,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc2,yc2,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc2,yc2,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn2:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 200
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc3,yc3,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc3,yc3,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc3,yc3,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn2:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 200
          goto 400
c
c         左
200       continue
          ifg=3
          x1=1d0
          y1=1d0
          x2=0d0
          y2=2d0
          x3=0d0
          y3=0d0
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc1,yc1,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc1,yc1,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc1,yc1,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn3:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 300
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc2,yc2,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc2,yc2,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc2,yc2,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn3:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 300
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc3,yc3,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc3,yc3,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc3,yc3,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn3:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 300
          goto 400
c
c         下
300       continue
          ifg=4
          x1=1d0
          y1=1d0
          x2=0d0
          y2=0d0
          x3=2d0
          y3=0d0
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc1,yc1,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc1,yc1,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc1,yc1,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn4:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 500
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc2,yc2,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc2,yc2,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc2,yc2,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn4:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 500
          ar =artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)
          dn1=artri(xc3,yc3,0d0,x2 ,y2 ,0d0,x3 ,y3 ,0d0)/ar
          dn2=artri(x1 ,y1 ,0d0,xc3,yc3,0d0,x3 ,y3 ,0d0)/ar
          dn3=artri(x1 ,y1 ,0d0,x2 ,y2 ,0d0,xc3,yc3,0d0)/ar
          dn=dn1+dn2+dn3
c          write(*,*)"dn4:",sngl(dn)
          if(dn.gt.1d0+1d-10)goto 500
          goto 400
c
c         飛ばす
500       continue
          icube=1
          return
c
c         判定する
400       continue
          write(*,*)"ifg:",ifg
c          write(*,*)"xc1:",sngl(xc1),sngl(yc1)
c          write(*,*)"xc2:",sngl(xc2),sngl(yc2)
c          write(*,*)"xc3:",sngl(xc3),sngl(yc3)
          icube=1
          if(ifg.eq.1 .or. ifg.eq.3)then   !縦
            if(xc1.eq.xc2 .or. xc2.eq.xc3 .or. xc3.eq.xc1)icube=0
          end if
          if(ifg.eq.2 .or. ifg.eq.4)then   !横
            if(yc1.eq.yc2 .or. yc2.eq.yc3 .or. yc3.eq.yc1)icube=0
          end if
          write(*,*)"icube:",icube
c          pause
c
        return
        end
c---------------------------------------------------------------------------
c
c       格子のバランスをとる
c
c---------------------------------------------------------------------------


c---------------------------------------------------------------------------
c
c       外側の格子を削除する
c         内側の節点の任意の方向に面が存在する
c         ある方向に存在する面の最も近いの面で
c         線分(節点-面の重心)と内向き法線の角度は<=90°
c         線分(節点-面の重心)と他の面が交差していないか調べる
c---------------------------------------------------------------------------
        subroutine delouternp(ns,nels,xs,ys,zs,ms,nt,xt,yt,zt,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3),mt(ind)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xt(ind),yt(ind),zt(ind)
        dimension xu(ind),yu(ind),zu(ind)
        data PI/3.1415926535897932384626433832795/
c
        do kt=1,nt
          x0=xt(kt)
          y0=yt(kt)
          z0=zt(kt)
c
c         最も近い面を探す
          nmin=-999
          dmin=999d0
          do kels=1,nels
            xg=(xs(ms(kels,1))+xs(ms(kels,2))+xs(ms(kels,3)))/3d0
            yg=(ys(ms(kels,1))+ys(ms(kels,2))+ys(ms(kels,3)))/3d0
            zg=(zs(ms(kels,1))+zs(ms(kels,2))+zs(ms(kels,3)))/3d0
            dd=dsqrt((x0-xg)**2+(y0-yg)**2+(z0-zg)**2)
            if(dd.lt.dmin)nmin=kels
            if(dd.lt.dmin)dmin=dd
          enddo
c
c         外積:面は内向きなので内向きの法線
          v1x=xs(ms(nmin,2))-xs(ms(nmin,1))
          v1y=ys(ms(nmin,2))-ys(ms(nmin,1))
          v1z=zs(ms(nmin,2))-zs(ms(nmin,1))
          v2x=xs(ms(nmin,3))-xs(ms(nmin,1))
          v2y=ys(ms(nmin,3))-ys(ms(nmin,1))
          v2z=zs(ms(nmin,3))-zs(ms(nmin,1))
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
c
c         線分:面の重心→節点
          xg=(xs(ms(nmin,1))+xs(ms(nmin,2))+xs(ms(nmin,3)))/3d0
          yg=(ys(ms(nmin,1))+ys(ms(nmin,2))+ys(ms(nmin,3)))/3d0
          zg=(zs(ms(nmin,1))+zs(ms(nmin,2))+zs(ms(nmin,3)))/3d0
          vx=x0-xg
          vy=y0-yg
          vz=z0-zg
c
c         角度
          mt(kt)=0
          sht=finpr(dnx,dny,dnz,vx,vy,vz)
          if(sht.lt.90d0/180d0*PI)mt(kt)=1
        enddo
c
c       表面に近い節点を除く(面ができなくなる)
c
        do 10 kt=1,nt
          x=xt(kt)
          y=yt(kt)
          z=zt(kt)
          do ks=1,ns
            dd=dsqrt((xs(ks)-x)**2+(ys(ks)-y)**2+(zs(ks)-z)**2)
            if(dd.lt.1d-1)mt(kt)=1
            if(dd.lt.1d-1)goto 10
          enddo
10      continue
c
c       節点を並べなおす
        nu=0
        do 100 kt=1,nt
          if(mt(kt).eq.1)goto 100
          nu=nu+1
          xu(nu)=xt(kt)
          yu(nu)=yt(kt)
          zu(nu)=zt(kt)
100     continue
c
c       移し変え
        write(*,*)"nt,nu",nt,nu
        nt=nu
        do kt=1,nt
          xt(kt)=xu(kt)
          yt(kt)=yu(kt)
          zt(kt)=zu(kt)
        enddo
      return
      end
