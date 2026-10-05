
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
        subroutine octtree(ns,xs,ys,zs,nt,xt,yt,zt,dl0,ng,ind)
        implicit double precision(a-h,o-z)
        dimension mt(ind)
        dimension xs(ind),ys(ind),zs(ind)          !表面の節点
        dimension xt(ind),yt(ind),zt(ind),dl(ind)  !格子
        dimension x0(ind),y0(ind),z0(ind),d0(ind)  !格子:分割分
        dimension xu(ind),yu(ind),zu(ind),du(ind)  !コピー用
c*******************************************************************
c       必要以上に分割しない
c*******************************************************************
c        nt=1        !格子の数
c        xt(1)=0d0   !格子の座標
c        yt(1)=0d0
c        zt(1)=0d0
c        dl(1)=dl0   !格子の幅
c*******************************************************************
c       任意の幅に分割
c*******************************************************************
        nt=0
        do 700 i=0,ng   !ng:分割数
        do 700 j=0,ng
        do 700 k=0,ng
          nt=nt+1
          xt(nt)=dl0/dble(ng)*dble(i)
          yt(nt)=dl0/dble(ng)*dble(j)
          zt(nt)=dl0/dble(ng)*dble(k)
          dl(nt)=dl0/dble(ng)
700     continue
c*******************************************************************
c
c       表面の節点を2つ以上含む格子を分割
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
        subroutine initsuf(ns,xs,ys,zs,nels,ms)
        include "header.h"
        dimension xs(ind),ys(ind),zs(ind)
        dimension ms(ine,3),ms0(ine,3) !ms0(nel0,3)    :今回登録した要素
        dimension le(ine,2),le0(ine,2) !le(ie,2)       :辺
        dimension mn(ind)
        data PI/3.1415926535897932384626433832795/
c
c       内側の基準点
c        x0=0d0
c        x0=0.75d0
c        y0=0d0
c        z0=0d0
c       外側の基準点
        x0=-99d0
        y0=-99d0
        z0=-99d0
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

        ie0=0   !今回登録した辺の数
        do 300 ke=1,ie
          k1=le(ke,1)
          k2=le(ke,2)
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
            if(isame(k1,k2,ks,nels,ms,ine).eq.1)goto 200           !既にある面は飛ばす
            if(iacute(k1,k2,ks,nels,ms,ns,xs,ys,zs).eq.1)goto 200  !鋭角(acute angle)になる面は飛ばす
            if(inele(k1,k2,ks,nels,ms).eq.1)goto 200               !1辺の要素の数は最大2個
            call edge(k1,k2,ks,xs,ys,zs,sht,ind)                   !角度の計算
            if(sht.gt.smax)then
              kmax=ks
              smax=sht
c              call pltsetline(k1,k2,ks,5d0,4)
c              write(*,*)"ks",ks,sngl(sht/PI*180d0),"/",sngl(smax/PI*180d0)
c              pause
            endif
200       continue
          if(smax.eq.-999d0)write(*,*)"節点が見つかりません,initsuf"
          if(smax.eq.-999d0)stop
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
          pause
300     continue
        if(ie0.eq.0)return     !辺がない
c
        ie=ie0
        do 500 ke=1,ie0
        do 500 j=1,2
500     le(ke,j)=le0(ke,j)
c
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
        if( shita .le. (89d0/180d0*PI) )then   !節点入れ換え
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
c
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
        x1=xs(kt)
        y1=ys(kt)
        z1=zs(kt)
        x2=xs(ks)
        y2=ys(ks)
        z2=zs(ks)
        x0=(xs(k1)+xs(k2))/2d0
        y0=(ys(k1)+ys(k2))/2d0
        z0=(zs(k1)+zs(k2))/2d0
c
        v1x=x1-x0
        v1y=y1-y0
        v1z=z1-z0
        v2x=x2-x0
        v2y=y2-y0
        v2z=z2-z0
c
        sht=finpr(v1x,v1y,v1z,v2x,v2y,v2z)        
c
c       面の角度が<90°の時は飛ばす
        iacute=0
        if(sht.lt.90d0/180d0*PI)iacute=1
c
        return
        end
c---------------------------------------------------------------------------
c
c       1辺の要素の数を探す.最大2個
c
c---------------------------------------------------------------------------
        function inele(k1,k2,ks,nels,ms)
        include "header.h"
        dimension ms(ine,3)
c
        ifg=0
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
c        do 10 kt=1,nt
c          x=xt(kt)
c          y=yt(kt)
c          z=zt(kt)
c          do ks=1,ns
c            dd=dsqrt((xs(ks)-x)**2+(ys(ks)-y)**2+(zs(ks)-z)**2)
c            if(dd.lt.1d-1)mt(kt)=1
c            if(dd.lt.1d-1)goto 10
c          enddo
c10      continue
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
c---------------------------------------------------------------------------
c
c       面の節点の補間と削除
c       表面を均一にする
c---------------------------------------------------------------------------
        subroutine resuf(ns,xs,ys,zs,nels,ms,dl)
        include "header.h"
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind),r(3)
        data PI/3.1415926535897932384626433832795/
c
1000    continue
        do 100 i=1,nels
c
c         節点の座標
          x1=xs(ms(i,1))
          x2=xs(ms(i,2))
          x3=xs(ms(i,3))
          y1=ys(ms(i,1))
          y2=ys(ms(i,2))
          y3=ys(ms(i,3))
          z1=zs(ms(i,1))
          z2=zs(ms(i,2))
          z3=zs(ms(i,3))
c
c         面の形状
          dr=difdr(x1,y1,z1,x2,y2,z2,x3,y3,z3)    !面の形状係数
          if(dr.lt.0.5d0)goto 100
c
c         辺の長さを算出
          dl1=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
          dl2=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
          dl3=dsqrt((x1-z2)**2+(y1-y2)**2+(z1-z2)**2)
          if(dl1.ge.2d0*dl .or. dl2.ge.2d0*dl .or. dl3.ge.2d0*dl)then

            write(*,*)"befoer div"
            call pltsufms(i,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            pause
            call divsuf(i,ns,xs,ys,zs,nels,ms,dl)  !面の分割
            write(*,*)"after div"
            call pltsufms(nels,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            pause

            goto 1000
          endif
c
c         θi(i=1,2,3)の算出
          v1x=x2-x1
          v1y=y2-y1
          v1z=z2-z1
          v2x=x3-x2
          v2y=y3-y2
          v2z=z3-z2
          v3x=x1-x3
          v3y=y1-y3
          v3z=z1-z3
          sht10=finpr(v1x,v1y,v1z,-v3x,-v3y,-v3z)
          sht20=finpr(v2x,v2y,v2z,-v1x,-v1y,-v1z)
          sht30=finpr(v3x,v3y,v3z,-v2x,-v2y,-v2z)
c
c         角度を並べ替える
          sht1=sht10
          if(sht1.lt.sht20)sht1=sht20
          if(sht1.lt.sht30)sht1=sht30
          sht3=sht10
          if(sht3.gt.sht20)sht3=sht20
          if(sht3.gt.sht30)sht3=sht30
          sht2=sht10
          if(sht2.eq.sht30 .or. sht2.eq.sht10)sht2=sht20
          if(sht2.eq.sht10. or. sht2.eq.sht20)sht2=sht30
c
c         節点を探す
          if(sht1.eq.sht10)kmax=1
          if(sht1.eq.sht20)kmax=2
          if(sht1.eq.sht30)kmax=3
          if(sht3.eq.sht10)kmin=1
          if(sht3.eq.sht20)kmin=2
          if(sht3.eq.sht30)kmin=3
c
          write(*,*)"i:",i,"/",nels
          write(*,*)"sht:",sngl(sht1/PI*180d0),sngl(sht2/PI*180d0)
     &                                         ,sngl(sht3/PI*180d0)
          pause
c
c         面を削除する
c         θ1>120∩θ2<30∩θ<30
          if(sht1.gt.120d0/180d0*PI .and. sht2.lt.30d0/180d0*PI
     &                              .and. sht3.lt.30d0/180d0*PI)then

            call pltsufms(i,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            write(*,*)"befoer del1"
            pause

            call delsuf1(i,kmax,ns,xs,ys,zs,nels,ms)

            call pltsufms(nels,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            write(*,*)"after del1"
            pause

            goto 1000
          endif
c         θ1>70∩θ2>70∩θ<30
          if(sht1.gt.70d0/180d0*PI .and. sht2.gt.70d0/180d0*PI
     &                              .and. sht3.lt.30d0/180d0*PI)then

            call pltsufms(i,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            write(*,*)"befoer del2"
            pause

            call delsuf2(i,kmin,ns,xs,ys,zs,nels,ms)

            call pltsufms(i,ns,xs,ys,zs,nels,ms)   !面iと節点を共有する面を描く
            write(*,*)"after del2"
            pause

            goto 1000
          endif
100     continue
c
        return
        end
c
c       面の分割
c
        subroutine divsuf(kels,ns,xs,ys,zs,nels,ms,dl)
        include "header.h"
        dimension ms(ine,3),ms0(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
c
        x1=xs(ms(kels,1))
        x2=xs(ms(kels,2))
        x3=xs(ms(kels,3))
        y1=ys(ms(kels,1))
        y2=ys(ms(kels,2))
        y3=ys(ms(kels,3))
        z1=zs(ms(kels,1))
        z2=zs(ms(kels,2))
        z3=zs(ms(kels,3))
        dl1=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        dl2=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
        dl3=dsqrt((x1-z2)**2+(y1-y2)**2+(z1-z2)**2)
        if(dl1.ge.2d0*dl)then
          l1=ms(kels,3)    !共有する辺の節点
          l2=ms(kels,2)    !共有する辺の節点
          mm=ms(kels,1)    !辺以外の節点
        elseif(dl2.ge.2d0*dl)then
          l1=ms(kels,1)    !共有する辺の節点
          l2=ms(kels,3)    !共有する辺の節点
          mm=ms(kels,2)    !辺以外の節点
        elseif(dl3.ge.2d0*dl)then
          l1=ms(kels,2)    !共有する辺の節点
          l2=ms(kels,1)    !共有する辺の節点
          mm=ms(kels,3)    !辺以外の節点
        endif
c
c       辺を共有する面を探す
        do 100 i=1,nels
          if(i.eq.kels)goto 100
          if(ms(i,1).eq.l1 .and. ms(i,2).eq.l2)jels=i     !辺1と辺を共有する面
          if(ms(i,2).eq.l1 .and. ms(i,3).eq.l2)jels=i
          if(ms(i,3).eq.l1 .and. ms(i,1).eq.l2)jels=i
          if(ms(i,1).eq.l1 .and. ms(i,2).eq.l2)nn=ms(i,3) !向かい側の節点
          if(ms(i,2).eq.l1 .and. ms(i,3).eq.l2)nn=ms(i,1)
          if(ms(i,3).eq.l1 .and. ms(i,1).eq.l2)nn=ms(i,2)
100     continue
c
c       節点を追加する
        ns=ns+1
        xs(ns)=(xs(l1)+xs(l2))/2d0
        ys(ns)=(ys(l1)+ys(l2))/2d0
        zs(ns)=(zs(l1)+zs(l2))/2d0
c
c       要素を追加する
        nels=nels+1
        ms(nels,1)=mm
        ms(nels,2)=ns
        ms(nels,3)=l1
        nels=nels+1
        ms(nels,1)=mm
        ms(nels,2)=l2
        ms(nels,3)=ns
        nels=nels+1
        ms(nels,1)=ns
        ms(nels,2)=l2
        ms(nels,3)=nn
        nels=nels+1
        ms(nels,1)=ns
        ms(nels,2)=nn
        ms(nels,3)=l1
c
c       要素を削除する
        nel0=0
        do 200 i=1,nels
          if(i.eq.kels .or. i.eq.jels)goto 200
          nel0=nel0+1
          do 300 j=1,3
300       ms0(nel0,j)=ms(i,j)
200     continue
c
        nels=nel0
        do 400 i=1,nels
        do 400 j=1,3
400     ms(i,j)=ms0(i,j)
c
        return
        end
c
c       面の削除
c
        subroutine delsuf1(kels,kmax,ns,xs,ys,zs,nels,ms)
        include "header.h"
        dimension ms(ine,3),ms0(ind,3)
        dimension xs(ind),ys(ind),zs(ind)
c
        if(kmax.eq.1)then
          l1=ms(kels,3)    !共有する辺の節点
          l2=ms(kels,2)    !共有する辺の節点
          mm=ms(kels,1)    !辺以外の節点
        elseif(kmax.eq.2)then
          l1=ms(kels,1)    !共有する辺の節点
          l2=ms(kels,3)    !共有する辺の節点
          mm=ms(kels,2)    !辺以外の節点
        elseif(kmax.eq.3)then
          l1=ms(kels,2)    !共有する辺の節点
          l2=ms(kels,1)    !共有する辺の節点
          mm=ms(kels,3)    !辺以外の節点
        endif
c
c       辺を共有する面を探す
        do 100 i=1,nels
          if(i.eq.kels)goto 100
          if(ms(i,1).eq.l1 .and. ms(i,2).eq.l2)jels=i     !辺1と辺を共有する面
          if(ms(i,2).eq.l1 .and. ms(i,3).eq.l2)jels=i
          if(ms(i,3).eq.l1 .and. ms(i,1).eq.l2)jels=i
          if(ms(i,1).eq.l1 .and. ms(i,2).eq.l2)nn=ms(i,3) !向かい側の節点
          if(ms(i,2).eq.l1 .and. ms(i,3).eq.l2)nn=ms(i,1)
          if(ms(i,3).eq.l1 .and. ms(i,1).eq.l2)nn=ms(i,2)
100     continue
c
c       節点の位置を変える
        x=(xs(l1)+xs(l2))/2d0
        y=(ys(l1)+ys(l2))/2d0
        z=(zs(l1)+zs(l2))/2d0
        xs(mm)=x
        ys(mm)=y
        zs(mm)=z
c
c       要素を追加する
        nels=nels+1
        ms(nels,1)=mm
        ms(nels,2)=l2
        ms(nels,3)=nn
        nels=nels+1
        ms(nels,1)=mm
        ms(nels,2)=nn
        ms(nels,3)=l1
c
c       要素を削除する
        nel0=0
        do 200 i=1,nels
          if(i.eq.kels .or. i.eq.jels)goto 200
          nel0=nel0+1
          do 300 j=1,3
300       ms0(nel0,j)=ms(i,j)
200     continue
c
        nels=nel0
        do 400 i=1,nels
        do 400 j=1,3
400     ms(i,j)=ms0(i,j)
c
        return
        end
c
c       面の削除
c
        subroutine delsuf2(kels,kmin,ns,xs,ys,zs,nels,ms)
        include "header.h"
        dimension ms(ine,3),ms0(ind,3),mn(ind)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xs0(ind),ys0(ind),zs0(ind)
c
        if(kmin.eq.1)then
          l1=ms(kels,3)    !共有する辺の節点
          l2=ms(kels,2)    !共有する辺の節点
        elseif(kmin.eq.2)then
          l1=ms(kels,1)    !共有する辺の節点
          l2=ms(kels,3)    !共有する辺の節点
        elseif(kmin.eq.3)then
          l1=ms(kels,2)    !共有する辺の節点
          l2=ms(kels,1)    !共有する辺の節点
        endif
c
c       辺を共有する面を探す
        do 100 i=1,nels
          if(i.eq.kels)goto 100
          if(ms(i,1).eq.l1 .and. ms(i,2).eq.l2)jels=i     !辺1と辺を共有する面
          if(ms(i,2).eq.l1 .and. ms(i,3).eq.l2)jels=i
          if(ms(i,3).eq.l1 .and. ms(i,1).eq.l2)jels=i
100     continue
c
c       要素を削除する
        nel0=0
        do 200 i=1,nels
          if(i.eq.kels .or. i.eq.jels)goto 200
          nel0=nel0+1
          do 300 j=1,3
300       ms0(nel0,j)=ms(i,j)
200     continue
c
        nels=nel0
        do 400 i=1,nels
        do 400 j=1,3
400     ms(i,j)=ms0(i,j)
c
c       節点l2の位置を変える
        x=(xs(l1)+xs(l2))/2d0
        y=(ys(l1)+ys(l2))/2d0
        z=(zs(l1)+zs(l2))/2d0
        xs(l2)=x
        ys(l2)=y
        zs(l2)=z
c
c       節点並べ替え
        ns0=0
        do 500 ks=1,ns
          if(ks.eq.l1)goto 500
          ns0=ns0+1
          xs0(ns0)=xs(ks)
          ys0(ns0)=ys(ks)
          zs0(ns0)=zs(ks)
          mn(ks)=ns0
500     continue
c
c       節点コピー
        ns=ns0
        do 600 ks=1,ns
          xs(ks)=xs0(ks)
          ys(ks)=ys0(ks)
          zs(ks)=zs0(ks)
600     continue
c
c       要素の節点並べ替え
        do 700 i=1,nels
        do 700 j=1,3
          if(ms(i,j).eq.l1)ms(i,j)=l2   !節点l1を節点l2に変える
          ms(i,j)=mn(ms(i,j))
700     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       面の節点の位置を１回最適化する
c
c---------------------------------------------------------------------------
        subroutine nsoptimize(ns,nels,xs,ys,zs,ms)
        include "header.h"
        dimension ms(ine,3)
        dimension ms0(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        dimension nelb(ine),neb(ine)
        data ceft/1d0/    !緩和係数
        data drstd/5d-1/   !基準値
c 
        do 100 ks=1,ns
          call nssur(ks,nels,ms,nb,nelb,neb)   !ksを含む面 ks=ms(nelb(nb),neb(nb))
          drmax=0d0
          drav=0d0
          do kb=1,nb
            x1=xs(ms(nelb(kb),1))
            x2=xs(ms(nelb(kb),2))
            x3=xs(ms(nelb(kb),3))
            y1=ys(ms(nelb(kb),1))
            y2=ys(ms(nelb(kb),2))
            y3=ys(ms(nelb(kb),3))
            z1=zs(ms(nelb(kb),1))
            z2=zs(ms(nelb(kb),2))
            z3=zs(ms(nelb(kb),3))
            dr=difdr(x1,y1,z1,x2,y2,z2,x3,y3,z3)    !面の形状係数
c            dr=difdr3(x1,y1,z1,x2,y2,z2,x3,y3,z3)    !面の形状係数
            drav=drav+dr/dble(nb)
            if(dr.gt.drmax)drmax=dr
          enddo
          if(drmax.lt.drstd)goto 100
c
c         ksの位置を最適化する

          nel0=nb
          do 200 kb=1,nb
          do 200 j=1,3
200       ms0(kb,j)=ms(nelb(kb),j)
          write(*,*)"ks:",ks,"/",ns
          write(*,*)"befoer optimize"
          write(*,*)"drmax:",sngl(drmax),",drav:",sngl(drav)
          call pltsetnp(ns,xs,ys,zs,ind)
          call pltsetme(nel0,ms0,3d0,2,ind,ine)
          pause

c          call decidens(ms,nb,nelb,neb,xs,ys,zs,x0,y0,z0)
          call decidens2(ms,nb,nelb,neb,xs,ys,zs,x0,y0,z0,ifg,drmax)
          if(ifg.eq.0)write(*,*)"no optimized"
          write(*,*)"x:",sngl(xs(ks)),sngl(ys(ks)),sngl(zs(ks))
          if(ifg.eq.1)then
            xs(ks)=xs(ks)*(1d0-ceft)+x0*ceft
            ys(ks)=ys(ks)*(1d0-ceft)+y0*ceft
            zs(ks)=zs(ks)*(1d0-ceft)+z0*ceft
          end if
          write(*,*)"x:",sngl(xs(ks)),sngl(ys(ks)),sngl(zs(ks))

          write(*,*)"after optimize"
          do kb=1,nb
            if(ks.ne.ms(nelb(kb),neb(kb)))stop
          enddo
          drmax=0d0
          drav=0d0
          do kb=1,nb
            x1=xs(ms(nelb(kb),1))
            x2=xs(ms(nelb(kb),2))
            x3=xs(ms(nelb(kb),3))
            y1=ys(ms(nelb(kb),1))
            y2=ys(ms(nelb(kb),2))
            y3=ys(ms(nelb(kb),3))
            z1=zs(ms(nelb(kb),1))
            z2=zs(ms(nelb(kb),2))
            z3=zs(ms(nelb(kb),3))
            dr=difdr(x1,y1,z1,x2,y2,z2,x3,y3,z3)    !面の形状係数
c            dr=difdr3(x1,y1,z1,x2,y2,z2,x3,y3,z3)    !面の形状係数
            drav=drav+dr/dble(nb)
            if(dr.gt.drmax)drmax=dr
          enddo
          write(*,*)"drmax:",sngl(drmax),",drav:",sngl(drav)
          call pltsetnp(ns,xs,ys,zs,ind)
          call pltsetme(nel0,ms0,3d0,2,ind,ine)
          pause
100     continue
c
        return
        end
c
c       節点ksを含む面
c
        subroutine nssur(ks,nels,ms,nb,nelb,neb)
        include "header.h"
        dimension ms(ine,3)
        dimension nelb(ine),neb(ine)
c
        nb=0
        do 100 i=1,nels
        do 100 j=1,3
          if(ks.ne.ms(i,j))goto 100
          nb=nb+1
          nelb(nb)=i    !面番号
          neb(nb)=j     !面内番号
100     continue
        return
        end
c
c       面の形状係数を計算する
c
        function difdr(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        xg=(x1+x2+x3)/3d0
        yg=(y1+y2+y3)/3d0
        zg=(z1+z2+z3)/3d0
        r1=dsqrt((x1-xg)**2+(y1-yg)**2+(z1-zg)**2)
        r2=dsqrt((x2-xg)**2+(y2-yg)**2+(z2-zg)**2)
        r3=dsqrt((x3-xg)**2+(y3-yg)**2+(z3-zg)**2)
c
        a=r1+r2+r3
        b=r1**2+r2**2+r3**2
        difdr=dsqrt(9d0*b-3d0*a**2)/a
        return
        end
c
c       面の形状係数を計算する(無次元化しない)
c
        function difdr2(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        xg=(x1+x2+x3)/3d0
        yg=(y1+y2+y3)/3d0
        zg=(z1+z2+z3)/3d0
        r1=dsqrt((x1-xg)**2+(y1-yg)**2+(z1-zg)**2)
        r2=dsqrt((x2-xg)**2+(y2-yg)**2+(z2-zg)**2)
        r3=dsqrt((x3-xg)**2+(y3-yg)**2+(z3-zg)**2)
c
        a=r1+r2+r3
        b=r1**2+r2**2+r3**2
c        difdr=dsqrt(9d0*b-3d0*a**2)/a
        difdr2=dsqrt(9d0*b-3d0*a**2)
        return
        end
c
c       表面の節点ksの最適化した位置を探す
c
        subroutine decidens(ms,nb,nelb,neb,xs,ys,zs,x0,y0,z0)
        include "header.h"
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        dimension nelb(ine),neb(ine)
c
        xmax=-1d0
        ymax=-1d0
        zmax=-1d0
        xmin= 1d0
        ymin= 1d0
        zmin= 1d0
        do 100 kb=1,nb
        do 100 j=1,3
          if(xmax.lt.xs(ms(nelb(kb),j)))xmax=xs(ms(nelb(kb),j))
          if(ymax.lt.ys(ms(nelb(kb),j)))ymax=ys(ms(nelb(kb),j))
          if(zmax.lt.zs(ms(nelb(kb),j)))zmax=zs(ms(nelb(kb),j))
          if(xmin.gt.xs(ms(nelb(kb),j)))xmin=xs(ms(nelb(kb),j))
          if(ymin.gt.ys(ms(nelb(kb),j)))ymin=ys(ms(nelb(kb),j))
          if(zmin.gt.zs(ms(nelb(kb),j)))zmin=zs(ms(nelb(kb),j))
100     continue
        xmax=xmax+1d0
        ymax=ymax+1d0
        zmax=zmax+1d0
        xmin=xmin-1d0
        ymin=ymin-1d0
        zmin=zmin-1d0
c
c       10分割を5回繰り返す
        nx=10
        ny=10
        nz=10
        do nn=1,10
          drmin=1d3
          do 200 kx=0,nx-1
          do 200 ky=0,ny-1
          do 200 kz=0,nz-1
            x=xmin+(xmax-xmin)*(dble(kx)+5d-1)/dble(nx)
            y=ymin+(ymax-ymin)*(dble(ky)+5d-1)/dble(ny)
            z=zmin+(zmax-zmin)*(dble(kz)+5d-1)/dble(nz)
            drmax=-1d3
            do kb=1,nb
              kelb=nelb(kb)
              keb=neb(kb)
              x1=xs(ms(kelb,1))
              x2=xs(ms(kelb,2))
              x3=xs(ms(kelb,3))
              y1=ys(ms(kelb,1))
              y2=ys(ms(kelb,2))
              y3=ys(ms(kelb,3))
              z1=zs(ms(kelb,1))
              z2=zs(ms(kelb,2))
              z3=zs(ms(kelb,3))
              if(keb.eq.1)then
                x1=x
                y1=y
                z1=z
              elseif(keb.eq.2)then
                x2=x
                y2=y
                z2=z
              elseif(keb.eq.3)then
                x3=x
                y3=y
                z3=z
              endif
c              dr=difdr2(x1,y1,z1,x2,y2,z2,x3,y3,z3)        !面の形状係数
              dr=difdr3(x1,y1,z1,x2,y2,z2,x3,y3,z3)        !面の形状係数
              if(dr.gt.drmax)drmax=dr
            enddo
            if(drmax.lt.drmin)then
              drmin=drmax
              nxmin=kx
              nymin=ky
              nzmin=kz
c              write(*,*)"nn,kx:",nn,kx,ky,kz,sngl(drmax)
c              pause
            endif
200       continue
          xmin=xmin+(xmax-xmin)*dble(nxmin)/dble(nx)
          ymin=ymin+(ymax-ymin)*dble(nymin)/dble(ny)
          zmin=zmin+(zmax-zmin)*dble(nzmin)/dble(nz)
          xmax=xmin+(xmax-xmin)*dble(nxmin+1)/dble(nx)
          ymax=ymin+(ymax-ymin)*dble(nymin+1)/dble(ny)
          zmax=zmin+(zmax-zmin)*dble(nzmin+1)/dble(nz)
        enddo
        x0=xmin+(xmax-xmin)*dble(nxmin)/dble(nx)
        y0=ymin+(ymax-ymin)*dble(nymin)/dble(ny)
        z0=zmin+(zmax-zmin)*dble(nzmin)/dble(nz)
c
        return
        end
c
c       面の形状係数を計算する(外心から計算)
c
        function difdr3(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        xg=(x1+x2+x3)/3d0
        yg=(y1+y2+y3)/3d0
        zg=(z1+z2+z3)/3d0
        call circm(x1,y1,z1,x2,y2,z2,x3,y3,z3,xc,yc,zc)
c        difdr3=dsqrt((xc-xg)**2d0+(yc-yg)**2d0+(zc-zg)**2d0)
        difdr3=dsqrt((x1-xg)**2d0+(y1-yg)**2d0+(z1-zg)**2d0)
c
        return
        end
c
c       直交座標上の3角形の外心を探す
c
        subroutine circm(x1,y1,z1,x2,y2,z2,x3,y3,z3,xc,yc,zc)
        implicit double precision(a-h,o-z)
c
        ck1=(x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1)
        ck2=(x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)
        ck=ck1/ck2
        x4=x1+ck*(x2-x1)
        y4=y1+ck*(y2-y1)
        z4=z1+ck*(z2-z1)
c
        vx=dsqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        vy=dsqrt((x3-x4)**2+(y3-y4)**2+(z3-z4)**2)
c
        x10=0d0
        y10=0d0
        x20=1d0
        y20=0d0
        x30=ck
        y30=1d0*vy/vx
c
        ckx1=y30*x20**2+y30*y20**2
     &      +y20*x10**2+y20*y10**2
     &      +y10*x30**2+y10*y30**2
     &      -y10*x20**2-y10*y20**2
     &      -y20*x30**2-y20*y30**2
     &      -y30*x10**2-y30*y10**2
        ckx2=2d0*(x10*y20+x20*y30+x30*y10-x30*y20-x20*y10-x10*y30)
        ckx=ckx1/ckx2
        cky1=x30*x20**2+x30*y20**2
     &      +x20*x10**2+x20*y10**2
     &      +x10*x30**2+x10*y30**2
     &      -x10*x20**2-x10*y20**2
     &      -x20*x30**2-x20*y30**2
     &      -x30*x10**2-x30*y10**2
        cky2=2d0*(x10*y20+x20*y30+x30*y10-x30*y20-x20*y10-x10*y30)
        cky=-cky1/cky2
c
c        dr1=dsqrt((x10-ckx)**2d0+(y10-cky)**2d0)
c        dr2=dsqrt((x20-ckx)**2d0+(y20-cky)**2d0)
c        dr3=dsqrt((x30-ckx)**2d0+(y30-cky)**2d0)
c        write(*,*)"ck:",sngl(ck),sngl(ckx),sngl(cky)
c        write(*,*)"dr0:",sngl(dr1),sngl(dr2),sngl(dr3)
c
        xc=x1+ckx*(x2-x1)+cky*(x3-x4)*vx/vy
        yc=y1+ckx*(y2-y1)+cky*(y3-y4)*vx/vy
        zc=z1+ckx*(z2-z1)+cky*(z3-z4)*vx/vy
c
c        write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
c        write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
c        write(*,*)"x3:",sngl(x3),sngl(y3),sngl(z3)
c        write(*,*)"xc:",sngl(xc),sngl(yc),sngl(zc)
        dr1=dsqrt((x1-xc)**2d0+(y1-yc)**2d0+(z1-zc)**2d0)
        dr2=dsqrt((x2-xc)**2d0+(y2-yc)**2d0+(z2-zc)**2d0)
        dr3=dsqrt((x3-xc)**2d0+(y3-yc)**2d0+(z3-zc)**2d0)
        if(dabs(dr1-dr2).gt.1d-4 .or. dabs(dr2-dr3).gt.1d-4)then
          write(*,*)"in circm"
          write(*,*)"dr:",sngl(dr1),sngl(dr2),sngl(dr3)
          stop
        endif
c
        return
        end
c
c       表面の節点ksの最適化した位置を探す
c
        subroutine decidens2(ms,nb,nelb,neb,xs,ys,zs,x0,y0,z0,ifg,drmax)
        include "header.h"
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        dimension nelb(ine),neb(ine)
c
        ifg=0   !0:節点が見つからなかった,1:節点が見つかった
c
        xmax=-1d0
        ymax=-1d0
        zmax=-1d0
        xmin= 1d0
        ymin= 1d0
        zmin= 1d0
        do 100 kb=1,nb
        do 100 j=1,3
          if(xmax.lt.xs(ms(nelb(kb),j)))xmax=xs(ms(nelb(kb),j))
          if(ymax.lt.ys(ms(nelb(kb),j)))ymax=ys(ms(nelb(kb),j))
          if(zmax.lt.zs(ms(nelb(kb),j)))zmax=zs(ms(nelb(kb),j))
          if(xmin.gt.xs(ms(nelb(kb),j)))xmin=xs(ms(nelb(kb),j))
          if(ymin.gt.ys(ms(nelb(kb),j)))ymin=ys(ms(nelb(kb),j))
          if(zmin.gt.zs(ms(nelb(kb),j)))zmin=zs(ms(nelb(kb),j))
100     continue
        xmax0=xmax+1d0
        ymax0=ymax+1d0
        zmax0=zmax+1d0
        xmin0=xmin-1d0
        ymin0=ymin-1d0
        zmin0=zmin-1d0
c
c       面積
        ar0=0d0
        do kb=1,nb
          kelb=nelb(kb)
          x1=xs(ms(kelb,1))
          x2=xs(ms(kelb,2))
          x3=xs(ms(kelb,3))
          y1=ys(ms(kelb,1))
          y2=ys(ms(kelb,2))
          y3=ys(ms(kelb,3))
          z1=zs(ms(kelb,1))
          z2=zs(ms(kelb,2))
          z3=zs(ms(kelb,3))
          ar=artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
          if(ar.lt.0d0)write(*,*)"in decidens2 ar:",sngl(ar)
          if(ar.lt.0d0)stop
          ar0=ar0+ar
        enddo
c
        drmin=1d3
        do while(1.eq.1)
          xmax=xmax0
          ymax=ymax0
          zmax=zmax0
          xmin=xmin0
          ymin=ymin0
          zmin=zmin0
c         -------------------------------------------------------------------------------
c         10分割を5回繰り返す
          ifg2=0   !0:節点が見つからなかった,1:節点が見つかった
          nx=10
          ny=10
          nz=10
          do nn=1,10
            do 200 kx=0,nx-1
            do 200 ky=0,ny-1
            do 200 kz=0,nz-1
              x=xmin+(xmax-xmin)*(dble(kx)+5d-1)/dble(nx)
              y=ymin+(ymax-ymin)*(dble(ky)+5d-1)/dble(ny)
              z=zmin+(zmax-zmin)*(dble(kz)+5d-1)/dble(nz)
              dr=-1d3
              ar=0d0
              do kb=1,nb
                kelb=nelb(kb)
                keb=neb(kb)
                x1=xs(ms(kelb,1))
                x2=xs(ms(kelb,2))
                x3=xs(ms(kelb,3))
                y1=ys(ms(kelb,1))
                y2=ys(ms(kelb,2))
                y3=ys(ms(kelb,3))
                z1=zs(ms(kelb,1))
                z2=zs(ms(kelb,2))
                z3=zs(ms(kelb,3))
                if(keb.eq.1)then
                  x1=x
                  y1=y
                  z1=z
                elseif(keb.eq.2)then
                  x2=x
                  y2=y
                  z2=z
                elseif(keb.eq.3)then
                  x3=x
                  y3=y
                  z3=z
                endif
                ar=ar+artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
                dr2=difdr(x1,y1,z1,x2,y2,z2,x3,y3,z3)        !面の形状係数
c                dr2=difdr2(x1,y1,z1,x2,y2,z2,x3,y3,z3)        !面の形状係数
c                dr2=difdr3(x1,y1,z1,x2,y2,z2,x3,y3,z3)        !面の形状係数
                if(dr2.gt.dr)dr=dr2
              enddo
              if(ar.gt.ar0)goto 200
              if(dr.ge.drmin)goto 200
              ifg2=1
              armin=ar
              drmin=dr
              nxmin=kx
              nymin=ky
              nzmin=kz
c              write(*,*)"nn,kx:",nn,kx,ky,kz,sngl(drmax)
c              write(*,*)"ar:",sngl(ar),"/",sngl(ar0)
c              pause
200         continue
            xmin=xmin+(xmax-xmin)*dble(nxmin)/dble(nx)
            ymin=ymin+(ymax-ymin)*dble(nymin)/dble(ny)
            zmin=zmin+(zmax-zmin)*dble(nzmin)/dble(nz)
            xmax=xmin+(xmax-xmin)*dble(nxmin+1)/dble(nx)
            ymax=ymin+(ymax-ymin)*dble(nymin+1)/dble(ny)
            zmax=zmin+(zmax-zmin)*dble(nzmin+1)/dble(nz)
          enddo
          if(ifg2.eq.0)exit
c         -------------------------------------------------------------------------------
c
          if(drmin.ge.drmax)exit
          ifg=1
          ar0=armin
          x0=xmin+(xmax-xmin)*dble(nxmin)/dble(nx)
          y0=ymin+(ymax-ymin)*dble(nymin)/dble(ny)
          z0=zmin+(zmax-zmin)*dble(nzmin)/dble(nz)
        enddo
c

        ar=0d0
        do kb=1,nb
          kelb=nelb(kb)
          keb=neb(kb)
          x1=xs(ms(kelb,1))
          x2=xs(ms(kelb,2))
          x3=xs(ms(kelb,3))
          y1=ys(ms(kelb,1))
          y2=ys(ms(kelb,2))
          y3=ys(ms(kelb,3))
          z1=zs(ms(kelb,1))
          z2=zs(ms(kelb,2))
          z3=zs(ms(kelb,3))
          if(keb.eq.1)then
            x1=x0
            y1=y0
            z1=z0
          elseif(keb.eq.2)then
            x2=x0
            y2=y0
            z2=z0
          elseif(keb.eq.3)then
            x3=x0
            y3=y0
            z3=z0
          endif
          ar=ar+artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        enddo
        write(*,*)"ar:",sngl(ar),"/",sngl(ar0)
c
        return
        end



c---------------------------------------------------------------------------
c
c       面付近に節点を発生させる
c
c---------------------------------------------------------------------------
        subroutine npsuf(ns,nels,xs,ys,zs,ms,nt,xt,yt,zt)
        include "header.h"
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xt(ind),yt(ind),zt(ind)
        dimension xk(ind),yk(ind),zk(ind),dl(ind)
        dimension xg(ind),yg(ind),zg(ind)
        dimension mk(ind),mt(ind)
c
c       面付近に格子を発生させる
        nk=0
        do kels=1,nels
          x1=xs(ms(kels,1))
          x2=xs(ms(kels,2))
          x3=xs(ms(kels,3))
          y1=ys(ms(kels,1))
          y2=ys(ms(kels,2))
          y3=ys(ms(kels,3))
          z1=zs(ms(kels,1))
          z2=zs(ms(kels,2))
          z3=zs(ms(kels,3))
c
c         面からの距離
          dl12=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
          dl23=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
          dl31=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
          dlav=(dl12+dl23+dl31)/3d0
c
c         外積
          v1x=x2-x1
          v1y=y2-y1
          v1z=z2-z1
          v2x=x3-x1
          v2y=y3-y1
          v2z=z3-z1
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz) !外積:|dn|=1
c
c         内心
          call innerpoint(x1,y1,z1,x2,y2,z2,x3,y3,z3,xin,yin,zin)
c
c         辺を共有する面との角度を比較する

c
c         格子を発生させる
          nk=nk+1
          xk(nk)=xin-dnx*dlav
          yk(nk)=yin-dny*dlav
          zk(nk)=zin-dnz*dlav
          dl(nk)=dlav
        enddo
c------------------------------------------------------------------------------
c
c       nt優先
c
c------------------------------------------------------------------------------
c
c       格子nkを削除する
c
        do 100 i=1,nk
100     mk(i)=0     !0:削除しない,1:削除する
c
c       nkと比較
        do 200 i=1,nk-1
          if(mk(i).eq.1)goto 200
          do 300 j=i+1,nk
            dr=dsqrt((xk(i)-xk(j))**2+(yk(i)-yk(j))**2+(zk(i)-zk(j))**2)
            dlav=(dl(i)+dl(j))/2d0
            if(dr.lt.dlav)mk(j)=1
300       continue
200     continue
c
c       nsと比較
        do 400 i=1,nk
          if(mk(i).eq.1)goto 400
          do 500 j=1,ns
            dr=dsqrt((xk(i)-xs(j))**2+(yk(i)-ys(j))**2+(zk(i)-zs(j))**2)
            if(dr.lt.0.9d0*dl(i))mk(i)=1
500       continue
400     continue
c
c       格子ntを削除する
c
        do 800 i=1,nt
800     mt(i)=0     !0:削除しない,1:削除する
c
c       ntと比較
        do 600 i=1,nt
          do 700 j=1,nk
            if(mk(j).eq.1)goto 700
            dr=dsqrt((xt(i)-xk(j))**2+(yt(i)-yk(j))**2+(zt(i)-zk(j))**2)
            if(dr.lt.dl(j))mt(i)=1
700       continue
600     continue
c
c       並べ替え
c
        ng=0
        k=0
        do 900 i=1,nk
          if(mk(i).eq.1)goto 900
          ng=ng+1

          k=k+1
          write(*,*)"pass1 i,nk:",k,nk

          xg(ng)=xk(i)
          yg(ng)=yk(i)
          zg(ng)=zk(i)
900     continue
        k=0
        do 1000 i=1,nt
          if(mt(i).eq.1)goto 1000
          ng=ng+1

          k=k+1
          write(*,*)"pass2 i,nt:",k,nt

          xg(ng)=xt(i)
          yg(ng)=yt(i)
          zg(ng)=zt(i)
1000    continue
c
        nt=ng
        do i=1,nt
          xt(i)=xg(i)
          yt(i)=yg(i)
          zt(i)=zg(i)
        enddo
        write(*,*)"npsuf nt:",nt
c
        return
        end
c----------------------------------------------------------------------
c
c       3角形の内心を算出する
c
c----------------------------------------------------------------------
        subroutine innerpoint(x1,y1,z1,x2,y2,z2,x3,y3,z3,xin,yin,zin)
        implicit double precision(a-h,o-z)
c
        p12=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)  !|p12|=|p21|
        p21=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        p23=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)  !|p23|=|p32|
        p32=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        p31=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)  !|p31|=|p13|
        p13=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
c
        a1=(x1-x2)/p12+(x3-x2)/p32
        a2=(y1-y2)/p12+(y3-y2)/p32
        a3=(z1-z2)/p12+(z3-z2)/p32
        b1=-(x2-x3)/p23-(x1-x3)/p13
        b2=-(y2-y3)/p23-(y1-y3)/p13
        b3=-(z2-z3)/p23-(z1-z3)/p13
c
        if(dabs(x1).lt.1d-8 .and. dabs(x2).lt.1d-8
     &                      .and. dabs(x3).lt.1d-8)then    !x平面上
          xi2=y2
          xi3=y3
          ai=a2
          bi=b2
        else
          xi2=x2
          xi3=x3
          ai=a1
          bi=b1
        endif
c
        d2=-bi*(x3+y3+z3-x2-y2-z2)+(b1+b2+b3)*(-xi2+xi3)
        d3= ai*(x3+y3+z3-x2-y2-z2)-(a1+a2+a3)*(-xi2+xi3)
        d=ai*(b1+b2+b3)-(a1+a2+a3)*bi
        dl2=d2/d
        dl3=d3/d
c
        xin=x2+dl2*((x1-x2)/p12+(x3-x2)/p32)
        yin=y2+dl2*((y1-y2)/p12+(y3-y2)/p32)
        zin=z2+dl2*((z1-z2)/p12+(z3-z2)/p32)
c
        xin0=x3+dl3*((x2-x3)/p23+(x1-x3)/p13)
        yin0=y3+dl3*((y2-y3)/p23+(y1-y3)/p13)
        zin0=z3+dl3*((z2-z3)/p23+(z1-z3)/p13)
c
        if(dabs(xin-xin0).gt.1d-5 .or. dabs(yin-yin0).gt.1d-5
     &                            .or. dabs(zin-zin0).gt.1d-5)then
          write(*,*)"innerpoint"
          write(*,*)"xin",sngl(xin),sngl(yin),sngl(zin)
          write(*,*)"xin",sngl(xin0),sngl(yin0),sngl(zin0)
          stop
        endif
c
c        write(*,*)"xin",sngl(xin),sngl(yin),sngl(zin)
c        write(*,*)"xin",sngl(xin0),sngl(yin0),sngl(zin0)
c       検算
c        ar=artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
c        dr=2d0*ar/(p12+p23+p31)
c        write(*,*)"dr",sngl(dr)
c
        return
        end
