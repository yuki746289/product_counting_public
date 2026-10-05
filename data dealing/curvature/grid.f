c---------------------------------------------------------------------------
c
c       四角:内部の節点を含む
c
c---------------------------------------------------------------------------
        subroutine grid(mp,mmp,np0,x0,y0,z0,xg,yg,zg,ind)
        implicit double precision(a-h,o-z)
        dimension mmp(ind),lp(ind)
        dimension x0(ind),y0(ind),z0(ind)
c
        xmin=-1d0
        xmax= 1d0
        ymin=-1d0
        ymax= 1d0
        zmin=-1d0
        zmax= 1d0
c
        dd=0.5d0
        il=int((xmax-xmin)/dd)
        jl=int((ymax-ymin)/dd)
        kl=int((zmax-zmin)/dd)
        np0=0
        do 100 i=1,il+1
          do 100 j=1,jl+1
            do 100 k=1,kl+1
              np0=np0+1
              x0(np0)=xmin+dd*int(i-1)
              y0(np0)=ymin+dd*int(j-1)
              z0(np0)=zmin+dd*int(k-1)
100     continue
c
c       mpの順番:mmp(np0-1)
c
        do 200 kp0=1,np0
200     lp(kp0)=0   !0:数えてない,1:数えた
c
        do 300 i=1,np0  !mmp(np0)
          dmin=999d0
          do 400 kp0=1,np0
            if(lp(kp0).eq.1)goto 400
            xp=x0(kp0)
            yp=y0(kp0)
            zp=z0(kp0)
            dd=dsqrt((xp-xg)**2+(yp-yg)**2+(zp-zg)**2)
            if(dd.lt.dmin)mp=kp0
            if(dd.lt.dmin)dmin=dd
400       continue
          mmp(i)=mp
          lp(mp)=1
300     continue
        mp=mmp(1)
c
        return
        end
c---------------------------------------------------------------------------
c
c       球
c
c---------------------------------------------------------------------------
        subroutine grid1(ns,xs,ys,zs,dl0,ind)
        implicit double precision(a-h,o-z)
        dimension xs(ind),ys(ind),zs(ind)
        data PI/3.14159265358979323/ !84626433832795
c
        dl0=2.2d0   !領域の幅
        dr=1d0      !球の直径
        dphi=10d0/180d0*PI  !縦回転
        dpsi=10d0/180d0*PI  !横回転
        ns=0
c
        r=dble(i)*dr
c       軸の下
        ns=ns+1
        xs(ns)=0d0
        ys(ns)=-r
        zs(ns)=0d0
c       軸の上
        ns=ns+1
        xs(ns)=0d0
        ys(ns)=r
        zs(ns)=0d0
c       側面
        do 10 kphi=0,35    !縦回転
        do 10 kpsi=0,17    !横回転
          phi=dble(j)*dphi    !縦回転
          psi=dble(k)*dpsi    !横回転
          if(kphi.eq.9 .or. kphi.eq.27)goto 10 !軸の上と下は重なるので飛ばす
          ns=ns+1
          xs(ns)=dr*dcos(phi)*dcos(psi)
          ys(ns)=dr*dsin(phi)
          zs(ns)=dr*dcos(phi)*dsin(psi)
10     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       トーラス
c
c---------------------------------------------------------------------------
        subroutine grid2(mp,mmp,np0,x0,y0,z0,ind)
        implicit double precision(a-h,o-z)
        dimension x0(ind),y0(ind),z0(ind)
        dimension mmp(ind)
        data PI/3.1415926535897932384626433832795/
c
        rc=5d0  !ドーナツの中心の半径
        rd=2d0  !ドーナツの半径
        shita=20d0/180d0*PI
c
c       中心の節点
        np0=0
        do 100 i=0,17
          psi=shita*dble(i)
          np0=np0+1
          x0(np0)=rc*dcos(psi)
          y0(np0)=0d0
          z0(np0)=rc*dsin(psi)
100     continue
c
c       外側
        do 200 i=0,17
          phi=shita*dble(i)
          do 200 j=0,17
            psi=shita*dble(j)
            np0=np0+1
            x0(np0)=rd*dcos(phi)*dcos(psi)+rc*dcos(psi)
            y0(np0)=rd*dsin(phi)
            z0(np0)=rd*dcos(phi)*dsin(psi)+rc*dsin(psi)
200     continue
c
c       節点の順番
        do 300 kp0=1,np0
300     mmp(kp0)=kp0
        mp=mmp(1)
c
        return
        end
c---------------------------------------------------------------------------
c
c       四角:表面の節点のみ
c
c---------------------------------------------------------------------------
        subroutine grid3(ns,xs,ys,zs,dl0,mgr,ind)
        implicit double precision(a-h,o-z)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xs(ind),ys(ind),zs(ind)
        dimension mgr(ind)
c
        dl0=2.2d0   !領域の幅
c
        xmin=-1d0
        xmax= 1d0
        ymin=-1d0
        ymax= 1d0
        zmin=-1d0
        zmax= 1d0
c
        dd=0.5d0
        il=int((xmax-xmin)/dd)
        jl=int((ymax-ymin)/dd)
        kl=int((zmax-zmin)/dd)
        np=0
        do 100 i=1,il+1
        do 100 j=1,jl+1
        do 100 k=1,kl+1
          np=np+1
          xx(np)=xmin+dd*int(i-1)
          yy(np)=ymin+dd*int(j-1)
          zz(np)=zmin+dd*int(k-1)
100     continue
c
        ns=0
        do 200 kp=1,np
          ifg=0
          if(xx(kp).eq.-1d0 .or. xx(kp).eq.1d0)ifg=1
          if(yy(kp).eq.-1d0 .or. yy(kp).eq.1d0)ifg=1
          if(zz(kp).eq.-1d0 .or. zz(kp).eq.1d0)ifg=1
          if(ifg.eq.0)goto 200
          ns=ns+1
          xs(ns)=xx(kp)
          ys(ns)=yy(kp)
          zs(ns)=zz(kp)
200     continue
c
c       節点グループ
c         面要素を発生させるする時に同じ節点グループから節点候補を選ぶ
        do 300 ks=1,ns
300     mgr(ks)=-1
        do 400 ks=1,ns
          ifg=0
          if(xs(ks).eq.-1d0)ifg=ifg+1
          if(xs(ks).eq. 1d0)ifg=ifg+1
          if(ys(ks).eq.-1d0)ifg=ifg+1
          if(ys(ks).eq. 1d0)ifg=ifg+1
          if(zs(ks).eq.-1d0)ifg=ifg+1
          if(zs(ks).eq. 1d0)ifg=ifg+1
c
          if(ifg.gt.1)mgr(ks)=0    !境界
          if(ifg.gt.1)goto 400
          if(xs(ks).eq.-1d0)mgr(ks)=1
          if(xs(ks).eq. 1d0)mgr(ks)=2
          if(ys(ks).eq.-1d0)mgr(ks)=3
          if(ys(ks).eq. 1d0)mgr(ks)=4
          if(zs(ks).eq.-1d0)mgr(ks)=5
          if(zs(ks).eq. 1d0)mgr(ks)=6
400     continue
        do ks=1,ns
        write(*,*)"mgr:",ks,mgr(ks)
        enddo
c
        return
        end
c---------------------------------------------------------------------------
c
c       四角:表面の節点のみ
c
c---------------------------------------------------------------------------
        subroutine grid4(ns,xs,ys,zs,dl0,mgr,ind)
        implicit double precision(a-h,o-z)
        dimension xt(ind),yt(ind),zt(ind)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xs0(ind),ys0(ind),zs0(ind)
        dimension mgr(ind)
        dimension mm(ind)
c
        dl0=2d0     !領域の幅
        n=5         !分割数
c
c       x=[0,1]
c       y=[0,1]
c       z=0
        ns=0
        do i=0,n*2
          x=dble(i)*dl0/dble(n*2)
          if(mod(i,2) .eq. 0)then
            do j=0,n
              y=dble(j)*dl0/dble(n)
              ns=ns+1
              xs(ns)=x
              ys(ns)=y
              zs(ns)=0d0
            enddo
          else
            do j=0,n-1
              y=(dble(j)+5d-1)*dl0/dble(n)
              ns=ns+1
              xs(ns)=x
              ys(ns)=y
              zs(ns)=0d0
            enddo
          endif
        enddo
c
c       コピー
        ns0=ns
        do ks=1,ns
          xs0(ks)=xs(ks)
          ys0(ks)=ys(ks)
          zs0(ks)=zs(ks)
        enddo
c
c       x=[0,1]
c       y=[0,1]
c       z=1d0
        do ks=1,ns0
          ns=ns+1
          xs(ns)=xs0(ks)
          ys(ns)=ys0(ks)
          zs(ns)=dl0
        enddo
c
c       x=0d0,dl
c       y=[0,1]
c       z=[0,1]
        do ks=1,ns0
          ns=ns+1
          xs(ns)=0d0
          ys(ns)=ys0(ks)
          zs(ns)=xs0(ks)
        enddo
c
        do ks=1,ns0
          ns=ns+1
          xs(ns)=dl0
          ys(ns)=ys0(ks)
          zs(ns)=xs0(ks)
        enddo
c
c       x=[0,1]
c       y=0d0,dl
c       z=[0,1]
        do ks=1,ns0
          ns=ns+1
          xs(ns)=xs0(ks)
          ys(ns)=0d0
          zs(ns)=ys0(ks)
        enddo
c
        do ks=1,ns0
          ns=ns+1
          xs(ns)=xs0(ks)
          ys(ns)=dl0
          zs(ns)=ys0(ks)
        enddo
c
c       重なった節点を除く
        nt=0
        do 10 ks=1,ns
          do kt=1,nt
            if(xt(kt).eq.xs(ks) .and. yt(kt).eq.ys(ks) .and.
     &         zt(kt).eq.zs(ks))goto 10
          enddo
          nt=nt+1
          xt(nt)=xs(ks)
          yt(nt)=ys(ks)
          zt(nt)=zs(ks)
          mm(nt)=0
10      continue
c
c       節点番号を割り振り直す(重心から遠い順)
        xc=1d0
        yc=1d0
        zc=1d0
        ns=0
        do while(1.eq.1)
          dmax=0d0
          do 100 kt=1,nt
            if(mm(kt).eq.1)goto 100
            dd=dsqrt((xt(kt)-xc)**2+(yt(kt)-yc)**2+(zt(kt)-zc)**2)
            if(dd.gt.dmax)n0=kt
            if(dd.gt.dmax)dmax=dd
100       continue
          ns=ns+1
          xs(ns)=xt(n0)
          ys(ns)=yt(n0)
          zs(ns)=zt(n0)
          mm(n0)=1
          if(ns.eq.nt)exit
        end do
c
c       節点グループ
c         面要素を発生させるする時に同じ節点グループから節点候補を選ぶ
        do 300 ks=1,ns
300     mgr(ks)=-1
        do 400 ks=1,ns
          ifg=0
          if(xs(ks).eq.0d0)ifg=ifg+1
          if(xs(ks).eq.2d0)ifg=ifg+1
          if(ys(ks).eq.0d0)ifg=ifg+1
          if(ys(ks).eq.2d0)ifg=ifg+1
          if(zs(ks).eq.0d0)ifg=ifg+1
          if(zs(ks).eq.2d0)ifg=ifg+1
c
          if(ifg.gt.1)mgr(ks)=0    !境界
          if(ifg.gt.1)goto 400
          if(xs(ks).eq.0d0)mgr(ks)=1
          if(xs(ks).eq.2d0)mgr(ks)=2
          if(ys(ks).eq.0d0)mgr(ks)=3
          if(ys(ks).eq.2d0)mgr(ks)=4
          if(zs(ks).eq.0d0)mgr(ks)=5
          if(zs(ks).eq.2d0)mgr(ks)=6
400     continue
c        do ks=1,ns
c      write(*,*)"mgr:",ks,mgr(ks),sngl(xs(ks)),sngl(ys(ks)),sngl(zs(ks))
c        enddo
c
        return
        end
