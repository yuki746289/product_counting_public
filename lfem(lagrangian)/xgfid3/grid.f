c---------------------------------------------------------------------------
c
c       四角
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
