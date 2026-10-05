        implicit double precision(a-h,o-z)
        parameter(ind=9999)
        dimension xx(ind),yy(ind),zz(ind)
        data PI/3.1415926535897932384626433832795/
c
        k1=1
        k2=2
        k3=3
        kk=4
c
c       p1(0,0,-1)
        xx(k1)=0d0
        yy(k1)=0d0
        zz(k1)=-1d0
c
c       p2(0,0,1)
        xx(k2)=0d0
        yy(k2)=0d0
        zz(k2)=1d0
c
c       p3(-1,0,0)
        xx(k3)=-1d0
        yy(k3)=0d0
        zz(k3)=0d0
c
c       p4(1,-1,0)
        xx(kk)=1d0
        yy(kk)=1d0
        zz(kk)=0d0
c
c       中点を計算する
        xg=(xx(k1)+xx(k2))/2d0
        yg=(yy(k1)+yy(k2))/2d0
        zg=(zz(k1)+zz(k2))/2d0
c
c       基準のベクトル
        vx=xx(kk)-xg
        vy=yy(kk)-yg
        vz=zz(kk)-zg
c
c       接線ベクトル
        vtx=xx(k3)-xg
        vty=yy(k3)-yg
        vtz=zz(k3)-zg
c
c       法線ベクトル
        v1x=xx(k2)-xx(k1)
        v1y=yy(k2)-yy(k1)
        v1z=zz(k2)-zz(k1)
        v2x=xx(k3)-xx(k1)
        v2y=yy(k3)-yy(k1)
        v2z=zz(k3)-zz(k1)
        call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
c
c       内積から角度を計算する
        phi=finpr(vx,vy,vz,dnx,dny,dnz) !法線との角度
        phi=PI-phi
        psi=finpr(vx,vy,vz,vtx,vty,vtz) !接線との角度
        if(phi.ge.0d0    .and. phi.le.PI/2d0)shita=psi
        if(phi.ge.PI/2d0 .and. phi.le.PI    )shita=2d0*PI-psi
c
        if(shita.lt.0d0 .or. shita.gt.2d0*PI)then
          write(*,*)"fshita=",shita/PI*180d0
          stop
        endif
        write(*,*)"shita=",shita/PI*180d0
c
        stop
        end
c
c------------------------------------------------------------
c
c       外積の計算
c       |n→|=1
c------------------------------------------------------------
        subroutine cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
        implicit double precision(a-h,o-z)
        dnx=v1y*v2z-v1z*v2y
        dny=v1z*v2x-v1x*v2z
        dnz=v1x*v2y-v1y*v2x
        dn=dsqrt(dnx**2+dny**2+dnz**2)
        if(dn.lt.1d-5)write(*,*)"外積:",dn
        if(dn.lt.1d-5)stop
        dnx=dnx/dn
        dny=dny/dn
        dnz=dnz/dn
c       検算
c        dd1=v1x*dnx+v1y*dny+v1z*dnz
c        dd2=v2x*dnx+v2y*dny+v2z*dnz
c        write(*,*)"dd1,dd2",dd1,dd2
c
        return
        end
c
c       内積から角度を計算する
        function finpr(px,py,pz,qx,qy,qz)
        implicit double precision(a-h,o-z)
        data PI/3.14159265358979323/ !84626433832795
c
        d1=dsqrt(px*px+py*py+pz*pz);
        d2=dsqrt(qx*qx+qy*qy+qz*qz);
        dd=px*qx+py*qy+pz*qz;
        d12=dd/d1/d2
c
        if(dabs(d12-(-1d0)).lt.1d-15)then       !180°:1d-16の誤差がでる
          shita=PI
        elseif(dabs(d12-1d0).lt.1d-15)then     !0°
          shita=0d0
        else
          shita=dacos(dd/d1/d2)                 !0°<shita<180°
        endif
c
        if(shita.lt.0d0 .or. shita.gt.PI)then
          write(*,*)'finpr:shita',shita
          stop
        endif
c
        finpr=shita
        return
        end
