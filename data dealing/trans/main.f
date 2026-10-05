        implicit double precision(a-h,o-z)
        dimension xx(4999),yy(4999),zz(4999)
        dimension xn(4999),yn(4999),zn(4999)
        dimension dd(4),dn(4)
        data PI/3.1415926535897932384626433832795/
c       -----------
c       x軸から反時計周りの角度φ
c       角度:dacos(x0/r0)
c       符号:dasin(z0/r0)/dabs(dasin(z0/r0))
c       -----------
c
c       φ=0°
        x0=1d0
        y0=1d0
        z0=0d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(0d0),sngl(dacos(x0/r0)/PI*180d0)
     &                             ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=45°
        x0=1d0
        y0=1d0
        z0=1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(45d0),sngl(dacos(x0/r0)/PI*180d0)
     &                             ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=90°
        x0=0d0
        y0=1d0
        z0=1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(90d0),sngl(dacos(x0/r0)/PI*180d0)
     &                             ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=135°
        x0=-1d0
        y0=1d0
        z0=1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(135d0),sngl(dacos(x0/r0)/PI*180d0)
     &                              ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=180°
        x0=-1d0
        y0=1d0
        z0=0d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(180d0),sngl(dacos(x0/r0)/PI*180d0)
     &                              ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=-135°
        x0=-1d0
        y0=1d0
        z0=-1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(-135d0),sngl(dacos(x0/r0)/PI*180d0)
     &                               ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=-90°
        x0=0d0
        y0=1d0
        z0=-1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(-90d0),sngl(dacos(x0/r0)/PI*180d0)
     &                               ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c       φ=-45°
        x0=1d0
        y0=1d0
        z0=-1d0
        r=dsqrt(x0**2+y0**2+z0**2)
        r0=dsqrt(x0**2+z0**2)
        phi=dasin(y0/r)
        write(*,*)"psi",sngl(-45d0),sngl(dacos(x0/r0)/PI*180d0)
     &                              ,sngl(dasin(z0/r0)/PI*180d0)
     &                             ,sngl(x0),sngl(z0)
c
c
c-------------------------------------------------------------------------------
c
c
        write(*,*)"----"
c       φ=90°
        y0=1d0
        r=1d0
        phi=dasin(y0/r)
        write(*,*)"phi",sngl(90d0),sngl(dasin(y0/r)/PI*180d0),sngl(y0)
c       φ=45°
        y0=1d0/dsqrt(2d0)
        r=1d0
        phi=dasin(y0/r)
        write(*,*)"phi",sngl(45d0),sngl(dasin(y0/r)/PI*180d0),sngl(y0)
c       φ=0°
        y0=0d0
        r=1d0
        phi=dasin(y0/r)
        write(*,*)"phi",sngl(0d0),sngl(dasin(y0/r)/PI*180d0),sngl(y0)
c       φ=-45°
        y0=-1d0/dsqrt(2d0)
        r=1d0
        phi=dasin(y0/r)
        write(*,*)"phi",sngl(-45d0),sngl(dasin(y0/r)/PI*180d0),sngl(y0)
c       φ=-90°
        y0=-1d0
        r=1d0
        phi=dasin(y0/r)
        write(*,*)"phi",sngl(-90d0),sngl(dasin(y0/r)/PI*180d0),sngl(y0)
c
c
c-------------------------------------------------------------------------------
c
c
        write(*,*)"-----"
        x0=1d0/dsqrt(2d0)
        y0=1d0
        z0=1d0/dsqrt(2d0)
        dphi=20d0/180d0*PI     !-90°<=φ<=90°
        dpsi=10d0/180d0*PI
c
        call rot12(x0,y0,z0,x1,y1,z1,dphi,dpsi)
        call rot21(x1,y1,z1,x2,y2,z2,dphi,dpsi)
c
        write(*,*)"x0",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"x1",sngl(x1),sngl(y1),sngl(z1)
        write(*,*)"x2",sngl(x2),sngl(y2),sngl(z2)
c
c
        x0= 1d-9
        y0= 1d0
        z0= 0d0
        dphi=-10d0/180d0*PI     !-90°<=φ<=90°
        dpsi=0d0/180d0*PI
c
        call rot12(x0,y0,z0,x1,y1,z1,dphi,dpsi)
        call rot21(x1,y1,z1,x2,y2,z2,dphi,dpsi)
c
        write(*,*)"x0",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"x1",sngl(x1),sngl(y1),sngl(z1)
        write(*,*)"x2",sngl(x2),sngl(y2),sngl(z2)
c
c       面の回転
        x0= 0d0
        y0= 0d0
        z0= 1d0
        phi=45d0/180d0*PI
        psi=45d0/180d0*PI
        eta=0d0/180d0*PI
c
        write(*,*)"x0",0,sngl(x0),sngl(y0),sngl(z0)
        do i=1,36
          call rot12r(x0,y0,z0,x1,y1,z1,phi,psi,eta)
c          call rot12(x0,y0,z0,x2,y2,z2,phi,psi)
          write(*,*)"x1",i,sngl(x1),sngl(y1),sngl(z1)
c          write(*,*)"x2",i,sngl(x2),sngl(y2),sngl(z2)
c          phi=phi+10d0/180d0*PI
c          psi=psi+10d0/180d0*PI
          eta=eta+10d0/180d0*PI
c          pause
        enddo
c
        x=0d0
        y=0d0
        z=1d0
        phi=dasin(1d0/dsqrt(3d0))
        psi=PI*3d0/4d0
        eta=0d0/180d0*PI
        write(*,*)"xyz:",sngl(x),sngl(y),sngl(z)
        write(*,*)"dgr:",sngl(phi/PI*180d0)
     &                  ,sngl(psi/PI*180d0),sngl(eta/PI*180d0)
        call rotzyx(-x,-y,-z,x1,y1,z1,-eta,-psi,-phi)
        write(*,*)"xyz:",sngl(x1),sngl(y1),sngl(z1)
c
        do n=1,49999
          phi=dble(rand())*360d0/180d0*PI
          psi=dble(rand())*360d0/180d0*PI
          eta=dble(rand())*360d0/180d0*PI
          write(*,*)"dgr:",n,sngl(phi/PI*180d0),sngl(psi/PI*180d0)
     &                       ,sngl(eta/PI*180d0)
          do j=1,4
            xx(j)=rand()  !x=[0,1]
            yy(j)=rand()
            zz(j)=rand()
            dd(j)=dsqrt(xx(j)**2+yy(j)**2+zz(j)**2)
            write(*,*)"n:",n,j,sngl(xx(j)),sngl(yy(j)),sngl(zz(j))
            call rotxyz(xx(j),yy(j),zz(j),xn(j),yn(j),zn(j),phi,psi,eta)
            dn(j)=dsqrt(xn(j)**2+yn(j)**2+zn(j)**2)
          enddo
c
          vl1=calvl(xx(1),yy(1),zz(1),xx(2),yy(2),zz(2)
     &             ,xx(3),yy(3),zz(3),xx(4),yy(4),zz(4))
          vl2=calvl(xn(1),yn(1),zn(1),xn(2),yn(2),zn(2)
     &             ,xn(3),yn(3),zn(3),xn(4),yn(4),zn(4))
          write(*,*)"vl:",vl1,vl2
          if(dabs(vl1-vl2).gt.1d-7)stop
          do j=1,4
            if(dabs(dd(j)-dn(j)).gt.1d-7)stop
          enddo
c
c         検算
          do j=1,4
            call rotzyx(xn(j),yn(j),zn(j),x3,y3,z3,-eta,-psi,-phi)
            dx=x3-xx(j)
            dy=y3-yy(j)
            dz=z3-zz(j)
            if(dabs(dx).gt.1d-7 .or. dabs(dy).gt.1d-7 .or. 
     &                               dabs(dz).gt.1d-7)then
              write(*,*)"phi,psi,eta:",sngl(phi/PI*180d0)
     &             ,sngl(psi/PI*180d0),sngl(eta/PI*180d0)
              write(*,*)"dx:",sngl(dx),sngl(dy),sngl(dz)
              stop
            endif
          enddo
        enddo
c
        stop
        end
c------------------------------------------------------------
c
c       座標をずらす
c
c       x軸から反時計周りの角度φ
c       角度:dacos(x0/r0)
c       符号:dasin(z0/r0)/dabs(dasin(z0/r0))
c------------------------------------------------------------
        subroutine trans(x0,y0,z0,x1,y1,z1,dphi,dpsi)
        implicit double precision(a-h,o-z)
        data PI/3.1415926535897932384626433832795/
c
c       x=rcos(φ)cos(ψ)
c       y=rsin(φ)
c       z=rcos(φ)sin(ψ)
        r=dsqrt(x0**2+y0**2+z0**2)  !2乗の代わりにdatan(y/x)を使う。r0=x/sin(φ),r0=y/cos(φ) 0ではない方で計算する
        r0=dsqrt(x0**2+z0**2)
        phi0=dasin(y0/r)
        if(dasin(z0/r0).ne.0d0)
     &    psi0=dacos(x0/r0) * dasin(z0/r0)/dabs(dasin(z0/r0))
        if(dasin(z0/r0).eq.0d0)psi0=dacos(x0/r0)
c
        phi=phi0+dphi   !-90°<=phi<=90°
        psi=psi0+dpsi   !  0°<=psi< 360°×n
        write(*,*)"phi,psi",phi/PI*180d0,psi/PI*180d0,"befr"
        if(phi.gt. PI/2d0)then
          phi= PI/2d0-(phi-PI/2d0)  ! 90°-( 93°-90°)= 87°
          psi=psi+PI
          dphi=-dphi
        endif
        if(phi.lt.-PI/2d0)then
          phi=-PI/2d0-(phi+PI/2d0)  !-90°-(-93°+90°)=-87°
          psi=psi+PI
          dphi=-dphi
        endif
        write(*,*)"phi,psi",phi/PI*180d0,psi/PI*180d0,"afte"
c
        x1=r*dcos(phi)*dcos(psi)
        y1=r*dsin(phi)
        z1=r*dcos(phi)*dsin(psi)
c
        return
        end
c------------------------------------------------------------
c
c       軸の回転(y,z軸の回転:点の移動用)
c
c------------------------------------------------------------
        subroutine rot12(x0,y0,z0,x1,y1,z1,phi,psi)
        implicit double precision(a-h,o-z)
c
c       x1=x0cosφcosψ+y0(-sinφcosψ)+z0(-sinψ)
c       y1=x0sinφ     +y0cosφ        +z0(0)
c       z1=x0cosφsinψ+y0(-sinφsinψ)+z0cosψ
        x1=x0*dcos(phi)*dcos(psi) -y0*dsin(phi)*dcos(psi) -z0*dsin(psi)
        y1=x0*dsin(phi)           +y0*dcos(phi)           +z0*(0d0)
        z1=x0*dcos(phi)*dsin(psi) -y0*dsin(phi)*dsin(psi) +z0*dcos(psi)
c
        return
        end
c
        subroutine rot21(x0,y0,z0,x1,y1,z1,phi,psi)
        implicit double precision(a-h,o-z)
c
c       x1= x0cosψcosφ +y0sinφ +z0sinψcosφ
c       y1=-x0cosψsinφ +y0cosφ -z0sinψsinφ
c       z1=-x0sinψ      +y0(0)   +z0cosψ
        x1= x0*dcos(psi)*dcos(phi) +y0*dsin(phi) +z0*dsin(psi)*dcos(phi)
        y1=-x0*dcos(psi)*dsin(phi) +y0*dcos(phi) -z0*dsin(psi)*dsin(phi)
        z1=-x0*dsin(psi)           +y0*(0d0)     +z0*dcos(psi)
c
        return
        end
c------------------------------------------------------------
c
c       軸の回転(z→y→x軸の回転:面の移動用)
c
c------------------------------------------------------------
        subroutine rot12r(x0,y0,z0,x1,y1,z1,phi,psi,eta)
        implicit double precision(a-h,o-z)
c
        vxx=dcos(phi)*dcos(psi)
        vyx=-dsin(phi)*dcos(eta)*dcos(psi)-dcos(phi)*dsin(psi)*dsin(eta)
        vzx=-dcos(eta)*dsin(psi)
        vxy=dsin(phi)
        vyy=dcos(phi)*dcos(eta)
        vzy=dsin(eta)
        vxz=dcos(phi)*dsin(psi)
        vyz=dcos(phi)*dcos(psi)*dsin(eta)-dsin(phi)*dcos(eta)*dsin(psi)
        vzz=dcos(eta)*dcos(psi)
c
        return
        end
c------------------------------------------------------------
c
c       軸の回転(x→y→z軸の回転:面の移動用)
c
c------------------------------------------------------------
        subroutine rotxyz(x0,y0,z0,x1,y1,z1,sht,phi,psi)
        implicit double precision(a-h,o-z)
c
        vxx= dcos(phi)*dcos(psi)
        vxy= dcos(sht)*dsin(psi)-dsin(sht)*dsin(phi)*dcos(psi)
        vxz=-dsin(sht)*dsin(psi)-dcos(sht)*dsin(phi)*dcos(psi)
        vyx=-dcos(phi)*dsin(psi)
        vyy= dcos(sht)*dcos(psi)+dsin(sht)*dsin(phi)*dsin(psi)
        vyz=-dsin(sht)*dcos(psi)+dcos(sht)*dsin(phi)*dsin(psi)
        vzx= dsin(phi)
        vzy= dsin(sht)*dcos(phi)
        vzz= dcos(sht)*dcos(phi)
c
        x1=x0*vxx+y0*vyx+z0*vzx
        y1=x0*vxy+y0*vyy+z0*vzy
        z1=x0*vxz+y0*vyz+z0*vzz
c
        return
        end
c------------------------------------------------------------
c
c       軸の回転(z→y→x軸の回転:面の移動用)
c
c------------------------------------------------------------
        subroutine rotzyx(x0,y0,z0,x1,y1,z1,sht,phi,psi)
        implicit double precision(a-h,o-z)
c
        vxx= dcos(sht)*dcos(phi)
        vxy= dsin(sht)*dcos(phi)
        vxz=-dsin(phi)

        vyx=-dcos(sht)*dsin(phi)*dsin(psi)-dsin(sht)*dcos(psi)
        vyy=-dsin(sht)*dsin(phi)*dsin(psi)+dcos(sht)*dcos(psi)
        vyz=-dcos(phi)*dsin(psi)

        vzx= dcos(sht)*dsin(phi)*dcos(psi)-dsin(sht)*dsin(psi)
        vzy= dsin(sht)*dsin(phi)*dcos(psi)+dcos(sht)*dsin(psi)
        vzz= dcos(phi)*dcos(psi)
c
        x1=x0*vxx+y0*vyx+z0*vzx
        y1=x0*vxy+y0*vyy+z0*vzy
        z1=x0*vxz+y0*vyz+z0*vzz
c
        return
        end
c------------------------------------------------------------
c
c       3角形の面積
c
c------------------------------------------------------------
        function artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        vx=(y1-y3)*(z2-z3)-(z1-z3)*(y2-y3)
        vy=(z1-z3)*(x2-x3)-(x1-x3)*(z2-z3)
        vz=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)
        artri=dsqrt(vx**2+vy**2+vz**2)/2d0
        return
        end
c------------------------------------------------------------
c
c       4面体の体積
c
c------------------------------------------------------------
        function calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        implicit double precision(a-h,o-z)
c
        a11=1d0
        a21=1d0
        a31=1d0
        a41=1d0
        a12=x1
        a22=x2
        a32=x3
        a42=x4
        a13=y1
        a23=y2
        a33=y3
        a43=y4
        a14=z1
        a24=z2
        a34=z3
        a44=z4
c
        det=a11*a22*a33*a44-a11*a22*a34*a43-a11*a23*a32*a44
     &      +a11*a23*a34*a42
     &      +a11*a24*a32*a43-a11*a24*a33*a42-a12*a21*a33*a44
     &      +a12*a21*a34*a43
     &      +a12*a23*a31*a44-a12*a23*a34*a41-a12*a24*a31*a43
     &      +a12*a24*a33*a41
     &      +a13*a21*a32*a44-a13*a21*a34*a42-a13*a22*a31*a44
     &      +a13*a22*a34*a41
     &      +a13*a24*a31*a42-a13*a24*a32*a41-a14*a21*a32*a43
     &      +a14*a21*a33*a42
     &      +a14*a22*a31*a43-a14*a22*a33*a41-a14*a23*a31*a42
     &      +a14*a23*a32*a41
        calvl=dabs(det)/6d0
c        calvl=det/6d0
c
        return
        end
