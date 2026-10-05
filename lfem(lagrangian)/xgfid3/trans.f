c------------------------------------------------------------
c
c       軸の回転
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
c       検算
c
c------------------------------------------------------------
        subroutine check(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs,ys,zs,eps,ic)
        implicit double precision(a-h,o-z)
c
        dd1=dsqrt((x1-xs)**2+(y1-ys)**2+(z1-zs)**2)
        dd2=dsqrt((x2-xs)**2+(y2-ys)**2+(z2-zs)**2)
        dd3=dsqrt((x3-xs)**2+(y3-ys)**2+(z3-zs)**2)
c
        dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &    +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &    +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c
        ic=0
        if(dd.lt.eps)ic=1
c        if(dd.ge.eps)write(*,*)"check,dd",dd
c        if(dd.ge.eps)stop
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
        vxx= dcos(phi)*dcos(psi)
        vyx=-dsin(phi)*dcos(psi)*dcos(eta)-dsin(psi)*dsin(eta)
        vzx= dsin(phi)*dcos(psi)*dsin(eta)-dsin(psi)*dcos(eta)
        vxy= dsin(phi)
        vyy= dcos(phi)*dcos(eta)
        vzy=-dcos(phi)*dsin(eta)
        vxz= dcos(phi)*dsin(psi)
        vyz=-dsin(phi)*dsin(psi)*dcos(eta)+dcos(psi)*dsin(eta)
        vzz= dsin(phi)*dsin(psi)*dsin(eta)+dcos(psi)*dcos(eta)
c
        x1=x0*vxx+y0*vyx+z0*vzx
        y1=x0*vxy+y0*vyy+z0*vzy
        z1=x0*vxz+y0*vyz+z0*vzz
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
