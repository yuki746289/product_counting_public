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