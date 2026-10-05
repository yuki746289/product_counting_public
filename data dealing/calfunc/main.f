        implicit double precision(a-h,o-z)
        parameter(ind=9999,ine=9999)
        dimension xk(ind),yk(ind),zk(ind)
c
c
        x0=0d0
        y0=0d0
        z0=0d0
c
c        z=2.56*x**2+3.42*y**2
        nk=1
        xk(nk)=2.32d0
        yk(nk)=1.24d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
        nk=nk+1
        xk(nk)=4.32d0
        yk(nk)=3.24d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
        nk=nk+1
        xk(nk)=-1.32d0
        yk(nk)=2.94d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
        nk=nk+1
        xk(nk)=2.82d0
        yk(nk)=-1.53d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
        nk=nk+1
        xk(nk)=-3.62d0
        yk(nk)=-1.77d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
        nk=nk+1
        xk(nk)=-1.32d0
        yk(nk)=0.24d0
        zk(nk)=2.56d0*xk(nk)**2+3.42d0*yk(nk)**2
c
c         ∑を計算する
          sgx=0d0       !∑x
          sgy=0d0       !∑y
          sgz=0d0       !∑z
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx2y2=0d0    !∑x^2y^2
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          do kk=1,nk
            sgx=sgx+xk(kk)  !∑x
            sgy=sgy+yk(kk)  !∑y
            sgz=sgz+zk(kk)  !∑z
            sgx4=sgx4+xk(kk)**4               !∑x^4
            sgy4=sgy4+yk(kk)**4               !∑y^4
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
            sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
          enddo
c
c         z=ax^2+by^2の係数a,bを計算する
          a1=sgy4*sgx2z-sgx2y2*sgy2z
          a2=sgy4*sgx4-sgx2y2*sgx2y2
          a=a1/a2
          b1=sgx2y2*sgx2z-sgx4*sgy2z
          b2=sgx2y2*sgx2y2-sgx4*sgy4
          b=b1/b2
          write(*,*)"a=",a,"b=",b
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
          p=2d0*a*x0    !∂z/∂x
          q=2d0*b*y0    !∂z/∂y
          r=2d0*a       !∂2z/∂x2
          s=0d0         !∂2z/∂x∂y
          t=2d0*b       !∂z/∂y2
          h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
          h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
          h=h1/h2   !平均曲率Η
          write(*,*)"h=",h
c
          stop
          end
