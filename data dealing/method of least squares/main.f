        implicit double precision(a-h,o-z)
        parameter(ind=9999,ine=9999)
        dimension xk(ind),yk(ind),zk(ind)
c
          nk=4
          xk(1)=0d0
          yk(1)=0d0
          zk(1)=0d0
          xk(2)=0d0
          yk(2)=5d-1
          zk(2)=1d0
          xk(3)=0d0
          yk(3)=-5d-1
          zk(3)=1d0
          xk(4)=-dsqrt(3d0)
          yk(4)=0d0
          zk(4)=1d0
c
          nk=9
          xk(1)= 0d0
          yk(1)= 0d0
          zk(1)= 0d0
          xk(2)= 0d0
          yk(2)= 1d0
          zk(2)= 1d0
          xk(3)= 0d0
          yk(3)=-1d0
          zk(3)= 1d0
c
          xk(4)= 1d0
          yk(4)= 0d0
          zk(4)= 0d0
          xk(5)= 1d0
          yk(5)= 1d0
          zk(5)= 1d0
          xk(6)= 1d0
          yk(6)=-1d0
          zk(6)= 1d0
c
          xk(7)=-1d0
          yk(7)= 0d0
          zk(7)= 0d0
          xk(8)=-1d0
          yk(8)= 1d0
          zk(8)= 1d0
          xk(9)=-1d0
          yk(9)=-1d0
          zk(9)= 1d0
c
          nk=7
          xk(1)=0.0000000d0
          yk(1)=0.0000000d0
          zk(1)=0.0000000d0
          xk(2)=0.2000000d0
          yk(2)=0.2000000d0
          zk(2)=0.0000000d0
          xk(3)=0.4000000d0
          yk(3)=0.0000000d0
          zk(3)=0.0000000d0
          xk(4)=0.2000000d0
          yk(4)=0.0000000d0
          zk(4)=0.2000000d0
          xk(5)=0.0000000d0
          yk(5)=0.4000000d0
          zk(5)=0.0000000d0
          xk(6)=0.0000000d0
          yk(6)=0.0000000d0
          zk(6)=0.4000000d0
          xk(7)=0.0000000d0
          yk(7)=0.2000000d0
          zk(7)=0.2000000d0
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
c
c         誤差を計算する
          ee=0d0
          do kk=1,nk
            ee=ee+dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
          x0=xk(1)
          y0=yk(1)
          z0=zk(1)
c
          p=2d0*a*x0    !∂z/∂x
          q=2d0*b*y0    !∂z/∂y
          r=2d0*a       !∂2z/∂x2
          s=0d0         !∂2z/∂x∂y
          t=2d0*b       !∂z/∂y2
c         平均曲率:H
          h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
          h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
          h=h1/h2   !平均曲率Η
c         Gaussの曲率:K
          gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c         単位法ベクトルの方向余弦:ex,ey,ez →外向き
          e=dsqrt(1d0+p**2+q**2)
          ex=-p/e
          ey=-q/e
          ez=1d0/e
c
          write(*,*)"a,b,ee:",sngl(a),sngl(b),sngl(ee)
          write(*,*)"h:",sngl(h)
          write(*,*)"ex:",sngl(ex),sngl(ey),sngl(ez)
c
        stop
        end
