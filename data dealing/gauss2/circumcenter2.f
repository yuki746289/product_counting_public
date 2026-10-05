c----------------------------------------------------------------------
c
c       4面体の外心
c
c----------------------------------------------------------------------
        subroutine certet2(x1,y1,z1,x2,y2,z2,x3,y3
     &                              ,z3,x4,y4,z4,x10,y10,z10,eps)
        implicit double precision(a-h,o-z)
        data PI/3.1415926535897932384626433832795/
        phi1=30d0/180d0*PI
        psi1=30d0/180d0*PI
        phi2=45d0/180d0*PI
        psi2=45d0/180d0*PI
c
c         3角形の外心の計算
          write(*,*)"面1"
          call cercm2(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs4,ys4,zs4,eps,kch)

          write(*,*)"面2"
          call cercm2(x2,y2,z2,x3,y3,z3,x4,y4,z4,xs1,ys1,zs1,eps,kch)

          write(*,*)"面3"
          call cercm2(x3,y3,z3,x4,y4,z4,x1,y1,z1,xs2,ys2,zs2,eps,kch)

          write(*,*)"面4"
          do i=1,3
            if(i.eq.1)
     &       call cercm2(x4,y4,z4,x1,y1,z1,x2,y2,z2,xs3,ys3,zs3,eps,kch)
            if(i.eq.2 .and. kch.eq.0)then
              call rot12(x4,y4,z4,x40,y40,z40,phi1,psi1)
              call rot12(x1,y1,z1,x10,y10,z10,phi1,psi1)
              call rot12(x2,y2,z2,x20,y20,z20,phi1,psi1)
              call cercm2(x40,y40,z40,x10,y10,z10,x20,y20,z20
     &                                          ,xs30,ys30,zs30,eps,kch)
              call rot21(x40 ,y40 ,z40 ,x4 ,y4 ,z4 ,phi1,psi1)
              call rot21(x10 ,y10 ,z10 ,x1 ,y1 ,z1 ,phi1,psi1)
              call rot21(x20 ,y20 ,z20 ,x2 ,y2 ,z2 ,phi1,psi1)
              call rot21(xs30,ys30,zs30,xs3,ys3,zs3,phi1,psi1)
            endif
            if(i.eq.3 .and. kch.eq.0)then
              call rot12(x4,y4,z4,x40,y40,z40,phi2,psi2)
              call rot12(x1,y1,z1,x10,y10,z10,phi2,psi2)
              call rot12(x2,y2,z2,x20,y20,z20,phi2,psi2)
              call cercm2(x40,y40,z40,x10,y10,z10,x20,y20,z20
     &                                          ,xs30,ys30,zs30,eps,kch)
              call rot21(x40 ,y40 ,z40 ,x4 ,y4 ,z4 ,phi2,psi2)
              call rot21(x10 ,y10 ,z10 ,x1 ,y1 ,z1 ,phi2,psi2)
              call rot21(x20 ,y20 ,z20 ,x2 ,y2 ,z2 ,phi2,psi2)
              call rot21(xs30,ys30,zs30,xs3,ys3,zs3,phi2,psi2)
            endif
          enddo
c
c         4面体の外心の計算
          v1x=x4-x2
          v1y=y4-y2
          v1z=z4-z2
          v2x=x3-x2
          v2y=y3-y2
          v2z=z3-z2
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dn1x,dn1y,dn1z) !外積:|dn1|=1
          v1x=x3-x1
          v1y=y3-y1
          v1z=z3-z1
          v2x=x4-x1
          v2y=y4-y1
          v2z=z4-z1
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dn2x,dn2y,dn2z) !外積:|dn2|=1
          v1x=x4-x1
          v1y=y4-y1
          v1z=z4-z1
          v2x=x2-x1
          v2y=y2-y1
          v2z=z2-z1
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dn3x,dn3y,dn3z) !外積:|dn3|=1
          v1x=x2-x1
          v1y=y2-y1
          v1z=z2-z1
          v2x=x3-x1
          v2y=y3-y1
          v2z=z3-z1
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dn4x,dn4y,dn4z) !外積:|dn4|=1
c
          a1=dn3x
          a2=dn3y
          a3=dn3z
          b1=-dn4x
          b2=-dn4y
          b3=-dn4z
          c1=xs4-xs3
          c2=ys4-ys3
          c3=zs4-zs3
c          write(*,*)"p1",sngl(x1),sngl(y1),sngl(z1)
c          write(*,*)"p2",sngl(x2),sngl(y2),sngl(z2)
c          write(*,*)"p3",sngl(x3),sngl(y3),sngl(z3)
c          write(*,*)"p4",sngl(x4),sngl(y4),sngl(z4)
c          write(*,*)"s1",sngl(xs1),sngl(ys1),sngl(zs1)
c          write(*,*)"s2",sngl(xs2),sngl(ys2),sngl(zs2)
c          write(*,*)"s3",sngl(xs3),sngl(ys3),sngl(zs3)
c          write(*,*)"s4",sngl(xs4),sngl(ys4),sngl(zs4)
c          write(*,*)"dn1",sngl(dn1x),sngl(dn1y),sngl(dn1z)
c          write(*,*)"dn2",sngl(dn2x),sngl(dn2y),sngl(dn2z)
c          write(*,*)"dn3",sngl(dn3x),sngl(dn3y),sngl(dn3z)
c          write(*,*)"dn4",sngl(dn4x),sngl(dn4y),sngl(dn4z)
c          write(*,*)sngl(a1),"k1+",sngl(b1),"k2=",sngl(c1)
c          write(*,*)sngl(a2),"k1+",sngl(b2),"k2=",sngl(c2)
c          write(*,*)sngl(a3),"k1+",sngl(b3),"k2=",sngl(c3)
          if(dabs(a3*b1-a1*b3).gt.eps)then
            dk1=(b1*c3-b3*c1)/(a3*b1-a1*b3)
            dk2=(a3*c1-a1*c3)/(a3*b1-a1*b3)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                         .and. dabs(g3).lt.eps)goto 300
          endif
          if(dabs(a2*b3-a3*b2).gt.eps)then
            dk1=(b3*c2-b2*c3)/(a2*b3-a3*b2)
            dk2=(a2*c3-a3*c2)/(a2*b3-a3*b2)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                         .and. dabs(g3).lt.eps)goto 300
          endif
          if(dabs(a1*b2-a2*b1).gt.eps)then
            dk1=(b2*c1-b1*c2)/(a1*b2-a2*b1)
            dk2=(a1*c2-a2*c1)/(a1*b2-a2*b1)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                         .and. dabs(g3).lt.eps)goto 300
          endif
          write(*,*)"誤差 : ",g1,g2,g3,g4
          stop
300       continue
          x10=xs3+dk1*dn3x
          y10=ys3+dk1*dn3y
          z10=zs3+dk1*dn3z
          x20=xs4+dk2*dn4x
          y20=ys4+dk2*dn4y
          z20=zs4+dk2*dn4z
c          write(*,*)"p30",sngl(dk1),sngl(x10),sngl(y10),sngl(z10)
c          write(*,*)"p40",sngl(dk2),sngl(x20),sngl(y20),sngl(z20)
          g1=dsqrt((x1-x10)**2+(y1-y10)**2+(z1-z10)**2)
          g2=dsqrt((x2-x10)**2+(y2-y10)**2+(z2-z10)**2)
          g3=dsqrt((x3-x10)**2+(y3-y10)**2+(z3-z10)**2)
          g4=dsqrt((x4-x10)**2+(y4-y10)**2+(z4-z10)**2)
          g=dabs(g1-(g1+g2+g3+g4)/4d0)+dabs(g2-(g1+g2+g3+g4)/4d0)
     &     +dabs(g3-(g1+g2+g3+g4)/4d0)+dabs(g4-(g1+g2+g3+g4)/4d0)
          write(*,*)"g",g
          if(g.lt.eps)return    !誤差:1倍
          write(*,*)"g",g
          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3),sngl(g4)
          stop
c
          return
          end




c
c
c-------サブルーチン------------------------------------------------------
c
c





        subroutine cercm2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,eps,kch)
        implicit double precision(a-h,o-z)
        parameter(n=3,m=4)
        dimension a(n,m),iwork(n)
c
          kch=0
c
c         面の式:z=ax+by+c
          call eq31r(x1,y1,z1,x2,y2,z2,x3,y3,z3,eps,ich,h1,h2,h3)
          if(ich.eq.0)write(*,*)"ストップ,ich",ich
          if(ich.eq.0)return
c
          write(*,*)"ich,h:",ich,sngl(h1),sngl(h2),sngl(h3)
c----------------------------------------------------------------------
c
c         外心
c
c----------------------------------------------------------------------
          !------------------------------------------------------------
          !       点:x=x0
          !------------------------------------------------------------
          if(ich.eq.1)then  !x=x0
            x0=h1
c
            cc=y1*z3+y3*z2+y2*z1-y1*z2-y2*z3-y3*z1
            z0=(y1*y3**2+y3*y2**2+y2*y1**2-y1*y2**2-y2*y3**2-y3*y1**2
     &         +y1*z3**2+y3*z2**2+y2*z1**2-y1*z2**2-y2*z3**2-y3*z1**2)
     &         /2d0
            y0=(z1*y2**2+z2*y3**2+z3*y1**2-z1*y3**2-z3*y2**2-z2*y1**2
     &         +z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0
            z0=z0/cc
            y0=y0/cc
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=1",dd
            stop
          endif
          !------------------------------------------------------------
          !       点:y=y0
          !------------------------------------------------------------
          if(ich.eq.2)then  !y=y0
            y0=h1
c
            cc=x1*z3+x3*z2+x2*z1-x1*z2-x2*z3-x3*z1
            z0=(x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &         +x1*z3**2+x3*z2**2+x2*z1**2-x1*z2**2-x2*z3**2-x3*z1**2)
     &         /2d0
            x0=(z1*x2**2+z2*x3**2+z3*x1**2-z1*x3**2-z3*x2**2-z2*x1**2
     &         +z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0
            z0=z0/cc
            x0=x0/cc
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=2",dd
            stop
          endif
          !------------------------------------------------------------
          !       点:z=z0
          !------------------------------------------------------------
          if(ich.eq.3)then  !z=z0
            z0=h1
c
            cc=x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1
            y0=(x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &         +x1*y3**2+x3*y2**2+x2*y1**2-x1*y2**2-x2*y3**2-x3*y1**2)
     &         /2d0
            x0=(y1*x2**2+y2*x3**2+y3*x1**2-y1*x3**2-y3*x2**2-y2*x1**2
     &         +y1*y2**2+y2*y3**2+y3*y1**2-y1*y3**2-y3*y2**2-y2*y1**2)
     &         /2d0
            y0=y0/cc
            x0=x0/cc
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=3",dd
            stop
          endif
          !------------------------------------------------------------
          !       直線:y=ax+b
          !------------------------------------------------------------
          if(ich.eq.11)then
            aa=h1
            bb=h2
c
            cc=(x1*z3+x3*z2+x2*z1-x1*z2-x2*z3-x3*z1)
     &        +(y1*z3+y3*z2+y2*z1-y1*z2-y2*z3-y3*z1)*aa
            z0=(x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &         +x1*y3**2+x3*y2**2+x2*y1**2-x1*y2**2-x2*y3**2-x3*y1**2
     &         +x1*z3**2+x3*z2**2+x2*z1**2-x1*z2**2-x2*z3**2-x3*z1**2)
     &         /2d0
     &         +(y1*x3**2+y3*x2**2+y2*x1**2-y1*x2**2-y2*x3**2-y3*x1**2
     &         + y1*y3**2+y3*y2**2+y2*y1**2-y1*y2**2-y2*y3**2-y3*y1**2
     &         + y1*z3**2+y3*z2**2+y2*z1**2-y1*z2**2-y2*z3**2-y3*z1**2)
     &         /2d0*aa
     &         +(x1*y2+x2*y3+x3*y1-x1*y3-x3*y2-x2*y1)*bb
            x0=(z1*x2**2+z2*x3**2+z3*x1**2-z1*x3**2-z3*x2**2-z2*x1**2
     &         +z1*y2**2+z2*y3**2+z3*y1**2-z1*y3**2-z3*y2**2-z2*y1**2
     &         +z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0
     &         +(z1*y3+z3*y2+z2*y1-z1*y2-z2*y3-z3*y1)*bb
            z0=z0/cc
            x0=x0/cc
            y0=aa*x0+bb
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=11",dd
            stop
          endif
          !------------------------------------------------------------
          !       直線:z=ay+b
          !------------------------------------------------------------
          if(ich.eq.12)then
            aa=h1
            bb=h2
c
            cc=(x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1)
     &        +(x1*z3+x3*z2+x2*z1-x1*z2-x2*z3-x3*z1)*aa
            x0=(y1*x2**2+y2*x3**2+y3*x1**2-y1*x3**2-y3*x2**2-y2*x1**2
     &         +y1*y2**2+y2*y3**2+y3*y1**2-y1*y3**2-y3*y2**2-y2*y1**2
     &         +y1*z2**2+y2*z3**2+y3*z1**2-y1*z3**2-y3*z2**2-y2*z1**2)
     &         /2d0
     &         +(z1*x2**2+z2*x3**2+z3*x1**2-z1*x3**2-z3*x2**2-z2*x1**2
     &         + z1*y2**2+z2*y3**2+z3*y1**2-z1*y3**2-z3*y2**2-z2*y1**2
     &         + z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0*aa
     &         +(y1*z3+y3*z2+y2*z1-y1*z2-y2*z3-y3*z1)*bb
            y0=(x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &         +x1*y3**2+x3*y2**2+x2*y1**2-x1*y2**2-x2*y3**2-x3*y1**2
     &         +x1*z3**2+x3*z2**2+x2*z1**2-x1*z2**2-x2*z3**2-x3*z1**2)
     &         /2d0
     &         +(x1*z2+x2*z3+x3*z1-x1*z3-x3*z2-x2*z1)*bb
            x0=x0/cc
            y0=y0/cc
            z0=aa*y0+bb
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=12",dd
            stop
          endif
          !------------------------------------------------------------
          !       直線:x=az+b
          !------------------------------------------------------------
          if(ich.eq.13)then
            aa=h1
            bb=h2
c
            cc=(y1*z3+y3*z2+y2*z1-y1*z2-y2*z3-y3*z1)
     &        +(y1*x3+y3*x2+y2*x1-y1*x2-y2*x3-y3*x1)*aa
            y0=(z1*x2**2+z2*x3**2+z3*x1**2-z1*x3**2-z3*x2**2-z2*x1**2
     &         +z1*y2**2+z2*y3**2+z3*y1**2-z1*y3**2-z3*y2**2-z2*y1**2
     &         +z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0
     &         +(x1*x2**2+x2*x3**2+x3*x1**2-x1*x3**2-x3*x2**2-x2*x1**2
     &         + x1*y2**2+x2*y3**2+x3*y1**2-x1*y3**2-x3*y2**2-x2*y1**2
     &         + x1*z2**2+x2*z3**2+x3*z1**2-x1*z3**2-x3*z2**2-x2*z1**2)
     &         /2d0*aa
     &         +(z1*x3+z3*x2+z2*x1-z1*x2-z2*x3-z3*x1)*bb
            z0=(y1*x3**2+y3*x2**2+y2*x1**2-y1*x2**2-y2*x3**2-y3*x1**2
     &         +y1*y3**2+y3*y2**2+y2*y1**2-y1*y2**2-y2*y3**2-y3*y1**2
     &         +y1*z3**2+y3*z2**2+y2*z1**2-y1*z2**2-y2*z3**2-y3*z1**2)
     &         /2d0
     &         +(y1*x2+y2*x3+y3*x1-y1*x3-y3*x2-y2*x1)*bb
            y0=y0/cc
            z0=z0/cc
            x0=aa*z0+bb
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=13",dd
            stop
          endif
          !------------------------------------------------------------
          !       面:z=ax+by+c
          !------------------------------------------------------------
          if(ich.eq.20)then
            aa=h1   !z=ax+by+cの係数
            bb=h2
            cc=h3
c
            cf=(x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1)
     &        +(z1*y3+z3*y2+z2*y1-z1*y2-z2*y3-z3*y1)*aa
     &        +(x1*z3+x3*z2+x2*z1-x1*z2-x2*z3-x3*z1)*bb
            x0=(y1*x2**2+y2*x3**2+y3*x1**2-y1*x3**2-y3*x2**2-y2*x1**2
     &         +y1*y2**2+y2*y3**2+y3*y1**2-y1*y3**2-y3*y2**2-y2*y1**2
     &         +y1*z2**2+y2*z3**2+y3*z1**2-y1*z3**2-y3*z2**2-y2*z1**2)
     &         /2d0
     &         +(z1*x2**2+z2*x3**2+z3*x1**2-z1*x3**2-z3*x2**2-z2*x1**2
     &         + z1*y2**2+z2*y3**2+z3*y1**2-z1*y3**2-z3*y2**2-z2*y1**2
     &         + z1*z2**2+z2*z3**2+z3*z1**2-z1*z3**2-z3*z2**2-z2*z1**2)
     &         /2d0*bb
     &         -(y1*z2+y2*z3+y3*z1-y1*z3-y3*z2-y2*z1)*cc
            y0=(x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &         +x1*y3**2+x3*y2**2+x2*y1**2-x1*y2**2-x2*y3**2-x3*y1**2
     &         +x1*z3**2+x3*z2**2+x2*z1**2-x1*z2**2-x2*z3**2-x3*z1**2)
     &         /2d0
     &         +(z1*x3**2+z3*x2**2+z2*x1**2-z1*x2**2-z2*x3**2-z3*x1**2
     &         + z1*y3**2+z3*y2**2+z2*y1**2-z1*y2**2-z2*y3**2-z3*y1**2
     &         + z1*z3**2+z3*z2**2+z2*z1**2-z1*z2**2-z2*z3**2-z3*z1**2)
     &         /2d0*aa
     &         +(x1*z2+x2*z3+x3*z1-x1*z3-x3*z2-x2*z1)*cc
            x0=x0/cf
            y0=y0/cf
            z0=aa*x0+bb*y0+cc
c           検算
            dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
            dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
            dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
            dd=dabs(dabs(dd1)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd2)-dabs(dd1+dd2+dd3)/3d0)
     &        +dabs(dabs(dd3)-dabs(dd1+dd2+dd3)/3d0)
c            write(*,*)"dd",dd
            if(dd.lt.eps)goto 100
            write(*,*)"ich=20",dd
            stop
          endif
c
          return   !kch=0
c
100       continue
          kch=1
          return
          end
c------------------------------------------------------------
c
c       連立方程式を解く
c       ax+by+cz=d
c------------------------------------------------------------
        subroutine eq31r(x1,y1,z1,x2,y2,z2,x3,y3,z3,eps,ich,h1,h2,h3)
        implicit double precision(a-h,o-z)
        parameter(n=3,m=4)
        dimension a(n,m),iwork(n)
        ich=0
c
        h1=-99d0
        h2=-99d0
        h3=-99d0
c        write(*,*)"eq31r"
c------------------------------------------------------------
c
c       点
c
c------------------------------------------------------------
        if(dabs(x1-x2).lt.eps .and. dabs(x1-x3).lt.eps)then
c          write(*,*)"x : ",x1,x2,x3
          ich=1
          h1=(x1+x2+x3)/3d0
          return
        else if(dabs(y1-y2).lt.eps .and. dabs(y1-y3).lt.eps)then
c          write(*,*)"y : ",y1,y2,y3
          ich=2
          h1=(y1+y2+y3)/3d0
          return
        else if(dabs(z1-z2).lt.eps .and. dabs(z1-z3).lt.eps)then
c          write(*,*)"z : ",z1,z2,z3
          ich=3
          h1=(z1+z2+z3)/3d0
          return
        endif
c------------------------------------------------------------
c
c       直線:y=ax+b
c
c------------------------------------------------------------
c        write(*,*)"直線:y=ax+b"
        if(dabs(x1-x3).gt.eps)then !4-1)
          a0=(y1-y3)/(x1-x3)
          b0=(x3*y1-x1*y3)/(x3-x1)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=11
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(x3-x2).gt.eps)then !4-2)
          a0=(y3-y2)/(x3-x2)
          b0=(x2*y3-x3*y2)/(x2-x3)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=11
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(x2-x1).gt.eps)then !4-3)
          a0=(y2-y1)/(x2-x1)
          b0=(x1*y2-x2*y1)/(x1-x2)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=11
            h1=a0
            h2=b0
            return
          endif
        endif
c------------------------------------------------------------
c
c       直線:z=ay+b
c
c------------------------------------------------------------
c        write(*,*)"直線:z=ay+b"
        if(dabs(y1-y3).gt.eps)then
          a0=(z1-z3)/(y1-y3)
          b0=(y3*z1-y1*z3)/(y3-y1)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=12
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(y3-y2).gt.eps)then
          a0=(z3-z2)/(y3-y2)
          b0=(y2*z3-y3*z2)/(y2-y3)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=12
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(y2-y1).gt.eps)then
          a0=(z2-z1)/(y2-y1)
          b0=(y1*z2-y2*z1)/(y1-y2)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=12
            h1=a0
            h2=b0
            return
          endif
        endif
c------------------------------------------------------------
c
c       直線:x=az+b
c
c------------------------------------------------------------
c        write(*,*)"直線:x=az+b"
        if(dabs(z1-z3).gt.eps)then
          a0=(x1-x3)/(z1-z3)
          b0=(z3*x1-z1*x3)/(z3-z1)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=13
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(z3-z2).gt.eps)then
          a0=(x3-x2)/(z3-z2)
          b0=(z2*x3-z3*x2)/(z2-z3)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=13
            h1=a0
            h2=b0
            return
          endif
        endif
        if(dabs(z2-z1).gt.eps)then
          a0=(x2-x1)/(z2-z1)
          b0=(z1*x2-z2*x1)/(z1-z2)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                       .and. dabs(g3).lt.eps)then
            ich=13
            h1=a0
            h2=b0
            return
          endif
        endif
c------------------------------------------------------------
c
c       面(斜め45°):z=a(x+y)+c
c
c------------------------------------------------------------
c        if(dabs((x1+y1)-(x3+y3)).gt.eps)then
c          a2=(z1-z3)/((x1+y1)-(x3+y3))
c          c2=((x1+y1)*z3-(x3+y3)*z1)/((x1+y1)-(x3+y3))
c          g1=a2*(x1+y1)+c2-z1
c          g2=a2*(x2+y2)+c2-z2
c          g3=a2*(x3+y3)+c2-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c        if(dabs((x3+y3)-(x2+y2)).gt.eps)then
c          a1=(z3-z3)/((x3+y3)-(x2+y2))
c          c1=((x3+y3)*z2-(x2+y2)*z3)/((x3+y3)-(x2+y2))
c          g1=a1*(x1+y1)+c1-z1
c          g2=a1*(x2+y2)+c1-z2
c          g3=a1*(x3+y3)+c1-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c        if(dabs((x2+y2)-(x1+y1)).gt.eps)then
c          a1=(z2-z1)/((x2+y2)-(x1+y1))
c          c1=((x2+y2)*z1-(x1+y1)*z2)/((x2+y2)-(x1+y1))
c          g1=a1*(x1+y1)+c1-z1
c          g2=a1*(x2+y2)+c1-z2
c          g3=a1*(x3+y3)+c1-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c------------------------------------------------------------
c
c       面(斜め45°):z=a(x-y)+c
c
c------------------------------------------------------------
c        if(dabs((x1-y1)-(x3-y3)).gt.eps)then
c          a2=(z1-z3)/((x1-y1)-(x3-y3))
c          c2=((x1-y1)*z3-(x3-y3)*z1)/((x1-y1)-(x3-y3))
c          g1=a2*(x1-y1)+c2-z1
c          g2=a2*(x2-y2)+c2-z2
c          g3=a2*(x3-y3)+c2-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c        if(dabs((x3-y3)-(x2-y2)).gt.eps)then
c          a1=(z3-z3)/((x3-y3)-(x2-y2))
c          c1=((x3-y3)*z2-(x2-y2)*z3)/((x3-y3)-(x2-y2))
c          g1=a1*(x1-y1)+c1-z1
c          g2=a1*(x2-y2)+c1-z2
c          g3=a1*(x3-y3)+c1-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c        if(dabs((x2-y2)-(x1-y1)).gt.eps)then
c          a1=(z2-z1)/((x2-y2)-(x1-y1))
c          c1=((x2-y2)*z1-(x1-y1)*z2)/((x2-y2)-(x1-y1))
c          g1=a1*(x1-y1)+c1-z1
c          g2=a1*(x2-y2)+c1-z2
c          g3=a1*(x3-y3)+c1-z3
c          if(dabs(g1).gt.eps .and. dabs(g2).gt.eps
c     &                       .and. dabs(g3).gt.eps)goto 100
c        endif
c------------------------------------------------------------
c
c       面:z=ax+by+c
c
c------------------------------------------------------------
        aa=(z1*y3+z3*y2+z2*y1-z1*y2-z2*y3-z3*y1)
     &    /(x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1)
        bb=(z1*x3+z3*x2+z2*x1-z1*x2-z2*x3-z3*x1)
     &    /(y1*x3+y3*x2+y2*x1-y1*x2-y2*x3-y3*x1)
        cc=(x2*y1*z3+x1*y3*z2+x3*y2*z1-x1*y2*z3-x2*y3*z1-x3*y1*z2)
     &    /(x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1)
c
        g1=x1*aa+y1*bb+cc - z1
        g2=x2*aa+y2*bb+cc - z2
        g3=x3*aa+y3*bb+cc - z3
        g=dabs(g1)+dabs(g2)+dabs(g3)
c        write(*,*)"a,b,c",sngl(aa),sngl(bb),sngl(cc)
        write(*,*)"g",sngl(g1),sngl(g2),sngl(g3)
        if(dabs(g).lt.eps)then
          ich=20
          h1=aa
          h2=bb
          h3=cc
          return
        endif
c
        write(*,*)"p1",sngl(x1),sngl(y1),sngl(z1)
        write(*,*)"p2",sngl(x2),sngl(y2),sngl(z2)
        write(*,*)"p3",sngl(x3),sngl(y3),sngl(z3)
c
        return
        end