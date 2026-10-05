c----------------------------------------------------------------------
c
c       4面体の外心
c
c       ～～～訂正～～～
c       2007/04 : 式修正
c----------------------------------------------------------------------
        subroutine certet(x1,y1,z1,x2,y2,z2,x3,y3,z3
     &                              ,x4,y4,z4,x10,y10,z10,r2,eps)    !g1:外接円の半径
        implicit double precision(a-h,o-z)
        data PI/3.1415926535897932384626433832795/
        phi1=30d0/180d0*PI
        psi1=30d0/180d0*PI
        phi2=45d0/180d0*PI
        psi2=45d0/180d0*PI
c
c         3角形の外心の計算
c
c         面1
c          write(*,*)"面1"
          call cercm2(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs4,ys4,zs4,eps,kch)
          do 10 i=1,2   !軸回転させる
c            write(*,*)"i1",i
            if(kch.eq.1)exit
c
            if(i.eq.1)phi=phi1
            if(i.eq.1)psi=psi1
            if(i.eq.2)phi=phi2
            if(i.eq.2)psi=psi2
            call rot12(x1,y1,z1,x10,y10,z10,phi,psi)
            call rot12(x2,y2,z2,x20,y20,z20,phi,psi)
            call rot12(x3,y3,z3,x30,y30,z30,phi,psi)
            call cercm2(x10,y10,z10,x20,y20,z20,x30,y30,z30
     &                                          ,xs40,ys40,zs40,eps,kch)
            call rot21(x10 ,y10 ,z10 ,x1 ,y1 ,z1 ,phi,psi)
            call rot21(x20 ,y20 ,z20 ,x2 ,y2 ,z2 ,phi,psi)
            call rot21(x30 ,y30 ,z30 ,x3 ,y3 ,z3 ,phi,psi)
            call rot21(xs40,ys40,zs40,xs4,ys4,zs4,phi,psi)
            call check(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs4,ys4,zs4,eps,ic) !検算
            if(ic.eq.1)exit
10        continue
          if(kch.eq.0 .and. ic.eq.0)write(*,*)"面1,kch,ic",kch,ic
          if(kch.eq.0 .and. ic.eq.0)stop
c
c         面2
c          write(*,*)"面2"
          call cercm2(x2,y2,z2,x3,y3,z3,x4,y4,z4,xs1,ys1,zs1,eps,kch)
          do 20 i=1,2
c            write(*,*)"i2",i
            if(kch.eq.1)exit
c
            if(i.eq.1)phi=phi1
            if(i.eq.1)psi=psi1
            if(i.eq.2)phi=phi2
            if(i.eq.2)psi=psi2
            call rot12(x2,y2,z2,x20,y20,z20,phi,psi)
            call rot12(x3,y3,z3,x30,y30,z30,phi,psi)
            call rot12(x4,y4,z4,x40,y40,z40,phi,psi)
            call cercm2(x20,y20,z20,x30,y30,z30,x40,y40,z40
     &                                          ,xs10,ys10,zs10,eps,kch)
            call rot21(x20 ,y20 ,z20 ,x2 ,y2 ,z2 ,phi,psi)
            call rot21(x30 ,y30 ,z30 ,x3 ,y3 ,z3 ,phi,psi)
            call rot21(x40 ,y40 ,z40 ,x4 ,y4 ,z4 ,phi,psi)
            call rot21(xs10,ys10,zs10,xs1,ys1,zs1,phi,psi)
            call check(x2,y2,z2,x3,y3,z3,x4,y4,z4,xs1,ys1,zs1,eps,ic)   !検算
            if(ic.eq.1)exit
20        continue
          if(kch.eq.0 .and. ic.eq.0)write(*,*)"面2,kch,ic",kch,ic
          if(kch.eq.0 .and. ic.eq.0)stop
c
c         面3
c          write(*,*)"面3"
          call cercm2(x3,y3,z3,x4,y4,z4,x1,y1,z1,xs2,ys2,zs2,eps,kch)
          do 30 i=1,2
c            write(*,*)"i3",i
            if(kch.eq.1)exit
c
            if(i.eq.1)phi=phi1
            if(i.eq.1)psi=psi1
            if(i.eq.2)phi=phi2
            if(i.eq.2)psi=psi2
            call rot12(x3,y3,z3,x30,y30,z30,phi,psi)
            call rot12(x4,y4,z4,x40,y40,z40,phi,psi)
            call rot12(x1,y1,z1,x10,y10,z10,phi,psi)
            call cercm2(x30,y30,z30,x40,y40,z40,x10,y10,z10
     &                                          ,xs20,ys20,zs20,eps,kch)
            call rot21(x30 ,y30 ,z30 ,x3 ,y3 ,z3 ,phi,psi)
            call rot21(x40 ,y40 ,z40 ,x4 ,y4 ,z4 ,phi,psi)
            call rot21(x10 ,y10 ,z10 ,x1 ,y1 ,z1 ,phi,psi)
            call rot21(xs20,ys20,zs20,xs2,ys2,zs2,phi,psi)
            call check(x3,y3,z3,x4,y4,z4,x1,y1,z1,xs2,ys2,zs2,eps,ic)   !検算
            if(ic.eq.1)exit
30        continue
          if(kch.eq.0 .and. ic.eq.0)write(*,*)"面3,kch,ic",kch,ic
          if(kch.eq.0 .and. ic.eq.0)stop
c
c         面4
c          write(*,*)"面4"
          call cercm2(x4,y4,z4,x1,y1,z1,x2,y2,z2,xs3,ys3,zs3,eps,kch)
          do 40 i=1,2
c            write(*,*)"i4",i
            if(kch.eq.1)exit
c
            if(i.eq.1)phi=phi1
            if(i.eq.1)psi=psi1
            if(i.eq.2)phi=phi2
            if(i.eq.2)psi=psi2
            call rot12(x4,y4,z4,x40,y40,z40,phi,psi)
            call rot12(x1,y1,z1,x10,y10,z10,phi,psi)
            call rot12(x2,y2,z2,x20,y20,z20,phi,psi)
            call cercm2(x40,y40,z40,x10,y10,z10,x20,y20,z20
     &                                          ,xs30,ys30,zs30,eps,kch)
            call rot21(x40 ,y40 ,z40 ,x4 ,y4 ,z4 ,phi,psi)
            call rot21(x10 ,y10 ,z10 ,x1 ,y1 ,z1 ,phi,psi)
            call rot21(x20 ,y20 ,z20 ,x2 ,y2 ,z2 ,phi,psi)
            call rot21(xs30,ys30,zs30,xs3,ys3,zs3,phi,psi)
            call check(x4,y4,z4,x1,y1,z1,x2,y2,z2,xs3,ys3,zs3,eps,ic)   !検算
            if(ic.eq.1)exit
40        continue
          if(kch.eq.0 .and. ic.eq.0)write(*,*)"面4,kch,ic",kch,ic
          if(kch.eq.0 .and. ic.eq.0)stop
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
c
          if(dabs(a3*b1-a1*b3).gt.eps)then
            dk1=(b1*c3-b3*c1)/(a3*b1-a1*b3)
            dk2=(a3*c1-a1*c3)/(a3*b1-a1*b3)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            g=sngl(g1)+sngl(g2)+sngl(g3)
c            write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
            if(dabs(g).lt.eps)goto 300
          endif
          if(dabs(a2*b3-a3*b2).gt.eps)then
            dk1=(b3*c2-b2*c3)/(a2*b3-a3*b2)
            dk2=(a2*c3-a3*c2)/(a2*b3-a3*b2)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            g=sngl(g1)+sngl(g2)+sngl(g3)
c            write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
            if(dabs(g).lt.eps)goto 300
          endif
          if(dabs(a1*b2-a2*b1).gt.eps)then
            dk1=(b2*c1-b1*c2)/(a1*b2-a2*b1)
            dk2=(a1*c2-a2*c1)/(a1*b2-a2*b1)
            g1=a1*dk1+b1*dk2-c1
            g2=a2*dk1+b2*dk2-c2
            g3=a3*dk1+b3*dk2-c3
            g=sngl(g1)+sngl(g2)+sngl(g3)
c            write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
            if(dabs(g).lt.eps)goto 300
          endif
          write(*,*)"gg",g
          stop
c
300       continue
          x10=xs3+dk1*dn3x
          y10=ys3+dk1*dn3y
          z10=zs3+dk1*dn3z
          x20=xs4+dk2*dn4x
          y20=ys4+dk2*dn4y
          z20=zs4+dk2*dn4z
c          write(*,*)"p10",sngl(dk1),sngl(x10),sngl(y10),sngl(z10)
c          write(*,*)"p20",sngl(dk2),sngl(x20),sngl(y20),sngl(z20)
          g1=(x1-x10)**2+(y1-y10)**2+(z1-z10)**2
          g2=(x2-x10)**2+(y2-y10)**2+(z2-z10)**2
          g3=(x3-x10)**2+(y3-y10)**2+(z3-z10)**2
          g4=(x4-x10)**2+(y4-y10)**2+(z4-z10)**2
          g=dabs(g1-(g1+g2+g3+g4)/4d0)+dabs(g2-(g1+g2+g3+g4)/4d0)
     &     +dabs(g3-(g1+g2+g3+g4)/4d0)+dabs(g4-(g1+g2+g3+g4)/4d0)
          r2=(g1+g2+g3+g4)/4d0
c          write(*,*)"g",g
          if(g.lt.eps*1d1)return    !誤差:10倍
c
          g1=(x1-x20)**2+(y1-y20)**2+(z1-z20)**2
          g2=(x2-x20)**2+(y2-y20)**2+(z2-z20)**2
          g3=(x3-x20)**2+(y3-y20)**2+(z3-z20)**2
          g4=(x4-x20)**2+(y4-y20)**2+(z4-z20)**2
          g=dabs(g1-(g1+g2+g3+g4)/4d0)+dabs(g2-(g1+g2+g3+g4)/4d0)
     &     +dabs(g3-(g1+g2+g3+g4)/4d0)+dabs(g4-(g1+g2+g3+g4)/4d0)
          r2=(g1+g2+g3+g4)/4d0
c          write(*,*)"g",g
          if(g.lt.eps*1d1)return    !誤差:10倍

          write(*,*)"g",g
c          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3),sngl(g4)
          stop
c
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
c          write(*,*)"ich,h:",ich,sngl(h1),sngl(h2),sngl(h3)
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
            return   !kch=0
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
            return   !kch=0
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
            return   !kch=0
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
            return   !kch=0
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
            return   !kch=0
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
            return   !kch=0
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
            if(dd.lt.eps)goto 100
            write(*,*)"ich=20",dd
            return   !kch=0
          endif
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
c       誤差なくす
c
c------------------------------------------------------------
c        x10=x1
c        y10=y1
c        z10=z1
c        x20=x2
c        y20=y2
c        z20=z2
c        x30=x3
c        y30=y3
c        z30=z3
c        x1=x1-mod(x1,1d-7)
c        y1=y1-mod(y1,1d-7)
c        z1=z1-mod(z1,1d-7)
c        x2=x2-mod(x2,1d-7)
c        y2=y2-mod(y2,1d-7)
c        z2=z2-mod(z2,1d-7)
c        x3=x3-mod(x3,1d-7)
c        y3=y3-mod(y3,1d-7)
c        z3=z3-mod(z3,1d-7)
c------------------------------------------------------------
c
c       点
c
c------------------------------------------------------------
c        write(*,*)"点"
        if(dabs(x1-x2).lt.eps .and. dabs(x1-x3).lt.eps)then
c          write(*,*)"x : ",x1,x2,x3
          ich=1
          h1=(x1+x2+x3)/3d0
c          call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
          return
        else if(dabs(y1-y2).lt.eps .and. dabs(y1-y3).lt.eps)then
c          write(*,*)"y : ",y1,y2,y3
          ich=2
          h1=(y1+y2+y3)/3d0
c          call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
          return
        else if(dabs(z1-z2).lt.eps .and. dabs(z1-z3).lt.eps)then
c          write(*,*)"z : ",z1,z2,z3
          ich=3
          h1=(z1+z2+z3)/3d0
c          call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
          return
        endif
c------------------------------------------------------------
c
c       直線:y=ax+b
c
c------------------------------------------------------------
c        write(*,*)"直線:y=ax+b"
c        write(*,*)"x",sngl(dabs(x1-x3)),sngl(dabs(x3-x2))
c     &                                  ,sngl(dabs(x2-x1))
        if(dabs(x1-x3).gt.eps)then !4-1)
          a0=(y1-y3)/(x1-x3)
          b0=(x3*y1-x1*y3)/(x3-x1)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=11
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(x3-x2).gt.eps)then !4-2)
          a0=(y3-y2)/(x3-x2)
          b0=(x2*y3-x3*y2)/(x2-x3)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=11
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(x2-x1).gt.eps)then !4-3)
          a0=(y2-y1)/(x2-x1)
          b0=(x1*y2-x2*y1)/(x1-x2)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=11
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
c------------------------------------------------------------
c
c       直線:z=ay+b
c
c------------------------------------------------------------
c        write(*,*)"直線:z=ay+b"
c        write(*,*)"y",sngl(dabs(y1-y3)),sngl(dabs(y3-y2))
c     &                                  ,sngl(dabs(y2-y1))
        if(dabs(y1-y3).gt.eps)then
          a0=(z1-z3)/(y1-y3)
          b0=(y3*z1-y1*z3)/(y3-y1)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=12
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(y3-y2).gt.eps)then
          a0=(z3-z2)/(y3-y2)
          b0=(y2*z3-y3*z2)/(y2-y3)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=12
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(y2-y1).gt.eps)then
          a0=(z2-z1)/(y2-y1)
          b0=(y1*z2-y2*z1)/(y1-y2)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=12
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
c------------------------------------------------------------
c
c       直線:x=az+b
c
c------------------------------------------------------------
c        write(*,*)"直線:x=az+b"
c        write(*,*)"z",sngl(dabs(z1-z3)),sngl(dabs(z3-z2))
c     &                                  ,sngl(dabs(z2-z1))
        if(dabs(z1-z3).gt.eps)then
          a0=(x1-x3)/(z1-z3)
          b0=(z3*x1-z1*x3)/(z3-z1)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=13
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(z3-z2).gt.eps)then
          a0=(x3-x2)/(z3-z2)
          b0=(z2*x3-z3*x2)/(z2-z3)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=13
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
            return
          endif
        endif
        if(dabs(z2-z1).gt.eps)then
          a0=(x2-x1)/(z2-z1)
          b0=(z1*x2-z2*x1)/(z1-z2)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
          g=dabs(g1)+dabs(g2)+dabs(g3)
c          write(*,*)"a0,b0",sngl(a0),sngl(b0)
c          write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
          if(dabs(g).lt.eps)then !誤差1倍
            ich=13
            h1=a0
            h2=b0
c            call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
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
c        write(*,*)"面:z=ax+by+c"
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
c        write(*,*)"eps < g",sngl(eps),sngl(g1),sngl(g2),sngl(g3)
c        write(*,*)"a,b,c",sngl(aa),sngl(bb),sngl(cc)
        if(dabs(g).lt.eps)then !誤差1倍
          ich=20
          h1=aa
          h2=bb
          h3=cc
c          call copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10,z10,x20,y20,z20
c     &                                                     ,x30,y30,z30)
          return
        endif
c
        write(*,*)"p1:",x1,y1,z1
        write(*,*)"p2:",x2,y2,z2
        write(*,*)"p3:",x3,y3,z3
c
        return
        end
c
        subroutine copy(x1,y1,z1,x2,y2,z2,x3,y3,z3,x10,y10
     &                           ,z10,x20,y20,z20,x30,y30,z30)
        implicit double precision(a-h,o-z)
c
        x1=x10
        y1=y10
        z1=z10
        x2=x20
        y2=y20
        z2=z20
        x3=x30
        y3=y30
        z3=z30
        return
        end
c------------------------------------------------------------
c
c       ガウスの消去法(ピボット選択)
c
c------------------------------------------------------------
        subroutine sweep(a,n,m,eps,iwork,ill)
        implicit double precision(a-h,o-z)
        dimension a(n,m),iwork(n)
        ill=0
c
        do 10 i=1,n
10      iwork(i)=i
c
        do 17 k=1,n
c
c         行列(k-n, k-n)の係数aij の最大値aij,max
c         i=ip, j=iq
          dmax=dabs(a(k,k))
          ip=k
          iq=k
          do 11 j=k,n
          do 11 i=k,n
            if( dmax.ge.dabs(a(i,j)) )goto 11
            dmax=dabs(a(i,j))
            ip=i
            iq=j
11        continue
c
c         係数が小さかった時
          if(dmax.le.eps) goto 20
c
c         k列とiq列の入れ換え
c
          do 12 i=1,n
            w=a(i,k)
            a(i,k)=a(i,iq)
12          a(i,iq)=w
c
c         k行とip行の入れ換え
c
          do 13 j=k,m
            w=a(k,j)
            a(k,j)=a(ip,j)
13          a(ip,j)=w
c
c         iwork(入れ換え後の列)=入れ換え前の列
c         iwork(入れ換え前の列)=入れ換え後の列
c
          i=iwork(k)
          iwork(k)=iwork(iq)
          iwork(iq)=i
c
c         akj=akj / akk  (akk=1d0になる)
          do 14 j=k+1,m
14        a(k,j)=a(k,j)/a(k,k)
c
c         注目している行を引いていく
c         aij=aij - aik*(akj/akk)
          do 16 i=1,n   !行
            if(i.eq.k)goto 16  !注目している行は飛ばす(「他の行」から「注目している行」を引く)
            do 15 j=k+1,m       !列 : 行列の右上3角形だけ考える。左下3角形の係数は0になる
15          a(i,j)=a(i,j)-a(i,k)*a(k,j)     !aij=aij - aik*(akj/akk)
16        continue
17      continue
c
c       列を元に戻す
        do 19 j=n+1,m    !j=4+1,5
        do 18 i=1,n
          iw=iwork(i)
18        a(iw,n)=a(i,j)
          do 19 i=1,n
19      a(i,j)=a(i,n)    !a(i,5)=a(i,4)
        return
c
20      write(6,100)
        ill=1
        return
100     format(1H ,3X,'MATRIX IS ILL')
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