        implicit double precision(a-h,o-z)
c        goto 100
c
c       4面体の座標
        x1=0d0
        y1=0d0
        z1=0d0

        x2=1d0
        y2=0d0
        z2=0d0

        x3=0d0
        y3=1d0
        z3=0d0

        x4=0d0
        y4=0d0
        z4=1d0
c
c        call circumct(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c        call circumct2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
        call circumct3(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c
c       表示
        d1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        d2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
        d3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
        d4=dsqrt((x4-x0)**2+(y4-y0)**2+(z4-z0)**2)
        write(*,*)"x0,y0,z0:",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"d1",d1
        write(*,*)"d2",d2
        write(*,*)"d3",d3
        write(*,*)"d4",d4
c
c       4面体の座標
        x1=0.12d0
        y1=0.16d0
        z1=0.11d0
        x2=1.85d0
        y2=0.13d0
        z2=0.15d0
        x3=0.85d0
        y3=0.12d0
        z3=1.10d0
        x4=0.64d0
        y4=1.19d0
        z4=0.74d0
c
c        call circumct(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c        call circumct2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
        call circumct3(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c
c       表示
        d1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        d2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
        d3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
        d4=dsqrt((x4-x0)**2+(y4-y0)**2+(z4-z0)**2)
        write(*,*)"x0,y0,z0:",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"d1",d1
        write(*,*)"d2",d2
        write(*,*)"d3",d3
        write(*,*)"d4",d4
c
c       4面体の座標
        x1=-10d0
        y1=-10d0
        z1=-10d0
        x2= 10d0
        y2=-10d0
        z2=-10d0
        x3=  0d0
        y3=  0d0
        z3= 10d0
        x4=  0d0
        y4= -1d0
        z4= -1d0
c
c        call circumct(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c        call circumct2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
        call circumct3(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c
c       表示
        d1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        d2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
        d3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
        d4=dsqrt((x4-x0)**2+(y4-y0)**2+(z4-z0)**2)
        write(*,*)"x0,y0,z0:",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"d1",d1
        write(*,*)"d2",d2
        write(*,*)"d3",d3
        write(*,*)"d4",d4
c
c       4面体の座標
100     x1= 0d0
        y1=-1d0
        z1=-1d0
        x2= 0d0
        y2=10d0
        z2= 0d0
        x3= 0d0
        y3= 0d0
        z3=10d0
        x4=-5d-1
        y4=-1d0
        z4=-1d0
c
c        call circumct(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c        call circumct2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
        call circumct3(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0)
c
c       表示
        d1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
        d2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
        d3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
        d4=dsqrt((x4-x0)**2+(y4-y0)**2+(z4-z0)**2)
        write(*,*)"x0,y0,z0:",sngl(x0),sngl(y0),sngl(z0)
        write(*,*)"d1",d1
        write(*,*)"d2",d2
        write(*,*)"d3",d3
        write(*,*)"d4",d4
c
        stop
        end
c--------------------------------------------------------------
c
c       4面体の外心を計算する
c
c--------------------------------------------------------------
        subroutine circumct(x1,y1,z1,x2,y2,z2
     &                     ,x3,y3,z3,x4,y4,z4,x0,y0,z0)
        implicit double precision(a-h,o-z)
c
c       面の外心(circumcenter)の計算
        call sfcr3(x4,y4,z4,x1,y1,z1,x2,y2,z2,xs1,ys1,zs1)  !sfcr([点1],[点2],[点3],[外心])
        call sfcr3(x4,y4,z4,x2,y2,z2,x3,y3,z3,xs2,ys2,zs2)
        call sfcr3(x4,y4,z4,x3,y3,z3,x1,y1,z1,xs3,ys3,zs3)
c
c       4面体の外心の計算
        call vlcr2(x4,y4,z4,xs1,ys1,zs1
     &                     ,xs2,ys2,zs2,xs3,ys3,zs3,x0,y0,z0)

        return
        end
c--------------------------------------------------------------
c
c       面の外心から4面体の外心を計算する部分
c
c--------------------------------------------------------------
        subroutine vlcr2(x4,y4,z4,xs1,ys1,zs1
     &                     ,xs2,ys2,zs2,xs3,ys3,zs3,x0,y0,z0)
        implicit double precision(a-h,o-z)
c
        xt1=xs1-x4
        yt1=ys1-y4
        zt1=zs1-z4

        xt2=xs2-x4
        yt2=ys2-y4
        zt2=zs2-z4

        xt3=xs3-x4
        yt3=ys3-y4
        zt3=zs3-z4
        write(*,*)"xt1",sngl(xt1),sngl(yt1),sngl(zt1)
        write(*,*)"xt2",sngl(xt2),sngl(yt2),sngl(zt2)
        write(*,*)"xt3",sngl(xt3),sngl(yt3),sngl(zt3)
c
        a1=xt1
        b1=yt1
        c1=zt1
        d1=xt1**2+yt1**2+zt1**2

        a2=xt2
        b2=yt2
        c2=zt2
        d2=xt2**2+yt2**2+zt2**2

        a3=xt3
        b3=yt3
        c3=zt3
        d3=xt3**2+yt3**2+zt3**2

        ee1=(yt1-yt3)*yt2*yt3
        ee2=xt2*(yt3*yt3-yt1*yt1)
        ee=a1*b2*c3+a2*b3*c1+a3*b1*c2-a2*b1*c3-a1*b3*c2-a3*b2*c1
        x0=b1*c2*d3+b2*c3*d1+b3*c1*d2-d2*b1*c3-d1*b3*c2-d3*b2*c1
        y0=c1*a2*d3+c2*a3*d1+c3*a1*d2-d2*c1*a3-d1*c3*a2-d3*c2*a1
        z0=a1*b2*d3+a2*b3*d1+a3*b1*d2-d2*a1*b3-d1*a3*b2-d3*a2*b1
        x0=x0/ee
        y0=y0/ee
        z0=z0/ee
c
c        call sfcr3(xt1,yt1,zt1,xt2,yt2,zt2,xt3,yt3,zt3,x0,y0,z0)
c
        x0=x0+x4
        y0=y0+y4
        z0=z0+z4

        return
        end
c--------------------------------------------------------------
c
c       3角形の外心を計算する
c
c--------------------------------------------------------------
        subroutine sfcr3(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs,ys,zs)  !sfcr([点1],[点2],[点3],[外心])
        implicit double precision(a-h,o-z)
c
c       p1(x1,y1,z1)を原点に持って来る
        xd1=x1-x1
        yd1=y1-y1
        zd1=z1-z1
        xd2=x2-x1
        yd2=y2-y1
        zd2=z2-z1
        xd3=x3-x1
        yd3=y3-y1
        zd3=z3-z1
c
c       面の方程式(ax+by+cz=0)と内積の条件2つ
        a1=yd2*zd3-yd3*zd2
        b1=zd2*xd3-zd3*xd2
        c1=xd2*yd3-xd3*yd2
        d1=0d0
        a2=xd2
        b2=yd2
        c2=zd2
        d2=(xd2*xd2+yd2*yd2+zd2*zd2)/2d0
        a3=xd3
        b3=yd3
        c3=zd3
        d3=(xd3*xd3+yd3*yd3+zd3*zd3)/2d0
c
c       連立方程式の解
        ee=a1*b2*c3+a2*b3*c1+a3*b1*c2-a2*b1*c3-a1*b3*c2-a3*b2*c1
        xs=b1*c2*d3+b2*c3*d1+b3*c1*d2-d2*b1*c3-d1*b3*c2-d3*b2*c1
        ys=c1*a2*d3+c2*a3*d1+c3*a1*d2-d2*c1*a3-d1*c3*a2-d3*c2*a1
        zs=a1*b2*d3+a2*b3*d1+a3*b1*d2-d2*a1*b3-d1*a3*b2-d3*a2*b1
        xs=xs/ee
        ys=ys/ee
        zs=zs/ee
c
c       p1(x1,y1,z1)を元の位置に戻す
        xs=xs+x1
        ys=ys+y1
        zs=zs+z1
c
        f1=dsqrt((x1-xs)**2+(y1-ys)**2+(z1-zs)**2)
        f2=dsqrt((x2-xs)**2+(y2-ys)**2+(z2-zs)**2)
        f3=dsqrt((x3-xs)**2+(y3-ys)**2+(z3-zs)**2)
c        write(*,*)"x1",sngl(x1),sngl(y1),sngl(z1)
c        write(*,*)"x2",sngl(x2),sngl(y2),sngl(z2)
c        write(*,*)"x3",sngl(x3),sngl(y3),sngl(z3)
        write(*,*)"f1,f2,f3",sngl(f1),sngl(f2),sngl(f3)

        return
        end
c--------------------------------------------------------------
c
c       面の外心から4面体の外心を計算する部分
c
c--------------------------------------------------------------
        subroutine vlcr(x4,y4,z4,xs1,ys1,zs1
     &                     ,xs2,ys2,zs2,xs3,ys3,zs3,x0,y0,z0)
        implicit double precision(a-h,o-z)
c
        cc1=(xs1-x4)*xs1+(ys1-y4)*ys1+(zs1-z4)*zs1
        cc2=(xs2-x4)*xs2+(ys2-y4)*ys2+(zs2-z4)*zs2
        cc3=(xs3-x4)*xs3+(ys3-y4)*ys3+(zs3-z4)*zs3
        c11=(xs1-x4)/cc1
        c12=(ys1-y4)/cc1
        c13=(zs1-z4)/cc1
        c21=(xs2-x4)/cc2
        c22=(ys2-y4)/cc2
        c23=(zs2-z4)/cc2
        c31=(xs3-x4)/cc3
        c32=(ys3-y4)/cc3
        c33=(zs3-z4)/cc3
        call eq3(c11,c12,c13,c21,c22,c23,c31,c32,c33,x0,y0,z0)

        return
        end
c--------------------------------------------------------------
c
c       3角形の外心を計算する
c
c--------------------------------------------------------------
        subroutine sfcr2(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs,ys,zs)  !sfcr([点1],[点2],[点3],[外心])
        implicit double precision(a-h,o-z)

c
c       外心の計算
        a1=x2-x1
        b1=y2-y1
        c1=z2-z1
        d1=(x2**2-x1**2+y2**2-y1**2+z2**2-z1**2)/2d0

        a2=x3-x2
        b2=y3-y2
        c2=z3-z2
        d2=(x3**2-x2**2+y3**2-y2**2+z3**2-z2**2)/2d0

        a3=x1-x3
        b3=y1-y3
        c3=z1-z3
        d3=(x1**2-x3**2+y1**2-y3**2+z1**2-z3**2)/2d0

        xyz=a1*b2*c3+a2*b3*c1+a3*b1*c2-a2*b1*c3-a1*b3*c2-a3*b2*c1
        xs=d1*(b1*c2+b2*c3+b3*c1-b1*c3-b3*c2-b2*c1)/xyz
        ys=d2*(c1*a2+c2*a3+c3*a1-c1*a3-c3*a2-c2*a1)/xyz
        zs=d3*(a1*b2+a2*b3+a3*b1-a1*b3-a3*b2-a2*b1)/xyz

        return
        end
c--------------------------------------------------------------
c
c       3角形の外心を計算する
c
c--------------------------------------------------------------
        subroutine sfcr(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs,ys,zs)  !sfcr([点1],[点2],[点3],[外心])
        implicit double precision(a-h,o-z)
c
c       x1=x2=x3=constの時
c
        if(dabs(x1-x2).lt.1d-10 .and. dabs(x2-x3).lt.1d-10
     &                          .and. dabs(x3-x1).lt.1d-10)then
          xs=(x1+x2+x3)/3d0
          call sfcr2d(y1,z1,y2,z2,y3,z3,ys,zs)
c
c       y1=y2=y3=constの時
c
        else if(dabs(y1-y2).lt.1d-10 .and. dabs(y2-y3).lt.1d-10
     &                                .and. dabs(y3-y1).lt.1d-10)then
          ys=(y1+y2+y3)/3d0
          call sfcr2d(z1,x1,z2,x2,z3,x3,zs,xs)
c
c       z1=z2=z3=constの時
c
        else if(dabs(z1-z2).lt.1d-10 .and. dabs(z2-z3).lt.1d-10
     &                                .and. dabs(z3-z1).lt.1d-10)then
        zs=(z1+z2+z3)/3d0
        call sfcr2d(x1,y1,x2,y2,x3,y3,xs,ys)
c
c       それ以外
c
        else
c
c         面:ax+by+cz=1の係数a,b,cの計算
          call eq3(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c)
c
c         外心の計算
          x21=(x2+x1)/2d0
          y21=(y2+y1)/2d0
          z21=(z2+z1)/2d0
          x31=(x3+x1)/2d0
          y31=(y3+y1)/2d0
          z31=(z3+z1)/2d0
c
          cc2=(x2-x1)*x21+(y2-y1)*y21+(z2-z1)*z21
          cc3=(x3-x1)*x31+(y3-y1)*y31+(z3-z1)*z31
          c11=a
          c12=b
          c13=c
          c21=(x2-x1)/cc2
          c22=(y2-y1)/cc2
          c23=(z2-z1)/cc2
          c31=(x3-x1)/cc3
          c32=(y3-y1)/cc3
          c33=(z3-z1)/cc3
          call eq3(c11,c12,c13,c21,c22,c23,c31,c32,c33,xs,ys,zs)
        endif

        return
        end
c--------------------------------------------------------------
c
c       2次元で3角形の外心を計算する
c
c--------------------------------------------------------------
        subroutine sfcr2d(x1,y1,x2,y2,x3,y3,xs,ys)  !sfcr([点1],[点2],[外心])
        implicit double precision(a-h,o-z)

        x1=y1*y3**2+y3*y2**2+y2*y1**2-y1*y2**2-y2*y3**2-y3*y1**2
     &    +y1*x3**2+y3*x2**2+y2*x1**2-y1*x2**2-y2*x3**2-y3*x1**2
        x2=2d0*(x1*y3+x3*y2+x2*y1-x1*y2-x2*y3-x3*y1)
        xs=x1/x2
c
        y1=x1*x3**2+x3*x2**2+x2*x1**2-x1*x2**2-x2*x3**2-x3*x1**2
     &    +x1*y3**2+x3*y2**2+x2*y1**2-x1*y2**2-x2*y3**2-x3*y1**2
        y2=2d0*(y1*x3+y3*x2+y2*x1-y1*x2-y2*x3-y3*x1)
        ys=y1/y2

        return
        end
c--------------------------------------------------------------
c
c       3元1次方程式の係数の計算
c         ax1+by1+cz1=1
c         ax2+by2+cz2=1
c         ax3+by3+cz3=1
c
c--------------------------------------------------------------
        subroutine eq3(x1,y1,z1,x2,y2,z2,x3,y3,z3,a,b,c)
        implicit double precision(a-h,o-z)

        d=x1*y2*z3+x2*y3*z1+x3*y1*z2-x2*y1*z3-x1*y3*z2-x3*y2*z1
        a=(y1*z2+y2*z3+y3*z1-y1*z3-y3*z2-y2*z1)/d
        b=(z1*x2+z2*x3+z3*x1-z1*x3-z3*x2-z2*x1)/d
        c=(x1*y2+x2*y3+x3*y1-x1*y3-x3*y2-x2*y1)/d

        return
        end