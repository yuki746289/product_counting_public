c----------------------------------------------------------------------
c
c       4面体の外心
c
c----------------------------------------------------------------------
        subroutine certet(x1,y1,z1,x2,y2,z2,x3,y3
     &                              ,z3,x4,y4,z4,x10,y10,z10,eps)
        implicit double precision(a-h,o-z)
c
c         3角形の外心の計算
          call cercm(x1,y1,z1,x2,y2,z2,x3,y3,z3,xs4,ys4,zs4,eps)
          call cercm(x2,y2,z2,x3,y3,z3,x4,y4,z4,xs1,ys1,zs1,eps)
          call cercm(x3,y3,z3,x4,y4,z4,x1,y1,z1,xs2,ys2,zs2,eps)
          call cercm(x4,y4,z4,x1,y1,z1,x2,y2,z2,xs3,ys3,zs3,eps)
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
          write(*,*)"g",sngl(g1),sngl(g2),sngl(g3),sngl(g4)
          if(dabs(g1-g2).lt.eps .and. dabs(g2-g3).lt.eps
     &       .and. dabs(g3-g4).lt.eps .and. dabs(g4-g1).lt.eps)return
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





        subroutine cercm(x1,y1,z1,x2,y2,z2,x3,y3,z3,x0,y0,z0,eps)
        implicit double precision(a-h,o-z)
        parameter(n=3,m=4)
        dimension a(n,m),iwork(n)
c
          x13=(x1+x3)/2d0
          x32=(x3+x2)/2d0
          x21=(x2+x1)/2d0
          y13=(y1+y3)/2d0
          y32=(y3+y2)/2d0
          y21=(y2+y1)/2d0
          z13=(z1+z3)/2d0
          z32=(z3+z2)/2d0
          z21=(z2+z1)/2d0
c
c         面の式:z=ax+by+c
          call eq31(x1,y1,z1,x2,y2,z2,x3,y3,z3,eps,ich,h1,h2,h3)
          write(*,*)"ich,h:",ich,sngl(h1),sngl(h2),sngl(h3)
c          write(*,*)"p1",sngl(x1),sngl(y1),sngl(z1)
c          write(*,*)"p2",sngl(x2),sngl(y2),sngl(z2)
c          write(*,*)"p3",sngl(x3),sngl(y3),sngl(z3)
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
            a1=y1-y3
            a2=y3-y2
            a3=y2-y1
            b1=z1-z3
            b2=z3-z2
            b3=z2-z1
            c1=(y1-y3)*y13+(z1-z3)*z13
            c2=(y3-y2)*y32+(z3-z2)*z32
            c3=(y2-y1)*y21+(z2-z1)*z21
            if(dabs(a3*b1-a1*b3).gt.eps)then
              z0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              y0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              z0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              y0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              z0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              y0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       点:y=y0
          !------------------------------------------------------------
          if(ich.eq.2)then  !y=y0
            y0=h1
c
            a1=x1-x3
            a2=x3-x2
            a3=x2-x1
            b1=z1-z3
            b2=z3-z2
            b3=z2-z1
            c1=(x1-x3)*x13+(z1-z3)*z13
            c2=(x3-x2)*x32+(z3-z2)*z32
            c3=(x2-x1)*x21+(z2-z1)*z21
            if(dabs(a3*b1-a1*b3).gt.eps)then
              z0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              x0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              z0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              x0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              z0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              x0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       点:z=z0
          !------------------------------------------------------------
          if(ich.eq.3)then  !z=z0
            z0=h1
c
            a1=x1-x3
            a2=x3-x2
            a3=x2-x1
            b1=y1-y3
            b2=y3-y2
            b3=y2-y1
            c1=(x1-x3)*x13+(y1-y3)*y13
            c2=(x3-x2)*x32+(y3-y2)*y32
            c3=(x2-x1)*x21+(y2-y1)*y21
            if(dabs(a3*b1-a1*b3).gt.eps)then
              y0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              x0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              y0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              x0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              y0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              x0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       直線:y=ax+b
          !------------------------------------------------------------
          if(ich.eq.11)then
            aa=h1
            bb=h2
c
            a1=(x1-x3)+(y1-y3)*aa
            a2=(x3-x2)+(y3-y2)*aa
            a3=(x2-x1)+(y2-y1)*aa
            b1=z1-z3
            b2=z3-z2
            b3=z2-z1
            c1=(x1-x3)*x13+(y1-y3)*y13+(z1-z3)*z13-(y1-y3)*bb
            c2=(x3-x2)*x32+(y3-y2)*y32+(z3-z2)*z32-(y3-y2)*bb
            c3=(x2-x1)*x21+(y2-y1)*y21+(z2-z1)*z21-(y2-y1)*bb
            if(dabs(a3*b1-a1*b3).gt.eps)then
              z0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              x0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              y0=aa*x0+bb
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              z0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              x0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              y0=aa*x0+bb
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              z0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              x0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              y0=aa*x0+bb
              g1=a1*x0+b1*z0-c1
              g2=a2*x0+b2*z0-c2
              g3=a3*x0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       直線:z=ay+b
          !------------------------------------------------------------
          if(ich.eq.12)then
            aa=h1
            bb=h2
c
            a1=(x1-x3)
            a2=(x3-x2)
            a3=(x2-x1)
            b1=(y1-y3)+(z1-z3)*aa
            b2=(y3-y2)+(z3-z2)*aa
            b3=(y2-y1)+(z2-z1)*aa
            c1=(x1-x3)*x13+(y1-y3)*y13+(z1-z3)*z13-(z1-z3)*bb
            c2=(x3-x2)*x32+(y3-y2)*y32+(z3-z2)*z32-(z3-z2)*bb
            c3=(x2-x1)*x21+(y2-y1)*y21+(z2-z1)*z21-(z2-z1)*bb
            if(dabs(a3*b1-a1*b3).gt.eps)then
              y0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              x0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              z0=aa*y0+bb
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              y0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              x0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              z0=aa*y0+bb
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              y0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              x0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              z0=aa*y0+bb
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       直線:x=az+b
          !------------------------------------------------------------
          if(ich.eq.13)then
            aa=h1
            bb=h2
c
            a1=(y1-y3)
            a2=(y3-y2)
            a3=(y2-y1)
            b1=(x1-x3)*aa+(z1-z3)
            b2=(x3-x2)*aa+(z3-z2)
            b3=(x2-x1)*aa+(z2-z1)
            c1=(x1-x3)*x13+(y1-y3)*y13+(z1-z3)*z13-(x1-x3)*bb
            c2=(x3-x2)*x32+(y3-y2)*y32+(z3-z2)*z32-(x3-x2)*bb
            c3=(x2-x1)*x21+(y2-y1)*y21+(z2-z1)*z21-(x2-x1)*bb
            if(dabs(a3*b1-a1*b3).gt.eps)then
              z0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              y0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              x0=aa*z0+bb
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              z0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              y0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              x0=aa*z0+bb
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              z0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              y0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              x0=aa*z0+bb
              g1=a1*y0+b1*z0-c1
              g2=a2*y0+b2*z0-c2
              g3=a3*y0+b3*z0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
          !------------------------------------------------------------
          !       面:z=ax+by+c
          !------------------------------------------------------------
          if(ich.eq.20)then
            aa=h1   !z=ax+by+cの係数
            bb=h2
            cc=h3
            a1=(x1-x3)+(z1-z3)*aa
            a2=(x3-x2)+(z3-z2)*aa
            a3=(x2-x1)+(z2-z1)*aa
            b1=(y1-y3)+(z1-z3)*bb
            b2=(y3-y2)+(z3-z2)*bb
            b3=(y2-y1)+(z2-z1)*bb
            c1=(x1-x3)*x13+(y1-y3)*y13+(z1-z3)*z13-(z1-z3)*cc
            c2=(x3-x2)*x32+(y3-y2)*y32+(z3-z2)*z32-(z3-z2)*cc
            c3=(x2-x1)*x21+(y2-y1)*y21+(z2-z1)*z21-(z2-z1)*cc
c
            if(dabs(a3*b1-a1*b3).gt.eps)then
              y0=(a3*c1-a1*c3)/(a3*b1-a1*b3)
              x0=(b1*c3-b3*c1)/(a3*b1-a1*b3)
              z0=aa*x0+bb*y0+cc
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a2*b3-a3*b2).gt.eps)then
              y0=(a2*c3-a3*c2)/(a2*b3-a3*b2)
              x0=(b3*c2-b2*c3)/(a2*b3-a3*b2)
              z0=aa*x0+bb*y0+cc
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            if(dabs(a1*b2-a2*b1).gt.eps)then
              y0=(a1*c2-a2*c1)/(a1*b2-a2*b1)
              x0=(b2*c1-b1*c2)/(a1*b2-a2*b1)
              z0=aa*x0+bb*y0+cc
              g1=a1*x0+b1*y0-c1
              g2=a2*x0+b2*y0-c2
              g3=a3*x0+b3*y0-c3
              if(dabs(g1).lt.eps .and. dabs(g2).lt.eps
     &                           .and. dabs(g3).lt.eps)goto 100
            endif
            stop
          endif
100       continue
          dd1=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
          dd2=dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2)
          dd3=dsqrt((x3-x0)**2+(y3-y0)**2+(z3-z0)**2)
c          write(*,*)"dd",dd1,dd2,dd3
          if(dabs(dd1-dd2).lt.eps .and. dabs(dd2-dd3).lt.eps
     &                            .and. dabs(dd3-dd1).lt.eps)return
          write(*,*)"dd",dd1,dd2,dd3
          pause
c
          return
          end
c------------------------------------------------------------
c
c       連立方程式を解く
c       ax+by+cz=d
c------------------------------------------------------------
        subroutine eq31(x1,y1,z1,x2,y2,z2,x3,y3,z3,eps,ich,h1,h2,h3)
        implicit double precision(a-h,o-z)
        parameter(n=3,m=4)
        dimension a(n,m),iwork(n)
        h1=-99d0
        h2=-99d0
        h3=-99d0
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
        if(dabs(x1-x3).gt.eps)then !4-1)
          a0=(y1-y3)/(x1-x3)
          b0=(x3*y1-x1*y3)/(x3-x1)
          g1=a0*x1+b0-y1
          g2=a0*x2+b0-y2
          g3=a0*x3+b0-y3
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
        if(dabs(y1-y3).gt.eps)then
          a0=(z1-z3)/(y1-y3)
          b0=(y3*z1-y1*z3)/(y3-y1)
          g1=a0*y1+b0-z1
          g2=a0*y2+b0-z2
          g3=a0*y3+b0-z3
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
        if(dabs(z1-z3).gt.eps)then
          a0=(x1-x3)/(z1-z3)
          b0=(z3*x1-z1*x3)/(z3-z1)
          g1=a0*z1+b0-x1
          g2=a0*z2+b0-x2
          g3=a0*z3+b0-x3
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
        a(1,1)=x1
        a(1,2)=y1
        a(1,3)=1d0
        a(2,1)=x2
        a(2,2)=y2
        a(2,3)=1d0
        a(3,1)=x3
        a(3,2)=y3
        a(3,3)=1d0
c
        a(1,4)=z1
        a(2,4)=z2
        a(3,4)=z3
c
        call sweep(a,n,m,eps,iwork,ill)  !n×n行列, m=n+1 (右辺の値)
        if(ill.eq.1)write(*,*)"ill"
        if(ill.eq.1)write(*,*)"p1:",sngl(x1),sngl(y1),sngl(z1)
        if(ill.eq.1)write(*,*)"p2:",sngl(x2),sngl(y2),sngl(z2)
        if(ill.eq.1)write(*,*)"p2:",sngl(x3),sngl(y3),sngl(z3)
        if(ill.eq.1)stop
c        do 200 i=1,n
c200       write(6,10)i,(a(i,j),j=n+1,m)
c10      format(1H ,10X,'X',I1,'=',4F8.4)
        g1=x1*a(1,m)+y1*a(2,m)+1d0*a(3,m) - z1
        g2=x2*a(1,m)+y2*a(2,m)+1d0*a(3,m) - z2
        g3=x3*a(1,m)+y3*a(2,m)+1d0*a(3,m) - z3
        g=dabs(g1)+dabs(g2)+dabs(g3)
        if(dabs(g).lt.eps)then
          ich=20
          h1=a(1,m)
          h2=a(2,m)
          h3=a(3,m)
          return
        endif
c
        write(*,*)"誤差 : ",g1,g2,g3
        pause
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
        if(dn.lt.1d-10)write(*,*)"外積:",dn
        if(dn.lt.1d-10)stop
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