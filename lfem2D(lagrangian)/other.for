c-------------------------------------------------------------------------------
c
c      曲率半径と中心を計算する
c
c-------------------------------------------------------------------------------
       subroutine radii(xx,yy,n1,n2,n3,rd,pxx,pyy)
       include "head.for"
       dimension xx(ind),yy(ind)
c
       x1= xx(n1)
       x2= xx(n2)
       x3= xx(n3)
c
       y1= yy(n1)
       y2= yy(n2)
       y3= yy(n3)
c
       call circum(x1,x2,x3,y1,y2,y3,rd,x,y)
c
       if(rd.eq.0d0)then
         rd=1d0
         pxx=0d0
         pyy=0d0
       else
         pxx=x2-x
         pyy=y2-y
c
         pxx=pxx/dsqrt(pxx**2)
         pyy=pyy/dsqrt(pyy**2)
       endif
c
       return
       end
c-------------------------------------------------------------------------------
c
c       3角形の面積
c
c-------------------------------------------------------------------------------
        function artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        vx=(y1-y3)*(z2-z3)-(z1-z3)*(y2-y3)
        vy=(z1-z3)*(x2-x3)-(x1-x3)*(z2-z3)
        vz=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)
        artri=dsqrt(vx**2+vy**2+vz**2)/2d0
        return
        end
c-------------------------------------------------------------------------------
c
c          外接円の中心と半径を計算する
c
c          引数
c          x1,x2,x3,y1,y2,y3
c
c          戻り値
c          rr :外接円の半径
c          x,y:外接円の中心半径
c
c-------------------------------------------------------------------------------
       subroutine circum(x1,x2,x3,y1,y2,y3,rr,x,y)
       implicit double precision(a-h,o-z)
c
c      辺の長さ
       a=dsqrt((x2-x3)**2+(y2-y3)**2)
       b=dsqrt((x3-x1)**2+(y3-y1)**2)
       c=dsqrt((x1-x2)**2+(y1-y2)**2)
c
c      面積
       ar=artri(x1,y1,0d0,x2,y2,0d0,x3,y3,0d0)
       if(ar.lt.1d-7)then
         rr=0d0
         return
       endif
c
c      外接円の中心
       x=a**2*(b**2+c**2-a**2)*x1
     &  +b**2*(c**2+a**2-b**2)*x2
     &  +c**2*(a**2+b**2-c**2)*x3
       x=x/16d0/ar**2
       y=a**2*(b**2+c**2-a**2)*y1
     &  +b**2*(c**2+a**2-b**2)*y2
     &  +c**2*(a**2+b**2-c**2)*y3
       y=y/16d0/ar**2
c
c      外接円の半径
       r1=dsqrt((x1-x)**2+(y1-y)**2)
       r2=dsqrt((x2-x)**2+(y2-y)**2)
       r3=dsqrt((x3-x)**2+(y3-y)**2)
       rr=(r1+r2+r3)/3d0
c
c      検算
       if(dabs(r1-rr).gt.1d-5 .or.
     &    dabs(r2-rr).gt.1d-5 .or. 
     &    dabs(r3-rr).gt.1d-5)then
          write(*,*)"circum"
          write(*,*)"radian is different"
          write(*,*)"r:",sngl(r1),sngl(r2),sngl(r3),sngl(rr)
          stop
        endif
c
       return
       end
