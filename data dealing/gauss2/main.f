        implicit double precision(a-h,o-z)
c
        eps=1d-10
c
        x1= 0d0
        y1= 1.1354d-8
        z1=-1.134987d0
        x2= 0d0
        y2= 1.8454d-8
        z2= 1d0
        x3= 1d0
        y3= 1.7861d-9
        z3= 0d0
        x4= 1.1268d-9
        y4= 1.7431d0
        z4= 0d0
        call certet2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0,eps)
        pause
c
        i=0
        do 200 i1=-2,2
          x1=dble(i1)
        do 200 i2=-2,2
          y1=dble(i2)
        do 200 i3=-2,2
          z1=dble(i3)
        do 200 j1=-2,2
          x2=dble(j1)
        do 200 j2=-2,2
          y2=dble(j2)
        do 200 j3=-2,2
          z2=dble(j3)
        do 200 k1=-2,2
          x3=dble(k1)
        do 200 k2=-2,2
          y3=dble(k2)
        do 200 k3=-2,2
          z3=dble(k3)
        do 200 l1=-2,2
          x4=dble(l1)
        do 200 l2=-2,2
          y4=dble(l2)
        do 200 l3=-2,2
          z4=dble(l3)
c
c         面積の計算
c          ar1=artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
c          ar2=artri(x2,y2,z2,x3,y3,z2,x4,y4,z4)
c          ar3=artri(x3,y3,z3,x4,y4,z4,x1,y1,z1)
c          ar4=artri(x4,y4,z4,x1,y1,z1,x2,y2,z2)
c          if(dabs(ar1).lt.eps .or. dabs(ar2).lt.eps .or.
c     &       dabs(ar3).lt.eps .or. dabs(ar4).lt.eps)goto 200
c
c         体積の計算
          vl=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          if(dabs(vl).lt.eps)goto 200
c          write(*,*)"vl:",sngl(vl)
          i=i+1
          write(*,*)i,i1,i2,i3
c
c         外心の計算
c          call certet(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0,eps)
          call certet2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x0,y0,z0,eps)
200     continue
c
        stop
        end