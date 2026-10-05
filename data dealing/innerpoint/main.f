c       3角形の内心を計算する
        implicit double precision(a-h,o-z)
        dimension ms(999,3)
        dimension xs(999),ys(999),zs(999)
c
        n1=1
        n2=2
        n3=3
        ms(1,1)=n1
        ms(1,2)=n2
        ms(1,3)=n3
c
        xs(n1)=0d0
        ys(n1)=1d0
        zs(n1)=0d0
        xs(n2)=1d0
        ys(n2)=0d0
        zs(n2)=0d0
        xs(n3)=-1d0
        ys(n3)=0d0
        zs(n3)=0d0
c
        x1=xs(n1)
        x2=xs(n2)
        x3=xs(n3)
        y1=ys(n1)
        y2=ys(n2)
        y3=ys(n3)
        z1=zs(n1)
        z2=zs(n2)
        z3=zs(n3)
c
        call innerpoint(x1,y1,z1,x2,y2,z2,x3,y3,z3,xin,yin,zin)
c
        stop
        end
c----------------------------------------------------------------------
c
c       3角形の内心を算出する
c
c----------------------------------------------------------------------
        subroutine innerpoint(x1,y1,z1,x2,y2,z2,x3,y3,z3,xin,yin,zin)
        implicit double precision(a-h,o-z)
c
        p12=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)  !|p12|=|p21|
        p21=dsqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        p23=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)  !|p23|=|p32|
        p32=dsqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        p31=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)  !|p31|=|p13|
        p13=dsqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
c
        a1=(x1-x2)/p12+(x3-x2)/p32
        a2=(y1-y2)/p12+(y3-y2)/p32
        a3=(z1-z2)/p12+(z3-z2)/p32
        b1=-(x2-x3)/p23-(x1-x3)/p13
        b2=-(y2-y3)/p23-(y1-y3)/p13
        b3=-(z2-z3)/p23-(z1-z3)/p13
c
        if(dabs(x1).lt.1d-8 .and. dabs(x2).lt.1d-8
     &                      .and. dabs(x3).lt.1d-8)then    !x平面上
          xi2=y2
          xi3=y3
          ai=a2
          bi=b2
        else
          xi2=x2
          xi3=x3
          ai=a1
          bi=b1
        endif
c
        d2=-bi*(x3+y3+z3-x2-y2-z2)+(b1+b2+b3)*(-xi2+xi3)
        d3= ai*(x3+y3+z3-x2-y2-z2)-(a1+a2+a3)*(-xi2+xi3)
        d=ai*(b1+b2+b3)-(a1+a2+a3)*bi
        dl2=d2/d
        dl3=d3/d
c
        xin=x2+dl2*((x1-x2)/p12+(x3-x2)/p32)
        yin=y2+dl2*((y1-y2)/p12+(y3-y2)/p32)
        zin=z2+dl2*((z1-z2)/p12+(z3-z2)/p32)
c
        xin0=x3+dl3*((x2-x3)/p23+(x1-x3)/p13)
        yin0=y3+dl3*((y2-y3)/p23+(y1-y3)/p13)
        zin0=z3+dl3*((z2-z3)/p23+(z1-z3)/p13)
c
        write(*,*)"xin",sngl(xin),sngl(yin),sngl(zin)
        write(*,*)"xin",sngl(xin0),sngl(yin0),sngl(zin0)
c       検算
        ar=artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        dr=2d0*ar/(p12+p23+p31)
        write(*,*)"dr",sngl(dr)
c
        stop
        end
c----------------------------------------------------------
c
c       3角形の面積
c
c----------------------------------------------------------
        function artri(x1,y1,z1,x2,y2,z2,x3,y3,z3)
        implicit double precision(a-h,o-z)
c
        vx=(y1-y3)*(z2-z3)-(z1-z3)*(y2-y3)
        vy=(z1-z3)*(x2-x3)-(x1-x3)*(z2-z3)
        vz=(x1-x3)*(y2-y3)-(y1-y3)*(x2-x3)
        artri=dsqrt(vx**2+vy**2+vz**2)/2d0
        return
        end
