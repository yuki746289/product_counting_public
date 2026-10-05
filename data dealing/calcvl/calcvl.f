        implicit double precision(a-h,o-z)
        parameter(ind=99999,ine=99999)
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
c
        im=1
        ne(1,1)=2
        ne(1,2)=3
        ne(1,3)=4
        ne(1,4)=1
c
c       ß“_1
        xx(1)=1d0
        yy(1)=1d0
        zz(1)=1d0
c       ß“_2
        xx(2)=3d0
        yy(2)=1d0
        zz(2)=1d0
c       ß“_3
        xx(3)=1d0
        yy(3)=1d0
        zz(3)=2d0
c       ß“_4
        xx(4)=1d0
        yy(4)=5d0
        zz(4)=1d0
c
        a11=1d0
        a21=1d0
        a31=1d0
        a41=1d0
        a12=xx(ne(im,1))
        a22=xx(ne(im,2))
        a32=xx(ne(im,3))
        a42=xx(ne(im,4))
        a13=yy(ne(im,1))
        a23=yy(ne(im,2))
        a33=yy(ne(im,3))
        a43=yy(ne(im,4))
        a14=zz(ne(im,1))
        a24=zz(ne(im,2))
        a34=zz(ne(im,3))
        a44=zz(ne(im,4))
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
        vl=dabs(det)/6d0
c        if(vl.le.0d0)then
          write(*,*)'vl',sngl(vl),det/6d0
          write(*,*)'x',sngl(xx(1)),sngl(xx(2)),sngl(xx(3)),sngl(xx(4))
          write(*,*)'y',sngl(yy(1)),sngl(yy(2)),sngl(yy(3)),sngl(yy(4))
          write(*,*)'z',sngl(zz(1)),sngl(zz(2)),sngl(zz(3)),sngl(zz(4))
c          pause
c        endif
        write(*,*)'vl-th',2d0*1d0/2d0*4d0/3d0   !—˜_’l
c
        stop
        end