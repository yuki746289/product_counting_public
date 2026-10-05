       implicit double precison(a-h,o-z)
       parameter(ind=9999,ibw=999)
       dimension xx(ind),yy(ind),zz(ind) 
       xr=4d0
       yr=2d0
       zr=3d0
       n=6+1       !節点の数
c      x座標
       xx(1)=4d0
       xx(2)=4d0
       xx(3)=4d0
       xx(4)=4d0
       xx(5)=3.9d0
       xx(6)=4.1d0
       xx(7)=xr
c      y座標
       yy(1)=2d0
       yy(2)=2d0
       yy(3)=1.9d0
       yy(4)=2.1d0
       yy(5)=2d0
       yy(6)=2d0
       yy(7)=yr
c      z座標
       zz(1)=2.9d0
       zz(2)=3.1d0
       zz(3)=3d0
       zz(4)=3d0
       zz(5)=3d0
       zz(6)=3d0
       zz(7)=zr
       call boundfs(n,xx,yy,zz,xr,yr,zr
     &                          ,x0,y0,z0,dnx,dny,dnz,ind,ibw)
       rr=dsqrt((xr-x0)**2+(yr-y0)**2+(zr-z0)**2)    !曲率半径
       write(*,*)"the counting values"
       write(*,*)"x0,y0,z0,rr",sngl(x0),sngl(y0),sngl(z0),sngl(rr)
       write(*,*)"dnx,dny,dnz",sngl(dnx),sngl(dny),sngl(dnz)
c      the values of theory
       x0t=2d0
       y0t=1d0
       z0t=1d0
       rr=dsqrt((xr-x0t)**2+(yr-y0t)**2+(zr-z0t)**2)    !曲率半径
       zk=1d0/rr          !曲率
       dnx=(xr-x0t)/rr     !法線ベクトルの方向余弦のx成分
       dny=(yr-y0t)/rr     !法線ベクトルの方向余弦のy成分
       dnz=(zr-z0t)/rr     !法線ベクトルの方向余弦のz成分
       write(*,*)"the theory values"
       write(*,*)"x0t,y0t,z0t,rr",sngl(x0t),sngl(y0t),sngl(z0t),sngl(rr)
       write(*,*)"dnx,dny,dnz",sngl(dnx),sngl(dny),sngl(dnz)
c
       stop
       end
c
c       サブルーチン
c
       subroutine boundfs(n,xx,yy,zz,xr,yr,zr
     &                          ,x0,y0,z0,dnx,dny,dnz,ind,ibw)
       implicit double precision(a-h,o-z)
       dimension sax(3,ibw),say(3,ibw),saz(3,ibw)
       dimension sfx(3),sfy(3),sfz(3)
       dimension xx(ind),yy(ind),zz(ind)
c
       do 1000 k=0,4
         do 1100 i=1,n
           x(k+1)=2d0*(xx(i)**k+xx(i)**(k+1)+xx(i)**(k+2))
           y(k+1)=2d0*(yy(i)**k+yy(i)**(k+1)+yy(i)**(k+2))
           z(k+1)=2d0*(zz(i)**k+zz(i)**(k+1)+zz(i)**(k+2))
1100     continue
1000   continue
c
       nbw=3     !バンド幅
       do 1200 i=1,3
         sfx(i)=-x(i+1)
         sfy(i)=-y(i+1)
         sfz(i)=-z(i+1)
         do 1300 j=1,3
           sax(i,j-i+nbw)=x((i-1)+(j-1)+1)
           say(i,j-i+nbw)=y((i-1)+(j-1)+1)
           saz(i,j-i+nbw)=z((i-1)+(j-1)+1)
1300     continue
1200   continue
c
       call gauss(3,nbw,sax,sfx,ind,ibw)
       call gauss(3,nbw,say,sfy,ind,ibw)
       call gauss(3,nbw,saz,sfz,ind,ibw)
c
       x0=-sfx(1)/2d0     !a1=sfx(1)
       y0=-sfy(1)/2d0     !b1=sfy(1)
       z0=-sfz(1)/2d0     !c1=sfz(1)
       rr=dsqrt((xr-x0)**2+(yr-y0)**2+(zr-z0)**2)    !曲率半径
       zk=1d0/rr          !曲率
       dnx=(xr-x0)/rr     !法線ベクトルの方向余弦のx成分
       dny=(yr-y0)/rr     !法線ベクトルの方向余弦のy成分
       dnz=(zr-z0)/rr     !法線ベクトルの方向余弦のz成分
       return
       end
c
       subroutine gauss(n,nbw,a,x,ind,ibw)
       implicit double precision (a-h,o-z)
       dimension a(ind,ibw),x(ind)
c
       do 2 i=1,n
       if(a(i,nbw).le.0) then
        write(*,*) '[gauss]',i,a(i,nbw)
        stop
       endif
2      continue
c
       nbwm=nbw-1
       do 10 i=1,n
         aa=a(i,nbw)
         x(i)=x(i)/aa
         do 20 j=i,min(n,i+nbwm)
20       a(i,nbw+j-i)=a(i,nbw+j-i)/aa
         do 40 k=i+1,min(n,i+nbwm)
           cc=a(k,nbw+i-k)
           do 30 j=i+1,min(n,i+nbwm)
30         a(k,nbw+j-k)=a(k,nbw+j-k)-cc*a(i,nbw+j-i)
40       x(k)=x(k)-cc*x(i)
10     continue
       do 70 i=n-1,1,-1
       do 70 k=i+1,min(n,i+nbwm)
70     x(i)=x(i)-a(i,nbw+k-i)*x(k)
       return
       end