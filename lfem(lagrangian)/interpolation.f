c------------------------------------------------------------------------------------
c
c       変数をコピーする
c
c------------------------------------------------------------------------------------
        subroutine copydata(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt,idx
     &                ,npn,nen,neln,xn,yn,zn,un,vn,wn,pn,cn,tn,idn)
        include "header.h"
        dimension ne(ine,4), nen(ine,4),idx(0:ine),idn(0:ine)
        dimension xx(ind),yy(ind),zz(ind),uu(ind),vv(ind),ww(ind)
        dimension xn(ind),yn(ind),zn(ind),un(ind),vn(ind),wn(ind)
        dimension cc(ind),cn(ind),tt(ind),tn(ind),pp(ind),pn(ind)
c
        npn=np
        neln=nele
        do 10 i=1,nele
          idn(i)=idx(i)
          do 10 j=1,4
          nen(i,j)=ne(i,j)
10      continue
c
        do 20 kp=1,np
          xn(kp)=xx(kp)
          yn(kp)=yy(kp)
          zn(kp)=zz(kp)
          un(kp)=uu(kp)
          vn(kp)=vv(kp)
          wn(kp)=ww(kp)
          pn(kp)=pp(kp)
c          tn(kp)=tt(kp)
c          cn(kp)=cc(kp)
20      continue
c
        return
        end
c------------------------------------------------------------------------------------
c
c       速度、圧力を補完する
c
c------------------------------------------------------------------------------------
        subroutine interpolation(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt
     &                ,idx,npn,nen,neln,xn,yn,zn,un,vn,wn,pn,cn,tn,idn)
        include "header.h"
        dimension ne(ine,4),nen(ine,4),idx(0:ine),idn(0:ine)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xn(ind),yn(ind),zn(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind),cc(ind),tt(ind)
        dimension un(ind),vn(ind),wn(ind),pn(ind),cn(ind),tn(ind)
c
        do 10 kp=1,np
          dnmin=99d0
          do 20 keln=1,neln
            k1=nen(keln,1)
            k2=nen(keln,2)
            k3=nen(keln,3)
            k4=nen(keln,4)
c           形状関数を計算
            call caldn(xn(k1),yn(k1),zn(k1),xn(k2),yn(k2),zn(k2)
     &                ,xn(k3),yn(k3),zn(k3),xn(k4),yn(k4),zn(k4)
     &                ,xx(kp),yy(kp),zz(kp),dn)
c            write(*,*)"keln:",keln,sngl(dn)
            if(dn.lt.dnmin)kmin=keln
            if(dn.lt.dnmin)dnmin=dn
20        continue
          if(dnmin.lt.0.9d0 .or. dnmin.gt.1.2d0)then
            write(*,*)"in interpolation dn:",dnmin
            stop
          endif
          k1=nen(kmin,1)
          k2=nen(kmin,2)
          k3=nen(kmin,3)
          k4=nen(kmin,4)
          xp=xx(kp)
          yp=yy(kp)
          zp=zz(kp)
          x1=xn(k1)
          x2=xn(k2)
          x3=xn(k3)
          x4=xn(k4)
          y1=yn(k1)
          y2=yn(k2)
          y3=yn(k3)
          y4=yn(k4)
          z1=zn(k1)
          z2=zn(k2)
          z3=zn(k3)
          z4=zn(k4)
          v1=calvl(xp,yp,zp,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          v2=calvl(x1,y1,z1,xp,yp,zp,x3,y3,z3,x4,y4,z4)
          v3=calvl(x1,y1,z1,x2,y2,z2,xp,yp,zp,x4,y4,z4)
          v4=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp)
          vt=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          dn1=v1/vt
          dn2=v2/vt
          dn3=v3/vt
          dn4=v4/vt
          uu(kp)=un(k1)*dn1+un(k2)*dn2+un(k3)*dn3+un(k4)*dn4
          vv(kp)=vn(k1)*dn1+vn(k2)*dn2+vn(k3)*dn3+vn(k4)*dn4
          ww(kp)=wn(k1)*dn1+wn(k2)*dn2+wn(k3)*dn3+wn(k4)*dn4
          pp(kp)=pn(k1)*dn1+pn(k2)*dn2+pn(k3)*dn3+pn(k4)*dn4
10      continue
c
        do 30 kele=1,nele
30      idx(kele)=idn(1)
c
        return
        end
c------------------------------------------------------------------------------------
c
c       節点kpの速度、圧力を補完する
c
c------------------------------------------------------------------------------------
        subroutine interpolation2(kp,xp,yp,zp
     &                   ,np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
        include "header.h"
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind),cc(ind),tt(ind)
        data ceft/2d-1/   !緩和係数
c
          dnmin=99d0
          do 20 kele=1,nele
            k1=ne(kele,1)
            k2=ne(kele,2)
            k3=ne(kele,3)
            k4=ne(kele,4)
c           形状関数を計算
            call caldn(xx(k1),yy(k1),zz(k1),xx(k2),yy(k2),zz(k2)
     &                ,xx(k3),yy(k3),zz(k3),xx(k4),yy(k4),zz(k4)
     &                ,xp,yp,zp,dn)
c            write(*,*)"kele:",kele,sngl(dn)
            if(dn.lt.dnmin)kmin=kele
            if(dn.lt.dnmin)dnmin=dn
20        continue
c          if(dnmin.lt.0.9d0 .or. dnmin.gt.1.2d0)then
          if(dnmin.lt.0.9d0 .or. dnmin.gt.5d0)then
            write(*,*)"in interpolation dn:",dnmin
            stop
          endif
          k1=ne(kmin,1)
          k2=ne(kmin,2)
          k3=ne(kmin,3)
          k4=ne(kmin,4)
          x1=xx(k1)
          x2=xx(k2)
          x3=xx(k3)
          x4=xx(k4)
          y1=yy(k1)
          y2=yy(k2)
          y3=yy(k3)
          y4=yy(k4)
          z1=zz(k1)
          z2=zz(k2)
          z3=zz(k3)
          z4=zz(k4)
          v1=calvl(xp,yp,zp,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          v2=calvl(x1,y1,z1,xp,yp,zp,x3,y3,z3,x4,y4,z4)
          v3=calvl(x1,y1,z1,x2,y2,z2,xp,yp,zp,x4,y4,z4)
          v4=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp)
          vt=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          dn1=v1/vt
          dn2=v2/vt
          dn3=v3/vt
          dn4=v4/vt
          u0=uu(k1)*dn1+uu(k2)*dn2+uu(k3)*dn3+uu(k4)*dn4
          v0=vv(k1)*dn1+vv(k2)*dn2+vv(k3)*dn3+vv(k4)*dn4
          w0=ww(k1)*dn1+ww(k2)*dn2+ww(k3)*dn3+ww(k4)*dn4
          p0=pp(k1)*dn1+pp(k2)*dn2+pp(k3)*dn3+pp(k4)*dn4
c
          uu(kp)=(1d0-ceft)*uu(kp)+ceft*u0
          vv(kp)=(1d0-ceft)*vv(kp)+ceft*v0
          ww(kp)=(1d0-ceft)*ww(kp)+ceft*w0
          pp(kp)=(1d0-ceft)*pp(kp)+ceft*p0
c
          xx(kp)=xp
          yy(kp)=yp
          zz(kp)=zp
c
        return
        end
