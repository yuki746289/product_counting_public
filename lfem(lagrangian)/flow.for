       subroutine setmat(sa,sf,np,nbw)    !nbw=4*nbw
       include "header.h"
       dimension sa(4*ind,ibw),sf(4*ind)
       do 1000 i=1,np
         sf(4*i-3)=0d0
         sf(4*i-2)=0d0
         sf(4*i-1)=0d0
         sf(4*i-0)=0d0
         do 1100 j=1,2*nbw
           sa(4*i-3,j)=0d0
           sa(4*i-2,j)=0d0
           sa(4*i-1,j)=0d0
           sa(4*i-0,j)=0d0
1100     continue
1000   continue
       return
       end
c------------------------------------------------------------
c      直接法でマトリックスを解く場合
c------------------------------------------------------------
       subroutine flow(nel,nbw,sa,sf,ne,xx,yy,zz    !ne(ine,4):流動領域の要素のみ格納
     &                     ,uvwp0,tt,idx,ind,ine,ibw)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)
       dimension sa(4*ind,ibw),sf(4*ind),uvwp0(4*ind)
       dimension x(4),y(4),z(4),b(4),c(4),d(4),zkesu(4,4)
       dimension tt(ind)    !物性値を求める時に用いる
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
       common /phyc/dens(1001),visc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
       do 200 im=1,nel     !要素imについて
         id=idx(im)
         usum=0d0
         vsum=0d0
         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           x(i)=xx(ne(im,i))
           y(i)=yy(ne(im,i))
           z(i)=zz(ne(im,i))
           usum=usum+uvwp0(4*ne(im,i)-3)
           vsum=vsum+uvwp0(4*ne(im,i)-2)
           wsum=wsum+uvwp0(4*ne(im,i)-1)
210      continue
c
c        係数行列の成分を算出する
c
         b(1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
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
         if(vl.le.0d0)then
           write(*,*)'vl',sngl(vl)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
c        無次元数および物性値を求める
c
         rey=re(id)     !1/Re
         svc=sv(id)**2  !(1/Ma)^2
         ggy=gv(id)     !1/Fr
         vis=(calvis(tt(ne(im,1)))+calvis(tt(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
     &       +calvis(tt(ne(im,3)))+calvis(tt(ne(im,4))))/4d0/visc(id)
         den=(calden(tt(ne(im,1)))+calden(tt(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
     &       +calden(tt(ne(im,3)))+calden(tt(ne(im,4))))/4d0/dens(id)
         vis=1d0    !無次元粘度(温度変化分)
         den=1d0    !無次元密度(温度変化分)
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c         write(*,*)'vis,den',sngl(vis),sngl(den)

         do 280 i=1,4   !重力考慮:y軸方向
280      sf(ne(im,i)*4-2)=sf(ne(im,i)*4-2)+ggy*vl/4d0
c
         do 290 i=1,4     !要素imの行
           iu=ne(im,i)*4-3
           iv=ne(im,i)*4-2
           iw=ne(im,i)*4-1
           ip=ne(im,i)*4-0
           do 290 j=1,4   !要素imの列
             ju=ne(im,j)*4-3
             jv=ne(im,j)*4-2
             jw=ne(im,j)*4-1
             jp=ne(im,j)*4-0
c             
             cmm=vl/20d0*zkesu(i,j)     !係数行列:[C]
             sxx=b(i)*b(j)/36d0/vl      !係数行列:[Sxx]
             syy=c(i)*c(j)/36d0/vl      !係数行列:[Syy]
             szz=d(i)*d(j)/36d0/vl      !係数行列:[Szz]
             sxy=b(i)*c(j)/36d0/vl      !係数行列:[Sxy]
             syx=c(i)*b(j)/36d0/vl      !係数行列:[Syx]
             syz=c(i)*d(j)/36d0/vl      !係数行列:[Syz]
             szy=d(i)*c(j)/36d0/vl      !係数行列:[Szy]
             szx=d(i)*b(j)/36d0/vl      !係数行列:[Szx]
             sxz=b(i)*d(j)/36d0/vl      !係数行列:[Sxz]
             chxij=b(j)/24d0            !係数行列:[Hx]
             chxji=b(i)/24d0            !係数行列:[Hx]T
             chyij=c(j)/24d0            !係数行列:[Hy]
             chyji=c(i)/24d0            !係数行列:[Hy]T
             chzij=d(j)/24d0            !係数行列:[Hz]
             chzji=d(i)/24d0            !係数行列:[Hz]T

c             cxx=(uvwp0(4*ne(im,i)-3)+usum)*b(j)/120d0  !係数行列:[Cxx]←対流項
c             cyy=(uvwp0(4*ne(im,i)-2)+vsum)*c(j)/120d0  !係数行列:[Cyy]←対流項
c             czz=(uvwp0(4*ne(im,i)-1)+wsum)*d(j)/120d0  !係数行列:[Czz]←対流項
             cxx=0d0
             cyy=0d0
             czz=0d0
c             
             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
             a44=cmm                    !連続の式:[P]の時間の偏微分項
c
             b11=(2d0*sxx+syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式X成分:[U]の距離の偏微分項
             b12=syx*rey*vis/den        !運動量収支式X成分:[V]の距離の偏微分項
             b13=szx*rey*vis/den        !運動量収支式X成分:[W]の距離の偏微分項
             b14=-chxji/den             !運動量収支式X成分:[P]の距離の偏微分項
c             
             b21=sxy*rey*vis/den        !運動量収支式Y成分:[U]の距離の偏微分項
             b22=(sxx+2d0*syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Y成分:[V]の距離の偏微分項
             b23=szy*rey*vis/den        !運動量収支式Y成分:[W]の距離の偏微分項
             b24=-chyji/den             !運動量収支式Y成分:[P]の距離の偏微分項
c             
             b31=sxz*rey*vis/den        !運動量収支式Z成分:[U]の距離の偏微分項
             b32=syz*rey*vis/den        !運動量収支式Z成分:[V]の距離の偏微分項
             b33=(sxx+syy+2d0*szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Z成分:[W]の距離の偏微分項
             b34=-chzji/den             !運動量収支式Z成分:[P]の距離の偏微分項
c             
             b41=svc*chxij               !連続の式:[U]の距離の偏微分項
             b42=svc*chyij               !連続の式:[V]の距離の偏微分項
             b43=svc*chzij               !連続の式:[W]の距離の偏微分項
             b44=0d0                    !連続の式:[P]の距離の偏微分項
c             
c            運動量収支式X成分
             sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11   ![U]の項のj成分
             sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12   ![V]の項のj成分
             sa(iu,nbw+jw-iu)=sa(iu,nbw+jw-iu)       +b13   ![W]の項のj成分
             sa(iu,nbw+jp-iu)=sa(iu,nbw+jp-iu)       +b14
             sf(iu)=sf(iu)+a11/dt*uvwp0(ju)
c             
c            運動量収支式Y成分
             sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21   ![U]の項のj成分
             sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22   ![V]の項のj成分
             sa(iv,nbw+jw-iv)=sa(iv,nbw+jw-iv)       +b23   ![W]の項のj成分
             sa(iv,nbw+jp-iv)=sa(iv,nbw+jp-iv)       +b24
             sf(iv)=sf(iv)+a22/dt*uvwp0(jv)
c
c            運動量収支式Z成分
             sa(iw,nbw+ju-iw)=sa(iw,nbw+ju-iw)       +b31   ![U]の項のj成分
             sa(iw,nbw+jv-iw)=sa(iw,nbw+jv-iw)       +b32   ![V]の項のj成分
             sa(iw,nbw+jw-iw)=sa(iw,nbw+jw-iw)+a33/dt+b33   ![W]の項のj成分
             sa(iw,nbw+jp-iw)=sa(iw,nbw+jp-iw)       +b34
             sf(iw)=sf(iw)+a33/dt*uvwp0(jw)
c
c            連続の式
             sa(ip,nbw+ju-ip)=sa(ip,nbw+ju-ip)       +b41
             sa(ip,nbw+jv-ip)=sa(ip,nbw+jv-ip)       +b42
             sa(ip,nbw+jw-ip)=sa(ip,nbw+jw-ip)       +b43
             sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a44/dt+b44
             sf(ip)=sf(ip)+a44/dt*uvwp0(jp)
290      continue
200    continue
c
       return
       end
c------------------------------------------------------------
c
c      SIMPLE法でマトリクスを解く場合
c
c------------------------------------------------------------
c
c      p0(np)を用いてN-S式からuu(np),vv(np),ww(np)のマトリクスsa(3*ind,ibw),sf(3*ind)を組み立てる
c
       subroutine flowns(nel,nbw,sa,sf,ne,xx,yy,zz    !ne(ine,4):新_新
     &                  ,uvwp,uvwp0,t0,idx)             !nbw=3*nbw
       include "header.h"
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)  !新_新
       dimension uvwp(4*ind),uvwp0(4*ind)
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension x(4),y(4),z(4),b(4),c(4),d(4),zkesu(4,4)
       dimension t0(ind)    !物性値を求める時に用いる
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
c       common /phyc/dens(1001),visc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
c       open(1,file='sa.res')
c       rewind(1)
c
       do 200 im=1,nel     !nel:新_新
         id=idx(im)
         usum=0d0
         vsum=0d0
         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           x(i)=xx(ne(im,i))
           y(i)=yy(ne(im,i))
           z(i)=zz(ne(im,i))
c-----------------------------------------------------------------
c
           usum=usum+uvwp(4*ne(im,i)-3)    !uvwp(4*新_新-i)
           vsum=vsum+uvwp(4*ne(im,i)-2)
           wsum=wsum+uvwp(4*ne(im,i)-1)
c
c-----------------------------------------------------------------
210      continue
c
c        係数行列の成分を算出する
c
         b(1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
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
         if(vl.le.0d0)then
           write(*,*)'im',im
           write(*,*)'vl',sngl(vl)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
c        無次元数および物性値を求める
c
         rey=re(id)     !1/Re
c         svc=sv(id)**2  !(1/Ma)^2
         ggy=gv(id)     !1/Fr
         vis=(calvis(t0(ne(im,1)))+calvis(t0(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
     &       +calvis(t0(ne(im,3)))+calvis(t0(ne(im,4))))
     &       /4d0/calvis(-999d0)
         den=(calden(t0(ne(im,1)))+calden(t0(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
     &       +calden(t0(ne(im,3)))+calden(t0(ne(im,4))))
     &       /4d0/calden(-999d0)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,id',im,id
c         write(*,*)'dt,vl',sngl(dt),sngl(vl)
c         write(*,*)'usum,vsum,wsum',sngl(usum),sngl(vsum),sngl(wsum)
c         write(*,*)'b(i)',sngl(b(1)),sngl(b(2)),sngl(b(3)),sngl(b(4))
c         write(*,*)'c(i)',sngl(c(1)),sngl(c(2)),sngl(c(3)),sngl(c(4))
c         write(*,*)'d(i)',sngl(d(1)),sngl(d(2)),sngl(d(3)),sngl(d(4))
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c         write(*,*)'vis,den',sngl(vis),sngl(den)
c         write(*,*)'---------------------------------------------------'
c         pause
c
         do 280 i=1,4   !重力考慮:-y軸方向
280      sf(ne(im,i)*3-1)=sf(ne(im,i)*3-1)+ggy*vl/4d0
c
         do 290 i=1,4     !要素imの行
           iu=ne(im,i)*3-2
           iv=ne(im,i)*3-1
           iw=ne(im,i)*3-0
           do 290 j=1,4   !要素imの列
             ju=ne(im,j)*3-2
             jv=ne(im,j)*3-1
             jw=ne(im,j)*3-0
c             
             cmm=vl/20d0*zkesu(i,j)     !係数行列:[C]
             sxx=b(i)*b(j)/36d0/vl      !係数行列:[Sxx]
             syy=c(i)*c(j)/36d0/vl      !係数行列:[Syy]
             szz=d(i)*d(j)/36d0/vl      !係数行列:[Szz]
             sxy=b(i)*c(j)/36d0/vl      !係数行列:[Sxy]
             syx=c(i)*b(j)/36d0/vl      !係数行列:[Syx]
             syz=c(i)*d(j)/36d0/vl      !係数行列:[Syz]
             szy=d(i)*c(j)/36d0/vl      !係数行列:[Szy]
             szx=d(i)*b(j)/36d0/vl      !係数行列:[Szx]
             sxz=b(i)*d(j)/36d0/vl      !係数行列:[Sxz]
c             chxij=b(j)/24d0            !係数行列:[Hx]
             chxji=b(i)/24d0            !係数行列:[Hx]T
c             chyij=c(j)/24d0            !係数行列:[Hy]
             chyji=c(i)/24d0            !係数行列:[Hy]T
c             chzij=d(j)/24d0            !係数行列:[Hz]
             chzji=d(i)/24d0            !係数行列:[Hz]T
c-----------------------------------------------------------------
c
             cxx=(uvwp(4*ne(im,i)-3)+usum)*b(j)/120d0  !係数行列:[Cxx]←対流項
             cyy=(uvwp(4*ne(im,i)-2)+vsum)*c(j)/120d0  !係数行列:[Cyy]←対流項
             czz=(uvwp(4*ne(im,i)-1)+wsum)*d(j)/120d0  !係数行列:[Czz]←対流項
c
c-----------------------------------------------------------------
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im',im,'i',i,'j',j
c         write(*,*)'cmm',sngl(cmm)
c         write(*,*)'sxx,syy,szz',sngl(sxx),sngl(syy),sngl(szz)
c         write(*,*)'sxy,syx',sngl(sxy),sngl(syx)
c         write(*,*)'syz,szy',sngl(syz),sngl(szy)
c         write(*,*)'szx,sxz',sngl(szx),sngl(sxz)
c       write(*,*)'chxji,chyji,chzji',sngl(chxji),sngl(chyji),sngl(chzji)
c         write(*,*)'cxx,cyy,czz',sngl(cxx),sngl(cyy),sngl(czz)
c         write(*,*)'---------------------------------------------------'
c         pause
c
             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
c
             b11=(2d0*sxx+syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式X成分:[U]の距離の偏微分項
             b12=syx*rey*vis/den                             !運動量収支式X成分:[V]の距離の偏微分項
             b13=szx*rey*vis/den                             !運動量収支式X成分:[W]の距離の偏微分項
             b14=-chxji/den                                  !運動量収支式X成分:[P]の距離の偏微分項
c
             b21=sxy*rey*vis/den                             !運動量収支式Y成分:[U]の距離の偏微分項
             b22=(sxx+2d0*syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Y成分:[V]の距離の偏微分項
             b23=szy*rey*vis/den                             !運動量収支式Y成分:[W]の距離の偏微分項
             b24=-chyji/den                                  !運動量収支式Y成分:[P]の距離の偏微分項
c
             b31=sxz*rey*vis/den                             !運動量収支式Z成分:[U]の距離の偏微分項
             b32=syz*rey*vis/den                             !運動量収支式Z成分:[V]の距離の偏微分項
             b33=(sxx+syy+2d0*szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Z成分:[W]の距離の偏微分項
             b34=-chzji/den                                  !運動量収支式Z成分:[P]の距離の偏微分項
c
c         write(1,*)'---------------------------------------------------'
c         write(1,*)'im',im,'i',i,'j',j
c         write(1,*)'aii/dt',sngl(a11/dt),sngl(a22/dt),sngl(a33/dt)
c         write(1,*)'b1i',sngl(b11),sngl(b12),sngl(b13),sngl(b14)
c         write(1,*)'b2i',sngl(b21),sngl(b22),sngl(b23),sngl(b24)
c         write(1,*)'b3i',sngl(b31),sngl(b32),sngl(b33),sngl(b34)
c         write(1,*)'---------------------------------------------------'
c         pause
c
c            運動量収支式X成分:iu行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11   ![U]の項のj成分
             sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12   ![V]の項のj成分
             sa(iu,nbw+jw-iu)=sa(iu,nbw+jw-iu)       +b13   ![W]の項のj成分
             sf(iu)=sf(iu)+a11/dt*uvwp0(4*ne(im,j)-3)
     &                    -b14*uvwp(4*ne(im,j)-0)       !圧力項b14{P}τ+⊿τを右辺に移項
c
c            運動量収支式Y成分:iv行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21   ![U]の項のj成分
             sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22   ![V]の項のj成分
             sa(iv,nbw+jw-iv)=sa(iv,nbw+jw-iv)       +b23   ![W]の項のj成分
             sf(iv)=sf(iv)+a22/dt*uvwp0(4*ne(im,j)-2)
     &                    -b24*uvwp(4*ne(im,j)-0)       !圧力項b24{P}τ+⊿τを右辺に移項
c
c            運動量収支式Z成分:iw行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iw,nbw+ju-iw)=sa(iw,nbw+ju-iw)       +b31   ![U]の項のj成分
             sa(iw,nbw+jv-iw)=sa(iw,nbw+jv-iw)       +b32   ![V]の項のj成分
             sa(iw,nbw+jw-iw)=sa(iw,nbw+jw-iw)+a33/dt+b33   ![W]の項のj成分
             sf(iw)=sf(iw)+a33/dt*uvwp0(4*ne(im,j)-1)
     &                    -b34*uvwp(4*ne(im,j)-0)       !圧力項b34{P}τ+⊿τを右辺に移項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)im,j,sngl(uvwp0(4*ne(im,j)-3))
c     &              ,sngl(uvwp0(4*ne(im,j)-2)),sngl(uvwp0(4*ne(im,j)-1))
c         write(*,*)im,j,sngl(uvwp(4*ne(im,j)-3))
c     &              ,sngl(uvwp(4*ne(im,j)-2)),sngl(uvwp(4*ne(im,j)-1))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iu,nbw+ju-iu))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iu,nbw+jv-iu))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iu,nbw+jw-iu))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iv,nbw+ju-iv))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iv,nbw+jv-iv))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iv,nbw+jw-iv))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iw,nbw+ju-iw))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iw,nbw+jv-iw))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iw,nbw+jw-iw))
c         write(*,*)'---------------------------------------------------'
c         pause
290      continue
200    continue
c
c       close(1)
c       do i=1,4*3
c       do j=1,4*3
c       write(*,*)'sa',i,j,sngl(sa(i,j-i+nbw))
c       pause
c       enddo
c       enddo
c       pause
c
       return
       end
c
c      u0(np),v0(np),w0(np)を用いて連続の式からpp(np)を計算する.sa(ind,ibw),sf(ind)
c
       subroutine flowcn(nel,nbw,sa,sf,ne,xx,yy,zz    !ne(ine,4):流動領域の要素のみ格納
     &                 ,uvwp,uvwp0,t0,idx)              !nbw=1*nbw
       include "header.h"
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)
       dimension uvwp(4*ind),uvwp0(4*ind)
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension x(4),y(4),z(4),b(4),c(4),d(4),zkesu(4,4)
       dimension t0(ind)    !物性値を求める時に用いる
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
c       common /phyc/dens(1001),visc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
       do 200 im=1,nel     !要素imについて
         id=idx(im)
c         usum=0d0
c         vsum=0d0
c         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           x(i)=xx(ne(im,i))
           y(i)=yy(ne(im,i))
           z(i)=zz(ne(im,i))
c           usum=usum+uvwp0(4*ne(im,i)-3)
c           vsum=vsum+uvwp0(4*ne(im,i)-2)
c           wsum=wsum+uvwp0(4*ne(im,i)-1)
210      continue
c
c        係数行列の成分を算出する
c
         b(1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
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
         if(vl.le.0d0)then
           write(*,*)'vl',sngl(vl)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
c        無次元数および物性値を求める
c
c         rey=re(id)     !1/Re
         svc=sv(id)**2  !(1/Ma)^2
c         ggy=gv(id)     !1/Fr
c         vis=(calvis(t0(ne(im,1)))+calvis(t0(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
c     &       +calvis(t0(ne(im,3)))+calvis(t0(ne(im,4))))
c     &       /4d0/calvis(-999d0)
c         den=(calden(t0(ne(im,1)))+calden(t0(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
c     &       +calden(t0(ne(im,3)))+calden(t0(ne(im,4))))
c     &       /4d0/calden(-999d0)
c         write(*,*)'svc',sngl(svc)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,id',im,id
c         write(*,*)'dt,vl',sngl(dt),sngl(vl)
c         write(*,*)'usum,vsum,wsum',sngl(usum),sngl(vsum),sngl(wsum)
c         write(*,*)'b(i)',sngl(b(1)),sngl(b(2)),sngl(b(3)),sngl(b(4))
c         write(*,*)'c(i)',sngl(c(1)),sngl(c(2)),sngl(c(3)),sngl(c(4))
c         write(*,*)'d(i)',sngl(d(1)),sngl(d(2)),sngl(d(3)),sngl(d(4))
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c         write(*,*)'vis,den',sngl(vis),sngl(den)
c         write(*,*)'---------------------------------------------------'
c         pause
c
         do 290 i=1,4     !要素imの行
           ip=ne(im,i)
           do 290 j=1,4   !要素imの列
             jp=ne(im,j)
c             
             cmm=vl/20d0*zkesu(i,j)     !係数行列:[C]
c             sxx=b(i)*b(j)/36d0/vl      !係数行列:[Sxx]
c             syy=c(i)*c(j)/36d0/vl      !係数行列:[Syy]
c             szz=d(i)*d(j)/36d0/vl      !係数行列:[Szz]
c             sxy=b(i)*c(j)/36d0/vl      !係数行列:[Sxy]
c             syx=c(i)*b(j)/36d0/vl      !係数行列:[Syx]
c             syz=c(i)*d(j)/36d0/vl      !係数行列:[Syz]
c             szy=d(i)*c(j)/36d0/vl      !係数行列:[Szy]
c             szx=d(i)*b(j)/36d0/vl      !係数行列:[Szx]
c             sxz=b(i)*d(j)/36d0/vl      !係数行列:[Sxz]
             chxij=b(j)/24d0            !係数行列:[Hx]
c             chxji=b(i)/24d0            !係数行列:[Hx]T
             chyij=c(j)/24d0            !係数行列:[Hy]
c             chyji=c(i)/24d0            !係数行列:[Hy]T
             chzij=d(j)/24d0            !係数行列:[Hz]
c             chzji=d(i)/24d0            !係数行列:[Hz]T
c             cxx=(uvwp0(4*ne(im,i)-3)+usum)*b(j)/120d0  !係数行列:[Cxx]←対流項
c             cyy=(uvwp0(4*ne(im,i)-2)+vsum)*c(j)/120d0  !係数行列:[Cyy]←対流項
c             czz=(uvwp0(4*ne(im,i)-1)+wsum)*d(j)/120d0  !係数行列:[Czz]←対流項
c             
             a44=cmm                    !連続の式:[P]の時間の偏微分項
c
             b41=svc*chxij              !連続の式:[U]の距離の偏微分項
             b42=svc*chyij              !連続の式:[V]の距離の偏微分項
             b43=svc*chzij              !連続の式:[W]の距離の偏微分項
             b44=0d0                    !連続の式:[P]の距離の偏微分項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,i,j',im,i,j
c         write(*,*)'cmm',sngl(cmm)
c       write(*,*)'chxij,chyij,chzij',sngl(chxij),sngl(chyij),sngl(chzij)
c         write(*,*)'a44',sngl(a44)
c         write(*,*)'b4i',sngl(b41),sngl(b42),sngl(b43),sngl(b44)
c         write(*,*)'uvwp0',sngl(uvwp0(4*ne(im,j)-3)),
c     &                      sngl(uvwp0(4*ne(im,j)-2)),
c     &                      sngl(uvwp0(4*ne(im,j)-1)),
c     &                           uvwp0(4*ne(im,j)-0)
c         write(*,*)'uvwp',sngl(uvwp(4*ne(im,j)-3)),
c     &                     sngl(uvwp(4*ne(im,j)-2)),
c     &                     sngl(uvwp(4*ne(im,j)-1)),
c     &                          uvwp(4*ne(im,j)-0)
c         write(*,*)'---------------------------------------------------'
c         pause
c
c            連続の式
             sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a44/dt+b44
             sf(ip)=sf(ip)+a44/dt*uvwp0(4*jp-0)-b41*uvwp(4*ne(im,j)-3)     !{U},{V},{W}の項を右辺に移項
     &                                         -b42*uvwp(4*ne(im,j)-2)
     &                                         -b43*uvwp(4*ne(im,j)-1)
290      continue
200    continue
c
       return
       end
c------------------------------------------------------------
c
c      SIMPLE法でマトリクスを解く場合
c      →オイラー座標(形状関数の値:一定)
c------------------------------------------------------------
c
c      p0(np)を用いてN-S式からuu(np),vv(np),ww(np)のマトリクスsa(3*ind,ibw),sf(3*ind)を組み立てる
c
       subroutine flowns0(nel,nbw,sa,sf,ne,xx,yy,zz    !ne(ine,4):新_新
     &                  ,uvwp,uvwp0,t0,idx,imt)          !nbw=3*nbw
       include "header.h"
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)  !新_新
       dimension uvwp(4*ind),uvwp0(4*ind)
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension x(4),y(4),z(4),b(ine,4),c(ine,4),d(ine,4),zkesu(4,4)
       dimension vl(ine)
       dimension t0(ind)    !物性値を求める時に用いる
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
c       common /phyc/dens(1001),visc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
c       open(1,file='sa.res')
c       rewind(1)
c
       do 200 im=1,nel     !nel:新_新
         id=idx(im)
         usum=0d0
         vsum=0d0
         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           if(imt.eq.1)x(i)=xx(ne(im,i))
           if(imt.eq.1)y(i)=yy(ne(im,i))
           if(imt.eq.1)z(i)=zz(ne(im,i))
c-----------------------------------------------------------------
c
           usum=usum+uvwp(4*ne(im,i)-3)    !uvwp(4*新_新-i)
           vsum=vsum+uvwp(4*ne(im,i)-2)
           wsum=wsum+uvwp(4*ne(im,i)-1)
c
c-----------------------------------------------------------------
210      continue
c*****************************************************************
c
c        初回のみ計算
c
c*****************************************************************
         if(imt.eq.1)then     !0:形状関数を計算しない,1:形状関数を計算する
c
c        係数行列の成分を算出する
c
         b(im,1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(im,2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(im,3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(im,4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(im,1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(im,2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(im,3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(im,4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(im,1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(im,2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(im,3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(im,4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
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
         vl(im)=dabs(det)/6d0
         if(vl(im).le.0d0)then
           write(*,*)'im',im
           write(*,*)'vl',sngl(vl(im))
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
         endif
c*****************************************************************
c
c        無次元数および物性値を求める
c
         rey=re(id)     !1/Re
c         svc=sv(id)**2  !(1/Ma)^2
         ggy=gv(id)     !1/Fr
         vis=(calvis(t0(ne(im,1)))+calvis(t0(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
     &       +calvis(t0(ne(im,3)))+calvis(t0(ne(im,4))))
     &       /4d0/calvis(-999d0)
         den=(calden(t0(ne(im,1)))+calden(t0(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
     &       +calden(t0(ne(im,3)))+calden(t0(ne(im,4))))
     &       /4d0/calden(-999d0)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,id',im,id
c         write(*,*)'dt,vl(im)',sngl(dt),sngl(vl(im))
c         write(*,*)'usum,vsum,wsum',sngl(usum),sngl(vsum),sngl(wsum)
c         write(*,*)'b(im,i)',sngl(b(im,1)),sngl(b(im,2))
c     &                       ,sngl(b(im,3)),sngl(b(im,4))
c         write(*,*)'c(im,i)',sngl(c(im,1)),sngl(c(im,2))
c     &                       ,sngl(c(im,3)),sngl(c(im,4))
c         write(*,*)'d(im,i)',sngl(d(im,1)),sngl(d(im,2))
c     &                       ,sngl(d(im,3)),sngl(d(im,4))
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c         write(*,*)'vis,den',sngl(vis),sngl(den)
c         write(*,*)'---------------------------------------------------'
c         pause
c
         do 280 i=1,4   !重力考慮:-y軸方向
280      sf(ne(im,i)*3-1)=sf(ne(im,i)*3-1)+ggy*vl(im)/4d0
c
         do 290 i=1,4     !要素imの行
           iu=ne(im,i)*3-2
           iv=ne(im,i)*3-1
           iw=ne(im,i)*3-0
           do 290 j=1,4   !要素imの列
             ju=ne(im,j)*3-2
             jv=ne(im,j)*3-1
             jw=ne(im,j)*3-0
c             
             cmm=vl(im)/20d0*zkesu(i,j)           !係数行列:[C]
             sxx=b(im,i)*b(im,j)/36d0/vl(im)      !係数行列:[Sxx]
             syy=c(im,i)*c(im,j)/36d0/vl(im)      !係数行列:[Syy]
             szz=d(im,i)*d(im,j)/36d0/vl(im)      !係数行列:[Szz]
             sxy=b(im,i)*c(im,j)/36d0/vl(im)      !係数行列:[Sxy]
             syx=c(im,i)*b(im,j)/36d0/vl(im)      !係数行列:[Syx]
             syz=c(im,i)*d(im,j)/36d0/vl(im)      !係数行列:[Syz]
             szy=d(im,i)*c(im,j)/36d0/vl(im)      !係数行列:[Szy]
             szx=d(im,i)*b(im,j)/36d0/vl(im)      !係数行列:[Szx]
             sxz=b(im,i)*d(im,j)/36d0/vl(im)      !係数行列:[Sxz]
c             chxij=b(im,j)/24d0                   !係数行列:[Hx]
             chxji=b(im,i)/24d0                   !係数行列:[Hx]T
c             chyij=c(im,j)/24d0                   !係数行列:[Hy]
             chyji=c(im,i)/24d0                   !係数行列:[Hy]T
c             chzij=d(im,j)/24d0                   !係数行列:[Hz]
             chzji=d(im,i)/24d0                   !係数行列:[Hz]T
c-----------------------------------------------------------------
c
             cxx=(uvwp(4*ne(im,i)-3)+usum)*b(im,j)/120d0  !係数行列:[Cxx]←対流項
             cyy=(uvwp(4*ne(im,i)-2)+vsum)*c(im,j)/120d0  !係数行列:[Cyy]←対流項
             czz=(uvwp(4*ne(im,i)-1)+wsum)*d(im,j)/120d0  !係数行列:[Czz]←対流項
c
c-----------------------------------------------------------------
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im',im,'i',i,'j',j
c         write(*,*)'cmm',sngl(cmm)
c         write(*,*)'sxx,syy,szz',sngl(sxx),sngl(syy),sngl(szz)
c         write(*,*)'sxy,syx',sngl(sxy),sngl(syx)
c         write(*,*)'syz,szy',sngl(syz),sngl(szy)
c         write(*,*)'szx,sxz',sngl(szx),sngl(sxz)
c       write(*,*)'chxji,chyji,chzji',sngl(chxji),sngl(chyji),sngl(chzji)
c         write(*,*)'cxx,cyy,czz',sngl(cxx),sngl(cyy),sngl(czz)
c         write(*,*)'---------------------------------------------------'
c         pause
c
             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
c
             b11=(2d0*sxx+syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式X成分:[U]の距離の偏微分項
             b12=syx*rey*vis/den                             !運動量収支式X成分:[V]の距離の偏微分項
             b13=szx*rey*vis/den                             !運動量収支式X成分:[W]の距離の偏微分項
             b14=-chxji/den                                  !運動量収支式X成分:[P]の距離の偏微分項
c
             b21=sxy*rey*vis/den                             !運動量収支式Y成分:[U]の距離の偏微分項
             b22=(sxx+2d0*syy+szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Y成分:[V]の距離の偏微分項
             b23=szy*rey*vis/den                             !運動量収支式Y成分:[W]の距離の偏微分項
             b24=-chyji/den                                  !運動量収支式Y成分:[P]の距離の偏微分項
c
             b31=sxz*rey*vis/den                             !運動量収支式Z成分:[U]の距離の偏微分項
             b32=syz*rey*vis/den                             !運動量収支式Z成分:[V]の距離の偏微分項
             b33=(sxx+syy+2d0*szz)*rey*vis/den+cxx+cyy+czz   !運動量収支式Z成分:[W]の距離の偏微分項
             b34=-chzji/den                                  !運動量収支式Z成分:[P]の距離の偏微分項
c
c         write(1,*)'---------------------------------------------------'
c         write(1,*)'im',im,'i',i,'j',j
c         write(1,*)'aii/dt',sngl(a11/dt),sngl(a22/dt),sngl(a33/dt)
c         write(1,*)'b1i',sngl(b11),sngl(b12),sngl(b13),sngl(b14)
c         write(1,*)'b2i',sngl(b21),sngl(b22),sngl(b23),sngl(b24)
c         write(1,*)'b3i',sngl(b31),sngl(b32),sngl(b33),sngl(b34)
c         write(1,*)'---------------------------------------------------'
c         pause
c
c            運動量収支式X成分:iu行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11   ![U]の項のj成分
             sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12   ![V]の項のj成分
             sa(iu,nbw+jw-iu)=sa(iu,nbw+jw-iu)       +b13   ![W]の項のj成分
             sf(iu)=sf(iu)+a11/dt*uvwp0(4*ne(im,j)-3)
     &                    -b14*uvwp(4*ne(im,j)-0)       !圧力項b14{P}τ+⊿τを右辺に移項
c
c            運動量収支式Y成分:iv行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21   ![U]の項のj成分
             sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22   ![V]の項のj成分
             sa(iv,nbw+jw-iv)=sa(iv,nbw+jw-iv)       +b23   ![W]の項のj成分
             sf(iv)=sf(iv)+a22/dt*uvwp0(4*ne(im,j)-2)
     &                    -b24*uvwp(4*ne(im,j)-0)       !圧力項b24{P}τ+⊿τを右辺に移項
c
c            運動量収支式Z成分:iw行の節点ne(im,j)の速度成分ju,jv,jw
             sa(iw,nbw+ju-iw)=sa(iw,nbw+ju-iw)       +b31   ![U]の項のj成分
             sa(iw,nbw+jv-iw)=sa(iw,nbw+jv-iw)       +b32   ![V]の項のj成分
             sa(iw,nbw+jw-iw)=sa(iw,nbw+jw-iw)+a33/dt+b33   ![W]の項のj成分
             sf(iw)=sf(iw)+a33/dt*uvwp0(4*ne(im,j)-1)
     &                    -b34*uvwp(4*ne(im,j)-0)       !圧力項b34{P}τ+⊿τを右辺に移項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)im,j,sngl(uvwp0(4*ne(im,j)-3))
c     &              ,sngl(uvwp0(4*ne(im,j)-2)),sngl(uvwp0(4*ne(im,j)-1))
c         write(*,*)im,j,sngl(uvwp(4*ne(im,j)-3))
c     &              ,sngl(uvwp(4*ne(im,j)-2)),sngl(uvwp(4*ne(im,j)-1))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iu,nbw+ju-iu))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iu,nbw+jv-iu))
c         write(*,*)'sa(3*',ne(im,i),'-2,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iu,nbw+jw-iu))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iv,nbw+ju-iv))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iv,nbw+jv-iv))
c         write(*,*)'sa(3*',ne(im,i),'-1,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iv,nbw+jw-iv))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-2)',
c     &                                         sngl(sa(iw,nbw+ju-iw))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-1)',
c     &                                         sngl(sa(iw,nbw+jv-iw))
c         write(*,*)'sa(3*',ne(im,i),'-0,3*',ne(im,j),'-0)',
c     &                                         sngl(sa(iw,nbw+jw-iw))
c         write(*,*)'---------------------------------------------------'
c         pause
290      continue
200    continue
c
c       close(1)
c       do i=1,4*3
c       do j=1,4*3
c       write(*,*)'sa',i,j,sngl(sa(i,j-i+nbw))
c       pause
c       enddo
c       enddo
c       pause
c
       return
       end
c
c      u0(np),v0(np),w0(np)を用いて連続の式からpp(np)を計算する.sa(ind,ibw),sf(ind)
c
       subroutine flowcn0(nel,nbw,sa,sf,ne,xx,yy,zz    !ne(ine,4):流動領域の要素のみ格納
     &                 ,uvwp,uvwp0,t0,idx,imt)              !nbw=1*nbw
       include "header.h"
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)
       dimension uvwp(4*ind),uvwp0(4*ind)
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension x(4),y(4),z(4),b(4),c(4),d(4),zkesu(4,4)
       dimension t0(ind)    !物性値を求める時に用いる
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
c       common /phyc/dens(1001),visc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
       do 200 im=1,nel     !要素imについて
         id=idx(im)
c         usum=0d0
c         vsum=0d0
c         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           x(i)=xx(ne(im,i))
           y(i)=yy(ne(im,i))
           z(i)=zz(ne(im,i))
c           usum=usum+uvwp0(4*ne(im,i)-3)
c           vsum=vsum+uvwp0(4*ne(im,i)-2)
c           wsum=wsum+uvwp0(4*ne(im,i)-1)
210      continue
c
c        係数行列の成分を算出する
c
         b(1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
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
         if(vl.le.0d0)then
           write(*,*)'vl',sngl(vl)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
c        無次元数および物性値を求める
c
c         rey=re(id)     !1/Re
         svc=sv(id)**2  !(1/Ma)^2
c         ggy=gv(id)     !1/Fr
c         vis=(calvis(t0(ne(im,1)))+calvis(t0(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
c     &       +calvis(t0(ne(im,3)))+calvis(t0(ne(im,4))))
c     &       /4d0/calvis(-999d0)
c         den=(calden(t0(ne(im,1)))+calden(t0(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
c     &       +calden(t0(ne(im,3)))+calden(t0(ne(im,4))))
c     &       /4d0/calden(-999d0)
c         write(*,*)'svc',sngl(svc)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,id',im,id
c         write(*,*)'dt,vl',sngl(dt),sngl(vl)
c         write(*,*)'usum,vsum,wsum',sngl(usum),sngl(vsum),sngl(wsum)
c         write(*,*)'b(i)',sngl(b(1)),sngl(b(2)),sngl(b(3)),sngl(b(4))
c         write(*,*)'c(i)',sngl(c(1)),sngl(c(2)),sngl(c(3)),sngl(c(4))
c         write(*,*)'d(i)',sngl(d(1)),sngl(d(2)),sngl(d(3)),sngl(d(4))
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c         write(*,*)'vis,den',sngl(vis),sngl(den)
c         write(*,*)'---------------------------------------------------'
c         pause
c
         do 290 i=1,4     !要素imの行
           ip=ne(im,i)
           do 290 j=1,4   !要素imの列
             jp=ne(im,j)
c             
             cmm=vl/20d0*zkesu(i,j)     !係数行列:[C]
c             sxx=b(i)*b(j)/36d0/vl      !係数行列:[Sxx]
c             syy=c(i)*c(j)/36d0/vl      !係数行列:[Syy]
c             szz=d(i)*d(j)/36d0/vl      !係数行列:[Szz]
c             sxy=b(i)*c(j)/36d0/vl      !係数行列:[Sxy]
c             syx=c(i)*b(j)/36d0/vl      !係数行列:[Syx]
c             syz=c(i)*d(j)/36d0/vl      !係数行列:[Syz]
c             szy=d(i)*c(j)/36d0/vl      !係数行列:[Szy]
c             szx=d(i)*b(j)/36d0/vl      !係数行列:[Szx]
c             sxz=b(i)*d(j)/36d0/vl      !係数行列:[Sxz]
             chxij=b(j)/24d0            !係数行列:[Hx]
c             chxji=b(i)/24d0            !係数行列:[Hx]T
             chyij=c(j)/24d0            !係数行列:[Hy]
c             chyji=c(i)/24d0            !係数行列:[Hy]T
             chzij=d(j)/24d0            !係数行列:[Hz]
c             chzji=d(i)/24d0            !係数行列:[Hz]T
c             cxx=(uvwp0(4*ne(im,i)-3)+usum)*b(j)/120d0  !係数行列:[Cxx]←対流項
c             cyy=(uvwp0(4*ne(im,i)-2)+vsum)*c(j)/120d0  !係数行列:[Cyy]←対流項
c             czz=(uvwp0(4*ne(im,i)-1)+wsum)*d(j)/120d0  !係数行列:[Czz]←対流項
c             
             a44=cmm                    !連続の式:[P]の時間の偏微分項
c
             b41=svc*chxij              !連続の式:[U]の距離の偏微分項
             b42=svc*chyij              !連続の式:[V]の距離の偏微分項
             b43=svc*chzij              !連続の式:[W]の距離の偏微分項
             b44=0d0                    !連続の式:[P]の距離の偏微分項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,i,j',im,i,j
c         write(*,*)'cmm',sngl(cmm)
c       write(*,*)'chxij,chyij,chzij',sngl(chxij),sngl(chyij),sngl(chzij)
c         write(*,*)'a44',sngl(a44)
c         write(*,*)'b4i',sngl(b41),sngl(b42),sngl(b43),sngl(b44)
c         write(*,*)'uvwp0',sngl(uvwp0(4*ne(im,j)-3)),
c     &                      sngl(uvwp0(4*ne(im,j)-2)),
c     &                      sngl(uvwp0(4*ne(im,j)-1)),
c     &                           uvwp0(4*ne(im,j)-0)
c         write(*,*)'uvwp',sngl(uvwp(4*ne(im,j)-3)),
c     &                     sngl(uvwp(4*ne(im,j)-2)),
c     &                     sngl(uvwp(4*ne(im,j)-1)),
c     &                          uvwp(4*ne(im,j)-0)
c         write(*,*)'---------------------------------------------------'
c         pause
c
c            連続の式
             sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a44/dt+b44
             sf(ip)=sf(ip)+a44/dt*uvwp0(4*jp-0)-b41*uvwp(4*ne(im,j)-3)     !{U},{V},{W}の項を右辺に移項
     &                                         -b42*uvwp(4*ne(im,j)-2)
     &                                         -b43*uvwp(4*ne(im,j)-1)
290      continue
200    continue
c
       return
       end
