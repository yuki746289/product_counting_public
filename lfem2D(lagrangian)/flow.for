c-----------------------------------------------------------------------
c
c      マトリックスの初期化
c
c-----------------------------------------------------------------------
       subroutine setmat(nw,nbw,sa,sf,dsf)
       include "head.for"
       dimension sa(ind,ibw),sf(ind),dsf(ind)
c
       if(nw*3.gt.ind .or. nbw*6.gt.ibw)then
         write(*,*)"setmat"
         write(*,*)"ind:",nw*3,ind
         write(*,*)"ibw:",nbw*6,ibw
         stop
       endif
c
       do i=1,nw*3
         sf(i)=0d0
         dsf(i)=0d0
         do j=1,nbw*6
           sa(i,j)=0d0
         enddo
       enddo
c
       return
       end
c-----------------------------------------------------------------------
c
c      マトリックスの作成
c
c-----------------------------------------------------------------------
       subroutine flow(nelw,nbw,new,idw,xn,yn,uvp0,dt,vre,vma,vfr,sa,sf)
       include "head.for"
       dimension new(ine,3),idw(0:ine)
       dimension xn(ind),yn(ind)
       dimension uvp0(3*ind)
       dimension vre(-5:ica),vma(-5:ica),vfr(-5:ica)
       dimension sa(ind,ibw),sf(ind)
c
       dimension x(3),y(3),z(3)
       dimension b(3),c(3),d(3)
       dimension zkesu(3,3)
       data zkesu/2d0,1d0,1d0,
     &            1d0,2d0,1d0,
     &            1d0,1d0,2d0/
c
       do 200 im=1,nelw     !要素imについて
         id=idw(im)
         do 210 i=1,3      !要素imを構成する節点の座標
           x(i)=xn(new(im,i))
           y(i)=yn(new(im,i))
210      continue
c
c        係数行列の成分を算出する
         b(1)=y(2)-y(3)
         b(2)=y(3)-y(1)
         b(3)=y(1)-y(2)
         c(1)=x(3)-x(2)
         c(2)=x(1)-x(3)
         c(3)=x(2)-x(1)
         ar=((x(1)-x(2))*(y(1)-y(3))-(x(1)-x(3))*(y(1)-y(2)))/2d0
         if(ar.le.0d0)then
           write(*,*)"flow"
           write(*,*)'ar',sngl(ar)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3))
           pause
         endif
c
c        無次元数および物性値を求める
         rey=1d0/vre(id)     !1/Re
         svc=1d0/vma(id)**2  !1/Ma^2
c         ggy=1d0/vfr(id)**2  !1/Fr^2
c         write(*,*)'rey,svc,ggy',sngl(rey),sngl(svc),sngl(ggy)
c
c        重力考慮:y軸方向
c         do 280 i=1,3   !重力考慮:y軸方向
c280      sf(new(im,i)*3-1)=sf(new(im,i)*3-1)+ggy*ar/3d0
c
c        マトリックス作成
         do 290 i=1,3     !要素imの行
           iu=new(im,i)*3-2
           iv=new(im,i)*3-1
           ip=new(im,i)*3-0
           do 290 j=1,3   !要素imの列
             ju=new(im,j)*3-2
             jv=new(im,j)*3-1
             jp=new(im,j)*3-0
c             
             cmm=ar/12d0*zkesu(i,j)     !係数行列:[C]
             sxx=b(i)*b(j)/4d0/ar       !係数行列:[Sxx]
             syy=c(i)*c(j)/4d0/ar       !係数行列:[Syy]
             sxy=b(i)*c(j)/4d0/ar       !係数行列:[Sxy]
             syx=c(i)*b(j)/4d0/ar       !係数行列:[Syx]
             chxij=b(j)/6d0             !係数行列:[Hx]
             chxji=b(i)/6d0             !係数行列:[Hx]T
             chyij=c(j)/6d0             !係数行列:[Hy]
             chyji=c(i)/6d0             !係数行列:[Hy]T
c             
             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
c
             b11=(2d0*sxx+syy)*rey  !運動量収支式X成分:[U]の距離の偏微分項
             b12=syx*rey            !運動量収支式X成分:[V]の距離の偏微分項
             b13=-chxji             !運動量収支式X成分:[P]の距離の偏微分項
c             
             b21=sxy*rey            !運動量収支式Y成分:[U]の距離の偏微分項
             b22=(sxx+2d0*syy)*rey  !運動量収支式Y成分:[V]の距離の偏微分項
             b23=-chyji             !運動量収支式Y成分:[P]の距離の偏微分項
c             
             b31=svc*chxij          !連続の式:[U]の距離の偏微分項
             b32=svc*chyij          !連続の式:[V]の距離の偏微分項
             b33=0d0                !連続の式:[P]の距離の偏微分項
c             
c            運動量収支式X成分
             sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11   ![U]の項のj成分
             sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12   ![V]の項のj成分
             sa(iu,nbw+jp-iu)=sa(iu,nbw+jp-iu)       +b13
             sf(iu)=sf(iu)+a11/dt*uvp0(ju)
c             
c            運動量収支式Y成分
             sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21   ![U]の項のj成分
             sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22   ![V]の項のj成分
             sa(iv,nbw+jp-iv)=sa(iv,nbw+jp-iv)       +b23
             sf(iv)=sf(iv)+a22/dt*uvp0(jv)
c
c            連続の式
             sa(ip,nbw+ju-ip)=sa(ip,nbw+ju-ip)       +b31
             sa(ip,nbw+jv-ip)=sa(ip,nbw+jv-ip)       +b32
             sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a33/dt+b33
             sf(ip)=sf(ip)+a33/dt*uvp0(jp)
290      continue
200    continue
c
       return
       end
