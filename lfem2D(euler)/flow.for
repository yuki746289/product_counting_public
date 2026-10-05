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
c-----------------------------------------------------------------------
c
c      マトリックスの作成:連続の式
c
c-----------------------------------------------------------------------
       subroutine flowcn(nelw,nbw,new,idw,xn,yn,uvp0,uvp,dt
     &                                     ,vre,vma,vfr,sa,sf)
       include "head.for"
       dimension new(ine,3),idw(0:ine)
       dimension xn(ind),yn(ind)
       dimension uvp0(3*ind),uvp(3*ind)
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
           ip=new(im,i)
           do 290 j=1,3   !要素imの列
             jp=new(im,j)
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
c             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
c             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
c
c             b11=(2d0*sxx+syy)*rey  !運動量収支式X成分:[U]の距離の偏微分項
c             b12=syx*rey            !運動量収支式X成分:[V]の距離の偏微分項
c             b13=-chxji             !運動量収支式X成分:[P]の距離の偏微分項
c             
c             b21=sxy*rey            !運動量収支式Y成分:[U]の距離の偏微分項
c             b22=(sxx+2d0*syy)*rey  !運動量収支式Y成分:[V]の距離の偏微分項
c             b23=-chyji             !運動量収支式Y成分:[P]の距離の偏微分項
c             
             b31=svc*chxij          !連続の式:[U]の距離の偏微分項
             b32=svc*chyij          !連続の式:[V]の距離の偏微分項
             b33=0d0                !連続の式:[P]の距離の偏微分項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)im,j,sngl(uvp0(1*new(im,j)-0))
c         write(*,*)im,j,sngl(uvp(1*new(im,j)-0))
c         write(*,*)'sa(1*',new(im,i),'-0,1*',new(im,j),'-0)',
c     &                                         sngl(sa(ip,nbw+jp-ip))
c         write(*,*)'sf(1*',new(im,i),'-0)',sngl(sf(ip))
c         write(*,*)'---------------------------------------------------'
c         if(new(im,i).eq.109)pause
c         if(new(im,j).eq.109)pause
c         pause
c
c            連続の式
             sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a33/dt+b33
             sf(ip)=sf(ip)+a33/dt*uvp0(3*jp-0)-b31*uvp(3*new(im,j)-2)     !{U},{V},{W}の項を右辺に移項
     &                                        -b32*uvp(3*new(im,j)-1)
c             write(*,*)"sf:",sngl(a33/dt*uvp0(3*jp-0))
c     &                       ,sngl(b31*uvp(3*new(im,j)-2))
c     &                       ,sngl(b32*uvp(3*new(im,j)-1))
290      continue
200    continue
c
       return
       end
c-----------------------------------------------------------------------
c
c      マトリックスの作成:運動量収支式
c
c-----------------------------------------------------------------------
       subroutine flowns0(nelw,nbw,new,idw,xn,yn,uvp0,uvp
     &                                       ,dt,vre,vma,vfr,sa,sf)
       include "head.for"
       dimension new(ine,3),idw(0:ine)
       dimension xn(ind),yn(ind)
       dimension uvp0(3*ind),uvp(3*ind)
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
         usum=0d0
         vsum=0d0
         do 210 i=1,3      !要素imを構成する節点の座標
           x(i)=xn(new(im,i))
           y(i)=yn(new(im,i))
c
           usum=usum+uvp(3*new(im,i)-2)    !uvwp(4*新_新-i)
           vsum=vsum+uvp(3*new(im,i)-1)
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
           iu=new(im,i)*2-1
           iv=new(im,i)*2-0
           do 290 j=1,3   !要素imの列
             ju=new(im,j)*2-1
             jv=new(im,j)*2-0
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
             cxx=(uvp(3*new(im,i)-2)+usum)*b(j)/30d0  !係数行列:[Cxx]←対流項
             cyy=(uvp(3*new(im,i)-1)+vsum)*c(j)/30d0  !係数行列:[Cyy]←対流項
c             
             a11=cmm                    !運動量収支式X成分:[U]の時間の偏微分項
             a22=cmm                    !運動量収支式Y成分:[V]の時間の偏微分項
             a33=cmm                    !運動量収支式Z成分:[W]の時間の偏微分項
c
             b11=(2d0*sxx+syy)*rey+cxx+cyy  !運動量収支式X成分:[U]の距離の偏微分項
             b12=syx*rey            !運動量収支式X成分:[V]の距離の偏微分項
             b13=-chxji             !運動量収支式X成分:[P]の距離の偏微分項
c             
             b21=sxy*rey            !運動量収支式Y成分:[U]の距離の偏微分項
             b22=(sxx+2d0*syy)*rey+cxx+cyy  !運動量収支式Y成分:[V]の距離の偏微分項
             b23=-chyji             !運動量収支式Y成分:[P]の距離の偏微分項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im',im,'i',i,'j',j
c         write(*,*)'aii/dt',sngl(a11/dt),sngl(a22/dt),sngl(a33/dt)
c         write(*,*)'b1i',sngl(b11),sngl(b12),sngl(b13)
c         write(*,*)'b2i',sngl(b21),sngl(b22),sngl(b23)
c         write(*,*)'b3i',sngl(b31),sngl(b32),sngl(b33)
c         write(*,*)'---------------------------------------------------'
c         pause
c
c            運動量収支式X成分
             sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11   ![U]の項のj成分
             sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12   ![V]の項のj成分
c             sf(iu)=sf(iu)+a11/dt*uvp0(ju)-b13*uvp(3*new(im,j)-0)
             sf(iu)=sf(iu)+a11/dt*uvp0(3*new(im,j)-2)
     &                    -b13*uvp(3*new(im,j)-0)
c
c            運動量収支式Y成分
             sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21   ![U]の項のj成分
             sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22   ![V]の項のj成分
c             sf(iv)=sf(iv)+a22/dt*uvp0(jv)-b23*uvp(3*new(im,j)-0)
             sf(iv)=sf(iv)+a22/dt*uvp0(3*new(im,j)-1)
     &                    -b23*uvp(3*new(im,j)-0)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)im,j,sngl(uvp0(3*new(im,j)-2))
c     &                  ,sngl(uvp0(3*new(im,j)-1))
c     &                  ,sngl(uvp0(3*new(im,j)-0))
c         write(*,*)im,j,sngl(uvp(3*new(im,j)-2))
c     &                  ,sngl(uvp(3*new(im,j)-1))
c     &                  ,sngl(uvp(3*new(im,j)-0))
c         write(*,*)'sa(2*',new(im,i),'-1,2*',new(im,j),'-1)',
c     &                                         sngl(sa(iu,nbw+ju-iu))
c         write(*,*)'sa(2*',new(im,i),'-1,2*',new(im,j),'-0)',
c     &                                         sngl(sa(iu,nbw+jv-iu))
c         write(*,*)'sa(2*',new(im,i),'-0,2*',new(im,j),'-1)',
c     &                                         sngl(sa(iv,nbw+ju-iv))
c         write(*,*)'sa(2*',new(im,i),'-0,2*',new(im,j),'-0)',
c     &                                         sngl(sa(iv,nbw+jv-iv))
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im:',im,nelw
c         write(*,*)'sf(2*',new(im,i),'-1)',sngl(sf(iu))
c         write(*,*)'sf(2*',new(im,i),'-0)',sngl(sf(iv))
c         pause
c         write(*,*)"a:",sngl(a11/dt*uvp0(3*new(im,j)-2))
c     &                  ,sngl(-b13*uvp(3*new(im,j)-0))
290      continue
200    continue
c
       return
       end
