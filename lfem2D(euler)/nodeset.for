c-------------------------------------------------------------------------
c
c       グローバル変数
c       np           :節点数
c       xx(np)       :x座標
c       yy(np)       :y座標
c       zz(np)       :z座標
c       rr(np)       :節点半径
c       nb(np)        0:内部の節点,1:境界の節点
c       ns(np)       :周囲の節点数
c       ms(np,ns(np)):周囲の節点番号 移動ベクトル用
c       nt(np)       :周囲の節点数
c       mt(np,nt(np)):周囲の節点番号 節点の追加・削除用
c       vx,vy,vz     :移動ベクトル
c
c       ローカル変数
c       x0(ind),y0(ind),z0(ind),rr0(ind)
c       nb0(ind)
c       ns0(ind)
c       np0
c       kp
c       ic
c       i,j,k
c       dd
c       nd    :削除する節点数
c       md(nd):削除する節点番号
c       d,r
c
c      ******************************* 
c      z0[np0]=zz[i];			//2次元
c      z0[np0]=zz[i]+rr[i]/2;	//3次元
c
c      追加・削除の節点数
c      境界節点の節点半径
c      ******************************* 
c
c-------------------------------------------------------------------------
       subroutine setnp(np,xx,yy,zz)
       include "head.for"
c
c      外部変数
       dimension xx(ind),yy(ind), zz(ind)
       dimension rr(ind)
       dimension nb(ind)
       dimension ns(ind),ms(ind,499)
       dimension nt(ind),mt(ind,499)
c      ローカル変数
       dimension x0(ind),y0(ind),z0(ind),rr0(ind)
       dimension nb0(ind),ns0(ind)
       dimension md(ind)     !削除する節点
c
       data r0/0.02/   !半径の指標:
       data re/0.2/    !勾配の指標:平均半径の2倍
       data ve/0.1/    !節点の移動量:1/10
c
c       初期化
       ic=1;
c
       !境界の節点半径rrを決める
       !戻り値:rr(np)
       call setr(np,xx,yy,zz,nb,rr)
c
c      内部の節点:2d
       np=np+1;
       xx(np)=0.5;
       yy(np)=0.5;
       zz(np)=0.0;
       rr(np)=r0;
       nb(np)=0;
c
c       描画
       call plt_node(ic,np,xx,yy,zz,rr,nb,ind)

100    continue

       if(ic.eq.499)goto 200
       ic=ic+1
       write(*,*)"---------------------------------"
       write(*,*)"ic:",ic,"np:",np


       !周囲の節点:移動ベクトル用
       !戻り値:ns(np),ms(np,ns(np))
       call surround(np,xx,yy,zz,ns,ms,re,ind)


       !rrの補正:内部の節点
       !戻り値:rr(np)
       call revr(np,xx,yy,zz,nb,ns,ms,rr,ind)


       !周囲の節点:節点追加・削除用
       call surround2(np,xx,yy,zz,rr,nt,mt,ind)

         np0=0  !増やす節点
         nd=0  !削除する節点
         do i=1,np
         if(nb(i).eq.0)then

           !削除
c           if(nt(i).gt.3)then
           if(nt(i).ge.5)then
c           if(nt(i).ge.10)then
             d=100.0
             do j=1,nt(i)
               k=mt(i,j)
               dd=dsqrt((xx(i)-xx(k))**2
     &                 +(yy(i)-yy(k))**2
     &                 +(zz(i)-zz(k))**2)
               dd=dd/(rr(i)+rr(k))
               if(dd.lt.d)then
                 d=dd
               endif
             enddo
             if(d.le.0.6)then
               write(*,*)"del_0",i
               nd=nd+1
               md(nd)=i
             endif
           endif


           !追加
c           if(nt(i).lt.2)then
           if(nt(i).lt.5)then
c           if(nt(i).lt.10)then
             d=100.0
             do j=1,nt(i)
               k=mt(i,j)
               dd=dsqrt((xx(i)-xx(k))**2
     &                 +(yy(i)-yy(k))**2
     &                 +(zz(i)-zz(k))**2)
               dd=dd/(rr(i)+rr(k))
               if(dd.lt.d)d=dd
             enddo
             if(d.ge.0.8)then
               np0=np0+1
               x0(np0)=xx(i)-rr(i)/2.0  !2次元
               y0(np0)=yy(i)-rr(i)/2.0
               z0(np0)=zz(i)
               rr0(np0)=rr(i)
               write(*,*)"add_0:",i,nb(i),x0(np0),y0(np0),z0(np0)
             endif
           endif
       !    write(*,*)"ar:",i,r
         endif
         enddo


         !節点を増やす
         do kp=1,np0
           xx(np+kp)=x0(kp)
           yy(np+kp)=y0(kp)
           zz(np+kp)=z0(kp)
           rr(np+kp)=rr0(kp)
           nb(np+kp)=0
           ns(np+kp)=0
         enddo
         np=np+np0

         !節点削除
         np0=np
         do kp=1,np
           x0(kp)=xx(kp)
           y0(kp)=yy(kp)
           z0(kp)=zz(kp)
           rr0(kp)=rr(kp)
           nb0(kp)=nb(kp)
           ns0(kp)=ns(kp)
         enddo

         if(nd.gt.0)nd=1    !削除は1つずつ
         np=0
         do kp=1,np0
           k=0
           do j=1,nd
             if(kp.eq.md(j))k=1 !削除する節点はとばす
           enddo
           if(k.eq.0)then
             np=np+1
             xx(np)=x0(kp)
             yy(np)=y0(kp)
             zz(np)=z0(kp)
             rr(np)=rr0(kp)
             nb(np)=nb0(kp)
             ns(np)=ns0(kp)
           endif
         enddo


         !更新/////////////////////////////////

         !周囲の節点
         call surround(np,xx,yy,zz,ns,ms,re,ind)

         !rrの補正
         !戻り値:rr(np)
         call revr(np,xx,yy,zz,nb,ns,ms,rr,ind)

         !/////////////////////////////////////


         !描画
         call plt_node(ic,np,xx,yy,zz,rr,nb,ind)


         !節点の移動
         do kp=1,np

         !内部の節点を移動する
         if(nb(kp).eq.0)then

           write(*,*)"xx:",kp,sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
     &                        ,sngl(rr(kp))

           !節点の移動方向を決める
           call cnavra(kp,xx,yy,zz,rr,ns,ms,vx,vy,vz,ve,nb,ind)

           xx(kp)=xx(kp)+vx
           yy(kp)=yy(kp)+vy
           zz(kp)=zz(kp)+vz
         endif
         enddo

         !状態評価

         goto 100

200      continue
c         
         return
         end
c-------------------------------------------------------------------------
c
c       節点の移動ベクトルを計算する
c
c       kp    :節点
c       xx(np):x座標
c       yy(np):y座標
c       zz(np):z座標
c       rr(np):節点半径
c       ns(np):周囲の節点数
c       ms(np,ns(np)):周囲の節点番号
c       vx,vy,vz     :移動ベクトル
c       ve           :節点の基準移動量
c       nb(np)        0:内部の節点,1:境界の節点
c       ind          :節点数の上限
c
c-------------------------------------------------------------------------
       subroutine cnavra(kp,xx,yy,zz,rr,ns,ms,vx,vy,vz,ve,nb,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind),rr(ind)
       dimension ns(ind),ms(ind,499)
       dimension nb(ind)

c         int kp
c         double x,y,z
c         double cx,cy,cz
c         double dx,dy,dz
c         double cc
c         double dd
c         int i,k

         x=xx(kp)
         y=yy(kp)
         z=zz(kp)

         cx=0.0
         cy=0.0
         cz=0.0

         !方向ベクトル
         do i=1,ns(kp)
           k=ms(kp,i)

           dx=x-xx(k)  !斥力方向
           dy=y-yy(k)
           dz=z-zz(k)

c           write(*,*)"dx:",kp,k,dx,dy,dz

           dd=dsqrt(dx**2+dy**2+dz**2)

           if(dd.le.0.000001)then
             dx=1.0
             dy=0.0
             dz=0.0
           else
             dx=dx/dd
             dy=dy/dd
             dz=dz/dd
           endif

           cc=dd/(rr(kp)+rr(k))

c           write(*,*)"cc:",kp,k,dd,rr(kp)+rr(k),cc,w(cc)

           cx=cx+w(cc,nb(k))*dx
           cy=cy+w(cc,nb(k))*dy
           cz=cz+w(cc,nb(k))*dz
         enddo

         cc=dsqrt(cx**2+cy**2+cz**2)

         if(cc.gt.0.0)then
           vx=ve*rr(kp)*cx/cc
           vy=ve*rr(kp)*cy/cc
           vz=ve*rr(kp)*cz/cc
         else
           vx=0.0
           vy=0.0
           vz=0.0
         endif

c         write(*,*)"vx:",vx,vy,vz,cc
c         pause

       return
       end
c-----------------------------------------------------------------------
c
c      重み関数
c        境界節点の斥力は大きくする
c
c-----------------------------------------------------------------------
       function w(cc,nb)
       implicit double precision(a-h,o-z)

         !斥力
         if(cc.le.1.0)then
       !   w=-log10(cc+0.0001)
           if(nb.eq.0)w=-100.0*cc+100.0     !内部の節点
           if(nb.eq.1)w=-1000.0*cc+1000.0   !境界の節点

         !引力
         elseif(cc.ge.1.0)then
       !   w=-log10(cc+0.0001)/10.0
           w=-0.01*cc+0.01
         endif

       end
c-----------------------------------------------------------------------
c
c      境界の節点半径rrを決める
c
c-----------------------------------------------------------------------
       subroutine setr(np,xx,yy,zz,nb,rr)
       include "head.for"
       dimension xx(ind),yy(ind),zz(ind),rr(ind)
       dimension nb(ind)
c
       do 100 i=1,np
         dmin=1000d0
         do 200 j=1,np
           if(i.eq.j)goto 200
           dd=dsqrt((xx(i)-xx(j))**2+(yy(i)-yy(j))**2+(zz(i)-zz(j))**2)
           if(dd.lt.dmin)dmin=dd
200      continue
         rr(i)=dmin/2d0
         nb(i)=1
c         write(*,*)"rr:",i,rr(i)
c         pause
100    continue
c
       return
       end
c-----------------------------------------------------------------------
c
c      周囲の節点:節点半径rrの補正用
c
c-----------------------------------------------------------------------
       subroutine surround(np,xx,yy,zz,ns,ms,re,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind)
       dimension ns(ind),ms(ind,499)
c
       do i=1,np
           ns(i)=0
           do j=1,np
           if(i.ne.j)then
             dd=dsqrt((xx(i)-xx(j))**2
     &               +(yy(i)-yy(j))**2
     &               +(zz(i)-zz(j))**2)
             if(dd.le.re)then
               ns(i)=ns(i)+1  !周囲の節点の数
               ms(i,ns(i))=j  !周囲の節点の番号
             endif
           endif
         enddo
       enddo
c
       return
       end
c-----------------------------------------------------------------------
c
c      周囲の節点:節点追加・削除用
c
c-----------------------------------------------------------------------
       subroutine surround2(np,xx,yy,zz,rr,nt,mt,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind),rr(ind)
       dimension nt(ind),mt(ind,499)
c
       do i=1,np
         nt(i)=0
         do j=1,np
         if(i.ne.j)then
           dd=dsqrt((xx(i)-xx(j))**2+(yy(i)-yy(j))**2+(zz(i)-zz(j))**2)
           if(dd.le.rr(i)*2.1)then
             nt(i)=nt(i)+1  !周囲の節点の数
             mt(i,nt(i))=j  !周囲の節点の番号
           endif
         endif
         enddo
       enddo
c
       return
       end
c-----------------------------------------------------------------------
c
c      内部節点のrr補正
c
c-----------------------------------------------------------------------
       subroutine revr(np,xx,yy,zz,nb,ns,ms,rr,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind),rr(ind)
       dimension ns(ind),ms(ind,499)
       dimension nb(ind)
c
       do i=1,np
       if(nb(i).eq.0 .and. ns(i).gt.0)then
         r = 0.0
         d = 0.0
         do j=1,ns(i)
           k=ms(i,j)
           dd=dsqrt((xx(i)-xx(k))**2+(yy(i)-yy(k))**2+(zz(i)-zz(k))**2)
           dd=1.0/(dd+0.001)
           d=d+dd
           r=r+rr(k)*dd
       !    write(*,*)"rr:",i,j,dd,rr(k)
         enddo
         dd=1.0/(0.0+0.001)
         d=d+dd
         r=r+rr(i)*dd
         rr(i)=r/d
       !  write(*,*)"r0:",i,r,d,r/d,rr(i)
       endif
       enddo
c
       return
       end
