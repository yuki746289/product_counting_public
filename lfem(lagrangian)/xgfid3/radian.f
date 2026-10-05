c----------------------------------------------------------------------
c
c       曲率を計算する
c
c----------------------------------------------------------------------
       subroutine calrad(np,nele,ne,xx,yy,zz,nb,nelb,neb,hh,vnx,vny,vnz)
        include "header.h"
        dimension ne(ine,4),nec(ine,4),mk(ind)
        dimension nelb(ine),neb(ine),melb(ine),meb(ine)
        dimension xt(ind),yt(ind),zt(ind)
        dimension xx(ind),yy(ind),zz(ind)
        dimension hh(ind),vnx(ind),vny(ind),vnz(ind)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension jm(4,3)
        dimension mm(ind)
c-------表示用-------------------------------------------------
        dimension mc(ine,3)
c--------------------------------------------------------------
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
        data PI/3.1415926535897932384626433832795/
c
c        write(*,*)"np in calrad"
c        call pltrset()
c        call pltsetnp(np,xx,yy,zz,ind)
c        call pltsetne(nele,ne,2d0,2,ine)
c        pause
c
        do j=1,3
          jm(1,j)=jm1(j)
          jm(2,j)=jm2(j)
          jm(3,j)=jm3(j)
          jm(4,j)=jm4(j)
        enddo
c
c       初期化
        do 600 kp=1,np
          hh(kp)=0d0
          vnx(kp)=0d0
          vny(kp)=0d0
          vnz(kp)=0d0
600     continue
c
c       隣接する要素
        do 10 kele=1,nele
        do 10 j=1,4
          nec(kele,j)=0
10      continue
c
        do 100 i=1,nele-1
        do 200 j=1,4
          if(nec(i,j).ne.0)goto 200
          m1=ne(i,jm(j,1))
          m2=ne(i,jm(j,2))
          m3=ne(i,jm(j,3))
          do k=i+1,nele
          do l=1,4
            n1=ne(k,jm(l,1))
            n2=ne(k,jm(l,2))
            n3=ne(k,jm(l,3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m2.eq.n3 .and. m3.eq.n2 .and. m1.eq.n1) .or.
     &         (m3.eq.n3 .and. m1.eq.n2 .and. m2.eq.n1))then
              nec(i,j)=k
              nec(k,l)=i
              goto 200
            endif
          enddo
          enddo
200     continue
100     continue
c
c       表面の面と節点の数
        nb=0    !表面の面の数
        mp=0    !表面の節点の数
        do 300 i=1,nele
        do 300 j=1,4
          if(nec(i,j).ne.0)goto 300
          nb=nb+1       !表面の数
          nelb(nb)=i    !表面の要素の番号
          neb(nb)=j     !表面の面の番号
c
          do 500 k=1,3
            do kp=1,mp
              if(ne(i,jm(j,k)).eq.mm(kp))goto 500
            enddo
            mp=mp+1
            mm(mp)=ne(i,jm(j,k))
500       continue
300     continue
c        write(*,*)"nb",nb
c        write(*,*)"mp",mp
c
c       曲率を計算する
c        do 1000 kp=1,np
        do 1000 km=1,mp
          kp=mm(km)
c          write(*,*)"kp",kp,sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c          pause
c
c         kpを含む表面
          mb=0      !表面の数
          do kb=1,nb
            m1=ne(nelb(kb),jm(neb(kb),1))
            m2=ne(nelb(kb),jm(neb(kb),2))
            m3=ne(nelb(kb),jm(neb(kb),3))
            if(kp.eq.m1 .or. kp.eq. m2 .or. kp.eq.m3)then
              mb=mb+1
              melb(mb)=nelb(kb)
              meb(mb)=neb(kb)
              !表示用
              mc(mb,1)=m1
              mc(mb,2)=m2
              mc(mb,3)=m3
            endif
          enddo
          if(mb.eq.0)goto 1000  !表面が無いときは飛ばす
c          call pltsetme(mb,mc,3d0,2,ind,ine)
c          write(*,*)"mb",mb
c
c         kpの周囲の節点を探す
          nk=1  !kpの周囲の節点(kp含む)
          mk(nk)=kp
          xt(nk)=xx(kp)
          yt(nk)=yy(kp)
          zt(nk)=zz(kp)
          do 400 kb=1,mb
          do 400 j=1,3
            m=ne(melb(kb),jm(meb(kb),j))
            do kk=1,nk
              if(m.eq.mk(kk))goto 400
            enddo
            nk=nk+1
            mk(nk)=m
            xt(nk)=xx(m)
            yt(nk)=yy(m)
            zt(nk)=zz(m)
400       continue
c          call nssurf(kp,nels,ms,nk,mk)
c          write(*,*)"nk",nk
c          do i=1,nk
c            write(*,*)"nk",i,sngl(xt(i)),sngl(yt(i)),sngl(zt(i))
c          enddo
c
c         法線ベクトルの和
          dnx=0d0
          dny=0d0
          dnz=0d0
          do kb=1,mb
            m1=ne(melb(kb),jm(meb(kb),1))
            m2=ne(melb(kb),jm(meb(kb),2))
            m3=ne(melb(kb),jm(meb(kb),3))
            v1x=xx(m2)-xx(m1)
            v1y=yy(m2)-yy(m1)
            v1z=zz(m2)-zz(m1)
            v2x=xx(m3)-xx(m1)
            v2y=yy(m3)-yy(m1)
            v2z=zz(m3)-zz(m1)
            call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx0,dny0,dnz0)
            dnx=dnx+dnx0    !x軸方向の方向余弦
            dny=dny+dny0
            dnz=dnz+dnz0
          enddo
          dn=dsqrt(dnx**2+dny**2+dnz**2)
          dnx=dnx/dn
          dny=dny/dn
          dnz=dnz/dn
c          write(*,*)"dn",sngl(dnx),sngl(dny),sngl(dnz) !座標変換して角度を算出する
c
c         曲率の計算:z=ax^2+y^2
           call calcsphere(nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
!          call calcx2y2(kp,nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
c          call calcx2y22(nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
c          call calcx2y23(nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
c          call curvature(nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
c          call curvature2(nk,xt,yt,zt,dnx1,dny1,dnz1,ee,h)
c
c         法線ベクトルの向きを決める
          h=-h
          shita=0d0
          shita=finpr(dnx,dny,dnz,dnx1,dny1,dnz1)
          shita=PI-shita
          if(shita.gt.90d0/180d0*PI)then h=-h    !曲率の面の向き
c          dn=dsqrt(dnx1**2+dny1**2+dnz1**2)
c          write(*,*)"shita:",sngl(shita/PI*180d0)
c          write(*,*)"dn:",sngl(dnx1),sngl(dny1),sngl(dnz1),sngl(dn)
c          write(*,*)"h:",h
c          call pltnt(xx(kp)+dnx1,yy(kp)+dny1,zz(kp)+dnz1,5d0,2)
c
c         曲率を登録する
          hh(kp)=h
          vnx(kp)=dnx1
          vny(kp)=dny1
          vnz(kp)=dnz1
c          write(*,*)"h,ee",h,ee
c          pause
1000    continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       節点を法線方向に座標変換
c
c---------------------------------------------------------------------------
        subroutine transcord(nk,mk,xx,yy,zz,dnx,dny,dnz,xk,yk,zk,ind)
        implicit double precision(a-h,o-z)
        dimension mk(ind),xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
c
        x0=xx(mk(1))
        y0=yy(mk(1))
        z0=zz(mk(1))
        s=(dnx**2+dnx**2)*(dnx**2+dnz**2)
     &   +(dnz-dny)*dnx*dny*dsqrt(dnx**2+dnz**2)+(dnx+dnz)*dny**2*dnz
        do kk=1,nk
          x=xx(mk(kk))
          y=yy(mk(kk))
          z=zz(mk(kk))
          xk(kk)=(dnz*x-dnx*z)*(dnx**2+dnz**2)
     &          +(dny*x-dnx*y)*dnx*dsqrt(dnx**2+dnz**2)
     &          +(dnz*y-dny*z)*dny*dnz
          yk(kk)=(z-x)*dnx*dny-(x+z)*dny*dnz+dnx**2*y+dnz**2*y
          zk(kk)=(dnx*x+dnz*z)*(dnx**2+dnz**2)
     &          +(dnz*y-dny*x)*dnx*dsqrt(dnx**2+dnz**2)
     &          +(dnx*y+dny*z)*dny*dnz
          xk(kk)=xk(kk)/s-x0
          yk(kk)=yk(kk)/s-y0
          zk(kk)=zk(kk)/s-z0
        enddo
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:z=ax^2+y^2
c
c---------------------------------------------------------------------------
        subroutine calcx2y2(kp,nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        dimension mk(ind)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
        common /radian/phi2(ind),psi2(ind),eta2(ind)
        data PI/3.1415926535897932384626433832795/
c
c       前回のタイムステップの角度を参照
        if(phi2(kp).ne.-999d0 .and. psi2(kp).ne.-999d0
     &                        .and. eta2(kp).ne.-999d0)then
          phi0=phi2(kp)
          psi0=psi2(kp)
          eta0=eta2(kp)
          write(*,*)"kp:",kp,sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
     &                       ,sngl(eta0/PI*180d0)
          pause
c          goto 1000
        end if
c
c       5°ずつ探す
        e0=1d10
        emax0=1d10
        do 200 i=0,180,5
        do 200 j=0,180,5
        do 200 k=0,90,5
          phi=dble(i)/180d0*PI
          psi=dble(j)/180d0*PI
          eta=dble(k)/180d0*PI
c
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         ∑を計算する
          sgx=0d0       !∑x
          sgy=0d0       !∑y
          sgz=0d0       !∑z
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx2y2=0d0    !∑x^2y^2
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          do kk=1,nk
            sgx=sgx+xk(kk)  !∑x
            sgy=sgy+yk(kk)  !∑y
            sgz=sgz+zk(kk)  !∑z
            sgx4=sgx4+xk(kk)**4               !∑x^4
            sgy4=sgy4+yk(kk)**4               !∑y^4
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
            sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
          enddo
c
c         z=ax^2+by^2の係数a,bを計算する
          a1=sgy4*sgx2z-sgx2y2*sgy2z
          a2=sgy4*sgx4-sgx2y2*sgx2y2
          a=0d0
          if(dabs(a2).gt.1d-10)a=a1/a2
          b1=sgx2y2*sgx2z-sgx4*sgy2z
          b2=sgx2y2*sgx2y2-sgx4*sgy4
          b=0d0
          if(dabs(a2).gt.1d-10)b=b1/b2
c          if(a*b.lt.0d0)goto 200
c
c         誤差を計算する
          ee=0d0
          emax=0d0
          do kk=1,nk
            e1=dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
            ee=ee+e1
            if(emax.lt.e1)emax=e1
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
c          if(ee.lt.e0)then
          if(emax.lt.emax0)then
            a0=a
            b0=b
            e0=ee
            emax0=emax
            ii=i
            jj=j
            kk=k
            phi0=phi
            psi0=psi
            eta0=eta
c
            x0=xk(1)
            y0=yk(1)
            z0=zk(1)
c
            p=2d0*a*x0    !∂z/∂x
            q=2d0*b*y0    !∂z/∂y
            r=2d0*a       !∂2z/∂x2
            s=0d0         !∂2z/∂x∂y
            t=2d0*b       !∂z/∂y2
c           平均曲率:H
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2   !平均曲率Η
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
c
c            do kk=1,nk
c              write(*,*)"x:",kk,sngl(xk(kk)),sngl(yk(kk)),sngl(zk(kk))
c            enddo
c        write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(e0),sngl(emax0)
c        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
c     &                                              ,sngl(eta0/PI*180d0)
c        e=dsqrt(ex**2+ey**2+ez**2)
c          write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
c          write(*,*)"dn:",sngl(dnx),sngl(dny),sngl(dnz)
          endif
200     continue
c
c       前回のタイムステップの角度を参照
1000    continue
c
c       1°ずつ探す
        phi1=phi0
        psi1=psi0
        eta1=eta0
        do 100 i=-5,5
        do 100 j=-5,5
        do 100 k=-5,5
          phi=phi1-dble(i)/180d0*PI
          psi=psi1-dble(j)/180d0*PI
          eta=eta1-dble(k)/180d0*PI
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         ∑を計算する
          sgx=0d0       !∑x
          sgy=0d0       !∑y
          sgz=0d0       !∑z
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx2y2=0d0    !∑x^2y^2
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          do kk=1,nk
            sgx=sgx+xk(kk)  !∑x
            sgy=sgy+yk(kk)  !∑y
            sgz=sgz+zk(kk)  !∑z
            sgx4=sgx4+xk(kk)**4               !∑x^4
            sgy4=sgy4+yk(kk)**4               !∑y^4
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
            sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
          enddo
c
c         z=ax^2+by^2の係数a,bを計算する
          a1=sgy4*sgx2z-sgx2y2*sgy2z
          a2=sgy4*sgx4-sgx2y2*sgx2y2
          a=0d0
          if(dabs(a2).gt.1d-10)a=a1/a2
          b1=sgx2y2*sgx2z-sgx4*sgy2z
          b2=sgx2y2*sgx2y2-sgx4*sgy4
          b=0d0
          if(dabs(a2).gt.1d-10)b=b1/b2

c          if(a*b.lt.0d0)goto 200
c
c         誤差を計算する
          ee=0d0
          emax=0d0
          do kk=1,nk
            e1=dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
            ee=ee+e1
            if(emax.lt.e1)emax=e1
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
c          if(ee.lt.e0)then
          if(emax.lt.emax0)then
            a0=a
            b0=b
            e0=ee
            emax0=emax
            phi0=phi
            psi0=psi
            eta0=eta
c
            x0=xk(1)
            y0=yk(1)
            z0=zk(1)
c
            p=2d0*a*x0    !∂z/∂x
            q=2d0*b*y0    !∂z/∂y
            r=2d0*a       !∂2z/∂x2
            s=0d0         !∂2z/∂x∂y
            t=2d0*b       !∂z/∂y2
c           平均曲率:H
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2   !平均曲率Η
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
          endif
100     continue
c
c       0.1°ずつ探す
        phi1=phi0
        psi1=psi0
        eta1=eta0
        do 300 i=-5,5
        do 300 j=-5,5
        do 300 k=-5,5
          phi=phi1-dble(i)*1d-1/180d0*PI
          psi=psi1-dble(j)*1d-1/180d0*PI
          eta=eta1-dble(k)*1d-1/180d0*PI
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         ∑を計算する
          sgx=0d0       !∑x
          sgy=0d0       !∑y
          sgz=0d0       !∑z
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx2y2=0d0    !∑x^2y^2
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          do kk=1,nk
            sgx=sgx+xk(kk)  !∑x
            sgy=sgy+yk(kk)  !∑y
            sgz=sgz+zk(kk)  !∑z
            sgx4=sgx4+xk(kk)**4               !∑x^4
            sgy4=sgy4+yk(kk)**4               !∑y^4
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
            sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
          enddo
c
c         z=ax^2+by^2の係数a,bを計算する
          a1=sgy4*sgx2z-sgx2y2*sgy2z
          a2=sgy4*sgx4-sgx2y2*sgx2y2
          a=0d0
          if(dabs(a2).gt.1d-10)a=a1/a2
          b1=sgx2y2*sgx2z-sgx4*sgy2z
          b2=sgx2y2*sgx2y2-sgx4*sgy4
          b=0d0
          if(dabs(a2).gt.1d-10)b=b1/b2

c          if(a*b.lt.0d0)goto 200
c
c         誤差を計算する
          ee=0d0
          emax=0d0
          do kk=1,nk
            e1=dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
            ee=ee+e1
            if(emax.lt.e1)emax=e1
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
c          if(ee.lt.e0)then
          if(emax.lt.emax0)then
            a0=a
            b0=b
            e0=ee
            emax0=emax
            phi0=phi
            psi0=psi
            eta0=eta
c
            x0=xk(1)
            y0=yk(1)
            z0=zk(1)
c
            p=2d0*a*x0    !∂z/∂x
            q=2d0*b*y0    !∂z/∂y
            r=2d0*a       !∂2z/∂x2
            s=0d0         !∂2z/∂x∂y
            t=2d0*b       !∂z/∂y2
c           平均曲率:H
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2   !平均曲率Η
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
          endif
300     continue
c
        phi2(kp)=phi0
        psi2(kp)=psi0
        eta2(kp)=eta0
c
c        write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(e0),sngl(emax0)
c        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
c     &                                              ,sngl(eta0/PI*180d0)
c        e=dsqrt(ex**2+ey**2+ez**2)
c        write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:z=ax^2+y^2
c
c---------------------------------------------------------------------------
        subroutine calcx2y22(nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        dimension mk(ind)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
        data PI/3.1415926535897932384626433832795/
c
        e0=1d10
        emax0=1d10
        phi0=0d0
        psi0=0d0
        eta0=0d0
c        phi0=dble(320)/180d0*PI
c        psi0=dble(220)/180d0*PI
c        eta0=dble(180)/180d0*PI

        do 200 i=1,1000
          phi=phi0
          psi=psi0
          eta=eta0
          do 300 j=1,6
            if(j.eq.1)phi=phi0+1d-2/180d0*PI
            if(j.eq.2)phi=phi0-1d-2/180d0*PI
            if(j.eq.3)psi=psi0+1d-2/180d0*PI
            if(j.eq.4)psi=psi0-1d-2/180d0*PI
            if(j.eq.5)eta=eta0+1d-2/180d0*PI
            if(j.eq.6)eta=eta0-1d-2/180d0*PI
            do l=1,nk
              x0=xx(l)-xx(1)
              y0=yy(l)-yy(1)
              z0=zz(l)-zz(1)
              call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
              xk(l)=x1
              yk(l)=y1
              zk(l)=z1
            enddo
c
c           ∑を計算する
            sgx=0d0       !∑x
            sgy=0d0       !∑y
            sgz=0d0       !∑z
            sgx4=0d0      !∑x^4
            sgy4=0d0      !∑y^4
            sgx2y2=0d0    !∑x^2y^2
            sgx2z=0d0     !∑x^2z
            sgy2z=0d0     !∑y^2z
            do kk=1,nk
              sgx=sgx+xk(kk)  !∑x
              sgy=sgy+yk(kk)  !∑y
              sgz=sgz+zk(kk)  !∑z
              sgx4=sgx4+xk(kk)**4               !∑x^4
              sgy4=sgy4+yk(kk)**4               !∑y^4
              sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
              sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
              sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
            enddo
c
c           z=ax^2+by^2の係数a,bを計算する
            a1=sgy4*sgx2z-sgx2y2*sgy2z
            a2=sgy4*sgx4-sgx2y2*sgx2y2
            a=a1/a2
            b1=sgx2y2*sgx2z-sgx4*sgy2z
            b2=sgx2y2*sgx2y2-sgx4*sgy4
            b=b1/b2
c            if(a*b.lt.0d0)goto 200
c
c           誤差を計算する
            ee=0d0
            emax=0d0
            do kk=1,nk
              e1=dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
              ee=ee+e1
              if(emax.lt.e1)emax=e1
            enddo
            if(emax.lt.emax0)goto 400
c            if(ee.lt.e0)goto 400
300       continue
          exit !ループを抜ける
c
c         平均曲率Η=(κ1+κ2)/2を計算する
400       continue
          a0=a
          b0=b
          e0=ee
          emax0=emax
          phi0=phi
          psi0=psi
          eta0=eta
c
          x0=xk(1)
          y0=yk(1)
          z0=zk(1)
c
          p=2d0*a*x0    !∂z/∂x
          q=2d0*b*y0    !∂z/∂y
          r=2d0*a       !∂2z/∂x2
          s=0d0         !∂2z/∂x∂y
          t=2d0*b       !∂z/∂y2
c         平均曲率:H
          h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
          h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
          h=h1/h2   !平均曲率Η
c         Gaussの曲率:K
          gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c         単位法ベクトルの方向余弦:ex,ey,ez →外向き
          e=dsqrt(1d0+p**2+q**2)
          ex=-p/e
          ey=-q/e
          ez=1d0/e
          call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
c
          do kk=1,nk
            write(*,*)"x:",kk,sngl(xk(kk)),sngl(yk(kk)),sngl(zk(kk))
          enddo
          write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(e0),sngl(emax0)
          write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0)
     &         ,sngl(psi0/PI*180d0),sngl(eta0/PI*180d0)
          e=dsqrt(ex**2+ey**2+ez**2)
          write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
          write(*,*)"dn:",sngl(dnx),sngl(dny),sngl(dnz)
200     continue
c
        write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(e0),sngl(emax0)
        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
     &                                              ,sngl(eta0/PI*180d0)
        e=dsqrt(ex**2+ey**2+ez**2)
        write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:z=ax^2+y^2
c
c---------------------------------------------------------------------------
        subroutine calcx2y23(nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        dimension mk(ind)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
        data PI/3.1415926535897932384626433832795/
c
c       基準のベクトル
        px=0d0
        py=0d0
        pz=0d0
        do l=2,nk
          vx=xx(l)-xx(1)
          vy=yy(l)-yy(1)
          vz=zz(l)-zz(1)
          vv=dsqrt(vx**2+vy**2+vz**2)
          vx=vx/vv
          vy=vy/vv
          vz=vz/vv
          px=px+vx
          py=py+vy
          pz=pz+vz
        enddo
        phi1=finpr(px,py,pz,1d0,0d0,0d0) !x軸との角度
        psi1=finpr(px,py,pz,0d0,1d0,0d0) !y軸との角度
        eta1=finpr(px,py,pz,0d0,0d0,1d0) !z軸との角度
        write(*,*)"phi,psi,eta:",sngl(phi1/PI*180d0),sngl(psi1/PI*180d0)
     &                                              ,sngl(eta1/PI*180d0)
c
c       1°ずつ探す
        e0=1d10
        emax0=1d10
        do 200 i=0,180,3
        do 200 j=0,180,3
        do 200 k=0,90,3
c          phi=phi1+dble(i)/180d0*PI
c          psi=psi1+dble(j)/180d0*PI
c          eta=eta1+dble(k)/180d0*PI

          phi=dble(i)/180d0*PI
          psi=dble(j)/180d0*PI
          eta=dble(k)/180d0*PI
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         ∑を計算する
          sgx=0d0       !∑x
          sgy=0d0       !∑y
          sgz=0d0       !∑z
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx2y2=0d0    !∑x^2y^2
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          do kk=1,nk
            sgx=sgx+xk(kk)  !∑x
            sgy=sgy+yk(kk)  !∑y
            sgz=sgz+zk(kk)  !∑z
            sgx4=sgx4+xk(kk)**4               !∑x^4
            sgy4=sgy4+yk(kk)**4               !∑y^4
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2 !∑x^2y^2
            sgx2z=sgx2z+xk(kk)**2*zk(kk)      !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)      !∑y^2z
          enddo
c
c         z=ax^2+by^2の係数a,bを計算する
          a1=sgy4*sgx2z-sgx2y2*sgy2z
          a2=sgy4*sgx4-sgx2y2*sgx2y2
          a=0d0
          if(dabs(a2).gt.1d-10)a=a1/a2
          b1=sgx2y2*sgx2z-sgx4*sgy2z
          b2=sgx2y2*sgx2y2-sgx4*sgy4
          b=0d0
          if(dabs(a2).gt.1d-10)b=b1/b2
c          if(a*b.lt.0d0)goto 200
c
c         誤差を計算する
          ee=0d0
          emax=0d0
          do kk=1,nk
            e1=dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2)
            ee=ee+e1
            if(emax.lt.e1)emax=e1
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
c          if(ee.lt.e0)then
          if(emax.lt.emax0)then
            a0=a
            b0=b
            e0=ee
            emax0=emax
            phi0=phi
            psi0=psi
            eta0=eta
c
            x0=xk(1)
            y0=yk(1)
            z0=zk(1)
c
            p=2d0*a*x0    !∂z/∂x
            q=2d0*b*y0    !∂z/∂y
            r=2d0*a       !∂2z/∂x2
            s=0d0         !∂2z/∂x∂y
            t=2d0*b       !∂z/∂y2
c           平均曲率:H
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2   !平均曲率Η
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
c
c            do kk=1,nk
c              write(*,*)"x:",kk,sngl(xk(kk)),sngl(yk(kk)),sngl(zk(kk))
c            enddo
c        write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(emax0)
c        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
c     &                                              ,sngl(eta0/PI*180d0)
c          e=dsqrt(ex**2+ey**2+ez**2)
c          write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
c          write(*,*)"dn:",sngl(dnx),sngl(dny),sngl(dnz)
c            call pltrset()
c            call pltnt(0d0,1d0,0d0,5d0,2)  !赤
c            call pltsetnp(nk,xk,zk,yk,ind)
          endif
200     continue
c
        write(*,*)"a,b,ee:",sngl(a0),sngl(b0),sngl(emax0)
        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
     &                                              ,sngl(eta0/PI*180d0)
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:z=ax^2+by^2+cxy+dx+ey
c
c---------------------------------------------------------------------------
        subroutine curvature(nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        parameter(n=5,m=6)
        dimension sa(n,m),sf(n)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
        data PI/3.1415926535897932384626433832795/
c
c       10°ずつ探す
        e0=1d10
        do 200 i=-9,9
        do 200 j=0,35
        do 200 k=0,35
          phi=dble(i*10)/180d0*PI
          psi=dble(j*10)/180d0*PI
          eta=dble(k*10)/180d0*PI
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         z=ax^2+by^2+cxy+dx+eyの係数を計算する
          call calccoefficient(nk,xk,yk,zk,a,b,c,d,e,ill)
          if(ill.eq.1)goto 200
c          write(*,*)"ct:",sngl(a),sngl(b),sngl(c),sngl(d),sngl(e)
c
c         誤差を計算する:z=ax^2+by^2+cxy+dx+ey
          ee=0d0
          do kk=1,nk
            ee=ee+dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2
     &                -c*xk(kk)*yk(kk)-d*xk(kk)-e*yk(kk))
          enddo
c
c         平均曲率Η=(κ1+κ2)/2を計算する
          if(ee.lt.e0)then
            e0=ee
            ii=i
            jj=j
            kk=k
            phi0=(ii*10)/180d0*PI
            psi0=(jj*10)/180d0*PI
            eta0=(kk*10)/180d0*PI
            x0=xx(1)
            y0=yy(1)
            z0=zz(1)
c
            p=2d0*a*x0+c*y0+d   !∂z/∂x
            q=2d0*b*y0+c*x0+e   !∂z/∂y
            r=2d0*a             !∂2z/∂x2
            s=c                 !∂2z/∂x∂y
            t=2d0*b             !∂z/∂y2
c           平均曲率:H
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
          endif
200     continue
c
c       1°ずつ探す
        do 100 i=-10,10
        do 100 j=-10,10
        do 100 k=-10,10
          phi=(dble(ii)*10+dble(i))/180d0*PI
          psi=(dble(jj)*10+dble(j))/180d0*PI
          eta=(dble(kk)*10+dble(k))/180d0*PI
          do l=1,nk
            x0=xx(l)-xx(1)
            y0=yy(l)-yy(1)
            z0=zz(l)-zz(1)
            call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
            xk(l)=x1
            yk(l)=y1
            zk(l)=z1
          enddo
c
c         z=ax^2+by^2+cxy+dx+eyの係数を計算する
          call calccoefficient(nk,xk,yk,zk,a,b,c,d,e,ill)
          if(ill.eq.1)goto 100
c          write(*,*)"ct:",sngl(a),sngl(b),sngl(c),sngl(d),sngl(e)
c
c         誤差を計算する:z=ax^2+by^2+cxy+dx+ey
          ee=0d0
          do kk=1,nk
            ee=ee+dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2
     &                -c*xk(kk)*yk(kk)-d*xk(kk)-e*yk(kk))
          enddo
c
c         z=ax^2+by^2の平均曲率Η=(κ1+κ2)/2を計算する
          if(ee.lt.e0)then
            e0=ee
            phi0=phi
            psi0=psi
            eta0=eta
            x0=xx(1)
            y0=yy(1)
            z0=zz(1)
c
            p=2d0*a*x0+c*y0+d   !∂z/∂x
            q=2d0*b*y0+c*x0+e   !∂z/∂y
            r=2d0*a             !∂2z/∂x2
            s=c                 !∂2z/∂x∂y
            t=2d0*b             !∂z/∂y2
            h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
            h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
            h=h1/h2   !平均曲率Η
c           Gaussの曲率:K
            gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c           単位法ベクトルの方向余弦:ex,ey,ez →外向き
            e=dsqrt(1d0+p**2+q**2)
            ex=-p/e
            ey=-q/e
            ez=1d0/e
            call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
          endif
100     continue
c
        write(*,*)"ee:",sngl(e0)
        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
     &                                              ,sngl(eta0/PI*180d0)
c        write(*,*)"h:",h
        e=dsqrt(ex**2+ey**2+ez**2)
        write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
        do k=1,nk
          x0=xx(k)-xx(1)
          y0=yy(k)-yy(1)
          z0=zz(k)-zz(1)
          call rotxyz(x0,y0,z0,x1,y1,z1,phi0,psi0,eta0)
          xk(k)=x1
          yk(k)=y1
          zk(k)=z1
        enddo
c        call pltrset()
c        call pltnt(0d0,0d0,1d0,5d0,2)  !赤
c        call pltsetnp(nk,xk,yk,zk,ind)
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:z=ax^2+by^2+cxy+dx+ey
c
c---------------------------------------------------------------------------
        subroutine curvature2(nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        parameter(n=5,m=6)
        dimension sa(n,m),sf(n)
        dimension xx(ind),yy(ind),zz(ind)
        dimension xk(ind),yk(ind),zk(ind)
        data PI/3.1415926535897932384626433832795/
c
        e0=1d10
        phi0=0d0
        psi0=0d0
        eta0=0d0
        do 200 i=1,1000
          phi=phi0
          psi=psi0
          eta=eta0
          do 300 j=1,6
            if(j.eq.1)phi=phi0+1d0/180d0*PI
            if(j.eq.2)phi=phi0-1d0/180d0*PI
            if(j.eq.3)psi=psi0+1d0/180d0*PI
            if(j.eq.4)psi=psi0-1d0/180d0*PI
            if(j.eq.5)eta=eta0+1d0/180d0*PI
            if(j.eq.6)eta=eta0-1d0/180d0*PI
            do l=1,nk
              x0=xx(l)-xx(1)
              y0=yy(l)-yy(1)
              z0=zz(l)-zz(1)
              call rotxyz(x0,y0,z0,x1,y1,z1,phi,psi,eta)
              xk(l)=x1
              yk(l)=y1
              zk(l)=z1
            enddo
c
c           z=ax^2+by^2+cxy+dx+eyの係数を計算する
            call calccoefficient(nk,xk,yk,zk,a,b,c,d,e,ill)
            if(ill.eq.1)goto 300
c            write(*,*)"ct:",sngl(a),sngl(b),sngl(c),sngl(d),sngl(e)
c
c           誤差を計算する:z=ax^2+by^2+cxy+dx+ey
            ee=0d0
            do kk=1,nk
              ee=ee+dabs(zk(kk)-a*xk(kk)**2-b*yk(kk)**2
     &                  -c*xk(kk)*yk(kk)-d*xk(kk)-e*yk(kk))
            enddo
            if(ee.lt.e0)goto 400
300       continue
          exit !ループを抜ける
c
c         平均曲率Η=(κ1+κ2)/2を計算する
400       continue
          e0=ee
          phi0=phi
          psi0=psi
          eta0=eta
          x0=xx(1)
          y0=yy(1)
          z0=zz(1)
c
          p=2d0*a*x0+c*y0+d   !∂z/∂x
          q=2d0*b*y0+c*x0+e   !∂z/∂y
          r=2d0*a             !∂2z/∂x2
          s=c                 !∂2z/∂x∂y
          t=2d0*b             !∂z/∂y2
c         平均曲率:H
          h1=r*(1d0+q**2)-2d0*p*q*s+t*(1d0+p**2)
          h2=2d0*(1d0+p**2+q**2)**(3d0/2d0)
          h=h1/h2
c         Gaussの曲率:K
          gk=(r*t-s**2)/(1d0+p**2+q**2)**2
c         単位法ベクトルの方向余弦:ex,ey,ez →外向き
          e=dsqrt(1d0+p**2+q**2)
          ex=-p/e
          ey=-q/e
          ez=1d0/e
          call rotzyx(-ex,-ey,-ez,dnx,dny,dnz,-eta0,-psi0,-phi0)  !法線ベクトル
200     continue
c
        write(*,*)"ee:",sngl(e0)
        write(*,*)"phi,psi,eta:",sngl(phi0/PI*180d0),sngl(psi0/PI*180d0)
     &                                              ,sngl(eta0/PI*180d0)
c        write(*,*)"h:",h
        e=dsqrt(ex**2+ey**2+ez**2)
        write(*,*)"ex:",sngl(-ex),sngl(-ey),sngl(-ez),sngl(e)
        do k=1,nk
          x0=xx(k)-xx(1)
          y0=yy(k)-yy(1)
          z0=zz(k)-zz(1)
          call rotxyz(x0,y0,z0,x1,y1,z1,phi0,psi0,eta0)
          xk(k)=x1
          yk(k)=y1
          zk(k)=z1
        enddo
        call pltrset()
        call pltnt(0d0,0d0,1d0,5d0,2)  !赤
        call pltsetnp(nk,xk,yk,zk,ind)
c
        return
        end
c---------------------------------------------------------------------------
c
c       z=ax^2+by^2+cxy+dx+eyの係数を計算する
c
c---------------------------------------------------------------------------
        subroutine calccoefficient(nk,xk,yk,zk,a,b,c,d,e,ill)
        include "header.h"
        parameter(n=5,m=6)
        dimension sa(n,m),sf(n)
        dimension xk(ind),yk(ind),zk(ind)
c
c         ∑を計算する
          sgx4=0d0      !∑x^4
          sgy4=0d0      !∑y^4
          sgx3y1=0d0    !∑x^3y^1
          sgx2y2=0d0    !∑x^2y^2
          sgx1y3=0d0    !∑x^1y^3
c
          sgx3=0d0      !∑x^3
          sgy3=0d0      !∑y^3
          sgx2y1=0d0    !∑x^2y^1
          sgx1y2=0d0    !∑x^1y^2
c
          sgx2=0d0      !∑x^2
          sgy2=0d0      !∑y^2
          sgx1y1=0d0    !∑x^1y^1
c
          sgx2z=0d0     !∑x^2z
          sgy2z=0d0     !∑y^2z
          sgxyz=0d0     !∑xyz
c
          sgxz=0d0      !∑xz
          sgyz=0d0      !∑yz
          do kk=1,nk
            sgx4=sgx4+xk(kk)**4                 !∑x^4
            sgy4=sgy4+yk(kk)**4                 !∑y^4
            sgx3y1=sgx3y1+xk(kk)**3*yk(kk)**1   !∑x^3y^1
            sgx2y2=sgx2y2+xk(kk)**2*yk(kk)**2   !∑x^2y^2
            sgx1y3=sgx1y3+xk(kk)**1*yk(kk)**3   !∑x^1y^3
c
            sgx3=sgx3+xk(kk)**3                 !∑x^3
            sgy3=sgy3+yk(kk)**3                 !∑y^3
            sgx2y1=sgx2y1+xk(kk)**2*yk(kk)**1   !∑x^2y^1
            sgx1y2=sgx1y2+xk(kk)**1*yk(kk)**2   !∑x^1y^2
c
            sgx2=sgx2+xk(kk)**2                 !∑x^2
            sgy2=sgy2+yk(kk)**2                 !∑y^2
            sgx1y1=sgx1y1+xk(kk)**1*yk(kk)**1   !∑x^1y^1
c
            sgx2z=sgx2z+xk(kk)**2*zk(kk)        !∑x^2z
            sgy2z=sgy2z+yk(kk)**2*zk(kk)        !∑y^2z
            sgxyz=sgxyz+xk(kk)*yk(kk)*zk(kk)    !∑xyz
c
            sgxz=sgxz+xk(kk)*zk(kk)             !∑xz
            sgyz=sgyz+yk(kk)*zk(kk)             !∑yz
          enddo
c
c         z=ax^2+by^2+cxy+dx+eyの係数a,b,c,d,eを計算する
c         左辺
          sa(1,1)=sgx4
          sa(1,2)=sgx2y2
          sa(1,3)=sgx3y1
          sa(1,4)=sgx3
          sa(1,5)=sgx2y1
          sa(2,1)=sgx2y2
          sa(2,2)=sgy4
          sa(2,3)=sgx1y3
          sa(2,4)=sgx1y2
          sa(2,5)=sgy3
          sa(3,1)=sgx3y1
          sa(3,2)=sgx1y3
          sa(3,3)=sgx2y2
          sa(3,4)=sgx2y1
          sa(3,5)=sgx1y2
          sa(4,1)=sgx3
          sa(4,2)=sgx1y2
          sa(4,3)=sgx2y1
          sa(4,4)=sgx2
          sa(4,5)=sgx1y1
          sa(5,1)=sgx2y1
          sa(5,2)=sgy3
          sa(5,3)=sgx1y2
          sa(5,4)=sgx1y1
          sa(5,5)=sgy2
c         右辺
          sa(1,6)=sgx2z
          sa(2,6)=sgy2z
          sa(3,6)=sgxyz
          sa(4,6)=sgxz
          sa(5,6)=sgyz
c
          eps=1d-6  !0とみなす値
c
          call gauss2(sa,n,m,eps,sf,ill)  !n:行,m:列
          if(ill.eq.1)return
c
          a=sa(1,6)
          b=sa(2,6)
          c=sa(3,6)
          d=sa(4,6)
          e=sa(5,6)
        return
        end
c---------------------------------------------------------------------------
c
c         kpの周囲の節点を探す
c
c---------------------------------------------------------------------------
        subroutine nssurf(kp,nels,ms,xs,ys,zs,nk,xt,yt,zt)
        include "header.h"
        dimension ms(ine,3)
        dimension mk(ind),mk0(ind)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xt(ind),yt(ind),zt(ind)
c
c       kpの周囲の節点(kp含む)
        nk=1
        mk(nk)=kp
        xt(nk)=xs(kp)
        yt(nk)=ys(kp)
        zt(nk)=zs(kp)
        do 100 i=1,nels
          if(ms(i,1).ne.kp .and. ms(i,2).ne.kp
     &                     .and. ms(i,3).ne.kp)goto 100
          do 200 j=1,3
            do kk=1,nk
              if(ms(i,j).eq.mk(kk))goto 200
            enddo
            nk=nk+1
            mk(nk)=ms(i,j)
            xt(nk)=xs(mk(nk))
            yt(nk)=ys(mk(nk))
            zt(nk)=zs(mk(nk))
200       continue
100     continue
c
c       mk(nk)の周囲の点
        nk0=nk
        do 300 kk=1,nk
300     mk0(kk)=mk(kk)
c
        do 400 i=1,nels
          ifg=0
          do 600 k0=1,nk0
            if(ms(i,1).eq.mk0(k0) .or. ms(i,2).eq.mk0(k0)
     &                       .or. ms(i,3).eq.mk0(k0))ifg=1
600       continue
          if(ifg.eq.0)goto 400
          do 500 j=1,3
            do kk=1,nk
              if(ms(i,j).eq.mk(kk))goto 500
            enddo
c            nk=nk+1
c            mk(nk)=ms(i,j)
c            xt(nk)=xs(mk(nk))
c            yt(nk)=ys(mk(nk))
c            zt(nk)=zs(mk(nk))
500       continue
400     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c         曲率の計算:(x-xc)^2+(y-yc)^2+(z-zc)^2=r^2
c
c---------------------------------------------------------------------------
        subroutine calcsphere(nk,xx,yy,zz,dnx,dny,dnz,ee,h)
        include "header.h"
        dimension xx(ind),yy(ind),zz(ind)
c
        x1=xx(1)    !面上の節点
        y1=yy(1)    !面上の節点
        z1=zz(1)    !面上の節点
        a1=0d0
        a2=0d0
        a3=0d0
        b1=0d0
        b2=0d0
        b3=0d0
        c1=0d0
        c2=0d0
        c3=0d0
        d1=0d0
        d2=0d0
        d3=0d0
c
c      係数を計算する
        do kk=1,nk
          xi=xx(kk)
          yi=yy(kk)
          zi=zz(kk)
c
          a1=a1+(x1-xi)*(x1-xi)
          a2=a2+(y1-yi)*(x1-xi)
          a3=a3+(z1-zi)*(x1-xi)
c
          b1=b1+(x1-xi)*(y1-yi)
          b2=b2+(y1-yi)*(y1-yi)
          b3=b3+(z1-zi)*(y1-yi)
c
          c1=c1+(x1-xi)*(z1-zi)
          c2=c2+(y1-yi)*(z1-zi)
          c3=c3+(z1-zi)*(z1-zi)
c
          d1=d1-(x1-xi)*(xi**2-x1**2+yi**2-y1**2+zi**2-z1**2)/2d0
          d2=d2-(y1-yi)*(xi**2-x1**2+yi**2-y1**2+zi**2-z1**2)/2d0
          d3=d3-(z1-zi)*(xi**2-x1**2+yi**2-y1**2+zi**2-z1**2)/2d0
        enddo
c
c       球の中心を計算する
        xyz=a1*b2*c3+a2*b3*c1+a3*b1*c2
     &     -a3*b2*c1-a2*b1*c3-a1*b3*c2
        if(xyz.eq.0d0)then
          h=0d0
          dnx=0d0
          dny=0d0
          dnz=0d0
          return
        endif
        xc=(b1*c2*d3+b2*c3*d1+b3*c1*d2-b3*c2*d1-b2*c1*d3-b1*c3*d2)/xyz
        yc=(a1*d2*c3+a2*d3*c1+a3*d1*c2-a3*d2*c1-a2*d1*c3-a1*d3*c2)/xyz
        zc=(a1*b2*d3+a2*b3*d1+a3*b1*d2-a3*b2*d1-a2*b1*d3-a1*b3*d2)/xyz
c
c       平均曲率と法線の方向余弦を計算する
        r=dsqrt((x1-xc)**2+(y1-yc)**2+(z1-zc)**2)
        h=1d0/r
        dnx=xc-x1
        dny=yc-y1
        dnz=zc-z1
        dn=dsqrt(dnx**2+dny**2+dnz**2)
c        dnx=dnx/dn !内向きを正にする
c        dny=dny/dn !内向きを正にする
c        dnz=dnz/dn !内向きを正にする
c
c       誤差を計算する
        ee=0d0
        do kk=1,nk
          xi=xx(kk)
          yi=yy(kk)
          zi=zz(kk)
          e=(xi-xc)**2+(yi-yc)**2+(zi-zc)**2-r**2
          ee=ee+dabs(e)
        enddo
c
        return
        end
