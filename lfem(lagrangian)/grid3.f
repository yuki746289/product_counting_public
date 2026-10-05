c---------------------------------------------------------------------------
c
c       4面体の節点位置を調節する
c
c---------------------------------------------------------------------------
        subroutine remesh2(np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
        include "header.h"
        dimension ne(ine,4),ne0(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension x0(4),y0(4),z0(4)
        dimension nelb(ine),neb(ine)
        dimension nec(ine,4)
        dimension mm(ind),mk1(ind),mk2(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind),cc(ind),tt(ind)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension jm(4,3)
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
        data ceft1/10d-1/   !緩和係数:内部の節点
        data ceft2/10d-1/   !緩和係数:表面の節点
        data drstd/6d-1/   !基準値
c
        do j=1,3
          jm(1,j)=jm1(j)
          jm(2,j)=jm2(j)
          jm(3,j)=jm3(j)
          jm(4,j)=jm4(j)
        enddo
c
c       隣接する要素
        do 10 kele=1,nele
        do 10 j=1,4
10      nec(kele,j)=0
c
        do 200 i=1,nele
        do 200 j=1,4
          if(nec(i,j).ne.0)goto 200
          m1=ne(i,jm(j,1))
          m2=ne(i,jm(j,2))
          m3=ne(i,jm(j,3))
          do 100 k=i+1,nele
          do 100 l=1,4
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
100       continue
200     continue
c
c       表面の節点と内部の節点に分ける
        do 600 kp=1,np
600     mm(kp)=0    !0:内部の節点,1:表面の節点
        do 700 i=1,nele
        do 700 j=1,4
          if(nec(i,j).ne.0)goto 700
          mm(ne(i,jm(j,1)))=1
          mm(ne(i,jm(j,2)))=1
          mm(ne(i,jm(j,3)))=1
700     continue
c
c       節点を並べ替える:内部の節点
        nk1=0
        do 800 kp=1,np
          if(mm(kp).eq.1)goto 800
          nk1=nk1+1
          mk1(nk1)=kp
800     continue
c
c       節点を並べ替える:表面の節点
        nk2=0
        do 900 kp=1,np
          if(mm(kp).eq.0)goto 900
          nk2=nk2+1
          mk2(nk2)=kp
900     continue
c
c       -------------------------------------------------------------
c
c       節点の位置を最適化する:内部の節点
        write(*,*)"on inner nodes"
        icount=0    !内部の節点を修正した回数:10回最適化する
500     continue
        do 300 kk=1,nk1
          kp=mk1(kk)
c          write(*,*)"kp:",kp,"/",np
c         kpを含む要素
          call npsurf(kp,nele,ne,nb,nelb,neb)
c
          drmax=0d0
          do kb=1,nb
            do i=1,4
              x0(i)=xx(ne(nelb(kb),i))
              y0(i)=yy(ne(nelb(kb),i))
              z0(i)=zz(ne(nelb(kb),i))
            enddo
c            dr=difdr3d(x0,y0,z0)
            call difdr4d(x0,y0,z0,dr,dr1,dr2)
            if(dr.gt.drmax)drmax=dr
          enddo
c          if(drmax.lt.drstd)goto 300   !基準値未満の節点は考慮しない
c
c          write(*,*)"befoer optimize:inner"
c          write(*,*)"drmax:",sngl(drmax)
c          write(*,*)"xx:",sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c          do 400 kb=1,nb
c          do 400 j=1,4
c400       ne0(kb,j)=ne(nelb(kb),j)
c          call pltsetnp(np,xx,yy,zz,ind)
c          call pltsetne(nb,ne0,2d0,2,ine)
c          pause
c         最適な節点を探す
          call decidenp(kp,np,nele,ne,nb,nelb,neb,xx,yy,zz
     &                                        ,xc,yc,zc,drmin,1d0)
          if(drmin.gt.drmax)goto 300   !最適化されない場合は飛ばす
c         節点kpの速度、圧力、座標を補完する
          xp=xx(kp)*(1d0-ceft1)+xc*ceft1
          yp=yy(kp)*(1d0-ceft1)+yc*ceft1
          zp=zz(kp)*(1d0-ceft1)+zc*ceft1
          call interpolation2(kp,xp,yp,zp
     &                   ,np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
c
c          write(*,*)"after optimize"
c          drmax=0d0
c          do kb=1,nb
c            do i=1,4
c              x0(i)=xx(ne(nelb(kb),i))
c              y0(i)=yy(ne(nelb(kb),i))
c              z0(i)=zz(ne(nelb(kb),i))
c            enddo
cc            dr=difdr3d(x0,y0,z0)
c            call difdr4d(x0,y0,z0,dr,dr1,dr2)
c            if(dr.gt.drmax)drmax=dr
c          enddo
c          write(*,*)"drmax:",kb,sngl(drmax)
c          write(*,*)"xx:",sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c          call pltsetnp(np,xx,yy,zz,ind)
c          call pltsetne(nb,ne0,2d0,2,ine)
c          pause
300     continue
c
c       -------------------------------------------------------------
c
c       形状係数の最大値を計算する
c        drmax2=-1d0
c        do 1300 kp=1,np
c         kpを含む要素
c          call npsurf(kp,nele,ne,nb,nelb,neb)
c
c          drmax=0d0
c          do kb=1,nb    !節点kpを含む要素の数
c            do i=1,4
c              x0(i)=xx(ne(nelb(kb),i))
c              y0(i)=yy(ne(nelb(kb),i))
c              z0(i)=zz(ne(nelb(kb),i))
c            enddo
cc            dr=difdr3d(x0,y0,z0)
c            call difdr4d(x0,y0,z0,dr,dr1,dr2)
c            if(dr.gt.drmax)drmax=dr
c          enddo
c          if(drmax.gt.drmax2)drmax2=drmax
c1300    continue
        drmax2=elshape2(np,nele,ne,xx,yy,zz,ind,ine)
        if(drmax2.gt.drstd)then
          icount=icount+1   !内部の節点を修正した回数
          if(icount.le.10)goto 500   !最初から最適化し直す
        else
          return   !最適化終了
        endif
c
c       -------------------------------------------------------------
c
c       節点の位置を最適化する:表面の節点
        write(*,*)"on surface nodes"
1200    continue
        do 1000 kk=1,nk2
          kp=mk2(kk)
c         kpを含む要素
          call npsurf(kp,nele,ne,nb,nelb,neb)
c
          drmax=0d0
          do kb=1,nb
            do i=1,4
              x0(i)=xx(ne(nelb(kb),i))
              y0(i)=yy(ne(nelb(kb),i))
              z0(i)=zz(ne(nelb(kb),i))
            enddo
c            dr=difdr3d(x0,y0,z0)
            call difdr4d(x0,y0,z0,dr,dr1,dr2)
            if(dr.gt.drmax)drmax=dr
          enddo
          if(drmax.lt.drstd)goto 1000  !基準値未満の節点は考慮しない
c
c          write(*,*)"befoer optimize:surface"
c          write(*,*)"drmax:",sngl(drmax)
c          write(*,*)"xx:",sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c          do 1100 kb=1,nb
c          do 1100 j=1,4
c1100       ne0(kb,j)=ne(nelb(kb),j)
c          call pltsetnp(np,xx,yy,zz,ind)
c          call pltsetne(nb,ne0,2d0,2,ine)
c          pause
c         最適な節点を探す
          call decidenp(kp,np,nele,ne,nb,nelb,neb,xx,yy,zz
     &                                        ,xc,yc,zc,drmin,1.5d0)
          if(drmin.gt.drmax)goto 1000   !最適化されない場合は飛ばす
c         節点kpの速度、圧力、座標を補完する
          xp=xx(kp)*(1d0-ceft2)+xc*ceft2
          yp=yy(kp)*(1d0-ceft2)+yc*ceft2
          zp=zz(kp)*(1d0-ceft2)+zc*ceft2
          call interpolation2(kp,xp,yp,zp
     &                   ,np,ne,nele,xx,yy,zz,uu,vv,ww,pp,cc,tt)
c
c          write(*,*)"after optimize"
c          drmax=0d0
c          do kb=1,nb
c            do i=1,4
c              x0(i)=xx(ne(nelb(kb),i))
c              y0(i)=yy(ne(nelb(kb),i))
c              z0(i)=zz(ne(nelb(kb),i))
c            enddo
cc            dr=difdr3d(x0,y0,z0)
c            call difdr4d(x0,y0,z0,dr,dr1,dr2)
c            if(dr.gt.drmax)drmax=dr
c          enddo
c          write(*,*)"drmax:",kb,sngl(drmax)
c          write(*,*)"xx:",sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c          call pltsetnp(np,xx,yy,zz,ind)
c          call pltsetne(nb,ne0,2d0,2,ine)
c          pause
1000     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       kpを含む要素
c
c---------------------------------------------------------------------------
        subroutine npsurf(kp,nele,ne,nb,nelb,neb)
        include "header.h"
        dimension ne(ine,4)
        dimension nelb(ine),neb(ine)
c
        nb=0
        do 100 i=1,nele
        do 100 j=1,4
          if(ne(i,j).ne.kp)goto 100
          nb=nb+1
          nelb(nb)=i
          neb(nb)=j
100     continue
c
        return
        end
c---------------------------------------------------------------------------
c
c       最適な節点を探す
c
c---------------------------------------------------------------------------
        subroutine decidenp(kp,np,nele,ne,nb,nelb,neb,xx,yy,zz,xc,yc,zc
     &    ,drmin,regct)
        include "header.h"
        dimension ne(ine,4)
        dimension nelb(ine),neb(ine)
        dimension xx(ind),yy(ind),zz(ind)
        dimension x0(4),y0(4),z0(4)
c       -------表示用-------
        dimension ne0(ine,4)
c
        xmax=-1d0
        ymax=-1d0
        zmax=-1d0
        xmin= 1d0
        ymin= 1d0
        zmin= 1d0
        do 100 kb=1,nb
        do 100 j=1,4
          if(xmax.lt.xx(ne(nelb(kb),j)))xmax=xx(ne(nelb(kb),j))
          if(ymax.lt.yy(ne(nelb(kb),j)))ymax=yy(ne(nelb(kb),j))
          if(zmax.lt.zz(ne(nelb(kb),j)))zmax=zz(ne(nelb(kb),j))
          if(xmin.gt.xx(ne(nelb(kb),j)))xmin=xx(ne(nelb(kb),j))
          if(ymin.gt.yy(ne(nelb(kb),j)))ymin=yy(ne(nelb(kb),j))
          if(zmin.gt.zz(ne(nelb(kb),j)))zmin=zz(ne(nelb(kb),j))
100     continue
        xmax=(xmax+xmin)/2d0+(xmax-xmin)/2d0*regct
        ymax=(ymax+ymin)/2d0+(ymax-ymin)/2d0*regct
        zmax=(zmax+zmin)/2d0+(zmax-zmin)/2d0*regct
        xmim=(xmax+xmin)/2d0-(xmax-xmin)/2d0*regct
        ymim=(ymax+ymin)/2d0-(ymax-ymin)/2d0*regct
        zmim=(zmax+zmin)/2d0-(zmax-zmin)/2d0*regct

        xmax0=xmax
        ymax0=ymax
        zmax0=zmax
        xmim0=xmin
        ymim0=ymin
        zmim0=zmin
c        write(*,*)"xmax",sngl(xmax),sngl(ymax),sngl(zmax)
c        write(*,*)"xmin",sngl(xmin),sngl(ymin),sngl(zmin)
c
c       10分割を5回繰り返す
        nx=10
        ny=10
        nz=10
        drmin=1d3
        do nn=1,10
          iflag=0
          do 200 kx=0,nx-1
          do 200 ky=0,ny-1
          do 200 kz=0,nz-1
            x=xmin+(xmax-xmin)*(dble(kx)+5d-1)/dble(nx)
            y=ymin+(ymax-ymin)*(dble(ky)+5d-1)/dble(ny)
            z=zmin+(zmax-zmin)*(dble(kz)+5d-1)/dble(nz)
            dr=-1d3
            do kb=1,nb
              do i=1,4
                x0(i)=xx(ne(nelb(kb),i))
                y0(i)=yy(ne(nelb(kb),i))
                z0(i)=zz(ne(nelb(kb),i))
              enddo
              x0(neb(kb))=x
              y0(neb(kb))=y
              z0(neb(kb))=z
              vl=calvl2(x0(1),y0(1),z0(1),x0(2),y0(2),z0(2)
     &                 ,x0(3),y0(3),z0(3),x0(4),y0(4),z0(4))
              if(vl.lt.1d-7)goto 200
c              dr2=difdr3d(x0,y0,z0)
              call difdr4d(x0,y0,z0,dr2,dummy1,dummy2)
              if(dr2.gt.dr)dr=dr2
            enddo
c            write(*,*)"dr:",dr,drmin
            if(dr.ge.drmin)goto 200
            drmin=dr
            nxmin=kx
            nymin=ky
            nzmin=kz
            xc=x
            yc=y
            zc=z
            iflag=1
200       continue
          if(iflag.eq.0)then
            if(nn.ne.1)goto 300 !ループを抜ける
            write(*,*)"the appropriate point was not found in decidenp."
            write(*,*)"nn,drmin:",nn,drmin
c
            xmax=xmax0
            ymax=ymax0
            zmax=zmax0
            xmim=xmin0
            ymim=ymin0
            zmim=zmin0
            write(*,*)"xmax",sngl(xmax),sngl(ymax),sngl(zmax)
            write(*,*)"xmin",sngl(xmin),sngl(ymin),sngl(zmin)
            write(*,*)"d:",sngl(xmax-xmin),sngl(ymax-ymin)
     &                     ,sngl(zmax-zmin)
c
            do kb=1,nb
              do i=1,4
                x0(i)=xx(ne(nelb(kb),i))
                y0(i)=yy(ne(nelb(kb),i))
                z0(i)=zz(ne(nelb(kb),i))
              enddo
              vl=calvl2(x0(1),y0(1),z0(1),x0(2),y0(2),z0(2)
     &                 ,x0(3),y0(3),z0(3),x0(4),y0(4),z0(4))  !符号を考慮
              call difdr4d(x0,y0,z0,dr,dummy1,dummy2)
              write(*,*)"kb:",kb,"/",nb,sngl(vl),sngl(dr)
            enddo
c
            do 500 kb=1,nb
            do 500 j=1,4
500         ne0(kb,j)=ne(nelb(kb),j)
            call pltsetnp(np,xx,yy,zz,ind)
            call pltsetne(nb,ne0,2d0,2,ine)
            pause
c
            write(*,*)"nx,dr:",nx0,ny0,nz0,sngl(dr)
            stop
          endif
          xmin=xmin+(xmax-xmin)*dble(nxmin)/dble(nx)
          ymin=ymin+(ymax-ymin)*dble(nymin)/dble(ny)
          zmin=zmin+(zmax-zmin)*dble(nzmin)/dble(nz)
          xmax=xmin+(xmax-xmin)*dble(nxmin+1)/dble(nx)
          ymax=ymin+(ymax-ymin)*dble(nymin+1)/dble(ny)
          zmax=zmin+(zmax-zmin)*dble(nzmin+1)/dble(nz)
        enddo
300     continue
        if(drmin.eq.1d3)write(*,*)"the appropriate node was not found."
        if(drmin.eq.1d3)stop
c
        return
        end
c---------------------------------------------------------------------------
c
c       要素の形状係数を計算する
c
c---------------------------------------------------------------------------
        function difdr3d(x0,y0,z0)
        implicit double precision(a-h,o-z)
        dimension x0(4),y0(4),z0(4)
c
        xg=(x0(1)+x0(2)+x0(3)+x0(4))/3d0
        yg=(y0(1)+y0(2)+y0(3)+y0(4))/3d0
        zg=(z0(1)+z0(2)+z0(3)+z0(4))/3d0
        r1=dsqrt((x0(1)-xg)**2+(y0(1)-yg)**2+(z0(1)-zg)**2)
        r2=dsqrt((x0(2)-xg)**2+(y0(2)-yg)**2+(z0(2)-zg)**2)
        r3=dsqrt((x0(3)-xg)**2+(y0(3)-yg)**2+(z0(3)-zg)**2)
        r4=dsqrt((x0(4)-xg)**2+(y0(4)-yg)**2+(z0(4)-zg)**2)
c
        a=r1+r2+r3+r4
        b=r1**2+r2**2+r3**2+r4**2
        difdr3d=2d0/a*dsqrt(4d0*b-a**2)

c        a=r(1)+r(2)+r(3)+r(4)
c        b=r(1)**2+r(2)**2+r(3)**2+r(4)**2
c        dr=2d0*dsqrt(4d0*b-a**2)/a
c
        return
        end
c---------------------------------------------------------------------------
c
c       要素の形状係数を計算する
c
c---------------------------------------------------------------------------
        subroutine difdr4d(x0,y0,z0,dr,dr1,dr2)
        implicit double precision(a-h,o-z)
        dimension x0(4),y0(4),z0(4)
c
        xg=(x0(1)+x0(2)+x0(3)+x0(4))/4d0
        yg=(y0(1)+y0(2)+y0(3)+y0(4))/4d0
        zg=(z0(1)+z0(2)+z0(3)+z0(4))/4d0
c
        call circumct3(x0,y0,z0,xc,yc,zc,rr)
c
        dr1=dsqrt((xc-xg)**2+(yc-yg)**2+(zc-zg)**2)
        dr2=rr
        dr=dr1/dr2  !dirdr4=[0,1)
c
        return
        end
c---------------------------------------------------------------------------
c
c       4面体の外心を計算する
c
c---------------------------------------------------------------------------
        subroutine circumct3(x0,y0,z0,xc,yc,zc,rr)
        implicit double precision(a-h,o-z)
        dimension x0(4),y0(4),z0(4)
c
        x1=x0(1)
        x2=x0(2)
        x3=x0(3)
        x4=x0(4)
        y1=y0(1)
        y2=y0(2)
        y3=y0(3)
        y4=y0(4)
        z1=z0(1)
        z2=z0(2)
        z3=z0(3)
        z4=z0(4)
c
        a1=2d0*(x1-x2)
        a2=2d0*(x1-x3)
        a3=2d0*(x1-x4)
        b1=2d0*(y1-y2)
        b2=2d0*(y1-y3)
        b3=2d0*(y1-y4)
        c1=2d0*(z1-z2)
        c2=2d0*(z1-z3)
        c3=2d0*(z1-z4)
        d1=x1**2-x2**2+y1**2-y2**2+z1**2-z2**2
        d2=x1**2-x3**2+y1**2-y3**2+z1**2-z3**2
        d3=x1**2-x4**2+y1**2-y4**2+z1**2-z4**2
c       外心の座標
        xyz=a1*b2*c3+a2*b3*c1+a3*b1*c2-a3*b2*c1-a2*b1*c3-a1*b3*c2
        xc=(b1*c2*d3+b2*c3*d1+b3*c1*d2-b3*c2*d1-b2*c1*d3-b1*c3*d2)/xyz
        yc=(a1*d2*c3+a2*d3*c1+a3*d1*c2-a3*d2*c1-a2*d1*c3-a1*d3*c2)/xyz
        zc=(a1*b2*d3+a2*b3*d1+a3*b1*d2-a3*b2*d1-a2*b1*d3-a1*b3*d2)/xyz
c       外接球の半径
        r1=dsqrt((x1-xc)**2+(y1-yc)**2+(z1-zc)**2)
        r2=dsqrt((x2-xc)**2+(y2-yc)**2+(z2-zc)**2)
        r3=dsqrt((x3-xc)**2+(y3-yc)**2+(z3-zc)**2)
        r4=dsqrt((x4-xc)**2+(y4-yc)**2+(z4-zc)**2)
c       検算
        ifg=0
        if(dabs(r1-r2).gt.1d-7)ifg=1
        if(dabs(r2-r3).gt.1d-7)ifg=1
        if(dabs(r3-r4).gt.1d-7)ifg=1
        if(dabs(r4-r1).gt.1d-7)ifg=1
        if(ifg.eq.1)write(*,*)"in circumct3 ifg:",ifg
        if(ifg.eq.1)then
          write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
          write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
          write(*,*)"x3:",sngl(x3),sngl(y3),sngl(z3)
          write(*,*)"x4:",sngl(x4),sngl(y4),sngl(z4)
          write(*,*)"r1:",r1
          write(*,*)"r2:",r2
          write(*,*)"r3:",r3
          write(*,*)"r4:",r4
          write(*,*)"xc:",sngl(xc),sngl(yc),sngl(zc)
        endif
        if(ifg.eq.1)stop
        rr=(r1+r2+r3+r4)/4d0
c
        return
        end
c---------------------------------------------------------------------------
c
c       全4面体要素中から最小体積の要素を探して、その最小体積を返す
c
c---------------------------------------------------------------------------
        function calvlmin(np,ne,nele,xx,yy,zz)
        include "header.h"
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
c
        vlmin=999d0
        do i=1,nele
          x1=xx(ne(i,1))
          x2=xx(ne(i,2))
          x3=xx(ne(i,3))
          x4=xx(ne(i,4))
          y1=yy(ne(i,1))
          y2=yy(ne(i,2))
          y3=yy(ne(i,3))
          y4=yy(ne(i,4))
          z1=zz(ne(i,1))
          z2=zz(ne(i,2))
          z3=zz(ne(i,3))
          z4=zz(ne(i,4))
          vl=calvl2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          if(vl.lt.vlmin)vlmin=vl
        enddo
        calvlmin=vlmin
c
        return
        end
c---------------------------------------------------------------------------
c
c       要素の変形:界面に2つの面が存在する要素を分割する
c
c---------------------------------------------------------------------------
        subroutine divne(np,nele,ne,xx,yy,zz,uu,vv,ww,pp)
        include "header.h"
        dimension ne(ine,4)
        dimension nec(ine,4)
        dimension me0(ine),mele(ine)
        dimension xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        dimension ml(ine,2)
c
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension jm(4,3)
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
c
        do j=1,3
          jm(1,j)=jm1(j)
          jm(2,j)=jm2(j)
          jm(3,j)=jm3(j)
          jm(4,j)=jm4(j)
        enddo
c
c       隣接する要素
        do 10 kele=1,nele
        do 10 j=1,4
10      nec(kele,j)=0
c
        do 200 i=1,nele
        do 200 j=1,4
          if(nec(i,j).ne.0)goto 200
          m1=ne(i,jm(j,1))
          m2=ne(i,jm(j,2))
          m3=ne(i,jm(j,3))
          do 100 k=i+1,nele
          do 100 l=1,4
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
100       continue
200     continue
c
c       界面に2つの面が存在する要素の内部の辺を探す
        nel0=0
        do 300 kele=1,nele
          n=0
          if(nec(kele,1).eq.0)n=n+1
          if(nec(kele,2).eq.0)n=n+1
          if(nec(kele,3).eq.0)n=n+1
          if(nec(kele,4).eq.0)n=n+1
          if(n.ne.2)goto 300
c
          if(nec(kele,1).eq.0 .and. nec(kele,2).eq.0)then
            l1=ne(kele,1)
            l2=ne(kele,2)
          endif
          if(nec(kele,1).eq.0 .and. nec(kele,3).eq.0)then
            l1=ne(kele,1)
            l2=ne(kele,3)
          endif
          if(nec(kele,1).eq.0 .and. nec(kele,4).eq.0)then
            l1=ne(kele,1)
            l2=ne(kele,4)
          endif
          if(nec(kele,2).eq.0 .and. nec(kele,3).eq.0)then
            l1=ne(kele,2)
            l2=ne(kele,3)
          endif
          if(nec(kele,2).eq.0 .and. nec(kele,4).eq.0)then
            l1=ne(kele,2)
            l2=ne(kele,4)
          endif
          if(nec(kele,3).eq.0 .and. nec(kele,4).eq.0)then
            l1=ne(kele,3)
            l2=ne(kele,4)
          endif
c
c         すでに登録されていないか確認する
          do kel0=1,nel0
            if(l1.eq.ml(kel0,1) .and. l2.eq.ml(kel0,2))goto 300
            if(l1.eq.ml(kel0,2) .and. l2.eq.ml(kel0,1))goto 300
          enddo
c
c         登録する
          nel0=nel0+1
          ml(nel0,1)=l1
          ml(nel0,2)=l2
300     continue
c
c       要素を分割する
        do kel0=1,nel0
          l1=ml(kel0,1)     !分割する辺の節点
          l2=ml(kel0,2)     !分割する辺の節点
c
c         辺が存在する要素を探す
          iflag=1
          nl=0
          do kele=1,nele
            k1=ne(kele,1)
            k2=ne(kele,2)
            k3=ne(kele,3)
            k4=ne(kele,4)
            if(k1.eq.l1 .and. k2.eq.l2)iflag=1
            if(k1.eq.l1 .and. k3.eq.l2)iflag=1
            if(k1.eq.l1 .and. k4.eq.l2)iflag=1
            if(k2.eq.l1 .and. k3.eq.l2)iflag=1
            if(k2.eq.l1 .and. k4.eq.l2)iflag=1
            if(k3.eq.l1 .and. k4.eq.l2)iflag=1
            if(k1.eq.l2 .and. k2.eq.l1)iflag=1
            if(k1.eq.l2 .and. k3.eq.l1)iflag=1
            if(k1.eq.l2 .and. k4.eq.l1)iflag=1
            if(k2.eq.l2 .and. k3.eq.l1)iflag=1
            if(k2.eq.l2 .and. k4.eq.l1)iflag=1
            if(k3.eq.l2 .and. k4.eq.l1)iflag=1
            if(iflag.eq.1)nl=nl+1
            if(iflag.eq.1)mele(nl)=kele
          enddo
c
c         要素を分割する
          call divne2(np,nele,ne,xx,yy,zz,uu,vv,ww,pp,nl,mele,l1,l2)
        enddo
c
        return
        end
c
c         要素を分割する
c
        subroutine divne2(np,nele,ne,xx,yy,zz,uu,vv,ww,pp,nl,mele,l1,l2)
        include "header.h"
        dimension ne(ine,4),ne0(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        dimension mele(ine)
c
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension jm(4,3)
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
c
        do kl=1,nl
          kele=mele(kl)
          k1=ne(kele,1)
          k2=ne(kele,2)
          k3=ne(kele,3)
          k4=ne(kele,4)
c
c         節点を登録する
          x=(xx(l1)+xx(l2))/2d0
          y=(yy(l1)+yy(l2))/2d0
          z=(zz(l1)+zz(l2))/2d0
          np=np+1
          xx(np)=x
          yy(np)=y
          zz(np)=z
c
c         要素を分割する
c         辺1-2
          if(k1.eq.l1 .and. k2.eq.l2 .or. k1.eq.l2 .and. k2.eq.l1)then
            !面1
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm1(1))
            ne(nele,3)=ne(kele,jm1(2))
            ne(nele,4)=ne(kele,jm1(3))
            !面2
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm2(1))
            ne(nele,3)=ne(kele,jm2(2))
            ne(nele,4)=ne(kele,jm2(3))
          endif
c
c         辺1-3
          if(k1.eq.l1 .and. k3.eq.l2 .or. k1.eq.l2 .and. k3.eq.l1)then
            !面1
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm1(1))
            ne(nele,3)=ne(kele,jm1(2))
            ne(nele,4)=ne(kele,jm1(3))
            !面3
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm3(1))
            ne(nele,3)=ne(kele,jm3(2))
            ne(nele,4)=ne(kele,jm3(3))
          endif
c
c         辺1-4
          if(k1.eq.l1 .and. k4.eq.l2 .or. k1.eq.l2 .and. k4.eq.l1)then
            !面1
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm1(1))
            ne(nele,3)=ne(kele,jm1(2))
            ne(nele,4)=ne(kele,jm1(3))
            !面4
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm4(1))
            ne(nele,3)=ne(kele,jm4(2))
            ne(nele,4)=ne(kele,jm4(3))
          endif
c
c         辺2-3
          if(k2.eq.l1 .and. k3.eq.l2 .or. k3.eq.l2 .and. k2.eq.l1)then
            !面2
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm2(1))
            ne(nele,3)=ne(kele,jm2(2))
            ne(nele,4)=ne(kele,jm2(3))
            !面3
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm3(1))
            ne(nele,3)=ne(kele,jm3(2))
            ne(nele,4)=ne(kele,jm3(3))
          endif
c
c         辺2-4
          if(k2.eq.l1 .and. k4.eq.l2 .or. k4.eq.l2 .and. k2.eq.l1)then
            !面2
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm2(1))
            ne(nele,3)=ne(kele,jm2(2))
            ne(nele,4)=ne(kele,jm2(3))
            !面4
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm4(1))
            ne(nele,3)=ne(kele,jm4(2))
            ne(nele,4)=ne(kele,jm4(3))
          endif
c
c         辺3-4
          if(k3.eq.l1 .and. k4.eq.l2 .or. k4.eq.l2 .and. k3.eq.l1)then
            !面3
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm3(1))
            ne(nele,3)=ne(kele,jm3(2))
            ne(nele,4)=ne(kele,jm3(3))
            !面4
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=ne(kele,jm4(1))
            ne(nele,3)=ne(kele,jm4(2))
            ne(nele,4)=ne(kele,jm4(3))
          endif
c
c         要素keleを削除する
          nel0=nele
          do 100 i=1,nele
            if(i.eq.kele)goto 100
            nel0=nel0+1
            do 200 j=1,4
200         ne0(nel0,j)=ne(i,j)
100       continue
c
          nele=nel0
          do 300 i=1,nel0
          do 300 j=1,4
300       ne(i,j)=ne0(i,j)
c
        enddo
c
        return
        end
