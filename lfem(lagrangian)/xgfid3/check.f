c---------------------------------------------------------------------------
c
c      作り直す要素を描く
c
c---------------------------------------------------------------------------
       subroutine pltne0(na,ma,ne,size,icolor,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),ma(ine)
c
       call fcolor(icolor,red,green,blue)
c
c      データを書く
c       open(1,file='setpne.res',position='REWIND')
       open(1,file='setne.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)na
       do 300 ka=1,na
300    write(1,*)(ne(ma(ka),j),j=1,4)
       close(1)
c
       return
       end
c---------------------------------------------------------------------------
c
c      要素を1つ描く
c
c---------------------------------------------------------------------------
       subroutine pltne1(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
     &                                    ,size,icolor,ind,ine)
       implicit double precision(a-h,o-z)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pne.res',position='REWIND')
       open(1,file='pne.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)4,1
       write(1,*)1,2,3,4
       write(1,*)sngl(x1),sngl(y1),sngl(z1)
       write(1,*)sngl(x2),sngl(y2),sngl(z2)
       write(1,*)sngl(x3),sngl(y3),sngl(z3)
       write(1,*)sngl(x4),sngl(y4),sngl(z4)
c
       return
       end
c---------------------------------------------------------------------------
c
c      体積確認
c
c---------------------------------------------------------------------------

        subroutine chkvol(nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
c
        do kele=1,nele
          x1=xx(ne(kele,1))
          y1=yy(ne(kele,1))
          z1=zz(ne(kele,1))
          x2=xx(ne(kele,2))
          y2=yy(ne(kele,2))
          z2=zz(ne(kele,2))
          x3=xx(ne(kele,3))
          y3=yy(ne(kele,3))
          z3=zz(ne(kele,3))
          x4=xx(ne(kele,4))
          y4=yy(ne(kele,4))
          z4=zz(ne(kele,4))
          vl=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
          if(vl.lt.1d-10)then
            write(*,*)"kele,vl",kele,sngl(vl)
            write(*,*)x1,y1,z1
            write(*,*)x2,y2,z2
            write(*,*)x3,y3,z3
            write(*,*)x4,y4,z4
            stop
          endif
        enddo
c
        return
        end
c---------------------------------------------------------------------------
c
c       np=1,2,3,4,5,6,7,8を除く
c
c---------------------------------------------------------------------------
        subroutine edit(np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
c
c       最後の8個の節点を最初に持って来る
c        do 500 kele=1,nele
c        do 500 j=1,4
c          kp=ne(kele,j)
c          if(kp.eq.np)ne(kele,j)=1
c          if(kp.eq.np-1)ne(kele,j)=2
c          if(kp.eq.np-2)ne(kele,j)=3
c          if(kp.eq.np-3)ne(kele,j)=4
c          if(kp.eq.np-4)ne(kele,j)=5    !
c          if(kp.eq.np-5)ne(kele,j)=6
c          if(kp.eq.np-6)ne(kele,j)=7
c          if(kp.eq.np-7)ne(kele,j)=8
c500     continue
c
c        do i=1,8
c          xx(i)=xx(np-i-1)
c          yy(i)=yy(np-i-1)
c          zz(i)=zz(np-i-1)
c        enddo
c        np=np-8
c
        nmax=-99
        nmin=99
        do 200 kele=1,nele
        do 200 j=1,4
          kp=ne(kele,j)
          if(kp.gt.nmax)nmax=kp
          if(kp.lt.nmin)nmin=kp
200     continue
        write(*,*)"nmax,nmin:",nmax,",",nmin
c
        do 100 kele=1,nele
        do 100 j=1,4
          ne(kele,j)=ne(kele,j)-8
          if(ne(kele,j).le.0)write(*,*)"0以下の節点があります:",kp-8
          if(ne(kele,j).le.0)pause
100     continue
c
        do kp=9,np
          xx(kp-8)=xx(kp)
          yy(kp-8)=yy(kp)
          zz(kp-8)=zz(kp)
        enddo
        np=np-8
c
        return
        end
c---------------------------------------------------------------------------
c
c      節点を外接円内に含む要素を探す
c      mb(na)
c---------------------------------------------------------------------------
        subroutine search(nele,ne,xx,yy,zz,xp,yp,zp,eps,nd,md,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),mb(ine),mc(ine),mc0(ine),md(ine)
        dimension xx(ind),yy(ind),zz(ind)
c
c       外接円を含む要素nel0を探す
        nel0=0
        dmin=999d0
        do kele=1,nele
          x1=xx(ne(kele,1))
          x2=xx(ne(kele,2))
          x3=xx(ne(kele,3))
          x4=xx(ne(kele,4))
          y1=yy(ne(kele,1))
          y2=yy(ne(kele,2))
          y3=yy(ne(kele,3))
          y4=yy(ne(kele,4))
          z1=zz(ne(kele,1))
          z2=zz(ne(kele,2))
          z3=zz(ne(kele,3))
          z4=zz(ne(kele,4))
c         形状関数を計算する(1次:Ni=Li)
          call caldn(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,dn)
          if(dn.lt.dmin)dmin=dn
          if(dn.lt.1d0+1d-10)nel0=kele
          if(dn.lt.1d0+1d-10)exit
        enddo
        if(nel0.eq.0)write(*,*)"nel0,dn",nel0,dmin
        if(nel0.eq.0)write(*,*)"xp",xp,yp,zp
        if(nel0.eq.0)stop
c
c       ボロノイ線図で節点を外接円内に含む要素を探す
        nd=1    !トータル:md(nd)
        md(nd)=nel0
c
        nc=1    !最初
        mc(nc)=nel0
        do while(1==1)

        nc0=0   !今回確定する要素の数
        do 100 kc=1,nc  !前回確定した要素の数:外接円内に節点を含む
          kel0=mc(kc)
c
          nb=0  !前回確定した要素の1つと隣接する要素の数
          do 200 kele=1,nele
            do kc0=1,nc0
              if(kele.eq.mc0(kc0))goto 200 !今回確定した要素
            enddo
            do kd=1,nd
              if(kele.eq.md(kd))goto 200   !前回までに確定した要素
            enddo
c
            kk=0
            do i=1,4
              i1=ne(kel0,i)
              do j=1,4
                j1=ne(kele,j)
                if(i1.eq.j1)kk=kk+1
              enddo
            enddo
            if(kk.ne.3)goto 200 !同じ節点が3つ:同じ面がある
c           外接円の確認:kele
            call certet(xx(ne(kele,1)),yy(ne(kele,1)),zz(ne(kele,1))
     &                 ,xx(ne(kele,2)),yy(ne(kele,2)),zz(ne(kele,2))
     &                 ,xx(ne(kele,3)),yy(ne(kele,3)),zz(ne(kele,3))
     &                 ,xx(ne(kele,4)),yy(ne(kele,4)),zz(ne(kele,4))
     &                 ,xc,yc,zc,r0,eps)    !r0:半径の2乗,誤差±1d-9
            rr=(xp-xc)**2+(yp-yc)**2+(zp-zc)**2  !外心-新節点の2乗
c            if(rr.gt.r0+1d-10)goto 200     !球領域を作成する時は、中心点を含む要素は全て候補になる(外接球等しい)
            if(rr.gt.r0-1d-12)goto 200
            write(*,*)"rr,r0",rr,r0
            nb=nb+1
            mb(nb)=kele
            if(nb.eq.4)exit !最大4つ
200       continue
c
          do kb=1,nb
            mc0(nc0+kb)=mb(kb)
          enddo
          nc0=nc0+nb    !今回確定した要素の数
100     continue
        if(nc0.eq.0)return   !終了
c       次考える要素
        nc=nc0
        do kc=1,nc0
          mc(kc)=mc0(kc)
        enddo
c       トータル
        do kc=1,nc0
          md(nd+kc)=mc0(kc)
        enddo
        nd=nd+nc0
c
        end do
c
        end
c
c       形状関数を計算する
       subroutine caldn(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,xp,yp,zp,dn)
        implicit double precision(a-h,o-z)
c
        v1=calvl(xp,yp,zp,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        v2=calvl(x1,y1,z1,xp,yp,zp,x3,y3,z3,x4,y4,z4)
        v3=calvl(x1,y1,z1,x2,y2,z2,xp,yp,zp,x4,y4,z4)
        v4=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,xp,yp,zp)
        vl=calvl(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4)
        if(vl.lt.1d-10)write(*,*)"caldn:vl",vl
        if(vl.lt.1d-10)stop
c
        dn1=v1/vl   !L1=N1
        dn2=v2/vl   !L2=N2
        dn3=v3/vl   !L3=N3
        dn4=v4/vl   !L4=N4
c
        dn=dn1+dn2+dn3+dn4
c
        return
        end
c---------------------------------------------------------------------------
c
c       出来上がった要素の表面の数ns(12になっているか)
c
c---------------------------------------------------------------------------
        subroutine region(np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension js1(3),js2(3),js3(3),js4(3)
        data js1/1,2,3/
        data js2/2,3,4/
        data js3/3,4,1/
        data js4/4,1,2/
c
        ns=0    !外側の面の数(4つ)
c
        do 100 kele=1,nele
c
c         面1
          is=0
          do 200 jele=1,nele
            if(jele.eq.kele)goto 200
            kk=0
            do k=1,3    !面1:1,2,3
            do j=1,4
                if(ne(kele,js1(k)).eq.ne(jele,j))kk=kk+1    !隣接する面が見つかった(keleの面1)
            enddo
            enddo
            if(kk.eq.3)is=1
            if(kk.eq.3)exit
200       continue
          if(is.eq.0)ns=ns+1    !隣接する面が見つからなかった数(外側の面の数)
c
c         面2
          is=0
          do 400 jele=1,nele
            if(jele.eq.kele)goto 400
            kk=0
            do k=1,3    !面2:2,3,4
            do j=1,4
                if(ne(kele,js2(k)).eq.ne(jele,j))kk=kk+1
            enddo
            enddo
            if(kk.eq.3)is=1
            if(kk.eq.3)exit
400       continue
          if(is.eq.0)ns=ns+1
c
c         面3
          is=0
          do 600 jele=1,nele
            if(jele.eq.kele)goto 600
            kk=0
            do k=1,3    !面3:3,4,1
            do j=1,4
                if(ne(kele,js3(k)).eq.ne(jele,j))kk=kk+1
            enddo
            enddo
            if(kk.eq.3)is=1
            if(kk.eq.3)exit
600       continue
          if(is.eq.0)ns=ns+1
c
c         面4
          is=0
          do 800 jele=1,nele
            if(jele.eq.kele)goto 800
            kk=0
            do k=1,3    !面4:4,1,2
            do j=1,4
                if(ne(kele,js4(k)).eq.ne(jele,j))kk=kk+1
            enddo
            enddo
            if(kk.eq.3)is=1
            if(kk.eq.3)exit
800       continue
          if(is.eq.0)ns=ns+1
100     continue
c
        write(*,*)"region,ns",ns
        if(ns.ne.12)stop
c
        return
        end
c---------------------------------------------------------------------------
c
c         ねじれた要素を探す
c         面の法線と線分(重心-節点)の内積が<89°の面を含む要素は含めない
c
c---------------------------------------------------------------------------
        subroutine chkne(xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3,shita)
        implicit double precision(a-h,o-z)
        data PI/3.1415926535897932384626433832795/
c
c       面 : p1,p2,p3, 節点候補 : pp(xp,yp,zp)
c       線分:重心-節点候補
        xg=(x1+x2+x3)/3d0
        yg=(y1+y2+y3)/3d0
        zg=(z1+z2+z3)/3d0
        vx=xp-xg
        vy=yp-yg
        vz=zp-zg
c       法線:dnx,dny,dnz
        v1x=x2-x1
        v1y=y2-y1
        v1z=z2-z1
        v2x=x3-x1
        v2y=y3-y1
        v2z=z3-z1
        call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz) !外積:|dn|=1
c
c       内積から角度を計算する
        shita=finpr(vx,vy,vz,dnx,dny,dnz)
        shita=PI-shita
c
        return
        end
c
c       内積から角度を計算する
        function finpr(px,py,pz,qx,qy,qz)
        implicit double precision(a-h,o-z)
        data PI/3.14159265358979323/ !84626433832795
c
        d1=dsqrt(px*px+py*py+pz*pz);
        d2=dsqrt(qx*qx+qy*qy+qz*qz);
        dd=px*qx+py*qy+pz*qz;
        d12=dd/d1/d2
c
        if(dabs(d12-(-1d0)).lt.1d-15)then       !180°:1d-16の誤差がでる
          shita=PI
        elseif(dabs(d12-1d0).lt.1d-15)then     !0°
          shita=0d0
        else
          shita=dacos(dd/d1/d2)                 !0°<shita<180°
        endif
c
        if(shita.lt.0d0 .or. shita.gt.PI)then
          write(*,*)'finpr:shita',shita
          stop
        endif
c
        finpr=shita
        return
        end
c------------------------------------------------------------------------
c
c       面に隣接する要素の番号:nec(nele,4)
c
c------------------------------------------------------------------------
        subroutine snec(nele,ne,xx,yy,zz,nec,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),nec(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension jm(4,3)
        data jm/4,4,4,1,2,3,1,3,3,1,2,2/
c
c       同じ面を探す
        do 400 iele=1,nele-1
        do 400 is=1,4   !面の番号
          i1=ne(iele,jm(is,1))  !面の節点
          i2=ne(iele,jm(is,2))
          i3=ne(iele,jm(is,3))
          do 500 jele=iele+1,nele
          do 500 js=1,4
            j1=ne(jele,jm(js,1))
            j2=ne(jele,jm(js,2))
            j3=ne(jele,jm(js,3))
            if((i1.eq.j1 .and. i2.eq.j2 .and. i3.eq.j3) .or.
     &         (i1.eq.j2 .and. i2.eq.j3 .and. i3.eq.j1) .or.
     &         (i1.eq.j3 .and. i2.eq.j1 .and. i3.eq.j2))then
              write(*,*)"snec,is,js"
              write(*,*)"iele",iele,i1,i2,i3
              write(*,*)"jele",jele,j1,j2,j3
            write(*,*)"iele",ne(iele,1),ne(iele,2),ne(iele,3),ne(iele,4)
            write(*,*)"jele",ne(jele,1),ne(jele,2),ne(jele,3),ne(jele,4)
              stop
            endif
500       continue
400     continue
c
c       隣接する要素
        do 10 i=1,nele
        do 10 j=1,4
10      nec(i,j)=0
c
        do 100 iele=1,nele-1
        do 100 is=1,4   !面の番号
          i1=ne(iele,jm(is,1))  !面の節点
          i2=ne(iele,jm(is,2))
          i3=ne(iele,jm(is,3))
          do 200 jele=iele+1,nele
          do 200 js=1,4
            j1=ne(jele,jm(js,1))
            j2=ne(jele,jm(js,2))
            j3=ne(jele,jm(js,3))
            if((i1.eq.j3 .and. i2.eq.j2 .and. i3.eq.j1) .or.
     &         (i1.eq.j1 .and. i2.eq.j3 .and. i3.eq.j2) .or.
     &         (i1.eq.j2 .and. i2.eq.j1 .and. i3.eq.j3))then
              nec(iele,is)=jele
              nec(jele,js)=iele
            endif
200       continue
100     continue
c
c       外側の面(節点:1～8, 面の数=12)
        nels=0
        do 300 kele=1,nele  !要素の番号
        do 300 ks=1,4       !面の番号
          if(nec(kele,ks).eq.0)then
            nels=nels+1
            k1=ne(kele,jm(ks,1))
            k2=ne(kele,jm(ks,2))
            k3=ne(kele,jm(ks,3))
            if(k1.ge.9 .or. k2.ge.9 .or. k3.ge.9)then
              write(*,*)"snec,kele",kele,k1,k2,k3
              stop
            endif
          endif
300     continue
        if(nels.ne.12)write(*,*)"snec,nels",nels
        if(nels.ne.12)stop
        write(*,*)"The mesh are maked correctly"
c
c       重心が他の要素内にない(要素が重なっている場合)
c        do iele=1,nele
c          i1=ne(iele,1)
c          i2=ne(iele,2)
c          i3=ne(iele,3)
c          i4=ne(iele,4)
c          xg=(xx(i1)+xx(i2)+xx(i3)+xx(i4))/4d0
c          yg=(yy(i1)+yy(i2)+yy(i3)+yy(i4))/4d0
c          zg=(zz(i1)+zz(i2)+zz(i3)+zz(i4))/4d0
c          do 600 jele=1,nele
c            if(jele.eq.iele)goto 600
c            j1=ne(jele,1)
c            j2=ne(jele,2)
c            j3=ne(jele,3)
c            j4=ne(jele,4)
cc           形状関数を計算
c
c            if(dn.gt.1d0+1d-14)then
c              write(*,*)"snec,dn",dn
c              stop
c            endif
c600       continue
c        enddo
c
        return
        end
c------------------------------------------------------------------------
c
c       埋没点を探す
c
c------------------------------------------------------------------------
       subroutine delnp(np,nele,ne,xx,yy,zz,np0,mmp,x0,y0,z0,ic,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),net(ine,4),mmp(ind),mg(ind),mm(ind)
        dimension xx(ind),yy(ind),zz(ind),x0(ind),y0(ind),z0(ind)
        dimension xt(ind),yt(ind),zt(ind)
c
c       mg(np)  0:埋没点,1:節点
c
        do 100 kp=1,np
100     mg(kp)=0
        do 200 i=1,nele
        do 200 j=1,4
200     mg(ne(i,j))=1
        do kp=1,np
          if(mg(kp).eq.0)goto 10
        enddo
        ic=0  !埋没点ない
        return
c
c       埋没点があった
c
10      ic=1
        np0=0   !埋没点 :mmp(np0),x0(mmp(np0)),y0(mmp(np0)),z0(mmp(np0))
        nt=0    !節点   :net(nelt,4),xt(nt),yt(nt),zt(nt)
        do kp=1,np
          if(mg(kp).eq.0)then  !埋没点
            np0=np0+1
            mmp(np0)=np0
            x0(np0)=xx(kp)
            y0(np0)=yy(kp)
            z0(np0)=zz(kp)
          elseif(mg(kp).eq.1)then  !節点:並べ替え
            nt=nt+1
            xt(nt)=xx(kp)
            yt(nt)=yy(kp)
            zt(nt)=zz(kp)
            mm(kp)=nt
          endif
        enddo
c
        do 400 i=1,nele
        do 400 j=1,4
400     net(i,j)=mm(ne(i,j))
c
c       元に戻す
c
        np=nt
        do kt=1,nt
          xx(kt)=xt(kt)
          yy(kt)=yt(kt)
          zz(kt)=zt(kt)
        enddo
        do 500 i=1,nele
        do 500 j=1,4
500     ne(i,j)=net(i,j)
c
        return
        end
c------------------------------------------------------------------------
c
c       領域の表面を含む要素が存在するか調べる
c
c       要素が表面をまたいでいなければよい
c       最初の表面どおりに要素の表面ができているとは限らない。
c
c------------------------------------------------------------------------
        subroutine checkms(ns0,nels,ms,xs,ys,zs
     &                     ,np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3),ne(ine,4)
        dimension mm(ind),mn(ine)
        dimension xx(ind),yy(ind),zz(ind),xs(ind),ys(ind),zs(ind)
        dimension ne0(ine,4),net(ine,4)
        data PI/3.1415926535897932384626433832795/
c       ========================================================================
c
c       表面の節点を座標から探す
c
c       ========================================================================
        nn=0
        do ks=1,ns0
          do kp=1,np
            if(dabs(xs(ks)-xx(kp)).lt.1d-14 .and.
     &         dabs(ys(ks)-yy(kp)).lt.1d-14 .and.
     &         dabs(zs(ks)-zz(kp)).lt.1d-14 )then
              nn=nn+1
              mm(ks)=kp
              exit
            endif
          enddo
        enddo
c       ns0があっていることを確認する
        write(*,*)"ns0,nn",ns0,nn
        pause
c       ========================================================================
c
c       表面の節点を含む要素を描かせる。
c       ip(4) -1:外側の節点,0:表面の節点,1:内側の節点
c       ========================================================================
        nel0=0
        do kele=1,nele
          kout=0    !内側にある要素内の節点の数(0～4)
          kin=0     !外側にある要素内の節点の数(0～4)
          do 200 j=1,4
            k=ne(kele,j)
c           -----------------------------------------------------------
c           表面上の節点かどうか判断する
c           -----------------------------------------------------------
            do ks=1,ns0
              if(k.eq.mm(ks))goto 200  !表面の節点は飛ばす
            enddo
c           -----------------------------------------------------------
c           節点が表面の内側にあるか外側にあるか調べる
c           -----------------------------------------------------------
            dmin=999d0
            mels=0
            do kels=1,nels  !最寄の面を探す
              xg=(xs(ms(kels,1))+xs(ms(kels,2))+xs(ms(kels,3)))/3d0
              yg=(ys(ms(kels,1))+ys(ms(kels,2))+ys(ms(kels,3)))/3d0
              zg=(zs(ms(kels,1))+zs(ms(kels,2))+zs(ms(kels,3)))/3d0
              dd=dsqrt((xg-xx(k))**2+(yg-yy(k))**2+(zg-zz(k))**2)
              if(dd.lt.dmin)mels=kels
              if(dd.lt.dmin)dmin=dd
            enddo
            if(mels.eq.0)stop
c
c           表面の法線ベクトルを計算する
            v1x=xs(ms(mels,2))-xs(ms(mels,1))
            v1y=ys(ms(mels,2))-ys(ms(mels,1))
            v1z=zs(ms(mels,2))-zs(ms(mels,1))
            v2x=xs(ms(mels,3))-xs(ms(mels,1))
            v2y=ys(ms(mels,3))-ys(ms(mels,1))
            v2z=zs(ms(mels,3))-zs(ms(mels,1))
            call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
c
c           ベクトル(表面の重心→節点k)を計算する
            xg=(xs(ms(mels,1))+xs(ms(mels,2))+xs(ms(mels,3)))/3d0
            yg=(ys(ms(mels,1))+ys(ms(mels,2))+ys(ms(mels,3)))/3d0
            zg=(zs(ms(mels,1))+zs(ms(mels,2))+zs(ms(mels,3)))/3d0
            vx=xx(k)-xg
            vy=yy(k)-yg
            vz=zz(k)-zg
c
c           内積から角度を計算する(表面の法線ベクトルは内向き)
            shita=finpr(vx,vy,vz,dnx,dny,dnz)
            shita=PI-shita
            if(shita.gt.90d0/180d0*PI)kout=kout+1   !外側にある節点
            if(shita.le.90d0/180d0*PI)kin=kin+1     !内側にある節点
200       continue
          if(kin.ge.1 .and. kout.ge.1)nel0=nel0+1
          if(kin.ge.1 .and. kout.ge.1)ne0(nel0,1)=ne(kele,1)
          if(kin.ge.1 .and. kout.ge.1)ne0(nel0,2)=ne(kele,2)
          if(kin.ge.1 .and. kout.ge.1)ne0(nel0,3)=ne(kele,3)
          if(kin.ge.1 .and. kout.ge.1)ne0(nel0,4)=ne(kele,4)
        enddo
c
c       表面の節点がばらついていて、間に内部の節点が挟まっていないか確認する(表面の節点の様に振舞う)
        write(*,*)"表面をまたぐ要素の数 nel0=",nel0
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nel0,ne0,2d0,2,ine)
        pause
c       ========================================================================
c
c       外側の要素を削除する(表面を持つ要素の重心から領域の内外を判定)
c       重心に最も近い表面:ms(nels,3),xs(ns0),ys(ns0),zs(ns0)
c
c       ========================================================================
        nel0=0  !削除する要素の数
        do 100 kele=1,nele      !表面を含む要素の数
c
c         np=1～8を持つ要素を削除する
          do j=1,4
            kp=ne(kele,j)
            if(kp.eq.1 .or. kp.eq.2 .or. kp.eq.3 .or. kp.eq.4 .or.
     &         kp.eq.5 .or. kp.eq.6 .or. kp.eq.7 .or. kp.eq.8)then
              nel0=nel0+1
              mn(nel0)=kele
              goto 100
            endif
          enddo
c
c         外側の要素を削除する
          k1=ne(kele,1)
          k2=ne(kele,2)
          k3=ne(kele,3)
          k4=ne(kele,4)
          k=0
          do ks=1,ns0
            if(k1.eq.mm(ks))k=k+1
            if(k2.eq.mm(ks))k=k+1
            if(k3.eq.mm(ks))k=k+1
            if(k4.eq.mm(ks))k=k+1
          enddo
c          if(k.ne.3)goto 100   !表面を持たない要素は飛ばす
          if(k.le.2)goto 100   !表面を持たない要素は飛ばす
c
c         最寄の表面との角度を計算する
          xg=(xx(k1)+xx(k2)+xx(k3)+xx(k4))/4d0
          yg=(yy(k1)+yy(k2)+yy(k3)+yy(k4))/4d0
          zg=(zz(k1)+zz(k2)+zz(k3)+zz(k4))/4d0
          dmin=999d0
          mels=0
          do kels=1,nels
            m1=ms(kels,1)
            m2=ms(kels,2)
            m3=ms(kels,3)
            xgs=(xs(m1)+xs(m2)+xs(m3))/3d0
            ygs=(ys(m1)+ys(m2)+ys(m3))/3d0
            zgs=(zs(m1)+zs(m2)+zs(m3))/3d0
            dd=dsqrt((xg-xgs)**2+(yg-ygs)**2+(zg-zgs)**2)
            if(dd.lt.dmin)mels=kels
            if(dd.lt.dmin)dmin=dd
          enddo
c
c         表面の法線ベクトルを計算する
          v1x=xs(ms(mels,2))-xs(ms(mels,1))
          v1y=ys(ms(mels,2))-ys(ms(mels,1))
          v1z=zs(ms(mels,2))-zs(ms(mels,1))
          v2x=xs(ms(mels,3))-xs(ms(mels,1))
          v2y=ys(ms(mels,3))-ys(ms(mels,1))
          v2z=zs(ms(mels,3))-zs(ms(mels,1))
          call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
c
c         ベクトル(表面の重心→要素の重心)を計算する
          xgs=(xs(ms(mels,1))+xs(ms(mels,2))+xs(ms(mels,3)))/3d0
          ygs=(ys(ms(mels,1))+ys(ms(mels,2))+ys(ms(mels,3)))/3d0
          zgs=(zs(ms(mels,1))+zs(ms(mels,2))+zs(ms(mels,3)))/3d0
          vx=xg-xgs
          vy=yg-ygs
          vz=zg-zgs
c
c         内積から角度を計算する(表面の法線ベクトルは内向き)
          shita=finpr(vx,vy,vz,dnx,dny,dnz)
          shita=PI-shita
          if(shita.gt.90d0/180d0*PI)nel0=nel0+1   !外側にある節点
          if(shita.gt.90d0/180d0*PI)mn(nel0)=kele !外側にある節点
100     continue
        write(*,*)"削除する要素の数:",nel0
        pause
c       ========================================================================
c
c       並べなおし
c
c       ========================================================================
        nelt=0
        do 300 kele=1,nele
          do kel0=1,nel0
            if(kele.eq.mn(kel0))goto 300   !削除する面を飛ばす
          enddo
          nelt=nelt+1
          do 400 j=1,4
400       net(nelt,j)=ne(kele,j)
300     continue
c
c       最後の8個の節点を最初に持って来る
        do 500 kelt=1,nelt
        do 500 j=1,4
          kp=net(kelt,j)
          if(kp.eq.np)net(kelt,j)=1
          if(kp.eq.np-1)net(kelt,j)=2
          if(kp.eq.np-2)net(kelt,j)=3
          if(kp.eq.np-3)net(kelt,j)=4
          if(kp.eq.np-4)net(kelt,j)=5
          if(kp.eq.np-5)net(kelt,j)=6
          if(kp.eq.np-6)net(kelt,j)=7
          if(kp.eq.np-7)net(kelt,j)=8
500     continue
c
        do i=1,8
          xx(i)=xx(np-i-1)
          yy(i)=yy(np-i-1)
          zz(i)=zz(np-i-1)
        enddo
        np=np-8
c
        nele=nelt
        do 600 kele=1,nele
        do 600 j=1,4
600     ne(kele,j)=net(kele,j)
c
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nele,ne,2d0,2,ine)
        pause
c
        return
        end





c------------------------------------------------------------------------
c
c       領域の表面を含む要素が存在するか調べる
c
c       要素が表面をまたいでいなければよい
c       最初の表面どおりに要素の表面ができているとは限らない。
c
c------------------------------------------------------------------------
        subroutine checkms2(ns0,nels,ms,xs,ys,zs
     &                     ,np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3),ne(ine,4)
        dimension mm(ind),mn(ine),m0(ind)
        dimension xx(ind),yy(ind),zz(ind),xs(ind),ys(ind),zs(ind)
        dimension ne0(ine,4),net(ine,4)
        dimension jm(4,3)
        dimension nelb(ine),neb(ine)
        dimension le(ine,3),le0(ine,3)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension mc(ine,3)    !表示用
        dimension mk(ind)
        data PI/3.1415926535897932384626433832795/
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
c       ========================================================================
c
c       表示
c
c       ========================================================================
        call pltrset()
        call pltsetnp(ns0,xs,ys,zs,ind)
c       ========================================================================
c
c       表面の節点を座標から探す
c
c       ========================================================================
        nn=0
        do ks=1,ns0
          do kp=1,np
            if(dabs(xs(ks)-xx(kp)).lt.1d-14 .and.
     &         dabs(ys(ks)-yy(kp)).lt.1d-14 .and.
     &         dabs(zs(ks)-zz(kp)).lt.1d-14 )then
              nn=nn+1
              mm(ks)=kp
              m0(kp)=ks
              exit
            endif
          enddo
        enddo
c       ns0があっていることを確認する
        if(ns0.ne.nn)write(*,*)"ns0,nn",ns0,nn
        if(ns0.ne.nn)stop
c       ========================================================================
c
c       表面ms(nels,3)の作成
c
c       ========================================================================
        call initsuf(ns0,xs,ys,zs,nels,ms)
        pause
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
c       ========================================================================
c
c       最初の表面を見つける
c
c       ========================================================================
        do kels=1,nels
          m1=mm(ms(kels,3))     !外側の要素を探すので反対向きの面
          m2=mm(ms(kels,2))
          m3=mm(ms(kels,1))
          do kele=1,nele
            !面1
            k1=ne(kele,jm1(1))
            k2=ne(kele,jm1(2))
            k3=ne(kele,jm1(3))
            if((m1.eq.k1 .and. m2.eq.k2 .and. m3.eq.k3) .or.
     &         (m1.eq.k2 .and. m2.eq.k3 .and. m3.eq.k1) .or.
     &         (m1.eq.k3 .and. m2.eq.k1 .and. m3.eq.k2))then
              nb=1              !要素(面)の数
              nelb(nb)=kele     !要素の番号
              neb(nb)=1         !要素の面の番号(1～4)
              goto 10
            endif
            !面2
            k1=ne(kele,jm2(1))
            k2=ne(kele,jm2(2))
            k3=ne(kele,jm2(3))
            if((m1.eq.k1 .and. m2.eq.k2 .and. m3.eq.k3) .or.
     &         (m1.eq.k2 .and. m2.eq.k3 .and. m3.eq.k1) .or.
     &         (m1.eq.k3 .and. m2.eq.k1 .and. m3.eq.k2))then
              nb=1
              nelb(nb)=kele
              neb(nb)=2
              goto 10
            endif
            !面3
            k1=ne(kele,jm3(1))
            k2=ne(kele,jm3(2))
            k3=ne(kele,jm3(3))
            if((m1.eq.k1 .and. m2.eq.k2 .and. m3.eq.k3) .or.
     &         (m1.eq.k2 .and. m2.eq.k3 .and. m3.eq.k1) .or.
     &         (m1.eq.k3 .and. m2.eq.k1 .and. m3.eq.k2))then
              nb=1
              nelb(nb)=kele
              neb(nb)=3
              goto 10
            endif
            !面4
            k1=ne(kele,jm4(1))
            k2=ne(kele,jm4(2))
            k3=ne(kele,jm4(3))
            if((m1.eq.k1 .and. m2.eq.k2 .and. m3.eq.k3) .or.
     &         (m1.eq.k2 .and. m2.eq.k3 .and. m3.eq.k1) .or.
     &         (m1.eq.k3 .and. m2.eq.k1 .and. m3.eq.k2))then
              nb=1
              nelb(nb)=kele
              neb(nb)=4
              goto 10
            endif
          enddo
        enddo
10     continue
c       ========================================================================
c
c       面を探す
c       nb          :要素(面)の数
c       nelb(nb)    :表面を持つ要素の番号
c       neb(nb)     :表面部分の要素の面の番号
c       ========================================================================
c
c       最初の面
        n1=ne(nelb(nb),jm(neb(nb),1))
        n2=ne(nelb(nb),jm(neb(nb),2))
        n3=ne(nelb(nb),jm(neb(nb),3))
c
c       面を探す辺
        ie=3
        le(1,1)=n2  !要素を作る辺の節点1
        le(1,2)=n1  !要素を作る辺の節点2
        le(1,3)=n3  !元の要素の節点(新しい要素との角度計算に使う)

        le(2,1)=n3  !要素を作る辺の節点1
        le(2,2)=n2  !要素を作る辺の節点2
        le(2,3)=n1  !元の要素の節点(新しい要素との角度計算に使う)

        le(3,1)=n1  !要素を作る辺の節点1
        le(3,2)=n3  !要素を作る辺の節点2
        le(3,3)=n2  !元の要素の節点(新しい要素との角度計算に使う)
c
c       面の法線は内向き

        nelc=0  !表示用
100     continue

        ie0=0   !今回登録した辺の数
        do 300 ke=1,ie
          k1=le(ke,1)
          k2=le(ke,2)
          kk=le(ke,3)   !反対側の節点
          call pltnt(xx(k1),yy(k1),zz(k1),10d0,1)  !1:白, 2:赤, 3:緑, 4:青, 5:黄色
c
c         既に登録した要素(面)
          do kb=1,nb
            m1=ne(nelb(kb),jm(neb(kb),1))
            m2=ne(nelb(kb),jm(neb(kb),2))
            m3=ne(nelb(kb),jm(neb(kb),3))
            if(k1.eq.m1 .and. k2.eq.m2)goto 300
            if(k1.eq.m2 .and. k2.eq.m3)goto 300
            if(k1.eq.m3 .and. k2.eq.m1)goto 300
          enddo
c
c         表面を探す
          shita=999d0/180d0*PI
          kelb=0
          call surns(m0(k1),m0(k2),m0(kk),nels,ms,nk,mk,ind,ine)  !周囲の節点
          do 400 kele=1,nele
          do 400 j=1,4  !要素の面の番号(1～4)
            m1=ne(kele,jm(j,1))
            m2=ne(kele,jm(j,2))
            m3=ne(kele,jm(j,3))
            do 1000 ks=1,ns0
              k=mm(ks)
c            do 1000 ks=1,nk
c              k=mm(mk(ks))
              if(k.eq.kk)goto 1000
              if((k1.eq.m1 .and. k2.eq.m2 .and. k.eq.m3).or.
     &           (k1.eq.m2 .and. k2.eq.m3 .and. k.eq.m1).or.
     &           (k1.eq.m3 .and. k2.eq.m1 .and. k.eq.m2))then
c                既にある面
c                do kb=1,nb
c                  m1=ne(nelb(kb),jm(neb(kb),1))
c                  m2= ne(nelb(kb),jm(neb(kb),2))
c                  m3=ne(nelb(kb),jm(neb(kb),3))
c                  if(k3.eq.m1 .and. k2.eq.m2)goto 1000
c                  if(k3.eq.m2 .and. k2.eq.m3)goto 1000
c                  if(k3.eq.m3 .and. k2.eq.m1)goto 1000
c                  if(k1.eq.m1 .and. k3.eq.m2)goto 1000
c                  if(k1.eq.m2 .and. k3.eq.m3)goto 1000
c                  if(k1.eq.m3 .and. k3.eq.m1)goto 1000
c                enddo
c                角度
c                sht=fshita(k1,k2,k,kk,xx,yy,zz,ind)
c                write(*,*)"k:",k1,k2,k
c                write(*,*)"sht=",sht/PI*180d0
c                if(sht.lt.shita)kelb=kele       !外側にめくれる側の面を探す
c                if(sht.lt.shita)keb=j
c                if(sht.lt.shita)k3=k
c                if(sht.lt.shita)shita=sht
c                外側にない面
                if(ifns(kele,j,ne,m1,m2,m3,kk,nels,ns0,ms
     &                         ,mm,m0,ind,ine).eq.0)goto 1000
                kelb=kele
                keb=j

              endif
1000        continue
400       continue
          if(kelb.eq.0)write(*,*)"表面がありません"
          if(kelb.eq.0)write(*,*)"k1,k2:",k1,k2
          if(kelb.eq.0)write(*,*)"kk,shita:",kk,shita/PI*180d0
          if(kelb.eq.0)stop
          write(*,*)"shita=",shita/PI*180d0
c
c         要素の面(表面部分)を登録
          nb=nb+1
          nelb(nb)=kelb
          neb(nb)=keb
c         次の辺を登録
          ie0=ie0+1
          le0(ie0,1)=k3 !辺の節点1
          le0(ie0,2)=k2 !辺の節点2
          le0(ie0,3)=k1 !辺以外の節点(角度の計算に使う)
          ie0=ie0+1
          le0(ie0,1)=k1 !辺の節点1
          le0(ie0,2)=k3 !辺の節点2
          le0(ie0,3)=k2 !辺以外の節点(角度の計算に使う)
c
          nelc=nelc+1
          mc(nelc,1)=ne(nelb(nb),jm(neb(nb),1))
          mc(nelc,2)=ne(nelb(nb),jm(neb(nb),2))
          mc(nelc,3)=ne(nelb(nb),jm(neb(nb),3))
          write(*,*)"nelc=",nelc
          do kelc=1,nelc
          write(*,*)"kelc:",mc(kelc,1),mc(kelc,2),mc(kelc,3)
          enddo
          call pltsetme(nelc,mc,3d0,2,ind,ine)
          pause
300     continue
        if(ie0.eq.0)return     !辺がない
c
        ie=ie0
        do 500 ke=1,ie0
        do 500 j=1,3
500     le(ke,j)=le0(ke,j)
c
        goto 100
c
        return
        end
c------------------------------------------------------------------------
c
c       面の角度を計算する
c
c------------------------------------------------------------------------
        function fshita(k1,k2,k3,kk,xx,yy,zz,ind)
        implicit double precision(a-h,o-z)
        dimension xx(ind),yy(ind),zz(ind)
        data PI/3.1415926535897932384626433832795/
c
c       垂線を計算する
        dk1=(xx(kk)-xx(k1))*(xx(k2)-xx(k1))
     &     +(yy(kk)-yy(k1))*(yy(k2)-yy(k1))
     &     +(zz(kk)-zz(k1))*(zz(k2)-zz(k1))
        dk2=(xx(k2)-xx(k1))**2+(yy(k2)-yy(k1))**2+(zz(k2)-zz(k1))**2
        dk=dk1/dk2
        xk=xx(k1)+dk*(xx(k2)-xx(k1))
        yk=yy(k1)+dk*(yy(k2)-yy(k1))
        zk=zz(k1)+dk*(zz(k2)-zz(k1))
c
c       基準のベクトル
        vx=xx(kk)-xk
        vy=yy(kk)-yk
        vz=zz(kk)-zk


c       -------------------------------------

c
c       垂線を計算する
        dk1=(xx(k3)-xx(k1))*(xx(k2)-xx(k1))
     &     +(yy(k3)-yy(k1))*(yy(k2)-yy(k1))
     &     +(zz(k3)-zz(k1))*(zz(k2)-zz(k1))
        dk2=(xx(k2)-xx(k1))**2+(yy(k2)-yy(k1))**2+(zz(k2)-zz(k1))**2
        dk=dk1/dk2
        xk=xx(k1)+dk*(xx(k2)-xx(k1))
        yk=yy(k1)+dk*(yy(k2)-yy(k1))
        zk=zz(k1)+dk*(zz(k2)-zz(k1))
c
c       接線ベクトル
        vtx=xx(k3)-xk
        vty=yy(k3)-yk
        vtz=zz(k3)-zk
c
c       法線ベクトル
        v1x=xx(k2)-xx(k1)
        v1y=yy(k2)-yy(k1)
        v1z=zz(k2)-zz(k1)
        v2x=xx(k3)-xx(k1)
        v2y=yy(k3)-yy(k1)
        v2z=zz(k3)-zz(k1)
        call cprod(v1x,v1y,v1z,v2x,v2y,v2z,dnx,dny,dnz)
c
c       内積から角度を計算する
        phi=finpr(vx,vy,vz,dnx,dny,dnz) !法線との角度
        phi=PI-phi
        psi=finpr(vx,vy,vz,vtx,vty,vtz) !接線との角度
        if(phi.ge.0d0    .and. phi.le.PI/2d0)shita=psi
        if(phi.ge.PI/2d0 .and. phi.le.PI    )shita=2d0*PI-psi
c
        if(shita.lt.0d0 .or. shita.gt.2d0*PI)then
          write(*,*)"fshita=",shita/PI*180d0
          stop
        endif
c
        fshita=shita
        return
        end
c------------------------------------------------------------------------
c
c       周囲の節点
c
c------------------------------------------------------------------------
        subroutine surns(ks1,ks2,ks3,nels,ms,nk,mk,ind,ine)
        implicit double precision(a-h,o-z)
        dimension mk(ind),ms(ine,3),mk0(ind)
c
        nk=0
        do 10 kels=1,nels
          m1=ms(kels,1)
          m2=ms(kels,2)
          m3=ms(kels,3)
          if((m1.ne.ks1 .and. m2.ne.ks1 .and. m3.ne.ks1) .and.
     &       (m1.ne.ks2 .and. m2.ne.ks2 .and. m3.ne.ks2) )goto 10
c         点1
          if(m1.eq.ks1 .or. m1.eq.ks2 .or. m1.eq.ks3)goto 20
          do kk=1,nk
            if(m1.eq.mk(kk))goto 20
          enddo
          nk=nk+1
          mk(nk)=m1
20        continue
c         点2
          if(m2.eq.ks1 .or. m2.eq.ks2 .or. m2.eq.ks3)goto 30
          do kk=1,nk
            if(m2.eq.mk(kk))goto 30
          enddo
          nk=nk+1
          mk(nk)=m2
30        continue
c         点3
          if(m3.eq.ks1 .or. m3.eq.ks2 .or. m3.eq.ks3)goto 40
          do kk=1,nk
            if(m3.eq.mk(kk))goto 40
          enddo
          nk=nk+1
          mk(nk)=m3
40        continue
10      continue
c
c
        write(*,*)"nk=",nk
        nk0=nk
        do 90 i=1,nk
90      mk0(i)=mk(i)
c
        do 50 kels=1,nels
          m1=ms(kels,1)
          m2=ms(kels,2)
          m3=ms(kels,3)
          ifg=0
          do i=1,nk0
            if(m1.eq.mk0(i))ifg=1
            if(m2.eq.mk0(i))ifg=1
            if(m3.eq.mk0(i))ifg=1
          enddo
          if(m1.eq.ks1 .or. m2.eq.ks1 .or. m3.eq.ks1)ifg=1
          if(m1.eq.ks2 .or. m2.eq.ks2 .or. m3.eq.ks2)ifg=1
          if(ifg.eq.0)goto 50

c         点1
          if(m1.eq.ks1 .or. m1.eq.ks2 .or. m1.eq.ks3)goto 60
          do kk=1,nk
            if(m1.eq.mk(kk))goto 60
          enddo
          nk=nk+1
          mk(nk)=m1
60        continue
c         点2
          if(m2.eq.ks1 .or. m2.eq.ks2 .or. m2.eq.ks3)goto 70
          do kk=1,nk
            if(m2.eq.mk(kk))goto 70
          enddo
          nk=nk+1
          mk(nk)=m2
70        continue
c         点3
          if(m3.eq.ks1 .or. m3.eq.ks2 .or. m3.eq.ks3)goto 80
          do kk=1,nk
            if(m3.eq.mk(kk))goto 80
          enddo
          nk=nk+1
          mk(nk)=m3
80        continue
50      continue
c
        write(*,*)"nk=",nk
c
        return
        end
c------------------------------------------------------------------------
c
c       外側の面
c
c------------------------------------------------------------------------
        function ifns(kele,j,ne,m1,m2,m3,kk,nels,ns0,ms,mm,m0,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),ms(ine,3)
        dimension mm(ind),m0(ind),mk(ind)
        dimension jm(4,3)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
c
        n1=mm(m1)
        n2=mm(m2)
        n3=mm(m3)
        if(j.eq.1)n4=ne(kele,1)
        if(j.eq.2)n4=ne(kele,2)
        if(j.eq.3)n4=ne(kele,3)
        if(j.eq.4)n4=ne(kele,4)
c
        call surns(m0(k1),m0(k2),m0(kk),nels,ms,nk,mk,ind,ine)  !周囲の節点
c
        ifns=0
        do kk=1,nk
          if( n4.eq.mm(mk(kk)) )ifns=1
        enddo
c
        return
        end















c------------------------------------------------------------------------
c
c       領域の表面を含む要素が存在するか調べる
c
c       要素が表面をまたいでいなければよい
c       最初の表面どおりに要素の表面ができているとは限らない。
c
c------------------------------------------------------------------------









        subroutine checkms3(ns,nels,ms,xs,ys,zs
     &                     ,np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ms(ine,3),ne(ine,4)
        dimension mm(ind),m0(ind)
        dimension xx(ind),yy(ind),zz(ind),xs(ind),ys(ind),zs(ind)
        dimension nelb(ine)
        dimension mc(ind)    !表示用
        dimension mk1(ind),mk2(ind),mk3(ind),mk4(ind)
        dimension ne0(ine,4)
c
c       ========================================================================
c
c       表示
c
c       ========================================================================
        call pltrset()
        call pltsetnp(ns,xs,ys,zs,ind)
c       ========================================================================
c
c       外側の節点
c
c       ========================================================================
        nc=8
        do 100 kc=1,8
100     mc(kc)=kc
c       ========================================================================
c
c       表面の節点を座標から探す
c
c       ========================================================================
        nn=0
        do ks=1,ns
          do kp=1,np
            if(dabs(xs(ks)-xx(kp)).lt.1d-14 .and.
     &         dabs(ys(ks)-yy(kp)).lt.1d-14 .and.
     &         dabs(zs(ks)-zz(kp)).lt.1d-14 )then
              nn=nn+1
              mm(ks)=kp
              m0(kp)=ks
              exit
            endif
          enddo
        enddo
c       nsがあっていることを確認する
        if(ns.ne.nn)write(*,*)"ns,nn",ns,nn
        if(ns.ne.nn)stop
c       ========================================================================
c
c       表面ms(nels,3)の作成
c
c       ========================================================================
        call initsuf(ns,xs,ys,zs,nels,ms)
        write(*,*)"np,ns:",np,ns
        pause
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
c       ========================================================================
c
c       最初の表面を見つける
c
c       ========================================================================
        nb=0    !削除する要素の数
        do 200 kele=1,nele
          n1=ne(kele,1)
          n2=ne(kele,2)
          n3=ne(kele,3)
          n4=ne(kele,4)
c
c         外側の節点を含む要素ははずす
          do kc=1,nc
            if(n1.eq.mc(kc) .or. n2.eq.mc(kc)
     &                      .or. n3.eq.mc(kc) .or. n4.eq.mc(kc))then
              nb=nb+1
              nelb(nb)=kele
              goto 200
            endif
          enddo
c
c         外側付近の要素をはずす
          n=0
          do ks=1,ns
            if(n1.eq.mm(ks))n=n+1
            if(n2.eq.mm(ks))n=n+1
            if(n3.eq.mm(ks))n=n+1
            if(n4.eq.mm(ks))n=n+1
          enddo
          if(n.ne.4)goto 200
c          write(*,*)"kele=",kele,"/",nele
          call surnp(m0(n1),nk1,mk1,nels,ms,ind,ine)
          call surnp(m0(n2),nk2,mk2,nels,ms,ind,ine)
          call surnp(m0(n3),nk3,mk3,nels,ms,ind,ine)
          call surnp(m0(n4),nk4,mk4,nels,ms,ind,ine)
c          write(*,*)"nk1=",nk1
c          write(*,*)"nk2=",nk2
c          write(*,*)"nk3=",nk3
c          write(*,*)"nk4=",nk4
c          pause
c
c         節点1
          iflag=0
          do k=1,nk2
            if(m0(n1).eq.mk2(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk3
            if(m0(n1).eq.mk3(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk4
            if(m0(n1).eq.mk4(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
c
c         節点2
          iflag=0
          do k=1,nk3
            if(m0(n2).eq.mk3(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk4
            if(m0(n2).eq.mk4(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk1
            if(m0(n2).eq.mk1(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
c
c         節点3
          iflag=0
          do k=1,nk4
            if(m0(n3).eq.mk4(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk1
            if(m0(n3).eq.mk1(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk2
            if(m0(n3).eq.mk2(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
c
c         節点4
          iflag=0
          do k=1,nk1
            if(m0(n4).eq.mk1(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk2
            if(m0(n4).eq.mk2(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
          iflag=0
          do k=1,nk3
            if(m0(n4).eq.mk3(k))iflag=1
          enddo
          if(iflag.eq.0)goto 300
c
          goto 200
300       continue
          nb=nb+1
          nelb(nb)=kele
200     continue
c
        write(*,*)"nb=",nb,"/",nele
        nel0=0
        do 400 kele=1,nele
          do kb=1,nb
            if(kele.eq.nelb(kb))goto 400
          enddo
          nel0=nel0+1
          do 500 j=1,4
500       ne0(nel0,j)=ne(kele,j)
400     continue
c
        nele=nel0
        do 600 kel0=1,nel0
        do 600 j=1,4
600     ne(kel0,j)=ne0(kel0,j)
c
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nele,ne,2d0,2,ine)
c        write(*,*)"end checkms3"

        do kele=1,nele
        do j=1,4
          k=ne(kele,j)
          if(ne(kele,j).le.8)write(*,*)"外側の節点があります:",k
          if(ne(kele,j).le.8)pause
        enddo
        enddo
        write(*,*)"外側の節点はありません"
c
        return
        end
c

        subroutine surnp(n0,nk,mk,nels,ms,ind,ine)
        implicit double precision(a-h,o-z)
        dimension mk(ine),ms(ine,3)
c
        nk=0
        do 100 kels=1,nels
          m1=ms(kels,1)
          m2=ms(kels,2)
          m3=ms(kels,3)
          if(n0.eq.m1)then
            do k=1,nk
              if(mk(k).eq.m2)goto 200
            enddo
            nk=nk+1
            mk(nk)=m2
200         continue
            do k=1,nk
              if(mk(k).eq.m3)goto 300
            enddo
            nk=nk+1
            mk(nk)=m3
300         continue
c
          elseif(n0.eq.m2)then
            do k=1,nk
              if(mk(k).eq.m3)goto 400
            enddo
            nk=nk+1
            mk(nk)=m3
400         continue
            do k=1,nk
              if(mk(k).eq.m1)goto 500
            enddo
            nk=nk+1
            mk(nk)=m1
500         continue
c
          elseif(n0.eq.m3)then
            do k=1,nk
              if(mk(k).eq.m1)goto 600
            enddo
            nk=nk+1
            mk(nk)=m1
600         continue
            do k=1,nk
              if(mk(k).eq.m2)goto 700
            enddo
            nk=nk+1
            mk(nk)=m2
700         continue
          endif
100     continue
c
        return
        end
c------------------------------------------------------------------------
c
c       要素形状の確認をする
c
c------------------------------------------------------------------------
        function elshape(np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4)
        dimension xx(ind),yy(ind),zz(ind),r(4)
        dimension ne0(ine,4)
c
        write(*,*)"in elshape np,nele:",np,nele
        drmax=0d0
        do 100 i=1,nele
          xg=(xx(ne(i,1))+xx(ne(i,2))+xx(ne(i,3))+xx(ne(i,4)))/4d0
          yg=(yy(ne(i,1))+yy(ne(i,2))+yy(ne(i,3))+yy(ne(i,4)))/4d0
          zg=(zz(ne(i,1))+zz(ne(i,2))+zz(ne(i,3))+zz(ne(i,4)))/4d0
          do 200 j=1,4
            x=xx(ne(i,j))
            y=yy(ne(i,j))
            z=zz(ne(i,j))
            r(j)=dsqrt((x-xg)**2+(y-yg)**2+(z-zg)**2)
200       continue
          a=r(1)+r(2)+r(3)+r(4)
          b=r(1)**2+r(2)**2+r(3)**2+r(4)**2
          dr=2d0*dsqrt(4d0*b-a**2)/a
          if(drmax.lt.dr)kele=i
          if(drmax.lt.dr)drmax=dr
100     continue
        write(*,*)"elshape:",kele,sngl(drmax)
        elshape=drmax
        write(10,*)'elshape:',elshape
c
        nel0=1
        ne0(1,1)=ne(kele,1)
        ne0(1,2)=ne(kele,2)
        ne0(1,3)=ne(kele,3)
        ne0(1,4)=ne(kele,4)
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nel0,ne0,2d0,2,ine)
c
        return
        end
c------------------------------------------------------------------------
c
c       要素形状の確認をする
c
c------------------------------------------------------------------------
        function elshape2(np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),ne0(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension x0(4),y0(4),z0(4)
c
c       形状係数計算
        drmax=0d0
        do 100 i=1,nele
          do j=1,4
            x0(j)=xx(ne(i,j))
            y0(j)=yy(ne(i,j))
            z0(j)=zz(ne(i,j))
          enddo
c          dr=difdr4d(x0,y0,z0,dr0,rr0)
          call difdr4d(x0,y0,z0,dr,dr1,dr2)
          if(drmax.lt.dr)kele=i
          if(drmax.lt.dr)dr10=dr1
          if(drmax.lt.dr)dr20=dr2
          if(drmax.lt.dr)drmax=dr
c
          vl=calvl2(x0(1),y0(1),z0(1),x0(2),y0(2),z0(2)
     &             ,x0(3),y0(3),z0(3),x0(4),y0(4),z0(4))
          if(vl.lt.1d-7)then
            kele=i
            dr10=1d0
            dr20=1d0
            drmax=1d0
          endif
100     continue
        write(*,*)"elshape2:"
     &             ,kele,sngl(drmax),"=",sngl(dr10),"/",sngl(dr20)
        elshape2=drmax
        write(10,*)'elshape2:',elshape2
c
c       描画
        nel0=1
        ne0(1,1)=ne(kele,1)
        ne0(1,2)=ne(kele,2)
        ne0(1,3)=ne(kele,3)
        ne0(1,4)=ne(kele,4)
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nel0,ne0,2d0,2,ine)
c
        return
        end
c------------------------------------------------------------------------
c
c       表面を探す
c
c------------------------------------------------------------------------
        subroutine findms(nele,ne,nels,ms)
        include "header.h"
        dimension ne(ine,4),ms(ine,3)
        dimension nec(ine,4)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        dimension jm(4,3)
        dimension mm(ind)
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
        nels=0    !表面の面の数
        do 300 i=1,nele
        do 300 j=1,4
          if(nec(i,j).ne.0)goto 300
          nels=nels+1
          ms(nels,1)=ne(i,jm(j,1))
          ms(nels,2)=ne(i,jm(j,2))
          ms(nels,3)=ne(i,jm(j,3))
300     continue
c
        return
        end
c------------------------------------------------------------------------
c
c       面形状の確認をする
c
c------------------------------------------------------------------------
        function sfshape(nels,ms,xx,yy,zz)
        include "header.h"
        dimension ms(ine,3)
        dimension xx(ind),yy(ind),zz(ind),r(4)
c
        drmax=0d0
        do 100 i=1,nels
          xg=(xx(ms(i,1))+xx(ms(i,2))+xx(ms(i,3)))/3d0
          yg=(yy(ms(i,1))+yy(ms(i,2))+yy(ms(i,3)))/3d0
          zg=(zz(ms(i,1))+zz(ms(i,2))+zz(ms(i,3)))/3d0
          do 200 j=1,3
            x=xx(ms(i,j))
            y=yy(ms(i,j))
            z=zz(ms(i,j))
            r(j)=dsqrt((x-xg)**2+(y-yg)**2+(z-zg)**2)
200       continue
          a=r(1)+r(2)+r(3)
          b=r(1)**2+r(2)**2+r(3)**2
          dr=dsqrt(9d0*b-3d0*a**2)/a
c          write(*,*)"a,b,dr",sngl(a),sngl(b),sngl(dr)
          if(drmax.lt.dr)kele=i
          if(drmax.lt.dr)drmax=dr
100     continue
        write(*,*)"sfshape:",kele,sngl(drmax)
        sfshape=drmax
        write(10,*)'sfshape:',sfshape
c
        return
        end
