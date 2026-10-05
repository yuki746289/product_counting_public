        subroutine remesh(np,ne,nele,xx,yy,zz)
        include "header.h"
        dimension ne(ine,4)                    !要素の節点
        dimension xx(ind),yy(ind),zz(ind)
        dimension x0(ind),y0(ind),z0(ind)      !節点の候補
        dimension mmp(ind)                     !mmp(np0)   節点を加える順番
        dimension ma(ine)                      !ma(na)     多面体の要素の番号
        dimension me(ine,3)                    !me(ns,3)   多面体の全部の面
        dimension mk(ine,3)                    !mk(nk,3)   多面体の表面の面
        dimension im(ine)                      !im(ns)     0:面を残す,1:面を削除する
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
c---    確認用
        dimension mb(ine),nec(ine,4),ne0(ine,4),net(ine,4)
c---    表面、格子用
        dimension ms(ine,3)
        dimension xs(ind),ys(ind),zs(ind)
        dimension xt(ind),yt(ind),zt(ind)
c---
        dimension ms0(ine,3)
c---
        dimension ip(ind),jp(ind)
        dimension jm(4,3)
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
c*********************************
        eps=1d-10    !誤差:外心
c*********************************
c
c
c--------------------------------------------------------------------------
c
c       節点を作る
c       1)表面を作る:節点を均一に分布させる
c       2)八分木で内部に節点を作る

c       描画のリセット
        call pltrset()
c
c       面の節点
        dl0=2.2d0   !dl0:領域の幅
        ng=10       !ng:分割数
        dl=dl0/dble(ng)     !格子幅
c
c       表面の面と節点の数
        ns=0    !表面の節点の数
        nels=0  !表面の面の数
        do 100 i=1,nele
        do 100 j=1,4
          m1=ne(i,jm(j,1))
          m2=ne(i,jm(j,2))
          m3=ne(i,jm(j,3))
          do 200 k=1,nele
          do 200 l=1,4
            if(i.eq.k .and. j.eq.l)goto 200
            n1=ne(k,jm(l,1))
            n2=ne(k,jm(l,2))
            n3=ne(k,jm(l,3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m2.eq.n3 .and. m3.eq.n2 .and. m1.eq.n1) .or.
     &         (m3.eq.n3 .and. m1.eq.n2 .and. m2.eq.n1))then
              goto 100
            endif
200       continue
c         節点
          do 300 m=1,3
            do ks=1,ns
              if(ne(i,jm(j,m)).eq.ip(ks))goto 300
            enddo
            ns=ns+1
            ip(ns)=ne(i,jm(j,m))    !ip(新)=古
            jp(ne(i,jm(j,m)))=ns    !jp(古)=新
300       continue
c         面
          nels=nels+1
          ms(nels,1)=jp(m1)
          ms(nels,2)=jp(m2)
          ms(nels,3)=jp(m3)
100     continue
c
        do ks=1,ns
          xs(ks)=xx(ip(ks))
          ys(ks)=yy(ip(ks))
          zs(ks)=zz(ip(ks))
        enddo
        call pltsetnp(ns,xs,ys,zs,ind)
        write(*,*)"end pltsetnp ns:",ns
c
c       面を発生させる
c        call initsuf(ns,xs,ys,zs,nels,ms)
        call pltsetme(nels,ms,3d0,2,ind,ine)
        pause
c
c       面の節点の補間と削除
c        call resuf(ns,xs,ys,zs,nels,ms,dl)
c
c       面の節点位置の修正
        call nsoptimize(ns,nels,xs,ys,zs,ms)
c
c       4面体要素の節点位置の修正
c        call npoptimize(ns,nels,xs,ys,zs,ms)
c
c       格子を作る
        call octtree(ns,xs,ys,zs,nt,xt,yt,zt,dl0,ng,ind)   !ng:分割数
c
c       外側の格子を除く
        call delouternp(ns,nels,xs,ys,zs,ms,nt,xt,yt,zt,ind,ine)
c
c       面付近に節点を発生させる
c        call npsuf(ns,nels,xs,ys,zs,ms,nt,xt,yt,zt)
c
c       外側の格子を除く
c        call delouternp(ns,nels,xs,ys,zs,ms,nt,xt,yt,zt,ind,ine)
c
c       内側の節点
        np0=0
        do kt=1,nt
          np0=np0+1
          mmp(np0)=np0
          x0(np0)=xt(kt)
          y0(np0)=yt(kt)
          z0(np0)=zt(kt)
        enddo
c
c       表面の節点
        ns0=ns
        do ks=1,ns
          np0=np0+1
          mmp(np0)=np0
          x0(np0)=xs(ks)
          y0(np0)=ys(ks)
          z0(np0)=zs(ks)
        enddo
        dx=dl0
        dy=dl0
        dz=dl0
c
c       並べ替え : 原点に近い節点 → その節点にに近い節点 → その節点にに近い節点
c        call sort(nt,xt,yt,zt,mmp,ind)
c        call sort(np0,x0,y0,z0,mmp,ind)
        call sort2(np0,x0,y0,z0,mmp,ind)
c
c       描画のリセット
        call pltrset()
        call pltsetnp(np0,x0,y0,z0,ind)
c
c
c--------------------------------------------------------------------------
c
c       囲む要素を作る

        np=8
        xx(1)=-dx
        yy(1)=-dy
        zz(1)=-dz
        xx(2)= dx
        yy(2)=-dy
        zz(2)=-dz
        xx(3)= dx
        yy(3)=-dy
        zz(3)= dz
        xx(4)=-dx
        yy(4)=-dy
        zz(4)= dz
        xx(5)=-dx
        yy(5)= dy
        zz(5)=-dz
        xx(6)= dx
        yy(6)= dy
        zz(6)=-dz
        xx(7)= dx
        yy(7)= dy
        zz(7)= dz
        xx(8)=-dx
        yy(8)= dy
        zz(8)= dz
c
        nele=6
        ne(1,1)=5
        ne(1,2)=1
        ne(1,3)=2
        ne(1,4)=4
        ne(2,1)=8
        ne(2,2)=5
        ne(2,3)=2
        ne(2,4)=4
        ne(3,1)=2
        ne(3,2)=5
        ne(3,3)=8
        ne(3,4)=6
        ne(4,1)=7
        ne(4,2)=4
        ne(4,3)=2
        ne(4,4)=3
        ne(5,1)=2
        ne(5,2)=8
        ne(5,3)=7
        ne(5,4)=6
        ne(6,1)=7
        ne(6,2)=8
        ne(6,3)=2
        ne(6,4)=4
        call region(np,nele,ne,xx,yy,zz,ind,ine)   !外の面:ns=12
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nele,ne,2d0,2,ine)
c-------------------------------------------------------------------
c
c       要素を作るループ
c
c-------------------------------------------------------------------
        idel=0  !ループの回数(埋没点)
2000    continue   !埋没点があった時
        idel=idel+1
c
        do 2100 iroop=1,np0
          write(*,*)"--------"
          write(*,*)"iroop",iroop,"/",np0
          write(*,*)"np,nele",np,nele
c
c         隣接する要素(チェック用、時間の掛かる原因)
c          call snec(nele,ne,xx,yy,zz,nec,ind,ine)
c
c         節点の候補を決める
          mp=mmp(iroop)
c
c         節点の候補を加える
          np=np+1
          xx(np)=x0(mp)
          yy(np)=y0(mp)
          zz(np)=z0(mp)
          call pltsetnp(np,xx,yy,zz,ind)
          call pltnt(xx(np),yy(np),zz(np),10d0,1)  !1:白, 2:赤, 3:緑, 4:青, 5:黄色
c
c--------------------------------------------------------------
c
c         節点の候補を外接円内に含む要素を探す
c         ma(na):多面体の要素
          xp=x0(mp)
          yp=y0(mp)
          zp=z0(mp)
          call search(nele,ne,xx,yy,zz,xp,yp,zp,eps,na,ma,ind,ine)
          call pltne0(na,ma,ne,2d0,2,ine)
c
c--------------------------------------------------------------
c
c
c         多面体の表面を探す
c         mk(nk,3):多面体の表面の節点
c
c         面を探す

1200      continue     !ma(na)を増やした

          ns=0        !4面体の面の数
          do ka=1,na
            kele=ma(ka)
            ns=ns+1       !面1
            me(ns,1)=ne(kele,jm1(1))
            me(ns,2)=ne(kele,jm1(2))
            me(ns,3)=ne(kele,jm1(3))
            ns=ns+1       !面2
            me(ns,1)=ne(kele,jm2(1))
            me(ns,2)=ne(kele,jm2(2))
            me(ns,3)=ne(kele,jm2(3))
            ns=ns+1       !面3
            me(ns,1)=ne(kele,jm3(1))
            me(ns,2)=ne(kele,jm3(2))
            me(ns,3)=ne(kele,jm3(3))
            ns=ns+1       !面4
            me(ns,1)=ne(kele,jm4(1))
            me(ns,2)=ne(kele,jm4(2))
            me(ns,3)=ne(kele,jm4(3))
          enddo
c
c         削除する面を探す:重なる面
          do 700 ks=1,ns
700       im(ks)=0    !0:残す,1:削除する
c
          do 800 ks=1,ns-1
            m1=me(ks,1)
            m2=me(ks,2)
            m3=me(ks,3)
            do 800 js=ks+1,ns
              n1=me(js,1)
              n2=me(js,2)
              n3=me(js,3)
              if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &           (m1.eq.n1 .and. m2.eq.n3 .and. m3.eq.n2) .or.
     &           (m1.eq.n2 .and. m2.eq.n1 .and. m3.eq.n3))then  !(1,2,3)=(3,2,1),(1,3,2),(2,1,3)
                im(ks)=1
                im(js)=1
              elseif(((m1.eq.n1).or.(m1.eq.n2).or.(m1.eq.n3)) .and.    !面の番号は反時計回りに付けられているか確認する
     &                ((m2.eq.n1).or.(m2.eq.n2).or.(m2.eq.n3)) .and.
     &                ((m3.eq.n1).or.(m3.eq.n2).or.(m3.eq.n3)))then
                write(*,*)"pass"
                write(*,*)"m",m1,m2,m3
                write(*,*)"n",n1,n2,n3
                do ka=1,na
                  k=ma(ka)
                  write(*,*)"na",ka,k,ne(k,1),ne(k,2),ne(k,3),ne(k,4)
                enddo
                stop
              endif
800       continue
c
          nk=0  !多面体の表面の数
          do 1100 ks=1,ns
            if(im(ks).eq.1)goto 1100
            nk=nk+1
            mk(nk,1)=me(ks,1)
            mk(nk,2)=me(ks,2)
            mk(nk,3)=me(ks,3)
1100      continue
          if(nk.gt.300)write(*,*)"iroop,na,nk,nele",iroop,na,nk,nele
          if(nk.gt.300)then
            write(*,*)"xp",xp,yp,zp
            stop    !同一の要素を探す(同じ節点)
          endif
c          if(nk.gt.100 .and. idel.eq.1)goto 2100   !埋没点にする
c
c--------------------------------------------------------------
c
c
c         多面体の要素の確認
c         vl<0d0+epsの要素:表面を共有する要素:ma(na+1)=kele

          do 1400 kk=1,nk
            k1=np
            k2=mk(kk,1)
            k3=mk(kk,2)
            k4=mk(kk,3)

c             面(k4,k3,k2)に隣接する要素を探す
              if(k4.le.8 .and. k3.le.8 .and. k2.le.8)goto 1400 !外側の要素
              do kele=1,nele
                if(ne(kele,jm1(1)).eq.k4 .and. ne(kele,jm1(2)).eq.k3 
     &                             .and. ne(kele,jm1(3)).eq.k2)goto 1300
                if(ne(kele,jm1(1)).eq.k2 .and. ne(kele,jm1(2)).eq.k4
     &                             .and. ne(kele,jm1(3)).eq.k3)goto 1300
                if(ne(kele,jm1(1)).eq.k3 .and. ne(kele,jm1(2)).eq.k2
     &                             .and. ne(kele,jm1(3)).eq.k4)goto 1300
                if(ne(kele,jm2(1)).eq.k4 .and. ne(kele,jm2(2)).eq.k3
     &                             .and. ne(kele,jm2(3)).eq.k2)goto 1300
                if(ne(kele,jm2(1)).eq.k2 .and. ne(kele,jm2(2)).eq.k4
     &                             .and. ne(kele,jm2(3)).eq.k3)goto 1300
                if(ne(kele,jm2(1)).eq.k3 .and. ne(kele,jm2(2)).eq.k2
     &                             .and. ne(kele,jm2(3)).eq.k4)goto 1300
                if(ne(kele,jm3(1)).eq.k4 .and. ne(kele,jm3(2)).eq.k3
     &                             .and. ne(kele,jm3(3)).eq.k2)goto 1300
                if(ne(kele,jm3(1)).eq.k2 .and. ne(kele,jm3(2)).eq.k4
     &                             .and. ne(kele,jm3(3)).eq.k3)goto 1300
                if(ne(kele,jm3(1)).eq.k3 .and. ne(kele,jm3(2)).eq.k2
     &                             .and. ne(kele,jm3(3)).eq.k4)goto 1300
                if(ne(kele,jm4(1)).eq.k4 .and. ne(kele,jm4(2)).eq.k3
     &                             .and. ne(kele,jm4(3)).eq.k2)goto 1300
                if(ne(kele,jm4(1)).eq.k2 .and. ne(kele,jm4(2)).eq.k4
     &                             .and. ne(kele,jm4(3)).eq.k3)goto 1300
                if(ne(kele,jm4(1)).eq.k3 .and. ne(kele,jm4(2)).eq.k2
     &                             .and. ne(kele,jm4(3)).eq.k4)goto 1300
              enddo
              write(*,*)"no kele",kele,k4,k3,k2,k1
              do ka=1,ka
                k=ma(ka)
                write(*,*)"na",k,ne(k,1),ne(k,2),ne(k,3),ne(k,4)
              enddo
              stop
1300          continue
c           vl<1d-10の要素を探す
            vl=calvl(xx(k1),yy(k1),zz(k1),xx(k2),yy(k2),zz(k2)
     &              ,xx(k3),yy(k3),zz(k3),xx(k4),yy(k4),zz(k4))
            if(vl.lt.1d-10)then
              na=na+1
              ma(na)=kele
              write(*,*)"iroop,nk,kk",iroop,kk,vl
              goto 1200
            endif
c           ねじれた要素(時計回り)を探す
            call chkne(xx(k1),yy(k1),zz(k1),xx(k2),yy(k2),zz(k2)
     &                 ,xx(k3),yy(k3),zz(k3),xx(k4),yy(k4),zz(k4),shita)
            if( shita .gt. (89d0/180d0*PI) )then   !ねじれた要素が見つかった
              na=na+1
              ma(na)=kele
           write(*,*)"iroop,nk,kk",iroop,nk,kk,kele,sngl(shita/PI*180d0)
              goto 1200
            endif
1400      continue
c
          do kk=1,nk
            k1=np
            k2=mk(kk,1)
            k3=mk(kk,2)
            k4=mk(kk,3)
            vl=calvl(xx(k1),yy(k1),zz(k1),xx(k2),yy(k2),zz(k2)
     &              ,xx(k3),yy(k3),zz(k3),xx(k4),yy(k4),zz(k4))
            if(vl.lt.1d-10)write(*,*)"iroop,kk,vl",iroop,kk,vl
            if(vl.lt.1d-10)stop
          enddo
c
c--------------------------------------------------------------
c
c
c         多面体の要素を作る
c         nele=nele+nk
c
c         元の要素保存しておく
          nel0=0
          do ka=1,na
            nel0=nel0+1
            do 1 j=1,4
1           ne0(nel0,j)=ne(ma(ka),j)
          enddo
c
c         残す面me(ns,3)と加える節点mpで要素を作る
          nb=0  !確認用
          do 900 kk=1,nk
            nele=nele+1
            ne(nele,1)=np
            ne(nele,2)=mk(kk,1)
            ne(nele,3)=mk(kk,2)
            ne(nele,4)=mk(kk,3)
c---        確認用
            nb=nb+1
            mb(nb)=nele
900       continue
c          call pltne0(nb,mb,ne,2d0,2,ine)
c
c         作り直す要素ma(na)を削除する
c         最後の要素は持ってこれない(削除要素で削除要素を穴埋めすると削除要素が残る:最後の要素が削除要素だった時)
          nelt=0
          do 1500 kele=1,nele
            do ka=1,na  !削除する要素
              if(kele.eq.ma(ka))goto 1500
            enddo
            nelt=nelt+1
            do 1600 j=1,4
1600        net(nelt,j)=ne(kele,j)
1500      continue
c
          nele=nelt
          do 1700 i=1,nele
          do 1700 j=1,4
1700      ne(i,j)=net(i,j)
c
c         元の要素が削除されているか
          write(*,*)"na,nk",na,nk
          do 2 k0=1,nel0
          do 2 k=1,nele
            if(ne0(k0,1).eq.ne(k,1) .and. ne0(k0,2).eq.ne(k,2) .and.
     &         ne0(k0,3).eq.ne(k,3) .and. ne0(k0,4).eq.ne(k,4))then
              write(*,*)"要素が削除されていません"
              write(*,*)"np",ne(k,1),ne(k,2),ne(k,3),ne(k,4)
              stop
            endif
2         continue
c
c         体積確認
          call chkvol(nele,ne,xx,yy,zz,ind,ine)
c
          call pltsetne(nele,ne,2d0,2,ine)
2100    continue
        write(*,*)"----ループ終了:",idel,"回目----"
c
c       埋没点を探す(多面体が大きくなった時に、多面体の中に節点が残る)
        call delnp(np,nele,ne,xx,yy,zz,np0,mmp,x0,y0,z0,ic,ind,ine)
        if(ic.eq.1)then   !埋没点があった
          write(*,*)"埋没点:np0",np0
          if(idel.lt.10)goto 2000   !idel:埋没点があったので、要素を切り直した回数
        endif
        write(*,*)"np0,np:",np0,np
c
c       外側の面12個
        call region(np,nele,ne,xx,yy,zz,ind,ine)
        call snec(nele,ne,xx,yy,zz,nec,ind,ine)
c
c       領域の表面を含む要素が存在するか調べる
c        call checkms(ns0,nels,ms,xs,ys,zs,np,nele,ne,xx,yy,zz,ind,ine)
c        call checkms2(ns0,nels,ms,xs,ys,zs,np,nele,ne,xx,yy,zz,ind,ine)
        call checkms3(ns0,nels,ms,xs,ys,zs,np,nele,ne,xx,yy,zz,ind,ine)
c
c       np=1～8削除
c        write(*,*)"before edit"
c        call pltrset()
c        call pltsetnp(np,xx,yy,zz,ind)
c        call pltsetne(nele,ne,2d0,2,ine)
c        pause
        call edit(np,nele,ne,xx,yy,zz,ind,ine)
c
c       表示
        call pltrset()
        call pltsetnp(np,xx,yy,zz,ind)
        call pltsetne(nele,ne,2d0,2,ine)
        if(ic.eq.0)write(*,*)"埋没点なし"
        if(ic.eq.1)write(*,*)"埋没点:idel,np0",idel,np0
c
        return
        end
