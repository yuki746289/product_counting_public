c-----------------------------------------------------------------------
c
c      デローニ分割で要素を発生させる
c
c      引数
c      np    :節点数
c      xx(np):x座標
c      yy(np):y座標
c
c      戻り値
c      nele      :要素数
c      ne(nele,3):要素内の節点番号
c
c      内部変数
c      nb(np)      0:残す節点,1:後で削除する節点
c      nl         :節点iを外接円内に含む要素数
c      ml(nl)     :節点iを外接円内に含む要素番号
c      nm         :辺の数
c      mm(nm,2)   :辺の2つの頂点の番号
c      lm(nm)      0:共有しない辺,1:共有する辺  ※0の辺で要素を作成する
c      ne0(nel0,3):要素削除用
c      ind        :節点数の上限
c      ine        :要素数の上限
c      nk         :埋没点の数
c      mk(nk)     :埋没点の節点番号
c      np0        :新しい節点数
c      x0(np0)    :新しいx座標
c      y0(np0)    :新しいy座標
c      z0(np0)    :新しいz座標
c      m0(np)     :新しい節点番号
c
c-----------------------------------------------------------------------
       subroutine mesh(np,xx,yy,nele,ne)
       include "head.for"
       dimension ne(ine,3)
       dimension xx(ind),yy(ind)
c
c      内部変数
       dimension nb(ind)
       dimension ml(ine)
       dimension mm(ind,2),lm(ind)
       dimension ne0(ine,3)
       dimension mk(ind)
       dimension m0(ind)
       dimension x0(ind),y0(ind),z0(ind)
c
c      範囲を決める
       xmin= 999d0
       xmax=-999d0
       ymin= 999d0
       ymax=-999d0
       do i=1,np
         nb(i)=0    !残す節点
         if(xx(i).lt.xmin)xmin=xx(i)
         if(xx(i).gt.xmax)xmax=xx(i)
         if(yy(i).lt.ymin)ymin=yy(i)
         if(yy(i).gt.ymax)ymax=yy(i)
       enddo
c
       xmin=xmin-(xmax-xmin)*2d0
       xmax=xmax+(xmax-xmin)*2d0
       ymin=ymin-(ymax-ymin)*2d0
       ymax=ymax+(ymax-ymin)*2d0
       write(*,*)sngl(xmin),sngl(ymin),sngl(xmax),sngl(ymax)
c
c      節点の設定
       np=np+1
       xx(np)=xmin
       yy(np)=ymin
       nb(np)=1  !後で削除する節点

       np=np+1
       xx(np)=xmax
       yy(np)=ymin
       nb(np)=1

       np=np+1
       xx(np)=xmax
       yy(np)=ymax
       nb(np)=1

       np=np+1
       xx(np)=xmin
       yy(np)=ymax
       nb(np)=1
c
c      要素の設定
       nele=1
       ne(nele,1)=np-3
       ne(nele,2)=np-2
       ne(nele,3)=np-1

       nele=2
       ne(nele,1)=np-3
       ne(nele,2)=np-1
       ne(nele,3)=np-0
c
c       call plt_mesh(0,np,nele,ne,xx,yy,ind,ine)
c
c      デローニ分割を行う
       do 100 i=1,np-4
         write(*,*)"----"
         write(*,*)"loop:",i,np-4,nele
c
c        削除しない節点
         nb(i)=0
c
c        節点を外接円内に含む要素
c        nl,ml(nl)
         nl=0
         do kele=1,nele
           x1=xx(ne(kele,1))
           x2=xx(ne(kele,2))
           x3=xx(ne(kele,3))
           y1=yy(ne(kele,1))
           y2=yy(ne(kele,2))
           y3=yy(ne(kele,3))
c          外接円の中心と半径を計算する
           call circum(x1,x2,x3,y1,y2,y3,rr,x,y)
c
           r=dsqrt((x-xx(i))**2+(y-yy(i))**2)
           if(r.lt.rr-1d-12)then
             nl=nl+1
             ml(nl)=kele
           endif
         enddo
c
c         do kl=1,nl
c           write(*,*)"ne:",kl,ml(kl)
c     &                     ,ne(ml(kl),1),ne(ml(kl),2),ne(ml(kl),3)
c         enddo
c         pause
c
c        共有しない辺を探す
c        nm,mm(nm,2),lm(nm)
         nm=0
         do kl=1,nl !外接円内に節点を含む要素
           kele=ml(kl)
c          辺1
           nm=nm+1
           mm(nm,1)=ne(kele,1)
           mm(nm,2)=ne(kele,2)
           lm(nm)=0
c          辺2
           nm=nm+1
           mm(nm,1)=ne(kele,2)
           mm(nm,2)=ne(kele,3)
           lm(nm)=0
c          辺3
           nm=nm+1
           mm(nm,1)=ne(kele,3)
           mm(nm,2)=ne(kele,1)
           lm(nm)=0
c
c          重複確認
           do 200 jl=1,nl
             if(kl.eq.jl)goto 200
             jele=ml(jl)
             j1=ne(jele,1)
             j2=ne(jele,2)
             j3=ne(jele,3)           
c            辺1
             if(mm(nm-2,1).eq.j1 .and. mm(nm-2,2).eq.j2)lm(nm-2)=1
             if(mm(nm-2,1).eq.j2 .and. mm(nm-2,2).eq.j3)lm(nm-2)=1
             if(mm(nm-2,1).eq.j3 .and. mm(nm-2,2).eq.j1)lm(nm-2)=1
             if(mm(nm-2,2).eq.j1 .and. mm(nm-2,1).eq.j2)lm(nm-2)=1
             if(mm(nm-2,2).eq.j2 .and. mm(nm-2,1).eq.j3)lm(nm-2)=1
             if(mm(nm-2,2).eq.j3 .and. mm(nm-2,1).eq.j1)lm(nm-2)=1
c            辺2
             if(mm(nm-1,1).eq.j1 .and. mm(nm-1,2).eq.j2)lm(nm-1)=1
             if(mm(nm-1,1).eq.j2 .and. mm(nm-1,2).eq.j3)lm(nm-1)=1
             if(mm(nm-1,1).eq.j3 .and. mm(nm-1,2).eq.j1)lm(nm-1)=1
             if(mm(nm-1,2).eq.j1 .and. mm(nm-1,1).eq.j2)lm(nm-1)=1
             if(mm(nm-1,2).eq.j2 .and. mm(nm-1,1).eq.j3)lm(nm-1)=1
             if(mm(nm-1,2).eq.j3 .and. mm(nm-1,1).eq.j1)lm(nm-1)=1
c            辺3
             if(mm(nm-0,1).eq.j1 .and. mm(nm-0,2).eq.j2)lm(nm-0)=1
             if(mm(nm-0,1).eq.j2 .and. mm(nm-0,2).eq.j3)lm(nm-0)=1
             if(mm(nm-0,1).eq.j3 .and. mm(nm-0,2).eq.j1)lm(nm-0)=1
             if(mm(nm-0,2).eq.j1 .and. mm(nm-0,1).eq.j2)lm(nm-0)=1
             if(mm(nm-0,2).eq.j2 .and. mm(nm-0,1).eq.j3)lm(nm-0)=1
             if(mm(nm-0,2).eq.j3 .and. mm(nm-0,1).eq.j1)lm(nm-0)=1
200        continue
         enddo
c
c         do km=1,nm
c           if(lm(km).eq.0)write(*,*)"km:",km,mm(km,1),mm(km,2)
c         enddo
c         pause
c
c        共有しない辺で要素を作成
         do 300 km=1,nm !共有しない辺
           if(lm(km).eq.1)goto 300
           nele=nele+1
           ne(nele,1)=i
           ne(nele,2)=mm(km,1)
           ne(nele,3)=mm(km,2)
300      continue
c
c        節点を外接円内に含む要素を削除する
c        nl,ml(nl)
         nel0=0
         do 400 kele=1,nele
           do kl=1,nl
             if(ml(kl).eq.kele)goto 400
           enddo
           nel0=nel0+1
           ne0(nel0,1)=ne(kele,1)
           ne0(nel0,2)=ne(kele,2)
           ne0(nel0,3)=ne(kele,3)
c           write(*,*)"ne0:",nel0,ne0(nel0,1),ne0(nel0,2),ne0(nel0,3)
400      continue
c
         nele=nel0
         do kele=1,nele
           ne(kele,1)=ne0(kele,1)
           ne(kele,2)=ne0(kele,2)
           ne(kele,3)=ne0(kele,3)
         enddo
c
c       call plt_mesh(i,np,nele,ne,xx,yy,ind,ine)
c
100    continue    !メインループ
c
c----------------------------------------------------------------------
c
c      最初の要素・節点を削除
       nel0=0
       do 700 kele=1,nele
         if(nb(ne(kele,1)).eq.1)goto 700
         if(nb(ne(kele,2)).eq.1)goto 700
         if(nb(ne(kele,3)).eq.1)goto 700
         nel0=nel0+1
         ne0(nel0,1)=ne(kele,1)
         ne0(nel0,2)=ne(kele,2)
         ne0(nel0,3)=ne(kele,3)
700    continue
c
       nele=nel0
       do kele=1,nele
         ne(kele,1)=ne0(kele,1)
         ne(kele,2)=ne0(kele,2)
         ne(kele,3)=ne0(kele,3)
       enddo
c
c----------------------------------------------------------------------
c
c      埋没点の確認
c      一旦要素の節点になった後、多角形内に残った節点を探す
       nk=0
       do 500 i=1,np
         do kele=1,nele
           if(i.eq.ne(kele,1))goto 500
           if(i.eq.ne(kele,2))goto 500
           if(i.eq.ne(kele,3))goto 500
         enddo
         nk=nk+1
         mk(nk)=i
500    continue
       write(*,*)"埋没点:",nk
c
c      埋没点を削除する
       np0=0
       do 600 i=1,np
         do kk=1,nk
           if(i.eq.mk(kk))goto 600
         enddo
         np0=np0+1
         x0(np0)=xx(i)
         y0(np0)=yy(i)
         m0(i)=np0
600    continue
c
       np=np0
       do i=1,np
         xx(i)=x0(i)
         yy(i)=y0(i)
       enddo
c
       do kele=1,nele
         ne(kele,1)=m0(ne(kele,1))
         ne(kele,2)=m0(ne(kele,2))
         ne(kele,3)=m0(ne(kele,3))
       enddo
c
c       call plt_mesh(1,np,nele,ne,xx,yy,ind,ine)
c
c       do kele=1,nele
c         do j=1,3
c           write(*,*)"ne:",kele,j,ne(kele,j)
c         enddo
c         pause
c       enddo
c
       return
       end
