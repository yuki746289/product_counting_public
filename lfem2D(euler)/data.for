c--------------------------------------------------------------------
c
c      計算条件を読み込む
c
c      iyy          0:最初から計算,1:途中から計算
c      dt          :タイムステップ
c      wdt1        :描画データ保存時間
c      wdt2        :計算データ保存時間
c      timax       :計算終了時間
c      vi          :代表速度
c      ri          :代表長さ
c      vre(idx)    :Re数
c      vma(idx)    :Ma数
c      vfr(idx)    :Fr数
c      vwe(idx,idx):We数
c
c      np          :節点数
c      nele        :要素数
c      ne(nele,3)  :要素番号
c      idx(nele)   :index番号
c      xx(np)      :x座標
c      yy(np)      :y座標
c
c--------------------------------------------------------------------
       subroutine data(iyy,dt,wdt1,wdt2,tmax,vi,ri,vre,vma,vfr,vwe
     &                  ,np,nele,ne,idx,xx,yy)
       include "head.for"
       dimension ne(ine,3)
       dimension idx(0:ine)
       dimension xx(ind),yy(ind),zz(ind)
       dimension vre(-5:ica),vma(-5:ica),vfr(-5:ica)
       dimension vwe(-5:ica,-5:ica)    !※※※※※※※
c
c      計算条件
       open(5,file="calc.inf")
       rewind(1)
c
c      開始時間
c      0:最初から計算, 1:途中から計算
       read(5,*) iyy
c
c      時間関係
       read(5,*) dt,wdt1,wdt2,tmax
c
c      代用速度・代用長さ
       read(5,*)vi,ri
c       
c      index領域の物性値
       read(5,*)n
       read(5,*)(i,vre(i),vma(i),vfr(i),i=1,n)
c       
c      index領域の界面の物性値
       read(5,*)n
       do i=1,n
         read(5,*)j,k,vwe(j,k)
         vwe(k,j)=vwe(j,k)
       enddo
       close(1)
c
c      グリッドファイル
       open(5,file="mesh.inp")
       rewind(1)
       read(5,*)np,nele
       read(5,*)(i,xx(i),yy(i),i=1,np)
       read(5,*)(i,ne(i,1),ne(i,2),ne(i,3),idx(i),i=1,nele)
       close(1)
c
        write(*,*)'np,nele',np,nele
        if(np*3.gt.ind .or. nele.gt.ine) then
          write(*,*) 'ine,ind=',ine,ind
          stop
        endif
c
       return
       end
