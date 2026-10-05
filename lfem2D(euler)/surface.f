c-----------------------------------------------------------
c
c       領域、界面を探す
c
c       引数
c       nele           :要素の数
c       ne(nele,3)     :要素内の節点番号
c       idx(nele)      :要素のindex番号
c       ine            :要素数
c       ica            :領域数
c       ibd            :境界の数
c
c       戻り値
c       nr             :領域の数
c       nc(nr)         :領域nrの要素の数
c       nelc(nr,nc(nr)):領域nrの要素の番号
c       ndx(nr)        :領域のindex番号
c       neled(nr)      :領域nrの最後の要素の番号
c
c       na              :境界の合計
c       mr(na)          :境界naの領域番号
c       nb(na)          :境界の接点数
c       nelb(na,nb(na)) :境界の要素の番号   
c       neb(na,nb(na))) :境界の要素内の番号
c       ncom(na,nb(na)) :境界のindex番号
c       dcom(na,nb(na))  -1d0:凹, 1d0:凸
c
c-----------------------------------------------------------
        subroutine sufele(nele,ne,idx,nr,nc,nelc
     &                    ,ndx,neled,na,mr,nb,nelb,neb,ncom)
        include 'head.for'
c
c       外部変数
        dimension ne(ine,3),idx(0:ine)
        dimension nc(ica),nelc(ica,-1:ibd)
        dimension ndx(0:ica)
        dimension neled(0:ica)
        dimension mr(ica)
        dimension nb(ica)
        dimension nelb(ica,-1:ibd)
        dimension neb(ica,-1:ibd)
        dimension ncom(ica,-1:ibd)
c
c       内部変数
        dimension nec(ine,3)   !nec(nele,3):隣の要素
        dimension l0(ine),ll(ine)
        dimension ik(0:ine)      !0:数えていない,1:数えた
        dimension ich(ine,3)
        dimension js(3)
        data js/2,3,1/
c
c       初期化
        do 100 i=1,nele
        do 100 j=1,3
100       nec(i,j)=0
c
c       隣の要素を探す
        do 200 i=1,nele
          do 300 j=1,3
            if(nec(i,j).ne.0)goto 300
            kp1=ne(i,j)
            kp2=ne(i,js(j))
            do 400 k=1,nele
            do 400 l=1,3
              jp1=ne(k,l)
              jp2=ne(k,js(l))
              if(kp1.eq.jp2 .and. kp2.eq.jp1)nec(i,j)=k
              if(kp1.eq.jp2 .and. kp2.eq.jp1)nec(k,l)=i
400         continue
300       continue
200     continue
c
c-----------------------------------------------------------
c       nr             :領域の数
c       nc(nr)         :領域nrの要素の数
c       nelc(nr,nc(nr)):領域nrの要素の番号
c       ndx(nr)        :領域のindex番号

        do 2100 kele=0,nele
2100    ik(kele)=0      !0:数えていない,1:数えた
        kele=1          !次の領域の元になる要素
c
        nr=0
600     continue       !領域が無くなるまで繰り返す

        nr=nr+1
        nc(nr)=1                !領域nrの要素の数.kele分
        nelc(nr,nc(nr))=kele    !kele:領域nrの最初の要素
        ndx(nr)=idx(kele)
        ik(kele)=1

        m0=1
        l0(1)=kele
        do 700 ic=1,ine
          mm=0
          do 800 i=1,m0
            le=l0(i)        !領域nrの要素
c           --------
            ke1=nec(le,1)   !隣の要素
            if(ke1.eq.0 .or. ik(ke1).eq.1)goto 2300    !隣に要素がない,既に数えた要素
            if(idx(ke1).eq.ndx(nr))then
              nc(nr)=nc(nr)+1       !領域nrの要素の数
              nelc(nr,nc(nr))=ke1   !領域nrの要素の番号
              ik(ke1)=1             !1:数えた
              mm=mm+1               !周りの要素の数
              ll(mm)=ke1            !周りの要素の番号
            endif
c           --------
2300        ke2=nec(le,2)   !隣の要素
            if(ke2.eq.0 .or. ik(ke2).eq.1)goto 2400    !隣に要素がない,既に数えた要素
            if(idx(ke2).eq.ndx(nr))then
              nc(nr)=nc(nr)+1       !領域nrの要素の数
              nelc(nr,nc(nr))=ke2   !領域nrの要素の番号
              ik(ke2)=1             !1:数えた
              mm=mm+1               !周りの要素の数
              ll(mm)=ke2            !周りの要素の番号
            endif
c           --------
2400        ke3=nec(le,3)   !隣の要素
            if(ke3.eq.0 .or. ik(ke3).eq.1)goto 800     !隣に要素がない,既に数えた要素
            if(idx(ke3).eq.ndx(nr))then
              nc(nr)=nc(nr)+1       !領域nrの要素の数
              nelc(nr,nc(nr))=ke3   !領域nrの要素の番号
              ik(ke3)=1             !1:数えた
              mm=mm+1               !周りの要素の数
              ll(mm)=ke3            !周りの要素の番号
            endif
800       continue
          if(mm.eq.0)goto 1000 !ループを抜ける
          m0=mm
          do 900 j=1,m0
900       l0(j)=ll(j)
700     continue
1000    continue
c
c       次の領域の要素keleを決める      ※※※※※※時間のかかるループ※※※※※
        kele=0
        do 1100 i=1,nele
          do 1200 kr=1,nr
          do 1200 kc=1,nc(kr)
            if(i.eq.nelc(kr,kc))goto 1100  !数えた要素はとばす
1200      continue
          kele=i
          goto 2200
1100    continue

2200    if(kele.ne.0)goto 600
c
c-----------------------------------------------------------
c       neled(nr)       :領域nrの最後の要素の番号

        do 500 i=0,ica
500     neled(i)=0
        do kr=1,nr
        neled(kr)=neled(kr-1)+nc(kr)
        enddo
c
c-----------------------------------------------------------
c       na              :境界の合計
c       mr(na)          :境界naの領域番号
c       nb(na)          :境界の接点数
c       nelb(na,nb(na)) :境界の要素の番号   
c       neb(na,nb(na))  :境界の要素内の番号
c       ncom(na,nb(na)  :境界のindex番号

        na=0
        do 2000 kr=1,nr
c         チェック用
          do 1300 i=1,nele
          do 1300 j=1,3
1300      ich(i,j)=0      !0:チェックしてない,1:チェックした境界要素
c
          do ib=1,ibd     !同じ領域の境界の数
c        
c           最初の境界の要素:kele
            do 1400 i=1,nc(kr)     !領域krの要素の数
              kele=nelc(kr,i)
              do 1500 j=1,3
                if(ich(kele,j).eq.1)goto 1500     !チェック
                ich(kele,j)=1
                if(idx(kele).ne.idx(nec(kele,j)))goto 1600     !境界が見つかった
1500          continue
1400        continue
            goto 2000 !次の領域:境界が見つからなかった
c
c           境界の領域番号
1600        na=na+1
            mr(na)=kr
            nb(na)=1
            nelb(na,1)=kele
            neb(na,1)=j
            ncom(na,1)=idx(nec(kele,j))
c
c           1周する
            do 1700 ic=1,ibd
              kp=ne(nelb(na,nb(na)),js(neb(na,nb(na))))       !次の境界の節点
              if(kp.eq.ne(nelb(na,1),neb(na,1)))goto 1900    !1周した
c             次の境界の節点を探す
              do 1800 i=1,nc(kr)
              do 1800 j=1,3
                jp=ne(nelc(kr,i),j)
                if(kp.eq.jp)then
                  if(idx(nec(nelc(kr,i),j)).eq.idx(nelc(kr,i)))goto 1800   !境界ではない要素はとばす
                  nb(na)=nb(na)+1
                  nelb(na,nb(na))=nelc(kr,i)
                  neb(na,nb(na))=j
                  ncom(na,nb(na))=idx(nec(nelc(kr,i),j))      !idx(0)=0   sub init()
                  ich(nelc(kr,i),j)=1    !チェックした境界
                  goto 1700
                endif
1800          continue
1700        continue

1900        continue !1周した
            nelb(na,-1      )=nelb(na,nb(na)-1)   !2つ前
            nelb(na, 0      )=nelb(na,nb(na)-0)   !1つ前
            nelb(na,nb(na)+1)=nelb(na,1)          !1つ先
            nelb(na,nb(na)+2)=nelb(na,2)          !2つ先
            neb(na,-1      )=neb(na,nb(na)-1)
            neb(na, 0      )=neb(na,nb(na)-0)
            neb(na,nb(na)+1)=neb(na,1)
            neb(na,nb(na)+2)=neb(na,2)
            ncom(na,-1     )=ncom(na,nb(na)-1)
            ncom(na, 0     )=ncom(na,nb(na)-0)
            ncom(na,nb(na)+1)=ncom(na,1)
            ncom(na,nb(na)+2)=ncom(na,2)
          enddo
2000    continue
c
        return
        end
