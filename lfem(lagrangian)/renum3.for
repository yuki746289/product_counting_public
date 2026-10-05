c       subroutine renum3(np,nele,ne,jwww,jmmm,nbmax,ind,ine,ibw)
       subroutine renum3(np,nele,ne,jwww,jmmm,nbmax)
       include "header.h"
       dimension ne(ine,4)           !各要素を構成する節点の番号
       dimension mb(ind)             !節点番号kpの節点と辺を共有している節点の数
c       integer*8 nsurn(ind,100)     !ある節点の周囲に存在する節点を格納した配列
       dimension iwn(ind),iws(ind),iwp(ind)    !辺を共有する節点の番号:nsurn(iwn(iw),iws(iw))=iwp(iw)
       dimension js(4,-1:3)          !要素内のある節点と辺を共有する節点の要素内番号
       dimension jww(ind),jmm(ind)   !jww(古)=新,jmm(新)=古
       dimension jwww(ind),jmmm(ind) !戻り値:jwww(古)=新,jmmm(新)=古
       dimension ln0(ind),ln(ind)    !ln0(n0),ln(nn):前段または考慮している段に存在している節点の番号
c
       dimension iwf(ind)            !iwf(kp):考慮している節点が格納されている配列iwn(kw)におけるkwの最小値
c       write(*,*)'in renum3,np',np       
c      要素内の節点番号1の節点
       js(1,-1)=3
       js(1,0)=4
       js(1,1)=2
       js(1,2)=3
       js(1,3)=4
c      要素内の節点番号2の節点
       js(2,-1)=3
       js(2,0)=4
       js(2,1)=1
       js(2,2)=3
       js(2,3)=4
c      要素内の節点番号3の節点
       js(3,-1)=2
       js(3,0)=4
       js(3,1)=1
       js(3,2)=2
       js(3,3)=4
c      要素内の節点番号4の節点
       js(4,-1)=2
       js(4,0)=3
       js(4,1)=1
       js(4,2)=2
       js(4,3)=3
c       write(*,*)'decided the arrangement of js(4,-1:3)'
c
c     iwf(np)を初期化する
c
      do kp=1,ind
        iwf(kp)=1     !iwf(kp):考慮している節点が格納されている配列iwn(kw)におけるkwの最小値
      enddo
c
c      mb(np),nsurn(np,mb(np))を求める
c
       iw=0
       do 1000 kp=1,np  !kp:考慮している節点
c         write(*,*)'kp,np',kp,np
c         pause
         mb(kp)=0
         do 1100 i=1,nele
           do 1200 j=1,4
             if(ne(i,j).ne.kp)goto 1200
             do 1300 k=1,3
c----------------既に考慮されている節点は考慮しない-----------------------------               
               do 1400 l=1,mb(kp)
c                 if(ne(i,js(j,k)).eq.nsurn(kp,l))goto 1300
                 if(ne(i,js(j,k)).eq.nsurn(kp,l,iwn,iws,iwp,iw,iwf,ind))
     &                                                    goto 1300
1400           continue
c--------------------------------------------------------------------------------
c               write(*,*)'iw,k,',iw,k,'ne(',i,',',j,')',kp
               iw=iw+1
               mb(kp)=mb(kp)+1
               iwn(iw)=kp
               iws(iw)=mb(kp)
               iwp(iw)=ne(i,js(j,k))
c               nsurn(kp,mb(kp))=ne(i,js(j,k))
c               write(*,*)'iw,iwn,iws,iwp',iw,iwn(iw),iws(iw),iwp(iw)
1300         continue
1200       continue
1100     continue
1000   continue
c
c     iwf(np)を求める
c
      do 10 kp=1,np
        do 20 kw=1,iw
          if(kp.ne.iwn(kw))goto 20
          iwf(kp)=kw    !iwf(kp):考慮している節点が格納されている配列iwn(kw)におけるkwの最小値
          goto 10
20      continue
10    continue
c      write(*,*)'decided the arrngements of mb(np) and nsurn(np,mb(np))'
c      write(*,*)'iw',iw
c      pause
c
c      新しい節点番号および古い節点番号を格納した配列jwww(ind),jmmm(ind)を求める
c      jwww(古)=新,jmmm(新)=古
       nbmax=9999999         !バンド幅
       mmax=9999999
       do 1500 nrt=1,np     !root node
c         write(*,*)"the number of the considering root node",nrt
c
         do 1800 i=1,np
         jww(i)=0
         jmm(i)=0
1800     continue
c
         n0=1
         nmax=-99999
         nmt=0
         mm=1
         ln0(1)=nrt
         jww(nrt)=1
         jmm(1)=nrt
2000     continue
         nn=0
         do 1600 k0=1,n0
           do 1700 kb=1,mb(ln0(k0))
             do j=1,mm     !既に考慮されている節点は考慮しない
               if(nsurn(ln0(k0),kb,iwn,iws,iwp,iw,iwf,ind).eq.jmm(j))
     &                                               goto 1700
             enddo
             mm=mm+1                                          !mm:新しい節点の番号
             jww(nsurn(ln0(k0),kb,iwn,iws,iwp,iw,iwf,ind))=mm !jww(古)=新
             jmm(mm)=nsurn(ln0(k0),kb,iwn,iws,iwp,iw,iwf,ind) !jmm(新)=古
             nn=nn+1                                          !nn:現在考慮している段に存在する節点の数
             ln(nn)=nsurn(ln0(k0),kb,iwn,iws,iwp,iw,iwf,ind)  !ln(nn):現在考慮している段に存在する節点の番号
c             write(*,*)'----------------------------------------------'
c             write(*,*)'nrt,n0,k0,mb(',ln0(k0),'),kb',
c     &                          nrt,n0,k0,mb(ln0(k0)),kb
c             write(*,*)'mm,nn',mm,nn
c             write(*,*)'----------------------------------------------'
c             pause
1700       continue
1600     continue
         nmt=nn+n0                   !前段と現在考慮している段に存在している節点の合計
         if(nmt.gt.nmax)nmax=nmt     !nmax:現在考慮しているtreeにおける隣接する節点の番号の差最大値
c         write(*,*)'nrt nmax,mmax',nrt,nmax,mmax
         if(nmax.gt.mmax)goto 1500  !現在考慮しているtreeの最大値nmaxが前回までの最大値mmaxより大きい場合は現在の考慮しているtreeを考慮しない
         n0=nn                       !n0:前段に存在している節点の数
         do i=1,nn
           ln0(i)=ln(i)              !ln0(n0):前段に存在している節点の番号
         enddo
         if(nn.ne.0)goto 2000       !全節点を考慮するまで
c         
         mmax=nmax
c
c        バンド幅を求める
c         write(*,*)'calc band'
         nbm=-99999
         do kele=1,nele
           do i=1,4
             do j=1,4
               kbm=iabs(jww(ne(kele,i))-jww(ne(kele,j)))
               if(kbm.gt.nbm)nbm=kbm     !nbm:現在考慮しているtreeのバンド幅
             enddo
           enddo
         enddo
c         write(*,*)'nbm',nbm
c         pause
c
         if(nbm.lt.nbmax)then
           nbmax=nbm
           do i=1,np
             jwww(i)=jww(i)
             jmmm(i)=jmm(i)
           enddo 
c           write(*,*)'root',nrt,'bnd',nbmax+1
         endif
1500   continue
       nbmax=nbmax+1     !バンド幅
       write(*,*)'---- band interval',nbmax,' ----'
c
       return
       end
c--------------------------------------------------------------------------
c
c       サブルーチン
c
c--------------------------------------------------------------------------
c      辺を共有する節点の番号:nsurn(iwn(iw),iws(iw))=iwp(iw)
       function nsurn(kp,l,iwn,iws,iwp,iw,iwf,ind)
       implicit double precision(a-h,o-z)
       dimension iwn(ind),iws(ind),iwp(ind),iwf(ind)
       do 100 kw=iwf(kp),iw   !kp:考慮している節点,l:考慮している節点の周囲に存在するl番目の節点
         if(iwn(kw).eq.kp .and. iws(kw).eq.l)then
           nsurn=iwp(kw)
           goto 200
         endif
100    continue
       write(*,*)"a node surrounding the considering one didn't find"
       stop
200    continue
       return
       end