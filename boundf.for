c       
c-----------------------------------------------------------------------
c
c      直接法でマトリクスを解く場合      sa(4*ind),ibw),sf(4*ind)
c                                          nbw=4*nbw
c-----------------------------------------------------------------------       
c
c      マトリクスを解く前に第一種の境界条件を代入する
c
       subroutine boundf1(np,nele,ne,idx,sa,sf,nbw,ind,ine,ibw)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sa(ind*4,ibw),sf(ind*4),idx(0:ine),js(4,-1:3)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/pmass/sc(1001)
       common/pheat/pr(1001),znu(1001,1001)
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c
c      第一種の境界条件を与える
c
       do 1000 i=1,nb(1)
         kp1=ne(nelb(1,i),js(neb(1,i),1))  !nelb(1,nb(1)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))  !neb(1,nb(1)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))
         do 1100 mfai=0,3
           if(mfai.eq.3)then
             fai=ub(i)
           elseif(mfai.eq.2)then
             fai=vb(i)
           elseif(mfai.eq.1)then
             fai=wb(i)
           elseif(mfai.eq.0)then
             fai=pb(i)
           endif
           if(fai.eq.-999d0)goto 1100
           aii1=sa(4*kp1-mfai,nbw)     !行列の対角成分を保存
           aii2=sa(4*kp2-mfai,nbw)     !行列の対角成分を保存
           aii3=sa(4*kp3-mfai,nbw)     !行列の対角成分を保存
           nbwm=nbw-4
c          節点kp1の行列を0にする
           do 1200 k=max(1,4*kp1-mfai-nbwm),min(4*np,4*kp1-mfai+nbwm)
             sa(k,(4*kp1-mfai)-k+nbw)=0d0
             sa(4*kp1-mfai,k-(4*kp1-mfai)+nbw)=0d0
1200       continue
c          節点kp2の行列を0にする
           do 1300 k=max(1,4*kp2-mfai-nbwm),min(4*np,4*kp2-mfai+nbwm)
             sa(k,(4*kp2-mfai)-k+nbw)=0d0
             sa(4*kp2-mfai,k-(4*kp2-mfai)+nbw)=0d0
1300       continue
c          節点kp3の行列を0にする
           do 1400 k=max(1,4*kp3-mfai-nbwm),min(4*np,4*kp3-mfai+nbwm)
             sa(k,(4*kp3-mfai)-k+nbw)=0d0
             sa(4*kp3-mfai,k-(4*kp3-mfai)+nbw)=0d0
1400       continue
           sa(4*kp1-mfai,nbw)=aii1     !行列の対角成分に値を戻す
           sa(4*kp2-mfai,nbw)=aii2     !行列の対角成分に値を戻す
           sa(4*kp3-mfai,nbw)=aii3     !行列の対角成分に値を戻す
           sf(4*kp1-mfai)=aii1*fai     !列ベクトルに値代入
           sf(4*kp2-mfai)=aii2*fai     !列ベクトルに値代入
           sf(4*kp3-mfai)=aii3*fai     !列ベクトルに値代入
1100     continue
1000   continue
c
c      生成項  φ[-/m3/s]
c
c       do 1700 i=1,nb(9)
c         kele=nelb(9,i)
cc-------------------------------------------------         
c         vl=((x(2)*z(4)-x(4)*z(2))*(y(3)-y(1))     !要素keleの体積を計算
c     &      -(x(3)*z(1)-x(1)*z(3))*(y(2)-y(4)))/6d0
cc-------------------------------------------------
c         tav(tt(ne(kele,1))+tt(ne(kele,2))
c     &      +tt(ne(kele,3))+tt(ne(kele,4)))/4d0
cc-------------------------------------------------     
c         dens=calden(tav)
c         do 1500 mfai=0,3
c           if(mfai.eq.3)then
c             fai=su(i)
c           elseif(mfai.eq.2)then
c             fai=sv(i)
c           elseif(mfai.eq.1)then
c             fai=sw(i)
c           elseif(mfai.eq.0)then
c             fai=sp(i)
c           endif
c           if(fai.eq.-999d0)goto 1500
c           do 1600j=1,4
c             kp=n3(kele,j)
c             sf(4*kp-mfai)=x0*fai*vl/4d0/v0**2                    !N-S式x,y,z成分の生成項:x0φV/4/v0^2
c             if(mfai.eq.0)sf(4*kp-mfai)=dens*v0*fai*vl            !連続の式の生成項:ρ0v0φV/4/x0/Ma^2
c     &                                  /4d0/x0/sv(idx(kele))**2
c1600       continue
c1500     continue
c1700   continue
       return
       end
c
c      マトリクスを解いた後に第一種の境界条件を代入する
c
       subroutine boundf1r(ne,sf,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sf(4*ind)     !引数
       dimension js(4,-1:3)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c
       do 1000 i=1,nb(1)
         kp1=ne(nelb(1,i),js(neb(1,i),1))     !nelb(1,nb(1)):第一種の境界条件を代入した要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))     !neb(1,nb(1)):第一種の境界条件を代入した要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))
         do 1100 mfai=0,3
           if(mfai.eq.3)then
             fai=ub(i)
           elseif(mfai.eq.2)then
             fai=vb(i)
           elseif(mfai.eq.1)then
             fai=wb(i)
           elseif(mfai.eq.0)then
             fai=pb(i)
           endif
           if(fai.eq.-999d0)goto 1100
           sf(4*kp1-mfai)=fai
           sf(4*kp2-mfai)=fai
           sf(4*kp3-mfai)=fai
1100     continue
1000   continue
       return
       end
c--------------------------------------------------------------------------
c
c      SIMPLE法において第一種の境界条件を代入する場合
c
c--------------------------------------------------------------------------
c
c      N-S式のマトリクスを解く前に第一種の境界条件を代入する     sa(3*ind,ibw),sf(3*ind)
c                                                                  nbw=3*nbw
       subroutine boundf1ns(np,ne,idx,sa,sf,nbw,ip,jww,nn,mn)
       parameter(ind=19999,ine=99999,ibw=1999)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sa(4*ind,ibw),sf(4*ind),idx(0:ine),js(4,-1:3)
       dimension ip(ind),jww(ind),mn(ind),mm(ind)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/pmass/sc(1001)
       common/pheat/pr(1001),znu(1001,1001)
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c----------------------------------------------------------------------
c
c      第一種の境界条件を与える
c
c----------------------------------------------------------------------
       do i=1,nn            !nn:境界条件が代入される節点の数
         mm(i)=0            !0:境界条件が代入されていない,1:既に境界条件が代入されている
       enddo               !mn(nn):境界条件が代入される節点の番号
c
       do 1000 i=1,nb(1)
         kp1=ne(nelb(1,i),js(neb(1,i),1))  !nelb(1,nb(1)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))  !neb(1,nb(1)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))  !kp1,kp2,kp3:古
c
         ifg1=0
         ifg2=0
         ifg3=0
         do j=1,nn      !境界条件が重複している節点を検索する
c          節点kp1
           if(kp1.eq.mn(j) .and. mm(j).eq.0)then       !kp1に境界条件が代入されていない場合:mm(j)=1とする
             mm(j)=1
           elseif(kp1.eq.mn(j) .and. mm(j).eq.1)then   !kp1に既に境界条件が代入されている場合は考慮しない
             ifg1=1
           endif
c          節点kp2
           if(kp2.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp2.eq.mn(j) .and. mm(j).eq.1)then
             ifg2=1
           endif
c          節点kp3
           if(kp3.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp3.eq.mn(j) .and. mm(j).eq.1)then
             ifg3=1
           endif
         enddo
c
         kp1=jww(ip(kp1))       !jww(ip(古)):新_新
         kp2=jww(ip(kp2))
         kp3=jww(ip(kp3))
c
         do 1100 mfai=0,2
           if(mfai.eq.2)then
             fai=ub(i)
           elseif(mfai.eq.1)then
             fai=vb(i)
           elseif(mfai.eq.0)then
             fai=wb(i)
           endif
           if(fai.eq.-999d0)goto 1100    !次の物理量を考慮
c-----------------------------------------------
           if(fai.lt.0d0)then
             write(*,*)'fai in boundf1ns',fai
             pause
           endif
c-----------------------------------------------
           nbwm=nbw-3
c
c----------マトリクスに節点kp1の境界条件を入力する---------------------------------------
c
           if(ifg1.eq.1)goto 1210   !境界条件の重複を避ける
c
           aii1=sa(3*kp1-mfai,nbw)     !行列の対角成分を保存
           do 1200 k=max(1,3*kp1-mfai-nbwm),min(3*np,3*kp1-mfai+nbwm)
             sf(k)=sf(k)-sa(k,(3*kp1-mfai)-k+nbw)*fai   !(3*kp1-mfai)列の右辺への移項
             sa(k,(3*kp1-mfai)-k+nbw)=0d0               !(3*kp1-mfai)列を0にする
             sa(3*kp1-mfai,k-(3*kp1-mfai)+nbw)=0d0      !(3*kp1-mfai)行を0にする
1200       continue
           sa(3*kp1-mfai,nbw)=aii1     !行列の対角成分に値を戻す
           sf(3*kp1-mfai)=aii1*fai     !列ベクトルに値代入
1210       continue
c
c----------マトリクスに節点kp2の境界条件を入力する---------------------------------------
c
           if(ifg2.eq.1)goto 1310   !境界条件の重複を避ける
c
           aii2=sa(3*kp2-mfai,nbw)     !行列の対角成分を保存
           do 1300 k=max(1,3*kp2-mfai-nbwm),min(3*np,3*kp2-mfai+nbwm)
             sf(k)=sf(k)-sa(k,(3*kp2-mfai)-k+nbw)*fai   !(3*kp2-mfai)列の右辺への移項
             sa(k,(3*kp2-mfai)-k+nbw)=0d0               !(3*kp2-mfai)列を0にする
             sa(3*kp2-mfai,k-(3*kp2-mfai)+nbw)=0d0      !(3*kp2-mfai)行を0にする
1300       continue
           sa(3*kp2-mfai,nbw)=aii2     !行列の対角成分に値を戻す
           sf(3*kp2-mfai)=aii2*fai     !列ベクトルに値代入
1310       continue
c
c----------マトリクスに節点kp3の境界条件を代入する---------------------------------------
c
           if(ifg3.eq.1)goto 1410   !境界条件の重複を避ける
c
           aii3=sa(3*kp3-mfai,nbw)     !行列の対角成分を保存
           do 1400 k=max(1,3*kp3-mfai-nbwm),min(3*np,3*kp3-mfai+nbwm)
             sf(k)=sf(k)-sa(k,(3*kp3-mfai)-k+nbw)*fai   !(3*kp3-mfai)列の右辺への移項
             sa(k,(3*kp3-mfai)-k+nbw)=0d0               !(3*kp3-mfai)列を0にする
             sa(3*kp3-mfai,k-(3*kp3-mfai)+nbw)=0d0      !(3*kp3-mfai)行を0にする
1400       continue
           sa(3*kp3-mfai,nbw)=aii3     !行列の対角成分に値を戻す
           sf(3*kp3-mfai)=aii3*fai     !列ベクトルに値代入
1410       continue
c
c-----------------------------------------------------------------------
c
1100     continue
1000   continue
c
c      生成項  φ[-/m3/s]
c
c       do 1700 i=1,nb(9)
c         kele=nelb(9,i)
cc-------------------------------------------------         
c         vl=((x(2)*z(4)-x(4)*z(2))*(y(3)-y(1))     !要素keleの体積を計算
c     &      -(x(3)*z(1)-x(1)*z(3))*(y(2)-y(4)))/6d0
cc-------------------------------------------------
c         tav(tt(ne(kele,1))+tt(ne(kele,2))
c     &      +tt(ne(kele,3))+tt(ne(kele,4)))/4d0
cc-------------------------------------------------     
c         dens=calden(tav)
c         do 1500 mfai=0,3
c           if(mfai.eq.3)then
c             fai=su(i)
c           elseif(mfai.eq.2)then
c             fai=sv(i)
c           elseif(mfai.eq.1)then
c             fai=sw(i)
c           elseif(mfai.eq.0)then
c             fai=sp(i)
c           endif
c           if(fai.eq.-999d0)goto 1500
c           do 1600j=1,4
c             kp=n3(kele,j)
c             sf(4*kp-mfai)=x0*fai*vl/4d0/v0**2                    !N-S式x,y,z成分の生成項:x0φV/4/v0^2
c             if(mfai.eq.0)sf(4*kp-mfai)=dens*v0*fai*vl            !連続の式の生成項:ρ0v0φV/4/x0/Ma^2
c     &                                  /4d0/x0/sv(idx(kele))**2
c1600       continue
c1500     continue
c1700   continue
       return
       end
c
c      連続の式のマトリクスを解く前に第一種の境界条件を代入する    sa(ind,ibw),sf(ind)
c                                                                    nbw=1*nbw
       subroutine boundf1cn(np,ne,idx,sa,sf,nbw,ip,jww,nn,mn)
       parameter(ind=19999,ine=99999,ibw=1999)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sa(4*ind,ibw),sf(4*ind),idx(0:ine),js(4,-1:3)
       dimension ip(ind),jww(ind),mn(ind),mm(ind)   !境界条件の重複を避ける(左辺に移項の際，考慮している行のfaiのみ値が更新される．他の行のsaは0)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/pmass/sc(1001)
       common/pheat/pr(1001),znu(1001,1001)
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c----------------------------------------------------------------------
c
c      第一種の境界条件を与える
c
c----------------------------------------------------------------------
       do i=1,nn
         mm(i)=0
       enddo
c
       do 1000 i=1,nb(1)      
         kp1=ne(nelb(1,i),js(neb(1,i),1))  !nelb(1,nb(1)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))  !neb(1,nb(1)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))
c
         ifg1=0
         ifg2=0
         ifg3=0
         do j=1,nn      !境界条件が重複している節点を検索する
c          節点kp1
           if(kp1.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp1.eq.mn(j) .and. mm(j).eq.1)then
             ifg1=1
           endif
c          節点kp2
           if(kp2.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp2.eq.mn(j) .and. mm(j).eq.1)then
             ifg2=1
           endif
c          節点kp3
           if(kp3.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp3.eq.mn(j) .and. mm(j).eq.1)then
             ifg3=1
           endif
         enddo
c
         kp1=jww(ip(kp1))       !jww(ip(古)):新_新
         kp2=jww(ip(kp2))
         kp3=jww(ip(kp3))
c
         fai=pb(i)
         if(fai.eq.-999d0)goto 1000
c-----------------------------------------------
         if(fai.lt.0d0)then
           write(*,*)'fai in boundf1cn',fai
           pause
         endif
c-----------------------------------------------
         nbwm=nbw-1
c
c----------マトリクスに節点kp1の境界条件を入力する---------------------------------------
c
         if(ifg1.eq.1)goto 1110   !境界条件の重複を避ける
c
         aii1=sa(kp1,nbw)     !行列の対角成分を保存
         do 1100 k=max(1,kp1-nbwm),min(np,kp1+nbwm)
           sf(k)=sf(k)-sa(k,kp1-k+nbw)*fai   !kp1列の右辺への移項
           sa(k,kp1-k+nbw)=0d0               !kp1列を0にする
           sa(kp1,k-kp1+nbw)=0d0             !kp1行を0にする
1100     continue
         sa(kp1,nbw)=aii1     !行列の対角成分に値を戻す
         sf(kp1)=aii1*fai     !列ベクトルに値代入
1110     continue
c
c----------マトリクスに節点kp2の境界条件を入力する---------------------------------------
c
         if(ifg2.eq.1)goto 1210   !境界条件の重複を避ける
c
         aii2=sa(kp2,nbw)     !行列の対角成分を保存
         do 1200 k=max(1,kp2-nbwm),min(np,kp2+nbwm)
           sf(k)=sf(k)-sa(k,kp2-k+nbw)*fai   !kp2列の右辺への移項
           sa(k,kp2-k+nbw)=0d0               !kp2列を0にする
           sa(kp2,k-kp2+nbw)=0d0             !kp2行を0にする
1200     continue
         sa(kp2,nbw)=aii2     !行列の対角成分に値を戻す
         sf(kp2)=aii2*fai     !列ベクトルに値代入
1210     continue
c
c----------マトリクスに節点kp3の境界条件を入力する---------------------------------------
c
         if(ifg3.eq.1)goto 1310   !境界条件の重複を避ける
c
         aii3=sa(kp3,nbw)     !行列の対角成分を保存
         do 1300 k=max(1,kp3-nbwm),min(np,kp3+nbwm)
           sf(k)=sf(k)-sa(k,kp3-k+nbw)*fai   !kp3列の右辺への移項
           sa(k,kp3-k+nbw)=0d0               !kp3列を0にする
           sa(kp3,k-kp3+nbw)=0d0             !kp3行を0にする
1300     continue
         sa(kp3,nbw)=aii3     !行列の対角成分に値を戻す
         sf(kp3)=aii3*fai     !列ベクトルに値代入
1310     continue
c
c-----------------------------------------------------------------------
c
1000   continue
c
c      生成項  φ[-/m3/s]
c
c       do 1400 i=1,nb(9)
c         kele=nelb(9,i)
cc-------------------------------------------------         
c         vl=((x(2)*z(4)-x(4)*z(2))*(y(3)-y(1))     !要素keleの体積を計算
c     &      -(x(3)*z(1)-x(1)*z(3))*(y(2)-y(4)))/6d0
cc-------------------------------------------------
c         tav(tt(ne(kele,1))+tt(ne(kele,2))
c     &      +tt(ne(kele,3))+tt(ne(kele,4)))/4d0
cc-------------------------------------------------     
c         dens=calden(tav)
c         do 1500 mfai=0,3
c           if(mfai.eq.3)then
c             fai=su(i)
c           elseif(mfai.eq.2)then
c             fai=sv(i)
c           elseif(mfai.eq.1)then
c             fai=sw(i)
c           elseif(mfai.eq.0)then
c             fai=sp(i)
c           endif
c           if(fai.eq.-999d0)goto 1500
c           do 1600j=1,4
c             kp=n3(kele,j)
c             sf(4*kp-mfai)=x0*fai*vl/4d0/v0**2                    !N-S式x,y,z成分の生成項:x0φV/4/v0^2
c             if(mfai.eq.0)sf(4*kp-mfai)=dens*v0*fai*vl            !連続の式の生成項:ρ0v0φV/4/x0/Ma^2
c     &                                  /4d0/x0/sv(idx(kele))**2
c1600       continue
c1500     continue
c1400   continue
       return
       end
c
c      N-S式のマトリクスを解いた後に第一種の境界条件を代入する   sf(3*ind)
c
       subroutine boundf1nsr(ne,sf,ip,jww,nn,mn,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sf(4*ind)     !引数
       dimension ip(ind),jww(ind),mn(ind),mm(ind)
       dimension js(4,-1:3)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c----------------------------------------------------------------------
c
c      第一種の境界条件を与える
c
c----------------------------------------------------------------------
       do i=1,nn
         mm(i)=0
       enddo
c
       do 1000 i=1,nb(1)      
         kp1=ne(nelb(1,i),js(neb(1,i),1))  !nelb(1,nb(1)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))  !neb(1,nb(1)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))
c
         ifg1=0
         ifg2=0
         ifg3=0
         do j=1,nn      !境界条件が重複している節点を検索する
c          節点kp1
           if(kp1.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp1.eq.mn(j) .and. mm(j).eq.1)then
             ifg1=1
           endif
c          節点kp2
           if(kp2.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp2.eq.mn(j) .and. mm(j).eq.1)then
             ifg2=1
           endif
c          節点kp3
           if(kp3.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp3.eq.mn(j) .and. mm(j).eq.1)then
             ifg3=1
           endif
         enddo
c
         kp1=jww(ip(kp1))       !jww(ip(古)):新_新
         kp2=jww(ip(kp2))
         kp3=jww(ip(kp3))
c
         do 1100 mfai=0,2
           if(mfai.eq.2)then
             fai=ub(i)
           elseif(mfai.eq.1)then
             fai=vb(i)
           elseif(mfai.eq.0)then
             fai=wb(i)
           endif
           if(fai.eq.-999d0)goto 1100
c--------------------------------------------------------------------------
c          列ベクトルに節点kp1の境界条件を入力する
c--------------------------------------------------------------------------
           if(ifg1.eq.1)goto 1200   !境界条件の重複を避ける
c
           sf(3*kp1-mfai)=fai
1200       continue
c--------------------------------------------------------------------------
c          列ベクトルに節点kp1の境界条件を入力する
c--------------------------------------------------------------------------
           if(ifg2.eq.1)goto 1300   !境界条件の重複を避ける
c
           sf(3*kp2-mfai)=fai
1300       continue
c--------------------------------------------------------------------------
c          列ベクトルに節点kp1の境界条件を入力する
c--------------------------------------------------------------------------
           if(ifg3.eq.1)goto 1400   !境界条件の重複を避ける
c
           sf(3*kp3-mfai)=fai
1400       continue
c--------------------------------------------------------------------------
1100     continue
1000   continue
       return
       end
c
c      連続の式のマトリクスを解いた後に第一種の境界条件を代入する
c
       subroutine boundf1cnr(ne,sf,ip,jww,nn,mn,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sf(4*ind)     !引数
       dimension ip(ind),jww(ind),mn(ind),mm(ind)
       dimension js(4,-1:3)
       common/prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c----------------------------------------------------------------------
c
c      第一種の境界条件を与える
c
c----------------------------------------------------------------------
       do i=1,nn
         mm(i)=0
       enddo
c
       do 1000 i=1,nb(1)      
         kp1=ne(nelb(1,i),js(neb(1,i),1))  !nelb(1,nb(1)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(1,i),js(neb(1,i),2))  !neb(1,nb(1)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(1,i),js(neb(1,i),3))
c
         ifg1=0
         ifg2=0
         ifg3=0
         do j=1,nn      !境界条件が重複している節点を検索する
c          節点kp1
           if(kp1.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp1.eq.mn(j) .and. mm(j).eq.1)then
             ifg1=1
           endif
c          節点kp2
           if(kp2.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp2.eq.mn(j) .and. mm(j).eq.1)then
             ifg2=1
           endif
c          節点kp3
           if(kp3.eq.mn(j) .and. mm(j).eq.0)then
             mm(j)=1
           elseif(kp3.eq.mn(j) .and. mm(j).eq.1)then
             ifg3=1
           endif
         enddo
c
         kp1=jww(ip(kp1))       !jww(ip(古)):新_新
         kp2=jww(ip(kp2))
         kp3=jww(ip(kp3))
c
         fai=pb(i)
         if(fai.eq.-999d0)goto 1000
c--------------------------------------------------------------------------
c          列ベクトルに節点kp1の境界条件を入力する
c--------------------------------------------------------------------------
         if(ifg1.eq.1)goto 1100   !境界条件の重複を避ける
c
         sf(kp1)=fai
1100     continue
c--------------------------------------------------------------------------
c          列ベクトルに節点kp2の境界条件を入力する
c--------------------------------------------------------------------------
         if(ifg2.eq.1)goto 1200   !境界条件の重複を避ける
c
         sf(kp2)=fai
1200     continue
c--------------------------------------------------------------------------
c          列ベクトルに節点kp3の境界条件を入力する
c--------------------------------------------------------------------------
         if(ifg3.eq.1)goto 1300   !境界条件の重複を避ける
c
         sf(kp3)=fai
1300     continue
c--------------------------------------------------------------------------
1000   continue
       return
       end