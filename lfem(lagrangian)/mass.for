c-----------------------------------------------------------------------------
c
c      物質移動に関するマトリクスを発生させる
c
c-----------------------------------------------------------------------------
       subroutine mtmass(nel,nbw,sa,sf,ne,xx,yy,zz    !nea(ine,4):物質移動領域の要素のみ格納
     &               ,uu,vv,ww,c0,t0,jmc,jpc,idx,ind,ine,ibw)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)
       dimension uu(ind),vv(ind),ww(ind)   !uu(古),vv(古),ww(古)
       dimension jmc(ind),jpc(ind)         !jpc(jmc(新_新)):古
       dimension sa(4*ind,ibw),sf(4*ind)
       dimension x(4),y(4),z(4),b(4),c(4),d(4),zkesu(4,4)
       dimension c0(ind),t0(ind)                   !c0(新_新),t0(新_新)
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
       common /pmass/sc(1001)
       common /time/dt,tmax,tsv,ttmp    !dtを用いる
       data zkesu/2d0,1d0,1d0,1d0,
     &            1d0,2d0,1d0,1d0,
     &            1d0,1d0,2d0,1d0,
     &            1d0,1d0,1d0,2d0/
c
       do 200 im=1,nel     !要素imについて:新_新
         id=idx(im)
         usum=0d0
         vsum=0d0
         wsum=0d0
         do 210 i=1,4      !要素imを構成する節点の座標
           x(i)=xx(ne(im,i))
           y(i)=yy(ne(im,i))
           z(i)=zz(ne(im,i))
           usum=usum+uu(jpc(jmc(ne(im,i))))  !jpc(jmc(新_新)):古
           vsum=vsum+vv(jpc(jmc(ne(im,i))))
           wsum=wsum+ww(jpc(jmc(ne(im,i))))
210      continue
c
c        係数行列の成分を算出する
c
         b(1)=-y(3)*z(4)+z(3)*y(4)+y(2)*z(4)
     &        -y(2)*z(3)-z(2)*y(4)+z(2)*y(3)
         b(2)= y(3)*z(4)-z(3)*y(4)-y(1)*z(4)
     &        +y(1)*z(3)+z(1)*y(4)-z(1)*y(3)
         b(3)=-y(2)*z(4)+z(2)*y(4)+y(1)*z(4)
     &        -y(1)*z(2)-z(1)*y(4)+z(1)*y(2)
         b(4)= y(2)*z(3)-z(2)*y(3)-y(1)*z(3)
     &        +y(1)*z(2)+z(1)*y(3)-z(1)*y(2)
         c(1)= x(3)*z(4)-z(3)*x(4)-x(2)*z(4)
     &        +x(2)*z(3)+z(2)*x(4)-z(2)*x(3)
         c(2)=-x(3)*z(4)+z(3)*x(4)+x(1)*z(4)
     &        -x(1)*z(3)-z(1)*x(4)+z(1)*x(3)
         c(3)= x(2)*z(4)-z(2)*x(4)-x(1)*z(4)
     &        +x(1)*z(2)+z(1)*x(4)-z(1)*x(2)
         c(4)=-x(2)*z(3)+z(2)*x(3)+x(1)*z(3)
     &        -x(1)*z(2)-z(1)*x(3)+z(1)*x(2)
         d(1)=-x(3)*y(4)+y(3)*x(4)+x(2)*y(4)
     &        -x(2)*y(3)-y(2)*x(4)+y(2)*x(3)
         d(2)= x(3)*y(4)-y(3)*x(4)-x(1)*y(4)
     &        +x(1)*y(3)+y(1)*x(4)-y(1)*x(3)
         d(3)=-x(2)*y(4)+y(2)*x(4)+x(1)*y(4)
     &        -x(1)*y(2)-y(1)*x(4)+y(1)*x(2)
         d(4)= x(2)*y(3)-y(2)*x(3)-x(1)*y(3)
     &        +x(1)*y(2)+y(1)*x(3)-y(1)*x(2)
c
c        体積Vを算出する
c
         a11=1d0
         a21=1d0
         a31=1d0
         a41=1d0
         a12=x(1)
         a22=x(2)
         a32=x(3)
         a42=x(4)
         a13=y(1)
         a23=y(2)
         a33=y(3)
         a43=y(4)
         a14=z(1)
         a24=z(2)
         a34=z(3)
         a44=z(4)
c
         det=a11*a22*a33*a44-a11*a22*a34*a43-a11*a23*a32*a44
     &      +a11*a23*a34*a42
     &      +a11*a24*a32*a43-a11*a24*a33*a42-a12*a21*a33*a44
     &      +a12*a21*a34*a43
     &      +a12*a23*a31*a44-a12*a23*a34*a41-a12*a24*a31*a43
     &      +a12*a24*a33*a41
     &      +a13*a21*a32*a44-a13*a21*a34*a42-a13*a22*a31*a44
     &      +a13*a22*a34*a41
     &      +a13*a24*a31*a42-a13*a24*a32*a41-a14*a21*a32*a43
     &      +a14*a21*a33*a42
     &      +a14*a22*a31*a43-a14*a22*a33*a41-a14*a23*a31*a42
     &      +a14*a23*a32*a41
         vl=dabs(det)/6d0
         if(vl.le.0d0)then
           write(*,*)'vl',sngl(vl)
           write(*,*)'xx',sngl(x(1)),sngl(x(2)),sngl(x(3)),sngl(x(4))
           write(*,*)'yy',sngl(y(1)),sngl(y(2)),sngl(y(3)),sngl(y(4))
           write(*,*)'zz',sngl(z(1)),sngl(z(2)),sngl(z(3)),sngl(z(4))
           pause
         endif
c
c        無次元数および物性値を求める
c
         sch=sc(id)     !1/Pe
         rey=re(id)     !1/Re
         pe=sch*rey     !1/Pe
c         vis=(calvis(t0(ne(im,1)))+calvis(t0(ne(im,2)))     !無次元粘度:μ-=μ(T)/μ0
c     &       +calvis(t0(ne(im,3)))+calvis(t0(ne(im,4))))
c     &       /4d0/calvis(-999d0)
c         den=(calden(t0(ne(im,1)))+calden(t0(ne(im,2)))     !無次元密度:ρ-=ρ(T)/ρ0
c     &       +calden(t0(ne(im,3)))+calden(t0(ne(im,4))))
c     &       /4d0/calden(-999d0)
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'im,id',im,id
c         write(*,*)'dt,vl',sngl(dt),sngl(vl)
c         write(*,*)'usum,vsum,wsum',sngl(usum),sngl(vsum),sngl(wsum)
c         write(*,*)'b(i)',sngl(b(1)),sngl(b(2)),sngl(b(3)),sngl(b(4))
c         write(*,*)'c(i)',sngl(c(1)),sngl(c(2)),sngl(c(3)),sngl(c(4))
c         write(*,*)'d(i)',sngl(d(1)),sngl(d(2)),sngl(d(3)),sngl(d(4))
c         write(*,*)'sch,rey',sngl(sch),sngl(rey)
c         write(*,*)'vis,den',sngl(vis),sngl(den)
c         write(*,*)'---------------------------------------------------'
c         pause
c
         do 290 i=1,4     !要素imの行
           ic=ne(im,i)
           do 290 j=1,4   !要素imの列
             jc=ne(im,j)
c             
             cmm=vl/20d0*zkesu(i,j)     !係数行列:[C]
             sxx=b(i)*b(j)/36d0/vl      !係数行列:[Sxx]
             syy=c(i)*c(j)/36d0/vl      !係数行列:[Syy]
             szz=d(i)*d(j)/36d0/vl      !係数行列:[Szz]
c             sxy=b(i)*c(j)/36d0/vl      !係数行列:[Sxy]
c             syx=c(i)*b(j)/36d0/vl      !係数行列:[Syx]
c             syz=c(i)*d(j)/36d0/vl      !係数行列:[Syz]
c             szy=d(i)*c(j)/36d0/vl      !係数行列:[Szy]
c             szx=d(i)*b(j)/36d0/vl      !係数行列:[Szx]
c             sxz=b(i)*d(j)/36d0/vl      !係数行列:[Sxz]
c             chxij=b(j)/24d0            !係数行列:[Hx]
c             chxji=b(i)/24d0            !係数行列:[Hx]T
c             chyij=c(j)/24d0            !係数行列:[Hy]
c             chyji=c(i)/24d0            !係数行列:[Hy]T
c             chzij=d(j)/24d0            !係数行列:[Hz]
c             chzji=d(i)/24d0            !係数行列:[Hz]T
             cxx=(uu(jpc(jmc(ne(im,i))))+usum)*b(j)/120d0  !係数行列:[Cxx]←対流項
             cyy=(vv(jpc(jmc(ne(im,i))))+vsum)*c(j)/120d0  !係数行列:[Cyy]←対流項
             czz=(ww(jpc(jmc(ne(im,i))))+wsum)*d(j)/120d0  !係数行列:[Czz]←対流項
c             
             a11=cmm                    !時間の偏微分項
c
             b11=cxx+cyy+czz            !対流項
             b12=pe*(sxx+syy+szz)       !拡散項
c
c         write(*,*)'---------------------------------------------------'
c         write(*,*)'a11',sngl(a11)
c         write(*,*)'b1i',sngl(b11),sngl(b12)
c         write(*,*)'---------------------------------------------------'
c         pause
c
c            物質移動収支式
             sa(ic,nbw+jc-ic)=sa(ic,nbw+jc-ic)+a11/dt+b11+b12
             sf(ic)=sf(ic)+a11/dt*c0(ne(im,j))
290      continue
200    continue
c
       return
       end
c-----------------------------------------------------------------------------
c
c      物質移動収支式のマトリクスを解く前に第一種の境界条件を代入する    sa(ind,ibw),sf(ind)
c
c-----------------------------------------------------------------------------
       subroutine boundm(np,ne,idx,sa,sf,nbw,ipc,jwc,nn,mn)      !np:新_新
       parameter(ind=19999,ine=99999,ibw=1999)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sa(4*ind,ibw),sf(4*ind),idx(0:ine),js(4,-1:3)   !ne(kele,4):古
       dimension ipc(ind),jwc(ind),mn(ind),mm(ind)   !境界条件の重複を避ける(左辺に移項の際，考慮している行のfaiのみ値が更新される．他の行のsaは0)
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
       common/pmass/sc(1001)
       common/pheat/pr(1001),znu(1001,1001)
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
       common/bnd1/cb(49999),tb(49999)     !第一種の境界条件
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
       do i=1,nn        !nn:境界条件を代入する節点の数
         mm(i)=0        !0:境界条件が代入されていない，1:境界条件が代入されている
       enddo           !mn(nn):境界条件を代入する節点の番号
c
       do 1000 i=1,nb(2)      
         kp1=ne(nelb(2,i),js(neb(2,i),1))  !nelb(2,nb(2)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(2,i),js(neb(2,i),2))  !neb(2,nb(2)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(2,i),js(neb(2,i),3))  !kp1,kp2,kp3:古
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
         kp1=jwc(ipc(kp1))       !jwc(ipc(古)):新_新
         kp2=jwc(ipc(kp2))
         kp3=jwc(ipc(kp3))
         fai=cb(i)
         if(fai.eq.-999d0)goto 1000    !次の面を考慮
c-----------------------------------------------
         if(fai.lt.0d0)then
           write(*,*)'fai in boundm',fai
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
c       do 1400 i=1,nb(10)
c         kele=nelb(10,i)
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
c-----------------------------------------------------------------------------
c
c      物質移動収支式のマトリクスを解いた後に第一種の境界条件を代入する    sa(ind,ibw),sf(ind)
c
c-----------------------------------------------------------------------------
       subroutine boundmr(ne,sf,ipc,jwc,nn,mn,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),sf(4*ind)     !ne(kele,4):古
       dimension ipc(ind),jwc(ind),mn(ind),mm(ind)     !mm(np):判定
       dimension js(4,-1:3)
       common /prop/re(1001),st(-1:1001,-1:1001),sv(1001),gv(1001)
     &                                                        ,vv0,xx0
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
       common/bnd1/cb(49999),tb(49999)     !第一種の境界条件
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
       do 1000 i=1,nb(2)      
         kp1=ne(nelb(2,i),js(neb(2,i),1))  !nelb(2,nb(2)):第一種の境界条件を代入する要素の番号
         kp2=ne(nelb(2,i),js(neb(2,i),2))  !neb(2,nb(2)):第一種の境界条件を代入する要素の面の番号
         kp3=ne(nelb(2,i),js(neb(2,i),3))
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
         kp1=jwc(ipc(kp1))       !jwt(ipc(古)):新_新
         kp2=jwc(ipc(kp2))
         kp3=jwc(ipc(kp3))
c
         fai=cb(i)
         if(fai.eq.-999d0)goto 1000    !次の面を考慮
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