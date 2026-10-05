c
c     マトリックス[sat],[sft]を組み立てる
c
      subroutine temp(nel,nbw,ne,xx,yy,sa,sf,dt,t0
     &                ,pp,p0,c0,idx,jwt,ibw,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine)
      dimension sa(ind,ibw),sf(ind)      
      dimension t0(ind),c0(ind),pp(ind),p0(ind)  !t0:温度,c0:濃度，pp,p0:圧力
      dimension x(3),y(3),b(3),c(3),li(9),lj(9)
      dimension jwt(ind)
      dimension c1(3,3),c2(3,3),c3(3,3)
c      common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
c      common /doper/ vi,di,sta,str
      common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
      data li/1,1,1,2,2,2,3,3,3/ ,lj/1,2,3,1,2,3,1,2,3/
      data c1/6,2,2,2,2,1,2,1,2/
      data c2/2,2,1,2,6,2,1,2,2/
      data c3/2,1,2,1,2,2,2,2,6/

c
c     物性値(有次元)
      alpha=vratio(1)              !膨張率α[1/K]
      capa0=dcp(-999d0,-999d0)     !代表の定容比熱Cp[J/K･kg]
      visc0=viscos(-999d0,-999d0)  !代表の粘度μ0[Pa･s]
      dens0=density(-999d0,-999d0) !代表の密度ρ0[kg/m^3]
      di=1d-3                      !代表長さr0[m]
      zk0=fk(-999d0,-999d0)        !代表の熱伝導度k0[J/s･m･K] 
c     無次元数
      pe=capa0*visc0/zk0           !ペクレ数:Pe=Cp0*μ0/k0
c
      do 200 im=1,nel
        if(idx(im).ne.1)goto 200      !インデックス番号1の領域のみ考慮
        do 210 i=1,3
          x(i)=xx(ne(im,i))
          y(i)=yy(ne(im,i))
210     continue
        b(1)=y(2)-y(3)
        b(2)=y(3)-y(1)
        b(3)=y(1)-y(2)
        c(1)=x(3)-x(2)
        c(2)=x(1)-x(3)
        c(3)=x(2)-x(1)
        ar=((x(1)-x(2))*(y(1)-y(3))-(x(1)-x(3))*(y(1)-y(2)))/2d0 !要素面積
        rav=(x(1)+x(2)+x(3))/3d0                                 !平均半径
        if(ar.le.0) then
          write(*,*)"the value of element's area is negative"
        endif
c       無次元熱伝導度K=k/k0 [-]
        zk=(fk(t0(ne(im,1)),c0(ne(im,1)))
     &     +fk(t0(ne(im,2)),c0(ne(im,2)))
     &     +fk(t0(ne(im,3)),c0(ne(im,3))))/(zk0*3d0)
c       無次元密度ρ-=ρ/ρ0 [-]
        dens=(density(t0(ne(im,1)),c0(ne(im,1)))
     &       +density(t0(ne(im,2)),c0(ne(im,2)))
     &       +density(t0(ne(im,3)),c0(ne(im,3))))/(dens0*3d0)
c       無次元比熱Cp-=Cp/Cp0 [-]
        capa=(dcp(t0(ne(im,1)),c0(ne(im,1)))
     &       +dcp(t0(ne(im,2)),c0(ne(im,2)))
     &       +dcp(t0(ne(im,3)),c0(ne(im,3))))/(capa0*3d0)
        do 290 l=1,9
          il=li(l)
          jl=lj(l)
          ii=jwt(ne(im,il))
          jj=jwt(ne(im,jl))
          a11=(c1(il,jl)*x(1)+c2(il,jl)*x(2)+c3(il,jl)*x(3))*ar/60d0	   ![B]
          b11=1d0/pe*zk/capa/dens*(b(il)*b(jl)+c(il)*c(jl))*rav/4d0/ar   !1/Pe*K/Cp/ρ*([Srr]+[Szz])
          c11=alpha*visc0**2/(capa0*dens0**2*di**2)/capa/dens*a11        !αi*μ0^2/(Cp0*ρ0*r0^2)/Cp/ρ*[B]
          sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11/dt+b11
          sf(ii)=sf(ii)+a11/dt*t0(ne(im,jl))
c     &                 +c11*(pp(ne(im,jl))-p0(ne(im,jl)))/dt
          if(sa(ii,jj-ii+nbw).eq.0) write(*,*)" (temp) sa=0",ii,jj
290     continue
200     continue
c
      return
      end
c
c     自由表面に境界条件を与える(第3種，蒸発線熱)
c
      subroutine boundtmp(nel,nbw,ne,xx,yy,sa,sf,
     &                    t0,c0,idx,jwt,ibw,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,3),nec(ine,3),xx(ind),yy(ind),idx(0:ine)
      dimension sa(ind,ibw),sf(ind)
      dimension t0(ind),c0(ind)		    !t0:温度T,c0:濃度CA
      dimension li(9),lj(9)
      dimension js2(3),jwt(ind)
c-----物性値---------------------------------------------------
      dimension flx(1001)
c--------------------------------------------------------------
c      common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
c      common /doper/ vi,di,sta,str
      common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
      data li/1,1,1,2,2,2,3,3,3/ ,lj/1,2,3,1,2,3,1,2,3/
      data js2/2,3,1/
c
      call getnec(nel,ne,nec,ine)   !nec取得
c
c     物性値(有次元)
      dl0=hvap(-999d0,-999d0)       !蒸発線熱L[J/kg]
      zhh=zheat(1)                  !境膜伝熱係数h[W/m^2･K]
	capa0=dcp(-999d0,-999d0)      !定容比熱Cp[J/K･kg]
      di=1d-3                       !代表長さr0[m]
      visc0=viscos(-999d0,-999d0)   !代表の粘度μ0[Pa･s]
      dens0=density(-999d0,-999d0)  !代表の密度ρ0[kg/m^3]
      zk0=fk(-999d0,-999d0)         !代表の熱伝導度k0[J/s･m･K] 
      call evpr2(ne,flx,t0,c0,ine,ind) !蒸発速度m[kg/m^2･s]を求める
c     境界温度
      ttw=tplate(1)                 !平板の温度Tw[℃]
      tb=tair(1)/tplate(1)          !無次元の気相温度Θb[-]
c     無次元数
      pe=capa0*visc0/zk0            !ペクレ数:Pe=Cp0*μ0/k0
      dnu=zhh*di/zk0                !ヌッセルト数:Nu=h*r0/k0
c     自由表面の境界条件
      open(1,file="boundtmp.res")
      open(2,file="boundxy.res")
      rewind(1)
      rewind(2)
      lr=2   !領域番号
      do 400 i=1,nb(lr)
        if(idx(nec(nelb(lr,i),neb(lr,i))).ne.0)goto 400    !自由表面のみ考慮
        im=nelb(lr,i)
        rav=(xx(ne(nelb(lr,i),neb(lr,i)))                  !rav:平均半径
     &      +xx(ne(nelb(lr,i+1),neb(lr,i+1))))/2d0
        rl=dsqrt((xx(ne(nelb(lr,i),neb(lr,i)))             !rl:1辺の長さ
     &           -xx(ne(nelb(lr,i+1),neb(lr,i+1))))**2
     &          +(yy(ne(nelb(lr,i),neb(lr,i)))
     &           -yy(ne(nelb(lr,i+1),neb(lr,i+1))))**2)
c       無次元熱伝導度K=k/k0 [-]
        zk=(fk(t0(ne(im,1)),c0(ne(im,1)))
     &     +fk(t0(ne(im,2)),c0(ne(im,2)))
     &     +fk(t0(ne(im,3)),c0(ne(im,3))))/(zk0*3d0)
c       無次元密度ρ-=ρ/ρ0 [-]
        dens=(density(t0(ne(im,1)),c0(ne(im,1)))
     &       +density(t0(ne(im,2)),c0(ne(im,2)))
     &       +density(t0(ne(im,3)),c0(ne(im,3))))/(dens0*3d0)
c       無次元比熱Cp-=Cp/Cp0 [-]
        capa=(dcp(t0(ne(im,1)),c0(ne(im,1)))
     &       +dcp(t0(ne(im,2)),c0(ne(im,2)))
     &       +dcp(t0(ne(im,3)),c0(ne(im,3))))/(capa0*3d0)
        do 450 l=1,9
          il=li(l)     !マトリックスの(il,jl)成分
          jl=lj(l)
          ii=jwt(ne(im,il))
          jj=jwt(ne(im,jl))
          a11=dnu/pe/capa/dens         !Nu/Pe/Cp/ρ
          if(il.eq.neb(lr,i) .and. jl.eq.neb(lr,i))then
            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rl/12d0
     &        *(3d0*xx(ne(nelb(lr,i),neb(lr,i)))
     &         +xx(ne(nelb(lr,i),js2(neb(lr,i)))))
c            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rl*rav*2d0/6d0
          else if(il.eq.js2(neb(lr,i)) .and. jl.eq.js2(neb(lr,i)))then
            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rl/12d0
     &        *(xx(ne(nelb(lr,i),neb(lr,i)))
     &         +3d0*xx(ne(nelb(lr,i),js2(neb(lr,i)))))
c            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rav*2d0*rl/6d0
          else if(il.eq.neb(lr,i) .and. jl.eq.js2(neb(lr,i)) .or. 
     &            jl.eq.neb(lr,i) .and. il.eq.js2(neb(lr,i)))then
            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rl/12d0
     &        *(xx(ne(nelb(lr,i),neb(lr,i)))
     &         +xx(ne(nelb(lr,i),js2(neb(lr,i)))))
c            sa(ii,jj-ii+nbw)=sa(ii,jj-ii+nbw)+a11*rav*rl/6d0
          endif
450     continue
c       書き込み
        write(1,*)"-------------------------"
        write(1,*)ne(nelb(lr,i),neb(lr,i))
        write(1,*)"rav,rl",sngl(rav),sngl(rl)
        write(1,*)"zk,dens,capa",sngl(zk),sngl(dens),sngl(capa)
        write(1,*)"tb,ttw,di",sngl(tb),sngl(ttw),sngl(di)
        write(1,*)"pe,dnu",sngl(pe),sngl(dnu)
        do l=1,9
          il=li(l)
          jl=lj(l)
          ii=jwt(ne(im,il))
          jj=jwt(ne(im,jl))
          if(il.eq.neb(lr,i) .and. jl.eq.neb(lr,i))then
            write(1,*)"--",l,ne(im,il),ne(im,jl),sngl(a11*rl/12d0
     &        *(3d0*xx(ne(nelb(lr,i),neb(lr,i)))
     &         +xx(ne(nelb(lr,i),js2(neb(lr,i))))))
          else if(il.eq.js2(neb(lr,i)) .and. jl.eq.js2(neb(lr,i)))then
            write(1,*)"--",l,ne(im,il),ne(im,jl),sngl(a11*rl/12d0
     &        *(xx(ne(nelb(lr,i),neb(lr,i)))
     &         +3d0*xx(ne(nelb(lr,i),js2(neb(lr,i))))))
          else if(il.eq.neb(lr,i) .and. jl.eq.js2(neb(lr,i)) .or. 
     &            jl.eq.neb(lr,i) .and. il.eq.js2(neb(lr,i)))then
            write(1,*)"--",l,ne(im,il),ne(im,jl),sngl(a11*rl/12d0
     &        *(xx(ne(nelb(lr,i),neb(lr,i)))
     &         +xx(ne(nelb(lr,i),js2(neb(lr,i))))))
          endif
          write(1,*)"sa",l,ne(im,il),ne(im,jl),sngl(sa(ii,jj-ii+nbw))
        enddo
c
c       無次元比熱L-=L/L0 [-]
        dl=(hvap(t0(ne(im,1)),c0(ne(im,1)))
     &     +hvap(t0(ne(im,2)),c0(ne(im,2)))
     &     +hvap(t0(ne(im,3)),c0(ne(im,3))))/(dl0*3d0)
c       無次元蒸発速度M･=r0*m･/μ0 [-]
        zm=di*(flx(i)+flx(i+1))/2d0/visc0
c
        d11=1d0/pe/capa/dens*(dnu*tb-(visc0*dl0/zk0/ttw)*zk*zm*dl)  !1/Pe/Cp/ρ*(Nu*Θb-(μ0*L0/k0/Tw)*K*M･*L)
        sf(jwt(ne(nelb(lr,i),neb(lr,i))))=
     &          sf(jwt(ne(nelb(lr,i),neb(lr,i))))+d11*rl/6d0
     &                         *(2d0*xx(ne(nelb(lr,i),neb(lr,i)))
     &                         +xx(ne(nelb(lr,i+1),neb(lr,i+1))))
c        sf(jwt(ne(nelb(lr,i),neb(lr,i))))=
c     &           sf(jwt(ne(nelb(lr,i),neb(lr,i))))+d11*rav*rl/2d0
        sf(jwt(ne(nelb(lr,i+1),neb(lr,i+1))))=
     &          sf(jwt(ne(nelb(lr,i+1),neb(lr,i+1))))+d11*rl/6d0
     &                         *(xx(ne(nelb(lr,i),neb(lr,i)))
     &                         +2d0*xx(ne(nelb(lr,i+1),neb(lr,i+1))))
c        sf(jwt(ne(nelb(lr,i+1),neb(lr,i+1))))=
c     &           sf(jwt(ne(nelb(lr,i+1),neb(lr,i+1))))+d11*rav*rl/2d0
c       書き込み
        write(1,*)"--",ne(nelb(lr,i),neb(lr,i)),
     &            sngl(xx(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(yy(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(d11*rl/6d0*(2d0*xx(ne(nelb(lr,i),neb(lr,i)))
     &                           +xx(ne(nelb(lr,i+1),neb(lr,i+1)))))
        write(1,*)"--",ne(nelb(lr,i),neb(lr,i)),
     &            sngl(xx(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(yy(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(d11*rl/6d0*(xx(ne(nelb(lr,i),neb(lr,i)))
     &                           +2d0*xx(ne(nelb(lr,i+1),neb(lr,i+1)))))
c
        write(1,*)"sf",ne(nelb(lr,i),neb(lr,i)),
     &            sngl(xx(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(yy(ne(nelb(lr,i),neb(lr,i)))),
     &            sngl(sf(jwt(ne(nelb(lr,i),neb(lr,i)))))
        write(1,*)"sf",ne(nelb(lr,i),js2(neb(lr,i))),
     &            sngl(xx(ne(nelb(lr,i),js2(neb(lr,i))))),
     &            sngl(yy(ne(nelb(lr,i),js2(neb(lr,i))))),
     &            sngl(sf(jwt(ne(nelb(lr,i),js2(neb(lr,i))))))
c
        write(2,*)"sf",ne(nelb(lr,i),neb(lr,i)),",",
     &            sngl(xx(ne(nelb(lr,i),neb(lr,i)))),",",
     &            sngl(yy(ne(nelb(lr,i),neb(lr,i))))
        write(2,*)"sf",ne(nelb(lr,i),js2(neb(lr,i))),",",
     &            sngl(xx(ne(nelb(lr,i),js2(neb(lr,i))))),",",
     &            sngl(yy(ne(nelb(lr,i),js2(neb(lr,i)))))
c
400     continue
        close(2)
        close(1)
      return
      end
c
c     ある1要素における熱伝導度の偏微分∂k/∂R,∂k/∂Zを求める(3角形1次要素)
c
      subroutine dfk(zk,x,y,dkr,dkz)
      implicit double precision(a-h,o-z)
      dimension x(3),y(3),zk(3)
c
      ar=((x(1)-x(2))*(y(1)-y(3))-(x(1)-x(3))*(y(1)-y(2)))/2d0
      dkr=1d0/(2d0*ar)
     &    *((y(2)-y(3))*zk(1)+(y(3)-y(1))*zk(2)+(y(1)-y(2))*zk(3))
      dkz=1d0/(2d0*ar)
     &    *((x(3)-x(2))*zk(1)+(x(1)-x(3))*zk(2)+(x(2)-x(1))*zk(3))
      return
      end
c
c     平板との境界面に第1種の境界条件を与える(Θ=1)
c     マトリックスを解く前
      subroutine boundtmp2(nel,npt,nbw,ne,sa,sf,idx,jwt,ibw,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,3),nec(ine,3),idx(0:ine)
      dimension sa(ind,ibw),sf(ind),jwt(ind)
      dimension nbb(10),kb(2001,10),tb(2001,10)   !内部の配列
      common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
      call getnec(nel,ne,nec,ine)
c
      lr=2   !領域番号
      nbb(lr)=0
      do 500 i=1,nb(lr)
        if(idx(nec(nelb(lr,i),neb(lr,i))).eq.0)goto 510    !終了
        nbb(lr)=nbb(lr)+1                            !境界上の節点の数
        kb(nbb(lr),lr)=jwt(ne(nelb(lr,i),neb(lr,i))) !境界上の節点の番号
        tb(nbb(lr),lr)=1d0                           !境界上の節点の温度:Tw/Tw[-]
500   continue
510   continue
        nbb(lr)=nbb(lr)+1                            !境界上の節点の数
        kb(nbb(lr),lr)=jwt(ne(nelb(lr,i),neb(lr,i))) !境界上の節点の番号
        tb(nbb(lr),lr)=1d0                           !境界上の節点の温度:Tw/Tw[-]
c
      nbwm=nbw-1
      lr=2
      do 530 i=1,nbb(lr)
        ii=kb(i,lr)
        aa=sa(ii,nbw)
        do 540 j=max(1,ii-nbwm),min(ii+nbwm,npt)
        sf(j)=sf(j)-sa(j,nbw+ii-j)*tb(i,lr)
        sa(ii,nbw+j-ii)=0d0
 540    sa(j,nbw+ii-j)=0d0
        sa(ii,nbw)=aa
        sf(ii)=aa*tb(i,lr)
 530  continue
c
      return
      end
c
c     平板との境界面に第1種の境界条件を与える(Θ=1)
c     マトリックスを解いた後
      subroutine boundtmp3(nel,ne,tp,idx,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,3),nec(ine,3),tp(ind),idx(0:ine)
      common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
      call getnec(nel,ne,nec,ine)
c
      lr=2   !領域番号
      do 500 i=1,nb(lr)
        if(idx(nec(nelb(lr,i),neb(lr,i))).eq.0)goto 510    !終了
        tp(ne(nelb(lr,i),neb(lr,i)))=1d0            !境界上の節点の温度:Tw/Tw[-]
500   continue
510   continue
        tp(ne(nelb(lr,i),neb(lr,i)))=1d0
      return
      end
c
c     メッシュ再発生後，温度を補間する
c
        subroutine grdfct(nel,np,ne,xx,yy,tt,
     &                    neln,npn,nen,xxn,yyn,ttn,irx,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),irx(0:ine)
        dimension nen(ine,3),xxn(ind),yyn(ind)
        dimension tt(ind),ttn(ind)
        dimension al(3),al0(3)
        dimension tb(1001)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /ard/ arr(8001)
c
      do 80 i=1,nel
        arr(i)=(xx(ne(i,1))-xx(ne(i,3)))*
     &       (yy(ne(i,2))-yy(ne(i,3)))
     &      -(xx(ne(i,2))-xx(ne(i,3)))*
     &       (yy(ne(i,1))-yy(ne(i,3)))
 80   continue
c
      do 100 n=1,npn		  !メッシュ再発生後の節点
      x=xxn(n)
      y=yyn(n)
      almin=1d10
      do 140 i=1,nel		 !メッシュ再発生前の要素
        if(arr(i).lt.1d-30) goto 140
        ifg=0
c1       alm=-1d0
        alm1=1000d0
        alm2=-1000d0
        al(1)=((x-xx(ne(i,3)))*(yy(ne(i,2))-yy(ne(i,3)))
     &      -(xx(ne(i,2))-xx(ne(i,3)))*(y-yy(ne(i,3))))
     &      /arr(i)
        al(2)=(-(x-xx(ne(i,3)))*(yy(ne(i,1))-yy(ne(i,3)))
     &       +(xx(ne(i,1))-xx(ne(i,3)))*(y-yy(ne(i,3))))
     &      /arr(i)
        al(3)=1d0-al(1)-al(2)
        do 152 j=1,3
        alm1=dmin1(alm1,al(j))    !最小面積
152     alm2=dmax1(alm2,al(j))    !最大面積
        do 150  j=1,3             !-5>L1,L2,L3 または 5<L1,L2,L3の時は考慮しない(時間短縮)
        if(al(j).lt.-5d0 .or. al(j).gt.5d0) goto 140
150     continue
        if(alm1.gt.-1d-12 .and. alm2.lt.1d0+1d-12) goto 200  !節点が要素iの中に存在する場合
        alm=dmax1(-alm1,alm2-1d0)   !考慮している要素のずれの最大値
        if(alm.lt.almin) then
          almin=alm
          imin=i
          do 145 j=1,3
 145      al0(j)=al(j)
        endif
 140  continue
      i=imin        !節点に最も近い位置に存在する要素の番号
      do 155 j=1,3
 155  al(j)=al0(j)
c
 200  ttn(n)=0d0
      ttn(n)=al(1)*tt(ne(i,1))+al(2)*tt(ne(i,2))+al(3)*tt(ne(i,3))  !T=N1*T1+N2*T2+N3*T3
c
      open(1,file="grdfct.res",access='append')
      write(1,*)n,",",i,",",sngl(arr(i)),",",sngl(al(1))
     &      ,",",sngl(al(2)),",",sngl(al(3)),",",sngl(ttn(n))
      close(1)
c
      if(ttn(n).lt.0d0) ttn(n)=0d0
 100  continue
c
      nel=neln
      np=npn
      do 300 i=1,nel
      do 300 j=1,3
 300  ne(i,j)=nen(i,j)
      do 310 i=1,np
      xx(i)=xxn(i)
      yy(i)=yyn(i)
 310  tt(i)=ttn(i)
c     他の領域との境界上の節点温度補間
      lr=2
      do i=1,nb(lr)            !領域lrの境界上の温度格納
        tb(i)=tt(ne(nelb(lr,i),neb(lr,i)))
      enddo
      do i=1,nel
        if(irx(i).eq.1)then   !領域1の温度0
        do j=1,3
        tt(ne(i,j))=0d0
        ttn(ne(i,j))=0d0
        enddo
        endif
      enddo
      do i=1,nb(lr)
        tt(ne(nelb(lr,i),neb(2,i)))=tb(i)
        ttn(ne(nelb(lr,i),neb(2,i)))=tb(i)
      enddo
c
      return
      end
c
c     接する要素を格納した配列nec(ine,3)を求める
c
      subroutine getnec(nelt,ne,nec,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,3),nec(ine,3),js2(3)
      data js2/2,3,1/
c
        do 100 i=1,nelt
        do 100 j=1,3
 100    nec(i,j)=0
c
        do 120 i=1,nelt
          do 130 j=1,nelt
            if (j.eq.i) goto 130
            do 140 k1=1,3
              if (nec(j,k1).ne.0) goto 140
              k2=js2(k1)
              do 150 jc=1,3
              if ((ne(i,jc).eq.ne(j,k2)).and.
     &                  (ne(i,js2(jc)).eq.ne(j,k1))) then
                 goto 190
              elseif ((ne(i,jc).eq.ne(j,k1)).and.          !
     &                  (ne(i,js2(jc)).eq.ne(j,k2))) then
                 write(*,*)"sufele/node number given on the",
     &                       " opposite direction in a element"
                 pause
              endif
 150          continue
 140        continue
            goto 130
 190        nec(i,jc)=j
            nec(j,k1)=i
            if(nec(i,1).ne.0.and.nec(i,2).ne.0.and.nec(i,3).ne.0)
     &          goto 120
 130      continue
 120     continue
      return
      end
c
c     全節点温度をファイルに保存する
c
      subroutine avs2dmtmp(nelt,npt,ne,idx,xx,yy,tp,cc,
     &                  ti,rfile,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),idx(0:ine),xx(ind),yy(ind)
        dimension tp(ind),cc(ind)
        dimension flx(1001),nec(ine,3)
        character lpp*8,rfile*30,vfile*30
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c
c       avs fileの作成
c
        write(*,*)"in avs2dmtmp"    !
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
 100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'-2d-tmp'//lpp//'.inp',vfile)
        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        idmin=0
        idmax=-9999
        do 200 i=1,nr
        idmin=min(idmin,ndx(i))
 200    idmax=max(idmax,ndx(i))
c       idmin=min(0,-idmin)
c
        open(3,file=vfile)
         rewind(3)
c        write(3,*) '#',1            !added '#' 
c        write(3,*) '#data_geom'
c        write(3,*) '#step1'
c        write(3,*) npt,nelt,
        write(3,*)  npt,nelt,1,0,0    !
        do 300 i=1,npt
 300    write(3,'(i4,3e12.4)') i,-sngl(xx(i)),0e0,sngl(yy(i))

        do 400 i=1,nelt
 400    write(3,'(2i4,a4,3i5)') i,idx(i),'tri',(ne(i,j),j=1,3)
c       write(3,*) 4,0				!vanished 
        write(3,*) 1,1
        write(3,*) 'temperature,'
        do 500 i=1,npt
 500    write(3,'(i4,10e12.4)')
     &          i,sngl(tp(i))

        close(3)
c
c       液滴の気相との界面における蒸発速度分布を求める
c
        call evpr2(ne,flx,tp,cc,ine,ind)
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 600 i=1,8
 600    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr('vap-rate',lpp//'.res',vfile)
        write(*,*) 'vpa-rate/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
        write(1,*)"time",ti
        call getnec(nelt,ne,nec,ine)  !nec(ine,3)取得(隣り合う要素の番号を得る)
        lr=2
        do 700 i=1,nb(lr)
          if(idx(nec(nelb(lr,i),neb(lr,i))).ne.0)goto 700    !自由表面のみ考慮
          k=ne(nelb(lr,i),neb(lr,i))
          write(1,*)k,",",sngl(xx(k)),",",sngl(yy(k)),",",sngl(flx(i))
c
          if(idx(nec(nelb(lr,i+1),neb(lr,i+1))).ne.0)then    !最後の節点
          k=ne(nelb(lr,i+1),neb(lr,i+1))
          write(1,*)k,",",sngl(xx(k)),",",sngl(yy(k)),",",sngl(flx(i))
          endif
 700    continue
        close(1)        
        write(*,*)"logging ended in avs2dmtmp"   !
c
        return
        end



