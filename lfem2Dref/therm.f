c
c      蒸発速度を求める(3角形2次要素)
c
       subroutine evpr(ne2,tp,cc0,flx,ine,ind)
	 implicit double precision(a-h,o-z)
	 dimension ne2(ine,6),flx(1001),tp(ind),cc0(ind)
       common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                   ncomw(21,1001),ncopw(21,1001)
	 lr=2
       ee=1d-2                 !ε[-]:蒸発係数
       pi=3.1415926535897932d0 !384626433832795:倍精度の仮数部は10進法で約17桁
       r=8.31451d0             !R[J/K･mol]:気体定数
       zm=fwm(1)               !M[kg/mol]:分子量(1:水)
       tv=tair(1)              !Tv[K]:気相の温度(1:条件1)
       pv=fpv(tv)              !Pv[Pa]:Tvにおける飽和蒸気圧
       cini=cntini(1)          !初期濃度cA0 [mol/m3]
c
       do 300 i=1,nbw(lr),2    !1要素づつ
         if(ncomw(lr,i).ne.0)goto 300     !自由表面のみ考慮
	   flx(i)=0d0
         flx(i+1)=0d0
         flx(i+2)=0d0
c
         k1=ne2(nelbw(lr,i),nebw(lr,i))      !2次要素内の節点番号1-3の節点の番号は1次要素の節点番号と同じ
         k3=ne2(nelbw(lr,i+2),nebw(lr,i+2))
	   tf1=tp(k1)*tplate(1)  !Tf[K]:液表面温度
         tf3=tp(k3)*tplate(1)
         tf2=(tf1+tf3)/2d0
         pf1=fpv(tf1)          !Pf[Pa]:Tvにおける飽和蒸気圧
         pf2=fpv(tf2)
         pf3=fpv(tf3)
c         cf1=cini*(cc0(k1)+1d0)  !cA[mol/m3]:初期濃度
c         cf3=cini*(cc0(k3)+1d0)
c         cf2=(cf1+cf2)/2d0
	   fl1=ee*dsqrt(zm/(2d0*pi*r))*
     &          (pf1/dsqrt(tf1)-pv/dsqrt(tv))
	   fl2=ee*dsqrt(zm/(2d0*pi*r))*
     &          (pf2/dsqrt(tf2)-pv/dsqrt(tv))
	   fl3=ee*dsqrt(zm/(2d0*pi*r))*
     &          (pf3/dsqrt(tf3)-pv/dsqrt(tv))
	   fl1=dmax1(fl1,0d0)
	   fl2=dmax1(fl2,0d0)
	   fl3=dmax1(fl3,0d0)
c
	   if(fl1.eq.0d0)then
	     flx(i)=fl1
	   else
	     flx(i)=flx(i)+(fl1-flx(i))*0.6d0
	   endif
c
	   if(fl2.eq.0d0)then
	     flx(i+1)=fl2
	   else
	     flx(i+1)=flx(i+1)+(fl2-flx(i+1))*0.6d0
	   endif
c
	   if(fl3.eq.0d0)then
	     flx(i+2)=fl3
	   else
	     flx(i+2)=flx(i+2)+(fl3-flx(i+2))*0.6d0
	   endif
c	   flx(i)=2.5d-2         !実験から求められた蒸発速度の100倍
 300	 continue
	 return
	 end
c
c      蒸発速度を求める(3角形1次要素)
c
       subroutine evpr2(ne,flx,tp0,cc0,ine,ind)
	 implicit double precision(a-h,o-z)
	 dimension ne(ine,3),flx(1001),tp0(ind),cc0(ind)
       common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
	 lr=2                    !領域2
       ee=1d-2                 !ε[-]:蒸発係数
       pi=3.1415926535897932d0 !384626433832795:倍精度の仮数部は10進法で約17桁
       r=8.31451d0             !R[J/K･mol]:気体定数
       zm=fwm(1)               !M[kg/mol]:分子量(1:水)
       tv=tair(1)              !Tv[K]:気相温度(1:条件1)
       pv=fpv(tv)              !Pv[Pa]:Tfにおける飽和蒸気圧
       cini=cntini(1)          !初期濃度cA0 [mol/m3]
c
       do 300 i=1,nb(lr)
         if(ncom(lr,i).ne.0)goto 300     !自由表面のみ考慮
	   flx(i)=0d0
         flx(i+1)=0d0
c
	   k1=ne(nelb(lr,i),neb(lr,i))
	   k2=ne(nelb(lr,i+1),neb(lr,i+1))
	   tf1=tp0(k1)*tplate(1)    !Tv[K]:液表面温度
	   tf2=tp0(k2)*tplate(1)    !Tv[K]:液表面温度
         pf1=fpv(tf1)             !Pv[Pa]:Tvにおける飽和蒸気圧
         pf2=fpv(tf2)             !Pv[Pa]:Tvにおける飽和蒸気圧
c         cf1=cini*(cc0(k1)+1d0)     !cA[mol/m3]:初期濃度
c         cf3=cini*(cc0(k3)+1d0)     !cA[mol/m3]:初期濃度
	   fl1=ee*dsqrt(zm/(2d0*pi*r))*(pf1/dsqrt(tf1)-pv/dsqrt(tv))
	   fl2=ee*dsqrt(zm/(2d0*pi*r))*(pf2/dsqrt(tf2)-pv/dsqrt(tv))
	   fl1=dmax1(fl1,0d0)
	   fl2=dmax1(fl2,0d0)
c
	   if(fl1.eq.0d0)then
	     flx(i)=fl1
	   else
	     flx(i)=flx(i)+(fl1-flx(i))*0.6d0
	   endif
c
	   if(fl2.eq.0d0)then
	     flx(i+1)=fl2
	   else
	     flx(i+1)=flx(i+1)+(fl2-flx(i+1))*0.6d0
	   endif
c	   flx(i)=2.5d-2       !実験から求められた蒸発速度の100倍
c	   flx(i)=1d-4
c	   flx(i+1)=1d-4
300	 continue
	 return
	 end
c
c      熱伝導度k[J/s･m･K]を求める
c
       function fk(tt,cc)
       implicit double precision(a-h,o-z)
c       tw=tplate(1)            !Tw[K]:平板温度
c       cini=cntini(1)          !初期濃度cA0[mol/m3]
c       tmp=tt*tw               !温度T[K]
c       cnt=cini*(cc+1d0)       !濃度cA[mol/m3]
c       fk=fk(tmp,cnt)とする
       fk0=0.602d0    !代表の熱伝導度k0
       fk=0.602d0     !熱伝導度kを温度ttと濃度ccの関数で与える
       if(tt.eq.-999d0 .and. cc.eq.-999d0)fk=fk0
       return
       end
c
c      蒸発潜熱L[J/kg]を求める
c
       function hvap(tt,cc)
       implicit double precision(a-h,o-z)
c       tw=tplate(1)            !Tw[K]:平板温度
c       cini=cntini(1)          !初期濃度cA0[mol/m3]
c       tmp=tt*tw               !温度T[K]
c       cnt=cini*(cc+1d0)       !濃度cA[mol/m3]
c      hvap=hvap(tmp,cnt)とする
       hvap0=44.016d3/18.016d-3   !1:水[J/kg]
       hvap=44.016d3/18.016d-3   !1:水[J/kg]
       if(tt.eq.-999d0 .and. cc.eq.-999d0)hvap=hvap0
c       if(n.eq.2)  hvap=44.016d3             !2:水[J/mol]
       return
       end
c
c      膨張率α[1/K]を求める
c
       function vratio(n)
       implicit double precision(a-h,o-z)
       if(n.eq.1) vtario=2.1d-4  !1:水
       return
       end
c
c      定容比熱Cp[J/K･kg]を求める
c
       function dcp(tt,cc)
	 implicit double precision(a-h,o-z)
c       tw=tplate(1)            !Tw[K]:平板温度
c       cini=cntini(1)          !初期濃度cA0[mol/m3]
c       tmp=tt*tw               !温度T[K]
c       cnt=cini*(cc+1d0)       !濃度cA[mol/m3]
c       dcp=dcp(tmp,cnt)とする
       dcp0=75.29d0/18.016d-3        !代表の比熱:水 [J/K･mol]/[g/mol]
       dcp=75.29d0/18.016d-3    !比熱Cpを温度ttと濃度ccの関数で与える
       if(tt.eq.-999d0 .and. cc.eq.-999d0)dcp=dcp0
       return
	 end
c
c      境膜伝熱係数h[W/m^2･K]を求める
c
       function zheat(n)
       implicit double precision(a-h,o-z)
       if(n.eq.1) zheat=(3d0+10d0)/2d0   !1:空気(自然対流)
       return
       end
c
c      粘度μ[Pa･s]を求める
c
       function viscos(tt,cc)
       implicit double precision(a-h,o-z)
c       tw=tplate(1)            !Tw[K]:平板温度
c       cini=cntini(1)          !初期濃度cA0[mol/m3]
c       tmp=tt*tw               !温度T[K]
c       cnt=cini*(cc+1d0)       !濃度cA[mol/m3]
c       viscos=viscos(tmp,cnt)とする
       viscos0=100.2d-5   !代表の粘度:水(20℃)
       viscos=100.2d-5    !粘度μを温度ttと濃度ccの関数で与える
       if(tt.eq.-999d0 .and. cc.eq.-999d0)viscos=viscos0
       return
       end
c
c      密度ρ[kg/m^3]を求める
c
       function density(tt,cc)
       implicit double precision(a-h,o-z)
c       tw=tplate(1)            !Tw[K]:平板温度
c       cini=cntini(1)          !初期濃度cA0[mol/m3]
c       tmp=tt*tw               !温度T[K]
c       cnt=cini*(cc+1d0)       !濃度cA[mol/m3]
c       density=density(tmp,cnt)とする
       density0=998.2d0     !代表の密度:水(20℃)
       density=998.2d0      !密度ρを温度ttと濃度ccの関数で与える
       if(tt.eq.-999d0 .and. cc.eq.-999d0)density=density0
       return
       end
c
c      分子量M[kg/mol]を求める
c
       function fwm(n)
	 implicit double precision(a-h,o-z)
       if(n.eq.1)fwm=18.016d-3         !1:水
c       if(n.eq.2)fwm=319.85d-3         !Methyleue Blue hydrate
       return
	 end
c
c      飽和蒸気圧p[Pa]を求める
c
       function fpv(t)
	 implicit double precision(a-h,o-z)
c      固定点(p*,t*)として，三重点を取る
       pstar=6.11d2       !p*[Pa]
       tstar=273.16d0     !t*[K]
c
       zhvap=44.016d3     !2:水[J/mol]
c
       r=8.31451d0        !気体定数:R[J/mol･K]
c
       c=zhvap/r*(1d0/t-1d0/tstar)
       fpv=pstar*dexp(-c)
	 return
	 end
c
       function fsound(id)
	 implicit double precision(a-h,o-z)
	 if(id.eq.1)then     !in liquid
	   fsound=1150d0
	   return
	 elseif(id.eq.3)then  !in gaseous
	   fsound=340d0
	   return
	 endif
	   write(*,*)"[fsound] id=",id
	 stop
	 end
c-----------------------------------------------
c                    境界条件
c-----------------------------------------------
c
c      液滴内初期温度[℃]を与える
c
       function tini(n)
       implicit double precision(a-h,o-z)
       if(n.eq.1)tini=25d0+273d0      !条件1
       return
       end
c
c      平板温度[℃]を与える
c
       function tplate(n)
       implicit double precision(a-h,o-z)
c       if(n.eq.1)tplate=25d0+273d0    !条件1
       if(n.eq.1)tplate=30d0+273d0    !条件1
       return
       end
c
c      気相温度[℃]を与える
c
       function tair(n)
       implicit double precision(a-h,o-z)
c       if(n.eq.1)tair=25d0+273d0      !条件1
       if(n.eq.1)tair=10d0+273d0      !条件1
       return
       end
c
c      液滴内初期濃度cA0 [mol/m3]
c      (蒸発速度を求める際にのみ用いる)
       function cntini(n)
       implicit double precision(a-h,o-z)
       if(n.eq.1)cntini=0d0        !条件1
       return
       end


