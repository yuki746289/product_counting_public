       implicit double precision(a-h,o-z)
       parameter(ind=99999,ine=99999,ibw=9999)
       dimension ne(ine,4),new(ine,4),nen(ine,4)
       dimension xx(ind),yy(ind),zz(ind)
       dimension xn(ind),yn(ind),zn(ind)
       dimension idx(0:ine),idn(0:ine)
       dimension ip(ind),jp(ind),ie(ine),je(ine)
       dimension jww(ind),jmm(ind)
       dimension pp(ind),uu(ind),vv(ind),ww(ind)
       dimension p0(ind),u0(ind),v0(ind),w0(ind)
       dimension pn0(ind),un0(ind),vn0(ind),wn0(ind)
       dimension sa(4*ind,ibw),sf(4*ind)
c
       do i=0,nele
       idx(i)=0     !idx(0:nele)
       enddo
c       
       nelew=0
       npw=0
       do 1000 i=1,np
         ip(i)=0
         jp(i)=0
1000   continue
c
c      流動領域の節点に新しい番号を割り振る
c
       do 1100 kele=1,nele  !古
         if(idx(kele).eq.1 .or. idx(kele).eq.2 .or.  !index番号が流動領域の場合
     &      idx(kele).eq.3 .or. idx(kele).eq.4)then
           nelew=nelew+1
           do 1200 i=1,4
             if(ip(ne(kele,i)).eq.0)then
               npw=npw+1
               new(nelew,i)=npw
               ip(ne(kele,i))=npw     !ip(古)=新
               jp(npw)=ne(kele,i)     !jp(新)=古
               ie(kele)=nelew         !ie(古)=新
               jp(nelew)=kele         !jp(新)=古
             elseif(ip(ne(kele,i)).ne.0)then
               new(nelew,i)=ip(ne(kele,i))
             endif
1200       continue
         endif
1100   continue
c
c      バンド幅最適化
       call renum3(npw,nelew,new,jww,jmm,nbw,ind,ine,ibw)  !jww(新)=新_新
c                                                          !jmm(新_新)=新
       do 1300 kele=0,nelew  !新_新(要素は新=新_新)
         idn(kele)=idx(je(kele))    !je(新_新)=古
         do 1300 i=1,4              !idn(0:nelew)
         nen(kele,i)=jww(new(kele,i))
1300   continue
       do 1400 kp=1,npw     !新_新
         xn(kp)=xx(jp(jmm(kp)))     !jp(新)=古
         yn(kp)=yy(jp(jmm(kp)))     !jmm(新_新)=新
         zn(kp)=zz(jp(jmm(kp)))
         pn0(kp)=p0(jp(jmm(kp)))    !jp(jmm(新_新))=古
         un0(kp)=u0(jp(jmm(kp)))
         vn0(kp)=v0(jp(jmm(kp)))
         wn0(kp)=w0(jp(jmm(kp)))
         tn0(kp)=t0(jp(jmm(kp)))    !物性値を計算するために用いる
1400  continue
      npn=npw
      neln=nelw
c
c     圧力,速度を計算する
c
      do 1500 i=1,npn
        uvwp0(4*i-3)=u0(i)
        uvwp0(4*i-2)=v0(i)
        uvwp0(4*i-1)=w0(i)
        uvwp0(4*i-0)=p0(i)
1500  continue
      call flow(nelew,4*nbw,sa,sf,nen,xn,yn,zn   !nen(ine,4):流動領域の要素のみ格納
     &                     ,uvwp0,tn0,idn,ind,ine,ibw)
      call boundf1(npw,nelew,nen,idn,sa,sf,4*nbw,ind,ine,ibw)
      call gauss(4*npw,4*nbw,sa,sf,ind,ibw)
      call boundf(nen,sf,ind,ine)
c
c     求めた値を配列に代入する
c
      do 1500 kp=1,npw     !新_新
        uu(jp(jmm(kp)))=sf(4*kp-3)     !jp(新)=古
        vv(jp(jmm(kp)))=sf(4*kp-2)     !jmm(新_新)=新
        ww(jp(jmm(kp)))=sf(4*kp-1)
        pp(jp(jmm(kp)))=sf(4*kp-0)     !jp(jmm(新_新))=古
1500  continue