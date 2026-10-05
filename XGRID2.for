c----------------------------------------------------------------
c
c      グリッドファイル,'.1','.2',作成関数
c      最大節点数:9999999, 最大要素数:9999999
c
c----------------------------------------------------------------
c      引数
c      lrg:分割数を変更する領域の1つの領域番号
c      ldx,ldy,ldz:領域lrgの分割数
c      nrg:領域数
c      mp(nrg,8):領域の各頂点の節点番号
c      x(np),y(np),z(np):各頂点の座標
c      mr(nrg):index番号
       include 'addchr.for'
c	 include 'renum3.for'
       implicit double precision(a-h,o-z)
       parameter(ind=99999,ine=99999,ibw=1999,irg=9999)
c       
       dimension xl(12,1001),yl(12,1001),zl(12,1001)         !辺上に存在する節点の座標
       dimension xc(1001),yc(1001),zc(1001)                  !計算した辺上に存在する節点の座標
       dimension xt(3),yt(3),zt(3)                  !辺上に存在する節点の座標
c
       dimension xs1(1001,1001),ys1(1001,1001),zs1(1001,1001)     !面上に存在する節点の座標
       dimension xs2(1001,1001),ys2(1001,1001),zs2(1001,1001)     !面上に存在する節点の座標
       dimension xs3(1001,1001),ys3(1001,1001),zs3(1001,1001)     !面上に存在する節点の座標
       dimension xs4(1001,1001),ys4(1001,1001),zs4(1001,1001)     !面上に存在する節点の座標
       dimension xs5(1001,1001),ys5(1001,1001),zs5(1001,1001)     !面上に存在する節点の座標
       dimension xs6(1001,1001),ys6(1001,1001),zs6(1001,1001)     !面上に存在する節点の座標
c
       integer*8 iwr(ind),iwx(ind),iwy(ind),iwz(ind)              !iwr(nn),iwx(nn),iwy(nn),iwz(nn)
       dimension xr(ind),yr(ind),zr(ind)                          !各領域内を構成する6面体の座標
c
       dimension ifline(irg,12)                    !0:直線,1:曲線
       dimension mdx(1001),mdy(1001),mdz(1001)     !各領域の分割数
       dimension mk(1001,6)                        !領域の各面を共有している領域の番号
       dimension mp(1001,8)                        !領域の各頂点の節点番号
       dimension xp(1001,12,3),yp(1001,12,3),zp(1001,12,3)  !辺を構成する節点の座標
       dimension xv(1001,8),yv(1001,8),zv(1001,8)  !各6面体の頂点座標
c
       dimension mm(1001,12)                       !節点番号の割り当てを示す
                                                    !0:割り振られている,1:割り振られていない
       dimension ms(1001,6)                        !6面体の面上に対する節点番号の割り当てを判別する
                                                    !0:割り当てられていない,1:割り当てられている
       dimension lli(ind),llj(ind)                 !辺の片端の節点番号
       dimension llk(ind)                          !辺に対する節点番号の割り当てを示す.
                                                    !0:割り当てられている,1:割り当てられていない
       dimension lm(ind)                           !割り当てる節点番号を格納した配列
c
       integer*8 lwr(ind),lwx(ind),lwy(ind),lwz(ind),lwp(ind)     !lwr(np),lwx(np),lwy(np),lwz(np)
       dimension xx(ind),yy(ind),zz(ind)           !各節点の座標
       dimension mr(irg),md(ine)                   !各領域および6面体のindex番号
       dimension ne(ine,8)                         !各6面体を構成する節点の番号
c
       dimension idx(ine)                          !各要素のindex番号
       dimension new(ine,4)                        !各4面体を構成する節点の番号
       dimension jww(ind),jmm(ind)                 !jww(古)=新,jmm(古)=新
c
       character*32 rfile,pfile,gfile1,gfile2      !file names
c
c
       dimension x(4),y(4),z(4)
c------bd.infファイルの作成---------------------------------------
       dimension ubt(ind),vbt(ind),wbt(ind),pbt(ind),cbt(ind),tbt(ind) !第一種の境界条件
	 dimension ugt(ind),vgt(ind),wgt(ind),cgt(ind),tgt(ind)			 !第ニ種の境界条件
	 dimension cht(ind),tht(ind)									 !第三種の境界条件
	 dimension scvt(ind),scmt(ind),sctt(ind)						 !生成項
c
       dimension ub(ine,6),vb(ine,6),wb(ine,6),pb(ine,6)               !第一種の境界条件
	 dimension cb(ine,6),tb(ine,6)									 !第一種の境界条件
	 dimension ug(ine,6),vg(ine,6),wg(ine,6),cg(ine,6),tg(ine,6)	 !第ニ種の境界条件
	 dimension ch(ine,6),th(ine,6)									 !第三種の境界条件
	 dimension scv(ine,6),scm(ine,6),sct(ine,6)						 !生成項
c
	 dimension mbr(ind),mbf(ind),mbn(ind),faib(ind)					 !faib(nn):物理量φ
	 dimension kbr(17,49999),kbf(17,49999)
	 dimension n(17)												 !境界条件の種類
c
       dimension kb(11),kelb(11,49999),keb(11,49999)
	 dimension nb(11),nelb(11,49999),neb(11,49999)
	 dimension mg(ine)			!mg(nele):6面体要素neleが存在する領域の番号krg
	 dimension mgr(ine)			!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
	 dimension kmp(ind)			!kmp(kn):X-Z平面状に存在する節点の番号を格納した配列
	 dimension nf(4)
c-----------------------------------------------------------------
c
c      infファイルを読み込む
c
       open(1,file='input.txt')
       rewind(1)
       read(1,*)rfile
       close(1)
       call addchr(rfile,'.inf',pfile)
       write(*,*)'the input file'
       write(*,*)pfile
       pause
c
       open(1,file=pfile)
       rewind(1)
       read(1,*)npr    !頂点の数
       read(1,*)nrg    !領域数
       do 100 krg=1,nrg
         read(1,*)mr(krg)  !領域krgのインデックス番号
         read(1,*)mdx(krg),mdy(krg),mdz(krg)           !各軸方向の分割数
         do 200 il=1,12
           read(1,*)ifline(krg,il)                     !辺の種類(0:直線,1:曲線)
           nside=2                             !直線
           if(ifline(krg,il).eq.1)nside=3      !曲線
           read(1,*)(xp(krg,il,k),k=1,nside)  !各頂点のx座標
           read(1,*)(yp(krg,il,k),k=1,nside)  !各頂点のy座標
           read(1,*)(zp(krg,il,k),k=1,nside)  !各頂点のz座標
c--------------------------------------------------------------------------------
c           write(*,*)'krg il',krg,il
c           write(*,*)ifline(krg,il)            !辺の種類(0:直線,1:曲線)
c           write(*,*)(xp(krg,il,k),k=1,nside)  !各頂点のx座標
c           write(*,*)(yp(krg,il,k),k=1,nside)  !各頂点のy座標
c           write(*,*)(zp(krg,il,k),k=1,nside)  !各頂点のz座標
c           pause
c--------------------------------------------------------------------------------
200      continue
100    continue
       do 300 krg=1,nrg
300    read(1,*)(mk(krg,i),i=1,6)  !各面が接する領域の番号
       do 400 krg=1,nrg
400    read(1,*)(mp(krg,i),i=1,8)  !各頂点の節点番号
       close(1)
c
c      各6面体の頂点座標配列xv(nrg,8),yv(nrg,8),zv(nrg,8)を求める
c
       do 500 krg=1,nrg
         !節点1
         xv(krg,1)=xp(krg,1,1)
         yv(krg,1)=yp(krg,1,1)
         zv(krg,1)=zp(krg,1,1)
         !節点2
         xv(krg,2)=xp(krg,2,1)
         yv(krg,2)=yp(krg,2,1)
         zv(krg,2)=zp(krg,2,1)
         !節点3
         xv(krg,3)=xp(krg,3,1)
         yv(krg,3)=yp(krg,3,1)
         zv(krg,3)=zp(krg,3,1)
         !節点4
         xv(krg,4)=xp(krg,4,1)
         yv(krg,4)=yp(krg,4,1)
         zv(krg,4)=zp(krg,4,1)
         !節点5
         xv(krg,5)=xp(krg,9,1)
         yv(krg,5)=yp(krg,9,1)
         zv(krg,5)=zp(krg,9,1)
         !節点6
         xv(krg,6)=xp(krg,10,1)
         yv(krg,6)=yp(krg,10,1)
         zv(krg,6)=zp(krg,10,1)
         !節点7
         xv(krg,7)=xp(krg,11,1)
         yv(krg,7)=yp(krg,11,1)
         zv(krg,7)=zp(krg,11,1)
         !節点8
         xv(krg,8)=xp(krg,12,1)
         yv(krg,8)=yp(krg,12,1)
         zv(krg,8)=zp(krg,12,1)
500    continue

c --------------
c       do krg=1,nrg
c       write(*,*)'krg',krg
c       write(*,*)'x',sngl(xv(krg,1)),sngl(xv(krg,2)),sngl(xv(krg,3)),
c     &                sngl(xv(krg,4)),sngl(xv(krg,5)),sngl(xv(krg,6)),
c     &                sngl(xv(krg,7)),sngl(xv(krg,8))
c       write(*,*)'y',sngl(yv(krg,1)),sngl(yv(krg,2)),sngl(yv(krg,3)),
c     &                sngl(yv(krg,4)),sngl(yv(krg,5)),sngl(yv(krg,6)),
c     &                sngl(yv(krg,7)),sngl(yv(krg,8))
c       write(*,*)'z',sngl(zv(krg,1)),sngl(zv(krg,2)),sngl(zv(krg,3)),
c     &                sngl(zv(krg,4)),sngl(zv(krg,5)),sngl(zv(krg,6)),
c     &                sngl(zv(krg,7)),sngl(zv(krg,8))
c       enddo
c -------------

c
c      各領域を構成する節点の座標を求める
c
       nn=0     !nn:節点数
       do 1100 krg=1,nrg
         ndx=mdx(krg)
         ndy=mdy(krg)
         ndz=mdz(krg)
c         
c        辺上の節点座標を求める
c
         do 1200 il=1,12   !辺番号
           nm=nmd(il,ndx,ndy,ndz)     !考慮している辺の分割数
           if(ifline(krg,il).eq.0)then   !直線
c             write(*,*)'in ifline 0'
             xt(1)=xp(krg,il,1)     !辺の始点
             yt(1)=yp(krg,il,1)     !辺の始点
             zt(1)=zp(krg,il,1)     !辺の始点
             xt(2)=xp(krg,il,2)     !辺の終点
             yt(2)=yp(krg,il,2)     !辺の終点
             zt(2)=zp(krg,il,2)     !辺の終点
             call straight(xt,yt,zt,nm,xc,yc,zc)
             do 1300 km=1,nm+1   !考慮している辺上に存在する各節点の座標
               xl(il,km)=xc(km)
               yl(il,km)=yc(km)
               zl(il,km)=zc(km)
c               write(*,*)'krg,nm+1',krg,nm+1
c               write(*,*)'xl(',il,',',km,')',sngl(xl(il,km))
c               write(*,*)'yl(',il,',',km,')',sngl(yl(il,km))
c               write(*,*)'zl(',il,',',km,')',sngl(zl(il,km))
c               pause
1300         continue
           elseif(ifline(krg,il).eq.1)then  !曲線
             do i=1,3
               xt(i)=xp(krg,il,i)
               yt(i)=yp(krg,il,i)
               zt(i)=zp(krg,il,i)
             enddo
             call curve(xt,yt,zt,nm,xc,yc,zc)
             do 1400 km=1,nm+1
               xl(il,km)=xc(km)
               yl(il,km)=yc(km)
               zl(il,km)=zc(km)
1400         continue
           endif
1200     continue
c         write(*,*)"ended deciding the nodes' positions of
c     & each side of hexahedral elements"
c
c        面上の節点を求める
c
         do 1500 kx=1,ndx+1
         do 1500 kz=1,ndz+1
c
c          面1:X-Z平面
c
           xs1(kx,kz)=xl(4,ndz+2-kz)+(xl(2,kz)-xl(4,ndz+2-kz))
     &                              *dble(kx-1)/dble(ndx)
           zs1(kx,kz)=zl(1,kx)+(zl(3,ndx+2-kx)-zl(1,kx))
     &                        *dble(kz-1)/dble(ndz)
           ysx=yl(4,ndz+2-kz)+(yl(2,kz)-yl(4,ndz+2-kz))
     &                       *dble(kx-1)/dble(ndx)
           ysz=yl(1,kx)+(yl(3,ndx+2-kx)-yl(1,kx))
     &                 *dble(kz-1)/dble(ndz)
           ys1(kx,kz)=(ysx+ysz)/2d0
c
c          面6:X-Z平面
c
           xs6(kx,kz)=xl(12,ndz+2-kz)+(xl(10,kz)-xl(12,ndz+2-kz))
     &                               *dble(kx-1)/dble(ndx)
           zs6(kx,kz)=zl(9,kx)+(zl(11,ndx+2-kx)-zl(9,kx))
     &                        *dble(kz-1)/dble(ndz)
           ysx=yl(12,ndz+2-kz)+(yl(10,kz)-yl(12,ndz+2-kz))
     &                        *dble(kx-1)/dble(ndx)
           ysz=yl(9,kx)+(yl(11,ndx+2-kx)-yl(9,kx))
     &                 *dble(kz-1)/dble(ndz)
           ys6(kx,kz)=(ysx+ysz)/2d0
1500     continue
c
         do 1600 ky=1,ndy+1
         do 1600 kz=1,ndz+1
c
c          面2:Y-Z平面
c
           ys2(ky,kz)=yl(2,kz)+(yl(10,kz)-yl(2,kz))
     &                        *dble(ky-1)/dble(ndy)
           zs2(ky,kz)=zl(5,ky)+(zl(6,ky)-zl(5,ky))
     &                        *dble(kz-1)/dble(ndz)
           xsy=xl(2,kz)+(xl(10,kz)-xl(2,kz))
     &                 *dble(ky-1)/dble(ndy)
           xsz=xl(5,ky)+(xl(6,ky)-xl(5,ky))
     &                 *dble(kz-1)/dble(ndz)
           xs2(ky,kz)=(xsy+xsz)/2d0
c
c          面4:Y-Z平面
c
           ys4(ky,kz)=yl(4,ndz+2-kz)+(yl(12,ndz+2-kz)-yl(4,ndz+2-kz))
     &                              *dble(ky-1)/dble(ndy)
           zs4(ky,kz)=zl(8,ky)+(zl(7,ky)-zl(8,ky))
     &                        *dble(kz-1)/dble(ndz)
           xsy=xl(4,ndz+2-kz)+(xl(12,ndz+2-kz)-xl(4,ndz+2-kz))
     &                       *dble(ky-1)/dble(ndy)
           xsz=xl(8,ky)+(xl(7,ky)-xl(8,ky))
     &                 *dble(kz-1)/dble(ndz)
           xs4(ky,kz)=(xsy+xsz)/2d0
1600     continue
c
         do 1700 kx=1,ndx+1
         do 1700 ky=1,ndy+1
c
c          面3:X-Y平面
c
           xs3(kx,ky)=xl(7,ky)+(xl(6,ky)-xl(7,ky))
     &                        *dble(kx-1)/dble(ndx)
           ys3(kx,ky)=yl(3,ndx+2-kx)+(yl(11,ndx+2-kx)-yl(3,ndx+2-kx))
     &                              *dble(ky-1)/dble(ndy)
           zsx=zl(7,ky)+(zl(6,ky)-zl(7,ky))
     &                 *dble(kx-1)/dble(ndx)
           zsy=zl(3,ndx+2-kx)+(zl(11,ndx+2-kx)-zl(3,ndx+2-kx))
     &                       *dble(ky-1)/dble(ndy)
           zs3(kx,ky)=(zsx+zsy)/2d0
c
c          面5:X-Y平面
c
           xs5(kx,ky)=xl(8,ky)+(xl(5,ky)-xl(8,ky))
     &                        *dble(kx-1)/dble(ndx)
           ys5(kx,ky)=yl(1,kx)+(yl(9,kx)-yl(1,kx))
     &                        *dble(ky-1)/dble(ndy)
           zsx=zl(8,ky)+(zl(5,ky)-zl(8,ky))
     &                 *dble(kx-1)/dble(ndx)
           zsy=zl(1,kx)+(zl(9,kx)-zl(1,kx))
     &                 *dble(ky-1)/dble(ndy)
           zs5(kx,ky)=(zsx+zsy)/2d0
1700     continue
c         write(*,*)"ended deciding the nodes' positions on each surface
c     & of hexahedral elements"
c
c        6面体内部の各節点座標を求める
c
         do 1800 kx=1,ndx+1
         do 1800 ky=1,ndy+1
         do 1800 kz=1,ndz+1
c           write(*,*)'krg,kx,ky,kz',krg,kx,ky,kz
           nn=nn+1      !nn:節点番号は重なる
           iwr(nn)=krg           
           iwx(nn)=kx
           iwy(nn)=ky
           iwz(nn)=kz
c
c          X座標:面4-2
c
           xr(nn)=xs4(ky,kz)+(xs2(ky,kz)-xs4(ky,kz))
     &                                     *dble(kx-1)/dble(ndx)
c           write(*,*)'xr',sngl(xs2(ky,kz)),sngl(xs4(ky,kz)),sngl(xr(nn))
c
c          Y座標:面1-6
c
           yr(nn)=ys1(kx,kz)+(ys6(kx,kz)-ys1(kx,kz))
     &                                     *dble(ky-1)/dble(ndy)
c           write(*,*)'yr',sngl(ys1(kx,kz)),sngl(ys6(kx,kz)),sngl(yr(nn))     
c
c          Z座標:面5-3
c
           zr(nn)=zs5(kx,ky)+(zs3(kx,ky)-zs5(kx,ky))
     &                                     *dble(kz-1)/dble(ndz)
c           write(*,*)'zr',sngl(zs5(kx,ky)),sngl(zs3(kx,ky)),sngl(zr(nn))
c           write(*,*)'nn',nn
c           write(*,*)'xr,yr,zy',sngl(xr(nn)),sngl(yr(nn)),sngl(zr(nn))
c           pause
1800    continue
1100  continue
c      write(*,*)"ended deciding the nodes' positions in each hexahedral
c     & element"
c      pause
c
c      各節点に番号を割り振る
c      nrg,mk(nrg,8),xr(nn),yr(nn),zr(nn),iwr(nn),iwx(nn),iwy(nn),iwz(nn) →
c      np,xx(np),yy(np),zz(np),lwr(np),lwx(np),lwy(np),lwz(np)
c       np=0
c       ii=0
c       do 2300 krg=1,nrg            !突起のある領域は節点番号が重なる
c         ndx=mdx(krg)
c         if(mk(krg,2).eq.0)ndx=ndx+1        !奥のY-Z平面を考慮
c         do 2400 kx=1,ndx
c           ndy=mdy(krg)
c           if(mk(krg,6).eq.0)ndy=ndy+1      !奥のZ-X平面を考慮
c           do 2500 ky=1,ndy
c             ndz=mdz(krg)
c             if(mk(krg,3).eq.0)ndz=ndz+1    !奥のX-Y平面を考慮
c             do 2600 kz=1,ndz
c               ii=ii+1      !隣り合う領域が共有する面上の節点を考慮
c               np=np+1      !節点番号は重ならない
c               lwr(ii)=krg
c               lwx(ii)=kx
c               lwy(ii)=ky
c               lwz(ii)=kz
c               lwp(ii)=np
c               write(*,*)'krg,kx,ky,kz,ii',krg,kx,ky,kz,ii
cc
c               do i=1,nn
c                 if(krg.eq.iwr(i) .and. kx.eq.iwx(i) .and. 
c     &               ky.eq.iwy(i) .and. kz.eq.iwz(i))then
c                   xx(np)=xr(i)
c                   yy(np)=yr(i)
c                   zz(np)=zr(i)
c                 endif
c               enddo
cc
c2600         continue
c2500       continue
c2400     continue
c2300   continue
c---------------------------------------------------------------------------------
c---------------------------------------------------------------------------------
c------6面体の頂点に節点番号を割り振る----------------------------------------------------
       np=npr
       ii=0
       do kp=1,np
         do krg=1,nrg   !領域番号
         do kk=1,8      !頂点番号
           if(mp(krg,kk).eq.kp)then    !領域番号krg,頂点番号kkの6面体
             ii=ii+1
             lwr(ii)=krg
             lwp(ii)=kp
             if(kk.eq.1)then       !頂点1
               lwx(ii)=1
               lwy(ii)=1
               lwz(ii)=1
             elseif(kk.eq.2)then   !頂点2
               lwx(ii)=mdx(krg)+1
               lwy(ii)=1
               lwz(ii)=1
             elseif(kk.eq.3)then   !頂点3
               lwx(ii)=mdx(krg)+1
               lwy(ii)=1
               lwz(ii)=mdz(krg)+1
             elseif(kk.eq.4)then   !頂点4
               lwx(ii)=1
               lwy(ii)=1
               lwz(ii)=mdz(krg)+1
             elseif(kk.eq.5)then   !頂点5
               lwx(ii)=1
               lwy(ii)=mdy(krg)+1
               lwz(ii)=1
             elseif(kk.eq.6)then   !頂点6
               lwx(ii)=mdx(krg)+1
               lwy(ii)=mdy(krg)+1
               lwz(ii)=1
             elseif(kk.eq.7)then   !頂点7
               lwx(ii)=mdx(krg)+1
               lwy(ii)=mdy(krg)+1
               lwz(ii)=mdz(krg)+1
             elseif(kk.eq.8)then   !頂点8
               lwx(ii)=1
               lwy(ii)=mdy(krg)+1
               lwz(ii)=mdz(krg)+1
             endif
           endif
         enddo
         enddo
       enddo  
c      
c       open(1,file='vertex.res')
c       rewind(1)
c       do kr=1,nrg
c         nnn=0
c         do i=1,ii
c           if(lwr(i).eq.kr)then
c             write(1,*)i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
c             nnn=nnn+1
c           endif
c         enddo
c         write(1,*)'kr,nnn',kr,nnn
c       enddo
c       close(1)
c
       write(*,*)'only vertexes np,ii',np,ii  !ii=64
c      mm(krg,12):節点番号の割り当てを示す．0:割り振られている,1:割り振られていない
c      lln:辺の数
c      lli(lln):辺の片端の節点番号
c      llj(lln):辺の片端の節点番号
c      llk(lln):辺に対する節点番号の割り当てを示す.0:割り当てられている,1:割り当てられていない
c
c------初期化------------------------------------------------------------------------------
c
       lln=0                    !lln=54:辺の数
       do 1810 krg=1,nrg
       do 1820 kl=1,12
         kds=icals(mp,krg,kl)
         kde=icale(mp,krg,kl)
         do kkn=1,lln
           if(lli(kkn).eq.kds .and. llj(kkn).eq.kde)goto 1820  !既に考慮されている
           if(lli(kkn).eq.kde .and. llj(kkn).eq.kds)goto 1820  !既に考慮されている
         enddo
         lln=lln+1      !考慮されていない辺
         lli(lln)=kds
         llj(lln)=kde
         llk(lln)=0     !節点番号が割り当てられていない
         mm(krg,kl)=0   !節点番号が割り当てられていない
1820   continue
1810   continue
       write(*,*)'the number of sides lln',lln
c
c------各6面体の辺上に節点番号を割り当てる---------------------------------------------
c
c      ii:各6面体を構成する節点の数(重複)
       npn=npr  !節点の数
       do 1840 krg=1,nrg
       do 1850 kl=1,12
         ndx=mdx(krg)
         ndy=mdy(krg)
         ndz=mdz(krg)
         kds=icals(mp,krg,kl)       !考慮すべき辺の始点の節点番号
         kde=icale(mp,krg,kl)       !考慮すべき辺の終点の節点番号
         do kkn=1,lln
           if(lli(kkn).eq.kds .and. llj(kkn).eq.kde .and. llk(kkn).eq.1)  !既に節点番号が割り当てられた辺
     &                                                      goto 1850     !考慮しない
           if(lli(kkn).eq.kde .and. llj(kkn).eq.kds .and. llk(kkn).eq.1)  !既に節点番号が割り当てられた辺
     &                                                      goto 1850     !考慮しない
         enddo
c         if(mm(lrg,ll).eq.1)goto 1850      !考慮している領域の辺上に節点が既に割り当てられている場合は，考慮しない
c
c         pause
c         write(*,*)'pass1,krg,kl,kds,kde',krg,kl,kds,kde
c
c--------割り当てる節点番号を格納した配列lm(nn)を求める------------------------------------
         if(kl.eq.1 .or. kl.eq.3 .or. kl.eq.9 .or. kl.eq.11)then        !X軸に平行な辺
           do kk=2,ndx
             npn=npn+1
             lm(kk)=npn
           enddo
         elseif(kl.eq.5 .or. kl.eq.6 .or. kl.eq.7 .or. kl.eq.8)then    !Y軸に平行な辺
           do kk=2,ndy
             npn=npn+1
             lm(kk)=npn
           enddo
         elseif(kl.eq.2 .or. kl.eq.10 .or. kl.eq.12 .or. kl.eq.4)then  !Z軸に平行な辺
           do kk=2,ndz
             npn=npn+1
             lm(kk)=npn
           enddo
         endif
c--------------------------------------------------------------------------------------------
         do 1860 lrg=1,nrg
         do 1870 ll=1,12
           lds=icals(mp,lrg,ll)
           lde=icale(mp,lrg,ll)
           if((kds.eq.lds .and. kde.eq.lde) .or.   !考慮している節点番号で構成されてている辺
     &        (kds.eq.lde .and. kde.eq.lds))then
             do kkn=1,lln
               if(lli(kkn).eq.kds .and. llj(kkn).eq.kde)llk(kkn)=1
               if(lli(kkn).eq.kde .and. llj(kkn).eq.kds)llk(kkn)=1
             enddo
c
c             write(*,*)'pass2,lrg,ll,lds,lde,ii',lrg,ll,lds,lde,ii
             if(ll.eq.1)then       !辺1
               do kk=2,ndx
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=kk
                 lwy(ii)=1
                 lwz(ii)=1
                 lwp(ii)=lm(kk)     !X軸に正の向き
                 mm(lrg,ll)=1       !考慮した6面体の辺
               enddo
             elseif(ll.eq.2)then   !辺2
               do kk=2,ndz
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=mdx(lrg)+1
                 lwy(ii)=1
                 lwz(ii)=kk
                 lwp(ii)=lm(kk)     !Z軸に正の向き
                 mm(lrg,ll)=1       !考慮した6面体の辺
               enddo
             elseif(ll.eq.3)then   !辺3
               do kk=2,ndx
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=ndx+2-kk
                 lwy(ii)=1
                 lwz(ii)=mdz(lrg)+1
                 lwp(ii)=lm(ndx+2-kk) !X軸に負の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.4)then   !辺4
               do kk=2,ndz
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=1
                 lwy(ii)=1
                 lwz(ii)=ndz+2-kk
                 lwp(ii)=lm(ndz+2-kk) !Z軸に負の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.5)then   !辺5
               do kk=2,ndy
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=mdx(lrg)+1
                 lwy(ii)=kk
                 lwz(ii)=1
                 lwp(ii)=lm(kk)     !Y軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.6)then   !辺6
               do kk=2,ndy
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=mdx(lrg)+1
                 lwy(ii)=kk
                 lwz(ii)=mdz(lrg)+1
                 lwp(ii)=lm(kk)     !Y軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.7)then   !辺7
               do kk=2,ndy
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=1
                 lwy(ii)=kk
                 lwz(ii)=mdz(lrg)+1
                 lwp(ii)=lm(kk)     !Y軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.8)then   !辺8
               do kk=2,ndy
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=1
                 lwy(ii)=kk
                 lwz(ii)=1
                 lwp(ii)=lm(kk)     !Y軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.9)then   !辺9
               do kk=2,ndx
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=kk
                 lwy(ii)=mdy(lrg)+1
                 lwz(ii)=1
                 lwp(ii)=lm(kk)     !X軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.10)then  !辺10
               do kk=2,ndz
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=mdx(lrg)+1
                 lwy(ii)=mdy(lrg)+1
                 lwz(ii)=kk
                 lwp(ii)=lm(kk)     !Z軸に正の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.11)then  !辺11
               do kk=2,ndx
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=ndx+2-kk
                 lwy(ii)=mdy(lrg)+1
                 lwz(ii)=mdz(lrg)+1
                 lwp(ii)=lm(ndx+2-kk) !X軸に負の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             elseif(ll.eq.12)then  !辺12
               do kk=2,ndz
                 ii=ii+1
                 lwr(ii)=lrg
                 lwx(ii)=1
                 lwy(ii)=mdy(lrg)+1
                 lwz(ii)=ndz+2-kk
                 lwp(ii)=lm(ndz+2-kk) !Z軸に負の向き
                 mm(lrg,ll)=1       !考慮した領域番号lrg,辺番号llの6面体
               enddo
             endif
           endif
1870     continue
1860     continue
1850   continue
1840   continue
c       write(*,*)'put away the number of nodes on all vertex'
       write(*,*)'vertexes and sides npn,ii',npn,ii
c
c       open(1,file='sides.res')
c       rewind(1)
c       do kr=1,nrg
c         nnn=0
c         do i=1,ii
c           if(lwr(i).eq.kr)then
c             write(1,*)lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
c             nnn=nnn+1
c           endif
c         enddo
c         write(1,*)'kr,nnn',kr,nnn
c       enddo
c       close(1)
c      
c------6面体の面上に節点番号を割り振る--------------------------------------------------
c      ii:各6面体を構成する節点の数(重複)
c      npn:現在割り当てられた節点番号の数
c      ms(nrg,6):6面体の面上に対する節点番号の割り当てを示す.0:割り当てられていない,1:割り当てられている
c      初期化
       do 1880 krg=1,nrg
       do 1880 ks=1,6
c       write(*,*)'mk(',krg,',',ks,')',mk(krg,ks)
       ms(krg,ks)=0
1880   continue

c
       do 1890 krg=1,nrg
       do 1900 ks=1,6
         if(ms(krg,ks).eq.1)goto 1900  !既に節点番号が割り当てられている面は考慮しない
         ndx=mdx(krg)
         ndy=mdy(krg)
         ndz=mdz(krg)
         lrg=mk(krg,ks)     !lrg:面を共有する領域の番号
         ls=isf(ks)         !ls:面を共有する領域の面番号
         ifs=1
         if(lrg.eq.0)ifs=0  !面を共有する領域が存在しない
c         write(*,*)'pass krg,ks,lrg,ls,npn',krg,ks,lrg,ls,npn
c         write(*,*)'ndx,ndy,ndz,ifs',ndx,ndy,ndz,ifs
c         pause
c
         if(ks.eq.1)then   !ks:面1,ls:面6←ZX平面
           do kz=2,ndz
           do kx=2,ndx
             npn=npn+1           
             ii=ii+1    !面1
             lwr(ii)=krg
             lwx(ii)=kx
             lwy(ii)=1
             lwz(ii)=kz
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号krg,面番号ksの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1    !面6
               lwr(ii)=lrg
               lwx(ii)=kx
               lwy(ii)=mdy(lrg)+1
               lwz(ii)=kz
               lwp(ii)=npn
               ms(lrg,ls)=1   !考慮した領域番号lrg,面番号lsの6面体
             endif
           enddo
           enddo
         elseif(ks.eq.2)then   !ks:面2,ls:面4←YZ平面
           do ky=2,ndy
           do kz=2,ndz
             npn=npn+1
             ii=ii+1    !面2
             lwr(ii)=krg
             lwx(ii)=mdx(krg)+1
             lwy(ii)=ky
             lwz(ii)=kz
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号krg,面番号ksの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1    !面4
               lwr(ii)=lrg
               lwx(ii)=1
               lwy(ii)=ky
               lwz(ii)=kz
               lwp(ii)=npn
               ms(lrg,ls)=1   !考慮した領域番号lrg,面番号lsの6面体
             endif
           enddo
           enddo
         elseif(ks.eq.3)then   !ks:面3,ls:面5←XY平面
           do kx=2,ndx
           do ky=2,ndy
             npn=npn+1
             ii=ii+1    !面3
             lwr(ii)=krg
             lwx(ii)=kx
             lwy(ii)=ky
             lwz(ii)=mdz(krg)+1
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号krg,面番号ksの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1    !面5
               lwr(ii)=lrg
               lwx(ii)=kx
               lwy(ii)=ky
               lwz(ii)=1
               lwp(ii)=npn
               ms(lrg,ls)=1   !考慮した領域番号lrg,面番号lsの6面体
             endif
           enddo
           enddo
         elseif(ks.eq.4)then   !ks:面4,ls:面2←YZ平面
           do ky=2,ndy
           do kz=2,ndz
             npn=npn+1
             ii=ii+1    !面4
             lwr(ii)=krg
             lwx(ii)=1
             lwy(ii)=ky
             lwz(ii)=kz
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号lrg,面番号lsの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1    !面2
               lwr(ii)=lrg
               lwx(ii)=mdx(lrg)+1
               lwy(ii)=ky
               lwz(ii)=kz
               lwp(ii)=npn
               ms(lrg,ls)=1
             endif
           enddo
           enddo
         elseif(ks.eq.5)then   !ks:面5,ls:面3←XY平面
           do kx=2,ndx
           do ky=2,ndy
             npn=npn+1
             ii=ii+1    !面5
             lwr(ii)=krg
             lwx(ii)=kx
             lwy(ii)=ky
             lwz(ii)=1
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号krg,面番号ksの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1    !面3
               lwr(ii)=lrg
               lwx(ii)=kx
               lwy(ii)=ky
               lwz(ii)=mdz(lrg)+1
               lwp(ii)=npn
               ms(lrg,ls)=1   !考慮した領域番号lrg,面番号lsの6面体
             endif
           enddo
           enddo
         elseif(ks.eq.6)then   !ks:面6,ls:面1←ZX平面
           do kz=2,ndz
           do kx=2,ndx
             npn=npn+1
             ii=ii+1    !面6
             lwr(ii)=krg
             lwx(ii)=kx
             lwy(ii)=mdy(krg)+1
             lwz(ii)=kz
             lwp(ii)=npn
             ms(krg,ks)=1   !考慮した領域番号krg,面番号ksの6面体
c
             if(ifs.eq.1)then !面を共有する領域が存在する
               ii=ii+1
               lwr(ii)=lrg
               lwx(ii)=kx
               lwy(ii)=1
               lwz(ii)=kz
               lwp(ii)=npn
               ms(lrg,ls)=1   !考慮した領域番号lrg,面番号lsの6面体
             endif
           enddo
           enddo
         endif
1900   continue
1890   continue
c       write(*,*)'put away the number of nodes on all surface'
       write(*,*)'vertexes, sides and surfaces npn,ii',npn,ii
c------各6面体内部に節点番号を割り当てる--------------------------------------------------
c      ii:各6面体を構成する節点の数(重複)
c      npn:現在割り当てられた節点番号の数
       do 1910 krg=1,nrg
       do 1910 kx=2,mdx(krg)
       do 1910 ky=2,mdy(krg)
       do 1910 kz=2,mdz(krg)
         npn=npn+1
         ii=ii+1
         lwr(ii)=krg    !領域番号krgの6面体
         lwx(ii)=kx
         lwy(ii)=ky
         lwz(ii)=kz
         lwp(ii)=npn
c         write(*,*)ii,lwr(ii),lwx(ii),lwy(ii),lwz(ii),lwp(ii)
1910   continue
       np=npn
       write(*,*)'all nodes npn,ii',npn,ii
c       write(*,*)"ended allotting the number of nodes to each vertexes"
c       pause
c---------------------------------------------------------------------------------
c---------------------------------------------------------------------------------       
c
c      各6面体要素の節点番号を格納した配列ne(nele,8)を求める
c
       nele=0
       do 2700 krg=1,nrg
         do 2800 kx=1,mdx(krg)
           do 2900 ky=1,mdy(krg)
             do 3000 kz=1,mdz(krg)
               nele=nele+1
               do i=1,ii
                 if(krg  .eq.lwr(i) .and. kx.eq.lwx(i) .and.
     &               ky  .eq.lwy(i) .and. kz.eq.lwz(i))then
                   ne(nele,1)=lwp(i)
c                   write(*,*)'1,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx+1.eq.lwx(i) .and.
     &               ky  .eq.lwy(i) .and. kz  .eq.lwz(i))then
                   ne(nele,2)=lwp(i)
c                   write(*,*)'2,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx+1.eq.lwx(i) .and.
     &               ky  .eq.lwy(i) .and. kz+1.eq.lwz(i))then
                   ne(nele,3)=lwp(i)
c                   write(*,*)'3,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx  .eq.lwx(i) .and.
     &               ky  .eq.lwy(i) .and. kz+1.eq.lwz(i))then
                   ne(nele,4)=lwp(i)
c                   write(*,*)'4,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx  .eq.lwx(i) .and.
     &               ky+1.eq.lwy(i) .and. kz  .eq.lwz(i))then
                   ne(nele,5)=lwp(i)
c                   write(*,*)'5,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx+1.eq.lwx(i) .and.
     &               ky+1.eq.lwy(i) .and. kz  .eq.lwz(i))then
                   ne(nele,6)=lwp(i)
c                   write(*,*)'6,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx+1.eq.lwx(i) .and.
     &               ky+1.eq.lwy(i) .and. kz+1.eq.lwz(i))then
                   ne(nele,7)=lwp(i)
c                   write(*,*)'7,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 elseif(krg  .eq.lwr(i) .and. kx  .eq.lwx(i) .and.
     &               ky+1.eq.lwy(i) .and. kz+1.eq.lwz(i))then
                   ne(nele,8)=lwp(i)
c                   write(*,*)'8,i,lwr,lwx,lwy,lwx,lwp',
c     &             i,lwr(i),lwx(i),lwy(i),lwz(i),lwp(i)
                 endif
               enddo
               md(nele)=mr(krg)
c               write(*,*)'---------------------------------------------'
c               write(*,*)'krg,kx,ky,kz',krg,kx,ky,kz
c               write(*,*)'ne(',nele,',',1,')',ne(nele,1)
c               write(*,*)'ne(',nele,',',2,')',ne(nele,2)
c               write(*,*)'ne(',nele,',',3,')',ne(nele,3)
c               write(*,*)'ne(',nele,',',4,')',ne(nele,4)
c               write(*,*)'ne(',nele,',',5,')',ne(nele,5)
c               write(*,*)'ne(',nele,',',6,')',ne(nele,6)
c               write(*,*)'ne(',nele,',',7,')',ne(nele,7)
c               write(*,*)'ne(',nele,',',8,')',ne(nele,8)
c--------------------------------------------------------------------------
                mg(nele)=krg		!mg(nele):6面体要素neleが存在する領域の番号krg
c--------------------------------------------------------------------------
3000         continue
2900       continue
2800     continue
2700   continue
       write(*,*)'nele',nele
c       write(*,*)"ended putting away the number of nodes in the
c     & arrangements of the hexahedral element,ne(ind,8)"
c
c      4面体に分割する:6面体1つ→4面体6つ
c
       nelw=0
       do 3100 kele=1,nele     !6面体
         nelw=nelw+1
         new(nelw,1)=ne(kele,6)
         new(nelw,2)=ne(kele,5)
         new(nelw,3)=ne(kele,8)
         new(nelw,4)=ne(kele,1)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
         nelw=nelw+1
         new(nelw,1)=ne(kele,4)
         new(nelw,2)=ne(kele,1)
         new(nelw,3)=ne(kele,2)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
         nelw=nelw+1
         new(nelw,1)=ne(kele,8)
         new(nelw,2)=ne(kele,1)
         new(nelw,3)=ne(kele,4)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
         nelw=nelw+1
         new(nelw,1)=ne(kele,7)
         new(nelw,2)=ne(kele,6)
         new(nelw,3)=ne(kele,8)
         new(nelw,4)=ne(kele,3)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
         nelw=nelw+1
         new(nelw,1)=ne(kele,3)
         new(nelw,2)=ne(kele,4)
         new(nelw,3)=ne(kele,2)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
         nelw=nelw+1
         new(nelw,1)=ne(kele,3)
         new(nelw,2)=ne(kele,8)
         new(nelw,3)=ne(kele,4)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
	   mgr(nelw)=mg(kele)	!mgr(nelw):4面体要素nelwが存在する領域の番号mg(nele)
c-------------------------------------------------------------------------------
3100   continue
       write(*,*)'nelw',nelw
c       write(*,*)"ended putting away the number of nodes in the
c     & arrangements of the tetrahedral element,new(ine,4)"
c
c      x,y,z座標を格納する配列xx(np),yy(np),zz(np)を求める
c
       do kp=1,np
         do kk=1,ii
           if(lwp(kk).eq.kp)then
             do jj=1,nn
               if(lwr(kk).eq.iwr(jj) .and. lwx(kk).eq.iwx(jj) .and.
     &            lwy(kk).eq.iwy(jj) .and. lwz(kk).eq.iwz(jj))then
                 xx(kp)=xr(jj)
                 yy(kp)=yr(jj)
                 zz(kp)=zr(jj)
c                 write(*,*)lwr(kk),lwx(kk),lwy(kk),lwz(kk)
c                 write(*,*)'kp,jj',kp,jj,sngl(xr(jj)),
c     &                                    sngl(yr(jj)),sngl(zr(jj))
c                 pause
               endif
             enddo
           endif
         enddo
       enddo
c
c       open(1,file='xx.res')
c       rewind(1)
c       write(1,*)'nelw,npw',nelw,np
c       do kele=1,nelw
c       write(1,*)kele,new(kele,1),new(kele,2),new(kele,3),new(kele,4),
c     &          xx(new(kele,1)),yy(new(kele,1)),zz(new(kele,1)),
c     &          xx(new(kele,2)),yy(new(kele,2)),zz(new(kele,2)),
c     &          xx(new(kele,3)),yy(new(kele,3)),zz(new(kele,3)),
c     &          xx(new(kele,4)),yy(new(kele,4)),zz(new(kele,4)),
c     &          idx(kele)
c       enddo
c       close(1)
c
c        体積Vを算出する
         do im=1,nelw
         do i=1,4
           x(i)=xx(new(im,i))
           y(i)=yy(new(im,i))
           z(i)=zz(new(im,i))
         enddo
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
           write(*,*)'im,vl',im,vl
           pause
         endif
         enddo
c
c
c      バンド幅を最適化
c
c------初期化-------------------------------------------------------------
       npw=np
       nbw=0
       do i=1,ind
         jww(i)=0
         jmm(i)=0
       enddo
c--------------------------------------------------------------------------
c       call renum3(npw,nelw,new,jww,jmm,nbw)
c***************************************************
        do i=1,npw
          jww(i)=i
          jmm(i)=i
        enddo
c
c        バンド幅を求める
c         write(*,*)'calc band'
         nbm=-99999
         do kele=1,nelw
           do i=1,4
             do j=1,4
               kbm=iabs(jww(new(kele,i))-jww(new(kele,j)))
               if(kbm.gt.nbm)nbm=kbm     !nbm:現在考慮しているtreeのバンド幅
             enddo
           enddo
         enddo
         nbw=nbm+1
         write(*,*)'band',nbw
c***************************************************
       write(*,*)"ended suiting the number of nodes"       
c
c      *.1および*,2ファイルの作成
c
       call addchr(rfile,'.1',gfile1)
       open(1,file=gfile1)
       rewind(1)
       write(1,*)npw
       do 3300 i=1,npw
3300   write(1,*)i,jww(i),xx(i),yy(i),zz(i)
       close(1)
c
       call addchr(rfile,'.2',gfile2)
       open(2,file=gfile2)
       rewind(2)
       write(2,*)nelw,npw,nbw
       do 3400 i=1,nelw
3400   write(2,*)i,(jww(new(i,j)),j=1,4)
     &             ,(xx(new(i,j)),yy(new(i,j)),zz(new(i,j)),j=1,4)
     &             ,idx(i)
       close(2)
c       write(*,*)'ended creating the file'
       write(*,*)'created ',gfile1
       write(*,*)'created ',gfile2
c       pause
c***************************************************************************
c
c      bdpv.infから境界条件のファイルbd.infを作成する
c
c***************************************************************************
       write(*,*)'**** about bd files ****'
       open(1,file='bdpv.inf')
       rewind(1)
       read(1,*)nn     !nn:境界条件の数
	 do 3500 i=1,nn
3500   read(1,*)mbr(i),mbf(i),mbn(i),faib(i)
       close(1)
       write(*,*)'readed bdpv.inf'
c----------------------------------------------------------
c
c      n(17),kbr(17,n(17)),kbf(17,n(17))の作成
c      n(17):境界条件の数
c      kbr(17,n(17)):領域番号
c      kbf(17,n(17)):面番号
c----------------------------------------------------------
c
c      初期化
c
       do 3600 i=1,17
3600 	 n(i)=0
c
c      全領域から検索
c
       do 3700 i=1,nn
	   if(mbn(i).eq.1)then		!ub(nb(1))
	     n(1)=n(1)+1
	     kbr(1,n(1))=mbr(i)
	     kbf(1,n(1))=mbf(i)
	     ubt(n(1))=faib(i)
	   elseif(mbn(i).eq.2)then	!vb(nb(2))
	     n(2)=n(2)+1
	     kbr(2,n(2))=mbr(i)
	     kbf(2,n(2))=mbf(i)
	     vbt(n(2))=faib(i)
	   elseif(mbn(i).eq.3)then	!wb(nb(3))
	     n(3)=n(3)+1
	     kbr(3,n(3))=mbr(i)
	     kbf(3,n(3))=mbf(i)
	     wbt(n(3))=faib(i)
	   elseif(mbn(i).eq.4)then	!pb(nb(4))
	     n(4)=n(4)+1
	     kbr(4,n(4))=mbr(i)
	     kbf(4,n(4))=mbf(i)
	     pbt(n(4))=faib(i)
	   elseif(mbn(i).eq.5)then	!自由表面
	     n(5)=n(5)+1
	     kbr(5,n(5))=mbr(i)
	     kbf(5,n(5))=mbf(i)
	   elseif(mbn(i).eq.6)then	!cb(nb(6))
	     n(6)=n(6)+1
	     kbr(6,n(6))=mbr(i)
	     kbf(6,n(6))=mbf(i)
	     cbt(n(6))=faib(i)
	   elseif(mbn(i).eq.7)then	!tb(nb(7))
	     n(7)=n(7)+1
	     kbr(7,n(7))=mbr(i)
	     kbf(7,n(7))=mbf(i)
	     tbt(n(7))=faib(i)
	   elseif(mbn(i).eq.8)then	!ug(n(8))
	     n(8)=n(8)+1
	     kbr(8,n(8))=mbr(i)
	     kbf(8,n(8))=mbf(i)
	     ugt(n(8))=faib(i)
	   elseif(mbn(i).eq.9)then	!vg(n(9))
	     n(9)=n(9)+1
	     kbr(9,n(9))=mbr(i)
	     kbf(9,n(9))=mbf(i)
	     vgt(n(9))=faib(i)
	   elseif(mbn(i).eq.10)then	!wg(n(10))
	     n(10)=n(10)+1
	     kbr(10,n(10))=mbr(i)
	     kbf(10,n(10))=mbf(i)
	     wgt(n(10))=faib(i)
	   elseif(mbn(i).eq.11)then	!cg(n(11))
	     n(11)=n(11)+1
	     kbr(11,n(11))=mbr(i)
	     kbf(11,n(11))=mbf(i)
	     cgt(n(11))=faib(i)
	   elseif(mbn(i).eq.12)then	!tg(n(12))
	     n(12)=n(12)+1
	     kbr(12,n(12))=mbr(i)
	     kbf(12,n(12))=mbf(i)
	     tgt(n(12))=faib(i)
	   elseif(mbn(i).eq.13)then	!ch(n(13))
	     n(13)=n(13)+1
	     kbr(13,n(13))=mbr(i)
	     kbf(13,n(13))=mbf(i)
	     cht(nb(13))=faib(i)
	   elseif(mbn(i).eq.14)then	!th(n(14))
	     n(14)=n(14)+1
	     kbr(14,n(14))=mbr(i)
	     kbf(14,n(14))=mbf(i)
	     tht(n(14))=faib(i)
	   elseif(mbn(i).eq.15)then	!scv(n(15))
	     n(15)=n(15)+1
	     kbr(15,n(15))=mbr(i)
	     kbf(15,n(15))=mbf(i)
	     scvt(n(15))=faib(i)
	   elseif(mbn(i).eq.16)then	!scm(n(16))
	     n(16)=n(16)+1
	     kbr(16,n(16))=mbr(i)
	     kbf(16,n(16))=mbf(i)
	     scmt(n(16))=faib(i)
	   elseif(mbn(i).eq.17)then	!sct(n(17))
	     n(17)=n(17)+1
	     kbr(17,n(17))=mbr(i)
	     kbf(17,n(17))=mbf(i)
	     sctt(n(17))=faib(i)
	   endif
3700   continue
c       write(*,*)'made n(17),kbr(17,n(17)),kbf(17,n(17))'
c
c	 do i=1,17
c	 write(*,*)'--------------------------------------'
c	 write(*,*)'n(',i,')',n(i)
c	 do j=1,n(i)
c	 write(*,*)j,kbr(i,j),kbf(i,j)
c	 enddo
c	 pause
c	 enddo
c
c----------------------------------------------------------
c
c      kb(11),kelb(11,kb(11)),keb(11,kb(11))の作成
c
c----------------------------------------------------------
c
c	 初期化
c
       do 3800 krg=1,nrg
	 do 3800 ks=1,6
	   ub(krg,ks)=-999d0
	   vb(krg,ks)=-999d0
	   wb(krg,ks)=-999d0
	   pb(krg,ks)=-999d0
	   cb(krg,ks)=-999d0
	   tb(krg,ks)=-999d0
	   ug(krg,ks)=-999d0
	   vg(krg,ks)=-999d0
	   wg(krg,ks)=-999d0
	   cg(krg,ks)=-999d0
	   tg(krg,ks)=-999d0
	   ch(krg,ks)=-999d0
	   th(krg,ks)=-999d0
	   scv(krg,ks)=-999d0
	   scm(krg,ks)=-999d0
	   sct(krg,ks)=-999d0
3800   continue
c
       do i=1,11
	   kb(i)=0
	 enddo
c
c	 全領域から検索
c
       do 3900 krg=1,nrg
	 do 3900 ks=1,6
c
c        kb(1),kelb(1,kb(1)),keb(1,kb(1))について
c        ub(kb(1)),vb(kb(1)),wb(kb(1)),pb(kb(1))
	   ifg1=0
	   do i=1,n(1)	!ubt(n(1))
	     if(kbr(1,i).eq.krg .and. kbf(1,i).eq.ks)then
	       if(ifg1.eq.0)then
	         kb(1)=kb(1)+1
	         kelb(1,kb(1))=krg
	         keb(1,kb(1))=ks
		     ifg1=1
	       endif
c	       ub(kb(1))=ubt(i)
		   ub(kelb(1,kb(1)),keb(1,kb(1)))=ubt(i)	!ub(krg,ks)
	     endif
	   enddo
c
	   do i=1,n(2)	!vbt(n(2))
	     if(kbr(2,i).eq.krg .and. kbf(2,i).eq.ks)then
	       if(ifg1.eq.0)then
	         kb(1)=kb(1)+1
	         kelb(1,kb(1))=krg
	         keb(1,kb(1))=ks
	         ifg1=1
	       endif
c	       vb(kb(1))=vbt(i)
		   vb(kelb(1,kb(1)),keb(1,kb(1)))=vbt(i)	!vb(krg,ks)
	     endif
	   enddo
c
	   do i=1,n(3)	!wbt(n(3))
	     if(kbr(3,i).eq.krg .and. kbf(3,i).eq.ks)then
	       if(ifg1.eq.0)then
	         kb(1)=kb(1)+1
	         kelb(1,kb(1))=krg
	         keb(1,kb(1))=ks
	         ifg1=1
	       endif
c	       wb(kb(1))=wbt(i)
		   wb(kelb(1,kb(1)),keb(1,kb(1)))=wbt(i)	!wb(krg,ks)
	     endif
	   enddo
c
	   do i=1,n(4)	!pbt(n(4))
	     if(kbr(4,i).eq.krg .and. kbf(4,i).eq.ks)then
	       if(ifg1.eq.0)then
	         kb(1)=kb(1)+1
	         kelb(1,kb(1))=krg
	         keb(1,kb(1))=ks
	         ifg1=1
	       endif
c	       pb(kb(1))=pbt(i)
		   pb(kelb(1,kb(1)),keb(1,kb(1)))=pbt(i)	!pb(krg,ks)
	     endif
	   enddo
	   ifg1=0
c
c        kb(2),kelb(2,kb(2)),keb(2,kb(2))について
c        cb(kb(2))
	   ifg2=0
c
	   do i=1,n(6)	!cbt(n(6))
	     if(kbr(6,i).eq.krg .and. kbf(6,i).eq.ks)then
	       if(ifg2.eq.0)then
	         kb(2)=kb(2)+1
	         kelb(2,kb(2))=krg
	         keb(2,kb(2))=ks
	         ifg2=1
	       endif
c	       cb(kb(2))=cbt(i)
		   cb(kelb(2,kb(2)),keb(2,kb(2)))=cbt(i)	!cb(krg,ks)
	     endif
	   enddo
	   ifg2=0
c
c        kb(3),kelb(3,kb(3)),keb(3,kb(3))について
c        tb(kb(3))
	   ifg3=0
c
	   do i=1,n(7)	!tbt(n(7))
	     if(kbr(7,i).eq.krg .and. kbf(7,i).eq.ks)then
	       if(ifg3.eq.0)then
	         kb(3)=kb(3)+1
	         kelb(3,kb(3))=krg
	         keb(3,kb(3))=ks
	         ifg3=1
	       endif
c	       tb(kb(3))=tbt(i)
		   tb(kelb(3,kb(3)),keb(3,kb(3)))=tbt(i)	!tb(krg,ks)
	     endif
	   enddo
	   ifg3=0
c
c        kb(4),kelb(4,kb(4)),keb(4,kb(4))について
c        ug(kb(4)),vg(kb(4)),wg(kb(4))
	   ifg4=0
c
	   do i=1,n(8)	!ugt(n(8))
	     if(kbr(8,i).eq.krg .and. kbf(8,i).eq.ks)then
	       if(ifg4.eq.0)then
	         kb(4)=kb(4)+1
	         kelb(4,kb(4))=krg
	         keb(4,kb(4))=ks
	         ifg4=1
	       endif
c	       ug(kb(4))=ugt(i)
		   ug(kelb(4,kb(4)),keb(4,kb(4)))=ugt(i)	!ug(krg,ks)
	     endif
	   enddo
c
	   do i=1,n(9)	!vgt(n(9))
	     if(kbr(9,i).eq.krg .and. kbf(9,i).eq.ks)then
	       if(ifg4.eq.0)then
	         kb(4)=kb(4)+1
	         kelb(4,kb(4))=krg
	         keb(4,kb(4))=ks
	         ifg4=1
	       endif
c	       vg(kb(4))=vgt(i)
		   vg(kelb(4,kb(4)),keb(4,kb(4)))=vgt(i)	!vg(krg,ks)
	     endif
	   enddo
c
	   do i=1,n(10)	!wgt(n(10))
	     if(kbr(10,i).eq.krg .and. kbf(10,i).eq.ks)then
	       if(ifg4.eq.0)then
	         kb(4)=kb(4)+1
	         kelb(4,kb(4))=krg
	         keb(4,kb(4))=ks
	         ifg4=1
	       endif
c	       wg(kb(4))=wgt(i)
		   wg(kelb(4,kb(4)),keb(4,kb(4)))=wgt(i)	!wg(krg,ks)
	     endif
	   enddo
	   ifg4=0
c
c        kb(5),kelb(5,kb(5)),keb(5,kb(5))について
c        cg(kb(5))
	   ifg5=0
c
	   do i=1,n(11)	!cgt(n(11))
	     if(kbr(11,i).eq.krg .and. kbf(11,i).eq.ks)then
	       if(ifg5.eq.0)then
	         kb(5)=kb(5)+1
	         kelb(5,kb(5))=krg
	         keb(5,kb(5))=ks
	         ifg5=1
	       endif
c	       cg(kb(5))=cgt(i)
		   cg(kelb(5,kb(5)),keb(5,kb(5)))=cgt(i)	!cg(krg,ks)
	     endif
	   enddo
	   ifg5=0
c
c        kb(6),kelb(6,kb(6)),keb(6,kb(6))について
c        tg(kb(6))
	   ifg6=0
c
	   do i=1,n(12)	!tgt(n(12))
	     if(kbr(12,i).eq.krg .and. kbf(12,i).eq.ks)then
	       if(ifg6.eq.0)then
	         kb(6)=kb(6)+1
	         kelb(6,kb(6))=krg
	         keb(6,kb(6))=ks
	         ifg6=1
	       endif
c	       tg(kb(6))=tgt(i)
		   tg(kelb(6,kb(6)),keb(6,kb(6)))=tgt(i)	!tg(krg,ks)
	     endif
	   enddo
	   ifg6=0
c
c        kb(7),kelb(7,kb(7)),keb(7,kb(7))について
c        ch(kb(7))
	   ifg7=0
c
	   do i=1,n(13)	!cht(n(13))
	     if(kbr(13,i).eq.krg .and. kbf(13,i).eq.ks)then
	       if(ifg7.eq.0)then
	         kb(7)=kb(7)+1
	         kelb(7,kb(7))=krg
	         keb(7,kb(7))=ks
	         ifg7=1
	       endif
c	       ch(kb(7))=cht(i)
		   ch(kelb(7,kb(7)),keb(7,kb(7)))=cht(i)	!ch(krg,ks)
	     endif
	   enddo
	   ifg7=0
c
c        kb(8),kelb(8,kb(8)),keb(8,kb(8))について
c        th(kb(8))
	   ifg8=0
c
	   do i=1,n(14)	!tht(n(14))
	     if(kbr(14,i).eq.krg .and. kbf(14,i).eq.ks)then
	       if(ifg8.eq.0)then
	         kb(8)=kb(8)+1
	         kelb(8,kb(8))=krg
	         keb(8,kb(8))=ks
	         ifg8=1
	       endif
c	       th(kb(8))=tht(i)
		   th(kelb(8,kb(8)),keb(8,kb(8)))=tht(i)	!th(krg,ks)
	     endif
	   enddo
	   ifg8=0
c
c        kb(9),kelb(9,kb(9)),keb(9,kb(9))について
c        scv(kb(9))
	   ifg9=0
c
	   do i=1,n(15)	!scvt(n(15))
	     if(kbr(15,i).eq.krg .and. kbf(15,i).eq.ks)then
	       if(ifg9.eq.0)then
	         kb(9)=kb(9)+1
	         kelb(9,kb(9))=krg
	         keb(9,kb(9))=ks
	         ifg9=1
	       endif
c	       scv(kb(9))=scvt(i)
		   scv(kelb(9,kb(9)),keb(9,kb(9)))=scvt(i)	!scv(krg,ks)
	     endif
	   enddo
	   ifg9=0
c
c        kb(10),kelb(10,kb(10)),keb(10,kb(10))について
c        scm(kb(10))
	   ifg10=0
c
	   do i=1,n(16)	!scmt(n(16))
	     if(kbr(16,i).eq.krg .and. kbf(16,i).eq.ks)then
	       if(ifg10.eq.0)then
	         kb(10)=kb(10)+1
	         kelb(10,kb(10))=krg
	         keb(10,kb(10))=ks
	         ifg10=1
	       endif
c	       scm(kb(10))=scmt(i)
		   scm(kelb(10,kb(10)),keb(10,kb(10)))=scmt(i)	!scm(krg,ks)
	     endif
	   enddo
	   ifg10=0
c
c        kb(11),kelb(11,kb(11)),keb(11,kb(11))について
c        sct(kb(11))
	   ifg11=0
c
	   do i=1,n(17)	!sctt(n(17))
	     if(kbr(17,i).eq.krg .and. kbf(17,i).eq.ks)then
	       if(ifg11.eq.0)then
	         kb(11)=kb(11)+1
	         kelb(11,kb(11))=krg
	         keb(11,kb(11))=ks
	         ifg11=1
	       endif
c	       sct(kb(11))=sctt(i)
		   sct(kelb(11,kb(11)),keb(11,kb(11)))=sctt(i)	!sct(krg,ks)
	     endif
	   enddo
	   ifg11=0
c
3900   continue		!krg,ks
c       write(*,*)'made kb(11),kelb(11,kb(11)),keb(11,kb(11))'
c
c       do i=1,11
c	   write(*,*)'----------------------------------------------'
c	   write(*,*)'kb(',i,')',kb(i)
c	   do j=1,kb(i)
c           write(*,*)j,kelb(i,j),keb(i,j)
c	   enddo
c	   pause
c	 enddo
c----------------------------------------------------------
c
c      nb(11),nelb(11,nb(11)),neb(11,nb(11))の作成
c      kelb(11,kb(11)):境界条件を代入する6面体の番号
c      keb(11,kb(11)):境界条件を代入する6面体の面番号
c      mgr(nelw):4面体要素nelwが存在する領域の番号
c----------------------------------------------------------
c
c      初期化
c
       numb1=0		!境界条件の数
	 numb2=0
	 numb3=0
	 numb4=0
	 numb5=0
	 numb6=0
	 numb7=0
	 numb8=0
	 numb9=0
	 numb10=0
	 numb11=0
       open(1,file='bd.inf')
	 rewind(1)
c
       do 4000 kk=1,11			!境界条件の種類
       do 4100 lb=1,kb(kk)		!kb(11):境界条件の数←境界条件は領域の面ごと
         krg=kelb(kk,lb)	    !kelb(11,kb(11)):6面体の領域の番号
	   ks=keb(kk,lb)			!keb(11,kb(11)):6面体の領域の面の番号
c
c        平面上に存在する節点の番号を格納した配列kmp(kn)の作成
c
	   kn=0
c
c        面1:Z-X平面←6面体の領域
c
         if(ks.eq.1)then
	     ky=1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwy(i).eq.ky)then
	         kn=kn+1
	         kmp(kn)=lwp(i)		!kmp(kn):X-Z平面上に存在する節点の番号を格納した配列
	       endif
	     enddo
c
c        面2:Y-Z平面←6面体の領域
c
	   elseif(ks.eq.2)then
	     kx=mdx(krg)+1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwx(i).eq.kx)then
	         kn=kn+1
	         kmp(kn)=lwp(i)
	       endif
	     enddo
c        面3:X-Y平面←6面体の領域
         elseif(ks.eq.3)then
	     kz=mdz(krg)+1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwz(i).eq.kz)then
	         kn=kn+1
	         kmp(kn)=lwp(i)
	       endif
	     enddo
c        面4:Y-Z平面←6面体の領域
         elseif(ks.eq.4)then
	     kx=1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwx(i).eq.kx)then
	         kn=kn+1
	         kmp(kn)=lwp(i)
	       endif
	     enddo
c        面5:X-Y平面←6面体の領域
         elseif(ks.eq.5)then
	     kz=1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwz(i).eq.kz)then
	         kn=kn+1
	         kmp(kn)=lwp(i)
	       endif
	     enddo
c        面6:Z-X平面←6面体の領域
         elseif(ks.eq.6)then
	     ky=mdy(krg)+1
	     do i=1,ii
	       if(lwr(i).eq.krg .and. lwy(i).eq.ky)then
	         kn=kn+1
	         kmp(kn)=lwp(i)
	       endif
	     enddo
c
	   endif
c	   write(*,*)'krg,ks,kn',krg,ks,kn
c	   do i=1,kn
c	   write(*,*)'kmp(',i,')',kmp(i)
c	   enddo
c	   pause
c
c        nb(kk),nelb(kk,nb(kk)),neb(kk,nb(kk))を決定する
c        領域番号krg,面番号ksについて
         do i=1,11
  	     nb(i)=0
      	 enddo
c
	   do 4200 kelw=1,nelw
           nf(1)=0		!0:節点が平面上に存在しない,1:節点が平面上に存在する
	     nf(2)=0
	     nf(3)=0
	     nf(4)=0
           do k=1,kn
             if(new(kelw,1).eq.kmp(k))nf(1)=1
	       if(new(kelw,2).eq.kmp(k))nf(2)=1
	       if(new(kelw,3).eq.kmp(k))nf(3)=1
	       if(new(kelw,4).eq.kmp(k))nf(4)=1
           enddo
c          面1:j=1,2,3←4面体要素
           if(nf(1).eq.1 .and. nf(2).eq.1 .and. nf(3).eq.1)then
c	       write(*,*)'surface number 1'
	       nb(kk)=nb(kk)+1
	       nelb(kk,nb(kk))=kelw
	       neb(kk,nb(kk))=1
c          面2:j=2,3,4←4面体要素
           elseif(nf(2).eq.1 .and. nf(3).eq.1 .and. nf(4).eq.1)then
c	       write(*,*)'surface number 2'
	       nb(kk)=nb(kk)+1
	       nelb(kk,nb(kk))=kelw
	       neb(kk,nb(kk))=2
c          面3:j=3,4,1←4面体要素
           elseif(nf(3).eq.1 .and. nf(4).eq.1 .and. nf(1).eq.1)then
c	       write(*,*)'surface number 3'
	       nb(kk)=nb(kk)+1
	       nelb(kk,nb(kk))=kelw
	       neb(kk,nb(kk))=3
c          面4:j=4,1,2←4面体要素
           elseif(nf(4).eq.1 .and. nf(1).eq.1 .and. nf(2).eq.1)then
c	       write(*,*)'surface number 4'
	       nb(kk)=nb(kk)+1
	       nelb(kk,nb(kk))=kelw
	       neb(kk,nb(kk))=4
	     endif
4200     continue
c         write(*,*)'made nb(11),nelb(11,nb(11)),neb(11,nb(11))'
c
c        bd.infファイルを作成する
c        krg=kelb(kk,lb)			!kelb(11,kb(11)):6面体の領域の番号
c     	 ks=keb(kk,lb)				!keb(11,kb(11)):6面体の領域の面の番号
         do 4300 i=1,nb(kk)			!nb(kk):領域番号krg,面番号ks上の4面体要素の数
	     if(kk.eq.1)then
c            ub(krg,ks),vb(krg,ks),wb(krg,ks),pb(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),ub(krg,ks),vb(krg,ks)
     &                                      ,wb(krg,ks),pb(krg,ks)
 	       numb1=numb1+1	!境界条件の数
	     elseif(kk.eq.2)then
c            cb(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),cb(krg,ks)
	       numb2=numb2+1	!境界条件の数
           elseif(kk.eq.3)then
c            tb(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),tb(krg,ks)
	       numb3=numb3+1	!境界条件の数
	     elseif(kk.eq.4)then
c            ug(krg,ks),vg(krg,ks),wg(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),ug(krg,ks),vg(krg,ks)
     &                                                 ,wg(krg,ks)
	       numb4=numb4+1	!境界条件の数
	     elseif(kk.eq.5)then
c            cg(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),cg(krg,ks)
	       numb5=numb5+1	!境界条件の数
           elseif(kk.eq.6)then
c            tg(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),tg(krg,ks)
	       numb6=numb6+1	!境界条件の数
	     elseif(kk.eq.7)then
c            ch(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),ch(krg,ks)
	       numb7=numb7+1	!境界条件の数
           elseif(kk.eq.8)then
c	       th(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),th(krg,ks)
	       numb8=numb8+1	!境界条件の数
	     elseif(kk.eq.9)then
c            scv(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),scv(krg,ks)
	       numb9=numb9+1	!境界条件の数
           elseif(kk.eq.10)then
c            scm(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),scm(krg,ks)
	       numb10=numb10+1	!境界条件の数
	     elseif(kk.eq.11)then
c            sct(krg,ks)
             write(1,*)nelb(kk,i),neb(kk,i),sct(krg,ks)
	       numb11=numb11+1	!境界条件の数
	     endif
4300     continue
c
4100   continue		!kb(11):境界条件の数←境界条件は領域の面ごと
4000   continue		!kk:境界条件の種類
c
       close(1)
       write(*,*)'created bd.inf'
c
       open(1,file='bdnum.inf')
	 rewind(1)
	 write(1,*)numb1
	 write(1,*)numb2
	 write(1,*)numb3
	 write(1,*)numb4
	 write(1,*)numb5
	 write(1,*)numb6
	 write(1,*)numb7
	 write(1,*)numb8
	 write(1,*)numb9
	 write(1,*)numb10
	 write(1,*)numb11
	 close(1)
       write(*,*)'created bdnum.inf'
c
       stop
       end
c-------------------------------------------------------------------------
c
c      サブルーチン
c
c-------------------------------------------------------------------------
c      境界条件を配列に代入する
c       subroutine setbnd(km,ub,vb,wb,pb,cb,tb,ug,vg,wg,cg,tg,
c     &                         ch,th,scv,scm,sct,ind)
c
       subroutine straight(xt,yt,zt,nm,xc,yc,zc)
       implicit double precision(a-h,o-z)
       dimension xt(3),yt(3),zt(3)
       dimension xc(1001),yc(1001),zc(1001)  !分割した各節点の座標
c
       x1=xt(1)
       y1=yt(1)
       z1=zt(1)
       x2=xt(2)
       y2=yt(2)
       z2=zt(2)
c
       do 1000 km=1,nm+1
         xc(km)=x1+(x2-x1)*dble(km-1)/dble(nm)
         yc(km)=y1+(y2-y1)*dble(km-1)/dble(nm)
         zc(km)=z1+(z2-z1)*dble(km-1)/dble(nm)
1000   continue
       return
       end
c
       subroutine curve(xt,yt,zt,nm,xc,yc,zc)
       implicit double precision(a-h,o-z)
       dimension xt(3),yt(3),zt(3)
       dimension xc(1001),yc(1001),zc(1001)
c
       x1=xt(1)
       y1=yt(1)
       z1=zt(1)
       x2=xt(2)
       y2=yt(2)
       z2=zt(2)
       x3=xt(3)
       y3=yt(3)
       z3=zt(3)
c      z=aa*x+bb*y+cc
       aa=((z3-z1)*(y2-y1)-(z2-z1)*(y3-y1))
     &   /((y2-y1)*(x3-x1)-(y3-y1)*(x2-x1))
       bb=((z2-z1)*(x3-x1)-(z3-z1)*(x2-x1))
     &   /((y2-y1)*(x3-x1)-(y3-y1)*(x2-x1))
       cc=z1-aa*x1-bb*y1
c
       t1=(x3**2-x1**2+y3**2-y1**2+z3**2-z1**2-2d0*cc*(z3-z1))
     &   /2d0/(x3-x1+aa*(z3-z1))
       t2=(x2**2-x1**2+y2**2-y1**2+z2**2-z1**1-2d0*cc*(z2-z1))
     &   /2d0/(x2-x1+aa*(z2-z1))
       s1=(y3-y1+bb*(z3-z1))/(x3-x1+aa*(z3-z1))
       s2=(y2-y1+bb*(z2-z1))/(x2-x1+aa*(z2-z1))
       y0=(t1-t2)/(s1-s2)
       x0=1d0/2d0*(t2+t1-(s2+s1)*y0)
       z0=aa*x0+bb*y0+cc
c
       shita0=dacos(((x1-x0)*(x2-x0)+(y1-y0)*(y2-y0)+(z1-z0)*(z2-z0))
     &       /dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
     &       /dsqrt((x2-x0)**2+(y2-y0)**2+(z2-z0)**2))
       xmax=dmax1(x1,x2,x3)
       xmin=dmin1(x1,x2,x3)
       ymax=dmax1(y1,y2,y3)
       ymin=dmin1(y1,y2,y3)
       zmax=dmax1(z1,z2,z3)
       zmin=dmin1(z1,z2,z3)
       dx=1d-7   !xの増分
       dy=1d-7   !yの増分
       dz=1d-7   !zの増分
       dd=1d-5   !終了条件
       ro=dsqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
c
       shita=0d0
       do 1000 km=1,nm+1
         shita=shita+shita0*dble(km-1)/dble(nm)
         x=(xmax+xmin)/2d0
         y=(ymax+ymin)/2d0
         z=(zmax+zmin)/2d0
         do 2000 m=1,100000
           call calf(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,df,dd)
           call calg(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dg,dd,shita)
           call calh(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dh,dd,shita,shita0)
           if(df.lt.dd .and. dg.lt.dd .and. dh.lt.dd)goto 3000
2000     continue
3000     continue
         xc(km)=x
         yc(km)=y
         zc(km)=z
1000   continue
c
       return
       end
c
       subroutine calf(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,df,dd)
       implicit double precision(a-h,o-z)
       f0=funcf(x,y,z,r0,x0,y0,z0)
       do 1000 m=1,100000
c
c        X軸方向
c
         f1=funcf(x-dx,y,z,r0,x0,y0,z0)
         f2=funcf(x+dx,y,z,r0,x0,y0,z0)
         if(f1.lt.f0)then
           x=x-dx
           if(x.lt.xmin)x=x+dx
         elseif(f2.lt.f0)then
           x=x+dx
           if(x.gt.xmax)x=x-dx
         endif
c
c        Y軸方向
c
         f1=funcf(x,y-dy,z,r0,x0,y0,z0)
         f2=funcf(x,y+dy,z,r0,x0,y0,z0)
         if(f1.lt.f0)then
           y=y-dy
           if(y.lt.ymin)y=y+dy
         elseif(f2.lt.f0)then
           y=y+dy
           if(y.gt.ymax)y=y-dy
         endif
c
c        Z軸方向
c
         f1=funcf(x,y,z-dz,r0,x0,y0,z0)
         f2=funcf(x,y,z+dz,r0,x0,y0,z0)
         if(f1.lt.f0)then
           z=z-dz
           if(z.lt.zmin)z=z+dz
         elseif(f2.lt.f0)then
           z=z+dz
           if(z.gt.zmax)z=z-dz
         endif
         f0=funcf(x,y,z,r0,x0,y0,z0)
         if(f0.lt.dd)goto 2000  !終了判定
1000   continue
2000   continue
       df=f0
       return
       end
c
       subroutine calg(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dg,dd,shita)
       implicit double precision(a-h,o-z)
       g0=funcg(x,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
       do 1000 m=1,100000
c
c        X軸方向
c
         g1=funcg(x-dx,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
         g2=funcg(x+dx,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
         if(g1.lt.g0)then
           x=x-dx
           if(x.lt.xmin)x=x+dx
         elseif(g2.lt.g0)then
           x=x+dx
           if(x.gt.xmax)x=x-dx
         endif
c
c        Y軸方向
c
         g1=funcg(x,y-dy,z,r0,shita,x0,y0,z0,x1,y1,z1)
         g2=funcg(x,y+dy,z,r0,shita,x0,y0,z0,x1,y1,z1)
         if(g1.lt.g0)then
           y=y-dy
           if(y.lt.ymin)y=y+dy
         elseif(g2.lt.g0)then
           y=y+dy
           if(y.gt.ymax)y=y-dy
         endif
c
c        Z軸方向
c
         g1=funcg(x,y,z-dz,r0,shita,x0,y0,z0,x1,y1,z1)
         g2=funcg(x,y,z+dz,r0,shita,x0,y0,z0,x1,y1,z1)
         if(g1.lt.g0)then
           z=z-dz
           if(z.lt.zmin)z=z+dz
         elseif(g2.lt.g0)then
           z=z+dz
           if(z.gt.zmax)z=z-dz
         endif
         g0=funcg(x,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
         if(g0.lt.dd)goto 2000  !終了判定
1000   continue
2000   continue
       dg=g0
       return
       end
c
       subroutine calh(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dh,dd,shita,shita0)
       implicit double precision(a-h,o-z)
       h0=funch(x,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
       do 1000 m=1,100000
c
c        X軸方向
c
         h1=funch(x-dx,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         h2=funch(x+dx,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         if(h1.lt.h0)then
           x=x-dx
           if(x.lt.xmin)x=x+dx
         elseif(h2.lt.h0)then
           x=x+dx
           if(x.gt.xmax)x=x-dx
         endif
c
c        Y軸方向
c
         h1=funch(x,y-dy,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         h2=funch(x,y+dy,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         if(h1.lt.h0)then
           y=y-dy
           if(y.lt.ymin)y=y+dy
         elseif(h2.lt.h0)then
           y=y+dy
           if(y.gt.ymax)y=y-dy
         endif
c
c        Z軸方向
c
         h1=funch(x,y,z-dz,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         h2=funch(x,y,z+dz,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         if(h1.lt.h0)then
           z=z-dz
           if(z.lt.zmin)z=z+dz
         elseif(h2.lt.h0)then
           z=z+dz
           if(z.gt.zmax)z=z-dz
         endif
         h0=funch(x,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         if(h0.lt.dd)goto 2000  !終了判定
1000   continue
2000   continue
       dh=h0
       return
       end
c
c      function
c
       function nmd(il,ndx,ndy,ndz)
       implicit double precision(a-h,o-z)
       if(il.eq.1)then
         nmd=ndx
       elseif(il.eq.2)then
         nmd=ndz
       elseif(il.eq.3)then
         nmd=ndx
       elseif(il.eq.4)then
         nmd=ndz
       elseif(il.eq.5)then
         nmd=ndy
       elseif(il.eq.6)then
         nmd=ndy
       elseif(il.eq.7)then
         nmd=ndy
       elseif(il.eq.8)then
         nmd=ndy
       elseif(il.eq.9)then
         nmd=ndx
       elseif(il.eq.10)then
         nmd=ndz
       elseif(il.eq.11)then
         nmd=ndx
       elseif(il.eq.12)then
         nmd=ndz
       endif
       return
       end
c
       function funcf(x,y,z,r0,x0,y0,z0)
       implicit double precision(a-h,o-z)
       funcf=(x-x0)**2+(y-y0)**2+(z-z0)**2-r0
       return
       end
c
       function funcg(x,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
       implicit double precision(a-h,o-z)
       funcg=(x-x0)*(x1-x0)+(y-y0)*(y1-y0)+(z-z0)*(z1-z0)-r0*dcos(shita)
       return
       end
c
       function funch(x,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
       implicit double precision(a-h,o-z)
       funch=(x-x0)*(x2-x0)+(y-y0)*(y2-y0)+(z-z0)*(z2-z0)
     &      -r0*dcos(shita0-shita)
       return
       end
c      領域番号lrgの6面体の辺番号llの辺を構成する始点の節点番号を求める
       function icals(mp,lrg,ll)
       implicit double precision(a-h,o-z)
       dimension mp(1001,8)        !領域の各頂点の節点番号
       if(ll.eq.1)then
         icals=mp(lrg,1)
       elseif(ll.eq.2)then
         icals=mp(lrg,2)
       elseif(ll.eq.3)then
         icals=mp(lrg,3)
       elseif(ll.eq.4)then
         icals=mp(lrg,4)
       elseif(ll.eq.5)then
         icals=mp(lrg,2)
       elseif(ll.eq.6)then
         icals=mp(lrg,3)
       elseif(ll.eq.7)then
         icals=mp(lrg,4)
       elseif(ll.eq.8)then
         icals=mp(lrg,1)
       elseif(ll.eq.9)then
         icals=mp(lrg,5)
       elseif(ll.eq.10)then
         icals=mp(lrg,6)
       elseif(ll.eq.11)then
         icals=mp(lrg,7)
       elseif(ll.eq.12)then
         icals=mp(lrg,8)
       endif
       return
       end
c      領域番号lrgの6面体の辺番号llの辺を構成する終点の節点番号を求める
       function icale(mp,lrg,ll)
       implicit double precision(a-h,o-z)
       dimension mp(1001,8)        !領域の各頂点の節点番号
       if(ll.eq.1)then
         icale=mp(lrg,2)
       elseif(ll.eq.2)then
         icale=mp(lrg,3)
       elseif(ll.eq.3)then
         icale=mp(lrg,4)
       elseif(ll.eq.4)then
         icale=mp(lrg,1)
       elseif(ll.eq.5)then
         icale=mp(lrg,6)
       elseif(ll.eq.6)then
         icale=mp(lrg,7)
       elseif(ll.eq.7)then
         icale=mp(lrg,8)
       elseif(ll.eq.8)then
         icale=mp(lrg,5)
       elseif(ll.eq.9)then
         icale=mp(lrg,6)
       elseif(ll.eq.10)then
         icale=mp(lrg,7)
       elseif(ll.eq.11)then
         icale=mp(lrg,8)
       elseif(ll.eq.12)then
         icale=mp(lrg,5)
       endif
       return
       end
c      面を共有する6面体の面番号を求める
       function isf(ks)
       implicit double precision(a-h,o-z)
       if(ks.eq.1)isf=6
       if(ks.eq.2)isf=4
       if(ks.eq.3)isf=5
       if(ks.eq.4)isf=2
       if(ks.eq.5)isf=3
       if(ks.eq.6)isf=1
       return
       end