c----------------------------------------------------------------
c
c      グリッドファイル,'.1','.2',作成関数
c      1つの領域における最大分割数:99,最大領域数:100
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
       implicit double precision(a-h,o-z)
       parameter(ind=9999,ine=9999,ibw=999,irg=9999)
c       
       dimension xl(12,1001),yl(12,1001),zl(12,1001)         !辺上に存在する節点の座標
       dimension xc(1001),yc(1001),zc(1001)                  !計算した辺上に存在する節点の座標
       dimension xt(1001),yt(1001),zt(1001)                  !辺上に存在する節点の座標
c
       dimension xs1(1001,1001),ys1(1001,1001),zs1(1001,1001)     !面上に存在する節点の座標
       dimension xs2(1001,1001),ys2(1001,1001),zs2(1001,1001)     !面上に存在する節点の座標
       dimension xs3(1001,1001),ys3(1001,1001),zs3(1001,1001)     !面上に存在する節点の座標
       dimension xs4(1001,1001),ys4(1001,1001),zs4(1001,1001)     !面上に存在する節点の座標
       dimension xs5(1001,1001),ys5(1001,1001),zs5(1001,1001)     !面上に存在する節点の座標
       dimension xs6(1001,1001),ys6(1001,1001),zs6(1001,1001)     !面上に存在する節点の座標
c
       integer*8 krr(10,10,10)                                 !xr(krg,krr(kx,ky,kz))
c--------------------------------------------------------------------------------------------       
       integer*8 xr(100,1000),yr(100,1000),zr(100,1000)           !各領域内を構成する6面体の頂点座標
c--------------------------------------------------------------------------------------------       
c
       dimension ifline(irg,12)                    !0:直線,1:曲線
       dimension mdx(1001),mdy(1001),mdz(1001)     !各領域の分割数
       dimension mk(1001,6)                        !領域の各面を共有している領域の番号
       dimension mp(1001,8)                        !領域の各頂点の節点番号
       dimension mx(1001),my(1001),mz(1001)        !分割数を変更する各領域の番号
       dimension xp(1001,12,3),yp(1001,12,3),zp(1001,12,3)  !辺を構成する節点の座標
c
c--------------------------------------------------------------------------------------------
       integer*8 kprg(10,10,10),nprg(100,1000) !領域内のある位置における節点の番号:nprg(krg,kprg(kdx,kdy,kdz))
c--------------------------------------------------------------------------------------------       
       dimension mprg(ind),npx(1001,ind),npy(1001,ind),npz(1001,ind) !各節点の領域内の位置
       dimension xx(ind),yy(ind),zz(ind)           !各節点の座標
       dimension mr(irg),md(ine)                   !各領域および6面体のindex番号
       dimension ne(ine,8)                         !各6面体を構成する節点の番号
c
       dimension idx(ine)                          !各要素のindex番号
       dimension new(ine,4)                        !各4面体を構成する節点の番号
       dimension jww(ind),jmm(ind)                 !jww(古)=新,jmm(古)=新
c
       character*128 rfile,pfile,gfile1,gfile2      !file names
c
c      infファイルを読み込む
c
       open(1,file='input.txt')
       rewind(1)
       read(1,*)pfile
       close(1)
       write(*,*)'the name of inputted file,',pfile
       pause
c
       open(1,file=pfile)
       rewind(1)
       read(1,*)nrg    !領域数
       do 100 krg=1,nrg
         read(1,*)mr(krg)  !領域krgのインデックス番号
         read(1,*)mdx(krg),mdy(krg),mdz(krg)           !各軸方向の分割数
         do 200 il=1,12
           read(1,*)ifline(krg,il)                     !辺の種類(0:直線,1:曲線)
           nside=2
           if(ifline(krg,il).eq.1)nside=3
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
c      各領域を構成する節点の座標を求める
c
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
             xt(1)=xp(krg,il,1)
             yt(1)=yp(krg,il,1)
             zt(1)=zp(krg,il,1)
             xt(2)=xp(krg,il,2)
             yt(2)=yp(krg,il,2)
             zt(2)=zp(krg,il,2)
             call straight(xt,yt,zt,nm,xc,yc,zc)
             do 1300 km=1,nm+1   !考慮している辺上に存在する各節点の座標
               xl(il,km)=xc(km)
               yl(il,km)=yc(km)
               zl(il,km)=zc(km)

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
         write(*,*)"ended deciding the nodes' positions of
     & each side of hexahedral elements"
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
     &                 *dble(kz-1)/dble(ndz)
           xsz=xl(5,ky)+(zl(6,ky)-zl(5,ky))
     &                 *dble(kz-1)/dble(ndz)
           xs2(ky,kx)=(xsy+xsz)/2d0
c
c          面4:Y-Z平面
c
           ys4(ky,kz)=yl(4,ndz+2-kz)+(yl(12,ndz+2-kz)-yl(4,ndz+2-kz))
     &                              *dble(ky-1)/dble(ndy)
           zs4(ky,kz)=zl(8,ky)+(zl(7,ky)-zl(8,ky))
     &                        *dble(kz-1)/dble(ndz)
           xsy=xl(4,ndz+2-kz)+(xl(12,ndz+2-kz)-xl(4,ndx+2-kz))
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
         write(*,*)"ended deciding the nodes' positions on each surface
     & of hexahedral elements"
c
c        6面体内部の各節点座標を求める
c
         do 1800 kx=1,ndx+1
         do 1800 ky=1,ndy+1
         do 1800 kz=1,ndz+1
           write(*,*)'kx,ky,kz',kx,ky,kz
c-----------------------------------------------------------------------------
c
           krr(kx,ky,kz)=kx+10*ky+100*kz
c
c-----------------------------------------------------------------------------
c
c          X座標:面4-2
c
           xr(krg,krr(kx,ky,kz))=xs4(ky,kz)+(xs2(ky,kz)-xs4(ky,kz))
     &                                     *dble(kx-1)/dble(ndx)
           write(*,*)'xr',sngl(xr(krg,krr(kx,ky,kz)))
c
c          Y座標:面1-6
c
           yr(krg,krr(kx,ky,kz))=ys1(kx,kz)+(ys6(kx,kz)-ys1(kx,kz))
     &                                     *dble(ky-1)/dble(ndy)
           write(*,*)'yr',yr(krg,krr(kx,ky,kz))
     
c
c          Z座標:面5-3
c
           zr(krg,krr(kx,ky,kz))=zs5(kx,ky)+(zs3(kx,ky)-zs5(kx,ky))
     &                                     *dble(kz-1)/dble(ndz)
           write(*,*)'zr',sngl(zr(krg,krr(kx,ky,kz)))
           pause
1800     continue
1100  continue
      write(*,*)"ended deciding the nodes' positions in each hexahedral
     & element"
      pause
c
c      各節点に番号を割り振る
c      nrg,mk(nrg,8),xr(nrg,ndx,ndy,ndz),yr(nrg,ndx,ndy,ndz),zr(nrg,ndx,ndy,ndz) →
c      np,nprg(nrg,ndx,ndy,ndz),mprg(np),npx(nrg,np),npy(nrg,np),npz(nrg,np),xx(np),yy(np),zz(np)
       np=0
       do 2300 krg=1,nrg
         do 2400 kdx=1,mdx(krg)+1     !ndx=mdx(nrg)
           if(mk(krg,2).lt.krg .and. kdx.eq.mdx(krg)+1) goto 2400 !面上の節点
           if(mk(krg,4).lt.krg .and. kdx.eq.1)goto 2400           !面上の節点
           do 2500 kdy=1,mdy(krg)+1     !ndy=mdy(nrg)
             if(mk(krg,6).lt.krg .and. kdy.eq.mdy(krg)+1)goto 2500  !面上の節点
             if(mk(krg,1).lt.krg .and. kdy.eq.1)goto 2500           !面上の節点
             do 2600 kdz=1,mdz(krg)+1     !ndz=mdz(nrg)
               if(mk(krg,3).lt.krg .and. kdz.eq.mdz(krg)+1)goto 2600  !面上の節点
               if(mk(krg,5).lt.krg .and. kdz.eq.1)goto 2600           !面上の節点
               np=np+1
c---------------------------------------------------------------------------------
c               
               kprg(kdx,kdy,kdz)=kdx+10*kdy+100*kdz
c
c---------------------------------------------------------------------------------
               nprg(krg,kprg(kdx,kdy,kdz))=np
               mprg(np)=krg
               npx(krg,np)=kdx
               npy(krg,np)=kdy
               npz(krg,np)=kdz
c
               xx(np)=xr(krg,krr(kdx,kdy,kdz))
               yy(np)=yr(krg,krr(kdx,kdy,kdz))
               zz(np)=zr(krg,krr(kdx,kdy,kdz))
2600         continue
2500       continue
2400     continue
2300   continue 
      write(*,*)"ended allotting the number of nodes to each vertexes"
      pause
c
c      6面体を4面体に分割する
c
       nele=0
       do 2700 krg=1,nrg
         do 2800 kdx=1,mx(krg)
           do 2900 kdy=1,my(krg)
             do 3000 kdz=1,mz(krg)
               nele=nele+1
               ne(nele,1)=nprg(krg,kprg(kdx,kdy,kdz))
               ne(nele,2)=nprg(krg,kprg(kdx+1,kdy,kdz))
               ne(nele,3)=nprg(krg,kprg(kdx+1,kdy,kdz+1))
               ne(nele,4)=nprg(krg,kprg(kdx,kdy,kdz+1))
               ne(nele,5)=nprg(krg,kprg(kdx,kdy+1,kdz))
               ne(nele,6)=nprg(krg,kprg(kdx+1,kdy+1,kdz))
               ne(nele,7)=nprg(krg,kprg(kdx+1,kdy+1,kdz+1))
               ne(nele,8)=nprg(krg,kprg(kdx,kdy+1,kdz+1))
               md(nele)=mr(krg)
3000         continue
2900       continue
2800     continue
2700   continue
       write(*,*)"ended putting away the number of nodes in the
     & arrangements of the hexahedral element,ne(ind,8)"
       pause
c
c      4面体に分割する
c
       nelw=0
       do 3100 kele=1,nele     !6面体
         new(nelw,1)=ne(kele,6)
         new(nelw,2)=ne(kele,5)
         new(nelw,3)=ne(kele,8)
         new(nelw,4)=ne(kele,1)
         idx(nelw)=md(kele)
c
         nelw=nelw+1
         new(nelw,1)=ne(kele,4)
         new(nelw,2)=ne(kele,1)
         new(nelw,3)=ne(kele,2)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
         nelw=nelw+1
         new(nelw,1)=ne(kele,8)
         new(nelw,2)=ne(kele,1)
         new(nelw,3)=ne(kele,4)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
         nelw=nelw+1
         new(nelw,1)=ne(kele,7)
         new(nelw,2)=ne(kele,6)
         new(nelw,3)=ne(kele,8)
         new(nelw,4)=ne(kele,3)
         idx(nelw)=md(kele)
c
         nelw=nelw+1
         new(nelw,1)=ne(kele,3)
         new(nelw,2)=ne(kele,4)
         new(nelw,3)=ne(kele,2)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
c
         nelw=nelw+1
         new(nelw,1)=ne(kele,3)
         new(nelw,2)=ne(kele,8)
         new(nelw,3)=ne(kele,4)
         new(nelw,4)=ne(kele,6)
         idx(nelw)=md(kele)
3100   continue
       write(*,*)"ended putting away the number of nodes in the
     & arrangements of the tetrahedral element,new(ine,4)"
       pause
c
c      バンド幅を最適化
c
       npw=np
       call renum3(npw,nelw,new,jww,jmm,nbw,ind,ine,ibw)
       write(*,*)"ended suiting the number of nodes"
       pause       
c
c      *.1および*,2ファイルの作成
c
       call addchr(gfile1,rfile,'.1')
       open(1,file=grile1)
       rewind(1)
       write(1,*)npw
       do 3300 i=1,npw
3300   write(1,*)i,jww(i),xx(i),yy(i),zz(i)
       close(1)
c
       call addchr(gfile2,rfile,'.2')
       open(2,file=gfile2)
       rewind(1)
       write(1,*)nelw,npw,nbw
       do 3400 i=1,nelew
3400   write(1,*)i,(jww(new(i,j)),j=1,4)
     &             ,(xx(new(i,j)),yy(new(i,j)),zz(new(i,j)),j=1,4)
     &             ,idx(nelw)
       close(1)
       write(*,*)"ended creating the file,'.1','.2'"
       pause       
c
       stop
       end
c-------------------------------------------------------------------------
c
c      サブルーチン
c
c-------------------------------------------------------------------------
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
     &               ,r0,df)
           call calg(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dg,shita)
           call calh(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dh,shita,shita0)
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
     &               ,r0,df)
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
           z=d+dz
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
     &               ,r0,dg,shita)
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
           z=d+dz
           if(z.gt.zmax)z=z-dz
         endif
         g0=funcg(x,y,z,r0,shita,x0,y0,z0,x1,y1,z1)
         if(f0.lt.dd)goto 2000  !終了判定
1000   continue
2000   continue
       dg=g0
       return
       end
c
       subroutine calh(x,y,z,xmax,ymax,zmax,xmin,ymin,zmin
     &               ,dx,dy,dz,x1,x2,x3,y1,y2,y3,z1,z2,z3
     &               ,r0,dh,shita,shita0)
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
           z=d+dz
           if(z.gt.zmax)z=z-dz
         endif
         h0=funch(x,y,z,r0,shita,shita0,x0,y0,z0,x2,y2,z2)
         if(f0.lt.dd)goto 2000  !終了判定
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
