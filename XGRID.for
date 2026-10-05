c      引数
c      lrg:分割数を変更する領域の1つの領域番号
c      ldx,ldy,ldz:領域lrgの分割数
c      nrg:領域数
c      mp(nrg,8):領域の各頂点の節点番号
c      x(np),y(np),z(np):各頂点の座標
c      mr(nrg):index番号
       include 'addchr.for'       
       implicit double precision(a-h,o-z)
       parameter(ind=9999,ine=9999,ibw=999)
       character*20 rfile,pfile,gfile1,gfile2
       dimension mlrg(1001),mldx(1001),mldy(1001),mldz(1001) !分割数を変更する領域の数および分割数
       dimension mdx(1001),mdy(1001),mdz(1001)     !各領域の分割数
       dimension mk(1001,6)                        !領域の各面を共有している領域の番号
       dimension mx(1001),my(1001),mz(1001)        !分割数を変更する各領域の番号
       dimension kx(1001),ky(1001),kz(1001)        !見つかった分割数を変更する領域の番号
       dimension xr(1001,1001,1001,1001),yr(1001,1001,1001)
     &           ,zr(1001,1001,1001,1001)           !各領域内を構成する6面体の頂点座標
       dimension xp(1001,8),yp(1001,8),zp(1001,8)  !各領域の頂点座標
       dimension mp(1001,8)                        !領域の各頂点の節点番号
       dimension x(1001),y(1001),z(1001)           !各頂点の座標
       dimension mr(nrg),mdx(ine)                  !各領域および6面体のindex番号
       dimension ne(ine,8)                         !各6面体を構成する節点の番号
c
       dimension nprg(1001,1001,1001,1001)         !領域内のある位置における節点の番号
       dimension mprg(ind),npx(1001,ind),npy(1001,ind),npz(1001,ind) !各節点の領域内の位置
       dimension xx(ind),yy(ind),zz(ind)           !各節点の座標
       dimension idx(ine)                          !各要素のindex番号
       dimension new(ine,4)                        !各4面体を構成する節点の番号
c
c      infファイルを読み込む
c
       write(*,*)'inpur the name of informatoin file'
       read(*,*)rfile
       call addchr(pfile,rfile,'.inf')
       open(1,file=pfile)
       rewind(1)
       read(1,*)nrg,lp,iflag,nnn
       read(1,*)(x(i),i=1,lp)
       read(1,*)(y(i),i=1,lp)
       read(1,*)(z(i),i=1,lp)
       do 100 krg=1,nrg
100    read(1,*)lrg,(mk(lrg,i),i=1,6)
       do 200 krg=1,nrg
200    read(1,*)lrg,(mp(lrg,i),i=1,8),mr(lrg)
       read(1,*)nmlg,inid     !分割数を変更する領域の数およびデフォルトの分割数
       do 300 kmlg=1,nmlg
300    read(1,*)mlrg(kmlg),mldx(kmlg),mldy(kmlg),mldz(kmlg)
       close(1)
c
c      分割数の初期化
c
       do 1000 krg=1,nrg
         mdx(krg)=inid
         mdy(krg)=inid
         mdz(krg)=inid
1000   continue
c
c      分割数を変更する:mdx(nrg),mdy(nrg),mdz(nrg)を求める
c
       do 1100 kmlg=1,nmlg
         lrg=mlrg(kmlg)
         ldx=mldx(kmlg)
         ldy=mldy(kmlg)
         ldz=mldz(kmlg)
c       
c        x軸方向:nx,mx(nx),mdx(mx(nx))を求める
         nx=1      !分割数を変更する領域の数
         mx(1)=lrg !分割数を変更する各領域の番号
         do m=1,1000000
           nx0=nx
           do i=1,nx
             kx(4*(i-1)+1)=mk(mx(i),1)
             kx(4*(i-1)+2)=mk(mx(i),3)
             kx(4*(i-1)+3)=mk(mx(i),6)
             kx(4*(i-1)+4)=mk(mx(i),5)
           enddo
           do 1200 i=1,4*nx
             if(kx(i).eq.0)goto 1200  !領域番号0の領域と接している場合
             do 1300 j=1,nx
               if(kx(i).eq.mx(j))goto 1200  !既に考慮されている領域の場合
1300         continue
             nx=nx+1
             mx(nx)=kx(i)
1200       continue
           if(nx0.eq.nx)goto 1400     !終了条件:分割数を変更する領域の数nxが増えなかった場合
         enddo
c
1400     do i=1,nx
           mdx(mx(i))=ldx
         enddo
c
c        y軸方向:ny,my(ny),mdy(my(ny))を求める
         ny=1      !分割数を変更する領域の数
         my(1)=lrg !分割数を変更する各領域の番号
         do m=1,1000000
           ny0=ny
           do i=1,ny
             ky(4*(i-1)+1)=mk(my(i),2)
             ky(4*(i-1)+2)=mk(my(i),3)
             ky(4*(i-1)+3)=mk(my(i),4)
             ky(4*(i-1)+4)=mk(my(i),5)
           enddo
           do 1500 i=1,4*ny
             if(ky(i).eq.0)goto 1500  !領域番号0の領域と接している場合
             do 1600 j=1,ny
               if(ky(i).eq.my(j))goto 1500  !既に考慮されている領域の場合
1600         continue
             ny=ny+1
             my(ny)=ky(i)
1500       continue
           if(ny0.eq.ny)goto 1700     !終了条件:分割数を変更する領域の数nxが増えなかった場合
         enddo
c
1700     do i=1,ny
           mdy(my(i))=ldy
         enddo
c
c        z軸方向:nz,mz(nz),mdz(mz(nz))を求める
         nz=1      !分割数を変更する領域の数
         mz(1)=lrg !分割数を変更する各領域の番号
         do m=1,1000000
           nz0=nz
           do i=1,nz
             kz(4*(i-1)+1)=mk(mz(i),2)
             kz(4*(i-1)+2)=mk(mz(i),3)
             kz(4*(i-1)+3)=mk(mz(i),4)
             kz(4*(i-1)+4)=mk(mz(i),5)
           enddo
           do 1800 i=1,4*nz
             if(kz(i).eq.0)goto 1800  !領域番号0の領域と接している場合
             do 1900 j=1,nz
               if(kz(i).eq.mz(j))goto 1800  !既に考慮されている領域の場合
1900         continue
             nz=nz+1
             mz(nz)=kz(i)
1800       continue
           if(nz0.eq.nz)goto 2000     !終了条件:分割数を変更する領域の数nxが増えなかった場合
         enddo
c
2000     do i=1,nz
           mdz(mz(i))=ldz
         enddo
1100   continue
c         
c      各領域の頂点座標を求める
c      nrg,x(mp(nrg,8)),y(mp(nrg,8)),z(mp(nrg,8)) → xp(nrg,8),yp(nrg,8),zp(nrg,8)
       do 2100 krg=1,nrg
         do i=1,8
           xp(krg,i)=x(mp(krg,i))
           yp(krg,i)=y(mp(krg,i))
           zp(krg,i)=z(mp(krg,i))
         continue
2100   continue
c
c      各領域を構成する6面体の頂点座標を求める
c      mdx(nrg),mdy(nrg),mdz(nrg),xp(nrg,8),yp(nrg,8),zp(nrg,8) → 
c      xr(nrg,ndx+1,ndy+1,ndz+1),yr(nrg,ndx+1,ndy+1,ndz+1),zy(nrg,ndx+1,ndy+1,ndz+1)
       do 2200 krg=1,nrg
         ndx=mdx(krg)
         ndy=mdy(krg)
         ndz=mdz(krg)
c        辺1-5,4-8,2-6,3-7上の座標を求める
         do j=0,ndy
           xr(krg,1,j+1,1)=xp(krg,1)+(xp(krg,5)-xp(krg,1))
     &                              *dble(j)/dble(ndy)
           yr(krg,1,j+1,1)=yp(krg,1)+(yp(krg,5)-yp(krg,1))
     &                              *dble(j)/dble(ndy)
           zr(krg,1,j+1,1)=zp(krg,1)+(zp(krg,5)-zp(krg,1))
     &                              *dble(j)/dble(ndy)
           xr(krg,1,j+1,ndz+1)=xp(krg,4)+(xp(krg,8)-xp(krg,4))
     &                                  *dble(j)/dble(ndy)
           yr(krg,1,j+1,ndz+1)=yp(krg,4)+(yp(krg,8)-yp(krg,4))
     &                                  *dble(j)/dble(ndy)
           zr(krg,1,j+1,ndz+1)=zp(krg,4)+(zp(krg,8)-zp(krg,4))
     &                                  *dble(j)/dble(ndy)
           xr(krg,ndx+1,j+1,1)=xp(krg,2)+(xp(krg,6)-xp(krg,2))
     &                                  *dble(j)/dble(ndy)
           yr(krg,ndx+1,j+1,1)=yp(krg,2)+(yp(krg,6)-yp(krg,2))
     &                                  *dble(j)/dble(ndy)
           zr(krg,ndx+1,j+1,1)=zp(krg,2)+(zp(krg,6)-zp(krg,2))
     &                                  *dble(j)/dble(ndy)
           xr(krg,ndx+1,j+1,ndz+1)=xp(krg,3)+(xp(krg,7)-xp(krg,3))
     &                                      *dble(j)/dble(ndy)
           yr(krg,ndx+1,j+1,ndz+1)=yp(krg,3)+(yp(krg,7)-yp(krg,3))
     &                                      *dble(j)/dble(ndy)
           zr(krg,ndx+1,j+1,ndz+1)=zp(krg,3)+(zp(krg,7)-zp(krg,3))
     &                                      *dble(j)/dble(ndy)
         enddo
c        面2,4上の座標を求める
         do j=1,ndy+1
           do k=1,ndz+1
             xr(krg,1,j,k)=xr(krg,1,j,1)
     &                    +(xr(krg,1,j,ndz+1)-xr(krg,1,j,1))
     &                    *dble(k-1)/dble(ndz)
             yr(krg,1,j,k)=yr(krg,1,j,1)
     &                    +(yr(krg,1,j,ndz+1)-yr(krg,1,j,1))
     &                    *dble(k-1)/dble(ndz)
             zr(krg,1,j,k)=zr(krg,1,j,1)
     &                    +(zr(krg,1,j,ndz+1)-zr(krg,1,j,1))
     &                    *dble(k-1)/dble(ndz)
             xr(krg,ndx+1,j,k)=xr(krg,ndx+1,j,1)
     &                    +(xr(krg,ndx+1,j,ndz+1)-xr(krg,ndx+1,j,1))
     &                    *dble(k-1)/dble(ndz)
             yr(krg,ndx+1,j,k)=yr(krg,ndx+1,j,1)
     &                    +(yr(krg,ndx+1,j,ndz+1)-yr(krg,ndx+1,j,1))
     &                    *dble(k-1)/dble(ndz)
              zr(krg,ndx+1,j,k)=zr(krg,ndx+1,j,1)
     &                    +(zr(krg,ndx+1,j,ndz+1)-zr(krg,ndx+1,j,1))
     &                    *dble(k-1)/dble(ndz)
           enddo
         enddo
c        全体の座標を求める
         do i=1,ndx+1
           do j=1,ndy+1
             do k=1,ndz+1
               xr(krg,i,j,k)=xr(krg,1,j,k)
     &                      +(xr(krg,ndx+1,j,k)-xr(krg,1,j,k))
     &                      *dble(i-1)/dble(ndx)
               yr(krg,i,j,k)=yr(krg,1,j,k)
     &                      +(yr(krg,ndx+1,j,k)-yr(krg,1,j,k))
     &                      *dble(i-1)/dble(ndx)
               zr(krg,i,j,k)=zr(krg,1,j,k)
     &                      +(zr(krg,ndx+1,j,k)-zr(krg,1,j,k))
     &                      *dble(i-1)/dble(ndz)
             enddo
           enddo
         enddo
2200   continue
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
             if(mk(krg,1).lt.krg .nad. kdy.eq.1)goto 2500           !面上の節点
             do 2600 kdz=1,mdz(krg)+1     !ndz=mdz(nrg)
               if(mk(krg,3).lt.krg .and. kdz.eq.mdz(krg)+1)goto 2600  !面上の節点
               if(mk(krg,5).lt.krg .and. kdz.eq.1)goto 2600           !面上の節点
               np=np+1
               nprg(krg,kdx,kdy,kdz)=np
               mprg(np)=krg
               npx(krg,np)=kdx
               npy(krg,np)=kdy
               npz(krg,np)=kdz
c
               xx(np)=xr(krg,kdx,kdy,kdz)
               yy(np)=yr(krg,kdx,kdy,kdz)
               zz(np)=zr(krg,kdx,kdy,kdz)
2600         continue
2500       continue
2400     continue
2300   continue 
c
c      6面体を4面体に分割する
c
       nele=0
       do 2700 krg=1,nrg
         do 2800 kdx=1,mx(krg)
           do 2900 kdy=1,my(krg)
             do 3000 kdz=1,mz(krg)
               nele=nele+1
               ne(nele,1)=nprg(krg,kdx,kdy,kdz)
               ne(nele,2)=nprg(krg,kdx+1,kdy,kdz)
               ne(nele,3)=nprg(krg,kdx+1,kdy,kdz+1)
               ne(nele,4)=nprg(krg,kdx,kdy,kdz+1)
               ne(nele,5)=nprg(krg,kdx,kdy+1,kdz)
               ne(nele,6)=nprg(krg,kdx+1,kdy+1,kdz)
               ne(nele,7)=nprg(krg,kdx+1,kdy+1,kdz+1)
               ne(nele,8)=nprg(krg,kdx,kdy+1,kdz+1)
               mdx(nele)=mr(krg)
3000         continue
2900       continue
2800     continue
2700   continue
c
c      4面体に分割する
c
       nelew=0
       do 3100 kele=1,nele     !6面体
         new(nelew,1)=ne(kele,6)
         new(nelew,2)=ne(kele,5)
         new(nelew,3)=ne(kele,8)
         new(nelew,4)=ne(kele,1)
         idx(nelew)=mdx(kele)
c
         nelew=nelew+1
         new(nelew,1)=ne(kele,4)
         new(nelew,2)=ne(kele,1)
         new(nelew,3)=ne(kele,2)
         new(nelew,4)=ne(kele,6)
         idx(nelew)=mdx(kele)
c
         nelew=nelew+1
         new(nelew,1)=ne(kele,8)
         new(nelew,2)=ne(kele,1)
         new(nelew,3)=ne(kele,4)
         new(nelew,4)=ne(kele,6)
         idx(nelew)=mdx(kele)
c
         nelew=nelew+1
         new(nelew,1)=ne(kele,7)
         new(nelew,2)=ne(kele,6)
         new(nelew,3)=ne(kele,8)
         new(nelew,4)=ne(kele,3)
         idx(nelew)=mdx(kele)
c
         nelew=nelew+1
         new(nelew,1)=ne(kele,3)
         new(nelew,2)=ne(kele,4)
         new(nelew,3)=ne(kele,2)
         new(nelew,4)=ne(kele,6)
         idx(nelew)=mdx(kele)
c
         nelew=nelew+1
         new(nelew,1)=ne(kele,3)
         new(nelew,2)=ne(kele,8)
         new(nelew,3)=ne(kele,4)
         new(nelew,4)=ne(kele,6)
         idx(nelew)=mdx(kele)         
3200   continue
c
c      バンド幅を最適化
c
       call renum3(
c
c      *.1および*,2ファイルの作成
c
       call addchr(gfile1,rfile,'.1')
       open(1,file=grile1)
       rewind(1)
       write(1,*)np
       do 3300 i=1,np
3300   write(1,*)i,jww(i),xx(i),yy(i),zz(i)
       close(1)
c
       call addchd(gfile2,rfile,'.2')
       open(2,file=gfile2)
       rewind(1)
       write(1,*)nelew,np,nbw
       do 3400 i=1,nelew
3400   write(1,*)i,(jww(new(i,j)),j=1,4)
     &             ,(xx(new(i,j)),yy(new(i,j)),zz(new(i,j)),j=1,4)
3500 &             ,idx(nelew)
       close(1)
c
       stop
       end