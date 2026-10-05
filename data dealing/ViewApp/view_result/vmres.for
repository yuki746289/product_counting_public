c********************************************************************
c
c       VRML(virtual reality modeling langage)によるgrid fileの表示
c
c********************************************************************
        include 'addchr.for'
        implicit double precision(a-h,o-z)
        parameter(ind=19999,ine=99999,ibw=1999)
c-------inp file-------------------------------------------------------------        
        dimension ne(ine,4),xx(ind),yy(ind),zz(ind)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        character*30 rfile,vfile,cnone
c-------grid file------------------------------------------------------------
        dimension mr(ine)                           !領域のインデクス番号
        dimension ifline(ine,12)                    !0:直線,1:曲線
        dimension mdx(1001),mdy(1001),mdz(1001)     !各領域の分割数
        dimension mk(1001,6)                        !領域の各面を共有している領域の番号
        dimension mp(1001,8)                        !領域の各頂点の節点番号
        dimension xp(1001,12,3),yp(1001,12,3),zp(1001,12,3)  !辺を構成する節点の座標
        dimension xv(1001,8),yv(1001,8),zv(1001,8)  !各6面体の頂点座標
        character*32 pfile                          !file names
c-----------------------------------------------------------------------------
        data rt/3d-1/        !rt:スケール
c
c------------------------------------------------------------------------
c
c       inp fileの読み込み
c
c------------------------------------------------------------------------
       open(1,file='input.txt')
       rewind(1)
       read(1,*)rfile
       close(1)
       call addchr(rfile,'.inp',vfile)
       close(1)
       write(*,*)'readed',vfile
c
       open(1,file=vfile)
       rewind(1)
       write(1,*)np,nel,none,none,none
       do 10 i=1,np
10     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
       do 20 i=1,nel
20     write(1,*)i,1,cnone,(ne(i,j),j=1,4)
       write(1,*)none,none,none
       write(1,*)cnone,cnone
       write(1,*)cnone,cnone
       do 30 i=1,np
30     write(1,*)i,sngl(pp(i)),sngl(uu(i)),sngl(vv(i)),sngl(ww(i))
       close(1)
c------------------------------------------------------------------------
c
c       流れ場の表示
c
c------------------------------------------------------------------------
        do 800 mm=1,2
          if(mm.eq.1)then
            open(1,file='vgrid.inf')
          elseif(mm.eq.2)then
            open(1,file='vgrid.wrl')
          endif
c
          rewind(1)
          write(1,*)'#VRML V2.0 utf8'
          write(1,*)'DEF g Transform {'
          write(1,*)'translation 0 -2.0 0.0'
          write(1,*)'children ['
          write(1,*)'DEF tree Shape {'
          write(1,*)'appearance Appearance {'
          write(1,*)'material Material {'
          write(1,*)'diffuseColor 1 1 1 '
          write(1,*)'}'
          write(1,*)'}'
          write(1,*)'geometry IndexedLineSet {'
          write(1,*)'coord Coordinate {'
          write(1,*)'point ['
          do 400 i=1,np
            write(1,*)sngl(xx(i)),sngl(yy(i)),sngl(zz(i)),','
            write(1,*)sngl(xx(i)+uu(i)*rt),
     &                 sngl(yy(i)+vv(i)*rt),sngl(zz(i)+ww(i)*rt)
            if(i.ne.np)write(1,*)','
400       continue
          write(1,*)']'
          write(1,*)'}'
          write(1,*)'coordIndex [ '
          do 500 i=1,np
500       write(1,*)(2*i-1)-1,2*i-1,'-1'
          write(1,*)']'
          write(1,*)'	color Color {'
          write(1,*)'color [ '
          do 600 i=1,np
            vab=dsqrt(uu(i)**2+vv(i)**2+ww(i)**2)   !an absolute value:絶対値
            write(1,*)vab,'0.5 0.5'
            if(i.ne.np)write(1,*)','
600       continue
          write(1,*)']'
          write(1,*)'}'
          write(1,*)'colorIndex ['
          do 700 i=1,np
700       write(1,*)i-1
          write(1,*)']'
          write(1,*)'	colorPerVertex FALSE'
          write(1,*)'}'
          write(1,*)'}'
          write(1,*)']'
          write(1,*)'}'
c
c          close(1)
c
          if(mm.eq.1)write(*,*)'created vgrid.txt'
          if(mm.eq.2)write(*,*)'created vgrid.wrl'
c***
          goto 900
c***
800     continue
c
900     continue
c------------------------------------------------------------------------
c
c       grid fileの読み込み
c
c------------------------------------------------------------------------
       open(2,file='input2.txt')
       rewind(2)
       read(2,*)rfile
       close(2)
       call addchr(rfile,'.inf',pfile)
       write(*,*)'the input file'
       write(*,*)pfile
c
       open(2,file=pfile)
       rewind(2)
       read(2,*)npr    !頂点の数
       read(2,*)nrg    !領域数
       do 1000 krg=1,nrg
         read(2,*)mr(krg)  !領域krgのインデックス番号
         read(2,*)mdx(krg),mdy(krg),mdz(krg)           !各軸方向の分割数
         do 1100 il=1,12
           read(2,*)ifline(krg,il)                     !辺の種類(0:直線,1:曲線)
           nside=2                             !直線
           if(ifline(krg,il).eq.1)nside=3      !曲線
           read(2,*)(xp(krg,il,k),k=1,nside)  !各頂点のx座標
           read(2,*)(yp(krg,il,k),k=1,nside)  !各頂点のy座標
           read(2,*)(zp(krg,il,k),k=1,nside)  !各頂点のz座標
c--------------------------------------------------------------------------------
c           write(*,*)'krg il',krg,il
c           write(*,*)ifline(krg,il)            !辺の種類(0:直線,1:曲線)
c           write(*,*)(xp(krg,il,k),k=1,nside)  !各頂点のx座標
c           write(*,*)(yp(krg,il,k),k=1,nside)  !各頂点のy座標
c           write(*,*)(zp(krg,il,k),k=1,nside)  !各頂点のz座標
c           pause
c--------------------------------------------------------------------------------
1100      continue
1000    continue
       do 1200 krg=1,nrg
1200    read(2,*)(mk(krg,i),i=1,6)  !各面が接する領域の番号
       do 1300 krg=1,nrg
1300    read(2,*)(mp(krg,i),i=1,8)  !各頂点の節点番号
       close(2)
c
c      各6面体の頂点座標配列xv(nrg,8),yv(nrg,8),zv(nrg,8)を求める
c
       do 1400 krg=1,nrg
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
1400    continue
c------------------------------------------------------------------------
c
c       解析領域の表示
c
c------------------------------------------------------------------------
          write(1,*)'DEF g2 Transform {'
          write(1,*)'translation 0 -2.0 0.0'
          write(1,*)'children ['
          write(1,*)'DEF tree2 Shape {'
          write(1,*)'appearance Appearance {'
          write(1,*)'material Material {'
          write(1,*)'diffuseColor 1 1 1 '
          write(1,*)'}'
          write(1,*)'}'
          write(1,*)'geometry IndexedLineSet {'
          write(1,*)'coord Coordinate {'
          write(1,*)'point ['
          do 1500 krg=1,nrg
          do 1500 i=1,8
          write(1,*)sngl(xv(krg,i)),sngl(yv(krg,i)),sngl(zv(krg,i))
          if(krg.eq.nrg .and. i.eq.8)goto 1600
          write(1,*)','
1500      continue
1600      continue
          write(1,*)']'
          write(1,*)'}'
          write(1,*)'coordIndex [ '
          do 1700 krg=1,nrg
          write(1,*)8*(krg-1)+2-1,8*(krg-1)+3-1,8*(krg-1)+7-1,    !面:2
     &               8*(krg-1)+6-1,8*(krg-1)+2-1,'-1'
          write(1,*)8*(krg-1)+4-1,8*(krg-1)+3-1,8*(krg-1)+7-1,    !面:3
     &               8*(krg-1)+8-1,8*(krg-1)+4-1,'-1'
          write(1,*)8*(krg-1)+1-1,8*(krg-1)+4-1,8*(krg-1)+8-1,    !面:4
     &               8*(krg-1)+5-1,8*(krg-1)+1-1,'-1'
          write(1,*)8*(krg-1)+1-1,8*(krg-1)+2-1,8*(krg-1)+6-1,    !面:5
     &               8*(krg-1)+5-1,8*(krg-1)+1-1,'-1'
1700      continue
          write(1,*)']'
          write(1,*)'	color Color {'
          write(1,*)'color [ 1.0 0.0 0.0 ]'
          write(1,*)'}'
          write(1,*)'colorIndex ['
          do 1800 krg=1,nrg
1800       write(1,*)' 0 0 0 0'
          write(1,*)']'
          write(1,*)'	colorPerVertex FALSE'
          write(1,*)'}'
          write(1,*)'}'
          write(1,*)']'
          write(1,*)'}'
c
        close(1)
        stop
        end