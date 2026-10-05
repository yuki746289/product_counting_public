c********************************************************************
c
c       VRML(virtual reality modeling langage)によるgrid fileの表示
c
c********************************************************************
        include 'addchr.for'
        implicit double precision(a-h,o-z)
        parameter(ind=19999,ine=99999,ibw=1999)
        dimension ne(ine,4),xx(ind),yy(ind),zz(ind),idx(0:ine)
        character*30 rfile,gfile
c------------------------------------------------------------------------
c
c       grid fileの読み込み
c
c------------------------------------------------------------------------
       open(1,file='input.txt')
       rewind(1)
       read(1,*)rfile
       close(1)
       call addchr(rfile,'.2',gfile)
       close(1)
       write(*,*)'readed',gfile
c
       open(1,file=gfile)
       rewind(1)
       read(1,*)nel,np,nbw          !全要素数,全節点数,バンド幅
       do 100 i=1,nel
         read(1,*)kel,(ne(kel,j),j=1,4)
         read(1,*)(xx(ne(kel,j)),yy(ne(kel,j)),zz(ne(kel,j)),j=1,4)
         read(1,*)idx(kel)
100    enddo
       close(1)
c------------------------------------------------------------------------
c
c       VRML codeの作成
c
c------------------------------------------------------------------------
        do 500 mm=1,2
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
          do 200 i=1,np
            write(1,*)sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
            if(i.ne.np)write(1,*)','
200       continue
          write(1,*)']'
          write(1,*)'}'
          write(1,*)'coordIndex [ '
          do 300 i=1,nel
            write(1,*)ne(i,1)-1,ne(i,2)-1,ne(i,3)-1,ne(i,1)-1,'-1'
            write(1,*)ne(i,2)-1,ne(i,3)-1,ne(i,4)-1,ne(i,2)-1,'-1'
            write(1,*)ne(i,3)-1,ne(i,4)-1,ne(i,1)-1,ne(i,3)-1,'-1'
            write(1,*)ne(i,4)-1,ne(i,1)-1,ne(i,2)-1,ne(i,4)-1,'-1'
300       continue

          write(1,*)']'
          write(1,*)'	color Color {'
          write(1,*)'color [ 1.0 1.0 1.0, 1.0 0.0 0.0, 0.0 1.0 0.0, 
     &                        0.0 0.0 1.0, 0.75 0.75 0.75 ]'    !0:白,1:赤,2:緑,3:青,4:茶
          write(1,*)'}'
          write(1,*)'colorIndex ['
          do 400 i=1,nel
            if(idx(i).eq.1)then
              write(1,*)'0 0 0 0'
              write(1,*)'0 0 0 0'
              write(1,*)'0 0 0 0'
              write(1,*)'0 0 0 0'
            elseif(idx(i).eq.2)then
              write(1,*)'1 1 1 1'
              write(1,*)'1 1 1 1'
              write(1,*)'1 1 1 1'
              write(1,*)'1 1 1 1'
            elseif(idx(i).eq.3)then
              write(1,*)'2 2 2 2'
              write(1,*)'2 2 2 2'
              write(1,*)'2 2 2 2'
              write(1,*)'2 2 2 2'
            elseif(idx(i).eq.4)then
              write(1,*)'3 3 3 3'
              write(1,*)'3 3 3 3'
              write(1,*)'3 3 3 3'
              write(1,*)'3 3 3 3'
            else
              write(1,*)'4 4 4 4'
              write(1,*)'4 4 4 4'
              write(1,*)'4 4 4 4'
              write(1,*)'4 4 4 4'
            endif
400       continue
          write(1,*)']'
          write(1,*)'	colorPerVertex FALSE'
          write(1,*)'}'
          write(1,*)'}'
          write(1,*)']'
          write(1,*)'}'
c
          close(1)
c
          if(mm.eq.1)write(*,*)'created vgrid.txt'
          if(mm.eq.2)write(*,*)'created vgrid.wrl'
c***
          goto 600
c***
500     continue
c
600     continue
        stop
        end