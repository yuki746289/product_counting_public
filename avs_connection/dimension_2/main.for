       implicit double precision(a-h,o-z)
       parameter(ind=4999,ine=4999)
       dimension ne(ine,3)
       dimension xx(ind),yy(ind),zz(ind)
       dimension uu(ind),vv(ind),pp(ind)
c
       character*40 tfile
       character*12 lpp
c
       character*499 dummy
c
       data num/30/    !読み込むデータ数
       data ti/1d-7/
c
       open(10,file="./data_2d/data.inp")
       rewind(10)
c
       do n=1,num
c
c        読み込み
         write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
         do i=1,12
           if(lpp(i:i).eq.' ') lpp(i:i)='0'
         enddo
         call addchr("./data_2d/data",lpp//'.inp',tfile)
         open(1,file=tfile)
         rewind(1)
         read(1,*)dummy    !"# flow data" !コメント行
         read(1,*)dummy    !1             !ステップ数
         read(1,*)dummy    !"data"        !データの繰り返しタイプ
         read(1,*)dummy    !"step1 ",ti   !ステップ番号,コメント
         read(1,*)np,nele  !節点数,要素数
         do 200 i=1,np
200      read(1,*)i,xx(i),yy(i),zz(i)    !節点番号,x,y,z座標
         do 300 i=1,nele
300      read(1,*)i,dummy,dummy,(ne(i,j),j=1,3)    !要素番号,1,tri,ne(i,1-3)
         read(1,*)dummy    !"3 0"     !節点のデータ数,要素のデータ数
         read(1,*)dummy    !"3 1 1 1" !節点のデータ成分数,1(構成数,スカラーのみなので1),1,1,1
         read(1,*)dummy    !"u,[-]"   !節点データ成分1のラベル,単位
         read(1,*)dummy    !"v,[-]"
         read(1,*)dummy    !"p,[-]"
         do 400 i=1,np
400      read(1,*)i,uu(i),vv(i),pp(i)    !節点番号,節点データ1,2,3
         close(1)
c
c        書き込み
         if(n.eq.1)then
           write(10,*)"# flow data" !コメント行
           write(10,*)1             !ステップ数
           write(10,*)"data"        !データの繰り返しタイプ
           write(10,*)"step1 ",n   !ステップ番号,コメント
           write(10,*)np,nele       !節点数,要素数
           do 500 i=1,np
500        write(10,*)i,sngl(xx(i)),sngl(yy(i)),sngl(0d0)    !節点番号,x,y,z座標
           do 600 i=1,nele
600        write(10,*)i,1,"tri",(ne(i,j),j=1,3)    !要素番号,0,tri,ne(i,1-3)
         else
           write(10,*)"step",n," ",n   !ステップ番号,コメント
         endif
c
         write(10,*)"3 0"     !節点のデータ数,要素のデータ数
         write(10,*)"3 1 1 1" !節点のデータ成分数,1(構成数,スカラーのみなので1),1,1,1
         write(10,*)"u,[-]"   !節点データ成分1のラベル,単位
         write(10,*)"v,[-]"
         write(10,*)"p,[-]"
         do 700 i=1,np
700      write(10,*)i,sngl(uu(i)),sngl(vv(i)),sngl(pp(i))    !節点番号,節点データ1,2,3
       enddo
c
       close(10)
c
       stop
       end
