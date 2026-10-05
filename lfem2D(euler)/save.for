c-------------------------------------------------------------------
c
c      avs用のデータを保存する
c
c-------------------------------------------------------------------
       subroutine svdata(np,nele,ne,xx,yy,uu,vv,pp,idx,ti)
       include "head.for"
       dimension ne(ine,3)
       dimension xx(ind),yy(ind)
       dimension uu(ind),vv(ind),pp(ind)
       dimension idx(0:ine)
c
       character*40 tfile
       character*12 lpp
c
       write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
       do i=1,12
         if(lpp(i:i).eq.' ') lpp(i:i)='0'
	 enddo
       call addchr("./inp/data",lpp//'.inp',tfile)
       open(1,file=tfile)
       rewind(1)
       write(1,*)"# flow data" !コメント行
       write(1,*)1             !ステップ数
       write(1,*)"data"        !データの繰り返しタイプ
       write(1,*)"step1 ",ti   !ステップ番号,コメント
       write(1,*)np,nele       !節点数,要素数
       do 200 i=1,np
200    write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(0d0)    !節点番号,x,y,z座標
       do 300 i=1,nele
300    write(1,*)i,1,"tri",(ne(i,j),j=1,3)    !要素番号,0,tri,ne(i,1-3)
       write(1,*)"3 0"     !節点のデータ数,要素のデータ数
       write(1,*)"3 1 1 1" !節点のデータ成分数,1(構成数,スカラーのみなので1),1,1,1
       write(1,*)"u,[-]"   !節点データ成分1のラベル,単位
       write(1,*)"v,[-]"
       write(1,*)"p,[-]"
       do 400 i=1,np
400    write(1,*)i,sngl(uu(i)),sngl(vv(i)),sngl(pp(i))    !節点番号,節点データ1,2,3
       close(1)
c
       return
       end
