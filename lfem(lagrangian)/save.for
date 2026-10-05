c-------------------------------------------------------------------------------------
c
c       流動解析の結果を保存する
c
c-------------------------------------------------------------------------------------
        subroutine svdata(np,nele,ne,xx,yy,zz,uu,vv,ww,pp,pnx,pny,pnz
     &                                         ,idx,ti,rfile,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),xx(ine),yy(ind),zz(ind),idx(0:ine)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        dimension ut(ind),vt(ind),wt(ind),pt(ind)
        character lpp*12,rfile*40,vfile*40
c------記録用--------------------------------------------------------
       dimension pnx(ind),pny(ind),pnz(ind)
       dimension ptx(ind),pty(ind),ptz(ind)
c--------------------------------------------------------------------
c
        do i=1,nele
          if(idx(i).eq.1 .or. idx(i).eq.2 .or. 
     &       idx(i).eq.3 .or. idx(i).eq.4)then     !要素iが流動領域の場合
            do j=1,4
              ut(ne(i,j))=uu(ne(i,j))
              vt(ne(i,j))=vv(ne(i,j))
              wt(ne(i,j))=ww(ne(i,j))
              pt(ne(i,j))=pp(ne(i,j))
c
              ptx(ne(i,j))=pnx(ne(i,j))
              pty(ne(i,j))=pny(ne(i,j))
              ptz(ne(i,j))=pnz(ne(i,j))
            enddo
          else    !要素が流動領域ではない場合
            do j=1,4
              ut(ne(i,j))=0d0
              vt(ne(i,j))=0d0
              wt(ne(i,j))=0d0
              pt(ne(i,j))=0d0
c
              ptx(ne(i,j))=0d0
              pty(ne(i,j))=0d0
              ptz(ne(i,j))=0d0
            enddo
          endif
        enddo
c
c        write(lpp,'(f12.1)') ti*1d3  !input a value into lpp
c        write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
        write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
        do 100 i=1,12
100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,lpp//'.inp',vfile)
c        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
c        write(1,*)"# flow data" !コメント行
        write(1,*)1             !ステップ数
        write(1,*)"data"        !データの繰り返しタイプ
        write(1,*)"step1 ",ti   !ステップ番号,コメント
        write(1,*)np,nele       !節点数,要素数
        do 200 i=1,np
200     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))    !節点番号,x,y,z座標
        do 300 i=1,nele
300     write(1,*)i,1,"tet",(ne(i,j),j=1,4)    !要素番号,0,tet,ne(i,1-4)
        write(1,*)"7 0"          !節点のデータ数,要素のデータ数
        write(1,*)"7 1 1 1 1 1 1 1"    !節点のデータ成分数,1(構成数,スカラーのみなので1),1,1,1
        write(1,*)"u,[-]"    !節点データ成分1のラベル,単位
        write(1,*)"v,[-]"
        write(1,*)"w,[-]"
        write(1,*)"p,[-]"
        write(1,*)"sgx,[-]"    !表面張力
        write(1,*)"sgy,[-]"
        write(1,*)"sgz,[-]"
        do 400 i=1,np
400     write(1,*)i,sngl(ut(i)),sngl(vt(i)),sngl(wt(i)),sngl(pt(i))
     &              ,sngl(ptx(i)),sngl(pty(i)),sngl(ptz(i))   !節点番号,節点データ1,2,3,4
        close(1)
c
        return
        end
c-------------------------------------------------------------------------------------
c
c       物質移動解析の結果を保存する
c
c-------------------------------------------------------------------------------------
        subroutine svdatam(np,nele,ne,xx,yy,zz,cc,idx,ti,rfile,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),xx(ine),yy(ind),zz(ind),idx(0:ine)
        dimension cc(ind),ct(ind)
        character lpp*12,rfile*40,vfile*40
c
        do 10 i=1,nele
          if(idx(i).eq.2 .or. idx(i).eq.4 .or. 
     &       idx(i).eq.5 .or. idx(i).eq.6)then     !要素iが物質移動領域の場合
            do 20 j=1,4
20          ct(ne(i,j))=cc(ne(i,j))
          else    !要素が物質移動領域ではない場合
            do 30 j=1,4
30          ct(ne(i,j))=0d0
          endif
10      continue
c
c        write(lpp,'(f12.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
        do 100 i=1,12
100     if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'m'//lpp//'.inp',vfile)
c        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
        write(1,*)np,nele,1,0,0
        do 200 i=1,np
200     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
        do 300 i=1,nele
300     write(1,*)i,1,'tet',(ne(i,j),j=1,4)
        write(1,*)1,1    !Ca
        write(1,*)'cont','[-]'
        do 400 i=1,np
400     write(1,*)i,sngl(ct(i))
        close(1)
c
        return
        end
c-------------------------------------------------------------------------------------
c
c       熱移動解析の結果を保存する
c
c-------------------------------------------------------------------------------------
        subroutine svdatat(np,nele,ne,xx,yy,zz,tt,idx,ti,rfile,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),xx(ine),yy(ind),zz(ind),idx(0:ine)
        dimension tt(ind),tm(ind)
        character lpp*12,rfile*40,vfile*40
c
        do 10 i=1,nele
          if(idx(i).eq.3 .or. idx(i).eq.4 .or. 
     &       idx(i).eq.5 .or. idx(i).eq.7)then     !要素iが物質移動領域の場合
            do 20 j=1,4
20          tm(ne(i,j))=tt(ne(i,j))
          else    !要素が物質移動領域ではない場合
            do 30 j=1,4
30          tm(ne(i,j))=0d0
          endif
10      continue
c
c        write(lpp,'(f12.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
        do 100 i=1,12
100     if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'t'//lpp//'.inp',vfile)
c        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
        write(1,*)np,nele,1,0,0
        do 200 i=1,np
200     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
        do 300 i=1,nele
300     write(1,*)i,1,'tet',(ne(i,j),j=1,4)
        write(1,*)1,1    !Ca
        write(1,*)'temp','[-]'
        do 400 i=1,np
400     write(1,*)i,sngl(tm(i))
        close(1)
c
        return
        end