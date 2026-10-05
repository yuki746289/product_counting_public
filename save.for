c-------------------------------------------------------------------------------------
c
c       流動解析の結果を保存する
c
c-------------------------------------------------------------------------------------
        subroutine svdata(np,nele,ne,xx,yy,zz,uu,vv,ww,pp,idx,
     &                                         ti,rfile,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),xx(ine),yy(ind),zz(ind),idx(0:ine)
        dimension uu(ind),vv(ind),ww(ind),pp(ind)
        dimension ut(ind),vt(ind),wt(ind),pt(ind)
        character lpp*8,rfile*20,vfile*20
c
        do i=1,nele
          if(idx(i).eq.1 .or. idx(i).eq.2 .or. 
     &       idx(i).eq.3 .or. idx(i).eq.4)then     !要素iが流動領域の場合
            do j=1,4
              ut(ne(i,j))=uu(ne(i,j))
              vt(ne(i,j))=vv(ne(i,j))
              wt(ne(i,j))=ww(ne(i,j))
              pt(ne(i,j))=pp(ne(i,j))
            enddo
          else    !要素が流動領域ではない場合
            do j=1,4
              ut(ne(i,j))=0d0
              vt(ne(i,j))=0d0
              wt(ne(i,j))=0d0
              pt(ne(i,j))=0d0
            enddo
          endif
        enddo
c
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
100    if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,lpp//'.inp',vfile)
c        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
        do 200 i=1,nele
200     write(1,*)i,idx(i),(ne(i,j),j=1,4)
        do 300 i=1,np
300     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
        do 400 i=1,np
400     write(1,*)i,sngl(ut(i)),sngl(vt(i)),sngl(wt(i)),sngl(pt(i))
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
        character lpp*8,rfile*20,vfile*20
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
c        write(lpp,'(f8.1)') ti*1d3  !input a value into lpp
        write(lpp,'(f8.1)') ti*1d4  !input a value into lpp
        do 100 i=1,8
100     if(lpp(i:i).eq.' ') lpp(i:i)='0'
        call addchr(rfile,'m'//lpp//'.inp',vfile)
c        write(*,*) 'avsuvp/ti,file=',ti,vfile
c
        open(1,file=vfile)
        rewind(1)
        do 200 i=1,nele
200     write(1,*)i,idx(i),(ne(i,j),j=1,4)
        do 300 i=1,np
300     write(1,*)i,sngl(xx(i)),sngl(yy(i)),sngl(zz(i))
        do 400 i=1,np
400     write(1,*)i,sngl(ct(i))
c
        return
        end