	
	implicit double precision (a-h,o-z)
      	parameter(ind=9001,ine=6001)
      	dimension xx(ind),yy(ind),xi(ind),yi(ind)
	dimension ne(ine,3),id(ine),zz(6)
	character*30 rfile,lfile,sfile
c
c---------------------------------- MS/CQ Fortran only
	open(1,file='intplot.txt')
	read(1,'(a)') rfile
        read(1,'(a)') lfile
        read(1,'(a)') sfile
	close(1)

	open(11,file=rfile)
	open(12,file=lfile)
        open(13,file=sfile)
c----------------------------------

      read(11,*) np
      do 10 i=1,np
  10  read(11,*) j,k,xx(k),yy(k)
      close(11)

      read(12,*) nel,np,nbw
      do 20 i=1,nel
 20   read(12,*) k,(ne(k,j),j=1,3),zz,id(i)
      close(12)

      read(13,*) intot
      do 30 i=1,intot
 30   read(13,*) xi(i),yi(i)
c-------------------------------------
	write(*,*) 'np,nel,nbw'
        write(*,*) np,nel,nbw

	call plt(nel,ne,xx,yy,xi,yi,intot,id)
	stop
	end

	subroutine plt(nel,ne,xx,yy,xi,yi,intot,id)
      	parameter(ind=9001,ine=6001)
      	implicit double precision (a-h,o-z)
      	dimension ne(ine,3),id(ine)
      	real*8 xx(ind),yy(ind),xi(3),yi(3)
c

	open(8,file="pic.inf")
      	rewind(8)
      	read(8,*) xam,yam,xo,yo,xul,xum,yul,yum

      	call xint
      	call xviewp(-1,xl,xu,yl,yu)
      	call chrsize(3d0,-2d0,0d0)
      	call plot(xo,yo,-3)
      	call usrcord(xam,yam,xul,xum,yul,yum)

	do 10 i=1,nel
	if(id(i).eq.2) call newpen(2)
	if(id(i).ne.2) call newpen(1)
	do 20 j=1,3
	 jj=j+1
	 if(jj.eq.4) jj=1
        k=ne(i,j)
	l=ne(i,jj)
        call plot(xnml(xx(k)),ynml(yy(k)),3)
	call plot(xnml(xx(l)),ynml(yy(l)),2)
 20     enddo
 10	enddo


      	call newpen(3)

	do i=1,intot-1
	call plot(xnml(xi(i)),ynml(yi(i)),3)
	call plot(xnml(xi(i+1)),ynml(yi(i+1)),2)
	enddo

      	call doplot
      	call xend
      	return
      	end
c
