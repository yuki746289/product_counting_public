c-----------------------------------------------------------
c
c       —v‘f‚ð•`‚­
c
c-----------------------------------------------------------
        subroutine plt_mesh(i,np,nele,ne,xx,yy,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,3)
        dimension xx(ind),yy(ind)
        character*32 rfile
        character*32 fname        
c
        rfile=fname(i)
c       
        open(1,file=rfile)
        rewind(1)
        write(1,*)np,nele
        do 100 kele=1,nele
100       write(1,*)kele,ne(kele,1),ne(kele,2),ne(kele,3)
c
        do 200 kp=1,np
200       write(1,*)kp,sngl(xx(kp)),sngl(yy(kp))
        close(1)
c
        return
        end
