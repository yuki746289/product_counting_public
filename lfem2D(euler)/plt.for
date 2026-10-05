c-----------------------------------------------------------
c
c       —v‘f‚ð•`‚­
c
c-----------------------------------------------------------
        subroutine pltne(np,nele,ne,xx,yy)
        include "head.for"
        dimension ne(ine,3)
        dimension xx(ind),yy(ind)
c       
        open(1,file="pne.res")
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
c-----------------------------------------------------------
c
c       ‹«ŠE‚ð•`‚­
c
c-----------------------------------------------------------
        subroutine pltnb(na,nb,nelb,neb,ncom)
        include "head.for"
        dimension nb(ica)
        dimension nelb(ica,-1:ibd)
        dimension neb(ica,-1:ibd)
        dimension ncom(ica,-1:ibd)
c       
        open(1,file="pnb.res")
        rewind(1)
        write(1,*)na
        do ka=1,na
          write(1,*)nb(ka)
          do kb=1,nb(ka)
            write(1,*)ka,kb,nelb(ka,kb),neb(ka,kb),ncom(ka,kb)
          enddo
        enddo
        close(1)
c
        return
        end
