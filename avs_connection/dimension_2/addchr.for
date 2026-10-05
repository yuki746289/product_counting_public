        subroutine addchr(a,b,f)
        character*(*) a,b,f
        character*1 s
        s=' '
        ia=index(a,s)-1
        ib=index(b,s)-1
        if(ia+1.eq.0) ia=len(a)
        if(ib+1.eq.0) ib=len(b)
        f=a(:ia)//b(:ib)
        if(ia+ib.gt.len(f)) goto 10
        return
 10     write(*,*) 
        write(*,*) '## error in subroutine fchr ##'
        write(*,*) '   set func.=',f
        write(*,*)
        return
        end
