c******************************************************************
c
c       サブルーチン内の値の保持
c       →サブルーチンが呼ばれた後，サブルーチン内の値は保持されている．
c******************************************************************

        implicit double precision(a-h,o-z)
c
        do 10 i=1,100
10      call value(i)
c
        stop
        end
c
        subroutine value(i)
        implicit double precision(a-h,o-z)
        if(i.eq.1)a=1d0
        if(i.ne.1)a=a+1d0
        write(*,*)'i,a',i,a
c
        return
        end