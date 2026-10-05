c-----------------------------------------------------------------------------------
c
c      KEX‚ĚŹÁ‹Ž–@
c
c-----------------------------------------------------------------------------------
        subroutine gauss(n,nbw,a,x)
        include "head.for"
        dimension a(ind,ibw),x(ind)
c
        do 2 i=1,n
        if(a(i,nbw).le.0d0) then
         write(*,*) '[gauss]',i,a(i,nbw)
         stop
        endif
 2      continue
c
        nbwm=nbw-1
        do 10 i=1,n
          aa=a(i,nbw)
          x(i)=x(i)/aa
          do 20 j=i,min(n,i+nbwm)
 20       a(i,nbw+j-i)=a(i,nbw+j-i)/aa
          do 40 k=i+1,min(n,i+nbwm)
            cc=a(k,nbw+i-k)
            do 30 j=i+1,min(n,i+nbwm)
 30         a(k,nbw+j-k)=a(k,nbw+j-k)-cc*a(i,nbw+j-i)
 40       x(k)=x(k)-cc*x(i)
 10     continue

        do 70 i=n-1,1,-1
        do 70 k=i+1,min(n,i+nbwm)
 70     x(i)=x(i)-a(i,nbw+k-i)*x(k)
        return
        end
