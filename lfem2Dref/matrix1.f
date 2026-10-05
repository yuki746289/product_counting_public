        subroutine setmat(n,m,sa,sf,dsf,ind,ibw)
        implicit double precision (a-h,o-z)
        dimension sa(ind,ibw),sf(ind),dsf(ind)
        do 10 i=1,n
        sf(i)=0d0
        dsf(i)=0d0
        do 10 j=1,m
 10     sa(i,j)=0d0
      return
      end

        subroutine boundasym(np,nbw,sa,sf,dsf,nbd,kbd,bdv,ind,ibw)
        implicit double precision (a-h,o-z)
        dimension kbd(1001),bdv(1001),sa(ind,ibw),sf(ind),dsf(ind)
c
        do 20 i=1,np		   !応力の釣り合い条件
        sf(i)=sf(i)+dsf(i)
 20     continue
c
        nbwm=nbw-1           !境界上の節点速度0の条件
        do 30 i=1,nbd
          ii=kbd(i)
          if(ii.lt.0) goto 30
          aa=sa(ii,nbw)
                if(aa.eq.0) then
                write(*,*) ii,aa
                stop
                endif
          do 40 j=max(1,ii-nbwm),min(ii+nbwm,np)
          sf(j)=sf(j)-sa(j,nbw+ii-j)*bdv(i)
          sa(ii,nbw+j-ii)=0d0
 40       sa(j,nbw+ii-j)=0d0
          sa(ii,nbw)=aa
          sf(ii)=aa*bdv(i)
 30     continue
        return
        end
c
        subroutine gauss(n,nbw,a,x,ind,ibw)
        implicit double precision (a-h,o-z)
        dimension a(ind,ibw),x(ind)
c
        do 2 i=1,n
        if(a(i,nbw).le.0) then
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
