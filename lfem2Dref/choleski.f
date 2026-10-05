        subroutine choleski2(n,nbw,a,x,ind,ibw)    !np2,nbw2,sa2,sf2,ind,ibw
        implicit double precision (a-h,o-z)
        dimension a(ind,ibw),x(ind),b(ind)
c
c       write(*,*)"a(ind,ibw)=","a(",ind,",",ibw,")"
	  write(*,*)"nbw=",nbw,"n=",n,"in chokeski2"
        if(a(1,nbw).eq.0) then
          write(*,*) '[choleski]',"a(1,nbw)=",a(1,nbw)
          stop
        endif
c       
        nbwm=nbw-1
        do 400 i=1,n
 400	  b(i)=x(i)
c
        write(*,*)"calc a and b"
        do i=2,n
	    do j=i,min(n,i+nbwm)
	      do 500 k=1,i-1  !no consideration of the right side of band reagion
	        if(i-(k-1)+nbwm.gt.min(n,k+nbwm)-(k-1)+nbwm)goto 500
	        p=a(k,it(i,k,nbwm))/a(k,it(k,k,nbwm))
	        if(j.gt.min(n,k+nbwm))goto 500
              a(i,it(j,i,nbwm))=a(i,it(j,i,nbwm))-p*a(k,it(j,k,nbwm))
 500	      continue
	    enddo
	    do 600 k=1,i-1	  !no consideration of the right side of band reagion
	      if(i-(k-1)+nbwm.gt.min(n,k+nbwm)-(k-1)+nbwm)goto 600
	      p=a(k,it(i,k,nbwm))/a(k,it(k,k,nbwm))
            b(i)=b(i)-p*b(k)
 600	    continue
	  enddo
c
        x(n)=b(n)/a(n,it(n,n,nbwm))
	  do i=n-1,1,-1
	    s=b(i)
	    do k=i+1,min(n,i+nbwm)
	      s=s-a(i,it(k,i,nbwm))*x(k)
	    enddo
	    x(i)=s/a(i,it(i,i,nbwm))
	  enddo
        write(*,*)"choreski ended"
        return
        end
c
	  subroutine choleski(n,nbw,a,x,ind,ibw)
        implicit double precision (a-h,o-z)
        dimension a(ind,ibw),x(ind),b(ind),aa(ind,3000)
c
        write(*,*)"a(ind,ibw)=","a(",ind,",",ibw,")"
	  write(*,*)"nbw=",nbw,"n=",n
        if(a(1,nbw).eq.0) then
          write(*,*) '[choleski]',"a(1,nbw)=",a(1,nbw)
          stop
        endif
c
        open(1,file="a.res")
	  open(2,file="aa.res")
	  rewind(1)
	  rewind(2)
        nbwm=nbw-1
        do 100 i=1,n 
	  do 100 j=1,n
100     aa(i,j)=0d0
        do 200 i=1,n
	  do 200 j=i,min(n,i+nbwm)	    
        aa(i,j)=a(i,j-i+nbw)
	  write(2,*)"aa(",i,",",j,")=",aa(i,j)
	  write(1,*)"a(",i,",",j-i+nbw,")=",a(i,j-i+nbw)
200     continue
        pause
c
        do i=1,n
	    s=0d0
	    do j=1,n
	      s=s+abs(aa(i,j))-abs(aa(i,i))
	    enddo
	    if(s.ge.aa(i,i))then
	      write(*,*)"not positive definite"
	      stop
	    endif
	  enddo
c
        do i=1,ind
	    b(i)=x(i)
	  enddo
c
        p=dsqrt(aa(1,1))
	  aa(1,1)=dsqrt(aa(1,1))
	  write(*,*)"aa(1,1)=",aa(1,1)
        do j=2,n
	    aa(1,j)=aa(1,j)/p
	    write(*,*)"aa(1",",",j,")=",aa(1,j)
	  enddo
c
        do i=2,n
	    s=aa(i,i)
	    write(*,*)"a(",i,",",i,")=",aa(i,i)
	    do k=1,i-1
	      s=s-aa(k,i)*aa(k,i)
	    enddo
	    write(*,*)"s=",s
	    p=dsqrt(s)
	    aa(i,i)=dsqrt(s)
	    write(*,*)"aa(",i,",",i,")=",aa(i,i),"p=",p
	    do j=i+1,n
	      s=aa(i,j)
	      do k=1,i-1
	        s=s-aa(k,i)*aa(k,j)
	      enddo
	      write(*,*)"s=",s,"p=",p
	      aa(i,j)=s/p
	      write(*,*)"aa(",i,",",j,")=",aa(i,j)
	    enddo
	  enddo

c
        b(1)=b(1)/aa(1,1)
	  write(*,*)"b(1)=",b(1)
	  do i=2,n
	    s=b(i)
	    do k=1,i-1
	      s=s-aa(k,i)*b(k)
	    enddo
	    b(i)=s/aa(i,i)
	    write(*,*)"b(",i,")=",b(i)
	  enddo
c
        x(n)=b(n)/aa(n,n)
	  do i=n-1,1,-1
	    s=b(i)
	    do k=i+1,n
	      s=s-aa(i,k)*x(k)
	    enddo
	    x(i)=s/aa(i,i)
	  enddo
	  return
	  end
c
        function it(j,i,nbwm)
	    it=j-(i-1)+nbwm
	  return 
	  end