        subroutine renum2(index,np,ne,nel,jww,nband)
        implicit double precision (a-h,o-z)
        parameter(ind=28001,ine=18001)
c        integer*2 ne(ine,50)
        dimension ne(ine,3)
        integer*2 psurp(50,ind)
c used in some condition(integer*2)
        integer*2 nesune(ine,50)
        dimension jw(ind),jww(ind)
        dimension nl(ind),nl0(ind)
        dimension npsurp(ind)
c used in some condition(dimension)
        dimension pnesune(ine)
        dimension nebnd(ine)
        dimension pnearb(ind)
        dimension pbnd(ind)
c
c  get psurp(50,np),npsurp(np) arraies
        call psurroundp(index,np,nel,ne,psurp,npsurp,ine,ind)
c        write(*,*)"got psurp(nel,50),npsurp(nel)"
c for some conditions
c  get nesune(nel,50),pnesune(nel) arraies
        call nesurroundne(nel,ne,nesune,pnesune,ine,ind)
c        write(*,*)"got nesune(nel,50),pnesune(nel)"
c  get nebnd(nbnd),nbnd
        call nebound(nel,pnesune,nebnd,nbnd,ine,ind)
c        write(*,*)"got nebnd(nbnd) nbnd=",nbnd
c  get pnearb(npnearb),npnearb
c        call pnearbnd(index,ne,nel,psurp,npsurp,nebnd,
c     $                     nbnd,pnearb,npnearb,ine,ind)
c        write(*,*)"got pnearb(npnearb) npnearb=",npnearb
c  get pbnd(npbnd) array
c        call findbnd2el(index,np,nel,ne,pbnd,npbnd,ine,ind)
c        call findbnd3el(index,np,nel,ne,pbnd,npbnd,ine,ind)
c        write(*,*)"got pbnd(npbnd) array npbnd=",npbnd
c  find band
        if(index.eq.1)then      !index denotes the element order
        npe=3                   !npe denotes the number of nodes in one element
        elseif(index.eq.2)then
        npe=6
        elseif(index.eq.3)then
        npe=10
        else
        write(*,*)"element order=",index,"in renum2"
        stop
        endif
c        write(*,*)"element order=",index,"npe=",npe
c
        mp=9999
        mprt=-9999
	  nband=9999
        do 40 m=1,np
c
c some condition
c        if(npsurp(m).gt.3)go to 40           !only nodex combined with less than 3 nodes
c        !serch only nodes from middle upward
c        do ii=1,nbnd        !consider nodes existing near boundition(?) only
c           if(m.eq.pnearb(ii))go to 1300
c        enddo
c
c         do ii=1,npbnd       !consider boundary nodes
c            if(m.eq.pbnd(ii))go to 1300
c         enddo
c
c        do ii=1,nbnd
c        do j=1,6  !or 10
c            do k=1,npsurp(ne(nebnd(ii),j))      !consider nodes extistiong near boundition
c            if(m.eq.psurp(k,ne(nebnd(ii),j)))go to 1300
c            enddo
c           if(m.eq.ne(nebnd(ii),j))go to 1300  !consider nodes existing near boudition
c           if(m.eq.ne(nebnd(ii),j))then        !consider nodes existing boundition
c              if(npsurp(m).le.3)go to 1300
c              if(j.le.3 .and. npsurp(m).le.4)then   !for 2order element
c              go to 1300
c              elseif(j.gt.3 .and. npsurp(m).le.2)then
c              go to 1300
c              endif
c           endif
c        enddo
c        enddo
c        go to 40
c1300    continue
c
c
	do 80 i=1,np
 80	jw(i)=0
	jw(m)=1
	nv0=1
	nl0(1)=m
	nn=1
	do 100 nlev=2,np
	  nv=0
	  do 120 i=1,nv0
	    n=nl0(i)
	        do 180 l=1,npsurp(n)
 		  if(jw(psurp(l,n)).eq.0) then
		  nn=nn+1
		  jw(psurp(l,n))=nn
	          nv=nv+1
		  nl(nv)=psurp(l,n)
	          endif
 180	        continue
 120	  continue
          if(nv.gt.mprt)mprt=nv
          if(mprt.gt.mp)then   !ge or gt
          mprt=-9999
          go to 40
          endif
	  if(nv.eq.0) goto 2000
	  do 240 i=1,nv
 240	  nl0(i)=nl(i)
	  nv0=nv
 100	continue
c
 2000   mp=mprt 
c
c       find band
c
	nbwt=0
	do 200 i=1,nel
	do 200 j=1,npe
	do 200 k=1,npe
c 200	nbwt=max0(nbwt,iabs(jw(ne(i,j))-jw(ne(i,k))))
        if(iabs(jw(ne(i,j))-jw(ne(i,k))).gt.nbwt)
     $         nbwt=iabs(jw(ne(i,j))-jw(ne(i,k)))
 200    continue
	if(nbwt.lt.nband) then
        nband=nbwt      !the smallest band width at m root node
        nroot=m         !the old root node 
c        write(*,*)"R max",nband,"root node",m
	do 220 i=1,np
 220	jww(i)=jw(i)           !jww(np) means new(nel,j)
        endif
c
 40	continue
        nband=nband+1
        write(*,*)"the smallest band width(R+1)",nband
c        write(*,*)"old root node",nroot
c        write(*,*)"x point",xx(nroot),"y point",yy(nroot)
c
        return
        end

        subroutine psurroundp(index,np,nel,ne,psurp,npsurp,ine,ind)
	  implicit double precision (a-h,o-z)
c
c        integer*2 ne(ine,50)
        dimension ne(ine,3)
        integer*2 psurp(50,ind)
        integer*2 jn(50,50),jnn(50,50,50)
        dimension npsurp(ind)
c get jn(50,50)
c one order element
        jnn(1,1,1)=2
        jnn(1,1,2)=3
        jnn(1,2,1)=1
        jnn(1,2,2)=3
        jnn(1,3,1)=1
        jnn(1,3,2)=2
c two order element
        jnn(2,1,1)=4
        jnn(2,1,2)=6
        jnn(2,2,1)=4
        jnn(2,2,2)=5
        jnn(2,3,1)=5
        jnn(2,3,2)=6
        jnn(2,4,1)=1
        jnn(2,4,2)=2
        jnn(2,5,1)=2
        jnn(2,5,2)=3
        jnn(2,6,1)=1
        jnn(2,6,2)=3
c three order element
        jnn(3,1,1)=4
        jnn(3,1,2)=9
        jnn(3,1,3)=10
        jnn(3,2,1)=5
        jnn(3,2,2)=6
        jnn(3,2,3)=10
        jnn(3,3,1)=7
        jnn(3,3,2)=8
        jnn(3,3,3)=10
        jnn(3,4,1)=1
        jnn(3,4,2)=5
        jnn(3,4,3)=10
        jnn(3,5,1)=2
        jnn(3,5,2)=4
        jnn(3,5,3)=10
        jnn(3,6,1)=2
        jnn(3,6,2)=7
        jnn(3,6,3)=10
        jnn(3,7,1)=3
        jnn(3,7,2)=6
        jnn(3,7,3)=10
        jnn(3,8,1)=3
        jnn(3,8,2)=9
        jnn(3,8,3)=10
        jnn(3,9,1)=1
        jnn(3,9,2)=8
        jnn(3,9,3)=10
        jnn(3,10,1)=1
        jnn(3,10,2)=2
        jnn(3,10,3)=3
        jnn(3,10,4)=4
        jnn(3,10,5)=5
        jnn(3,10,6)=6
        jnn(3,10,7)=7
        jnn(3,10,8)=8
        jnn(3,10,9)=9
        do i=1,50
        do j=1,50
        jn(i,j)=jnn(index,i,j)
        enddo
        enddo
c get npe
        if(index.eq.1)then
        npe=3
        nsp=2
        elseif(index.eq.2)then
        npe=6
        nsp=2
        elseif(index.eq.3)then
        npe=10
        nsp=3
        else
        write(*,*)"index=",index,"in psrroundp"
        stop
        endif
c
        do 50 i=1,np	
        do 60 j=1,50
        psurp(j,i)=-999
 60     continue
        npsurp(i)=-999
 50     continue
c
        m=0
        do 700 kp=1,np
           do 600 i=1,nel
           do 300 j=1,npe
              if(kp.eq.ne(i,j))then
                do 400 l=1,nsp
                do 500 k=1,50
                if(psurp(k,kp).eq.ne(i,jn(j,l))) go to 400
 500            continue
                m=m+1
                psurp(m,kp)=ne(i,jn(j,l))
 400            continue
              endif
 300       continue
 600       continue
           npsurp(kp)=m
           m=0
 700    continue
c
        return
        end

       subroutine nesurroundne(nel,ne,nesune,pnesune,ine,ind)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50)
       dimension ne(ine,3)
       integer*2 nesune(ine,50)
       dimension pnesune(ine)
       dimension j1(3)
c
       j1(1)=2
       j1(2)=3
       j1(3)=1
c
       npe=3
       nnel=0
       do 800 i=1,nel
       do 900 j=1,npe
          do 1000 k=1,nel
          do 1100 l=1,npe
          if(ne(i,j).eq.ne(k,l) .and. i.ne.k)then
          do m=1,3
             if(ne(i,j1(j)).eq.ne(k,m))then
             nnel=nnel+1
             nesune(i,nnel)=k
             endif
          enddo
          endif
 1100     continue
 1000     continue
 900   continue
       pnesune(i)=nnel
       nnel=0
 800   continue
c
       return
       end

       subroutine nebound(nel,pnesune,nebnd,nbnd,ine,ind)
	 implicit double precision (a-h,o-z)
c
       dimension pnesune(ine)
       dimension nebnd(ine)
c
       nbnd=0
       do 1200 i=1,nel
       if(pnesune(i).le.2)then
       nbnd=nbnd+1
       nebnd(nbnd)=i
       endif
 1200  continue
c
       return
       end
c
       subroutine pnearbnd(index,ne,nel,psurp,npsurp,nebnd
     $              ,nbnd,pnearb,npnearb,ine,ind)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50)
       dimension ne(ine,3)
       integer*2 psurp(50,ind)
       dimension npsurp(ind)
       dimension nebnd(ine)
       dimension pnearb(ind)
c
       if(index.eq.1)then
       npe=3
       elseif(index.eq.2)then
       npe=6
       elseif(index.eq.3)then
       npe=10
       else
       write(*,*)"index=",index,"in pnearbnd"
       stop
       endif
       l=0
       do 1400 i=1,nbnd
       do 1500 j=1,npe
          do 1600 k=1,npsurp(ne(nebnd(i),j))
             if(l.ge.1)then
                do m=1,l
                if(psurp(k,ne(nebnd(i),j)).eq.pnearb(m))go to 1600
                enddo
             endif
             l=l+1
             pnearb(l)=psurp(k,ne(nebnd(i),j))
 1600     continue
 1500  continue
 1400  continue
       npnearb=l
c
       return
       end

       subroutine generate2el(np,ne,nel,xx,yy,ine,ind,
     $                 np2,ne2,xx2,yy2,nbwt)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50),ne2(ine,50)
	 integer*2 nesune(ine,50)
       dimension ne(ine,3)
c       integer*2 ne2(ine,50)
       dimension ne2(ine,6)
       dimension xx(ind),yy(ind)
       dimension xx2(ind),yy2(ind)
       dimension pnesune(ine)
       dimension jn1(3,2),jl(3,3),j1(3)
c
       jn1(1,1)=2
       jn1(1,2)=3
       jn1(2,1)=1
       jn1(2,2)=3
       jn1(3,1)=1
       jn1(3,2)=2
c
       jl(1,2)=4
       jl(1,3)=6
       jl(2,1)=4
       jl(2,3)=5
       jl(3,1)=6
       jl(3,2)=5
c
       j1(1)=2
       j1(2)=3
       j1(3)=1
c
       index=2
       nptmp=np
       do i=1,nel
       do j=1,6
         if(j.le.3)then
           ne2(i,j)=ne(i,j)
           xx2(ne(i,j))=xx(ne(i,j))
           yy2(ne(i,j))=yy(ne(i,j))
         elseif(j.ge.4)then
           ne2(i,j)=0
         endif
       enddo
       enddo
       do i=np+1,ind
         xx2(i)=0d0
         yy2(i)=0d0
       enddo
c
c generate two order mesh
       call nesurroundne(nel,ne,nesune,pnesune,ine,ind)
c       write(*,*)"got nesune(nel,50),pnesune(nel) arraies"
c
       npe=3
       do i=1,nel
       do j=1,npe
         m=j1(j)
         do k=1,pnesune(i)        !about inner element
           do l=1,npe
           if(ne2(i,j).eq.ne2(nesune(i,k),l))then
             do n=1,2
             if(ne2(i,m).eq.ne2(nesune(i,k),jn1(l,n)))then
               if(ne2(nesune(i,k),jl(l,jn1(l,n))).eq.0)then
               np=np+1
               ne2(i,jl(j,m))=np
               xx2(np)=(xx2(ne2(i,j))+xx2(ne2(i,m)))/2d0
               yy2(np)=(yy2(ne2(i,j))+yy2(ne2(i,m)))/2d0
               elseif(ne2(nesune(i,k),jl(l,jn1(l,n))).ne.0)then
               ne2(i,jl(j,m))=ne2(nesune(i,k),jl(l,jn1(l,n)))
               endif
             endif
             enddo
           endif
           enddo
         enddo
       enddo
       if(pnesune(i).le.2)then     !about boudary element
         do j=1,3
         m=j1(j)
         if(ne2(i,jl(j,m)).eq.0)then
           np=np+1
           ne2(i,jl(j,m))=np
           xx2(np)=(xx2(ne2(i,j))+xx2(ne2(i,m)))/2d0
           yy2(np)=(yy2(ne2(i,j))+yy2(ne2(i,m)))/2d0
         endif
         enddo
       endif
       enddo
c
       np2=np
       np=nptmp
c      find band width
	 nbwt=0
	 do 200 i=1,nel
	 do 200 j=1,npe
	 do 200 k=1,npe
c 200	 nbwt=max0(nbwt,iabs(jw(ne(i,j))-jw(ne(i,k))))
         if(iabs(ne(i,j)-ne(i,k)).gt.nbwt)
     $          nbwt=iabs(ne(i,j)-ne(i,k))
 200   continue
       nbwt=nbwt+1
c	  
       return 
       end



       subroutine generate3el(np,ne,nel,xx,yy,ine,ind,
     $                 np3,ne3,xx3,yy3,nbwt)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50),ne3(ine,50)
       dimension ne(ine,3)
       integer*2 ne3(ine,50)
       integer*2 nesune(ine,50)
       dimension xx(ind),yy(ind)
       dimension xx3(ind),yy3(ind)
       dimension pnesune(ine)
       dimension jn1(3,2),jl(3,3,2),j1(3)
c
       jn1(1,1)=2
       jn1(1,2)=3
       jn1(2,1)=1
       jn1(2,2)=3
       jn1(3,1)=1
       jn1(3,2)=2
c
       jl(1,2,1)=4
       jl(2,1,1)=4
       jl(1,2,2)=5
       jl(2,1,2)=5
       jl(2,3,1)=6
       jl(3,2,1)=6
       jl(2,3,2)=7
       jl(3,2,2)=7
       jl(3,1,1)=8
       jl(1,3,1)=8
       jl(3,1,2)=9
       jl(1,3,2)=9
c
       j1(1)=2
       j1(2)=3
       j1(3)=1
c
       index=3
       nptmp=np
       do i=1,nel
       do j=1,10
         if(j.le.3)then
           ne3(i,j)=ne(i,j)
           xx3(ne(i,j))=xx(ne(i,j))
           yy3(ne(i,j))=yy(ne(i,j))
         elseif(j.ge.4)then
           ne3(i,j)=0
         endif
       enddo
       enddo
       do i=nptmp,ind
         xx3(i)=0d0
         yy3(i)=0d0
       enddo
       nptmp=np
c
c generate three order mesh
       call nesurroundne(nel,ne,nesune,pnesune,ine,ind)
c       write(*,*)"got nesune(nel,50),pnesune(nel) arraies"
c
       npe=3
       do i=1,nel
       do j=1,npe
         m=j1(j)
         do k=1,pnesune(i)        !about inner element
           do l=1,npe
           if(ne3(i,j).eq.ne3(nesune(i,k),l))then
             do n=1,2
             if(ne3(i,m).eq.ne3(nesune(i,k),jn1(l,n)))then
               do nn=1,2
               if(ne3(nesune(i,k),jl(l,jn1(l,n),nn)).eq.0)then
                 np=np+1
                 ne3(i,jl(j,m,nn))=np
                 if(nn.eq.1)then
                   xx3(np)=(2*xx3(ne3(i,j))+xx3(ne3(i,m)))/3d0
                   yy3(np)=(2*yy3(ne3(i,j))+yy3(ne3(i,m)))/3d0
                 elseif(nn.eq.2)then
                   xx3(np)=(xx3(ne3(i,j))+2*xx3(ne3(i,m)))/3d0
                   yy3(np)=(yy3(ne3(i,j))+2*yy3(ne3(i,m)))/3d0
                 endif
               elseif(ne3(nesune(i,k),jl(l,jn1(l,n),nn)).ne.0)then
                 ne3(i,jl(j,m,nn))=
     $                ne3(nesune(i,k),jl(l,jn1(l,n),nn))
               endif
               enddo
             endif
             enddo
           endif
           enddo
         enddo
       enddo
       if(pnesune(i).le.2)then     !about boudary element
         do j=1,3
         m=j1(j)
         do nn=1,2
           if(ne3(i,jl(j,m,nn)).eq.0)then
             np=np+1
             ne3(i,jl(j,m,nn))=np
             if(nn.eq.1)then
               xx3(np)=(2*xx3(ne3(i,j))+xx3(ne3(i,m)))/3d0
               yy3(np)=(2*yy3(ne3(i,j))+yy3(ne3(i,m)))/3d0
             elseif(nn.eq.2)then
               xx3(np)=(xx3(ne3(i,j))+2*xx3(ne3(i,m)))/3d0
               yy3(np)=(yy3(ne3(i,j))+2*yy3(ne3(i,m)))/3d0
             endif
           endif
         enddo
         enddo
       endif
       np=np+1
       ne3(i,10)=np
       xx3(np)=(xx(ne3(i,1))+xx(ne3(i,2))+xx(ne3(i,3)))/3d0
       yy3(np)=(yy(ne3(i,1))+yy(ne3(i,2))+yy(ne3(i,3)))/3d0
       enddo
c
       np3=np
       np=nptmp
c      find band width
	 nbwt=0
	 do 200 i=1,nel
	 do 200 j=1,10
	 do 200 k=1,10
c 200	 nbwt=max0(nbwt,iabs(jw(ne(i,j))-jw(ne(i,k))))
         if(iabs(ne(i,j)-ne(i,k)).gt.nbwt)
     $          nbwt=iabs(ne(i,j)-ne(i,k))
 200   continue
       nbwt=nbwt+1
c
       return
       end

       subroutine findbnd2el(index,np,nel,ne,pbnd2,npbnd2,ine,ind)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50)
       dimension ne(ine,3)
       integer*2 psurp(50,ind)
       integer*2 nesune(ine,50)
       dimension pbnd2(ind)
       dimension npsurp(ind)
       dimension pnesune(ine)
       dimension nebnd(ine)
       dimension j1(3),jl(3,3),jj(3)
c
       j1(1)=2
       j1(2)=3
       j1(3)=1
       jl(1,2)=4
       jl(2,1)=4
       jl(2,3)=5
       jl(3,2)=5
       jl(3,1)=6
       jl(3,2)=6
c
       call psurroundp(index,np,nel,ne,psurp,npsurp,ine,ind)
       call nesurroundne(nel,ne,nesune,pnesune,ine,ind)
       call nebound(nel,pnesune,nebnd,nbnd,ine,ind)
c
       k=0
       do 1700 i=1,nbnd
       do 1800 j=1,3
         if(npsurp(ne(nebnd(i),j)).le.5 .and.
     $                npsurp(ne(i,j1(j))).le.5)then
           jj(1)=j
           jj(2)=j1(j)
           jj(3)=jl(j,j1(j))
           if(k.eq.0)then
             do m=1,3
             k=k+1
             pbnd2(k)=npsurp(ne(nebnd(i),jj(m)))
             enddo
           elseif(k.ge.1)then
           do 1900 m=1,3
             do l=1,k
             if(npsurp(ne(nebnd(i),jj(m))).eq.pbnd2(k))goto 1900 
             enddo
             k=k+1
             pbnd2(k)=npsurp(ne(nebnd(i),jj(m)))
 1900      continue
           endif
         endif
 1800  continue
 1700  continue      
       npbnd2=k
c
       return
       end

       subroutine findbnd3el(index,np,nel,ne,pbnd3,npbnd3,ine,ind)
	 implicit double precision (a-h,o-z)
c
c       integer*2 ne(ine,50)
       dimension ne(ine,3)
       integer*2 psurp(50,ind)
       integer*2 nesune(ine,50)
       dimension pbnd3(ind)
       dimension npsurp(ind)
       dimension pnesune(ine)
       dimension nebnd(ine)
       dimension j1(3),jl(3,3,2),jj(4)
c
       j1(1)=2
       j1(2)=3
       j1(3)=1
       jl(1,2,1)=4
       jl(2,1,1)=4
       jl(1,2,2)=5
       jl(2,1,2)=5
       jl(2,3,1)=6
       jl(3,2,1)=6
       jl(2,3,2)=7
       jl(3,2,2)=7
       jl(3,1,1)=8
       jl(1,3,1)=8
       jl(3,1,2)=9
       jl(1,3,2)=9
c
       call psurroundp(index,np,nel,ne,psurp,npsurp,ine,ind)
       call nesurroundne(nel,ne,nesune,pnesune,ine,ind)
       call nebound(nel,pnesune,nebnd,nbnd,ine,ind)
c
       k=0
       do 2100 i=1,nbnd
       do 2000 j=1,3
         if(npsurp(ne(nebnd(i),j)).le.9 .and.
     $        npsurp(ne(nebnd(i),j1(j))).le.9)then
           jj(1)=j
           jj(2)=j1(j)
           jj(3)=jl(j,j1(j),1)
           jj(4)=jl(j,j1(j),2)
           if(k.eq.0)then
              do m=1,4
              k=k+1
              pbnd3(k)=npsurp(ne(nebnd(i),jj(m)))
              enddo
           elseif(k.ge.1)then 
           do 2200 m=1,4
             do l=1,k
             if(npsurp(ne(nebnd(i),jj(m))).eq.pbnd3(k))goto 2200
             enddo
             k=k+1
             pbnd3(k)=npsurp(ne(nebnd(i),jj(m)))
 2200      continue
           endif
         endif
 2000  continue
 2100  continue
       npbnd3=k
c
       return
       end
c
c
       subroutine genet2el(np,nel,ne,xx,yy,np2,ne2,xx2,yy2,nbwt)
	 implicit double precision(a-h,o-z)
       parameter(ind=28001,ine=18001)
	 integer*2 mb(ine,ine)
	 dimension ne(ine,3),xx(ind),yy(ind)
	 dimension ne2(ine,6),xx2(ind),yy2(ind)
	 dimension j1(3)
	 data j1/2,3,1/
c
       do 10 i=1,nel
	 do 10 j=1,nel
10     mb(i,j)=-999
c
       do 20 i=1,nel
	 do 20 j=1,6
	   if(j.le.3)then
	   ne2(i,j)=ne(i,j)
	   else
	   ne2(i,j)=0
         endif
20     continue
c
       do 30 i=1,ind
	   if(i.le.np)then
	   xx2(i)=xx(i)
	   yy2(i)=yy(i)
	   else
	   xx2(i)=0d0
	   yy2(i)=0d0
         endif
30     continue
c
       np2=np
	 do 100 im=1,nel
	 do 200 j=1,3
	   do 300 im2=1,nel
	   do 400 k=1,3
	     if(ne2(im,j).eq.ne2(im2,j1(k)) .and. 
     &                     ne2(im,j1(j)).eq.ne2(im2,k))then
c	       if(mb(max(im,im2),min(im,im2)).eq.-999)then
	       if(mb(im,im2).eq.-999)then
	         np2=np2+1
	         ne2(im,j+3)=np2
	         ne2(im2,k+3)=np2
	         xx2(np2)=(xx2(ne2(im,j))+xx2(ne2(im,j1(j))))/2d0
	         yy2(np2)=(yy2(ne2(im,j))+yy2(ne2(im,j1(j))))/2d0
c	         mb(max(im,im2),min(im,im2))=0
	         mb(im,im2)=0
	         mb(im2,im)=0
	         goto 200
	       else
	         goto 200
	       endif
	     endif
400      continue
300      continue
         np2=np2+1   !«E
	   ne2(im,j+3)=np2
	   xx2(np2)=(xx2(ne2(im,j))+xx2(ne2(im,j1(j))))/2d0
	   yy2(np2)=(yy2(ne2(im,j))+yy2(ne2(im,j1(j))))/2d0
200    continue
100    continue
c
c      find band width
	 nbwt=0
	 do 500 i=1,nel
	 do 500 j=1,6
	 do 500 k=1,6
c 200	 nbwt=max0(nbwt,iabs(jw(ne(i,j))-jw(ne(i,k))))
         if(iabs(ne2(i,j)-ne2(i,k)).gt.nbwt)
     $          nbwt=iabs(ne2(i,j)-ne2(i,k))
 500   continue
       nbwt=nbwt+1

       return
	 end