      subroutine renum(nel,ne,np,jnt,ind,ine,nbw)
c      common /array/ jv(10000,4),mcon(40000),ncon(10000),
c     &               nodes,lments,idiff
      common /array/ jv(40000),mcon(40000),ncon(10000),    !
     &               nnod,nelm,ndif
      dimension ne(ine,3),jntr(ind),jnt(ind)
      integer*2 jv,mcon,ncon
c
co1/18      write(*,*) 'renum'
c
        write(*,*)"in renum"     !
        if(nel.gt.10000) then
                write(*,*) 'error in renum / 10000 < ine=',ine
c               call clos
                stop
        endif   
c            do 389 i=1,10000
c            do 389 j=1,4
c  389       jv(i,j)=0
            do 389 i=1,40000  !
  389       jv(i)=0
            do 400 i=1,nel
            do 400 j=1,3
c  400       jv(i,j)=ne(i,j)
  400       jv(10000*(j-1)+i)=ne(i,j)    !
c            nodes=np
c            lments=nel
            nnod=np      !
            nelm=nel	    !
            call conect
            call autrnm(jnt,nbw,ind)
c            do 430 i=1,np
c              jntr(jnt(i))=i
c  430       continue
      return
      end
c
      subroutine conect
      common /array/ jv(40000),mcon(40000),ncon(10000)
     &               ,nnod,nelm,ndif
      integer*2 jv,mcon,ncon
c
co1/18      write(*,*) 'nnod=',nnod     !
      ndif=nnod
      do 100 j=1,nnod
  100 ncon(j)=0
      do 400 j=1,nelm
        do 300 i=1,4
          jn=jv(10000*(i-1)+j)
          if(jn.eq.0) go to 400
            js=(jn-1)*8
            do 200 i2=1,4
              if(i2.eq.i) go to 200
                jj=jv(10000*(i2-1)+j)
                if(jj.eq.0) go to 300
                  m=ncon(jn)
                  if(m.eq.0) go to 160
                    do 150 i3=1,m
                      if(mcon(js+i3).eq.jj) go to 200
  150               continue
  160             ncon(jn)=ncon(jn)+1
                  mcon(js+ncon(jn))=jj
                  if(iabs(jn-jj).gt.ndif) ndif=iabs(jn-jj)
  200       continue
  300   continue
  400 continue
      return
      end
c
      subroutine autrnm(jnew,nbw,ind)
c      common/array/jv(40000),mcon(40000),ncon(10000),
c     &             nnod,nflm,ndif
      common/array/jv(40000),mcon(40000),ncon(10000),		!
     &              nnod,nelm,ndif
      dimension newj(10000),jold(10000),jnew(ind)
      integer*2 jv,mcon,ncon,newj,jold
      nj=nnod
      mm=ndif
      do 500 l=1,nnod
        do 100 j=1,nnod
          jold(j)=0
  100   newj(j)=0
c
        mx=0
        m=1
        newj(1)=l
        jold(l)=1
        n=1
c
  200   kk=ncon(newj(m))
          if(kk.eq.0) go to 310
            js=(newj(m)-1)*8
            do 300 jj=1,kk
            kn=mcon(js+jj)
            if(jold(kn).gt.0) go to 300
              n=n+1
              newj(n)=kn
              jold(kn)=n
              ldif=iabs(m-n)
              if(ldif.ge.mm) go to 500
                if(ldif.gt.mx) mx=ldif
  300       continue
            if (n.eq.nj) go to 400
  310     m=m+1
          go to 200
c*****
  400   mm=mx
        do 450 j=1,nnod
  450   jnew(j)=jold(j)
  500 continue
c
      nbw=mm+1
      return
      end
