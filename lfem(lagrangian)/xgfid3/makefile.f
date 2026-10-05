c------------------------------------------------------------------------------
c
c       ファイルを作成する
c
c------------------------------------------------------------------------------
        subroutine makefile(np,nele,ne,xx,yy,zz,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),nec(ine,4)
        dimension xx(ind),yy(ind),zz(ind)
        dimension nelb(ine),neb(ine)
c
c       grid.1
        open(1,file="grid.1")
        rewind(1)
        do 10 kp=1,np
10      write(1,*)kp,kp,xx(kp),yy(kp),zz(kp)
        close(1)
c
c       grid.2
        open(1,file="grid.2")
        rewind(1)
        write(1,*)nele,np,0
        do 20 i=1,nele
        write(1,*)i,(ne(i,j),j=1,4)
     &            ,xx(ne(i,1)),yy(ne(i,1)),zz(ne(i,1))
     &            ,xx(ne(i,2)),yy(ne(i,2)),zz(ne(i,2))
     &            ,xx(ne(i,3)),yy(ne(i,3)),zz(ne(i,3))
     &            ,xx(ne(i,4)),yy(ne(i,4)),zz(ne(i,4)),1
20      continue
        close(1)
c
c       bdnum.inf
        open(1,file="bdnum.inf")
        rewind(1)
        do 30 kb=1,11
30      write(1,*)0  !nb(1-11)
        close(1)
c
c       bd.inf
        call findsuf(nele,ne,nelb,neb,ine)  !表面を探す
        open(1,file="bd.inf")
        rewind(1)
        close(1)
c
c       初期条件
        open(1,file="icd.inf")
        rewind(1)
        write(1,*)0d0  !uini
        write(1,*)0d0  !vini
        write(1,*)0d0  !wini
        write(1,*)0d0  !pini
        close(1)
        open(1,file="icdc.inf")
        rewind(1)
        write(1,*)0d0  !cini
        close(1)
        open(1,file="icdh.inf")
        rewind(1)
        write(1,*)0d0  !hini
        close(1)
c
        return
        end
c
c       表面を探す
c
        subroutine findsuf(nele,ne,nelb,neb,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),nec(ine,4),nelb(ine),neb(ine)
        dimension js(4,-1:3)
c
        js(1,-1)=2
        js(1,0)=3
        js(1,1)=1
        js(1,2)=2
        js(1,3)=3
c
        js(2,-1)=3
        js(2,0)=2
        js(2,1)=4
        js(2,2)=3
        js(2,3)=2
c
        js(3,-1)=4
        js(3,0)=1
        js(3,1)=3
        js(3,2)=4
        js(3,3)=1
c
        js(4,-1)=1
        js(4,0)=4
        js(4,1)=2
        js(4,2)=1
        js(4,3)=4
c
c       隣接する要素
        do 30 kele=1,nele
        do 30 j=1,4
30      nec(kele,j)=0
c
        do 100 i=1,nele-1
        do 200 j=1,4
          m1=ne(i,js(j,1))
          m2=ne(i,js(j,2))
          m3=ne(i,js(j,3))
          do k=i+1,nele
          do l=1,4
            n1=ne(k,js(l,1))
            n2=ne(k,js(l,2))
            n3=ne(k,js(l,3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m2.eq.n3 .and. m3.eq.n2 .and. m1.eq.n1) .or.
     &         (m3.eq.n3 .and. m1.eq.n2 .and. m2.eq.n1))then
              nec(i,j)=k
              nec(k,l)=i
              goto 200
            endif
          enddo
          enddo
200     continue
100     continue
c
c       表面
        nb=0
        do 300 i=1,nele
        do 300 j=1,4
          if(nec(i,j).ne.0)goto 300
          nb=nb+1       !表面の数
          nelb(nb)=i    !表面の要素の番号
          neb(nb)=j     !表面の面の番号
300     continue
c        write(*,*)"nb",nb
c
        return
        end

