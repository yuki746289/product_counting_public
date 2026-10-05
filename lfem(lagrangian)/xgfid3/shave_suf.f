c---------------------------------------------------------------------------
c
c       面ができているか確認する
c       面がない場合は、表面に節点を追加して要素を作り直す
c---------------------------------------------------------------------------
        subroutine decorate(np,nele,ne,xx,yy,zz,ns,ms,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),ms(ine,3),mk(ine)
        dimension xx(ind),yy(ind),zz(ind)
c
        nk=0    !無い面
        do 100 ks=1,ns
c
          do kele=1,nele
            kk=0
            do j=1,4
              j1=ne(kele,j)
              do i=1,3
                i1=ms(ks,i)
                if(i1.eq.j1)kk=kk+1
              enddo
            enddo
            if(kk.eq.3)goto 100
          enddo
          nk=nk+1
          mk(nk)=ks
100     continue
        if(nk.eq.0)return
        write(*,*)"no suf"
        stop
c
c       節点を追加して要素を作り直す
c        do 200 kk=1,nk
c          k1=ms(mk(kk,1))
cc          k2=ms(mk(kk,2))
c          k3=ms(mk(kk,3))
c          xp=(xx(k1)+xx(k2)+xx(k3))/3d0
c          yp=(yy(k1)+yy(k2)+yy(k3))/3d0
c          zp=(zz(k1)+zz(k2)+zz(k3))/3d0

c
        return
        end
c---------------------------------------------------------------------------
c
c       外側の要素を削る
c
c---------------------------------------------------------------------------
        subroutine shavesuf(np,nele,ne,xx,yy,zz,ns,ms,ind,ine)
        implicit double precision(a-h,o-z)
        dimension ne(ine,4),ms(ine,3),msr(ine)
        dimension xx(ind),yy(ind),zz(ind)
        dimension jm1(3),jm2(3),jm3(3),jm4(3)  !4面体の面の番号
        data jm1/4,2,3/    !面1の節点の番号:ne(nele,jm1(1-3))
        data jm2/4,3,1/    !面2の節点の番号:ne(nele,jm2(1-3))
        data jm3/4,1,2/    !面3の節点の番号:ne(nele,jm3(1-3))
        data jm4/1,3,2/    !面4の節点の番号:ne(nele,jm4(1-3))
c
c       外側の要素を削る    ※ms(ns)は内向き
        nels=0  !外側の要素の数
        do 100 ks=1,ns
          m1=ms(ks,1)
          m2=ms(ks,2)
          m3=ms(ks,3)
          do kele=1,nele
c
c           面1
            n1=ne(kele,jm1(1))
            n2=ne(kele,jm1(2))
            n3=ne(kele,jm1(3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m1.eq.n1 .and. m2.eq.n3 .and. m3.eq.n2) .or.
     &         (m1.eq.n2 .and. m2.eq.n1 .and. m3.eq.n3))then  !(1,2,3)=(3,2,1),(1,3,2),(2,1,3)
              nels=nels+1
              msr(nels)=kele
              goto 100
            endif
c
c           面2
            n1=ne(kele,jm2(1))
            n2=ne(kele,jm2(2))
            n3=ne(kele,jm2(3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m1.eq.n1 .and. m2.eq.n3 .and. m3.eq.n2) .or.
     &         (m1.eq.n2 .and. m2.eq.n1 .and. m3.eq.n3))then  !(1,2,3)=(3,2,1),(1,3,2),(2,1,3)
              nels=nels+1
              msr(nels)=kele
              goto 100
            endif
c
c           面3
            n1=ne(kele,jm3(1))
            n2=ne(kele,jm3(2))
            n3=ne(kele,jm3(3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m1.eq.n1 .and. m2.eq.n3 .and. m3.eq.n2) .or.
     &         (m1.eq.n2 .and. m2.eq.n1 .and. m3.eq.n3))then  !(1,2,3)=(3,2,1),(1,3,2),(2,1,3)
              nels=nels+1
              msr(nels)=kele
              goto 100
            endif
c
c           面4
            n1=ne(kele,jm4(1))
            n2=ne(kele,jm4(2))
            n3=ne(kele,jm4(3))
            if((m1.eq.n3 .and. m2.eq.n2 .and. m3.eq.n1) .or.
     &         (m1.eq.n1 .and. m2.eq.n3 .and. m3.eq.n2) .or.
     &         (m1.eq.n2 .and. m2.eq.n1 .and. m3.eq.n3))then  !(1,2,3)=(3,2,1),(1,3,2),(2,1,3)
              nels=nels+1
              msr(nels)=kele
              goto 100
            endif
          enddo
          if(kele.eq.nele)write(*,*)"not ns"
          if(kele.eq.nele)stop
100     continue
c
c       外側の要素を削る
        do kels=1,nels
          kele=msr(kels)
          ne(kele,1)=ne(nele,1) !最後の要素を持って来る:最後の要素が外側の要素だったら削除されない
          ne(kele,2)=ne(nele,2)
          ne(kele,3)=ne(nele,3)
          ne(kele,4)=ne(nele,4)
          nele=nele-1
        enddo
c
        return
        end