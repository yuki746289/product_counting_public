        subroutine coalesce(nelt,ne,xx,yy,irx,idx,icoal,ind,ine)
        implicit double precision (a-h,o-z)
        parameter (iar=8001)
        dimension ne(ine,3),xx(ind),yy(ind),irx(0:ine),idx(0:ine)
        dimension ar(iar)
        dimension ifg(2),lc(2),jc(2),xc(2),yc(2),si(2),sk(2)
        dimension kelb(1001),keb(1001),kcop(1001)
        dimension js2(3)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /regcon/ncon(21),kcon(21,21)
        common /sufcood/xb(21,1001),yb(21,1001)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
        data js2/2,3,1/,epc/1d-6/,ep/-1d-12/
! nelb,nebは要素の辺ではなく節点を示すと考える．
c
c       write(*,*) 'coalesce/nelt=',nelt
        do 30 k=1,nelt
          ar(k)=(xx(ne(k,2))-xx(ne(k,1)))*(yy(ne(k,3))-yy(ne(k,1)))
     &      -(xx(ne(k,3))-xx(ne(k,1)))*(yy(ne(k,2))-yy(ne(k,1)))
 30     continue
c
        icoal=0
        do 60 i=1,nr
          if(.NOT.(btest(mswi(ndx(i)),0))) goto 60      !流動領域を対象
c       write(*,*) 'nr,idx=',i,ndx(i)
 70       irc=-9999
          icy=nb(i)
          do 80 jst=1,nb(i)                             !始点の探索
c           write(*,*) 'jst,ncop=',jst,ncop(i,jst)
            if(ncop(i,jst).ne.0) goto 90
            n=ne(nelb(i,jst),neb(i,jst))
            x=xx(n)
            y=yy(n)
            call checkcoal(nelt,ne,xx,yy,irx,n,irc,inder,ar,epc,
     &                  ine,ind,iar)
            if(inder.eq.0) goto 90
            irc=inder
 80       continue
                write(*,*) 'coalesce/ error at DO 80'
                stop
 90       ifg(1)=0
c         write(*,*) 'jst,nb(i)=',jst,nb(i)
          do 100 j=jst+1,nb(i)+1
            js=mcy(j,icy)
            if(ncop(i,js).ne.0) then
              ifg(1)=-1
            else
              n=ne(nelb(i,js),neb(i,js))
              call checkcoal(nelt,ne,xx,yy,irx,n,irc,inder,ar,epc,
     &                  ine,ind,iar)
              if(inder.eq.0) then
                ifg(1)=0
              else
                jc(1)=j-1
                k=inder
                irc=k
                if(ifg(1).eq.0) ifg(1)=1        !交差
                goto 110
              endif
            endif
 100      continue
          goto 60
c
 110      ifg(2)=0
          write(*,*) '1/jc,ifg,k=',jc(1),ifg(1),k
          do 120 j=jc(1)+1,jc(1)+nb(i)
            js=mcy(j,icy)
c               write(*,*) '2/ncop',j,js,ncop(i,js)
            if(ncop(i,js).ne.0) then
              ifg(2)=-1                         !接する
              jc(2)=j           !-1
              goto 140
            else
              n=ne(nelb(i,js),neb(i,js))
              call checkcoal(nelt,ne,xx,yy,irx,n,irc,inder,ar,epc,
     &                  ine,ind,iar)
c               write(*,*) '2/inder',j,inder
              if(inder.eq.0) then
                jc(2)=j-1
                ifg(2)=1                        !交差
                goto 140
              endif
            endif
 120      continue
                  write(*,*) 'coalesce/ error at DO 120'
                  stop
c
 140      kcy=nb(k)
          write(*,*) '2/jc,ifg,k=',jc(2),ifg(2),k
          lst=1
          do 160 m=2,1,-1
            if(ifg(m).eq.-1) then
              j=jc(m)
              js=mcy(j,icy)
              do 170 l=lst,lst+nb(k)
                ls=mcy(l,kcy)
                if(ne(nelb(i,js),neb(i,js)).eq.
     &                          ne(nelb(k,ls),neb(k,ls))) goto 200
 170          continue
                write(*,*) 'coalesce /error at DO 170 /m=',m
                stop
            else
              j=jc(m)
              js=mcy(j,icy)
              jsp=mcy(j+1,icy)
              xi1=xx(ne(nelb(i,js),neb(i,js)))
              yi1=yy(ne(nelb(i,js),neb(i,js)))
              xi2=xx(ne(nelb(i,jsp),neb(i,jsp)))
              yi2=yy(ne(nelb(i,jsp),neb(i,jsp)))
c       call pltb(ne,xx,yy,xi1,yi1,xi2,yi2,ind,ine)
c       pause
              do 180 l=lst,lst+nb(k)
                ls=mcy(l,kcy)
                lsp=mcy(l+1,kcy)
                xk1=xx(ne(nelb(k,ls),neb(k,ls)))
                yk1=yy(ne(nelb(k,ls),neb(k,ls)))
                xk2=xx(ne(nelb(k,lsp),neb(k,lsp)))
                yk2=yy(ne(nelb(k,lsp),neb(k,lsp)))
c       call pltb(ne,xx,yy,xk1,yk1,xk2,yk2,ind,ine)
c       pause
                call crss(xi1,yi1,xi2,yi2,xk1,yk1,xk2,yk2,
     &                  xc(m),yc(m),si(m),sk(3-m),inder,ep)
c               write(*,*) '3/l,inder',l,inder
 180            if(inder.eq.1) goto 200
            endif
 200      lc(3-m)=l
          lst=l
 160      continue
c
          write(*,*) 'jc1,jc2,lc1,lc2',jc(1),jc(2),lc(1),lc(2)
          icoal=1
          if(ifg(1).eq.1) then
            j1=jc(1)
            if(si(1).gt.1-1d-4) j1=j1+1
            xx(ne(nelb(i,mcy(j1,icy)),neb(i,mcy(j1,icy))))=xc(1)
            yy(ne(nelb(i,mcy(j1,icy)),neb(i,mcy(j1,icy))))=yc(1)
            l2=lc(2)+1
            if(sk(2).lt.1d-4) l2=l2-1
            xx(ne(nelb(k,mcy(l2,kcy)),neb(k,mcy(l2,kcy))))=xc(1)
            yy(ne(nelb(k,mcy(l2,kcy)),neb(k,mcy(l2,kcy))))=yc(1)
          else
           j1=jc(1)
           l2=lc(2)
          endif
          if(ifg(2).eq.1) then
            j2=jc(2)+1
            if(si(2).lt.1d-4) j2=j2-1
            xx(ne(nelb(i,mcy(j2,icy)),neb(i,mcy(j2,icy))))=xc(2)
            yy(ne(nelb(i,mcy(j2,icy)),neb(i,mcy(j2,icy))))=yc(2)
            l1=lc(1)
            if(sk(1).gt.1-1d-4) l1=l1+1
            xx(ne(nelb(k,mcy(l1,kcy)),neb(k,mcy(l1,kcy))))=xc(2)
            yy(ne(nelb(k,mcy(l1,kcy)),neb(k,mcy(l1,kcy))))=yc(2)
          else
           j2=jc(2)
           l1=lc(1)
          endif
c       call pltb(ne,xx,yy,xc(1),yc(1),xc(2),yc(2),ind,ine)
        write(*,*) 'j1,j2,l1,l2',j1,j2,l1,l2
c       pause
c
          jj=0
          do 250 j=jst,j1-1
          jj=jj+1
          kcop(jj)=ncop(i,j)
          keb(jj)=neb(i,j)
 250      kelb(jj)=nelb(i,j)
          if(l2.lt.l1) l2=l2+kcy
          do 260 l=l2,l1,-1
          jj=jj+1
          kcop(jj)=1
          ncop(k,l)=1
          keb(jj)=neb(k,mcy(l,kcy))
 260      kelb(jj)=nelb(k,mcy(l,kcy))
c         write(*,*) 'j2,nb(i)=',j2,nb(i)
          do 270 j=j2+1,jst+nb(i)-1
          jj=jj+1
          kcop(jj)=ncop(i,mcy(j,icy))
          keb(jj)=neb(i,mcy(j,icy))
 270      kelb(jj)=nelb(i,mcy(j,icy))
c
          nb(i)=jj
          kcop(jj+1)=kcop(1)
          keb(jj+1)=keb(1)
          kelb(jj+1)=kelb(1)
c
c               Adding/erasing common points of three region
          if(ifg(1).eq.1) then
            ncp=ncp+1
            kcp(ncp)=ne(nelb(k,l2),neb(k,l2))
            icp(ncp)=k
            jcp(ncp)=l2
c           xcp(ncp)=xx(kcp(ncp))
c           ycp(ncp)=yy(kcp(ncp))
            icl(ncp)=-999
          else
            do 272 j=1,ncp
 272        if(kcp(j).eq.ne(nelb(i,j1),neb(i,j1))) kcp(j)=-kcp(j)
          endif
          if(ifg(2).eq.1) then
            ncp=ncp+1
            kcp(ncp)=ne(nelb(k,l1),neb(k,l1))
            icp(ncp)=k
            jcp(ncp)=l1
c           xcp(ncp)=xx(kcp(ncp))
c           ycp(ncp)=yy(kcp(ncp))
            icl(ncp)=-999
          else
            do 274 j=1,ncp
 274        if(kcp(j).eq.ne(nelb(i,j2),neb(i,j2))) kcp(j)=-kcp(j)
          endif
c
          do 280 j=1,nb(i)+1
          ncop(i,j)=kcop(j)
          neb(i,j)=keb(j)
 280      nelb(i,j)=kelb(j)
c
          if(i.ne.k) then
            do 300 j=1,ncon(i)
 300        if(kcon(i,j).eq.k) goto 305
            ncon(i)=ncon(i)+1
            kcon(i,ncon(i))=k
            ncon(k)=ncon(k)+1
            kcon(k,ncon(k))=i
          endif
c
 305      do 310 j=1,nelt
            if(irx(j).ne.i) goto 310
            do 320 m=1,3
              n=ne(j,m)
              call checkcoal(nelt,ne,xx,yy,irx,n,k,inder,ar,epc,
     &                  ine,ind,iar)
              if(inder.ne.0) then
                irx(j)=0
                idx(j)=0
                goto 310
              endif
 320        continue
 310      continue
c
c       call pltb(ne,xx,yy,0d0,0d0,0d0,0d0,ind,ine)
c       pause
        goto 70
 60     continue
c
        if(icoal.eq.0) return
c
        call setar1
c
c       call pltb(ne,xx,yy,0d0,0d0,0d0,0d0,ind,ine)
c       pause
c
        return
        end
c
        function mcy(i,n)
        if(i.lt.1) then
          mcy=n-i
        elseif(i.gt.n) then
          mcy=i-n
        else
          mcy=i
        endif
        return
        end

        subroutine checkcoal(nel,ne,xx,yy,idx,n,idc,inder,ar,epc,
     &                  ine,ind,iar)
        implicit double precision (a-h,o-z)
        dimension xx(ind),yy(ind),ne(ine,3),idx(0:ine),ar(iar)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        integer*2 jp1(3),jp2(3)
        data jp1/2,3,1/,jp2/3,1,2/
          x=xx(n)
          y=yy(n)
          do 100 j=1,nel
          if(idc.gt.-900 .and. idx(j).ne.idc) goto 100
          if(ne(j,1).eq.n .or. ne(j,2).eq.n .or. ne(j,3).eq.n) goto 100
          do 120 l=1,3
            al=((xx(ne(j,jp1(l)))-x)*(yy(ne(j,jp2(l)))-y)
     &         -(xx(ne(j,jp2(l)))-x)*(yy(ne(j,jp1(l)))-y))/ar(j)
 120      if(al.lt.-epc .or. al.gt.1d0+epc) goto 100
          inder=idx(j)
          return
 100    continue
        inder=0
        return
        end
c
        subroutine crss(xi1,yi1,xi2,yi2,xk1,yk1,xk2,yk2,
     &                                          xs,ys,si,sk,inder,ep)
        implicit double precision (a-h,o-z)
        inder=0
        a2=yi2-yi1
        b2=-(xi2-xi1)
        c2=a2*xi1+b2*yi1
        if((xk2.eq.xi1.and.yk2.eq.yi1).or.
     &                  (xk2.eq.xi2.and.yk2.eq.yi2)) then
          xs=xk2
          ys=yk2
          goto 90
        endif
c
c       if((dabs(xk2-xi1).lt.1d-12 .and. dabs(yk2-yi1).lt.1d-12) .or.
c     &    (dabs(xk2-xi2).lt.1d-12 .and. dabs(yk2-yi2).lt.1d-12)) then
c          xs=xk2
c          ys=yk2
c          goto 90
c       endif
c 
        a1=yk2-yk1
        b1=-(xk2-xk1)
        c1=a1*xk1+b1*yk1
        if(a1*b2-a2*b1.eq.0d0) return
        xs=(c1*b2-c2*b1)/(a1*b2-a2*b1)
        if(xi1.eq.xi2 .and. (xs-xk1)*(xk2-xs).ge.ep) goto 20
        if(xk1.eq.xk2 .and. (xs-xi1)*(xi2-xs).ge.ep) goto 20
        if((xs-xi1)*(xi2-xs).lt.ep.or.(xs-xk1)*(xk2-xs).lt.ep) return
 20     ys=(c1*a2-c2*a1)/(a2*b1-a1*b2)
        if(yi1.eq.yi2 .and. (ys-yk1)*(yk2-ys).ge.ep) goto 90
        if(yk1.eq.yk2 .and. (ys-yi1)*(yi2-ys).ge.ep) goto 90
        if((ys-yi1)*(yi2-ys).lt.ep.or.(ys-yk1)*(yk2-ys).lt.ep) return
 90     inder=1
        si=dsqrt((xs-xi1)**2+(ys-yi1)**2)
     &                  /dsqrt((xi2-xi1)**2+(yi2-yi1)**2)
        sk=dsqrt((xs-xk1)**2+(ys-yk1)**2)
     &                  /dsqrt((xk2-xk1)**2+(yk2-yk1)**2)
        return
        end
