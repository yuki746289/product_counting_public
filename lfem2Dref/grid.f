         subroutine grd(nel,np,ne,xx,yy,idx,irx,jnt,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          new,xw,yw,idw,irw,jnw,
     &          nelrs,nelre,npc,nbwc,mcr,ncr,kcr,
     &          ar0,ind,ine,ibw,ica)


        implicit double precision (a-h,o-z)
!　ncomをしゅうきかんすうにするかどうか？
        parameter (ixg=301,iyg=501,ixy=151604)          !ixy=ixg*iyg
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine),irx(0:ine)
        dimension nen(ine,3),xn(ind),yn(ind),idn(0:ine),irn(0:ine)
        dimension new(ine,3),xw(ind),yw(ind),idw(0:ine),irw(0:ine)
        dimension jnt(ica,5,ind),jnw(ind),ar0(ine)
        dimension mcr(ica),ncr(ica,5),kcr(ica,5,21)
        dimension npc(ica,5),nbwc(ica,5),nelrs(21),nelre(21)
        dimension xg(ixg,iyg),yg(ixg,iyg),xgi(ixg),ygi(iyg)
        dimension kxs(21),kys(21),kcal(21)
        dimension xxb(1001),yyb(1001),nncop(0:1001)
        integer*2 kumn(ixy),kux(ixy),kuy(ixy),kuf(ixy)
        integer*2 mg(ixg,iyg),mk(ixg,iyg),kum(ixg,iyg),mf(ixg,iyg,2)
        integer*2 mmg(21,ixg,iyg),js2(3)
c        integer*2 kumnr(ixy),kufn(ixy)     !
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /sufcood/xb(21,1001),yb(21,1001)
        common /regcon/ncon(21),kcon(21,21)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
        data js2/2,3,1/
        eta=.55
c---- arが負の時その接点をskipするが。jcpをずらさなくてよいか？
        write(*,*)"in grd"      ! 
        do 10 i=1,nr
          jj=0
          do 20 j=1,nb(i)
            n=nelb(i,j)
            ar=(xx(ne(n,2))-xx(ne(n,1)))*(yy(ne(n,3))-yy(ne(n,1)))
     &        -(xx(ne(n,3))-xx(ne(n,1)))*(yy(ne(n,2))-yy(ne(n,1)))
            if(ar.le.0) goto 20
            jj=jj+1
            xb(i,jj)=xx(ne(nelb(i,j),neb(i,j)))
            yb(i,jj)=yy(ne(nelb(i,j),neb(i,j)))
 20       continue
          nb(i)=jj
          xb(i,jj+1)=xb(i,1)
 10     yb(i,jj+1)=yb(i,1)
        open(1,file="xb.res")   !
        rewind(1)
        do i=1,nr
        write(1,*)"nr=",i
        do j=1,nb(i)
        write(1,*)"xb(",i,",",j,")=",xb(i,j),
     &            "yb(",i,",",j,")=",yb(i,j)
        enddo
        enddo
        close(1)               !
        nset0=nset
        igerr=0
c
        write(*,*)"mesh started in grd"         !
 90     call mesh(nxg,nyg,xg,yg,xgi,ygi,kxs,kys,mg,xx,yy,
     &          ind,igerr,ixg,iyg)
        write(*,*) 'nset,nel,nxg,nyg=',nset,nel,nxg,nyg
        write(*,*)"mesh ended in grd"          !
c               call pltm(nxg,nyg,xg,yg,ixg,iyg)
c               pause
c
        neln=0
        npn=0
        num=0
        do 100 j=1,nyg
        do 100 i=1,nxg
        num=num+1
        kum(i,j)=num
        kux(num)=i
 100    kuy(num)=j
c        do i=1,num    !
c          kufn(i)=0
c        enddo
c
c       write(*,*) 'nset,nr=',nset,nr,kset(1),kset(2),kset(3)
c        open(1,file="n.res")      !
c        rewind(1)
        if(kset(1).ne.kset(2))then   !
          write(*,*)"kset(1)=",kset(1),"kset(2)=",kset(2)
          pause
          endif                                           !
        do 200 ls=1,nset
c               write(*,*) 'ls=',ls
          do 220 i=1,nxg
          do 220 j=1,nyg
          mg(i,j)=0
          mk(i,j)=0
          do 220 k=1,2
 220      mf(i,j,k)=0
          do 240 i=1,num
 240      kuf(i)=0

          igerr=0
          do 320 lr=1,nr
            if(kset(lr).ne.ls) goto 320
            nbp=nb(lr)+1
            do 330 i=1,nbp
            xxb(i)=xb(lr,i)
            yyb(i)=yb(lr,i)
 330        nncop(i)=ncop(lr,i)
            nncop(0)=ncop(lr,nb(lr))
            nncop(nbp)=ncop(lr,1)
            do 340 i=1,nxg
            do 340 j=1,nyg
 340        mg(i,j)=-iabs(mg(i,j))
            ii0=kxs(lr)
            jj0=kys(lr)
            idc=ndx(lr)
            irc=lr
            write(*,*)"adjust starter in grd irc=",irc,"ii0=",ii0,
     &               "jj0=",jj0   !
            call adjust(nxg,nyg,xg,yg,mk,mg,nncop,nbp,xxb,yyb,
     &          ii0,jj0,idc,ixg,iyg,igerr)
c            call adjust2(nxg,nyg,xg,yg,mk,mg,nncop,nbp,xxb,yyb,    !mk(ii,jj)が求まらない
c     &          ii0,jj0,idc,ixg,iyg,igerr)
c            call adjust3(nxg,nyg,xg,yg,mk,mg,nncop,nbp,xxb,yyb,
c     &          ii0,jj0,idc,ixg,iyg,igerr)

            open(1,file="mg.res")
	      rewind(1)
	      do i=1,nxg
	      do j=1,nyg
              if(mg(i,j).eq.1)then
			write(1,*)i,j,sngl(xg(i,j)),sngl(yg(i,j))
			endif
		  enddo
		  enddo
		  close(1)	      

            do 350 i=1,nxg
            do 350 j=1,nyg
 350        mmg(lr,i,j)=mg(i,j)
c
            if(igerr.ge.10) then
                  write(*,*) ' err/grd igerr=',igerr
                  stop
            elseif(igerr.ne.0) then
                  write(*,*) 'grd/adding a new mesh/igerr=',igerr
c                       call pltm(nxg,nyg,xg,yg,ixg,iyg)
c                       pause
                  goto 90
             endif
c
c               call pltg(nbp,xxb,yyb,nxg,nyg,xg,yg,mg,ixg,iyg)
c               pause
 320        continue
c            call pltg2(nxg,nyg,xg,yg,mmg,ixg,iyg)
c            pause
c
            nelw=0
            nelw0=0
            do 360 lr=1,nr
              write(*,*)"kset(",lr,")=",kset(lr),"ls=",ls    !
              if(kset(lr).ne.ls) goto 360
              do 380 i=1,nxg
              do 380 j=1,nyg
 380          mg(i,j)=mmg(lr,i,j)
              idc=ndx(lr)
              irc=lr
              eta=etas(idc)
              call ardiv(nel,ne,xx,yy,irx,nelw,new,idw,irw,idc,irc,
     &                  eta,nxg,nyg,xg,yg,mg,mk,kuf,kum,mf,
     &                  ine,ind,ixg,iyg,ixy)
              write(*,*)"ardiv ended lr=",lr     !
c              write(1,*)"nr=",lr  !
c              do knum=1,num       !
c              if(kuf(knum).eq.1)
c     &        write(1,*)"kum(",kux(knum),",",kuy(knum),")=",knum
c              enddo                !
              nelrs(lr)=nelw0+1+neln
              nelre(lr)=nelw+neln
              nelw0=nelw
              write(*,*) 'ls,lr,nelr(e-s)',ls,lr,nelrs(lr),nelre(lr)
              write(*,*)"nelw0=",nelw0,"nelw=",nelw,"neln=",
     &                                  neln,"npn=",npn   !
 360        continue
c
c            npw=0             ! 
c            do 400 i=1,num    !each region
c              if(kuf(i).eq.1) then
c                if(kufn(i).eq.0)then
c                  npw=npw+1
c                  kumn(i)=npw+npn     !each region
c                  kumnr(i)=npw+npn    !all region
c                  xw(npw+npn)=xg(kux(i),kuy(i))
c                  yw(npw+npn)=yg(kux(i),kuy(i))
c                  kuf(i)=-1       !each region
c                  kufn(i)=-1      !all region
c                elseif(kufn(i).eq.-1)then
c                  kumn(i)=kumnr(i)
c                  xw(kumn(i))=xg(kux(i),kuy(i))
c                  yw(kumn(i))=yg(kux(i),kuy(i))
c                  kuf(i)=-1
c                  write(*,*)"lr=",lr,"kumn(",i,")=",kumn(i)  !
c                endif
c              endif
c 400        continue
cc
c            do 420 i=1,nelw
c            idn(i+neln)=idw(i)
c            irn(i+neln)=irw(i)
c            do 425 j=1,3
c 425        nen(i+neln,j)=kumn(new(i,j))
c c420    continue
cc
c          do 440 i=npn+1,npn+npw
c          xn(i)=xw(i)
c 440      yn(i)=yw(i)
c
c          do 460 i=1,nelw
c          idn(i+neln)=idw(i)
c          irn(i+neln)=irw(i)
c          do 460 j=1,3
c          if(kufm(new(i,j)).eq.1)then
c            nen(i+neln,j)=new(i,j)+npn   !
c          else
c            write(*,*)"kufm(",new(i,j),")=",kufm(new(i,j))
c            pause
c          endif
c 460      continue
c
          npw=0
          do 400 i=1,num
              if(kuf(i).eq.1)then
              npw=npw+1
                        kumn(i)=npw
                        xw(npw)=xg(kux(i),kuy(i))
                        yw(npw)=yg(kux(i),kuy(i))
                        kuf(i)=-1
                  endif
400       continue
c
          do 420 i=1,nelw
            do 425 j=1,3
425       new(i,j)=kumn(new(i,j))
420       continue
c
          do 440 i=1,npw
            xn(i+npn)=xw(i)
440       yn(i+npn)=yw(i)
c
          do 460 i=1,nelw
            idn(i+neln)=idw(i)
            irn(i+neln)=irw(i)
          do 460 j=1,3
460       nen(i+neln,j)=new(i,j)+npn
c
          npn=npn+npw
          neln=neln+nelw
c       write(*,*) 'ls,npn,neln=',ls,npn,nelw
c       call plt(nelw,new,xw,yw,ind,ine)
c       pause
c
 200    continue
c        close(1)  !
        write(*,*)"loop 200 ended"   !
        open(1,file="nen.res")
        rewind(1)
        do i=1,neln
        do j=1,3
          write(1,*)"nen(",i,",",j,")=",nen(i,j)
        enddo
        enddo
        close(1)
        if(npn.gt.ind .or. neln.gt.ine) then
          write(*,*) 'grd/ increase ind or ine'
          write(*,*) '    npn*3,neln=',npn*3,neln
          stop
        endif
        write(*,*)"before thinmesh"
c        call plt(neln,nen,xn,yn,ind,ine)
c        pause
c---------------------------
         write(*,*)"thinmesh started in ger neln=",neln,"nr=",nr   !
         call thinmesh2(nen,npn,neln,nr,nelrs,nelre,xn,yn,idn,irn,
     &                       new,jnw,ind,ine)
c----------------------------
        write(*,*)"sufele started in grd nr=",nr     !
        call sufele(neln,nen,xn,yn,idn,xx,yy,ine,ind)
        write(*,*)"sufele ended in grd nr=",nr       !
        do 660 i=1,neln
        ar0(i)=
     &  (xn(nen(i,2))-xn(nen(i,1)))*(yn(nen(i,3))-yn(nen(i,1)))
     &  -(xn(nen(i,3))-xn(nen(i,1)))*(yn(nen(i,2))-yn(nen(i,1)))
        if(ar0(i).le.0) then
          write(*,*) 'grd/ar0(i)<0)'
            write(*,*)"ar0(",i,")=",ar0(i)
            write(*,*)"xn=",xn(nen(i,1)),xn(nen(i,2)),xn(nen(i,3))
            write(*,*)"yn=",yn(nen(i,1)),yn(nen(i,2)),yn(nen(i,3))
          stop
        endif
 660    continue
c----------------------------
        do 700 m=1,msw
        write(*,*) m
          call calreg(m-1,ncal,kcal)
          do 720 i=1,ixy
 720      kux(i)=0
          mcr(m)=0
          nelw=0
          npw=0
          npwold=0
          do 730 n=1,ncal
            ncr(m,n)=0
            do 740 i=1,nr
              if(kcal(i).ne.n) goto 740
              ii=ii+1
              mcr(m)=max(mcr(m),kcal(i))
              ncr(m,n)=ncr(m,n)+1
              kcr(m,n,ncr(m,n))=i
              write(*,*) 'm,n,i,nelr=',m,n,i,nelrs(i),nelre(i)
              do 750 j=nelrs(i),nelre(i)
                nelw=nelw+1
                do 750 k=1,3
                if(nen(j,k).eq.0) stop
                if(kux(nen(j,k)).eq.0) then
                  npw=npw+1
                  kux(nen(j,k))=npw
                  kuy(npw)=nen(j,k)
                endif
 750          new(nelw,k)=kux(nen(j,k))
 740        continue
c           do 780 i=1,ind
c 780       jnw(i)=0
            write(*,*)"renum started in grd"    !
            call renum(nelw,new,npw,jnw,ind,ine,nbwc(m,n))
c            call renum2(1,npn,nen,neln,jww,ine,ind,nbwc(m,n))
            write(*,*)"renum ended in grd"      !
            do 760 j=1,npw
 760        jnt(m,n,kuy(j))=jnw(j)
c 760       jnt(m,n,kuy(j)+npwold)=jnw(j)+npwold
c           npwold=npw
c       write(*,*) m,n,npw
            npc(m,n)=npw
                write(*,*) 'm,n,np,nel,nbw=',m,n,npw,nelw,nbwc(m,n)
                if(npw*3.gt.ind .or. nbwc(m,n)*6.gt.ibw) then
                  write(*,*) 'grd/ npw*3,nbw*6=',npw*3,nbwc(m,n)*6
                  stop
                endif
 730      continue
 700    continue
       open(1,file="xygrd.res")
         rewind(1)
         do i=1,npw
         write(1,*)"xn(",i,")=",xn(i),"yn(",i,")=",yn(i)
         enddo
         close(1)
         open(1,file="negrd.res")
         rewind(1)
         do i=1,nelw
         do j=1,3
         write(1,*)"new(",i,",",j,")=",new(i,j)
         enddo
         enddo
         close(1)
c               write(*,*) 'npw,nelw=',npw,nelw
c
c       call plt(neln,nen,xn,yn,ind,ine)
c       pause
c
        return
          do 9230 m=1,mcr(1)
            do 9240 i=1,ind
 9240       jnw(i)=jnt(1,m,i)
            ncrg=ncr(1,m)
            nelw=0
            do 9250 n=1,ncrg
              kcrg=kcr(1,m,n)
c       write(*,*) 'm,n,i,nelr=',m,n,kcrg,nelrs(kcrg),nelre(kcrg)
              do 9260 i=nelrs(kcrg),nelre(kcrg)
              nelw=nelw+1
              do 9260 j=1,3
              new(nelw,j)=jnw(nen(i,j))
              xw(new(nelw,j))=xn(nen(i,j))
              yw(new(nelw,j))=yn(nen(i,j))
 9260         continue
 9250       continue
cc        call plt(nelw,new,xw,yw,ind,ine)
        pause
 9230   continue
c
        return
      end
c
        subroutine adjust(nxg,nyg,xg,yg,mk,mg,ncop,nbp,xb,yb,ii0,jj0
     &                  ,idc,ixg,iyg,igerr)
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension xb(1001),yb(1001),sb(1001),ncop(0:1001)
        dimension ssp(1001)
        integer*2 mk(ixg,iyg),mg(ixg,iyg)
        integer*2 isp(1001),jsp(1001)
        integer*2 is0(8),is1(8),js0(8),js1(8),it(8),jt(8),kt(8)
        common /addmesh/ i0,j0,i1,j1
c       data is0/1,0,0,-1,-1, 0, 0, 1/, is1/1,1,-1,-1,-1,-1, 1, 1/
c       data js0/0,1,1, 0, 0,-1,-1, 0/ ,js1/1,1, 1, 1,-1,-1,-1,-1/
        data is0/1,1,0,-1,-1, -1, 0, 1/, is1/1,0,-1,-1,-1,0, 1, 1/
        data js0/0,1,1, 1, 0,-1,-1, -1/ ,js1/1,1, 1, 0,-1,-1,-1,0/
        data it/0,0,-1,-1,-1,-1,0,0/,jt/0,0,0,0,-1,-1,-1,-1/
        data kt/1,1,2,2,1,1,2,2/
c-------------
        write(*,*)"in adjust ii0=",ii0,"jj0=",jj0     !
c       write(*,*) 'adjust'
c       write(*,*) nxg,nyg,nbp,ii0,jj0,ixg,iyg
c       stop
c
        sb(1)=0
        dsav=0
        do 100 i=2,nbp
        sb(i)=sb(i-1)+dsqrt((xb(i)-xb(i-1))**2+(yb(i)-yb(i-1))**2)
 100    dsav=dsav+(sb(i)-sb(i-1))
        s0=sb(1)
        ib0=1
        dsav=dsav/dble(nbp-1)
c       dsav=1e-4
c       write(*,*) dsav
c
        i0=0
        j0=0
c-------------
        iist=ii0
        jjst=jj0
        ii=ii0
        jj=jj0
        mg(ii,jj)=1
        ksp=1
        ssp(ksp)=s
        isp(ksp)=ii
        jsp(ksp)=jj
cc
c       write(*,*) 'iist,jjst=',iist,jjst
        write(*,*)"ii=",ii,"jj=",jj,"iist=",iist,"jjst=",jjst     !
c-------------
c
        ep=1d-12
        do 210 lp=1,10000

          sm=9999
          do 220 lqq=ib0,nbp-1,6
          do 230 lq=lqq,min(lqq+5,nbp-1)
          do 230 k=1,8
            i0=ii+is0(k)
            j0=jj+js0(k)
            i1=ii+is1(k)
            j1=jj+js1(k)
            if(i0.lt.1 .or. j0.lt.1) goto 230
            if(i1.lt.1 .or. j1.lt.1) goto 230
            if(i0.gt.nxg .or. j0.gt.nyg) goto 230
            if(i1.gt.nxg .or. j1.gt.nyg) goto 230
            xg0=xg(i0,j0)
            yg0=yg(i0,j0)
            xg1=xg(i1,j1)
            yg1=yg(i1,j1)
            ib=lq
            call cross(xb,yb,xg0,yg0,xg1,yg1,xs,ys,s,ib,-ep,sb,inder)
c       if(ib.eq.37.and.inder.ne.0) then
c                write(*,*) i0,j0,i1,j1,sngl(s),sngl(sm)
c                write(*,*) xg0,yg0
c                write(*,*) xg1,yg1
c                write(*,*) xs,ys
c       endif
            if(inder.eq.0 .or. s.le.s0+dsav*1d-4 .or. s.gt.sm) goto 230
c            if(inder.eq.0 .or. s.le.s0+1d-12 .or. s.gt.sm) goto 230
            sm=s
            km=k
            xsm=xs
            ysm=ys
            lqm=lq
 230      continue
          if(sm.lt.9000) goto 240
 220      continue
                write(*,*) '##adjust/error ib0,ii,jj=',ib0,ii,jj
                write(*,*)"lp=",lp,"xg(",ii,",",jj,")=",xg(ii,jj),    !
     &                             "yg(",ii,",",jj,")=",yg(ii,jj),         
     &            "xg(",ii+1,",",jj,")=",xg(ii+1,jj),
     &            "yg(",ii+1,",",jj,")=",yg(ii+1,jj),
     &            "xg(",ii+1,",",jj+1,")=",xg(ii+1,jj+1),
     &            "yg(",ii+1,",",jj+1,")=",yg(ii+1,jj+1),
     &            "xg(",ii,",",jj+1,")=",xg(ii,jj+1),
     &            "yg(",ii,",",jj+1,")=",yg(ii,jj+1),
     &            "xg(",ii-1,",",jj+1,")=",xg(ii-1,jj+1),
     &            "yg(",ii-1,",",jj+1,")=",yg(ii-1,jj+1),
     &            "xg(",ii-1,",",jj,")=",xg(ii-1,jj),
     &            "yg(",ii-1,",",jj,")=",yg(ii-1,jj),
     &            "xg(",ii-1,",",jj-1,")=",xg(ii-1,jj-1),
     &            "yg(",ii-1,",",jj-1,")=",yg(ii-1,jj-1),
     &            "xg(",ii,",",jj-1,")=",xg(ii,jj-1),
     &            "yg(",ii,",",jj-1,")=",yg(ii,jj-1),
     &            "xg(",ii+1,",",jj-1,")=",xg(ii+1,jj-1),
     &            "yg(",ii+1,",",jj-1,")=",yg(ii+1,jj-1)
cc                call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
                stop
c -----
 240      continue
          s=sm
          xs=xsm
          ys=ysm
          k=km
          ib=lqm
          i0=ii+is0(k)
          j0=jj+js0(k)
          i1=ii+is1(k)
          j1=jj+js1(k)
          xg0=xg(i0,j0)
          yg0=yg(i0,j0)
          xg1=xg(i1,j1)
          yg1=yg(i1,j1)
          dg=dmax1(dabs(xg1-xg0),dabs(yg1-yg0))
c -----
          if(dabs(s-sb(nbp)).lt.ep) then                !end point
            is=iist
            js=jjst
            if((ii-is)*(jj-js).ne.0.and.mk(ii+it(k),jj+jt(k)).eq.0)
     &                               mk(ii+it(k),jj+jt(k))=kt(k)
            goto 290
          endif
          is=i1
          js=j1
          if(dabs(xs-xg0).lt.ep.and. dabs(ys-yg0).lt.ep) then
            is=i0
                  js=j0
          elseif(((xs-xg0)**2+(ys-yg0)**2).le.
     &                          0.01*((xs-xg1)**2+(ys-yg1)**2))then
            is=i0
            js=j0
          endif

c       write(*,*) '*',is,js
          if(mg(is,js).eq.1) then
            if(is.eq.i0 .and. js.eq.j0) then
              is=i1
              js=j1
            else
              is=i0
              js=j0
            endif
            if(mg(is,js).eq.1) then
                write(*,*) 'adjust/(err 1) already used (is,js)=',is,js
                write(*,*) 'nb,ib,ii,jj',nbp-1,ib,ii,jj
                write(*,*) 's,sb(nbp)=',s,sb(nbp)
                write(*,*) xb(1),yb(1)
                write(*,*) xb(ib),yb(ib)
                write(*,*) xb(nbp),yb(nbp)
cc                call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
                igerr=igerr+1
                return
            endif
          endif
          if(mg(is,js).lt.0) then
c               write(*,*) "!",ii,jj,is,js,ncop(ib),ncop(ib+1)
            if(ncop(ib).ne.0 .and.ncop(ib+1).ne.0) then
            elseif(dabs(xs-xb(ib)).lt.ep .and. dabs(ys-yb(ib)).lt.ep
     &                                  .and. ncop(ib).ne.0) then
            elseif(dabs(xs-xb(ib+1)).lt.ep .and. dabs(ys-yb(ib+1)).lt.ep
     &                                  .and. ncop(ib+1).ne.0) then
c               write(*,*) "$"
            else
c               write(*,*) "%"
              if(is.eq.i0 .and. js.eq.j0) then
                iss=i1
                jss=j1
              else
                iss=i0
                jss=j0
              endif
              if(mg(iss,jss).ne.-1) then
                is=iss
                js=jss
                if(mg(is,js).eq.1) then
                  write(*,*) 'adjust/(err 2)already used (is,js)=',is,js
                  write(*,*) 'ii,jj,i0,j0,i1,j1',ii,jj,i0,j0,i1,j1
cc                  call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
                  igerr=igerr+1
                  pause
                  return
                endif
              endif
            endif
c           write(*,*) "#",ii,jj,is,js
          endif
c
          if((ii-is)*(jj-js).ne.0.and.mk(ii+it(k),jj+jt(k)).eq.0)
     &                               mk(ii+it(k),jj+jt(k))=kt(k)
          xg(is,js)=xs
          yg(is,js)=ys
c         write(*,*) ib,is,js,mg(1,js),ncop(ib)
          mg(is,js)=1
          ii0=ii
          jj0=jj
          ii=is
          jj=js
          ib0=ib
          s0=s
c         call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
c         pause
          ksp=ksp+1
          ssp(ksp)=s
          isp(ksp)=ii
          jsp(ksp)=jj
 210  continue
c
 290    continue
c
c               call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
c               pause
c
        return
c
        isp(ksp+1)=isp(1)
        jsp(ksp+1)=jsp(1)
        isp(ksp+2)=isp(2)
        jsp(ksp+2)=jsp(2)
        isp(ksp+3)=isp(3)
        jsp(ksp+3)=jsp(3)
        do 320 k=1,ksp+1
          ax=xg(isp(k+1),jsp(k+1))-xg(isp(k),jsp(k))
          ay=yg(isp(k+1),jsp(k+1))-yg(isp(k),jsp(k))
          bx=xg(isp(k+2),jsp(k+2))-xg(isp(k+1),jsp(k+1))
          by=yg(isp(k+2),jsp(k+2))-yg(isp(k+1),jsp(k+1))
          cs=(ax*bx+ay*by)/dsqrt((ax*ax+ay*ay)*(bx*bx+by*by))
          if(cs.gt.0.99999999d0 .and. idc.gt.0) then
            xg(isp(k+1),jsp(k+1))=
     &                  (xg(isp(k),jsp(k))+xg(isp(k+2),jsp(k+2)))/2d0
            yg(isp(k+1),jsp(k+1))=
     &                  (yg(isp(k),jsp(k))+yg(isp(k+2),jsp(k+2)))/2d0
          endif
  320   continue
c       write(*,*) 'adjust end'
c       stop
c       write(*,*) 'idc=',idc
c       pause
c       if (iabs(idc).ne.2) pause                        ! use
c
c      do 295 j=1,jj
c 295  mg(1,j)=1
ccc      if(ifg.eq.0) goto 302 
ccc      if(lllppp.lt.2) then
ccc        do 300 i=1+9,ksp-5
ccc        inder=1
ccc        call suffc(ssp(i),fxs,fys,dmy1,dmy2,inder)
ccc        xg(isp(i),jsp(i))=fxs
ccc        yg(isp(i),jsp(i))=fys
ccc 300    continue
ccc        xtip=xx(ksf(nsf0))
ccc        call suffc(ssp(ksp),fxs,ytop,dmy1,dmy2,inder)
ccc      endif
 302  iia=ii
      jja=jj
c-------------------
      return
      end       
c
      subroutine cross(xb,yb,xg0,yg0,xg1,yg1,xs,ys,s,i,ep,sb,inder)
      implicit double precision (a-h,o-z)
      dimension xb(1001),yb(1001),sb(1001)
c
        inder=0         !交点は見つかっていない
        x0=xb(i)
        y0=yb(i)
        x1=xb(i+1)
        y1=yb(i+1)
c       if((x1.eq.xg0.and.y1.eq.yg0).or.
c     &                 (x1.eq.xg1.and.y1.eq.yg1)) then
        if((dabs(x1-xg0).lt.1d-12 .and. dabs(y1-yg0).lt.1d-12) .or.
     &     (dabs(x1-xg1).lt.1d-12 .and. dabs(y1-yg1).lt.1d-12)) then

          xs=x1
          ys=y1
          inder=1
c         write(*,*) 'cross'
          goto 90
        endif
c
        a2=yg1-yg0
        b2=-(xg1-xg0)
        c2=a2*xg0+b2*yg0
        a1=y1-y0
        b1=-(x1-x0)
        c1=a1*x0+b1*y0
        if(dabs(a1*b2-a2*b1).gt.1d-9) then
          xs=(c1*b2-c2*b1)/(a1*b2-a2*b1)
          ys=(c1*a2-c2*a1)/(a2*b1-a1*b2)
        else
          a2=xg1-xg0
          b2=-(yg1-yg0)
          c2=a2*yg0+b2*xg0
          a1=x1-x0
          b1=-(y1-y0)
          c1=a1*y0+b1*x0
          if(a1*b2-a2*b1.eq.0d0) return
          ys=(c1*b2-c2*b1)/(a1*b2-a2*b1)     !交点
          xs=(c1*a2-c2*a1)/(a2*b1-a1*b2)	   !交点
        endif
c
c       if(xg0.eq.xg1 .and. (xs-x0)*(x1-xs).ge.ep) goto 20
c       if(x0.eq.x1 .and. (xs-xg0)*(xg1-xs).ge.ep) goto 20
        if((xs-xg0)*(xg1-xs).lt.ep.or.(xs-x0)*(x1-xs).lt.ep) return
c
c 20    if(yg0.eq.yg1 .and. (ys-y0)*(y1-ys).ge.ep) goto 90
c       if(y0.eq.y1 .and. (ys-yg0)*(yg1-ys).ge.ep) goto 90
        if((ys-yg0)*(yg1-ys).lt.ep.or.(ys-y0)*(y1-ys).lt.ep) return
c
 90     s=sb(i)+dsqrt((xs-x0)**2+(ys-y0)**2)
c      ibn=i
c      write(*,*) '(cross)s',s
      inder=1
      return
      end
c
        subroutine meshadd(nxg,nyg,xg,yg,ixg,iyg)
        implicit double precision (a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        common /addmesh/ i0,j0,i1,j1
        if(i0.eq.i1) then
          do 100 i=1,nxg
            do 120 j=nyg,max(j0,j1),-1
            xg(i,j+1)=xg(i,j)
 120        yg(i,j+1)=yg(i,j)
            xg(i,max(j0,j1))=
     &          (xg(i,min(j0,j1))+xg(i,max(j0,j1)+1))/2d0
            yg(i,max(j0,j1))=
     &          (yg(i,min(j0,j1))+yg(i,max(j0,j1)+1))/2d0
 100      continue
          nyg=nyg+1
        else
          do 140 j=1,nyg
            do 160 i=nxg,max(i0,i1),-1
            xg(i+1,j)=xg(i,j)
 160        yg(i+1,j)=yg(i,j)
            xg(max(i0,i1),j)=
     &          (xg(min(i0,i1),j)+xg(max(i0,i1)+1,j))/2d0
            yg(max(i0,i1),j)=
     &          (yg(min(i0,i1),j)+yg(max(i0,i1)+1,j))/2d0
 140      continue
          nxg=nxg+1
        endif
        return
        end
c
      subroutine crossold(xb,yb,xg0,yg0,xg1,yg1,xs,ys,s,i,ep,sb,inder)
      implicit double precision (a-h,o-z)
      dimension xb(1001),yb(1001),sb(1001)
c
        inder=0         !交点は見つかっていない
        a2=yg1-yg0
        b2=-(xg1-xg0)
        c2=a2*xg0+b2*yg0
        x0=xb(i)
        y0=yb(i)
        x1=xb(i+1)
        y1=yb(i+1)
        if((x1.eq.xg0.and.y1.eq.yg0).or.
     &                  (x1.eq.xg1.and.y1.eq.yg1)) then
          xs=x1
          ys=y1
          goto 90
        endif
c
c       if((dabs(x1-xg0).lt.1d-12 .and. dabs(y1-yg0).lt.1d-12) .or.
c     &    (dabs(x1-xg1).lt.1d-12 .and. dabs(y1-yg1).lt.1d-12)) then
c          xs=x1
c          ys=y1
c          goto 90
c       endif
c 
        a1=y1-y0
        b1=-(x1-x0)
        c1=a1*x0+b1*y0
        if(a1*b2-a2*b1.eq.0d0) return
c
        xs=(c1*b2-c2*b1)/(a1*b2-a2*b1)
        if(xg0.eq.xg1 .and. (xs-x0)*(x1-xs).ge.ep) goto 20
        if(x0.eq.x1 .and. (xs-xg0)*(xg1-xs).ge.ep) goto 20
        if((xs-xg0)*(xg1-xs).lt.ep.or.(xs-x0)*(x1-xs).lt.ep) return
c
 20     ys=(c1*a2-c2*a1)/(a2*b1-a1*b2)
        if(yg0.eq.yg1 .and. (ys-y0)*(y1-ys).ge.ep) goto 90
        if(y0.eq.y1 .and. (ys-yg0)*(yg1-ys).ge.ep) goto 90
        if((ys-yg0)*(yg1-ys).lt.ep.or.(ys-y0)*(y1-ys).lt.ep) return
c
 90   s=sb(i)+dsqrt((xs-x0)**2+(ys-y0)**2)
c      ibn=i
c      write(*,*) '(cross)s',s
      inder=1
      return
      end
c
        subroutine ardiv(nelt,ne,xx,yy,irx,neln,nen,idn,irw,idc,irc,
     &          eta,nx,ny,xg,yg,mg,mk,kuf,kum,mf,ine,ind,ixg,iyg,ixy)
        implicit double precision (a-h,o-z)
        parameter (iar=8001)
        dimension ne(ine,3),nen(ine,3),xx(ind),yy(ind),irx(0:ine)
        dimension idn(0:ine),irw(0:ine),xg(ixg,iyg),yg(ixg,iyg)
        dimension xq(4),yq(4),as(2,2),ar(iar)      !,xy(501)   fcdで格子線を補正していない
        integer*2 kuf(ixy)
        integer*2 mg(ixg,iyg),mk(ixg,iyg)
        integer*2 kum(ixg,iyg),mf(ixg,iyg,2)
        integer*2 kdi(4),kdj(4),kdv(2,2,3)
        integer*2 kk(2),iz(5),jz(5),jp1(3),jp2(3)
        common /dgrd/ dx,dy,in0,in1,jn0,jn1
        common /gmovedata/ inmin,jnmin
        data kdi/0,1,1,0/, kdj/0,0,1,1/
        data iz/-1,1,0,0,0/,jz/0,0,-1,1,0/
        data (kdv(1,1,j),j=1,3)/1,3,4/,(kdv(1,2,j),j=1,3)/1,2,3/
        data (kdv(2,1,j),j=1,3)/1,2,4/,(kdv(2,2,j),j=1,3)/2,3,4/
        data jp1/2,3,1/,jp2/3,1,2/
        data epc /1d-15/
!ここでのエラーの多くは，adjustでのエラー
c
c----------
        write(*,*)"in ardiv nelt=",nelt     !
        open(1,file="irx.res")   !
        rewind(1)
        write(1,*)"irx=",irx,"nelt=",nelt,"in ardiv"
        do i=1,nelt
        write(1,*)"irx(",i,")=",irx(i)
        enddo
        close(1)
        if(iar.lt.nelt) then
                write(*,*) 'ardiv/increase iar > nel',iar,nelt
                stop
        endif
        do 8 k=1,nelt
          if(irx(k).ne.1 .and. irx(k).ne.2)then    !
            write(*,*)"irx(",k,")=",irx(k)
            pause
          endif
          if(irx(k).ne.irc) goto 8
          ar(k)=(xx(ne(k,2))-xx(ne(k,1)))*(yy(ne(k,3))-yy(ne(k,1)))
     &      -(xx(ne(k,3))-xx(ne(k,1)))*(yy(ne(k,2))-yy(ne(k,1)))
            if(ar(k).le.0d0) then
                write(*,*) 'ardiv / ar<0 / k,ar=',k,ar(k)
c               stop
            endif
 8      continue
c                       do 91 j=ny,1,-1
c 91                    write(*,'(i3,1x,60i1)') j,(mg(i,j),i=1,nx)
c                       write(*,*) 'mg',mg(1,1)
        iaa=0
        do 10 i=2,nx-1
        do 10 j=2,ny-1
c         if(mg(i,j).gt.0) goto 10
          if(mg(i,j).ne.0) goto 10
          x=xg(i,j)
          y=yg(i,j)
          do 12 k=1,nelt
            if(irx(k).ne.irc) goto 12
            do 14 l=1,3
              al=((xx(ne(k,jp1(l)))-x)*(yy(ne(k,jp2(l)))-y)
     &         -(xx(ne(k,jp2(l)))-x)*(yy(ne(k,jp1(l)))-y))/ar(k)
 14         if(al.lt.epc .or. al.gt.1d0-epc) goto 12
            iaa=iaa+1
            mg(i,j)=2
 12       continue
 10     continue
        write(*,*)"iaa=",iaa   !
c
        if(iaa.eq.0) then
          do 20 i=1,nx-1
          do 20 j=1,ny-1
            ll=0
            x=0d0
            y=0d0
            do 24 k=1,4
              if(mg(i+kdi(k),j+kdj(k)).eq.1) then
                ll=ll+1
                x=x+xg(i+kdi(k),j+kdj(k))/3d0
                y=y+yg(i+kdi(k),j+kdj(k))/3d0
              elseif(mg(i+kdi(k),j+kdj(k)).eq.0) then
                i0=i+kdi(k)
                j0=j+kdj(k)
              else
                ll=999
              endif
 24         continue
            if(ll.ne.3) goto 20
            do 26 k=1,nelt
              if(irx(k).ne.irc) goto 26
              do 28 l=1,3
                al=((xx(ne(k,jp1(l)))-x)*(yy(ne(k,jp2(l)))-y)
     &          -(xx(ne(k,jp2(l)))-x)*(yy(ne(k,jp1(l)))-y))/ar(k)
 28           if(al.lt.epc .or. al.gt.1d0-epc) goto 26
              mg(i0,j0)=2
              goto 30
 26         continue
 20       continue
          write(*,*) 'ardiv/no nod in region no.=',irc
          stop
c          stop   !
        endif
 30       continue
c                       write(*,*) '96'
c                       do 96 j=ny,1,-1
c 96                    write(*,'(60i1)') (mg(i,j),i=1,nx)
 40     do 45 i=2,nx-1
        do 50 j=2,ny-1
c         if(mg(i,j).gt.0) goto 50 
          if(mg(i,j).ne.0) goto 50 
          do 55  k=1,4
 55       if(mg(i+iz(k),j+jz(k)).eq.2) goto 60
          goto 50
 60       do 65  k=1,5
c 65      if(mg(i+iz(k),j+jz(k)).le.0) mg(i+iz(k),j+jz(k))=2
 65       if(mg(i+iz(k),j+jz(k)).eq.0) mg(i+iz(k),j+jz(k))=2
          goto 40
 50     continue
 45     continue
c       if(irc.eq.3) then
c                       write(*,*) '98'
c                       do 98 j=ny,1,-1
c 98                    write(*,'(i3,1x,60i1)') j,(mg(i,j),i=1,25)
c                       pause
c       endif
c  ------------------
        if(eta.lt.0) goto 280
c
        do 100 i=1,nx
          n1=0
          do 110 ii=1,11
            n1=n1+1
            ifg=0
            do 120 j=n1,ny
              if(ifg.eq.0) then
                if(mg(i,j).eq.1) ifg=1
              elseif(ifg.eq.1) then
                if(mg(i,j).eq.2) ifg=2
                n0=j-1
              elseif(ifg.eq.2) then
                if(mg(i,j).eq.1) goto 140
              endif
 120        continue
            goto 100
 140        n1=j
            if(n1.le.jnmin) goto 110
            n0=max(n0,jnmin)
            nm=(n0+n1)/2
            y0=yg(i,n0)
            ym=yg(i,nm)
            y1=yg(i,n1)
            k0=nm-n0+1
            k1=n1-nm+1
            k=n1-n0+1
c            if(k0.ge.3 .and. k1.ge.3) then
c              call fcd(k0,y0,(ym-y0)*(1d0-eta)+y0,ym,xy)
c              do 150 j=2,k0-1
c 150          yg(i,j+n0-1)=xy(j)
c              call fcd(k1,ym,(y1-ym)*eta+ym,y1,xy)
c              do 160 j=2,k1-1
c 160          yg(i,j+nm-1)=xy(j)
c           elseif(k.ge.3) then
c             call fcd(k,y0,(y1-y0)*0.5+y0,y1,xy)
c             do 170 j=2,k-1
c 170         yg(i,j+n0-1)=xy(j)
c            endif
 110      continue
 100    continue
c
        do 200 j=1,ny
          n1=0
          do 210 ii=1,11
            n1=n1+1
            ifg=0
            do 220 i=n1,nx
              if(ifg.eq.0) then
                if(mg(i,j).eq.1) ifg=1
              elseif(ifg.eq.1) then
                if(mg(i,j).eq.2) ifg=2
                n0=i-1
              elseif(ifg.eq.2) then
                if(mg(i,j).eq.1) goto 240
              endif
 220        continue
            goto 200
 240        n1=i
            if(n1.le.inmin) goto 210
            n0=max(n0,inmin)
            nm=(n0+n1)/2
            x0=xg(n0,j)
            xm=xg(nm,j)
            x1=xg(n1,j)
            k0=nm-n0+1
            k1=n1-nm+1
            k=n1-n0+1
c            if(k0.ge.3 .and. k1.ge.3) then
c              call fcd(k0,x0,(xm-x0)*(1d0-eta)+x0,xm,xy)
c              do 250 i=2,k0-1
c 250          xg(i+n0-1,j)=xy(i)
c              call fcd(k1,xm,(x1-xm)*eta+xm,x1,xy)
c              do 260 i=2,k1-1
c 260          xg(i+nm-1,j)=xy(i)
c           elseif(k.ge.3) then
c             call fcd(k,x0,(x1-x0)*0.5+x0,x1,xy)
c             do 270 i=2,k-1
c 270         xg(i+n0-1,j)=xy(i)
c            endif
 210      continue
 200    continue
c------------------
 280    continue
        do 400 jj=1,ny-1
        do 410 ii=1,nx-1
          if(mf(ii,jj,1).ne.0 .and. mf(ii,jj,2).ne.0) goto 410
          m=1
          k1=0
          k2=0
          do 420 i=1,4
          xq(i)=xg(ii+kdi(i),jj+kdj(i))
          yq(i)=yg(ii+kdi(i),jj+kdj(i))
          if(mg(ii+kdi(i),jj+kdj(i)).eq.1) k1=k1+1
 420      if(mg(ii+kdi(i),jj+kdj(i)).eq.2) k2=k2+1
          if((xq(1)-xq(3))**2+(yq(1)-yq(3))**2.gt.
     &            ((xq(2)-xq(4))**2+(yq(2)-yq(4))**2)*1.001) m=2
c--
          ars=0d0
          do 430 i=1,2
          do 435 j=1,3
          xq(j)=xg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))
 435      yq(j)=yg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))
 430      ars=ars+
     &       (xq(1)-xq(2))*(yq(1)-yq(3))-(xq(1)-xq(3))*(yq(1)-yq(2))
c
          if((k1+k2).eq.4) then
            m1=m
            m2=iabs(not(-(m1-1)))+1
            do 440 m=1,2
            do 440 i=1,2
            do 445 j=1,3
            xq(j)=xg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))
 445        yq(j)=yg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))
 440        as(m,i)=
     &        (xq(1)-xq(2))*(yq(1)-yq(3))-(xq(1)-xq(3))*(yq(1)-yq(2))
            if(as(m1,1)/ars.gt.0.01 .and. as(m1,2)/ars.gt.0.01) then
              m=m1
            elseif(as(m2,1)/ars.gt.0.01.and.as(m2,2)/ars.gt.0.01) then
              m=m2
            elseif(as(m1,1).gt.0d0 .and. as(m1,2).gt.0d0) then
              m=m1
            elseif(as(m2,1).gt.0d0 .and. as(m2,2).gt.0d0) then
              m=m2
            else
              write(*,*) 'err in ardiv/ii,jj,as='
     &                          ,ii,jj,((as(i,j)/ars,j=1,2),i=1,2)
              write(*,*) 'k1,k2=',k1,k2
                x=xg(ii,jj)
                y=yg(ii,jj)
cc                call pltardiv(nelt,ne,xx,yy,irx,irc,x,y,ind,ine,
cc     &                  nx,ny,xg,yg,ixg,iyg)
              stop
            endif
          endif
c
          kk(1)=0
          kk(2)=0
          if(k1.eq.0)then
            if(k2.eq.4) then
              kk(1)=1
              kk(2)=1
            endif
          elseif(k1.eq.1)then
            if(k2.eq.3) then
              kk(1)=1
              kk(2)=1
            endif
          elseif(k1.eq.2)then
            if(k2.eq.1)then
              m=mk(ii,jj)
                if(m.eq.0) then
                 write(*,*) 'err in ardiv/ii,jj,k1,k2,m=',ii,jj,k1,k2,m
                 stop
                endif
              do 450 i=1,2
              do 450 j=1,3
 450          if(mg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j))).eq.2)
     &                                                     goto 455
 455          kk(i)=1
            elseif(k2.eq.2)then
              kk(1)=1
              kk(2)=1
            endif
          elseif(k1.eq.3)then
            if(k2.eq.0)then
              m=mk(ii,jj)
              if(m.ne.0)then
                kk(1)=1
                kk(2)=1
                do 460 i=1,2
                do 460 j=1,3
 460            if(mg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j))).le.0) iv=i
                kk(iv)=0
                iv=iabs(not(-(iv-1)))+1              !special case
                if(mf(ii,jj,iv).ne.0) kk(iv)=0       !
              endif
            else                        !k2=1
              if(mk(ii,jj).eq.0) then
                kk(1)=1
                kk(2)=1
              else
                m=mk(ii,jj)
                kk(1)=0
                kk(2)=0
                do 467 i=1,2
                do 467 j=1,3
 467            if(mg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j))).eq.2) iv=i
                kk(iv)=1
                iv=iabs(not(-(iv-1)))+1
                if(mf(ii,jj,iv).ne.0) then
                  kk(iv)=0
                else
                  x=0d0
                  y=0d0
                  do 468 j=1,3
                  x=x+xg(ii+kdi(kdv(m,iv,j)),jj+kdj(kdv(m,iv,j)))/3d0
 468              y=y+yg(ii+kdi(kdv(m,iv,j)),jj+kdj(kdv(m,iv,j)))/3d0
                  call chekin(nelt,ne,xx,yy,irx,ar,irc,x,y,inder,
     &                      ine,ind,iar)
                  if(inder.eq.1) kk(iv)=1
                endif
              endif
            endif
          elseif(k1.eq.4)then
            iv=0
            as(1,1)=0d0
            as(1,2)=0d0
            do 470 i=1,2
            x=0d0
            y=0d0
            do 475 j=1,3
            x=x+xg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))/3d0
 475        y=y+yg(ii+kdi(kdv(m,i,j)),jj+kdj(kdv(m,i,j)))/3d0
            call chekin(nelt,ne,xx,yy,irx,ar,irc,x,y,inder,ine,ind,iar)
c               call pltardiv(nelt,ne,xx,yy,irx,irc,x,y,ind,ine,
c     &                 nx,ny,xg,yg,ixg,iyg)
c       pause
            if(inder.eq.1) then
              iv=iv+1
              as(1,i)=1
            endif
 470        continue
            if(mk(ii,jj).eq.0) then
              if(iv.eq.2) then
                kk(1)=1
                kk(2)=1
              elseif(iv.eq.1)then
                if (mf(ii,jj,1).eq.0 .and. mf(ii,jj,2).eq.0) then
                  kk(1)=1
                  kk(2)=1
                  write(*,*) 'ardiv(iv=1)/kk(1)=kk(2)=1/ii,jj=',ii,jj
                else
c                 x=0d0
c                 y=0d0
c                 do 478 j=1,4
c                 x=x+xg(ii+kdi(j),jj+kdj(j))/4d0
c478              y=y+yg(ii+kdi(j),jj+kdj(j))/4d0
c                 call chekin(nelt,ne,xx,yy,irx,ar,irc,x,y,
c     &                                         inder,ine,ind,iar)
c                 if(inder.eq.1) then
c                   kk(1)=1
c                   kk(2)=1
c                   write(*,*) 'ardiv(sq)/kk(1)=kk(2)=1/ii,jj=',ii,jj
c                 else
                    kk(1)=0
                    kk(2)=0
                    write(*,*) 'error!/ardiv ... iv=1',iv,ii,jj
                    stop
c                 endif
                endif
              endif
            else
              m=mk(ii,jj)
              if(mf(ii,jj,1).ne.0) then
                kk(1)=0
                kk(2)=1
              elseif(mf(ii,jj,2).ne.0)then
                kk(1)=1
                kk(2)=0
              else
                kk(1)=int(as(1,1)+.01)
                kk(2)=int(as(1,2)+.01)
              endif
            endif
          else
            write(*,*) 'err/ardiv ii,jj,k1,k2=',ii,jj,k1,k2
            stop
          endif
c
          if(kk(1).eq.1 .or. kk(2).eq.1) then
            if(mf(ii,jj,1).ne.0 .and. m.ne.mf(ii,jj,1)) then
              write(*,*) 'ardiv/1:ii,jj,m,mf=',ii,jj,m,mf(ii,jj,1)
              call clos
              stop
            endif
            if(mf(ii,jj,2).ne.0 .and. m.ne.mf(ii,jj,2)) then
              write(*,*) 'ardiv/2:ii,jj,m,mf=',ii,jj,m,mf(ii,jj,2)
              call clos
              stop
            endif
          endif
c--------
c         if(kk(1).eq.0 .and. kk(2).eq.0) mk(ii,jj)=0
          if(kk(1).eq.1 .or. kk(2).eq.1) mk(ii,jj)=m
          do 500 k=1,2
            if(kk(k).eq.1 .and. mf(ii,jj,k).eq.0) then
              neln=neln+1
              idn(neln)=idc
              irw(neln)=irc
              mf(ii,jj,k)=m
              do 510 j=1,3
                ki=ii+kdi(kdv(m,k,j))
                kj=jj+kdj(kdv(m,k,j))
                n=kum(ki,kj)
                nen(neln,j)=n
                if(kuf(n).eq.0) kuf(n)=1
 510          continue
            endif
 500      continue
c--------
 410    continue
 400    continue
c
      return
      end
c
       subroutine chekin(nelt,ne,xx,yy,idx,ar,idc,x,y,inder,ine,ind,iar)
       implicit double precision (a-h,o-z)
       dimension xx(ind),yy(ind),ne(ine,3),idx(0:ine),ar(iar)
       integer*2 jp1(3),jp2(3)
       data jp1/2,3,1/,jp2/3,1,2/
        do 100 k=1,nelt
          if(idx(k).ne.idc .or. ar(k).le.0d0) goto 100
          do 120 l=1,3
            al=((xx(ne(k,jp1(l)))-x)*(yy(ne(k,jp2(l)))-y)
     &         -(xx(ne(k,jp2(l)))-x)*(yy(ne(k,jp1(l)))-y))/ar(k)
 120      if(al.lt.-1d-15 .or. al.gt.(1d0+1d-15)) goto 100
          inder=1
          return
 100    continue
        inder=0
        return
        end
c
        subroutine fcd(nrws,p1,p2,p3,xy)
        implicit double precision (a-h,o-z)
        dimension xy(501)
        deta=2d0/dble(nrws-1)
        do 10 i=1,nrws
          eta=dble(i-1)*deta-1d0
          z1=-.5d0*eta*(1d0-eta)
          z2=(1d0+eta)*(1d0-eta)
          z3=.5d0*eta*(1d0+eta)
          xy(i)=p1*z1+p2*z2+p3*z3
 10     continue
        return
        end
c
        subroutine sufele(nelt,ne,xx,yy,idx,x0,y0,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),idx(0:ine),xx(ind),yy(ind),x0(ind),y0(ind)
        dimension kelb(1001),keb(1001)
        dimension kcpn(101),js2(3)
        dimension icpn(101),jcpn(101),icln(101),stan(101),strn(101)
        integer*2 ichk(8001)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /sufcood/xb(21,1001),yb(21,1001)
        common /regcon/ ncon(21),kcon(21,21)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
        common /dnec/ nec(8001,3)
        data js2/2,3,1/
c
        write(*,*)"in sufele"      !
c       write(*,*) 'sufele0',nelt
                if(nelt.gt.8001) then
                write(*,*) 'sufele/increase nec/nelt=',nelt
                stop
                endif
c
        do 100 i=1,nelt
        do 100 j=1,3
 100    nec(i,j)=0
c
        write(*,*)"make nec(nelt,3) array"     !
        do 120 i=1,nelt
          do 130 j=1,nelt
            if (j.eq.i) goto 130
            do 140 k1=1,3
              if (nec(j,k1).ne.0) goto 140
              k2=js2(k1)
              do 150 jc=1,3
              if ((ne(i,jc).eq.ne(j,k2)).and.
     &                  (ne(i,js2(jc)).eq.ne(j,k1))) then
                 goto 190
              elseif ((ne(i,jc).eq.ne(j,k1)).and.          !
     &                  (ne(i,js2(jc)).eq.ne(j,k2))) then
                 write(*,*)"sufele/node number given on the",
     &                       " opposite direction in a element"
                   write(*,*)"ne(",i,",",jc,")","ne(",j,",",k1,")=",
     &                                        ne(i,jc)
                 write(*,*)"ne(",i,",",js2(jc),")","ne(",j,",",k2,")=",
     &                                        ne(j,k2)
                 write(*,*)"idx(",i,")=",idx(i),"idx(",j,")=",idx(j)
                 write(*,*)"xx(",ne(j,k1),")=",xx(ne(j,k1))
                           write(*,*)"xx(",ne(j,k2),")=",xx(ne(j,k2))
                 write(*,*)"yy(",ne(j,k1),")=",yy(ne(j,k1))
                           write(*,*)"yy(",ne(j,k2),")=",yy(ne(j,k2))
                 pause
              endif
 150          continue
 140        continue
            goto 130
 190        nec(i,jc)=j
            nec(j,k1)=i
            if(nec(i,1).ne.0.and.nec(i,2).ne.0.and.nec(i,3).ne.0)
     &          goto 120
 130      continue
 120     continue
         open(1,file="nec1.res")     !
         rewind(1)
         nnn=0
         do i=1,nelt
         do j=1,3
           nnn=nnn+1
           write(1,*)"nec(",i,",",j,")=",nec(i,j)
         enddo
         enddo
         write(1,*)"n=",nnn
         close(1)
c
c Searching Boundary Elements and Arrange the Element Data
c
c       write(*,*) 'sufele1'
        nn=0
        do 200 ld=-9,9
c         write(*,*) 'ld=',ld
          kb=0
          do 210 i=1,nelt
            if(idx(i).ne.ld) goto 210
            do 220 j=1,3
              if(idx(nec(i,j)).ne.ld) then
                kb=kb+1
                kelb(kb)=i
                keb(kb)=j
              endif
 220        continue
 210      continue
          if(kb+2.gt.1001) then
            write(*,*) 'sufele/increase kelb, nelb et al./kb=',kb
            stop
          endif
          if(kb.eq.0) goto 200
          open(1,file="kb.res")       !
          rewind(1)
          write(1,*)"kb=",kb
          do i=1,kb
          write(1,*)"i=",i,"ne(",kelb(i),",",keb(i),")="
     &                       ,ne(kelb(i),keb(i))
          enddo
          close(1)
          open(1,file="ne.res")       !
          rewind(1)
          write(1,*)"nelt=",nelt
          do i=1,nelt
          do j=1,3
          write(1,*)"ne(",i,",",j,")=",ne(i,j),"x=",xx(ne(i,j)),
     &                  "y=",yy(ne(i,j))
          enddo
          enddo
          close(1)
                ! Search start point: on origin > on y-axis > on xaxis
 230      imin=0
c         smin=9999
c         ifg=0
          xmin=9999
          ymin=9999
c
          do 232 i=1,8000
 232      ichk(i)=0
          do 234 i=1,kb
            if(kelb(i).lt.0) goto 234
            ichk(ne(kelb(i),keb(i)))=ichk(ne(kelb(i),keb(i)))+1
            if(ichk(ne(kelb(i),keb(i))).ge.2) then
              imin=i
             write(*,*) 'sufele/boundary surfaces are crossed/index=',ld
cc              call pltsuf(ne,xx,yy,kb,kelb,keb,ind,ine)
              stop
c             goto 260
            endif
 234      continue
c
          do 240 i=1,kb
            if(kelb(i).lt.0) goto 240
            x1=xx(ne(kelb(i),keb(i)))
            y1=yy(ne(kelb(i),keb(i)))
            x2=xx(ne(kelb(i),js2(keb(i))))
            y2=yy(ne(kelb(i),js2(keb(i))))
            if(x1.le.xmin .and. y1.le.ymin) then
              xmin=x1
              ymin=y1
              imin=i
            endif
c
c           if(dabs(x1).lt.1d-12 .and. dabs(y1).lt.1d-12) then
c             imin=i
c             goto 260
c           elseif(dabs(y1).lt.1d-12 .and. dabs(x1).lt.smin) then
c             ifg=1
c             smin=dabs(x1)
c             imin=i
c           elseif(dabs(x1).lt.1d-12 .and. dabs(y1).lt.smin) then
c             ifg=1
c             smin=dabs(y1)
c             imin=i
c           elseif(ifg.eq.0) then
c             s1=x1*x1+y1*y1
c             if(s1.lt.smin) then
c               smin=s1
c               imin=i
c             endif
c           endif
 240      continue
c
 260    if(imin.eq.0) goto 200
          nn=nn+1
          if(nn.gt.21) then
            write(*,*) 'sufele/increase nb, ndx et al./nn=',nn
            stop
          endif
          kk=1
          nelb(nn,kk)=kelb(imin)
          neb(nn,kk)=keb(imin)
          nst=ne(kelb(imin),keb(imin))
          n0=ne(kelb(imin),js2(keb(imin)))
          kelb(imin)=-kelb(imin)
c
 270      do 280 i=1,kb
            if(kelb(i).lt.0) goto 280
            if(n0.eq.ne(kelb(i),keb(i))) then
              kk=kk+1
              nelb(nn,kk)=kelb(i)
              neb(nn,kk)=keb(i)
              n0=ne(kelb(i),js2(keb(i)))
              kelb(i)=-kelb(i)
              if(n0.eq.nst) goto 290
              goto 270
            endif
 280      continue
                write(*,*) 'sufele/find the next element n0,nn=',n0,nn
                stop
 290      nb(nn)=kk
          ndx(nn)=ld
          nelb(nn,0)=nelb(nn,kk)
          nelb(nn,-1)=nelb(nn,kk-1)
          neb(nn,0)=neb(nn,kk)
          neb(nn,-1)=neb(nn,kk-1)
          nelb(nn,kk+1)=nelb(nn,1)
          nelb(nn,kk+2)=nelb(nn,2)
          neb(nn,kk+1)=neb(nn,1)
          neb(nn,kk+2)=neb(nn,2)
          open(1,file="nelb.res")     !
          rewind(1)
          write(1,*)"ld=",ld,"nn=",nn
          do i=1,kk
          inn=ne(nelb(nn,i),neb(nn,i))
          write(1,*)"ne(",nelb(nn,i),",",neb(nn,i),")=",
     &       ne(nelb(nn,i),neb(nn,i)),"x=",xx(inn),"y=",yy(inn)
          enddo
          close(1)
          goto 230     !imin=0となる
 200    continue
        nr=nn
        write(*,*)"in sufele nr=",nr     !
c
c       連続する領域の検索
c
        do 310 i=1,nr
 310    ncon(i)=0
        do 320 ii=1,nr
c               write(*,*) ii,nb(ii)
          do 330 jj=ii+1,nr
            do 340 i=1,nb(ii)
              n1=ne(nelb(ii,i),neb(ii,i))
              n2=ne(nelb(ii,i),js2(neb(ii,i)))
              do 350 j=1,nb(jj)
                k1=ne(nelb(jj,j),neb(jj,j))
                k2=ne(nelb(jj,j),js2(neb(jj,j)))
                if(n1.eq.k2 .and. n2.eq.k1) then
                    ncon(ii)=ncon(ii)+1
                    kcon(ii,ncon(ii))=jj
                    ncon(jj)=ncon(jj)+1
                    kcon(jj,ncon(jj))=ii
                    goto 330
                endif
c       if(dabs(xx(n1)-xx(k2)).lt.1d-6.and.dabs(yy(n1)-yy(k2)).lt.1d-6
c     &         .and.
c     & dabs(xx(n2)-xx(k1)).lt.1d-6.and.dabs(yy(n2)-yy(k1)).lt.1d-6)then
c       write(*,*) ii,jj,ndx(ii),ndx(jj)
c       endif
 350          continue
 340        continue
 330      continue
 320    continue
c
        call setar1
c
c Contact region of boundary surface
c
        do 410 i=1,nr
        do 410 j=1,nb(i)+1
 410    ncom(i,j)=0
        do 420 ii=1,nr
          do 430 i=1,nb(ii)
c           if(ncom(ii,i).ne.0) goto 430
            n1=ne(nelb(ii,i),neb(ii,i))
            n2=ne(nelb(ii,i),js2(neb(ii,i)))
            do 440 jj=1,nr
              if(jj.eq.ii) goto 440
              do 450 j=1,nb(jj)
                if(ncom(jj,j).ne.0) goto 450
                k1=ne(nelb(jj,j),neb(jj,j))
                k2=ne(nelb(jj,j),js2(neb(jj,j)))
                if(n1.eq.k2 .and. n2.eq.k1) then
                  ncom(ii,i)=jj
                  ncom(jj,j)=ii
                endif
 450          continue
 440        continue
 430      continue
c         ncom(ii,nb(i)+1)=ncom(ii,1)
 420    continue
c
        open(1,file="ncom.res")     !
        do i=1,nb(2)
        write(1,*)"xx(",ne(nelb(2,i),neb(2,i)),")=",
     &                                     xx(ne(nelb(2,i),neb(2,i))),
     &   "yy(",ne(nelb(2,i),neb(2,i)),")=",yy(ne(nelb(2,i),neb(2,i))),
     &        "ncom(2,",i,")=",ncom(2,i)
        enddo
        close(1)
c
c
c Set flag on boundary node
c
        do 460 i=1,nr
        do 460 j=1,nb(i)+1
        ncop(i,j)=0
 460    if(ncom(i,mcy(j-1,nb(i))).ne.0 .or. 
     &                  ncom(i,mcy(j,nb(i))).ne.0) ncop(i,j)=1
c
c Common points of three region
c
        ncpn=0
        do 490 i=1,nr
          if(btest(not(mswi(ndx(i))),0)) goto 490
          do 500 j=1,nb(i)
c       write(*,*) j,ncom(i,mcy(j-1,nb(i))),ncom(i,j)
            if(ncom(i,mcy(j-1,nb(i))) .eq. ncom(i,j)) goto 500
c       write(*,*) "#"
            do 510 k=1,ncpn
 510        if(ne(nelb(i,j),neb(i,j)).eq.kcpn(k)) goto 500
            n=ne(nelb(i,j),neb(i,j))
c       write(*,*) '*',xx(n),yy(n)
            do 520 k=1,ncp
              if(kcp(k).lt.0 .or. icl(k).lt.-900) goto 520
c       write(*,*) k,x0(kcp(k)),y0(kcp(k)),icl(k)
              if(dabs(xx(n)-x0(kcp(k))).lt.1d-12 .and.
     &                  dabs(yy(n)-y0(kcp(k))).lt.1d-12) then
                ncpn=ncpn+1
                kcpn(ncpn)=n
                icpn(ncpn)=i
                jcpn(ncpn)=j
                icln(ncpn)=icl(k)
                stan(ncpn)=staa(k)
                strn(ncpn)=strr(k)
                goto 500
              endif
 520        continue
            id=ndx(i)
            jd=ndx(ncom(i,j))
            jdm=ndx(ncom(i,mcy(j-1,nb(i))))
            im=icoml(id,jd,jdm)
c       write(*,*) '* id,jd,jdm,im,j=',id,jd,jdm,im,j
            if(im.ne.0) then
c       write(*,*) '* ncpn,i,j=',ncpn,i,j
              ncpn=ncpn+1
              kcpn(ncpn)=n
              icpn(ncpn)=i
              jcpn(ncpn)=j
              icln(ncpn)=0
              stan(ncpn)=0
              strn(ncpn)=180
              goto 500
            endif
 500      continue
 490    continue

        do 540 i=1,nr
          do 540 j=1,nb(i)
            if(ncom(i,mcy(j-1,nb(i))) .eq. ncom(i,j)) goto 540
            do 560 k=1,ncpn
 560        if(ne(nelb(i,j),neb(i,j)).eq.kcpn(k)) goto 540
            n=ne(nelb(i,j),neb(i,j))
            ncpn=ncpn+1
            kcpn(ncpn)=n
            icpn(ncpn)=i
            jcpn(ncpn)=j
            icln(ncpn)=-999
            stan(ncpn)=0
            strn(ncpn)=0
 540    continue

        ncp=ncpn
        do 580 i=1,ncp
          icp(i)=icpn(i)
          jcp(i)=jcpn(i)
          kcp(i)=kcpn(i)
          icl(i)=icln(i)
          staa(i)=stan(i)
          strr(i)=strn(i)
 580    continue

        do 930 i=1,ncp
 930    write(*,*) 'i,icp,jcp,icl=',i,icp(i),jcp(i),icl(i)
c       pause
      return
      end
c
        subroutine setar1
        implicit double precision (a-h,o-z)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /regcon/ncon(21),kcon(21,21)
        common /sufcood/xb(21,1001),yb(21,1001)
        nset=0
        do 368 i=1,nr
 368    kset(i)=0
        do 370 i=1,nr
          if(kset(i).ne.0) goto 370
          nset=nset+1
          kset(i)=nset
          do 380 j=1,ncon(i)
          kset(kcon(i,j))=nset
          do 380 k=1,ncon(kcon(i,j))
          kset(kcon(kcon(i,j),k))=nset
          do 380 l=1,ncon(kcon(kcon(i,j),k))
          kset(kcon(kcon(kcon(i,j),k),l))=nset
          do 380 m=1,ncon(kcon(kcon(kcon(i,j),k),l))
          kset(kcon(kcon(kcon(kcon(i,j),k),l),m))=nset
          do 380 n=1,ncon(kcon(kcon(kcon(kcon(i,j),k),l),m))
 380      kset(kcon(kcon(kcon(kcon(kcon(i,j),k),l),m),n))=nset
 370    continue
c
        write(*,*) 'nset=',nset
       write(*,*) 'kset=',(kset(i),i=1,nr)
c
        do 410 i=1,nr
          if(kset(i).eq.0) then
          write(*,*) 'setar1/you need further iteration in DO 370'
          stop
          endif
 410    continue
        return
        end
c
        subroutine calreg(ic,ncal,kcal)
        implicit double precision (a-h,o-z)
        dimension kcal(21)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /regcon/ncon(21),kcon(21,21)
        common /reginf/ msw,mswi(-9:9),etas(-9:9),mfin(-9:9)
c
        ncal=0
        do 320 i=1,nr
320     kcal(i)=0
c
        do 370 i=1,nr
c      write(*,*) btest(mswi(ndx(i)),ic)
        if(btest(mswi(ndx(i)),ic) .and. kcal(i).eq.0) then
        ncal=ncal+1
        kcal(i)=ncal
        do 380 j=1,ncon(i)
        if(btest(mswi(ndx(kcon(i,j))),ic)) then
        kcal(kcon(i,j))=ncal
        do 390 k=1,ncon(kcon(i,j))
        if(btest(mswi(ndx(kcon(kcon(i,j),k))),ic)) then
          kcal(kcon(kcon(i,j),k))=ncal
        do 400 l=1,ncon(kcon(kcon(i,j),k))
        if(btest(mswi(ndx(kcon(kcon(kcon(i,j),k),l))),ic)) then
        kcal(kcon(kcon(kcon(i,j),k),l))=ncal
        do 410 m=1,ncon(kcon(kcon(kcon(i,j),k),l))
       if(btest(mswi(ndx(kcon(kcon(kcon(kcon(i,j),k),l),m)))
     &      ,ic)) then
        kcal(kcon(kcon(kcon(kcon(i,j),k),l),m))=ncal
        do 420 n=1,ncon(kcon(kcon(kcon(kcon(i,j),k),l),m))
        if(btest(mswi(ndx(kcon(kcon(kcon(kcon(kcon(i,j),k),l),m),n)))
     & ,ic)) 
     &  kcal(kcon(kcon(kcon(kcon(kcon(i,j),k),l),m),n))=ncal
420     continue
        endif
410     continue
        endif
400     continue
        endif
390     continue
        endif
380     continue
        endif
370     continue
c
c       write(*,*) 'ic,ncal=',ic,ncal
c       write(*,*) 'kcal=',(kcal(i),i=1,nr)
c
        return
        end

c
!       subroutine setar1sub(irc,nn)
!       common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
!     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
!       common /regcon/ncon(21),kcon(21,21)
!       common /reginf/ msw,mswi(-9:9),mfin(-9:9),etas(9)
!        ni=mswi(ndx(nn))
!       if(irc.eq.0) then
!         kset(nn)=nset
!         if(ni.ge.2) irc=nn
!       else
!         if(ni.lt.2) then
!            kset(nn)=nset
!         else
!           do 100 i=1,ncon(nn)
!             if(kcon(nn,i).eq.irc) then
!                kset(nn)=nset
!               return
!             endif
! 100       continue
!         endif
!       endif
!       return
!       end
c
        subroutine genirx(nelt,irx,idx,ine)
        dimension irx(0:ine),idx(0:ine)
c       integer*2 nec
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /dnec/ nec(8001,3)
c
        do 10 i=1,nelt
 10     irx(i)=0
c
        do 40 i=1,nr
        do 40 j=1,nb(i)
 40     irx(nelb(i,j))=i
c
 80     do 100 n1=1,nelt
          if(irx(n1).eq.0) goto 100
          irc=irx(n1)
          idc=idx(n1)
          do 105 i=1,3
            n2=nec(n1,i)
            if(irx(n2).eq.0 .and. idx(n2).eq.idc) then
              irx(n2)=irc
              do 110 j=1,3
                n3=nec(n2,j)
                if(irx(n3).eq.0 .and. idx(n3).eq.idc) then
                  irx(n3)=irc
                  do 120 k=1,3
                  n4=nec(n3,k)
 120              if(irx(n4).eq.0 .and.idx(n4).eq.idc) irx(n4)=irc
                endif
 110          continue
            endif
 105      continue
 100    continue
c
 210    ifg=0
        do 220 i=1,nelt
          if(irx(i).ne.0) goto 220
          ifg=1
c         write(*,*) i
          do 240 k=1,3
            n=nec(i,k)
            if(irx(n).ne.0 .and. idx(n).eq.idx(i)) then
              irx(i)=irx(n)
              goto 220
            endif
 240      continue
 220    continue
        if(ifg.eq.1) goto 80
c
        return
        end
c
      subroutine clos
      close(1)
      close(2)
      return
      end

        subroutine mesh(nxg,nyg,xg,yg,xgi,ygi,kxs,kys,mg,xx,yy,
     &                  ind,igerr,ixg,iyg)
        implicit double precision (a-h,o-z)
        dimension xx(ind),yy(ind)
        dimension xg(ixg,iyg),yg(ixg,iyg),xgi(ixg),ygi(iyg)
        dimension kxs(21),kys(21)
        integer*2 mg(ixg,iyg)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /sufcood/xb(21,1001),yb(21,1001)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
c
        write(*,*)"in mesh"
        open(1,file="mesh.res")         !
        rewind(1)
        do i=1,nr
        write(1,*)"nb(",i,")=",nb(i)
        do j=1,nb(i)
        write(1,*)"xb(",i,",",j,")=",xb(i,j),"yb(",i,",",j,")=",yb(i,j)
        enddo
        enddo
        close(1)
        if(igerr.eq.0) then
          write(*,*)"igerr=",igerr
          write(*,*)"gmesh started in mesh"        !
          call gmesh(nxg,nyg,xgi,ygi,ixg,iyg)
          write(*,*)"gmesh ended in mesh"          !
          do 100 i=1,nxg
          do 100 j=1,nyg
          xg(i,j)=xgi(i)
 100      yg(i,j)=ygi(j)
        else
          write(*,*)"meshadd started in mesh"     !
          call meshadd(nxg,nyg,xg,yg,ixg,iyg)
          write(*,*)"meshadd ended in mesh"       !
        endif
c
          do 120 i=1,nxg
          do 120 j=1,nyg
 120      mg(i,j)=0
        write(*,*)"It's no prob before fitting start points"   !
c
c       call pltm(nxg,nyg,xg,yg,ixg,iyg)
c
c Fitting start points of broundary with grid
c
        do 220 k=1,nr
          x=xb(k,1)
          y=yb(k,1)
          do 222 i=1,k-1
            if(dabs(x-xb(i,1)).lt.1d-12 .and. 
     &                          dabs(y-yb(i,1)).lt.1d-12) then
              kxs(k)=kxs(i)
              kys(k)=kys(i)
              goto 220
            endif
 222      continue
          smin=9999
          do 230 i=1,nxg
          do 230 j=1,nyg
            s=(x-xg(i,j))**2+(y-yg(i,j))**2
            if(s.lt.smin) then
              smin=s
              imin=i
              jmin=j
            endif
 230      continue
          if(mg(imin,jmin).eq.0) then
            xg(imin,jmin)=x
            yg(imin,jmin)=y
            mg(imin,jmin)=1
            kxs(k)=imin
            kys(k)=jmin
          else
            write(*,*) 'mesh/imin,jmin,mk=',imin,jmin,mg(imin,jmin)
            write(*,*) 'xg,yg=',xg(imin,jmin),yg(imin,jmin)
            stop
          endif
 220    continue
        write(*,*)"It's no prob before fitting grid"
c
c               call pltm(nxg,nyg,xg,yg,ixg,iyg)
c
c Fitting Grid to Characteristic Point
c
                                                !Common points
        write(*,*)"ncp=",ncp    !
        do 500 k=1,ncp
          if(kcp(k).le.0) goto 500
c         x=xb(icp(k),jcp(k))
c         y=yb(icp(k),jcp(k))
          x=xx(kcp(k))
          y=yy(kcp(k))
c       write(*,*) k,x,y
          smin=9999
          do 510 i=1,nxg
          do 510 j=1,nyg
            s=(x-xg(i,j))**2+(y-yg(i,j))**2
            if(s.lt.smin) then
              smin=s
              imin=i
              jmin=j
            endif
 510      continue
c       write(*,*) k,imin,jmin
          if(mg(imin,jmin).eq.0) then
            xg(imin,jmin)=x
            yg(imin,jmin)=y
            mg(imin,jmin)=1
c           write(*,*) 'k,imin,jmin=',k,imin,jmin
          endif
c           write(*,*) '!k,imin,jmin=',k,imin,jmin
 500    continue
c
c       call pltm(nxg,nyg,xg,yg,ixg,iyg)
      return
      end
c
        subroutine insertgi(nzg,zg,z,k,jg,izg,ig)
        implicit double precision (a-h,o-z)
        dimension zg(izg)
        integer*2 jg(2,ig)
c
        if(z.le.zg(1)) then
          dz=(zg(2)-zg(1))*0.1
          if(zg(1)-z.lt.1d-6) then
            zg(1)=z
            jg(k,1)=1
          elseif(zg(1)-z.lt.dz .and. jg(k,1).eq.0) then
            zg(1)=z
            jg(k,1)=1
          else
c           z1=zg(1)
c           z2=zg(2)
            do 10 j=nzg,1,-1
            jg(k,j+1)=jg(k,j)
 10         zg(j+1)=zg(j)
            zg(1)=z
            jg(k,1)=1
            nzg=nzg+1
c           write(*,*) 'insertgi/(1)k,z=',k,z,z1,z2
          endif
          return
        endif
c
        do 100 i=2,nzg
 100    if(zg(i).ge.z) goto 120
c
        z1=zg(i)
        zg(i)=z
        jg(k,i)=1
        nzg=nzg+1
c       write(*,*) 'insertgi/(2)k,z=',k,z,z1
        return
c
 120    dz=zg(i)-zg(i-1)
c
        if(zg(i)-z.lt.dz*0.1) then
           if(jg(k,i).eq.0) zg(i)=z
          jg(k,i)=1
        elseif(z-zg(i-1).lt.dz*0.1) then
          if(jg(k,i-1).eq.0) zg(i-1)=z
          jg(k,i-1)=1
        elseif(zg(i)-z.le.dz*0.5 .and. jg(k,i).eq.0) then
c         zgold=zg(i)
          zg(i)=z
          jg(k,i)=1
c         write(*,*) 'move/(1)k,z,i=',k,z,zgold
        elseif(z-zg(i-1).lt.dz*0.5 .and. jg(k,i-1).eq.0) then
c         zgold=zg(i-1)
          zg(i-1)=z
          jg(k,i-1)=1
c         write(*,*) 'move/(2)k,z,i=',k,z,zgold
        else
c         z1=zg(i-1)
c         z2=zg(i)
          do 140 j=nzg,i,-1
          jg(k,j+1)=jg(k,j)
 140      zg(j+1)=zg(j)
          zg(i)=z
           jg(k,i)=1
          nzg=nzg+1
c         write(*,*) 'insertgi/(3)k,z=',k,z,z1,z2
        endif
        return
        end
c
        subroutine checkelemnt(nel,ne,xx,yy,idx,ar0,hmax,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),idx(ine),ar0(ine)
        hmax=-1e10
        do 100 i=1,nel
c         if(idx(i).lt.0) goto 100
          ar=(xx(ne(i,2))-xx(ne(i,1)))*(yy(ne(i,3))-yy(ne(i,1)))
     &       -(xx(ne(i,3))-xx(ne(i,1)))*(yy(ne(i,2))-yy(ne(i,1)))
          if(ar.le.0d0) then
            hmax=9999
            return
          else
            hmax=dmax1(ar/ar0(i),ar0(i)/ar,hmax)
          endif
 100    continue
        return
        end
c
        subroutine grdfcn(nel,np,ne,xx,yy,idx,irx,jnw,
     &          neln,npn,nen,xn,yn,idn,irn,
     &          uu,pp,vv,uun,vvn,ppn,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine),irx(0:ine)
        dimension nen(ine,3),xn(ind),yn(ind),idn(0:ine),irn(0:ine)
        dimension jnw(ind)
        dimension uu(ind),vv(ind),pp(ind)
        dimension uun(ind),vvn(ind),ppn(ind)
        dimension al(3),al0(3),jp1(3),jp2(3)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /ard/ arr(8001)
        data jp1/2,3,1/,jp2/3,1,2/

c        write(*,*) 'grdfcn'
        do 20 i=1,nel
 20     arr(i)=(xx(ne(i,2))-xx(ne(i,1)))*(yy(ne(i,3))-yy(ne(i,1)))
     &      -(xx(ne(i,3))-xx(ne(i,1)))*(yy(ne(i,2))-yy(ne(i,1)))
        do 40 i=1,neln
        do 40 j=1,3
        uun(nen(i,j))=0d0
        vvn(nen(i,j))=0d0
        ppn(nen(i,j))=0d0
 40     jnw(nen(i,j))=0
c
        do 100 m=1,neln         !neln(nset,msw-1)
          idc=idn(m)
          if(idc.lt.0) goto 100
          do 120 k=1,3
            n=nen(m,k)
            if(jnw(n).eq.1) goto 120
            x=xn(n)
            y=yn(n)
            jnw(n)=1
            almin=1d10
            do 140 i=1,nel
             if(idx(i).ne.idc) goto 140
             alm1=1000
             alm2=-1000
             do 150 j=1,3
               al(j)=((xx(ne(i,jp1(j)))-x)*(yy(ne(i,jp2(j)))-y)
     &             -(xx(ne(i,jp2(j)))-x)*(yy(ne(i,jp1(j)))-y))/arr(i)
               if(al(j).lt.-5d0 .or. al(j).gt.5d0) goto 140
               alm1=dmin1(alm1,al(j))
               alm2=dmax1(alm2,al(j))
 150         continue
             if(alm1.ge.0 .and. alm1.le.1d0 .and. 
     &                  alm2.ge.0d0 .and. alm2.le.1d0) goto 200
             alm=dmax1(-alm1,alm2-1d0)
             if(alm.lt.almin) then
               almin=alm
               imin=i
               do 145 j=1,3
 145           al0(j)=al(j)
             endif
 140      continue
c
          i=imin
          do 155 j=1,3
 155      al(j)=al0(j)
c               write(*,*) 'm,imin,almin,idc=',m,i,sngl(almin),idc
 200      do 210 j=1,3
          uun(n)=uun(n)+al(j)*uu(ne(i,j))
          vvn(n)=vvn(n)+al(j)*vv(ne(i,j))
 210      ppn(n)=ppn(n)+al(j)*pp(ne(i,j))
c210      ttn(n)=ttn(n)+al(j)*tt(ne(i,j))
 120    continue
 100    continue
c
        nel=neln
        np=npn
        do 820 j=1,nel
        idx(j)=idn(j)
        irx(j)=irn(j)
        do 820 k=1,3
 820    ne(j,k)=nen(j,k)
        do 830 j=1,np
        uu(j)=uun(j)
        vv(j)=vvn(j)
        pp(j)=ppn(j)
        xx(j)=xn(j)
 830    yy(j)=yn(j)
 800    continue
      return
      end
c
c    added subroutines more
c
        subroutine sufelew
        implicit double precision (a-h,o-z)
          dimension j12(3)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
          common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                   ncomw(21,1001),ncopw(21,1001)
          data j12/4,5,6/
c
          do lr=1,nr
            nbw(lr)=0
            do k=1,nb(lr)
              nbw(lr)=nbw(lr)+1
              nelbw(lr,nbw(lr))=nelb(lr,k)
              nebw(lr,nbw(lr))=neb(lr,k)
              ncomw(lr,nbw(lr))=ncom(lr,k)
              ncopw(lr,nbw(lr))=ncop(lr,k)
              nbw(lr)=nbw(lr)+1
              nelbw(lr,nbw(lr))=nelb(lr,k)
              nebw(lr,nbw(lr))=j12(neb(lr,k))
              ncomw(lr,nbw(lr))=ncom(lr,k)
              ncopw(lr,nbw(lr))=ncop(lr,k)
            enddo
c
            nelbw(lr,-1)=nelbw(lr,nbw(lr)-1)
            nelbw(lr,0)=nelbw(lr,nbw(lr))
            nelbw(lr,nbw(lr)+1)=nelbw(lr,1)
            nelbw(lr,nbw(lr)+2)=nelbw(lr,2)
          nebw(lr,-1)=nebw(lr,nbw(lr)-1)
            nebw(lr,0)=nebw(lr,nbw(lr))
            nebw(lr,nbw(lr)+1)=nebw(lr,1)
            nebw(lr,nbw(lr)+2)=nebw(lr,2)         
c
          ncomw(lr,nbw(lr)+1)=ncomw(lr,1)
            ncomw(lr,nbw(lr)+2)=ncomw(lr,2)
            ncopw(lr,nbw(lr)+1)=ncopw(lr,1)
            ncopw(lr,nbw(lr)+2)=ncopw(lr,2)
          enddo
c	    
	    open(1,file="ncomw.res")
	    do lr=1,nr
	    write(1,*)"lr,nbw(lr)",lr,nbw(lr)
	    do i=1,nbw(lr)
	    write(1,*)"ncomw(",lr,",",i,")",ncomw(lr,i)
	    enddo
	    enddo
	    close(1)
c
          return
          end                      
c
        subroutine adjust2(nxg,nyg,xg,yg,mk,mg,ncop,nbp,xb,yb,ii0,jj0
     &                  ,idc,ixg,iyg,igerr)
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension xb(1001),yb(1001)
        dimension ncop(0:1001)   !使わない
        dimension ss(4),ik(4),jk(4)
        integer*2 mk(ixg,iyg),mg(ixg,iyg)
        integer*2 is0(8),is1(8),js0(8),js1(8),it(8),jt(8),kt(8)
        logical ifg,icross
        common /addmesh/ i0,j0,i1,j1
        data is0/1,1,0,-1,-1, -1, 0, 1/, is1/1,0,-1,-1,-1,0, 1, 1/
        data js0/0,1,1, 1, 0,-1,-1, -1/ ,js1/1,1, 1, 0,-1,-1,-1,0/
        data it/0,0,-1,-1,-1,-1,0,0/,jt/0,0,0,0,-1,-1,-1,-1/
        data kt/1,1,2,2,1,1,2,2/
        data ik/0,1,1,0/,jk/0,0,1,1/
c-------------
        write(*,*)"in adjust ii0=",ii0,"jj0=",jj0     !
c       write(*,*) 'adjust'
c       write(*,*) nxg,nyg,nbp,ii0,jj0,ixg,iyg
c       stop
c
        ib0=1
        ii=ii0
        jj=jj0
        mg(ii,jj)=1    !sub meshで節点座標が代入されている
        write(*,*)"lq,x1,y1",1,xb(1),yb(1)
        do 210 lp=1,10000
          do 230 lq=ib0,nbp-1
          do 230 k=1,8
            i0=ii+is0(k)
            j0=jj+js0(k)
            i1=ii+is1(k)
            j1=jj+js1(k)
            if(i0.lt.1 .or. j0.lt.1) goto 230
            if(i1.lt.1 .or. j1.lt.1) goto 230
            if(i0.gt.nxg .or. j0.gt.nyg) goto 230
            if(i1.gt.nxg .or. j1.gt.nyg) goto 230
            xg0=xg(i0,j0)
            yg0=yg(i0,j0)
            xg1=xg(i1,j1)
            yg1=yg(i1,j1)
            ifg=icross(xb,yb,xg0,yg0,xg1,yg1,lq)
            if(ifg)goto 240
 230      continue
          write(*,*)"ib0=",ib0,"npb-1=",nbp-1
            stop
c
 240      continue
          x1=xb(lq+1)     !考慮する境界上の節点
          y1=yb(lq+1)
          write(*,*)sngl(xg0),sngl(yg0),sngl(xg1),sngl(yg1)
          write(*,*)"lq,x0,y0",lq,xb(lq),yb(lq)
          write(*,*)"lq+1,x1,y1",lq+1,x1,y1
c          pause 
c
          ep=1d-12
            if(dabs(x1-xg(ii0,jj0)).lt.ep .and.    !終了条件
     &         dabs(y1-yg(ii0,jj0)).lt.ep)then     
            if(mk(ii+it(k),jj+jt(k)).eq.0)mk(ii+it(k),jj+jt(k))=kt(k)
            goto 4000   !終了
            endif
c 
            do 500 i=-100,100      !節点の存在するグリッドを検索
            do 600 j=-100,100
              if(ii+i.lt.1 .or. ii+i.gt.nxg)goto 500
              if(jj+j.lt.1. or. jj+j.gt.nyg)goto 600
              if(x1.ge.xg(ii+i,jj+j) .and. x1.le.xg(ii+i+1,jj+j+1))then
              if(y1.ge.yg(ii+i,jj+j) .and. 
     &           y1.le.yg(ii+i+1,jj+j+1))goto 2000
              endif
600       continue
500       continue
          write(*,*)"didn't find"
          stop
c       
2000      continue
          iw=ii
          jw=jj
          ii=ii+i     !仮
          jj=jj+j
          do i=1,4                        !各頂点から節点までの距離をss(4)に格納
          ss(i)=dsqrt((x1-xg(ii+ik(i),jj+jk(i)))**2
     &               +(y1-yg(ii+ik(i),jj+jk(i)))**2)   
          enddo
          do i=1,4                         !節点が既に使われていないかチェック
          kk=nsmall(i,ss,4)                 !i番目に小さな値が格納された配列の番号を代入
          if(mg(ii+ik(kk),jj+jk(kk)).ne.1)goto 3000
          enddo
          write(*,*)"all grid points are already considered"
          stop
c
3000      continue
          write(*,*)"xg,yg",xg(ii+ik(kk),jj+jk(kk))
     &                     ,yg(ii+ik(kk),jj+jk(kk)),"修正前"
          mg(ii+ik(kk),jj+jk(kk))=1          !最も節点に近い頂点(格子点)に節点の座標を代入
          xg(ii+ik(kk),jj+jk(kk))=x1
          yg(ii+ik(kk),jj+jk(kk))=y1
          write(*,*)"xg,yg",xg(ii+ik(kk),jj+jk(kk))
     &                     ,yg(ii+ik(kk),jj+jk(kk)),"修正後"
          call chek(ii,jj,xg,yg,kk,ixg,iyg)   !内角が180°以上の頂点がないかチェック．あったら修正
          if(mk(iw+it(k),jw+jt(k)).eq.0)
     &         mk(iw+it(k),jw+jt(k))=kt(k)
          ii=ii+ik(kk)
          jj=jj+jk(kk)
          ib0=lq+1
c          call pltg(npb,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
210       continue
c
4000    continue
c       call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
c       pause
        return
        end
c
        subroutine chek(ii,jj,xg,yg,kk,ixg,iyg)
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension ik(4),jk(4),ij(4,4)
        data ik/0,1,1,0/,jk/0,0,1,1/
c
c        write(*,*)"in check"
        do i=1,4
          k=0
          do j=1,4
          if(i.eq.2.and.j.eq.4)k=4
          if(i.eq.3.and.j.eq.3)k=4
          if(i.eq.4.and.j.eq.2)k=4
        ij(i,j)=i+j-1-k
          enddo
          enddo 
        x1=xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))   !グリッドの4つの頂点
          x2=xg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))
          x3=xg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))
          x4=xg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))
          y1=yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))
          y2=yg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))
          y3=yg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))
          y4=yg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))
c       対角線の傾きと切片
        if(dabs(x2-x4).lt.1d-12)then  !al=∞,am≠∞
c          write(*,*)"al=∞,am≠∞"
1100      am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=x2
          if(xx.gt.dmin1(x1,x3) .and. xx.lt.dmax1(x1,x3))return
          write(*,*)"x1,y1修正 al=∞,am≠∞"
          x1=x1-(x3-x1)/2d0    !x1,y1を外側にずらす
          if(dabs(y1-y4).gt.dabs(y1-y2))y1=y1+(y4-y1)/2d0
          if(dabs(y1-y4).lt.dabs(y1-y2))y1=y1+(y2-y1)/2d0
          xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
          yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=y1
          write(*,*)"al=∞,am≠∞"
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1100
        endif
c
        if(dabs(x1-x3).lt.1d-12)then    !al≠∞,am=∞
c          write(*,*)"al≠∞,am=∞"
          al=(y2-y4)/(x2-x4)
          bl=(x2*y4-x4*y2)/(x2-x4)
          xx=x1
          yy=al*x1+bl
1200      if(yy.gt.dmin1(y1,y3) .and. yy.lt.dmax1(y1,y3))return
          write(*,*)"x1,y1修正 al≠∞,am=∞"
          y1=y1-(y3-y1)/2d0           !y1を外側にずらす.x1=const
          xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
          yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=y1
          write(*,*)"al≠∞,am=∞"
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1200
        endif
c
c          write(*,*)"al≠∞,am≠∞"
          al=(y2-y4)/(x2-x4)     !al≠∞,am≠∞
          bl=(x2*y4-x4*y2)/(x2-x4)
1000      am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=-(bl-bm)/(al-am)    !交点
          yy=am*xx+bm
          if(xx.gt.dmin1(x1,x3) .and. xx.lt.dmax1(x1,x3)  .and.  !交点が存在するか
     &   xx.gt.dmin1(x2,x4) .and. xx.lt.dmax1(x2,x4))return
          write(*,*)"x1,y1修正 al≠∞,am≠∞"
          if(dabs(y1-y4).ge.dabs(y1-y2))x1=x1-(x3-x1)/2d0    !交点が存在しなかった場合，頂点を外側にずらす
          if(dabs(y1-y4).lt.dabs(y1-y2))y1=y1-(y3-y1)/2d0
          xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
          yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=y1
          write(*,*)"al≠∞,am≠∞"
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1000
c
          end
c
           function icross(xb,yb,xg0,yg0,xg1,yg1,lq)
          implicit double precision(a-h,o-z)
          dimension xb(1001),yb(1001)
          logical icross
c
          icross=.FALSE.
          x0=xb(lq)
          y0=yb(lq)
           x1=xb(lq+1)
          y1=yb(lq+1)
c         2つの線分(対角線)の傾きal,amと切片bl,bm
          if(dabs(xg1-xg0).lt.1d-12 .and. dabs(x1-x0).lt.1d-12)then   !al=∞,am=∞
            write(*,*)"al=∞,am=∞ in icross"
            if(dabs(xg1-x1).lt.1d-12)then
		    if(dmax1(y0,y1).gt.dmin1(yg0,yg1) .or. 
     &           dmin1(y0,y1).lt.dmax1(yg0,yg1))icross=.TRUE.
	      endif
            goto 1000
          endif
c
          if(dabs(xg1-xg0).lt.1d-12)then      !al=∞,am≠∞
            write(*,*)"al=∞,am≠∞ in icross"
            am=(y1-y0)/(x1-x0)
            bm=(x1*y0-x0*y1)/(x1-x0)
            xx=xg1
              yy=am*xg1+bm
            if(yy.ge.dmin1(yg0,yg1) .and. yy.le.dmax1(yg0,yg1) .and.
     &       xx.ge.dmin1(x0,x1) .and. xx.le.dmax1(x0,x1))
     &                                               icross=.TRUE.
            goto 1000
           endif
c
          if(dabs(x1-x0).lt.1d-12)then        !al≠∞,am=∞    
             write(*,*)"al≠∞,am=∞ in icross"
            al=(yg1-yg0)/(xg1-xg0)
            bl=(xg1*yg0-xg0*yg1)/(xg1-xg0)
             xx=x1
            yy=al*x1+bl
            if(yy.ge.dmin1(y0,y1) .and. yy.le.dmax1(y0,y1) .and.
     &     xx.ge.dmin1(xg0,xg1) .and. xx.le.dmax1(xg0,xg1))
     &                                            icross=.TRUE.
            goto 1000
          endif
c		 
          write(*,*)"al≠∞,am≠∞ in icross"
          al=(yg1-yg0)/(xg1-xg0)              !al≠∞,am≠∞
          bl=(xg1*yg0-xg0*yg1)/(xg1-xg0)
          am=(y1-y0)/(x1-x0)
          bm=(x1*y0-x0*y1)/(x1-x0)
          if(dabs(al-am).lt.1d-12)then  !al≠∞,am≠∞ & al=am    傾きが等しい
            if(dabs(bl-bm).lt.1d-12)icross=.TRUE.
            goto 1000
          endif
          xx=-(bl-bm)/(al-am)        !交点のx座標
          if(xx.ge.dmin1(x0,x1) .and. xx.le.dmax1(x0,x1) .and.   !交点が存在するか
     &     xx.ge.dmin1(xg0,xg1) .and. xx.le.dmax1(xg0,xg1))
     &                                            icross=.TRUE.
c
1000    continue  !終了
        end
c         
        function nsmall(k,ss,n)
          implicit double precision(a-h,o-z)
          dimension ss(n),num(n)
c
        do 1000 ii=1,n     !ii番目に小さい値が格納されている配列の番号を求める
          smin=9999
            do 2000 i=1,n
                do j=1,ii-1     !既に格納された値は考慮しない
              if(i.eq.num(j))goto 2000
            enddo
            if(ss(i).lt.smin)then
            smin=ss(i)
            num(ii)=i
            endif
2000        continue
1000      continue
        nsmall=num(k)
        end
c
        subroutine adjust3(nxg,nyg,xg,yg,mk,mg,ncop,nbp,xb,yb,ii0,jj0
     &                  ,idc,ixg,iyg,igerr)
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension xb(1001),yb(1001)
        dimension ncop(0:1001)   !使わない
        dimension ss(4),is(4),js(4),ik(4),jk(4)
	  dimension is2(8),js2(8)
c	  dimension k2(8),k3(8)    !sub chekに使う
        integer*2 mk(ixg,iyg),mg(ixg,iyg)
        integer*2 is0(8),is1(8),js0(8),js1(8),it(8),jt(8),kt(8)
        logical ifg
        common /addmesh/ i0,j0,i1,j1
        data is0/1,1,0,-1,-1, -1, 0, 1/, is1/1,0,-1,-1,-1,0, 1, 1/
        data js0/0,1,1, 1, 0,-1,-1, -1/ ,js1/1,1, 1, 0,-1,-1,-1,0/
        data it/0,0,-1,-1,-1,-1,0,0/,jt/0,0,0,0,-1,-1,-1,-1/
        data kt/1,1,2,2,1,1,2,2/
        data ik/0,1,1,0/,jk/0,0,1,1/
	  data is2/1,0,0,-1,-1,0,0,1/,js2/0,1,1,0,0,-1,-1,0/
c	  data k2/2,3,3,4,4,1,1,2/,k3/3,4,4,1,1,2,2,3/
c-------------
        write(*,*)"in adjust ii0=",ii0,"jj0=",jj0     !
c       write(*,*) 'adjust'
c       write(*,*) nxg,nyg,nbp,ii0,jj0,ixg,iyg
c       stop
c
        ib0=1
        ii=ii0
        jj=jj0
        mg(ii,jj)=1    !sub meshで節点座標が代入されている
        write(*,*)"lq,x1,y1",1,xb(1),yb(1)
        do 210 lp=1,10000
          do 235 lq=ib0,nbp-1
          do 230 k=1,8
            i0=ii+is0(k)
            j0=jj+js0(k)
            i1=ii+is1(k)
            j1=jj+js1(k)
            if(i0.lt.1 .or. j0.lt.1) goto 230
            if(i1.lt.1 .or. j1.lt.1) goto 230
            if(i0.gt.nxg .or. j0.gt.nyg) goto 230
            if(i1.gt.nxg .or. j1.gt.nyg) goto 230
            xg0=xg(i0,j0)
            yg0=yg(i0,j0)
            xg1=xg(i1,j1)
            yg1=yg(i1,j1)
            call igrid(xb,yb,xg0,yg0,xg1,yg1,lq,xc,yc,ifg)
            if(ifg)goto 240
 230      continue
          if(dabs(xb(lq+1)-xg(ii0,jj0)).lt.1d-12 .and.
     &       dabs(yb(lq+1)-yg(ii0,jj0)).lt.1d-12)then
	      write(*,*)"error in deciding ending point"
	      stop
	    endif
 235      continue
          write(*,*)"ib0=",ib0,"npb-1=",nbp-1
            stop
c
 240      continue
          x1=xb(lq+1)     !考慮する境界上の節点
          y1=yb(lq+1)
          write(*,*)sngl(xg0),sngl(yg0),sngl(xg1),sngl(yg1)
          write(*,*)"lq,x0,y0",lq,sngl(xb(lq)),sngl(yb(lq))
          write(*,*)"lq+1,x1,y1",lq+1,x1,y1
	    write(*,*)"xc,yc",xc,yc
c       節点の存在するグリッドを検索
        do 500 j=100,1,-1
        do 600 i=100,-100,-1
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 600
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 500
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &       (x1.ge.xg(ii+i,jj+j+1) .and. 
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &	   y1.le.yg(ii+i+1,jj+j+1)))goto 2000
          endif
600     continue
500     continue
        j=0
        do 550 i=100,0,-1
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 550
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 550
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &       (x1.ge.xg(ii+i,jj+j+1) .and. 
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &       y1.le.yg(ii+i+1,jj+j+1)))goto 2000
          endif
550     continue
        do 650 i=-100,-1
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 650
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 650
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &       (x1.ge.xg(ii+i,jj+j+1) .and. 
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &       y1.le.yg(ii+i+1,jj+j+1)))goto 2000
          endif
650     continue
        j=-1
        do 750 i=-100,-1
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 750
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 750
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &       (x1.ge.xg(ii+i,jj+j+1) .and. 
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &       y1.le.yg(ii+i+1,jj+j+1)))goto 2000
	    endif
750     continue
        do 700 j=-2,-100,-1      !節点の存在するグリッドを検索
        do 800 i=-100,100
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 800
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 700
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &       (x1.ge.xg(ii+i,jj+j+1) .and.
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &       y1.le.yg(ii+i+1,jj+j+1)))goto 2000
	    endif
800     continue
700     continue
        j=-1
        do 850 i=100,0,-1
          if(ii+i.lt.1 .or. ii+i+1.gt.nxg)goto 850
          if(jj+j.lt.1. or. jj+j+1.gt.nyg)goto 850
          if((x1.ge.xg(ii+i,jj+j) .and. 
     &       x1.le.xg(ii+i+1,jj+j)) .or.
     &	   (x1.ge.xg(ii+i,jj+j+1) .and.
     &       x1.le.xg(ii+i+1,jj+j+1)))then
          if((y1.ge.yg(ii+i,jj+j) .and. 
     &       y1.le.yg(ii+i,jj+j+1)) .or.
     &       (y1.ge.yg(ii+i+1,jj+j) .and.
     &       y1.le.yg(ii+i+1,jj+j+1)))goto 2000
	    endif
850     continue
        write(*,*)"didn't find"
        stop
c
2000    continue
        do l=1,4  !節点を含むグリッドがii,jjを含むグリッドの周囲に存在するか調べる
	    if(ii+i+ik(l).eq.i0 .and. jj+j+jk(l).eq.j0)then
	    do m=1,4
	      if(ii+i+ik(m).eq.i1 .and. jj+j+jk(m).eq.j1)goto 2500
	    enddo
	    endif
	  enddo
c
          write(*,*)"節点見つからない"
	    write(*,*)"k,ii,jj",k,ii,jj
	    write(*,*)"ii*,jj*",ii+i,jj+j
          if(k.eq.2 .or. k.eq.3 .or. k.eq.6 .or. k.eq.7)then    !ii,jjを含むグリッド
c
	      if(dabs(xg(i0,j0)-xc).lt.dabs(xg(i1,j1)-xc))then
	        if(mg(i0,j0).ne.1)then
	          write(*,*)"xg,yg",xg(i0,j0),yg(i0,j0),"修正前"
	          mg(i0,j0)=1
	          xg(i0,j0)=xc
	          yg(i0,j0)=yc
	          ki=i0
	          kj=j0
	        elseif(mg(i1,j1).ne.1)then
	          write(*,*)"xg,yg",xg(i1,j1),yg(i1,j1),"修正前"
	          mg(i1,j1)=1
	          xg(i1,j1)=xc
	          yg(i1,j1)=yc
	          ki=i1
	          kj=j1
	        else
	          write(*,*)"error1節点が見つからない"
	          pause
	        endif
	      else
	        if(mg(i1,j1).ne.1)then
	          write(*,*)"xg,yg",xg(i1,j1),yg(i1,j1),"修正前"
	          mg(i1,j1)=1
	          xg(i1,j1)=xc
	          yg(i1,j1)=yc
	          ki=i1
	          kj=j1
	        elseif(mg(i0,j0).ne.1)then
	          write(*,*)"xg,yg",xg(i0,j0),yg(i0,j0),"修正前"
	          mg(i0,j0)=1
	          xg(i0,j0)=xc
	          yg(i0,j0)=yc
	          ki=i0
	          kj=j0
	        else
	          write(*,*)"error2節点が見つからない"
	          pause
	        endif
	      endif
c
	    else
	      if(dabs(yg(i0,j0)-yc).lt.dabs(yg(i1,j1)-yc))then
	        if(mg(i0,j0).ne.1)then
	          write(*,*)"xg,yg",xg(i0,j0),yg(i0,j0),"修正前"
	          mg(i0,j0)=1
	          xg(i0,j0)=xc
	          yg(i0,j0)=yc
	          ki=i0
	          kj=j0
	        elseif(mg(i1,j1).ne.1)then
	          write(*,*)"xg,yg",xg(i1,j1),yg(i1,j1),"修正前"	          
	          mg(i1,j1)=1
	          xg(i1,j1)=xc
	          yg(i1,j1)=yc
	          ki=i1
	          kj=j1
	        else
	          write(*,*)"error3節点が見つからない"
	          pause
	        endif
	      else
	        if(mg(i1,j1).ne.1)then
	           write(*,*)"xg,yg",xg(i1,j1),yg(i1,j1),"修正前"
	           mg(i1,j1)=1
	           xg(i1,j1)=xc
	           yg(i1,j1)=yc
	           ki=i1
	           kj=j1
	         elseif(mg(i0,j0).ne.1)then
	           write(*,*)"xg,yg",xg(i0,j0),yg(i0,j0),"修正前"	          
	           mg(i0,j0)=1
	           xg(i0,j0)=xc
	           yg(i0,j0)=yc
	           ki=i0
			   kj=j0
			 else
			   write(*,*)"error4節点が見つからない"
	           pause
			 endif
	      endif
	    endif
	    mk(ii+it(k),jj+jt(k))=kt(k)
c
        write(*,*)"xg,yg",xg(ki,kj),yg(ki,kj),"修正後"
	  do 5000 i=1,4    !位置をずらした節点を含むグリッドをチェック
c	  do 5000 j=1,4
	  call chek3(ki-ik(i),kj-jk(i),xg,yg,1,ixg,iyg)
5000    continue
	  ii=ki
	  jj=kj
	  xb(lq)=xg(ii,jj)
	  yb(lq)=yg(ii,jj)
	  ib0=lq
c          call pltg(npb,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
	  pause
	  goto 210        !xb(lq),yb(lq)の位置を変えてやり直す
c
2500      continue        !ii,jjを含むグリッドの周囲のグリッドに節点が存在する
          write(*,*)"節点見つかる"
	    write(*,*)"k,ii,jj",k,ii,jj
	    write(*,*)"ii*,jj*",ii+i,jj+j
          iw=ii
          jw=jj
	    ii=ii+i
	    jj=jj+j
	    is(1)=i0
	    js(1)=j0
	    is(2)=i0+is2(k)
	    js(2)=j0+js2(k)
	    is(3)=i1+is2(k)
	    js(3)=j1+js2(k)
	    is(4)=i1
	    js(4)=j1
c
          if(dabs(x1-xg(ii0,jj0)).lt.1d-12 .and.    !終了条件
     &       dabs(y1-yg(ii0,jj0)).lt.1d-12)then
	      do i=1,4
	        if(is(i).eq.ii0 .and. js(i).eq.jj0)then
	        kk=i
	        goto 4500
	        endif
	      enddo
            write(*,*)"error in end"
	      stop
c
4500        continue
	      write(*,*)"終了条件 kk,xc,yc",kk,xc,yc
          write(*,*)"xg,yg",xg(is(kk),js(kk)),yg(is(kk),js(kk)),"修正前"
	      if(kk.eq.2)then
	        mg(is(1),js(1))=1
	        xg(is(1),js(1))=xc
	        yg(is(1),js(1))=yc
	        ki=is(1)    !ずらした節点
	        kj=js(1)
	      elseif(kk.eq.3)then
	        mg(is(4),js(4))=1
	        xg(is(4),js(4))=xc
	        yg(is(4),js(4))=yc
	        ki=is(4)     !ずらした節点
	        kj=js(4)
	      endif
          write(*,*)"xg,yg",xg(is(kk),js(kk)),yg(is(kk),js(kk)),"修正後"
c            call chek2(is,js,xg,yg,kk,ixg,iyg)    !内角が180°以上の頂点がないかチェック．あったら修正
c	      if(kk.eq.2)call chek(iw+it(k),jw+jt(k),xg,yg,k2(k),ixg,iyg)
c	      if(kk.eq.3)call chek(iw+it(k),jw+jt(k),xg,yg,k3(k),ixg,iyg)
            if(kk.eq.2 .or. kk.eq.3)then
              do 5500 i=1,4    !位置をずらした節点を含むグリッドをチェック
c	        do 5500 j=1,4
	        call chek3(ki-ik(i),kj-jk(i),xg,yg,1,ixg,iyg)
5500          continue
	      endif
            if(mk(iw+it(k),jw+jt(k)).eq.0)mk(iw+it(k),jw+jt(k))=kt(k)
            goto 4000   !終了
          endif
c
	    do i=1,4
          ss(i)=dsqrt((xg(is(i),js(i))-x1)**2+(yg(is(i),js(i))-y1)**2)
	    enddo
	    if((dabs(xg(is(1),js(1))-x1).lt.1d-12 .and.    !格子線上（垂直または水平と仮定）に節点があるとき
     &        dabs(xg(is(1),js(1))-xg(is(2),js(2))).lt.1d-12) .or.
     &       (dabs(yg(is(1),js(1))-y1).lt.1d-12 .and. 
     &        dabs(yg(is(1),js(1))-yg(is(2),js(2))).lt.1d-12))then
	      ss(3)=ss(3)+1d2
	      ss(4)=ss(4)+1d2
	    elseif((dabs(xg(is(2),js(2))-x1).lt.1d-12 .and. 
     &            dabs(xg(is(2),js(2))-xg(is(3),js(3))).lt.1d-12) .or.
     &           (dabs(yg(is(2),js(2))-y1).lt.1d-12 .and. 
     &            dabs(yg(is(2),js(2))-yg(is(3),js(3))).lt.1d-12))then
	      ss(4)=ss(4)+1d2
	      ss(1)=ss(1)+1d2
	    elseif((dabs(xg(is(3),js(3))-x1).lt.1d-12 .and. 
     &            dabs(xg(is(3),js(3))-xg(is(4),js(4))).lt.1d-12) .or.
     &           (dabs(yg(is(3),js(3))-y1).lt.1d-12 .and. 
     &            dabs(yg(is(3),js(3))-yg(is(4),js(4))).lt.1d-12))then
	      ss(1)=ss(1)+1d2
	      ss(2)=ss(2)+1d2
	    elseif((dabs(xg(is(4),js(4))-x1).lt.1d-12 .and. 
     &            dabs(xg(is(4),js(4))-xg(is(1),js(1))).lt.1d-12) .or.
     &           (dabs(yg(is(4),js(4))-y1).lt.1d-12 .and. 
     &            dabs(yg(is(4),js(4))-yg(is(1),js(1))).lt.1d-12))then
	      ss(2)=ss(2)+1d2
	      ss(3)=ss(3)+1d2
          endif
c
	    do i=1,4
	      kk=nsmall(i,ss,4)
	      if(mg(is(kk),js(kk)).ne.1)goto 3000
	    enddo
          write(*,*)"all grid points are already considered"
          stop
c
3000      continue
          write(*,*)sngl(xg(is(1),js(1))),sngl(yg(is(1),js(1)))
          write(*,*)sngl(xg(is(2),js(2))),sngl(yg(is(2),js(2)))
          write(*,*)sngl(xg(is(3),js(3))),sngl(yg(is(3),js(3)))
          write(*,*)sngl(xg(is(4),js(4))),sngl(yg(is(4),js(4)))
          write(*,*)"kk,xc,yc",kk,xc,yc
          write(*,*)"xg,yg",xg(is(kk),js(kk)),yg(is(kk),js(kk)),"修正前"
c
	    if(kk.eq.1 .or. kk.eq.4)then
            mg(is(kk),js(kk))=1      !最も節点に近い頂点(格子点)に節点の座標を代入
            xg(is(kk),js(kk))=x1
            yg(is(kk),js(kk))=y1
	      ki=is(kk)
	      kj=js(kk)
c
	    elseif(kk.eq.2)then
	      if(mg(is(1),js(1)).ne.1)then
              mg(is(kk),js(kk))=1
              xg(is(kk),js(kk))=x1
              yg(is(kk),js(kk))=y1
	        ki=is(kk)
	        kj=js(kk)
	        mg(is(1),js(1))=1
	        xg(is(1),js(1))=xc
	        yg(is(1),js(1))=yc
	        ki2=is(1)
	        kj2=js(1)
	      elseif(mg(is(4),js(4)).ne.1)then
              mg(is(3),js(3))=1
              xg(is(3),js(3))=x1
              yg(is(3),js(3))=y1
              ki=is(3)
	        kj=js(3)
	        mg(is(4),js(4))=1
	        xg(is(4),js(4))=xc
	        yg(is(4),js(4))=yc
	        ki2=is(4)
	        kj2=js(4)
           else
	        write(*,*)"error20"
	        stop
            endif
c
	    elseif(kk.eq.3)then
	      if(mg(is(4),js(4)).ne.1)then
              mg(is(kk),js(kk))=1
              xg(is(kk),js(kk))=x1
              yg(is(kk),js(kk))=y1
	        ki=is(kk)
	        kj=js(kk)
	        mg(is(4),js(4))=1
	        xg(is(4),js(4))=xc
	        yg(is(4),js(4))=yc
	        ki2=is(4)
	        kj2=js(4)
	      elseif(mg(is(1),js(1)).ne.1)then
              mg(is(2),js(2))=1
              xg(is(2),js(2))=x1
              yg(is(2),js(2))=y1
	        ki=is(2)
	        kj=js(2)
	        mg(is(1),js(1))=1
	        xg(is(1),js(1))=xc
	        yg(is(1),js(1))=yc
	        ki2=is(1)
	        kj2=js(1)
	      else
	        write(*,*)"error10"
	        stop
	      endif
	    endif
          write(*,*)"xg,yg",xg(is(kk),js(kk)),yg(is(kk),js(kk)),"修正後"
	    pause
c          call chek2(is,js,xg,yg,kk,ixg,iyg)    !内角が180°以上の頂点がないかチェック．あったら修正
c          if(kk.eq.2)call chek(iw+it(k),jw+jt(k),xg,yg,k2(k),ixg,iyg)
c          if(kk.eq.3)call chek(iw+it(k),jw+jt(k),xg,yg,k3(k),ixg,iyg)
          do 6000 i=1,4   !位置をずらした節点を含むグリッドをチェック
c	    do 6000 j=1,4
	    call chek3(ki-ik(i),kj-jk(i),xg,yg,1,ixg,iyg)
6000      continue
	    if(kk.eq.2 .or. kk.eq.3)then
	      do 6500 i=1,4
c	      do 6500 j=1,4
	      call chek3(ki2-ik(i),kj2-jk(i),xg,yg,1,ixg,iyg)
6500        continue
	    endif
          if(mk(iw+it(k),jw+jt(k)).eq.0)mk(iw+it(k),jw+jt(k))=kt(k)
          ii=is(kk)
          jj=js(kk)
          ib0=lq+1
c          call pltg(npb,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
          pause
210       continue
c
4000    continue
c       call pltg(nbp,xb,yb,nxg,nyg,xg,yg,mg,ixg,iyg)
       pause
        return
        end
c
        subroutine chek2(is,js,xg,yg,kk,ixg,iyg)
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension is(4),js(4),ij(4,4)
c
        do i=1,4
          k=0
          do j=1,4
          if(i.eq.2.and.j.eq.4)k=4
          if(i.eq.3.and.j.eq.3)k=4
          if(i.eq.4.and.j.eq.2)k=4
          ij(i,j)=i+j-1-k
          enddo
        enddo 
          x1=xg(is(ij(kk,1)),js(ij(kk,1)))   !グリッドの4つの頂点
          x2=xg(is(ij(kk,2)),js(ij(kk,2)))
          x3=xg(is(ij(kk,3)),js(ij(kk,3)))
          x4=xg(is(ij(kk,4)),js(ij(kk,4)))
          y1=yg(is(ij(kk,1)),js(ij(kk,1)))
          y2=yg(is(ij(kk,2)),js(ij(kk,2)))
          y3=yg(is(ij(kk,3)),js(ij(kk,3)))
          y4=yg(is(ij(kk,4)),js(ij(kk,4)))
c       対角線の傾きと切片
        if(dabs(x2-x4).lt.1d-12)then  !al=∞,am≠∞
          write(*,*)"al=∞,am≠∞"
1100      am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=x2
          if(xx.gt.dmin1(x1,x3) .and. xx.lt.dmax1(x1,x3))return
          write(*,*)"x1,y1修正 al=∞,am≠∞"
          x1=x1-(x3-x1)/2d0    !x1,y1を外側にずらす
          if(dabs(y1-y4).gt.dabs(y1-y2))y1=y1+(y4-y1)/2d0
          if(dabs(y1-y4).lt.dabs(y1-y2))y1=y1+(y2-y1)/2d0
          xg(is(ij(kk,1)),js(ij(kk,1)))=x1
          yg(is(ij(kk,1)),js(ij(kk,1)))=y1
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1100
        endif
c
        if(dabs(x1-x3).lt.1d-12)then    !al≠∞,am=∞
          write(*,*)"al≠∞,am=∞"
          al=(y2-y4)/(x2-x4)
          bl=(x2*y4-x4*y2)/(x2-x4)
          xx=x1
          yy=al*x1+bl
1200      if(yy.gt.dmin1(y1,y3) .and. yy.lt.dmax1(y1,y3))return
          write(*,*)"x1,y1修正 al≠∞,am=∞"
          y1=y1-(y3-y1)/2d0           !y1を外側にずらす.x1=const
          xg(is(ij(kk,1)),js(ij(kk,1)))=x1
          yg(is(ij(kk,1)),js(ij(kk,1)))=y1
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
	    write(*,*)"xx",xx
         pause
          goto 1200
        endif
c
          write(*,*)"al≠∞,am≠∞"
          al=(y2-y4)/(x2-x4)     !al≠∞,am≠∞
          bl=(x2*y4-x4*y2)/(x2-x4)
1000      am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=-(bl-bm)/(al-am)    !交点
          yy=am*xx+bm
          if(xx.gt.dmin1(x1,x3) .and. xx.lt.dmax1(x1,x3)  .and.  !交点が存在するか
     &       xx.gt.dmin1(x2,x4) .and. xx.lt.dmax1(x2,x4))return
          write(*,*)"x1,y1修正 al≠∞,am≠∞"
          if(dabs(y1-y4).ge.dabs(y1-y2))x1=x1-(x3-x1)/2d0    !交点が存在しなかった場合，頂点を外側にずらす
          if(dabs(y1-y4).lt.dabs(y1-y2))y1=y1-(y3-y1)/2d0
          xg(is(ij(kk,1)),js(ij(kk,1)))=x1
          yg(is(ij(kk,1)),js(ij(kk,1)))=y1
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
	    write(*,*)"xx",xx
          pause
          goto 1000
c
          end
c
          subroutine igrid(xb,yb,xg0,yg0,xg1,yg1,lq,xc,yc,icross)
          implicit double precision(a-h,o-z)
          dimension xb(1001),yb(1001)
          logical icross
c
          icross=.FALSE.
          x0=xb(lq)
          y0=yb(lq)
          x1=xb(lq+1)
          y1=yb(lq+1)
c         2つの線分(対角線)の傾きal,amと切片bl,bm
          if(dabs(xg1-xg0).lt.1d-12 .and. dabs(x1-x0).lt.1d-12)then !al=∞,am=∞
c            write(*,*)"al=∞,am=∞ in icross"
            if(dabs(xg1-x1).lt.1d-12 )then
		    if(dmax1(y0,y1).gt.dmin1(yg0,yg1) .or. 
     &           dmin1(y0,y1).lt.dmax1(yg0,yg1))then
              xx=(xg0+xg1)/2d0
			yy=(yg0+yg1)/2d0
			icross=.TRUE.
	        endif
	      endif
            goto 1000
          endif
c
          if(dabs(xg1-xg0).lt.1d-12)then      !al=∞,am≠∞
c            write(*,*)"al=∞,am≠∞ in icross"
            am=(y1-y0)/(x1-x0)
            bm=(x1*y0-x0*y1)/(x1-x0)
            xx=xg1
            yy=am*xg1+bm
            if(yy+1d-12.ge.dmin1(yg0,yg1) .and. 
     &         yy-1d-12.le.dmax1(yg0,yg1) .and.
     &         xx+1d-12.ge.dmin1(x0,x1) .and. 
     &         xx-1d-12.le.dmax1(x0,x1))icross=.TRUE.
            goto 1000
            endif
c
          if(dabs(x1-x0).lt.1d-12)then        !al≠∞,am=∞    
c           write(*,*)"al≠∞,am=∞ in icross"
            al=(yg1-yg0)/(xg1-xg0)
            bl=(xg1*yg0-xg0*yg1)/(xg1-xg0)
            xx=x1
            yy=al*x1+bl
            if(yy+1d-12.ge.dmin1(y0,y1) .and. 
     &         yy-1d-12.le.dmax1(y0,y1) .and.
     &         xx+1d-12.ge.dmin1(xg0,xg1) .and. 
     &         xx-1d-12.le.dmax1(xg0,xg1))icross=.TRUE.
            goto 1000
          endif
c		 
c          write(*,*)"al≠∞,am≠∞ in icross"
          al=(yg1-yg0)/(xg1-xg0)              !al≠∞,am≠∞
          bl=(xg1*yg0-xg0*yg1)/(xg1-xg0)
          am=(y1-y0)/(x1-x0)
          bm=(x1*y0-x0*y1)/(x1-x0)
          if(dabs(al-am).lt.1d-12)then  !al≠∞,am≠∞ & al=am    傾きが等しい
            if(dabs(bl-bm).lt.1d-12)icross=.TRUE.
            goto 1000
          endif
          xx=-(bl-bm)/(al-am)        !交点のx座標
	    yy=am*xx+bm
          if(xx+1d-12.ge.dmin1(x0,x1) .and. 
     &       xx-1d-12.le.dmax1(x0,x1) .and.   !交点が存在するか
     &       xx+1d-12.ge.dmin1(xg0,xg1) .and. 
     &       xx-1d-12.le.dmax1(xg0,xg1))icross=.TRUE.
c
1000    continue
        xc=xx
	  yc=yy
        return
        end
c
        subroutine chek3(ii,jj,xg,yg,kk,ixg,iyg)    !kkランダム
        implicit double precision(a-h,o-z)
        dimension xg(ixg,iyg),yg(ixg,iyg)
        dimension ik(4),jk(4),ij(4,4)
        data ik/0,1,1,0/,jk/0,0,1,1/
c
c        write(*,*)"in check"
        do i=1,4
          k=0
          do j=1,4
          if(i.eq.2.and.j.eq.4)k=4
          if(i.eq.3.and.j.eq.3)k=4
          if(i.eq.4.and.j.eq.2)k=4
          ij(i,j)=i+j-1-k
          enddo
          enddo 
        x1=xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))   !グリッドの4つの頂点
          x2=xg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))
          x3=xg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))
          x4=xg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))
          y1=yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))
          y2=yg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))
          y3=yg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))
          y4=yg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))
c       対角線の傾きと切片
        if(dabs(x1-x3).lt.1d-12 .and. dabs(x2-x4).lt.1d-12)then   !al=∞,am=∞
	    write(*,*)"グリッドが平らになっている"
	    write(*,*)"ii,jj",ii,jj
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
	    pause
	   endif
c	   
        if(dabs(x2-x4).lt.1d-12)then  !al=∞,am≠∞
c          write(*,*)"al=∞,am≠∞"
1100      am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=x2
          if(xx-1d-12.gt.dmin1(x1,x3) .and. 
     &       xx+1d-12.lt.dmax1(x1,x3))return
c
	    if(x1+1d-12.ge.dmin1(x3,xx) .and. 
     &       x1-1d-12.le.dmax1(x3,xx) .and.
     &       y1+1d-12.ge.dmin1(y3,yy) .and. 
     &       y1-1d-12.le.dmax1(y3,yy))then
		  x1=x1-(x3-x1)/2d0    !交点が存在しなかった場合，頂点を外側にずらす
            y1=y1-(y3-y1)/2d0
            xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
            yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=y1
            write(*,*)"x1,y1修正 al=∞,am≠∞"
	    elseif(x2+1d-12.ge.dmin1(x4,xx) .and.
     &           x2-1d-12.le.dmax1(x4,xx) .and.
     & 	       y2+1d-12.ge.dmin1(y4,yy) .and.
     &           y2-1d-12.le.dmax1(y4,yy))then
	      y2=y2-(y4-y2)/2d0
            yg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))=y2
            write(*,*)"y2修正 al=∞,am≠∞"
	    elseif(x3+1d-12.ge.dmin1(x1,xx) .and.
     &           x3-1d-12.le.dmax1(x1,xx) .and.
     &	       y3+1d-12.ge.dmin1(y1,yy) .and.
     &           y3-1d-12.le.dmax1(y1,yy))then
	      x3=x3-(x1-x3)/2d0
	      y3=y3-(y1-y3)/2d0
            xg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))=x3
            yg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))=y3
          write(*,*)"x3,y3修正 al=∞,am≠∞"
	    elseif(x4+1d-12.ge.dmin1(x2,xx) .and.
     &           x4-1d-12.le.dmax1(x2,xx) .and.
     &	       y4+1d-12.ge.dmin1(y2,yy) .and.
     &           y4-1d-12.le.dmax1(y2,yy))then
            y4=y4-(y2-y4)/2d0
            yg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))=y4
            write(*,*)"y4修正 al=∞,am≠∞"
	    endif
c
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1100
        endif
c
        if(dabs(x1-x3).lt.1d-12)then    !al≠∞,am=∞
c          write(*,*)"al≠∞,am=∞"
1200      al=(y2-y4)/(x2-x4)
          bl=(x2*y4-x4*y2)/(x2-x4)
          xx=x1
          yy=al*x1+bl
          if(yy-1d-12.gt.dmin1(y1,y3) .and. 
     &       yy+1d-12.lt.dmax1(y1,y3))return
c
	    if(x1+1d-12.ge.dmin1(x3,xx) .and. 
     &       x1-1d-12.le.dmax1(x3,xx) .and.
     &       y1+1d-12.ge.dmin1(y3,yy) .and. 
     &       y1-1d-12.le.dmax1(y3,yy))then
		  y1=y1-(y3-y1)/2d0    !交点が存在しなかった場合，頂点を外側にずらす
            yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
            write(*,*)"y1修正 al≠∞,am=∞"
	    elseif(x2+1d-12.ge.dmin1(x4,xx) .and.
     &           x2-1d-12.le.dmax1(x4,xx) .and.
     & 	       y2+1d-12.ge.dmin1(y4,yy) .and.
     &           y2-1d-12.le.dmax1(y4,yy))then
	      x2=x2-(x4-x2)/2d0
	      y2=y2-(y4-y2)/2d0
            xg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))=x2
            yg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))=y2
            write(*,*)"x2,y2修正 al≠∞,am=∞"
         elseif(x3+1d-12.ge.dmin1(x1,xx) .and.
     &           x3-1d-12.le.dmax1(x1,xx) .and.
     &	       y3+1d-12.ge.dmin1(y1,yy) .and.
     &           y3-1d-12.le.dmax1(y1,yy))then
	      y3=y3-(y1-y3)/2d0
            yg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))=x3
            write(*,*)"y3修正 al≠∞,am=∞"
	    elseif(x4+1d-12.ge.dmin1(x2,xx) .and.
     &           x4-1d-12.le.dmax1(x2,xx) .and.
     &	       y4+1d-12.ge.dmin1(y2,yy) .and.
     &           y4-1d-12.le.dmax1(y2,yy))then
	      x4=x4-(x2-x4)/2d0
            y4=y4-(y2-y4)/2d0
            xg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))=x4
            yg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))=y4
            write(*,*)"x4,y4修正 al≠∞,am=∞"
	    endif
c
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
          pause
          goto 1200
        endif
c
c          write(*,*)"al≠∞,am≠∞"
1000      al=(y2-y4)/(x2-x4)     !al≠∞,am≠∞
          bl=(x2*y4-x4*y2)/(x2-x4)
          am=(y1-y3)/(x1-x3)
          bm=(x1*y3-x3*y1)/(x1-x3)
          xx=-(bl-bm)/(al-am)    !交点
          yy=am*xx+bm
          if(xx-1d-12.gt.dmin1(x1,x3) .and. 
     &       xx+1d-12.lt.dmax1(x1,x3)  .and.  !交点が存在するか
     &       xx-1d-12.gt.dmin1(x2,x4) .and. 
     &       xx+1d-12.lt.dmax1(x2,x4))return
c
	    if(x1+1d-12.ge.dmin1(x3,xx) .and. 
     &       x1-1d-12.le.dmax1(x3,xx) .and.
     &       y1+1d-12.ge.dmin1(y3,yy) .and. 
     &       y1-1d-12.le.dmax1(y3,yy))then
		  x1=x1-(x3-x1)/2d0    !交点が存在しなかった場合，頂点を外側にずらす
            y1=y1-(y3-y1)/2d0
            xg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=x1
            yg(ii+ik(ij(kk,1)),jj+jk(ij(kk,1)))=y1
            write(*,*)"x1,y1修正 al≠∞,am≠∞"
	    elseif(x2+1d-12.ge.dmin1(x4,xx) .and.
     &           x2-1d-12.le.dmax1(x4,xx) .and.
     & 	       y2+1d-12.ge.dmin1(y4,yy) .and.
     &           y2-1d-12.le.dmax1(y4,yy))then
	      x2=x2-(x4-x2)/2d0
	      y2=y2-(y4-y2)/2d0
            xg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))=x2
            yg(ii+ik(ij(kk,2)),jj+jk(ij(kk,2)))=y2
            write(*,*)"x2,y2修正 al≠∞,am≠∞"
	    elseif(x3+1d-12.ge.dmin1(x1,xx) .and.
     &           x3-1d-12.le.dmax1(x1,xx) .and.
     &	       y3+1d-12.ge.dmin1(y1,yy) .and.
     &           y3-1d-12.le.dmax1(y1,yy))then
	      x3=x3-(x1-x3)/2d0
	      y3=y3-(y1-y3)/2d0
            xg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))=x3
            yg(ii+ik(ij(kk,3)),jj+jk(ij(kk,3)))=y3
            write(*,*)"x3,y3修正 al≠∞,am≠∞"
	    elseif(x4+1d-12.ge.dmin1(x2,xx) .and.
     &           x4-1d-12.le.dmax1(x2,xx) .and.
     &	       y4+1d-12.ge.dmin1(y2,yy) .and.
     &           y4-1d-12.le.dmax1(y2,yy))then
	      x4=x4-(x2-x4)/2d0
            y4=y4-(y2-y4)/2d0
            xg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))=x4
            yg(ii+ik(ij(kk,4)),jj+jk(ij(kk,4)))=y4
            write(*,*)"x4,y4修正 al≠∞,am≠∞"
	    endif
c
          write(*,*)"x1,y1",x1,y1
          write(*,*)"x2,y2",x2,y2
          write(*,*)"x3,y3",x3,y3
          write(*,*)"x4,y4",x4,y4
	    write(*,*)"xx,yy",xx,yy
          pause
          goto 1000
c
          end
