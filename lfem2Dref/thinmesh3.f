        subroutine thinmesh2(ne,np,nel,nr,nelrs,nelre,xx,yy,idx,irx,
     &                      new,jw,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine),irx(0:ine)
        dimension nelrs(21),nelre(21)
        dimension new(ine,3),jw(ind),mb(ind)
        dimension kelb(1001),keb(1001),ff(ind),nec(ine,3),kec(ine,3)
          dimension f0(50)     !
        common /thininf/ dthin(-9:9,-9:9),athin(-9:9,-9:9),
     &                  mthin(-9:9,-9:9)
c
        write(*,*)"in thinmesh" 
c       call plt(nel,ne,xx,yy,ind,ine)
c        call pltnt(nel,ne,xx,yy,idx,ind,ine)
c        pause
        do 100 lr=1,nr
          nels=nelrs(lr)
          nele=nelre(lr)
          id=idx(nels)
          ir=irx(nels)
c         do 120 i=nels,nele
c         do 120 j=1,3
c 120     nec(i,j)=0
c         do 140 i=1,np
c 140     ff(i)=-9999
c
          do 160 jr=0,nr
            jd=0
            if(jr.ne.0) jd=idx(nelrs(jr))
            mth=mthin(id,jd)
            df=dthin(id,jd)
            alph=athin(id,jd)-1d0
            write(*,*) lr,jr,mth
            if(mth.eq.0) goto 160

          do 120 i=1,nel
          do 120 j=1,3
          kec(i,j)=0
 120      nec(i,j)=0
          do 140 i=1,np
 140      ff(i)=-9999

              nels=nelrs(lr)
              nele=nelre(lr)
              call thinbound(ne,nel,nels,nele,idx,jd
     &                                  ,kb,kelb,keb,nec,kec,mb,ine,ind)
              write(*,*) 'kb=',kb
              if(kb.eq.0) goto 160
              call potential(ne,nels,nele,xx,yy,kb,kelb,keb,ff,ine,ind)
              do 180 m=1,mth
c              f0(m)=df*dble(m)*(1d0+alph*dble(m-1))
              f0(m)=df*dble(m)
 180          write(*,*)"f0(",m,")=",f0(m)
              call thindiv3(ne,nel,np,xx,yy,jw,mth,
     &            new,npw,nelw,ff,f0,nec,ine,ind,kb,kelb,keb,idx)
        write(*,*) 'np,npw=',np,npw
c       call plt(nelw,new,xx,yy,ind,ine)
c       pause
c       call plt(nel,ne,xx,yy,ind,ine)
c       pause
              do 200 i=nel,nele+1,-1
              idx(i+nelw)=idx(i)
              irx(i+nelw)=irx(i)
              do 200 j=1,3
 200          ne(i+nelw,j)=ne(i,j)
c
              do 220 i=1,nelw
              idx(nele+i)=id
              irx(nele+i)=ir
              do 220 j=1,3
              ne(i+nele,j)=new(i,j)
 220          nec(i+nele,j)=0
c
              do 240 i=lr+1,nr
 240          nelrs(i)=nelrs(i)+nelw
              do 260 i=lr,nr
 260          nelre(i)=nelre(i)+nelw
c
            nel=nel+nelw
            np=npw

c       call plt(nel,ne,xx,yy,ind,ine)
c        call pltnt(nel,ne,xx,yy,idx,ind,ine)
 160      continue
 100    continue
c        call plt(nel,ne,xx,yy,ind,ine)
        return
        end
c
        subroutine potential(ne,nels,nele,xx,yy,kb,kelb,keb,ff,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),ff(ind)
        dimension kelb(1001),keb(1001),js2(3)
        data js2/2,3,1/
c
        do 100 n=nels,nele
        do 110 m=1,3
          if(ff(ne(n,m)).ge.0d0) goto 110
          x=xx(ne(n,m))
          y=yy(ne(n,m))
          hmin=9999
          do 120 i=1,kb
            n1=ne(kelb(i),keb(i))
            n2=ne(kelb(i),js2(keb(i)))
            h=dmin1((x-xx(n1))**2+(y-yy(n1))**2,
     &                  (x-xx(n2))**2+(y-yy(n2))**2)
            if(h.lt.hmin) then
              imin=i
              hmin=h
            endif
 120      continue
          n1=ne(kelb(imin),keb(imin))
          n2=ne(kelb(imin),js2(keb(imin)))
          a=yy(n2)-yy(n1)
          b=-(xx(n2)-xx(n1))
          c=(xx(n2)-xx(n1))*yy(n1)-(yy(n2)-yy(n1))*xx(n1)
          ff(ne(n,m))=dabs(a*x+b*y+c)/sqrt(a*a+b*b)
 110    continue
 100    continue
        return
        end
c
        subroutine thindiv3(ne,nel,np,xx,yy,jw,mth,
     &      new,npw,nelw,ff,f0,nec,ine,ind,kb,kelb,keb,idx)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),jw(ind),ff(ind)
        dimension new(ine,3),idx(0:ine)
        dimension js2(3),js3(3)
        dimension nss(ine,3),nec(ine,3)
        dimension kelb(1001),keb(1001)
        dimension f0(50),na(50),nb(50),xa(50),ya(50),xb(50),yb(50)
          dimension nc(50),xc(50),yc(50)
          dimension k1(50),ik(50),j1(50),j2(50),j3(50)
        data js2/2,3,1/
        data js3/3,1,2/   
        data eps/-1d-3/
c
        write(*,*) 'in thindiv np=',np,"ind=",ind,"mth=",mth
        do 20 i=1,ind
 20       jw(i)=i
        do 40 i=1,ine
        do 40 j=1,3
 40     nss(i,j)=0
c
c       find start point(境界が繋がっている場合は適用できない)
        npw=np
        nelw=0
        do 100 i=1,kb
          nn=ne(kelb(i),js2(keb(i)))
          ichk=0
          do 200 j=1,kb
            nm=ne(kelb(j),keb(j))
            if(nn.eq.nm)ichk=1
 200     continue
         if(ichk.eq.0)goto 800
 100    continue
 800    lastkb=i
        imend=kelb(lastkb)
        kelb(kb+1)=kelb(lastkb)
        keb(kb+1)=js2(keb(lastkb))
c
        write(*,*)"last element im=",imend,"kb=",lastkb
        write(*,*)"xx=",xx(ne(kelb(lastkb),keb(lastkb)))
        write(*,*)"yy=",yy(ne(kelb(lastkb),keb(lastkb)))
        write(*,*)"xx=",xx(ne(kelb(kb+1),keb(kb+1)))
        write(*,*)"yy=",yy(ne(kelb(kb+1),keb(kb+1)))
c        pause
        ikb=lastkb
        do 300 it=1,5000
          do 400 i=1,kb
            if(ne(kelb(i),js2(keb(i))).eq.
     &           ne(kelb(ikb),keb(ikb)))then
              ikb=i
              goto 300
            endif
 400      continue
          goto 700
 300    continue
 700    continue
        imin=kelb(ikb)
        write(*,*)"imin=",imin,"kb=",ikb
c
c       about first element
          im=imin
          call findk1(ne,kelb,keb,kb,im,j1(1),k1(1),ine)
        j2(1)=js2(j1(1))
        j3(1)=js3(j1(1))
          n1=ne(im,j1(1))
        n2=ne(im,j2(1))
        n3=ne(im,j3(1))
        write(*,*)"first element im=",im,"k1(1)=",k1(1)
        write(*,*)"xx,yy",xx(n1),yy(n1)
        write(*,*)"xx,yy",xx(n2),yy(n2)
        write(*,*)"xx,yy",xx(n3),yy(n3)
        write(*,*)"ff0,ff",f0(1),f0(2),ff(n1),ff(n2),ff(n3)
          if(k1(1).eq.1)then    !the case of k1(1)=1 at first element
            write(*,*)"k1(1)",1
            do m=1,mth
            xb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(xx(n3)-xx(n1))+xx(n1)
            yb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(yy(n3)-yy(n1))+yy(n1)
            enddo
c
c
            if(ff(n1)-ff(n3).le.f0(mth))then
            write(*,*)"ff(",n1,")-ff(",n3,")=",ff(n1)-ff(n3)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1)
              ddtt=2d-1
              if(mth.gt.4)ddtt=1d-1
              if(mth.gt.9)then
                write(*,*)"mth is too many"
                stop
              endif
              do m=1,mth
                xb(m)=(1d0-ddtt*dble(m))*(xx(n3)-xx(n1))+xx(n1)
                yb(m)=(1d0-ddtt*dble(m))*(yy(n3)-yy(n1))+yy(n1)
              enddo
            endif
c
c
          do m=1,mth
          write(*,*)"xb,yb at fst",xb(m),yb(m)
            nelw=nelw+1
            npw=npw+1
            nb(m)=npw
            xx(nb(m))=xb(m)
            yy(nb(m))=yb(m)
            new(nelw,1)=n2
            if(m.eq.1)then
              new(nelw,2)=n3
              new(nelw,3)=nb(1)
            else
              new(nelw,2)=nb(m-1)
              new(nelw,3)=nb(m)
            endif
          enddo
          nss(im,j3(1))=nb(mth)
          write(*,*)"ended first element"
          write(*,*)"next element nec(",im,",",j3(1),")=",nec(im,j3(1))
          im=nec(im,j3(1))
c
        elseif(k1(1).eq.2)then   !the case of k1(1)=2 at first emelent
          write(*,*)"k1(1)",2
          do m=1,mth
            xb(m)=f0(m)/(ff(n1)-ff(n2))*(xx(n1)-xx(n2))+xx(n2)
            yb(m)=f0(m)/(ff(n1)-ff(n2))*(yy(n1)-yy(n2))+yy(n2)
          enddo
c
c
          if(ff(n2)-ff(n1).le.f0(mth))then
            write(*,*)"ff(",n2,")-ff(",n1,")=",ff(n2)-ff(n1)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1)
            ddtt=2d-1
            if(mth.gt.4)ddtt=1d-1
            if(mth.gt.9)then
              write(*,*)"mth is too many"
              stop
            endif
            do m=1,mth
              xb(m)=ddtt*dble(m)*(xx(n2)-xx(n1))+xx(n1)
              yb(m)=ddtt*dble(m)*(yy(n2)-yy(n1))+yy(n1)
            enddo
          endif
c
c
          write(*,*)"xb,yb at fst",xb(m),yb(m)
          do m=1,mth
            nelw=nelw+1
            npw=npw+1
            nb(m)=npw
            xx(nb(m))=xb(m)
            yy(nb(m))=yb(m)
            new(nelw,1)=n3
            if(m.eq.1)then
              new(nelw,2)=n1
              new(nelw,3)=nb(1)
            else
              new(nelw,2)=nb(m-1)
              new(nelw,3)=nb(m)
            endif
          enddo
          nss(im,j1(1))=nb(mth)
          write(*,*)"ended first element"
          write(*,*)"next element nec(",im,",",j1(1),")=",nec(im,j1(1))
          im=nec(im,j1(1))
        endif
        do m=1,mth
        na(m)=nb(m)
        xa(m)=xb(m)
        ya(m)=yb(m)
        enddo
        write(*,*)"xx,yy",xx(ne(im,j1(1))),yy(ne(im,j1(1)))
        write(*,*)"xx,yy",xx(ne(im,js2(j1(1)))),yy(ne(im,js2(j1(1))))
        write(*,*)"xx,yy",xx(ne(im,js3(j1(1)))),yy(ne(im,js3(j1(1))))
c        pause
c
c       larger elements than first element
c
        do 500 it=1,100000
c          write(*,*)"started im=",im
            call findk1(ne,kelb,keb,kb,im,j1(1),k1(1),ine)
            j2(1)=js2(j1(1))
            j3(1)=js3(j1(1))
c
          if(k1(1).eq.1)then
              im2=nec(im,j3(1))
              call findk1(ne,kelb,keb,kb,im2,j1(2),k1(2),ine)
              j2(2)=js2(j1(2))
              j3(2)=js3(j1(2))
            if(k1(2).eq.1)goto 1000  !deal with the case of k1(1)=1,k1(2)=1
            if(k1(2).eq.2)goto 2000  !deal with the case of k1(1)=1,k1(2)=2
            write(*,*)"wrong value in k1(2) in the case of k1(1)=1"
              stop
          elseif(k1(1).eq.2)then
              im2=nec(im,j1(1))
              call findk1(ne,kelb,keb,kb,im2,j1(2),k1(2),ine)     
              j2(2)=js2(j1(2))
              j3(2)=js3(j1(2))
            if(k1(2).eq.1)goto 3000  !deal with the case of k1(1)=2,k1(2)=1
            if(k1(2).eq.2)goto 4000  !deal with the case of k1(1)=2,k1(2)=2
            write(*,*)"wrong value in k1(2) in the case of k1(1)=2"
              stop
          else
            write(*,*)"wrong value in k1(1)"
            pause
          endif
c
 1000     continue     !the case of k1(1)=1,k1(2)=1
          write(*,*)"1000"
          n1=ne(im,j1(1))
          n2=ne(im,j2(1))
          n3=ne(im,j3(1))
          n4=ne(im2,j3(2))
          do m=1,mth
            xb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(xx(n3)-xx(n1))+xx(n1)
            yb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(yy(n3)-yy(n1))+yy(n1)
            xc(m)=(f0(m)-ff(n1))/(ff(n4)-ff(n1))*(xx(n4)-xx(n1))+xx(n1)
            yc(m)=(f0(m)-ff(n1))/(ff(n4)-ff(n1))*(yy(n4)-yy(n1))+yy(n1)
            enddo
c          if(ff(n1)-ff(n3).lt.ff(n1)-f0(mth) .or.
c     &       ff(n1)-ff(n4).lt.ff(n1)-f0(mth))then
          if(ff(n1)-ff(n3).le.f0(mth)-eps .or.
     &       ff(n1)-ff(n4).le.f0(mth)-eps)then
            write(*,*)"ff(",n1,")-ff(",n3,")=",ff(n1)-ff(n3)
            write(*,*)"ff(",n1,")-ff(",n4,")=",ff(n1)-ff(n3)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1),"k1(2)=",k1(2)
              ddtt=2d-1
              if(mth.gt.4)ddtt=1d-1
              if(mth.gt.9)then
                write(*,*)"mth is too many"
                stop
              endif
              do m=1,mth
                if(ff(n1)-ff(n3).le.f0(mth)-eps)then
                  xb(m)=(1d0-ddtt*dble(m))*(xx(n3)-xx(n1))+xx(n1)
                yb(m)=(1d0-ddtt*dble(m))*(yy(n3)-yy(n1))+yy(n1)
                endif
                if(ff(n1)-ff(n4).le.f0(mth)-eps)then
                  xc(m)=(1d0-ddtt*dble(m))*(xx(n4)-xx(n1))+xx(n1)
                  yc(m)=(1d0-ddtt*dble(m))*(yy(n4)-yy(n1))+yy(n1)
                endif
              enddo
            endif
c          if(idx(nec(im,j3(1))).eq.-1)goto 5000    !the end of element
c          if(idx(nec(im2,j3(2))).eq.-1)goto 6000    !the end of element
          if(im.eq.imend)goto 5000      !the end of element
          if(im2.eq.imend)goto 6000     !the end of element
          do m=1,mth
              npw=npw+1
              nb(m)=npw
              xx(nb(m))=xb(m)
              yy(nb(m))=yb(m)
            npw=npw+1
            nc(m)=npw
            xx(nc(m))=xc(m)
            yy(nc(m))=yc(m)
            enddo
          nss(im,j2(1))=na(mth)      !downer element
          nss(im,j3(1))=nb(mth)
            nss(im2,j2(2))=nb(mth)
            nss(im2,j3(2))=nc(mth)
          ik(1)=1               !upper element
          if((xa(mth)-xx(n4))**2+(ya(mth)-yy(n4))**2.gt.
     &               (xc(mth)-xx(n2))**2+(yc(mth)-yy(n2))**2)ik(1)=2
          if(ik(1).eq.1)then
              do m=1,mth
              if(m.eq.1)then
              nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=n3
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=na(1)
                new(nelw,2)=n3
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=n3
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=nb(1)
                new(nelw,2)=n4
                new(nelw,3)=nc(1)       
                  else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
              endif
              enddo
          elseif(ik(1).eq.2)then
              do m=1,mth
              if(m.eq.1)then
              nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=nb(1)
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=n3
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=n3
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=nb(1)
                new(nelw,2)=n4
                new(nelw,3)=nc(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=nc(m)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
              endif
              enddo
          endif
c
c         write(*,*)"end im=",im,"im2=",im2,"k1(1)=",k1(1),"k1(2)=",k1(2)
c         write(*,*)"next element nec(",im2,",",j3(2),")=",nec(im2,j3(2))
          im=nec(im2,j3(2))
            do m=1,mth
          na(m)=nc(m)
          xa(m)=xc(m)
          ya(m)=yc(m)
            enddo
          goto 500
c
 2000     continue     !the case of k1(1)=1,k1(2)=2
          write(*,*)"2000"
          n1=ne(im,j1(1))
          n2=ne(im,j2(1))
          n3=ne(im,j3(1))
            n4=ne(im2,j2(2))
            do m=1,mth
            xb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(xx(n3)-xx(n1))+xx(n1)
            yb(m)=(f0(m)-ff(n1))/(ff(n3)-ff(n1))*(yy(n3)-yy(n1))+yy(n1)
            xc(m)=f0(m)/(ff(n4)-ff(n3))*(xx(n4)-xx(n3))+xx(n3)
            yc(m)=f0(m)/(ff(n4)-ff(n3))*(yy(n4)-yy(n3))+yy(n3)
            enddo
c          if(ff(n1)-ff(n3).lt.ff(n1)-f0(mth) .or.
c     &       ff(n4)-ff(n3).lt.f0(mth))then
          if(ff(n1)-ff(n3).le.f0(mth)-eps .or.
     &       ff(n4)-ff(n3).le.f0(mth)-eps)then
            write(*,*)"ff(",n1,")-ff(",n3,")=",ff(n1)-ff(n3)
            write(*,*)"ff(",n4,")-ff(",n3,")=",ff(n4)-ff(n3)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1),"k1(2)=",k1(2)
              ddtt=2d-1
              if(mth.gt.4)ddtt=1d-1
              if(mth.gt.9)then
                write(*,*)"mth is too many"
                stop
              endif
              do m=1,mth
                if(ff(n1)-ff(n3).le.f0(mth)-eps)then
                  xb(m)=(1d0-ddtt*dble(m))*(xx(n3)-xx(n1))+xx(n1)
                yb(m)=(1d0-ddtt*dble(m))*(yy(n3)-yy(n1))+yy(n1)
                endif
                if(ff(n4)-ff(n3).le.f0(mth)-eps)then
                  xc(m)=ddtt*dble(m)*(xx(n4)-xx(n3))+xx(n3)
                  yc(m)=ddtt*dble(m)*(yy(n4)-yy(n3))+yy(n3) 
                endif
              enddo
            endif
c          if(idx(nec(im,j3(1))).eq.-1)goto 5000    !the end of element
c          if(idx(nec(im2,j1(2))).eq.-1)goto 6000    !the end of element
          if(im.eq.imend)goto 5000    !the end of element
          if(im2.eq.imend)goto 6000    !the end of element
            do m=1,mth
          npw=npw+1
          nc(m)=npw
          xx(nc(m))=xc(m)
          yy(nc(m))=yc(m)
            enddo
          nss(im,j2(1))=na(mth)      !downer element
          nss(im,j3(1))=nc(mth)
            nss(im2,j1(2))=nc(mth)
          ik(1)=1               !upper element
          if((xa(mth)-xx(n3))**2+(ya(mth)-yy(n3))**2.gt.
     &               (xc(mth)-xx(n2))**2+(yc(mth)-yy(n2))**2)ik(1)=2
          if(ik(1).eq.1)then
              do m=1,mth
              if(m.eq.1)then
              nelw=nelw+1
              new(nelw,1)=n2
              new(nelw,2)=n3
              new(nelw,3)=na(1)
              nelw=nelw+1
              new(nelw,1)=na(1)
              new(nelw,2)=n3
              new(nelw,3)=nc(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
              endif
              enddo
          elseif(ik(1).eq.2)then
              do m=1,mth
              if(m.eq.1)then
              nelw=nelw+1
              new(nelw,1)=n2
              new(nelw,2)=n3
              new(nelw,3)=nc(1)
              nelw=nelw+1
              new(nelw,1)=n2
              new(nelw,2)=nc(1)
              new(nelw,3)=na(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nc(m)
                new(nelw,3)=na(m)
              endif
              enddo
          endif
c
c         write(*,*)"end im=",im,"im2=",im2,"k1(1)=",k1(1),"k1(2)=",k1(2)
c         write(*,*)"next element nec(",im2,",",j1(2),")=",nec(im2,j1(2))
          im=nec(im2,j1(2))
            do m=1,mth
          na(m)=nc(m)
          xa(m)=xc(m)
          ya(m)=yc(m)
            enddo
          goto 500
c
 3000     continue     !the case of k1(1)=2,k1(2)=1         
          write(*,*)"3000"
          n1=ne(im,j1(1))
          n2=ne(im,j2(1))
          n3=ne(im,j3(1))
          n4=ne(im2,j3(2))
          do m=1,mth
            xb(m)=f0(m)/(ff(n2)-ff(n1))*(xx(n2)-xx(n1))+xx(n1)
            yb(m)=f0(m)/(ff(n2)-ff(n1))*(yy(n2)-yy(n1))+yy(n1)
            xc(m)=(f0(m)-ff(n2))/(ff(n4)-ff(n2))*(xx(n4)-xx(n2))+xx(n2)
            yc(m)=(f0(m)-ff(n2))/(ff(n4)-ff(n2))*(yy(n4)-yy(n2))+yy(n2)
c            xb(m)=(f0(m)-ff(n1))/(ff(n2)-ff(n1))*(xx(n2)-xx(n1))+xx(n1)
c            yb(m)=(f0(m)-ff(n1))/(ff(n2)-ff(n1))*(yy(n2)-yy(n1))+yy(n1)
          enddo
c          if(ff(n2)-ff(n1).lt.f0(mth) .or.
c     &       ff(n2)-ff(n4).lt.ff(n2)-f0(mth))then
          if(ff(n2)-ff(n1).le.f0(mth)-eps .or.
     &       ff(n2)-ff(n4).le.f0(mth)-eps)then
            write(*,*)"ff(",n2,")-ff(",n1,")=",ff(n2)-ff(n1)
            write(*,*)"ff(",n2,")-ff(",n4,")=",ff(n2)-ff(n4)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1),"k1(2)=",k1(2)
              ddtt=2d-1
              if(mth.gt.4)ddtt=1d-1
              if(mth.gt.9)then
                write(*,*)"mth is too many"
                stop
              endif
              do m=1,mth
              if(ff(n2)-ff(n1).le.f0(mth)-eps)then
                xb(m)=ddtt*dble(m)*(xx(n2)-xx(n1))+xx(n1)
              yb(m)=ddtt*dble(m)*(yy(n2)-yy(n1))+yy(n1)
              endif
              if(ff(n2)-ff(n4).le.f0(mth)-eps)then
                xc(m)=(1d0-ddtt*dble(m))*(xx(n4)-xx(n2))+xx(n2)
                yc(m)=(1d0-ddtt*dble(m))*(yy(n4)-yy(n2))+yy(n2) 
              endif
              enddo
          endif
       
         do m=1,mth
           write(*,*)"xa,ya",xa(m),ya(m)
           write(*,*)"xb,yb",xb(m),yb(m)
           write(*,*)"xc,yc",xc(m),yc(m)
         enddo
c        if(idx(nec(im,j1(1))).eq.-1)goto 5000
c        if(idx(nec(im2,j3(2))).eq.-1)goto 6000
          if(im.eq.imend)goto 5000     !the end of element
          if(im2.eq.imend)goto 6000    !the end of element
            do m=1,mth
          npw=npw+1
          nc(m)=npw
          xx(nc(m))=xc(m)
          yy(nc(m))=yc(m)
            enddo
          ik(1)=1
          if((xx(n4)-xa(mth))**2+(yy(n4)-ya(mth))**2.gt.
     &            (xx(n1)-xc(mth))**2+(yy(n1)-yc(mth))**2)ik(1)=2
          if(ik(1).eq.1)then       !upper elements
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=na(1)
                new(nelw,2)=n4
                new(nelw,3)=nc(1)
              else
                nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=nc(m-1)
              new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
              endif
              enddo
          elseif(ik(1).eq.2)then
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=nc(1)
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nc(1)
                new(nelw,3)=na(1)
              else
                    nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=nc(m-1)
              new(nelw,3)=nc(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nc(m)
                new(nelw,3)=na(m)
              endif
              enddo
          endif
c
          nss(im,j1(1))=na(mth)     !basic element
          nss(im2,j2(2))=na(mth)
            nss(im2,j3(2))=nc(mth)
c
c         write(*,*)"end im=",im,"im2=",im2,"k1(1)=",k1(1),"k1(2)=",k1(2)
c         write(*,*)"next element nec(",im2,",",j3(2),")=",nec(im2,j3(2))
          im=nec(im2,j3(2))
          do m=1,mth
          na(m)=nc(m)
          xa(m)=xc(m)
          ya(m)=yc(m)
            enddo
          goto 500
c
 4000     continue     !the case of k1(1)=2,k1(2)=2         
          write(*,*)"4000"
          n1=ne(im,j1(1))
          n2=ne(im,j2(1))
          n3=ne(im,j3(1))
          n4=ne(im2,j2(2))
          do m=1,mth
            xb(m)=f0(m)/(ff(n2)-ff(n1))*(xx(n2)-xx(n1))+xx(n1)
            yb(m)=f0(m)/(ff(n2)-ff(n1))*(yy(n2)-yy(n1))+yy(n1)
              xc(m)=f0(m)/(ff(n4)-ff(n1))*(xx(n4)-xx(n1))+xx(n1)
              yc(m)=f0(m)/(ff(n4)-ff(n1))*(yy(n4)-yy(n1))+yy(n1)
c            xb(m)=(f0(m)-ff(n1))/(ff(n2)-ff(n1))*(xx(n2)-xx(n1))+xx(n1)
c            yb(m)=(f0(m)-ff(n1))/(ff(n2)-ff(n1))*(yy(n2)-yy(n1))+yy(n1)
          enddo
          if(ff(n2)-ff(n1).le.f0(mth)-eps .or.
     &       ff(n4)-ff(n1).le.f0(mth)-eps)then
            write(*,*)"ff(",n2,")-ff(",n1,")=",ff(n2)-ff(n1)
            write(*,*)"ff(",n4,")-ff(",n1,")=",ff(n4)-ff(n1)
            write(*,*)"f0(",mth,")=",f0(mth)
            write(*,*)"k1(1)=",k1(1),"k1(2)=",k1(2)
            ddtt=2d-1
            if(mth.gt.4)ddtt=1d-1
            if(mth.gt.9)then
              write(*,*)"mth is too many"
              stop
            endif
            do m=1,mth
              if(ff(n2)-ff(n1).le.f0(mth)-eps)then
                xb(m)=ddtt*dble(m)*(xx(n2)-xx(n1))+xx(n1)
              yb(m)=ddtt*dble(m)*(yy(n2)-yy(n1))+yy(n1)
              endif
              if(ff(n4)-ff(n1).le.f0(mth)-eps)then
                xc(m)=ddtt*dble(m)*(xx(n4)-xx(n1))+xx(n1)
                yc(m)=ddtt*dble(m)*(yy(n4)-yy(n1))+yy(n1)
              endif
            enddo
          endif
c          if(idx(nec(im,j1(1))).eq.-1)goto 5000
c          if(idx(nec(im2,j1(2))).eq.-1)goto 6000
          if(im.eq.imend)goto 5000   !the end of element
          if(im2.eq.imend)goto 6000  !the end of element
            do m=1,mth
            npw=npw+1
            nc(m)=npw
            xx(nc(m))=xc(m)
            yy(nc(m))=yc(m)
            enddo
          ik(1)=1
          if((xc(1)-xa(mth))**2+(yc(1)-ya(mth))**2.gt.
     &            (xa(1)-xc(mth))**2+(ya(1)-yc(mth))**2)ik(1)=2
          if(ik(1).eq.1)then       !upper elements
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nc(1)
                new(nelw,3)=na(1)
              else
               nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=nc(m-1)
              new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nc(m-1)
                new(nelw,3)=nc(m)
              endif
              enddo
          elseif(ik(1).eq.2)then
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nc(1)
                new(nelw,3)=na(1)
              else
                nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=nc(m-1)
              new(nelw,3)=nc(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nc(m)
                new(nelw,3)=na(m)
              endif
              enddo
          endif
c
          nelw=nelw+1
            new(nelw,1)=na(mth)
            new(nelw,2)=nc(mth)
            new(nelw,3)=n2
          nss(im,j1(1))=na(mth)     !basic element
          nss(im2,j1(2))=nc(mth)
c
c         write(*,*)"end im=",im,"im2=",im2,"k1(1)=",k1(1),"k1(2)=",k1(2)
c         write(*,*)"next element nec(",im2,",",j1(2),")=",nec(im2,j1(2))
          im=nec(im2,j1(2))
          do m=1,mth
          na(m)=nc(m)
          xa(m)=xc(m)
          ya(m)=yc(m)
            enddo
          goto 500
c
 500      continue
c
 5000   continue        !the end of element
        write(*,*)"the latest element im=",im,"k1(1)=",k1(1)
        if(k1(1).eq.1)then
          nss(im,j2(1))=na(mth)
            do m=1,mth
            if(m.eq.1)then
            nelw=nelw+1
            new(nelw,1)=n2
            new(nelw,2)=n3
            new(nelw,3)=na(1)
            else
              nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=n3
              new(nelw,3)=na(m)
            endif
            enddo
            goto 7000
          elseif(k1(1).eq.2)then
            nss(im,j1(1))=na(mth)
            do m=1,mth
            if(m.eq.1)then
              nelw=nelw+1
              new(nelw,1)=n1
              new(nelw,2)=n2
              new(nelw,3)=na(1)
            else
              nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=n2
              new(nelw,3)=na(m)
            endif
            enddo
            goto 7000
          else
            write(*,*)"error at the latest element im=",im
          endif
c
 6000   continue
        write(*,*)"the latest element im2=",im2
          write(*,*)"k1(1)=",k1(1),"k1(2)=",k1(2)
          if(k1(1).eq.1)goto 6300    !deal with the case of k1(1)=1
          if(k1(1).eq.2)goto 6600        !deal with the case of k1(1)=2
          write(*,*)"error at the latest element im=",im
          stop
c
 6300     if(k1(2).eq.1)then      !the case of k1(1)=1,k1(2)=1
          do m=1,mth
              npw=npw+1
              nb(m)=npw
              xx(nb(m))=xb(m)
              yy(nb(m))=yb(m)
            enddo
          ik(1)=1
            if((xx(n3)-xa(mth))**2+(yy(n3)-ya(mth))**2.gt.
     &       (xx(n2)-xb(mth))**2+(yy(n2)-yb(mth))**2)ik(1)=2
            nss(im,j2(1))=na(mth)
            nss(im,j3(1))=nb(mth)
            nss(im2,j2(2))=nb(mth)
            if(ik(1).eq.1)then     !ik(1)=1
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=n3
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=na(1)
                new(nelw,2)=n3
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=n3
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            elseif(ik(1).eq.2)then              !ik(1)=2
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=nb(1)
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=n2
                new(nelw,2)=n3
                new(nelw,3)=nb(1)
                nelw=nelw+1
                new(nelw,1)=n3
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            endif
                goto 7000
          elseif(k1(2).eq.2)then     !the case of k1(1)=1,k1(2)=2
            nss(im,j2(1))=na(mth)
            do m=1,mth
            if(m.eq.1)then
              nelw=nelw+1
              new(nelw,1)=n2
              new(nelw,2)=n3
              new(nelw,3)=na(1)
            else
              nelw=nelw+1
              new(nelw,1)=na(m-1)
              new(nelw,2)=n3
              new(nelw,3)=na(m)
            endif
            enddo
            goto 7000
          else
            write(*,*)"error at the latest element im2=",im2
          endif
c
 6600     if(k1(2).eq.1)then         !the case of k1(1)=2,k1(2)=1
          do m=1,mth
              npw=npw+1
              nb(m)=npw
              xx(nb(m))=xb(m)
              yy(nb(m))=yb(m)
            enddo
          ik(1)=1
            if((xb(1)-xx(n3))**2+(yb(1)-yy(n3))**2.gt.
     &       (xa(1)-xx(n2))**2+(yb(1)-yy(n2))**2)ik(1)=2
            if(ik(1).eq.1)then     !ik(1)=1
              nss(im,j1(1))=nb(mth)
              nss(im2,j2(2))=nb(mth)
              nelw=nelw+1
              new(nelw,1)=na(mth)
              new(nelw,2)=nb(mth)
              new(nelw,3)=n3
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nb(1)
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            elseif(ik(1).eq.2)then                 !ik(1)=2
              nss(im,j1(1))=na(mth)
              nss(im2,j2(2))=nb(mth)
              nelw=nelw+1
              new(nelw,1)=na(mth)
              new(nelw,2)=nb(mth)
              new(nelw,3)=n2
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nb(1)
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            endif
            goto 7000
          elseif(k1(2).eq.2)then     !the case of k1(1)=2,k1(2)=2,ik(1)=2
            do m=1,mth
              npw=npw+1
              nb(m)=npw
              xx(nb(m))=xb(m)
              yy(nb(m))=yb(m)
            enddo
            ik(1)=1
            if((xb(1)-xx(n3))**2+(yb(1)-yy(n3))**2.gt.
     &       (xa(1)-xx(n1))**2+(ya(1)-yy(n2))**2)ik(1)=2
            if(ik(1).eq.1)then      !ik(1)=1
              nss(im,j1(1))=nb(mth)
              nss(im2,j1(2))=nb(mth)
                  nelw=nelw+1
                  new(nelw,1)=na(mth)
                  new(nelw,2)=nb(mth)
                  new(nelw,3)=n3
                  do m=1,mth
                  if(m.eq.1)then
                    nelw=nelw+1
                        new(nelw,1)=n1
                        new(nelw,2)=nb(1)
                        new(nelw,3)=na(1)
                        nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            elseif(ik(1).eq.2)then  !ik(1)=2
              nss(im,j1(1))=na(mth)
              nss(im2,j1(2))=nb(mth)
              nelw=nelw+1
              new(nelw,1)=na(mth)
              new(nelw,2)=nb(mth)
              new(nelw,3)=n2
              do m=1,mth
              if(m.eq.1)then
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=nb(1)
                new(nelw,3)=na(1)
                nelw=nelw+1
                new(nelw,1)=n1
                new(nelw,2)=n4
                new(nelw,3)=nb(1)
              else
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m)
                new(nelw,3)=na(m)
                nelw=nelw+1
                new(nelw,1)=na(m-1)
                new(nelw,2)=nb(m-1)
                new(nelw,3)=nb(m)
                nelw=nelw+1
                new(nelw,1)=nb(m-1)
                new(nelw,2)=n4
                new(nelw,3)=nb(m)
              endif
              enddo
            endif
            goto 7000
          else
            write(*,*)"error at the latest element im2=",im2
            stop
          endif
c
 7000   do 600 i=1,nel
        do 600 j=1,3
 600    if(nss(i,j).ne.0)ne(i,j)=nss(i,j)
c
        return
        end
c
        subroutine thinbound(ne,nel,nels,nele,idx,jd,
     &                          kb,kelb,keb,nec,kec,mb,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),nec(ine,3),kec(ine,3),idx(0:ine)
        dimension kelb(1001),keb(1001),mb(ind)
        dimension js2(3)
        data js2/2,3,1/
c
        write(*,*)"in thinbound"   !
c        do 120 i=nels,nele
        do 120 i=1,nel     !
        do 130 j=1,3
          if(nec(i,j).ne.0) goto 130
          do 140 n=nels,nele    !1,nel
            if(n.eq.i) goto 130
            do 150 m=1,3
              if(ne(i,j).eq.ne(n,js2(m)).and.
     &                ne(i,js2(j)).eq.ne(n,m)) then
              nec(i,j)=n
              kec(i,j)=m
              nec(n,m)=i
              kec(n,m)=j
              goto 130
              endif
 150        continue
 140      continue
 130    continue
 120    continue
c
        do 180 i=1,ind
 180    mb(i)=0
        kb=0
        do 200 i=nels,nele
        do 200 j=1,3
!         if(idx(nec(i,j)).eq.jd) then
          if(idx(nec(i,j)).eq.0) then
            kb=kb+1
            kelb(kb)=i
            keb(kb)=j
            mb(ne(i,j))=1
            mb(ne(i,js2(j)))=1
          endif
 200    continue
c
        return
        end
c
        subroutine findk1(ne,kelb,keb,kb,iim,jj1,kk1,ine)
          implicit double precision(a-h,o-z)
          dimension ne(ine,3)
          dimension kelb(1001),keb(1001)
          dimension jchk(3)
c
          m=0          !deciding j1,j2,j3 about im element
          do j=1,3
          jchk(j)=0
          enddo
          do i=1,kb+1
            if(ne(iim,1).eq.ne(kelb(i),keb(i)))then
              m=m+1
              jchk(1)=1
            elseif(ne(iim,2).eq.ne(kelb(i),keb(i)))then
              m=m+1
              jchk(2)=1
            elseif(ne(iim,3).eq.ne(kelb(i),keb(i)))then
              m=m+1
              jchk(3)=1
            endif
          enddo
c
          if(m.eq.2)then
            kk1=1
            do j=1,3
              if(jchk(j).eq.0)jj1=j
            enddo
          elseif(m.eq.1)then
            kk1=2
            do j=1,3
              if(jchk(j).eq.1)jj1=j
            enddo
          else
            write(*,*)"error in thindiv m=",m
            pause
          endif
c
          return
            end
