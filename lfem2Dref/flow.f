        subroutine flow(nel,nbw,sa,sf,ne,xx,yy,uvp,idx,dt,ind,ine,ibw)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine)
        dimension sa(ind,ibw),sf(ind),uvp(ind)
        dimension x(3),y(3),b(3),c(3),zkesu(3,3)
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        data zkesu/2d0,1d0,1d0,1d0,2d0,1d0,1d0,1d0,2d0/
c
c         open(1,file="neflow.res")
c         rewind(1)
c         do i=1,nel
c         do j=1,3
c         write(1,*)"ne(",i,",",j,")=",ne(i,j)
c         enddo
c         enddo
c         close(1)
c      write(*,*) 'flow/nel=',nel,dt,nbw
        do 200 im=1,nel
        id=idx(im)
          do 210 i=1,3
          x(i)=xx(ne(im,i))
 210      y(i)=yy(ne(im,i))
          b(1)=y(2)-y(3)
          b(2)=y(3)-y(1)
          b(3)=y(1)-y(2)
          c(1)=x(3)-x(2)
          c(2)=x(1)-x(3)
          c(3)=x(2)-x(1)
          ar=((x(1)-x(2))*(y(1)-y(3))-(x(1)-x(3))*(y(1)-y(2)))/2d0
                if(ar.le.0) then
                write(*,*) ' ar<=0 (flow) ar,im,id=',ar,im,idx(im)
                write(*,*)"xx(",ne(im,1),")=",x(1),"xx(",ne(im,2),")=",
     &                     x(2),"xx(",ne(im,3),")=",xx(3)
                write(*,*)"yy(",ne(im,2),")=",y(1),"yy(",ne(im,2),")=",
     &                 y(2),"yy(",ne(im,3),")=",y(3)
c               write(*,*) '   n=',ne(im,1),ne(im,2),ne(im,3)
c               write(*,*) '   x=',x(1),x(2),x(3)
c               write(*,*) '   y=',y(1),y(2),y(3)
c               call pltnt(nel,ne,xx,yy,idx,ind,ine)
cc                call plterr(nel,ne,xx,yy,x,y,ine,ind)
                pause
c               stop
                endif
         rav=(x(1)+x(2)+x(3))/3d0
c        tav=(tt(ne(im,1))+tt(ne(im,2))+tt(ne(im,3)))/3d0
c        uav=(uu(ne(im,1))+uu(ne(im,2))+uu(ne(im,3)))/3d0
c        vav=(vv(ne(im,1))+vv(ne(im,2))+vv(ne(im,3)))/3d0
c        pav=(pp(ne(im,1))+pp(ne(im,2))+pp(ne(im,3)))/3d0
         den=dens(id)         !ro
         vis=visc(id)			!1/Re
c     Ma数が上限を越えるのを防ぐ
         sv=den*svc(id)**2	!ro/(Ma**2)
c          sv1=den*svc(id)      !ro/Ma
c	    sv2=svc(id)          !Ma          sv=sv1*sv2
c	    write(*,*)"sv",sv
c
         ggy=den*gy(id)
c       write(*,*) den,vis,sv,ggy,ar,rav
         do 280 i=1,3
 280     sf(ne(im,i)*3-1)=sf(ne(im,i)*3-1)+ggy*ar/3d0*rav
c
        do 290 i=1,3
          iu=ne(im,i)*3-2
          iv=ne(im,i)*3-1
          ip=ne(im,i)*3
          do 290 j=1,3
            ju=ne(im,j)*3-2
            jv=ne(im,j)*3-1
            jp=ne(im,j)*3
            cmm=ar/12d0*rav*zkesu(i,j)
            sxx=b(i)*b(j)/4d0/ar*rav
            syy=c(i)*c(j)/4d0/ar*rav
            sxy=b(i)*c(j)/4d0/ar*rav
            syx=c(i)*b(j)/4d0/ar*rav
            chxij=b(j)/6d0*rav
            chxji=b(i)/6d0*rav
            chyij=c(j)/6d0*rav
            chyji=c(i)/6d0*rav
            a11=cmm*den
            a22=cmm*den
            a33=cmm
            b11=(2d0*sxx+syy)*vis+2d0*vis*cmm/rav**2
            b12=syx*vis
            b13=-chxji-cmm/rav
            b21=sxy*vis
            b22=(sxx+2d0*syy)*vis
            b23=-chyji
c      Ma数を1乗ずつ掛けて変数が上限を越えるのを防ぐ
c            write(*,*)"sv,chxij,sv,cmm,rav"
c	      write(*,*)sv,chxij,sv,cmm,rav
c	      write(*,*)sv*chxij,sv*cmm
c            b31=sv1*chxij*sv2+sv1*cmm/rav*sv2     !sv=sv1*sv2
c            b32=sv1*chyij*sv2
c	      write(*,*)"b31,b32",b31,b32
c	      pause
            b31=sv*chxij+sv*cmm/rav
            b32=sv*chyij
c
            b33=0d0
c---- 
            sa(iu,nbw+ju-iu)=sa(iu,nbw+ju-iu)+a11/dt+b11
            sa(iu,nbw+jv-iu)=sa(iu,nbw+jv-iu)       +b12
            sa(iu,nbw+jp-iu)=sa(iu,nbw+jp-iu)       +b13
            sf(iu)=sf(iu)+a11/dt*uvp(ju)
c---- 
            sa(iv,nbw+ju-iv)=sa(iv,nbw+ju-iv)       +b21
            sa(iv,nbw+jv-iv)=sa(iv,nbw+jv-iv)+a22/dt+b22
            sa(iv,nbw+jp-iv)=sa(iv,nbw+jp-iv)       +b23
            sf(iv)=sf(iv)+a22/dt*uvp(jv)
c----- 
            sa(ip,nbw+ju-ip)=sa(ip,nbw+ju-ip)       +b31
            sa(ip,nbw+jv-ip)=sa(ip,nbw+jv-ip)       +b32
            sa(ip,nbw+jp-ip)=sa(ip,nbw+jp-ip)+a33/dt+b33
            sf(ip)=sf(ip)+a33/dt*uvp(jp)
c----
 290    continue
 200    continue
        return
        end
c
        subroutine boundf1(nelf,ne,irw,jnw,ncrg,kcrg,
     &                                  nbd,kbd,bdv,uvp,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),uvp(ind),jnw(ind),irw(0:ine)
        dimension kcrg(21),kbd(1001),bdv(1001),jbcp(101)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /dbound/ kfbnd(-9:9,-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
c
c       write(*,*) 'boundf1/nelf,ncrg=',nelf,ncrg
c         nbd=0
          nbd=0        !
          do 100 k=1,ncrg
            i=kcrg(k)
            id=ndx(i)            !id:領域iのインデックス番号
              do 100 j=1,nb(i)   !i=2
              jd=ndx(ncom(i,j))   !jd:領域iの境界上の節点jが接している領域のインデックス番号
              ks=kfbnd(min(id,jd),max(id,jd))  !ks=0 at jd=0 or ks=1 at jd=-1(contacting with nr=1)
              if(ks.eq.1 .or. ks.eq.3) then
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j),neb(i,j)))*3-2
                bdv(nbd)=0d0
                if(nbd.ne.1 .and. kbd(nbd).eq.kbd(nbd-1)) nbd=nbd-1
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j+1),neb(i,j+1)))*3-2
                bdv(nbd)=0d0
              endif
              if(ks.eq.2 .or. ks.eq.3) then
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j),neb(i,j)))*3-1
                bdv(nbd)=0d0
                if(nbd.ne.1 .and. kbd(nbd).eq.kbd(nbd-1)) nbd=nbd-1
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j+1),neb(i,j+1)))*3-1
                bdv(nbd)=0d0
              endif
 100      continue
c
                if(nbd.gt.1001) then
                  write(*,*) 'boundf1/increase 1001'
                  stop
                endif
c
          do 120 i=1,ncp
            if(icl(i).lt.-900) goto 120
            do 140 j=1,nbd
            if(kbd(j).eq.jnw(kcp(i))*3-2) then
              if(icl(i).eq.0) then
                kbd(j)=iabs(kbd(j))
              else
                kbd(j)=-iabs(kbd(j))
                write(*,*)"kbd(",j,")=",kbd(j)       !
                pause
              endif
            endif
 140        continue
 120      continue
c
        call bounduvp(nbd,kbd,bdv,uvp,ind)
c
        return
        end
c
        subroutine bounduvp(nbd,kbd,bdv,uvp,ind)
        implicit double precision (a-h,o-z)
        dimension uvp(ind),kbd(1001),bdv(1001)
c
c       write(*,*) 'boundf1/nbd=',nbd,nelf
c
        do 200 i=1,nbd
 200    if(kbd(i).gt.0) uvp(kbd(i))=bdv(i)
c
        return
        end
c
        subroutine fst(ne,xx,yy,uu,vv,ine,ind)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),uu(ind),vv(ind)
        dimension uclm(101,6)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /dconhy/ stsa(21),zksa(21),zmsa(21),       !%yama(6.24)
     &          stsr(21),zksr(21),zmsr(21)
        save uclm
        data iav/5/
c
        do 100 n=1,ncp
          st(n)=-999
          if(icl(n).lt.-900) goto 100
          i=icp(n)
          j=jcp(n)
          id=ndx(i)
          jd=ndx(ncom(i,j))
          jdm=ndx(ncom(i,mcy(j-1,nb(i))))
          im=icoml(id,jd,jdm)
          stn=sten(min(id,jd),max(id,jd))
          if(stn.lt.1d-12) return
          kk=ne(nelb(i,j),neb(i,j))
          kp=ne(nelb(i,mcy(j+1,nb(i))),neb(i,mcy(j+1,nb(i))))
          km=ne(nelb(i,mcy(j-1,nb(i))),neb(i,mcy(j-1,nb(i))))
          ucl=uu(kk)
c----------------------
          ucl=0d0
          ii=0
          do 120 i=iav,2,-1
            uclm(n,i)=uclm(n,i-1)
            if(uclm(n,i).gt.-9000) then
            ii=ii+1
            ucl=ucl+uclm(n,i)
            endif
 120      continue
          uclm(n,1)=uu(kk)
          ucl=(ucl+uu(kk))/dble(ii+1)
c-----------------------
c         ucl=uu(kk)                    !!!!
c----------------------
        write(*,*) '            fst/ucl(av)=',ii,sngl(ucl)
c
          ancl(n)=atan((yy(kp)-yy(kk))/(xx(kk)-xx(kp)))
          if(ancl(n).lt.0d0) ancl(n)=ancl(n)+3.1415926d0
          ancl(n)=ancl(n)*180d0/3.1415926d0
c         write(*,*) ucl,st(n),ancl(n)
c
          if((icl(n).eq.1  .and. ucl.lt.0d0) .or.
     &                (icl(n).eq.-1 .and. ucl.gt.0d0)) then
            if(icl(n).eq.1) then
                staa(n)=dmax1(ancl(n),stsa(im))
                strr(n)=stsr(im)
            else if(icl(n).eq.-1) then
                staa(n)=stsa(im)
                strr(n)=dmin1(ancl(n),stsr(im))
            endif
            icla(n)=icl(n)
            icl(n)=0
c         endif
          elseif(icl(n).eq.0) then
            if(ancl(n).gt.staa(n)) then
              icla(n)=0
              icl(n)=1
              ucl=0d0
              st(n)=stsa(im)
              do 140 i=1,iav
 140          uclm(n,i)=-9999
            elseif(ancl(n).lt.strr(n)) then
              icla(n)=0
              icl(n)=-1
              ucl=0d0
              st(n)=stsr(im)
              do 160 i=1,iav
 160          uclm(n,i)=-9999
            else
              if(icla(n).eq.1.and.ancl(n).lt.staa(n)) staa(n)=stsa(im)
              if(icla(n).eq.-1.and.ancl(n).gt.strr(n)) strr(n)=stsr(im)
            endif
          endif
          if(icl(n).eq.0) then
           st(n)=-999
          else
           st(n)=999    !fconhy(im,icl(n),ucl)
          endif
        write(*,*) icla(n),icl(n),staa(n),strr(n)
 100    continue
        return
        end     
c
        subroutine boundfs(nelf,ne,xx,yy,tp0,uvp,jnw,cc0,
     &                                ncrg,kcrg,dsf,ind,ine)
c        subroutine boundfs(nelf,ne,xx,yy,uvp,jnw,ncrg,kcrg,dsf,
c     &                                      jbw,px2,py2,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),xx(ind),yy(ind),dsf(ind),jnw(ind),kcrg(21)
        dimension uvp(ind),px(ind),py(ind)
	  dimension flx(1001),tp0(ind),cc0(ind)    !蒸発の項の計算に用いる
c	  dimension jbw(ind),px2(ind),py2(ind)     !法線ベクトルの方向余弦の格納に用いる
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /doper/ vi,di,sta,str     !蒸発の項の計算に用いる
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
c
        do 20 n=1,ncrg
        i=kcrg(n)
        id=ndx(i)
        do 25 j=1,nb(i)
        nj=ne(nelb(i,j),neb(i,j))
        if(jnw(nj).gt.nelf) then
          write(*,*) jnw(nj),nelf
          stop
        endif
        px(jnw(nj))=0d0
 25     py(jnw(nj))=0d0
 20     continue
c
        do 100 n=1,ncrg
          i=kcrg(n)
          id=ndx(i)
          do 140 j=1,nb(i)
            nj=ne(nelb(i,j),neb(i,j))
            if((px(jnw(nj)).ne.0d0 .and. py(jnw(nj)).ne.0d0)) goto 140
            jd=ndx(ncom(i,j))
            jdm=ndx(ncom(i,mcy(j-1,nb(i))))
c           write(*,*) id,jdm,jd
            if(jd.ne.id .and. jdm.ne.id .and. jd.eq.jdm) then
              stn=sten(min(id,jd),max(id,jd))
              if(stn.lt.1d-12) goto 140
              n1=jnw(ne(nelb(i,j-1),neb(i,j-1)))
              n2=jnw(nj)
              n3=jnw(ne(nelb(i,j+1),neb(i,j+1)))
            elseif(jdm.eq.-1 .and. jd.ne.id) then
              stn=sten(min(id,jd),max(id,jd))
              if(stn.lt.1d-12) goto 140
              n1=-1
              n2=jnw(nj)
              n3=jnw(ne(nelb(i,j+1),neb(i,j+1)))
            elseif(jd.eq.-1 .and. jdm.ne.id) then
              stn=sten(min(id,jdm),max(id,jdm))
              if(stn.lt.1d-12) goto 140
              n1=jnw(ne(nelb(i,j-1),neb(i,j-1)))
              n2=jnw(nj)
              n3=-1
            elseif(jd.ne.id .and. jdm.ne.id. and. jdm.ne.jd .and.
     &             (jd.gt.0 .or. jdm.gt.0)) then
              m=ncom(i,j)
              if(m.eq.0) m=ncom(i,mcy(j-1,nb(i)))
              do 170 k=1,nb(m)
 170          if(ne(nelb(m,k),neb(m,k)).eq.nj) goto 180
                write(*,*) 'boundfs/error id,jdm,jd=',id,jdm,jd
                stop
 180          if(sten(min(id,jd),max(id,jd)).gt.1d-12) then
                n1=jnw(ne(nelb(m,k-1),neb(m,k-1)))
                n2=jnw(nj)
                n3=jnw(ne(nelb(i,j+1),neb(i,j+1)))
                stn=sten(min(id,jd),max(id,jd))
c               write(*,*) '#',n1,n2,n3
c    &             ,radii(xx,yy,n1,n2,n3,0d0,pxx,pyy,fxd,fyd,ind),stn
              elseif(sten(min(id,jdm),max(id,jdm)).gt.1d-12) then
                n1=jnw(ne(nelb(i,j-1),neb(i,j-1)))
                n2=jnw(nj)
                n3=jnw(ne(nelb(m,k+1),neb(m,k+1)))
                stn=sten(min(id,jdm),max(id,jdm))
c               write(*,*) '!',n1,n2,n3
c     &            ,radii(xx,yy,n1,n2,n3,0d0,pxx,pyy,fxd,fyd,ind),stn
              else
                goto 140
              endif
            else
              goto 140
            endif
            rd=radii(xx,yy,n1,n2,n3,0d0,pxx,pyy,fxd,fyd,ind)
            px(n2)=-2d0*stn*rd*pxx
            py(n2)=-2d0*stn*rd*pyy
c----------------2次要素の計算に用いる------------------------------------
c           px2(jbw(n2))=pxx        !法線ベクトルnの方向余弦nrを格納
c	      py2(jbw(n2))=pyy		  !法線ベクトルnの方向余弦nzを格納
c	      write(*,*)"pxx,pyy",pxx,pyy
c--------------------------------------------------------------------------
c       if(id.eq.1) write(*,*) j,jdm,jd,rd
c       if(jd.eq.-1 .and. jdm.eq.0 .and. id.eq.1) write(*,*) j,rd
 140      continue
 100    continue
c
        do 200 n=1,ncp
          if(st(n).lt.-900) goto 200
          i=icp(n)
          id=ndx(i)
          do 210 k=1,ncrg
 210      if(kcrg(n).eq.i) goto 220
          goto 200
 220      j=jcp(n)
          jd=ndx(ncom(i,j))
          jdm=ndx(ncom(i,mcy(j-1,nb(i))))
          stn=sten(min(id,jd),max(id,jd))
          nm1=jnw(ne(nelb(i,j-1),neb(i,j-1)))
          n1=jnw(ne(nelb(i,j),neb(i,j)))
          n2=jnw(ne(nelb(i,j+1),neb(i,j+1)))
          n3=jnw(ne(nelb(i,j+2),neb(i,j+2)))
          n4=jnw(ne(nelb(i,j+3),neb(i,j+3)))
c-----
          kk=ne(nelb(i,j),neb(i,j))
          ucl=uvp(jnw(kk)*3-2)
          im=icoml(id,jd,jdm)
          st(n)=fconhy(im,icl(n),ucl)
c       write(*,*) n,st(n)
c-----
          rd=radic(st(n),xx,yy,n1,n2,n3,n4,pxx,pyy,ind)
          px(n1)=-2d0*stn*rd*pxx
          py(n1)=-2d0*stn*rd*pyy
c----------------2次要素の計算に用いる------------------------------------
c          px2(jbw(n1))=pxx        !法線ベクトルnの方向余弦nrを格納
c	    py2(jbw(n1))=pyy		  !法線ベクトルnの方向余弦nzを格納
c          write(*,*)"pxx,pyy",pxx,pyy
c--------------------------------------------------------------------------
 200    continue
c
        do 300 n=1,ncrg
          i=kcrg(n)
          do 300 j=1,nb(i)
            j0=jnw(ne(nelb(i,j),neb(i,j)))
            j1=jnw(ne(nelb(i,j+1),neb(i,j+1)))
            rav=(xx(j0)+xx(j1))/2d0
            rl=dsqrt((xx(j1)-xx(j0))**2+(yy(j1)-yy(j0))**2)
            dsf(3*j0-2)=dsf(3*j0-2)+(2d0*px(j0)+px(j1))*rl*rav/6d0
            dsf(3*j1-2)=dsf(3*j1-2)+(px(j0)+2d0*px(j1))*rl*rav/6d0
            dsf(3*j0-1)=dsf(3*j0-1)+(2d0*py(j0)+py(j1))*rl*rav/6d0
 300      dsf(3*j1-1)=dsf(3*j1-1)+(py(j0)+2d0*py(j1))*rl*rav/6d0
c       write(*,*) 'boundfs/end'
c------------------------------------------------------------------
c----------------蒸発を考慮----------------------------------------
c------------------------------------------------------------------
          open(1,file="suf-flow.res")
	    rewind(1)
          call evpr2(ne,flx,tp0,cc0,ine,ind)
          lr=2    !流動領域
	    id=1    !流動領域
          visc0=viscos(-999d0,-999d0)  !代表の粘度μ0[Pa･s]   (therm.f参照)
          do 400 i=1,nb(lr)
	      if(ncom(lr,i).ne.0)goto 400     !自由表面のみ考慮
            j0=jnw(ne(nelb(lr,i),neb(lr,i)))
            j1=jnw(ne(nelb(lr,i+1),neb(lr,i+1)))
            rav=(xx(j0)+xx(j1))/2d0                          !中心軸からの平均距離
            rl=dsqrt((xx(j1)-xx(j0))**2+(yy(j1)-yy(j0))**2)  !辺の長さ
   	      sv=svc(id)**2*di/visc0                           !(c/(μ/ρr0))^2*r0/μ
            dsf(3*j0)=dsf(3*j0)-sv*flx(i)*rl/2d0*rav
            dsf(3*j1)=dsf(3*j1)-sv*flx(i+1)*rl/2d0*rav
 400      continue		         !400の行はcontinue にする
          close(1)
c------------------------------------------------------------------
c-------------------------------------------------------------------

        return
        end
c
        function radii(xx,yy,n1,n2,n3,cc,pxx,pyy,fxd,fyd,ind)
        implicit double precision (a-h,o-z)
        dimension s(3),x(3),y(3),xx(ind),yy(ind)
c
      if(n1.lt.0) then
        n1=n3
        x(1)=-xx(n3)
        y(1)= yy(n3)
      else
        x(1)= xx(n1)
        y(1)= yy(n1)
      end if
        x(2)= xx(n2)
        y(2)= yy(n2)
      if(n3.lt.0) then
        n3=n1
        x(3)=-xx(n1)
        y(3)= yy(n1)
      else
        x(3)= xx(n3)
        y(3)= yy(n3)
      end if
c
        s(1)=-dsqrt((x(2)-x(1))**2+(y(2)-y(1))**2)
        s(2)=0d0
        s(3)= dsqrt((x(3)-x(2))**2+(y(3)-y(2))**2)
c
        call lagrange(s,x,cc,fx,fxd,fxdd)
        call lagrange(s,y,cc,fy,fyd,fydd)
c
        if(dabs(fx).lt.1d-12) then
          radii=fxd*fydd/(fxd**2)**1.5
        else
          radii=(fx**2*(fxd*fydd-fxdd*fyd)+(fxd**2+fyd**2)*fx*fyd)
     &                            /(2d0*fx**2*(fxd**2+fyd**2)**1.5d0)
        endif
        pxx=fyd/dsqrt(fxd**2+fyd**2)
        pyy=-fxd/dsqrt(fxd**2+fyd**2)
        return
        end
c
        function radic(st,xx,yy,n1,n2,n3,n4,pxx,pyy,ind)
        implicit double precision (a-h,o-z)
        dimension s(3),x(3),y(3),xx(ind),yy(ind)
        data pi/3.1415926d0/
c
        r=radii(xx,yy,n1,n2,n3,0d0,px2,py2,fxd,fyd,ind)
        r=radii(xx,yy,n2,n3,n4,0d0,px3,py3,fxd,fyd,ind)
c
        stc=pi-st*pi/180d0
        sns=dsin(stc)
        css=dcos(stc)
        s(1)=-dsqrt((xx(n2)-xx(n1))**2+(yy(n2)-yy(n1))**2)
        s(2)=0d0
        s(3)= dsqrt((xx(n3)-xx(n2))**2+(yy(n3)-yy(n2))**2)
        x(1)=css
        x(2)=-py2
        x(3)=-py3
        y(1)=sns
        y(2)=px2
        y(3)=px3
        call lagrange(s,x,s(1),fx,dcs,fxdd)
        call lagrange(s,y,s(1),fy,dsn,fydd)
        if(xx(n1).ne.0d0) then
          radic=(dsn*css-dcs*sns+sns/xx(n1))/2d0
        else
          radic=0d0
        endif
        pxx=dsin(stc)
        pyy=-dcos(stc)
        return
        end
c
      subroutine lagrange(xl,yl,xi,ff,dydx,d2y)
      implicit double precision(a-h,o-z)
      dimension xl(3),yl(3)

      a1=yl(1)/((xl(1)-xl(2))*(xl(1)-xl(3)))
      a2=yl(2)/((xl(2)-xl(3))*(xl(2)-xl(1)))
      a3=yl(3)/((xl(3)-xl(1))*(xl(3)-xl(2)))
c
      ff  =a1*(xi-xl(2))*(xi-xl(3))
     &    +a2*(xi-xl(3))*(xi-xl(1))
     &    +a3*(xi-xl(1))*(xi-xl(2))
      dydx=a1*(2d0*xi-(xl(2)+xl(3)))
     &    +a2*(2d0*xi-(xl(3)+xl(1)))
     &    +a3*(2d0*xi-(xl(1)+xl(2)))
      d2y =2d0*(a1+a2+a3)
c
       return
       end
c
        subroutine shape(x,zn,znd,zndd)
        implicit double precision (a-h,o-z)
        dimension zn(4),znd(4),zndd(4)
        zn(1)=-x/2d0*(1d0-x)
          zn(2)=(1d0+x)*(1d0-x)
        zn(3)=x/2d0*(1d0+x)
         znd(1)=x-.5d0
        znd(2)=-2d0*x
        znd(3)=x+.5d0
        zndd(1)=1d0
        zndd(2)=-2d0
        zndd(3)=1d0
        return
        end

c       
        subroutine boundf2(nelf,ne,irw,jnw,ncrg,kcrg,
     &                                  nbd,kbd,bdv,uvp,ind,ine)
        implicit double precision (a-h,o-z)
        dimension ne(ine,3),uvp(ind),jnw(ind),irw(0:ine)
        dimension kcrg(21),kbd(1001),bdv(1001),jbcp(101)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /dbound/ kfbnd(-9:9,-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
c
c       write(*,*) 'boundf1/nelf,ncrg=',nelf,ncrg
c         nbd=0

c
c---------------------追加-----------------------------------
          do i=1,nb(2)
		  if(ncom(2,i).eq.0)goto 1000
		  enddo
1000      ii=i
c------------------------------------------------------------


          nbd=0        !
          do 100 k=1,ncrg
            i=kcrg(k)
            id=ndx(i)            !id:領域iのインデックス番号
              do 100 j=1,nb(i)   !i=2
              jd=ndx(ncom(i,j))   !jd:領域iの境界上の節点jが接している領域のインデックス番号
              ks=kfbnd(min(id,jd),max(id,jd))
c              if(ks.eq.1 .or. ks.eq.3) then     !u=0    
              if((ks.eq.1 .or. ks.eq.3)) then
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j),neb(i,j)))*3-2
                bdv(nbd)=0d0
                if(nbd.ne.1)then 
                  if(kbd(nbd).eq.kbd(nbd-1)) nbd=nbd-1
                endif
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j+1),neb(i,j+1)))*3-2
                bdv(nbd)=0d0
              endif
c              if(ks.eq.2 .or. ks.eq.3) then        !v=0
              if((ks.eq.2 .or. ks.eq.3) .and. j.lt.ii) then      !条件追加
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j),neb(i,j)))*3-1
                bdv(nbd)=0d0
                if(nbd.ne.1)then
                  if(kbd(nbd).eq.kbd(nbd-1)) nbd=nbd-1
                endif
                nbd=nbd+1
                kbd(nbd)=jnw(ne(nelb(i,j+1),neb(i,j+1)))*3-1
                bdv(nbd)=0d0
              endif
 100      continue
c
                if(nbd.gt.1001) then
                  write(*,*) 'boundf1/increase 1001'
                  stop
                endif
c
          do 120 i=1,ncp
            if(icl(i).lt.-900) goto 120
            do 140 j=1,nbd
            if(kbd(j).eq.jnw(kcp(i))*3-2) then
              if(icl(i).eq.0) then
                kbd(j)=iabs(kbd(j))
              else
                kbd(j)=-iabs(kbd(j))
                write(*,*)"kbd(",j,")=",kbd(j)       !
                pause
              endif
            endif
 140        continue
 120      continue
c
        call bounduvp(nbd,kbd,bdv,uvp,ind)
c
        return
        end
