c-----------------------------------------------------------------------------------
c
c      第1種、第2種の境界条件をマトリックスに代入する
c
c-----------------------------------------------------------------------------------
        subroutine boundasym(nw,nbw,nbd,mbd,vbd,dsf,sa,sf)
        include "head.for"
        dimension mbd(ibd),vbd(ibd)
        dimension dsf(ind)
        dimension sa(ind,ibw),sf(ind)
c
c       第2種境界条件(応力)
        do i=1,nw
          sf(i)=sf(i)+dsf(i)
        enddo
c
c       第1種境界条件
        nbwm=nbw-1           !境界上の節点速度0の条件
        do 30 i=1,nbd
          ii=mbd(i)
          if(ii.lt.0) goto 30
          aa=sa(ii,nbw)
          if(aa.eq.0) then
            write(*,*) ii,aa
            stop
          endif
          do 40 j=max(1,ii-nbwm),min(ii+nbwm,nw)
            sf(j)=sf(j)-sa(j,nbw+ii-j)*vbd(i)   !右辺に移項
            sa(ii,nbw+j-ii)=0d0
            sa(j,nbw+ii-j)=0d0
 40       continue
          sa(ii,nbw)=aa
          sf(ii)=aa*vbd(i)
 30     continue
        return
        end
c-----------------------------------------------------------------------------------
c
c      第1種境界条件の配列を作成する
c
c      引数
c      new(nelw,3)
c      la
c      mb(la)
c      melb(la,mb(la))
c      meb(la,mb(la))
c      mcom(la,mb(la))
c
c      戻り値
c      nbd      :境界の節点数
c      mbd(nbd) :境界の節点番号
c      vbd(nbd) :境界の速度
c
c-----------------------------------------------------------------------------------
       subroutine boundf1(new,la,mb,melb,meb,mcom,nbd,mbd,vbd)
       include "head.for"
       dimension new(ine,3)
       dimension idw(0:ine)
       dimension ndx(0:ica)
       dimension mb(ica)
       dimension melb(ica,-1:ibd)
       dimension meb(ica,-1:ibd)
       dimension mcom(ica,-1:ibd)
       dimension uvp(ind*3),uvp0(ind*3)
c
       dimension mbd(ibd),vbd(ibd)
c
       nbd=0
       do ka=1,la
         do j=1,mb(ka)
           jd=mcom(ka,j)    !接している要素のindex番号
c
c          U=0
           if(jd.eq.-1 .or. jd.eq.-3)then
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j),meb(ka,j))*3-2
             vbd(nbd)=0d0
             if(nbd.ne.1)then
               if(mbd(nbd).eq.mbd(nbd-1))nbd=nbd-1
             endif
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j+1),meb(ka,j+1))*3-2
             vbd(nbd)=0d0
           endif
c
c          V=0
           if(jd.eq.-2 .or. jd.eq.-3)then
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j),meb(ka,j))*3-1
             vbd(nbd)=0d0
             if(nbd.ne.1)then
               if(mbd(nbd).eq.mbd(nbd-1))nbd=nbd-1
             endif
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j+1),meb(ka,j+1))*3-1
             vbd(nbd)=0d0
           endif
         enddo
       enddo
c
       return
       end
c-----------------------------------------------------------------------------------
c
c      第2種境界条件の配列を作成する(自由表面)
c
c      引数
c
c      戻り値
c      dsf(nw*3)
c
c
c-----------------------------------------------------------------------------------
        subroutine boundfs(kr,ndx,new,xn,yn,vwe
     &                      ,la,mb,melb,meb,mcom,dsf)
        include "head.for"
        dimension new(ine,3)
        dimension xn(ind),yn(ind)
        dimension ndx(0:ica)
        dimension mb(ica)
        dimension melb(ica,-1:ibd)
        dimension meb(ica,-1:ibd)
        dimension mcom(ica,-1:ibd)
        dimension vwe(-5:ica,-5:ica)
        dimension dsf(ind)
c
        dimension px(ind),py(ind)
c
        id=ndx(kr)
        if(id.le.0)return
c
        do 100 ka=1,la
c
          do 140 j=1,mb(ka)
            m1=new(melb(ka,j-2),meb(ka,j-2))
            m2=new(melb(ka,j-1),meb(ka,j-1))
            m3=new(melb(ka,j+0),meb(ka,j+0))
            m4=new(melb(ka,j+1),meb(ka,j+1))
            m5=new(melb(ka,j+2),meb(ka,j+2))
c
            id1=mcom(ka,j-2)
            id2=mcom(ka,j-1)
            id3=mcom(ka,j-0)
            id4=mcom(ka,j+1)
c            id5=mcom(melb(ka,j+2),meb(ka,j+2))
c
c           中央
            if(id2.eq.id3 .and. id2.ge.0)then
              n1=m2
              n2=m3
              n3=m4
              sten=1d0/vwe(id,id2)
c           左側
            elseif(id1.eq.id2 .and. id1.ge.0)then
              n1=m1
              n2=m2
              n3=m3
              sten=1d0/vwe(id,id1)
c           右側
            elseif(id3.eq.id4 .and. id3.ge.0)then
              n1=m3
              n2=m4
              n3=m5
              sten=1d0/vwe(id,id3)
            else
              px(m3)=0d0
              py(m3)=0d0
              goto 140
            endif
c
c           曲率計算
            call radii(xn,yn,n1,n2,n3,rd,pxx,pyy)
            px(m3)=-2d0*sten/rd*pxx
            py(m3)=-2d0*sten/rd*pyy
c            write(*,*)"bd:",m3,sngl(rd),sngl(pxx),sngl(pyy)
c            pause
 140      continue
 100    continue
c
c       応力代入
        do 200 ka=1,la
        do 200 j=1,mb(ka)
          j0=new(melb(ka,j+0),meb(ka,j+0))
          j1=new(melb(ka,j+1),meb(ka,j+1))
c
          rl=dsqrt((xn(j1)-xn(j0))**2+(yn(j1)-yn(j0))**2)
          dsf(3*j0-2)=dsf(3*j0-2)+(2d0*px(j0)+    px(j1))*rl/6d0
          dsf(3*j1-2)=dsf(3*j1-2)+(    px(j0)+2d0*px(j1))*rl/6d0
          dsf(3*j0-1)=dsf(3*j0-1)+(2d0*py(j0)+    py(j1))*rl/6d0
          dsf(3*j1-1)=dsf(3*j1-1)+(    py(j0)+2d0*py(j1))*rl/6d0
c          write(*,*)"p:",j0,j1
c     &                   ,sngl((2d0*px(j0)+    px(j1))*rl/6d0)
c     &                   ,sngl((    px(j0)+2d0*px(j1))*rl/6d0)
c     &                   ,sngl((2d0*py(j0)+    py(j1))*rl/6d0)
c     &                   ,sngl((    py(j0)+2d0*py(j1))*rl/6d0)
c          write(*,*)"j0:",j0,sngl(px(j0)),sngl(py(j0))
c          write(*,*)"j1:",j1,sngl(px(j1)),sngl(py(j1))
c          pause
200     continue
c
        return
        end
