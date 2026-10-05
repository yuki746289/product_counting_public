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
c        do i=1,nw
c          sf(i)=sf(i)+dsf(i)
c        enddo
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
       subroutine boundf1(new,la,mb,melb,meb,mcom,nbd,mbd,vbd,mbd2,vbd2)
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
       dimension mbd2(ibd),vbd2(ibd)
c
       nbd=0
       do ka=1,la
         do j=1,mb(ka)
           jd=mcom(ka,j)    !接している要素のindex番号
c
c          U=1, V=0
           if(jd.eq.-1 .or. jd.eq.-3)then
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j),meb(ka,j))*2-1
             vbd(nbd)=1d0
c
             mbd2(nbd)=new(melb(ka,j),meb(ka,j))*2-0
             vbd2(nbd)=0d0
             if(nbd.ne.1)then
               if(mbd(nbd).eq.mbd(nbd-1))nbd=nbd-1
             endif
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j+1),meb(ka,j+1))*2-1
             vbd(nbd)=1d0
c
             mbd2(nbd)=new(melb(ka,j+1),meb(ka,j+1))*2-0
             vbd2(nbd)=0d0
           endif
c
c          V=0, U=0
           if(jd.eq.-2 .or. jd.eq.-3)then
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j),meb(ka,j))*2-0
             vbd(nbd)=0d0
c
             mbd2(nbd)=new(melb(ka,j),meb(ka,j))*2-0
             vbd2(nbd)=0d0
             if(nbd.ne.1)then
               if(mbd(nbd).eq.mbd(nbd-1))nbd=nbd-1
             endif
             nbd=nbd+1
             mbd(nbd)=new(melb(ka,j+1),meb(ka,j+1))*2-0
             vbd(nbd)=0d0
c
             mbd2(nbd)=new(melb(ka,j+1),meb(ka,j+1))*2-0
             vbd2(nbd)=0d0
           endif
         enddo
       enddo
c
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
       subroutine boundfp(new,la,mb,melb,meb,mcom,nbp,mbp,vbp)
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
       dimension mbp(ibd),vbp(ibd)
c
       nbp=0
       do ka=1,la
         do j=1,mb(ka)
           jd=mcom(ka,j)    !接している要素のindex番号
c
c          P=0
           if(jd.eq.-1 .or. jd.eq.-2 .or. jd.eq.-3)then
             nbp=nbp+1
             mbp(nbp)=new(melb(ka,j),meb(ka,j))     !saの対角成分
             vbp(nbp)=0d0
             if(nbp.ne.1)then
               if(mbp(nbp).eq.mbp(nbp-1))nbp=nbp-1
             endif
             nbp=nbp+1
             mbp(nbp)=new(melb(ka,j+1),meb(ka,j+1))
             vbp(nbp)=0d0
           endif
         enddo
       enddo
c
       return
       end
