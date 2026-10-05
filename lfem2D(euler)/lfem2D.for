       include "head.for"
       dimension ne(ine,3),net(ine,3),new(ine,3)
       dimension idx(0:ine),idt(0:ine),idw(0:ine)
       dimension xx(ind),yy(ind)
       dimension xt(ind),yt(ind)
       dimension xw(ind),yw(ind)
       dimension xn(ind),yn(ind)
c
c      物性値
       dimension vre(-5:ica),vma(-5:ica),vfr(-5:ica)
       dimension vwe(-5:ica,-5:ica)    !※※※※※※※
c
c       境界
        dimension nc(ica),nelc(ica,-1:ibd)
        dimension ndx(0:ica)
        dimension neled(0:ica)
        dimension mr(ica)
        dimension nb(ica)
        dimension nelb(ica,-1:ibd)
        dimension neb(ica,-1:ibd)
        dimension ncom(ica,-1:ibd)
c
        dimension mb(ica)
        dimension melb(ica,-1:ibd)
        dimension meb(ica,-1:ibd)
        dimension mcom(ica,-1:ibd)
c
c       バンド幅
        dimension ip(ind),jp(ind)
        dimension iq(ine),jq(ine)
        dimension jww(ind),jvv(ind)
c
c       物理量
        dimension uu(ind),vv(ind),pp(ind)
        dimension ut(ind),vt(ind),pt(ind)
        dimension uw(ind),vw(ind),pw(ind)
c
c       マトリックス
        dimension uvp(ind*3),uvp0(ind*3)
        dimension mbd(ibd),vbd(ibd)    !速度
        dimension mbp(ibd),vbp(ibd)    !圧力
        dimension mbd2(ibd),vbd2(ibd)    !速度
        dimension sa(ind,ibw),sf(ind),dsf(ind)
c
        character*40 tfile
        character*12 lpp
c
        data zkasu/1.0/    !緩和係数
        data emax/1d-7/  !収束判定:時間間隔よりも小さくする
c
c     ファイルを読み込む
       call data(iyy,dt,wdt1,wdt2,tmax,vi,ri,vre,vma,vfr,vwe
     &           ,np,nele,ne,idx,xx,yy)
c
c      初期化
       if(iyy.eq.0)then
         ti=0d0
         tsav=wdt1     !保存時間
         ttmp=wdt2     !保存時間
         do kp=1,np
           uu(kp)=0d0
           vv(kp)=0d0
           pp(kp)=0d0
         enddo
       else
         open(1,file="tmp.inp")
         rewind(1)
         read(1,*)ti,tsav,ttmp
         read(1,*)np,nele
         read(1,*)(ne(i,1),ne(i,2),ne(i,3),i=1,nele)
         read(1,*)(idx(i),i=0,nele)
         read(1,*)(xx(i),yy(i),i=1,np)
         read(1,*)(uu(i),vv(i),pp(i),i=1,np)
         close(1)
       endif
c
c      界面を探す
       call sufele(nele,ne,idx,nr,nc,nelc
     &                    ,ndx,neled,na,mr,nb,nelb,neb,ncom)
c
c        call pltne(np,nele,ne,xx,yy)
c        call pltnb(na,nb,nelb,neb,ncom)
c-----------------------------------------------------------------------------
c
c      計算ループ
c
c-----------------------------------------------------------------------------
       it=0
100    continue
       it=it+1
c
       write(*,*)"------------------------------"
       write(*,*)"ti:",sngl(ti)
c
c      領域ループ
       do 1000 kr=1,nr
         if(ndx(kr).le.0)goto 1000
c
c        領域krの要素
c        nt,nelt
c        net(nelt,3)
c        xt(nt),yt(nt)
c        ut(nt),vt(nt),pt(nt)
         nelt=0
         do kc=1,nc(kr)
           nelt=nelt+1
           net(nelt,1)=ne(nelc(kr,kc),1)
           net(nelt,2)=ne(nelc(kr,kc),2)
           net(nelt,3)=ne(nelc(kr,kc),3)
           idt(nelt)=idx(nelc(kr,kc))
           iq(nelt)=nelc(kr,kc)     !iq(新)=古
           jq(nelc(kr,kc))=nelt     !jq(古)=新
c           write(*,*)"kc:",kc,nc(kr),nelc(kr,kc)
c           write(*,*)"kc:",kc,nc(kr)
c     &                      ,net(nelt,1),net(nelt,2),net(nelt,3)
c           pause
         enddo
c
         nt=0
         do 200 i=1,np
           ich=0
           do kelt=1,nelt
             if(net(kelt,1).eq.i)ich=1
             if(net(kelt,2).eq.i)ich=1
             if(net(kelt,3).eq.i)ich=1
           enddo
           if(ich.eq.0)goto 200
           nt=nt+1
           xt(nt)=xx(i)
           yt(nt)=yy(i)
           ut(nt)=uu(i)
           vt(nt)=vv(i)
           pt(nt)=pp(i)
           ip(nt)=i     !ip(新)=古
           jp(i)=nt     !ip(古)=新
200      continue
c
         do kelt=1,nelt
           net(kelt,1)=jp(net(kelt,1))
           net(kelt,2)=jp(net(kelt,2))
           net(kelt,3)=jp(net(kelt,3))
c           write(*,*)"kelt:",kelt,nelt
c     &                      ,net(kelt,1),net(kelt,2),net(kelt,3)
c           pause
         enddo
c
c       バンド幅
c       nw,nelw
c       new(nelw,3)
c       xw(nw),yw(nw)
c       uw(nw),vw(nw),pw(nw)
c        call renum2(index,np,ne,nel,jww,ine,ind,nband)
c        do kt=1,nt
c          jvv(jww(kt))=kt    !jvv(新_新)=新
c        enddo
c
        do kt=1,nt
          jww(kt)=nt-kt+1    !jww(新)=新_新
          jvv(jww(kt))=kt    !jvv(新_新)=新
        enddo
c
        nbm=-99999
        do kele=1,nelt
        do i=1,3
        do j=1,3
          kbm=iabs(net(kele,i)-net(kele,j))
          if(kbm.gt.nbm)nbm=kbm     !nbm:現在考慮しているtreeのバンド幅
        enddo
        enddo
        enddo
        nbw=nbm+1
        write(*,*)'band',nbw
c
        nw=nt
        do kw=1,nw  !新_新
          xw(kw)=xt(jvv(kw))
          yw(kw)=yt(jvv(kw))
          uw(kw)=ut(jvv(kw))
          vw(kw)=vt(jvv(kw))
          pw(kw)=pt(jvv(kw))
        enddo
c
        nelw=nelt
        do kelw=1,nelw
          new(kelw,1)=jww(net(kelw,1))
          new(kelw,2)=jww(net(kelw,2))
          new(kelw,3)=jww(net(kelw,3))
          idw(kelw)=idt(kelw)
        enddo
c
c       境界の節点
c       la,mb(la),melb(la,mb(la)),meb(la,mb(la))
        la=0
        do 300 ka=1,na
          if(kr.ne.mr(ka))goto 300 !領域krの境界を考慮
          la=la+1
          mb(la)=nb(ka)
          do kb=-1,nb(ka)+2
            melb(la,kb)=jq(nelb(ka,kb))
            meb(la,kb)=neb(ka,kb)
            mcom(la,kb)=ncom(ka,kb)
          enddo
300     continue
c
c        call pltne(nw,nelw,new,xw,yw)
c        call pltnb(la,mb,melb,meb,mcom)
c
c       マトリックス
c       uvp0(nw*3)
c       xn(nw),yn(nw)
        do kw=1,nw
          uvp0(kw*3-2)=uw(kw)
          uvp0(kw*3-1)=vw(kw)
          uvp0(kw*3-0)=pw(kw)
          uvp(kw*3-2)=uw(kw)
          uvp(kw*3-1)=vw(kw)
          uvp(kw*3-0)=pw(kw)
          sf(kw)=pw(kw)
c
          xn(kw)=xw(kw)
          yn(kw)=yw(kw)
        enddo
c
c       第1種の境界条件
c       nbd,mbd(nbd),vbd(nbd)
        call boundf1(new,la,mb,melb,meb,mcom,nbd,mbd,vbd,mbd2,vbd2)
c        call boundfp(new,la,mb,melb,meb,mcom,nbp,mbp,vbp)
c
c       流動解析
c       nelw,new(nelw,3),idw(nelw)
c       nw,xw(nw),yw(nw)
c       uw(nw),vw(nw),pw(nw)
c       la,mb(la),melb(la,mb(la)),meb(la,mb(la))
        do 400 lk=1,26
c
c         N-S式のマトリクスを解く  sa(3*ind,ibw),sf(3*ind)
          do i=1,nw     !pn0(npn)の代入
            uvp(3*i-0)=uvp(3*i-0)+zkasu*(sf(i)-uvp(3*i-0)) !zkasu:緩和係数
c            write(*,*)"p uvp:",i,nw,sngl(uvp(3*i-2)),sngl(uvp(3*i-1))
c            pause
          enddo
c
          call setmat(nw,nbw,sa,sf,dsf)                        !マトリクスを初期化する
          call flowns0(nelw,nbw*2,new,idw,xn,yn,uvp0,uvp
     &                                 ,dt,vre,vma,vfr,sa,sf)
          call boundasym(nw*2,nbw*2,nbd,mbd,vbd,dsf,sa,sf)
          call boundasym(nw*2,nbw*2,nbd,mbd2,vbd2,dsf,sa,sf)
          call gauss(2*nw,2*nbw,sa,sf)             !ガウスの消去法でマトリクスを解く
c          call boundf1nsr(ne,sf,ip,jww,nn,mn,ind,ine)         !列ベクトルsfにu,v,pの境界条件を代入
c
c         速度u,v,pの残差を求める
          err1=-999d0
          umx=-999d0                                           !列ベクトルに用いる速度:uvp0(3*ind)
          vmx=-999d0                                           !対流項に用いる速度:uvp(3*ind)
          do i=1,nw
            umx=dmax1(umx,dabs(sf(2*i-1)))
            vmx=dmax1(vmx,dabs(sf(2*i-0)))
            if(err1.lt.dabs(uvp(3*i-2)-sf(2*i-1)) 
     &                +dabs(uvp(3*i-0)-sf(2*i-0)))then
              err1=dabs(uvp(3*i-2)-sf(2*i-1))
     &            +dabs(uvp(3*i-1)-sf(2*i-0))
c              ier=i
            endif
          enddo         

c         do im=1,nelw
c           do i=1,3     !要素imの行
c             ik=new(im,i)
c             do j=1,3   !要素imの列
c               jk=new(im,j)
c               write(*,*)'sa(1*',new(im,i),'-0,1*',new(im,j),'-0)'
c     &                                          ,sngl(sa(ik,nbw+jk-ik))
c             enddo
c           enddo
c         enddo
c         pause
c
c         do im=1,nelw
c           do i=1,3     !要素imの行
c             ik=new(im,i)
c             write(*,*)'sf(1*',new(im,i),'-0)',sngl(sf(ik))
c           enddo
c         enddo
c         pause

c
c         連続の式のマトリクスを解く  sa(ind,ibw),sf(ind)
          do i=1,nw
            uvp(3*i-2)=uvp(3*i-2)+zkasu*(sf(2*i-1)-uvp(3*i-2))    !速度のx成分
            uvp(3*i-1)=uvp(3*i-1)+zkasu*(sf(2*i-0)-uvp(3*i-1))    !速度のy成分
c            if(it.ge.2)then
c              write(*,*)"u uvp:",i,nw,sngl(uvp(3*i-2)),sngl(uvp(3*i-1))
c              write(*,*)"xx:",sngl(xn(i)),sngl(yn(i))
c              pause
c            endif
          enddo
          call setmat(nw,nbw,sa,sf,dsf)                        !マトリクスを初期化する
          call flowcn(nelw,nbw*1,new,idw,xn,yn,uvp0,uvp
     &                                ,dt,vre,vma,vfr,sa,sf)
c          call boundasym(nw*1,nbw*1,nbp,mbp,vbp,dsf,sa,sf)
c          call boundf1cn(npn,ne,idn,sa,sf,nbw,ip,jww,nn,mn)   !マトリクスsaにpの境界条件を代入
          call gauss(1*nw,1*nbw,sa,sf)             !ガウスの消去法でマトリクスを解く
c          call boundf1cnr(ne,sf,ip,jww,nn,mn,ind,ine)         !列ベクトルsfにpの境界条件を代入
c
c         圧力pの残差を求める
          err2=-999d0
          pmx=-999d0
          do 410 i=1,nw
            pmx=dmax1(pmx,dabs(sf(i)))
            if(err2.lt.dabs(uvp(3*i-0)-sf(i)))then
              err2=dabs(uvp(3*i-0)-sf(i))
c              ier=i
            endif
410      continue
c1900      err=dmax1(err,dabs(uvp0(3*i-0)-sf(i)))
          err=err1+err2
c
c         表示
          write(*,*)"--------------------------------------------------"
          write(*,*)'lk,err1 err2 ',lk,err1,err2
          write(*,*)'umx,vmx,pmx',sngl(umx),sngl(vmx),
     &                                 sngl(pmx)
c
          if(err.lt.emax)goto 500   !終了判定
400    continue    !収束ループ
c
       open(10,file="err.res",access='append')
       write(10,*)"ti",sngl(ti),"  lk",lk,"  err",sngl(err),
     &               "   No convergion"
       close(10)
       write(*,*) ' !!!No convergion in lk loop/err=',sngl(err)
500    continue
c
c      値コピー
       do kw=1,nw   !新_新
         uu(ip(jvv(kw)))=uvp(3*kw-2)
         vv(ip(jvv(kw)))=uvp(3*kw-1)
         pp(ip(jvv(kw)))=uvp(3*kw-0)
       enddo
c
       call pltne(np,nele,ne,xx,yy)
       call pltnb(na,nb,nelb,neb,ncom)
c-------------------------------------------------------------------------------
c
c      要素の再発生
c
c-------------------------------------------------------------------------------
c       if(icheck(nele,ne,xx,yy).eq.1)then
c         write(*,*)"mesh is regenerated."
cc        要素の再発生
c         call grid(np,nele,ne,idx,xx,yy,uu,vv,pp,nr
c     &                 ,nc,nelc,ndx,mr,na,nb,neb,nelb)
cc        界面を探す
c         call sufele(nele,ne,idx,nr,nc,nelc
c     &                      ,ndx,neled,na,mr,nb,nelb,neb,ncom)
cc
c         open(10,file="err.res",access='append')
c         write(10,*)"ti",sngl(ti),"  remesh"
c         close(10)
cc         pause
c       endif
c-------------------------------------------------------------------------------
c
c      データ保存
c
c-------------------------------------------------------------------------------
c
c      avs fileの保存
       if(ti.gt.(tsav-dt*1d-2) .and. ti.lt.(tsav+dt*1d-2))then
         write(*,*)'**** save avs file data ****'
         call svdata(np,nele,ne,xx,yy,uu,vv,pp,idx,ti)
         tsav=tsav+wdt1
       endif
c
c      tmp fileの保存
       if(ti.gt.(ttmp-dt*1d-2) .and. ti.lt.(ttmp+dt*1d-2))then
         write(*,*)'**** save tmp file data ****'
c         write(lpp,'(f12.1)') ti*1d4  !input a value into lpp
         write(lpp,'(f12.1)') ti*1d6  !input a value into lpp
         do i=1,12
           if(lpp(i:i).eq.' ') lpp(i:i)='0'
	   enddo
         call addchr("./tmp/tmp",lpp//'.inp',tfile)
         open(1,file=tfile)
         rewind(1)
         write(1,*)ti,tsav,ttmp
         write(1,*)np,nele
         write(1,*)(ne(i,1),ne(i,2),ne(i,3),i=1,nele)
         write(1,*)(idx(i),i=0,nele)
         write(1,*)(xx(i),yy(i),i=1,np)
         write(1,*)(uu(i),vv(i),pp(i),i=1,np)
         close(1)
         ttmp=ttmp+wdt2
       endif
c-------------------------------------------------------------------------------
c
c      タイムステップ
c
c-------------------------------------------------------------------------------
1000   continue !タイムステップ
c
       ti=ti+dt
       if(ti.lt.tmax)goto 100
c
       stop
       end
