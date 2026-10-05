c-----------------------------------------------------------
c
c      要素の再発生
c
c-----------------------------------------------------------
       subroutine grid(np,nele,ne,idx,xx,yy,uu,vv,pp,nr
     &                 ,nc,nelc,ndx,mr,na,nb,neb,nelb)
       include "head.for"
       dimension ne(ine,3),idx(0:ine)
       dimension xx(ind),yy(ind)
       dimension uu(ind),vv(ind),pp(ind)
c
       dimension nc(ica),nelc(ica,-1:ibd)
       dimension ndx(0:ica)
       dimension neled(0:ica)
       dimension mr(ica)
       dimension nb(ica)
       dimension nelb(ica,-1:ibd)
       dimension neb(ica,-1:ibd)
c
       dimension ne0(ine,3)
       dimension net(ine,3)
       dimension new(ine,3),idw(0:ine)
       dimension x0(ind),y0(ind),z0(ind)
       dimension u0(ind),v0(ind),p0(ind)
       dimension xw(ind),yw(ind),zw(ind)
       dimension uw(ind),vw(ind),pw(ind)
       dimension m0(ind)
       dimension nf(ica),xf(ica,ibd),yf(ica,ibd),mf(ica,ibd)
       dimension ng(ica),xg(ica,ibd),yg(ica,ibd),mg(ica,ibd)
       dimension ip(ind),jp(ind),mw(ind),mn(ind)
       dimension melt(ine)
c
c      領域ごとに発生させた節点
c      nf(nr),xf(nr,nf(nr)),yf(nr,nf(nr)),mf(nr,nf(nr))
c      ng(nr),xg(nr,nf(nr)),yg(nr,nf(nr))
       do 100 kr=1,nr
         la=0
         np0=0
         do 200 ka=1,na
           if(ndx(kr).ne.mr(kr))goto 200
           do kb=1,nb(ka)
             np0=np0+1
             x0(np0)=xx(ne(nelb(ka,kb),neb(ka,kb)))
             y0(np0)=yy(ne(nelb(ka,kb),neb(ka,kb)))
             z0(np0)=0d0
             ip(np0)=ne(nelb(ka,kb),neb(ka,kb))     !ip(新)=古
             jp(ne(nelb(ka,kb),neb(ka,kb)))=np0     !jp(古)=新
           enddo
200      continue
         n0=np0
c         write(*,*)"n0:",n0
c
c        内部の節点を発生させる
         call setnp(np0,x0,y0,z0)
c         write(*,*)"np0:",np0
c
         nf(kr)=0
         do kp=1,n0
           nf(kr)=nf(kr)+1
           xf(kr,nf(kr))=x0(kp)
           yf(kr,nf(kr))=y0(kp)
           mf(kr,nf(kr))=ip(kp)  !境界の節点番号
         enddo
c         write(*,*)"nf:",kr,nf(kr)
c
         ng(kr)=0
         do kp=n0+1,np0
           ng(kr)=ng(kr)+1
           xg(kr,ng(kr))=x0(kp)
           yg(kr,ng(kr))=y0(kp)
         enddo
c         write(*,*)"ng:",kr,ng(kr)
100    continue    !領域ループ
c
c      領域ごとに要素を発生させる
c      nw,xw(nw),yw(nw),uw(nw),vw(nw),pw(nw)
c      nelw,new(nelw,3)
c      mw(nw)
       nw=0
       nelw=0
       do 300 kr=1,nr
         np0=0
c
c        境界の節点
         do kf=1,nf(kr)
           np0=np0+1
           x0(np0)=xf(kr,kf)
           y0(np0)=yf(kr,kf)
           m0(np0)=mf(kr,kf)
         enddo
c
c        内部の節点
         do kg=1,ng(kr)
           np0=np0+1
           x0(np0)=xg(kr,kg)
           y0(np0)=yg(kr,kg)
           m0(np0)=0
         enddo
c         write(*,*)"np0:",kr,np0
c         pause
c
c        要素を発生させる
         call mesh(np0,x0,y0,nel0,ne0)
c         do kele=1,nel0
c           do j=1,3
c             write(*,*)"ne0:",kele,j,ne0(kele,j)
c           enddo
c           pause
c         enddo
c
c        境界の裏側についている要素を削除する
         nelt=0
         do 400 ka=1,na
           if(ndx(kr).ne.mr(kr))goto 400
           do kb=1,nb(ka)
             k1=jp(ne(nelb(ka,kb+0),neb(ka,kb+0)))
             k2=jp(ne(nelb(ka,kb+1),neb(ka,kb+1)))
             do kel0=1,nel0
               l1=ne0(kel0,1)
               l2=ne0(kel0,2)
               l3=ne0(kel0,3)
               if(k2.eq.l1 .and. k1.eq.l2 .or.
     &            k2.eq.l2 .and. k1.eq.l3 .or.
     &            k2.eq.l3 .and. k1.eq.l1)then
                 nelt=nelt+1
                 melt(nelt)=kel0
               endif
             enddo
           enddo
400      continue
c
         do 900 kel0=1,nel0
           l1=ne0(kel0,1)
           l2=ne0(kel0,2)
           l3=ne0(kel0,3)
c
c          3つの節点が境界の要素
           ik1=0
           ik2=0
           ik3=0
           do 700 ka=1,na
             if(ndx(kr).ne.mr(kr))goto 700
             do kb=1,nb(ka)
               kp=jp(ne(nelb(ka,kb+0),neb(ka,kb+0)))
               if(kp.eq.l1)ik1=1
               if(kp.eq.l2)ik2=1
               if(kp.eq.l3)ik3=1
             enddo
700        continue
           if(ik1.eq.0 .or. ik2.eq.0 .or. ik3.eq.0)goto 900
c
c          境界の節点の順番と合っていない要素
           do 800 ka=1,na
             if(ndx(kr).ne.mr(kr))goto 800
             do kb=1,nb(ka)
               k1=jp(ne(nelb(ka,kb+0),neb(ka,kb+0)))
               k2=jp(ne(nelb(ka,kb+1),neb(ka,kb+1)))
               if(k1.eq.l1 .and. k2.eq.l2)goto 900
               if(k1.eq.l2 .and. k2.eq.l3)goto 900
               if(k1.eq.l3 .and. k2.eq.l1)goto 900
             enddo
             nelt=nelt+1
             melt(nelt)=kel0
800        continue
900      continue
c
         write(*,*)"nel0:",nel0
         write(*,*)"nelt:",nelt
         do 600 kele=1,nel0
           do kelt=1,nelt
             if(melt(kelt).eq.kele)goto 600
           enddo
           nelw=nelw+1
           new(nelw,1)=ne0(kele,1)
           new(nelw,2)=ne0(kele,2)
           new(nelw,3)=ne0(kele,3)
           idw(nelw)=ndx(kr)
600      continue
c
c         nel0=nelw
c         do kel0=1,nel0
c           ne0(kel0,1)=new(kel0,1)
c           ne0(kel0,2)=new(kel0,2)
c           ne0(kel0,3)=new(kel0,3)
c         enddo
c
c        速度・圧力の補間
         nelt=0
         do kc=1,nc(kr)
           nelt=nelt+1
           net(nelt,1)=ne(nelc(kr,kc),1)
           net(nelt,2)=ne(nelc(kr,kc),2)
           net(nelt,3)=ne(nelc(kr,kc),3)
         enddo
         call refvl(nelt,net,xx,yy,uu,vv,pp,x0,y0,u0,v0,p0)
c
         do kp=1,np0
           nw=nw+1
           xw(nw)=x0(kp)
           yw(nw)=y0(kp)
           uw(nw)=u0(kp)
           vw(nw)=v0(kp)
           pw(nw)=p0(kp)
           mw(nw)=m0(kp)    !0:内部の節点,0以外:境界の節点
           ip(nw)=kp
           jp(kp)=nw
         enddo
300    continue
       call pltne(nw,nelw,new,xw,yw)
c       write(*,*)"nw:",nw
c       write(*,*)"nelw:",nelw
c       pause
c
c      境界の節点番号が重複しているので割り振り直す
c      np,ip(np),jp(kw)
       np=0
       do 500 kw=1,nw
c        内部の節点
         if(mw(kw).eq.0)then
           np=np+1
           ip(np)=kw    !ip(新)=古(複数)
           jp(kw)=np    !jp(古)=新
           mn(np)=mw(kw)
c        境界の節点
         else
           do kp=1,np
             if(mw(kw).eq.mn(kp))then
               ip(kp)=kw
               jp(kw)=kp
               goto 500
             endif
           enddo
           np=np+1
           ip(np)=kw
           jp(kw)=np
           mn(np)=mw(kw)
         endif
500    continue
c       write(*,*)"np:",np
c       pause
c
       do kw=1,nw
         xx(jp(kw))=xw(kw)
         yy(jp(kw))=yw(kw)
         uu(jp(kw))=uw(kw)
         vv(jp(kw))=vw(kw)
         pp(jp(kw))=pw(kw)
       enddo
c
       nele=nelw
       do kele=1,nele
         ne(kele,1)=jp(new(kele,1))
         ne(kele,2)=jp(new(kele,2))
         ne(kele,3)=jp(new(kele,3))
         idx(kele)=idw(kele)
       enddo
c
       return
       end
c-----------------------------------------------------------
c
c       要素分割するか調べる
c
c-----------------------------------------------------------
        function icheck(nele,ne,xx,yy)
        include 'head.for'
        dimension ne(ine,3)
        dimension xx(ind),yy(ind)
c        
        do kele=1,nele
          x1=xx(ne(kele,1))
          y1=yy(ne(kele,1))
          x2=xx(ne(kele,2))
          y2=yy(ne(kele,2))
          x3=xx(ne(kele,3))
          y3=yy(ne(kele,3))
c
          ar=artri(x1,y1,0d0,x2,y2,0d0,x3,y3,0d0)
          d1=dsqrt((x2-x3)**2+(y2-y3)**2)
          d2=dsqrt((x3-x1)**2+(y3-y1)**2)
          d3=dsqrt((x1-x2)**2+(y1-y2)**2)
c
          h1=ar/d1
          h2=ar/d2
          h3=ar/d3
c
          if(h1.gt.h2*2d0 .or. h1.gt.h3*2d0)then
            icheck=1
            write(*,*)"h",sngl(h1),sngl(2d0*h2/h1),sngl(2d0*h3/h1)
            write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
            write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
            write(*,*)"x3:",sngl(x3),sngl(y3),sngl(z3)
            pause
            return
          endif
c
          if(h2.gt.h3*2d0 .or. h2.gt.h1*2d0)then
            icheck=1
            write(*,*)"h",sngl(h2),sngl(2d0*h3/h2),sngl(2d0*h1/h2)
            write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
            write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
            write(*,*)"x3:",sngl(x3),sngl(y3),sngl(z3)
            pause
            return
          endif
c
          if(h3.gt.h1*2d0 .or. h3.gt.h2*2d0)then
            icheck=1
            write(*,*)"h",sngl(h3),sngl(2d0*h1/h3),sngl(2d0*h2/h3)
            write(*,*)"x1:",sngl(x1),sngl(y1),sngl(z1)
            write(*,*)"x2:",sngl(x2),sngl(y2),sngl(z2)
            write(*,*)"x3:",sngl(x3),sngl(y3),sngl(z3)
            pause
            return
          endif

        enddo
c
        icheck=0    !0:要素そのまま,1:要素再分割
        return
        end
c-----------------------------------------------------------
c
c       速度を補間する
c
c-----------------------------------------------------------
        subroutine refvl(nele,ne,xx,yy,uu,vv,pp,x0,y0,u0,v0,p0)
        include 'head.for'
        dimension ne(ine,3)
        dimension xx(ind),yy(ind)
        dimension x0(ind),y0(ind)
        dimension uu(ind),vv(ind),pp(ind)
        dimension u0(ind),v0(ind),p0(ind)
c
        do 100 kp0=1,np0
          xi=x0(kp0)
          yi=y0(kp0)
          do 200 kele=1,nele
            kp1=ne(kele,1)
            kp2=ne(kele,2)
            kp3=ne(kele,3)
c
            x1=xx(kp1)
            y1=yy(kp1)
            x2=xx(kp2)
            y2=yy(kp2)
            x3=xx(kp3)
            y3=yy(kp3)
c
            ar =artri(x1,y1,0d0,x2,y2,0d0,x3,y3,0d0)
            cn1=artri(xi,yi,0d0,x2,y2,0d0,x3,y3,0d0)/ar
            cn2=artri(x1,y1,0d0,xi,yi,0d0,x3,y3,0d0)/ar
            cn3=artri(x1,y1,0d0,x2,y2,0d0,xi,yi,0d0)/ar
            if((cn1+cn2+cn3).gt.1d0+1d-10)goto 200
c
            u0(kp0)=cn1*uu(kp1)+cn2*uu(kp2)+cn3*uu(kp3)
            v0(kp0)=cn1*vv(kp1)+cn2*vv(kp2)+cn3*vv(kp3)
            p0(kp0)=cn1*pp(kp1)+cn2*pp(kp2)+cn3*pp(kp3)
            goto 100
200      continue
100    continue

       return
       end
