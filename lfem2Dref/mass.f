      subroutine bunpu(np,tt,xx,yy,fo,ind)
	implicit double precision (a-h,o-z)
	dimension xx(ind),yy(ind),tt(ind),con(ind)
c
      write(*,*)"fo=",fo,"np=",np
      pi=3.1415926535897d0
      do 10 i=1,np
	  if(xx(i).lt.1d-1)goto 10
	  if(yy(i).gt.5d0 .and. yy(i).le.7d0) then
	    r=dsqrt((xx(i)-1d-1)**2+(yy(i)-6d0)**2)
	    if(r.lt.1d-5) r=1d-3
	  else if(yy(i).lt.5d0 .and. yy(i).ge.3d0) then
	    r=dsqrt((xx(i)-1d-1)**2+(yy(i)-4d0)**2)
	    if(r.lt.1d-5) r=1d-3
	  elseif(yy(i).eq.5d0)then
	    r=1d0
	  else
	    goto 10
	  endif
	  cono=0d0
	  do 20 n=1,20000
	    con(i)=cono+2d0/pi/r*((-1d0)**n)/dble(n)*
     &           dsin(dble(n)*pi*r)*
     &  	       dexp(-(dble(n**2))*(pi**2)*fo)
	    cono=con(i)
 20	  continue
	  tt(i)=cono+1d0
	  if(r.eq.1d-3) tt(i)=0d0
 10	continue
c
	return
	end
c
      subroutine btw(new,nuv,kuv,uvb,tt,ine,ind)
      implicit double precision (a-h,o-z)
      dimension new(ine,6),kuv(2001,10),uvb(2001,10),tt(ind),nuv(10)
	common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c--
      lr=2
      nuv(lr)=0
      do 10 i=1,nbw(lr)+1
	if(ncomw(lr,i).ne.0)goto 10
      k=new(nelbw(lr,i),nebw(lr,i))
      nuv(lr)=nuv(lr)+1
      kuv(nuv(lr),lr)=k
      uvb(nuv(lr),lr)=1.0d0
      tt(k)=uvb(nuv(lr),lr)
 10   continue
      return
      end
c
      subroutine calw(nel,nen,xxw,yyw,t0,nelrs,nelre,cov,ind,ine)
	implicit double precision (a-h,o-z)
	dimension new(ine,6),xxw(ind),yyw(ind)
	dimension t0(ind)
	dimension pn(6,7)
      dimension el(4,2),ww(7)
c	dimension p(6,6,7),a(6,6,7),b(6,6,7),c(6,6,7),d(6,6,7)
	dimension nelrs(21),nelre(21)
	dimension jp21(6),nen(ine,6)
	common /dprop/dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
	common /mass/dmas(9)
	data jp21/1,3,5,2,4,6/
cn	data (el(i,1),i=1,3)//
cn	data (el(i,2),i=1,3)/0.5d0,0.5d0,0d0/
c
      do 200 i=1,nel
	do 200 j=1,6
200	new(i,j)=nen(i,jp21(j))
      lr=2     !considering area
      id=1
	pi=3.1415926535897d0
      pe=(1d0/visc(id))*dmas(id)  !Pe=Re*Sc
	write(*,*)"pe,visc(id),dmas(id)",pe,visc(id),dmas(id)
	zlam=1d0/pe
	ww(1)=-27d0/96d0
	ww(2)=25d0/96d0
	ww(3)=25d0/96d0
	ww(4)=25d0/96d0
	el(1,1)=1d0/3d0
	el(2,1)=11d0/15d0
	el(3,1)=2d0/15d0
	el(4,1)=2d0/15d0
	el(1,2)=1d0/3d0
	el(2,2)=2d0/15d0
	el(3,2)=2d0/15d0
	el(4,2)=11d0/15d0
cn	ww(1)=1d0/6d0
cn	ww(2)=1d0/6d0
cn	ww(3)=1d0/6d0
	do 10 i=1,4
	pn(1,i)=el(i,1)*(2d0*el(i,1)-1d0)
	pn(2,i)=4d0*el(i,1)*el(i,2)
	pn(3,i)=el(i,2)*(2d0*el(i,2)-1d0)
	pn(4,i)=4d0*el(i,2)*(1d0-el(i,1)-el(i,2))
	pn(5,i)=(1d0-el(i,1)-el(i,2))*(1d0-2d0*el(i,1)-2d0*el(i,2))
 10	pn(6,i)=4d0*el(i,1)*(1d0-el(i,1)-el(i,2))
c
cn	do 10 i=1,4
cn	pn(1,i)=cc1*(2d0*cc1-1d0)
cn	pn(2,i)=4d0*cc1*cc2
cn	pn(3,i)=cc2*(2d0*cc2-1d0)
cn	pn(4,i)=4d0*cc2*(1d0-cc1-cc2)
cn	pn(5,i)=(1d0-cc1-cc2)*(1d0-2d0*cc1-2d0*cc2)
cn 10	pn(6,i)=4d0*cc1*(1d0-cc1-cc2)
c
cn	aa=0.445948490915965d0
cn	bb=0.091576213509771d0
cn	ww(2)=0.111690794839005d0
cn	ww(4)=0.111690794839005d0
cn	ww(6)=0.111690794839005d0
cn	ww(1)=0.054975871827661d0
cn	ww(3)=0.054975871827661d0
cn	ww(5)=0.054975871827661d0
cnn	www=1d0/2d0
cnn	cc1=1d0/3d0
cnn	cc2=1d0/3d0
cnn	pn(1,1)=cc1*(2d0*cc1-1d0)
cnn	pn(2,1)=4d0*cc1*cc2
cnn	pn(3,1)=cc2*(2d0*cc2-1d0)
cnn	pn(4,1)=4d0*cc2*(1d0-cc1-cc2)
cnn	pn(5,1)=(1d0-cc1-cc2)*(1d0-2d0*cc1-2d0*cc2)
cnn	pn(6,1)=4d0*cc1*(1d0-cc1-cc2)
cn	do 5 i=1,6
cn	if(i.eq.2.or.i.eq.6) cc1=aa
cn	if(i.eq.3.or.i.eq.5) cc1=bb
cn	if(i.eq.4) cc1=1d0-2d0*aa
cn	if(i.eq.1) cc1=1d0-2d0*bb
cn	if(i.eq.2.or.i.eq.4) cc2=aa
cn	if(i.eq.1.or.i.eq.5) cc2=bb
cn	if(i.eq.6) cc2=1d0-2d0*aa
cn	if(i.eq.3) cc2=1d0-2d0*bb
cn	pn(1,i)=cc1*(2d0*cc1-1d0)
cn	pn(2,i)=4d0*cc1*cc2
cn	pn(3,i)=cc2*(2d0*cc2-1d0)
cn	pn(4,i)=4d0*cc2*(1d0-cc1-cc2)
cn	pn(5,i)=(1d0-cc1-cc2)*(1d0-2d0*cc1-2d0*cc2)
cn	pn(6,i)=4d0*cc1*(1d0-cc1-cc2)
cn 5	continue
c
c--
	cov=0d0
c	do 100 im=1,nel
      do 100 im=nelrs(lr),nelre(lr)
      call jcbd(new,im,xxw,yyw,detj,ind,ine)
	rav=0d0
	covp=0d0
	do 110 j=1,6
 110	rav=rav+xxw(new(im,j))/6.d0
	do 150 i=1,4  !6
	do 170 j=1,6
 170	covp=covp+2d0*pi*rav*ww(i)*pn(j,i)*detj*t0(new(im,j))
 150	continue
	cov=cov+covp
cn	do 130 j=1,6
cn 130	covp=covp+2d0*pi*rav*www*pn(j,1)*detj*t0(new(im,j))
cn	cov=cov+covp
 100	continue
	return
	end

      subroutine sufar(ne,xx,yy,area,ind,ine)   !1 order element
      implicit double precision (a-h,o-z)
      dimension ne(ine,3),xx(ind),yy(ind)
      common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &       neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
c---
      lr=2
 	pi=3.1415926535897d0
	area=0d0
	sss1=0d0
	do 20 i=1,nb(lr)-1
	if(ncom(lr,i).ne.0)goto 20
	sarea=0d0
	sss=dsqrt((xx(ne(nelb(lr,i+1),neb(lr,i+1)))
     &          -xx(ne(nelb(lr,i),neb(lr,i))))**2
     &	     +(yy(ne(nelb(lr,i+1),neb(lr,i+1)))
     &          -yy(ne(nelb(lr,i),neb(lr,i))))**2)
	ssfr=sss/40d0
	do 21 j=0,40-1
	xi=(xx(ne(nelb(lr,i+1),neb(lr,i+1)))
     &   -xx(ne(nelb(lr,i),neb(lr,i))))/40d0
     &	   *dble(j)+xx(ne(nelb(lr,i),neb(lr,i)))
	xi1=(xx(ne(nelb(lr,i+1),neb(lr,i+1)))
     &   -xx(ne(nelb(lr,i),neb(lr,i))))/40d0
     &	   *dble(j+1)+xx(ne(nelb(lr,i),neb(lr,i)))
c	sfarea=2d0*pi*(xi+xi1-1d-1)/2d0*ssfr   !S=2π*R*L  (断面の半分を用いているので2を掛ける)
	sfarea=2d0*pi*(xi+xi1)/2d0*ssfr         !領域1の厚みを考慮する必要はない
	sarea=sarea+sfarea
 21	continue
	area=area+sarea
c
 20	continue
	return
	end
c
      subroutine bound(np,nbw,sa,sf,nb,kb,tb,ind,ibw)
      implicit double precision (a-h,o-z)
      dimension nb(10),kb(2001,10),tb(2001,10),sa(ind,ibw),sf(ind)
      nbwm=nbw-1
cn      do 10 jjj=1,2
cn      jj=jjj
      lr=2
      do 30 i=1,nb(lr)
      ii=kb(i,lr)
      aa=sa(ii,nbw)
      do 40 j=max(1,ii-nbwm),min(ii+nbwm,np)
      sf(j)=sf(j)-sa(j,nbw+ii-j)*tb(i,lr)
      sa(ii,nbw+j-ii)=0d0
 40   sa(j,nbw+ii-j)=0d0
      sa(ii,nbw)=aa
      sf(ii)=aa*tb(i,lr)
 30   continue
cn      if (ico.eq.1) goto 50
cn 10   continue 
 50   continue 
      return
      end
c
      subroutine noudot(nel,nbww,new,xxw,yyw,tp,
     &                     ne2,jww,jp,saw,sfw,dt,t0,cc0,ind,ine,ibw)
	implicit double precision (a-h,o-z)
	dimension new(ine,6),xxw(ind),yyw(ind)
	dimension saw(ind,ibw),sfw(ind),t0(ind)
	dimension li(36),lj(36),pn(6,7)
      dimension ww(7),vj(2,2),vjr(2,2),den(2,6,7)        !,el(3,2),ej(2,6)
	dimension p(6,6,7),a(2,6,7),b(6,2,7),c(6,6,7)
	dimension flx(1001)        !蒸発速度m・[kg/m2・s]
      dimension tp(ind)          !温度
      dimension cc0(ind)         !前回計算した濃度
	dimension jww(ind),jp(ind),ne2(ine,6) !蒸発の計算に用いる
      common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
	common /mass/ dmas(9)
      common /doper/ vi,di,sta,str   !蒸発に用いる
      common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                   ncomw(21,1001),ncopw(21,1001)


c      	data (ej(1,j),j=1,6)/3d0,2d0,0d0,-2d0,-3d0,0d0/
c      	data (ej(2,j),j=1,6)/0d0,2d0,3d0,0d0,-3d0,-2d0/
c	data (el(i,1),i=1,3)/0.5d0,0d0,0.5d0/
c	data (el(i,2),i=1,3)/0.5d0,0.5d0,0d0/
	data li/1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,
     &			4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6/
	data lj/1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6,
     &			1,2,3,4,5,6,1,2,3,4,5,6,1,2,3,4,5,6/
c
c	lr=2
	id=1
	write(*,*)"dt=",dt
      write(*,*)"visc(id),dmas(id)",visc(id),dmas(id)
	pe=(1d0/visc(id))*dmas(id)  !Pe=Re*Sc
	write(*,*)"pe=",pe
	zlam=1d0/pe
	alpha=0.05961587d0
	beta=0.47014206d0
	gamma=0.10128651d0
	delta=0.79742699d0
	ww(1)=0.11250d0
	ww(2)=0.66197075d0
	ww(3)=0.66197075d0
	ww(4)=0.66197075d0
	ww(5)=0.06296959d0
	ww(6)=0.06296959d0
	ww(7)=0.06296959d0
	do 5 i=1,7		      !積分点7
	if(i.eq.1) cc1=1d0/3d0
	if(i.eq.2) cc1=alpha
	if(i.eq.3.or.i.eq.4) cc1=beta
	if(i.eq.5.or.i.eq.7) cc1=gamma
	if(i.eq.6) cc1=delta
	if(i.eq.1) cc2=1d0/3d0
	if(i.eq.2.or.i.eq.3) cc2=beta
	if(i.eq.4) cc2=alpha
	if(i.eq.5.or.i.eq.6) cc2=gamma
	if(i.eq.7) cc2=delta
c
	den(1,1,i)=4d0*cc1-1d0
	den(1,4,i)=4d0*cc2
	den(1,2,i)=0d0
	den(1,5,i)=-4d0*cc2
	den(1,3,i)=4d0*cc1+4d0*cc2-3d0
	den(1,6,i)=4d0*(-2d0*cc1-cc2+1d0)
	den(2,1,i)=0d0
	den(2,4,i)=4d0*cc1
	den(2,2,i)=4d0*cc2-1d0
	den(2,5,i)=4d0*(-cc1-2d0*cc2+1d0)
	den(2,3,i)=4d0*cc1+4d0*cc2-3d0
	den(2,6,i)=-4d0*cc1
	pn(1,i)=cc1*(2d0*cc1-1d0)
	pn(4,i)=4d0*cc1*cc2
	pn(2,i)=cc2*(2d0*cc2-1d0)
	pn(5,i)=4d0*cc2*(1d0-cc1-cc2)
	pn(3,i)=(1d0-cc1-cc2)*(1d0-2d0*cc1-2d0*cc2)
	pn(6,i)=4d0*cc1*(1d0-cc1-cc2)
 5	continue
c
	do 10 im=1,nel
c      do 10 im=nelrs(lr),nelre(lr)
	vj(1,1)=0d0
	vj(1,2)=0d0
	vj(2,1)=0d0
	vj(2,2)=0d0
	vjr(1,1)=0d0
	vjr(1,2)=0d0
	vjr(2,1)=0d0
	vjr(2,2)=0d0
c	do 40 j=1,6
c	vj(1,1)=vj(1,1)+ej(1,j)*xxw(new(im,j))
c	vj(1,2)=vj(1,2)+ej(1,j)*yyw(new(im,j))
c	vj(2,1)=vj(2,1)+ej(2,j)*xxw(new(im,j))
c 40 	vj(2,2)=vj(2,2)+ej(2,j)*yyw(new(im,j))
cc	vj(1,1)=-xxw(new(im,1))-3d0*xxw(new(im,5))+4d0*xxw(new(im,6))   !要素内の節点番号変更
cc	vj(1,2)=-yyw(new(im,1))-3d0*yyw(new(im,5))+4d0*yyw(new(im,6))
cc	vj(2,1)=-xxw(new(im,3))-3d0*xxw(new(im,5))+4d0*xxw(new(im,4))
cc	vj(2,2)=-yyw(new(im,3))-3d0*yyw(new(im,5))+4d0*yyw(new(im,4))
c
	vj(1,1)=-xxw(new(im,1))-3d0*xxw(new(im,3))+4d0*xxw(new(im,6))
	vj(1,2)=-yyw(new(im,1))-3d0*yyw(new(im,3))+4d0*yyw(new(im,6))
	vj(2,1)=-xxw(new(im,2))-3d0*xxw(new(im,3))+4d0*xxw(new(im,5))
	vj(2,2)=-yyw(new(im,2))-3d0*yyw(new(im,3))+4d0*yyw(new(im,5))
c
c	vj(1,1)=xxw(new(im,1))-xxw(new(im,3))
c	vj(1,2)=yyw(new(im,1))-yyw(new(im,3))
c	vj(2,1)=xxw(new(im,2))-xxw(new(im,3))
c	vj(2,2)=yyw(new(im,2))-yyw(new(im,3))
c
	detj=vj(1,1)*vj(2,2)-vj(1,2)*vj(2,1)
c	detj=dabs(detj)
c	do 45 j=1,6
c	vjr(1,1)=vjr(1,1)+ej(2,j)*yyw(new(im,j))/detj
c	vjr(1,2)=vjr(1,2)-ej(1,j)*yyw(new(im,j))/detj
c	vjr(2,1)=vjr(2,1)-ej(2,j)*xxw(new(im,j))/detj
c 45	vjr(2,2)=vjr(2,2)+ej(1,j)*xxw(new(im,j))/detj
cc	vjr(1,1)=(-yyw(new(im,3))-3d0*yyw(new(im,5))   !要素内の節点番号変更
cc    &	+4d0*yyw(new(im,4)))
cc   &	/detj
cc	vjr(1,2)=-(-yyw(new(im,1))-3d0*yyw(new(im,5))
cc     &	+4d0*yyw(new(im,6)))
cc     &	/detj
cc	vjr(2,1)=-(-xxw(new(im,3))-3d0*xxw(new(im,5))
cc     &	+4d0*xxw(new(im,4)))
cc     &	/detj
cc	vjr(2,2)=(-xxw(new(im,1))-3d0*xxw(new(im,5))
cc     &	+4d0*xxw(new(im,6)))
cc     &	/detj
	vjr(1,1)=vj(2,2)/detj   !逆行列
	vjr(1,2)=-vj(1,2)/detj
	vjr(2,1)=-vj(2,1)/detj
	vjr(2,2)=vj(1,1)/detj
	detj=dabs(detj)
c
	do 50 n=1,7		  !積分点
	do 60 i=1,2		  !L1,L2
	do 70 j=1,6		  !節点
	a(i,j,n)=vjr(i,1)*den(1,j,n)+vjr(i,2)*den(2,j,n)
	b(j,i,n)=vjr(i,1)*den(1,j,n)+vjr(i,2)*den(2,j,n)
 70	continue
 60	continue
	do 80 i=1,6
	do 90 j=1,6
	c(i,j,n)=b(i,1,n)*a(1,j,n)+b(i,2,n)*a(2,j,n)
	p(i,j,n)=pn(i,n)*pn(j,n)
 90	continue
 80	continue
 50	continue

	open(40,file='pmat.dat')
	rewind(40)
	do 105 ki=1,7
	do 105 kj=1,6
	do 105 kk=1,6
 105	write(40,*) kj,kk,ki,sngl(p(kj,kk,ki))
	close(40)
c	open(41,file='cmat.dat')
c	rewind(41)
c	do 60 ki=1,3
c	do 60 kj=1,6
c	do 60 kk=1,6
c 60	write(41,*) kj,kk,ki,sngl(c(kj,kk,ki))
c	close(41)
c	pause
c--
	rav=0d0
	do 110 j=1,6
 110	rav=rav+xxw(new(im,j))/6.d0
c	write(*,*) 'detj=',sngl(detj)
c      rav=rav-1d-1    !領域1の厚みを除く
	do 150 i=1,7	!積分点
	do 200 l=1,36	!節点
	ii=li(l)
	jj=lj(l)
	ik=new(im,ii)
	jk=new(im,jj)
	saw(ik,nbww+jk-ik)=saw(ik,nbww+jk-ik)+
     &				p(ii,jj,i)*detj*rav*ww(i)/dt+
     &				c(ii,jj,i)*detj*rav*ww(i)*zlam
	sfw(ik)=sfw(ik)+p(ii,jj,i)*detj*rav*ww(i)/dt*t0(jk)
 200	continue
 150	continue
 10	continue
c     蒸発を考慮
	write(*,*)"蒸発の境界条件"
	lr=2
      visc0=viscos(-999d0,-999d0)  !代表の粘度μ0[Pa･s]
      dens0=density(-999d0,-999d0)       !代表の密度ρ0[kg/m^3]
      call evpr(ne2,tp,cc0,flx,ine,ind)   !蒸発速度flx(nbw(lr))を求める(3角形2次要素)
c
	open(1,file="cont.res")
	rewind(1)
c
      do 300 i=1,nbw(2),2               !jp(古)=新,jww(新)=新_新
	  if(ncomw(2,i).ne.0)goto 300     !自由表面のみ考慮
	  j0=jww(jp(ne2(nelbw(2,i),nebw(2,i))))         !頂点の節点j0
	  j1=jww(jp(ne2(nelbw(2,i+1),nebw(2,i+1))))     !辺上の節点j1
	  j2=jww(jp(ne2(nelbw(2,i+2),nebw(2,i+2))))     !頂点の節点j2
        denst=(density(tp(ne2(nelbw(2,i),nebw(2,i))), !無次元密度ρ-=ρ/ρ0[-]
     &               cc0(ne2(nelbw(2,i),nebw(2,i))))
     &       +density(tp(ne2(nelbw(2,i+1),nebw(2,i+1))),
     &               cc0(ne2(nelbw(2,i+1),nebw(2,i+1)))))/(dens0*2d0)
	  rav=(xxw(j0)+xxw(j1)+xxw(j2))/3d0             !辺j0-j2の平均半径R
        rl=dsqrt((xxw(j0)-xxw(j2))**2+(yyw(j0)-yyw(j2))**2)   !辺j0-j2の長さl
c
        hi=-1d0/denst*di/visc0*flx(i)   !-1/ρ-*r0/μ*m･
c        hi=-flx(i)/vi/ro              !m/(v0*ro)
c
c	write(*,*)"j0,j1",j0,j1
c	write(*,*)sngl(xxw(j0)),sngl(yyw(j0))
c	write(*,*)sngl(xxw(j1)),sngl(yyw(j1))
c	write(*,*)"rav,rl,tav,t0,t1",sngl(rav),sngl(rl),
c     &                    sngl(tav),sngl(t0(j0)),sngl(t0(j1))
c	write(*,*)"flx(",i,"),hi",flx(i),hi
c
c       液相内に流入する(濃度が増加する)ので右辺第2項は正となる
	  sfw(j0)=sfw(j0)-hi*(t0(j0)+1d0)*rl*1d0/6d0*rav  !頂点の節点j0
	  sfw(j1)=sfw(j1)-hi*(t0(j1)+1d0)*rl*4d0/6d0*rav  !辺上の節点j1
	  sfw(j2)=sfw(j2)-hi*(t0(j2)+1d0)*rl*1d0/6d0*rav  !頂点の節点j2
	  write(1,*)j0,sngl(xxw(j0)),sngl(yyw(j0)),sngl(t0(j0)),
     &               sngl(-hi*(t0(j0)+1d0)*rl*1d0/6d0*rav)
        write(1,*)j1,sngl(xxw(j1)),sngl(yyw(j1)),sngl(t0(j1)),
     &               sngl(-hi*(t0(j1)+1d0)*rl*4d0/6d0*rav)
        write(1,*)j2,sngl(xxw(j2)),sngl(yyw(j2)),sngl(t0(j2)),
     &               sngl(-hi*(t0(j2)+1d0)*rl*1d0/6d0*rav)
c	  
300   continue       !300の行はcontinueにする
      close(1)
c
	return
	end
c	
      subroutine riron(ti,pe,fo,flth,amth)
	implicit double precision (a-h,o-z)
c	dimension con(2000),con1(2000),con2(2000)
c
	pi=3.1415926535897d0
	rr=1.0d0
	foo=ti/pe+fo
	con=0d0
	con1=0d0
	con2=0d0
	cono=0d0
	cono1=0d0
	cono2=0d0
	do 150 n=1,1000
	con=cono+6d0/(pi**2)/dble(n**2)*
     &			dexp(-foo*(dble(n**2))*(pi**2))
	cono=con
	con1=cono1-2d0/pi/(rr*2)*((-1d0)**n)/dble(n)*
     &		dsin(dble(n)*pi*rr)*
     &		dexp(-foo*((dble(n))**2)*(pi**2))
	con2=cono2+2d0/rr*((-1d0)**n)*
     &		dcos(dble(n)*pi*rr)*
     &		dexp(-foo*((dble(n))**2)*(pi**2))
	cono1=con1
	cono2=con2
 150	continue
	flth=con1+con2
	amth=-con+1d0
	return
	end
c
      subroutine grdfc5(nel,np,nbw1,ne,xx,yy,tt,
     &      neln,npn,nbwn,nen,xxn,yyn,ttn,irx,ind,ine)
      implicit double precision (a-h,o-z)
      dimension ne(ine,6),nen(ine,6)
      dimension xx(ind),yy(ind),tt(ind)
      dimension xxn(ind),yyn(ind),ttn(ind)
	dimension irx(0:ine),tb(ind)
c      dimension arr(6001),al(3),al0(3),jp1(3),jp2(3),jp3(3),jp4(3)
c      data jp1/2,3,1/,jp2/3,1,2/,jp3/1,3,5/,jp4/3,5,1/
      dimension arr(6001),al(3),al0(3),jp3(3)
	common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
      data jp3/1,2,3/
         if(ine.gt.8001) then
           write(*,*) 'change dim in grdfc / ine=',ine
           call clos
           stop
         endif
      do 80 i=1,nel
c      arr(i)=(xx(ne(i,jp3(2)))-xx(ne(i,jp3(1))))*
c     &	     (yy(ne(i,jp3(3)))-yy(ne(i,jp3(1))))
c     &      -(xx(ne(i,jp3(3)))-xx(ne(i,jp3(1))))*
c     &	     (yy(ne(i,jp3(2)))-yy(ne(i,jp3(1))))
      arr(i)=(xx(ne(i,jp3(1)))-xx(ne(i,jp3(3))))*
     &       (yy(ne(i,jp3(2)))-yy(ne(i,jp3(3))))
     &      -(xx(ne(i,jp3(2)))-xx(ne(i,jp3(3))))*
     &       (yy(ne(i,jp3(1)))-yy(ne(i,jp3(3))))
 80   continue
c
c	write(*,*) 'a-loop80'
      do 100 n=1,npn
c	write(*,*) 'in loop100-1'
      x=xxn(n)
      y=yyn(n)
      almin=1d10
      do 140 i=1,nel
      if(arr(i).lt.1d-30) goto 140
        ifg=0
c1        alm=-1d0
	  alm1=1000d0
      alm2=-1000d0
c1        do 150 j=1,3
c        al(j)=((xx(ne(i,jp3(jp1(j))))-x)*(yy(ne(i,jp3(jp2(j))))-y)
c     &        -(xx(ne(i,jp3(jp2(j))))-x)*(yy(ne(i,jp3(jp1(j))))-y))
c     &		/arr(i)
	al(1)=((x-xx(ne(i,jp3(3))))*(yy(ne(i,jp3(2)))-yy(ne(i,jp3(3))))
     &      -(xx(ne(i,jp3(2)))-xx(ne(i,jp3(3))))*(y-yy(ne(i,jp3(3)))))
     &      /arr(i)
	al(2)=(-(x-xx(ne(i,jp3(3))))*(yy(ne(i,jp3(1)))-yy(ne(i,jp3(3))))
     &       +(xx(ne(i,jp3(1)))-xx(ne(i,jp3(3))))*(y-yy(ne(i,jp3(3)))))
     &      /arr(i)
	al(3)=1d0-al(1)-al(2)
c      if(al(j).lt.-5d0 .or. al(j).gt.5d0) goto 140
	do 152 j=1,3
	alm1=dmin1(alm1,al(j))
 152	alm2=dmax1(alm2,al(j))
c1        if(al(j).lt.0d0) alm=dmax1(alm,-al(j))
c1        if(al(j).gt.1d0) alm=dmax1(alm,al(j)-1d0)
c1 150  continue
c1       if(alm.lt.0d0) goto 200
c       if(alm1.ge.0 .and. alm2.le.1d0) goto 200
       if(alm1.gt.-1d-12 .and. alm2.lt.1d0+1d-12) goto 200
       alm=dmax1(-alm1,alm2-1d0)
       if(alm.lt.almin) then
         almin=alm
         imin=i
         do 145 j=1,3
 145     al0(j)=al(j)
       endif
 140  continue
c
      open(1,file="grdfc5.res",access='append')
      write(1,*)"didn't pass"
      close(1)
c
      i=imin
      do 155 j=1,3
 155  al(j)=al0(j)
c      write(*,*) 'n,i,almin=',n,i,almin
 200  ttn(n)=0d0
c      do 210 jj=1,6
c 210  ttn(n)=ttn(n)+al(1)*(2d0*al(1)-1d0)*tt(ne(i,1))
c	write(*,*) 'in loop100-2'
cc      ttn(n)=al(1)*(2d0*al(1)-1d0)*tt(ne(i,1))       !要素内節点番号変更
cc    &	      +4d0*al(1)*al(2)*tt(ne(i,2))
cc     &	      +al(2)*(2d0*al(2)-1d0)*tt(ne(i,3))
cc    &	      +4d0*al(2)*(1d0-al(1)-al(2))*tt(ne(i,4))
cc   &	   +(1d0-al(1)-al(2))*(2d0*(1d0-al(1)-al(2))-1d0)*tt(ne(i,5))
cc     &	      +4d0*al(1)*(1d0-al(1)-al(2))*tt(ne(i,6))
      ttn(n)=al(1)*(2d0*al(1)-1d0)*tt(ne(i,1))
     &	      +4d0*al(1)*al(2)*tt(ne(i,4))
     &	      +al(2)*(2d0*al(2)-1d0)*tt(ne(i,2))
     &	      +4d0*al(2)*(1d0-al(1)-al(2))*tt(ne(i,5))
     &	   +(1d0-al(1)-al(2))*(2d0*(1d0-al(1)-al(2))-1d0)*tt(ne(i,3))
     &	      +4d0*al(1)*(1d0-al(1)-al(2))*tt(ne(i,6))
c      ttn(n)=al(1)*tt(ne(i,1))+al(2)*tt(ne(i,2))+al(3)*tt(ne(i,3))
c
      open(1,file="grdfc5.res",access='append')
      write(1,*)n,",",i,",",sngl(arr(i)),",",sngl(al(1))
     &      ,",",sngl(al(2)),",",sngl(al(3)),",",sngl(ttn(n))
      close(1)
c

      if(ttn(n).lt.0d0) ttn(n)=0d0
 100  continue
c	write(*,*) 'a loop100'
c
      nel=neln
      np=npn
      nbw1=nbwn
      do 300 i=1,nel
      do 300 j=1,6
 300  ne(i,j)=nen(i,j)
      do 310 i=1,np
      xx(i)=xxn(i)
      yy(i)=yyn(i)
c      uu(i)=uun(i)
c      vv(i)=vvn(i)
c      pp(i)=ppn(i)
 310  tt(i)=ttn(i)
c
      do i=1,nbw(2)
	tb(i)=tt(ne(nelbw(2,i),nebw(2,i)))
	enddo
	do i=1,nel
	if(irx(i).eq.1)then
	do j=1,6
	tt(ne(i,j))=0d0
	ttn(ne(i,j))=0d0
	enddo
	endif
	enddo
	do i=1,nbw(2)
	tt(ne(nelbw(2,i),nebw(2,i)))=tb(i)
	ttn(ne(nelbw(2,i),nebw(2,i)))=tb(i)
	enddo
c
cn      do 500 i=1,nel
cn      arr2=(xx(ne(i,2))-xx(ne(i,1)))*(yy(ne(i,3))-yy(ne(i,1)))
cn     &         -(xx(ne(i,3))-xx(ne(i,1)))*(yy(ne(i,2))-yy(ne(i,1)))
cn      do 500 j=1,3
cn        j1=j+1
cn        if(j1.eq.4) j1=1
cn        zl=dsqrt((xx(ne(i,j))-xx(ne(i,j1)))**2
cn     &                           +(yy(ne(i,j))-yy(ne(i,j1)))**2)
cn        hl=arr2/zl
cn        zh(i,j)=zl/hl
cn 500  continue
c
c      call axnod(np,nel,ne,xx,yy,id,ind,ine,ico,ido)
c      call sufnod3(np,nel,ne,id,xx,yy,ico,ido,ind,ine)
c	call nelsuf(nel,ne,ind,ine)
      return
      end
c
	subroutine jcbd(new,im,xxw,yyw,detj,ind,ine)
      implicit double precision (a-h,o-z)
      dimension new(ine,6),xxw(ind),yyw(ind)
c      dimension ej(2,6)
c      data (ej(1,j),j=1,6)/3d0,2d0,0d0,-2d0,-3d0,0d0/
c      data (ej(2,j),j=1,6)/0d0,2d0,3d0,0d0,-3d0,-2d0/
c...
c	do 10 im=1,neldn
        vj11=0d0
	  vj12=0d0
        vj21=0d0
	  vj22=0d0
	  detj=0d0
c
cc      vj11=-xxw(new(im,1))-3d0*xxw(new(im,5))+4d0*xxw(new(im,6))  !要素内の節点番号変更
cc	  vj12=-yyw(new(im,1))-3d0*yyw(new(im,5))+4d0*yyw(new(im,6))
cc	  vj21=-xxw(new(im,3))-3d0*xxw(new(im,5))+4d0*xxw(new(im,4))
cc	  vj22=-yyw(new(im,3))-3d0*yyw(new(im,5))+4d0*yyw(new(im,4))
	  vj11=-xxw(new(im,1))-3d0*xxw(new(im,3))+4d0*xxw(new(im,6))
	  vj12=-yyw(new(im,1))-3d0*yyw(new(im,3))+4d0*yyw(new(im,6))
	  vj21=-xxw(new(im,2))-3d0*xxw(new(im,3))+4d0*xxw(new(im,5))
	  vj22=-yyw(new(im,2))-3d0*yyw(new(im,3))+4d0*yyw(new(im,5))
	  detj=vj11*vj22-vj12*vj21
	  detj=dabs(detj)
	  return
      end
c
        subroutine getnek(nelrs,nelre,ne2,nek,nelk,npk,
     &               ip,jp,ine,ind)
	  implicit double precision (a-h,o-z)
	  dimension ne2(ine,6),nelrs(21),nelre(21)
	  dimension nek(ine,6),ip(ind),jp(ind)
c
	  lr=2
	  nelk=0
	  do 10 i=nelrs(lr),nelre(lr)
	  nelk=nelk+1
	  do 10 j=1,6
10      nek(nelk,j)=ne2(i,j)
c
        npk=0
	  do 20 i=1,nelk
	  do 30 j=1,6
        do k=1,npk
	      if(nek(i,j).eq.ip(k))goto 30
	    enddo
	    npk=npk+1
	    ip(npk)=nek(i,j)
		jp(nek(i,j))=npk
30       continue
20       continue
c
        do i=1,nelk
	  do j=1,6
	  nek(i,j)=jp(nek(i,j))
	  enddo
	  enddo
        return
	  end
c
      subroutine btw2(new,nuv,kuv,uvb,tt,jww,jp,ine,ind)
      implicit double precision (a-h,o-z)
      dimension new(ine,6),kuv(2001,10),uvb(2001,10),tt(ind),nuv(10)
      dimension jww(ind),jp(ind)
      common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c--
      lr=2
      nuv(lr)=0
      do 10 i=1,nbw(lr)+1
	  if(ncomw(lr,i).ne.0)goto 10
        k=jww(jp(new(nelbw(lr,i),nebw(lr,i))))
        nuv(lr)=nuv(lr)+1
        kuv(nuv(lr),lr)=k
        uvb(nuv(lr),lr)=1.0d0
        tt(k)=uvb(nuv(lr),lr)
 10   continue
      return
      end

      subroutine calw2(nel,new,xxw,yyw,t0,nelrs,nelre,cov,ind,ine)
	implicit double precision (a-h,o-z)  !初期濃度から増分した濃度の全モル数を求める
	dimension new(ine,6),xxw(ind),yyw(ind),t0(ind)
	dimension pn(6,7),el(4,2),ww(7)
	dimension nelrs(21),nelre(21)
	common /dprop/dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
	common /mass/dmas(9)
c
      lr=2     !considering area
c      id=1
c      pe=(1d0/visc(id))*dmas(id)  !Pe=Re*Sc
c	write(*,*)"pe,visc(id),dmas(id)",pe,visc(id),dmas(id)
c	zlam=1d0/pe
	pi=3.1415926535897d0     !倍精度の仮数部は約14桁
	ww(1)=-27d0/96d0  !重み係数
	ww(2)=-27d0/96d0
	ww(3)=25d0/96d0
	ww(4)=25d0/96d0
	el(1,1)=1d0/3d0	  !積分点
	el(2,1)=11d0/15d0
	el(3,1)=2d0/15d0
	el(4,1)=2d0/15d0
	el(1,2)=1d0/3d0
	el(2,2)=2d0/15d0
	el(3,2)=2d0/15d0
	el(4,2)=11d0/15d0
	do 10 i=1,4       !形状関数
	pn(1,i)=el(i,1)*(2d0*el(i,1)-1d0)
	pn(4,i)=4d0*el(i,1)*el(i,2)
	pn(2,i)=el(i,2)*(2d0*el(i,2)-1d0)
	pn(5,i)=4d0*el(i,2)*(1d0-el(i,1)-el(i,2))
	pn(3,i)=(1d0-el(i,1)-el(i,2))*(1d0-2d0*el(i,1)-2d0*el(i,2))
 10	pn(6,i)=4d0*el(i,1)*(1d0-el(i,1)-el(i,2))

c--
	cov=0d0
c	do 100 im=1,nel
      do 100 im=nelrs(lr),nelre(lr)	   !領域2の全要素
	call jcbd(new,im,xxw,yyw,detj,ind,ine)	  !ヤコビアンマトリックス
	rav=0d0
	covp=0d0
	do 110 j=1,6						
 110	rav=rav+xxw(new(im,j))/6d0
c      rav=rav-1d-1      !領域1の厚み除く
	do 150 i=1,4	  !積分点
	do 170 j=1,6	  !形状関数
 170	covp=covp+2d0*pi*rav*ww(i)*pn(i,j)*detj*t0(new(im,j))  !2π*Rnel*Snel
 150	continue
	cov=cov+covp
 100	continue
	return
	end
c
      subroutine vdrop(nel,ne,xx,yy,idx,vvv,ind,ine)
      implicit double precision (a-h,o-z)
      dimension ne(ine,3),xx(ind),yy(ind),idx(0:ine)
	dimension ar(ine)
c
      pi=3.1415926535897d0
	vvv=0d0
      do 100 k=1,nel	   !index1の領域の全要素
	if(idx(k).ne.1)goto 100
          ar(k)=(xx(ne(k,2))-xx(ne(k,1)))*(yy(ne(k,3))-yy(ne(k,1)))     !ar:面積*2
     &      -(xx(ne(k,3))-xx(ne(k,1)))*(yy(ne(k,2))-yy(ne(k,1)))
	      arr=ar(k)/2d0   !arr:面積
	    rav=(xx(ne(k,1))+xx(ne(k,2))+xx(ne(k,3)))/3d0
	    vvv=vvv+2d0*pi*rav*arr         !2π*Rnel*Snel
 100	continue
      return
	end
c
      subroutine calw3(nel,new,xxw,yyw,t0,nelrs,nelre,cov,ind,ine)
	implicit double precision (a-h,o-z)  !液滴内の全モル数を求める
	dimension new(ine,6),xxw(ind),yyw(ind),t0(ind)
	dimension pn(6,7),el(4,2),ww(7)
	dimension nelrs(21),nelre(21)
	common /dprop/dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
	common /mass/dmas(9)
c
      lr=2     !considering area
c      id=1
c      pe=(1d0/visc(id))*dmas(id)  !Pe=Re*Sc
c	write(*,*)"pe,visc(id),dmas(id)",pe,visc(id),dmas(id)
c	zlam=1d0/pe
	pi=3.1415926535897d0     !倍精度の仮数部は約14桁
	ww(1)=-27d0/96d0  !重み係数
	ww(2)=-27d0/96d0
	ww(3)=25d0/96d0
	ww(4)=25d0/96d0
	el(1,1)=1d0/3d0	  !積分点(4点)
	el(2,1)=11d0/15d0
	el(3,1)=2d0/15d0
	el(4,1)=2d0/15d0
	el(1,2)=1d0/3d0
	el(2,2)=2d0/15d0
	el(3,2)=2d0/15d0
	el(4,2)=11d0/15d0
	do 10 i=1,4       !形状関数
	pn(1,i)=el(i,1)*(2d0*el(i,1)-1d0)
	pn(4,i)=4d0*el(i,1)*el(i,2)
	pn(2,i)=el(i,2)*(2d0*el(i,2)-1d0)
	pn(5,i)=4d0*el(i,2)*(1d0-el(i,1)-el(i,2))
	pn(3,i)=(1d0-el(i,1)-el(i,2))*(1d0-2d0*el(i,1)-2d0*el(i,2))
 10	pn(6,i)=4d0*el(i,1)*(1d0-el(i,1)-el(i,2))

c--
	cov=0d0
c	do 100 im=1,nel
      do 100 im=nelrs(lr),nelre(lr)	   !領域2の全要素
	call jcbd(new,im,xxw,yyw,detj,ind,ine)	  !ヤコビアンマトリックス
	rav=0d0
	covp=0d0
	do 110 j=1,6						
 110	rav=rav+xxw(new(im,j))/6d0
c      rav=rav-1d-1      !領域1の厚み除く
	do 150 i=1,4	  !積分点
	do 170 j=1,6	  !形状関数
 170	covp=covp+2d0*pi*rav*ww(i)*pn(i,j)*detj*(1d0+t0(new(im,j)))
 150	continue
	cov=cov+covp
 100	continue
	return
	end
c
c
      subroutine btw22(new,nuv,kuv,uvb,tt,ine,ind)
      implicit double precision (a-h,o-z)
      dimension new(ine,6),kuv(2001,10),uvb(2001,10),tt(ind),nuv(10)
      common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c--
      lr=2
      nuv(lr)=0
      do 10 i=1,nbw(lr)+1
	  if(ncomw(lr,i).ne.0)goto 10
        k=new(nelbw(lr,i),nebw(lr,i))
        nuv(lr)=nuv(lr)+1
        kuv(nuv(lr),lr)=k
        uvb(nuv(lr),lr)=1d0
        tt(k)=uvb(nuv(lr),lr)
 10   continue
      return
      end
c
      subroutine btwvap(new,nuv,kuv,uvb,ttk,tt,jww,jp,ine,ind)
      implicit double precision (a-h,o-z)
      dimension new(ine,6),kuv(2001,10),uvb(2001,10),nuv(10)
	dimension tt(ind),ttk(ind)
      dimension jww(ind),jp(ind)
      common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c--
      lr=2
      nuv(lr)=0
      do 10 i=1,nbw(lr)
	  if(ncomw(lr,i).ne.0)goto 10
        j=new(nelbw(lr,i),nebw(lr,i))			      !np_old
        k=jww(jp(new(nelbw(lr,i),nebw(lr,i))))	  !np_new
        nuv(lr)=nuv(lr)+1
        kuv(nuv(lr),lr)=k
        uvb(nuv(lr),lr)=tt(j)
        ttk(k)=uvb(nuv(lr),lr)
 10   continue
      return
      end
c
c------------------2次要素の自由表面上の境界の濃度を求める-----------------
        subroutine boundfsw(new,uu,vv,tt,px2,py2,ind,ine)
        implicit double precision (a-h,o-z)
	  dimension new(ine,6),tt(ind),uu(ind),vv(ind)
        dimension px2(ind),py2(ind)      !法線ベクトルの方向余弦nr,nz
	  dimension flx(1001)              !蒸発速度
	  dimension jj(0:10),hi(0:10),vf(0:10)
        common /region/ nr,nb(21),ndx(0:21),nelb(21,-1:1001),
     &         neb(21,-1:1001),ncom(21,1001),ncop(21,1001),nset,kset(21)
        common /dprop/ dens(9),visc(9),svc(9),gy(9),sten(-9:9,-9:9)
        common /cpoint/ ncp,kcp(101),icp(101),jcp(101),
     &       st(101),ancl(101),icl(101),icla(101),staa(101),strr(101)
     &       ,icoml(-9:9,-9:9,-9:9)
        common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c	  data pi/3.14159265358979d0/      !32384626433832795   倍精度の仮数部は約17桁
c
c------------------自由表面上の濃度を蒸発速度，法線方向の移動速度--------------
c------------------分子量，初期濃度から求める----------------------------------
        ca0=1d-2         !溶質[mol/m3]の濃度
c      call evpr(new,tt,jww,jp,flx,ine,ind)
	  write(*,*)"end evpr"
	    lr=2   !領域
          do 300 j=1,nbw(lr),2
	      write(*,*)"nbw(2),j",nbw(lr),j
	      if(ncomw(lr,j).ne.0)goto 300    !自由界面のみ考慮
            jj(0)=new(nelbw(lr,j),nebw(lr,j))     !ne2(ine,6)
            jj(1)=new(nelbw(lr,j+2),nebw(lr,j+2))
            jj(2)=new(nelbw(lr,j),nebw(lr,j)+3)    !辺j-j+2上の節点(j=4,5,6)
c------------------------1辺の境界濃度を求める---------------------
            px2(jj(2))=(px2(jj(0))+px2(jj(1)))/2d0
	      py2(jj(2))=(py2(jj(0))+py2(jj(1)))/2d0
            uu(jj(2))=(uu(jj(0))+uu(jj(1)))/2d0
            vv(jj(2))=(vv(jj(0))+vv(jj(1)))/2d0
            write(*,*)"jj(0-2)",(jj(i),i=0,2)
	      write(*,*)"uu(0-2)",sngl(uu(jj(0)))
     &                         ,sngl(uu(jj(1))),sngl(uu(jj(2)))
	      write(*,*)"vv(0-2)",sngl(vv(jj(0)))
     &                         ,sngl(vv(jj(1))),sngl(vv(jj(2)))
	      write(*,*)"px2(0-2)",sngl(px2(jj(0)))
     &                          ,sngl(px2(jj(1))),sngl(px2(jj(2)))
	      write(*,*)"py2(0-2)",sngl(py2(jj(0)))
     &                          ,sngl(py2(jj(1))),sngl(py2(jj(2)))
	      do 800 k=0,2
800         hi(k)=0d0
            do 500 i=0,1     !自由表面における法線方向の移動速度を求める
c           角度θを求める
            if(dabs(uu(jj(i))).lt.1d-12 			!u=0,v=0
     &   .and. dabs(vv(jj(i))).lt.1d-12)goto 400
c
            vsin=vv(jj(i))/dsqrt(uu(jj(i))**2+vv(jj(i))**2)
	      vcos=uu(jj(i))/dsqrt(uu(jj(i))**2+vv(jj(i))**2)
	      if(dabs(vsin).lt.1d-12)shita=dacos(vcos)   !v=0
	      if(dabs(vcos).lt.1d-12)shita=dasin(vsin)	 !u=0
	      if(dabs(vsin).ge.1d-12 .and. dabs(vcos).ge.1d-12)
     &                                   shita=datan(vsin/vcos)
	      write(*,*)"shita",shita
c		  角度αを求める
            if(dabs(px2(jj(i))).lt.1d-12 			!px2=0,py2=0
     &   .and. dabs(py2(jj(i))).lt.1d-12)goto 400
c
            psin=py2(jj(i))/dsqrt(px2(jj(i))**2+py2(jj(i))**2)
	      pcos=py2(jj(i))/dsqrt(px2(jj(i))**2+py2(jj(i))**2)
	      if(dabs(psin).lt.1d-12)alpha=dacos(pcos)   !py2=0
	      if(dabs(pcos).lt.1d-12)alpha=dasin(psin)	 !px2=0
	      if(dabs(psin).ge.1d-12 .and. dabs(pcos).ge.1d-12)
     &                                   alpha=datan(psin/pcos)
	      write(*,*)"alpha",alpha
c
	      velcy=dsqrt(uu(jj(i))**2+vv(jj(j))**2)
            vf(i)=velcy*dcos(shita-alpha)    !自由表面における法線方向の移動速度を求める
	      write(*,*)"velcy,vf(",i,")",velcy,vf(i)
500		  continue
c
	      vf(2)=(vf(0)+vf(1))/2d0               !辺上の速度を求める
	      write(*,*)"vf",(vf(i),i=0,2)
	      do 600 i=0,2
	      hi(i)=0d0
	      if(dabs(vf(i)).lt.1d-12)goto 600
            hi(i)=flx(j+i)/vf(i)/fwm(0)/ca0     !蒸発の項を求める(分子量は水の値)
600         continue
c
400         do 700 i=0,2
            write(*,*)"vf(",i,")","hi(",i,")","tt(",jj(i),")",
     &                 sngl(vf(i)),sngl(hi(i)),sngl(tt(jj(i)))
700	      tt(jj(i))=tt(jj(i))+hi(i)                   !1d0→Ca1に修正
c
300       continue
c
        write(*,*)"pass 300"
        open(1,file="them.res")
	  rewind(1)
	  lr=2
	  write(1,*)"nbw(2)",nbw(2)
	  do i=1,nbw(2)
	  k=new(nelbw(lr,i),nebw(lr,i))
	  write(1,*)k,sngl(tt(k))
	  enddo
	  close(1)
       write(*,*) 'boundfsw/end'
c
        return
        end
c
      subroutine bunpu2(np2,tt,ind)
	implicit double precision (a-h,o-z)
	dimension tt(ind)
c
c      ca0=1d-2
      write(*,*)"np2=",np2     !2次要素
      do 10 i=1,np2
c	  if(xx(i)+1d-12.lt.0d0)goto 10
c	  if(yy(i)+1d-12.lt.0d0)goto 10
	  tt(i)=0d0       !CA=(cA-cA0)/cA0    cA=cA0  at t=0
 10	continue
c
	return
	end
c
        subroutine boundfsw2(nel,ne,xx,yy,new,
     &                        tt,irx,time1,time2,ind,ine)
        implicit double precision (a-h,o-z)
	  dimension ne(ine,3),xx(ind),yy(ind),ar(ine)
	  dimension new(ine,6),tt(ind),irx(0:ine)
	  dimension flx(1001)              !蒸発速度
        common /doper/ vi,di,sta,str
        common /regionw/ nbw(21),nelbw(21,-1:1001),nebw(21,-1:1001),
     &                 ncomw(21,1001),ncopw(21,1001)
c	  data pi/3.14159265358979d0/      !32384626433832795   倍精度の仮数部は約17桁
c
        write(*,*)"in boundfsw2"
        arr=0d0
        lr=2
        do 8 k=1,nel       !液滴の断面積を求める
          if(irx(k).ne.lr) goto 8
          ar(k)=(xx(ne(k,2))-xx(ne(k,1)))*(yy(ne(k,3))-yy(ne(k,1)))
     &      -(xx(ne(k,3))-xx(ne(k,1)))*(yy(ne(k,2))-yy(ne(k,1)))
	    arr=arr+ar(k)
          if(ar(k).le.0d0) then
              write(*,*) 'ardiv / ar<0 / k,ar=',k,ar(k)
c             stop
          endif
 8      continue
        arr=arr*di**2     !有次元化[m2]
c
        ca0=1d-2         !溶質[mol/m3]の濃度
	  dt=(time2-time1)*di/vi       !経過時間[s]
c      call evpr(ne,tt,jww,jp,flx,ine,ind)
	    lr=2   !領域
          do 300 j=1,nbw(lr)
	      if(ncomw(lr,j).ne.0)goto 300    !自由界面のみ考慮
            j0=new(nelbw(lr,j),nebw(lr,j))     !頂点の節点
            j1=new(nelbw(lr,j+1),nebw(lr,j+1)) !辺上の節点
	      rl=dsqrt((xx(j0)-xx(j1))**2+(yy(j0)-yy(j1))**2)*di   !辺の長さ[m]
	      hl=1d-2/(1d0-1d-2)*flx(i)/fwm(0)*rl*dt   !l*w*dt=0.01/(100-0.01)*m/M*l*dt [mol/m]
	      tt(j0)=tt(j0)+hl/2d0/arr/ca0               !CA=CA1+l*w*dt/(3A*cA0)   [-]
		  tt(j1)=tt(j1)+hl/2d0/arr/ca0
	      write(*,*)"j,tt(",j0,")","tt(",j1,")",j,tt(j0),tt(j1)
	      pause
300       continue							         !分子量は水の値
c
        open(1,file="them.res")
	  rewind(1)
	  lr=2
	  write(1,*)"nbw(2)",nbw(2)
	  do i=1,nbw(2)
	  k=new(nelbw(lr,i),nebw(lr,i))
	  write(1,*)k,sngl(tt(k))
	  enddo
	  close(1)
        write(*,*)"boundfsw2/end"
c
        return
        end