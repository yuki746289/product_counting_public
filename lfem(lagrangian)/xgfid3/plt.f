c---------------------------------------------------------
c
c     節点を登録する
c
c---------------------------------------------------------
       subroutine pltsetnp(np,xx,yy,zz,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind)
c
       open(1,file='setnp.res')
       rewind(1)
       write(1,*)np
       do 200 kp=1,np
200    write(1,*)sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
c
       return
       end
c---------------------------------------------------------
c
c     表示する要素を追加する
c
c---------------------------------------------------------
       subroutine pltsetne(nele,ne,size,icolor,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),ne0(ine,4)
c
       call fcolor(icolor,red,green,blue)
c
c      データを書く
       open(1,file='setne.res')
       rewind(1)
c
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)nele
       do 300 i=1,nele      !追加のデータ
300    write(1,*)(ne(i,j),j=1,4)
       close(1)
c
       return
       end
c---------------------------------------------------------
c
c     表示する面を追加する
c
c---------------------------------------------------------
       subroutine pltsetme(nels,ms,size,icolor,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ms(ine,3)
c
       call fcolor(icolor,red,green,blue)
c
c      データを書く
       open(1,file='setme.res')
       rewind(1)
c
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)nels
       do 300 i=1,nels      !追加のデータ
300    write(1,*)ms(i,1),ms(i,2),ms(i,3)
       close(1)
c
       return
       end
c---------------------------------------------------------
c
c     表示する線を追加する
c
c---------------------------------------------------------
       subroutine pltsetline(k1,k2,kmax,size,icolor)
       implicit double precision(a-h,o-z)
c
       call fcolor(icolor,red,green,blue)
c
c      データを書く
       open(1,file='setline.res')
       rewind(1)
       write(1,*)k1,k2,kmax,size,red,green,blue
       close(1)
c
       return
       end
c---------------------------------------------------------
c
c     リセットする
c
c---------------------------------------------------------
       subroutine pltrset()
       implicit double precision(a-h,o-z)
c
c      要素のデータをリセットする
       open(1,file='setne.res')
       rewind(1)
       write(1,*)0e0
       write(1,*)0e0,0e0,0e0
       write(1,*)0
       close(1)
c
c      面のデータをリセットする
       open(1,file='setme.res')
       rewind(1)
       write(1,*)0e0
       write(1,*)0e0,0e0,0e0
       write(1,*)0
       close(1)
c
c      データを書く
       open(1,file='setline.res')
       rewind(1)
       write(1,*)1     !1:描画する,0:描画しない
       write(1,*)1,1,1,1d0,1d0,1d0,1d0  !k1,k2,kmax,size,red,green,blue
       close(1)
c
       return
       end
c---------------------------------------------------------
c
c      面kelsと節点を共有する面を描く
c
c---------------------------------------------------------
       subroutine pltsufms(kels,ns,xs,ys,zs,nels,ms)
       include "header.h"
       dimension ms(ine,3),ms0(ine,3)
       dimension xs(ind),ys(ind),zs(ind)
c
       call pltrset()
       call pltsetnp(ns,xs,ys,zs,ind)
c
       ks1=ms(kels,1)
       ks2=ms(kels,2)
       ks3=ms(kels,3)
c
       nel0=1
       ms0(nel0,1)=ks1
       ms0(nel0,2)=ks2
       ms0(nel0,3)=ks3
c
       do 100 i=1,nels
         if(i.eq.kels)goto 100
         k1=ms(i,1)
         k2=ms(i,2)
         k3=ms(i,3)
         if(k1.eq.ks1 .or. k1.eq.ks2 .or. k1.eq.ks3 .or.
     &      k2.eq.ks1 .or. k2.eq.ks2 .or. k2.eq.ks3 .or.
     &      k3.eq.ks1 .or. k3.eq.ks2 .or. k3.eq.ks3)then
           nel0=nel0+1
           ms0(nel0,1)=k1
           ms0(nel0,2)=k2
           ms0(nel0,3)=k3
           goto 100
         endif
100    continue
c
       call pltsetme(nel0,ms0,3d0,2,ind,ine)
c
       return
       end



c---------------------------------------------------------
c
c     要素を登録する
c
c---------------------------------------------------------
       subroutine pltne(np,nele,ne,xx,yy,zz,size,icolor,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4)
       dimension xx(ind),yy(ind),zz(ind)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pne.res',position='REWIND')
       open(1,file='pne.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)np,nele
       do 100 i=1,nele
100    write(1,*)(ne(i,j),j=1,4)
       do 200 kp=1,np
200    write(1,*)sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))

       return
       end
c
c------------------------------------
c
       subroutine pltnet(np,na,ma,ne,xx,yy,zz,size,icolor,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),ma(ine)
       dimension xx(ind),yy(ind),zz(ind)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pnet.res',position='REWIND')
       open(1,file='pnet.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)np,na
       do 100 ka=1,na
100    write(1,*)(ne(ma(ka),j),j=1,4)
       do 200 kp=1,np
200    write(1,*)sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
       close(1)
c
       return
       end
c
c------------------------------------
c
       subroutine pltms(ns,nels,ms,xs,ys,zs,size,icolor,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ms(ine,3)
       dimension xs(ind),ys(ind),zs(ind)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pms.res',position='REWIND')
       open(1,file='pms.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)ns,nels
       do 100 i=1,nels
100    write(1,*)(ms(i,j),j=1,3)
       do 200 ks=1,ns
200    write(1,*)sngl(xs(ks)),sngl(ys(ks)),sngl(zs(ks))

       return
       end
c
c------------------------------------
c
       subroutine pltnp(np,xx,yy,zz,size,icolor,ind)
       implicit double precision(a-h,o-z)
       dimension xx(ind),yy(ind),zz(ind)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pnp.res',position='REWIND')
       open(1,file='pnp.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
        write(1,*)np
       do 100 kp=1,np
100    write(1,*)sngl(xx(kp)),sngl(yy(kp)),sngl(zz(kp))
       close(1)
c
       return
       end
c
c------------------------------------
c
       subroutine pltnt(x,y,z,size,icolor)
       implicit double precision(a-h,o-z)
c
       call fcolor(icolor,red,green,blue)
c
c       open(1,file='pnt.res',position='REWIND')
       open(1,file='pnt.res')
       rewind(1)
       write(1,*)size
       write(1,*)red,green,blue
       write(1,*)sngl(x),sngl(y),sngl(z)
       close(1)
c
       return
       end
c
c------------------------------------
c
       subroutine fcolor(icolor,red,green,blue)
       implicit double precision(a-h,o-z)
c      1:白
c      2:赤
c      3:緑
c      4:青
c      5:黄色

       if(icolor.eq.1)then  !白
         red  =1d0
         green=1d0
         blue =1d0
       elseif(icolor.eq.2)then  !赤
         red  =1d0
         green=0d0
         blue =0d0
       elseif(icolor.eq.3)then  !緑
         red  =0d0
         green=1d0
         blue =0d0
       elseif(icolor.eq.4)then  !青
         red  =0d0
         green=0d0
         blue =1d0
       elseif(icolor.eq.5)then  !黄色
         red  =0d0
        green=1d0
          blue =1d0
       endif

       return
       end