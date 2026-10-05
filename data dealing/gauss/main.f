        implicit double precision(a-h,o-z)
        parameter(n=4,m=5)
        dimension a(n,m),iwork(n)
c
        a(1,1)= 3d0
        a(1,2)= 2d0
        a(1,3)= 7d0
        a(1,4)= 1d0
        a(2,1)= 1d0
        a(2,2)= 5d0
        a(2,3)= 1d0
        a(2,4)=-1d0
        a(3,1)= 4d0
        a(3,2)= 1d0
        a(3,3)= 3d0
        a(3,4)=-2d0
        a(4,1)= 1d0
        a(4,2)= 6d0
        a(4,3)= 4d0
        a(4,4)= 3d0
c
        a(1,5)= 8d0
        a(2,5)= 5d0
        a(3,5)= 7d0
        a(4,5)=13d0
c
        eps=1d-5
c
        call gauss(a,n,m,eps,iwork,ill)  !n×n行列, m=n+1 (右辺の値)
        if(ill.eq.1)goto 100
        do 200 i=1,n
200       write(6,10)i,(a(i,j),j=n+1,m)
10      format(1H ,10X,'X',I1,'=',4F8.4)

        write(*,*)3d0*3.5d0+2d0*1d0+7d0*(-1d0)+1d0*2.5d0- 8d0
        write(*,*)1d0*3.5d0+5d0*1d0+1d0*(-1d0)-1d0*2.5d0- 5d0
        write(*,*)4d0*3.5d0+1d0*1d0+3d0*(-1d0)-2d0*2.5d0- 7d0
        write(*,*)1d0*3.5d0+6d0*1d0+4d0*(-1d0)+3d0*2.5d0-13d0

c
100     write(*,*)"ill"
        stop
        end
c------------------------------------------------------------
c
c       ガウスの消去法
c
c------------------------------------------------------------
        subroutine gauss(a,n,m,eps,iwork,ill)
        implicit double precision(a-h,o-z)
        dimension a(n,m),iwork(n)
        ill=0
c
        do 10 i=1,n
10      iwork(i)=i
c
        do 17 k=1,n
c
c         行列(k-n, k-n)の係数aij の最大値aij,max
c         i=ip, j=iq
          dmax=dabs(a(k,k))
          ip=k
          iq=k
          do 11 j=k,n
          do 11 i=k,n
            if( dmax.ge.dabs(a(i,j)) )goto 11
            dmax=dabs(a(i,j))
            ip=i
            iq=j
11        continue
c
c         係数が小さかった時
          if(dmax.le.eps) goto 20
c
c         k列とiq列の入れ換え
c
          do 12 i=1,n
            w=a(i,k)
            a(i,k)=a(i,iq)
12          a(i,iq)=w
c
c         k行とip行の入れ換え
c
          do 13 j=k,m
            w=a(k,j)
            a(k,j)=a(ip,j)
13          a(ip,j)=w
c
c         iwork(入れ換え後の列)=入れ換え前の列
c         iwork(入れ換え前の列)=入れ換え後の列
c
          i=iwork(k)
          iwork(k)=iwork(iq)
          iwork(iq)=i
c
c         akj=akj / akk  (akk=1d0になる)
          do 14 j=k+1,m
14        a(k,j)=a(k,j)/a(k,k)
c
c         注目している行を引いていく
c         aij=aij - aik*(akj/akk)
          do 16 i=1,n   !行
            if(i.eq.k)goto 16  !注目している行は飛ばす(「他の行」から「注目している行」を引く)
            do 15 j=k+1,m       !列 : 行列の右上3角形だけ考える。左下3角形の係数は0になる
15          a(i,j)=a(i,j)-a(i,k)*a(k,j)     !aij=aij - aik*(akj/akk)
16       continue
17     continue
c
c      列を元に戻す
       do 19 j=n+1,m    !j=4+1,5
       do 18 i=1,n
         iw=iwork(i)
18       a(iw,n)=a(i,j)
         do 19 i=1,n
19     a(i,j)=a(i,n)    !a(i,5)=a(i,4)
       return
c
20     write(6,100)
       ill=1
       return
100    format(1H ,3X,'MATRIX IS ILL')
       end