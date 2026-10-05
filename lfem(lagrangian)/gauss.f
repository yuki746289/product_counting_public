        subroutine gauss(n,nbw,a,x,ind,ibw)
        implicit double precision (a-h,o-z)
        dimension a(ind,ibw),x(ind)
c
        do 2 i=1,n
        if(a(i,nbw).le.0) then
         write(*,*) '[gauss]',i,a(i,nbw)
         stop
        endif
 2      continue
c
        nbwm=nbw-1
        do 10 i=1,n
          aa=a(i,nbw)
          x(i)=x(i)/aa
          do 20 j=i,min(n,i+nbwm)
 20       a(i,nbw+j-i)=a(i,nbw+j-i)/aa
          do 40 k=i+1,min(n,i+nbwm)
            cc=a(k,nbw+i-k)
            do 30 j=i+1,min(n,i+nbwm)
 30         a(k,nbw+j-k)=a(k,nbw+j-k)-cc*a(i,nbw+j-i)
 40       x(k)=x(k)-cc*x(i)
 10     continue
        do 70 i=n-1,1,-1
        do 70 k=i+1,min(n,i+nbwm)
 70     x(i)=x(i)-a(i,nbw+k-i)*x(k)
        return
        end

c------------------------------------------------------------
c
c       ガウスの消去法
c
c------------------------------------------------------------
        subroutine gauss2(a,n,m,eps,iwork,ill)
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
20     continue
c       write(6,100)
       ill=1
c       write(*,*)"drmax:",k,"/",n,":",sngl(drmax)
       return
100    format(1H ,3X,'MATRIX IS ILL')
       end
