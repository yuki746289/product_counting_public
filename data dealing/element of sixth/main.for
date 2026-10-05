        implicit double precision(a-h,o-z)
c
        call calc1 !局座標による評価
        call calc2 !体積座標による評価
        stop
        end
c
c       calc1
c
        subroutine calc1
        implicit double precision(a-h,o-z)
        dimension gzi(8),eth(8),zth(8)
        gzi(1)=-1
        gzi(2)= 1
        gzi(3)= 1
        gzi(4)=-1
        gzi(5)=-1
        gzi(6)= 1
        gzi(7)= 1
        gzi(8)=-1
        eth(1)=-1
        eth(2)=-1
        eth(3)= 1
        eth(4)= 1
        eth(5)=-1
        eth(6)=-1
        eth(7)= 1
        eth(8)= 1
        zth(1)=-1
        zth(2)=-1
        zth(3)=-1
        zth(4)=-1
        zth(5)= 1
        zth(6)= 1
        zth(7)= 1
        zth(8)= 1
c
        open(1,file="grid1.res")
        rewind(1)
        gzai=-1.5d0
        etha=-1.5d0
        ztha=-1.5d0
        do 100 i=1,10
        do 200 j=1,10
        do 300 k=1,10
          zn1=(1d0+gzi(1)*gzai)*(1d0+eth(1)*etha)*(1d0+zth(1)*ztha)/8d0
          zn2=(1d0+gzi(2)*gzai)*(1d0+eth(2)*etha)*(1d0+zth(2)*ztha)/8d0
          zn3=(1d0+gzi(3)*gzai)*(1d0+eth(3)*etha)*(1d0+zth(3)*ztha)/8d0
          zn4=(1d0+gzi(4)*gzai)*(1d0+eth(4)*etha)*(1d0+zth(4)*ztha)/8d0
          zn5=(1d0+gzi(5)*gzai)*(1d0+eth(5)*etha)*(1d0+zth(5)*ztha)/8d0
          zn6=(1d0+gzi(6)*gzai)*(1d0+eth(6)*etha)*(1d0+zth(6)*ztha)/8d0
          zn7=(1d0+gzi(7)*gzai)*(1d0+eth(7)*etha)*(1d0+zth(7)*ztha)/8d0
          zn8=(1d0+gzi(8)*gzai)*(1d0+eth(8)*etha)*(1d0+zth(8)*ztha)/8d0
          zzn=zn1+zn2+zn3+zn4+zn5+zn6+zn7+zn8
          ztha=ztha+0.3d0
          write(1,*)sngl(gzai),sngl(etha),sngl(ztha),sngl(zzn)
          write(*,*)sngl(gzai),sngl(etha),sngl(ztha),sngl(zzn)
c          pause
300     continue
        ztha=-1.5d0
        etha=etha+0.3d0
200     continue
        etha=-1.5d0
        gzai=gzai+0.3d0
100     continue
        close(1)
        return
        end
c
c       calc2
c
        subroutine calc2
        implicit double precision(a-h,o-z)
        open(1,file="grid2.res")
        rewind(1)
        xl=0d0
        xh=1d0
        yl=0d0
        yh=1d0
        zl=0d0
        zh=1d0
        x=-1.5d0
        y=-1.5d0
        z=-1.5d0
        do 100 i=1,10
        do 200 j=1,10
        do 300 k=1,10
c         (1)6面体要素の体積の算出
          vl=1d0
c         (2)形状関数の算出
          zn1=dabs(xh-xl)*dabs(y -yl)*dabs(zh-zl)/3d0/vl  !面xz×高さy(low)
          zn2=dabs(xh-xl)*dabs(yh-yl)*dabs(z -zl)/3d0/vl  !面xy×高さz(low)
          zn3=dabs(xh- x)*dabs(yh-yl)*dabs(zh-zl)/3d0/vl  !面yz×高さx(high)
          zn4=dabs(xh-xl)*dabs(yh-yl)*dabs(zh- z)/3d0/vl  !面xy×高さz(high)
          zn5=dabs(x -xl)*dabs(yh-yl)*dabs(zh-zl)/3d0/vl  !面yz×高さx(low)
          zn6=dabs(xh-xl)*dabs(yh- y)*dabs(zh-zl)/3d0/vl  !面xz×高さy(high)
          zzn=zn1+zn2+zn3+zn4+zn5+zn6
          z=z+0.3d0
          write(1,*)sngl(x),sngl(y),sngl(z),sngl(zzn)
          write(*,*)sngl(x),sngl(y),sngl(z),sngl(zzn)
c          pause
300     continue
        z=-1.5d0
        y=y+0.3d0
200     continue
        y=-1.5d0
        x=x+0.3d0
100     continue
        close(1)
        stop
        end