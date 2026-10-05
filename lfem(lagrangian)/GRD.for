       subroutine findsuf(np,ne,xx,yy,zz,ngx,ngy,ngz,
     &                     xg,yg,zg,ms,mx,my,mz,ind,ine)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind)
       dimension xg(10001),yg(10001),zg(10001)
       dimension ms(ind),mx(ind),my(ind),mz(ind)
       dimension lx(1001),ly(1001),lz(1001),lp(ind)
       dimension kxx(1001),kyy(1001),kz(1001)
c
       nm0=1
       lx(nm0)=
       ly(nm0)=
       lz(nm0)=
       nx0=
       ny0=
       nz0=
       do 1000 km=1,nm0
         call adjust(km,lx,ly,lz,lp,xg,yg,zg,nx0,ny0,nz0,
     &               xx,yy,zz,kn,kx,ky,kz,kp,ind,ine)
         
       