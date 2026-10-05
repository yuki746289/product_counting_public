c--------------------------------------------------------------
c
c      計算条件を読み込む
c
c--------------------------------------------------------------
       subroutine data(np,nele,ne,xx,yy,zz,nbw,iyy,ivs,idx,rfile)
       implicit double precision(a-h,o-z)
       parameter(ind=19999,ine=99999,ibw=1999)             
       dimension ne(ine,4),xx(ind),yy(ind),zz(ind)
       dimension idx(0:ine)
       dimension idtmp(1001),jdtmp(1001)
c------contidion file---------------------------------------------------------------------------------
       dimension ug(ind),vg(ind),wg(ind),cg(ind),tg(ind)   !第二種の境界条件
       dimension ch(ind),th(ind)                           !第三種の境界条件
       dimension scm(ind),sct(ind)                         !生成項φ[1/m3/s]
c-----------------------------------------------------------------------------------------------------
       character*20 gfile,gfile1,gfile2,ifile,bfile,bfile1,rfile
       common /time/dt,timax,tsv,ttmp
       common /prop/re(1001),st(1001,1001),sv(1001),gv(1001),vv0,xx0
       common /pmass/sc(1001)
       common /pheat/pr(1001),znu(1001,1001)
       common /nxyz/nx,ny,nz
c------contidion file---------------------------------------------------------------------------------       
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
       common/bnd1/cb(49999),tb(49999)     !第一種の境界条件
       common/init/uini,vini,wini,pini
       common/initch/cini,hini
c-----------------------------------------------------------------------------------------------------            
c
       open(1,file='input.txt')
       rewind(1)
       read(1,*)rfile
       close(1)
       call addchr(rfile,'.inf',ifile)
c
       write(*,*)'information file   ',ifile
       open(1,file=ifile)
       rewind(1)
       read(1,*)iyy   !0:最初から,1:続きから
       read(1,*)ivs   !0:描画しない,1:描画する
       if(iyy.eq.0)write(*,*)'calc from the beginning iyy',iyy
       if(iyy.eq.1)write(*,*)'calc from the middle    iyy',iyy
       if(ivs.eq.0)write(*,*)'no display              ivs',ivs
       if(ivs.eq.1)write(*,*)'display periodically    ivs',ivs
c
       read(1,*)dt,tmax,tsv,ttmp
       read(1,*)nx,ny,nz
       read(1,*)v0,x0
       read(1,*)gfile
       call addchr(gfile,'.1',gfile1)
       call addchr(gfile,'.2',gfile2)
c      *.1ファイルの読み込み
       open(2,file=gfile1)
       rewind(2)
       read(2,*)np
       do 100 i=1,np
100    read(2,*)none,kp,xx(kp),yy(kp),zz(kp)
       close(2)
c      *.2ファイルの読み込み
       open(2,file=gfile2)
       rewind(2)
       read(2,*)nele,np,nbw
       write(*,*)"grid file   ",gfile1,gfile2
       do 200 i=1,nele
200    read(2,*)kele,ne(kele,1),ne(kele,2),ne(kele,3),ne(kele,4),
     &          xx(ne(kele,1)),yy(ne(kele,1)),zz(ne(kele,1)),
     &          xx(ne(kele,2)),yy(ne(kele,2)),zz(ne(kele,2)),
     &          xx(ne(kele,3)),yy(ne(kele,3)),zz(ne(kele,3)),
     &          xx(ne(kele,4)),yy(ne(kele,4)),zz(ne(kele,4)),idx(kele)
       close(2)
c
c      1/Re
       read(1,*)m1                               !m1:index番号の数
       if(m1.ne.0)read(1,*)(idtmp(i),i=1,m1)     !idtmp(m1):index番号
       if(m1.ne.0)read(1,*)(re(idtmp(i)),i=1,m1) !re(idtmp(m1)):各index番号の領域の1/Re
c       
c      1/We
       read(1,*)m2                               !m2:index番号の数
       if(m2.ne.0)read(1,*)(idtmp(i),i=1,m2)     !idtmp(m2):index番号
       if(m2.ne.0)read(1,*)(jdtmp(i),i=1,m2)     !jdtmp(m2):index番号
       if(m2.ne.0)read(1,*)(st(idtmp(i),jdtmp(i)),i=1,m2)  !st(idtmp(m2),jdtmp(m2)):各index番号の領域の1/We
c
c      1/Ma
       read(1,*)m3                               !m3:index番号の数
       if(m3.ne.0)read(1,*)(idtmp(i),i=1,m3)     !idtmp(m3):index番号
       if(m3.ne.0)read(1,*)(sv(idtmp(i)),i=1,m3) !sv(idtmp(m3)):各index番号の領域の1/Ma
c
c      1/Fr
       read(1,*)m4                               !m4:index番号の数
       if(m4.ne.0)read(1,*)(idtmp(i),i=1,m4)     !idtmp(m4):index番号
       if(m4.ne.0)read(1,*)(gv(idtmp(i)),i=1,m4) !gv(idtmp(m4)):各index番号の領域の1/Fr
c
c      1/Sc
       read(1,*)m5                               !m5:index番号の数
       if(m5.ne.0)read(1,*)(idtmp(i),i=1,m5)     !idtmp(m5):index番号
       if(m5.ne.0)read(1,*)(sc(idtmp(i)),i=1,m5) !sc(idtmp(m5)):各index番号の領域の1/Sc
c
c      1/Pr
       read(1,*)m6                               !m6:index番号の数
       if(m6.ne.0)read(1,*)(idtmp(i),i=1,m6)     !idtmp(m6):index番号
       if(m6.ne.0)read(1,*)(pr(idtmp(i)),i=1,m6) !pr(idtmp(m6)):各index番号の領域の1/Pr
c
c      1/Nu
       read(1,*)m7                               !m7:index番号の数
       if(m7.ne.0)read(1,*)(idtmp(i),i=1,m7)     !idtmp(m7):index番号
       if(m7.ne.0)read(1,*)(jdtmp(i),i=1,m7)     !jdtmp(m7):index番号
       if(m7.ne.0)read(1,*)(znu(idtmp(i),jdtmp(i)),i=1,m7) !znu(idtmp(m7),jdtmp(m7)):各index番号の領域の1/Nu
c
       close(1)
c       write(*,*)"readed information file"
c-------------------------------------------------------------
c      Boundary Condition fileの読み込み
c-------------------------------------------------------------
       write(*,*)'boundary condition file              bd.inf'
       write(*,*)'boundary condition file              bdnum.inf'
       open(1,file='bd.inf')
       open(2,file='bdnum.inf')
       rewind(1)
       rewind(2)
c
c      第一種の境界条件
c       
c      velocity
       read(2,*)nb(1)
       do 2000 i=1,nb(1)
2000   read(1,*)nelb(1,i),neb(1,i),ub(i),vb(i),wb(i),pb(i)
c      mass
       read(2,*)nb(2)
       do 2100 i=1,nb(2)
2100   read(1,*)nelb(2,i),neb(2,i),cb(i)    !cb(nb(2)):濃度
c      temperature
       read(2,*)nb(3)
       do 2200 i=1,nb(3)
2200   read(1,*)nelb(3,i),neb(3,i),tb(i)    !tb(nb(3)):温度
c
c      第二種の境界条件
c
c      velocity gradient
       read(2,*)nb(4)
       do 2300 i=1,nb(4)
2300   read(1,*)nelb(4,i),neb(4,i),ug(i),vg(i),wg(i)  !せん断応力τ[N/m2]
c      mass gradient
       read(2,*)nb(5)
       do 2400 i=1,nb(5)
2400   read(1,*)nelb(5,i),neb(5,i),cg(i)     !cg(nb(5)):物質流束J[mol/m2/s]
c      temperature gradient
       read(2,*)nb(6)
       do 2500 i=1,nb(6)
2500   read(1,*)nelb(6,i),neb(6,i),tg(i)     !tg(nb(5)):熱流束q[J/m2/s]
c
c      第三種の境界条件
c
c      concentration
       read(2,*)nb(7)
       do 2600 i=1,nb(7)
2600   read(1,*)nelb(7,i),neb(7,i),ch(i)     !ch(nb(7)):バルク濃度Ch[mol/m3]
c      temperature
       read(2,*)nb(8)
       do 2700 i=1,nb(7)
2700   read(1,*)nelb(8,i),neb(8,i),th(i)     !th(nb(8)):バルク温度Th[K]
c
c      生成項
c
c      velocity(連続の式）
       read(2,*)nb(9)
       do 2800 i=1,nb(9)
2800   read(1,*)nelb(9,i),neb(9,i),scv(i)     !scv(nb(9)): 質量生成速度[kg/m3/s]
c      concentration
       read(2,*)nb(10)
       do 2900 i=1,nb(10)
2900   read(1,*)nelb(10,i),neb(10,i),scm(i)   !scm(nb(10)):溶質成分の生成速度[mol/m3/s]
c      temperature
       read(2,*)nb(11)
       do 3000 i=1,nb(11)
3000  read(1,*)nelb(11,i),neb(11,i),sct(i)    !sct(nb(11)):熱量の生成速度[J/m3/s]
      close(1)
      close(2)
c-------------------------------------------------------------
c      Initial Condition fileの読み込み
c-------------------------------------------------------------
       write(*,*)'initial  condition file on the flow  icd.inf'
       open(1,file='icd.inf')
       rewind(1)
       read(1,*)uini
       read(1,*)vini
       read(1,*)wini
       read(1,*)pini
       close(1)
c
       write(*,*)'initial  condition file on the mass  icdc.inf'
       open(1,file='icdc.inf')
       rewind(1)
       read(1,*)cini
       close(1)
c
       write(*,*)'initial  condition file on the heat  icdh.inf'
       open(1,file='icdh.inf')
       rewind(1)
       read(1,*)hini
       close(1)
c
       return
       end
c--------------------------------------------------------------
c
c      境界条件を修正する
c
c--------------------------------------------------------------
       subroutine bdmodify(kr,nn,ne,mn,ine,ind)
       implicit double precision(a-h,o-z)
       dimension ne(ine,4)
       dimension js(4,-1:3),mn(ind)
       common/bndvwl/nelb(100,9999),neb(100,9999),nb(100),
     &                pb(49999),ub(49999),vb(49999),wb(49999),scv(49999)
c
       js(1,-1)=2
       js(1,0)=3
       js(1,1)=1
       js(1,2)=2
       js(1,3)=3
c
       js(2,-1)=3
       js(2,0)=2
       js(2,1)=4
       js(2,2)=3
       js(2,3)=2
c
       js(3,-1)=4
       js(3,0)=1
       js(3,1)=3
       js(3,2)=4
       js(3,3)=1
c
       js(4,-1)=1
       js(4,0)=4
       js(4,1)=2
       js(4,2)=1
       js(4,3)=4
c
       nn=0
       do 10 i=1,nb(kr)
         kp1=ne(nelb(kr,i),js(neb(kr,i),1))     !nelb(1,nb(1)):第一種の境界条件を代入した要素の番号
         kp2=ne(nelb(kr,i),js(neb(kr,i),2))     !neb(1,nb(1)):第一種の境界条件を代入した要素の面の番号
         kp3=ne(nelb(kr,i),js(neb(kr,i),3))
         do j=1,nn
           if(kp1.eq.mn(j))goto 20             !nn:境界条件を代入する節点の数
         enddo                                 !mm(nn):境界条件を代入する節点の番号
         nn=nn+1
         mn(nn)=kp1
20       do j=1,nn
           if(kp2.eq.mn(j))goto 30
         enddo
         nn=nn+1
         mn(nn)=kp2
30       do j=1,nn
           if(kp3.eq.mn(j))goto 10
         enddo
         nn=nn+1
         mn(nn)=kp3
10     continue
       return
       end
c--------------------------------------------------------------
c
c      初期条件を入力する
c
c--------------------------------------------------------------
c
c     流動領域
c
      subroutine initf(np,nele,ne,u0,v0,w0,p0,uu,vv,ww,pp,idx,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,4)
      dimension u0(ind),v0(ind),w0(ind),p0(ind)
      dimension uu(ind),vv(ind),ww(ind),pp(ind)
      dimension idx(0:ine)
      common/init/uini,vini,wini,pini
c
      do 10 i=1,np
        u0(i)=0d0
        v0(i)=0d0
        w0(i)=0d0
        p0(i)=0d0
        uu(i)=0d0
        vv(i)=0d0
        ww(i)=0d0
        pp(i)=0d0
10    continue
c
      do 100 kele=1,nele
        if(idx(kele).eq.1 .or. idx(kele).eq.2 .or.  !index番号が流動領域の場合
     &     idx(kele).eq.3 .or. idx(kele).eq.4)then
          do 200 j=1,4
            kp=ne(kele,j)
            u0(kp)=uini
            v0(kp)=vini
            w0(kp)=wini
            p0(kp)=pini
200       continue
        endif
100   continue
      return
      end
c
c     物質移動領域
c
      subroutine initc(np,nele,ne,c0,cc,idx,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,4),c0(ind),cc(ind)
      dimension idx(0:ine)
      common/initch/cini,hini
c
      do 10 i=1,np
        c0(i)=0d0
        cc(i)=0d0
10    continue
c
      do 100 kele=1,nele
        if(idx(kele).eq.2 .or. idx(kele).eq.4 .or.  !index番号が物質移動領域の場合
     &     idx(kele).eq.5 .or. idx(kele).eq.6)then
          do 200 j=1,4
            kp=ne(kele,j)
            c0(kp)=cini
200       continue
        endif
100   continue
      return
      end
c
c     熱移動領域
c
      subroutine initt(np,nele,ne,t0,tt,idx,ind,ine)
      implicit double precision(a-h,o-z)
      dimension ne(ine,4),t0(ind),tt(ind)
      dimension idx(0:ine)
      common/initch/cini,hini
c
      do 10 i=1,np
        t0(i)=0d0
        tt(i)=0d0
10    continue
c
      do 100 kele=1,nele
        if(idx(kele).eq.3 .or. idx(kele).eq.4 .or.  !index番号が熱移動領域の場合
     &     idx(kele).eq.5 .or. idx(kele).eq.7)then
          do 200 j=1,4
            kp=ne(kele,j)
            t0(kp)=hini
200       continue
        endif
100   continue
      return
      end