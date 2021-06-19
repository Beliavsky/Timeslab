c&FFTS
c&FFTS
        subroutine ffts(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   FFT(xr,xi,nfreqs,nout,isign,zr,zi)  or
c
c   FFT(x,nfreqs,nout,isign,zr,zi)
c
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,i10*1
        integer*2 ickl,icki,ickse
c
c
        if(ickl(args(1),npxr,nxr).eq.1) go to 99
        nxi=nxr
        if(nargs.eq.7) then
                if(ickl(args(2),npxi,nxi).eq.1) go to 99
        endif
        m1=3
        if(nargs.eq.6) m1=2
        if(icki(args(m1),m1,nfreqs,1,1).eq.1) go to 99
        if(nfreqs.gt.min0(nxr,nxi)) then
                call error(args,m1,2)
                go to 99
                endif
        m1=m1+1
        if(icki(args(m1),m1,nout,0,0).eq.1) go to 99
        if(nout.lt.1.or.nout.gt.nfreqs) then
                call error(args,m1,2)
                go to 99
                endif
        m1=m1+1
        if(icki(args(m1),m1,isign,0,0).eq.1) go to 99
        if(isign.ne.1.and.isign.ne.-1) then
                call error(args,m1,2)
                go to 99
                endif
c
c
        call movxy(wk,array(nstart(npxr)),4*nfreqs)
        if(nargs.eq.7) call movxy(wk(nfreqs+1),array(nstart(npxi)),
     1  4*nfreqs)
        if(nargs.eq.6) call movct(wk(nfreqs+1),4*nfreqs,char(0))
c
c
        if(ickse(2*nfreqs).eq.1) go to 99
        call fft(wk,wk(nfreqs+1),nfreqs,nfreqs,nfreqs,isign)
c
c
        if(imflag.ne.1) go to 71
        imflag=0
        write(*,5)
   5    format(' Using the Brute Force DFT: It might take a long time')
        ns1=nstart(npxr)
        iopt=1
        ns2=ns1
        if(nargs.eq.7) then
                iopt=2
                ns2=nstart(npxi)
        endif
        n1=nfreqs+1
        n2=n1+nfreqs
        n3=n2+nout
c
c
        if(ickse(2*nfreqs+2*nout).eq.1) go to 99
        call fft1(array(ns1),array(ns2),nfreqs,iopt,isign,nout,
     1  wk,wk(n1),wk(n2),wk(n3))
c
c
        lab='Real Part of FFT'
        m1=m1+1
        call ckadda(args(m1),nout,lab,n2,iref)
        if(iref.eq.1) go to 99
        lab='Imaginary Part of FFT'
        m1=m1+1
        call ckadda(args(m1),nout,lab,n3,iref)
        go to 99
c
c
  71    continue
        lab='Real Part of FFT'
        m1=m1+1
        call ckadda(args(m1),nout,lab,1,iref)
        if(iref.eq.1) go to 99
        lab='Imaginary Part of FFT'
        m1=m1+1
        call ckadda(args(m1),nout,lab,nfreqs+1,iref)
c
c
  99    continue
        return
        end
c&FFT1
        subroutine fft1(a,b,n,iopt,isign,m,cc,ss,zr,zi)
c******************************************************************
c
c   Subroutine to do brute force n-point DFT. a and b are
c   real and imaginary parts of input, zr and zi are first m points
c   of output. If iopt=1 (2), then b is assumed (not assumed)
c   to be zeros. isign is 1 or -1 for DFT or IDFT. cc and ss are
c   work vectors of length n.
c
c*******************************************************************
c
        dimension a(1),b(1),cc(1),ss(1),zr(1),zi(1)
c
c
        on=float(n)
        oi=float(isign)
        twopi=8.*atan(1.0)
        c1=twopi/on
        n1=(n/2)+1
        cc(1)=1.0
        ss(1)=0.0
c
c
        jj=n+1
        do 10 j=2,n1
        jj=jj-1
        arg=float(j-1)*c1
        cc(j)=cos(arg)
        ss(j)=oi*sin(arg)
        cc(jj)=cc(j)
  10    ss(jj)=-ss(j)
        if(2*(n/2).eq.n) ss(n1)=-ss(n1)
c
c
        mm=m
        iopt1=0
        if(iopt.eq.1.and.m.gt.n1) then
                mm=n1
                iopt1=1
        endif
        iptev=0
        if(2*(n/2).eq.n) iptev=1
c
c
        if(iptev.eq.1) ajj=a(n1)
        if(iptev.eq.1.and.iopt.eq.2) bjj=b(n1)
        jj=n+1
        do 15 j=2,n1
        jj=jj-1
        ajp=a(j)+a(jj)
        ajm=a(j)-a(jj)
        a(jj)=ajm
  15    a(j)=ajp
        if(iopt.eq.1) go to 17
        jj=n+1
        do 16 j=2,n1
        jj=jj-1
        bjp=b(j)+b(jj)
        bjm=b(j)-b(jj)
        b(jj)=bjm
  16    b(j)=bjp
  17    continue
c
c
        cc1=0.0
        if(iopt.eq.2) cc1=b(1)
c
c
        do 40 i=1,mm
        c1=a(1)
        c2=cc1
        ind=1
        im1=i-1
        jj=n+1
c
        if(iopt.eq.1) then
                do 30 j=2,n1
                jj=jj-1
                ind=ind+im1
                if(ind.gt.n) ind=ind-n
                cci=cc(ind)
                ssi=ss(ind)
                c1=c1+a(j)*cci
                c2=c2+a(jj)*ssi
  30            continue
                if(iptev.eq.1) c1=c1-ajj*cci
        else
                do 31 j=2,n1
                jj=jj-1
                ind=ind+im1
                if(ind.gt.n) ind=ind-n
                cci=cc(ind)
                ssi=ss(ind)
                c1=c1+a(j)*cci-b(jj)*ssi
                c2=c2+a(jj)*ssi+b(j)*cci
  31            continue
                if(iptev.eq.1) then
                        c1=c1-ajj*cci
                        c2=c2-bjj*cci
                endif
        endif
        zr(i)=c1
        zi(i)=c2
  40    continue
c
c   Unscramble a's and b's:
c
        jj=n+1
        do 45 j=2,n1
        jj=jj-1
        ajp=(a(j)+a(jj))/2.
        ajm=(a(j)-a(jj))/2.
        a(jj)=ajm
  45    a(j)=ajp
        if(iptev.eq.1) a(n1)=ajj
        if(iopt.eq.1) go to 47
        jj=n+1
        do 46 j=2,n1
        jj=jj-1
        bjp=(b(j)+b(jj))/2.
        bjm=(b(j)-b(jj))/2.
        b(jj)=bjm
  46    b(j)=bjp
        if(iptev.eq.1) b(n1)=bjj
  47    continue
c
c

        if(iopt1.eq.1) then
                jj=n-n1+2
                do 50 j=n1+1,m
                jj=jj-1
                zr(j)=zr(jj)
  50            zi(j)=-zi(jj)
        endif
c
c
        return
        end
c&REPL
        subroutine repl(args,nargs,vname)
c******************************************************************
c
c   Subroutine to handle the command
c
c   y=REPLACE(x,xi,rep,nr)
c
c  or
c
c   y=REPLACE(x,n,iopt,a[,b,c])   iopt=1 (below), 2(above), 3(both)
c                                      4 (in between)
c                                      5 (b if le a, c if >a)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
        integer*2 ickl,icki,ickr
c
c
        if(ickl(args(1),npx,nx).eq.1) go to 99
        if(nargs.eq.4) then
                call ckint(args(2),n)
                if(n.eq.-999) go to 5
        endif
c
c   handle second case:
c
        if(icki(args(2),2,n,5,nx).eq.1) go to 99
        if(icki(args(3),3,iopt,1,1).eq.1) go to 99
        if(iopt.gt.5) then
                call error(args,3,2)
                go to 99
        endif
        if(iopt.eq.3.and.nargs.ne.5) then
                write(*,901)
 901            format(' Wrong number of arguments in REPLACE')
                imflag=1
                go to 99
        endif
        if(iopt.ge.4.and.nargs.ne.6) then
                write(*,901)
                imflag=1
                go to 99
        endif
        if(ickr(args(4),4,a,0,0.).eq.1) go to 99
        if(iopt.gt.2) then
                if(ickr(args(5),5,b,0,0.).eq.1) go to 99
        endif
        if(iopt.ge.4) then
                if(ickr(args(6),6,c,0,0.).eq.1) go to 99
        endif
        call movxy(wk,array(nstart(npx)),4*n)
        if(iopt.eq.1) then
                do 1 i=1,n
   1            if(wk(i).lt.a) wk(i)=a
                go to 7
        endif
        if(iopt.eq.2) then
                do 2 i=1,n
   2            if(wk(i).gt.a) wk(i)=a
                go to 7
        endif
        if(iopt.eq.3) then
                do 3 i=1,n
                if(wk(i).lt.a) wk(i)=a
   3            if(wk(i).gt.b) wk(i)=b
                go to 7
        endif
        if(iopt.eq.4) then
                do 4 i=1,n
   4            if(a.le.wk(i).and.wk(i).le.b) wk(i)=c
                go to 7
        endif
        if(iopt.eq.5) then
                do 6 i=1,n
                cc=wk(i)
                if(cc.le.a) wk(i)=b
   6            if(cc.gt.a) wk(i)=c
                go to 7
        endif
   7    lab='Truncated Array'
        call ckadda(vname,n,lab,1,iref)
        go to 99
c
c   handle first case:
c
   5    if(ickl(args(2),npxi,nxi).eq.1) go to 99
        if(nxi.gt.nx) then
                call error(args,2,2)
                go to 99
        endif
        if(icki(args(4),4,nr,1,1).eq.1) go to 99
        if(nr.gt.min0(nx,nxi)) then
                call error(args,4,2)
                go to 99
        endif
        rep=0.
        call cklist(args(3),xlist,numarr,nprep)
        if(nprep.ne.-999) then
                iptr=1
                nrep=nend(nprep)-nstart(nprep)+1
                if(nr.gt.nrep) then
                        call error(args,4,2)
                        go to 99
                endif
                go to 10
        endif
        iptr=0
        if(ickr(args(3),3,rep,0,0.).eq.1) go to 99
c
c
  10    nsxi=nstart(npxi)-1
        if(iptr.eq.1) nxr=nstart(nprep)-1
        call movxy(wk,array(nstart(npx)),4*nx)
        do 20 i=1,nr
        ii=array(nsxi+i)
        if(ii.lt.1.or.ii.gt.nx) then
                write(*,15)
  15            format(' Index outside range of array')
                imflag=1
                go to 99
        endif
        wk(ii)=rep
        if(iptr.eq.1) wk(ii)=array(nxr+i)
  20    continue
c
c
        lab='Array from REPLACE'
        call ckadda(vname,nx,lab,1,iref)
c
c
  99    continue
        return
        end
c&PARCOR
        subroutine parcor(args,nargs,vname)
c******************************************************************
c
c   Subroutine to handle the command
c
c   part=PARCORR(x,n,m) or
c
c   part=PARCORR(x,n,m,rvar)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
        double precision inprod
        integer*2 ickl,ickse
c
c
        if(ickl(args(1),npx,nx).eq.1) go to 99
        call ckint(args(2),n)
        if(n.lt.1.or.n.gt.nx) then
                call error(args,2,2)
                go to 99
        endif
        call ckint(args(3),m)
        if(m.lt.1.or.m.gt.n-1) then
                call error(args,3,2)
                go to 99
        endif
c
c
        npm=n+m
        nws=2*(n+m)+2*m
        if(ickse(nws).eq.1) go to 99
        n4=4*n
        npm2=2*npm
        npm24=4*npm2
        call movct(wk,npm24,char(0))
        call movxy(wk,array(nstart(npx)),n4)
        call movxy(wk(npm+2),array(nstart(npx)),n4)
        r0=inprod(wk(1),wk(1),n)/float(n)
        rvar=r0
        do 120 i=1,m
c
c   get partial :
c
        top=2.*inprod(wk(1),wk(npm+1),npm)
        bot=inprod(wk(1),wk(1),npm)+inprod(wk(npm+1),wk(npm+1),npm)
        part=top/bot
        if(part.le.-1.) part=-.9999
        if(part.ge.1.) part=.9999
        wk(npm2+i)=part
        rvar=rvar*(1.-part*part)
        wk(npm2+m+i)=rvar
c
c   shift :
c
        do 160 j=1,npm-1
        i1=npm-j+1
        i2=i1-1
        ni1=npm+i1
        ni2=npm+i2
        temp=wk(ni1)
        wk(ni1)=wk(ni2)-part*wk(i2)
 160    wk(i1)=wk(i1)-part*temp
        wk(1)=wk(1)-part*wk(npm+1)
        wk(npm+1)=0.0
c
c   end recursion :
c
 120    continue
c
c
        lab='Partial Autocorrelations'
        call ckadda(vname,m,lab,npm2+1,iref)
        if(iref.eq.1.or.nargs.eq.3) go to 99
        lab='Residual Variances'
        call ckadda(args(4),m,lab,npm2+m+1,iref)
c
c
  99    continue
        return
        end
c&DNSITY
        subroutine dnsity(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   DENSITY(x,nobs,rbins,npts,kern,y,fy)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        integer*2 ickl,icki,ickse
        pi=4.*atan(1.0)
c
c
        if(ickl(args(1),npx,nx).eq.1) go to 99
        call ckint(args(2),nobs)
        if(nobs.lt.2.or.nobs.gt.nx) then
                call error(args,2,2)
                go to 99
        endif
        call ckreal(args(3),rbins)
        if(rbins.lt.0..or.rbins.gt.float(nobs)) then
                call error(args,3,2)
                go to 99
        endif
        if(icki(args(4),4,npts,0,0).eq.1) go to 99
        if(iabs(npts).lt.2) then
                call error(args,4,2)
                go to 99
        endif
        iptdet=0
        if(npts.lt.0) then
                iptdet=1
                npts=-npts
                if(ickl(args(6),npy,ny).eq.1) go to 99
                if(npts.gt.ny) then
                        call error(args,6,1)
                        go to 99
                endif
                call movxy(wk,array(nstart(npy)),4*npts)
                do 1 i=1,npts-1
                if(wk(i).gt.wk(i+1)) then
                        call error(args,6,2)
                        go to 99
                endif
   1            continue
        endif
        call ckint(args(5),iopt)
        if(iopt.lt.1.or.iopt.gt.7) then
                call error(args,5,2)
                go to 99
        endif
        if(ickse(npts).eq.1) go to 99
        if(ickse(2*npts+nobs).eq.1) go to 99
c
c
        call movxy(wk,array(nstart(npx)),4*nobs)
        call qsort(nobs,wk)
        xmin=wk(1)
        xmax=wk(nobs)
        rx=xmax-xmin
        if(rx.lt.1.e-20) then
                write(*,10)
  10            format(' Range of Data less than 1.e-20')
                imflag=1
                go to 99
        endif
        if(rbins.eq.0.) then
                call xbars2(wk,nobs,xbar,s2)
                ot=1./3.
                if(iopt.gt.1) ot=1./5.
                rbins=(rx*(float(nobs)**ot)/(3.5*sqrt(s2)))+1.
        endif
        bw=(rx/rbins)/2.
        onh=float(nobs)*bw
        yinc=rx/float(npts-1)
        nobsp=nobs+npts
c
c
        j1=1
        do 100 i=1,npts
        if(iptdet.eq.1) y=array(nstart(npy)+i-1)
        if(iptdet.eq.0) y=xmin+float(i-1)*yinc
        yl=y-bw
        yu=y+bw
                do 20 j=j1,nobs
                if(wk(j).ge.yl.and.wk(j).le.yu) go to 30
  20            continue
                nin=0
                go to 60
  30            nin=1
                j1=j
                if(j.eq.nobs) go to 60
                do 40 k=j+1,nobs
                if(wk(k).gt.yu) go to 50
  40            nin=nin+1
  50    continue
  60    call smooth(wk(j1),nin,y,bw,onh,iopt,pi,fy)
        wk(nobs+i)=y
        wk(nobsp+i)=fy
 100    continue
c
c
        lab='y'
        call ckadda(args(6),npts,lab,nobs+1,iref)
        if(iref.eq.1) go to 99
        lab='f(y)'
        call ckadda(args(7),npts,lab,nobsp+1,iref)
c
c
  99    continue
        return
        end
c&SMOOTH
        subroutine smooth(x,nin,y,bw,onh,iopt,pi,fy)
c*******************************************************************
c
c   Subroutine to find fy using a kernel function in dnsity
c
c********************************************************************
c
        dimension x(1)
c
        fy=0.
        if(nin.eq.0) go to 99
        if(iopt.eq.1) then
                fy=.5*float(nin)/onh
                go to 99
        endif
        if(iopt.eq.2) then
                do 10 i=1,nin
                        u=(y-x(i))/bw
  10            fy=fy+(1.+cos(pi*u))/2.
                fy=fy/onh
                go to 99
        endif
        if(iopt.eq.3) then
                do 20 i=1,nin
                        u=(y-x(i))/bw
  20            fy=fy+.75*(1.-u*u)
                fy=fy/onh
                go to 99
        endif
        if(iopt.eq.4) then
                do 30 i=1,nin
                        u=(y-x(i))/bw
  30            fy=fy+.9375*(1.-u*u)**2
                fy=fy/onh
                go to 99
        endif
        if(iopt.eq.5) then
                do 40 i=1,nin
                        u=(y-x(i))/bw
  40            fy=fy+(1.-abs(u))
                fy=fy/onh
                go to 99
        endif
        if(iopt.eq.6) then
                do 50 i=1,nin
                        u=(y-x(i))/bw
  50            fy=fy+1.012339*exp(-pi*u*u)
                fy=fy/onh
                go to 99
        endif
        if(iopt.eq.7) then
                do 60 i=1,nin
                        u=(y-x(i))/bw
                        au=abs(u)
                        if(au.le..5) fy=fy+1.333-8.*(u*u-au**3)
                        if(au.gt..5) fy=fy+2.6667*(1.-au)**3
  60            continue
                fy=fy/onh
                go to 99
        endif
c
c
  99    continue
        return
        end
