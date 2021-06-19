c&TOEPL
c&TOEPL
	subroutine toepl(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   A=TOEPL(r,r0,n)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,ickr,icki
c
c
	if(icki(args(3),3,n,1,1).eq.1) go to 99
	if(n.gt.100) then
		call error(args,3,2)
		go to 99
	endif
	if(ickr(args(2),2,r0,0,0.).eq.1) go to 99
	wk(1)=r0
	if(n.eq.1) go to 20
	if(ickl(args(1),npr,nr).eq.1) go to 99
	if(n-1.gt.nr) then
		call error(args,3,2)
		go to 99
		endif
	do 5 i=1,n
   5	wk((i-1)*n+i)=r0
	n1=nstart(npr)-1
	do 10 i=1,n
	do 10 j=1,n
	if(i.eq.j) go to 10
	wk((i-1)*n+j)=array(n1+iabs(i-j))
  10	continue
  20	continue
	lab='Toeplitz matrix'
	call ckadda(vname,n*n,lab,1,iref)
c
c
  99	continue
	return
	end
c&FILT
	subroutine filt(args,nargs,vname)
c******************************************************************
c
c   Subroutine to handle the command
c
c   y=FILT(x,beta,beta0,n,m)
c
c   y(i)=beta0*x(i+m)+Sum(j=1,m) of beta(j)*x(i+m-j), i=1,...,n-m
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(ickl(args(2),npb,nb).eq.1) go to 99
        if(ickr(args(3),3,beta0,0,0.).eq.1) go to 99
	if(icki(args(4),4,n,0,0).eq.1) go to 99
        if(icki(args(5),5,m,0,0).eq.1) go to 99
	ipte=0
	if(n.lt.1.or.n.gt.nx) ipte=4
	if(m.lt.1.or.m.gt.nb.or.m.ge.n) ipte=5
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
	if(ickse(2*n).eq.1) go to 99
c
c
	nmm=n-m
	npm=n+m
	call movxy(wk,array(nstart(npx)),4*n)
	call movxy(wk(n+1),array(nstart(npb)),4*m)
	do 10 i=1,nmm
		ipm=i+m
		c=beta0*wk(ipm)
		do 5 j=1,m
   5		c=c+wk(n+j)*wk(ipm-j)
  10	wk(npm+i)=c
c
c
	lab='Filtered Time Series'
	call ckadda(vname,nmm,lab,npm+1,iref)
c
c
  99	continue
	return
	end
c&DOUBLE
	subroutine double(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   DOUBLE(fxr,fxi,nfreqs)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickse
c
c
	if(icki(args(3),3,nfreqs,0,0).eq.1) go to 99
	if(ickl(args(1),npr,nr).eq.1) go to 99
	if(ickl(args(2),npi,ni).eq.1) go to 99
	n1=(nfreqs/2)+1
	if(n1.ne.nr.or.n1.ne.ni) then
		call error(args,3,2)
		go to 99
	endif
	if(ickse(2*nfreqs).eq.1) go to 99
	call movxy(wk,array(nstart(npr)),4*n1)
	call movxy(wk(nfreqs+1),array(nstart(npi)),4*n1)
	do 10 i=nfreqs,n1+1,-1
	wk(i)=wk(nfreqs-i+2)
  10	wk(nfreqs+i)=-wk(2*nfreqs-i+2)
c
c
	lab=label(npr)
	call ckadda(args(1),nfreqs,lab,1,iref)
	if(iref.eq.1) go to 99
	lab=label(npi)
	call ckadda(args(2),nfreqs,lab,nfreqs+1,iref)
c
c
  99	continue
	return
	end
c&REVRSE
	subroutine revrse(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to handle the command
c
c   y=reverse(x,n)
c
c*****************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(icki(args(2),2,n,0,0).eq.1) go to 99
	if(n.lt.1.or.n.gt.nx) then
		call error(args,2,2)
		go to 99
	endif
c
c
	call movxy(wk,array(nstart(npx)),4*n)
	nb2=n/2
	do 10 i=1,nb2
	n1=n-i+1
	c=wk(i)
	wk(i)=wk(n1)
  10    wk(n1)=c
	lab='array in reverse order'
	call ckadda(vname,n,lab,1,iref)
c
c
  99	continue
	return
	end
c&MCD
	subroutine mcd(args,nargs)
c******************************************************************
c
c   Subroutine to handle the command
c
c   MCHOL(a,n,l,d,ier)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),n)
	if(n.lt.1.or.n*n.gt.na) then
		call error(args,2,2)
		go to 99
	endif
c
c
	n1=1
	n2=n1+n*n
	n3=n2+n
	n4=n3+n*n
	if(ickse(n4).eq.1) go to 99
c
c
	call mchol(array(nstart(npa)),n,n,wk(n1),wk(n2),ier)
	call ckaddi(args(5),ier,iref)
	if(ier.eq.1) then
		write(*,10)
  10		format(' Matrix not positive definite in MCHOL')
		imflag=1
		go to 99
	endif
	if(iref.eq.1) go to 99
	lab='Lower Triangular Cholesky Factor'
	if(n.eq.1) go to 30
	do 20 i=2,n
	im1=i-1
  20	call movct(wk(im1*n+1),4*im1,char(0))
  30	continue
	call ckadda(args(3),n*n,lab,n1,iref)
	if(iref.eq.1) go to 99
	lab='Diagonal Cholesky Factor'
	call movct(wk(n3),4*n*n,char(0))
	do 40 i=1,n
	im1=i-1
  40	wk(n3+im1*n+im1)=wk(n2+im1)
	call ckadda(args(4),n*n,lab,n3,iref)
c
c
  99	continue
	return
	end
c&SWEEP
	subroutine sweep(args,nargs,vname)
c******************************************************************
c
c   Subroutine to handle the command
c
c   B = SWEEP(A,n,k1,k2,ier)  or  B = SWEEP(A,n,x,m,ier)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),n)
	if(n.lt.1.or.n*n.gt.na) then
		call error(args,2,2)
		go to 99
	endif
	call ckint(args(3),k1)
	if(k1.eq.-999) go to 20
	call ckint(args(4),k2)
	if(k1.lt.1.or.k1.gt.k2.or.k2.gt.n) then
		call error(args,3,2)
	 	go to 99
	endif
	do 10 i=k1,k2
  10	wk(i-k1+1)=i
	m=k2-k1+1
	go to 50
c
c
  20	if(ickl(args(3),npx,nx).eq.1) go to 99
	call ckint(args(4),m)
	if(m.lt.1.or.m.gt.nx) then
		call error(args,4,2)
		go to 99
	endif
	ns1=nstart(npx)-1
	do 40 i=1,m
	ii=array(ns1+i)
		if(ii.lt.1.or.ii.gt.n) then
			write(*,30)
  30			format(' Illegal Diagonal Number')
			imflag=1
			go to 99
		endif
  40	wk(i)=ii
c
c
  50	continue
	if(ickse(m+n*n).eq.1) go to 99
	call movxy(wk(m+1),array(nstart(npa)),4*n*n)
	do 70 i=1,m
c
c   find and swap first remaining nonzero diagonal:
c
		do 65 j=i,m
		ii=wk(j)
		xi=wk(m+(ii-1)*n+ii)
		if(abs(xi).gt.1.e-25) go to 66
  65		continue
		ier=1
		go to 56
  66	continue
	c=wk(j)
	wk(j)=wk(i)
	wk(i)=c
	ii=wk(i)
  55	call swpk12(wk(m+1),n,n,ii,ii,ier)
  56	if(ier.eq.1) then
		write(*,60)
  60		format(' All remaining diagonals zero in SWEEP')
		imflag=1
		go to 80
	endif
  70	continue
c
c
  80	call ckaddi(args(5),ier,iref)
	if(iref.eq.1.or.ier.eq.1) go to 99
	lab='Swept Matrix'
	call ckadda(vname,n*n,lab,m+1,iref)
c
c
  99	continue
	return
	end
c&GS
	subroutine gs(args,nargs)
c*****************************************************************
c
c   Subroutine to handle the command
c
c   GS(x,n,m,q,r,ier)
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickse
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(icki(args(2),2,n,1,1).eq.1) go to 99
	if(icki(args(3),3,m,1,1).eq.1) go to 99
	if(n*m.gt.nx.or.n.lt.m) then
		call error(args,2,2)
		go to 99
	endif
c
c
	n1=1
	n2=n1+n*m
	n3=n2+m*m
	if(ickse(n3).eq.1) go to 99
	call grmsmt(array(nstart(npx)),n,m,wk(n1),wk(n2),ier)
	call ckaddi(args(6),ier,iref)
	if(iref.eq.1.or.ier.eq.1) go to 99
	lab='Orthogonal Matrix in GS'
	call ckadda(args(4),n*m,lab,n1,iref)
	if(iref.eq.1) go to 99
	lab='Upper Triangular Matrix in GS'
	call ckadda(args(5),m*m,lab,n2,iref)
c
c
  99	continue
	return
	end
c&GRMSMT
	subroutine grmsmt(x,n,m,q,r,ier)
c*******************************************************************
c
c   Subroutine to calculate the GS decomposition X=QR of the nxm
c   matrix X.
c
c********************************************************************
c
	dimension x(n,m),q(n,m),r(m,m)
	double precision inprod
c
c
	ier=0
	call movxy(q,x,4*n*m)
	call movct(r,4*m*m,char(0))
	do 10 i=1,m
  10	r(i,i)=1.
	if(m.eq.1) go to 99
	do 50 i=1,m
	di=inprod(q(1,i),q(1,i),n)
	if(di.lt.1.e-25) then
		write(*,15)
  15		format(' X singular in GS')
		ier=1
		go to 99
	endif
	if(i.eq.m) go to 50
		do 20 j=i+1,m
		r(i,j)=inprod(q(1,j),q(1,i),n)/di
  20		call accum(q(1,j),q(1,i),n,r(i,j),q(1,j))
  50	continue
c
c
  99	continue
	return
	end
c&ACCUM
	subroutine accum(y,x,n,a,z)
c******************************************************************
c
c   Subroutine to do the accumulation
c
c   z=y-ax
c
c*******************************************************************
c
	dimension y(n),x(n),z(n)
c
c
	do 10 i=1,n
  10	z(i)=y(i)-a*x(i)
c
c
	return
	end
c&MABS
	subroutine mabs(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   y=ABS(x,n)  or  y=abs(x)
c
c********************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl
c
c
	if(nargs.eq.1) then
		call ckreal(args(1),x)
		if(x.eq.999.) then
                        if(ickl(args(1),npx,nx).eq.1) go to 99
                        n=nx  
			go to 10
		endif
		y=abs(x)
		call ckaddr(vname,y,iref)
		go to 99
	endif
c
c
	if(nargs.eq.2) then
		if(ickl(args(1),npx,nx).eq.1) go to 99
		call ckint(args(2),n)
		if(n.lt.1.or.n.gt.nx) then
			call error(args,2,2)
			go to 99
		endif
        endif
c
c
  10	ns1=nstart(npx)-1
	do 20 i=1,n
  20	wk(i)=abs(array(ns1+i))
	lab='Absolute Value of an Array'
	call ckadda(vname,n,lab,1,iref)
c
c
  99	continue
	return
	end
