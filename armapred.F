c&ARMAPRED
c&MXPRED
	subroutine mxpred(args,nargs)
c********************************************************************
c
c   Subroutine to handle the command
c
c   ARMAPRED(x,n,alpha,beta,np,nq,rv,ntf,ntl,nhf,nhl,xpreds,ier[,stds])
c
c********************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	call ckint(args(2),n)
	if(n.lt.1.or.n.gt.nx) then
		call error(args,2,2)
		go to 99
	endif
c
c
	if(icki(args(5),5,np,1,0).eq.1) go to 99
	if(np.gt.0) then
		if(ickl(args(3),npa,na).eq.1) go to 99
		if(np.gt.na) then
			call error(args,5,2)
			go to 99
		endif
	endif
c
c
	if(icki(args(6),6,nq,1,0).eq.1) go to 99
	if(nq.gt.0) then
		if(ickl(args(4),npb,nb).eq.1) go to 99
		if(nq.gt.nb) then
			call error(args,6,2)
			go to 99
		endif
	endif
c
c
    	if(np.le.0) np=0
	if(nq.le.0) nq=0
	if(ickr(args(7),7,rvar,2,0.).eq.1) go to 99
c
c
	call ckint(args(8),ntf)
	call ckint(args(9),ntl)
	call ckint(args(10),nhf)
	call ckint(args(11),nhl)
	ipte=0
	if(ntf.lt.1.or.ntf.gt.ntl) ipte=8
	if(ntl.lt.1.or.ntf.gt.n) ipte=9
	if(nhf.lt.1.or.nhf.gt.nhl) ipte=10
	if(nhl.lt.1) ipte=11
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
	mxpqp1=max0(np,nq)+1
	nr=max0(np,nq+1)
	nr2=nr*nr
	npreds=(ntl-ntf+1)*(nhl-nhf+1)
	n1=1
	n2=n1+n
	n3=n2+nr
	n4=n3+nr
	n5=n4+mxpqp1
	n6=n5+nr2
	n7=n6+nr2
	n8=n7+nr2
	n9=n8+nr
	n10=n9+nr2
	n11=n10+mxpqp1
	n12=n11+nr2
	n13=n12+mxpqp1
	n14=n13+nr2
	n15=n14+npreds
	call movxy(wk,array(nstart(npx)),4*n)
	if(np.gt.0) call movxy(wk(n2),array(nstart(npa)),4*np)
	if(nq.gt.0) call movxy(wk(n3),array(nstart(npb)),4*nq)
	iptstd=0
	if(nargs.eq.14) iptstd=1
	call mxpd(wk(n1),wk(n2),wk(n3),np,nq,rvar,ntf,ntl,nhf,nhl,
     1	iptstd,nr,wk(n4),wk(n5),wk(n6),wk(n7),wk(n8),wk(n9),
     1	wk(n10),wk(n11),wk(n12),wk(n13),wk(n14),wk(n15))
	lab='Predictions'
	call ckadda(args(12),npreds,lab,n14,iref)
	if(iref.eq.1) go to 99
	if(nargs.eq.14) then
		lab='Prediction Standard Errors'
		call ckadda(args(14),npreds,lab,n15,iref)
	endif
c
c
  99	continue
	return
	end
c
c
c&MXPD
	subroutine mxpd(x,alpha,beta,np,nq,rvar,ntf,ntl,nhf,nhl,
     1	iptstd,nr,xt,st,w,a,xh,sh,xt1,st1,xt2,st2,xp,std)
c********************************************************************
c
c
c********************************************************************
c
	dimension x(1),alpha(1),beta(1),
     1	xt(1),st(nr,nr),w(nr,nr),a(nr,nr),xh(nr),
     1	sh(nr,nr),xt1(1),st1(nr,nr),xp(1),std(1),xt2(1),
     1	st2(nr,nr)
c
c   initial conditions:
c
	mxpq=max0(np,nq)
	mxpqp1=mxpq+1
	if(np.gt.0) call vneg(alpha,np)
    	if(nq.gt.0) call vneg(beta,nq)
c
c   Form Initial State Vector (xh) and Covariance Matrix (sh):
c
	call wilson(alpha,np,beta,nq,xt,mxpqp1,xt1,mxpqp1,
     1	xt2,mxpq,ier)
	if(ier.eq.1) go to 99
	if(np.gt.0) call vneg(alpha,np)
    	if(nq.gt.0) call vneg(beta,nq)
c
c
	do 50 i=1,nr
	xh(i)=0.
	do 50 j=1,i
    	sh(i,j)=xt(i-j+1)*rvar
  50	sh(j,i)=sh(i,j)
	if(nr.gt.1) then
		call mxma(alpha,beta,np,nq,nr,xt1)
		do 51 j=2,nr
		do 51 k=j,nr
		c=sh(j,k)
			do 52 l=0,j-2
			c1=1.
			if(l.gt.0) c1=xt1(l)
			c2=1.
			l1=l+k-j
			if(l1.gt.0) c2=xt1(l1)
  52			c=c-c1*c2*rvar
		sh(j,k)=c
  51		sh(k,j)=c
	endif
c
c   form a and w:
c
	xt1(1)=1.
	if(nr.gt.1) then
		do 55 i=2,nr
		c=0.0
		if(i-1.le.nq) c=beta(i-1)
		if(np.gt.0) then
			do 54 j=1,min0(np,i-1)
  54			c=c-alpha(j)*xt1(i-j)
		endif
  55		xt1(i)=c
	endif
	do 56 i=1,nr
	do 56 j=1,nr
  56	w(i,j)=rvar*xt1(i)*xt1(j)
	do 57 i=1,nr
	do 57 j=1,nr
  57	a(i,j)=0.0
	if(np.gt.0) then
		do 58 i=1,np
  58		a(nr,nr-i+1)=-alpha(i)
	endif
	if(nr.gt.1) then
		do 59 i=1,nr-1
  59		a(i,i+1)=1.
	endif
c
c   iterate:
c
	npr=0
	nr2=nr*nr
	nr24=4*nr2
	nr4=4*nr
c
c
	do 100 it=0,ntl
c
c  x tilda and sigma tilda:
c
	call mnvar(xh,a,sh,w,nr,nr,xt,st)
c
c   Forecasts:
c
	if(it.ge.ntf.and.it.le.ntl) then
        	call movxy(xt1,xt,nr4)
	        call movxy(st1,st,nr24)
 	        do 60 ih=1,nhl
		        if(ih.ge.nhf.and.ih.le.nhl) then
			        npr=npr+1
			        xp(npr)=xt1(1)
			        if(iptstd.eq.1) std(npr)=sqrt(st1(1,1))
		        endif
		        if(ih.eq.nhl) go to 60
c
c
		        call mnvar(xt1,a,st1,w,nr,nr,xt2,st2)
		        call movxy(xt1,xt2,nr4)
		        call movxy(st1,st2,nr24)
  60	        continue
	endif
	if(it.eq.ntl) go to 100
c
c   Update x hat and sigma hat:
c
	st11=st(1,1)
	do 65 i=1,nr
  65	xt1(i)=st(1,i)/st11
c
c   xhat:
c
	z=x(it+1)-xt(1)
	do 70 i=1,nr
  70	xh(i)=xt(i)+xt1(i)*z
c
c  sigma hat:
c
c	do 66 i=1,nr
c	do 66 j=1,nr
c  66	st1(i,j)=0.
c	do 67 i=1,nr
c    	st1(i,i)=1.
c  67	st1(i,1)=st1(i,1)-xt1(i)
c	do 68 i=1,nr
c	do 68 j=1,nr
c	c=0.
c		do 69 k=1,nr
c  69		c=c+st1(i,k)*st(k,j)
c  68	sh(i,j)=c
        do 68 i=1,nr
        do 68 j=1,nr
  68    sh(i,j)=st(i,j)-xt1(i)*st(1,j)
c
c
 100	continue
c
c
c
  99	continue
	return
	end
c&MNVAR
	subroutine mnvar(x,a,sig,w,ndim,n,y,b)
c*******************************************************************
c
c   Subroutine to calculate
c
c   y=a*x     and     b=a*sig*a' + w
c
c   for nx1 x and nxn a and sig  in Kalman Filter for ARMAPRED
c
c*******************************************************************
c
	dimension x(n),a(ndim,ndim),sig(ndim,ndim),b(ndim,ndim),y(n)
	dimension w(ndim,ndim)
c
c
	do 20 i=1,n
		if(i.lt.n) then
			c=x(i+1)
			go to 20
		endif
c
		c=0.0
		do 10 j=1,n
  10		c=c+a(i,j)*x(j)
  20	y(i)=c
c
c
	if(n.eq.1) go to 40
	do 30 i=1,n-1
	do 30 j=1,i
	c=w(i,j)+sig(i+1,j+1)
	b(i,j)=c
  30	b(j,i)=c
c
c   i=n:
c
  40	continue
	do 70 j=1,n
		c=w(n,j)
		do 60 k=1,n
	 	if(j.lt.n) then
			c1=sig(k,j+1)
			go to 60
	 	endif
			c1=0.
			do 50 l=1,n
  50			c1=c1+sig(k,l)*a(j,l)
  60		c=c+a(n,k)*c1
  	b(n,j)=c
  70	b(j,n)=c
c
c
	return
	end
c&VNEG
	subroutine vneg(x,n)
c********************************************************************
c
c   Subroutine to negate a vector
c
c*********************************************************************
c
	dimension x(n)
	do 10 i=1,n
  10	x(i)=-x(i)
	return
	end
c&MXMA
	subroutine mxma(alpha,beta,np,nq,n,d)
c********************************************************************
c
c   Subroutine to find the first n coefficients of the MA(ì)
c   representation of an ARMA(alpha,beta,np,nq) process.
c
c********************************************************************
c
	dimension alpha(1),beta(1),d(n)
c
	do 10 i=1,n
  10	d(i)=0.0
	if(nq.gt.0) then
		do 20 i=1,nq
  20		d(i)=beta(i)
	endif
c
c
	if(np.gt.0) then
		do 40 i=1,n
		c=d(i)
		if(np.eq.0) go to 40
			do 30 j=1,min0(i,np)
			c1=1.
			if(i-j.gt.0) c1=d(i-j)
  30			c=c-alpha(j)*c1
  40		d(i)=c
	endif
c
c
	return
	end
