c&AR
c&ARCORR
        subroutine arcorr(args,nargs,vname)
c******************************************************************
c
c   Subroutine to process ARCORR.
c
c   rho=ARCORR(alpha,nord,rvar,ncorr,r0,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
        call ckint(args(2),nord)
        if(nord.lt.1.or.nord.gt.na) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(icki(args(4),4,ncorr,1,0).eq.1) go to 99
	if(ickse(nord+ncorr).eq.1) go to 99
c
c   calculate and store :
c
	call movxy(wk,array(nstart(np)),4*nord)
	n1=nord+1
	call arcv(wk,nord,rvar,ncorr,wk(n1),r0,ier)
        call ckaddi(args(6),ier,iref)
	if(ier.ne.0) then
		write(*,45) ier
  45    format(' process NONSTATIONARY, partial ',i3,' outside (-1,1)')
		imflag=1
		go to 99
	endif
	if(iref.eq.1) go to 99
        call ckaddr(args(5),r0,iref)
        if(iref.eq.1) go to 99
        if(ncorr.eq.0) go to 99
        do 50 i=1,ncorr
  50    wk(n1+i-1)=wk(n1+i-1)/r0
        lab='autoregressive correlations'
        call ckadda(vname,ncorr,lab,n1,iref)
c
c
  99    continue
        return
        end
c&ARCV
        subroutine arcv(alpha,nord,rvar,ncov,r,r0,ier)
c********************************************************************
c
c   Subroutine to calculate AR covariances r0,r(1),...,r(ncov)
c   from an AR(alpha,nord,rvar) process. If ncov=0 then only r0 is
c   calculated.
c
c   If partial(i) is greater than 1 in absolute value, then
c   arcv returns with the error flag ier set equal to i. If ier=0
c   on return, then there were no errors in arcv.
c
c
c*********************************************************************
c
        dimension alpha(nord),r(1)
c
c   Find R0:
c
	call arprt(alpha,nord,ier)
	if(ier.ne.0) go to 99
	r0=rvar
	do 10 i=1,nord
  10	r0=r0/(1.-alpha(i)*alpha(i))
	if(ncov.eq.0) go to 99
c
c   now work the other way:
c
	r(1)=-alpha(1)*r0
	if(ncov.eq.1) go to 99
	if(nord.eq.1) go to 50
        do 20 i=2,min0(nord,ncov)
        	c=alpha(i)
         	nl=i/2
         	do 15 j=1,nl
         		imj=i-j
         		temp=alpha(j)+c*alpha(imj)
         		alpha(imj)=alpha(imj)+c*alpha(j)
  15			alpha(j)=temp
		c=-alpha(i)*r0
		do 17 j=1,i-1
  17		c=c-alpha(j)*r(i-j)
	r(i)=c
  20    continue
	if(ncov.le.nord) go to 99
  50	do 60 i=nord+1,ncov
  		c=0.
		do 55 j=1,nord
  55		c=c-alpha(j)*r(i-j)
  60	r(i)=c
c
c
  99    continue
        return
        end
c&ARDT
        subroutine ardt(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to generate data from an AR
c
c   x=ARDT(alpha,nord,rvar,seed,n,ier,r0)
c
c   or
c
c   x=ARDT(alpha,nord,n,f) : x(i)=f(i)-sum(j=1,nord) of alpha(j)*x(i-j)
c   i=nord+1,...,n and x(i)=f(i),i=1,...,nord.
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,i10*1
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
        call ckint(args(2),nord)
        if(nord.gt.na.or.nord.le.0) then
		call error(args,2,2)
		go to 99
	endif
c
c   handle the difference equation case here:
c
        if(nargs.eq.4) then
		if(icki(args(3),3,n,2,nord).eq.1) go to 99
		if(ickl(args(4),nf,nfl).eq.1) go to 99
                if(nfl.lt.n) then
    	        	call error(args,3,2)
                	go to 99
	        endif
		if(ickse(n+nord).eq.1) go to 99
	        call movxy(wk,array(nstart(nf)),4*n)
	        n1=nstart(np)-1
	        do 702 i=nord+1,n
        	c=wk(i)
                	do 703 j=1,nord
 703 		        c=c-array(n1+j)*wk(i-j)
 702 	        wk(i)=c
	        lab='Values of Difference Equation'
	        go to 704
	endif
c
c   Non-difference equation case:
c
    	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(ickr(args(4),4,rseed,1,0.).eq.1) go to 99
        if(rseed.ne.0.) i4seed=rseed
	if(icki(args(5),5,n,1,nord).eq.1) go to 99
c
c   generate data :
c
        na=nstart(np)
	nws=n+nord*nord+nord
	if(ickse(nws).eq.1) go to 99
        call ardt1(array(na),nord,rvar,n,i4seed,nord,wk(n+1),
     1  wk(n+nord*nord+1),r0,wk,ier)
c
c   put stuff in right place :
c
        call ckaddi(args(6),ier,iref)
        if(iref.eq.1) go to 99
        if(ier.eq.0) go to 70
        write(*,61) ier
  61    format(' process NONSTATIONARY, partial ',i3,' outside (-1,1)')
        go to 99
  70    continue
        call ckaddr(args(7),r0,iref)
        if(iref.eq.1) go to 99
        lab='autoregressive process'
 704    call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&ARDT1
        subroutine ardt1(alpha,nord,rvar,n,dseed,ndim,a,d,r0,x,ier)
c********************************************************************
c
c   Subroutine to generate a realization x(1),...,x(n) from an
c   AR(nord,alpha,rvar) process using double precision integer seed
c   dseed. The (nordxnord) matrix a (which has row dimension ndim in
c   the calling program) contains the coefficients for orders
c   1,...,nord on output and r0 is the true variance of the process
c   on output. The vector d contains the true residual variances for
c   orders 0,...,nord-1 on output.
c
c   If partial(i) is greater than 1 in absolute value, then
c   ardt1 returns with the error flag ier set equal to i. If ier=0
c   on return, then there were no errors in ardt1.
c
c*********************************************************************
c
        dimension alpha(nord),a(ndim,ndim),d(nord),x(n)
        integer*4 dseed
c
c   generate n(0,1) noise for first nord errors and n(0,sqrt(rvar))
c   noise for the last n-nord errors. We'll fix the first nord later:
c
        call whtsim(n,dseed,x)
c
c   find lower order alpha's and rvar's (we use these to get first
c   nord x's) :
c
	call arall(alpha,nord,rvar,ndim,a,d,r0,ier)
	if(ier.ne.0) go to 99
c
c
c
c   generate data :
c
c   first do x(1),...,x(nord) :
c
	do 20 i=1,nord
  20	x(i)=sqrt(d(i))*x(i)
	if(nord.eq.1) go to 80
        do 50 i=2,nord
        im1=i-1
        c=x(i)
                do 45 j=1,im1
  45            c=c-a(im1,j)*x(i-j)
  50    x(i)=c
  80    continue
c
c   now do x(nord+1),...,x(n) :
c
	sqrv=sqrt(rvar)
	do 70 i=nord+1,n
  70	x(i)=sqrv*x(i)
        do 90 i=nord+1,n
        c=x(i)
                do 85 j=1,nord
  85            c=c-alpha(j)*x(i-j)
  90    x(i)=c
c
c   finish :
c
  99    continue
        return
        end
c&ARALL
	subroutine ARALL(alpha,nord,rvar,ndim,a,d,r0,ier)
c********************************************************************
c
c
c*********************************************************************
c
	dimension alpha(nord),a(ndim,ndim),d(nord)
c
c
	ier=0
        part=alpha(nord)
        if(abs(part).lt.1.0) go to 15
        ier=nord
        go to 99
  15    continue
        do 20 i=1,nord
  20    a(nord,i)=alpha(i)
        d(nord)=rvar/(1.-part*part)
        if(nord.eq.1) go to 40
        do 30 i=nord-1,1,-1
        ip1=i+1
        part=a(ip1,ip1)
        if(abs(part).lt.1.) go to 21
        ier=ip1
        go to 99
  21    continue
        part2=part*part
        opart2=1.-part2
                do 25 j=1,i
  25            a(i,j)=(a(ip1,j)-part*a(ip1,ip1-j))/opart2
        if(abs(a(i,i)).lt.1.) go to 26
        ier=i
        go to 99
  26    d(i)=d(ip1)/(1.-a(i,i)*a(i,i))
  30    continue
  40    continue
        r0=d(1)
c
c
  99	continue
	return
	end
c&ARPART
        subroutine arpart(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to process the ARPART command.
c
c   part=ARPART(alpha,nord,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
        call ckint(args(2),nord)
        if(nord.lt.1.or.nord.gt.na) then
		call error(args,2,2)
		go to 99
	endif
        ns=nstart(np)-1
        do 30 i=1,nord
  30    wk(i)=array(ns+i)
        call arprt(wk,nord,ier)
	call ckaddi(args(3),ier,iref)
	if(ier.ne.0) then
		write(*,35)
  35		format(' Process not stationary')
		imflag=1
	endif
	if(ier.ne.0.or.iref.eq.1) go to 99
        lab='partial autocorrelations'
	do 40 i=1,nord
  40	wk(i)=-wk(i)
        call ckadda(vname,nord,lab,1,iref)
c
c
  99    continue
        return
        end
c&ARPRT
        subroutine arprt(alpha,nord,ier)
c******************************************************************
c
c   Subroutine to find partials from coefficients.
c
c******************************************************************
c
        dimension alpha(nord)
	ier=0
        if(nord.eq.1) then
		if(abs(alpha(1)).ge.1.) ier=1
		go to 99
	endif
        do 10 i=1,nord
        mmi=nord-i
        mm=mmi+1
        c=alpha(mm)
	if(abs(c).ge.1.) then
		ier=i
		go to 99
	endif
        if(i.eq.nord) go to 10
        c1=1.-c*c
        nps=mmi/2
        if(mmi.eq.1) go to 30
                do 20 j=1,nps
                jj=mm-j
                alpha(jj)=(alpha(jj)-c*alpha(j))/c1
  20            alpha(j)=alpha(j)-c*alpha(jj)
        if(2*nps.eq.mmi) go to 10
  30    alpha(nps+1)=alpha(mm-nps-1)/(1.+c)
  10    continue
  99	continue
c
c
        return
        end
c&ARFILT
	subroutine arfilt(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to handle the command
c
c   y=ARFILT(alpha,nord,rvar,x,n,m,ier)
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,ickr,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),nord)
	if(nord.lt.1.or.nord.gt.na) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(ickl(args(4),npx,nx).eq.1) go to 99
	if(icki(args(5),5,n,2,0).eq.1) go to 99
	if(icki(args(6),6,m,2,0).eq.1) go to 99
	if(float(nx)/float(n).lt.float(m)) then
		call error(args,5,2)
		go to 99
	endif
	nws=2*n*m+nord*nord+2*nord
	if(ickse(nws).eq.1) go to 99
c
c
	call movxy(wk,array(nstart(npx)),4*n*m)
	n1=n*m+1
	n2=n1+n*m
	n3=n2+nord*nord
	call arflt1(array(nstart(npa)),nord,rvar,wk,wk(n1),n,m,
     1	wk(n2),wk(n3),ier)
	if(ier.eq.1) then
		write(*,20)
  20		format(' Process nonstationary in ARFILT')
 	endif
        call ckaddi(args(7),ier,iref)
	if(ier.eq.1.or.iref.eq.1) go to 99
	lab='Data filtered by an AR filter'
	call ckadda(vname,n*m,lab,n1,iref)
c
c
  99	continue
	return
	end
c&ARFLT1
	subroutine arflt1(alpha,nord,rvar,x,y,n,m,a,d,ier)
c****************************************************************
c
c
c****************************************************************
c
	dimension a(nord,nord),alpha(nord),x(n,m),y(n,m),d(nord)
	call arall(alpha,nord,rvar,nord,a,d,r0,ier)
	if(ier.ne.0) go to 99
	sqrv=sqrt(rvar)
	do 10 i=1,nord
  10	d(i)=sqrt(d(i))
	do 50 i=1,m
	y(1,i)=x(1,i)/d(1)
	do 40 j=2,nord
		jm1=j-1
		c=x(j,i)
		do 30 k=1,jm1
  30		c=c+a(jm1,k)*x(j-k,i)
  40	y(j,i)=c/d(j)
	do 20 j=nord+1,n
		c=x(j,i)
		do 15 k=1,nord
  15		c=c+a(nord,k)*x(j-k,i)
  20	y(j,i)=c/sqrv
  50	continue
c
c
  99	continue
	return
	end
