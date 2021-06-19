c&T3AF
c&BART
        subroutine bart(args,nargs,vname)
c****************************************************************
c
c   Subroutine to do Bartletts test for white noise.
c
c   pval=BARTTEST(cper,nfreqs,n,B)  or pval=BARTTEST(xmax)
c
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr
c
c
	if(nargs.eq.1) then
		if(ickr(args(1),1,xmax,2,0.0).eq.1) go to 99
		call bartf(xmax,pval)
		call ckaddr(vname,pval,iref)
		go to 99
	endif
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
        call ckint(args(2),nfreqs)
        n1=(nfreqs/2)+1
        if(n1.ne.na) then
	        call error(args,2,2)
        	go to 99
	endif
	if(icki(args(3),3,n,2,0).eq.1) go to 99
c
c
        nn=nstart(np)-1
        do 40 i=1,n1
  40    wk(i)=abs(array(nn+i)-(float(i)/float(n1)))
        call max(wk,n1,xmax,ind)
        fac=sqrt(2./float(n))
	xmax=xmax/fac
	call bartf(xmax,pval)
c
c
	call ckaddr(args(4),xmax,iref)
	if(iref.eq.1) go to 99
        call ckaddr(vname,pval,iref)
c
c
  99    continue
        return
        end
c&BARTF
	subroutine bartf(xmax,pval)
c***************************************************************
c
c****************************************************************
c
	pval=1.
	do 50 i=1,50
	c1=(-1.)**i
	oi=i
	c2=exp(-2.*xmax*xmax*oi*oi)
	if(c2.lt.1.e-20) go to 55
	c=c1*c2
  50	pval=pval+2.*c
  55	continue
	pval=1.-pval
c
c
	return
	end
c&CORR
        subroutine corr(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to find periodogram and correlogram.
c
c   corr(x,nobs,ncov,nfreqs,iopt,r0,per)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickse
	double precision inprod
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
        call ckint(args(2),nobs)
	call ckint(args(3),ncov)
	call ckint(args(4),nfreqs)
	call ckint(args(5),iopt)
	ipte=0
	if(nobs.le.0.or.nobs.gt.nx) ipte=2
	if(ncov.lt.0.or.ncov.ge.nobs) ipte=3
	if(nfreqs.lt.0) ipte=4
	if(nfreqs.ne.0.and.nfreqs.lt.nobs+ncov) ipte=4
	if(iopt.ne.1.and.iopt.ne.2) ipte=5
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c   see if theres room in wk (need nobs+2*nfreqs+ncov) :
c
	if(ickse(nfreqs).eq.1) go to 99
	if(ickse(2*nfreqs).eq.1) go to 99
	nws=2*nfreqs+nobs+ncov
	if(ickse(nws).eq.1) go to 99
c
	if(nfreqs.eq.0) then
		call movxy(wk,array(nstart(npx)),4*nobs)
		if(iopt.eq.1) call means(wk,nobs,1,xbar,1)
		on=float(nobs)
		r0=inprod(wk,wk,nobs)/on
		onr0=on*r0
		call ckaddr(args(6),r0,iref)
		if(iref.eq.1.or.ncov.eq.0) go to 99
		do 50 i=1,ncov
  50		wk(nobs+i)=inprod(wk,wk(i+1),nobs-i)/(onr0)
		lab='correlogram'
		call ckadda(vname,ncov,lab,nobs+1,iref)
		go to 99
	 endif
c
c
	if(nfreqs.ne.0) then
		call movct(wk,8*nfreqs,char(0))
		call movxy(wk,array(nstart(npx)),4*nobs)
		if(iopt.eq.1) call means(wk,nobs,1,xbar,1)
		on=float(nobs)
		on1=float(nfreqs)
		r0=inprod(wk,wk,nobs)/on
		call ckaddr(args(6),r0,iref)
		if(iref.eq.1) go to 99
		call fft(wk,wk(nfreqs+1),nfreqs,nfreqs,nfreqs,1)
		if(imflag.eq.1) go to 99
		n1=(nfreqs/2)+1
		do 60 i=1,nfreqs
  60		wk(i)=(wk(i)*wk(i)+wk(nfreqs+i)*wk(nfreqs+i))/on
		lab='Periodogram'
		call ckadda(args(7),n1,lab,1,iref)
		if(iref.eq.1.or.ncov.eq.0) go to 99
c
c
	        call movct(wk(nfreqs+1),4*nfreqs,zero)
        	call fft(wk,wk(nfreqs+1),nfreqs,nfreqs,nfreqs,-1)
	        div=on1*r0
	        do 66 i=1,ncov
  66    	wk(i)=wk(i+1)/div
		lab='Correlogram'
	        call ckadda(vname,ncov,lab,1,iref)
		go to 99
	endif
c
c
  99    continue
        return
        end
c&CUM
	subroutine cum(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to find cumulative sums:
c
c   y=cum(x,n,iopt)  iopt=1 means sums, 2 means averages
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	call ckint(args(2),n)
        if(n.lt.1.or.n.gt.nx) then
		call error(args,2,2)
		go to 99
	endif
	call ckint(args(3),iopt)
	if(iopt.ne.1.and.iopt.ne.2) then
		call error(args,3,2)
		go to 99
	endif
c
c
	call movxy(wk,array(nstart(npx)),4*n)
	if(n.gt.1) then
		do 10 i=2,n
  10		wk(i)=wk(i-1)+wk(i)
	endif
	lab='Cumulative Sums'
	if(iopt.eq.2) then
		do 30 i=1,n
  30		wk(i)=wk(i)/float(i)
		lab='Cumulative Averages'
	endif
c
c
	call ckadda(vname,n,lab,1,iref)
c
c
  99	continue
	return
	end
c&CUMSP
        subroutine cumsp(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to calculate the cumulative spectra of the first
c   (nfreqs/2)+1 values of a spectral array.
c
c   x=cumsp(spec,nfreqs)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),nfreqs)
        n1=(nfreqs/2)+1
        if(n1.ne.nn) then
	        call error(args,2,2)
        	go to 99
	endif
c
c
        c=0.0
        nn=nstart(np)-1
        do 30 i=1,n1
  30    c=c+array(nn+i)
        c1=0.0
        do 40 i=1,n1
        c1=c1+array(nn+i)
  40    wk(i)=c1/c
        lab='cumulative spectra'
        call ckadda(vname,n1,lab,1,iref)
c
c
  99    continue
        return
        end
c&DIFF
        subroutine diff(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to take the nd difference of x(1),...,x(n) and
c   put it in the array whose name is vname. args(1),args(2),
c   and args(3) are assumed to be nd,n, and the array name for x.
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl
c
c
	if(ickl(args(3),np1,nn).eq.1) go to 99
	ipte=0
	call ckint(args(1),n)
	if(n.lt.1.or.n.gt.nn) ipte=2
	call ckint(args(2),nd)
	if(nd.lt.1.or.nd.ge.n) ipte=3
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
        n1=nstart(np1)-1
        n2=n-nd
        do 50 i=1,n2
  50    wk(i)=array(n1+i+nd)-array(n1+i)
c
c
        lab='differenced array'
        call ckadda(vname,n2,lab,1,iref)
  99    continue
        return
        end
c&EXPVEC
        subroutine expvec(args,nargs,vname)
c********************************************************************
c
c   Subroutine to form y(i)=exp(x(i)), i=1,...,n
c
c   y=EXP(x,n)   or y=EXP(x)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,i10*1
	integer*2 ickl,ickr
c
c   handle y=exp(x):
c
	if(nargs.eq.1) then
		call cklist(args(1),xlist,numarr,np)
		if(np.gt.0) then
             		n=nend(np)-nstart(np)+1
			go to 5
		endif
		if(ickr(args(1),1,x,0,0.).eq.1) go to 99
		x=amax1(-87.5,x)
		if(x.gt.87.5) then
 		  	write(*,25)
  25			format(' An element is too big ')
 			imflag=1
 			go to 99
 		endif
		c=exp(x)
		call ckaddr(vname,c,iref)
		go to 99
	endif
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
        call ckint(args(2),n)
        if(n.le.0.or.n.gt.na) then
	        call error(args,2,2)
        	go to 99
	endif
c
c
   5    n1=nstart(np)-1
        do 30 i=1,n
	x=array(n1+i)
 	x=amax1(-87.5,x)
 	if(x.gt.87.5) then
 		write(*,25)
 		imflag=1
 		go to 99
 	endif
  30    wk(i)=exp(array(n1+i))
c
c
        lab='exponentiated array'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&EXTEND
        subroutine extend(args,nargs,vname)
c****************************************************************
c
c   Subroutine to extend a series by antidifferencing.
c
c   y=EXTEND(x,n,m,nd1[,nd2])
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickse
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
        call ckint(args(2),n)
        if(n.lt.1.or.n.gt.nx) then
	        call error(args,2,2)
        	go to 99
	endif
	if(icki(args(3),3,m,1,1).eq.1) go to 99
	if(ickse(m).eq.1) go to 99
	if(ickse(n+m).eq.1) go to 99
        call ckint(args(4),nd1)
	if(nd1.lt.1.or.nd1.ge.n) then
		call error(args,4,2)
		go to 99
	endif
	nd2=0
	if(nargs.eq.5) then
		if(icki(args(5),5,nd2,1,1).eq.1) go to 99
		if(nd1+nd2.ge.n) then
			call error(args,5,2)
			go to 99
		endif
	endif
c
c
	call movxy(wk,array(nstart(npx)),4*n)
	if(nd2.eq.0) then
		c=0.
		do 10 i=1,n-nd1
  10		c=c+wk(i+nd1)-wk(i)
		c=c/float(n-nd1)
		do 20 i=1,m
		npi=n+i
  20		wk(npi)=wk(npi-nd1)+c
		go to 50
	endif
c
c
	nd12=nd1+nd2
	c=0.
	do 30 i=1,n-nd12
	ipnd12=i+nd12
  30	c=c+wk(ipnd12)-wk(ipnd12-nd1)-wk(ipnd12-nd2)+wk(i)
	c=c/float(n-nd12)
        do 40 i=1,m
        npi=n+i
  40    wk(npi)=wk(npi-nd1)+wk(npi-nd2)-wk(npi-nd12)+c
c
c
  50    lab='extended array'
        call ckadda(vname,n+m,lab,1,iref)
c
c
  99    continue
        return
        end
c&EXTRCT
        subroutine extrct(args,nargs,vname)
c***************************************************************
c
c   Subroutine to form a vector by taking a piece out of another
c   array.
c
c   x=EXTRACT(y,nfrst,nlast). if nfrst=nlast, then x is a variable
c   not an array.
c
c***************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,vn1*15
	integer*2 ickl,ickr
c
c   handle B = EXTRACT(A,n,r1,c1,nr,nc)
c

	if(nargs.eq.6) then
		if(ickl(args(1),npa,na).eq.1) go to 99
		ipte=2
		call ckint(args(2),n)
		if(n.lt.1.or.n.gt.na) go to 901
		m=na/n
		ipte=3
		call ckint(args(3),nr1)
		if(nr1.lt.1.or.nr1.gt.n) go to 901
	 	ipte=4
		call ckint(args(4),nc1)
		if(nc1.lt.1.or.nc1.gt.m) go to 901
c
c
		call cklist(args(5),xlist,numarr,npx)
 		if(npx.eq.-999) go to 200
 		if(ickl(args(5),npx,nx).eq.1) go to 99
		call ckint(args(6),nc)
		if(nc.lt.1.or.nc.gt.nx) then
			call error(args,6,2)
			go to 99
		endif
		ns1=nstart(npx)-1
		do 100 i=1,nc
 100		iwk(i)=array(ns1+i)
		nr=n-nr1+1
		go to 110
c
c
200		continue
		ipte=5
		call ckint(args(5),nr)
		if(nr.lt.1.or.nr1+nr-1.gt.n) go to 901
		ipte=6
		call ckint(args(6),nc)
		if(nc.lt.1.or.nc1+nc-1.gt.m) go to 901
		do 105 i=1,nc
 105		iwk(i)=nc1+i-1
c
c
 110	continue
		ii=0
		ns1=nstart(npa)-1
		do 70 ncc=1,nc
		ncol=iwk(ncc)
		if(ncol.lt.1.or.ncol.gt.m) then
			write(*,175)
 175			format(' Illegal column in EXTRACT')
			imflag=1
			go to 99
		endif
		j1=(ncol-1)*n
		do 70 nrow=nr1,nr1+nr-1
		jj=j1+nrow
		ii=ii+1
		wk(ii)=array(ns1+jj)
  70		continue
		lab='Extracted Matrix'
		call ckadda(vname,nr*nc,lab,1,iref)
		go to 99
c
c
 901		call error(args,ipte,2)
		go to 99
	endif
c
c   handle z=extract(x,n,xr,exp,nextr)
c
	if(nargs.eq.3) go to 20
	if(nargs.eq.5) then

		call ckint(args(4),n)
		if(n.ne.-999) go to 25
		if(ickl(args(1),npx,nx).eq.1) go to 99
		if(ickr(args(3),3,xr,3,999.).eq.1) go to 99
		ipte=0
		call ckint(args(2),n)
		if(n.lt.1.or.n.gt.nx) ipte=2
		vn1=args(4)
		call ltoup(vn1,15)
		iopt=0
		if(vn1.eq.'LT') iopt=1
		if(vn1.eq.'LE') iopt=2
		if(vn1.eq.'EQ') iopt=3
		if(vn1.eq.'NE') iopt=4
		if(vn1.eq.'GE') iopt=5
		if(vn1.eq.'GT') iopt=6
		if(iopt.eq.0) ipte=4
		if(ipte.ne.0) then
			call error(args,ipte,2)
			go to 99
		endif
		nsx=nstart(npx)-1
		nextr=0
		if(iopt.eq.1) then
			do 11 i=1,n
			c=array(nsx+i)
    			if(c.ge.xr) go to 11
			nextr=nextr+1
			wk(nextr)=i
  11			continue
		endif
		if(iopt.eq.2) then
			do 12 i=1,n
			c=array(nsx+i)
     			if(c.gt.xr) go to 12
			nextr=nextr+1
      			wk(nextr)=i
  12			continue
		endif
		if(iopt.eq.3) then
			do 13 i=1,n
			c=array(nsx+i)
     			if(c.ne.xr) go to 13
			nextr=nextr+1
			wk(nextr)=i
  13			continue
		endif
		if(iopt.eq.4) then
			do 14 i=1,n
			c=array(nsx+i)
     			if(c.eq.xr) go to 14
			nextr=nextr+1
			wk(nextr)=i
  14			continue
		endif
		if(iopt.eq.5) then
			do 15 i=1,n
			c=array(nsx+i)
     			if(c.lt.xr) go to 15
			nextr=nextr+1
			wk(nextr)=i
  15			continue
		endif
		if(iopt.eq.6) then
			do 16 i=1,n
			c=array(nsx+i)
     			if(c.le.xr) go to 16
			nextr=nextr+1
			wk(nextr)=i
  16			continue
		endif
		call ckaddi(args(5),nextr,iref)
		if(nextr.eq.0.or.iref.eq.1) go to 99
		lab='Extracted Array'
		call ckadda(vname,nextr,lab,1,iref)
		go to 99
	endif
c
c   handle z=extract(x,xind,n)
c
  20	continue
	call cklist(args(2),xlist,numarr,npi)
	if(nargs.eq.3.and.npi.ne.-999) then
		if(ickl(args(1),npx,nx).eq.1) go to 99
		ni=nend(npi)-nstart(npi)+1
		call ckint(args(3),n)
		if(n.lt.1.or.n.gt.ni) then
			call error(args,3,2)
			go to 99
		endif
		nsx=nstart(npx)-1
		nsi=nstart(npi)-1
		do 5 i=1,n
		ii=array(nsi+i)
		if(ii.lt.1.or.ii.gt.nx) then
			write(*,2)
   2			format(' an index outside range of data')
			imflag=1
			go to 99
		endif
   5		wk(i)=array(nsx+ii)
		lab='Extracted Array'
		call ckadda(vname,n,lab,1,iref)
		go to 99
	endif
c
c
  25	continue
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),nfrst)
        call ckint(args(3),nlast)
	ipte=0
	if(nfrst.lt.1.or.nfrst.gt.nn) ipte=2
	if(nlast.lt.nfrst.or.nlast.gt.nn) ipte=3
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
        inc=1
        if(nargs.eq.5) then
	        call ckint(args(4),inc)
        	if(inc.le.0) then
			call error(args,4,2)
			go to 99
	 	endif
	endif
c
c
        n1=nstart(np)-1
        j=0
        do 30 i=nfrst,nlast,inc
        j=j+1
  30    wk(j)=array(n1+i)
        n=j
c
c
	if(nargs.eq.3.and.n.eq.1) then
		call ckaddr(vname,wk(1),iref)
		go to 99
	endif
		lab='extracted array'
		call ckadda(vname,n,lab,1,iref)
        	if(iref.ne.1.and.nargs.eq.5) call ckaddi(args(5),n,iref)
c
c
  99    continue
        return
        end
c&FORMLN
        subroutine formln(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to form a line x(i)=a+b*i, i=1,...,n.
c
c   y=LINE(n,a,b)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 icki,ickr,ickse
c
c
	if(icki(args(1),1,n,1,1).eq.1) go to 99
	if(ickse(n).eq.1) go to 99
	if(nargs.eq.1) then
		a=0.0
		b=1.0
		go to 10
	endif
	if(ickr(args(2),2,a,0,0.).eq.1) go to 99
	if(ickr(args(3),3,b,0,0.).eq.1) go to 99
	if(nargs.eq.4) then
		if(a.gt.b) then
			call error(args,3,2)
			go to 99
		endif
		if(n.eq.1) b=0.
		if(n.gt.1) b=(b-a)/float(n-1)
		do 20 i=1,n
  20		wk(i)=a+float(i-1)*b
		go to 35
	endif
c
c
  10	continue
        do 30 i=1,n
  30    wk(i)=a+b*float(i)
  35    lab='straight line'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
