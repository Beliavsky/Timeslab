c&MEANVAR
c&DIVSDS
        subroutine divsds(args,nargs,vname)
c***************************************************************
c
c   Subroutine to calculate and subtract periodic means.
c
c   y=divsds(x,n,nper,sds) or y=divsds(x,n,nper,sds,iopt)
c   where iopt=0 means just calculate standard deviations, iopt=1
c   means calculate and divide, and iopt=2 means multiply
c   inputted standard deviations.
c
c****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,ickr,ickse
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),n)
        call ckint(args(3),nper)
	ipte=0
        if(n.le.0.or.n.gt.nn) ipte=2
        if(nper.le.0.or.nper.gt.n) ipte=3
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
	call movxy(wk,array(nstart(np)),4*n)
        iopt=1
        if(nargs.eq.4) go to 60
        call ckint(args(5),iopt)
        if(iopt.lt.0.or.iopt.gt.2) then
		call error(args,5,2)
		go to 99
	endif
        if(iopt.eq.0.or.iopt.eq.1) go to 60
        if(nper.eq.1) then
		if(ickr(args(4),4,wk(n+1),2,0.).eq.1) go to 99
		go to 60
	endif
     	if(ickl(args(4),np1,na).eq.1) go to 99
        if(nper.gt.na) then
		call error(args,3,2)
		go to 99
	endif
        np2=nstart(np1)-1
        do 55 i=1,nper
        wk(n+i)=array(np2+i)
  	if(wk(n+i).le.0.) then
		write(*,51)
  51		format(' An element of sds is not positive')
		imflag=1
		go to 99
		endif
  55	continue
  60    continue
	nws=n+nper
	if(ickse(nws).eq.1) go to 99
c
c
c
        call stddev(wk,n,nper,wk(n+1),iopt,ier)
	if(ier.eq.1) then
		write(*,51)
		imflag=1
		go to 99
	endif
c
c   add arrays :
c
        iref=0
        if(iopt.eq.2) go to 85
        if(nper.gt.1) go to 70
        call ckaddr(args(4),wk(n+1),iref)
        go to 80
  70    continue
        lab='standard deviations'
        call ckadda(args(4),nper,lab,n+1,iref)
  80    if(iref.eq.1.or.iopt.eq.0) go to 99
        lab='series divided by periodic sds'
  85    if(iopt.eq.2) lab='series multiplied by periodic sds'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&MEANS
        subroutine means(x,n,nper,xbar,iopt)
c*******************************************************************
c
c   If iopt=0 or 1, this subroutine calculates the nper-periodic
c   means xbar(1),...,xbar(nper) of data x(1),...,x(n).
c   If iopt=1 these means are subtracted. If iopt=2,
c   inputted means are added.
c
c*******************************************************************
c
        dimension x(n),xbar(nper)
        if(iopt.eq.2) go to 20
        do 10 i=1,nper
                c=0.0
                nc=0.0
                do 5 j=i,n,nper
                nc=nc+1
   5            c=c+x(j)
  10    xbar(i)=c/float(nc)
        if(iopt.eq.0) go to 99
c
c
  20    do 30 i=1,nper
        c=xbar(i)
        if(iopt.eq.2) c=-c
                do 25 j=i,n,nper
  25            x(j)=x(j)-c
  30    continue
c
c
  99    continue
        return
        end
c&STDDEV
        subroutine stddev(x,n,nper,sd,iopt,ier)
c*******************************************************************
c
c   If iopt=0 or 1, this subroutine calculates the nper-periodic
c   standard deviations sd(1),...,sd(nper) of data x(1),...,x(n).
c   if iopt=1 these standard deviations are divided. If iopt=2,
c   inputted standard deviations are multiplied.
c
c*******************************************************************
c
        dimension x(n),sd(nper)
	ier=0
        if(iopt.eq.2) go to 20
	call means(x,n,nper,sd,0)
        do 10 i=1,nper
	xb=sd(i)
                c=0.0
                nc=0
                do 5 j=i,n,nper
                nc=nc+1
   5            c=c+(x(j)-xb)**2
  10    sd(i)=sqrt(c/float(nc))
        if(iopt.eq.0) go to 99
c
c
  20    do 30 i=1,nper
        c=sd(i)
	if(c.le.0.0) then
		ier=1
		go to 99
	endif
        if(iopt.eq.2) c=1./c
                do 25 j=i,n,nper
  25            x(j)=x(j)/c
  30    continue
c
c
  99    continue
        return
        end
c&SUBMNS
        subroutine submns(args,nargs,vname)
c***************************************************************
c
c   Subroutine to calculate and subtract periodic means.
c
c   y=submns(x,n,nper,xbar) or y=submns(x,n,nper,xbar,iopt)
c   where iopt=0 means only xbar is calculated, 1 means xbar
c   is subtracted, and 2 means inputted means are added.
c
c****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,ickr,ickse
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),n)
        call ckint(args(3),nper)
	ipte=0
        if(n.le.0.or.n.gt.nn) ipte=2
        if(nper.le.0.or.nper.gt.n) ipte=3
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
	call movxy(wk,array(nstart(np)),4*n)
        iopt=1
        if(nargs.eq.4) go to 60
        call ckint(args(5),iopt)
        if(iopt.lt.0.or.iopt.gt.2) then
		call error(args,5,2)
		go to 99
	endif
        if(iopt.eq.0.or.iopt.eq.1) go to 60
	if(nper.eq.1) then
		if(ickr(args(4),4,wk(n+1),0,0.).eq.1) go to 99
		go to 60
	endif
	if(ickl(args(4),np1,na).eq.1) go to 99
        if(nper.gt.na) then
		call error(args,3,2)
		go to 99
	endif
        np2=nstart(np1)-1
        do 55 i=1,nper
  55    wk(n+i)=array(np2+i)
  60    continue
	if(ickse(n+nper).eq.1) go to 99
c
c
        call means(wk,n,nper,wk(n+1),iopt)
c
c   add arrays :
c
        iref=0
        if(iopt.eq.2) go to 85
        if(nper.gt.1) go to 70
        call ckaddr(args(4),wk(n+1),iref)
        go to 80
  70    continue
        lab='mean'
        call ckadda(args(4),nper,lab,n+1,iref)
  80    if(iref.eq.1.or.iopt.eq.0) go to 99
        lab='series minus means'
  85    if(iopt.eq.2) lab='series plus means'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end


