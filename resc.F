c&RESCN
c&MAXMIN
	subroutine maxmin(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command MAXMIN(x,n,xmax,imax,xmin,imin)
c
c********************************************************************
c
#include 'tslabc'
	character args(nargs)*15
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
	call max(array(nstart(npx)),n,xmax,imax)
	call min(array(nstart(npx)),n,xmin,imin)
	call ckaddr(args(3),xmax,iref)
	if(iref.eq.1) go to 99
	call ckaddi(args(4),imax,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(5),xmin,iref)
	if(iref.eq.1) go to 99
	call ckaddi(args(6),imin,iref)
c
c
  99	continue
	return
	end
c&MAX
	subroutine max(x,n,xmax,ind)
	dimension x(n)
	xmax=x(1)
	ind=1
	if(n.eq.1) return
	do 10 i=2,n
	if(x(i).gt.xmax) then
		xmax=x(i)
		ind=i
		endif
  10	continue
	return
	end
c&MIN
	subroutine min(x,n,xmin,ind)
	dimension x(n)
	xmin=x(1)
	ind=1
	if(n.eq.1) return
	do 10 i=2,n
	if(x(i).lt.xmin) then
		xmin=x(i)
		ind=i
		endif
  10	continue
	return
	end
