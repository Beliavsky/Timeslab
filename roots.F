c&ROOTS
c&RTPOLY
	subroutine rtpoly(args,nargs)
c******************************************************************
c
c   Subroutine to handle the command
c
c   ROOTSPOLY(rr,ri,np,a)
c
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl
c
c
	if(ickl(args(1),nprr,nrr).eq.1) go to 99
	if(ickl(args(2),npri,nri).eq.1) go to 99
	call ckint(args(3),np)
	if(np.lt.1.or.np.gt.min0(nrr,nri)) then
		call error(args,3,2)
		go to 99
	endif
c
c
	n1=nstart(nprr)
	n2=nstart(npri)
	np4=4*np
	call movxy(wk,array(n1),np4)
	call movxy(wk(np+1),array(n2),np4)
	do 5 i=1,np
	if(abs(wk(i)).lt.1.e-20.and.abs(wk(np+i)).lt.1.e-20) then
		write(*,6)
   6		format(' One of the roots is too small')
		imflag=1
		go to 99
	endif
   5	continue
	n1=np+1
	n2=n1+np
	n3=n2+np
	n4=n3+2*(np+1)
	n5=n4+2*(np+1)
	call rtar(wk,wk(n1),np,wk(n2),wk(n3),wk(n4),wk(n5),iwk)
c
c
	lab='Polynomial Coefficients'
	call ckadda(args(4),np,lab,n2,iref)
c
c
  99	continue
	return
	end
c&RTAR
	subroutine rtar(rr,ri,n,alpha,r,c,cf,j)
c****************************************************************
c
c   subroutine to find the coefficients of a polynomial given it zeros
c
c********************************************************************
c
	complex r(1),c(1),pr,sum,cone
	dimension rr(n),ri(n),alpha(n),cf(1),j(1)
c
c
	cone=cmplx(1.0,0.)
	nn=n+1
	do 10 i=1,n
  10	r(i)=cone/cmplx(rr(i),ri(i))
	cf(nn)=1.
	do 14 m=1,n
	sum=cmplx(0.0,0.0)
	l=1
	j(1)=1
	go to 2
   1	j(l)=j(l)+1
   2	if(l-m) 3,5,50
   3	mm=m-1
	do 4 i=l,mm
	ii=i+1
   4	j(ii)=j(i)+1
   5	pr=cmplx(1.0,0.0)
	do 7 i=1,m
	ick=j(i)
   7	pr=-pr*r(ick)
	sum=sum+pr
	do 6 i=1,m
	l=m-i+1
	if(j(l)-n+m-l) 1,6,50
   6	continue
	mp=n-m+1
  14	cf(mp)=real(sum)
	do 15 i=1,n
  	nn=n-i+1
  15	alpha(i)=cf(nn)
	return
  50	ier=1
	return
	end
c&POLYRT
	subroutine polyrt(args,nargs)
c******************************************************************
c
c   Subroutine to handle the command
c
c   POLYROOTS(a,a0,np,maxit,rootr,rooti,ier)
c
c   roots of polynomial a0 + a(1)x + a(2)x**2 + ... + a(np)x**np
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,ickr,icki
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(ickr(args(2),2,a0,0,0.).eq.1) go to 99
	if(icki(args(3),3,np,1,1).eq.1) go to 99
	if(np.gt.na.or.np.gt.1000) then
		call error(args,3,2)
		go to 99
	endif
	if(icki(args(4),4,maxit,1,1).eq.1) go to 99
	if(maxit.gt.1000) then
		call error(args,4,2)
		go to 99
	endif
c
c
	wk(1)=a0
	call movxy(wk(2),array(nstart(npa)),4*np)
	nn=2*(np+1)
	n1=np+2
	n2=n1+nn
	n3=n2+nn
	n4=n3+nn
	n5=n4+np
	call proots(wk,np,maxit,wk(n1),wk(n2),wk(n3),wk(n4),wk(n5),ier)
c
c
	call ckaddi(args(7),ier,iref)
	if(ier.eq.1) then
		write(*,10)
  10		format(' Convergence not reached in POLYROOT')
		imflag=1
		go to 99
	endif
	if(iref.eq.1) go to 99
	lab='Real Part of Polynomial Roots'
	call ckadda(args(5),np,lab,n4,iref)
	if(iref.eq.1) go to 99
	lab='Imaginary Part of Polynomial Roots'
	call ckadda(args(6),np,lab,n5,iref)
c
c
  99	continue
	return
	end
c&PROOTS
	subroutine proots(ra,np,maxit,a,a1,roots,rr,ri,ier)
c*****************************************************************
c
c   calculate roots of a(1)+a(2)z+a(3)z**2+...+a(np+1)z**np
c
c*****************************************************************
c
	dimension ra(1),rr(1),ri(1)
	complex a(1),a1(1),roots(1),x,b,c
c
c
	eps=1.e-6
	do 1 i=1,np+1
	a(i)=ra(i)
   1	a1(i)=a(i)
	do 3 j=np,1,-1
		x=cmplx(0.,0.)
		call root1(a1,j,x,eps,0,maxit,ier)
		if(ier.eq.1) go to 99
	if(abs(aimag(x)).le.2.*eps**2*abs(real(x))) x=cmplx(real(x),0.)
		roots(j)=x
		b=a1(j+1)
		do 2 jj=j,1,-1
			c=a1(jj)
			a1(jj)=b
   2			b=x*b+c
   3	continue
	do 4 j=1,np
    	call root1(a,np,roots(j),eps,1,maxit,ier)
   4	if(ier.eq.1) go to 99
	do 6 j=2,np
		x=roots(j)
		do 5 i=j-1,1,-1
			if(real(roots(i)).le.real(x)) go to 10
			roots(i+1)=roots(i)
   5		continue
  	i=0
  10	roots(i+1)=x
   6	continue
	do 20 i=1,np
	rr(i)=real(roots(i))
  20	ri(i)=aimag(roots(i))
  99	continue
	return
	end
c&ROOT1
	subroutine root1(a,np,x,eps,iopt,maxit,ier)
c********************************************************************
c
c   calculate one root.
c
c********************************************************************
c
	complex a(1),x,dx,x1,b,d,f,g,h,sq,gp,gm,g2,zero
c
c
	ier=0
	zero=cmplx(0.,0.)
	if(iopt.eq.1) then
		dxold=cabs(x)
		npol=0
	endif
	do 2 it=1,maxit
		b=a(np+1)
		d=zero
		f=zero
		do 1 j=np,1,-1
			f=x*f+d
			d=x*d+b
   1		b=x*b+a(j)
		if(cabs(b).le.1.e-15) then
			dx=zero
		else if(cabs(d).le.1.e-15.and.cabs(f).le.1.e-15) then
			dx=cmplx(cabs(b/a(np+1))**(1./np),0.)
		else
			g=d/b
			g2=g*g
			h=g2-2.*f/b
			sq=csqrt((np-1)*(np*h-g2))
			gp=g+sq
			gm=g-sq
			if(cabs(gp).lt.cabs(gm)) gp=gm
			dx=np/gp
		endif
		x1=x-dx
		if(x.eq.x1) return
		x=x1
		if(iopt.eq.1) then
			npol=npol+1
			cdx=cabs(dx)
			if(npol.gt.9.and.cdx.ge.dxold) return
			dxold=cdx
		else
			if(cabs(dx).le.eps*cabs(x)) return
		endif
   2	continue
	ier=1
	return
	end
