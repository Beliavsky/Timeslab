c&DTMXLK
c&DTARMA
	subroutine dtarma(args,nargs)
c********************************************************************
c
c   Subroutine to handle the command
c
c   dtarma(y,n,np,nq,maxit,eps,alpha,beta,rvar,am2ll,ier,var)
c
c********************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,ickr,icki,ickse
c
c
	if(ickl(args(1),npy,ny).eq.1) go to 99
	call ckint(args(2),n)
	if(n.lt.1.or.n.gt.ny) then
		call error(args,2,2)
		go to 99
	endif
	if(icki(args(3),3,np,1,0).eq.1) go to 99
	if(icki(args(4),4,nq,1,0).eq.1) go to 99
	if(np.eq.0.and.nq.eq.0) then
		write(*,10)
  10		format(' One of np,nq must be positive')
		go to 99
	endif
	if(icki(args(5),5,maxit,1,1).eq.1) go to 99
	if(ickr(args(6),6,eps,2,0.).eq.1) go to 99
	if(np.gt.0) then
		if(ickl(args(7),npa,na).eq.1) go to 99
		if(np.gt.na) then
			call error(args,3,2)
			go to 99
		endif
	endif
	if(nq.gt.0) then
		if(ickl(args(8),npb,nb).eq.1) go to 99
		if(nq.gt.nb) then
			call error(args,4,2)
			go to 99
		endif
	endif
c
c
	if(np.gt.0) then
	call movxy(wk,array(nstart(npa)),4*np)
	call arprt(wk,np,ier)
		if(ier.ne.0) then
			write(*,11)
  11			format(' Starting Values Not Stationary')
			call ckaddi(args(11),ier,iref)
			imflag=1
			go to 99
		endif
	endif
	if(nq.gt.0) then
		call movxy(wk(np+1),array(nstart(npb)),4*nq)
		call arprt(wk(np+1),nq,ier)
		if(ier.ne.0) then
			write(*,12)
  12			format(' Starting Values Not Invertible')
			call ckaddi(args(11),ier,iref)
			imflag=1
			go to 99
		endif
	endif
c	write(*,13) (wk(j),j=1,np+nq)
c 13	format(1x,5f14.6)
	if(np.gt.0) then
		do 15 i=1,np
  15		wk(i)=-alog((1.-wk(i))/(1.+wk(i)))
	endif
	if(nq.gt.0) then
		do 20 i=1,nq
  20		wk(np+i)=-alog((1.-wk(np+i))/(1.+wk(np+i)))
	endif
	cws1=float(np+nq)
	cws2=float(n)
	cws3=float(max0(np,nq+1)+1)
	cws=3.*cws1+2.*cws2+5.*cws3+cws1*(cws1+1.)
	if(cws.gt.32767.) then
		if(ickse(maxnum+1).eq.1) go to 99
	endif
	nws=3*(np+nq)+2*n+5*(max0(np,nq+1)+1)+(np+nq)*(np+nq+1)
	if(ickse(nws).eq.1) go to 99
	nppnq=np+nq
	n1=nppnq+1
	n2=n1+nppnq
	n3=n2+nppnq
	n4=n3+n
	n5=n4+n
	mxpq=max0(np,nq+1)+1
	n6=n5+mxpq
	n7=n6+mxpq
	n8=n7+mxpq
	n9=n8+mxpq
	n10=n9+mxpq
	n11=n10+(np+nq)*(np+nq+1)
	do 30 i=1,nppnq
  30	wk(n2+i-1)=1.
	call movxy(wk(n3),array(nstart(npy)),4*n)
	call nelmin(nppnq,wk,wk(n1),ynewlo,eps,wk(n2),10,maxit,wk(n3),
     1	n,wk(n4),wk(n5),wk(n6),wk(n7),np,nq,wk(n8),wk(n9),np+nq,wk(n10),
     1	rv,am2ll,ier,std)
	call ckaddi(args(11),ier,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(10),am2ll,iref)
	if(iref.eq.1) go to 99
	if(np.gt.0) then
		lab='Autoregressive Coefficients'
		call ckadda(args(7),np,lab,n8,iref)
		if(iref.eq.1) go to 99
	endif
	if(nq.gt.0) then
		lab='Moving Average Coefficients'
		call ckadda(args(8),nq,lab,n9,iref)
		if(iref.eq.1) go to 99
	endif
	call ckaddr(args(9),rv,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(10),am2ll,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(12),var,iref)
c
c
  99	continue
	return
	end
c&NELMIN
	subroutine nelmin(n,start,min,ynewlo,reqmin,step,konvge,icount,
     1	dat,ny,e,vk,vl,vw,np,nq,alpha,beta,nppnq,p,rv,am2ll,ier,var)
c*********************************************************************
c
c   Algorithm AS 47 Applied Statistics (JRSS C)
c   (1971) Vol.20, No. 3
c
c   The Nelder-Mead Simplex Minimization Procedure
c
c
c     Purpose:  To find the minimum value of a user-specified function
c
c
c     Formal Parameters:
c
c       n: input  : the number of variables being minimized over
c   start: input  : real*8 array containing starting values
c     min: output : real*8 array containing minimizing values
c  ynewlo: output : real*8 variable containing minimum value of function
c  reqmin: input  : the terminating limit for variance of function vals
c    step: input  : real*8 array determining the size and shape of the
c                   initial simplex.  The relative magnitudes of its n
c                   elements should reflect the units of the n variables
c  konvge: input  : the convergence check is performed every konvge
c                   iterations
c  icount: input  : maximum number of function evaluations
c          output : function evaluations performed + 1000*number of
c                   restarts.  Negative icount value identifies input
c                   parameter faults(s). -1 means input reqmin.le.0,
c                   -10 means n.gt.20, -100 means konvge.le.0.
c
c     All variables and arrays are to be declared in the calling
c     program as double precision.
c
c     Auxiliary Algorithm:  The double precision function FN(A) finds
c     the function value at point A which is double precision with
c     n elements.
c
c     Reference Nelder,J.A. and Mead,R. (1965). A Simplex Method
c     For Function Minimization. Computer Journal, Vol. 7,308-313.
c
c********************************************************************
c
	Real*4 start(n),min(n),ynewlo,reqmin,step(n),
     1	p(nppnq,1),pstar(20),p2star(20),pbar(20),y(20),
     2	dn,dnn,z,ylo,rcoeff,ystar,ecoeff,y2star,ccoeff,
     3	curmin,del,armalk
	dimension dat(ny),alpha(1),beta(1)
	dimension e(ny),vw(1),vk(1),vl(1)
	double precision sum,summ,cmin
c
	data rcoeff,ecoeff,ccoeff/1.,2.,.5/
c
c
	ier=0
	kcount=icount
	icount=0
	if(reqmin.le.0.) icount=icount-1
	if(n.gt.20) icount=icount-10
	if(konvge.le.0) icount=icount-100
	if(icount.lt.0) return
c
	jcount=konvge
	dn=float(n)
	nn=n+1
	dnn=float(nn)
	del=1.
c
c   construction of initial simplex:
c
1001	do 1 i=1,n
   1	p(i,nn)=start(i)
	z=armalk(start,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	y(nn)=z
	sum=dble(z)
	summ=dble(z*z)
	do 2 j=1,n
	start(j)=start(j)+step(j)*del
		do 3 i=1,n
   3		p(i,j)=start(i)
	z=armalk(start,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	y(j)=z
	sum=sum+dble(z)
	summ=summ+dble(z)*dble(z)
   2	start(j)=start(j)-step(j)*del
c
c   simplex construction complete
c
c   find highest and lowest y values. ynewlo (=y(ihi)) indicates
c   the vertex of the simplex to be replaced.
c
1000	ylo=y(1)
	ynewlo=ylo
	ilo=1
	ihi=1
	do 5 i=2,nn
	if(y(i).ge.ylo) go to 4
	ylo=y(i)
	ilo=i
   4	if(y(i).le.ynewlo) go to 5
	ynewlo=y(i)
	ihi=i
   5	continue
	sum=sum-dble(ynewlo)
	summ=summ-dble(ynewlo)*dble(ynewlo)
c
c   calculate pbar. The centroid of the simplex vertices
c   excepting that with y value ynewlo.
c
	do 7 i=1,n
	z=0.
		do 6 j=1,nn
   6		z=z+p(i,j)
	z=z-p(i,ihi)
   7	pbar(i)=z/dn
c
c   reflection through the centroid
c
	do 8 i=1,n
   8	pstar(i)=(1.+rcoeff)*pbar(i)-rcoeff*p(i,ihi)
	ystar=armalk(pstar,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	icount=icount+1
	if(ystar.ge.ylo) go to 12
c
c   successful reflection, so extension
c
	do 9 i=1,n
   9	p2star(i)=ecoeff*pstar(i)+(1.-ecoeff)*pbar(i)
	y2star=armalk(p2star,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	icount=icount+1
c
c   retain extension or contraction
c
	if(y2star.ge.ylo) go to 19
  10	do 11 i=1,n
  11	p(i,ihi)=p2star(i)
	y(ihi)=y2star
	go to 900
c  no extension
  12	l=0
	do 13 i=1,nn
	if(y(i).gt.ystar) l=l+1
  13	continue
	if(l.gt.1) go to 19
	if(l.eq.0) go to 15
c
c  contraction on the reflection side of the centroid
c
	do 14 i=1,n
  14	p(i,ihi)=pstar(i)
	y(ihi)=ystar
c
c   contraction on the y(ihi) side of the centroid
c
  15	do 16 i=1,n
  16	p2star(i)=ccoeff*p(i,ihi)+(1.-ccoeff)*pbar(i)
	y2star=armalk(p2star,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	icount=icount+1
	if(y2star.le.y(ihi)) go to 10
c
c   contract whole simplex
c
	sum=0.d0
	summ=0.d0
	do 18 j=1,nn
		do 17 i=1,n
		p(i,j)=(p(i,j)+p(i,ilo))*0.5
  17		min(i)=p(i,j)
	y(j)=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	sum=sum+dble(y(j))
  18	summ=summ+dble(y(j))*dble(y(j))
	icount=icount+nn
	go to 901
c   retain reflection
  19	do 20 i=1,n
  20	p(i,ihi)=pstar(i)
	y(ihi)=ystar
 900	sum=sum+dble(y(ihi))
	summ=summ+dble(y(ihi))*dble(y(ihi))
 901	jcount=jcount-1
	if(jcount.ne.0) go to 1000
c
c   check to see if minimum reached
c
	if(icount.gt.kcount) go to 22
	jcount=konvge
	cmin=(summ-(sum*sum)/dble(dnn))/dble(dn)
	curmin=cmin
c
c   curmin is the variance of the n+1 fn values at the vertices
c
	if(curmin.ge.reqmin) go to 1000
c
c   factorial tests to check that ynewlo is a local minimum
c
  22	do 23 i=1,n
  23	min(i)=p(i,ihi)
	ynewlo=y(ihi)
	if(icount.gt.kcount) then
		ier=1
		return
	endif
	do 24 i=1,n
	del=step(i)*.001
	min(i)=min(i)+del
	z=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	if(z.lt.ynewlo) go to 25
	min(i)=min(i)-del-del
	z=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	if(z.lt.ynewlo) go to 25
  24	min(i)=min(i)+del
	am2ll=z
	var=curmin
	return
c
c   restart procedure
c
  25	do 26 i=1,n
  26	start(i)=min(i)
	del=.001
	icount=icount+1
	go to 1001
	end
c&ARMALK
	function armalk(theta,np,nq,n,y,alpha,beta,e,vw,vl,vk,rv,ier)
c********************************************************************
c
c   Subroutine to calculate -2log likelihood for ARMA(np,nq) process
c
c*********************************************************************
c
	dimension theta(1),y(n),alpha(1),beta(1),
     1	e(n),vw(1),vl(1),vk(1)
c
c   Find -alpha and -beta corresponding to input transformed partials:
c
	if(np.eq.0) go to 15
	do 5 i=1,np
	e1=exp(-theta(i))
   5	alpha(i)=(1.-e1)/(1.+e1)
c	write(*,6) (alpha(i),i=1,np)
c  6	format(1x,5f14.6)
	call partar(alpha,np,alpha)
	do 10 i=1,np
  10	alpha(i)=-alpha(i)
  15	if(nq.eq.0) go to 30
	do 20 i=1,nq
   	e1=exp(-theta(np+i))
  20	beta(i)=(1.-e1)/(1.+e1)
c	write(*,6) (beta(i),i=1,nq)
	call partar(beta,nq,beta)
	do 25 i=1,nq
  25	beta(i)=-beta(i)
  30	continue
c
c   call evaluator:
c
	tol=.0001
	call lkhood(alpha,np,beta,nq,y,e,n,sumsq,fact,vw,vl,vk,tol,ier)
c	write(*,31) ier
c 31	format(' ier = ',i5)
	if(ier.gt.0) return
	if(np.gt.0) then
		do 35 i=1,np
  35		alpha(i)=-alpha(i)
	endif
	if(nq.gt.0) then
 		do 40 i=1,nq
  40		beta(i)=-beta(i)
	endif
	rv=sumsq/float(n)
	al2pi=alog(8.*atan(1.0))
	c=float(n)*(1.+al2pi+alog(rv)+alog(fact))
c	write(*,123) c
c123	format(1x,f14.6)
	armalk=c
	return
	end
c&LKHOOD
	subroutine lkhood(alpha,np,beta,nq,y,e,n,sumsq,fact,vw,vl,
     1	vk,tol,ier)
c********************************************************************
c
c   Subroutine to find the 2 terms in the exact ARMA likelihood for
c   -alpha and -beta as in Melard (Applied Statistics, 1984, 104)
c
c   vw and vl are max(p,q+1)+1 vk is max(p,q+1)
c
c********************************************************************
c
	dimension alpha(1),beta(1),y(n),e(n),vw(1),vl(1),vk(1)
	data epsil1/1.e-10/
c
c
	mr=max0(np,nq+1)
	mrp1=mr+1
	fact=0.
	detman=1.0
	detcar=0.
	sumsq=0.
	mxpq=max0(np,nq)
	mxpqp1=mxpq+1
	mqp1=nq+1
	mpp1=np+1
c
c
	call wilson(alpha,np,beta,nq,vw,mxpqp1,vl,mxpqp1,vk,mxpq,ier)
	if(ier.gt.0) return
	vk(1)=vw(1)
	if(mr.eq.1) go to 15
	do 14 k=2,mr
		vk(k)=0.
		if(k.gt.np) go to 12
		do 11 j=k,np
  11			vk(k)=vk(k)+alpha(j)*vw(j+2-k)
  12		if(k.gt.mqp1) go to 14
		do 13 j=k,mqp1
  13			vk(k)=vk(k)-beta(j-1)*vl(j+1-k)
  14	continue
c
c
  15	r=vk(1)
	vl(mr)=0.
	do 16 j=1,mr
		vw(j)=0.
		if(j.ne.mr) vl(j)=vk(j+1)
		if(j.le.np) vl(j)=vl(j)+alpha(j)*r
  16    vk(j)=vl(j)
c
c
	last=mpp1-nq
	loop=np
	jfrom=mpp1
	vw(mpp1)=0.
	vl(mxpqp1)=0.
c
c
	if(n.le.0) go to 50
	do 29 i=1,n
c
c
		if(i.ne.last) go to 17
		loop=min0(np,nq)
		jfrom=loop+1
c
c
		if(nq.le.0) go to 30
  17		if(r.le.epsil1) go to 40
		if(abs(r-1.).lt.tol.and.i.gt.mxpq) go to 30
c
c
		detman=detman*r
  19		if(abs(detman).lt.1.) go to 20
		detman=detman*.0625
		detcar=detcar+4.
		go to 19
  20		if(abs(detman).ge..0625) go to 21
		detman=detman*16.
		detcar=detcar-4.
		go to 20
  21		vw1=vw(1)
		a=y(i)-vw1
		e(i)=a/sqrt(r)
		aor=a/r
		sumsq=sumsq+a*aor
		vl1=vl(1)
		alf=vl1/r
		r=r-alf*vl1
		if(loop.eq.0) go to 23
c
c
		do 22 j=1,loop
			flj=vl(j+1)+alpha(j)*vl1
			vw(j)=vw(j+1)+alpha(j)*vw1+aor*vk(j)
			vl(j)=flj-alf*vk(j)
  22 		vk(j)=vk(j)-alf*flj
  23		if(jfrom.gt.nq) go to 25
		do 24 j=jfrom,nq
			vw(j)=vw(j+1)+aor*vk(j)
			vl(j)=vl(j+1)-alf*vk(j)
  24		vk(j)=vk(j)-alf*vl(j+1)
  25		if(jfrom.gt.np) go to 27
		do 26 j=jfrom,np
  26		vw(j)=vw(j+1)+alpha(j)*y(i)
  27		continue
  29	continue
	go to 39
c
c
  30	nexti=i
	ier=-nexti
	do 31 i=nexti,n
  31	e(i)=y(i)
        if(np.eq.0) go to 34
        do 33 i=nexti,n
		do 32 j=1,np
  32		e(i)=e(i)-alpha(j)*y(i-j)
  33	continue
  34	if(nq.eq.0) go to 37
	do 36 i=nexti,n
		do 35 j=1,nq
  35		e(i)=e(i)+beta(j)*e(i-j)
  36	continue
c
c
  37	do 38 i=nexti,n
  38	sumsq=sumsq+e(i)*e(i)
  39	fn=n
	fact=detman**(1./fn)*2.**(detcar/fn)
	return
  40	ier=8
	return
  50	ier=9
	return
	end
c&WILSON
	subroutine wilson(alpha,np,beta,nq,acf,ma,cvli,mxpqp1,alph,
     1	mxpq,ier)
c********************************************************************
c
c   subroutine to find mx covariances
c
c********************************************************************
c
	dimension alpha(1),beta(1),acf(ma),cvli(mxpqp1),alph(mxpq)
	data epsil2/1.e-10/
c
c
	ier=0
	if(np.lt.0.or.nq.lt.0) ier=1
	if(mxpq.ne.max0(np,nq)) ier=2
	if(mxpqp1.ne.mxpq+1) ier=3
	if(ma.lt.mxpqp1) ier=4
	if(ier.gt.0) return
c
c
	acf(1)=1.
	cvli(1)=1.
	if(ma.eq.1) return
	do 1 i=2,ma
   1	acf(i)=0.
	if(mxpqp1.eq.1) return
	do 2 i=2,mxpqp1
   2	cvli(i)=0.
	do 9 k=1,mxpq
   9	alph(k)=0.
c
c
	if(nq.eq.0) go to 18
	do 13 k=1,nq
		cvli(k+1)=-beta(k)
		acf(k+1)=-beta(k)
		kc=nq-k
		if(kc.eq.0) go to 12
		do 11 j=1,kc
  11		acf(k+1)=acf(k+1)+beta(j)*beta(j+k)
  12		acf(1)=acf(1)+beta(k)*beta(k)
  13	continue
c
c
  18	if(np.eq.0) return
	do 19 k=1,np
	alph(k)=alpha(k)
  19	cvli(k)=alpha(k)
c
c
	do 29 k=1,mxpq
		kc=mxpq-k
		if(kc.ge.np) go to 24
		div=1.-alph(kc+1)*alph(kc+1)
		if(div.le.epsil2) go to 70
		if(kc.eq.0) go to 29
		do 23 j=1,kc
  23		alph(j)=(cvli(j)+alph(kc+1)*cvli(kc+1-j))/div
  24		if(kc.ge.nq) go to 26
		j1=max0(kc+1-np,1)
		do 25 j=j1,kc
  25		acf(j+1)=acf(j+1)+acf(kc+2)*alph(kc+1-j)
  26		if(kc.ge.np) go to 29
		do 27 j=1,kc
  27		cvli(j)=alph(j)
  29	continue
c
c
	acf(1)=.5*acf(1)
	do 33 k=1,mxpq
	if(k.gt.np) go to 33
	div=1.-alph(k)*alph(k)
		do 31 j=1,k+1
  31		cvli(j)=(acf(j)+alph(k)*acf(k+2-j))/div
		do 32 j=1,k+1
  32		acf(j)=cvli(j)
  33	continue
c
c
	do 43 i=1,ma
	miim1p=min0(i-1,np)
	if(miim1p.eq.0) go to 43
		do 42 j=1,miim1p
  42		acf(i)=acf(i)+alpha(j)*acf(i-j)
  43	continue
	acf(1)=acf(1)*2.
c
c
	cvli(1)=1.
	if(nq.le.0) go to 60
	do 53 k=1,nq
	cvli(k+1)=-beta(k)
	if(np.eq.0) go to 53
	mikp=min0(k,np)
		do 52 j=1,mikp
  52		cvli(k+1)=cvli(k+1)+alpha(j)*cvli(k+1-j)
  53	continue
c
  60	return
c
c
  70	ier=5
	return
	end
