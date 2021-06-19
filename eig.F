c&EIG
c&EIG
	subroutine eig(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   val=EIG(a,n,maxit,ier[,vec])
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname
	integer*2 ickl,icki,ickdg
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(icki(args(2),2,n,5,56).eq.1) go to 99
	ns=n*n
	if(n.lt.1.or.n*n.gt.na) then
		call error(args,2,2)
		go to 99
	endif
	if(icki(args(3),3,maxit,1,1).eq.1) go to 99
	if(maxit.gt.500) then
		call error(args,3,2)
		go to 99
	endif
c
c
	nbytes=4*n*n
	call movxy(wk,array(nstart(npa)),nbytes)
	call max(wk,ns,xmax,imax)
	call min(wk,ns,xmin,imin)
	xsc=xmax
	if(abs(xmin).gt.abs(xmax)) xsc=xmin
	if(xsc.eq.0.) then
		write(*,102)
 102		format(' Matrix all zeros')
		imflag=1
		go to 99
	endif
	do 3 i=1,ns
   3	wk(i)=wk(i)/xsc
	iptvec=0
	if(nargs.eq.5) iptvec=1
	ns=n*n
	ns1=ns+1
	ns2=ns1+n
	ns3=ns2+2*n*n
	ns4=ns3+2*n
	call eig2(wk,wk(ns1),wk(ns2),wk(ns3),wk(ns4),n,iptvec,
     1	maxit,ier)
	if(ier.eq.1) then
		write(*,2)
   2		format(' Nonconvergence in EIG')
		imflag=1
	endif
	call ckaddi(args(4),ier,iref)
	if(iref.eq.1.or.ier.ne.0) go to 99
	lab='Eigenvalues'
	do 4 i=1,n
   4	wk(ns+i)=wk(ns+i)*xsc
	call ckadda(vname,n,lab,ns+1,iref)
	if(iref.eq.1.or.iptvec.eq.0) go to 99
	lab='Eigenvectors'
	call ckadda(args(5),ns,lab,1,iref)
c
c
  99	continue
	return
	end
c&EIG2
	subroutine eig2(a,vals,wk1,wk2,wk3,n,iptvec,maxit,ier)
c*******************************************************************
c
c
c
c
c******************************************************************
c
	dimension a(n,n),vals(n)
	double precision wk1(n,n),wk2(n),wk3(n)
	double precision h,scale,f,g,hh,dd,r,s,c,p,b
c
c
c   Tridiagonalize:
c

	if(n.lt.1) then
		ier=2
		go to 99
	endif
	if(n.eq.1) then
		vals(1)=a(1,1)
		ier=0
		go to 99
	endif
	do 5 i=1,n
	do 5 j=1,n
   5	wk1(i,j)=a(i,j)
	do 18 i=n,2,-1
		l=i-1
		h=0.d0
		scale=0.d0
		if(l.gt.1) then
			do 11 k=1,l
     				scale=scale+dabs(wk1(i,k))
  11			continue
		if(scale.eq.0.d0) then
			wk2(i)=wk1(i,l)
		else
			do 12 k=1,l
				wk1(i,k)=wk1(i,k)/scale
				h=h+wk1(i,k)**2
  12 			continue
			f=wk1(i,l)
			g=-dsign(dsqrt(h),f)
			wk2(i)=scale*g
			h=h-f*g
			wk1(i,l)=f-g
			f=0.d0
			do 15 j=1,l
			if(iptvec.eq.1) wk1(j,i)=wk1(i,j)/h
			g=0.d0
			do 13 k=1,j
				g=g+wk1(j,k)*wk1(i,k)
  13			continue
			if(l.gt.j) then
				do 14 k=j+1,l
				g=g+wk1(k,j)*wk1(i,k)
  14				continue
			endif
			wk2(j)=g/h
			f=f+wk2(j)*wk1(i,j)
  15			continue
		hh=f/(h+h)
			do 17 j=1,l
			f=wk1(i,j)
			g=wk2(j)-hh*f
			wk2(j)=g
				do 16 k=1,j
				wk1(j,k)=wk1(j,k)-f*wk2(k)-g*wk1(i,k)
  16				continue
  17			continue
		endif
	else
		wk2(i)=wk1(i,l)
	endif
	wk3(i)=h
  18	continue
	if(iptvec.eq.1) wk3(1)=0.d0
	wk2(1)=0.d0
	do 23 i=1,n
	if(iptvec.eq.1) then
		l=i-1
		if(wk3(i).ne.0.d0) then
			do 21 j=1,l
			g=0.d0
				do 19 k=1,l
				g=g+wk1(i,k)*wk1(k,j)
  19				continue
				do 20 k=1,l
				wk1(k,j)=wk1(k,j)-g*wk1(k,i)
  20				continue
  21			continue
		endif
	endif
	wk3(i)=wk1(i,i)
	if(iptvec.eq.1) then
		wk1(i,i)=1.d0
		if(l.ge.1) then
			do 22 j=1,l
			wk1(i,j)=0.d0
			wk1(j,i)=0.d0
  22			continue
		endif
	endif
  23	continue
c
c
c
	do 31 i=2,n
  31	wk2(i-1)=wk2(i)
	wk2(n)=0.d0
	do 35 l=1,n
		iter=0
   1		do 32 m=l,n-1
			dd=dabs(wk3(m))+dabs(wk3(m+1))
			if(dabs(wk2(m))+dd.eq.dd) go to 2
  32		continue
   		m=n
   2		if(m.ne.l) then
			if(iter.eq.maxit) then
				ier=1
				go to 99
 		endif
			iter=iter+1
			g=(wk3(l+1)-wk3(l))/(2.*wk2(l))
			r=dsqrt(g**2+1.d0)
			g=wk3(m)-wk3(l)+wk2(l)/(g+dsign(r,g))
			s=1.d0
			c=1.d0
			p=0.d0
			do 34 i=m-1,l,-1
				f=s*wk2(i)
				b=c*wk2(i)
				if(dabs(f).ge.dabs(g)) then
					c=g/f
	 				r=dsqrt(c**2+1.d0)
					wk2(i+1)=f*r
					s=1.d0/r
					c=c*s
				else
					s=f/g
					r=dsqrt(s**2+1.d0)
					wk2(i+1)=g*r
					c=1.d0/r
 					s=s*c
				endif
				g=wk3(i+1)-p
				r=(wk3(i)-g)*s+2.d0*c*b
				p=s*r
				wk3(i+1)=g+p
				g=c*r-b
				if(iptvec.eq.1) then
					do 33 k=1,n
					f=wk1(k,i+1)
					wk1(k,i+1)=s*wk1(k,i)+c*f
					wk1(k,i)=c*wk1(k,i)-s*f
  33					continue
				endif
  34			continue
		wk3(l)=wk3(l)-p
		wk2(l)=g
		wk2(m)=0.d0
		go to 1
		endif
  35	continue
c
c
c
	do 40 i=1,n
	do 40 j=1,n
  40	a(i,j)=wk1(i,j)
	do 45 i=1,n
  45	vals(i)=wk3(i)
c
c   sort:
c
	do 53 i=1,n-1
 	k=i
		cp=vals(i)
		do 51 j=i+1,n
			if(vals(j).ge.cp) then
				k=j
				cp=vals(j)
			endif
  51		continue
		if(k.ne.i) then
			vals(k)=vals(i)
			vals(i)=cp
			if(iptvec.eq.1) then
				do 52 j=1,n
					cp=a(j,i)
					a(j,i)=a(j,k)
					a(j,k)=cp
  52				continue
			endif
		endif
  53	continue
	ier=0
c
c
  99	return
	end
