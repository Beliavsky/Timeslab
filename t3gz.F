c&T3GZ
c&TRANS
        subroutine trans(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to process the command y=TRANS(x,n,m)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(icki(args(2),2,n,1,1).eq.1) go to 99
	if(icki(args(3),3,m,5,nx/n).eq.1) go to 99
        ntm=n*m
c
        ns1=nstart(npx)-1
        do 10 j=1,m
        jj=(j-1)*n
        do 10 i=1,n
        ji=jj+i
  10    wk((i-1)*m+j)=array(ns1+ji)
c
        lab='Transpose of matrix'
        call ckadda(vname,ntm,lab,1,iref)
c
c
  99    continue
        return
        end
c&INPOLY
        subroutine inpoly(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to invert a polynomial.
c
c   beta=INVPOLY(alpha,np,nq)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
        call ckint(args(2),np)
        if(np.lt.1.or.np.gt.na) then
	        call error(args,2,2)
        	go to 99
	endif
	if(icki(args(3),3,nq,1,1).eq.1) go to 99
	if(ickse(nq).eq.1) go to 99
c
c
        nsa=nstart(npa)-1
        wk(1)=-array(nsa+1)
        if(nq.eq.1) go to 60
        do 50 i=2,nq
        c=0.0
        if(i.le.np) c=-array(nsa+i)
        m=min0(np,i-1)
        do 45 j=1,m
  45    c=c-array(nsa+j)*wk(i-j)
  50    wk(i)=c
  60    continue
c
c
        lab='reciprocal of polynomial'
        call ckadda(vname,nq,lab,1,iref)
c
  99    continue
        return
        end
c&LOGE
        subroutine loge(args,nargs,vname)
c****************************************************************
c
c   Subroutine to form the loge of an array.
c
c   y=loge(x,n)   or y=loge(x)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,i10*1
	integer*2 ickl,ickr,ickse
c
c
	if(nargs.eq.1) then
                call cklist(args(1),xlist,numarr,np)
                if(np.gt.0) then
                        n=nend(np)-nstart(np)+1
                        go to 10
                endif 
		call ckreal(args(1),x)
		if(x.eq.999..or.x.le.0.) then
			call error(args,1,2)
			go to 99
		endif
		c=alog(x)
		call ckaddr(vname,c,iref)
		go to 99
	endif
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),n)
        if(n.le.0.or.n.gt.nn) then
		call error(args,2,2)
		go to 99
	endif
  10	if(ickse(n).eq.1) go to 99
c
c
        n1=nstart(np)-1
        do 30 i=1,n
        c=array(n1+i)
	if(c.le.0.0) then
		write(*,25)
  25		format(' Nonpositive Value in LOGE')
		imflag=1
		go to 99
	endif
        wk(i)=alog(c)
  30    continue
c
c
        lab='logged array'
        call ckadda(vname,n,lab,1,iref)
  99    continue
        return
        end
c&MTPOLY
        subroutine mtpoly(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to multiply two polynomials.
c
c   gam=MULTPOLY(alpha,beta,np,nq)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(ickl(args(2),npb,nb).eq.1) go to 99
        call ckint(args(3),np)
        call ckint(args(4),nq)
	ipte=0
	if(np.lt.1.or.np.gt.na) ipte=3
	if(nq.lt.1.or.nq.gt.nb) ipte=4
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
        nppnq=np+nq
	call movct(wk,4*nppnq,char(0))
        nsa=nstart(npa)-1
        nsb=nstart(npb)-1
        do 65 i=1,np
  65    wk(i)=wk(i)+array(nsa+i)
        do 70 i=1,nq
  70    wk(i)=wk(i)+array(nsb+i)
        do 75 i=1,np
        ii=nsa+i
        do 75 j=1,nq
  75    wk(i+j)=wk(i+j)+array(ii)*array(nsb+j)
c
c   store :
c
        lab='product of two polynomials'
        call ckadda(vname,nppnq,lab,1,iref)
  99    continue
        return
        end
c&PARTAL
        subroutine partal(args,nargs,vname)
c********************************************************************
c
c   subroutine to call partar.
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),nord)
        if(nord.le.0.or.nord.gt.nn) then
	        call error(args,2,2)
        	go to 99
	endif
c
c
	ns1=nstart(np)-1
	do 10 i=1,nord
  10	wk(nord+i)=-array(ns1+i)
        call partar(wk(nord+1),nord,wk)
c
c
        lab='autoregressive coefficients'
        call ckadda(vname,nord,lab,1,iref)
c
c
  99    continue
        return
        end
c&PARTAR
      SUBROUTINE PARTAR(PART,NP,ALPHA)
C***********************************************************
C
C   SUBROUTINE TO DETERMINE  AUTOREGRESSIVE COEFFICIENTS
C   ALPHA(1),...,ALPHA(NP) GIVEN PARTIAL AUTOCORRELATIONS
C   PART(1),...,PART(NP)
C
C   INPUT :
C         NP,PART(1),...,PART(NP)
C
C   OUTPUT :
C         ALPHA(1),...,ALPHA(NP)
C
C   SUBROUTINES CALLED : NONE
C
C***********************************************************
C
      DIMENSION PART(NP),ALPHA(NP)
C
C   NP=1 :
C
      ALPHA(1)=PART(1)
      IF(NP.EQ.1)GO TO 99
C
C   NP.GT.1 :
C
C
C   FIX ORDER INP IN RECURSION :
C
      DO 20 INP=2,NP
C
C   FIND ALPHA(1),...,ALPHA(INP-1) :
C
      ALPHA(INP)=PART(INP)
      C=ALPHA(INP)
         NL=INP/2
         DO 10 J=1,NL
         INPMJ=INP-J
         TEMP=ALPHA(J)+C*ALPHA(INPMJ)
         ALPHA(INPMJ)=ALPHA(INPMJ)+C*ALPHA(J)
   10    ALPHA(J)=TEMP
   20 CONTINUE
C
C
   99 CONTINUE
      RETURN
      END
c&QTEST
        subroutine qtest(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to do portmanteau  Q test.
c
c   pval=QTEST(rho,m,np,nq,n,q)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,i10*1
        double precision dpval,dq,pdf
	integer*2 ickl,icki
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),m)
        if(m.le.1.or.m.gt.nn) then
	        call error(args,2,2)
        	go to 99
	endif
	if(icki(args(3),3,norda,1,0).eq.1) go to 99
	if(icki(args(4),4,nordb,1,0).eq.1) go to 99
        if(m-norda-nordb.le.0) then
	        write(*,4)
   4    	format(' m-np-nq is less than or equal to 0')
		imflag=1
        	go to 99
	endif
	if(icki(args(5),5,n,2,m).eq.1) go to 99
c
c   calculate Q:
c
        n1=nstart(np)-1
        q=0.0
        do 50 i=1,m
  50    q=q+array(n1+i)*array(n1+i)/(float(n-i))
        q=float(n)*float(n+2)*q
        dq=q
        call chicdf(dq,m-norda-nordb,dpval,pdf)
        pval=dpval
        pval=1.-pval
c
c   put q and pval where they belong :
c
        call ckaddr(args(6),q,iref)
        if(iref.eq.1) go to 99
        call ckaddr(vname,pval,iref)
c
c
  99    continue
        return
        end
c&TRIG
        subroutine trig(args,nargs,vname,iopt)
c******************************************************************
c
c   Subroutine to form an array whose ith element is either
c  amp*cos(if iopt=1) or amp*sin(if iopt=2) of twopi*(i-1)/per,i=1,...,n
c   where n,amp, and per are specified by args(1),args(2),args(3).
c
c*******************************************************************
c
#include 'tslabc'
        character*15 args(nargs),vname
	integer*2 icki,ickr
c
c   handle y=cos(x) or y=sin(x):
c
	if(nargs.eq.1) then
                call cklist(args(1),xlist,numarr,npx)
                if(npx.gt.0) then
                        n=nend(npx)-nstart(npx)+1
                        ns1=nstart(npx)-1
                        if(iopt.eq.1) then
                              	do 1 i=1,n
   1				wk(i)=cos(array(ns1+i))
				lab='Cosine of an array'
			endif
                        if(iopt.eq.2) then
                              	do 2 i=1,n
   2				wk(i)=sin(array(ns1+i))
				lab='Sine of an array'
			endif
                        go to 5
                endif
		if(ickr(args(1),1,x,0,0.).eq.1) go to 99
		if(iopt.eq.1) c=cos(x)
		if(iopt.eq.2) c=sin(x)
		call ckaddr(vname,c,iref)
		go to 99
	endif
c
c   handle 3 argument case:
c
        if(nargs.eq.3) then
                if(icki(args(1),1,n,2,0).eq.1) go to 99
        	if(ickr(args(2),2,amp,0,0.).eq.1) go to 99
        	if(ickr(args(3),3,per,2,0.).eq.1) go to 99
c
                twopi=8.*atan(1.0)
                if(iopt.eq.1) then
                	do 70 i=1,n
  70            	wk(i)=amp*cos((twopi*float(i-1))/per)
			lab='Cosine Curve'
		endif
                if(iopt.eq.2) then
    		        do 90 i=1,n
  90            	wk(i)=amp*sin((twopi*float(i-1))/per)
			lab='Sine Curve'
		endif
                go to 5
        endif
c
c   handle y=cos(x,n) or y=sin(x,n):
c
	if(nargs.eq.2) then
		if(ickl(args(1),npx,nx).eq.1) go to 99
		call ckint(args(2),n)
		if(n.lt.1.or.n.gt.nx) then
			call error(args,2,2)
			go to 99
		endif
   		n1=nstart(npx)-1
		if(iopt.eq.1) then
			do 10 i=1,n
  10			wk(i)=cos(array(n1+i))
			lab='Cosine of an array'
		endif
		if(iopt.eq.2) then
			do 15 i=1,n
  15			wk(i)=sin(array(n1+i))
			lab='Sine of an array'
		endif
                go to 5
	endif

c
c
   5    call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&WHTSIM
      subroutine whtsim(n,dseed,x)
c*****************************************************************
c
c   subroutine to generate n iid N(0,1) rv's x(1),...,x(n). dseed
c   is an integer*4 seed on input and output.
c
c*******************************************************************
c
      dimension x(n)
        integer*4 dseed
c
c   if n is odd do last one separately :
c
        nn=n-1
        if(2*(n/2).eq.n) nn=n
        if(n.eq.nn) go to 20
   5    call unif2(dseed,1,r1)
        call unif2(dseed,1,r2)
        r1=2.*(r1-.5)
        r2=2.*(r2-.5)
        s2=r1*r1+r2*r2
        if(s2.gt.1.) go to 5
        y=sqrt(-2.*(alog(s2)/s2))
        x(n)=r1*y
        if(n.eq.1) go to 99
c
c   now do the rest :
c
  20    call unif2(dseed,nn,x)
        nnb2=nn/2
        do 30 i=1,nnb2
        i2=2*i
        i1=i2-1
  12    r1=2.*(x(i1)-.5)
        r2=2.*(x(i2)-.5)
        s2=r1*r1+r2*r2
        if(s2.le.1.) go to 15
        call unif2(dseed,2,x(i1))
        go to 12
  15    y=sqrt(-2.*(alog(s2)/s2))
        x(i1)=y*r1
        x(i2)=y*r2
  30    continue
c
c
  99    continue
        return
        end
c&WN
        subroutine wn(args,nargs,vname)
c******************************************************************
c
c   Subroutine to form a white noise series.
c
c******************************************************************
c
#include 'tslabc'
        character*15 args(nargs),vname
	integer*2 icki,ickr,ickse
c
c
	if(ickr(args(1),1,rseed,1,0.).eq.1) go to 99
	if(rseed.ne.0.) i4seed=rseed
	if(icki(args(2),2,n,1,1).eq.1) go to 99
	if(ickse(n).eq.1) go to 99
c
c    normal data :
c
        if(nargs.eq.3) then
	        call ckint(args(3),iptdis)
		if(iptdis.lt.1.or.iptdis.gt.8) then
			call error(args,3,2)
			go to 99
 	endif
        	if(iptdis.ne.1) go to 100
	endif
        call whtsim(n,i4seed,wk)
        lab='normal white noise series'
        go to 51
c
c   uniform :
c
 100    if(iptdis.eq.2) then
	        call unif2(i4seed,n,wk)
        	lab='uniform white noise series'
	        go to 51
	endif
c
c   exponential :
c
        if(iptdis.eq.3) then
	        call unif2(i4seed,n,wk)
        	do 105 i=1,n
 105    	wk(i)=-alog(wk(i))
        	lab='exponential white noise series'
        	go to 51
	endif
c
c   logistic :
c
        if(iptdis.eq.4) then
	        call unif2(i4seed,n,wk)
        	do 125 i=1,n
 125    	wk(i)=alog(wk(i)/(1.-wk(i)))
        	lab='logistic white noise series'
        	go to 51
	endif
c
c   cauchy :
c
        if(iptdis.eq.5) then
	        call unif2(i4seed,n,wk)
        	pi=4.*atan(1.0)
        	do 135 i=1,n
 135    	wk(i)=tan(pi*(wk(i)-.5))
        	lab='cauchy white noise series'
        	go to 51
	endif
c
c   extreme value :
c
        if(iptdis.eq.6) then
	        call unif2(i4seed,n,wk)
        	do 145 i=1,n
 145    	wk(i)=alog(-alog(1.-wk(i)))
        	lab='extreme value white noise series'
        	go to 51
	endif
c
c   lognormal :
c
        if(iptdis.eq.7) then
	        call whtsim(n,i4seed,wk)
        	do 155 i=1,n
 155    	wk(i)=exp(wk(i))
        	lab='lognormal white noise series'
        	go to 51
	endif
c
c   double exponential :
c
        if(iptdis.eq.8) then
	        call unif2(i4seed,n,wk)
        	do 165 i=1,n
        	if(wk(i).le..5) wk(i)=alog(2.*wk(i))
 165    	if(wk(i).gt..5) wk(i)=-alog(2.*(1.-wk(i)))
        	lab='double exponential white noise series'
        	go to 51
	endif
        write(*,101)
 101    format(' unimplemented random number generator')
	imflag=1
        go to 99
c
c   add array :
c
  51    call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&SORT
        subroutine sort(args,nargs,vname)
c******************************************************************
c
c   Subroutine to sort an array.
c
c   y=SORT(x,n)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
	if(nargs.eq.1) then
		n=nend(np)-nstart(np)+1
		go to 5
	endif
        call ckint(args(2),n)
        if(n.le.0.or.n.gt.na) then
	        call error(args,2,2)
        	go to 99
	endif
c
c
   5	call movxy(wk,array(nstart(np)),4*n)
        call qsort(n,wk)
c
c
        lab='sorted array'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
       subroutine qsort(n,x)
c********************************************************
c
c   Quicksort a real array x of length n.
c
c********************************************************
c
       dimension x(n),ll(20),lr(20)
c
c   ns is # of pieces left to partition
c
c   ll and lr are left and right borders of pieces
c
       ns=1
       ll(1)=1
       lr(1)=n
c
c   After splitting we come here. Stop when no pieces left.
c
   5   if(ns.eq.0) go to 99
c
       i=ll(ns)
       j=lr(ns)
       nl=j-i+1
c
c   Sort pieces of size 1 or 2:
c
       if(nl.le.2) then
              ns=ns-1
              if(nl.eq.2.and.x(i).gt.x(j)) call swap(x(i),x(j))
              go to 5
       endif
c
c   Sort first, ``middle'' and last elements:
c
       nm=(i+j)/2
       if(x(i).gt.x(nm)) call swap(x(i),x(nm))
       if(x(nm).gt.x(j)) call swap(x(nm),x(j))
       if(x(i).gt.x(nm)) call swap(x(i),x(nm))
c
c   If piece of size 3, it's now sorted
c
       if(nl.eq.3) then
              ns=ns-1
              go to 5
       endif
c
c   Put middle (target) into 1st and keep a copy for comparisons
c
       ax=x(nm)
       call swap(x(i),x(nm))
c
c   Look for 1st one from left (starting with 2nd)  > target 
c   and first from right < target. If no such pair we end 
c   up at 20.
c
       itemp=i
       i=i+1
c
  10   if(x(i).le.ax) then
              if(i.eq.j) go to 20
              i=i+1
              go to 10
       endif
c
  15   if(x(j).ge.ax) then
              if(i.eq.j) go to 20
              j=j-1
              go to 15
       endif
c
c   Can only get to here if we found a pair to swap
c
       call swap(x(i),x(j))
c
c   Only go back to look for another pair if there's a chance 
c   to get one
c
       if(j-i.le.1) then
              go to 20
       else
              i=i+1
              j=j-1
              go to 10
       endif
c
c   Get to here when no more pairs to swap. i might be 1 space 
c   too far to right 
c
  20   if(x(i).gt.ax) i=i-1
       call swap(x(itemp),x(i))
c
c   Now put shorter piece at end of ll and lr (this guarantees 
c   small number of pieces at any one time). A piece might be 
c   empty, but that is taken care of above.
c
       nleft=i-ll(ns)
       nright=lr(ns)-i
       ns=ns+1
       if(nleft.le.nright) then
              ll(ns)=ll(ns-1)
              ll(ns-1)=i+1
              lr(ns)=i-1
       else
              lr(ns)=lr(ns-1)
              lr(ns-1)=i-1
              ll(ns)=i+1
       endif
       go to 5
c
  99   continue
       return
       end
       subroutine swap(a,b)
c*****************************************************
c
c   Swap two reals.
c
c*****************************************************
c
       c=b
       b=a
       a=c
c
       return
       end
