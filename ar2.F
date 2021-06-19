c&AR2
c&ARCR2
	subroutine arcr2(args,nargs)
c****************************************************************
c
c   Subroutine to handle the command
c
c   ARCORR2(a,p,sigma,m,rx0,ry0,rhoxy0,rhox,rhoy,rhoxy,rhoyx,ier)
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),np)
	if(np.lt.1.or.na.lt.4*np) then
		call error(args,2,2)
		go to 99
 	endif
	if(ickl(args(3),nps,ns).eq.1) go to 99
	if(ns.lt.4) then
		call error(args,3,2)
		go to 99
	endif
	call ckint(args(4),m)
	if(m.lt.0.or.m.gt.500) then
		call error(args,4,2)
		go to 99
	endif
c
c
	nws=5*(m+1)+9*np+21
	if(ickse(nws).eq.1) go to 99
	n1=4*np+1
	n2=n1+4
	n3=n2+4*np
	n4=n3+4*(m+1)
	n5=n4+4
	n6=n5+4
	n7=n6+4
	n8=n7+4
	n9=n8+np+1
	mp1=m+1
	n10=n9+mp1
	n11=n10+mp1
	n12=n11+mp1
	call movxy(wk,array(nstart(npa)),16*np)
	call movxy(wk(n1),array(nstart(nps)),16)
c
c
	call arcvmi(wk,wk(n1),np,m,100,.001,wk(n2),wk(n3),wk(n4),ier,
     1	wk(n5),wk(n6),wk(n7),wk(n8),rx0,ry0,rhoxy0,wk(n9),wk(n10),
     1	wk(n11),wk(n12))
c
c
	call ckaddi(args(12),ier,iref)
	if(iref.eq.1) go to 99
	if(ier.eq.1) then
		write(*,1)
   1		format(' Process not stationary')
		imflag=1
		go to 99
	endif
	call ckaddr(args(5),rx0,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(6),ry0,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(7),rhoxy0,iref)
	if(iref.eq.1.or.m.eq.0) go to 99
	lab='Autocorrelations of First Series'
	call ckadda(args(8),m,lab,n9,iref)
	if(iref.eq.1) go to 99
	lab='Autocorrelations of Second Series'
	call ckadda(args(9),m,lab,n10,iref)
	if(iref.eq.1) go to 99
	lab='Cross-correlations of Positive Lags'
	call ckadda(args(10),m,lab,n11,iref)
	if(iref.eq.1) go to 99
	lab='Cross-correlations of Negative Lags'
	call ckadda(args(11),m,lab,n12,iref)
c
c
  99	continue
	return
	end
c&ARCVMI
      SUBROUTINE ARCVMI(A,SIG,NP,M,MAXIT,DEL,B,R,R0,ier,wk,wk1,wk2,c1,
     1	rx0,ry0,rhoxy0,rhox,rhoy,rhoxy,rhoyx)
C*******************************************************
C
C   SUBROUTINE TO CALCULATE THE 1ST M+1 AUTOCOVARIANCES
C   (R(0),...,R(M)) OF AN NPTH ORDER AUTOREGRESSIVE
C   PROCESS OF DIMENSION 2.
C
C*******************************************************
C
      DIMENSION A(2,2,1),SIG(2,2),R(2,2,1),R0(2,2),B(2,2,1),
     1	WK(2,2),WK1(2,2),WK2(2,2),c1(1),rhox(1),rhoy(1),rhoxy(1),
     1	rhoyx(1)
C
      CALL ARMAM(A,NP,NP,B,wk)
      CALL MACVM(B,SIG,NP,R,R0)
	npp1=np+1
c	write(*,3) r0(1,1),r0(2,2),r0(1,2),r0(2,1)
c  3	format(10x,4e15.7)
C
C
	do 200 it=1,maxit
C
C   FIND NEW B :
C
      DO 20 IM=1,2
      DO 20 IN=1,2
         C=0.
         DO 25 IJ=1,NP
         NP1MIJ=NP+1-IJ
            DO 30 IL=1,2
  30        C=C+A(IM,IL,IJ)*B(IL,IN,NP1MIJ)
  25     CONTINUE
  20  WK(IM,IN)=-C
C
C   CALCULATE SIG*WKTRANSPOSE :
C
      DO 35 I=1,2
      DO 35 J=1,2
         C=0.
         DO 40 K=1,2
  40     C=C+SIG(I,K)*WK(J,K)
  35  WK1(I,J)=C
C
C   CALCULATE NEW R'S AND CHECK CONVERGENCE :
C
      DO 45 I=1,2
      DO 45 J=1,2
         C=0.
         DO 50 K=1,2
  50     C=C+WK(I,K)*WK1(K,J)
      WK2(I,J)=C
  45  R0(I,J)=R0(I,J)+C
         CALL FNORM(WK2,C1(1))
      DO 55 IV=1,NP
      NN=NP+1-IV
         DO 60 I=1,2
         DO 60 J=1,2
            C=0.
            DO 65 K=1,2
  65        C=C+B(I,K,NN)*WK1(K,J)
         WK2(I,J)=C
  60     R(I,J,IV)=R(I,J,IV)+C
  55  CALL FNORM(WK2,C1(IV+1))
c	write(*,2) c1(1),c1(2)
c  2	format(1x,2e15.7)
c
c
      CALL MAX(C1,NPP1,ACC,IND)
	if(abs(acc).gt.1.e10) then
		ier=1
		go to 99
	endif
	if(abs(acc).lt.del) go to 250
C
C   SHIFT B'S AND GO TO NEXT ITERATION :
C
	npm1=np-1
	if(np.gt.1) then
		do 110 iv=1,np-1
		ivp1=iv+1
		do 110 i=1,2
		do 110 j=1,2
 110		b(i,j,iv)=b(i,j,ivp1)
	endif
		do 120 i=1,2
		do 120 j=1,2
 120		b(i,j,np)=wk(i,j)
 200	continue
C
C   NONCONVERGENCE :
C
	ier=1
	go to 99
 250	continue
	ier=0
c
c   find the rest of R's if necessary:
c
	if(m.le.np) go to 99
	do 300 iv=npp1,m
	do 270 ll=1,2
	do 270 mm=1,2
		c=0.0
		do 260 j=1,np
		ivj=iv-j
 260		c=c-r(ll,1,ivj)*a(mm,1,j)-r(ll,2,ivj)*a(mm,2,j)
 270	r(ll,mm,iv)=c
 300	continue
	rx0=r0(1,1)
	ry0=r0(2,2)
	dv=sqrt(rx0*ry0)
	rhoxy0=r0(1,2)/dv
	if(m.eq.0) go to 99
	do 310 iv=1,m
	rhox(iv)=r(1,1,iv)/rx0
	rhoy(iv)=r(2,2,iv)/ry0
	rhoxy(iv)=r(1,2,iv)/dv
 310	rhoyx(iv)=r(2,1,iv)/dv
  99	continue
	return
	end
c&FNORM
	subroutine fnorm(a,c)
c****************************************************************
c
c
c****************************************************************
	dimension a(2,2)
	c=sqrt(a(1,1)**2+a(2,2)**2+a(1,2)**2+a(2,1)**2)
	return
	end
c&MACVM
      SUBROUTINE MACVM(B,SIG,NQ,R,R0)
C********************************************************
C
C   SUBROUTINE TO CALCULATE THE 1ST NQ+1 AUTOCOVARIANCES
C   R0,R(1),...,R(NQ) OF A MOVING AVERAGE PROCESS OF ORDER NQ
C
C********************************************************
C
      DIMENSION B(2,2,1),SIG(2,2),R(2,2,1),r0(2,2)
      NQP1=NQ+1
      DO 10 IVP1=1,NQP1
      IV=IVP1-1
      NQMVP1=NQ-IV+1
         DO 15 IM=1,2
         DO 15 IN=1,2
            C=0.
            DO 20 IJP1=1,NQMVP1
            IJ=IJP1-1
               DO 25 IL=1,2
                  DO 30 IR=1,2
                  C1=0.
                  C3=0.
                  IF(IJ.EQ.0) GO TO 35
                  C1=B(IM,IR,IJ)
                  GO TO 40
  35              IF(IM.EQ.IR) C1=1.
  40              IF((IJ+IV).EQ.0) GO TO 45
                  C3=B(IN,IL,IJ+IV)
                  GO TO 30
  45              IF(IN.EQ.IL) C3=1.
  30              C=C+C1*SIG(IR,IL)*C3
  25           CONTINUE
  20        CONTINUE
         IF(IV.EQ.0) GO TO 50
         R(IM,IN,IV)=C
         GO TO 15
  50     R0(IM,IN)=C
  15     CONTINUE
  10  CONTINUE
C
      RETURN
      END
c&ARMAM
      SUBROUTINE ARMAM(A,NP,NQ,B,wk1)
C********************************************************
C
C   SUBROUTINE TO FIND THE COEFFICIENT MATRICES B(1),...,B(NQ)
C   OF A MOVING AVERAGE REPRESENTATION OF ORDER NQ (OF
C   DIMENSION 2) OF AN AUTOREGRESSIVE PROCESS OF ORDER
C   NP WITH COEFFICIENT MATRICES A(1),...,A(NP).
C
C   INPUT :
C           NP,NQ
C           A : THREE DIMENSIONAL ARRAY CONTAINING
C           AUTOREGRESSIVE COEFFICIENTS BY ROW,COLUMN, LAG
C
C   OUTPUT :
C           B : THREE DIMENSIONAL ARRAY CONTAINING
C           MOVING AVERAGE COEFFICIENTS BY ROW,COL,LAG
C
C   SUBROUTINES CALLED : NONE
C
C********************************************************
C
      DIMENSION A(2,2,np),B(2,2,nq),wk1(2,2)
C
      DO 10 IK=1,NQ
         JJ=MIN0(NP,IK)
	 call movct(wk1,16,char(0))
         DO 30 IJ=1,JJ
         IF(IJ.EQ.IK) GO TO 40
         II=IK-IJ
            DO 35 I=1,2
            DO 35 J=1,2
               C=0.
               DO 36 L=1,2
  36           C=C+A(I,L,IJ)*B(L,J,II)
  35        WK1(I,J)=WK1(I,J)-C
            GO TO 30
  40     DO 41 I=1,2
         DO 41 J=1,2
  41     WK1(I,J)=WK1(I,J)-A(I,J,IJ)
  30     CONTINUE
         DO 45 I=1,2
         DO 45 J=1,2
  45     B(I,J,IK)=WK1(I,J)
  10  CONTINUE
C
      RETURN
      END
c&ARDT2
	subroutine ardt2(args,nargs,vname)
c********************************************************************
c
c   Subroutine to handle the command
c
c   X=ARDT2(A,p,SIGMA,seed,n,ier)
c
c   or
c
c   X=ARDT2(A,p,n,f,ier)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),np)
	if(np.lt.1.or.na.lt.4*np) then
		call error(args,2,2)
		go to 99
 	endif
c
c   handle difference equation case:
c
	if(nargs.eq.5) then
		if(icki(args(3),3,n,2,np).eq.1) go to 99
		if(ickl(args(4),npf,nf).eq.1) go to 99
		if(nf.lt.2*n) then
			call error(args,4,2)
			go to 99
		endif
		if(4.*float(n+np).gt.32767.) then
			if(ickse(maxnum+1).eq.1) go to 99
		endif
		nws=4*(n+np)
		if(ickse(nws).eq.1) go to 99
		n1=2*n+1
		n2=n1+4*np
		call movxy(wk,array(nstart(npf)),8*n)
		call movxy(wk(n1),array(nstart(npa)),16*np)
		call deq2(wk,wk(n1),np,n,wk(n2),ier)
		call ckaddi(args(5),ier,iref)
		if(ier.eq.1) then
			write(*,11)
			imflag=1
			go to 99
		endif
		if(iref.eq.1) go to 99
		lab='Values for Bivariate Difference Equation'
		call ckadda(vname,2*n,lab,n2,iref)
		go to 99
	endif
c
c   handle non-difference equation case:
c
	if(ickl(args(3),nps,ns).eq.1) go to 99
	if(ns.lt.4) then
		call error(args,3,2)
		go to 99
	endif
	if(ickr(args(4),4,rseed,1,0.).eq.1) go to 99
	if(rseed.ne.0.0) i4seed=rseed
	if(icki(args(5),5,n,2,np).eq.1) go to 99
		if(4.*float(n+np+1).gt.32767.) then
			if(ickse(maxnum+1).eq.1) go to 99
		endif
		nws=4*(n+np+1)
		if(ickse(nws).eq.1) go to 99
c
c
	n1=4*np+1
	n2=n1+4
	n3=n2+4*n
	call movxy(wk,array(nstart(npa)),16*np)
	call movxy(wk(n1),array(nstart(nps)),16)
	call ardt2c(wk,np,wk(n1),wk(n2),i4seed,n,wk(n3),ier)
	call ckaddi(args(6),ier,iref)
	if(ier.eq.1) then
		write(*,11)
  11		format(' Singular SIGMA or Nonstationary Process')
		imflag=1
		go to 99
	endif
	if(iref.eq.1) go to 99
	lab='Bivariate Autoregressive Data'
	call ckadda(vname,2*n,lab,n3,iref)
c
c
  99	continue
	return
	end
c&DEQ2
	subroutine deq2(f,a,np,n,x,ier)
c*******************************************************************
c
c   x(i)=f(i), i=1,...,np
c
c   x(np+i)=f(np+i)-Sum(j=1,np) a(j)x(np+i-j), i=1,...,n-np
c
c*******************************************************************
c
	dimension f(2,n),a(2,2,np),x(2,n)
c
c
	ier=0
	do 10 i=1,np
	do 10 j=1,2
  10	x(j,i)=f(j,i)
c
	do 50 i=1,n-np
	nppi=np+i
	do 20 k=1,2
	c=f(k,nppi)
	do 30 j=1,np
	do 30 l=1,2
  30	c=c-a(k,l,j)*x(l,nppi-j)
	if(abs(c).gt.1.e10) then
		ier=1
		go to 99
	endif
  20 	x(k,nppi)=c
  50	continue
  99	continue
c
	return
	end
c&ARDT2C
	subroutine ardt2c(a,np,sigma,e,iseed,n,x,ier)
c*******************************************************************
c
c   Subroutine to generate bivariate AR data
c
c*******************************************************************
c
	dimension a(2,2,np),sigma(2,2),e(2,1),x(2,1)
	integer*4 iseed
c
c
	ier=0
	det=sigma(1,1)*sigma(2,2)-sigma(1,2)**2
	if(sigma(1,1).le.0.0.or.det.le.0.0) then
		ier=1
		go to 99
	endif
	d=sigma(1,2)/sigma(1,1)
	c1=sqrt(sigma(1,1))
	c2=c1*d
	c3=sqrt(sigma(2,2)-d*d*sigma(1,1))
	call whtsim(4*n,iseed,e)
	n2=2*n
	do 20 i=1,n2
	e(2,i)=c3*e(2,i)+c2*e(1,i)
  20	e(1,i)=c1*e(1,i)
c
c
	call movct(x,8*np,char(0))
c
c
	do 50 i=1,n2-np
	nppi=np+i
	do 25 k=1,2
	c=e(k,nppi)
	do 35 j=1,np
	do 35 l=1,2
  35	c=c-a(k,l,j)*x(l,nppi-j)
	if(abs(c).gt.1.e10) then
		ier=1
		go to 99
	endif
  25 	x(k,nppi)=c
  50	continue
c
c
	do 70 i=1,n
	npi=n+i
	do 70 j=1,2
  70	x(j,i)=x(j,npi)
c
c
  99	continue
	return
	end
c&CORAR2
	subroutine corar2(args,nargs,vname)
c********************************************************************
c
c   Subroutine to handle the command
c
c   A=CORRAR2(rx0,ry0,rhoxy0,rhox,rhoy,rhoxy,rhoyx,m,iopto,n,
c                 nord,sigma,cat,ier)
c
c   or
c
c   A=CORRAR2(rx0,ry0,rhoxy0,rhox,rhoy,rhoxy,rhoyx,p,sigma,ier)
c
c********************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr
c
c
	if(ickr(args(1),1,rx0,2,0.).eq.1) go to 99
	if(ickr(args(2),2,ry0,2,0.).eq.1) go to 99
	call ckreal(args(3),rhoxy0)
	if(rhoxy0.eq.999..or.abs(rhoxy0).ge.1.) then
		call error(args,3,2)
		go to 99
	endif
	if(ickl(args(4),npx,nx).eq.1) go to 99
	if(ickl(args(5),npy,ny).eq.1) go to 99
	if(ickl(args(6),npxy,nxy).eq.1) go to 99
	if(ickl(args(7),npyx,nyx).eq.1) go to 99
	nn=min0(nx,ny,nxy,nyx)
	if(icki(args(8),8,m,0,0).eq.1) go to 99
	if(m.lt.1.or.m.gt.nn.or.m.gt.200) then
		call error(args,8,2)
		go to 99
        endif
c
c   trap nargs=10 case:
c
	if(nargs.eq.10) then
		iopto=1
		n=100
		go to 5
	endif
	call ckint(args(9),iopto)
	if(iopto.lt.1.or.iopto.gt.2) then
		call error(args,9,2)
		go to 99
	endif
	if(icki(args(10),10,n,2,m).eq.1) go to 99
c
c
   5	continue
	sxy=sqrt(rx0*ry0)
	wk(1)=rx0
	wk(2)=rhoxy0*sxy
	wk(3)=wk(2)
	wk(4)=ry0
	do 10 i=1,m
	n1=4*i
		wk(n1+1)=array(nstart(npx)+i-1)*rx0
		wk(n1+2)=array(nstart(npyx)+i-1)*sxy
		wk(n1+3)=array(nstart(npxy)+i-1)*sxy
		wk(n1+4)=array(nstart(npy)+i-1)*ry0
  10	continue
c
c
	m4=4*m
	n1=5
	n2=1
	n3=4*(m+1)+1
	n4=n3+m4
	n5=n4+m4
	n6=n5+m4
	n7=n6+m4
	n8=n7+m4
	n9=n8+m4
	n10=n9+m4
	n11=n10+m
	n12=n11+m4
	n13=n12+m4
c
c
	call cvarm2(wk(n1),wk(n2),m,iopto,n,wk(n3),wk(n4),wk(n5),
     1	wk(n6),wk(n7),wk(n8),wk(n9),wk(n10),nord,wk(n11),wk(n12),
     1	wk(n13),ier)
c
c
	call ckaddi(args(nargs),ier,iref)
	if(ier.eq.1) then
		write(*,91)
  91		format(' Singular Matrix in CORRAR2')
		imflag=1
		go to 99
	endif
	if(iref.eq.1) go to 99
c
c   trap nargs=10 case:
c
	if(nargs.eq.10) then
		lab='Autoregressive Coefficient Matrices'
		call ckadda(vname,4*m,lab,n12,iref)
		if(iref.eq.1) go to 99
		lab='Error Covariance Matrix'
		call ckadda(args(9),4,lab,n11+4*(m-1),iref)
		go to 99
	endif
	call ckaddi(args(11),nord,iref)
	if(iref.eq.1.or.nord.eq.0) go to 99
	lab='Autoregressive Coefficient Matrices'
	call ckadda(vname,4*nord,lab,n12,iref)
	if(iref.eq.1) go to 99
	lab='Autoregressive Error Variance Matrix'
	call ckadda(args(12),4,lab,n11+4*(nord-1),iref)
	if(iref.eq.1) go to 99
	lab='Bivariate CAT Criterion'
	call ckadda(args(13),m,lab,n13,iref)
c
c
  99	continue
	return
	end
c&CVARM2
	subroutine cvarm2(r,r0,m,iopto,n,a1,a2,ab1,ab2,del,delb,vb,tr,
     1	nord,v,a,cat,ier)
c********************************************************************
c
c   This is a direct application of Whittle's (1963 Biometrika)
c   multiple Levinson to d=2.  It is specially designed for TIMESLAB.
c
c
c********************************************************************
c
	dimension r(2,2,m),r0(2,2),a1(2,2,m),a2(2,2,m),ab1(2,2,m),
     1	ab2(2,2,m),del(2,2,m),delb(2,2,m),vb(2,2,m),tr(m),v(2,2,m),
     1	a(2,2,m),cat(m)
c
c
c  Initialize:
c
	ier=0
	nord=m
	call movxy(v,r0,16)
	call movxy(vb,r0,16)
	call movxy(del,r,16)
	call movxy(delb,del,16)
	do 5 i=1,2
	do 5 j=1,i
	c=del(i,j,1)
	del(i,j,1)=del(j,i,1)
   5	del(j,i,1)=c
	call solv2(vb,del,a1,ier)
	if(ier.eq.1) go to 99
	call solv2(v,delb,ab1,ier)
      	if(ier.eq.1) go to 99
	if(iopto.eq.2) then
		cd=r0(1,1)*r0(2,2)-r0(1,2)*r0(1,2)
		if(cd.le.0.) then
			ier=1
			go to 99
		endif
		catmin=-(r0(1,1)+r0(2,2))/cd
		nord=0
		on=n
		on2=on*on
	endif
c
c
c  Start Iteration:
c
c
	do 100 np=1,m
	if(np.eq.1) go to 40
	npm1=np-1
c
c   partials:
c
	call solv2(vb(1,1,npm1),del(1,1,npm1),a2(1,1,np),ier)
	if(ier.eq.1) go to 99
	call solv2(v(1,1,npm1),delb(1,1,npm1),ab2(1,1,np),ier)
	if(ier.eq.1) go to 99
c
c   the rest:
c
		do 30 j=1,npm1
		np1=npm1-j+1
			do 20 k=1,2
			do 20 l=1,2
			c1=a1(k,l,j)
			c2=ab1(k,l,j)
				do 10 i=1,2
 				c1=c1+a2(k,i,np)*ab1(i,l,np1)
  10				c2=c2+ab2(k,i,np)*a1(i,l,np1)
			a2(k,l,j)=c1
  20			ab2(k,l,j)=c2
  30		continue
	nmv=16*np
	call movxy(a1,a2,nmv)
	call movxy(ab1,ab2,nmv)
c
c   V's:
c
  40	continue
		do 60 k=1,2
		do 60 l=1,2
		c1=r0(k,l)
		c2=c1
			do 50 j=1,np
			do 50 i=1,2
			c1=c1+a1(k,i,j)*r(i,l,j)
  50			c2=c2+ab1(k,i,j)*r(l,i,j)
  		v(k,l,np)=c1
  60		vb(k,l,np)=c2
c
c   Do CAT if iopto=2:
c
	if(iopto.eq.2) then
		cd=v(1,1,np)*v(2,2,np)-v(1,2,np)*v(1,2,np)
		if(cd.le.0.) then
			ier=1
			go to 99
		endif
		tr(np)=(v(1,1,np)+v(2,2,np))/cd
		onp=np
		c2=((on-onp*2)*tr(np))/on
		c1=0.
			do 65 j=1,np
  65			c1=c1+(on-2.*float(j))*tr(j)
		c1=(2.*c1)/on2
		cat(np)=c1-c2
		if(cat(np).lt.catmin) then
			catmin=cat(np)
			nord=np
			call movxy(a,a1,16*np)
		endif
	endif
	if(np.eq.m) go to 100
c
c   Del's:
c
		do 80 k=1,2
		do 80 l=1,2
		c1=r(l,k,np+1)
		c2=r(k,l,np+1)
			do 70 j=1,np
			nn=np-j+1
			do 70 i=1,2
			c1=c1+a1(k,i,j)*r(l,i,nn)
  70			c2=c2+ab1(k,i,j)*r(i,l,nn)
		del(k,l,np)=c1
  80		delb(k,l,np)=c2
c
c
c
 100	continue
	if(iopto.eq.1) call movxy(a,a1,16*m)
c
c
  99	continue
	return
	end
c&SOLV2
	subroutine solv2(a,b,c,ier)
c*******************************************************************
c
c   Subroutine to find c=-b*a inverse.
c
c*******************************************************************
c
	dimension a(2,2),b(2,2),c(2,2),ai(2,2)
c
	ier=0
	da=a(1,1)*a(2,2)-a(1,2)*a(2,1)
	if(da.le.0.) then
		ier=1
		go to 99
	endif
	ai(1,1)=-a(2,2)/da
	ai(2,2)=-a(1,1)/da
	ai(1,2)=a(1,2)/da
	ai(2,1)=a(2,1)/da
	do 20 i=1,2
	do 20 j=1,2
		c1=0.
		do 10 k=1,2
  10 		c1=c1+b(i,k)*ai(k,j)
  20	c(i,j)=c1
c
c
  99	continue
	return
	end
c&ARSP2
	subroutine arsp2(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   ARSP2(a,nord,sigma,nfreqs,fxx,fyy,fxyr,fxyi)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickr
c
c
	if(icki(args(2),2,nord,2,0).eq.1) go to 99
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(na.lt.4*nord) then
		call error(args,2,2)
		go to 99
	endif
	if(ickl(args(3),npsg,nsg).eq.1) go to 99
	if(nsg.lt.4) then
		call error(args,3,2)
		go to 99
	endif
	call ckint(args(4),nfreqs)
	if(nfreqs.lt.2.or.nfreqs.gt.2048) then
		call error(args,4,2)
		go to 99
	endif
c
c
	nn=(nfreqs/2)+1
	frstw=0.
	delw=8.*atan(1.0)/float(nfreqs)
	call movxy(wk,array(nstart(npsg)),16)
	call movxy(wk(5),array(nstart(npa)),16*nord)
	n1=4*nord+5
	n2=n1+nn
	n3=n2+nn
	n4=n3+nn
c
	call arsp21(wk(5),wk,2,nord,nn,frstw,delw,wk(n1),wk(n2),
     1	wk(n3),wk(n4))
c
c
	lab='Autospectra of First Series'
	call ckadda(args(5),nn,lab,n1,iref)
	if(iref.eq.1) go to 99
	lab='Autospectra of Second Series'
	call ckadda(args(6),nn,lab,n4,iref)
	if(iref.eq.1) go to 99
	lab='Real Part Of Cross-Spectra'
	call ckadda(args(7),nn,lab,n2,iref)
	if(iref.eq.1) go to 99
	lab='Imaginary Part of Cross-Spectra'
	call ckadda(args(8),nn,lab,n3,iref)
c
c
  99	continue
	return
	end
c&ARSP21
      SUBROUTINE ARSP21(A,SIGMA,NDMAX,NORD,nn,frstw,delw,
     1	f11,f12r,f12i,f22)
C**********************************************************
C
C   SUBROUTINE TO CALCULATE THE 2 DIMENSIONAL SPECTRAL
C   DENSITY OF 2 DIMENSIONAL AUTOREGRESSIVE PROCESS WITH
C   PARAMETERS A(1),...,A(NORD) AND SIGMA AT THE
C   FREQUENCIES frstw+(k-1)*delw, k=1,...,nn.
C
C   INPUT :
C           NORD,NN
C           A(1),...,A(NORD) : STORED IN 3 DIMENSIONAL
C           ARRAY A BY ROW,COLUMN,INDEX
C           SIGMA
C           NDMAX : DIMENSION OF VARIOUS ARRAYS IN CALLING
C           PROGRAM (SEE DIMENSIONS BELOW)
C
C   OUTPUT :
C           F11,F22 : AUTOSPECTRA OF CHANNELS 1,2 RESPECTIVELY
C           F12R,F12I : CO AND QUADRATURE SPECTRA
C
C   SUBROUTINES CALLED : NONE
C
C**********************************************************
C
      DIMENSION A(NDMAX,NDMAX,NORD),SIGMA(NDMAX,NDMAX),
     1F11(NN),F12R(NN),F12I(NN),F22(NN),SI(2,2)
      COMPLEX C1,C2,C3,C5,C6,G(2,2),FFI(2,2),ARG,C2PI
C
C   FIND SIGMA INVERSE :
C
      C=SIGMA(1,1)*SIGMA(2,2)-SIGMA(1,2)**2
      SI(1,1)=SIGMA(2,2)/C
      SI(2,2)=SIGMA(1,1)/C
      SI(1,2)=-SIGMA(1,2)/C
      SI(2,1)=SI(1,2)
C
C   FIND TRANSFORM OF A (I.E. G) :
C
      C5=CMPLX(0.0,0.0)
      C6=CMPLX(1.0,0.0)
C
      DO 15 K=1,NN
      FAC1=frstw+float(k-1)*delw
         DO 20 L1=1,2
         DO 20 L2=1,2
            C1=C5
            IF(L1.EQ.L2) C1=C6
      IF(NORD.EQ.0) GO TO 20
            DO 25 J=1,NORD
            OJ=J
            ARG=CMPLX(0.0,FAC1*OJ)
  25        C1=C1+CMPLX(A(L1,L2,J),0.0)*CEXP(ARG)
  20     G(L1,L2)=C1
C
C   FIND FINV = GSTAR*SIGMAINVERSE*G :
C
         DO 30 I=1,2
         DO 30 J=1,I
            C1=C5
            DO 35 L=1,2
               C2=C5
               DO 40 IR=1,2
  40           C2=C2+CONJG(G(IR,I))*CMPLX(SI(IR,L),0.0)
  35        C1=C1+C2*G(L,J)
         FFI(I,J)=C1
         FFI(J,I)=CONJG(C1)
  30     CONTINUE
C
C   FIND F :
C
      C1=FFI(1,1)*FFI(2,2)-FFI(1,2)*FFI(2,1)
      C3=-FFI(1,2)/C1
      F11(K)=REAL(FFI(2,2)/C1)
      F22(K)=REAL(FFI(1,1)/C1)
      F12R(K)=REAL(C3)
      F12I(K)=AIMAG(C3)
C
C
C
  15  CONTINUE
C
C
C
      RETURN
      END
c&MCHOL
      SUBROUTINE MCHOL(R,MDIM,M,AL,D,IER)
C*******************************************************
C
C   SUBROUTINE TO OBTAIN THE MODIFIED CHOLESKY DECOMPOSITION
C   OF A SYMMETRIC POSITIVE DEFINITE
C   MXM MATRIX R=L*D*(LTRANSPOSE).
C
C   INPUT :
C           MDIM : DIMENSION OF R AND AL IN MAIN PROGRAM
C           M,R(MXM MATIX)
C
C   OUTPUT :
C           AL : L IS IN LOWER TRIANGULAR PORTION OF AL,
C           LTRANSPOSE IS IN UPPER TRIANGULAR POSITION.
C           D : VECTOR OF DIAGONAL ELEMENTS
C           IER : ERROR INDICATOR (IER=1 IS ERROR RETURN,
C           IER = 0 IS NORMAL RETURN
C
C   SUBROUTINES CALLED : NONE
C
C*********************************************************
C
      DIMENSION R(MDIM,MDIM),AL(MDIM,MDIM),D(M)
      IER=0
C
      DO 30 K=1,M
      AL(1,K)=R(1,K)
  30  AL(K,1)=AL(1,K)
C
      DO 2 K=1,M
C
      D(K)=R(K,K)
      AL(K,K)=1.
      IF(K.EQ.1) GO TO 10
         KM1=K-1
         DO 3 J=1,KM1
   3     D(K)=D(K)-AL(J,K)*AL(K,J)
  10  IF(D(K).LT.1.E-20) GO TO 20
      IF(K.EQ.1) GO TO 15
      IF(K.EQ.M) GO TO 6
      GO TO 16
  15  DO 17 J=1,M
  17  AL(J,1)=AL(J,1)/D(1)
      AL(1,1)=1.
      GO TO 2
  16     KP1=K+1
         DO 4 J=KP1,M
         C=R(K,J)
            DO 5 I=1,KM1
   5        C=C-AL(I,K)*AL(J,I)
         AL(K,J)=C
   4     AL(J,K)=AL(K,J)/D(K)
      GO TO 2
  20     DO 11 I=1,M
         AL(K,I)=0.
  11     AL(I,K)=0.
      IER=1
      D(K)=0.
      WRITE(*,12) K
  12  FORMAT(10X,'COLUMN',I4,'IS DEPENDENT ON PREVIOUS ONES')
C
   2  CONTINUE
C
   6  DO 7 I=1,M
         DO 7 J=I,M
   7  AL(I,J)=AL(J,I)
C
      RETURN
      END
