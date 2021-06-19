c&ARSPPEAK
c&ARPEAK
	subroutine arpeak(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   ARSPPEAK(alpha,np,rvar,n,ier,peak,stderr[,start])
c
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,ickr,icki
c
c
	ier=0
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),np)
	if(np.lt.1.or.np.gt.na.or.np.gt.55) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(icki(args(4),4,n,2,0).eq.1) go to 99
	ipts=0
	if(nargs.eq.8) then
		ipts=1
		call ckreal(args(8),start)
		if(start.le.0..or.start.ge..5) then
			call error(args,8,2)
			go to 99
		endif
	endif
c
c   find smallest relative minimum of 1/f:
c
	if(ipts.eq.1) go to 31
c
c
	call movct(wk,2048,char(0))
	wk(1)=1.
	call movxy(wk(2),array(nstart(npa)),4*np)
	nf=256
	call fft(wk,wk(257),nf,nf,nf,1)
	do 10 i=1,129
	ii=i+nf
  10	wk(i)=wk(i)*wk(i)+wk(ii)*wk(ii)
	start=0.0
	fm=1.e20
	do 20 i=2,128
	if(wk(i-1).gt.wk(i).and.wk(i).lt.wk(i+1)) then
		if(wk(i).lt.fm) then
			start=float(i-1)/256.
			fm=wk(i)
		endif
	endif
  20	continue
	if(start.eq.0.) then
		write(*,30)
  30		format(' AR spectra has no relative maxima')
		ier=1
		call ckaddi(args(5),ier,iref)
		if(iref.eq.1) go to 99
		go to 99
	endif
c
c   Find Peak Frequency:
c
  31	continue
	call macv(array(nstart(npa)),np,1.,wk,ri0)
	twopi=8.*atan(1.0)
	it=1
	omega=start
  35	fp=0.
	fpp=0.
	do 40 i=1,np
	oi=i
	fp=fp+oi*wk(i)*sin(twopi*oi*omega)
  40	fpp=fpp+oi*oi*wk(i)*cos(twopi*oi*omega)
	if(abs(fpp).lt.1.e-20) then
		write(*,45)
  45		format(' Zero 2nd Derivative Encountered')
		ier=2
		go to 99
	endif
	del=fp/(fpp*twopi)
	omega1=omega-del
	write(*,46) omega1
  46	format(1x,f14.6)
	if(abs(del/omega).lt.1.e-5) then
		if(omega1.lt.1.e-5.or.omega1.gt..49999) then
			write(*,50)
  50			format(' Convergence to 0 or .5')
			ier=3
		endif
		call ckaddr(args(6),omega1,iref)
		if(iref.eq.1) go to 99
		start=omega1
		go to 60
	endif
	it=it+1
	omega=omega1
	if(it.gt.100) then
		write(*,55)
  55		format(' Nonconvergence')
		ier=4
		go to 99
	endif
	go to 35
c
c   Asymptotic Standard Error:
c
  60	continue
	np2=np*np
	n1=1
	n2=n1+np
	n3=n2+np2
	n4=n3+np
	n5=n4+np2
	n6=n5+np2
	call movxy(wk,array(nstart(npa)),4*np)
	call freqcv(start,np,wk(n1),n,np,wk(n2),wk(n3),wk(n4),
     1	wk(n5),wk(n6),var)
	sd=sqrt(var)
	call ckaddr(args(7),sd,iref)
	if(iref.eq.1) go to 99
	call ckaddi(args(5),ier,iref)
c
c
  99	continue
	return
	end
c&FREQCV
      SUBROUTINE freqcv(freq,NP,ALPHA,NOBS,NDIM,COVRI,
     1RI,A,C,B,VAR)
C*****************************************************************
C
C   SUBROUTINE TO FIND THE ASYMPTOTIC VARIANCE OF A frequency
C   ESTIMATOR freq FOR AN AR(NP) PROCESS HAVING COEFFS ALPHA.
C
C   INPUT :
C            freq,NP,ALPHA(1),ALPHA(NP)
c            nobs: number of observations
C            NDIM : ROW DIMENSION OF VARIOUS ARRAYS IN CALLING PROG
C
C   OUTPUT :
C            VAR : ASYMPTOTIC VARIANCE
C            COVRI,RI
C
C********************************************************************
C
      DIMENSION ALPHA(NP),COVRI(NDIM,NDIM),RI(NP),A(NDIM,NDIM),
     1C(NDIM,NDIM),B(NP)
C
C
      CALL SCHUR(ALPHA,NDIM,NP,A)
      TWOPI=8.*ATAN(1.0)
      DO 10 IV=1,NP
      DO 10 IJ=1,NP
      CC=0.0
      IVPIJ=IV+IJ
      IJMIV=IJ-IV
      IF(IVPIJ.LE.NP) CC=CC+ALPHA(IVPIJ)
      IF(IJMIV.EQ.0) CC=CC+1.0
      IF(IJMIV.GT.0) CC=CC+ALPHA(IJMIV)
  10  C(IV,IJ)=CC
      DO 40 J=1,NP
      DO 40 K=1,J
         CC=0.0
         DO 30 L=1,NP
         C1=C(J,L)
         C2=0.0
            DO 20 M=1,NP
  20        C2=C2+A(L,M)*C(K,M)
  30     CC=CC+C1*C2
      COVRI(J,K)=CC
  40  COVRI(K,J)=CC
      DO 60 IV=1,NP
      CC=ALPHA(IV)
      IF(IV.EQ.NP) GO TO 60
         NPMIV=NP-IV
         DO 50 J=1,NPMIV
  50     CC=CC+ALPHA(J)*ALPHA(J+IV)
  60  RI(IV)=CC
C
C   FIND VECTOR B :
C
      DO 90 IV=1,NP
  90  B(IV)=FLOAT(IV)*SIN(IV*freq*twopi)
C
C   FIND VAR :
C
c
c   find -h''/4pi
c
	cc=0.0
	do 95 i=1,np
  95	cc=cc+float(i*i)*ri(i)*cos(twopi*float(i)*freq)
	cc=twopi*cc
      VAR=0.0
      DO 100 J=1,NP
      DO 100 K=1,NP
 100  VAR=VAR+B(J)*COVRI(J,K)*B(K)
      VAR=VAR/(cc*cc)
	var=var/float(nobs)

C
C
      RETURN
      END
