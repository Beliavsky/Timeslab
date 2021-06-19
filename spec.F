c&SPEC
c&WINDOW
        subroutine window(args,nargs,vname)
c****************************************************************
c
c   Subroutine to form the window periodogram spectral estimator.
c
c   window(rho,r0,ntrunc,nfreqs,ioptw)
c
c****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
	dimension al2(8)
	al2(1)=2.
	al2(2)=2./3.
	al2(3)=.795
	al2(4)=.539
	al2(5)=.586
	al2(6)=1.
	al2(7)=1.2
	al2(8)=1.666
c
c
	if(ickl(args(1),np,ncov).eq.1) go to 99
	if(ickr(args(2),2,r0,2,0.).eq.1) go to 99
        call ckint(args(3),ntrunc)
        if(ntrunc.le.0.or.ntrunc.gt.ncov) then
		call error(args,3,2)
		go to 99
	endif
	if(icki(args(4),4,nfreqs,1,ncov).eq.1) go to 99
	if(icki(args(5),4,ioptw,1,1).eq.1) go to 99
	if(ioptw.le.5) go to 20
        if(ioptw.gt.8.or.nargs.eq.5) then
		call error(args,5,2)
		go to 99
	endif
	call ckint(args(6),n)
	if(n.lt.1.or.n.gt.ncov.or.nfreqs.lt.n) then
		call error(args,6,2)
		go to 99
	endif
  20	continue
c
c
	if(ickse(nfreqs).eq.1) go to 99
	if(ickse(ntrunc+3*nfreqs).eq.1) go to 99
        nn=nstart(np)-1
	nnn=ntrunc
	if(ioptw.gt.5) nnn=n-1
        do 60 i=1,nnn
  60    wk(i)=r0*array(nn+i)
        call cvspw(wk,r0,ioptw,ntrunc,nfreqs,n,wk(nnn+1),
     1  wk(nnn+nfreqs+1),wk(nnn+2*nfreqs+1))
	if(imflag.eq.1) go to 99
        n1=(nfreqs/2)+1
c
c
        lab='window spectral estimator'
        call ckadda(vname,n1,lab,nnn+2*nfreqs+1,iref)
	if(nargs.eq.7) then
	if(icki(args(6),6,n,1,1).eq.1) go to 99
	c=al2(ioptw)*float(ntrunc)/float(n)
	c=exp(2.0*sqrt(c))
	call ckaddr(args(7),c,iref)
	endif
  99    continue
        return
        end
c&CVSPW
      SUBROUTINE CVSPW(R,R0,IOPT,NTRUNC,NFREQ,n,WORK1,WORK2,SPEC)
C******************************************************************
c
C   GIVEN R0,R(1),...,R(NTRUNC), THE FIRST NTRUNC+1
C   AUTOCOVARIANCES OF A TIME SERIES X(.), SUBROUTINE CVSPW
C   CALCULATES SPEC(1),...,SPEC(NFREQ), A SMOOTHED ESTIMATE
C   OF THE SPECTRAL DENSITY OF X(.) AT THE NFREQ POINTS
C   0,TWOPI/NFREQ,...,TWOPI*(NFREQ-1)/NFREQ.  THE WEIGHTING
C   FUNCTION IS SPECIFIED BY THE USER BY IOPT (SEE BELOW)
C
C   METHOD : THE CHOSEN WEIGHTING FUNCTION WORK1(1),...
C   ,WORK1(NTRUNC) IS CALCULATED FOR A SPECIFIED TRUNCATION
C   POINT NTRUNC.  THEN A VECTOR OF LENGTH NFREQ IS
C   DETERMINED AS (R0,WORK1(1)*R(1),...,WORK1(NTRUNC)*
C   R(NTRUNC),0,...,0).
C   THEN SUBROUTINE FFT IS USED TO FIND THE COSINE
C   TRANSFORM OF THIS VECTOR. THIS TRANSFORM IS DIVIDED BY
C   TWOPI GIVING THE DESIRED SPECTRUM.
C
C   INPUT :
C   NTRUNC; NFREQ; R0,R(1),...,R(NTRUNC)
C
C   IOPT :
C
C   AS GIVEN BY T.W.ANDERSON : THE STAT ANAL OF TIME
C   SERIES, THE FOLLOWING OPTIONS OF WEIGHTING FUNCTIONS
C   ARE AVAILABLE :
C
C   IOPT=1 : truncated periodogram
C   IOPT=2 : BARTLETT
C   IOPT=3 : Tukey
C   IOPT=4 : Parzen
C   IOPT=5 : Bohman
C
C   OUTPUT :
C   SPEC(1),...,SPEC(NFREQ)
C
C   AUXILLIARY : WORK1,WORK2
C
C   SUBROUTINES CALLED : FFT
C
C*******************************************************************
      DIMENSION R(1),WORK1(1),WORK2(NFREQ),SPEC(NFREQ)
C
C   INITIALIZE :
C
      PI=4.*ATAN(1.0)
      TWOPI=2.*PI
      DO 1 I=1,NFREQ
      SPEC(I)=0.
   1  WORK2(I)=0.
      SPEC(1)=R0
      C1=NTRUNC
      C2=NFREQ
C
C   CALCULATE SPECIFIED WEIGHTS :
C
      IF(IOPT.EQ.1) GO TO 2
      IF(IOPT.EQ.2) GO TO 4
      IF(IOPT.EQ.3) GO TO 6
      IF(IOPT.EQ.4) GO TO 8
      IF(IOPT.EQ.5) GO TO 10
      IF(IOPT.EQ.6) GO TO 21
      IF(IOPT.EQ.7) GO TO 23
      IF(IOPT.EQ.8) GO TO 25
      GO TO 99
C
C
C   Truncated Periodogram:
C
   2  DO 3 J=1,NTRUNC
   3  WORK1(J)=2.
      GO TO 14
C
C   BARTLETT :
C
   4  DO 5 J=1,NTRUNC
   5  WORK1(J)=2.*(1.-(FLOAT(J)/C1))
      GO TO 14
C
C   Tukey:
C
   6  DO 7 J=1,NTRUNC
   7  WORK1(J)=1.08+.92*COS((PI*FLOAT(J))/C1)
      GO TO 14
C
C   PARZEN :
C
   8  M1=NTRUNC/2
      DO 9 J=1,M1
      J1=J+M1
      O1=FLOAT(J)/C1
      O2=O1**2
      O3=O1**3
      WORK1(J)=2.*(1.-6.*O2+6.*O3)
   9  WORK1(J1)=4.*(1.-(FLOAT(J1)/C1))**3
	go to 14
c
c   Bohman:
c
  10	do 11 j=1,ntrunc
	uj=float(j)/c1
	piuj=pi*uj
    	work1(j)=(1.-uj)*cos(piuj)+(sin(piuj)/pi)
  11	work1(j)=2.*work1(j)
	go to 14
c
c   Daniell:
c
  21	do 22 j=1,n-1
	uj=float(j)/c1
	piuj=pi*uj
  22	work1(j)=2.*sin(piuj)/piuj
	go to 14
c
c   Bartlett-Priestley:
c
  23	do 24 j=1,n-1
	uj=float(j)/c1
	piuj=pi*uj
	sj=sin(piuj)
	cj=cos(piuj)
	piuj2=piuj*piuj
  24	work1(j)=6.*((sj/piuj)-cj)/piuj2
	go to 14
c
c   Parzen-Cogburn-Davis:
c
  25	do 26 j=1,n-1
	uj=float(j)/c1
  26	work1(j)=2./(1.+uj**4.)
C
C   CALL FFT :
C
  14	nn=ntrunc
	if(iopt.gt.5) nn=n-1
      DO 15 J=1,nn
  15  SPEC(J+1)=WORK1(J)*R(J)
      CALL FFT(SPEC,WORK2,NFREQ,NFREQ,NFREQ,1)
  99  CONTINUE
      RETURN
      END
c&CROSSP
        subroutine crossp(args,nargs)
c*******************************************************************
c
c   Subroutine to process the command
c
c   CROSSP(rhoxy,rhoyx,rx0,ry0,rhoxy0,ntrunc,nfreqs,ioptw,fxyr,fxyi)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),npxy,nxy).eq.1) go to 99
	if(ickl(args(2),npyx,nyx).eq.1) go to 99
 	if(ickr(args(3),3,rx0,2,0.).eq.1) go to 99
 	if(ickr(args(4),4,ry0,2,0.).eq.1) go to 99
 	if(ickr(args(5),5,rhoxy0,0,0.).eq.1) go to 99
        call ckint(args(6),ntrunc)
        if(ntrunc.le.0.or.ntrunc.gt.nxy.or.ntrunc.gt.nyx) then
                call error(args,6,2)
                go to 99
        endif
	if(icki(args(7),7,nfreqs,2,ntrunc).eq.1) go to 99
	if(ickse(nfreqs).eq.1) go to 99
	if(ickse(3*nfreqs+ntrunc).eq.1) go to 99
        call ckint(args(8),ioptw)
        if(ioptw.lt.1.or.ioptw.gt.5) then
                call error(args,8,2)
                go to 99
        endif
        if(ioptw.eq.4.and.2*(ntrunc/2).ne.ntrunc) then
                call error(args,6,2)
                go to 99
        endif
c
c
        srxy0=sqrt(rx0*ry0)
        pi=4.*atan(1.0)
        c1=float(ntrunc)
        c2=float(nfreqs)
        one=1.0
        two=2.0
        go to(5,10,15,20,25) ioptw
   5    do 6 i=1,ntrunc
   6    wk(i)=one
        go to 35
  10    do 11 i=1,ntrunc
  11    wk(i)=one-(float(i)/c1)
        go to 35
  15    do 16 i=1,ntrunc
  16    wk(i)=.54+.46*cos((pi*float(i))/c1)
        go to 35
  20    m1=ntrunc/2
        do 21 i=1,m1
        i1=i+m1
        o1=float(i)/c1
        o2=o1*o1
        o3=o1*o2
        wk(i)=one*(one-6.*(o2-o3))
  21    wk(i1)=two*(one-(float(i1)/c1))**3
	go to 35
  25	do 26 i=1,ntrunc
	o1=float(i)/c1
	piu=pi*o1
  26	wk(i)=(1.-o1)*cos(piu)+(sin(piu)/pi)
  35    continue
        nsxy=nstart(npxy)-1
        nsyx=nstart(npyx)-1
        call movct(wk(ntrunc+1),8*nfreqs,char(0))
        wk(ntrunc+1)=rhoxy0
        n1=(nfreqs/2)+1
        n2=ntrunc+nfreqs+1
        do 40 i=1,ntrunc
  40    wk(ntrunc+1+i)=wk(i)*(array(nsxy+i)+array(nsyx+i))
        call fft(wk(ntrunc+1),wk(n2),nfreqs,nfreqs,nfreqs,1)
	if(imflag.eq.1) go to 99
        do 50 i=ntrunc+1,ntrunc+nfreqs
  50    wk(i)=wk(i)*srxy0
	n3=ntrunc+2*nfreqs+1
	call movxy(wk(n3),wk(ntrunc+1),4*nfreqs)
c
c
        call movct(wk(ntrunc+1),8*nfreqs,char(0))
        wk(ntrunc+1)=0.0
        do 60 i=1,ntrunc
  60    wk(ntrunc+1+i)=wk(i)*(array(nsxy+i)-array(nsyx+i))
        call fft(wk(ntrunc+1),wk(n2),nfreqs,nfreqs,nfreqs,1)
	if(imflag.eq.1) go to 99
        do 70 i=ntrunc+nfreqs+1,ntrunc+2*nfreqs
  70    wk(i)=-wk(i)*srxy0
        lab='real part of cross-spectra'
        call ckadda(args(9),n1,lab,n3,iref)
        if(iref.eq.1) go to 99
        lab='imaginary part of cross-spectra'
        call ckadda(args(10),n1,lab,ntrunc+nfreqs+1,iref)
c
c
  99    continue
        return
        end
c&ARMASP
        subroutine armasp(args,nargs,vname,iptmx)
c*******************************************************************
c
c   Subroutine to calculate arma spectra
c
c   spec=ARMASP(alpha,beta,p,q,rvar,Q)
c
c   spec=ARSP(alpha,p,rvar,Q)
c
c   spec=MASP(beta,q,rvar,Q)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	nordq=0
	if(ickl(args(1),np1,n1).eq.1) go to 99
	if(ickr(args(nargs-1),nargs-1,rvar,2,0.0).eq.1) go to 99
	if(iptmx.eq.1) call ckint(args(2),nordp)
	if(iptmx.eq.2) call ckint(args(2),nordp)
        if(iptmx.eq.3) call ckint(args(3),nordp)
        if(nordp.le.0.or.nordp.gt.n1) then
		call error(args,3,2)
		go to 99
	endif
	if(iptmx.eq.3) then
		if(ickl(args(2),np2,n2).eq.1) go to 99
		call ckint(args(4),nordq)
		if(nordq.le.0.or.nordq.gt.n2) then
			call error(args,4,2)
			go to 99
		endif
	endif
	mm=max0(nordp,nordq)
	if(icki(args(nargs),nargs,nfreqs,2,mm).eq.1) go to 99
	if(ickse(nfreqs).eq.1) go to 99
	nws=2*nfreqs
	if(iptmx.eq.3) nws=3*nfreqs
	if(ickse(nws).eq.1) go to 99
c
c
        n2=2*nfreqs
	call movct(wk,8*nfreqs,char(0))
        wk(1)=1.0
	call movxy(wk(2),array(nstart(np1)),4*nordp)
c
        call fft(wk,wk(nfreqs+1),nfreqs,nfreqs,nfreqs,1)
	if(imflag.eq.1) go to 99
c
c
        n1=(nfreqs/2)+1
	if(iptmx.eq.2) then
		do 200 i=1,n1
 200		wk(i)=rvar*(wk(i)*wk(i)+wk(nfreqs+i)*wk(nfreqs+i))
		lab='MA Spectra'
		go to 97
	endif
	if(iptmx.eq.1) then
		do 201 i=1,n1
      		c=wk(i)*wk(i)+wk(nfreqs+i)*wk(nfreqs+i)
		if(c.lt.1.e-20) then
			write(*,65)
			imflag=1
			go to 99
		endif
 201		wk(i)=rvar/c
		lab='AR Spectra'
		go to 97
	endif

c
        do 82 i=1,n1
  82    wk(n2+i)=wk(i)*wk(i)+wk(nfreqs+i)*wk(nfreqs+i)
	call movct(wk,8*nfreqs,char(0))
        wk(1)=1.0
	call movxy(wk(2),array(nstart(np2)),4*nordq)
c
        call fft(wk,wk(nfreqs+1),nfreqs,nfreqs,nfreqs,1)
	if(imflag.eq.1) go to 99
c
        do 85 i=1,n1
	c1=rvar*(wk(i)*wk(i)+wk(nfreqs+i)*wk(nfreqs+i))
	c=wk(n2+i)
	if(c.lt.1.e-20) then
		write(*,65)
  65		format(' Spectra has infinite value')
		imflag=1
		go to 99
	endif
  85    wk(i)=c1/c
c
c
        lab='ARMA spectra'
  97    call ckadda(vname,n1,lab,1,iref)
c
c
  99    continue
        return
        end
