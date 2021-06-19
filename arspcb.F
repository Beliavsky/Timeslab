c&ARSPCB
c&ARSPCB
	subroutine arspcb(args,nargs)
c*****************************************************************
c
c   Subroutine to handle the command
c
c     ARSPCB(alpha,np,rvar,n,nfreqs,p,fl,fu)
c
c****************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,ickr,icki,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	call ckint(args(2),np)
	if(np.lt.1.or.np.gt.na.or.np.gt.50) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(icki(args(4),4,n,2,0).eq.1) go to 99
	if(icki(args(5),5,nfreqs,2,0).eq.1) go to 99
	call ckreal(args(6),p)
	if(p.le.0..or.p.ge.1.) then
		call error(args,6,2)
		go to 99
	endif
	if(ickse(nfreqs).eq.1) go to 99
	nws=5*(np+1)+3*(np+1)**2+3*nfreqs
	if(ickse(nws).eq.1) go to 99
c
c
	nn=(nfreqs/2)+1
	ndim=np+1
	ndim2=ndim*ndim
	n1=1
	n2=n1+np
	n3=n2+np
	n4=n3+np
	n5=n4+ndim2
	n6=n5+ndim2
	n7=n6+ndim2
	n8=n7+nn
	n9=n8+ndim
	n10=n9+ndim
	n11=n10+nfreqs
	call movxy(wk,array(nstart(npa)),4*np)
	call calcs(wk(n1),rvar,np,wk(n2),r0,wk(n3),gamm0,nfreqs,
     1	n,p,ndim,wk(n4),wk(n5),wk(n6),wk(n7),wk(n8),wk(n9),wk(n10),
     1	wk(n11))
	lab='Lower Confidence Band'
	call ckadda(args(7),nn,lab,n10,iref)
	if(iref.eq.1) go to 99
	lab='Upper Confidence Band'
	call ckadda(args(8),nn,lab,n11,iref)
c
c
  99	continue
	return
	end
c&CALCS
      SUBROUTINE CALCS(ALPHA,SIG2,NORD,R,R0,GAMM,GAMM0,NFREQS,
     1N,P,NDIM,SIGB,BB,SIGG,S,x,d,fl,fu)
C**************************************************************
C
C   SUBROUTINE TO CALCULATE THE QUANTITY S AT THE NFREQS EQUALLY
C   SPACED FREQUENCIES BETWEEN 0 AND PI INCLUSIVE AS IN
C   "SIMULTANEOUS CONFIDENCE BANDS FOR AUTOREGRESSIVE SPECTRA".
C   ON OUTPUT SIGB,BB,SIGG ARE SIGMA(BETA),B(BETA),SIGMA(GAMMA)
C   RESPECTIVELY. ALSO NDIM.GT.NORD. ALSO S IS CALCULATED FOR 100P
C   PERCENT BANDS.
C
C***************************************************************
C
c   s is (nfreqs/2)+1, d and x are nord+1
      DIMENSION ALPHA(NORD),R(NORD),SIGB(NDIM,NDIM),BB(NDIM,NDIM),
     1SIGG(NDIM,NDIM),S(1),X(1),GAMM(NORD),D(1),fl(1),fu(1)
	double precision dchi2,dp
C
C
	obrv=1./sig2
      CALL MACV(ALPHA,nord,obrv,GAMM,GAMM0)
C
C   FIND SIGB :
C
      CALL SCHUR(ALPHA,NDIM,NORD,SIGB)
      NP1=NORD+1
      DO 30 I=1,NORD
      SIGB(NP1,I)=0.0
  30  SIGB(I,NP1)=0.0
      SIGB(NP1,NP1)=2./(SIG2*SIG2)
C
C   FIND BB :
C
      BB(1,NP1)=SIG2*GAMM0
      DO 40 L=2,NP1
  40  BB(L,NP1)=SIG2*GAMM(L-1)
      DO 50 L=1,NP1
      DO 50 M=1,NORD
  50  BB(L,M)=0.0
      DO 70 J=1,NORD
      DO 60 IVP1=1,NP1
      IV=IVP1-1
      IF(J.GT.NORD-IV) GO TO 55
      BB(IVP1,J)=BB(IVP1,J)+ALPHA(J+IV)/SIG2
  55  IF(J.LT.IV) GO TO 60
      C1=1./SIG2
      IF(J.NE.IV) C1=C1*ALPHA(J-IV)
      BB(IVP1,J)=BB(IVP1,J)+C1
  60  CONTINUE
  70  CONTINUE
C
C   FIND S :
C
	dp=dble(p)
	call chiqnt(dp,np1,dchi2,ier)
	chi2=dchi2
      FAC=CHI2/FLOAT(N)
      CALL MCHOL(SIGB,NDIM,NP1,SIGG,D,IER)
      DO 311 I=1,NP1
 311  D(I)=SQRT(D(I))
      DO 312 I=1,NP1
      DO 312 J=I,NP1
 312  SIGG(I,J)=D(I)*SIGG(I,J)
      DO 320 I=1,NP1
      DO 320 J=1,NP1
      C=0.0
      DO 315 L=I,NP1
 315  C=C+SIGG(I,L)*BB(J,L)
 320  SIGB(I,J)=C
      X(1)=1.
      ON=FLOAT(NFREQS)
	n1=(nfreqs/2)+1
      DO 120 IW=1,N1
	w=float(i-1)/on
      C=0.0
      DO 100 J=2,NP1
      OJ=FLOAT(J-1)*twopi
 100  X(J)=2.*COS(OJ*W)
      C1=0.0
      DO 313 I=1,NP1
      C=0.0
      DO 314 J=1,NP1
 314  C=C+SIGB(I,J)*X(J)
 313  C1=C1+C*C
 120  S(IW)=SQRT(FAC*C1)
	do 400 i=1,nfreqs
	fl(i)=0.
 400	fu(i)=0.
	fl(1)=1.
	do 401 i=1,nord
 401	fl(i+1)=alpha(i)
	call fft(fl,fu,nfreqs,nfreqs,nfreqs,1)
	do 402 i=1,n1
	fl(i)=(fl(i)*fl(i)+fu(i)*fu(i))/sig2
	fu(i)=10000.
	if(fl(i)-s(i).gt.0) fu(i)=1./(fl(i)-s(i))
 402	fl(i)=1./(fl(i)+s(i))
C
C
C
      RETURN
      END
