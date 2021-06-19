c&COEFFSD
c&ARMASD
	subroutine armasd(args,nargs)
c****************************************************************
c
c   Subroutine to handle the command
c   COEFFSD(coeff,nord,n,sdc) or
c   COEFFSD(alpha,beta,nordp,nordq,n,sdalph,sdbeta,ier)
c
c****************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 icki,ickl,ickse
c
c   nargs=4:
c
	if(nargs.eq.4) then
		if(ickl(args(1),npc,nc).eq.1) go to 99
		call ckint(args(2),nord)
		if(nord.le.0.or.nord.gt.nc) then
			call error(args,2,2)
			go to 99
		endif
		call schur(array(nstart(npc)),nord,nord,wk)
		call ckint(args(3),n)
		if(n.le.0) then
			call error(args,3,2)
			go to 99
		endif
		on=n
		do 5 i=1,nord
  5		wk(nord*nord+i)=sqrt(wk((i-1)*nord+i)/on)
		lab='Standard Deviation of Coefficients'
		call ckadda(args(4),nord,lab,nord*nord+1,iref)
		go to 99
		endif
c
c   nargs=8:
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(ickl(args(2),npb,nb).eq.1) go to 99
	call ckint(args(3),nordp)
	if(nordp.le.0.or.nordp.gt.na) then
		call error(args,3,2)
		go to 99
	endif
	call ckint(args(4),nordq)
	if(nordq.le.0.or.nordq.gt.nb) then
		call error(args,4,2)
		go to 99
	endif
	if(icki(args(5),5,n,2,0).eq.1) go to 99
c
c
	mdim=2*max0(nordp,nordq)+1
	n1=mdim*mdim+1
	n2=n1+mdim
	n3=n2+mdim*mdim
	n4=n3+nordp+nordq
	if(ickse(n4).eq.1) go to 99
c
c
	call mxcsd(array(nstart(npa)),array(nstart(npb)),nordp,nordq,
     1	n,mdim,wk,wk(n1),wk(n2),wk(n3),ier)
	call ckaddi(args(8),ier,iref)
	if(ier.eq.1.or.iref.eq.1) go to 99
	lab='STD DEV of AR Coefficients'
	call ckadda(args(6),nordp,lab,n3,iref)
	if(iref.eq.1) go to 99
	lab='STD DEV of MA Coefficients'
	call ckadda(args(7),nordq,lab,n3+nordp,iref)
c
c
  99	continue
	return
	end
c&MXCSD
      SUBROUTINE MXCSD(ALPH,BETA,NP,NQ,N,MDIM,WK1,R,
     1AINF,SD,IER)
C**********************************************************
C
C   SUBROUTINE TO find the standard deviation of the estimated
c   coefficients (based on N observations) of an ARMA(NP,NQ,ALPH,BETA)
c   process.
C
C   INPUT :
C           NP,NQ,ALPH,BETA
C           MDIM : DIMENSION OF VARIOUS ARRAYS IN
C           CALLING PROGRAM (SEE DIMENSIONS BELOW)
C           MDIM.GE.(2*MAX(NP,NQ))+1
C
C   OUTPUT :
C           AINF
C           IER : 0 IS NORMAL RETURN, 1 MEANS THEREIS
C           A SINGULAR MATRIX IN THE PROCEDURE
c           sd(1),...,sd(np+nq)
C
C   SUBROUTINES CALLED : MXCVIN,SCHUR,swpk12,decomp,solv
C
C**********************************************************
C
      DIMENSION ALPH(NP),BETA(NQ),WK1(mDIM,mDIM),
     1AINF(MDIM,MDIM),R(1),sd(1)
C
C   IAB :
C
	if(np*nq.eq.0) go to 99
	ier=0
	on=n
      M=MAX0(NP,NQ)
      MM=2*M+1
      NPPNQ=NP+NQ
      CALL MXCVIN(ALPH,BETA,NP,NQ,MM,MDIM,AINF,R,IER)
	if(ier.eq.1) go to 99
      DO 1 I=1,NPPNQ
      DO 1 J=1,NPPNQ
   1  AINF(I,J)=0.
      DO 2 I=1,NP
      DO 2 J=1,NQ
      AINF(I,NP+J)=-R(I-J+M+1)
   2  CONTINUE
C
C   IAA :
C
      CALL SCHUR(ALPH,mDIM,NP,WK1)
	call swpk12(wk1,mdim,np,1,np,ier)
	if(ier.eq.1) go to 99
      DO 3 I=1,NP
      DO 3 J=I,NP
   3  AINF(I,J)=wk1(i,j)
C
C   IBB :
C
      CALL SCHUR(beta,mDIM,NP,WK1)
	call swpk12(wk1,mdim,nq,1,nq,ier)
	if(ier.eq.1) go to 99
      DO 4 I=1,NQ
      II=NP+I
      DO 4 J=I,NQ
   4  AINF(II,NP+J)=wk1(i,j)
C
      DO 5 I=1,NPPNQ
      DO 5 J=I,NPPNQ
   5  AINF(J,I)=AINF(I,J)
	call swpk12(ainf,mdim,nppnq,1,nppnq,ier)
	on=n
	do 10 i=1,nppnq
  10	sd(i)=sqrt(ainf(i,i)/on)
  99	continue
C
      RETURN
      END
c&MXCVIN
      SUBROUTINE MXCVIN(ALPH,BETA,NP,NQ,MM,NDIM,A,R,IER)
C**********************************************************
C
C   SUBROUTINE TO CALCULATE THE CROSS COVARIANCES FOR
C   LAGS (-MAX(NP,NQ),...,MAX(NP,NQ)), (STORED IN
C   R(1),...,R(MM),MM=(2*MAX(NP,NQ))+1) OF THE TWO
C   DIMENSIONAL PROCESS (X1(T),X2(T)) WHERE X1(.) AND
C   X2(.) ARE DEFINED BY :
C
C      SUM(J=0,NP) ALPH(J)*X1(T-J)=E(T)
C      SUM(J=0,NQ)BETA(J)*X2(T-J)=E(T)
C
C   INPUT :
C           NP,NQ,MM=2*MAX(NP,NQ)+1
C           ALPH,BETA
C           NDIM : DIMENSION OF A IN CALLING PROGRAM
C           (NDIM.GE.MM)
C
C   OUTPUT :
C           R(1),...,R(MM) (R(J) IS R12(J-MAX(NP,NQ)-1)
C           IER : 0 IS NORMAL RETURN, 1 MEANS A IS SINGULAR
C
C   SUBROUTINES CALLED : DECOMP,SOLV
C
C**********************************************************
C
      DIMENSION ALPH(NP),BETA(NQ),A(NDIM,NDIM),R(MM),
     1IP(100)
      DIMENSION LL(2)
      DATA NOUT/0/
C
C
C
      M=MAX0(NP,NQ)
      MP1=M+1
      MP2=M+2
      DO 1 I=1,MM
      DO 1 J=1,MM
   1  A(I,J)=0.
C
C
      IF(NP.LT.NQ) GO TO 50
C
C   FOR NP.GE.NQ :
C
      A(1,1)=1.
      DO 5 I=1,NP
   5  A(1,I+1)=ALPH(I)
      DO 6 J=2,MP1
      JM1=J-1
      A(J,1)=A(JM1,MM)
         DO 7 K=2,MM
   7     A(J,K)=A(JM1,K-1)
   6  CONTINUE
      M1=M-NQ+1
      DO 8 J=1,NQ
   8  A(MP2,M1+J)=BETA(NQ-J+1)
      A(MP2,M1+NQ+1)=1.
      MP3=M+3
      DO 9 J=MP3,MM
      JM1=J-1
      A(J,1)=A(JM1,MM)
         DO 10 K=2,MM
  10     A(J,K)=A(JM1,K-1)
   9  CONTINUE
      GO TO 100
C
C   FOR NP.LT.NQ :
C
  50  A(1,1)=1.
      DO 51 I=1,NP
  51  A(1,I+1)=ALPH(I)
      IF(M.EQ.1) GO TO 70
      DO 52 J=2,M
      JM1=J-1
      A(J,1)=A(JM1,MM)
         DO 53 K=2,MM
  53     A(J,K)=A(JM1,K-1)
  52  CONTINUE
  70  DO 71 J=1,NQ
  71  A(MP1,J)=BETA(NQ-J+1)
      A(MP1,NQ+1)=1.
      DO 72 J=MP2,MM
      JM1=J-1
      A(J,1)=A(JM1,MM)
         DO 73 K=2,MM
  73     A(J,K)=A(JM1,K-1)
  72  CONTINUE
C
C   SOLVE SYSTEM OF EQUATIONS :
C
 100  DO 101 I=1,MM
 101  R(I)=0.
      R(MP1)=1.
C
C
      CALL DECOMP(MM,NDIM,A,IP)
      IF(IP(MM).EQ.0) GO TO 125
      IER=0
      CALL SOLV(MM,NDIM,A,R,IP)
C
C   R CONTAINS R12(-M),...,R12(M)
C
      RETURN
C
 125  IER=1
      WRITE(NOUT,126)
 126  FORMAT(10X,25HSINGULAR MATRIX IN MXCVIN)
      RETURN
      END
c&SCHUR
      SUBROUTINE SCHUR(ALPH,NDIM,NP,A)
C*******************************************************
C
C   SUBROUTINE TO FORM THE SCHUR MATRIX OF THE COEFFICIENTS
C   ALPH(1),...,ALPH(NP) (SEE PAGANO : WHEN IS AN
C   AUTOREGRESSIVE SCHEME STATIONARY )
C
C   INPUT :
C           NDIM : DIMENSION OF MATRIX A IN MAIN PROGRAM
C           NP, ALPH(1),...,ALPH(NP)
C
C   OUTPUT :
C           A
C
C   SUBROUTINES CALLED : NONE
C
C*******************************************************
C
      DIMENSION A(NDIM,NDIM),ALPH(NP)
      DO 1 J=1,NP
      DO 1 K=J,NP
      C=0.
         DO 6 L=1,J
         IF(L.EQ.J) GO TO 2
         D1=ALPH(J-L)
         GO TO 3
   2     D1=1.
   3     IF(L.EQ.K) GO TO 4
         D2=ALPH(K-L)
         GO TO 6
   4     D2=1.
   6     C=C+D1*D2-ALPH(NP+L-J)*ALPH(NP+L-K)
      A(J,K)=C
   1  A(K,J)=A(J,K)
      RETURN
      END

