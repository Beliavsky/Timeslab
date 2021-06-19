c&ARMASL
c&ARMASL
	subroutine armasl(args,nargs)
c******************************************************************
c
c   Subroutine to handle the command
c
c   ARMASEL(x,n,mr,ma,k1,k2,kopt,p,nordp,nordq,alpha,beta,rvar,ier)
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki,ickr,ickse
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	call ckint(args(2),nn)
	n=iabs(nn)
	if(icki(args(3),3,mr,2,0).eq.1) go to 99
	if(n.le.mr.or.n.gt.nx) then
		call error(args,2,2)
		go to 99
		endif
	if(icki(args(4),4,ma,0,0).eq.1) go to 99
	if(ma.le.0.or.ma.ge.mr) then
		call error(args,4,2)
		go to 99
		endif
	m1=mr-ma
	if(m1.gt.40) then
		if(ickse(maxnum+1).eq.1) go to 99
	endif
	call ckint(args(5),k1)
	if(k1.lt.0.or.k1.gt.m1) then
		call error(args,5,2)
		go to 99
		endif
	call ckint(args(6),k2)
        if(k2.lt.0.or.k2.gt.m1) then
		call error(args,6,2)
		go to 99
        	endif
	if(icki(args(7),7,kopt,0,0).eq.1) go to 99
		iptin=0
		if(kopt.lt.0) then
		iptin=1
			if(-kopt.ne.k1+k2) then
			call error(args,7,2)
			go to 99
			endif
		if(k1.ne.0) then
			if(ickl(args(11),npk1,nk1).eq.1) go to 99
			do 10 i=1,k1
  10			iwk(i)=array(nstart(npk1)+i-1)
			endif
		if(k2.ne.0) then
			if(ickl(args(12),npk2,nk2).eq.1) go to 99
			do 20 i=1,k2
  20			iwk(k1+i)=-array(nstart(npk2)+i-1)
			endif
		k1=-k1
		k2=-k2
		kopt=-kopt
		go to 30
		endif
	if(kopt.gt.k1+k2) then
		call error(args,7,2)
		go to 99
		endif
  30	call ckreal(args(8),p)
	if(p.le.0..or.p.ge.1.) then
		call error(args,8,2)
		go to 99
		endif
c
c
c
	call dtcvc(array(nstart(npx)),n,mr,wk,ry0)
	call cvar(wk,ry0,1,n,ma,ma,wk(2*mr+1),wk(mr+1),rva)
	ndim=2*(mr-ma)+1
	n1=2*mr+1+ndim*ndim
	call cvmxsl(wk,ry0,wk(mr+1),mr,ma,ndim,k1,k2,kopt,nn,p,
     1	wk(2*mr+1),kp,iwk,wk(n1),iwk(ndim+1),rvar,ier,wk(n1+400),
     1	wk(n1+500),wk(n1+600),wk(n1+700),iwk(ndim+200),iptrec)
	call ckaddi(args(14),ier,iref)
	if(ier.eq.1.or.iref.eq.1) then
		imflag=1
		go to 99
		endif
	call slmx(kp,iwk(ndim+1),wk(n1),wk(n1+100),wk(n1+200),nordp,
     1	nordq)
	call ckaddi(args(9),nordp,iref)
	if(iref.eq.1) go to 99
	call ckaddi(args(10),nordq,iref)
	if(iref.eq.1) go to 99
	if(nordp.gt.0) then
		lab='AR part of Select ARMA'
		call ckadda(args(11),nordp,lab,n1+100,iref)
		if(iref.eq.1) go to 99
		endif
	if(nordq.gt.0) then
		lab='MA part of Select ARMA'
		call ckadda(args(12),nordq,lab,n1+200,iref)
		if(iref.eq.1) go to 99
		endif
	call ckaddr(args(13),rvar,iref)
c
c
  99	continue
     	return
	end
c&DTCVC
	subroutine dtcvc(x,n,m,rx,rx0)
c*****************************************************************
c
c   Subroutine to find the sample variance rx0 and first m sample
c   autocovariances rx(1),...,rx(m) of x(1),...,x(n)
c
c******************************************************************
c
	dimension x(n),rx(m)
	double precision inprod
c
	on=n
	rx0=inprod(x,x,n)/on
	do 10 i=1,m
  10	rx(i)=inprod(x,x(i+1),n-i)/on
c
	return
	end
c&SLMX
      SUBROUTINE SLMX(KP,IND,COF,ALPHA,BETA,NORDAR,NORDMA)
C****************************************************************
C
C   SUBROUTINE TO TRANSFORM MIXED SELECT SCHEME PARAMETERS
C   KP,IND,COF TO FULL MIXED SCHEME PARAMETERS NORDAR,NORDMA,
C   ALPHA,BETA.
C
C   INPUT :
C            KP : NUMBER OF LAGS IN SELECT SCHEME
C            IND(1),...,IND(KP) : LAGS (AR LAGS +, MA LAGS -)
C            COF(1),...,COF(KP) : COEFFS CORRESPONDING TO LAGS
C
C   OUTPUT :
C            NORDAR,NORDMA : AR AND MA PART ORDERS
C            ALPHA(1),...,ALPHA(NORDAR) : AR COEFFS
C            BETA(1),...,BETA(NORDMA) : MA COEFFS
C
C***************************************************************
C
      DIMENSION IND(KP),COF(KP),ALPHA(1),BETA(1)
C
C   FIND NORDAR,NORDMA :
C
      NORDAR=0
      NORDMA=0
      IF(KP.EQ.0) GO TO 99
      DO 20 I=1,KP
      II=IND(I)
      IF(II.GT.0) GO TO 10
      IF(-II.GT.NORDMA) NORDMA=-II
      GO TO 20
  10  CONTINUE
      IF(II.GT.NORDAR) NORDAR=II
  20  CONTINUE
C
C   ZERO OUT ALPHA AND BETA :
C
      IF(NORDAR.EQ.0) GO TO 40
      DO 30 I=1,NORDAR
      ALPHA(I)=0.0
  30  CONTINUE
  40  CONTINUE
      IF(NORDMA.EQ.0) GO TO 60
      DO 50 I=1,NORDMA
      BETA(I)=0.0
  50  CONTINUE
  60  CONTINUE
C
C   PUT COEFFS INTO ALPHA,BETA :
C
      DO 80 I=1,KP
      II=IND(I)
      IF(II.GT.0) GO TO 70
      IJ=-II
      BETA(IJ)=COF(I)
      GO TO 80
  70  CONTINUE
      ALPHA(II)=COF(I)
  80  CONTINUE
C
C   FINISH :
C
  99  CONTINUE
      RETURN
      END
c&SELREG
      SUBROUTINE SELREG(A,NDIM,K,NAME,KTOP,L,p,KP,COF,IND,IER,iptrec)
C*******************************************************
C
C   SUBROUTINE TO PERFORM SUBSET REGRESSION
C
C   INPUT :
C           A : COVARIANCE MATRIX TO BE ANALYZED
C           NDIM : DIMENSION OF A IN CALLING PROGRAM
C           K : NUMBER OF INDEPENDENT VARIABLES IN FULL
C           MODEL
C           NAME : NAME(I)=VARIABLE NUMBER OF VARIABLE ON
C           ITH DIAGONAL OF A.
C           KTOP : INDICATOR (0 MEANS SELREG AUTOMATICALLY
C           CHOOSES BEST MODEL, KTOP = M FORCES BEST M
C           VARIABLES INTO MODEL)
C           KTOP = -1 MEANS TERMINATE AFTER FIRST DELETION
C           IF RES VAR .LT. 0.100  OR TERMINATE IF RES VAR
C           .LT.0.001
C   IF KTOP.GE.99, CHISQ VALUES FOR P=.90 ARE USED AND
C   KTOP-100 IS USED TO DETERMINE OTHER OPTIONS.
C   IF KTOP.GE.199, THEN CHISQ VALUES FOR P=.80 ARE USED,
C   AS IS THE VALUE OF KTOP-200
C           L : SAMPLE SIZE (IF L IS NEGATIVE, NO PRINTING
C           IS DONE)
C
C   OUTPUT :
C           KP : NUMBER OF PREDICTOR VARIABLES IN CHOSEN
C           MODEL
C           IND : VECTOR TELLING WHICH VARIABLES (ACCORDING
C           TO NAME) HAVE BEEN CHOSEN
C           COF : VECTOR OF COEFFICIENTS FOR CHOSEN MODEL
C           IER : ERROR INDICATOR (1 MEANS ALL DIAGONAL
C           ELEMENTS OF A.LT.1.E-5 AT A CYCLE, 0 IS
C           NORMAL RETURN)
C
C   SUBROUTINES CALLED : NONE
C
C*******************************************************
C
      DIMENSION A(NDIM,NDIM),NAME(K),COF(K),IND(K),
     1CB(100),IND1(100)
	double precision du,df
C
      IOPTP=1
      LLOPT=L
      IF(L.GT.0) GO TO 50
      L=-L
      IOPTP=0
  50  CONTINUE
      KTOP1=KTOP
      IF(KTOP.GE.99) GO TO 2
	du=.95d0
      GO TO 4
   2  IF(KTOP.GE.199) GO TO 3
	du=.9d0
      KTOP=KTOP-100
      GO TO 4
   3    du=.8d0
      KTOP=KTOP-200
   4  CONTINUE
	du=p
C
      DO 10 I=1,K
      CB(I)=0.
  10  IND1(I)=0
      IER=0
      IF(IOPTP.EQ.0) GO TO 51
      WRITE(*,41) K,KTOP
	if(iptrec.eq.1) write(3,41,err=701) k,ktop
 701	continue
  41  FORMAT(//,10X,'REGRESSION ESTIMATION STAGEWISE SUMMARY'/10X,
     1'NUMBER OF VARIABLES IN FULL MODEL = ',I3/10X,
     1'MAXIMUM NUMBER OF COEFFICIENTS = ',I3/)
	if(iptrec.eq.1) write(3,1,err=702)
 702	continue
      WRITE(*,1)
   1  FORMAT(1X,75(1H-))
  51  CONTINUE
      N=K+1
C
C
C
      F2=6.635
      K1=0
      TOL=1.E-5
      KC=0
 100  KP=0
C
C
C   STAGE 1 : ANALYZE MATRIX A
C
C
      RTOT=0.
      VMIN=100000.
      VMAX=0.
      NMIN=0
      NMAX=0
C
C   KTOL COUNTS NUMBER OF DIAGONALS.GT.TOL
C
      KTOL=0
      DO 110 IKT=1,K
      V=A(IKT,N)*A(N,IKT)
C
C   CHECK IF VARIABLE IKT IS IN MODEL
C
      IF(V.LT.0.) GO TO 105
      CB(IKT)=0.
      IND1(IKT)=0
      GO TO 106
 105  CB(IKT)=A(IKT,N)
      IND1(IKT)=1
      KP=KP+1
C
C   CHECK DIAGONAL A(IKT,IKT) :
C
 106  IF(A(IKT,IKT).LT.TOL) GO TO 110
      KTOL=KTOL+1
      V=V/A(IKT,IKT)
      IF(V.LT.0.) GO TO 60
      RTOT=RTOT+V
      GO TO 80
C
C   FIND VMIN,NMIN :
C
  60  IF((ABS(V)-ABS(VMIN)).GT.0.) GO TO 110
      VMIN=V
      NMIN=IKT
      GO TO 110
C
C   FIND VMAX,NMAX :
C
  80  IF((V-VMAX).LE.0.) GO TO 110
      VMAX=V
      NMAX=IKT
 110  CONTINUE
C
C   ARE ALL DIAGONALS.LT.TOL ?
C
      IF(KTOL.EQ.0) GO TO 999
C
C   CALCULATE AND PRINT CRITERION :
C
      PHI=L-KP-1
      KT1=K-KP
      IF(KT1.EQ.0) GO TO 120
	call chiqnt(du,kt1,df,ierr)
	f1=df
      GO TO 125
 120  F1=1.E-10
 125  T1=(ABS(VMIN)*PHI)/(F2*A(N,N))
      T2=(RTOT*PHI)/(F1*A(N,N))
C
C
      IF(NMIN.EQ.0) GO TO 130
      KD1=NAME(NMIN)
      GO TO 140
 130  KD1=0
 140  IF(NMAX.EQ.0) GO TO 150
      KD2=NAME(NMAX)
      GO TO 160
 150  KD2=0
 160  IF(K1.EQ.0) GO TO 170
      KD3=NAME(K1)
      GO TO 180
 170  KD3=0
C
C   AKAIKE :
C
 180  PHI0=L+2*KP
      PHI1=PHI0+2.
      PHI2=PHI0-2.
      PHI3=2.*FLOAT(KP)/FLOAT(L)
      T3=ALOG(A(N,N))+PHI3
      T4=(VMAX*PHI1)/(2.*A(N,N))
      T5=(VMIN*PHI2)/(2.*A(N,N))
C
C
C
C   PRINT COEFFICIENTS AND VARIABLES
C
      IF(IOPTP.EQ.0) GO TO 52
      WRITE(*,190)
	if(iptrec.eq.1) write(3,190,err=703)
 703	continue
 190  FORMAT(10X,12HVAR IN MODEL,10X,5HCOEFF)
  52  CONTINUE
      IF(KP.EQ.0) GO TO 184
      I1=1
      DO 191 I=1,K
      IF(IND1(I).EQ.0) GO TO 191
      I2=NAME(I)
      IND(I1)=I2
      COF(I1)=CB(I)
      I1=I1+1
      IF(IOPTP.EQ.0) GO TO 191
      WRITE(*,192) I2,CB(I)
	if(iptrec.eq.1) write(3,192,err=704) i2,cb(i)
 704	continue
 192  FORMAT(10X,I12,5X,F10.5)
 191  CONTINUE
 184  IF(IOPTP.EQ.0) GO TO 55
	if(iptrec.eq.1) write(3,185,err=705)
 705	continue
      WRITE(*,185)
 185	format(1x,3x,5hC ADD,3x,5hC DEL,3x,5HV ADD,3x,5hV DEL,
     1	3x,5hV LST,3x,5hNPRED,2x,6hNO CYC,5x,4hRVAR)
	write(nout,186) t2,t1,kd2,kd1,kd3,kp,kc,a(n,n)
	if(iptrec.eq.1) write(3,186,err=710) t2,t1,kd2,kd1,kd3,kp,
     1	kc,a(n,n)
 710	continue
 186	format(1x,f8.3,f8.3,5i8,f9.4)
      WRITE(*,1)
	if(iptrec.eq.1) write(3,1,err=706)
 706	continue
  55  CONTINUE
C
C
C   STAGE 2 : CHOSE PIVOT OR TERMINATE
C
C
 200  IF((KTOP.GT.0).AND.(KP.LT.KTOP).AND.(NMAX.EQ.0)) GO TO 919
      IF((KTOP.GT.0).AND.(KP.LT.KTOP).AND.(NMAX.NE.0)) GO TO 300
      IF((KTOP.GT.0).AND.(KP.EQ.KTOP)) GO TO 280
C
C   AUTOMATIC :
C
C   CHECK FOR DELETE :
C
C
C   CHECK KTOP=-1 CONDITIONS :
C
      IF(KTOP.NE.-1) GO TO 224
      IF((KP.LT.KC).AND.(A(N,N).LT.0.1)) GO TO 280
      IF(A(N,N).LT.0.001) GO TO 225
C
C   CHECK DELETE CONDITIONS FOR KTOP.NE.-1 :
C
 224  IF((ABS(VMIN).GE.F2*A(N,N)/PHI).OR.(K1.EQ.NMIN)) GO TO 230
C
C   YES
C
 225  K1=NMIN
      GO TO 240
C
C   NO
C
C
C   CHECK ADD :
C
 230  IF(RTOT.LE.F1*(A(N,N)/PHI)) GO TO 280
C
C   YES
C
 300  K1=NMAX
 240  CONTINUE
C
C
C   STAGE 3 : PIVOT
C
C
      PIVOT=1./A(K1,K1)
      DO 350 JKT=1,N
 350  A(K1,JKT)=A(K1,JKT)*PIVOT
      DO 370 IKT=1,N
      IF(IKT.EQ.K1) GO TO 370
      TEMP=A(IKT,K1)
         DO 360 JKT=1,N
 360     A(IKT,JKT)=A(IKT,JKT)-TEMP*A(K1,JKT)
      A(IKT,K1)=-TEMP*PIVOT
 370  CONTINUE
      A(K1,K1)=PIVOT
C
C   RETURN  TO STAGE 1
C
      KC=KC+1
      IF(KC.LT.100) GO TO 100
C
C   TERMINATIONS :
C
      WRITE(*,913)
	if(iptrec.eq.1) write(3,913,err=707)
 707	continue
 913  FORMAT(/,10X,'TERMINATION OF SELREG SINCE KC.GE.100')
      GO TO 280
C
C
  919 WRITE(*,920)
	if(iptrec.eq.1) write(3,920,err=708)
 708	continue
 920  FORMAT(/,10X,'TERMINATION SINCE ALL REMAINING VARIABLES HAVE 0 PAR
     1TIAL VARIANCE')
      GO TO 280
C
C
  999 WRITE(*,1000)
	if(iptrec.eq.1) write(3,1000,err=709)
 709	continue
 1000 FORMAT(/,10X,'TERMINATION SINCE ALL DIAGONALS.LT.TOL')
      IER=1
 280  KTOP=KTOP1
      L=LLOPT
      RETURN
C
C
      END
c&CVMXSL
      SUBROUTINE CVMXSL(RY,RY0,ALPH,MR,MA,NDIM,K1,K2,KTOP,
     1L,p,A,KP,NAME,COF,IND,rvar,IER,rey,ree,wk,rye,nwk,iptrec)
C*******************************************************
C
C   SUBROUTINE TO PERFORM SUBSET MIXED SCHEME ESTIMATION
C        (G(L)Y(T)=H(L)E(T))
C
C   INPUT :
C           MR,RY0,RY(1),...,RY(MR) : 1ST MR+1 AUTOCOVARIANCES
C           OF Y
C      NOTE : IF MR IS NEGATIVE, THE COVARIANCES OBTAINED
C      FROM ARMXCV ARE NOT PRINTED
C           NDIM : DIMENSION OF MATRIX A IN CALLING
C           PROGRAM
C           MA,ALPH(1),...,ALPH(M) : COEFFICIENTS OF MATH
C           ORDER AUTOREGRESSIVE REPRESENTATION OF Y
C
C   K1,K2 :
C
C           NUMBER OF LAGS OF Y AND E RESPECTIVELY TO BE
C           INCLUDED IN THE FULL MODEL.
C           IF(K1.LT.MR-MA), THE MOST IMPORTANT LAGS OF Y
C           ARE CHOSEN FOR INCLUSION BY THE MAGNITUDE OF
C           THEIR AUTOCOVARIANCE. THE SAME PROCESS IS
C           CARRIED OUT FOR LAGS OF E IF(K2.LT.MR-MA).
C           IF(K1.EQ.0) THE RESULT IS A PURE MOVING
C           AVERAGE STEPWISE REGRESSION. IF(K2.EQ.0), THE
C           RESULT IS A PURE AUTOREGRESSIVE STEPWISE
C           REGRESSION.
C           IF(K1.EQ.-J1), OR (K2.EQ.-J2) OR BOTH, THE LAGS
C           SPECIFIED BY USING THE VECTOR NAME AS INPUT
C           (LAGS OF Y BEING POSITIVE AND IN THE FIRST
C           J1 POSITIONS, THE LAGS OF E BEING NEGATIVE
C           AND IN THE NEXT J2 POSITIONS) ARE FORCED INTO THE MODEL.
C
C           KTOP : INDICATOR (IF KTOP = 0, THE BEST MODEL
C           IS CHOSEN AUTOMATICALLY, IF KTOP = J, THE J
C           MOST IMPORTANT LAGS ARE FORCED INTO THE MODEL)
C           L : SAMPLE SIZE
C
C   OUTPUT :
C           KP : NUMBER OF PREDICTORS IN CHOSEN MODEL
C           NAME(1),...,NAME(K1+K2) : VECTOR OF LAGS IN FULL
C           MODEL
C
C   NOTE : IN VECTOR NAME, YLAGS ARE POSITIVE, E LAGS ARE
C   NEGATIVE.
C
C           COF : COEFFICIENTS OF CHOSEN LAGS
C           IND : LAGS CHOSEN (ACCORDING TO NAME)
C           IER (1 MEANS ALL DIAGONALS ARE LT 1.E-5 FOR
C           SOME CYCLE, 0 IS NORMAL RETURN)
C
C   SUBROUTINES CALLED : SELREG,ARMXCV,MAX
C
C*******************************************************
C
      DIMENSION RY(1),ALPH(MA),A(NDIM,NDIM),NAME(1),COF(1),
     1IND(1),REY(1),REE(1),WK(1),NWK(1),rye(1)
C
C   FIND AUTOCOVARIANCES :
C
      MOPT=MR
      IOPT=0
      IF(MR.GT.0) GO TO 110
      IOPT=1
      MR=-MR
 110  CONTINUE
      IF(MR-MA.GT.1) GO TO 103
      WRITE(*,104)
 104  FORMAT(10X,'MR-MA.LT.1 IN CVMXSL')
      GO TO 99
 103  CONTINUE
C
      CALL ARMXCV(RY,RY0,ALPH,MR,MA,MR-MA,REY,REY0,RYE,REE,REE0,
     1IOPT)
C
C
C   DETERMINE VECTOR NAME :
C
C
      KK=IABS(K1)+IABS(K2)
      IF(KK.EQ.0) RETURN
      DO 1 I=1,KK
      DO 1 J=1,KK
   1  A(I,J)=0.
      M1=MR-MA
C
C   K1 :
C
      IF(K1.GT.0) GO TO 10
C
C   K1.LE.0 :
C
      K1=IABS(K1)
      GO TO 20
  10  IF(K1.EQ.M1) GO TO 15
C
C   0.LT.K1.LT.M1 :
C
      DO 11 I=1,M1
      WK(I)=ABS(RY(I))
  11  NWK(I)=0
      DO 12 I=1,K1
      CALL MAX(WK,M1,WKMAX,MIND)
      NWK(MIND)=1
  12  WK(MIND)=0.
      I1=1
      DO 13 I=1,M1
      IF(NWK(I).EQ.0) GO TO 13
      NAME(K1-I1+1)=I
      I1=I1+1
  13  CONTINUE
      GO TO 20
C
C   K1.EQ.M1 :
C
  15  DO 16 I=1,K1
  16  NAME(I)=K1-I+1
C
C   K2 :
C
  20  IF(K2.GT.0) GO TO 30
C
C   K2.LE.0 :
C
      K2=IABS(K2)
      GO TO 40
  30  IF(K2.EQ.M1) GO TO 35
C
C   0.LT.K2.LT.M1 :
C
      DO 31 I=1,M1
      WK(I)=ABS(REY(I))
  31  NWK(I)=0
      DO 32 I=1,K2
      CALL MAX(WK,M1,WKMAX,MIND)
      NWK(MIND)=1
  32  WK(MIND)=0.
      I1=1
      DO 33 I=1,M1
      IF(NWK(I).EQ.0) GO TO 33
      NAME(K1+K2-I1+1)=-I
      I1=I1+1
  33  CONTINUE
      GO TO 40
C
C   K2.EQ.M1 :
C
  35  DO 36 I=1,K2
  36  NAME(K1+I)=I-K2-1
C
C
C   FORM MATRIX A :
C
C
  40  IF((K1.NE.0).AND.(K2.NE.0)) GO TO 50
      IF(K2.EQ.0) GO TO 45
C
C   PURE MOVING AVERAGE :
C
      DO 41 I=1,K2
      I1=NAME(I)
      DO 41 J=1,I
      I2=IABS(I1-NAME(J))
      IF(I2.EQ.0) GO TO 42
      A(I,J)=REE(I2)
      GO TO 41
  42  A(I,J)=REE0
  41  A(J,I)=A(I,J)
      GO TO 60
C
C   PURE AUTOREGRESSION :
C
  45  DO 46 I=1,K1
      I1=NAME(I)
      DO 46 J=1,I
      I2=IABS(I1-NAME(J))
      IF(I2.EQ.0) GO TO 47
      A(I,J)=RY(I2)
      GO TO 46
  47  A(I,J)=RY0
  46  A(J,I)=A(I,J)
      GO TO 60
C
C   MIXED :
C
  50  DO 54 I=1,KK
      I1=NAME(I)
      DO 54 J=1,I
      J1=NAME(J)
      I2=I1*J1
      I3=IABS(I1)
      I4=IABS(J1)
      I3=IABS(I3-I4)
      IF(I2.GT.0) GO TO 56
C   REY :
      IF(I3.EQ.0) GO TO 55
      A(I,J)=REY(I3)
      IF((I1+J1).GT.0) A(I,J)=RYE(I3)
      GO TO 54
  55  A(I,J)=REY0
      GO TO 54
  56  IF(I1.GT.0) GO TO 58
C   REE :
      IF(I3.EQ.0) GO TO 57
      A(I,J)=REE(I3)
      GO TO 54
  57  A(I,J)=REE0
      GO TO 54
C   RYY :
  58  IF(I3.EQ.0) GO TO 59
      A(I,J)=RY(I3)
      GO TO 54
  59  A(I,J)=RY0
C
C
  54  A(J,I)=A(I,J)
C
  60  N=KK+1
      A(N,N)=RY0
      DO 61 I=1,KK
      I1=NAME(I)
      IF(I1.GT.0) GO TO 63
      I2=-I1
      A(I,N)=REY(I2)
      GO TO 61
  63  A(I,N)=RY(I1)
  61  A(N,I)=A(I,N)
C
C   CALL SELREG :
C
C
      CALL SELREG(A,NDIM,KK,NAME,KTOP,L,p,KP,COF,IND,IER,iptrec)
	rvar=a(kk+1,kk+1)
      IF(KP.EQ.0) GO TO 99
      DO 70 I=1,KP
      IF(IND(I).GT.0) COF(I)=-COF(I)
  70  CONTINUE
  99  CONTINUE
      MR=MOPT
C
C
C
      RETURN
      END
c&ARMXCV
      SUBROUTINE ARMXCV(RYY,RYY0,ALPH,MR,MA,M,REY,REY0,RYE,REE,REE0,
     1IOPT)
C*******************************************************
C
C   SUBROUTINE TO FIND AUTOCOVARIANCES REE(0),...,REE(M)
C   AND CROSS-COVARIANCES REY(0),...,REY(M) OF AN AUTO-
C   REGRESSIVE PROCESS OF ORDER MA WITH COEFFICIENTS
C   ALPH(1),...,ALPH(MA), AND AUTOCOVARIANCES RYY(0),...
C   ,RYY(MR)
C
C      NOTE : M.LE.MR-MA
C
C   INPUT :
C           MR,MA,RYY0,RYY(1),...,RYY(MR),ALPH(1),...,
C           ALPH(MA),M
C           IOPT : 0 MEANS RESULTS ARE PRINTED, 1 THEY ARENT.
C
C   OUTPUT :
C           REY0,REY(1),...,REY(M),REE0,REE(1),REE(M)
C
C   SUBROUTINES CALLED : NONE
C
C*******************************************************
C
      DIMENSION RYY(MR),ALPH(MA),REY(M),REE(M)
      DIMENSION RYE(M)
C
C
      IF(M.LE.MR-MA) GO TO 10
      WRITE(*,1)
   1  FORMAT(10X,'M MUST BE .LE. MR-MA IN ARMXCV')
      RETURN
C
C   FIND REY0,REY(1),...,REY(M) :
C
  10  MAP1=MA+1
      MP1=M+1
      C=RYY0
      DO 20 I=1,MP1
      IF(I.NE.1) C=RYY(I-1)
         DO 25 J=2,MAP1
  25     C=C+ALPH(J-1)*RYY(I+J-2)
      IF(I.EQ.1) GO TO 26
      REY(I-1)=C
      GO TO 20
  26  REY0=C
  20  CONTINUE
      DO 210 I=1,M
         C=RYY(I)
         DO 211 J=1,MA
         J1=IABS(J-I)
         IF(J1.EQ.0) GO TO 212
         C=C+ALPH(J)*RYY(J1)
         GO TO 211
 212     C=C+ALPH(J)*RYY0
 211     CONTINUE
 210  RYE(I)=C
C
C   FIND REE0,REE(1),...,REE(M) :
C
      REE0=REY0
      DO 30 I=1,M
      I1=MIN0(I,MA)
         C=REY(I)
         DO 35 J=1,I1
         IF(I.EQ.J) GO TO 36
         C=C+ALPH(J)*REY(I-J)
         GO TO 35
  36     C=C+ALPH(J)*REY0
  35     CONTINUE
  30  REE(I)=C
C
C
	iopt=1
      IF(IOPT.EQ.1) RETURN
      WRITE(*,50)
  50  FORMAT(///,9X,2X,3HLAG,10X,6HRYY(V),10X,6HREY(V),
     110X,6HRYE(V),10X,6HREE(V),10X,
     16HPVH(V)/10X,89(1H-))
      II=0
      WRITE(*,51) II,RYY0,REY0,REY0,REE0,REY0
  51  FORMAT(10X,I5,2X,5(F14.7,2X))
      PVH=REY0
      DO 52 I=1,M
      PVH=PVH+(REY(I)**2)/REY0
  52  WRITE(*,51) I,RYY(I),REY(I),RYE(I),REE(I),PVH
      RETURN
      END
