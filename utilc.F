c&UTILC
c&CHICDF
      SUBROUTINE CHICDF(X,N,F,PDF)
      DOUBLE PRECISION C,E,F,G,GAMA,PDF,S,S1,X
      S=1.0D0
      C=1.0D0
      E=N
    1 S1=S
      E=E+2.0D0
      C=C*X/E
      S=S+C
      IF(DABS(S1-S).GT.1.0D-08) GO TO 1
      E=N/2.0D0
      G=GAMA(E)
      F=DEXP(E*DLOG(X/2.0D0)+DLOG(S)-(DLOG(E)+G+X/2.0D0))
      PDF=DEXP((E-1.0D0)*DLOG(X)-(E*DLOG(2.0D0)+G+X/2.0D0))
      IF(PDF.LT.1.0D-04) PDF=1.0D-04
      RETURN
      END
c&DECOMP
      SUBROUTINE DECOMP(N,NDIM,A,IP)
C...........................................................
C
C   MATRIX TRIANGULARIZATION BY GAUSSIAN ELIMINATION.
C
C   INPUT....
C     N = ORDER OF MATRIX
C     NDIM = DECLARED DIMENSION OF ARRAY A
C     A = MATRIX TO BE TRIANGULARIZED
C
C   OUTPUT....
C     A(I,J), I.LE.J = UPPER TRIANGULAR FACTOR, U.
C     A(I,J),I.GT.J = MULTIPLIERS = LOWER TRIANGULAR
C                     FACTOR, I-L
C     IP(K), K.LT.N = INDEX OF K-TH PIVOT ROW
C     IP(N) = (-1)**(NUMBER OF INTERCHANGES) OR 0
C
C   USE SUBROUTINE  SOLVE  TO OBTAIN SOLUTION OF LINEAR SYSTEM.
C   DETERM(A) = IP(N)*A(1,1)*A(2,2)*...*A(N,N)
C   IF IP(N)=0, A IS SINGULAR,  SOLVE  WILL DIVIDE BY ZERO
C   INTERCHANGES FINISHED IN U, ONLY PARTLY IN L.
C
C   REFERENCE...ALGORITHM 423 'LINEAR EQUATION SOLVER'
C               BY CLEVE B. MOLER, CACM, APRIL 1972 P. 274
C
C...........................................................
      DIMENSION A(NDIM,NDIM),IP(NDIM)
      IP(N)=1
      DO 6 K=1,N
      IF(K.EQ.N) GO TO 5
      KP1=K+1
      M=K
      DO 1 I=KP1,N
      IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I
    1 CONTINUE
      IP(K)=M
      IF(M.NE.K) IP(N)=-IP(N)
      T=A(M,K)
      A(M,K)=A(K,K)
      A(K,K)=T
      IF(T.EQ.0.) GO TO 5
      DO 2 I=KP1,N
    2 A(I,K)=-A(I,K)/T
      DO 4 J=KP1,N
      T=A(M,J)
      A(M,J)=A(K,J)
      A(K,J)=T
      IF(T.EQ.0.) GO TO 4
      DO 3 I=KP1,N
    3 A(I,J)=A(I,J)+A(I,K)*T
    4 CONTINUE
    5 IF(A(K,K).EQ.0.) IP(N)=0
    6 CONTINUE
      RETURN
      END
c&GAMA
      FUNCTION GAMA(X)
c*******************************************************************
c
c
c********************************************************************
      DOUBLE PRECISION GAMA,S,T,T2,T4,X
      S=0.0D0
      T=X
      IF(X.GT.10.0D0) GO TO 1
      T=X+10.0D0
      DO 2 I=1,10
      T=T-1.0D0
    2 S=S+DLOG(T)
      T=X+10.0D0
    1 T2=T*T
      T4=T2*T2
      GAMA=(9.189385332D-01)+(T-0.5D0)*DLOG(T)-T+(1.0D0/(12.0D0*T))-(1.0
     1D0/(360.0D0*T*T2))+(1.0D0/(1260.0D0*T*T4))-(1.0D0/(1680.0D0*T*T2*T
     24))-S
      RETURN
      END
c&SOLV
      SUBROUTINE SOLV(N,NDIM,A,B,IP)
C******************************************************************
C
C   SOLUTION OF LINEAR SYSTEM  A*X=B
C
C   DO NOT USE IF  DECOMP  HAS SET IP(N)=0
C
C   INPUT....
C     N = ORDER OF MATRIX A
C     NDIM = DECLARED DIMENSION OF ARRAY A
C     A = TRIANGULARIZED MATRIX OBTAINED FROM SUBROUTINE |DECOMP|
C     B = RIGHT HAND SIDE VECTOR
C     IP = PIVOT VECTOR OBTAINED FROM |DECOMP|
C
C   OUTPUT....
C     B = SOLUTION VECTOR, X.
C
C   REFERENCE...ALGORITHM 423 'LINEAR EQUATION SOLVER'
C               BY CLEVE B. MOLER, CACM, APRIL 1972 P. 274
C
C******************************************************************
      DIMENSION A(NDIM,NDIM),B(NDIM),IP(NDIM)
      IF(N.EQ.1) GO TO 9
      NM1=N-1
      DO 7 K=1,NM1
      KP1=K+1
      M=IP(K)
      T=B(M)
      B(M)=B(K)
      B(K)=T
      DO 7 I=KP1,N
    7 B(I)=B(I)+A(I,K)*T
      DO 8 KB=1,NM1
      KM1=N-KB
      K=KM1+1
      B(K)=B(K)/A(K,K)
      T=-B(K)
      DO 8 I=1,KM1
    8 B(I)=B(I)+A(I,K)*T
    9 B(1)=B(1)/A(1,1)
      RETURN
      END
c&TIMEF
        function timef()
c******************************************************************
c
c   Function to find the elapsed time since midnight in seconds.
c
c********************************************************************
c
      integer*4 timea,it
      it=timea()
      ctpsec=1193180./65536.
      timef=float(it)/ctpsec
      return
      end
c&TCHOL
      SUBROUTINE TCHOL(R,MDIM,M,AL,D,IER)
C*******************************************************
C
C   SUBROUTINE TO OBTAIN THE MODIFIED CHOLESKY DECOMPOSITION
C   OF A SYMMETRIC POSITIVE DEFINITE TOEPLITZ
C   MXM MATRIX A= TOEPL(1.0,R(1),........R(M-1)  )
C
C   INPUT :
C           MDIM : DIMENSION OF AL IN MAIN PROGRAM
C           M: # OF COVARIANCES + 1
C           R: VECTOR CONTAINING 2ND THRU MTH ELEMENTS OF 1ST ROW OF A
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
      DIMENSION R(M),AL(MDIM,MDIM),D(M)
C
      IER=0
      K=1
      D(1)=1.
      AL(1,1)=1.
      IF(M.EQ.1) RETURN
C
      DO 1 I=2,M
      AL(1,I)=R(I-1)
   1  AL(I,1)=AL(1,I)
C
      DO 2 K=2,M
      C=1.
      KM1=K-1
         DO 3 J=1,KM1
   3     C=C-AL(J,K)*AL(K,J)
      D(K)=C
      AL(K,K)=1.
      IF(D(K).LT.1.E-10) GO TO 10
      IF(K.EQ.M) GO TO 6
      KP1=K+1
         DO 4 J=KP1,M
         C=R(J-K)
            DO 5 I=1,KM1
   5        C=C-AL(I,K)*AL(J,I)
         AL(K,J)=C
   4     AL(J,K)=AL(K,J)/D(K)
         GO TO 2
  10     DO 11 J=KP1,M
         AL(K,J)=0.
  11     AL(J,K)=0.
      D(K)=0.
      IER=1
      WRITE(*,12) K
  12  FORMAT(10X,'COLUMN',I4,'IS DEPENDENT ON THE PREVIOUS ONES')
   2  CONTINUE
   6  DO 7 I=1,M
         DO 7 J=I,M
   7  AL(I,J)=AL(J,I)
C
      RETURN
      END
c&XBARS2
      subroutine xbars2(x,n,xbar,s2)
c*****************************************************************
c
c   subroutine to calculate the sample mean (xbar) and variance (s2)
c   of the numbers x(1),...,x(n).
c
c********************************************************************
c
      dimension x(n)
      double precision c,xx
      c=0.0d0
      do 10 i=1,n
  10  c=c+x(i)
      on=n
      on1=on-1
      xbar=c/dble(on)
      c=0.0d0
      do 20 i=1,n
      xx=dble(x(i)-xbar)
  20  c=c+xx*xx
      s2=c/dble(on1)
      return
      end
