c&DIST
c&DIST
        subroutine dist(args,nargs,vname)
c******************************************************************
c
c   Subroutine to process the command x=DIST(name,iopt,n,x) or
c   x=DIST(name,iopt,n,x,df) or x=DIST(name,iopt,n,x,df1,df2).
c
c   iopt=1 means pdf, 2 means cdf, 3 means quantile.
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15,lck1*15,lck2(15)*1,lck*1
        double precision z,z1,z2,z3,z4,z5,f,s2pi,gama
	integer*2 icki,ickr,ickse
        equivalence (lck1,lck2)
c
c   Handle Gamma Function Case:
c
	if(nargs.eq.3) then
		if(icki(args(1),1,n,1,1).eq.1) go to 99
		if(ickr(args(2),2,a,2,0.).eq.1) go to 99
		if(ickr(args(3),3,b,1,a).eq.1) go to 99
		del=0.0
		if(n.gt.1) del=(b-a)/float(n-1)
		if(ickse(n).eq.1) go to 99
		do 5 i=1,n
		z=dble(a+float(i-1)*del)
		z1=gama(z)
		z1=dexp(z1)
   5		wk(i)=z1
		lab='Gamma Function'
		if(n.eq.1) call ckaddr(vname,wk(1),iref)
		if(n.gt.1) call ckadda(vname,n,lab,1,iref)
		go to 99
	endif
c
c   Non Gamma Function Case:
c
        lck1=args(1)
        call ltoup(lck2,1)
        lck=lck2(1)
        if(lck.ne.'Z'.and.lck.ne.'T'.and.lck.ne.'F'.and.lck.ne.'C') then
	        call error(args,1,2)
        	go to 99
	endif
        call ckint(args(2),iopt)
        if(iopt.lt.1.or.iopt.gt.3) then
                call error(args,2,2)
                go to 99
        endif
	if(icki(args(3),3,n,1,1).eq.1) go to 99
	if(ickse(n).eq.1) go to 99
	if(n.eq.1) then
		if(ickr(args(4),4,x,0,0.).eq.1) go to 99
		wk(1)=x
		if(iopt.eq.3) then
			if(x.le.0..or.x.ge.1.) then
				call error(args,4,2)
				go to 99
			endif
		endif
	endif
	if(n.gt.1) then
		if(ickl(args(4),npx,nx).eq.1) go to 99
		if(n.lt.nx) then
			call error(args,3,2)
			go to 99
		endif
		call movxy(wk,array(nstart(npx)),4*n)
		if(iopt.eq.3) then
			do 1 i=1,n
			if(wk(i).le.0..or.wk(i).ge.1.) then
				call error(args,4,2)
				go to 99
			endif
   1			continue
		endif
	endif
	ipte=0
	if(lck.eq.'Z'.and.nargs.ne.4) ipte=1
	if(lck.eq.'T'.and.nargs.ne.5) ipte=1
	if(lck.eq.'C'.and.nargs.ne.5) ipte=1
	if(lck.eq.'F'.and.nargs.ne.6) ipte=1
	if(ipte.eq.1) then
		write(*,10)
  10		format(' Wrong Number of Arguments in DIST')
		imflag=1
		go to 99
	endif
	if(nargs.eq.5.or.nargs.eq.6) then
		if(icki(args(5),5,nd1,2,0).eq.1) go to 99
	endif
	if(nargs.eq.6) then
		if(icki(args(6),6,nd2,2,0).eq.1) go to 99
	endif
c
c   Z pdf:
c
	if(lck.eq.'Z') then
		if(iopt.eq.1) then
		        s2pi=dble(sqrt(8.*atan(1.)))
		        do 31 i=1,n
		        z=dble(wk(i))
  31		        wk(i)=dexp(-z*z/2.d0)/s2pi
    		        lab='Normal pdf'
		        go to 90
		endif
c
c   Z cdf:
c
		if(iopt.eq.2) then
		        do 33 i=1,n
		        z=dble(wk(i))
		        call norcdf(z,f,z1)
  33		        wk(i)=f
		        lab='Normal cdf'
		        go to 90
		endif
c
c   Z quantile:
c
		if(iopt.eq.3) then
		        do 35 i=1,n
		        z=dble(wk(i))
  		        call norqnt(z,f,iref)
  35 		        wk(i)=f
		        lab='Normal quantile'
		        go to 90
		endif
	endif
c
c   t pdf:
c
	if(lck.eq.'T') then
	        if(iopt.eq.1) then
        		z1=dble(float(nd1)/2.)
		        z2=.5d0
		        z1=gama(z1+z2)-gama(z1)-gama(z2)
		        on=nd1
		        z3=dble(on)
		        ak=float(nd1+1)/2.
		        z4=dble(ak)
		        z1=z1-.5d0*dlog(z3)
		        do 41 i=1,n
		        z2=dble(wk(i))
		        z=z1-z4*dlog(1.d0+((z2*z2)/z3))
  41		        wk(i)=dexp(z)
  		        lab='t density'
			go to 90
		endif
c
c   t cdf:
c
		if(iopt.eq.2) then
		        do 43 i=1,n
		        z=dble(wk(i))
		        call tdist(z,nd1,f,z1,z2,z3,z4)
  43 		        wk(i)=f
		        lab='t cdf'
		        go to 90
		endif
c
c   t quantile :
c
		if(iopt.eq.3) then
		        do 45 i=1,n
		        z=dble(wk(i))
		        call invst(z,nd1,f)
  45 		        wk(i)=f
		        lab='t quantile'
		        go to 90
		endif
	endif
c
c   chi square pdf:
c
        if(lck.eq.'C') then
	        if(iopt.eq.1) then
        		z1=dble(float(nd1)/2.)
		        f=z1-1.d0
		        z=-z1*dlog(2.d0)-gama(z1)
		        do 51 i=1,n
		        z2=dble(wk(i))
		        if(z2.gt.0.) go to 151
		        wk(i)=0.
		        go to 51
 151		        continue
		        z3=z2/2.d0
		        z4=z-z3+f*dlog(z2)
		        wk(i)=dexp(z4)
  51 			continue
		        lab='chi-square pdf'
		        go to 90
		endif
c
c   chi square cdf:
c
		if(iopt.eq.2) then
		        do 53 i=1,n
		        z=dble(wk(i))
		        call chicdf(z,nd1,f,z1)
  53 		        wk(i)=f
		        lab='chi-square cdf'
		        go to 90
		endif
c
c   chi square quantile:
c
		if(iopt.eq.3) then
		        do 55 i=1,n
		        z=dble(wk(i))
		        call chiqnt(z,nd1,f,iref)
  55		        wk(i)=f
  		        lab='chi-square quantiles'
		        go to 90
		endif
	endif
c
c   F pdf:
c
	if(lck.eq.'F') then
		if(iopt.eq.1) then
		        z1=dble(float(nd1)/2.)
		        z2=dble(float(nd2)/2.)
		        z=gama(z1+z2)-gama(z1)-gama(z2)
		        z3=dble(float(nd1)/float(nd2))
		        z=z-z1*dlog(z3)
		        z4=z1+z2
		        z5=z1-1.d0
		        do 61 i=1,n
		        z1=dble(wk(i))
		        if(z1.gt.0.d0) go to 161
		        wk(i)=0.
		        go to 61
 161 		        continue
		        f=z+z5*dlog(z1)
		        f=f-z4*dlog(1.d0+z3*z1)
		        wk(i)=dexp(f)
  61 			continue
		        lab='F pdf'
		        go to 90
		endif
c
c   F cdf:
c
		if(iopt.eq.2) then
		        do 63 i=1,n
		        z=dble(wk(i))
		        call fdist(z,nd1,nd2,f,z1)
  63		        wk(i)=f
		        lab='F cdf'
		        go to 90
		endif
c
c   F quantile:
c
		if(iopt.eq.3) then
		        do 65 i=1,n
		        z=dble(wk(i))
		        call fqntl(z,nd1,nd2,f,iref)
  65			wk(i)=f
		        lab='F quantile'
		endif
	endif
c
c
c
  90    if(n.eq.1) call ckaddr(vname,wk(1),iref)
        if(n.gt.1) call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&INVST
      SUBROUTINE INVST(CDF,N,X)
C-------------------------------------------------------------------
C
C   'INVST' COMPUTES THE QUANTILE FUNCTION OF THE STUDENT T DISTRIBUTION
C   WITH N DEGREES OF FREEDOM.
C
C   INPUT:   CDF (DOUBLE PRECISION)
C            N = DEGREES OF FREEDOM (INTEGER)
C
C   OUTPUT:  X = QUANTILE (DOUBLE PRECISION)
C
C   SUBPROGRAMS CALLED:  TDIST(X,N,CDF,Q,CEN,Q2,Q1)
C                        GAMA(X)
C
C   REFERENCE: ABRAMOWITZ & STEGUN, EQNS: 26.2.22, 26.7.3, 26.7.4,
C                                         26.7.5
C
C-------------------------------------------------------------------
      DOUBLE PRECISION A,C,CDF,CEN,E,GAMA,P,PDF,Q,Q1,Q2,T,X,X1,X2
      IF(CDF.LE.0.0D0.OR.CDF.GE.1.0D0) RETURN
      P=CDF
      IF(CDF.GT.0.5D0) P=1.0D0-CDF
      IF(N.GT.1) GO TO 201
      X1=DTAN(1.57079632679D0*(1.0D0-2.0D0*P))
      GO TO 200
  201 IF(N.GT.2) GO TO 202
      X1=(1.0D0-2.0D0*P)/DSQRT(2.0D0*P*(1.0D0-P))
      GO TO 200
  202 IF(N.GT.3) GO TO 203
      X1=1.57079632679D0
      A=X1*2.0D0*(1.0D0-2.0D0*P)
  204 X2=X1
      X1=(X2*DCOS(X2)-DSIN(X2)+A)/(1.0D0+DCOS(X2))
      IF(DABS(X1-X2).GT.1.0D-07) GO TO 204
      X1=1.7320508076D0*DTAN(X1/2.0D0)
      GO TO 200
  203 T=DSQRT((-2.0D0)*DLOG(P))
      E=N
      X2=T-(2.30753D0+0.27061D0*T)/((T*0.04481D0+0.99229D0)*T+1.0D0)
      X1=X2+((X2*X2+1.0D0)*X2)/(4.0D0*E)
  205 X2=X1
      T= (-0.5D0)*DLOG(3.141592654D0*E)-((E+1.0D0)/2.0D0)*DLOG(1.0D0+X2*
     1X2/E)+GAMA((E+1.0D0)/2.0D0)-GAMA(E/2.0D0)
      PDF=DEXP(T)
      IF(PDF.LT.1.0D-03) PDF=1.0D-03
      CALL TDIST(X2,N,C,Q,CEN,Q2,Q1)
      X1=X2+(Q-P)/PDF
       IF(DABS(X1-X2).GT.1.0D-06) GO TO 205
  200 X=X1
      IF(CDF.LT.0.5D0) X=X1*(-1.0D0)
      RETURN
      END
c&TDIST
      SUBROUTINE TDIST(X,N,CDF,H,CEN,H2,H1)
C------------------------------------------------------------------
C
C   'TDIST' COMPUTES THE DISTRIBUTION FUNCTION FOR THE STUDENT
C   T-DISTRIBUTION WITH N DEGREES OF FREEDOM.
C
C   INPUT:  X = POINT AT WHICH THE CDF IS TO BE EVALUATED
C               (DOUBLE PRECISION)
C           N = DEGREES OF FREEDOM (INTEGER)
C
C   OUTPUT:  CDF, H, CEN, H2, H1 (DOUBLE PRECISION EACH)
C
C   SUBPROGRAMS CALLED:  NONE
C
C   REFERENCE:  ABRAMOWITZ & STEGUN, EQNS:  26.7.3 AND 26.7.4
C
C   NOTE:  NO IFAULT CODES USED.
C
C-------------------------------------------------------------------
      DOUBLE PRECISION C,CDF,CEN,C2,DL,DU,E,H,H1,H2,S,SUM,T,THETA,X
      E=N
      T=DABS(X)
      THETA=DATAN(T/DSQRT(E))
      S=DSIN(THETA)
      C=DCOS(THETA)
      C2=C*C
      L=N-2
      IF(MOD(N,2).EQ.1) GO TO 200
      SUM=1.0D0
      IF(N.EQ.2) GO TO 203
      DU=-1.0D0
      DL=0.0D0
      E=1.0D0
      DO 11 I=2,L,2
      DU=DU+2.0D0
      DL=DL+2.0D0
      E=E*C2*(DU/DL)
   11 SUM=SUM+E
  203 CEN=S*SUM
      GO TO 202
  200 SUM=0.0D0
      IF(N.EQ.1) GO TO 201
      SUM=S*C
      IF(N.EQ.3) GO TO 201
      SUM=C
      E=C
      DU=0.0D0
      DL=1.0D0
      DO 10 I=3,L,2
      DU=DU+2.0D0
      DL=DL+2.0D0
      E=E*C2*(DU/DL)
   10 SUM=SUM+E
      SUM=S*SUM
  201 CEN=(THETA+SUM)*.6366197723D0
  202 H2=1.0D0-CEN
      H1=H2/2.0D0
      H=H1
      IF(X.LT.0.0D0) H=H1+CEN
      CDF=1.0D0-H
      RETURN
      END
c&NORQNT
      SUBROUTINE NORQNT(U,Q,IFAULT)
C-------------------------------------------------------------------
C
C   'NORQNT' COMPUTES THE QUANTILE FUNCTION OF THE STANDARD NORMAL
C   DISTRIBUTION.
C
C   INPUT:  U = POINT AT WHICH THE QUANTILE FUNCTION IS TO BE
C               EVAULUATED. (DOUBLE PRECISION)
C
C   OUTPUT:  Q (DOUBLE PRECISION)
C            IFAULT (INTEGER)
C
C   SUBPROGRAMS CALLED:  NORCDF
C
C   REFERENCE:  ABRAMOWITZ & STEGUN:  EQN. 26.2.22
C
C-------------------------------------------------------------------
      DOUBLE PRECISION D,F,P,PDF,Q,T,TOL,U,X1,X2
      IFAULT=2
      IF((U.LE.0.0D0).OR.(U.GE.1.0D0)) GO TO 99
      IFAULT=1
      TOL=1.0D-05
      IF(U-0.5D0) 201,202,203
  201 P=U
      GO TO 204
  202 Q=0.0D0
      GO TO 99
  203 P=1.0D0-U
  204 T=DSQRT(-2.0D0*DLOG(P))
      X1=T-(2.30753D0+0.27061D0*T)/((T*0.04481D0+0.99229D0)*T+1.0D0)
      IF(U.LT.0.5D0) X1=(-1.0D0)*X1
      ITER=0
  205 X2=X1
      CALL NORCDF(X2,F,PDF)
      IF(PDF.LT.1.0D-03) PDF=1.0D-03
      D=(F-U)/PDF
      X1=X2-D
      ITER=ITER+1
      IF(DABS(D).LT.TOL) GO TO 206
      IF(ITER.LE.20) GO TO 205
      IFAULT=IFAULT*3
  206 Q=X1
   99 RETURN
      END
c&NORCDF
      SUBROUTINE NORCDF(X,F,PDF)
c***************************************************************
c
c*****************************************************************
c
      DOUBLE PRECISION D1,D2,D3,F,PDF,T,W,X,Y
      Y=DABS(X)
      PDF=3.989422804D-01*DEXP(-0.5D0*Y*Y)
      IF(Y.LE.10.0D0) GO TO 200
      IF(X) 201,99,202
  201 F=0.0D0
      GO TO 99
  202 F=1.0D0
      GO TO 99
  200 IF(Y.LE.3.0D0) GO TO 203
      M=-10.0D0+(85.0D0+255.0D0/Y)/Y
      W=Y*Y
      T=W
      DO 10 I=1,M
      J=M-I+1
      D1=4*J-1
      D2=2*J
      D3=2*J+1
   10 W=T+D1-(D2*D3/W)
      F=PDF*((1.0D0/W)-1.0D0)/Y+1.0D0
      GO TO 204
  203 M=10.0D0*(1.0D0+Y)
      W=Y
      DO 11 I=1,M
      D1=2*(M-I+1)+1
   11 W=(W*Y/D1+1.0D0)*Y
      F=PDF*W+0.5D0
  204 IF(X.LT.0.0D0) F=1.0D0-F
   99 RETURN
      END
c&CHIQNT
      SUBROUTINE CHIQNT(U,N,Q,IFAULT)
C-------------------------------------------------------------------
C
C   'CHIQNT' COMPUTES THE QUANTILE FUNCTION OF THE CHI-SQUARE
C   DISTRIBUTION WITH N DEGREES OF FREEDOM.
C
C   INPUT:  U = POINT AT WHICH THE QUANTILE FUNCTION IS TO BE
C               EVALUATED (DOUBLE PRECISION)
C           N = NUMBER OF DEGREES OF FREEDOM
C               (INTEGER)
C
C   OUTPUT:  Q (DOUBLE PRECISION)
C            IFAULT (INTEGER)
C
C   SUBPROGRAMS CALLED:  CHICDF
C                           GAMA
C
C   REFERENCE:  TEXT
C               ABRAMOWITZ & STEGUN, EQNS: 6.1.41, 26.2.22, 26.4.6,
C                                          26.4.17
C
C-------------------------------------------------------------------
      DOUBLE PRECISION BD,E,F,P,PDF,Q,T,U,X,X1,X2,X3,Y1,Y2,Y3
      IF((0.0D0.LT.U).AND.(U.LT.1.0D0)) GO TO 200
      IFAULT=2
      GO TO 99
  200 IFAULT=1
      P=1.0D0-U
      IF(U.LT.0.5D0) P=U
      T=DSQRT(-2.0D0*DLOG(P))
      X=T-(0.27061D0*T+2.30753)/((0.04481*T+0.99229D0)*T+1.0D0)
      IF(U.LT.0.5D0) X=-1.0D0*X
      E=N
      T=2.0D0/(9.0D0*E)
      X2=DABS(E*(1.0D0-T+X*DSQRT(T))**3)
      ITER=0
  201 X1=X2
      ITER=ITER+1
      CALL CHICDF(X1,N,F,PDF)
      X2=X1+(U-F)/PDF
      IF(X2) 202,202,203
  203 IF(ITER.GT.20) GO TO 98
      BD=DABS(X1-X2)/DABS(X1)
      IF(BD.GT.1.0D-05) GO TO 201
      Q=X2
      GO TO 99
  202 IFAULT=IFAULT*3
      X2=X1
      Y2=F
      X1=0.0D0
      Y1=X1
  204 ITER=ITER+1
      IF(ITER.GT.30) GO TO 98
      IF(Y1-Y2) 205,206,205
  206 IFAULT=IFAULT*7
      GO TO 99
  205 X3=X1+(X1-X2)*(U-Y1)/(Y1-Y2)
      IF(X3) 207,210,208
  207 X3=0.0D0
  210 Y3=X3
      GO TO 209
  208 CALL CHICDF(X3,N,Y3,PDF)
  209 X1=X2
      Y1=Y2
      X2=X3
      Y2=Y3
      BD=DABS(X1-X2)/DABS(X1+1.0D-05)
      IF(BD.GT.1.0D-05) GO TO 204
      Q=X2
      GO TO 99
   98 IFAULT=IFAULT*5
   99 RETURN
      END
c&FDIST
      SUBROUTINE FDIST(X,N,M,CDF,H)
C-------------------------------------------------------------------
C
C   'FDIST' COMPUTES THE DISTRIBUTION FUNCTION OF THE F-DISTRIBUTION
C   WITH N AND M DEGREES OF FREEDOM.
C
C   INPUT:  X = POINT AT WHICH THE CDF IS TO BE EVALUATED
C               (DOUBLE PRECISION)
C           N = DEGREES OF FREEDOM IN THE NUMERATOR (INTEGER)
C           M = DEGREES OF FREEDOM IN THE DENOMINATOR (INTEGER)
C
C   OUTPUT:  CDF, H (DOUBLE PRECISION EACH)
C
C   SUBPROGRAMS CALLED: SUMZ(T,K,L,S)
C
C   REFERENCE:  ABRAMOWITZ & STEGUN, EQNS: 26.6.4, 26.6.5, 26.6.8
C
C   NOTE: NO IFAULT CODES
C
C-----------------------------------------------------------------
      DOUBLE PRECISION A,B,BETA,C,CDF,C2,D,DL,DU,E,EM,EN,H,S,SM,S2,T,
     1THETA,X,Y
      CDF=0.0D0
      H=1.0D0
      IF(X.LT.0.0D0) GO TO 99
      EN=N
      EM=M
      IF(MOD(M*N,2).EQ.1) GO TO 200
      Y=EM/(EM+EN*X)
      IF(MOD(N,2).EQ.1) GO TO 201
      IF(MOD(M,2).EQ.1) GO TO 202
      K=MIN0(N,M)
      L=MAX0(N,M)
      IF(K.EQ.N) GO TO 203
      GO TO 204
  201 K=M
      L=N
  204 T=Y
      CALL SUMZ(T,K,L,S)
      E=L
      E=E/2.0D0
      H=1.0D0-((1.0D0-Y)**E)*S
      CDF=1.0D0-H
      GO TO 99
  202 K=N
      L=M
  203 T=1.0D0-Y
      CALL SUMZ(T,K,L,S)
      E=L
      E=E/2.0D0
      H=(Y**E)*S
      CDF=1.0D0-H
      GO TO 99
  200 THETA=DATAN(DSQRT(EN*X/EM))
      S=DSIN(THETA)
      C=DCOS(THETA)
      S2=S*S
      C2=C*C
      SM=0.0D0
      IF(M.EQ.1) GO TO 205
      SM=S*C
      IF(M.EQ.3) GO TO 205
      SM=C
      B=C
      DU=0.0D0
      DL=1.0D0
      L=M-2
      DO 10 I=3,L,2
      DU=DU+2.0D0
      DL=DL+2.0D0
      B=(B*DU*C2)/DL
   10 SM=SM+B
      SM=S*SM
  205 A=(THETA+SM)*0.6366197723D0
      BETA=0.0D0
      D=1.0D0
      IF(N.EQ.1) GO TO 206
      IF(M.EQ.1) GO TO 208
      L=(M-1)/2
      DU=L+1
      DL=EM/2.0D0
      DO 11 I=1,L
      DU=DU-1.0D0
      DL=DL-1.0D0
   11 D=D*DU/DL
  208 D=S*(0.6366197723D0)*(C**M)*D
      SM=1.0D0
      IF(N.EQ.3) GO TO 207
      L=N-3
      DU=M-1
      DL=1.0D0
      B=1.0D0
      DO 12 I=2,L,2
      DU=DU+2.0D0
      DL=DL+2.0D0
      B=(B*DU*S2)/DL
   12 SM=SM+B
  207 BETA=D*SM
  206 H=1.0D0-A+BETA
      CDF=1.0D0-H
   99 RETURN
      END
c&SUMZ
      SUBROUTINE SUMZ(T,K,L,S)
c******************************************************************
c
c*******************************************************************
      DOUBLE PRECISION A,DL,DU,S,T
      S=1.0D0
      J=(K-2)/2
      IF(J.EQ.0) RETURN
      DU=L-2
      DL=0.0D0
      A=1.0D0
      DO 1 I=1,J
      DU=DU+2.0D0
      DL=DL+2.0D0
      A=(DU*A*T)/DL
    1 S=S+A
      RETURN
      END
c&FQNTL
      SUBROUTINE FQNTL(U,N,M,Q,IFAULT)
C-------------------------------------------------------------------
C
C   'FQNTL' COMPUTES THE QUANTILE FUNCTION OF THE F DISTRIBUTION
C   WITH N AND M DEGREES OF FREEDOM.
C
C   INPUT:  U = POINT AT WHICH THE QUANTILE FUNCTION IS TO BE
C               EVALUATED (DOUBLE PRECISION)
C           N = NUMBER OF DEGREES OF FREEDOM IN THE NUMERATOR
C               (INTEGER)
C           M = NUMBER OF DEGREES OF FREEDOM IN THE DENOMINATOR
C               (INTEGER)
C
C   OUTPUT:  Q (DOUBLE PRECISION)
C            IFAULT (INTEGER)
C
C   SUBPROGRAMS CALLED: FDIST
C                          SUMZ
C
C   REFERENCE:  ABRAMOWITZ & STEGUN, EQNS: 26.5.22, 26.2.23
C
C--------------------------------------------------------------------
      DOUBLE PRECISION A,B,D,E,H,P,Q,T,TOL,U,W,X,X1,X2,Y,Y1,Y2
      IFAULT=2
      IF((U.LE.0.0D0).OR.(U.GE.1.0D0)) GO TO 99
      IFAULT=1
      TOL=1.0D-05
      P=1.0D0-U
      IF(U.LE.0.5D0) P=U
      T=DSQRT(-2.0D0*DLOG(P))
      Y=T-(0.27061D0*T+2.30753D0)/((0.04481D0*T+0.9922D0)*T+1.0D0)
      IF(U.LT.0.5D0) Y=-1.0D0*Y
      IF((N.EQ.1).OR.(M.EQ.1)) GO TO 200
      A=M-1
      B=N-1
      H=2.0D0/((1.0D0/A)+(1.0D0/B))
      D=(Y*Y-3.0D0)/6.0D0
      W=(Y*DSQRT(H+D)/H)-((1.0D0/B)-(1.0D0/A))*(D+(5.0D0/6.0D0)-(2.0D0/(
     &3.0D0*H)))
      X1=DEXP(2.0D0*W)
      X2=1.01D0*X1
  205 CALL FDIST(X1,N,M,Y1,h)
      CALL FDIST(X2,N,M,Y2,h)
      ITER=1
  203 X=(U-Y1)*(X1-X2)/(Y1-Y2)+X1
      IF(DABS((X-X1)/X1).LT.TOL) GO TO 201
      ITER=ITER+1
      IF(ITER.GT.20) GO TO 202
      CALL FDIST(X,N,M,Y,h)
      X2=X1
      Y2=Y1
      X1=X
      Y1=Y
      GO TO 203
  201 Q=X
      GO TO 99
  202 IFAULT=IFAULT*3
      GO TO 99
  200 IF(N.EQ.1) GO TO 204
      IFAULT=IFAULT*5
      E=N
      X1=((DTAN(U*1.570796327D0+1.772453851D0))**2)/E
      X2=1.1D0*X1
      GO TO 205
  204 IFAULT=IFAULT*7
      E=M
      Y=DABS(Y)
      T=Y+(Y*(Y*Y+1.0D0)/(4.0D0*E))
      X1=T*T
      X2=1.1D0*X1
      GO TO 205
   99 RETURN
      END
