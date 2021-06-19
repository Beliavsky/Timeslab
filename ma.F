c&MA
c&MACV
      SUBROUTINE MACV(BETA,nq,sig,R,R0)
C*******************************************************
C
C   SUBROUTINE TO CALCULATE THE AUTOCOVARIANCES R0,R(1),
C   ...,R(NQ) FOR A MOVING AVERAGE PROCESS OF ORDER NQ
C   WITH PARAMETERS BETA(1),...,BETA(NQ), AND SIG (RES VAR)
C
C   INPUT :
C           NQ,BETA(1),...,BETA(NQ),SIG
C
C   OUTPUT :
C           R0,R(1),...,R(NQ)
C
C   SUBROUTINES CALLED : NONE
C
C*******************************************************
C
      DIMENSION BETA(NQ),R(NQ)
C
      C=1.
      DO 1 I=1,NQ
   1  C=C+BETA(I)*BETA(I)
      R0=C*SIG
C
      DO 2 I=1,NQ
         C=BETA(I)
         IF(I.EQ.NQ) GO TO 2
         NQMI=NQ-I
         DO 3 J=1,NQMI
   3     C=C+BETA(J)*BETA(J+I)
   2  R(I)=C*SIG
C
      RETURN
      END
c&MACORR
        subroutine macorr(args,nargs,vname)
c******************************************************************
c
c   Subroutine to find correlations of ma
c
c   rho=macorr(beta,nq,rvar,ncorr,r0)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),np,nb).eq.1) go to 99
        call ckint(args(2),nq)
        if(nq.le.0.or.nq.gt.nb) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(icki(args(4),4,ncorr,1,nq).eq.1) go to 99
	if(ickse(ncorr).eq.1) go to 99
c
c
        ns=nstart(np)
        call movct(wk,4*ncorr,char(0))
        call macv(array(ns),nq,rvar,wk,r0)
        do 50 i=1,nq
  50    wk(i)=wk(i)/r0
        call ckaddr(args(5),r0,iref)
        if(iref.eq.1) go to 99
        lab='MA correlations'
        call ckadda(vname,ncorr,lab,1,iref)
  99    continue
        return
        end
c&MADT
        subroutine madt(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to generate data from an ma
c
c   x=MADT(beta,nord,rvar,seed,n)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),np,nb).eq.1) go to 99
        call ckint(args(2),nord)
        if(nord.le.0.or.nord.gt.nb) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,rvar,2,0.).eq.1) go to 99
	if(ickr(args(4),4,rseed,1,0.).eq.1) go to 99
        if(rseed.ne.0.) i4seed=rseed
	if(icki(args(5),5,n,2,nord).eq.1) go to 99
	if(ickse(n).eq.1) go to 99
c
c
        call madt1(array(nstart(np)),rvar,nord,n,i4seed,wk)
c
c
        lab='moving average process'
        call ckadda(vname,n,lab,1,iref)
c
c
  99    continue
        return
        end
c&MADT1
      SUBROUTINE MADT1(BETA,rvar,NQ,N,dseed,x)
C*******************************************************
C
C   SUBROUTINE TO OBTAIN A SAMPLE X(1),...,X(N) FROM A
C   MOVING AVERAGE PROCESS OF ORDER NQ WITH PARAMETERS
C   BETA(1),...,BETA(NQ), AND rvar (RES VAR)
C
C   INPUT :
C           N,NQ,BETA(1),...,BETA(NQ),rvar
C           dseed : integer*4 SEED FOR WHTSIM
C
C   OUTPUT :
C           X(1),...,X(N)
C
C   SUBROUTINES CALLED : WHTSIM,unif2
C
C*******************************************************
C
      DIMENSION BETA(NQ),X(1)
        integer*4 dseed
C
C   SIMULATE WHITE NOISE :
C
      CALL WHTSIM(N+nq,dseed,X)
        sig=sqrt(rvar)
        do 10 i=1,n+nq
  10    x(i)=sig*x(i)
        do 20 i=1,n
        ii=i+nq
        c=x(ii)
        do 15 j=1,nq
  15    c=c+beta(j)*x(ii-j)
  20    x(i)=c
C
      RETURN
      END
