c&ARMA
c&ARMADT
        subroutine armadt(args,nargs,vname)
c******************************************************************
c
c   Subroutine to simulate ARMA data.
c
c   x=ARMADT(alpha,beta,np,nq,rvar,iseed,n,ncorr,rho,r0,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),na,nla).eq.1) go to 99
	if(ickl(args(2),nb,nlb).eq.1) go to 99
        call ckint(args(3),np)
        call ckint(args(4),nq)
	ipte=0
	if(np.lt.1.or.np.gt.nla) ipte=3
	if(nq.lt.1.or.nq.gt.nlb) ipte=4
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
	if(ickr(args(5),5,rvar,2,0.).eq.1) go to 99
	if(ickr(args(6),6,rseed,1,0.).eq.1) go to 99
	if(icki(args(7),7,n,2,max0(np,nq)).eq.1) go to 99
	if(icki(args(8),8,ncorr,1,max0(np,nq)).eq.1) go to 99
	if(rseed.ne.0.) i4seed=rseed
c
c
        ndim=max0(np,nq)+1
	cws=float(np+nq)+float(ndim)**2+5.*float(ncorr)+2.*float(n)
	if(cws.gt.32767.) then
		if(ickse(maxnum+1).eq.1) go to 99
	endif
	nws=np+nq+ndim*ndim+5*ncorr+2*n
	if(ickse(nws).eq.1) go to 99
	nd1=max0(ndim,ncorr+1)
        nsa=nstart(na)
        nsb=nstart(nb)
        ndim1=ndim*ndim+1
        ndim2=ndim1+ncorr
        ndim3=ndim2+ncorr
        ndim4=ndim3+n
	n5=ndim4+n
	n6=n5+nd1
	n7=n6+nd1
        call mxdt(array(nsa),array(nsb),rvar,np,nq,n,ndim,ncorr,i4seed,
     1  wk,wk(ndim1),iwk,wk(ndim2),r0,wk(ndim3),wk(ndim4),ier,wk(n5),
     1	wk(n6),wk(n7))
c
c
        call ckaddi(args(11),ier,iref)
        if(ier.ne.0.or.iref.eq.1) go to 99
        do 60 i=1,ncorr
  60    wk(ndim2-1+i)=wk(ndim2-1+i)/r0
        lab='ARMA correlations'
        call ckadda(args(9),ncorr,lab,ndim2,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(10),r0,iref)
        if(iref.eq.1) go to 99
        lab='ARMA data'
        call ckadda(vname,n,lab,ndim4,iref)
  99    continue
        return
        end
c&MXDT
      SUBROUTINE MXDT(ALPH,BETA,rvar,NP,NQ,N,NDIM,NR,dseed,
     1AL,D,ip,R,R0,E,X,IER,acf,cvli,alpha)
C*******************************************************
C
C   SUBROUTINE TO SIMULATE DATA X(1),...,X(N) FROM A
C   MIXED PROCESS OF ORDER (NP,NQ) WITH PARAMETERS
C   ALPH(1),...,ALPH(NP),BETA(1),...,BETA(NQ), AND rvar
C   (RES VAR)
C
C   INPUT :
C           NP,NQ,ALPH(1),...,ALPH(NP),BETA(1),...,BETA(NQ),
C           rvar,N
C           NDIM : DIMENSION OF AL IN MAIN PROGRAM
C           dseed : integer*4 SEED FOR WHITE NOISE GENERATOR
C           NR : NUMBER OF AUTOCOVARIANCES TO BE FOUND
C           (NR.GE.MAX(NP,NQ))
C
C   OUTPUT :
C           R0,R(1),...,R(NR)
C           X(1),...,X(N)
C           IER : ERROR INDICATOR (1 MEANS A MATIX IS SINGULAR,
C           0 IS NORMAL RETURN)
C           E(1),...,E(N) : SIMULATED DATA FOR MOVING
C           AVERAGE PART
C
C   SUBROUTINES CALLED :MXCV,DECOMP,SOLV,WHTSIM,RANDU,CVDT
C                        TCHOL,MADT
C
C*******************************************************
C
      DIMENSION ALPH(NP),BETA(NQ),AL(NDIM,NDIM),R(NR),D(NR),
     1E(N),X(N),ip(nr),acf(1),cvli(1),alpha(1)
        integer*4 dseed
C
C   FIND X(1),...,X(NP) :
C
	do 310 i=1,np
 310	alph(i)=-alph(i)
	do 320 i=1,nq
 320	beta(i)=-beta(i)
	ma=nr+1
	mxpq=max0(np,nq)
	mxpqp1=mxpq+1
	call wilson(alph,np,beta,nq,acf,ma,cvli,
     1	mxpqp1,alpha,mxpq,ier)
	do 330 i=1,np
 330	alph(i)=-alph(i)
	do 340 j=1,nq
 340	beta(j)=-beta(j)
	do 350 i=1,nr
 350	r(i)=acf(i+1)*rvar
	r0=acf(1)*rvar
	if(ier.ne.0) return
c     CALL MXCV(ALPH,BETA,rvar,NP,NQ,NDIM,NR,AL,d,ip,R,R0,IER)
c     IF(IER.EQ.1) RETURN
      CALL CVDT(R,R0,NDIM,NP,dseed,AL,D,e,X,IER)
      IF(IER.EQ.1) RETURN
C
C   FIND X(NP+1),...,X(N) :
C
      CALL MADT1(BETA,rvar,NQ,N,dseed,e)
      NPP1=NP+1
      DO 1 I=NPP1,N
         C=E(I)
         DO 2 J=1,NP
   2     C=C-ALPH(J)*X(I-J)
   1  X(I)=C
C
      RETURN
      END
c&ARMACR
        subroutine armacr(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to calculate arma correlations
c
c   corr=ARMACORR(alpha,beta,nordp,nordq,rvar,ncorr,r0,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	if(ickl(args(1),np1,na).eq.1) go to 99
	if(ickl(args(2),np2,nb).eq.1) go to 99
        call ckint(args(3),nordp)
        call ckint(args(4),nordq)
	ipte=0
	if(nordp.lt.1.or.nordp.gt.na) ipte=3
	if(nordq.lt.1.or.nordq.gt.nb) ipte=4
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
	if(ickr(args(5),5,rvar,2,0.).eq.1) go to 99
	if(icki(args(6),6,ncorr,1,max0(nordp,nordq)).eq.1) go to 99
c
c
	ns1=nstart(np1)
	ns2=nstart(np2)
	if(ickse(ncorr).eq.1) go to 99
	nws=nordp+nordq+ncorr+1+2*(max0(nordp,nordq)+1)
	if(ickse(nws).eq.1) go to 99
	call movxy(wk,array(ns1),4*nordp)
	call movxy(wk(nordp+1),array(ns2),4*nordq)
	do 310 i=1,nordp
 310	wk(i)=-wk(i)
	do 320 i=1,nordq
 320	wk(nordp+i)=-wk(nordp+i)
	mxpq=max0(nordp,nordq)
	mxpqp1=mxpq+1
	ma=ncorr+1
	n1=nordp+nordq+1
	n2=n1+ma
	n3=n2+mxpqp1
	call wilson(wk,nordp,wk(nordp+1),nordq,wk(n1),ma,wk(n2),
     1	mxpqp1,wk(n3),mxpq,ier)
        call ckaddi(args(8),ier,iref)
        if(ier.ne.0.or.iref.eq.1) go to 99
	r0=wk(n1)
        do 70 i=1,ncorr
  70    wk(n1+i)=wk(n1+i)/r0
	r0=r0*rvar
        call ckaddr(args(7),r0,iref)
        if(iref.eq.1) go to 99
        lab='ARMA correlations'
        call ckadda(vname,ncorr,lab,n1+1,iref)
c
  99    continue
        return
        end
