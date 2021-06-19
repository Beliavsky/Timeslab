c&CV
c&CORRMA
        subroutine corrma(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to do Wilsons algorithm
c
c   beta=corrma(rho,nq,r0,maxit,del,ier,rvar)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),np,nr).eq.1) go to 99
        call ckint(args(2),nq)
        if(nq.le.0.or.nq.gt.nr.or.nq.gt.1000) then
		call error(args,2,2)
		go to 99
	endif
	if(ickr(args(3),3,r0,2,0.).eq.1) go to 99
	if(icki(args(4),4,maxit,1,1).eq.1) go to 99
	if(ickr(args(5),5,del,2,0.).eq.1) go to 99
c
c
        n1=nstart(np)-1
        do 50 i=1,nq
  50    wk(i)=array(n1+i)*r0
	n1=nq+1
	n2=n1+nq
	n3=n2+nq
	n4=n3+nq
	n5=n4+nq
	n6=n5+1
	n7=n6+nq
	n8=n7+nq
	n9=n8+nq
	call cvmawl(wk,r0,nq,maxit,del,wk(n1),wk(n2),wk(n3),wk(n4),
     1	wk(n5),wk(n6),wk(n7),wk(n8),wk(n9),rvar,ier)
        call ckaddi(args(6),ier,iref)
	if(ier.ne.0) imflag=1
        if(iref.eq.1.or.ier.ne.0) go to 99
        call ckaddr(args(7),rvar,iref)
        if(iref.eq.1) go to 99
        lab='MA coefficients'
        call ckadda(vname,nq,lab,n9,iref)
  99    continue
        return
        end
c&CVMAWL
	subroutine cvmawl(r,r0,nq,maxit,del,t,f,fn,fp,g,d,dn,gam,
     1	beta,rvar,ier)
c******************************************************************
c
c   Subroutine to do Wilson's algorithm for MA processes.
c
c******************************************************************
c
	dimension r(nq),t(nq),f(nq),fn(nq),fp(nq),g(nq),d(nq),dn(nq),
     1	gam(nq),beta(nq)
c
c   starting values:
c
	t0=sqrt(r0)
	do 10 j=1,nq
  10	t(j)=r(j)/t0
c
c   start iterations (t plays the role of tau, g plays role of r):
c
	do 150 it=1,maxit
c
c
c	write(*,11) t0,(t(kk),kk=1,nq)
c 11	format(1x,6f12.6)
		do 20 j=1,nq
  20		f(j)=t(j)/t0
		g0=t0*t0
		do 22 i=1,nq
  22		g0=g0+t(i)*t(i)
		do 24 i=1,nq
			c=t0*t(i)
			if(i.eq.nq) go to 26
			do 25 j=1,nq-i
  25			c=c+t(j)*t(j+i)
  26			g(i)=c
  24		continue
		g0=(g0+r0)/t0
		do 30 j=1,nq
  30		g(j)=(g(j)+r(j))/t0
c
c   work backward:
c
		do 70 k=nq,1,-1
		gam(k)=g(k)
		fp(k)=f(k)
		fpk=fp(k)
		fpk2=1.-fpk*fpk
		if(fpk2.le.0.) go to 299
		if(k.eq.1) go to 70
			do 40 j=1,k-1
  40			fn(j)=(f(j)-fpk*f(k-j))/fpk2
			do 50 j=1,k-1
  50			g(j)=g(j)-gam(k)*fn(k-j)
			do 60 j=1,k-1
  60			f(j)=fn(j)
  70		continue
c
c   work forward:
c
		d0=g0/2.
		do 110 k=0,nq-1
		d(k+1)=gam(k+1)
		fpk=fp(k+1)
		fpk2=1.-fpk*fpk
		dn0=(d0-fpk*d(k+1))/fpk2
		dn(k+1)=(d(k+1)-fpk*d0)/fpk2
		if(k.eq.0) go to 90
			do 80 j=1,k
  80			dn(j)=(d(j)-fpk*d(k+1-j))/fpk2
  90		d0=dn0
			do 100 j=1,k+1
 100			d(j)=dn(j)
 110		continue
c
c   check convergence:
c
		eps=abs(t0-d0)/abs(d0)
		do 120 j=1,nq
		e1=abs(t(j)-d(j))/abs(d(j))
		if(e1.gt.eps) eps=e1
 120		continue
		t0=d0
		do 130 j=1,nq
 130		t(j)=d(j)
c
c   yes:
c
		if(eps.lt.del) then
			ier=0
			rvar=t0*t0
			do 140 j=1,nq
 140			beta(j)=t(j)/t0
			go to 99
		endif
c
c   no:
c
 150	continue
c
c   nonconvergence in maxit iterations:
c
	ier=1
	write(*,160) eps
 160 	format(' Nonconvergence in Wilson''s algorithm, eps = ',e12.5)
	go to 99
c
c   partial outside (-1,1):
c
 299	write(*,170)
 170	format(' Nonfeasability encountered in Wilson''s algorithm')
	ier=2
c
c
  99	continue
	return
	end
c&CRARMA
        subroutine crarma(args,nargs,vname)
c********************************************************************
c
c   Subroutine to get ARMA parameters from correlations.
c
c   rvar=CORRARMA(rho,r0,np,nq,maxit,del,ier,alpha,beta)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),nrho,nr).eq.1) go to 99
	if(ickr(args(2),2,r0,2,0.).eq.1) go to 99
	if(icki(args(3),3,np,1,1).eq.1) go to 99
	if(icki(args(4),4,nq,1,1).eq.1) go to 99
	if(np+nq.gt.90) then
		call error(args,4,2)
		go to 99
	endif
	if(nr.lt.np+nq) then
		call error(args,4,2)
		go to 99
	endif
	if(icki(args(5),5,maxit,1,1).eq.1) go to 99
	if(ickr(args(6),6,del,2,0.).eq.1) go to 99
c
c
        ndim=np+nq
	n1=ndim+1
	n2=n1+ndim*ndim
	n3=n2+ndim
	n4=n3+ndim
	n5=n4+ndim
	n6=n5+ndim
	n7=n6+ndim
	n8=n7+ndim
	n9=n8+ndim
	n10=n9+ndim
	n12=n11+ndim
	n13=n12+ndim
	call movxy(wk,array(nstart(nrho)),4*ndim)
        do 70 i=1,ndim
  70    wk(i)=wk(i)*r0
	call cvmx1(wk,r0,np,nq,ndim,maxit,del,iwk,wk(n1),wk(n2),
     1	wk(n3),wk(n4),wk(n5),wk(n6),rvar,ier,
     1	wk(n7),wk(n8),wk(n9),wk(n10),wk(n11),wk(n12),wk(n13))
        call ckaddi(args(7),ier,iref)
        if(iref.eq.1) go to 99
        call ckaddr(vname,rvar,iref)
        if(iref.eq.1) go to 99
        lab='AR part of ARMA coefficients'
        call ckadda(args(8),np,lab,n5,iref)
        if(iref.eq.1) go to 99
        lab='MA part of ARMA coefficients'
        call ckadda(args(9),nq,lab,n6,iref)
  99    continue
        return
        end
c&CVMX1
        subroutine cvmx1(r,r0,np,nq,ndim,maxit,del,ip,al,ra,ry,g,
     1  alpha,beta,rvar,ier,t,f,fn,fp,d,dn,gam)
c*******************************************************************
c
c   Subroutine to find arma paramters from covariances.
c
c   ndim.ge.np+nq
c
c******************************************************************
c
        dimension r(ndim),al(ndim,ndim),ra(ndim),ry(ndim),alpha(ndim),
     1  beta(ndim),ip(ndim),g(ndim),t(ndim),f(ndim),fn(ndim),d(ndim),
     1	dn(ndim),gam(ndim),fp(ndim)
c
c   find alpha :
c
        call hiyw(r,r0,np,nq,ndim,ip,al,alpha,ier)
        if(ier.eq.1) then
		ier=3
		write(*,1)
   1		format(' Singular matrix in high order Yule-Walker')
		go to 99
	endif
c
c   find ra :
c
        call macv(alpha,np,1.0,ra,ra0)
c
c   find ry :
c
        ry0=ra0*r0
        do 12 i=1,np
  12    ry0=ry0+2.*ra(i)*r(i)
        do 13 i=1,nq
                c=ra0*r(i)
                do 14 j=1,np
                if(i-j.eq.0) go to 16
                n1=iabs(i-j)
                c1=r(n1)
                go to 14
  16            c1=r0
  14            c=c+ra(j)*(c1+r(i+j))
  13    ry(i)=c
c
c   find corresponding betas :
c
	call cvmawl(ry,ry0,nq,maxit,del,t,f,fn,fp,g,d,dn,gam,
     1  beta,rvar,ier)
  99	continue
        return
        end
c&HIYW
        subroutine hiyw(r,r0,np,nq,ndim,ip,al,alpha,ier)
c******************************************************************
c
c   Subroutine to solve high order Yule-Walker equations. ndim.ge.np+nq
c
c**********************************************************************
c
        dimension r(ndim),al(ndim,ndim),alpha(np),ip(ndim)
c
        ier=0
        if(np.eq.1) go to 10
c
        do 1 i=1,np
        alpha(i)=-r(nq+i)
                do 1 j=i,np
                if(nq+i-j.eq.0) go to 3
                n1=iabs(nq+i-j)
                al(i,j)=r(n1)
                go to 1
   3            al(i,j)=r0
   1    al(j,i)=r(nq+j-i)
c
        call decomp(np,ndim,al,ip)
        if(ip(np).eq.0) go to 99
        call solv(np,ndim,al,alpha,ip)
        return
  10    alpha(1)=-r(nq+1)/r(nq)
        return
  99    write(*,100)
 100    format(' singularity in HIYW')
        ier=1
        return
        end
c&CVDT
      SUBROUTINE CVDT(R,R0,MDIM,M,dseed,AL,D,z,X,IER)
C*******************************************************
C
C   SUBROUTINE TO PRODUCE A SAMPLE X(1),...,X(M) HAVING
C   COVARIANCE MATRIX TOEPL(R0,R(1),...,R(M)).
C
C   INPUT :
C           MDIM : DIMENSION OF MATRIX AL IN MAIN PROGRAM
C           M,R0,R(1),...,R(M)
C           IR : integer*4 SEED FOR WHTSIM
C
C   OUTPUT :
C           X(1),...,X(M)
C           IER : (1 MEANS COV MATRIX NOT POSITIVE DEF.,
C           0 IS NORMAL RETURN)
C
C   AUXILLIARY :
C           AL,D : FACTORS IN MODIFIED CHOLESKY DECOMP
C           OF TOEPL(1.,CORR(1),...,CORR(M-1))
C
C   SUBROUTINES CALLED : TCHOL,WHTSIM
C
C*******************************************************
C
      DIMENSION R(M),AL(MDIM,MDIM),D(M),X(M),Z(m)
        integer*4 dseed
C
C   FIND SQUARE ROOT OF TOEPL(R0,R(1),...,R(M-1)) :
      DO 1 I=1,M
   1  R(I)=R(I)/R0
      CALL TCHOL(R,MDIM,M,AL,D,IER)
      DO 5 I=1,M
   5  R(I)=R(I)*R0
      IF(IER.EQ.1) GO TO 99
C
C   FIND SAMPLE WITH IDENTITY COV MATRIX :
C
      CALL WHTSIM(M,dseed,Z)
C
C   MULTIPLY BY SQUARE ROOT :
C
      CC=SQRT(R0)
      DO 2 I=1,M
   2  Z(I)=Z(I)*SQRT(D(I))
      DO 3 I=1,M
      C=0.
         DO 4 J=1,I
   4     C=C+AL(I,J)*Z(J)
   3  X(I)=C*CC
      RETURN
  99  WRITE(*,100)
 100  FORMAT(10X,'TOEPL(R0,R(1),...,R(M-1)) NOT POS DEF IN CVDT')
C
      RETURN
      END
c&CORR2
        subroutine corr2(args,nargs)
c*****************************************************************
c
c   Subroutine to process the command CORR2(x,y,n,m,iopt,rx0,ry0,
c   rhoxy0,rhox,rhoy,rhoxy,rhoyx)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        double precision inprod
	integer*2 ickl,ickse
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(ickl(args(2),npy,ny).eq.1) go to 99
	ipte=0
        call ckint(args(3),n)
        call ckint(args(4),m)
        call ckint(args(5),iopt)
        if(n.lt.1.or.n.gt.nx.or.n.gt.ny) ipte=3
        if(m.lt.1.or.m.ge.n) ipte=4
        if(iopt.ne.1.and.iopt.ne.2) ipte=5
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
	if(ickse(2*n+4*m).eq.1) go to 99
        call movxy(wk,array(nstart(npx)),4*n)
        call movxy(wk(n+1),array(nstart(npy)),4*n)
        if(iopt.eq.1) then
                call means(wk,n,1,xbar,1)
                call means(wk(n+1),n,1,xbar,1)
        endif
        on=float(n)
c
c   divide variances by n later
c
        rx0=inprod(wk,wk,n)
        ry0=inprod(wk(n+1),wk(n+1),n)
        srxy=sqrt(rx0*ry0)
        rxy0=inprod(wk,wk(n+1),n)
        n2=2*n
        n3=n2+m
        n4=n2+2*m
        n5=n2+3*m
        do 10 i=1,m
        nmi=n-i
        wk(n2+i)=inprod(wk,wk(1+i),nmi)/rx0
        wk(n3+i)=inprod(wk(n+1),wk(n+1+i),nmi)/ry0
        wk(n4+i)=inprod(wk,wk(n+1+i),nmi)/srxy
  10    wk(n5+i)=inprod(wk(n+1),wk(1+i),nmi)/srxy
        rx0=rx0/on
        ry0=ry0/on
        rxy0=rxy0/on
        rhoxy0=rxy0/sqrt(rx0*ry0)
c
c   put in right places :
c
        call ckaddr(args(6),rx0,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(7),ry0,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(8),rhoxy0,iref)
        if(iref.eq.1) go to 99
        lab='autocorrelation of first series'
        call ckadda(args(9),m,lab,n2+1,iref)
        if(iref.eq.1) go to 99
        lab='autocorrelation of second series'
        call ckadda(args(10),m,lab,n3+1,iref)
        if(iref.eq.1) go to 99
        lab='cross-correlation with + lags'
        call ckadda(args(11),m,lab,n4+1,iref)
        if(iref.eq.1) go to 99
        lab='cross-correlation with - lags'
        call ckadda(args(12),m,lab,n5+1,iref)
c
c
  99    continue
        return
        end
