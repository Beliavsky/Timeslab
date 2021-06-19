c&DTFORE
c&DTFORE
        subroutine dtfore(args,nargs)
c********************************************************************
c
c   Subroutine to process the command
c
c   DTFORE(y,n,mlag1,mlag2,npreds,lag1,lag2,alpha1,alpha2,ypred)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        double precision inprod
	integer*2 ickl,icki,ickse
c
	if(ickl(args(1),npy,nn).eq.1) go to 99
c
        call ckint(args(2),n)
        if(n.lt.1.or.n.gt.nn) then
                call error(args,2,2)
                go to 99
        endif
c
        call ckint(args(3),mlag1)
        if(iabs(mlag1).ge.n) then
                call error(args,3,2)
                go to 99
        endif
        call ckint(args(4),mlag2)
        if(iabs(mlag1)+iabs(mlag2).ge.n) then
                call error(args,4,2)
                go to 99
        endif
c
	if(icki(args(5),5,npreds,1,0).eq.1) go to 99
	nws=4*n+5*iabs(ml2)+2*iabs(ml1)+2*npreds
	if(ickse(nws).eq.1) go to 99
c
c   find alpha1 :
c
        call movxy(wk,array(nstart(npy)),4*n)
        if(mlag1.ne.0) go to 210
        lag1=0
        alpha1=0.0
        go to 220
 210    if(mlag1.gt.0) go to 215
        lag1=-mlag1
        go to 217
 215    continue
        do 10 i=1,mlag1
        c1=-inprod(wk,wk(i+1),n-i)
        c2=inprod(wk,wk,n-i)
  10    wk(n+i)=c1/c2
        do 11 i=1,mlag1
  11    wk(n+i)=abs(wk(n+i))
        call max(wk(n+1),mlag1,xmax,lag1)
 217    c1=-inprod(wk,wk(lag1+1),n-lag1)
        c2=inprod(wk,wk,n-lag1)
        alpha1=c1/c2
 220    continue
c
c   residuals and cat AR for residuals :
c
        n1=n-lag1
        do 20 i=1,n1
  20    wk(n+i)=wk(lag1+i)+alpha1*wk(i)
        on=n1
        n2=n+n1+1
        r0=inprod(wk(n+1),wk(n+1),n1)/on
        if(mlag2.eq.0) then
                lag2=0
                n3=n2
                n4=n3
                go to 310
        endif
        iopto=0
        if(mlag2.lt.0) then
                mlag2=-mlag2
                iopto=1
        endif
        n3=n2+mlag2
        n4=n3+mlag2
        do 30 i=1,mlag2
  30    wk(n2+i-1)=inprod(wk(n+1),wk(n+1+i),n1-i)/on
        call cvar(wk(n2),r0,iopto,n1,mlag2,lag2,wk(n3),wk(n4),rvar)
 310    continue
        lagt=lag1+lag2
        n5=n4+mlag2
        n6=n5+lagt
        if(lagt.ne.0) go to 320
        call movct(wk(n6),4*(n+npreds),char(0))
        go to 100
 320    continue
        call movct(wk(n5),4*lagt,char(0))
        wk(n5+lag1-1)=alpha1
        if(lag2.eq.0) go to 50
        do 40 i=0,lag2-1
        wk(n5+i)=wk(n5+i)+wk(n4+i)
  40    wk(n5+i+lag1)=wk(n5+i+lag1)+alpha1*wk(n4+i)
c
c   ypred (1st lagt values are original y's, then come 1 step predictors
c   then come 1 thru npreds step ahead predictors) :
c
c   n6 is 1st element of ypred while n6+n is location of 1st future
c   predictor :
c
  50    n6=n5+lagt
        do 60 i=1,lagt
  60    wk(n6+i-1)=wk(i)
 330    continue
        do 70 i=lagt+1,n
        c=0.0
                do 80 j=1,lagt
  80            c=c-wk(n5+j-1)*wk(i-j)
  70    wk(n6+i-1)=c
        if(npreds.eq.0) go to 100
c
c   n7 is first location past ypred :
c
        n7=n6+n+npreds
c
c   we'll get future predictors by putting last lagt elements of y
c   into wk(n7),...,wk(n7+lagt-1) and then run the filter off the
c   end of these locations
c
        do 85 i=1,lagt
  85    wk(n7+i-1)=wk(n-lagt+i)
c
c   n8+1 is location of first future predictor
c
        n8=n7+lagt-1
        do 90 i=1,npreds
        c=0.0
                do 87 j=1,lagt
  87            c=c-wk(n5+j-1)*wk(n8+i-j)
  90    wk(n8+i)=c
        call movxy(wk(n6+n),wk(n8+1),4*npreds)
 100    continue
c
c   put things in right place :
c
        call ckaddi(args(6),lag1,iref)
        if(iref.eq.1) go to 99
        if(lag1.ne.0) call ckaddr(args(8),alpha1,iref)
        if(iref.eq.1) go to 99
        call ckaddi(args(7),lag2,iref)
        if(iref.eq.1) go to 99
        if(lag2.eq.0) go to 110
        lab='AR coefficients'
        call ckadda(args(9),lag2,lab,n4,iref)
        if(iref.eq.1) go to 99
 110    continue
        lab='fitted values and forecasts'
        call ckadda(args(10),n+npreds,lab,n6,iref)
c
c
  99    continue
        return
        end
c&CORRAR
        subroutine corrar(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to process the command CORRAR.
c
c   alpha=CORRAR(rho,r0,nord,rvar)
c
c   alpha=CORRAR(rho,r0,maxord,n,nord,rvar,cat)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr,ickse
c
c
	iopt=1
	n=1
	if(ickl(args(1),np,na).eq.1) go to 99
	if(ickr(args(2),2,r0,2,0.).eq.1) go to 99
        call ckint(args(3),nord)
        if(nord.le.0.or.nord.gt.na) then
		call error(args,3,2)
		go to 99
	endif
	if(nargs.eq.7) then
		if(icki(args(4),4,n,2,0).eq.1) go to 99
		iopt=0
        endif
c
c   do calculations and store output :
c
	nws=3*nord
	if(ickse(nws).eq.1) go to 99
        n1=nstart(np)-1
        do 40 i=1,nord
  40    wk(i)=array(n1+i)*r0
        call cvar(wk,r0,iopt,n,nord,nn,wk(nord+1),wk(2*nord+1),rvar)
	if(iopt.eq.0) then
		call ckaddi(args(5),nn,iref)
		if(iref.eq.1) go to 99
		lab='CAT Criterion'
		call ckadda(args(7),nord,lab,nord+1,iref)
		if(nn.eq.0.or.iref.eq.1) go to 99
		call ckaddr(args(6),rvar,iref)
		lab='autoregressive coefficients'
		call ckadda(vname,nn,lab,2*nord+1,iref)
		go to 99
	endif
        lab='autoregressive coefficients'
        call ckadda(vname,nord,lab,2*nord+1,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(4),rvar,iref)
c
c
  99    continue
        return
        end
c&CVAR
      SUBROUTINE CVAR(r,r0,iopt,n,m,nord,wk,alpha,rvar)
C*************************************************                      
C                                                                       
c   For iopt=1 :
c   Subroutine to find alpha and rvar for an AR(alpha,m,rvar) given
c   covariances r0,r(1),...,r(m) (n is ignord and nord=m)
C                                                                       
c   For iopt=0 :
c   Find cat order nord.le.m as well as coefficients and rvar.
c
C*************************************************                      
C                                                                       
      DIMENSION r(m),alpha(m),wk(m),wk1(100),wk2(100)
        if(iopt.eq.1) nord=m
        if(iopt.eq.1) go to 5
        cat=-1.
        nord=0
        catmin=-1.
        on=n
        on1=1.-on
        on2=on*on
   5    continue
        rvar=r0
        do 6 i=1,m
   6    r(i)=r(i)/r0
c
c
        do 80 i=1,m
        c1=r(i)
        c2=1.
        if(i.eq.1) go to 20
                im1=i-1
                do 10 j=1,im1
                c1=c1+wk(j)*r(i-j)
  10            c2=c2+wk(j)*r(j)
  20    continue
        c=-c1/c2
        alpha(i)=c
        rvar=rvar*(1.-c*c)
        if(i.eq.1) go to 40
                do 30 j=1,im1
  30            alpha(j)=wk(j)+c*wk(i-j)
  40    continue
                do 50 j=1,i
  50            wk(j)=alpha(j)
        if(iopt.eq.1) go to 80
                rvar1=rvar/r0
                oi=float(i)
                x1=((on-oi)*on1)/(on2*(1.-c*c))
                x2=(on-oi+1.)/on
                cat=cat+((x1+x2)/rvar1)
		wk2(i)=cat
                if(cat.ge.catmin) go to 80
                catmin=cat
                nord=i
                rvarc=rvar
                        do 55 j=1,i
  55                    wk1(j)=alpha(j)
  80    continue
	if(iopt.eq.0) call movxy(wk,wk2,4*m)
c
c
        do 99 i=1,m
  99    r(i)=r(i)*r0
        if(iopt.eq.1.or.nord.eq.0) return
        do 100 i=1,nord
 100    alpha(i)=wk1(i)
        return
        end
