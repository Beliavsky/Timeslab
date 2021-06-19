c&DTAR
c&DTAR
        subroutine dtar(args,nargs,vname)
c***************************************************************
c
c   Subroutine to find AR parameters for a time series.
c
c   alpha=dtar(x,n,m,iopto,ioptl,ioptm,nord,r0,rvar)
c
c****************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
        double precision inprod
	integer*2 ickl,ickse
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),n)
	call ckint(args(3),m)
	call ckint(args(4),iopto)
	call ckint(args(5),ioptl)
	call ckint(args(6),ioptm)
	ipte=0
	if(ioptm.lt.1.or.ioptm.gt.6) ipte=6
	if(ioptl.lt.1.or.ioptl.gt.2) ipte=5
	if(iopto.lt.1.or.iopto.gt.2) ipte=4
	if(m.le.0.or.m.gt.n) ipte=3
	if(n.lt.1.or.n.gt.nn) ipte=2
	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif
c
c
	nws=2*(n+m)+3*m
	if(ickse(nws).eq.1) go to 99
        npm=n+m
        npm2=2*npm
        do 100 i=1,npm
        wk(i)=0.0
 100    wk(npm+i)=0.0
        n1=nstart(np)-1
        do 110 i=1,n
        c=array(n1+i)
        wk(i)=c
 110    wk(npm+i+1)=c
        r0=inprod(wk(1),wk(1),n)/float(n)
        rvar=r0
        rvar1=rvar
        brv=1.0
        urv=1.0
        cat=-1.0
        on=n
        oneon=1./on
        on1=1.-on
        on2=on*on
        catmin=-(1.+oneon)
        nord=0
        npm3=npm2+m
        npm4=mpm3+m
        if(ioptl.eq.1) write(*,111)
	if(iptrec.eq.1) write(3,111,err=701)
 701	continue
 111    format(5x,'i',5x,'part(i)',6x,'brv(i)',6x,'urv(i)',6x,'cat(i)',
     1  /1x,55(1h-))
c
c   start recursion :
c
        do 120 i=1,m
c
c   get partial :
c
        if(ioptm.eq.2.or.ioptm.eq.4.or.ioptm.eq.6) go to 130
c
c   yule-walker :
c
        top=-2.*inprod(wk(1),wk(npm+1),npm)
        bot=inprod(wk(1),wk(1),npm)+inprod(wk(npm+1),wk(npm+1),npm)
        part=top/bot
        go to 140
c
c   burg :
c
 130    top=-2.*inprod(wk(npm+i+1),wk(i+1),n-i)
        bot=inprod(wk(npm+i+1),wk(npm+i+1),n-i)+
     1      inprod(wk(i+1),wk(i+1),n-i)
        part=top/bot
 140    continue
        if(part.le.-1.) part=-.9999
        if(part.ge.1.) part=.9999
        wk(npm2+i)=part
        rvar1=rvar1*(1.-part*part)
c
c   brv,urv,cat :
c
        onpi=on+i
        onmi=on-i
        brv=brv*(1.-part*part)
        urv1=urv
        urv=(on/onmi)*brv
        cat=cat+(on1/(on*urv))+(1./urv1)
        wk(npm3+i)=cat
        if(cat.ge.catmin) go to 150
        catmin=cat
        rvar=r0*brv
        nord=i
 150    continue
c
c   shift :
c
        do 160 j=1,npm-1
        i1=npm-j+1
        i2=i1-1
        ni1=npm+i1
        ni2=npm+i2
        temp=wk(ni1)
        wk(ni1)=wk(ni2)+part*wk(i2)
 160    wk(i1)=wk(i1)+part*temp
        wk(1)=wk(1)+part*wk(npm+1)
        wk(npm+1)=0.0
c
c   write information if asked for :
c
        if(ioptl.eq.2) go to 120
        if(ioptb.eq.1) go to 162
        if(mod(i,22).eq.0) pause
 162    write(*,161) i,-part,brv,urv,cat
 161    format(1x,i5,4f12.6)
	if(iptrec.eq.1) write(3,161,err=702) i,-part,brv,urv,cat
 702	continue
c
c   end recursion :
c
 120    continue
c
c   put output arguments in proper place :
c
        if(iopto.eq.1) rvar=rvar1
        if(iopto.eq.1) nord=m
        write(*,175) nord
 175    format(' nord = ',i3)
	if(iptrec.eq.1) write(3,175,err=703) nord
 703	continue
        if(nord.eq.0) go to 217
        if(ioptm.le.2) go to 217
        if(ioptm.eq.3.or.ioptm.eq.4) go to 201
c
c   ioptm=5 or 6 :
c
        call subaic(wk(npm2+1),nord,n,r0,iwk,nord1,wk(npm4+1),rvar)
        call movxy(wk(npm2+1),wk(npm4+1),4*nord)
        nord=nord1
        write(*,176) nord
 176    format(' order after subset AIC = ',i2)
	if(iptrec.eq.1) write(3,176,err=704) nord
 704	continue
        if(nord.eq.0) go to 217
        write(*,214)
	if(iptrec.eq.1) write(3,214,err=705)
 705	continue
        do 177 i=1,nord
        if(wk(npm2+i).ne.0.) then
	write(*,215) i,wk(npm2+i)
	if(iptrec.eq.1) write(3,215,err=706) i,wk(npm2+i)
 706	continue
	endif
 177	continue
        go to 217
c
c   ioptm=3 or 4 :
c
 201    continue
        on=2./sqrt(float(n))
        write(*,214)
 214    format(' nonzero partials :')
	if(iptrec.eq.1) write(3,214,err=710)
 710	continue
        rvar=r0
	nord1=0
        do 213 i=1,nord
        if(abs(wk(npm2+i)).lt.on) wk(npm2+i)=0.
        if(wk(npm2+i).eq.0.) go to 213
        rvar=rvar*(1.-wk(npm2+i)*wk(npm2+i))
	nord1=i
        write(*,215) i,-wk(npm2+i)
 215    format(' part(',i2,') = ',f12.6)
	if(iptrec.eq.1) write(3,215,err=711) i,-wk(npm2+i)
 711	continue
 213    continue
	nord=nord1
 217    continue
        call ckaddi(args(7),nord,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(8),r0,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(9),rvar,iref)
        if(iref.eq.1) go to 99
        if(nord.eq.0) go to 219
        call partar(wk(npm2+1),nord,wk(1))
        lab='autoregressive coefficients'
        call ckadda(vname,nord,lab,1,iref)
 219    continue
        if(nargs.eq.9) go to 99
        if(iref.eq.1) go to 99
        lab='CAT'
        call ckadda(args(10),m,lab,npm3+1,iref)
c
c
  99    continue
        return
        end
c&SUBAIC
        subroutine subaic(part,nord,n,r0,iwk,nord1,parts,rvar)
c*****************************************************************
c
c   Subroutine to find the subset of the partial autocorrelations
c   that are significant according to the AIC criterion.
c
c   input :  nord,part(1),...,part(nord),n(sample size),r0
c
c   output:  nord1(highest lag partial), parts(1),...,parts(nord1):
c            new array of partials with 0's at nonsignificant lags.
c            rvar : residual variance of chosen model.
c
c*********************************************************************
c
        dimension part(nord),iwk(nord),parts(nord)
c
        do 10 i=1,nord
  10    parts(i)=abs(part(i))
        al1=0.0
        almin=0.0
        nsig=0
        nord1=0
        tovn=2./float(n)
        do 20 i=1,nord
        call max(parts,nord,pmax,indm)
        al1=al1+alog(1.-pmax*pmax)+tovn
        iwk(i)=indm
        parts(indm)=0.0
        if(al1.gt.almin)go to 20
        almin=al1
        nsig=i
        if(indm.gt.nord1) nord1=indm
 20     continue
        rvar=r0
        if(nsig.eq.0) go to 99
        do 30 i=1,nsig
        ii=iwk(i)
        rvar=rvar*(1.-part(ii)*part(ii))
  30    parts(ii)=part(ii)
  99    continue
c
c
        return
        end
