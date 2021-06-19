c&REG
c&REG
        subroutine reg(args,nargs)
c******************************************************************
c
c   Subroutine to do regression analysis
c
c   REG(y,x,n,m,beta,vibeta,rss,ni,t,res,logrssio[,vy])
c
c   or REG(y,x,n,m,beta,rss,t,res)
c
c
c   or REG(y,x,n,m,pval,kopt,kp,beta,ind,rss)
c
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki,ickr,ickse
c
c
c
	if(ickl(args(1),npy,ny).eq.1) go to 99
	if(ickl(args(2),npx,nx).eq.1) go to 99
	if(icki(args(3),3,n,0,0).eq.1) go to 99
	if(icki(args(4),4,m,2,0).eq.1) go to 99
	iptn=0
	if(n.lt.0) then
		iptn=1
		n=-n
        endif
	iptv=0
	if(n.le.m) iptv=1
c
c   REG(y,x,n,m,pval,kopt,kp,beta,ind,rss)
c
	if(nargs.eq.10) then
	cws=float(n)+float(n+3)*float(m)+float(m+1)**2
	if(cws.gt.32767.) then
		if(ickse(maxnum+1).eq.1) go to 99
	endif
	nws=n+(n+3)*m+(m+1)**2
	if(ickse(nws).eq.1) go to 99
		if(n.gt.ny.or.n.le.m.or.m.lt.1.or.n*m.gt.nx) then
			write(*,701)
			imflag=1
			go to 99
		endif
		call ckreal(args(5),pval)
		if(pval.le.0..or.p.ge.1.) then
			call error(args,5,2)
			go to 99
		endif
		call ckint(args(6),kopt)
		if(kopt.lt.0.or.kopt.gt.m) then
			call error(args,6,2)
			go to 99
		endif
		if(kopt.gt.0) then
			if(ickl(args(9),npi,ni).eq.1) go to 99
			if(kopt.gt.ni) then
				call error(args,6,2)
				go to 99
			endif
			ns1=nstart(npi)-1
			do 300 i=1,kopt
	 		ii=array(ns1+i)
			if(ii.lt.1.or.ii.gt.m) then
				write(*,305)
 305				format(' Illegal Diagonal in REG')
				imflag=1
				go to 99
			endif
 300			iwk(i)=ii
		endif
		if(kopt.eq.0) then
			do 310 i=1,m
 310			iwk(i)=i
		endif
		mp1=m+1
		n2=mp1*mp1+1
		call reg4(array(nstart(npy)),array(nstart(npx)),n,m,
     1		mp1,wk,kopt,pval,kp,wk(n2),iwk,rss,iwk(m+1),ier,iptrec)
		if(ier.eq.1) then
			imflag=1
			go to 99
		endif
		call ckaddi(args(7),kp,iref)
		if(iref.eq.1.or.kp.eq.0) go to 99
		lab='Stepwise Regression Coefficients'
		call ckadda(args(8),kp,lab,n2,iref)
		if(iref.eq.1) go to 99
		do 315 i=1,kp
 315		wk(i)=iwk(m+i)
		lab='Regressors Chosen'
		call ckadda(args(9),kp,lab,1,iref)
		if(iref.eq.1) go to 99
		call ckaddr(args(10),rss,iref)
		go to 99
	endif
c
c   do simple form:
c
	if(nargs.eq.8) then
	cws=2.*float(n+m)+float(n)*float(m)+float(m+1)**2
	if(cws.gt.32767.) then
		if(ickse(maxnum+1).eq.1) go to 99
	endif
	nws=2*(n+m)+n*m+(m+1)**2
	if(ickse(nws).eq.1) go to 99
		if(n.gt.ny.or.n.le.m.or.m.lt.1.or.n*m.gt.nx) then
			write(*,701)
 701			format(' error in n or m')
			imflag=1
			go to 99
			endif
		mp1=m+1
		n1=1
		n2=n1+mp1*mp1
		n3=n2+m
		n4=n3+m
		call reg3(array(nstart(npy)),array(nstart(npx)),n,m,mp1,
     1		wk(n1),wk(n2),rss,wk(n3),wk(n4),ier)
		if(ier.eq.1) then
			imflag=1
			go to 99
		endif
		lab='Regression Coefficients'
		call ckadda(args(5),m,lab,n2,iref)
		if(iref.eq.1) go to 99
		lab='t-values'
		call ckadda(args(7),m,lab,n3,iref)
		if(iref.eq.1) go to 99
		lab='Regression Residuals'
		call ckadda(args(8),n,lab,n4,iref)
		if(iref.eq.1) go to 99
		call ckaddr(args(6),rss,iref)
		go to 99
	endif
c
c
c
	if(iptv.eq.0) then
		if(float(m+1)**2+float(m)+float(n).gt.32767.) then
			if(ickse(maxnum+1).eq.1) go to 99
		endif
		nws=(m+1)**2+n+m
		if(ickse(nws).eq.1) go to 99
	endif
	if(iptv.eq.1) then
		if(float(n+m)**2+float(n+m).gt.32767.) then
			if(ickse(maxnum+1).eq.1) go to 99
		endif
		nws=(n+m)**2+n+m
		if(ickse(nws).eq.1) go to 99
	endif
        ntm=n*m
	ipte=0
	if(n.gt.ny) ipte=3
	if(ntm.gt.nx) ipte=4
	if(ipte.ne.0) then
               call error(args,ipte,2)
               go to 99
        endif
	if(ickl(args(5),npb,nb).eq.1) go to 99
        if(nb.lt.m) then
               call error(args,5,2)
               go to 99
        endif
	if(ickl(args(6),npqi,nqi).eq.1) go to 99
        if(nqi.lt.m*m) then
               call error(args,6,2)
               go to 99
        endif
	if(ickr(args(7),7,rss,0,0.).eq.1) go to 99
	if(icki(args(8),8,ni,0,0).eq.1) go to 99
        if(iptv.eq.1) then
                if(nargs.eq.11) then
                        write(*,5)
   5                    format(' Need a VY in REG')
			imflag=1
                        go to 99
                endif
		if(ickl(args(12),npv,nv).eq.1) go to 99
                if(nv.lt.n*n) then
                        call error(args,12,2)
                        go to 99
                endif
        endif
c
c   Not Using V:
c
        if(iptv.eq.1) go to 200
        n1=(m+1)*(m+1)+1
        n2=n1+m
        call reg1(array(nstart(npy)),array(nstart(npx)),
     1  array(nstart(npqi)),n,m+1,n,m,wk,array(nstart(npb)),
     1  wk(n1),rss,wk(n2),ni,alrss,ier,iptn)
        if(ier.eq.1) then
                write(*,101)
 101            format(' singular matrix in REG')
		imflag=1
                go to 99
                endif
c
c   Put stuff in right place:
c
        lab='Regression t values'
        call ckadda(args(9),m,lab,n1,iref)
        if(iref.eq.1) go to 99
        lab='Regression Residuals'
        call ckadda(args(10),n,lab,n2,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(7),rss,iref)
        if(iref.eq.1) go to 99
        call ckaddi(args(8),ni,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(11),alrss,iref)
        go to 99
c
c   Using V:
c
 200    continue
	n1=(n+m)*(n+m)+1
	n2=n1+n
	call reg2(array(nstart(npy)),array(nstart(npx)),n+m,n,m,ni,
     1	array(nstart(npb)),array(nstart(npqi)),array(nstart(npv)),rss,
     1	wk,wk(n1),wk(n2),alrss,ier,iptn)
	if(ier.eq.1) then
		write(*,101)
		imflag=1
		go to 99
		endif
        lab='Regression t values'
        call ckadda(args(9),m,lab,n2,iref)
        if(iref.eq.1) go to 99
        lab='Regression Residuals'
        call ckadda(args(10),n,lab,n1,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(7),rss,iref)
        if(iref.eq.1) go to 99
        call ckaddi(args(8),ni,iref)
	if(iref.eq.1) go to 99
	call ckaddr(args(11),alrss,iref)
c
c
  99    continue
        return
        end
c&REG4
	subroutine reg4(y,x,n,m,mp1,a,kopt,pval,kp,beta,indxi,rss,
     1	indxo,ier,iptrec)
c*********************************************************************
c
c
c**********************************************************************
c
	dimension y(n),x(n,m),a(mp1,mp1),beta(m),indxi(m),indxo(m)
	double precision inprod
c
c
	mm=m
	if(kopt.gt.0) mm=kopt
	mmp1=mm+1
	do 10 i=1,mm
	ii=indxi(i)
		a(i,mmp1)=inprod(x(1,ii),y,n)
		a(mmp1,i)=a(i,mmp1)
		do 10 j=i,mm
		jj=indxi(j)
		a(i,j)=inprod(x(1,ii),x(1,jj),n)
  10	a(j,i)=a(i,j)
	a(mmp1,mmp1)=inprod(y,y,n)
	call selreg(a,mp1,mm,indxi,ktop,n,pval,kp,beta,indxo,ier,iptrec)
	if(ier.eq.1) then
		write(*,20)
  20		format(' Singular Matrix in REG')
		go to 99
	endif
c
c
  99	continue
	return
	end
c&REG3
	subroutine reg3(y,x,n,m,mp1,xtx,beta,rss,t,res,ier)
c******************************************************************
c
c
c
c******************************************************************
c
	dimension y(n),x(n,m),xtx(mp1,mp1),beta(m),t(m),res(n)
	double precision inprod
c
c
	do 10 i=1,m
	xtx(mp1,i)=inprod(x(1,i),y,n)
	xtx(i,mp1)=xtx(mp1,i)
	do 10 j=1,i
	xtx(i,j)=inprod(x(1,i),x(1,j),n)
  10	xtx(j,i)=xtx(i,j)
	xtx(mp1,mp1)=inprod(y,y,n)
	call swpk12(xtx,mp1,mp1,1,m,ier)
	do 13 i=1,mp1
  13	if(xtx(i,i).le.0.) ier=1
	if(ier.eq.1) then
		write(*,15)
  15		format(' Singularity in REG')
		go to 99
	endif
	do 20 i=1,m
  20	beta(i)=xtx(i,mp1)
	rss=xtx(mp1,mp1)
	sde=sqrt(rss/float(n-m))
	do 30 i=1,m
  30	t(i)=beta(i)/(sde*sqrt(xtx(i,i)))
	do 40 i=1,n
	c=y(i)
		do 35 j=1,m
  35		c=c-x(i,j)*beta(j)
  40	res(i)=c
c
c
  99	continue
	return
	end
c&REG2
        subroutine reg2(y,x,mdim,n,m,ni,beta,vibeta,vy,rss,w,res,t,
     1  alrss,ier,iptn)
c******************************************************************
c
c
c******************************************************************
c
	dimension y(n),beta(m),vibeta(m,m),vy(n,n),w(mdim,mdim),res(n)
	dimension x(n,m),t(m)
c
c   Form Matrix to be Swept:
c
	iswn=-1
	if(iptn.eq.1) iswn=1
	ni=ni-iswn*n
        do 10 i=1,m
        do 10 j=1,m
  10    w(i,j)=vibeta(i,j)
        do 20 i=1,n
        do 20 j=1,n
  20    w(m+i,m+j)=vy(i,j)
        do 30 i=1,n
        do 30 j=1,m
        w(j,m+i)=x(i,j)
  30    w(m+i,j)=w(j,m+i)*iswn
c
c   Sweep on beta:
c
        call swpk12(w,mdim,m+n,1,m,ier)
        if(ier.eq.1) go to 99
c
c   Save Output vy:
c
        do 40 i=1,n
        do 40 j=1,n
  40    vy(i,j)=w(m+i,m+j)
c
c   Sweep on x:
c
        call swpk12(w,mdim,m+n,m+1,m+n,ier)
        if(ier.eq.1) go to 99
c
c   Find first residuals:
c
        do 50 i=1,n
                c=y(i)
                do 45 j=1,m
  45            c=c-x(i,j)*beta(j)
  50    res(i)=c
c
c   Find output beta:
c
        do 70 i=1,m
                c=beta(i)
                do 60 j=1,n
  60            c=c+iswn*w(i,m+j)*res(j)
  70    beta(i)=c
c
c   find rss out:
c
	alrss=alog(rss)
        do 80 i=1,n
        do 80 j=1,n
  80    rss=rss-iswn*res(i)*res(j)*w(m+i,m+j)
	alrss=alog(rss)-alrss
c
c   find t and residuals:
c
        ni1=ni-m
        s2rs=rss/float(ni1)
        do 85 i=1,m
  85    t(i)=beta(i)/sqrt(s2rs*w(i,i))
        do 87 i=1,n
                c=y(i)
                do 86 j=1,m
  86            c=c-beta(j)*x(i,j)
  87    res(i)=c

c
c   find output vibeta:
c
        call swpk12(w,mdim,m+n,1,m,ier)
        if(ier.eq.1) go to 99
        do 90 i=1,m
        do 90 j=1,m
  90    vibeta(i,j)=w(i,j)
c
c
  99    continue
        return
        end
c&REG1
        subroutine reg1(y,x,qi,ndim,mdim,n,m,xtxa,beta,t,rss,res,ni,
     1	alrss,ier,iptn)
c*******************************************************************
c
c   Subroutine to regress y on x.
c
c*******************************************************************
c
        dimension y(n),x(ndim,m),xtxa(mdim,mdim),beta(m),t(m),res(n)
        dimension qi(m,m)
        double precision inprod
c
c   Form augmented matrix :
c
	if(ni.eq.0) then
		rss=0.
		do 3 i=1,m
		do 3 j=1,m
   3		xtxa(i,j)=0.0
	endif
c
c
	if(ni.ne.0) then
		do 4 i=1,m
		do 4 j=1,m
   4		xtxa(i,j)=qi(i,j)
	endif
c
c
        do 7 i=1,n
                c=y(i)
                do 6 j=1,m
   6            c=c-x(i,j)*beta(j)
   7    res(i)=c
c
c
	iswn=1
	if(iptn.eq.1) iswn=-1
	oswn=iswn
        mp1=m+1
        do 10 i=1,m
        xtxa(i,mp1)=inprod(x(1,i),res,n)
        xtxa(mp1,i)=oswn*xtxa(i,mp1)
                do 10 j=i,m
                xtxa(i,j)=xtxa(i,j)+iswn*inprod(x(1,i),x(1,j),n)
  10    xtxa(j,i)=xtxa(i,j)
        xtxa(mp1,mp1)=inprod(res,res,n)
c
c
	rssin=rss+iswn*xtxa(mp1,mp1)
        do 15 i=1,m
        do 15 j=1,m
  15    qi(i,j)=xtxa(i,j)
c
c   Sweep:
c
        call swpk12(xtxa,mdim,mp1,1,m,ier)
        if(ier.eq.1) go to 99
c
c   calculate other stuff:
c
        do 20 i=1,m
  20    beta(i)=beta(i)+iswn*xtxa(i,mp1)
        do 50 i=1,n
                c=y(i)
                do 40 j=1,m
  40            c=c-beta(j)*x(i,j)
  50    res(i)=c
	rss=rss+iswn*xtxa(mp1,mp1)
	alrss=alog(rssin)-alog(rss)
	ni=ni+iswn*n
        s2rs=rss/float(ni-m)
        do 30 i=1,m
  30    t(i)=xtxa(i,mp1)/sqrt(s2rs*xtxa(i,i))
c
c
  99    continue
        return
        end
c&SWPK12
        subroutine swpk12(a,mdim,m,k1,k2,ier)
c*****************************************************************
c
c   Subroutine to sweep the mxm matrix a on its k1 st thru k2 th
c   diagonals.  ier is 1 if a is singular.
c
c******************************************************************
c
        dimension a(mdim,mdim)
	integer*2 iwk(100)
c
        ier=0
	do 1 k=k1,k2
   1	iwk(k)=k
        do 50 k=k1,k2
c
c   find first nonzero diagonal element and swap index with kth:
c
		do 70 kk=k,k2
		k3=iwk(kk)
  70		if(abs(a(k3,k3)).gt.1.e-25) go to 71
		go to 99
  71		kc=iwk(k)
		iwk(k)=iwk(kk)
		iwk(kk)=kc
	kk=iwk(k)
        d=1./a(kk,kk)
        a(kk,kk)=1.
                do 10 i=1,m
  10            a(kk,i)=d*a(kk,i)
                do 20 j=1,m
                if(j.eq.kk) go to 20
                a(j,kk)=-a(j,kk)*d
  20            continue
                do 40 j=1,m
                if(j.eq.kk) go to 40
                        do 30 i=1,m
                        if(i.eq.kk) go to 30
                        a(j,i)=a(j,i)+a(j,kk)*a(kk,i)/d
  30                    continue
  40            continue
  50    continue
        return
  99    ier=1
        return
        end
