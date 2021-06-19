c&SEAS
c&SEASES
        subroutine seases(args,nargs)
c********************************************************************
c
c   Subroutine to handle the command
c
c   seas(y,n,ords,coeffs,lags,it,eps,nback,rvar,ier,sds[,e])
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        integer*2 ickl,ickr,icki,ickse
c
c
        if(ickl(args(1),npy,ny).eq.1) go to 99
        call ckint(args(2),n)
        if(n.lt.1.or.n.gt.ny) then
                call error(args,2,2)
                go to 99
        endif
        if(ickl(args(3),npords,nords).eq.1) go to 99
        if(nords.lt.5) then
                call error(args,3,2)
                go to 99
        endif
        ns1=nstart(npords)-1
        nn=0
        do 10 i=1,5
        ii=array(ns1+i)
        nn=nn+ii
  10    iwk(i)=ii
c
c
        maxp=iwk(1)
        maxq=iwk(3)
        if(nn-iwk(5).eq.0) then
                write(*,15)
  15            format(' One of the orders must be positive')
                imflag=1
                go to 99
        endif
c
c   get coefficients:
c
        if(ickl(args(4),npcofs,ncofs).eq.1) go to 99
        if(ncofs.lt.nn) then
                call error(args,4,2)
                go to 99
        endif
c
c   get lags (if necessary):
c
        nlags1=iwk(2)+iwk(4)
        if(nlags1.gt.0) then
                if(ickl(args(5),nplags,nlags).eq.1) go to 99
                if(nlags.lt.nlags1) then
                        call error(args,5,2)
                        go to 99
                endif
                ns1=nstart(nplags)-1
                do 20 i=1,nlags1
  20            iwk(5+i)=array(ns1+i)
                npl=iwk(2)
                nql=iwk(4)
                if(npl.gt.0) maxp=maxp+iwk(5+npl)
                if(nql.gt.0) maxq=maxq+iwk(5+npl+nql)
        endif
c
c   get maxit:
c
        if(icki(args(6),6,maxit,3,0).eq.1) go to 99
        iptit=0
        if(maxit.lt.0) then
                iptit=1
                maxit=-maxit
                if(iptrec.eq.1) iptit=2
        endif
c
c   get eps and nt:
c
        if(ickr(args(7),7,eps,2,0.).eq.1) go to 99
        if(icki(args(8),8,nt,1,0).eq.1) go to 99
c
c   check workspace size:
c
        cws1=float(n+nt)
        cws2=float(nn+1)
        cws=3.*cws1+5.*float(nn)+3.*float(maxp)+float(maxq)+
     1  float(nn)*cws1+2.*cws2**2
        if(cws.gt.float(maxnum)) then
                if(ickse(maxnum+1).eq.1) go to 99
        endif
c
c   nn is the number of coefficients in the model
c
        ndim=nn+1
        n1=n+1
        n2=n1+nn
        n3=n2+maxp+1
        n4=n3+maxq+1
        n5=n4+n+nt
        n6=n5+n+nt
        n7=n6+maxp+maxq
        n8=n7+nt+1
        n9=n8+nn*(n+nt)
        n10=n9+(nn+1)*(nn+1)
        n11=n10+(nn+1)*(nn+1)
        n12=n11+nn
        n13=n12+nn
        n14=n13+nn
        m1=6
        m2=m1+nlags1+1
        m3=m2+nlags+1
        m4=m3+nlags1+1
        call movxy(wk(n1),array(nstart(npcofs)),4*nn)
        call movxy(wk,array(nstart(npy)),4*n)
c
c
        npnt=n+nt
        call marq(wk,n,ndim,maxit,eps,iwk,wk(n1),nt,npnt,nn,
     1  wk(n2),wk(n3),wk(n4),wk(n5),wk(n6),wk(n7),wk(n8),wk(n9),
     1  wk(n10),wk(n11),wk(n12),wk(n13),wk(n14),iwk(m1),iwk(m2),
     1  iwk(m3),iwk(m4),ier,ss,iptit)
c
c
        rvar=ss/float(n-nn)
        call ckaddi(args(10),ier,iref)
        if(iref.eq.1) go to 99
        call ckaddr(args(9),rvar,iref)
        if(iref.eq.1) go to 99
        lab='Coefficients From SEASEST'
        call ckadda(args(4),nn,lab,n1,iref)
        if(iref.eq.1) go to 99
        do 813 i=1,nn
 813    wk(n13+i-1)=sqrt(wk(n13+i-1)*rvar)
        lab='Standard Errors From SEASEST'
        call ckadda(args(11),nn,lab,n13,iref)
        if(iref.eq.1) go to 99
        if(nargs.eq.12) then
        lab='Errors From SEASEST'
        call ckadda(args(12),n,lab,n5-n+1,iref)
        if(iref.eq.1) go to 99
        endif
c
c
  99    continue
        return
        end
c&MARQ
        subroutine marq(y,n,ndim,maxit,eps,nords,coeffs,nt,npnt,ntot,
     1  alpha,beta,e,e1,a1,xx,x,a,as,g,h,d,beta1,lags,lp,lq,nl,ier,ss,
     1  iptit)
c******************************************************************
c
c******************************************************************
c
c
        double precision inprod
        dimension y(1),nords(1),coeffs(1),alpha(1),beta(1),e(1),e1(1),
     1  a1(1),xx(1),x(npnt,1),a(ndim,ndim),as(ndim,ndim),g(1),h(1),
     1  d(1),beta1(1),lags(1),lp(1),lq(1),nl(1)
c
c
c
c
        ppi=.1
        f2=2.
        delta=.05
        npnt=n+nt
        ntp=ntot+1
c
c
        do 100 it=1,maxit
c
c
        if(iptit.eq.1.or.iptit.eq.2) then
                write(*,300) it,(coeffs(i),i=1,min0(ntot,6))
                if(ntot.gt.6) write(*,301) (coeffs(i),i=7,ntot)
 301            format(6x,6f10.4)
                if(iptit.eq.2) then
                write(3,300,err=1900) it,(coeffs(i),i=1,min0(ntot,6))
             if(ntot.gt.6) write(3,301,err=1900) (coeffs(i),i=7,ntot)
 300            format(1x,i5,6f10.4)
1900            continue
                endif
        endif
c
        ss=seaslk(y,n,nords,coeffs,lags,lp,lq,alpha,beta,
     1  e,a1,nl,xx,nt)
c
        do 10 i=1,ntot
        coeffs(i)=coeffs(i)+delta
        ss1=seaslk(y,n,nords,coeffs,lags,lp,lq,alpha,beta,
     1             e1,a1,nl,xx,nt)
        coeffs(i)=coeffs(i)-delta
                do 20 j=1,npnt
  20            x(j,i)=(e(j)-e1(j))/delta
                do 30 j=1,i
                a(i,j)=inprod(x(1,i),x(1,j),npnt)
c                a(i,j)=inprod(x(nt+1,i),x(nt+1,j),n)
  30            a(j,i)=a(i,j)
c
c
                if(maxit.eq.1) then
                       ier=1
                       go to 200
                endif
c
c
                g(i)=inprod(x(1,i),e,npnt)
c                g(i)=inprod(x(nt+1,i),e(nt+1),n)
                d(i)=sqrt(a(i,i))
  10            continue
c
c
  25            do 40 i=1,ntot
                as(ntp,i)=g(i)/d(i)
                as(i,ntp)=as(ntp,i)
                        do 40 j=1,ntot
                        as(i,j)=a(i,j)/(d(i)*d(j))
  40            continue
                do 50 i=1,ntot
  50            as(i,i)=as(i,i)+ppi
                as(ntp,ntp)=1.
                call swpk12(as,ntp,ntp,1,ntot,ier)
c
c
                if(ier.ne.0) then
                        ier=2
                        go to 99
                endif
c
c
                do 60 i=1,ntot
                h(i)=as(i,ntp)/d(i)
  60            beta1(i)=coeffs(i)+h(i)
                ss1=seaslk(y,n,nords,beta1,lags,lp,lq,alpha,beta,
     1          e1,a1,nl,xx,nt)
c
c
                if(ss1.ge.ss) then
                        if(ppi.gt.200.) then
                                ier=3
                                go to 99
                        endif
                        ppi=ppi*f2
                        go to 25
                endif
c
c
                do 70 i=1,ntot
                cc=abs(h(i))
  70            if(cc.gt.eps) go to 75
                call movxy(coeffs,beta1,4*ntot)
                ier=0
                go to 200
  75            call movxy(coeffs,beta1,4*ntot)
                if(ppi.lt.1.e-20) then
                        ier=4
                        go to 99
                endif
                ppi=ppi/f2
c
c
 100    continue
        ier=1
 200    continue
  99    continue
        call swpk12(a,ndim,ntot,1,ntot,ier)
        do 110 i=1,ntot
 110    d(i)=a(i,i)
        if(nords(5).eq.1) then
                cc=coeffs(ntot)
                gg1=1.
                if(nords(1).gt.0) then
                        do 120 i=1,nords(1)
 120                    gg1=gg1+coeffs(i)
                endif
                gg2=1.
                if(nords(2).gt.0) then
                        do 130 i=1,nords(2)
 130                    gg2=gg2+coeffs(nords(1)+i)
                endif
                gg=gg1*gg2
                cc=cc*gg
                coeffs(ntot)=cc
                d(ntot)=d(ntot)*gg*gg
        endif
        return
        end
c&SEASLK
        real*4 function seaslk(x,n,nords,coeffs,lags,lp,lq,alpha,beta,
     1  e,a,nl,xx,nt)
c*********************************************************************
c
c   Function to evaluate a seasonal ARMA sum of squares.
c
c**********************************************************************
c
        dimension x(1),nords(1),coeffs(1),lags(1),lp(1),lq(1),alpha(1),
     1  beta(1),e(1),a(1),nl(1),xx(1)
c
c   Convert to subset ARMA representation:
c
        call convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,alpha,
     1  beta,a,nl,amu)
c
c   find eta:
c
        if(nt.eq.0) go to 715
        if(maxq.gt.0) then
                do 310 i=n-maxp+1,n-maxp+maxq
 310            e(i)=0.0
        endif
        do 360 i=n-maxp,1,-1
                mi=maxp+i
                mj=maxq+i
                c=x(i)-amu
                if(nlp.eq.0) go to 330
                do 320 j=1,nlp
                        jj=lp(j)
 320            c=c+alpha(j)*(x(i+jj)-amu)
 330            continue
                if(nlq.eq.0) go to 350
                do 340 k=1,nlq
                        kk=lq(k)
 340            c=c-beta(k)*e(i+kk)
 350            continue
 360    e(i)=c
c
c   now get x(0),x(-1),...
c
        do 400 i=0,-nt+1,-1
                c=0.0
                if(nlq.eq.0) go to 430
                do 420 k=1,nlq
                ii=i+lq(k)
                if(ii.lt.1) go to 420
                c=c+beta(k)*e(ii)
 420            continue
 430            if(nlp.eq.0) go to 450
                do 440 j=1,nlp
                ii=i+lp(j)
                if(ii.gt.0) c1=x(ii)-amu
                if(ii.le.0) c1=xx(nt+ii)-amu
 440            c=c-alpha(j)*c1
 450            continue
 400    xx(nt+i)=c+amu
c
c   now get e(-nt+1),...,e(n):
c
c
 715    ss=0.0
        ij=1-nt
        do 60 i=-nt+1,n
        if(i.gt.0) c=x(i)-amu
        if(i.le.0) c=xx(nt+i)-amu
                if(nlp.eq.0) go to 30
                do 20 j=1,nlp
                ii=i-lp(j)
                c1=0.0
                if(ii.gt.0) c1=x(ii)-amu
                if(ii.le.0.and.ii.ge.ij) c1=xx(nt+ii)-amu
  20            c=c+alpha(j)*c1
  30            continue
                if(nlq.eq.0) go to 50
                do 40 k=1,nlq
                ii=i-lq(k)
                if(ii.lt.ij) go to 40
                c=c-beta(k)*e(nt+ii)
  40            continue
  50            continue
        e(nt+i)=c
  60    ss=ss+c*c
c
c
        seaslk=ss
        return
        end
c&CONVT
        subroutine convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,
     1  alpha,beta,a,nl,amu)
c*******************************************************************
c
c   Subroutine to convert the input to SEASEST to that of a subset
c   ARMA model.
c
c   Input:
c      nords : p, P, q, Q, M 
c      coeffs: vector of length nords(1)+...+nords(5) with start values
c      lags  : vector of length nords(2)+nords(4) with subset lags
c     
c   Output:  
c      maxp, maxq : maximum AR and MA lag in subset version of model
c      nlp, nlq   : number of AR and MA lags in subset version of model
c      lp, lq     : arrays of length nlp and nlq with AR and MA lags
c                   in subset version of model
c      alpha, beta: arrays of length nlp and nlq with AR and MA coeffs
c                   in subset version of model
c      amu        : 0 (if M=1) or last element of coeffs (if M=1)
c
c   Work:
c      a, nl : arrays of length max(maxp,maxq)
c
c   Subprograms called: None
c
c*******************************************************************
c
        dimension nords(1),coeffs(1),lags(1),lp(1),lq(1),
     1  alpha(1),beta(1),a(1),nl(1)
c
c
        np=nords(1)
        npl=nords(2)
        nq=nords(3)
        nql=nords(4)
        amu=0.0
        if(nords(5).eq.1) amu=coeffs(np+npl+nq+nql+1)
        n1=0
        n2=n1+np
        n3=n2+npl
        n4=n3+nq
c
c
        maxp=np
        maxq=nq
        if(npl.gt.0) maxp=maxp+lags(npl)
        if(nql.gt.0) maxq=maxq+lags(npl+nql)
c
c
        nlp=0
        nlq=0
c
c
        if(maxp.gt.0) then
c
c
                do 10 i=1,maxp
                nl(i)=0
  10            a(i)=0.0
c
c
                if(np.gt.0) then
                        do 15 i=1,np
                        nl(i)=1
  15                    a(i)=coeffs(i)
                endif
c
c
                if(npl.gt.0) then
c
                        do 20 i=0,np
                                c1=1.
                                if(i.gt.0) c1=coeffs(i)
                        do 20 j=1,npl
                                jj=i+lags(j)
                                nl(jj)=1
  20                    a(jj)=a(jj)+c1*coeffs(n2+j)
                endif
c
c
                do 40 i=1,maxp
                        if(nl(i).gt.0) then
                                nlp=nlp+1
                                lp(nlp)=i
                                alpha(nlp)=a(i)
                        endif
  40            continue
c
c
        endif
c
c
        if(maxq.gt.0) then
c
                do 60 i=1,maxq
                nl(i)=0
  60            a(i)=0.0
c
                if(nq.gt.0) then
                        do 65 i=1,nq
                        nl(i)=1
  65                    a(i)=coeffs(n3+i)
                endif
c
                if(nql.gt.0) then
c
                        do 70 i=0,nq
                                c1=1.
                                if(i.gt.0) c1=coeffs(n3+i)
                        do 70 j=1,nql
                                jj=i+lags(npl+j)
                                nl(jj)=1
  70                    a(jj)=a(jj)+c1*coeffs(n4+j)
c
                endif
c
c
                do 90 i=1,maxq
                        if(nl(i).gt.0) then
                                nlq=nlq+1
                                lq(nlq)=i
                                beta(nlq)=a(i)
                        endif
  90            continue
c
c
        endif
c
c
        return
        end
c&SEASPR
        subroutine seaspr(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   SEASPRED(x,n,ords,coeffs,lags,rvar,tf,tl,hl,conf,xp,xpl,xpu,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        integer*2 ickl,icki,ickr,ickse
c
c   check that x exists (npx will be which array and nx is its length):
c
        if(ickl(args(1),npx,nx).eq.1) go to 99
c
c   check that n is defined and is at least 1 and not > nx:
c
        if(icki(args(2),2,n,1,1).eq.1) go to 99
        if(n.gt.nx) then
                call error(args,2,2)
                go to 99
        endif
c
c   move x[1],...,x[n] into wk(1),...,wk(n):
c
        call movxy(wk,array(nstart(npx)),4*n)
c
c   check that nords exists and has at least 8 elements:
c
        if(ickl(args(3),npords,nords).eq.1) go to 99
        if(nords.lt.8) then
                call error(args,3,2)
                go to 99
        endif
c
c   put ords into iwk(1),...,iwk(8) and get nn = # of coefficients:
c
        ns1=nstart(npords)-1
        nn=0
        do 10 i=1,8
        iwk(i)=array(ns1+i)
  10    if(i.le.5) nn=nn+iwk(i)
c
c   nd is # of 1st differences, ndd is # of Sth differences, and 
c   ns is S:
c
        nd=iwk(6)
        ndd=iwk(7)
        ns=iwk(8)
        if(nd.lt.0.or.ndd.lt.0.or.ns.lt.0) then
                call error(args,3,2)
                go to 99
        endif
c
c   check that coeffs exists and has at least nn+2 elements:
c
        if(ickl(args(4),npcofs,ncofs).eq.1) go to 99
        if(ncofs.lt.nn+2) then
                call error(args,4,2)
                go to 99
        endif
c
c   move parameters, m, and lambda into wk(n+1),...:
c
        call movxy(wk(n+1),array(nstart(npcofs)),4*(nn+2))
c
c   get maxp = p + d + D*S and maxq = q:
c
        maxp=iwk(1)
        maxp=maxp+nd+ndd*ns
        maxq=iwk(3)
c
c   read AR and MA lags and update maxp and maxq:
c
        nlags1=iwk(2)+iwk(4)
        if(nlags1.gt.0) then
                if(ickl(args(5),nplags,nlags).eq.1) go to 99
                if(nlags.lt.nlags1) then
                        call error(args,5,2)
                        go to 99
                endif
                ns1=nstart(nplags)-1
                do 20 i=1,nlags1
  20            iwk(8+i)=array(ns1+i)
                npl=iwk(2)
                nql=iwk(4)
                if(npl.gt.0) maxp=maxp+iwk(8+npl)
                if(nql.gt.0) maxq=maxq+iwk(8+npl+nql)
        endif
c
c   now npl = # of AR lags, nql = # of MA lags, and the lags are in
c   iwk(9),...
c
c   get rvar and make sure it's positive:
c
        if(ickr(args(6),6,rvar,2,0.).eq.1) go to 99
c
c   get ntf and make sure it is at least maxp+maxq:
c
        if(icki(args(7),7,ntf,1,1).eq.1) go to 99
        if(ntf.lt.maxp+maxq) then
                call error(args,7,2)
                go to 99
        endif
c
c   get ntl and make sure it is not > ntl or > n:
c
        if(icki(args(8),8,ntl,1,ntf).eq.1) go to 99
        if(ntl.gt.n) then
                call error(args,8,2)
                go to 99
        endif
c
c   get nhl and make sure it is at least 1:
c
        if(icki(args(9),9,nhl,1,1).eq.1) go to 99
c
c   get conf and make sure it is between 0 and 1:
c
        if(ickr(args(10),10,conf,2,0.0).eq.1) go to 99
        if(conf.ge.1.) then
                call error(args,10,2)
                go to 99
        endif
c
c   check that there is enough workspace:
c
        ncf=maxp+maxq+1
        cws=5.*float(ncf)+3.*float(n)+float(nhl)+
     1  3.*float(ntl-ntf+1)*float(nhl)
        if(cws.gt.float(maxnum)) then
                if(ickse(maxnum+1).eq.1) go to 99
        endif
c
c   partition the work array to become the arrays in sspr:
c
        ncf1=(ntl-ntf+1)*nhl
        n1=n+1
        n2=n1+nn+3
        n3=n2+ncf
        n4=n3+ncf
        n5=n4+ncf
        n6=n5+ncf
        n7=n6+n
        n8=n7+n+nhl
        n9=n8+ncf1
        n10=n9+ncf1
        m1=9
        m2=m1+ncf
        m3=m2+ncf
        m4=m3+ncf
c
c
        call sspr(wk,n,iwk,wk(n1),iwk(m1),conf,ntf,ntl,nhl,iwk(m2),
     1  iwk(m3),wk(n2),wk(n3),wk(n4),wk(n5),iwk(m4),wk(n6),wk(n7),
     1  wk(n8),wk(n9),wk(n10),ier,npds,rvar)
c
c
        call ckaddi(args(14),ier,iref)
        if(iref.eq.1.or.ier.eq.1) go to 99
        lab='Predictors from SEASPRED'
        call ckadda(args(11),ncf1,lab,n8,iref)
        if(iref.eq.1) go to 99
        lab='Lower Probability Limits'
        call ckadda(args(12),ncf1,lab,n9,iref)
        if(iref.eq.1) go to 99
        lab='Upper Probability Limits'
        call ckadda(args(13),ncf1,lab,n10,iref)
c
c
  99    continue
        return
        end
c&SSPR
        subroutine sspr(x,n,nords,coeffs,lags,conf,ntf,ntl,nhl,
     1  lp,lq,alpha,beta,a,a1,nl,e,xx,xp,xpl,xpu,ier,npds,rvar)
c*******************************************************************
c
c
c*******************************************************************
c
        real*4 fctr10
        double precision dp,dq
        dimension x(1),nords(1),coeffs(1),lags(1),lp(1),lq(1),
     1  alpha(1),beta(1),
     1  a(1),a1(1),nl(1),e(1),xx(1),xp(1),xpl(1),xpu(1)
c
c
        ier=0
c
c   express model (without differencing) as subset ARMA:
c
        call convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,
     1  alpha,beta,a,nl,amu)
c
c   adjust for 1st differences:
c
        nd=nords(6)
        ndd=nords(7)
        ns=nords(8)
        maxq2=maxq
        maxp2=maxp+nd+ndd*ns
c
        if(maxp2.eq.0) go to 26
        do 10 i=1,maxp2
                a(i)=0.0
  10    nl(i)=0
        do 20 i=0,nlp
                c1=1.0
                ii=0
                if(i.gt.0) then
                        ii=lp(i)
                        c1=alpha(i)
                endif
                do 20 j=0,nd
                        c2=((-1.)**j)*fctr10(nd,j)
                        do 20 k=0,ndd
                                c3=((-1.)**k)*fctr10(ndd,k)
                                ijks=ii+j+k*ns
                                if(ijks.gt.0) then
                                        nl(ijks)=1
                                        a(ijks)=a(ijks)+c1*c2*c3
                                endif
  20    continue
        nlp=0
        do 25 i=1,maxp2
                if(nl(i).eq.1) then
                        nlp=nlp+1
                        lp(nlp)=i
                        alpha(nlp)=a(i)
                endif
  25    continue
  26    continue
c
c   now nlp,nlq,alpha,beta,lp,lq are lags and coefficients of
c   full model
c
c   transform the series:
c
        ntot=nords(1)+nords(2)+nords(3)+nords(4)+nords(5)
        am=coeffs(ntot+1)
        alam=coeffs(ntot+2)
c
c   add the constant to make data positive:
c
        do 200 i=1,n
 200    x(i)=x(i)+am
c
c   power transform:
c
        iptlam=1
        if(abs(alam-1.).gt..00001) then
c 
c   check data + constant all positive:
c
                call min(x,n,xmin,imin)
                if(xmin.le.0.) then
                        write(*,5)
   5            format(' Data + constant not all positive in SEASPRED')
                        ier=1
                        go to 99
                endif
c
c   do log transform:
c
                if(abs(alam).lt..0001) then
                        iptlam=2
                        do 210 i=1,n
 210                    x(i)=alog(x(i))
                endif
c
c   other power transforms:
c
                if(abs(alam).ge..0001) then
                        iptlam=3
                        do 220 i=1,n
 220                    x(i)=x(i)**alam
                endif
        endif  
c
c   forecasts: first one possible is x(maxp2+maxq2+1):
c
c
c   find one step ahead forecast errors for x(maxp2+1),...,x(n):
c
        do 250 i=1,n
 250    e(i)=0.0
        do 260 i=maxp2+1,n
                c=x(i)-amu
                if(nlp.gt.0) then
                        do 270 j=1,nlp
                        if(i-lp(j).lt.1) go to 270
                        c=c+alpha(j)*x(i-lp(j))
 270                    continue
                endif
                if(nlq.gt.0) then
                        do 280 j=1,nlq
                        if(i-lq(j).lt.1) go to 280
                        c=c-beta(j)*e(i-lq(j))
 280                    continue
                endif
 260    e(i)=c
c
c   find forecasts:
c
        npds=1
        do 350 nt=ntf,ntl
        call movxy(xx,x,4*nt)
                do 360 nh=1,nhl
                ntph=nt+nh
                c=amu
                if(maxq2.gt.0) then
                        do 370 j=1,nlq
                        jj=ntph-lq(j)
                        if((jj.lt.1).or.(jj.gt.nt)) go to 370
                        c=c+beta(j)*e(ntph-lq(j))
 370                    continue
                endif
                if(maxp2.gt.0) then
                        do 380 j=1,nlp
 380                    c=c-alpha(j)*xx(ntph-lp(j))
                endif
 360            xx(ntph)=c
        call movxy(xp(npds),xx(nt+1),4*nhl)
        npds=npds+nhl
 350    continue
        npds=npds-1
c
c   find psi weights:
c
        if(maxp2.gt.0) then
                do 410 i=1,maxp2
 410            a(i)=0.0
                do 420 i=1,nlp
 420            a(lp(i))=alpha(i)
        endif
        if(maxq2.gt.0) then
                do 430 i=1,maxq2
 430            a1(i)=0.0
                do 440 i=1,nlq
 440            a1(lq(i))=beta(i)
        endif
        do 450 i=1,nhl
 450    xx(i)=0.0
        if(maxq2.gt.0) then
                do 451 i=1,maxq2
 451            xx(i)=a1(i)
        endif
        if(maxp2.gt.0) then
                do 453 i=1,nhl
                        c=xx(i)
                        do 452 j=1,min0(i,maxp2)
                                c1=1.0
                                if(i-j.gt.0) c1=xx(i-j)
 452                    c=c-a(j)*c1
 453            xx(i)=c
        endif
c
        e(1)=rvar
        if(nhl.gt.1) then
                do 470 i=2,nhl
                c=1.
                        do 480 j=1,i-1
 480                    c=c+xx(j)*xx(j)
 470            e(i)=rvar*c
        endif
c
c
        conf=(1.+conf)/2.
        dp=dble(conf)
        call norqnt(dp,dq,ifault)
        q=dq
        do 490 i=1,nhl
 490    e(i)=q*sqrt(e(i))
        do 500 nt=ntf,ntl
        ntt=(nt-ntf)*nhl
        do 500 nh=1,nhl
        nthh=ntt+nh
        xpl(nthh)=xp(nthh)-e(nh)
 500    xpu(nthh)=xp(nthh)+e(nh)
c
c   transform back:
c
        if(iptlam.eq.3) then
                call min(xpl,npds,xmin,imin)
                if(xmin.le.0.) then
                        write(*,505)
 505            format(' Can''t transform back because negative limit')
                        ier=1
                        go to 99
                endif
                alam1=1./alam
                do 510 i=1,npds
                xp(i)=xp(i)**alam1
                xpl(i)=xpl(i)**alam1
 510            xpu(i)=xpu(i)**alam1
        endif
        if(iptlam.eq.2) then
                do 520 i=1,npds
                xp(i)=exp(xp(i))
                xpl(i)=exp(xpl(i))
 520            xpu(i)=exp(xpu(i))
        endif
 600    do 610 i=1,npds
        xp(i)=xp(i)-am
        xpl(i)=xpl(i)-am
 610    xpu(i)=xpu(i)-am
c
c
  99    continue
        return
        end
c&FCTR10
        real*4 function fctr10(n,k)
c*************************************************************************
c
c   function to find the binomial coefficient { n choose k }.
c
c************************************************************************
c
        double precision c,c1,c2
c
        if(k.eq.0.or.k.eq.n) then
                fctr10=1.
                return
        endif
        kk=min0(k,n-k)
        c=1.d0
        do 10 i=1,kk
        c1=dble(float(n-i+1))
        c2=dble(float(i))
  10    c=c*c1/c2
        fctr10=c
        return
        end
