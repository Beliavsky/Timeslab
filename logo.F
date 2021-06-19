c&LOGO
c&LOGO
        subroutine logo
        write(*,10)
 10     format(/,1x,'TIMESLAB: A Time Series Analysis Laboratory')
        write(*,20)
 20     format(1x,'Unix Version 1.0, (C) Copyright 1997, H.J. Newton',/)
        return
        end

c&BINOM
        subroutine binom(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   x=BINOM(n,k)   or  x=binom(n,k,p)  or  x=binom(n,p)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
        integer*2 icki,ickse
c
c
        if(icki(args(1),1,n,1,1).eq.1) go to 99
c
c   handle x=binom(n,p) case:
c
        call ckreal(args(2),p)
        if(0.lt.p.and.p.lt.1.) then
                if(ickse(n+1).eq.1) go to 99
                alpq=alog(p)-alog(1.-p)
                wk(1)=float(n)*alog(1.-p)
                do 5 i=1,n
   5            wk(i+1)=wk(i)+alog(float(n-i+1))-alog(float(i))+alpq
                do 6 i=1,n+1
   6            wk(i)=exp(wk(i))
                lab='Binomial Probabilities'
                call ckadda(vname,n+1,lab,1,iref)
                go to 99
        endif
c
c
        if(icki(args(2),2,k,5,n).eq.1) go to 99
        if(k.lt.0) then
                call error(args,2,2)
                go to 99
        endif
        if(k.eq.0.or.k.eq.n) then
                x=0.
                go to 20
        endif
        kk=min(k,n-k)
        x=0.0
        do 10 i=1,kk
  10    x=x+alog(float(n-i+1))-alog(float(kk-i+1))
  20    continue
c
c
        if(nargs.eq.3) then
                call ckreal(args(3),p)
                if(p.lt.0..or.p.gt.1.) then
                        call error(args,3,2)
                        go to 99
                endif
                x=x+float(k)*alog(p)+float(n-k)*alog(1-p)
        endif
        x=exp(x)
c
c
        call ckaddr(vname,x,iref)
c
c
  99    continue
        return
        end
