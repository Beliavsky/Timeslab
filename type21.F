c&TYPE2
c&TYPE2
        subroutine type2(l,args,vname)
c*******************************************************************
c
c   Subroutine to handle type 2 commands.
c
c*******************************************************************
c
#include 'tslabc'
        character l(72)*1,vname*15,args(1)*15
	character ltmp(72)*1
	integer*2 icktyp,ickear,ickop,ickpop
c
c   ntmp is counter for temporary variables:
c
	ntmp=0
c
c   look for parentheses:
c
   1	continue
   	nparen=0
	do 61 nparr=1,lnb
  61	if(l(nparr).eq.'}') go to 62
	go to 5
  62 	do 63 nparl=nparr,1,-1
  63	if(l(nparl).eq.'{') go to 64
	go to 98
c
c   found a pair in characters nparl and nparr:
c
  64	nparen=1
	ntmp=ntmp+1
	lnbtmp=lnb
	nf=7
	if(ntmp.gt.9) nf=8
	call movxy(ltmp,l,lnb)
	lnb=nf-1+nparr-nparl-1
	call movxy(l(nf),ltmp(nparl+1),nparr-nparl-1)
	l(1)='t'
	l(2)='m'
	l(3)='p'
	l(4)='&'
	n1=ntmp/10
	n2=ntmp-10*n1
	if(n1.eq.0) then
		l(5)=char(48+n2)
		go to 65
	endif
	l(5)=char(48+n1)
	l(6)=char(48+n2)
  65	continue
	l(nf-1)='='
c
c   Find equal sign, see if it is preceded by a colon (necol=0 means
c   no, 1 means yes), and form output variable name (vname):
c
   5	do 10 neq=1,lnb
  10    if(l(neq).eq.'=') go to 20
  20	continue
	if(neq.lt.2) go to 98
	necol=0
	n1=neq-1
	if(l(neq-1).eq.':') then
		n1=neq-2
		necol=1
	endif
	if(n1.lt.1) go to 98
	call formnm(l,1,n1,vname)
c
c   see if output name is element of an array (noarr=1 means yes,
c   0 means no), and if so, abort if there is an error:
c
	noarr=0
	if(ickear(vname,x,np,i,ier).eq.1) then
		noarr=1
		if(ier.eq.1) go to 98
	endif
c
c   Handle xxxx='yyyyy':
c
	if(l(neq+1).eq.char(39).and.l(lnb).eq.char(39)) then
		if(neq+1.eq.lnb.or.lnb-neq.gt.15) go to 98
		call formnm(l,neq+1,lnb,args(1))
		call ckchar(args(1),args(2))
		if(args(2).eq.'***************') go to 98
		call ckaddc(vname,args(2),iref)
		go to 99
	endif
c
c   count the number of operations (nops=0 means none and thus an
c   assignment):
c
  25	nops=0
	do 30 i=neq+1,lnb
  30	if(ickop(l(i)).ne.0) nops=nops+1
c
c   handle the case when xxxx=-yyyy (insert '0' and go back to 25):
c
	if(nops.eq.1.and.l(neq+1).eq.'-') then
		lnb=lnb+1
		do 40 i=lnb,neq+1,-1
  40		l(i)=l(i-1)
		l(neq+1)='0'
		go to 25
	endif
c
c   handle assignment statements (other than xxxx='yyyy'):
c
	if(nops.eq.0) then
		if(lnb-neq.lt.1.or.lnb-neq.gt.15) go to 98
		call formnm(l,neq+1,lnb,args(1))
		n1=icktyp(args(1),args(2),n,x,0,nn,m1,iptd)
		if(iptd.eq.1) call clean1(args(1),ier)
		if(n1.eq.0) go to 98
		if(n1.eq.1) then
			call ckaddc(vname,args(2),iref)
			go to 99
		endif
 		if(necol.eq.0) then
 			if(n1.eq.2.or.n1.eq.3) m1=1
 			if(n1.eq.4.or.n1.eq.5) m1=2
 			if(n1.eq.6) m1=3
 		endif
		if(necol.eq.1) then
			if(n1.eq.2.or.n1.eq.3) then
				x=n
				m1=2
			endif
			if(n1.eq.4.or.n1.eq.5) then
				n=x
				m1=1
				if(noarr.eq.1) then
					x=n
					m1=2
				endif
			endif
			if(n1.eq.6) then
				do 50 i=1,nn
  50				wk(i)=int(wk(i))
				lab='Integer Part of an Array'
				m1=3
			endif
		endif
c
c   do storing:
c
     		if(m1.eq.1) call ckaddi(vname,n,iref)
	     	if(m1.eq.2) call ckaddr(vname,x,iref)
		if(m1.eq.3.and.noarr.eq.1) go to 98
		if(m1.eq.3) call ckadda(vname,nn,lab,1,iref)
		go to 99
	endif
c
c
c   Handle operations:
c
c   Find Operation:
c
	nmd=6
	if(ickpop(l,6,nlbd,nrbd,args(1),args(2),ier).eq.1) go to 200
	nmd=5
	if(ickpop(l,5,nlbd,nrbd,args(1),args(2),ier).eq.1) go to 200
	do 198 i=neq+1,lnb
	nmd=ickop(l(i))
	if(nmd.ne.3.and.nmd.ne.4) go to 198
    	if(ickpop(l,nmd,nlbd,nrbd,args(1),args(2),ier).eq.1) go to 200
 198	continue
	do 199 i=neq+1,lnb
	nmd=ickop(l(i))
	if(nmd.ne.1.and.nmd.ne.2) go to 199
    	if(ickpop(l,nmd,nlbd,nrbd,args(1),args(2),ier).eq.1) go to 200
 199	continue
c
c
	go to 98
 200	if(ier.eq.1) go to 98
c
c   Find left and right variable:
c
	n1=icktyp(args(1),args(3),nl,xl,0,nnl,m1,iptd1)
	n2=icktyp(args(2),args(4),nr,xr,nnl,nnr,m2,iptd2)
	if(n1.eq.0.or.n2.eq.0) go to 98
c
c   Get ready for calculation:
c
	nnl=1
	if(m1.eq.1) wk(1)=nl
	if(m1.eq.2) wk(1)=xl
	if(m1.eq.3) then
		call cklist(args(1),xlist,numarr,npl)
		nnl=nend(npl)-nstart(npl)+1
		call movxy(wk,array(nstart(npl)),4*nnl)
		if(iptd1.eq.1) call clean1(args(1),ier)
	endif
	nnr=1
	if(m2.eq.1) wk(nnl+1)=nr
	if(m2.eq.2) wk(nnl+1)=xr
	if(m2.eq.3) then
		call cklist(args(2),xlist,numarr,npr)
		nnr=nend(npr)-nstart(npr)+1
		call movxy(wk(nnl+1),array(nstart(npr)),4*nnr)
		if(iptd2.eq.1) call clean1(args(2),ier)
	endif
	if(nmd.eq.6) then
		if(m2.eq.1) nr=-nr
		do 205 i=1,nnr
 205		wk(nnl+i)=-wk(nnl+i)
		nmd=5
	endif
c
c   Determine type of result:
c
	if(m1.eq.3.or.m2.eq.3) iptadd=3
	if(m1.le.2.and.m2.le.2) iptadd=2
	if(m1.eq.1.and.m2.eq.1) iptadd=1
c
c
	nn=1
	npoint=1
	if(m1.eq.3.and.m2.eq.3) then
		nn=min0(nnl,nnr)
		if(nnr.lt.nnl) npoint=nnl+1
	endif
	if((iptadd.eq.3).and.(m1.lt.3.or.m2.lt.3)) then
		nn=max0(nnl,nnr)
		if(m1.lt.3) npoint=2
	endif
c
c
	nnlp1=nnl+1
	incl=m1/3
	incr=m2/3
	if(nmd.eq.1) then
		do 610 i=0,nn-1
 610		wk(npoint+i)=wk(1+i*incl)+wk(nnlp1+i*incr)
		go to 660
	endif
	if(nmd.eq.2) then
		do 620 i=0,nn-1
 620		wk(npoint+i)=wk(1+i*incl)-wk(nnlp1+i*incr)
		go to 660
	endif
  	if(nmd.eq.3) then
		do 630 i=0,nn-1
 630		wk(npoint+i)=wk(1+i*incl)*wk(nnlp1+i*incr)
		go to 660
	endif
  	if(nmd.eq.4) then
		do 640 i=0,nn-1
		c=wk(nnlp1+i*incr)
		if(c.eq.0.0) then
			write(*,201)
 201			format(' Division by zero')
			imflag=1
			go to 99
		endif
 640		wk(npoint+i)=wk(1+i*incl)/c
		go to 660
	endif
	if(nmd.eq.5) then
		do 650 i=0,nn-1
		cl=wk(1+i*incl)
		cr=wk(nnlp1+i*incr)
		if((cl.eq.0.0.and.cr.lt.0.0).or.
     1		(cl.lt.0.0.and.m2.ne.1)) then
			write(*,203)
 203			format(' Exponentiation error')
			imflag=1
			go to 99
		endif
		if(m2.eq.1) wk(npoint+i)=cl**nr
 650		if(m2.gt.1) wk(npoint+i)=cl**cr
		go to 660
	endif
c
c
 660	continue
	if(iptadd.eq.1.and.abs(wk(npoint)).gt.32767.) then
			write(*,975)
 975			format(' Integer overflow')
			imflag=1
			go to 99
	endif
	ntmp=ntmp+1
	call fnm(ntmp,args(5),'&')
	if(iptadd.eq.1) then
		c=wk(npoint)
		if(c.lt.0.0) nc=int(c-.00001)
		if(c.ge.0.0) nc=int(c+.00001)
		call ckaddi(args(5),nc,iref)
	endif
	if(iptadd.eq.2) call ckaddr(args(5),wk(npoint),iref)
	if(iptadd.eq.3) then
		lab='Arithmetic expression'
		call ckadda(args(5),nn,lab,npoint,iref)
	endif
	if(iref.eq.1) go to 99
	lnb1=lnb
	call mge(l,lnb1,args(5),nlbd,nrbd,ier)
	if(ier.eq.1) go to 98
	lnb=lnb1
	go to 25
c
c   error return:
c
  98	call err2
	go to 101
c
c   handle non-error return:
c
  99    if(nparen.eq.1) then
		lnb=lnbtmp
		call movxy(l,ltmp,lnb)
		lnb1=lnb
		call mge(l,lnb1,vname,nparl,nparr,ier)
		lnb=lnb1
		if(ier.eq.1) go to 98
		nparen=0
		go to 1
	endif
 101	if(ntmp.eq.0) return
	do 97 i=1,ntmp
	call fnm(i,args(1),'&')
  97	call clean1(args(1),ier)
        return
        end
c&CLEAN1
	subroutine clean1(arg,ier)
c********************************************************************
c
c   Subroutine to delete the variable arg.
c
c*******************************************************************
c
#include 'tslabc'
c
	character*15 arg
c
c
	ier=0
        call cklist(arg,xlist,numarr,np)
        if(np.ne.-999) then
		call clean(np)
		go to 99
	endif
	call cklist(arg,rlist,nrvars,np)
	if(np.ne.-999) then
		if(np.eq.nrvars) then
			nrvars=nrvars-1
			go to 99
		endif
                nrvars=nrvars-1
		do 10 j=np,nrvars
		rlist(j)=rlist(j+1)
  10		reals(j)=reals(j+1)
		go to 99
	endif
	call cklist(arg,ilist,nivars,np)
	if(np.ne.-999) then
		if(np.eq.nivars) then
			nivars=nivars-1
			go to 99
		endif
                nivars=nivars-1
		do 20 j=np,nivars
		ilist(j)=ilist(j+1)
  20		ints(j)=ints(j+1)
		go to 99
	endif
	call cklist(arg,chlist,nchars,np)
	if(np.ne.-999) then
		if(np.eq.nchars) then
			nchars=nchars-1
			go to 99
		endif
                nchars=nchars-1
		do 30 j=np,nchars
		chlist(j)=chlist(j+1)
  30		chars(j)=chars(j+1)
	go to 99
	endif
	ier=1
  99	return
	end
c&ICKPOP
	integer*2 function ickpop(l,iopt,nlbd,nrbd,argl,argr,ier)
c******************************************************************
c
c   Function to look for the first occurrence of a specified
c   operation according to:
c
c   iopt: 1(+), 2(-), 3(*), 4(/), 5(^), 6(^-)
c
c   and pull off the arguments on either side of it (argl,argr),
c   and the left and right boundaries of the whole expression
c   (nlbd,nrbd).
c
c*******************************************************************
c
#include 'tslabc'
	character l(1)*1,argl*15,argr*15,l1*1
	integer*2 ickop
c
c
	iopt1=iopt
	if(iopt.eq.6) iopt1=5
	ickpop=0
	do 10 i=1,lnb
    	if(ickop(l(i)).eq.iopt1) then
		if(iopt.eq.6.and.l(i+1).ne.'-') go to 10
		go to 20
	endif
  10	continue
	go to 99
c
c
  20	ickpop=1
	nlr=i-1
	nrl=i+1
	if(iopt.eq.6) nrl=i+2
c
c  make sure we don't have 2 operations in a row and that there
c  is something to the left and right:
c
	if(ickop(l(nlr)).ne.0.or.ickop(l(nrl)).ne.0) then
		ier=1
		go to 99
	endif
	if(nlr-1.lt.1.or.nrl.gt.lnb) then
		ier=1
		go to 99
	endif
c
c   go left until op or =
c
	do 30 nll=nlr-1,1,-1
  30	if(ickop(l(nll)).ne.0.or.l(nll).eq.'=') go to 40
c
c   neither happened:
c
	ier=1
	go to 99
c
c   go right until op or off the end:
c
  40	nll=nll+1
	nrr=lnb
	if(nrl.eq.lnb) go to 50
	do 45 nrr1=nrl+1,lnb
     	if(ickop(l(nrr1)).ne.0) then
		nrr=nrr1-1
		go to 50
	endif
  45	continue
	nrr=lnb
  50	continue
	nlbd=nll
	nrbd=nrr
c
c   make sure left and right stuff aren't too long or too short:
c
	nlg=nlr-nll+1
	nrg=nrr-nrl+1
	if(nlg.lt.1.or.nlg.gt.15.or.nrg.lt.1.or.nrg.gt.15) then
		ier=1
		go to 99
	endif
	call formnm(l,nll,nlr,argl)
	call formnm(l,nrl,nrr,argr)
c
c
  99	continue
	return
	end
c&MGE
	subroutine mge(l,lnb1,arg,nl,nr,ier)
c******************************************************************
c
c   Subroutine to replace the elements l(nl),...,l(nr) with the
c   nonblank part of the character*15 variable arg and update the
c   position (lnb1) of the last nonblank character of l.
c
c*******************************************************************
c
	character l(1)*1,arg*15,arg1*15,arg11(15)*1
	equivalence (arg1,arg11)
c
c
	arg1=arg
	ier=0
c
c
	do 10 nadd=1,15
  10	if(l(nadd).eq.char(32)) go to 20
  20	nadd=nadd-1
	ndel=nr-nl+1
	nchg=nadd-ndel
	if(lnb1+nchg.gt.72) then
		ier=1
		go to 99
	endif
c
c   shift l(nr+1),...,l(lnb1) right (if nchg>0) or left (if nchg<0)
c   nchg places (unless nr=lnb1):
c
	if(nr.eq.lnb1.or.nchg.eq.0) go to 50
	if(nchg.lt.0) then
		do 30 i=nr+1,lnb1
  30		l(i+nchg)=l(i)
		go to 50
	endif
	if(nchg.gt.0) then
		do 40 i=lnb1,nr+1,-1
  40		l(i+nchg)=l(i)
		go to 50
	endif
  50	do 60 i=1,nadd
  60	l(nl+i-1)=arg11(i)
	lnb1=lnb1+nchg
c
c
  99	continue
	return
	end
c&FNM
	subroutine fnm(n,arg,ch)
c*****************************************************************
c
c   Subroutine to form the character*15 variable tmp&n.
c
c******************************************************************
c
	character arg*15,arg1*15,arg11(15)*1,ch*1
	equivalence (arg1,arg11)
c
c
	call movct(arg1,15,char(32))
	arg11(1)='t'
	arg11(2)='m'
	arg11(3)='p'
	arg11(4)=ch
	n1=n/10
	n2=n-10*n1
	n3=1
	if(n1.eq.0) then
		arg11(5)=char(48+n2)
		go to 10
	endif
	arg11(5)=char(48+n1)
	arg11(6)=char(48+n2)
c
c
  10	arg=arg1
	return
	end
c&ICKOP
	integer*2 function ickop(l)
c*******************************************************************
c
c   Function to determine if l is not an operation (0) or else
c   which operation it is:
c
c   +:  1     -:  2     *:  3     /:  4     ^:  5
c
c*******************************************************************
c
	character l*1
	ickop=0
	if(l.eq.'+') ickop=1
	if(l.eq.'-') ickop=2
	if(l.eq.'*') ickop=3
	if(l.eq.'/') ickop=4
	if(l.eq.'^') ickop=5
	return
	end
c&ERR2
	subroutine err2
c********************************************************************
c
c   Subroutine to display error messages in subroutine type2.
c
c*********************************************************************
c
#include 'tslabc'
c
c
	write(*,10)
  10	format(' Illegal Syntax or Term in Type 2 Command')
	imflag=1
c
c
	return
	end
c&ICKEAR
	integer*2 function ickear(vname,x,np,i,ier)
c*****************************************************************
c
c   Function to determine if vname is an element of an array
c   (1 means yes, 0 no), and if so its value (x), which array (np),
c   and which element (i) it is. If any error is encountered,
c   ier is returned as a 1 (0 if no error).
c
c********************************************************************
c
#include 'tslabc'
	character vname*15,arg1*15,arg2*15
c
c
	ickear=0
	ier=0
	call findch(vname,1,15,'[',nocbl,ncbl)
	call findch(vname,1,15,']',nocbr,ncbr)
	if(nocbl+nocbr.eq.0) go to 99
c
c   there is a bracket:
c
	ickear=1
	ier=1
c
c   check for errors:
c
	if(nocbl.ne.1.or.nocbr.ne.1) go to 99
	if(ncbl.le.1.or.ncbl.ge.ncbr) go to 99
	call formnm(vname,1,ncbl-1,arg1)
	call cklist(arg1,xlist,numarr,np)
	if(np.eq.-999) go to 99
	call formnm(vname,ncbl+1,ncbr-1,arg2)
	call ckint(arg2,i)
	if(i.le.0.or.i.gt.nend(np)-nstart(np)+1) go to 99
c
c   got to here with no errors:
c
	x=array(nstart(np)+i-1)
	ier=0
  99	continue
	return
	end
c&ICKTYP
	integer*2 function icktyp(arg1,arg2,n,x,noff,nn,m,iptd)
c******************************************************************
c
c   Function to determine the type of arg1:
c
c      0: error
c      1: character variable
c      2: immediate integer
c      3: integer variable
c      4: immediate real
c      5: real variable or element of an array
c      6: array
c
c   and place its value into arg2 (1), n (2 or 3),
c   x (4 or 5), or wk(1),...,wk(nn) (if an array)
c
c*****************************************************************
c
#include 'tslabc'
	character arg1*15,arg2*15
	integer*2 ickear,ickse
c
c   we start out with icktyp=0 and only change it if we
c   find a type that fits:
c
	iptd=0
	nn=0
	icktyp=0
	m=0
c
c   try character variable:
c
	call cklist(arg1,chlist,nchars,np)
	if(np.ne.-999) then
		icktyp=1
		arg2=chars(np)
		go to 99
	endif
c
c   try immediate integer:
c
	read(arg1,'(bn,i8)',err=20) n
	icktyp=2
	m=1
	go to 99
c
c   try integer variable:
c
  20	continue
	call cklist(arg1,ilist,nivars,np)
	if(np.ne.-999) then
		icktyp=3
		m=1
		n=ints(np)
		go to 99
	endif
c
c   try real variable:
c
	call cklist(arg1,rlist,nrvars,np)
	if(np.ne.-999) then
		icktyp=5
		m=2
		x=reals(np)
		go to 99
	endif
c
c   try element of an array:
c
	if(ickear(arg1,x,np,i,ier).eq.0) go to 35
	if(ier.eq.1) go to 99
	icktyp=5
        m=2
	go to 99
  35	continue
c
c   try pi or exp:
c
	if(arg1.eq.'pi'.or.arg1.eq.'PI') then
		icktyp=5
		m=2
		x=4.*atan(1.0)
		go to 99
	endif
	if(arg1.eq.'exp'.or.arg1.eq.'EXP') then
		icktyp=5
		m=2
		x=exp(1.0)
		go to 99
	endif
c
c   try immediate real:
c
	read(arg1,'(f15.0)',err=40) x
	call findch(arg1,1,15,'.',nocd,ncd)
	if(nocd.eq.0) go to 40
	icktyp=4
	m=2
	go to 99
c
c   try array:
c
  40	continue
	call cklist(arg1,xlist,numarr,np)
	if(np.eq.-999) go to 99
	nsw=noff+nend(np)-nstart(np)+1
	if(ickse(nsw).eq.1) go to 99
	icktyp=6
	m=3
	nn=nend(np)-nstart(np)+1
	call movxy(wk(noff+1),array(nstart(np)),4*nn)
	lab=label(np)
 	call formnm(arg1,1,4,arg2)
 	if(arg2.eq.'tmp&') iptd=1
c
c
  99	continue
	return
	end
