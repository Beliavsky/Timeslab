c&CHECK
c&CKCHAR
        subroutine ckchar(vn1,vn2)
c************************************************************
c
c   Subroutine to check whether the character*15 variable is
c   an immediate character variable (i.e. starts and ends with
c   quotes) or included in the set of defined character variables.
c
c*************************************************************
c
#include 'tslabc'
        character vn1*15,vn2*15,v(15)*1,v2*15,v3*13,quote*1
	character vn5*15,vn6*6,vn7(6)*1,vn8(15)*1,vn9*15
	equivalence (vn6,vn7),(vn8,vn9)
        equivalence (v,v2),(v(2),v3)
        data quote / ''''/
c
c   check whether in quotes :
c
        v2 = vn1
        call lastnb(v,15,n)
        if(v(1).ne.quote.or.v(n).ne.quote) go to 10
        v(n)=' '
c
c   yes:
c
	call findch(v,1,n,'#',nocc,nc1)
	if(nocc.eq.2) then
		call findch(v,nc1+1,n,'#',nocc1,nc2)
		if(nc2.eq.nc1+1) go to 5
		call formnm(v,nc1+1,nc2-1,vn5)
		call ckint(vn5,np)
		if(np.eq.-999) go to 5
		write(vn6,6,err=5) np
   6		format(i6)
		call lastnb(vn7,6,n1)
		call movct(vn8,15,char(32))
		do 7 i=2,nc1-1
   7		vn8(i-1)=v(i)
		jj=nc1-2
		do 8 i=1,n1
		if(vn7(i).eq.' ') go to 8
		jj=jj+1
    		vn8(jj)=vn7(i)
   8		continue
		if(nc2.ne.15) then
			do 9 i=nc2+1,15
			if(v(i).eq.' ') go to 11
			jj=jj+1
			vn8(jj)=v(i)
   9			continue
  11			continue
		endif
		vn2=vn9
		go to 99
	endif
   5	continue
        if(n.le.2) then
        vn2='***************'
        else
        vn2 = v3
        endif
        go to 99
c
c   no :
c
  10    continue
        call cklist(vn1,chlist,nchars,np)
        if(np.eq.-999) vn2='***************'
        if(np.gt.0) vn2=chars(np)
c
c
  99    continue
        return
        end
c&CKREAL
        subroutine ckreal(arg,x)
c********************************************************************
c
c   Subroutine to see if the character*15 variable arg is a real number
c   or an element of a list of real variable names. If not then the
c   real variable x is returned as 999.. If yes, then x is returned
c   as either the value of the real number or the variable corresponding
c   to the variable name.
c
c********************************************************************
c
#include 'tslabc'
        character*15 arg,arg15,vn1,vn2
        character*1 arg1(15)
        equivalence (arg15,arg1)
        arg15=arg
c
c   if the first character in arg is  -,.,or 0-9 we assume it's
c   an immediate value:
c
        ii=ichar(arg1(1))
        if(ii.ge.48.and.ii.le.57) go to 50
        if(ii.eq.45.or.ii.eq.46) go to 50
c
c   maybe it's a variable name:
c
        call cklist(arg,rlist,nrvars,np)
        if(np.gt.0) then
		x=reals(np)
		go to 99
	endif
c
c   maybe it's an integer variable:
c
	call cklist(arg,ilist,nivars,np)
	if(np.gt.0) then
		x=ints(np)
		go to 99
	endif
c
c   maybe it's PI:
c
	if(arg.eq.'pi'.or.arg.eq.'PI') then
		x=4.*atan(1.0)
		go to 99
	endif
c
c   maybe it's exp:
c
	if(arg.eq.'exp'.or.arg.eq.'EXP') then
		x=exp(1.0)
		go to 99
	endif
c
c   maybe it's x[i]:
c
	call findch(arg,1,15,'[',noclb,nclb)
	if(noclb.ne.1) go to 10
	call findch(arg,1,15,']',nocrb,ncrb)
	if(nocrb.ne.1.or.ncrb.le.nclb) go to 10
	call formnm(arg,1,nclb-1,vn1)
	call formnm(arg,nclb+1,ncrb-1,vn2)
	call cklist(vn1,xlist,numarr,np)
	if(np.eq.-999) go to 10
	call ckint(vn2,i)
	if(i.eq.-999.or.i.lt.1.or.i.gt.nend(np)-nstart(np)+1) go to 10
c
c   yes:
c
	x=array(nstart(np)+i-1)
	go to 99
c
c   no:
c
  10	continue
        x=999.
        go to 99
c
c   assume it's a number(make sure there's a decimal point):
c
  50    continue
        do 31 i=1,15
  31    if(arg1(i).eq.'.') go to 32
ccc
ccc
	do 42 i=15,1,-1
    	if(ichar(arg1(i)).lt.48.or.ichar(arg1(i)).gt.57) go to 42
	go to 43
  42	continue
  43	if(i.ne.15) arg1(i+1)='.'
c	go to 60
ccc
ccc
  32    continue
        read(arg15,'(f15.0)',err=60) x
        go to 99
  60    x=999.
  99    continue
        return
        end
c&CKLIST
        subroutine cklist(name,nlist,n,np)
c*********************************************************************
c
c   Subroutine to check whether the character*15 variable name is
c   included in the character*15 array nlist(1),...,nlist(n). If so,
c   then its position is returned in np. Otherwise np is set to -999.
c   if n=0, then np is set to -999.
c
c*********************************************************************
c
        character*15 name,nlist(1)
        np=-999
        do 10 i=1,n
  10    if(name.eq.nlist(i)) go to 20
        return
  20    np=i
        return
        end
c&CKINT
        subroutine ckint(arg,n)
c*********************************************************************
c
c   Subroutine to see if the character*15 variable arg is an integer
c   number or an element of the list of integer variable names. If
c   not, the output integer n is set to -999. If yes then n is given
c   either the value corresponding to the variable name or else the
c   value of the number contained in arg.
c
c********************************************************************
c
#include 'tslabc'
        character*15 arg
        character*1 arg1
c
c   see if arg is a number :
c
        call movxy(arg1,arg,1)
        ii=ichar(arg1)
        if(48.le.ii.and.57.ge.ii) go to 50
        if(ii.eq.45) go to 50
c
c   no :
c
        call cklist(arg,ilist,nivars,np)
        if(np.gt.0) go to 20
        n=-999
        go to 99
  20    n=ints(np)
        go to 99
c
c   yes :
c
  50    continue
        read(arg,'(bn,i8)',err=60) n
        go to 99
  60    n=-999
  99    continue
        return
        end
c&CKADDA
        subroutine ckadda(vname,n,lab,m,iref)
c******************************************************************
c
c   Subroutine to check whether an array exists and add it if the
c   user wants. Iref is returned as a 1 if the user refuses to
c   overwrite and a 0 otherwise.
c
c********************************************************************
c
#include 'tslabc'
        character vname*15,i10*1
	integer*2 icka,ickvn,sswrqq
c
c
c
        iref=1
	n1=sswrqq()
	call ckexc(n1,ier)
	if(ier.ne.0) then
		imflag=1
		write(*,1) ier
   1		format(' Arithmetic Error ',i1)
		go to 99
	endif
	if(ickvn(vname).eq.1) go to 99
	if(icka(vname,1).eq.1) go to 99
        call cklist(vname,xlist,numarr,np)
        if(np.eq.-999) go to 20
        if(ioptov.eq.1) call clean(np)
        if(ioptov.eq.1) go to 20
        write(*,10) vname
  10    format(1x,a15,'already exists, overwrite it ?'$)
        read(*,15) i10
  15    format(a1)
        if(i10.ne.'Y'.and.i10.ne.'y') go to 99
        call clean(np)
c
c   see if array will fit :
c
  20	if(n.gt.num.or.numarr.ge.maxarr) then
		write(*,30)
  30	format(' Command Aborted-Too many arrays or too many elements')
 		imflag=1
		iref=1
		go to 99
	endif
c
c   add array :
c
        call addarr(n,vname,lab,m)
        iref=0
c
c
  99    continue
        return
        end
c&ICKA
	integer*2 function icka(vname,ntype)
c*******************************************************************
c
c   Function to determine if vname has been defined as a type other
c   than type ntype (1:array,2:integer,3:real,4:char)
c
c   icka=0 if not, 1 if so
c
c*******************************************************************
c
#include 'tslabc'
	character vname*15
c
	icka=1
	if(ntype.eq.1) go to 10
	call cklist(vname,xlist,numarr,np)
	if(np.ne.-999) go to 99
  10	if(ntype.eq.2) go to 20
	call cklist(vname,ilist,nivars,np)
	if(np.ne.-999) go to 99
  20	if(ntype.eq.3) go to 30
	call cklist(vname,rlist,nrvars,np)
	if(np.ne.-999) go to 99
  30	if(ntype.eq.4) go to 40
	call cklist(vname,chlist,nchars,np)
	if(np.ne.-999) go to 99
  40	icka=0
  99	if(icka.eq.1) then
		write(*,1) vname
   1	format(1x,
     1	'You already have a variable of a different type called ',a15)
		imflag=1
	endif
	return
	end
c&CKADDI
        subroutine ckaddi(vname,ival,iref)
c******************************************************************
c
c   Subroutine to check whether an integer exists and add it if the
c   user wants. Iref is returned as a 1 if the user refuses to
c   overwrite and a 0 otherwise.
c
c********************************************************************
c
#include 'tslabc'
        character vname*15,i10*1
	integer*2 icka,ickvn,sswrqq
c
c
c
	iref=1
	n1=sswrqq()
	call ckexc(n1,ier)
	if(ier.ne.0) then
		imflag=1
		write(*,1) ier
   1		format(' Arithmetic Error ',i1)
		go to 99
	endif
	if(ickvn(vname).eq.1) go to 99
	call findch(vname,1,15,'[',nocb,ncb)
	if(nocb.ne.0) then
		call ckaddr(vname,float(ival),iref)
		go to 99
	endif
c
c
        iref=1
	if(icka(vname,2).eq.1) go to 99
        call cklist(vname,ilist,nivars,np)
        if(np.eq.-999) go to 20
        if(ioptov.eq.1) ints(np)=ival
        if(ioptov.eq.1) go to 30
        write(*,10) vname
  10    format(1x,a15,'already exists, overwrite it ?'$)
        read(*,15) i10
  15    format(a1)
        if(i10.ne.'Y'.and.i10.ne.'y') go to 99
        ints(np)=ival
        go to 30
  20	if(nivars.eq.maxiv) then
		write(*,29)
  29		format(' Command Aborted-Too many integer variables')
		imflag=1
		iref=1
		go to 99
		endif
        nivars=nivars+1
        ints(nivars)=ival
        ilist(nivars)=vname
  30    continue
        iref=0
c
c
  99    continue
        return
        end
c&CKADDR
        subroutine ckaddr(vname,rval,iref)
c******************************************************************
c
c   Subroutine to check whether a real exists and add it if the
c   user wants. Iref is returned as a 1 if the user refuses to
c   overwrite and a 0 otherwise.
c
c********************************************************************
c
#include 'tslabc'
        character vname*15,i10*1,vn1*15,vn2*15
	integer*2 icka,ickvn,sswrqq
c
c   see if vname=x[i]:
c
	iref=1
	n1=sswrqq()
	call ckexc(n1,ier)
	if(ier.ne.0) then
		imflag=1
		write(*,1) ier
   1		format(' Arithmetic Error ',i1)
		go to 99
	endif
	if(ickvn(vname).eq.1) go to 99
	call findch(vname,1,15,'[',noclb,nclb)
	if(noclb.ne.1) go to 5
	call findch(vname,1,15,']',nocrb,ncrb)
	if(nocrb.ne.1.or.ncrb.le.nclb) go to 5
	call formnm(vname,1,nclb-1,vn1)
	call formnm(vname,nclb+1,ncrb-1,vn2)
	call cklist(vn1,xlist,numarr,np)
	if(np.eq.-999) go to 5
	call ckint(vn2,i)
	if(i.eq.-999.or.i.lt.1.or.i.gt.nend(np)-nstart(np)+1) go to 5
c
c   yes:
c
	array(nstart(np)+i-1)=rval
	iref=0
	go to 99
c
c   no:
c
   5	continue
        iref=1
	if(icka(vname,3).eq.1) go to 99
        call cklist(vname,rlist,nrvars,np)
        if(np.eq.-999) go to 20
        if(ioptov.eq.1) reals(np)=rval
        if(ioptov.eq.1) go to 30
        write(*,10) vname
  10    format(1x,a15,'already exists, overwrite it ?'$)
        read(*,15) i10
  15    format(a1)
        if(i10.ne.'Y'.and.i10.ne.'y') go to 99
        reals(np)=rval
        go to 30
  20	if(nrvars.eq.maxrv) then
		write(*,29)
  29		format(' Command Aborted-Too many real variables')
		imflag=1
		iref=1
		go to 99
		endif
        nrvars=nrvars+1
        reals(nrvars)=rval
        rlist(nrvars)=vname
  30    continue
        iref=0
c
c
  99    continue
        return
        end
c&CKADDC
        subroutine ckaddc(vname,cval,iref)
c******************************************************************
c
c   Subroutine to check whether a character exists and add it if the
c   user wants. Iref is returned as a 1 if the user refuses to
c   overwrite and a 0 otherwise.
c
c********************************************************************
c
#include 'tslabc'
        character vname*15,i10*1,cval*15
	integer*2 icka,ickvn
c
c
c
        iref=1
	if(ickvn(vname).eq.1) go to 99
	if(icka(vname,4).eq.1) go to 99
        call cklist(vname,chlist,nchars,np)
        if(np.eq.-999) go to 20
        if(ioptov.eq.1) chars(np)=cval
        if(ioptov.eq.1) go to 30
        write(*,10) vname
  10    format(1x,a15,'already exists, overwrite it ?'$)
        read(*,15) i10
  15    format(a1)
        if(i10.ne.'Y'.and.i10.ne.'y') go to 99
        chars(np)=cval
        go to 30
  20	if(nchars.eq.maxchr) then
		write(*,29)
  29		format(' Command Aborted-Too many character variables')
		imflag=1
		iref=1
		go to 99
		endif
        nchars=nchars+1
        chars(nchars)=cval
        chlist(nchars)=vname
  30    continue
        iref=0
c
c
  99    continue
        return
        end
c&ICKSE
	integer*2 function ickse(n)
c*******************************************************************
c
c   Function to check if ther's enough room in wk to do a command.
c
c  0:yes  1:no
c
c********************************************************************
c
#include 'tslabc'
c
	ickse=0
	if(n.gt.maxnum) then
		ickse=1
		imflag=1
		write(*,10)
  10		format(' Not enough work space to do this command')
        endif
	return
	end
c&ICKVN
	integer*2 function ickvn(vname)
c***************************************************************
c
c   Function to check whether the first character of vname is a Letter
c
c   1:no  0:yes
c
c****************************************************************
c
#include 'tslabc'
	character vname*15,v1(15)*1,vn1*15
	equivalence (vn1,v1)
c
c
	vn1=vname
	ickvn=1
	ii=ichar(v1(1))
	if(ii.ge.65.and.ii.le.90) ickvn=0
	if(ii.ge.97.and.ii.le.122) ickvn=0
	if(ickvn.eq.1) then
		write(*,10)
  10		format(' Variable names must start with a-z or A-Z')
		imflag=1
		go to 99
	endif
c
c
  99	continue
	return
	end
c&ICKL
 	integer*2 function ickl(name,npx,nx)
c*****************************************************************
c
c   function to check whether the array having name name exists and
c   if so to return its position (npx) and length nx
c
c   ickl = 1  if name not in list
c          0  if name is in list
c
c*******************************************************************
c
#include 'tslabc'
 	character name*15
c
c
 	ickl=0
 	call cklist(name,xlist,numarr,npx)
 	if(npx.eq.-999) then
 		call error(name,1,1)
 		ickl=1
 		go to 99
 		endif
	nx=nend(npx)-nstart(npx)+1
c
c
  99	continue
 	return
 	end
c&ICKI
	integer*2 function icki(name,narg,n,iopt,nck)
c********************************************************************
c
c   Function to convert a character*15 name to an integer n and check
c   whether it is in a certain range.
c
c   iopt: 0 only check for existence
c         1 check if n.lt.nck
c         2 check if n.le.nck
c         3 check if n.eq.nck
c         4 check if n.ge.nck
c         5 check if n.gt.nck
c
c   icki: 0 means everything ok
c         1 means check is true
c
c********************************************************************
c
#include 'tslabc'
	character name*15
c
	icki=0
	call ckint(name,n)
	if(n.eq.-999) then
		write(*,1) name
   1		format(1x,a15,' not a valid integer')
		icki=1
		imflag=1
		go to 100
	endif
	if(iopt.eq.0) go to 99
	go to(10,20,30,40,50) iopt
  10	if(n.lt.nck) icki=1
	go to 99
  20	if(n.le.nck) icki=1
	go to 99
  30	if(n.eq.nck) icki=1
	go to 99
  40	if(n.ge.nck) icki=1
	go to 99
  50	if(n.gt.nck) icki=1
  99	if(icki.eq.1) then
		write(*,2) name,narg
   2		format(1x,a15,' illegal value for argument ',i5)
		imflag=1
	endif
 100	continue
	return
	end
c&ICKR
	integer*2 function ickr(name,narg,x,iopt,xck)
c********************************************************************
c
c   Function to convert a character*15 name to a real x and check
c   whether it is in a certain range.
c
c   iopt: 0 only check for existence
c         1 check if x.lt.xck
c         2 check if x.le.xck
c         3 check if x.eq.xck
c         4 check if x.ge.xck
c         5 check if x.gt.xck
c
c   icki: 0 means everything ok
c         1 means check is true
c
c********************************************************************
c
#include 'tslabc'
	character name*15
c
	ickr=0
	call ckreal(name,x)
	if(x.eq.999.) then
		write(*,1) name
   1		format(1x,a15,' not a valid real number')
		ickr=1
		imflag=1
		go to 100
	endif
	if(iopt.eq.0) go to 99
	go to(10,20,30,40,50) iopt
  10	if(x.lt.xck) ickr=1
	go to 99
  20	if(x.le.xck) ickr=1
	go to 99
  30	if(x.eq.xck) ickr=1
	go to 99
  40	if(x.ge.xck) ickr=1
	go to 99
  50	if(x.gt.xck) ickr=1
  99	if(ickr.eq.1) then
		write(*,2) name,narg
   2		format(1x,a15,' illegal value for argument ',i5)
		imflag=1
	endif
 100	continue
	return
	end
