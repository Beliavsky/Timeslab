c&MACRO
c&MACRO
        subroutine macro(args,nargs,lun)
c*************************************************************
c
c   Subroutine to initiate a macro.
c
c   macro(fname)
c
c**************************************************************
c
#include 'tslabc'
        character args(nargs)*15,l15*15,l1(15)*1,arg1(15)*1,arg15*15
        character ename*50,ename1(50)*1
        logical logic
	equivalence (l1,l15),(arg1,arg15),(ename,ename1)
c
c   check whether macro file exists :
c
        iptenv=0
        call cklist(args(1),chlist,nchars,np)
        if(np.gt.0) args(1)=chars(np)
        inquire(file=args(1),exist=logic)
        if(logic) go to 10
	arg15=args(1)
	arg1(12)='.'
	arg1(13)='m'
	arg1(14)='a'
	arg1(15)='c'
	call rmblnk(arg1,15,larg)
	inquire(file=arg15,exist=logic)
	if(logic) then
		args(1)=arg15
		go to 10
	endif
        iptenv=1
        write(*,*) iptenv
        ename='/usr/local/lib/timeslab/'
c        call movxy(ename1(1),'/usr/local/lib/timeslab/',24)
        call movxy(ename1(25),args(1),15)
        call rmblnk(ename,50,larg)
        call movct(ename1(larg+1),50-larg,char(32))
c        write(*,*) ename
c        ename='/usr/local/lib/timeslab/hello.mac'
        write(*,*) ename
        inquire(file=ename,exist=logic)

        if(logic) go to 10
        call movxy(ename1(larg+1),'.mac',4)
        call rmblnk(ename,50,larg)
        write(*,*) ename
        inquire(file=ename,exist=logic)
        if(logic) go to 10        
        write(*,1) args(1),arg15
   1    format(1x,' No MACRO file having name 'a15,' or  ',a15)
	imflag=1
        go to 99
c
c
c
  10    continue
	if(lun.eq.105) then
        	write(*,11)
  11	        format(' Too many open MACRO files')
		imflag=1
	        go to 99
	        endif
	lun=lun+1
        write(*,*) lun
        if(iptenv.eq.0) open(lun,file=args(1))
        if(iptenv.eq.1) open(lun,file=ename)
        rewind lun
	if(nargs.eq.1) go to 99
  25	read(lun,26,end=30) l15
  26	format(a15)
	if(l1(1).eq.';') then
		call movxy(l1,l1(2),14)
		l1(15)=' '
		if(l15.eq.args(2)) go to 99
	endif
	go to 25
  30	write(*,31)
  31	format(' No entry point with this name')
	imflag=1
  99    continue
        return
        end
c&IFMAC
        subroutine ifmac(args,nargs,lun)
c********************************************************************
c
c   Subroutine to do an if statement in a macro.
c
c IF(var,neg,zero,pos)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,l72*72,l(72)*1,vn1*15,vn2*15,vn3*15,
     1	vn4*15,l15*15,l33*3
	equivalence (l,l72)
c
c   check lun:
c
	if(lun.eq.100) then
		write(*,1)
   1		format(' IF illegal except in MACRO')
		imflag=1
		go to 99
	endif
c
c   nargs=1 or 2 case:
c
	if(nargs.eq.1.or.nargs.eq.2) then
		vn4=args(1)
		call lastnb(vn4,15,lnba)
		call findch(vn4,1,lnba,'.',nocp,np1)
		if(nocp.ne.2.and.nocp.ne.3) then
			write(*,351)
 351			format(' Illegal If')
			imflag=1
			go to 99
		endif
		call findch(vn4,np1+1,lnba,'.',nocp,np2)
		call formnm(vn4,1,np1-1,vn1)
		call formnm(vn4,np2+1,lnba,vn2)
		call formnm(vn4,np1+1,np2-1,vn3)
		call ltoup(vn3,15)
		call ckint(vn1,n1)
		call ckint(vn2,n2)
		if(n1.ne.-999.and.n2.ne.-999) then
			x1=n1
			x2=n2
			go to 310
		endif
		call ckreal(vn1,x1)
		call ckreal(vn2,x2)
		if(x1.ne.999..and.x2.ne.999.) go to 310
		write(*,351)
		imflag=1
		go to 99
 310	if(vn3.eq.'LT'.or.vn3.eq.'LE'.or.vn3.eq.'EQ') go to 315
	if(vn3.eq.'GE'.or.vn3.eq.'GT'.or.vn3.eq.'NE') go to 315
	write(*,351)
	imflag=1
	go to 99
 315    if(vn3.eq.'LT'.and.x1.lt.x2) go to 320
	if(vn3.eq.'LE'.and.x1.le.x2) go to 320
	if(vn3.eq.'EQ'.and.x1.eq.x2) go to 320
	if(vn3.eq.'GE'.and.x1.ge.x2) go to 320
	if(vn3.eq.'GT'.and.x1.gt.x2) go to 320
	if(vn3.eq.'NE'.and.x1.ne.x2) go to 320
c
c   condition false, if nargs=2 we're finished, if nargs=1,
c   we must find the matching ENDIF:
c
	if(nargs.eq.2) go to 99
	nif1=1
	nend1=0
 316	read(lun,41,end=317) l72
	call rmblnk(l72,72,lnb72)
	call ltoup(l72,5)
	if(l72.eq.'ENDIF') then
		nend1=nend1+1
		if(nend1.eq.nif1) go to 99
		go to 316
	endif
	call movxy(l33,l72,3)
	if(l33.eq.'IF(') then
		call findch(l72,1,lnb72,',',nocc,nccoma)
		if(nocc.gt.0) go to 316
		nif1=nif1+1
		go to 316
	endif
	go to 316
 317	write(*,318)
 318	format(' No matching ENDIF')
	imflag=1
	go to 99
c
c   condition true, if nargs=1, we're finished, otherwise look
c   for target of IF:
c
 320	if(nargs.eq.1) go to 99
    	rewind lun
 325	read(lun,41,end=330) l72
	if(l(1).ne.';') go to 325
	call movxy(l15,l(2),15)
	if(l15.ne.args(2)) go to 325
	go to 99
 330	write(*,331)
 331	format(' Target of IF doesn''t exist')
	imflag=1
	go to 99
	endif
c
c
c   nargs=4:
c
c  see if var is in ilist or rlist :
c
        call ckint(args(1),np)
        if(np.eq.-999) go to 10
        if(np.eq.0) ioptsn=0
        if(np.lt.0) ioptsn=-1
        if(np.gt.0) ioptsn=1
        go to 25
  10    call ckreal(args(1),xp)
        if(xp.ne.-999.) go to 20
        call error(args,1,2)
        go to 99
  20    if(xp.eq.0.0) ioptsn=0
        if(xp.gt.0.0) ioptsn=1
        if(xp.lt.0.0) ioptsn=-1
c
c   check arguments 2-4:
c
  25    continue
        nck=3+ioptsn
        call ckint(args(nck),np)
        if(np.ne.-999) go to 30
        call error(args,nck,2)
        go to 99
c
c   reposition file:
c
  30    npm1=np-1
	if(npm1.eq.0) go to 99
	if(npm1.gt.0) then
		do 40 i=1,npm1
  45		read(lun,41,end=61) l72
  41		format(a72)
  40		continue
		go to 99
		endif
	do 51 i=1,1-np
  51	backspace lun
	go to 99
  61	lun=lun-1
c
c
  99    continue
        return
        end
