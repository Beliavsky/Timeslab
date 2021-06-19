c&TYPE3
c&TYPE3
        subroutine type3(l)
c****************************************************************
c
c   Subroutine to process type 3 commands :
c
c****************************************************************
c
#include 'tslabc'
        character l(72)*1,pname*15,args(21)*15,vname*15
c
c   find args and pname, convert pname to upper case, and go to
c   correct procedure :
c
        call parse(l,pname,args,nargs,vname,ierr)
	if(ierr.eq.1) then
		write(*,6701)
6701	format(' An argument or command name is too long or too short')
		imflag=1
		go to 99
	endif
        call ltoup(pname,15)
        call cklist(pname,pr3lst,numpr3,np)
        if(np.gt.0.and.np.le.numpr3) go to 3002
        write(*,3001)
3001    format(' undefined type 3 command ')
	imflag=1
        go to 99
3002    go to(3005,3010,3010,3020,3030,3040,3050,3060,3070,3080,3090,
     1  3100,3110,3120,3130,3140,3150,3160,3170,3180,3190,3200,3210,
     1  3220,3230,3240,3250,3260,3270,3280,3290,3300,3310,3320,3330,
     1  3340,3350,3360,3370,3380,3390,3400,3410,3420,3430,3440,
     1	3450,3460,3470,3480,3490,3500,3510,3520) np
c
c   procedure WN :
c
3005    if(nargs.ne.2.and.nargs.ne.3) go to 98
        call wn(args,nargs,vname)
        go to 99
3010    continue
c
c   procedure COS or procedure SIN :
c
        if(nargs.lt.1.or.nargs.gt.3) go to 98
        call trig(args,nargs,vname,np-1)
        go to 99
3020    continue
c
c   procedure DIFF :
c
        if(nargs.ne.3) go to 98
        call diff(args,nargs,vname)
        go to 99
3030    continue
c
c   procedure CORR :
c
        if(nargs.ne.7) go to 98
        call corr(args,nargs,vname)
        go to 99
3040    continue
c
c   procedure LOGE :
c
        if(nargs.ne.2.and.nargs.ne.1) go to 98
        call loge(args,nargs,vname)
        go to 99
3050    continue
c
c   procedure WINDOW :
c
        if(nargs.ne.5.and.nargs.ne.6.and.nargs.ne.7) go to 98
        call window(args,nargs,vname)
        go to 99
3060    continue
c
c   procedure SUBMNS :
c
        if(nargs.ne.4.and.nargs.ne.5) go to 98
        call submns(args,nargs,vname)
        go to 99
3070    continue
c
c   procedure DTAR :
c
        if(nargs.ne.9.and.nargs.ne.10) go to 98
        call dtar(args,nargs,vname)
        go to 99
3080    continue
c
c   procedure ARSP :
c
        if(nargs.ne.4) go to 98
        call armasp(args,nargs,vname,1)
        go to 99
3090    continue
c
c   procedure CUMSP :
c
        if(nargs.ne.2) go to 98
        call cumsp(args,nargs,vname)
        go to 99
3100    continue
c
c   procedure MASP :
c
        if(nargs.ne.4) go to 98
        call armasp(args,nargs,vname,2)
        go to 99
3110    continue
c
c   procedure ARMASP :
c
        if(nargs.ne.6) go to 98
        call armasp(args,nargs,vname,3)
        go to 99
3120    continue
c
c   procedure EXTRACT :
c
        if(nargs.ne.3.and.nargs.ne.5.and.nargs.ne.6) go to 98
        call extrct(args,nargs,vname)
        go to 99
3130    continue
c
c   procedure ARFILT:
c
        if(nargs.ne.7) go to 98
        call arfilt(args,nargs,vname)
        go to 99
3140    continue
c
c   procedure ARDT2 :
c
        if(nargs.ne.5.and.nargs.ne.6) go to 98
        call ardt2(args,nargs,vname)
        go to 99
3150    continue
c
c   procedure EXP :
c
        if(nargs.ne.2.and.nargs.ne.1) go to 98
        call expvec(args,nargs,vname)
        go to 99
3160    continue
c
c   procedure BARTTEST :
c
        if(nargs.ne.4.and.nargs.ne.1) go to 98
        call bart(args,nargs,vname)
        go to 99
3170    continue
c
c   procedure ARDT :
c
        if(nargs.ne.7.and.nargs.ne.4) go to 98
        call ardt(args,nargs,vname)
        go to 99
3180    continue
c
c   procedure QTEST :
c
        if(nargs.ne.6) go to 98
        call qtest(args,nargs,vname)
        go to 99
3190    continue
c
c   procedure EXTEND :
c
        if(nargs.ne.4.and.nargs.ne.5) go to 98
        call extend(args,nargs,vname)
        go to 99
3200    continue
c
c   procedure SORT :
c
        if(nargs.ne.1.and.nargs.ne.2) go to 98
        call sort(args,nargs,vname)
        go to 99
3210    continue
c
c   procedure PARTAR :
c
        if(nargs.ne.2) go to 98
        call partal(args,nargs,vname)
        go to 99
3220    continue
c
c   procedure LINE :
c
        if(nargs.ne.3.and.nargs.ne.1.and.nargs.ne.4) go to 98
        call formln(args,nargs,vname)
        go to 99
3230    continue
c
c   procedur ARMACORR :
c
        if(nargs.ne.8) go to 98
        call armacr(args,nargs,vname)
        go to 99
3240    continue
c
c   procedure CORRMA :
c
        if(nargs.ne.7) go to 98
        call corrma(args,nargs,vname)
        go to 99
3250    continue
c
c   procedure MACORR :
c
        if(nargs.ne.5) go to 98
        call macorr(args,nargs,vname)
        go to 99
3260    continue
c
c   procedure CORRARMA :
c
        if(nargs.ne.9) go to 98
        call crarma(args,nargs,vname)
        go to 99
3270    continue
c
c   procedure MADT :
c
        if(nargs.ne.5) go to 98
        call madt(args,nargs,vname)
        go to 99
3280    continue
c
c   procedure ARMADT :
c
        if(nargs.ne.11) go to 98
        call armadt(args,nargs,vname)
        go to 99
3290    continue
c
c   procedure MULTPOLY :
c
        if(nargs.ne.4) go to 98
        call mtpoly(args,nargs,vname)
        go to 99
3300    continue
c
c   procedure INVPOLY :
c
        if(nargs.ne.3) go to 98
        call inpoly(args,nargs,vname)
        go to 99
3310    continue
c
c   procedure CORRAR :
c
        if(nargs.ne.4.and.nargs.ne.7) go to 98
        call corrar(args,nargs,vname)
        go to 99
3320    continue
c
c   procedure ARCORR :
c
        if(nargs.ne.6) go to 98
        call arcorr(args,nargs,vname)
        go to 99
3330    continue
c
c   procedure DIVSDS :
c
        if(nargs.ne.4.and.nargs.ne.5) go to 98
        call divsds(args,nargs,vname)
        go to 99
3340    continue
c
c   procedure ARPART :
c
        if(nargs.ne.3) go to 98
        call arpart(args,nargs,vname)
        go to 99
3350    continue
c
c   procedure TRANS :
c
        if(nargs.ne.3) go to 98
        call trans(args,nargs,vname)
        go to 99
3360    continue
c
c   procedure DIST :
c
        if(nargs.ne.3.and.nargs.ne.4.and.nargs.ne.5.and.nargs.ne.6)
     1	go to 98
        call dist(args,nargs,vname)
        go to 99
3370    continue
c
c   procedure POLY :
c
        if(nargs.ne.4.and.nargs.ne.5) go to 98
        call poly(args,nargs,vname)
        go to 99
3380    continue
c
c   procedure DOT :
c
        if(nargs.ne.3) go to 98
        call dot(args,nargs,vname)
        go to 99
3390	continue
c
c   procedure MINV:
c
	if(nargs.ne.3) go to 98
	call minv(args,nargs,vname)
	go to 99
3400	continue
c
c   procedure MDEL:
c
	call mdel(args,nargs,vname)
	go to 99
3410	continue
c
c   procedure CORRAR2:
c
	if(nargs.ne.10.and.nargs.ne.14) go to 98
	call corar2(args,nargs,vname)
	go to 99
3420	continue
c
c    procedure MMULT:
c
	if(nargs.lt.3.or.nargs.gt.5) go to 98
	call mmult(args,nargs,vname)
	go to 99
3430	continue
c
c   procedure TOEPL:
c
	if(nargs.ne.3) go to 98
	call toepl(args,nargs,vname)
	go to 99
3440	continue
c
c   procedure FILT:
c
	if(nargs.ne.5) go to 98
	call filt(args,nargs,vname)
	go to 99
3450	continue
c
c   procedure REVERSE:
c
	if(nargs.ne.2) go to 98
	call revrse(args,nargs,vname)
	go to 99
3460	continue
c
c   procedure EIG:
c
	if(nargs.ne.4.and.nargs.ne.5) go to 98
	call eig(args,nargs,vname)
	go to 99
3470	continue
c
c   procedure REPLACE:
c
	if(nargs.lt.4.or.nargs.gt.6) go to 98
	call repl(args,nargs,vname)
	go to 99
3480	continue
c
c   procedure CUM:
c
	if(nargs.ne.3) go to 98
	call cum(args,nargs,vname)
	go to 99
3490	continue
c
c   procedure PARCORR:
c
	if(nargs.ne.3.and.nargs.ne.4) go to 98
	call parcor(args,nargs,vname)
	go to 99
3500	continue
c
c   procedure SWEEP:
c
	if(nargs.ne.5) go to 98
	call sweep(args,nargs,vname)
	go to 99
3510	continue
c
c   procedure ABS:
c
	if(nargs.ne.1.and.nargs.ne.2) go to 98
	call mabs(args,nargs,vname)
	go to 99
3520	continue
c
c   procedure BINOM:
c
	if(nargs.ne.2.and.nargs.ne.3) go to 98
	call binom(args,nargs,vname)
	go to 99
c
c
c    ERROR FINISH:
c
c
  98    write(*,101) pr3lst(np)
 101    format(' wrong number of arguments in procedure ',a15)
	imflag=1
c
c   finish :
c
  99    continue
        return
        end
