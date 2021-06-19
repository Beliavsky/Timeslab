c&T1NZ
c&READFL
        subroutine readfl(args,nargs)
c********************************************************************
c
c   Subroutine to read a file and put the information into the proper
c   places.
c
c*********************************************************************
c
#include 'tslabc'
        character args(nargs)*15,lform*30,l1(72)*1,nname*15
	character arg15*15,arg1(15)*1
        logical logic
	integer*2 icki
	equivalence (arg1,arg15)
c
c   check if file exits :
c
        inquire(file=args(1),exist=logic,err=11)
        if(logic) go to 20
        call cklist(args(1),chlist,nchars,np)
        if(np.eq.-999) go to 15
        args(1)=chars(np)
        inquire(file=args(1),exist=logic,err=11)
        if(logic) go to 20
  15    continue
	arg15=args(1)
	arg1(12)='.'
	arg1(13)='d'
	arg1(14)='a'
	arg1(15)='t'
	call rmblnk(arg1,15,larg)
	inquire(file=arg15,exist=logic)
	if(logic) then
		args(1)=arg15
		go to 20
	endif
c
c   no :
c
  11    write(*,10) args(1),arg15
  10    format(' No file named ',a15,' or ',a15)
	imflag=1
        go to 99
c
c   yes :
c
c
  20	nskip=0
     	if(nargs.eq.4) then
		if(icki(args(4),4,nskip,1,0).eq.1) go to 99
 	endif
        open(2,file=args(1),status='old')
	do 701 jd=1,nskip+1
        read(2,30,err=70) lab
  30    format(a40)
        read(2,40,err=70) l1
  40    format(72a1)
c
c   find n :
c
        ioptf=0
        do 41 i=1,72
  41    if(l1(i).ne.' ') go to 42
        go to 70
  42    do 43 j=i,72
  43    if(l1(j).eq.' ') go to 44
        j=72
        go to 55
  44    j=j-1
        call formnm(l1,i,j,nname)
        call ckint(nname,n)
        if(n.lt.1.or.n.gt.maxnum) go to 70
c
c   look for format statement :
c
        do 45 k=j+1,72
  45    if(l1(k).ne.' ') go to 46
        go to 55
  46    do 47 kk=k+1,72
  47    if(l1(kk).eq.' ') go to 48
        go to 70
  48    kk=kk-1
        if(l1(k).ne.'('.or.l1(kk).ne.')') go to 70
        if(kk-k+1.gt.30) go to 70
        ioptf=1
        call movct(lform,30,' ')
        call movxy(lform,l1(k),kk-k+1)
  55    continue
        if(ioptf.eq.0) read(2,*,end=70,err=70) (wk(i),i=1,n)
        if(ioptf.eq.1) read(2,lform,end=70,err=70) (wk(i),i=1,n)
        go to 701
  70    write(*,80)
  80    format(' file not set up correctly')
	imflag=1
        go to 99
c
c   add array :
c
 701	continue
        close(2)
        call ckadda(args(2),n,lab,1,iref)
        if(iref.eq.1) go to 99
c
c   add sample size :
c
        call ckaddi(args(3),n,iref)
  99    continue
        return
        end
c&SAVE
        subroutine save(args,nargs)
c*****************************************************************
c
c   Subroutine to save one of the defined arrays onto a disk file.
c
c   save(x,n,fname)
c
c*****************************************************************
c
#include 'tslabc'
        logical logic
        character args(nargs)*15,form*20
	integer*2 ickl
        form='(5e14.7)'
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),n)
        if(n.le.0.or.n.gt.nn) then
		call error(args,2,2)
		go to 99
	endif
        call cklist(args(3),chlist,nchars,np1)
        if(np1.ne.-999) args(3)=chars(np1)
        inquire(file=args(3),exist=logic,err=60)
	iptn=0
        if(logic) iptn=1
	if(iptn.eq.1) then
		open(2,file=args(3),status='old',err=60)
   5		read(2,6,end=7) i10
   6		format(a1)
		go to 5
   7		backspace 2
	endif
	if(iptn.eq.0) open(2,file=args(3),status='new',err=60)
c
c
        write(2,80,err=60) label(np)
  80    format(a40)
        write(2,90,err=60) n,form
  90    format(i5,4x,a20)
        i1=nstart(np)
        i2=nstart(np)+n-1
        write(2,form,err=60) (array(i),i=i1,i2)
        close(2)
        go to 99
  60    write(*,61)
  61    format(' Error in opening or writing to a file')
	imflag=1
  99    continue
        return
        end
c&SPEAKER
        subroutine speaker(args,nargs)
c*******************************************************************
c
c   Subroutine to make the speaker sound in frequency NFREQ for
c   NHUND hundredths of a second.
c
c   SPEAKER(nfreq,nhund)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 icki,ickl
c
c
	if(nargs.eq.3) then
	if(ickl(args(1),npf,nf).eq.1) go to 99
	if(ickl(args(2),nph,nh).eq.1) go to 99
	call ckint(args(3),n)
	if(n.lt.1.or.n.gt.min0(nf,nh)) then
		call error(args,3,2)
		go to 99
	endif
	call movxy(wk,array(nstart(npf)),4*n)
	call movxy(wk(n+1),array(nstart(nph)),4*n)
	do 5 i=1,n
    	c1=wk(i)
	c2=wk(n+i)
	if(c1.lt.21..or.c1.gt.32767..or.c2.lt.1..or.c2.gt.32767.) then
		write(*,6)
   6		format(' Illegal Input')
		imflag=1
		go to 99
	endif
   5	continue
	do 10 i=1,n
  10	call speak(int(wk(i)),int(wk(n+i)))
	go to 99
	endif
c
c
	if(icki(args(1),1,nfreq,1,21).eq.1) go to 99
	if(icki(args(2),2,nhund,1,1).eq.1) go to 99
c
c
        call speak(nfreq,nhund)
c
c
  99    continue
        return
        end
c&TIME
        subroutine time1(args,nargs)
c*******************************************************************
c
c   Subroutine to return the elapsed time since midnight in seconds.
c
c   TIME(t)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15

        t=secnds(0.0)
        call ckaddr(args(1),t,iref)
        return
        end
