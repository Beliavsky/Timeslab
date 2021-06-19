c&TYPE5
c&TYPE5
        subroutine type5(l,iend,lun)
c******************************************************************
c
c   subroutine to process type 5 commands.
c
c******************************************************************
#include 'tslabc'
        character l(72)*1,pname*15,p1name(1)*15
        character others(12)*15,pgstr*29
        character*1 pgc
        integer*4 pgbegin,pgcurs
        integer*4 system
        noth=9
        others(1)='TYPE1'
        others(2)='TYPE2'
        others(3)='TYPE3'
        others(4)='TYPE4'
        others(5)='TYPE5'
        others(6)='SEED'
        others(7)='GOTO'
        others(8)='WHILE'
        others(9)='END'

        iend=0
        call formnm(l,1,lnb,pname)
        p1name(1)=pname
        call ltoup(pname,15)
        call cklist(pname,pr5lst,numpr5,np)
        if(np.ge.1.and.np.le.numpr5) go to 100
c
c   Check if listing by just saying a variable name:
c
        call cklist(p1name(1),xlist,numarr,np1)
        call cklist(p1name(1),ilist,nivars,np2)
        call cklist(p1name(1),rlist,nrvars,np3)
        call cklist(p1name(1),chlist,nchars,np4)
        if(np1.gt.0.or.np2.gt.0.or.np3.gt.0.or.np4.gt.0) then
                nargs=1
                call listts(p1name,nargs)                 
                go to 99
        endif
c
c
        write(*,79)
  79    format(' unrecognizable type 5 command')
        imflag=1
        go to 99
 100    go to(5005,5006,5010,5020,5030,5040,5050,5060,5070,5080,5090,
     1  5100,5110,5120,5130,5140,5150,5160,5170,5180,5190,5200,5210,
     1  5220,5230,5240,5250,5260,5270,5280,5290,5300,5310,5320) np
c
c   quit :
c
5005    iend=1
        go to 99
c
c   help file :
c
5006    iend=1
        go to 99
c       iatth=kback*16+kfore
c       iatt=iatth
c       if(kfore.le.7) iatt=iatth+8
c       if(kfore.ge.8) iatth=iatt-8
c       call scrla(0,0,0,24,79,0,iatth)
c       call blnksc(1)
c       call wrtxtb('Available Topics:',17,0,0,24,iatt)
c       call wrtxtb(' For help about a topic, enter HELP(topic)',42,
c    1  0,17,24,iatt)
c       call wrtxtb('For an outline of what commands do what, enter',
c    1  47,0,0,23,iatt)
c       call wrtxtb(' HELP(outline)',14,0,47,23,iatt)
c       nrow=22
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,0,nrow,iatt)
c       call wrtxtb('Type 1 Commands:',18,0,31,nrow,iatt)
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,49,nrow,
c    1  iatt)
c       call srtlst(pr1lst,numpr1,iwk,iwk(numpr1+1),nrow,iatth)
c       nrow=nrow-1
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,0,nrow,iatt)
c       call wrtxtb('Type 3 Commands:',18,0,31,nrow,iatt)
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,49,nrow,
c    1  iatt)
c       call srtlst(pr3lst,numpr3,iwk,iwk(numpr3+1),nrow,iatth)
c       nrow=nrow-1
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,0,nrow,iatt)
c       call wrtxtb('Type 5 Commands:',18,0,31,nrow,iatt)
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,49,nrow,
c    1  iatt)
c       call srtlst(pr5lst,numpr5,iwk,iwk(numpr5+1),nrow,iatth)
c       nrow=nrow-1
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,0,nrow,iatt)
c       call wrtxtb('Other Topics:   ',18,0,31,nrow,iatt)
c       call wrtxtb('อออออออออออออออออออออออออออออออ',31,0,49,nrow,
c    1  iatt)
c       call srtlst(others,noth,iwk,iwk(noth+1),nrow,iatth)
c       call stcura(0,nrow-1,0)
c       call blnksc(0)
5010    continue
c
c   information about arrays defined so far :
c
        write(*,5016) numarr,num
5016    format(' Number of arrays and Free Elements: ',2i7)
        if(numarr.gt.0) then
        write(*,5015)
5015    format(1x,4hNAME,11x,6hLENGTH,2x,5hLABEL/1x,64(1h-))
        do 5019 i=1,numarr
5019    write(*,5011) xlist(i),nend(i)-nstart(i)+1,label(i)
5011    format(1x,a15,1x,i5,2x,a40)
        write(*,5021)
5021    format(1x,64(1h-))
        endif
        if(nivars.gt.0) write(*,501) (ilist(i),ints(i),i=1,nivars)
 501    format(1x,a15,' = ',i15,1x,a15,' = ',i15)
        if(nrvars.gt.0) write(*,502) (rlist(i),reals(i),i=1,nrvars)
 502    format(1x,a15,' = ',f15.8,1x,a15,' = ',f15.8)
        if(nchars.gt.0) write(*,503) (chlist(i),chars(i),i=1,nchars)
 503    format(1x,a15,' = ',a15,1x,a15,' = ',a15)
        if(iptrec.eq.0) go to 99
c
c   RECORD and INFO:
c
        write(3,5016,err=71) numarr,num
        if(numarr.gt.0) then
        write(3,5015,err=71)
        do 72 i=1,numarr
  72    write(3,5011,err=71) xlist(i),nend(i)-nstart(i)+1,label(i)
        write(3,5021,err=71)
        endif
        if(nivars.gt.0) then
        write(3,501,err=71) (ilist(i),ints(i),i=1,nivars)
        endif
        if(nrvars.gt.0) then
        write(3,502,err=71) (rlist(i),reals(i),i=1,nrvars)
        endif
        if(nchars.gt.0) then
        write(3,503,err=71) (chlist(i),chars(i),i=1,nchars)
        endif
        go to 99
  71    write(*,73)
  73    format(1x,'RECORD file full, I''m closing it')
        close(3)
        iptrec=0
        imflag=1
        go to 99
5020    continue
c
c   clean :
c
        call clean(0)
        go to 99
5030    continue
c
c   new label for some array :
c
        write(*,5032)
5032    format(' LABEL With no argument is no longer allowed')
        go to 99
5040    continue
c
c   overwrite on :
c
        ioptov=1
        go to 99
5050    continue
c
c   overwrite off :
c
        ioptov=0
        go to 99
5060    continue

c   batch operation on :
c
        ioptb=1
        ioptv1=ioptov
        ioptov=1
        go to 99
5070    continue
c
c   batch operation off :
c
        ioptb=0
        ioptov=ioptv1
        go to 99
5080    continue
c
c   plot screen on :
c
        ioptps=1
        go to 99
5090    continue
c
c   plot screen off :
c
        ioptps=0
        go to 99
5100    continue
c
c   speaker on :
c
        ioptsp=1
        go to 99
5110    continue
c
c   speaker off :
c
        ioptsp=0
        go to 99
5120    continue
c
c   echo :
c
        open(1,file='lpt1',err=5125)
        close(1)
        call echo
        go to 99
5125    write(*,5126)
5126    format(' printer not attached')
        imflag=1
        go to 99
5130    continue
c
c   dos :
c
        call dos(lun)
        go to 99
5140    continue
c
c   PROMPTON :
c
        ioptpr=1
        go to 99
5150    continue
c
c   promptoff :
c
        ioptpr=0
        go to 99
5160    continue
c
c   pause :
c
        write(*,5165)
5165    format(' PAUSE...STRIKE q TO BREAK, ANYTHING ELSE TO GO ON',$)
        call waitk(n10)
        if(n10.eq.113) ioptf3=1
        go to 99
5170    continue
c
c  page :
c
        open(1,file='lpt1',err=5175)
        write(1,*) char(12)
        close(1)
        go to 99
5175    write(*,5126)
        imflag=1
        go to 99
5180    continue
c
c   macro:
c
c       macopt=1
c       maccnt=max0(maccnt-1,0)
        go to 99
5190    continue
c
c   restart :
c
        iend=2
        go to 99
5200    continue
c
c   singleon:
c
        ioptsg=1
        go to 99
5210    continue
c
c   singleoff
c
        ioptsg=0
        go to 99
5220    continue
c
c   CLS:
c
        i=system('clear')
        go to 99
5230    continue
c
c   ABORTON:
c
        ioptf2=1
        go to 99
5240    continue
c
c   ABORTOFF:
c
        ioptf2=0
        go to 99
5250    continue
c
c   help1:
c
        if(iptpon.eq.1) then
                write(*,5251)
5251            format(' HELP only available in the text mode')
                imflag=1
                go to 99
        endif
        iatth=kback*16+kfore
        iatt=iatth
        if(kfore.le.7) iatt=iatth+8
        if(kfore.ge.8) iatth=iatt-8
        call srt1(pr1lst,pr3lst,pr5lst,others,numpr1,numpr3,numpr5,
     1  noth,wk,iwk,iwk(1001),iatth)
        go to 99
5260    continue
c 
c   PLOTON:
c
        if(pgbegin(0,'?',1,1).ne.1) stop
        go to 99
5270    continue
c
c   PLOTOFF:
c
        iptpon=0
        call pgclos
        go to 99
5280    continue
c
c   ENDIF:
c
        if(lun.eq.100) then
                write(*,5281)
5281            format(' ENDIF only allowed in macros')
        endif
        go to 99
 5290   continue
c
c   XON:
c
        if(pgbegin(0,'/xwin',1,1).ne.1) go to 99
        call pgask(ask)
        call pgpap(6.0,1.0)
        go to 99
 5300   continue
c
c   FIND:
c
        pgx=(pgxmin+pgxmax)/2.
        pgy=(pgymin+pgymax)/2.
        pgry=pgymax+.10*(pgymax-pgymin)
c        write(*,*) pgx,pgy
        write(pgstr,5301) pgx,pgy
 5301   format('x = ',f10.4,' y = ',f10.4)
        call pgptxt(pgxmax,pgry,0.0,1.0,pgstr)
 5305   i=pgcurs(pgx,pgy,pgc)
        call pgsci(0)
        call pgptxt(pgxmax,pgry,0.0,1.0,pgstr)
        write(pgstr,5301) pgx,pgy
        call pgsci(1)
        call pgptxt(pgxmax,pgry,0.0,1.0,pgstr)
c        write(*,*) i,pgx,pgy,pgc
        if((i.eq.0).or.(pgc.eq.'q')) go to 99
        go to 5305
 5310   continue
c
c  ASKON:
c
        ask=.true.
        call pgask(ask)
        go to 99
 5320   continue
c
c   ASKOFF:
c
        ask=.false.
        call pgask(ask)
        go to 99
c
c   finish :
c
  99    continue
        return
        end
c&SRTLST
        subroutine srtlst(list,n,il,il1,nrow,iatth)
c*******************************************************************
c
c   Subroutine to display a sorted list of the elements of the
c   character*15 array list(1),...,list(n).
c
c******************************************************************
c
        character list(n)*15,l89(100)*10
        dimension il(n)

c   get distinct items:

        nn=1
        list(n+1)=list(1)
        do 10 i=2,n
           do 5 j=1,nn
 5         if(list(i).eq.list(n+j)) go to 10
           nn=nn+1
           list(n+nn)=list(i)
 10     continue
        do 15 i=1,nn
 15     list(i)=list(n+i)
        n=nn

c   get il(i)=position of item i in sorted list:

        call sortl(list,n,il)

c   use this to get sorted list:

        do 20 i=1,n
 20     list(n+il(i))=list(i)
        do 25 i=1,n
 25     list(i)=list(n+i)

c   write it out:

c        write(*,35) (list(j),j=1,n)
c        write(*,36)
c 36     format(/)
c        write(*,35) list(1),list(2)

        nperl=8
        nlines=((n-1)/nperl)+1
        do 30 i=1,nlines
           k=0
           do 31 j=i,n,nlines
           k=k+1
 31        l89(k)=list(j)
 30     write(*,35) (l89(j),j=1,k)
 35     format(1x,8a10)

        return
        end
c&SRT1
        subroutine srt1(l1,l3,l5,oth,n1,n3,n5,no,lall,il,il1,iatth)
        character l1(n1)*15,l3(n3)*15,l5(n5)*15,oth(no)*15,lall(300)*15
        dimension il(1),il1(1)

        n=n1+n3+n5+no

        call movxy(lall,l1,15*n1)
        call movxy(lall(n1+1),l3,15*n3)
        call movxy(lall(n1+n3+1),l5,15*n5)
        call movxy(lall(n1+n3+n5+1),oth,15*no)

        write(*,10)
 10     format(/,' List of command names:',/)
        call srtlst(lall,n,il,il1,nrow,iatth)

        write(*,20)
 20     format(/,1x,'For help on a command, type HELP(name)')

        write(*,30)
 30     format(1x,
     1     'For a list of general help areas,type HELP(areas)',/)

        return
        end
      
      subroutine sortl(names,n,il)
c*********************************************************************
      dimension il(1)

      character*15 names(n),namei,namej
      character*1 namei1(15),namej1(15)

      equivalence (namei,namei1), (namej,namej1)

      do 10 i=1,n
         numlti=0
         namei=names(i)
         do 20 j=1,n
            namej=names(j)
            do 30 k=1,15
            if(ichar(namej1(k)).eq.ichar(namei1(k))) go to 30
            if(ichar(namej1(k)).gt.ichar(namei1(k))) go to 20
            if(ichar(namej1(k)).lt.ichar(namei1(k))) then
               numlti=numlti+1
               go to 20
            endif
 30         continue
 20      continue
 10   il(i)=numlti+1

      return
      end
