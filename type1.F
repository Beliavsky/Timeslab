c&TYPE1
c&TYPE1
        subroutine type1(l,lun)
#include 'tslabc'
        character l(72)*1,pname*15,args(21)*15,vname*15
        character args1(3)*15,i10*1
        character ck*15
        character lno*72,lno1(72)*1,l15*15,l151*15,l152*15
        character mpt*15,mpt1(15)*1
        integer*2 ickl,ickr
        integer pgbegin
        logical lgc
        equivalence (lno,lno1)
        equivalence (mpt,mpt1)
        dimension icnum(20)
c
c   parse command line, convert pname to upper case, and go to
c   appropriate place :
c
        call parse(l,pname,args,nargs,vname,ierr)
        if(ierr.eq.1) then
                write(*,6701)
6701    format(' An argument or command name is too long or too short')
                imflag=1
                go to 99
        endif
        call ltoup(pname,15)
        call cklist(pname,pr1lst,numpr1,np)
        if(np.gt.0.and.np.le.numpr1) go to 1002
        write(*,1001)
1001    format(' undefined type 1 procedure ')
        imflag=1
        go to 99
1002    continue
        go to(1005,1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,
     1  1110,1120,1130,1140,1150,1160,1170,1180,1190,1200,1210,1220,
     1  1230,1240,1250,1260,1270,1280,1290,1300,1310,1320,1330,1340,
     1  1350,1360,1370,1380,1390,1400,1410,1420,1430,1440,1450,
     1  1460,1470,1480,1490,1500,1510,1520,1530,1540,1550,1560,
     1  1570,1580,1590,1600,1610,1620,1630,1640,1650) np
c
c   procedure READ :
c
1005    continue
        if(nargs.ne.3.and.nargs.ne.4) go to 98
        call readfl(args,nargs)
        go to 99
c
c   procedure LIST :
c
1010    continue
        if(nargs.lt.1.or.nargs.gt.21) go to 98
        call listts(args,nargs)
        go to 99
1020    continue
c
c   procedure PLOT :
c
        iptt=0
        if(nargs.eq.1) then
                if(ickl(args(1),npx,nx).eq.1) go to 99
                write(args(2),1021) nx
1021            format(i15)
                call rmblnk(args(2),15,npx)
                nargs=2
        endif
        if(nargs.eq.6) then
                if(ickl(args(1),npx,nx).eq.1) go to 99
                iptt=1
                do 1022 i=7,2,-1
1022            args(i)=args(i-1)
                args(1)='t&&&&'
                do 1023 i=1,nx
1023            wk(i)=i
                lab='     '
                call ckadda(args(1),nx,lab,1,iref)
                if(iref.eq.1) go to 99
                nargs=7
        endif
        call plotts(args,nargs)
        if(iptt.eq.1) then
                call clean(numarr)
        endif
        go to 99
1030    continue
c
c   procedure SAVE :
c
        if(nargs.ne.3) go to 98
        call save(args,nargs)
        go to 99
1040    continue
c
c   procedure PLOTSP :
c
        ierd=0
        ngs=nargs/3
        if(3*ngs.ne.nargs) go to 98
        call loadpg(0.,0.5,-6.,6.)
        call pgenv(0,.5,-6.,6.,0,0)
        if(iptpon.eq.1.and.nplt.gt.0) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973    continue
        do 1046 i=1,ngs
        ii=3*(i-1)
        args1(1)=args(ii+1)
        args1(2)=args(ii+2)
        args1(3)=args(ii+3)
        call plotsp(args1,3,i,ierr)
        if(iptpon.eq.1) nplt=nplt+1
        if(ierr.eq.1) go to 99
1046    continue
c        call pltend(args,nargs,0.,0.5,-6.,6.)
c       if(ioptps.eq.1) call prtscn
c       if(ioptf1.eq.1) then
c               call rdscr(charsv,wk,ierd)
c               ioptf1=0
c               endif
c1047    if(ioptb.eq.1) go to 1048
c       icol=500
c       irow=150
c       call f5f9(icol,irow,5,iwk,0.0,0.5,-6.0,6.0,
c    1  npixx,npixy,nx1,ny1,wk,1,
c    1  kolor,ioptf3)
c1048    continue
c          if(iptpon.eq.1) then
c               call drrect(0,0,640,8,0,1)
c               call stcura(0,24,0)
c               go to 99
c       endif
c       call smodea(3)
c       call txtcls
c       if(ierd.eq.1) call error(args,1,5)
        go to 99
c
c
1050    continue
c
c   procedure HELP :
c
        if(nargs.ne.1.and.nargs.ne.2) go to 98
        if(iptpon.eq.1) then
                write(*,1051)
1051            format(' HELP only available in the text mode')
                imflag=1
                go to 99
        endif
        call help(args,nargs,kback,kfore)
        go to 99
1060    continue
c
c   procedure MACRO :
c
        if(nargs.ne.2.and.nargs.ne.1) go to 98
        call macro(args,nargs,lun)
        maccnt=0
        go to 99
1070    continue
c
c   procedure PRINT :
c
        if(nargs.lt.1.or.nargs.gt.21) go to 98
        mpt=args(nargs)
        if(mpt1(1).eq.char(39)) then
                write(*,1071)
1071            format(' No file name allowed in PRINT command')
                imflag=1
                go to 99
        endif
        nargs=nargs+1
        args(nargs)='''prn'''
        call listts(args,nargs)
        go to 99
1080    continue
c
c   procedure SPEAKER :
c
        if(nargs.ne.2.and.nargs.ne.3) go to 98
        call speaker(args,nargs)
        go to 99
1090    continue
c
c   procedure EDIT :
c

        go to 99
1100    continue
c
c   procedure PLOTCSP :
c
        if(nargs.ne.2.and.nargs.ne.4) go to 98
        call pltcsp(args,nargs)
        go to 99
c
c   procedure LISTSP :
c
1110    continue
        if(nargs.ne.2) go to 98
        call listsp(args,nargs)
        go to 99
1120    continue
c
c   procedure infqnt :
c
        if(nargs.ne.1.and.nargs.ne.2.and.nargs.ne.4.and.nargs.ne.5)
     1  go to 98
        call infqnt(args,nargs)
        go to 99
1130    continue
c
c   procedure IF :
c
        if(nargs.ne.4.and.nargs.ne.2.and.nargs.ne.1) go to 98
        call ifmac(args,nargs,lun)
        go to 99
1140    continue
c
c   procedure HIST :
c
        call hist(args,nargs)
        go to 99
1150    continue
c
c   procedure CLEAN(var) :
c
        do 1159 i=1,nargs
        ck=args(i)
        call ltoup(ck,15)
        if(ck.eq.'ALL') then
                nrvars=0
                nivars=0
                nchars=0
                numarr=0
                num=maxnum
                go to 1159
        endif
        if(ck.eq.'REALS') then
                nrvars=0
                go to 1159
        endif
        if(ck.eq.'INTS') then
                nivars=0
                go to 1159
        endif
        if(ck.eq.'CHARS') then
                nchars=0
                go to 1159
        endif
        if(ck.eq.'ARRAYS') then
                numarr=0
                num=maxnum
                go to 1159
        endif
        call clean1(args(i),ier)
        if(ier.eq.1) then
                call error(args,i,2)
                go to 99
        endif
1159    continue
        go to 99
1160    continue
c
c   procedure LENGTH :
c
        call cklist(args(1),xlist,numarr,np)
        if(np.eq.-999) then
                call error(args,1,1)
                go to 99
                endif
        write(*,1166) args(1),nend(np)-nstart(np)+1
1166    format(1x,a15,' has length ',i6)
        if(nargs.eq.2) call ckaddi(args(2),nend(np)-nstart(np)+1,iref)
        go to 99
1170    continue
c
c   procedure PLOT2 :
c
        if(nargs.eq.6.or.nargs.eq.7) go to 1175
        if(nargs.eq.10.or.nargs.eq.11) go to 1175
        go to 98
1175    continue
        call plot2(args,nargs)
        go to 99
1180    continue
c
c   procedure TIME :
c
        if(nargs.ne.1) go to 98
        call time1(args,nargs)
        go to 99
1190    continue
c
c   procedure DELAY :
c
        if(nargs.ne.1) go to 98
        call delay(args,nargs)
        go to 99
1200    continue
c
c   procedure COLOR :
c
        if(nargs.ne.1) go to 98
        call ckint(args(1),nk)
        if(nk.eq.-999) then
                call error(args,1,2)
                go to 99
                endif
        kolor=nk
        if(iptpon.eq.1) then
                call colora(0)
                call colora(kolor)
        endif
        go to 99
1210    continue
c
c   procedure CORR2 :
c
        if(nargs.ne.12) go to 98
        call corr2(args,nargs)
        go to 99
1220    continue
c
c   procedure CROSSP :
c
        if(nargs.ne.10) go to 98
        call crossp(args,nargs)
        go to 99
1230    continue
c
c   procedure POLAR :
c
        if(nargs.ne.5.and.nargs.ne.4) go to 98
        call polar(args,nargs)
        go to 99
1240    continue
c
c   procedure REG :
c
        if(nargs.ne.11.and.nargs.ne.12.and.nargs.ne.8.and.nargs.ne.10)
     1  go to 98
        call reg(args,nargs)
        go to 99
1250    continue
c
c   procedure DTFORE :
c
        if(nargs.ne.10) go to 98
        call dtfore(args,nargs)
        go to 99
1260    continue
c
c   procedure TEXTCOLOR :
c
        if(nargs.ne.2) go to 98
        call ckint(args(1),iback)
        if(iback.lt.0.or.iback.gt.7) then
                call error(args,1,2)
                go to 99
        endif
        call ckint(args(2),ifore)
        if(ifore.lt.0.or.ifore.gt.15.or.ifore.eq.iback) then
                call error(args,2,2)
                go to 99
        endif
        kback=iback
        kfore=ifore
        if(ifore.gt.7) ifore=ifore-8
        if(ifore.eq.kfore) write(*,*) char(27),char(91),char(48),
     1  char(109)
        if(ifore.ne.kfore) write(*,*) char(27),char(91),char(49),
     1  char(109)
        call txtcol(iback,ifore,1)
        go to 99
1270    continue
c
c   procedure RESCREEN:
c
        call rescn(args,nargs)
        go to 99
1280    continue
c
c   procedure MAXMIN:
c
        if(nargs.ne.6) go to 98
        call maxmin(args,nargs)
        go to 99
1290    continue
c
c   procedure ARSP2:
c
        if(nargs.ne.8) go to 98
        call arsp2(args,nargs)
        go to 99
1300    continue
c
c   procedure SAVESC:
c
        if(nargs.ne.1) go to 98
        call cklist(args(1),chlist,nchars,n1)
        charsv=args(1)
        if(n1.ne.-999) charsv=chars(n1)
        if(iptpon.eq.1) then
                call drrect(0,192,640,8,0,1)
                call rdscr(charsv,wk,ierd)
                go to 99
        endif
        ioptf1=1
        go to 99
1310    continue
c
c   procedure PRINTER:
c
        call ckint(args(1),ipt)
        if(ipt.eq.27) then
                do 1316 i=2,nargs
1316            call ckint(args(i),icnum(i-1))
                open(1,file='lpt1',err=1322)
                write(1,*) (char(icnum(i)),i=1,nargs-1)
                close(1)
                go to 99
        endif
        if(ioptls.ne.1.and.ioptls.ne.11.and.ioptls.ne.21) then
                write(*,1317)
1317    format(' iopt1 must be 27 if not using IBM graphics printer')
                imflag=1
                go to 99
                endif
        do 1319 i=1,nargs
        call ckint(args(i),ipt)
        if(ipt.eq.-999) then
                call error(args,i,2)
                go to 99
                endif
        if(ipt.ge.1.and.ipt.le.10) then
        open(1,file='lpt1',err=1322)
        if(ipt.eq.1) write(1,*) char(15),char(27),char(48)
        if(ipt.eq.2) write(1,*) char(18),char(27),char(51),char(36)
        if(ipt.eq.3) write(1,*) char(27),char(55)
        if(ipt.eq.4) write(1,*) char(27),char(54)
        if(ipt.eq.5) write(1,*) char(27),char(69)
        if(ipt.eq.6) write(1,*) char(27),char(70)
        if(ipt.eq.7) write(1,*) char(27),char(71)
        if(ipt.eq.8) write(1,*) char(27),char(72)
        if(ipt.eq.9) write(1,*) char(27),char(87),char(1)
        if(ipt.eq.10) write(1,*) char(27),char(87),char(0)
        close(1)
        go to 1319
        endif
        if(ipt.gt.30) then
                open(1,file='lpt1',err=1322)
                do 1315 j=31,ipt
1315            write(1,*) char(10)
                close(1)
                go to 1319
                endif
        call error(args,1,2)
1319    continue
        go to 99
1322    write(*,1323)
1323    format(' error in opening printer')
        imflag=1
        go to 99
1320    continue
c
c   procedure LISTM:
c
        if(nargs.ne.3.and.nargs.ne.4) go to 98
        call listm(args,nargs)
        go to 99
1330    continue
c
c   procedure NOTES:
c
        if(nargs.ne.1.and.nargs.ne.2) go to 98
        inquire(file=args(1),exist=lgc)
                if(.not.lgc) then
                        write(*,1332) args(1)
1332                    format(1x,a15,'Doesn''t exist')
                        imflag=1
                        go to 99
                        endif
                open(2,file=args(1))
                if(nargs.eq.1) then
                write(*,1331) args(1)
1331            format(/,' You have the following topics in ',a15)
1335            read(2,1336,end=99) lno
1336            format(a72)
                if(lno1(1).ne.' ') write(*,1337) (lno1(j),j=1,15)
1337            format(1x,15a1)
                go to 1335
                endif
c
c
                l15=args(2)
                call cklist(args(2),chlist,nchars,np)
                if(np.ne.-999) l15=chars(np)
                l151=l15
                call ltoup(l151,15)
                no=0
1338            read(2,1336,end=1333) lno
                        if(lno1(1).ne.' ') then
                        if(no.eq.1) then
                                close(2)
                                go to 99
                                endif
                        call movxy(l152,lno,15)
                        if(l152.ne.l15.and.l152.ne.l151) go to 1338
                        no=1
                        write(*,1339) lno
1339                    format(1x,a72)
                        go to 1338
                        endif
                        if(no.eq.0) go to 1338
                        write(*,1339) lno
                        go to 1338
1333            if(no.eq.0) then
                        write(*,1334) l15,l151,args(1)
1334                    format(1x,a15,',',a15,' not in ',a15)
                        endif
                        close(2)
                        go to 99
1340    continue
c
c   procedure ARMASEL:
c
        if(nargs.ne.14) go to 98
        call armasl(args,nargs)
                go to 99
1350    continue
c
c   procedure COEFFCSD:
c
        if(nargs.ne.4.and.nargs.ne.8) go to 98
        call armasd(args,nargs)
        go to 99
1360    continue
c
c   procedure FFT:
c
        if(nargs.ne.7.and.nargs.ne.6) go to 98
        call ffts(args,nargs)
        go to 99
1370    continue
c
c   procedure DOUBLE:
c
        if(nargs.ne.3) go to 98
        call double(args,nargs)
        go to 99
1380    continue
c
c   procedure DTARMA:
c
        if(nargs.ne.12) go to 98
        call dtarma(args,nargs)
        go to 99
1390    continue
c
c   procedure PRINTSEL:
c
        if(nargs.ne.1) go to 98
        call ckint(args(1),ipt)
        if(ipt.eq.2) ipt=22 
        if(ipt.eq.1.or.ipt.eq.3.or.(ipt.ge.21.and.ipt.le.24)) then
                if(ipt.gt.3) then
                open(1,file='LPT1')
                if(ipt.eq.24) write(1,*) char(27),char(42),char(116),
     1                  char(55),char(53),char(82)
                if(ipt.eq.23) write(1,*) char(27),char(42),char(116),
     1                  char(49),char(48),char(48),char(82)
                if(ipt.eq.22) write(1,*) char(27),char(42),char(116),
     1                  char(49),char(53),char(48),char(82)
                if(ipt.eq.21) write(1,*) char(27),char(42),char(116),
     1                  char(51),char(48),char(48),char(82)
                close(1)
                endif
                ioptls=ipt
                go to 99
        endif
        call error(args,1,2)
        imflag=1
        go to 99
1400    continue
c
c   procedure POLYROOTS:
c
        if(nargs.ne.7) go to 98
        call polyrt(args,nargs)
        go to 99
1410    continue
c
c   procedure ROOTSPOLY:
c
        if(nargs.ne.4) go to 98
        call rtpoly(args,nargs)
        go to 99
1420    continue
c
c   procedure ARSPCB:
c
        if(nargs.ne.8) go to 98
        call arspcb(args,nargs)
        go to 99
1430    continue
c
c   procedure ARSPPEAK:
c
        if(nargs.ne.7.and.nargs.ne.8) go to 98
        call arpeak(args,nargs)
        go to 99
1440    continue
c
c   procedure ARMAPRED:
c
        if(nargs.ne.13.and.nargs.ne.14) go to 98
        call mxpred(args,nargs)
        go to 99
1450    continue
c
c   procedure DENSITY:
c
        if(nargs.ne.7) go to 98
        call dnsity(args,nargs)
        go to 99
1460    continue
c
c   procedure MCHOL:
c
        if(nargs.ne.5) go to 98
        call mcd(args,nargs)
        go to 99
1470    continue
c
c   procedure GS:
c
        if(nargs.ne.6) go to 98
        call gs(args,nargs)
        go to 99
1480    continue
c
c   procedure PLOTK:
c
        if(nargs.ne.5.and.nargs.ne.9.and.nargs.ne.12) go to 98
        call plotk(args,nargs)
        go to 99
1490    continue
c
c   procedure SEASEST:
c
        if(nargs.ne.11.and.nargs.ne.12) go to 98
        call seases(args,nargs)
        go to 99
1500    continue
c
c   procedure SEASPRED:
c
        if(nargs.ne.14) go to 98
        call seaspr(args,nargs)
        go to 99
1510    continue
c
c   procedure RECORD:
c
        if(nargs.ne.1) go to 98
        if(args(1).eq.'CLOSE'.or.args(1).eq.'close') then
                iptrec=0
                close(3)
                go to 99
        endif
        if(iptrec.eq.1) then
                write(*,1491)
1491            format(' A RECORD file is already open')
                write(3,1491)
                imflag=1
                go to 99
        endif
        iptrec=1
        args(2)=args(1)
        call ltoup(args(2),15)
        if(args(2).eq.'LPT1'.or.args(2).eq.'PRN') then
                open(3,file='lpt1')
                go to 99
        endif
        call cklist(args(1),chlist,nchars,np)
        if(np.ne.-999) args(1)=chlist(np)
        inquire(file=args(1),exist=lgc)
        if(lgc) go to 1492
        open(3,file=args(1),status='new')
        go to 99
1492    open(3,file=args(1),status='old')
1493    read(3,1494,end=1495) i10
1494    format(a1)
        go to 1493
1495    backspace 3
        go to 99
1520    continue
c
c   procedure ARCORR2:
c
        if(nargs.ne.12) go to 98
        call arcr2(args,nargs)
        go to 99
1530    continue
c
c   procedure PLOTSIZE:
c
        if(nargs.eq.1) then
                call ckint(args(1),iptsz)
                if(iptsz.lt.0.or.iptsz.gt.8) then
                        call error(args,1,2)
                        go to 99
                endif
                if(iptsz.eq.0) then
                        npixx=480
                        npixy=120
                        nx1=55
                        ny1=20
                        go to 99
                endif
                if(iptsz.ge.1.and.iptsz.le.4) then
                        npixx=260
                        npixy=60
                        nx1=50
                        ny1=20
                        if(iptsz.eq.1.or.iptsz.eq.3) ny1=120
                        if(iptsz.eq.3.or.iptsz.eq.4) nx1=360
                        go to 99
                endif
                if(iptsz.eq.5.or.iptsz.eq.6) then
                        npixx=260
                        npixy=120
                        nx1=50
                        ny1=20
                        if(iptsz.eq.6) nx1=360
                        go to 99
                endif
                if(iptsz.eq.7.or.iptsz.eq.8) then
                        npixx=480
                        npixy=60
                        nx1=50
                        ny1=110
                        if(iptsz.eq.8) ny1=20
                        go to 99
                endif
        endif
c
c
        if(icki(args(1),1,npixx,1,20).eq.1) go to 99
        if(icki(args(2),2,npixy,1,20).eq.1) go to 99
        if(icki(args(3),3,nx1,1,20).eq.1) go to 99
        if(icki(args(4),4,ny1,1,20).eq.1) go to 99
        if(nargs.eq.4) go to 99
        if(icki(args(5),5,nticsx,2,2).eq.1) go to 99
        if(icki(args(6),6,nticsy,2,2).eq.1) go to 99
        if(nargs.eq.6) go to 99
        if(icki(args(7),7,nclabx,1,0).eq.1) go to 99
        if(icki(args(9),9,nclaby,1,0).eq.1) go to 99
        if(icki(args(8),8,ndecx,1,0).eq.1) go to 99
        if(icki(args(10),10,ndecy,1,0).eq.1) go to 99
        ipte=0
        if(nclabx.gt.0.and.ndecx.ge.nclabx) ipte=8
        if(nclaby.gt.0.and.ndecy.ge.nclaby) ipte=9
        if(ipte.ne.0) then
                call error(args,ipte,2)
                go to 99
        endif
        go to 99
1540    continue
c
c   procedure ERASE:
c
        if(iptpon.eq.0) go to 99
        if(nargs.ne.1.and.nargs.ne.4) go to 99
        if(nargs.eq.1) then
                call ckint(args(1),iptsz)
                if(iptsz.lt.0.or.iptsz.gt.8) go to 99
                if(iptsz.eq.0) call drrect(0,0,640,200,0,1)
                if(iptsz.eq.1) call drrect(0,100,320,100,0,1)
                if(iptsz.eq.2) call drrect(0,0,320,100,0,1)
                if(iptsz.eq.3) call drrect(320,100,320,100,0,1)
                if(iptsz.eq.4) call drrect(320,0,320,100,0,1)
                if(iptsz.eq.5) call drrect(0,0,320,200,0,1)
                if(iptsz.eq.6) call drrect(320,0,320,200,0,1)
                if(iptsz.eq.7) call drrect(0,100,640,100,0,1)
                if(iptsz.eq.8) call drrect(0,0,640,100,0,1)
                go to 99
        endif
        call ckint(args(1),n1)
        if(n1.lt.0.or.n1.gt.639) go to 99
        call ckint(args(2),n2)
        if(n2.lt.0.or.n2.gt.199) go to 99
        call ckint(args(3),n3)
        if(n3.lt.0.or.n1+n3.gt.640) go to 99
        call ckint(args(4),n4)
        if(n4.lt.0.or.n2+n4.gt.200) go to 99
        call drrect(n1,n2,n3,n4,0,1)
        go to 99
1550    continue
c
c   procedure GRMENU:
c
        if(nargs.ne.4) go to 98
        if(iptpon.ne.1) go to 99
        if(ickr(args(1),1,xmin,0,0.).eq.1) go to 99
        if(ickr(args(2),2,xmax,2,xmin).eq.1) go to 99
        if(ickr(args(3),3,ymin,0,0.).eq.1) go to 99
        if(ickr(args(4),4,ymax,2,ymin).eq.1) go to 99
        icol=500
        irow=150
        call f5f9(icol,irow,5,iwk,xmin,xmax,ymin,ymax,
     1  npixx,npixy,nx1,ny1,wk,1,kolor,ioptf3)
        call drrect(0,0,640,8,0,1)
        go to 99
1560    continue
c
c   procedure PLOT3D:
c
        if(nargs.ne.6) go to 98
        call plt3d(args,nargs)
        go to 99
1570    continue
c
c   procedure CONTOUR:
c
        if(nargs.ne.7) go to 98
        call cont(args,nargs)
        go to 99
 1580   continue
c
c   procedure DEVON:
c
        if(nargs.ne.1.and.nargs.ne.3) go to 98
        call cklist(args(1),chlist,nchars,np1)
        if(np1.ne.-999) args(1)=chars(np1)
        nrows=1
        ncols=1
        if(nargs.eq.3) then
           if(icki(args(2),2,nrows,5,10).eq.1) go to 99
           if(nrows.lt.1) then
              call error(args,2,2)
              go to 99
           endif
           if(icki(args(3),3,ncols,5,10).eq.1) go to 99
           if(ncols.lt.-10) then
              call error(args,3,3)
              go to 99
           endif
        endif
        write(*,*) nrows,ncols
        write(*,*) args(1)
        if(pgbegin(0,args(1),ncols,nrows).ne.1) go to 99
        call pgask(ask)
        call pgsch(float(nrows)/2.)
        go to 99
 1590   continue
c
c   procedure POINTS:
c
        if(nargs.ne.2.and.nargs.ne.3) go to 98
        if(ickl(args(1),npx,nx).eq.1) go to 99        
        if(ickl(args(2),npy,ny).eq.1) go to 99        
        n=min(nx,ny)
        if(nargs.eq.3) then
           if(icki(args(3),3,nn,5,n).eq.1) go to 99
           if(nn.lt.1) then
              call error(args,3,3)
              go to 99
           endif
           n=nn
        endif
        call pgpt(n,array(nstart(npx)),array(nstart(npy)),-3)
        go to 99
 1600   continue
c
c   procedure LINES:
c
        if(nargs.ne.2.and.nargs.ne.3) go to 98
        if(ickl(args(1),npx,nx).eq.1) go to 99        
        if(ickl(args(2),npy,ny).eq.1) go to 99        
        n=min(nx,ny)
        if(nargs.eq.3) then
           if(icki(args(3),3,nn,5,n).eq.1) go to 99
           if(nn.lt.1) then
              call error(args,3,3)
              go to 99
           endif
           n=nn
        endif
        call pgline(n,array(nstart(npx)),array(nstart(npy)))
        go to 99
 1610   continue
c
c   procedure AXES:
c
        if(nargs.ne.4) go to 98
        if(ickr(args(1),1,xmin,0,0.).eq.1) go to 99
        if(ickr(args(2),2,xmax,2,xmin).eq.1) go to 99
        if(ickr(args(3),3,ymin,0,0.).eq.1) go to 99
        if(ickr(args(4),4,ymax,2,ymin).eq.1) go to 99

        call loadpg(xmin,xmax,ymin,ymax)
        call pgenv(xmin,xmax,ymin,ymax,0,0)

        go to 99
 1620   continue
c
c   procedure LABELS:
c
        if(nargs.ne.3) go to 98
        if(ickl(args(1),npx,nx).eq.1) go to 99
        if(ickl(args(2),npy,ny).eq.1) go to 99
        if(ickl(args(3),npz,nz).eq.1) go to 99
        call pglab(label(npx),label(npy),label(npz))
        go to 99
 1630   continue
c
c  XON(n,m)
c
        if(nargs.ne.2) go to 98

        nrows=1
        ncols=1

        if(icki(args(1),1,nrows,5,10).eq.1) go to 99
        if(nrows.lt.1) then
           call error(args,2,2)
           go to 99
        endif

        if(icki(args(2),2,ncols,5,10).eq.1) go to 99
        if(ncols.lt.-10) then
           call error(args,3,3)
           go to 99
        endif

        if(pgbegin(0,'/xwin',ncols,nrows).ne.1) go to 99
        call pgask(ask)
c        call pgsch(float(nrows)/2.)
        go to 99
 1640   continue
c
c   FONT:
c
        if(ickr(args(1),1,ff,2,0.0).eq.1) go to 99
        call pgsch(ff)
        go to 99
 1650   continue
c
c   PANEL:
c
        if(nargs.ne.2) go to 98
        if(icki(args(1),1,npansx,5,10).eq.1) go to 99
        if(icki(args(2),2,npansy,5,10).eq.1) go to 99
        if(npansx.lt.1) then
           call error(args,1,1)
           go to 99
        endif
        if(npansy.lt.1) then
           call error(args,2,1)
           go to 99
        endif
        call pgpanl(npansx,npansy)
        go to 99
c
c   finish :
c
  98    write(*,101) pr1lst(np)
 101    format(' wrong number of arguments in procedure ',a15)
        imflag=1
  99    continue
        return
        end
