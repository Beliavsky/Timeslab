c&TIMESLAB
c$nofloatcalls
c$storage:2
c        interface to integer*2 function system [c]
c     1      (string[reference])
c        character*1 string
c        end
c&TIMESLAB
#include 'tslabc'
        character l(72)*1,args(21)*15,pname*15,fname*20,vname*15
        character vn1*15,vn2*15,vn3*15,lck3*3
        character i10*1,lck*6,lname*15,lck5*5
        character l72*72,lchr*1,l80(80)*1,l801(80)*1,v6*15,v7*20
        character v8(20)*1
        character ded10*10
c%%
        character ltmp1(72)*1
c%%
c######
        common /i/old
        character*1 old(72,20)
c######
        integer*4 system,status
        integer*4 time,pgbegin
        dimension noc(5),ncc(5)
        equivalence (l,l72),(v7,v8)
        logical logic
c
c******************************************************************
c
c   INITIALIZE:
c
c******************************************************************
c
c   lcwrqq is intrinsic function masking numerical errors from DOS
c
c   1    call lcwrqq(4927)
c        if(pgbegin(0,'?',1,1).ne.1) stop
       
 1      continue
c
c    *******************
c    *                 *
c    *   GRAPHICS      *
c    *                 *
c    *******************
c
c   kolor: graphics foreground color
c   npixx,npixy: length of x and y axis
c   (nx1,ny1): location of origin
c   nticsx,nticsy: number of tic marks on x and y axes
c   nclabx,ndecx: numeric labels on x axis tic marks in F nclabx.ndecx
c   nclaby,ndecy: numeric labels on y axis tic marks in F nclaby.ndecy
c   kolor: foreground color in graphics mode
c
        ask=.false.
        npixx=480
        npixy=120
        nx1=55
        ny1=30
        nticsx=10
        nticsy=10
        nclabx=8
        ndecx=2
        nclaby=6
        ndecy=2
        kolor=6
c
c   lun : logical unit number for macro being executed (100  if none)
c   lun1: stores lun when a  macro is interrupted
c
        lun=100
        lun1=100
c
c   these are number of arrays, integers, reals, and characters:
c
c   num is number of free elements
c
        numarr=0
        nivars=0
        nrvars=0
        nchars=0
        num=maxnum
c
c   zero out the array vector:
c
        call movct(array,4*num,char(0))
c
c   **********************************
c   *                                *
c   *   FLAGS and SWITCHES           *
c   *                                *
c   **********************************
c
c   ioptsp : speakeroff/speakeron (0/1)
c   ioptpr : promptoff/prompton (0/1)
c   ioptov : overoff/overon (0/1)
c   ioptov1: keeps track of ioptov during batchon
c   ioptb  : batchoff/batchon (0/1)
c   ioptps : psoff/pson (0/1)
c   ioptls : which screen dump routine (1: IBM graphics, 2:Laserjet,
c            3:Toshiba 321)
c   ioptf1 : 1 (0) if SAVESC has (has not) been issued and no screen
c            has yet been saved.
c   ioptf2 : abortoff/aborton (0/1)
c   ioptf3 : 1 (0) if a macro has just been (not been) broken
c            from the graphics menu
c   iptrec : 1 (0) if the RECORD command is (is not) `active'
c   ioptsg : singleoff/singleon (0/1)
c   iptpon : plotoff/ploton (0/1)
c   imflag : 1 (0) if previous command was (was not) an error
c   imfl   : keeps track of imflag
c
        ioptsp=0
        ioptpr=1
        ioptov=1
        ioptv1=1
        ioptb=0
        ioptps=0
        ioptls=1
        ioptf1=0
        ioptf2=1
        ioptf3=0
        iptrec=0
        ioptsg=0
        iptpon=0
        imflag=0
        imfl=0
c%%
        iptas=0
        iptbr=0
        iptbr1=0
        ntmp1=0
c%%
c#####
        n456=0
c#####
c
c   background and foreground colors in text mode:
c
        kback=1
        kfore=7
c
c   i4seed is the internal seed for random number generators.
c   timef returns the number of seconds since midnight
c
        i4seed=time()
c
c   txtcol(assembler) sets text mode colors and clears the screen
c
c   iattr is used in wrtxta (assembler)
c
c        call txtcol(kback,kfore,1)
c
c  *********************************
c  *                               *
c  *  Create Opening screen:       *
c  *                               *
c  *********************************
c
        call logo
        
c
c  *********************************************************
c  *                                                       *
c  *   Form second opening screen if `tslabw1' is present: *
c  *                                                       *
c  *********************************************************
c
        fname='tslabw1'
        inquire(file=fname,exist=logic)
        if(logic) call plt801(fname,kolor)
c
c  *****************************************
c  *                                       *
c  *   Put list of commands on the screen: *
c  *                                       *
c  *****************************************
c
        call txtcol(kback,kfore,1)
        call movct(l,72,char(32))
        call movxy(l,'help(news)',10)
        lnb=10
        call type1(l,lun)
c        call movxy(l,'help',4)
c        lnb=4
c        call type5(l,iend,lun)

c
c  *****************************************
c  *                                       *
c  *   Execute autoexec.mac if present:    *
c  *                                       *
c  *****************************************
c
        fname='autoexec.mac'
        inquire(file=fname,exist=logic)
        if(logic) then
                l72='macro(autoexec.mac)'
                lnb=19
                call type1(l,lun)
                go to 10
        endif
c
c   make speaker beep:
c
c        call speak(500,10)
c
c  **************************
c  *                        *
c  *  Process next command: *
c  *                        *
c  **************************
c
  10    continue
c
c   Before reading new command, take care of some things
c   that couldn't be handled during previous command:
c
c
c   1) A macro was broken from graphics menu:
c      --------------------------------------
c
        if(ioptf3.eq.1) then
                if(lun.ne.100) then
                        lun1=lun
                        lun=100
                endif
                ioptf3=0
        endif
c
c   2) There was an error in the previous command:
c      --------------------------------------------
c
        if(imflag.eq.1) then
        imfl=1
c
c               a) Macro running and prompt off:
c                  ------------------------------
c
        if(ioptpr.eq.0.and.lun.ne.100) write(*,972)
 972    format(' ERROR OCCURED IN :')
        if(ioptpr.eq.0.and.lun.ne.100) write(*,72) (l(i),i=1,lnb)
c
c               b) Abort off:
c                  -----------
c
                if(ioptf2.eq.0) then
                        imflag=0
                        go to 15
                endif
c
c               c) Abort on and macro running:
c                  ----------------------------
c
                if(lun.gt.100) backspace lun
                if(lun.gt.100) lun1=lun
                lun=100
                imflag=0
c%%
                if(iptbr.eq.1) then
                        do 7101 i=1,ntmp1
                        call fnm(i,vn1,'%')
7101                    call clean1(vn1,ier)
                        iptbr=0
                        ntmp1=0
                endif
c%%
        endif
  15    continue
c%%
        if(iptbr1.eq.1) then
                do 7201 i=1,ntmp1
                call fnm(i,vn1,'%')
7201            call clean1(vn1,ier)
                iptbr=0
                iptbr1=0
                ntmp1=0
        endif
        if(iptbr.eq.1) go to 6900
c%%
c
c  ******************************
c  *                            *
c  *   Read next line of macro: *
c  *                            *
c  ******************************
c
        if(lun.gt.100) then
                n=inkeya()
c
c   check if ctrl-prtsc was struck:
c
                if(n.eq.29184) call echo
c
c   check if F9 was struck:
c
                if(n.eq.17152) then
                        lun1=lun
                        lun=100
                        go to 10
                endif
c
c   Read:
c
                read(lun,51,end=52) l72
  51            format(a72)
                go to 1012
c
c   macro has ended:
c
  52            close(lun)
                lun=lun-1
                lun1=lun
                go to 10
c
c   go back to 10 if a blank line or nondisplayed comment:
c
1012            do 1013 i=1,72
1013            if(l(i).ne.' ') go to 1014
                go to 10
1014            if(l(i).eq.';'.and.l(i+1).ne.';') go to 10
c
c   go to 71 if this line is not blank nor a nondisplayed comment
c
                go to 71
        endif
c
c  ************************************
c  *                                  *
c  *   Read next line from keyboard:  *
c  *                                  *
c  ************************************
c
c
c   PLOTON in effect:
c
        if(iptpon.eq.1) then
                call drrect(0,192,640,8,0,1)
                call stcura(0,24,0)
                call wrtxt1('?',1,0,192,1)
                call rdname(8,192,72,l,1)
                if(imfl.eq.1) then
                        call drrect(0,184,640,8,0,1)
                        imfl=0
                endif
                go to 9879
        endif
c
c   PLOTOFF in effect:
c
        write(*,11)
  11    format(1x,2h? ,$)
        read(*,12) l72
 12     format(a72)
c        call histry(n456,l,kback,kfore)
        if(l(1).eq.'!') write(*,865)
 865    format(/)
9879    continue
c
c  *****************************************************
c  *                                                   *
c  *   Now have command from either macro or keyboard  *
c  *                                                   *
c  *****************************************************
c
c   find last nonblank character:
c
  71    continue
        call lastnb(l,72,lnb)
        if(lnb.eq.0) go to 10
c
c   Check for DOS command:
c
        if(l(1).eq.'!') then
                do 581 i=1,lnb-1
 581            l(i)=l(i+1)
                l(lnb)=char(0)
                status=system(l72)
c                i=system(l72)
                go to 10
        endif
c
c   remove blanks except between apostrophes and on comment lines:
c
        if(l(1).ne.';') call rmblnk(l,72,lnb)
c
c   if record active, put line out to file:
c
        if(iptrec.eq.1) then
                write(3,9876,err=9877) (l(i),i=1,lnb)
9876            format(72a1)
9877            continue
        endif
c
c   handle prompton during macro:
c
        if(lun.ne.100.and.ioptpr.eq.1) then
                if(iptpon.eq.0) write(*,72) (l(i),i=1,lnb)
  72            format(1x,'?',72a1)
                if(iptpon.eq.1) then
                        call drrect(0,192,640,8,0,1)
                        call stcura(0,24,0)
                        call wrtxt1('?',1,0,192,1)
                        call wrtxt1(l,lnb,8,192,1)
                endif
        endif
c
c   handle speakeron during macro:
c
        if(lun.ne.100.and.ioptsp.eq.1) call speak(500,10)
c
c   handle singleon during macro:
c
        if(lun.ne.100.and.ioptsg.eq.1) then
                call waitk(n10)
                if(n10.eq.256*67) then
                        ioptf3=1
                        backspace lun
                        go to 10
                endif
        endif
c
c
c
        if(l(1).eq.';') go to 10
c
c  *************************************************
c  *                                               *
c  *   Trap commands that don't fit command types: *
c  *                                               *
c  *************************************************
c
c   *************
c   *           *
c   *   GOTO    *
c   *           *
c   *************
c
        call movxy(lck5,l,5)
        call ltoup(lck5,5)
        if(lck5.eq.'GOTO(') then
c
c   make sure a macro is running:
c
                if(lun.eq.100) then
                        write(*,284)
 284                    format(' GOTO only allowed in MACRO''s')
                        imflag=1
                        go to 10
                endif
c
c   now rewind file and look for target of GOTO:
c
                rewind lun
                call movct(vn1,15,' ')
                call movxy(vn1,l(6),lnb-6)
 285            read(lun,51,end=286) l72
                if(l(1).ne.';') go to 285
                call movct(vn2,15,' ')
                do 288 i=2,72
                if(l(i).eq.' ') go to 289
 288            continue
 289            call movxy(vn2,l(2),min0(i-2,15))
c
c   found it:
c
                if(vn1.eq.vn2) go to 10
                go to 285
c
c   doesn't exist:
c
 286            write(*,287)
 287            format(' Target of GOTO doesn''t exist')
                imflag=1
                go to 10
        endif
c%%
c
c   see if there are any {   }:
c
6900    continue
        if(iptbr.eq.1) then
                call movct(l,72,char(32))
                call movxy(l,ltmp1,ibl-1)
                if(iptas.eq.1) then
                        call movxy(l(ibl),ltmp1(ibl+1),nceq)
                        call movxy(l(ibl+nceq),ltmp1(ibr+1),lnb1-ibr)
                        lnb=ibl-1+nceq+lnb1-ibr
                        iptas=0
                else
                        l(ibl)='t'
                        l(ibl+1)='m'
                        l(ibl+2)='p'
                        l(ibl+3)='%'
                        write(l(ibl+4),'(i1)') ntmp1
                        call movxy(l(ibl+5),ltmp1(ibr+1),lnb1-ibr)
                        lnb=lnb1+ibl-ibr+4
                endif
        endif
        do 6901 j=1,lnb
6901    if(l(j).eq.'}') go to 6902
        go to 6990
6902    do 6903 i=j,1,-1
6903    if(l(i).eq.'{') go to 6904
        go to 6990
6904    iptbr=1
        ibl=i
        ibr=j
        call movct(ltmp1,72,char(32))
        call movxy(ltmp1,l,lnb)
        lnb1=lnb
        call movct(l,72,char(32))
        call findch(ltmp1,ibl+1,ibr-1,'=',noceq,nceq)
        if(noceq.eq.1) then
                iptas=1
                nceq=nceq-ibl-1
                lnb=ibr-ibl-1
                call movxy(l,ltmp1(ibl+1),lnb)
        else
                l(1)='t'
                l(2)='m'
                l(3)='p'
                l(4)='%'
                ntmp1=ntmp1+1
                write(l(5),'(i1)') ntmp1
                l(6)='='
                lnb=6
                do 6905 k=i+1,j-1
                lnb=lnb+1
6905            l(6+k-i)=ltmp1(k)
        endif
        go to 7000
6990    if(iptbr.eq.1) iptbr1=1
7000    continue
        if(l(1).eq.'@') go to 10
c%%
c
c   *************
c   *           *
c   *  WHILE    *
c   *           *
c   *************
c
        call movxy(lck,l,6)
        call ltoup(lck,6)
        if(lck.eq.'WHILE(') then
c
c   make sure a macro is running:
c
                if(lun.eq.100) then
                        write(*,302)
 302                    format(' WHILE only legal in MACRO')
                        go to 10
                endif
c
c   find the arguments of the WHILE and check for errors:
c
                call findch(l,1,lnb,')',nocrp,nrp)
                if(nocrp.ne.1.or.nrp.ne.lnb) then
                        write(*,351)
 351                    format(' Illegal WHILE')
                        imflag=1
                        go to 10
                endif
                call findch(l,1,lnb,'.',nocp,np1)
                if(nocp.ne.2.and.nocp.ne.3) then
                        write(*,351)
                        imflag=1
                        go to 10
                endif
                call findch(l,np1+1,lnb,'.',nocp,np2)
c
c   vn1, vn2, vn3 are WHILE(vn1.vn3.vn2):
c
                call formnm(l,7,np1-1,vn1)
                call formnm(l,np2+1,lnb-1,vn2)
                call formnm(l,np1+1,np2-1,vn3)
                call ltoup(vn3,15)
c
c   check that vn1 and vn2 are defined:
c
                call ckint(vn1,n1)
                call ckint(vn2,n2)
                if(n1.ne.-999.and.n2.ne.-999) then
                        x1=n1
                        x2=n2
                        go to 310
                endif
                call ckreal(vn1,x1)
                call ckreal(vn2,x2)
                if(x1.eq.999..or.x2.eq.999.) then
                        write(*,351)
                        imflag=1
                        go to 10
                endif
c
c   make sure legal relation being used:
c
 310    if(vn3.eq.'LT'.or.vn3.eq.'LE'.or.vn3.eq.'EQ') go to 314
        if(vn3.eq.'NE'.or.vn3.eq.'GE'.or.vn3.eq.'GT') go to 314
c
c   illegal:
c
        write(*,351)
        imflag=1
        go to 10
c
c   if condition is true, just keep going:
c
 314    if(vn3.eq.'LT'.and.x1.lt.x2) go to 10
        if(vn3.eq.'LE'.and.x1.le.x2) go to 10
        if(vn3.eq.'EQ'.and.x1.eq.x2) go to 10
        if(vn3.eq.'NE'.and.x1.ne.x2) go to 10
        if(vn3.eq.'GE'.and.x1.ge.x2) go to 10
        if(vn3.eq.'GT'.and.x1.gt.x2) go to 10
c
c   condition false, look for matching END. Have to keep track
c   of how many other WHILE's there are:
c
c   nwh and nend1 are the number of while's and end's encountered:
c
c   keep reading until nwh=nend1 (error if end of file hit first)
c
        nwh=1
        nend1=0
 315    read(lun,51,end=320) l72
        call lastnb(l,72,lnb)
        if(lnb.eq.0) go to 315
        call rmblnk(l,72,lnb)
        call movxy(lck,l,6)
        call ltoup(lck,6)
        if(lck.eq.'WHILE(') then
                nwh=nwh+1
                go to 315
                endif
        if(lnb.ne.3) go to 315
        call movxy(lck3,l,3)
        call ltoup(lck3,3)
        if(lck3.eq.'END') then
                nend1=nend1+1
                if(nend1.eq.nwh) go to 10
                endif
        go to 315
 320    write(*,321)
 321    format(' Matching END missing')
        imflag=1
        go to 10
        endif
c
c   *************
c   *           *
c   *   END     *
c   *           *
c   *************
c
        if(lnb.ne.3) go to 399
                call movxy(lck3,l,3)
                call ltoup(lck3,3)
                if(lck3.eq.'END') then
c
c   check that a macro is running:
c
                        if(lun.eq.100) then
                                write(*,322)
 322                            format(' END only legal in MACROS')
                                go to 10
                        endif
c
c   read backward in the file until find the matching WHILE
c
c   only read back as many as 200 lines since we can't tell if
c   we've hit beginning of file
c
c   nwh and nend1 are number of whiles and ends encountered
c
                        nend1=1
                        nwh=0
                        ncnt=0
 325                    ncnt=ncnt+1
                        if(ncnt.gt.200) then
                                write(*,326)
 326     format(' END not within 200 lines of WHILE or unmatched END')
                                imflag=1
                                go to 10
                        endif
c
c   have to backspace twice to get to the line before
c
                        backspace lun
                        backspace lun
                        read(lun,51) l72
                        call lastnb(l,72,lnb)
                        if(lnb.eq.0) go to 325
                        call rmblnk(l,72,lnb)
                        call movxy(lck3,l,3)
                        call ltoup(lck3,3)
                        if(lnb.eq.3.and.lck3.eq.'END') then
                                nend1=nend1+1
                                go to 325
                        endif
                        call movxy(lck,l,6)
                        call ltoup(lck,6)
                        if(lck.eq.'WHILE(') then
                                nwh=nwh+1
                                if(nend1.eq.nwh) then
                                        backspace lun
                                        go to 10
                                endif
                        endif
                        go to 325
                endif
 399    continue
c
c   ****************
c   *              *
c   *  DEDICATION  *
c   *              *
c   ****************
c
        call movxy(ded10,l,10)
        call ltoup(ded10,10)
        if(ded10.eq.'DEDICATION') then
                call txtcls
                write(*,971)
 971            format(' Thank you Karah, Tim, Mom, Dad, and Linda')
                go to 10
        endif
c
c   ****************************
c   *                          *
c   *  MACRO with no argument  *
c   *                          *
c   ****************************
c
        call formnm(l,1,6,lck)
        call ltoup(lck,6)
        if(lck.eq.'MACRO ') then
                lun=lun1
                go to 10
        endif
c
c   **************************
c   *                        *
c   *  Some forms of LABEL   *
c   *                        *
c   **************************
c
        call movxy(lck,l,6)
        if(lck.ne.'LABEL('.and.lck.ne.'label(') go to 199
        call findch(l,1,lnb,'=',nocceq,neq)
        if(nocceq.eq.0) go to 199
        if(ichar(l(neq+1)).ne.39.or.ichar(l(lnb)).ne.39) go to 199
        if(l(neq-1).ne.')') go to 199
c
c   now know it's LABEL(xxxxxx)='yyyyyyyyy'
c
c   pretend it's a type 3 command and pull off arguments:
c
        call parse(l,pname,args,nargs,vname,ierr)
c
c   handle errors:
c
        if(ierr.eq.1.or.nargs.gt.3) then
                write(*,116)
                imflag=1
                go to 10
        endif
        if(nargs.gt.1.and.iptpon.eq.0) then
                write(*,751)
 751            format(' Multi-argument LABEL illegal in text mode')
                imflag=1
                go to 10
        endif
c
c   Handle left hand side:
c
c
c       a) one argument:
c          -------------
c
        if(nargs.eq.1) then
                ipt2=1
                if(ickl(args(1),np,npx).eq.1) go to 10
                go to 752
        endif
c
c       b) two arguments:
c          --------------
c
        if(nargs.eq.2) then
                ipt2=2
                call ltoup(args(1),15)
                if(args(1).ne.'X'.and.args(1).ne.'Y'.
     1          and.args(1).ne.'CAP') then
                        write(*,116)
                        imflag=1
                        go to 10
                endif
                if(icki(args(2),2,npt2,0,0).eq.1) go to 10
                if(args(1).eq.'X') then
                        ipth=1
                        nr1=ny1-12-8*(npt2+1)
                        nc1=nx1
                endif
                if(args(1).eq.'Y') then
                        ipth=0
                        nc1=nx1-8*(nclaby)-2-8*npt2
                        nr1=ny1+npixy
                endif
                if(args(1).eq.'CAP') then
                        ipth=1
                        nc1=nx1+2
                        nr1=ny1+npixy+8*npt2+2
                endif
                go to 752
        endif
c
c       c) three arguments:
c          ----------------
c
        if(nargs.eq.3) then
                ipt2=3
                call ckint(args(1),nc1)
                call ckint(args(2),nr1)
                call ckint(args(3),ipth)
                if(ipth.ne.0.and.ipth.ne.1) then
                        call error(args,3,2)
                        go to 10
                endif
        endif
c
c   Now pull off right hand side:
c
c   nl and nr are boundaries of region within apostrophes
c
 752    continue
        nl=neq+2
        nr=lnb-1
c
c   see if there are any substitutions (i.e. #xx# or @xx@) being
c   asked for:
c
        call findch(l,nl,nr,'#',nocps,ncps)
        call findch(l,nl,nr,'@',nocat,nat)
c
c  handle no substitution case:
c
        if(nocps.le.1.and.nocat.le.1) then
                if(nr-nl+1.le.40) then
                        if(ipt2.ne.1) then
                                call movct(l80,80,' ')
                                call movxy(l80,l(nl),nr-nl+1)
                                go to 754
                        endif
                        call movct(label(np),40,' ')
                        call movxy(label(np),l(nl),nr-nl+1)
                        go to 10
                endif
                write(*,115)
 115            format(' label must be less than 41 characters')
                imflag=1
                go to 10
        endif
c
c   handle substitution case:
c
c   we do #xx# first and then come back to 114 to do @xx@:
c
c   lchr is either # or @
c
        lchr='#'
c
c  move what's in apostrophe's into l80:   nlst is number of elements
c
        call movct(l80,80,' ')
        call movxy(l80,l(nl),nr-nl+1)
        nlst=nr-nl+1
c
c   start: we'll create new array called l801 containing
c          expanded vesrion of l80.  jj is number of elements in
c          new array
c
 114    jj=0
c
c   find number of occurrences of lchr and location of first one:
c
        call findch(l80,1,nlst,lchr,nocps,ncps)
        if(nocps.eq.0) go to 120
c
c   has to be even number of them
c
        if(2*(nocps/2).ne.nocps) then
                write(*,116)
 116            format(' error in LABEL')
                imflag=1
                go to 10
        endif
c
c   find second occurrence:
c
        call findch(l80,ncps+1,nlst,lchr,nocps,ncps1)
c
c   can't be next to each other:
c
        if(ncps.eq.ncps-1) then
                write(*,116)
                imflag=1
                go to 10
        endif
c
c   copy what's before first occurrence of lchr # into l801
c   (unless there's nothing):
c
        if(ncps.gt.1) then
                do 117 i=1,ncps-1
                jj=jj+1
                if(jj.gt.40) go to 131
 117            l801(jj)=l80(i)
        endif
c
c   pull of what's inside lchr's:   (check if errors too)
c
        call formnm(l80,ncps+1,ncps1-1,v6)
c
c       a) lchr=#:
c          --------
c
        if(lchr.eq.'#') then
                call ckint(v6,n1)
c
c   undefined
c
                if(n1.eq.-999) then
                        write(*,116)
                        imflag=1
                        go to 10
                endif
                write(v7,118) n1
 118            format(i20)
c
c  v7 (and thus v8) now contain integer to insert so insert it:
c
                do 119 i=1,20
                if(v8(i).eq.' ') go to 119
                jj=jj+1
                if(jj.gt.40) go to 131
                l801(jj)=v8(i)
 119            continue
        endif
c
c       b) lchr=@:
c          --------
c
        if(lchr.eq.'@') then
        call ckreal(v6,x1)
c
c   undefined
c
        if(x1.eq.999.) then
                write(*,116)
                imflag=1
                go to 10
        endif
c
c   insert:
c
        write(v7,122) x1
 122    format(f20.10)
                do 123 i=20,1,-1
                if(v8(i).ne.'0') go to 124
 123            continue
 124            do 125 j=1,i
                if(v8(j).eq.' ') go to 125
                jj=jj+1
                if(jj.gt.40) go to 131
                l801(jj)=v8(j)
 125    continue
        endif
c
c
c   we've now handled one pair of lchr's
c
c   if there are more, copy what's left in original stuff between
c   apostrophes onto the end of l801
c
        if(ncps1.ne.nlst) then
                do 121 i=ncps1+1,nlst
                jj=jj+1
                if(jj.gt.40) go to 131
 121            l801(jj)=l80(i)
        endif
c
c   now l801 contains either a version with all lchr's expanded
c   or a partially expanded version with the rest of the original
c   stuff appended
c
c   copy l801 back to l80 and go back to start again
c
c   if we're actually finished, the code above will know it
c
        call movct(l80,80,' ')
        call movxy(l80,l801,jj)
        nlst=jj
        go to 114
c
c   131 can only be reached if we have reached more than 40
c   characters in the expansion process
c
 131    write(*,115)
        imflag=1
        go to 10
c
c   we get to 120 when there are no more occurrences of the current
c   value of lchr
c
c
c   if lchr=#, we replace it with @ and start all over again
c
c
 120    continue
        if(lchr.eq.'#') then
                lchr='@'
                jj=0
                go to 114
        endif
c
c   now we have the expanded right hand side
c
 754    continue
c
c   ipt2 is 1 if all we're doing is replacing an array's label
c
c   ipt2 is 2 if we're going to put the characters on the screen
c
c
c   ipth is 1 (0) if label goes horizontally  (vertically)
c
c   put label on the screen:
c
        if(ipt2.ne.1) then
                call lastnb(l80,80,nc2)
c
c   horizontally:
c
                if(ipth.eq.1) then
                        if(nargs.eq.2.and.args(1).eq.'X') then
                                nc1=nc1+((npixx-8*nc2)/2)
                        endif
                        if(nc1.lt.0.or.nc1+8*nc2.gt.639.or.nr1.lt.0.or.
     1                  nr1.gt.192) then
                                write(*,116)
                                imflag=1
                                go to 10
                        endif
                        call wrtxt1(l80,nc2,nc1,nr1,1)
                        go to 10
                endif
c
c   vertically:
c
                if(ipth.eq.0) then
                        if(nargs.eq.2) nr1=nr1-((npixy-8*nc2)/2)
                        nr2=nr1-8*nc2
                        if(nr2.lt.0.or.nr1.gt.199.or.nc1.lt.0.or.
     1                  nc1.gt.632) then
                                write(*,116)
                                imflag=1
                                go to 10
                        endif
                        do 755 i=1,nc2
 755                    call wrtxt1(l80(i),1,nc1,nr1-8*i,1)
                        go to 10
                endif
        endif
c
c   replace an array's label:
c
        call movct(label(np),40,' ')
        call movxy(label(np),l80,40)
        go to 10
c
c  *******************
c  *                 *
c  *   END OF LABEL  *
c  *                 *
c  *******************
c
c
 199    continue
c
c  ****************************************************************
c  *                                                              *
c  *   Now find command type and pass it to appropriate routine:  *
c  *                                                              *
c  ****************************************************************
c
        call typchk(l,lnb,itype)
        if(itype.ne.0) go to 30
        write(*,25)
  25    format(' unrecognizable command type ')
        imflag=1
        go to 10
  30    continue
c
c   go to appropriate command type section :
c
        go to(1000,2000,3000,4000,5000) itype
c
c   type 1 command section :
c
1000    continue
        call type1(l,lun)
        go to 10
c
c   type 2 command :
c
2000    continue
        call type2(l,args,vname)
        go to 10
c
c
c
c
c   type 3 command :
c
3000    continue
        call type3(l)
        go to 10
c
c   type 4 command :
c
4000    continue
        call type4(l)
        go to 10
c
c   type 5 command :
c
5000    continue
        call type5(l,iend,lun)
        if(iend.eq.1) then
                if(iptpon.eq.1) call smodea(3)
                call txtcls
                go to 6000
        endif
        if(iend.eq.2) go to 1
        go to 10
c
c
c
6000    continue
        stop
        end
