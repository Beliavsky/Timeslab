c&T1AM
c&DELAY
        subroutine delay(args,nargs)
c********************************************************************
c
c   Subroutine to delay for x seconds.
c
c   DELAY(x)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickr
c
c
	if(ickr(args(1),1,x,2,0.).eq.1) go to 99
c
c
        t=secnds(0.0)
  20    t1=secnds(0.0)-t
        if(t1.le.x) go to 20
c
c
  99    continue
        return
        end
c&HELP
        subroutine help(args,nargs,kback,kfore)
c**************************************************************
c
c   subroutine to list part of help file
c
c***************************************************************
c
        character args(nargs)*15,l*80,l1(80)*1,pname*15,nch1*1
        character names(5)*15
	character fname*15,fn1*15,fn11(15)*1,fn12*50
      character name*15,name1*15,name2*65
      character*1 name3(15)



        dimension locnms(40)
        logical logic
        equivalence (l,l1),(fn1,fn11)
        equivalence (name1,name3)
c
c
        pname=args(1)
        call ltoup(pname,15)

        fn12='help'
        write(*,901) fn12
 901    format(1x,a50)
        inquire(file=fn12,exist=logic)
        if(logic) go to 8
        fn12='/usr/local/lib/timeslab/help'
        write(*,901) fn12
        inquire(file=fn12,exist=logic)
        if(logic) go to 8
        write(*,1)
   1    format(' Can''t find help file (help)')
	call error(args,1,4)
        go to 99
   8    open(1,file=fn12,status='old')
c
c

      
 103    read(1,31,end=98) name1,name2
 31     format(a15,a65)

        call ltoup(args(1),15)
        if(args(1).eq.name1) then
         write(*,32)
 32      format(/)
 40      write(*,37) name1,name2
 37      format(1x,a15,a65)
         read(1,31,end=99) name1,name2
         if(name3(1).ne.' ') go to 99
         go to 40
      endif

      go to 103

 98   write(*,101) args(1)
 101  format(1x,a15,'   not found')
  99    continue
        close(1)
        return
        end
c&HIST
        subroutine hist(args,nargs)
c****************************************************************
c
c   Subroutine to produce the histogram of a set of numbers.
c
c   HIST(x,n)  or  HIST(x,n,nbins,xmin,xmax,ymax)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki,ickr
        dimension ninbin(100),rf(100),xww(4),yww(4)
	ierd=0
c
c
        if(nargs.ne.1.and.nargs.ne.2.and.nargs.ne.3.and.nargs.ne.6)
     1	then
           write(*,501)
 501 	   format(' wrong number of arguments in HIST')
           imflag=1
	   go to 99
        endif
c
	if(ickl(args(1),np,nn).eq.1) go to 99

        if(nargs.ne.1) then
           if(icki(args(2),2,n,5,nn).eq.1) go to 99
           iptsol=1
           if(n.lt.0) then
		n=-n
		iptsol=0
           endif
        endif

	if(nargs.eq.1) then
		n=nend(np)-nstart(np)+1
		nargs=2
		iptsol=1
	endif

        iptscx=1
        iptbin=1
        iptscy=1

        if(nargs.eq.2) then
           call min(array(nstart(np)),n,xmin,ind)
           call max(array(nstart(np)),n,xmax,ind)
           rx=xmax-xmin
           call xbars2(x,n,xbar,s2)
           nbins=(rx/(3.5*sqrt(s2)))*float(n)**.33+1
           binwid=rx/float(nbins)           
        endif

        if(nargs.gt.2) then
           iptbin=0
           if(icki(args(3),3,nbins,1,0).eq.1) go to 99
           if(nargs.eq.3) go to 301
           iptscx=0
           call max(array(nstart(np)),n,xmax1,ind)
           call min(array(nstart(np)),n,xmin1,ind)
           if(ickr(args(4),4,xmin,0,0.).eq.1) go to 99
           if(xmin.gt.xmin1) xmin=xmin1
           if(ickr(args(5),5,xmax,0,0.).eq.1) go to 99
           if(xmax.lt.xmax1) xmax=xmax1
           if(ickr(args(6),6,ymax,0,0.).eq.1) go to 99
           iptscy=1
           if(ymax.gt.0.) iptscy=0

        endif
 301	if(iptpon.eq.1) go to 973
c        call smodea(6)
c        call colora(0)
c        call colora(kolor)
 973	continue
           binwid=(xmax-xmin)/float(nbins)
      do 40 i=1,nbins
      ninbin(i)=0
  40  continue
      n1=nstart(np)
      do 50 i=1,n
      nbin=((array(n1+i-1)-xmin)/binwid)+1
  50  ninbin(nbin)=ninbin(nbin)+1
c
c   find relative frequencies :
c
      do 70 i=1,nbins
  70  rf(i)=float(ninbin(i))/(n*binwid)
      call max(rf,nbins,ymax1,ind)
      if(iptscy.eq.1) ymax=ymax1
      if(iptscy.eq.0.and.ymax1.gt.ymax) ymax=ymax1

      write(*,*) xmin,xmax,ymax
      call pgenv(xmin,xmax,0.0,ymax,0,0)
      if(iptsol.eq.1) call pgsfs(2)
      if(iptsol.eq.0) call pgsfs(1)

      do 80 i=1,nbins
         xww(1)=xmin+(i-1)*binwid
         xww(2)=xww(1)
         xww(3)=xmin+i*binwid
         xww(4)=xww(3)
         yww(1)=0.0
         yww(2)=rf(i)
         yww(3)=rf(i)
         yww(4)=0.0
         call pgpoly(4,xww,yww)
 80   continue
c
c        n1=nstart(np)
c        call histgm(array(n1),n,0,0,iptbin,iptscx,iptscy,iptsol,nbins,
c     1  binwid,npixx,npixy,nx1,ny1,nclabx,ndecx,1,1,
c     1  xmin,xmax,ymax,wk,wk(101),iwk,nticsy,nclaby,ndecy)
c        call wrtxt5(label(np),40,nx1+2,ny1+npixy+2,1)
c	call pltend(args,nargs,xmin,xmax,0.,ymax)
  99    continue
        return
        end
c&HISTGM
      subroutine histgm(x,n,iptrcl,iptgp,iptbin,iptscx,iptscy,iptsol,
     1nbins,binwid,npixx,npixy,nx1,ny1,nclabx,ndecx,kolor1,kolor2,
     1xmin,xmax,ymax,rf,rf1,ninbin,nticsy,nclaby,ndecy)
c*******************************************************************
c
c   Subroutine to draw a histogram of size NPIXX by NPIXY pixels
c   with lower left corner screen location (NX1,NY1).
c
c   IPTRCL : 1 means this histogram goes on the most recently drawn
c            axes, 0 means new axes.
c
c   IPTGP  : 1 means data grouped. Thus N=number of bins, X(i)= 
c            frequency of ith bin, and the user is responsible for
c            setting NBINS=N, and setting the values of XMIN,XMAX,
c            and BINWID=(XMAX-XMIN)/NBINS.
c            0 means data ungrouped.
c
c   IPTBIN : 1 means HISTGM finds NBINS by
c
c       NBINS=(RX/(3.5*SIG))*cube root(n) +1
c
c     where RX=XMAX-XMIN, SIG=std dev of X. 0 means user sets NBINS.
c
c     IPTBIN=1 is overridden by IPTRCL=1 or IPTGP=1.
c
c   IPTSCX : 1 means HISTGM finds XMIN,XMAX, 0 means user sets them.
c
c     IPTSCX=0 automatically if IPTRCL=1 .
c
c   IPTSCY : 1 means HISTGM finds YMAX, 0 means user sets it.
c
c     IPTSCY=0 automatically if IPTRCL=1 .
c
c   IPTSOL : 1 means bars are solid, 0 means they're not (0 doesn't
c            work properly yet).
c
c   NBINS  : number of bins (input if IPTRCL=1 or IPTGP=1 or IPTBIN=0,
c            output otherwise).
c
c   BINWID : bin width (input if IPTRCL=1 or IPTGP=1, otherwise 
c            BINWID=(XMAX-XMIN)/NBINS).
c
c   NPIXX,NPIXY,NX1,NY1 : X axis is of length NPIXX pixels, Y axis is
c                         of length NPIXY, they intersect at (NX1,NY1).
c                         These shouldn't be changed from 1st call if
c                         IPTRCL=1. NX1.ge.34, NY1.ge.10.
c
c   NCLABX,NDECX : Labels on X axis are in F NCLABX.NDECX format.
c                  Labels on Y axis are in F 4.3 format.
c
c   KOLOR1,KOLOR2 : Color numbers in current palette of axes and bars
c                   (see GPLOT documentation for details).
c
c   XMIN,XMAX : The values of the left and right ends of the X axis.
c               c
c   YMAX : The value at the top of the Y axis.
c
c   RF(1),...,RF(NBINS) : The relative frequencies divided by the
c                         BINWID. These values determine the height
c                         of the bars relative to the values on the 
c                         Y axis.
c
c   IX(1),...,IX(N) : The bin numbers of the X's. These are only
c                     found if the X's are ungrouped.
c
c*******************************************************************
c
c   NOTES :
c
c      1) HISTGM is designed so that GPLOT can be called subsequently
c         to superimpose a smooth curve onto the histogram.
c
c      2) If IPTRCL=1, the new histogram is formed by forming bars
c         half as wide as the old bars in the left half of the old 
c         bars. Thus in medium resolution mode, make sure that KOLOR2
c         for the two calls are different.
c
c      3) If you're superimposing widely differing data sets, the
c         number of bins chosen will be very large and the bars may
c         be very narrow.
c
c*******************************************************************
c
      dimension x(n),rf(100),ninbin(100),rf1(100)
c
c   make sure xmin,xmax,binwid,nbins are defined (if data is grouped,
c   these must all be specified on inpu) :
c
      if(iptrcl.eq.0) go to 6
      do 5 i=1,nbins
   5  rf1(i)=rf(i)
   6  continue
      ndata=n
      if(iptgp.eq.0) go to 20
      ndata=0
      do 10 i=1,n
      ndata=ndata+x(i)
  10  ninbin(i)=x(i)
      go to 60
  20  if(iptscx.eq.0) go to 30
      call max(x,n,xmax,ind)
      call min(x,n,xmin,ind)
  30  continue
      rx=xmax-xmin
      if(iptbin.eq.0) go to 35
      call xbars2(x,n,xbar,s2)
      nbins=(rx/(3.5*sqrt(s2)))*float(n)**.33+1
  35  continue
      binwid=rx/float(nbins)
c
c   assign x's to bins :
c
      do 40 i=1,nbins
      ninbin(i)=0
  40  continue
      do 50 i=1,n
      nbin=((x(i)-xmin)/binwid)+1
  50  ninbin(nbin)=ninbin(nbin)+1
c
c   find relative frequencies :
c
  60  continue
      on=ndata
      do 70 i=1,nbins
  70  rf(i)=float(ninbin(i))/(on*binwid)
	call max(rf,nbins,ymax1,ind)
	if(iptscy.eq.1) ymax=ymax1
	 if(iptscy.eq.0.and.ymax1.gt.ymax) ymax=ymax1
      if(iptrcl.eq.1) go to 80
c      call pgenv(xmin,xmax,0.,ymax,0,0)
c      call xaxis(nx1,ny1,npixx,nbins,kolor1,xmin,xmax,0,nclabx,ndecx)
c      call yaxis(nx1,ny1,npixy,nticsy,kolor1,0.0,ymax,0,nclaby,ndecy)
  80  continue
      inc=npixx/nbins
      inc1=inc
      if(iptrcl.eq.1) inc1=inc/2
      onpixy=npixy
      do 90 i=1,nbins
      i1=1+(onpixy*rf(i)/ymax)
      i2=nx1+(i-1)*inc
      call drrect(i2,ny1+1,inc1,i1,kolor2,iptsol)
      if(iptsol.eq.0) go to 81
      call drrect(i2+inc1-1,ny1+1,1,i1,0,1)
  81  continue
      if(iptrcl.eq.0) go to 90
      i3=1+(onpixy*rf1(i)/ymax)
      if(i3.gt.i1) call drrect(i2,ny1+i1+1,inc1,i3-i1,0,1)
  90  continue
c
c
      return
      end
c&LISTSP
        subroutine listsp(args,nargs)
c*******************************************************************
c
c   Subroutine to list the first (nfreqs/2)+1 values of a spectral
c   array.
c
c   listsp(spec,nfreqs)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl
c
c   check arguments :
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),nfreqs)
        n1=(nfreqs/2)+1
        if(n1.gt.nn) then
        	call error(args,2,2)
		go to 99
	endif
c
c   produce list :
c
        write(*,25) label(np)
  25    format(1x,a40)
        write(*,26)
c 26    format(/,1x,'Freq',3x,'Value',1x,'Freq',3x,'Value',1x,
c    1  'Freq',3x,'Value',1x,'Freq',3x,'Value',1x,'Freq',3x,
c    1  'Value',1x,'Freq',3x,'Value'/1x,78(1h-))
  26    format(/,1x,'Freq',3x,'Value','|','Freq',3x,'Value','|',
     1  'Freq',3x,'Value','|','Freq',3x,'Value','|','Freq',3x,
     1  'Value','|','Freq',3x,'Value'/1x,78(1h-))
        nn=nstart(np)-1
        on=float(nfreqs)
c       do 30 i=1,n1
c       wk(i)=float(i-1)/on
c       wk(n1+i)=array(nn+i)
c 30    continue
        write(*,40) (float(i-1)/on,array(nn+i),i=1,n1)
c       write(*,40) (wk(i),wk(n1+i),i=1,n1)
  40    format(1x,f4.3,1x,f7.2,'|',f4.3,1x,f7.2,'|',f4.3,1x,f7.2,
     1  '|',f4.3,1x,f7.2,'|',f4.3,1x,f7.2,'|',f4.3,1x,f7.2)
c
c
  99    continue
        return
        end
c&LISTTS
        subroutine listts(args,nargs)
c**********************************************************************
c
c   Subroutine to handle the LIST command.
c
c**********************************************************************
c
#include 'tslabc'
        character*15 args(nargs)
        character lf*30,lf1(30)*1,aa(15)*1,aa1*15,l31*31,iprn*1
	integer*2 ickl,icki
	logical logic
        equivalence (lf,lf1),(aa,aa1)
        lf='(1x,i5,2h |,5g14.6)'
	iptful=0
	iun=0
	iptprn=0
	aa1=args(nargs)
	call lastnb(aa1,15,lnbprn)
	if(aa(1).eq.char(39).and.aa(lnbprn).eq.char(39)) then
		nargs=nargs-1
		iun=1
		aa(lnbprn)=' '
		lnbprn=lnbprn-2
		do 1 i=1,lnbprn
   1		aa(i)=aa(i+1)
		aa(lnbprn+1)=' '
		call ltoup(aa1,15)
		if(aa1.eq.'PRN'.or.aa1.eq.'LPT1') then
			open(1,file='prn')
			go to 55
		endif
		inquire(file=aa1,exist=logic)
		if(logic) then
			open(1,file=aa1,status='old')
			iptprn=1
   5			read(1,6,end=7) iprn
   6			format(a1)
			go to 5
   7			backspace 1
			go to 55
		endif
		open(1,file=aa1,status='new')
		iptprn=1
  55		continue
	endif
c
c   check if all arguments are scalars:
c
	do 200 i=1,nargs
	call ckreal(args(i),x1)
 	call cklist(args(i),ilist,nivars,np2)
 	call cklist(args(i),chlist,nchars,np3)
	if(x1.eq.999..and.np2.eq.-999.and.np3.eq.-999) go to 201
 200	continue
	go to 210
c
c   no:
c
 201 	call cklist(args(1),xlist,numarr,np)
	if(np.eq.-999.or.nargs.eq.3.or.nargs.gt.4) then
		write(*,202)
 202	format(' Variable Undefined or Arrays and Scalars are Mixed')
		imflag=1
		go to 99
	endif
	go to 220
c
c   yes:
c
 210	continue
	do 215 i=1,nargs
 	call cklist(args(i),ilist,nivars,np)
	if(np.ne.-999) then
		write(l31,'(a15,1h=,i15)') args(i),ints(np)
		call rmblnk(l31,31,l31b)
		write(iun,216,err=57) l31
		if(iptrec.eq.1) then
			write(3,216,err=801) l31
			go to 901
 801 			iptful=1
 901			continue
		endif
		go to 215
	endif
	call ckreal(args(i),x1)
	if(x1.ne.999.) then
		write(l31,'(a15,1h=,g14.6)') args(i),x1
		call rmblnk(l31,31,l31b)
		write(iun,216,err=57) l31
 216		format(1x,a31)
		if(iptrec.eq.1) then
			write(3,216,err=802) l31
			go to 902
 802 			iptful=1
 902			continue
		endif
		go to 215
	endif
 	call cklist(args(i),chlist,nchars,np)
	if(np.ne.-999) then
		write(l31,'(a15,1h=,a15)') args(i),chars(np)
		write(iun,216,err=57) l31
		if(iptrec.eq.1) then
			write(3,216,err=803) l31
			go to 903
 803 			iptful=1
 903			continue
		endif
		go to 215
	endif
 215	continue
	go to 99
c
c
c
 220	continue
c
c
c   check whether array exists :
c
	if(ickl(args(1),np,nn).eq.1) go to 99
c
c   check whether n is properly defined :
c
	if(nargs.eq.1) then
		n=nn
		nargs=2
		go to 221
	endif
	if(icki(args(2),2,n,1,0).eq.1) go to 99
        if(n.gt.nn) then
	        call error(args,2,2)
        	go to 99
	endif
c
c   do list :
c
 221    write(iun,60,err=57) label(np)
  60    format(1x,a40)
	if(iptrec.eq.1) then
		write(3,60,err=804) label(np)
		go to 904
 804 		iptful=1
 904		continue
	endif
        if(n.eq.0) go to 99
        npln=5
        if(nargs.eq.2) go to 75
c
c   Get format from 3rd and 4th arguments:
c
	if(icki(args(3),3,npln,1,1).eq.1) go to 99
        if(nargs.eq.3) then
		if(npln.gt.5) then
			write(*,789)
 789    format(' You must use a FORMAT if you want 6 or more per line')
		if(iptrec.eq.1) then
			write(3,789,err=805)
			go to 905
 805 			iptful=1
 905			continue
		endif
			imflag=1
			go to 99
		endif
	go to 75
	endif
        call cklist(args(4),chlist,nchars,npc)
        if(npc.ne.-999) args(4)=chars(npc)
        call movct(lf1(13),18,char(0))
        aa1=args(4)
        call lastnb(args(4),15,lastf)
        do 66 i=1,lastf
  66    lf1(12+i)=aa(i)
        lf1(lastf+12+1)=')'
  75    continue
        nn=((n-1)/npln)+1
        ii=nstart(np)-1
        do 70 i=1,nn
        i1=npln*(i-1)+1
        i2=min0(n,npln*i)
        write(iun,lf,err=57) npln*(i-1)+1,(array(j),j=ii+i1,ii+i2)
  	if(iptrec.eq.1) then
	write(3,lf,err=806) npln*(i-1)+1,(array(j),j=ii+i1,ii+i2)
		go to 70
 806 		iptful=1
	endif
  70	continue
        go to 99
  57    write(*,58)
  58	format(' Error in format statement or writing to a file')
	imflag=1
c
c
  99    continue
	if(iptprn.eq.1) close(1)
	if(iptful.eq.1) then
		write(*,807)
 807		format(1x,'RECORD file full, I''m closing it')
		close(3)
		iptrec=0
		imflag=1
	endif
        return
        end
