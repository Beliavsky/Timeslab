c&PLOT
c&PLOTTS
        subroutine plotts(args,nargs)
c*******************************************************************
c
c   Subroutine to plot an array.
c
c*******************************************************************
c
#include 'tslabc'
        character*15 args(nargs)
	integer*2 ickl,icki,ickr

	ierd=0
c
c   nargs=3 :
c
        ipttyp=20
        if(nargs.ne.3.and.nargs.ne.7) go to 201
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(ickl(args(2),npy,ny).eq.1) go to 99
	if(icki(args(3),3,n,0,0).eq.1) go to 99
        if(n.lt.0) ipttyp=34
        if(n.lt.0) n=-n
        if(n.le.0.or.n.gt.nx.or.n.gt.ny) then
	        call error(args,3,2)
        	go to 99
	endif
	if(ickse(n).eq.1) go to 99
	call movxy(wk,array(nstart(npx)),4*n)
        nn=nstart(npy)
        iptsc=0
        call min(wk,n,xmin,ind)
        call max(wk,n,xmax,ind)
        call min(array(nn),n,ymin,ind)
        call max(array(nn),n,ymax,ind)
        if(nargs.eq.3) go to 42
c
c   nargs=7 :
c
	if(ickr(args(4),4,xmin,0,0.).eq.1) go to 99
	if(ickr(args(5),5,xmax,0,0.).eq.1) go to 99
	if(ickr(args(6),6,ymin,0,0.).eq.1) go to 99
	if(ickr(args(7),7,ymax,0,0.).eq.1) go to 99
        iptsc=1
        nsx=nstart(npx)
        nsy=nstart(npy)
        call max(array(nsx),n,xmax1,ind)
        call min(array(nsx),n,xmin1,ind)
        call max(array(nsy),n,ymax1,ind)
        call min(array(nsy),n,ymin1,ind)
	if(xmax1.gt.xmax) xmax=xmax1
	if(xmin1.lt.xmin) xmin=xmin1
	if(ymax1.gt.ymax) ymax=ymax1
	if(ymin1.lt.ymin) ymin=ymin1
        go to 42
c
c   nargs=2 or 4 :
c
 201    continue
        if(nargs.ne.2.and.nargs.ne.4) go to 50
	if(ickl(args(1),np,na).eq.1) go to 99
	if(icki(args(2),2,n,0,0).eq.1) go to 99
        if(n.lt.0) ipttyp=34
        if(n.lt.0) n=-n
        if(n.gt.na) then
	     	call error(args,2,2)
        	go to 99
	endif
        nn=nstart(np)
        iptsc=0
        if(nargs.eq.2) go to 34
        iptsc=1
	xmin=1.
	xmax=n
	if(ickr(args(3),3,ymin,0,0.).eq.1) go to 99
	if(ickr(args(4),4,ymax,0,0.).eq.1) go to 99
	if(ymin.gt.ymax) then
		call error(args,4,2)
		go to 99
	endif
        call min(array(nn),n,ymin1,ind)
        call max(array(nn),n,ymax1,ind)
	if(ymin1.lt.ymin) ymin=ymin1
	if(ymax1.gt.ymax) ymax=ymax1
  34	continue
        do 40 i=1,n
  40    wk(i)=i
  42    continue
c
c
	if(iptpon.eq.1) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973	continue
c----------------
        if(nargs.eq.2) then
           xmin=1.
           xmax=n
           call min(array(nn),n,ymin,ind)
           call max(array(nn),n,ymax,ind)
        endif

        call loadpg(xmin,xmax,ymin,ymax)
        call pgenv(xmin,xmax,ymin,ymax,0,0)
        if(ipttyp.ne.20) call pgpt(n,wk,array(nn),-3)
        if(ipttyp.le.30) call pgline(n,wk,array(nn))
        if(nargs.eq.2.or.nargs.eq.4) call pglab('index','value',
     1   label(np))
        if(nargs.eq.3.or.nargs.eq.7) call pglab(label(npx),label(npy),
     1    'x-y plot')
        go to 99
c----------------
c	call gplot(wk,array(nn),n,npixx,npixy,nx1,ny1,0,nticsx,
c     1  nticsy,nclabx,ndecx,nclaby,ndecy,iptsc,xmin,xmax,ymin,ymax,
c     1	0,ipttyp,1,1)
	if(iptpon.eq.1) nplt=nplt+1
        if(nargs.eq.3.or.nargs.eq.7) then
		lab='    '
		if(lab.eq.label(npx)) call wrtxt5(label(npy),
     1		40,nx1+2,ny1+npixy+2,1)
		if(lab.ne.label(npx)) call wrtxt5(label(npy),
     1		40,nx1+2,ny1+npixy+10,1)
		if(lab.ne.label(npx)) then
        	call wrtxt1('VS  ',4,nx1+2,ny1+npixy+2,1)
        	call wrtxt5(label(npx),40,nx1+34,ny1+npixy+2,1)
 	endif
	else
        	call wrtxt5(label(np),40,nx1+2,ny1+npixy+2,1)
	endif
	call pltend(args,nargs,xmin,xmax,ymin,ymax)
	go to 99
  50    continue
        write(*,4)
   4    format(' only 2,3,4, or 7 argument plots implemented so far')
	imflag=1
  99    continue
        return
        end
c&PLTEND
	subroutine pltend(args,nargs,xmin,xmax,ymin,ymax)
c*****************************************************************
c
c
c*****************************************************************
c
#include 'tslabc'
	character args(nargs)*15
        if(ioptps.eq.1) call prtscn
	if(ioptf1.eq.1) then
		call rdscr(charsv,wk,ierd)
		ioptf1=0
		endif
        if(ioptb.eq.1) go to 41
	icol=500
	irow=150
	call f5f9(icol,irow,5,iwk,xmin,xmax,ymin,ymax,
     1  npixx,npixy,nx1,ny1,wk,1,kolor,ioptf3)
  41    continue
	if(iptpon.eq.1) then
		call drrect(0,0,640,8,0,1)
		call stcura(0,24,0)
		go to 99
	endif
        call smodea(3)
        call txtcls
	if(ierd.eq.1) call error(args,1,5)
  99	continue
	return
	end
c&WRTXT5
	subroutine wrtxt5(l,nchar,ncol,nrow,kolor)
	character l(1)*1
	do 10 i=nchar,1,-1
  10	if(l(i).ne.' ') go to 20
	nc=1
	go to 30
  20	nc=i
  30	call wrtxt1(l,nc,ncol,nrow,kolor)
	return
	end
c&F5F9
	subroutine f5f9(icol,irow,ncol,iwk,xmin,xmax,ymin,ymax,
     1	npixx,npixy,nx1,ny1,wk,iptcr,kolor,ioptf3)
c*******************************************************************
c
c   Subroutine to issue the Print screen and FIND messages on the
c   graphics screens and take appropriate action.
c
c*******************************************************************
c
	dimension iwk(1),wk(1)
	logical lgc
	character fnscn*15
c
c
c
	icol1=icol
	irow1=irow
 415	call drrect(0,0,640,8,1,1)
c   	call wrtxt1('F1:FIND  F3:PRT SCN  F5:SAVE SCN ',33,0,0,129)
c	call wrtxt1(' F9:BREAK MACRO ',16,264,0,129)
c	call wrtxt1(' ANYTHING ELSE:RETURN TO PROMPT',31,392,0,129)
    	call wrtxt3(' F1:FIND  F3:PrtSc  F5:SAVE SCN ',33,0,0,129)
	call wrtxt3(' F9:BREAK MACRO ',16,33,0,129)
	call wrtxt3(' ANY OTHER KEY:RETURN TO PROMPT',31,49,0,129)
        call waitk(n10)
	n10=n10/256
	if(n10.ne.63.and.n10.ne.59.and.n10.ne.61.and.n10.ne.67) go to 99
	if(n10.eq.67) then
		ioptf3=1
		go to 99
	endif
	if(n10.eq.63) then
		call drrect(0,0,640,8,0,1)
		call wrtxt1('Filename?',9,0,0,1)
		call rdname(80,0,40,fnscn,2)
		call drrect(0,0,640,8,0,1)
			call rdscr(fnscn,wk,ier)
			if(ier.eq.1) then
				call wait10
			endif
			go to 415
	endif
	if(n10.eq.61) then
		call drrect(0,0,640,8,0,1)
		call prtscn
		go to 415
		endif
	if(n10.eq.59) then
		call drrect(0,0,640,8,1,1)
	call wrtxt3(' Alt-F1: HELP, Alt-F10: Exit FIND',33,0,0,129)
	call drrect(icol1,irow1,ncol,ncol,0,1)
	icol1=icol
	irow1=irow
c        call findxy(icol1,irow1,ncol,iwk,xmin,xmax,ymin,ymax,
c     1  npixx,npixy,nx1,ny1,iptcr,kolor,wk)
		call drrect(0,0,640,8,0,1)
		go to 415
		endif
  99	continue
	return
	end
c&WAIT10
	subroutine wait10
	call drrect(0,0,640,8,0,1)
	call wrtxt1('Not enough room on disk to save screen',38,0,0,1)
	call wrtxt1(', Strike any key to continue',28,304,0,1)
	call waitk(n10)
	return
	end
c&RDNAME
	subroutine rdname(icol,irow,nchar,fname,iopt)
c****************************************************************
c
c   Subroutine to read and echo a character screen on a mode 6 screen
c   starting at pixel location (icol,irow).
c
c*****************************************************************
c
	character fname(nchar)*1
c
c
	if(iopt.eq.1) call movct(fname,nchar,char(32))
	if(iopt.eq.2) call movct(fname,nchar,char(0))
	ncol=icol
	nn=0
c  5	call waitk(n)
   5	call drrect(ncol,irow,8,2,1,1)
	t1=timef()
   6	t=timef()-t1
	if(t.le..1) go to 6
	call drrect(ncol,irow,8,2,0,1)
	t1=timef()
   7	t=timef()-t1
	if(t.le..1) go to 7
	n=inkeya()
	if(n.eq.0) go to 5
	n10=n/256
	if(n10.ne.0) go to 5
	if(n.eq.13) go to 99
	if(n.eq.8) then
		if(nn.eq.0) go to 5
		if(iopt.eq.1) fname(nn)=char(32)
		if(iopt.eq.2) fname(nn)=char(0)
		nn=nn-1
		ncol=ncol-8
		call drrect(ncol,irow,8,8,0,1)
		go to 5
		endif
c	if(n.lt.32.or.n.gt.127) go to 5
	nn=nn+1
	call wrtxt1(char(n),1,ncol,irow,1)
	fname(nn)=char(n)
	ncol=ncol+8
	go to 5
  99	fname(nn+1)=char(32)
	return
	end
c&PLOTSP
        subroutine plotsp(args,nargs,numplt,ierr)
c****************************************************************
c
c   Subroutine to plot the log of normalized spectra on a scale
c   of -6 to 6.
c
c   plotsp(spec,nfreqs,var)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
c
c   check arguments :
c
        ierr=1
        call cklist(args(1),xlist,numarr,np)
        if(np.gt.0) go to 10
        write(*,5) args(1)
 5      format(1x,a15,' undefined')
c        call wrtxt1(args(1),15,0,192,1)
c        call wrtxt1('undefined',9,120,192,1)
	imflag=1
        go to 99
  10    call ckint(args(2),nfreqs)
        nn=nend(np)-nstart(np)+1
        n1=(nfreqs/2)+1
        if(nn.ge.n1) go to 20
        write(*,15) args(2)
 15     format(1x,a15,' illegal value of nfreqs')
c        call wrtxt1(args(2),15,0,192,1)
c        call wrtxt1(' illegal value of nfreqs',25,120,192,1)
	imflag=1
        go to 99
  20    call ckreal(args(3),var)
        if(var.ne.999.) go to 40
        write(*,25) args(3)
 25     format(1x,a15,'illegal value of var')
c        call wrtxt1(args(3),15,0,192,1)
c        call wrtxt1(' illegal value of var in PLOTSP',31,120,192,1)
	imflag=1
        go to 99
  40    continue
c
c   prepare for plotting :
c
        on=nfreqs
        do 50 i=1,n1
  50    wk(i)=float(i-1)/on
        div=var
        aldiv=alog(div)
        n3=nstart(np)-1
        emin=div*exp(-6.)
        emax=div*exp(6.0)
        do 60 i=1,n1
        c=array(n3+i)
        if(c.gt.emax) go to 55
        if(c.lt.emin) go to 56
        c=alog(c)-aldiv
        go to 60
  55    c=6.
        go to 60
  56    c=-6.
  60    wk(n1+i)=c
c
c   plot :
c
        call pgline(n1,wk(1),wk(n1+1))
        call pglab('frequency','log spectra',label(np))
c        call gplot(wk(1),wk(n1+1),n1,
c     1  npixx,npixy,nx1,ny1,0,nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,
c     1  1,0.0,0.5,-6.0,6.0,0,20,1,1)
c        call wrtxt5(label(np),40,nx1+2,ny1+npixy+2+(numplt-1)*8,1)
        ierr=0
  99    continue
        return
        end
c&PLOT2
        subroutine plot2(args,nargs)
c*******************************************************************
c
c   Subroutine to plot two arrays on the same axes.
c
c   plot2(y1,y2,n,m,iopty1,iopty2) or
c   plot2(y1,y2,n,m,ioptx,iopty,x) or
c   plot2(y1,y2,n,m,iopty1,iopty2,ymin,ymax) or
c   plot2(y1,y2,n,m,iopty1,iopty2,x,ymin,ymax)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki
        dimension xww(2),yww(2)
	ierd=0
c
c   check input :
c
	if(ickl(args(1),npy1,ny11).eq.1) go to 99
	if(ickl(args(2),npy2,ny2).eq.1) go to 99
	if(icki(args(3),3,n,5,ny11).eq.1) go to 99
	if(icki(args(4),4,m,5,ny2).eq.1) go to 99

	if(icki(args(5),5,ioptx,0,0).eq.1) go to 99
	if(icki(args(6),6,iopty,0,0).eq.1) go to 99

	ipte=0

        if(ioptx.lt.1.or.ioptx.gt.4) ipte=5
	if(iopty.lt.1.or.iopty.gt.4) ipte=6

	if(ipte.ne.0) then
		call error(args,ipte,2)
		go to 99
	endif

        ijoin=0
        if(ioptx.eq.3.and.iopty.eq.3) then
                ijoin=1
                ioptx=1
                iopty=1
        endif

        nx=nstart(npy1)
        ny=nstart(npy2)
        call max(array(nx),n,vmax,ind)
        call max(array(ny),m,vv,ind)
        if(vv.gt.vmax) vmax=vv
        call min(array(nx),n,vmin,ind)
        call min(array(ny),m,vvl,ind)
        if(vvl.lt.vmin) vmin=vvl
        hmin=1
        hmax=max0(n,m)
        do 50 i=1,max0(n,m)
  50    wk(i)=i

	if(nargs.eq.10.or.nargs.eq.11) then
		call ckreal(args(nargs-1),vmin1)
		call ckreal(args(nargs),vmax1)
		if(vmin1.le.vmin) vmin=vmin1
		if(vmax1.ge.vmax) vmax=vmax1
		call ckreal(args(nargs-3),hmin1)
		call ckreal(args(nargs-2),hmax1)

		if(nargs.eq.11) then
                   hmin=hmin1
                   hmax=hmax1
                   go to 51
		endif
		if(hmin1.le.hmin) hmin=hmin1
		if(hmax1.ge.hmax) hmax=hmax1
	endif

  51	continue
        if(nargs.eq.6.or.nargs.eq.10) go to 57
	if(nargs.eq.7.or.nargs.eq.11) then
		if(ickl(args(7),npz,nz).eq.1) go to 99
                        if(nz.lt.max0(n,m)) then
                        call error(args,7,2)
                        go to 99
                endif
                call movxy(wk,array(nstart(npz)),4*max0(n,m))
                call max(wk,max0(n,m),hmax1,ind)
                call min(wk,max0(n,m),hmin1,ind)
		if(nargs.eq.7) then
			hmax=hmax1
			hmin=hmin1
 		go to 57
		endif
		if(hmax1.gt.hmax) hmax=hmax1
		if(hmin1.lt.hmin) hmin=hmin1
	endif
c
c
c
  57    continue
	if(iptpon.eq.1) go to 973
c        call smodea(6)
c        call colora(0)
c        call colora(kolor)
 973	continue


        iptx=20
        if(ioptx.eq.2) iptx=34
	if(ioptx.eq.4) iptx=14
        ipty=20
        if(iopty.eq.2) ipty=34
	if(iopty.eq.4) ipty=14
                if(ijoin.eq.1) then
                iptx=34
                ipty=34
                endif

        call pgenv(hmin,hmax,vmin,vmax,0,0)
        if(iptx.ne.20) call pgpt(n,wk,array(nx),-3)
        if(iptx.le.30) call pgline(n,wk,array(nx))
        if(ipty.ne.20) call pgpt(m,wk,array(ny),-3)
        if(ipty.le.30) call pgline(m,wk,array(ny))
        if(ijoin.eq.1) then
           do 71 i=1,min0(n,m)
           xww(1)=wk(i)
           xww(2)=wk(i)
           yww(1)=array(nx+i-1)
           yww(2)=array(ny+i-1)
 71        call pgline(2,xww,yww)
        endif
c        call gplot(wk,array(nx),n,
c     1  npixx,npixy,nx1,ny1,0,nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,
c     1  1,hmin,hmax,vmin,vmax,0,iptx,1,1)
c        call wrtxt5(label(npy1),40,nx1+2,ny1+npixy+2,1)
c        call gplot(wk,array(ny),m,
c     1  npixx,npixy,nx1,ny1,-1,nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,
c     1  1,hmin,hmax,vmin,vmax,0,ipty,1,1)
	if(iptpon.eq.1) nplt=nplt+1
c        call wrtxt5(label(npy2),40,nx1+2,ny1+npixy+10,1)
c
c                if(ijoin.eq.1) then
c                do 71 i=1,min0(n,m)
c  71            call linexy(wk(i),array(nx+i-1),wk(i),array(ny+i-1),
c     1          npixx,npixy,nx1,ny1,hmin,hmax,vmin,vmax,1)
c                endif
c
c	call pltend(args,nargs,hmin,hmax,vmin,vmax)
	go to 99
  99    continue
        return
        end
c&LINEXY
        subroutine linexy(x1,y1,x2,y2,npixx,npixy,nx1,ny1,xmin,xmax,
     1  ymin,ymax,kolor)
c********************************************************************
c
c   Subroutine to draw a line from the real point (x1,y1) to the real
c   point (x2,y2) on axes of length npixx,npixy having pixel
c   of intersection at (nx1,ny1).  The values of x and y at the start
c   and end of the axes are determined by xmin,xmax,ymin,ymax.  The
c   line is drawn in color number kolor.
c
c*******************************************************************
c
        rx=xmax-xmin
        ry=ymax-ymin
        if(rx.gt.1.e-15.and.ry.gt.1.e-15) go to 15
        call wrtxt1('Range of X or Y is Zero',23,0,192,1)
        go to 99
  15    continue
c
c   find screen location of x,y :
c
        onpixx=npixx
        onpixy=npixy
        i2=1+(onpixy*(y1-ymin)/ry)
        iy1=ny1+min0(npixy,i2)
        i1=1+(onpixx*(x1-xmin)/rx)
        ix1=nx1+min0(npixx,i1)
        i3=1+(onpixy*(y2-ymin)/ry)
        iy2=ny1+min0(npixy,i3)
        i4=1+(onpixx*(x2-xmin)/rx)
        ix2=nx1+min0(npixx,i4)
        call line(ix1,iy1,ix2,iy2,kolor)
  99    continue
        return
        end
c&PLTCSP
        subroutine pltcsp(args,nargs)
c*****************************************************************
c
c   Subroutine to plot the first (nfreqs/2)+1 values of a spectral
c   array.
c
c   PLOTCSP(spec,nfreqs)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki,ickr
	ierd=0
c
c
	if(ickl(args(1),np,nn).eq.1) go to 99
        call ckint(args(2),nfreqs)
        n1=(nfreqs/2)+1
        if(n1.ne.nn) then
		call error(args,2,2)
		go to 99
	endif
	iptsc=0
	if(nargs.eq.4) then
		iptsc=1
		if(ickr(args(3),3,ymin,0,0.).eq.1) go to 99
		if(ickr(args(4),4,ymax,1,ymin).eq.1) go to 99
	endif
c
c   plot :
c
        on=nfreqs
        n1=nstart(np)-1
        do 30 i=1,nn
        wk(i)=float(i-1)/on
  30    wk(nn+i)=array(n1+i)
        if(nargs.eq.2) then
           call min(wk(nn+1),nn,ymin,ind)
           call max(wk(nn+1),nn,ymax,ind)
        endif
	if(nargs.eq.4) then
		xmin=0.
		xmax=.5
                call min(wk(nn+1),nn,ymin1,ind)
                call max(wk(nn+1),nn,ymax1,ind)
		if(ymin1.lt.ymin) ymin=ymin1
		if(ymax1.gt.ymax) ymax=ymax1
	endif
	if(iptpon.eq.1.and.nplt.gt.0) go to 973
c        call smodea(6)
c        call colora(0)
c        call colora(kolor)
 973	continue
        call pgenv(0.0,0.5,ymin,ymax,0,0)
        call pgline(nn,wk,wk(nn+1))
        call pglab('frequency','value',label(np))
c        call gplot(wk(1),wk(nn+1),nn,
c     1  npixx,npixy,nx1,ny1,0,nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,
c     1  iptsc,xmin,xmax,ymin,ymax,0,20,1,1)
	if(iptpon.eq.1) nplt=nplt+1
c        call wrtxt5(label(np),40,nx1+2,ny1+npixy+2,1)
c	call pltend(args,nargs,xmin,xmax,ymin,ymax)
	go to 99
  99    continue
        return
        end
c&INFQNT
        subroutine infqnt(args,nargs)
c******************************************************************
c
c   Subroutine to produce an informative quantile plot.
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki
	ierd=0
c
c
	if(ickl(args(1),np,na).eq.1) go to 99
	if(nargs.eq.1) then
		n=nend(np)-nstart(np)+1
		go to 1701
	endif
	if(icki(args(2),2,n,5,na).eq.1) go to 99
1701	continue
c
c   prepare for plot :
c
	call movxy(wk,array(nstart(np)),4*n)
        call qsort(n,wk)
        i1=(n+1)/2
        i2=(n+2)/2
        xmed=(wk(i1)+wk(i2))/2.
        i1=.25*float(n)
        q1=(wk(i1)+wk(i1+1))/2.
        i1=.75*float(n)
        q3=(wk(i1)+wk(i1+1))/2.
	if(nargs.gt.2) then
		call ckaddr(args(3),xmed,iref)
		if(iref.eq.1) go to 99
		call ckaddr(args(4),q3-q1,iref)
		if(iref.eq.1) go to 99
	endif
	if(nargs.eq.5) go to 99
        tqr=2.*(q3-q1)
	if(tqr.eq.0.) tqr=1.
        on=n
        do 40 i=1,n
        wk(n+i)=(float(i)-.5)/on
        c=(wk(i)-xmed)/tqr
        wk(i)=amax1(c,-1.)
        wk(i)=amin1(wk(i),1.)
  40    continue
c
c   plot :
c
	if(iptpon.eq.1.and.nplt.gt.0) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973	continue
        call loadpg(0.,1.,-1.,1.,0,0)
        call pgenv(0.,1.,-1.,1.,0,0)
        call pgline(n,wk(n+1),wk)
c	call gplot(wk(n+1),wk,n,npixx,npixy,nx1,ny1,0,
c     1  nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,1,0.,1.,-1.,1.,0,20,1,1)
	if(iptpon.eq.1) nplt=nplt+1
	wk(1)=0.
	wk(2)=1.
	wk(3)=0.
	wk(4)=0.
        call pgline(2,wk,wk(3))
c        call gplot(wk,wk(3),2,
c     1  npixx,npixy,nx1,ny1,-1,8,5,5,2,5,2,1,
c     1  0.,1.,-1.,1.,0,20,1,1)
        wk(1)=0.
        wk(2)=1.
        wk(3)=-.5
        wk(4)=.5
        call pgline(2,wk,wk(3))
        wk(3)=-1.
        wk(4)=-1.
        call pgline(2,wk,wk(3))
        wk(3)=1.
        wk(4)=1.
        call pgline(2,wk,wk(3))
        call pglab('u','IQ(u)',label(np))
c        call gplot(wk,wk(3),2,
c     1  npixx,npixy,nx1,ny1,-1,8,5,5,2,5,2,1,
c     1  0.,1.,-1.,1.,0,20,1,1)
c	call line(nx1,ny1+npixy,nx1+npixx,ny1+npixy,1)
c	call line(nx1+npixx,ny1+npixy,nx1+npixx,ny1,1)
c        call wrtxt5(label(np),40,nx1+2,ny1+npixy+2,1)
c	call pltend(args,nargs,0.,1.,-1.,1.)
	go to 99
  99    continue
        return
        end
      subroutine loadpg(xmin,xmax,ymin,ymax)
#include 'tslabc'
      pgxmin=xmin
      pgxmax=xmax
      pgymin=ymin
      pgymax=ymax
      return
      end






