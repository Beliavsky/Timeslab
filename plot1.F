c&PLOT1
c&GPLOT
      SUBROUTINE GPLOT(X,Y,N,NPIXX,NPIXY,NX1,NY1,IPTAX,NTICSX,NTICSY,
     1NCLABX,NDECX,NCLABY,NDECY,IPTSC,XMIN,XMAX,YMIN,YMAX,IPTEQ,IPTTYP,
     1KOLOR1,KOLOR2)
c******************************************************************
c
c   
c        X,Y,N : vetors to be plotted and how many elements there are.
c
c        NPIXX,NPIXY : size of plotting region
c
c        NX1,NY1 :  lower left hand pixel of plotting region (including
c                   axes if drawn.
c
c        NTICSX,NTICSY : number of tic marks on the x and y axes
c                        respectively (ignored if iptax=-1)
c
c        NCLABX,NDECX  : labels on x axis are in F NCLABX.NDECX format
c        NCLABY,NDECY  : labels on y axis are in F NCLABY.NDECY format
c                        (all 4 of these arguments are ignored if
c                        IPTAX=-1 or +1)
c
c        IPTAX : -1: no axes, 1:axes, no label, 0:axes and label
c        IPTSC :  1: mins and maxs input, 0 they're calculated
c        IPTEQ :  1: if n=npixx and x's equally spaced then 1 pixel
c		     per column is turned on
c        IPTTYP : 20: line plot, 34:point plot(3x3 x), 14:points+lines
c
c        KOLOR1,KOLOR2: color of axes and lines, color of points
c
c        XMIN,XMAX  : values of X at screen locations (NX1,NY1) and
c                     (NX1+NPIXX,NY1) respectively, i.e. values at the
c                     ends of the graph in the horizontal direction.
c        YMIN,YMAX  : values of Y at screen locations (NX1,NY1) and
c                     (NX1,NY1+NPIXY).
c
c
c******************************************************************
c
      dimension x(n),y(n)
c
c   find max and min if iptsc=0:
c
	if(iptsc.eq.0) then
        	call max(x,n,xmax,ind)
        	call min(x,n,xmin,ind)
        	call max(y,n,ymax,ind)
        	call min(y,n,ymin,ind)
	endif
        rx=xmax-xmin
        ry=ymax-ymin
        if(rx.lt.1.e-15.or.ry.lt.1.e-15) then
        	call wrtxt1('Range of X or Y is Zero',23,0,192,1)
        	go to 99
	endif
c
c   draw axes if iptax.ne.-1 :
c
        if(iptax.ne.-1) then
		call xaxis(nx1,ny1,npixx,nticsx,kolor1,xmin,xmax,iptax,
     1  	nclabx,ndecx)
        	call yaxis(nx1,ny1,npixy,nticsy,kolor1,ymin,ymax,iptax,
     1  	nclaby,ndecy)
	endif
c
c   find screen location of x,y :
c
        onpixx=npixx
        onpixy=npixy
        do 30 i=1,n
        i2=1+(onpixy*(y(i)-ymin)/ry)
        iy=ny1+min0(npixy,i2)
        i1=1+(onpixx*(x(i)-xmin)/rx)
        ix=nx1+min0(npixx,i1)
        if(ipteq.eq.1.and.n.eq.npixx) ix=nx1+i
c
c   draw points unless ipttyp=20 :
c
        if(ipttyp.eq.20) go to 90
		if(ipttyp.eq.13.or.ipttyp.eq.33) then
		        call line(ix,iy,ix,iy,kolor2)
			go to 90
		endif
		if(ipttyp.eq.14.or.ipttyp.eq.34) then
		        call line(ix-1,iy-1,ix+1,iy+1,kolor2)
		        call line(ix-1,iy+1,ix+1,iy-1,kolor2)
			go to 90
		endif
		if(ipttyp.eq.15.or.ipttyp.eq.35) then
			call line(ix-1,iy+1,ix+1,iy+1,kolor2)
		        call line(ix+1,iy+1,ix+1,iy-1,kolor2)
			call line(ix+1,iy-1,ix-1,iy-1,kolor2)
		        call line(ix-1,iy-1,ix-1,iy+1,kolor2)
			go to 90
		endif
		if(ipttyp.eq.16.or.ipttyp.eq.36) then
			call line(ix-1,iy+1,ix-1,iy-1,kolor2)
		        call line(ix-1,iy-1,ix+1,iy-1,kolor2)
			call line(ix+1,iy-1,ix-1,iy+1,kolor2)
			go to 90
		endif
		if(ipttyp.eq.17.or.ipttyp.eq.37) then
			call line(ix-1,iy+1,ix+1,iy+1,kolor2)
		        call line(ix+1,iy+1,ix+1,iy-1,kolor2)
			call line(ix+1,iy-1,ix-1,iy+1,kolor2)
			go to 90
		endif
  90    continue
c
c   draw line plot unless ipttyp.gt.30 :
c
        if(ipttyp.gt.30) go to 30
        if(i.gt.1) call line(ixo,iyo,ix,iy,kolor1)
        ixo=ix
        iyo=iy
  30    continue
  99    continue
c
c   finish :
c
      return
      end
c&XAXIS
      SUBROUTINE XAXIS(NX1,NY1,NPIXS,NTICS,KOLOR,XMIN,XMAX,IPTLB,
     1NCLABX,NDECX)
c*****************************************************************
c
c   Subroutine to draw a horizontal axis from screen location
c   (NX1,NY1) to location (NX1+NPIXS,NY1) with NTICS equally spaced
c   tic marks (they are at locations NX1+j*INT(NPIXS/NTICS),
c   j=1,...,NTICS) with numerical labels (if IPTLB=0) placed at any
c   tic marks where they will fit. Between tic marks there are
c   NLENG=INT(NPIXS/NTICS) pixels. The labels are centered on the
c   tic marks in F NCLABX.NDECX format. Thus NTOT=8*NCLABX pixels
c   are required for each label.  If ntics is negative, then the
c   tic marks are below the axis.
c
c      If NTOT<NLENG, then all tic marks are labelled.
c      If NLENG.le.NTOT.lt.2*NLENG, the 2nd,4th,etc. are labelled.
c      IF 2*NLENG.le.NTOT.lt.3*NLENG, the 3rd,6th,etc. are labelled.
c      etc.
c
c   XMIN is the value of X corresponding to location (NX1,NY1).
c   XMAX is the value of X corresponding to location (NX1+NPIXS,NY1).
c
c   The tic marks extend 3 pixels above the axis.
c
c   The axis and labels are in color determined by KOLOR. (In
c   medium resolution this is 1,2,3 for green,red,yellow for
c   palette 0 or cyan,magenta,white for palette 1 . In high
c   resolution kolor is ignored and the axis and labels are in the
c   color specified in the calling program as the background color).
c
c********************************************************************
c
	character form*6,form1(6)*1,xchar*20,xchar1(20)*1
	equivalence (xchar,xchar1),(form,form1)
c
c   draw axis and tic marks:
c
	call line(nx1,ny1,nx1+npixs,ny1,kolor)
        ntics1=iabs(ntics)
        ipttic=3
        if(ntics.lt.0) ipttic=-3
	nbetw=npixs/ntics1
	do 10 i=1,ntics1
        ii=i*nbetw
        call line(nx1+ii,ny1,nx1+ii,ny1+ipttic,kolor)
  10    continue
c
c   draw labels if iptlb=0 :
c
        if(iptlb.ne.0.or.nclabx.eq.0) go to 99
	form='(f . )'
	write(form1(3),20,err=99) nclabx
	write(form1(5),20,err=99) ndecx
  20	format(i1)
c	write(form,20,err=21) nclabx,ndecx
c 20  	format('(f',i1,'.',i1,')')
c	go to 23
c 21	write(*,22) nclabx,ndecx,form
c 22	format(/,' nclabx = ',i5,' ndecx = ',i5,' form = ',a6)
c	go to 99
c 23	continue
c
c   set constants :
c
        nrow=ny1-12
        if(ntics.lt.0) nrow=ny1+4
        nleft=8*(nclabx-ndecx-1)+4
        nleng=8*nclabx
        ileng=0
c
c   fix tic mark number :
c
	rx=xmax-xmin
	x1=rx*float(nbetw)/float(npixs)
c
c
        do 30 i=0,ntics1
        ileng=ileng+nbetw
        if(i.ne.0.and.nleng.gt.ileng) go to 30
c
c   label fits so write it :
c
        ileng=0
	ncent=nx1+i*nbetw
      	ncol=ncent-nleft
	x=xmin+x1*float(i)
	write(xchar,form,err=99) x
c
c   put a 0 in if nothing to left of decimal point:
c
	do 27 j=1,nclabx-1
  27    if(xchar1(j).eq.' '.and.xchar1(j+1).eq.'.') xchar1(j)='0'
        call wrtxt1(xchar,nclabx,ncol,nrow,kolor)
c
c
  30    continue
c
c
  99	continue
        return
        end
c&YAXIS
      SUBROUTINE YAXIS(NX1,NY1,NPIXS,NTICS,KOLOR,YMIN,YMAX,IPTLB,
     1NCLABY,NDECY)
c*******************************************************************
c
c
c********************************************************************
c
	character form*6,form1(6)*1,ychar*20,ychar1(20)*1
	equivalence (ychar,ychar1),(form,form1)
c
c   draw axis :
c
      call line(nx1,ny1,nx1,ny1+npixs,kolor)
c
c   draw tic marks :
c
        ntics1=iabs(ntics)
        ipttic=3
        if(ntics.lt.0) ipttic=-3
      nbetw=npixs/ntics1
      do 10 i=1,ntics1
      ii=i*nbetw
      call line(nx1,ny1+ii,nx1+ipttic,ny1+ii,kolor)
  10  continue
c
c
c   draw labels if iptlb=0
c
c
      if(iptlb.ne.0.or.nclaby.eq.0) go to 99
c
c   form format statement :
c
	form='(f . )'
	write(form1(3),20,err=99) nclaby
	write(form1(5),20,err=99) ndecy
  20	format(i1)
c     write(form,20,err=99) nclaby,ndecy
c 20  format('(f',i1,'.',i1,')')
c
c   set constants :
c
      ncol=nx1-8*nclaby-2
        if(ntics.lt.0) ncol=nx1+2
      nleng=8
      ileng=0
c
c   fix tic mark number :
c
      do 30 i=0,ntics1
      ileng=ileng+nbetw
      if(i.ne.1.and.nleng.gt.ileng) go to 30
c
c   label fits, write it :
c
      ileng=0
      ncent=ny1+i*nbetw
      nrow=ncent-4
      y=ymin+(float(i*nbetw)/float(npixs))*(ymax-ymin)
      write(ychar,form,err=99) y
	do 28 j=1,nclaby-1
  28	if(ychar1(j).eq.' '.and.ychar1(j+1).eq.'.') ychar1(j)='0'
      call wrtxt1(ychar,nclaby,ncol,nrow,kolor)
  30  continue
c
c
  99  continue
      return
      end
c&PLOTK
	subroutine plotk(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   PLOTK(x,y,n,k,type[,xmin,xmax,ymin,ymax])
c
c*******************************************************************
#include 'tslabc'
	character args(nargs)*15,outnm*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(ickl(args(2),npy,ny).eq.1) go to 99
	if(icki(args(4),4,k,1,1).eq.1) go to 99
	call ckint(args(3),n)
	if(n.ne.-999) then
		if(n.lt.1.or.n*k.gt.min0(nx,ny)) then
			call error(args,3,2)
			go to 99
		endif
		do 1 i=1,k+1
   1		iwk(i)=(i-1)*n+1
		nk=n*k
	endif
	if(n.eq.-999) then
		call cklist(args(3),xlist,numarr,npn)
		if(npn.eq.-999) then
			call error(args,3,2)
			go to 99
		endif
		nn=nend(npn)-nstart(npn)+1
		if(nn.lt.k) then
			call error(args,4,2)
			go to 99
		endif
		nk=1
		ns1=nstart(npn)-1
		do 3 i=1,k
   		iwk(i)=nk
		nk=nk+array(ns1+i)
   3		continue
		iwk(k+1)=nk
		nk=nk-1
		if(nk.gt.min0(nx,ny)) then
			call error(args,3,2)
			go to 99
		endif
	endif
	if(k.eq.1) then
		call ckint(args(5),ntyp)
		if(ntyp.ne.-999) then
			iwk(1000+1)=ntyp
			go to 7
		endif
	endif
	if(ickl(args(5),nptyp,ntyp).eq.1) go to 99
	if(ntyp.lt.k) then
		call error(args,5,2)
		go to 99
	endif
c
c
	ns1=nstart(nptyp)-1
	do 6 i=1,k
   6	iwk(1000+i)=array(ns1+i)
   7	do 10 i=1,k
	nti=iwk(1000+i)
	if(nti.eq.2) then
		iwk(1000+i)=20
		go to 10
	endif
	if(nti.ge.10.and.nti.le.14) then
		iwk(1000+i)=nti+3
		go to 10
	endif
	if(nti.ge.30.and.nti.le.34) then
		iwk(1000+i)=nti+3
		go to 10
	endif
	call error(args,5,2)
	go to 99
c
c
  10	continue
c
c
	call max(array(nstart(npx)),nk,xmax,imax)
	call min(array(nstart(npx)),nk,xmin,imin)
	call max(array(nstart(npy)),nk,ymax,imax)
	call min(array(nstart(npy)),nk,ymin,imin)
c
c
	if(nargs.eq.9.or.nargs.eq.12) then
		if(ickr(args(6),6,xmin1,0,0.).eq.1) go to 99
		if(xmin1.le.xmin) xmin=xmin1
		if(ickr(args(7),7,xmax1,1,xmin1).eq.1) go to 99
		if(xmax1.gt.xmax) xmax=xmax1
		if(ickr(args(8),8,ymin1,0,0.).eq.1) go to 99
		if(ymin1.lt.ymin) ymin=ymin1
		if(ickr(args(9),9,ymax1,1,ymin1).eq.1) go to 99
		if(ymax1.gt.ymax) ymax=ymax1
	endif
c
c
        if(nargs.eq.12) then
                if(icki(args(10),10,nchp,1,300).eq.1) go to 99
                if(icki(args(11),11,nrhp,1,300).eq.1) go to 99
                if(nchp.gt.2700.or.nrhp.gt.2700) then
                        write(*,30)
  30                    format(' Too many rows or columns')
                        imflag=1
                        go to 99
                endif
c
                outnm=args(12)
                call cklist(args(12),chlist,nchars,np)
                if(np.ne.-999) outnm=chars(np)
c
                do 50 i=1,k
  50            iwk(i)=iwk(i+1)-iwk(i)
                call hpplot(array(nstart(npx)),array(nstart(npy)),
     1                      k,iwk,iwk(1001),nchp,nrhp,xmin,xmax,
     1                      ymin,ymax,nticsx,nticsy,20000,wk,outnm,
     1                      nclabx,ndecx,nclaby,ndecy,
     1                      label(npx),label(nptyp),label(npy),
     1                      wk(5001),iwk(1500),iwk(1700))
                go to 99
        endif
c
c
c
	if(iptpon.eq.1) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973	continue
	nsx1=nstart(npx)-1
	nsy1=nstart(npy)-1
	call wrtxt5(label(npy),40,nx1+2,ny1+npixy+10,1)
	lab='    '
	if(label(npx).eq.lab) go to 497
        call wrtxt5(label(npx),40,nx1+6,ny1+npixy+2,1)
 497	continue
c
c
        call loadpg(xmin,xmax,ymin,ymax)
        call pgenv(xmin,xmax,ymin,ymax,0,0)
	do 20 i=1,k
	   nsx=nsx1+iwk(i)
	   nsy=nsy1+iwk(i)
           ityp=iwk(1000+i)
           nn=iwk(i+1)-iwk(i)
           if(ityp.ne.20) call pgpt(nn,array(nsx),array(nsy),-3)
           if(ityp.le.30) call pgline(nn,array(nsx),array(nsy))
c        call gplot(array(nsx),array(nsy),iwk(i+1)-iwk(i),
c     1  npixx,npixy,nx1,ny1,
c     1	0,nticsx,nticsy,nclabx,ndecx,nclaby,ndecy,1,
c     1  xmin,xmax,ymin,ymax,0,iwk(1000+i),1,1)
  20	continue
        call pglab(' ',label(npy),label(npx))
c
c
c	call pltend(args,nargs,xmin,xmax,ymin,ymax)
  99	continue
	return
	end
