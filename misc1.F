c&MISC1
c&POLAR
        subroutine polar(args,nargs)
c*******************************************************************
c
c   Subroutine to process the command POLAR(xr,xi,n,amp,phase)
c
c   or POLAR(xr,xi,amp,phase)
c
c********************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl
c
c
	pib2=2.*atan(1.0)
c
c   nargs=4:
c
	if(nargs.eq.4) then
		if(ickr(args(1),1,xr,0,0.).eq.1) go to 99
		if(ickr(args(2),2,xi,0,0.).eq.1) go to 99
		amp=sqrt(xr*xr+xi*xi)
		if(xr.eq.0.) then
			phase=pib2
			if(xi.lt.0.) phase=-pib2
			go to 701
		endif
		phase=atan2(xi,xr)
 701		call ckaddr(args(3),amp,iref)
		if(iref.eq.1) go to 99
		call ckaddr(args(4),phase,iref)
		go to 99
	endif
c
c   nargs=5:
c
	if(ickl(args(1),nxr,nsxr).eq.1) go to 99
	if(ickl(args(2),nxi,nsxi).eq.1) go to 99
        call ckint(args(3),n)
        if(n.lt.1.or.n.gt.nsxr.or.n.gt.nsxi) then
                call error(args,3,2)
                go to 99
        endif
c
        nsxr=nstart(nxr)-1
        nsxi=nstart(nxi)-1
        do 10 i=1,n
        c1=array(nsxr+i)
        c2=array(nsxi+i)
        wk(i)=sqrt(c1*c1+c2*c2)
	if(abs(c1).lt.1.e-25) then
		wk(n+i)=pib2
		if(c2.lt.0.) wk(n+i)=-pib2
		go to 10
	endif
        wk(n+i)=atan2(c2,c1)
  10	continue
c
c
        lab='amplitude of a complex array'
        call ckadda(args(4),n,lab,1,iref)
        if(iref.eq.1) go to 99
        lab='phase of a complex array'
        call ckadda(args(5),n,lab,n+1,iref)
c
c
  99    continue
        return
        end
c&PLT801
        subroutine plt801(fname,ifore)
        character fname*20,l*80,l1(80)*1
	equivalence (l,l1)
        open(2,file=fname)
   5    read(2,10,end=39) l
  10    format(a)
	do 15 i=80,1,-1
  15	if(l1(i).ne.' ') go to 17
	i=1
  17    write(*,18) (l1(i),i=1,i)
 18     format(1x,80a1)
        go to 5
  39    continue
        close(2)
        return
        end
        subroutine dos(lun)
c****************************************************************
c
c   subroutine to handle dos commands
c
c*****************************************************************
c
        integer*2 system
        character comm*72,comm1(72)*1,comm4*4,prompt*1
        equivalence (comm,comm1)
        write(*,1)
   1    format(' to return to TIMESLAB prompt, enter EXIT or exit')
  10    call curdsk(ndsk)
        prompt=char(65+ndsk)
        if(lun.eq.100) write(*,2) prompt
   2    format(1x,'?',a1,'>'$)
	if(lun.eq.100) read(*,3) comm
	if(lun.ne.100) read(lun,3) comm
   3	format(a72)
	if(lun.ne.100) write(*,30) prompt,comm
  30	format(1x,'?',a1,'>',a72)
        call movxy(comm4,comm,4)
        call ltoup(comm4,4)
        if(comm4.eq.'EXIT') go to 99
        comm1(72)=char(0)
        i=system(comm)
	if(i.ne.0) then
		write(*,111)
 111	format(' Either can''t find COMMAND.COM or not enough memory')
		go to 99
	endif
        go to 10
c
c
  99    return
        end
c&PRTSCN
        subroutine prtscn
c******************************************************************
c
c
c******************************************************************
c
#include 'tslabc'
c
c
c	n4=ioptls
c        if(n4.eq.1) call grptsc(wk)
c	if(n4.gt.3) call lptsc
c	if(n4.eq.3) call tosh(wk)
c
c
        return
        end
c&POLY
        subroutine poly(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to process the command
c
c   y=POLY(alpha,np,n,a,b)  ,alpha=(A(0),...,A(np))
c
c   or y=POLY(alpha,np,n,x)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),na,nn).eq.1) go to 99
        call ckint(args(2),np)
        if(np.lt.1.or.np+1.gt.nn) then
                call error(args,2,2)
                go to 99
        endif
	if(icki(args(3),3,n,1,1).eq.1) go to 99
c
c
	if(nargs.eq.4) then
		if(ickl(args(4),npx,nx).eq.1) go to 99
		if(nx.lt.n) then
			call error(args,3,2)
			go to 99
		endif
		ns1=nstart(npx)-1
		go to 201
	endif
c
c
	if(ickr(args(4),4,a,0,0.).eq.1) go to 99
	if(ickr(args(5),5,b,0,0.).eq.1) go to 99
        if(a.gt.b) then
                call error(args,5,2)
                go to 99
        endif
c
c
        del=0.0
        if(n.gt.1) del=(b-a)/float(n-1)
 201	continue
        call movxy(wk,array(nstart(na)),4*(np+1))
        npp1=np+1
        do 20 i=1,n
	if(nargs.eq.4) xi=array(ns1+i)
        if(nargs.ne.4) xi=a+float(i-1)*del
                c=0.0
                do 10 j=1,npp1
  10            c=c*xi+wk(npp1-j+1)
  20    wk(npp1+i)=c
c
c
        if(n.gt.1) go to 30
        call ckaddr(vname,wk(npp1+1),iref)
        go to 99
  30    lab='polynomial'
        call ckadda(vname,n,lab,npp1+1,iref)
c
c
  99    continue
        return
        end
c&DOT
        subroutine dot(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to process the command
c
c   c=DOT(x,y,n)
c
c******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
        double precision inprod
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npx,nx).eq.1) go to 99
	if(ickl(args(2),npy,ny).eq.1) go to 99
	nxy=min0(nx,ny)
	if(icki(args(3),3,n,5,nxy).eq.1) go to 99
        c=inprod(array(nstart(npx)),array(nstart(npy)),n)
        call ckaddr(vname,c,iref)
c
c
  99    continue
        return
        end
c&MDEL
	subroutine mdel(args,nargs,vname)
c*****************************************************************
c
c   Subroutine to process the command
c
c     B=MDEL(A,n,m,nr1,...,nrk)
c
c****************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	dimension idel(20)
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(icki(args(2),2,n,2,0).eq.1) go to 99
	if(icki(args(3),3,m,2,0).eq.1) go to 99
        if(na.lt.n*m) then
		call error(args,1,2)
		go to 99
	endif
	do 10 i=1,nargs-3
	call ckint(args(3+i),idel(i))
	if(idel(i).le.0.or.idel(i).gt.n) then
		call error(args,3+i,2)
        	go to 99
		endif
  10	continue
c
c
	n1=nstart(npa)-1
	ii=0
	do 20 k=1,m
	do 30 i=1,n
		do 15 j=1,nargs-3
  15		if(i.eq.idel(j)) go to 30
	ii=ii+1
    	wk(ii)=array(n1+(k-1)*n+i)
  30	continue
  20	continue
	lab='Matrix with rows deleted'
	call ckadda(vname,ii,lab,1,iref)
c
c
  99	continue
	return
	end
c&ERROR
        subroutine error(args,nargs,iopt)
c******************************************************************
c
c
c******************************************************************
#include 'tslabc'
        character args(nargs)*15
	imflag=1
        go to(10,20,30,99,40) iopt
  10    write(*,15) args(nargs)
  15    format(1x,a15,' undefined ')
        go to 99
  20    write(*,25) args(nargs),nargs
  25    format(1x,a15,' illegal value for argument number ',i2)
        go to 99
  30    write(*,35)
  35    format(' range of array outside limits')
	go to 99
  40	write(*,45)
  45	format(' Not enough room on disk to save screen')
  99    continue
        return
        end
c&LISTM
	subroutine listm(args,nargs)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   LISTM(a,n,m[,k])
c
c*******************************************************************
c
#include 'tslabc'
	character args(nargs)*15
	integer*2 ickl,icki
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(icki(args(2),2,n,0,0).eq.1) go to 99
	if(icki(args(3),3,m,0,0).eq.1) go to 99
	k=1
	if(nargs.eq.4) then
		if(icki(args(4),4,k,0,0).eq.1) go to 99
		endif
	nn=n*m*k
	if(nn.gt.na) then
		call error(args,1,2)
		go to 99
		endif
c
c
	write(*,10) label(npa)
	if(iptrec.eq.1) write(3,10,err=701) label(npa)
 701	continue
  10	format(1x,a40)
	n1=nstart(npa)-1
	do 20 kk=1,k
	write(*,15)
  15	format(/)
	if(iptrec.eq.1) write(3,15,err=703)
 703	continue
	k1=(kk-1)*n*m
	do 20 nn=1,n
	if(iptrec.eq.1) then
	     	write(3,30,err=702) (array(n1+k1+(i-1)*n+nn),i=1,m)
 702		continue
	endif
  20	write(*,30) (array(n1+k1+(i-1)*n+nn),i=1,m)
  30	format(1x,6f12.6)
c
c
  99	continue
	return
	end
c&MMULT
	subroutine mmult(args,nargs,vname)
c******************************************************************
c
c   Subroutine to handle the commands
c
c   C=MMULT(A,B,n)  or  C=MMULT(A,B,n,m)  or  C=MMULT(A,B,n,m,k)
c
c      A,B square          C=A'B
c******************************************************************
c
#include 'tslabc'
	character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickse
c
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(ickl(args(2),npb,nb).eq.1) go to 99
	if(icki(args(3),3,n,2,0).eq.1) go to 99
c
c
	if(nargs.eq.3) then
		m=n
		k=n
		go to 20
	endif
	if(nargs.eq.4) then
		if(icki(args(4),4,m,2,0).eq.1) go to 99
		k=m
		go to 20
	endif
	if(nargs.eq.5) then
		if(icki(args(4),4,m,2,0).eq.1) go to 99
		if(icki(args(5),5,k,2,0).eq.1) go to 99
	endif
c
c
	if(na.lt.n*m.or.nb.lt.m*k) then
		write(*,10)
  10		format(' Matrix size wrong in MMULT')
		imflag=1
		go to 99
		endif
c
c
  20	if(nargs.eq.4) then
		if(ickse(m*m).eq.1) go to 99
		call mdot(array(nstart(npa)),array(nstart(npb)),m,m,
     1	n,wk)
		lab='Matrix Product'
		call ckadda(vname,m*m,lab,1,iref)
		go to 99
	endif
c
c
	if(ickse(n*m+n*k).eq.1) go to 99
	ii=0
	n1=nstart(npa)-1
	do 30 i=1,n
	do 30 j=1,m
	ii=ii+1
  30	wk(ii)=array(n1+(j-1)*n+i)
	call mdot(wk,array(nstart(npb)),n,k,m,wk(n*m+1))
	lab='Matrix Product'
	call ckadda(vname,n*k,lab,n*m+1,iref)
c
c
  99	continue
	return
	end
c&MDOT
	subroutine mdot(a,b,ncolsa,ncolsb,nrows,c)
c********************************************************************
c
c   Subroutine to find inner products odf corresponding columns of
c   a and b.
c
c********************************************************************
c
	dimension a(nrows,ncolsa),b(nrows,ncolsb),c(ncolsa,ncolsb)
	double precision inprod
c
c
	do 10 i=1,ncolsa
	do 10 j=1,ncolsb
  10	c(i,j)=inprod(a(1,i),b(1,j),nrows)
c
c
	return
	end
c&MINV
        subroutine minv(args,nargs,vname)
c*******************************************************************
c
c   Subroutine to handle the command
c
c   b=MINV(a,n,ier)
c
c*******************************************************************
c
#include 'tslabc'
        character args(nargs)*15,vname*15
	integer*2 ickl,icki,ickse
c
	if(ickl(args(1),npa,na).eq.1) go to 99
	if(icki(args(2),2,n,0,0).eq.1) go to 99
        if(na.lt.n*n) then
              call error(args,2,2)
              go to 99
        endif
	n1=2*n*n
	if(ickse(n1).eq.1) go to 99
c
c
        call movxy(wk,array(nstart(npa)),4*n*n)
	call decomp(n,n,wk,iwk)
	ier=0
	if(iwk(n).eq.0) then
		ier=1
	        write(*,20)
  20            format(' Singular Matrix in MINV')
		imflag=1
	endif
        call ckaddi(args(3),ier,iref)
	if(ier.eq.1.or.iref.eq.1) go to 99
	n1=n*n
	do 30 i=1,n
		n2=n1+(i-1)*n
		do 25 j=1,n
  25		wk(n2+j)=0.0
		wk(n2+i)=1.0
		call solv(n,n,wk,wk(n2+1),iwk)
  30	continue
c
c
        lab='Inverse Matrix'
        call ckadda(vname,n*n,lab,n1+1,iref)
c
c
  99    continue
        return
        end
