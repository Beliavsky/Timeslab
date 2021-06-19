c&HPPLOT
c&HPPLOT
        subroutine hpplot(x,y,k,n,nptype,ncols,nrows,
     1  xmin,xmax,ymin,ymax,nticsx,nticsy,nbfsz,buff,fname,
     1  nclabx,ndecx,nclaby,ndecy,labx,laby,cap,fbuff,nptr,nw)
c*********************************************************************
c
c
c
c
c
c*********************************************************************
c
        dimension x(1),y(1),n(1),nptype(1),nptr(1),nw(1)
        character buff(1)*1,fname*15,labx(40)*1,laby(40)*1,
     1   cap(40)*1,fbuff(1)*1
        integer*2 strlen,twidth
        logical logic
c
        nsz=6
c
c
        inquire(file='tslfnt.bin',exist=logic)
c
        if(logic) go to 5
        write(*,2)
   2    format(' Can''t find TSLFNT.BIN')
        return
c
   5    open(201,file='tslfnt.bin',form='binary')
        read(201) (nw(i),i=1,128)
        read(201) (nptr(i),i=1,129)
        read(201) (fbuff(i),i=1,8748)
        close(201)
c
c
        rx=xmax-xmin
        ry=ymax-ymin
c
        if(rx.lt.1.e-15.or.ry.lt.1.e-15) then
           write(*,7)
   7       format(' Range of x or y too small')
           return
        endif
c
        lengx=ncols-200
        lengy=nrows-200
        nx1=150
        ny1=150
c
c
        open(201,file=fname,form='binary',status='new')
        write(201,err=98) char(27),'*','t','3','0','0','R'
        write(201,err=98) char(27),'*','r','1','A'
c
c
        ncols=8*(ncols/8)
        nbtr=ncols/8
        nrp=nbfsz/nbtr
        npasss=((nrows-1)/nrp)+1
        nrtop=0
        nrbot=nrp-1
c
c
        do 100 npass=1,npasss
c
c
          write(*,60) npass,npasss
  60      format(' Pass ',i2,' of ',i2)
          write(*,61) nrtop,nrbot
  61      format(' Top and bottom rows: ',i5,2x,i5)
c
          call movct(buff,20000,char(0))
c
          lx=twidth(labx,40,nw)
          lc=twidth(cap,40,nw)
          lw=ny1+lengy-((lengy-36*strlen(laby,40))/2)-36
c
          call hptext(labx,strlen(labx,40),nx1+(lengx-lx)/2,ny1-100,
     1           buff,nrtop,nrbot,nrows,nbtr,fbuff,nptr,nw)
          call hptext(cap,strlen(cap,40),nx1+(lengx-lc)/2,ny1+lengy+20,
     1            buff,nrtop,nrbot,nrows,nbtr,fbuff,nptr,nw)
c
          do 17 ic=1,strlen(laby,40)
             call hptext(laby(ic),1,nx1-140,lw,buff,nrtop,nrbot,nrows,
     1            nbtr,fbuff,nptr,nw)
  17      lw=lw-36
c
c   ***  LOOP OVER OVERLAYS  ***
c
           noff=0
           do 50 nover=1,k
              ntype=nptype(nover)-3  
c
              if(nover.eq.1) then
c
                 call hpxax(nx1,ny1,lengx,nticsx,xmin,xmax,0,nclabx,
     1                ndecx,nrows,buff,nrtop,nrbot,nbtr,fbuff,nptr,nw)
                 call hpyax(nx1,ny1,lengy,nticsy,ymin,ymax,0,nclaby,
     1                ndecy,nrows,buff,nrtop,nrbot,nbtr,fbuff,nptr,nw)
c
              endif
c
c   ***  LOOP OVER (x,y) PAIRS  ***
c
              do 30 i=1,n(nover)   
c
              ix=nx1+
     1    max0(0,min0(int(.5+lengx*(x(noff+i)-xmin)/rx),lengx-1))
              iy=ny1+
     1    max0(0,min0(int(.5+lengy*(y(noff+i)-ymin)/ry),lengy-1))
c
              if(ntype.eq.10.or.ntype.eq.30) then
c
                 call hpdot(ix,iy,buff,nrtop,nrbot,nrows,nbtr)
                 go to 29
c
              endif
c
              if(ntype.eq.11.or.ntype.eq.31) then
c
                 call hpline(ix-nsz,iy-nsz,ix+nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix-nsz,iy+nsz,ix+nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 go to 29
c
              endif
c
              if(ntype.eq.12.or.ntype.eq.32) then    
c
                 call hpline(ix-nsz,iy+nsz,ix+nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix+nsz,iy+nsz,ix+nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix+nsz,iy-nsz,ix-nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix-nsz,iy-nsz,ix-nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 go to 29
c
              endif 
c
              if(ntype.eq.13.or.ntype.eq.33) then    
c
                 call hpline(ix-nsz,iy+nsz,ix-nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix-nsz,iy-nsz,ix+nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix+nsz,iy-nsz,ix-nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 go to 29
c
              endif
c
              if(ntype.eq.14.or.ntype.eq.34) then    
c
                 call hpline(ix-nsz,iy+nsz,ix+nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix+nsz,iy+nsz,ix+nsz,iy-nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 call hpline(ix+nsz,iy-nsz,ix-nsz,iy+nsz,buff,
     1                       nrtop,nrbot,nrows,nbtr)   
                 go to 29
c
              endif
c
  29          continue
              if(ntype.ge.30) go to 30
              if(i.gt.1) call hpline(ixo,iyo,ix,iy,buff,nrtop,
     1                   nrbot,nrows,nbtr)
              ixo=ix
              iyo=iy
c
c
  30          continue
c
  50       noff=noff+n(nover)
c
c   ***  SEND RESULTS OF THIS PASS TO FILE  ***
c
           do 40 nr=0,nrp-1
              no=nr*nbtr
              do 35 nbytes=nbtr,1,-1
  35             if(buff(nr*nbtr+nbytes).ne.char(0)) go to 37
              nbytes=1
  37          continue
              write(201,err=98) char(27),'*','b'
c
              nb1=nbytes/100
              nb2=(nbytes-100*nb1)/10
              nb3=nbytes-100*nb1-10*nb2
              if(nb1.ne.0) write(201,err=98) char(48+nb1),
     1          char(48+nb2),char(48+nb3)
              if(nb1.eq.0.and.nb2.ne.0) write(201,err=98) char(48+nb2),
     1          char(48+nb3)
              if(nb1.eq.0.and.nb2.eq.0) write(201,err=98) char(48+nb3)
c
              write(201,err=98) 'W'              
c
              do 38 i=1,nbytes
  38          write(201,err=98) buff(no+i)

c
c
  40       continue
        nrtop=nrtop+nrp
        nrbot=min0(nrows-1,nrbot+nrp)
 100    continue
c
        write(201,err=98) char(27),'*','r','B'
        return
c
c
  98    write(*,99) outnm
  99    format(' Disk error writing ',a15)
        close(201)
        return
        end
c&HPXAX
        subroutine hpxax(nx1,ny1,leng,ntics,xmin,xmax,iptax,nclab,
     1             ndec,nrows,buff,nrtop,nrbot,nbtr,fbuff,nptr,nw)
c*********************************************************************
c
c
c
c*********************************************************************
c
        character form*6,form1(6)*1,st*20
        character*1 buff(1),fbuff(1),st1(20)
        dimension nptr(1),nw(1)
        integer*2 twidth
        equivalence (form,form1),(st,st1)
c
        if(iptax.ne.1) then
           form='(f . )'
	         write(form1(3),10) nclab
	         write(form1(5),10) ndec
  10	     format(i1)
        endif
c
        call hpline(nx1,ny1,nx1+leng,ny1,buff,nrtop,nrbot,nrows,nbtr)
c
        rx=xmax-xmin
        do 30 i=0,ntics
c
           ix=nx1+(i*leng)/ntics
c
           call hpline(ix,ny1,ix,ny1+10,buff,nrtop,nrbot,nrows,nbtr)
c
           if(iptax.eq.1) go to 30
           xx=xmin+(i*rx)/ntics
           write(st,form) xx
c
           do 20 j=1,nclab-1
  20          if(st1(j).eq.' '.and.st1(j+1).eq.'.') st1(j)='0'
           do 21 j=nclab,1,-1
           if(st1(j).ne.'0'.and.st1(j).ne.'.') go to 22
           if(st1(j).eq.'.') then
                st1(j)=' '
                go to 22
           endif
  21       if(st1(j).eq.'0') st1(j)=' '
  22       continue
c
           ncc=ix-twidth(st,nclab,nw)/2
           nrr=ny1-50
           call hptext(st,nclab,ncc,nrr,buff,nrtop,nrbot,nrows,
     1                 nbtr,fbuff,nptr,nw)
c
  30    continue
c
c
        return
        end
c&HPYAX
        subroutine hpyax(nx1,ny1,leng,ntics,ymin,ymax,iptax,nclab,
     1             ndec,nrows,buff,nrtop,nrbot,nbtr,fbuff,nptr,nw)
c*********************************************************************
c
c
c
c*********************************************************************
c
        character form*6,form1(6)*1,st*20
        character*1 buff(1),fbuff(1),st1(20)
        dimension nptr(1),nw(1)
        integer*2 twidth
        equivalence (form,form1),(st,st1)
c
        if(iptax.ne.1) then
           form='(f . )'
	         write(form1(3),10) nclab
	         write(form1(5),10) ndec
  10	     format(i1)
        endif
c
        call hpline(nx1,ny1,nx1,ny1+leng,buff,nrtop,nrbot,nrows,nbtr)
c
        ry=ymax-ymin
        do 30 i=0,ntics
c
           iy=ny1+(i*leng)/ntics
c
           call hpline(nx1,iy,nx1+10,iy,buff,nrtop,nrbot,nrows,nbtr)
c
           if(iptax.eq.1) go to 30
           yy=ymin+(i*ry)/ntics
           write(st,form) yy
c
           do 20 j=1,nclab-1
  20          if(st1(j).eq.' '.and.st1(j+1).eq.'.') st1(j)='0'
           do 21 j=nclab,1,-1
           if(st1(j).ne.'0'.and.st1(j).ne.'.') go to 22
           if(st1(j).eq.'.') then
                st1(j)=' '
                go to 22
           endif
  21       if(st1(j).eq.'0'.or.st1(j).eq.'.') st1(j)=' '
  22       continue
c
           ncc=nx1-twidth(st,nclab,nw)-2
           nrr=iy-18
           call hptext(st,nclab,ncc,nrr,buff,nrtop,nrbot,nrows,
     1                 nbtr,fbuff,nptr,nw)
c
  30    continue
c
c
        return
        end
c&TWIDTH
        integer*2 function twidth(lab,maxt,nw)
c*********************************************************************
c
c
c
c*********************************************************************
c
        dimension nptr(1),nw(1)
        character*1 lab(1)
        integer*2 strlen
c
        twidth=0
        do 10 i=1,strlen(lab,maxt)
  10    twidth=twidth+nw(ichar(lab(i))+1)
c
        return
        end
c&HPTEXT
        subroutine hptext(text,maxt,ncol,nrow,buff,nrtop,nrbot,
     1    nrows,nbtr,fbuff,nptr,nw)
c*********************************************************************
c
c
c
c********************************************************************
c
        character*1 buff(1),fbuff(1),text(1)
        integer*2 strlen,ibtor,ibtand,getfnt,nptr(1),nw(1)
c
        nr=nrows-36-nrow
c
        do 100 i=0,35
        nc=ncol
c
        if(nr.lt.nrtop.or.nr.gt.nrbot) go to 100
c
           do 50 j=1,strlen(text,maxt)
c
           nch=ichar(text(j))
           ncw=nw(nch+1)
           ibyte=getfnt(nch,0,0,nbr,nptr,fbuff)
c
              do 30 jj=0,nbr-1
                 ibyte=getfnt(nch,i,jj,nbr,nptr,fbuff)
                 kk=min0(8,ncw-8*jj)-1
c
                 do 20 k=0,kk
                    if(ibtand(ibyte,nshftr(128,k)).ne.0) then
                       indbuf=(nr-nrtop)*nbtr+((nc+k)/8)+1
                       jbyte=ichar(buff(indbuf))
                       jbyte=ibtor(jbyte,nshftr(128,mod(nc+k,8)))
                       buff(indbuf)=char(jbyte) 
                    endif
  20             continue      
c
  30          nc=nc+kk+1
c
  50       continue
c
 100    nr=nr+1
c
        return
        end
c&GETFNT
        integer*2 function getfnt(nchar,nrow,ncol,nbr,nptr,fbuff)
c*********************************************************************
c
c
c
c*********************************************************************
c
        dimension nptr(129)
        character*1 fbuff(1)
c
c
        nbr=(nptr(nchar+2)-nptr(nchar+1))/36
        indbf=nptr(nchar+1)+nbr*nrow+ncol+1
        getfnt=ichar(fbuff(indbf))
c
        return
        end
c&STRLEN
        integer*2 function strlen(str,maxl)
c*********************************************************************
c
c
c
c*********************************************************************
c
        character*1 str(1)
c
        do 10 strlen=maxl,1,-1
  10    if(str(strlen).ne.' ') return
        strlen=1
c
        return
        end
c&HPLINE
        subroutine hpline(ncol1,nrow1,ncol2,nrow2,buff,nrtop,nrbot,
     1                    nrows,nbtr)
c*********************************************************************
c
c
c
c********************************************************************
c
        character*1 buff(1)
        nrmin=min0(nrow1,nrow2)
        nrmax=max0(nrow1,nrow2)
        if((nrmax.lt.nrows-1-nrbot).or.(nrmin.gt.nrows-1-nrtop)) return
c
c
        if(nrow1.eq.nrow2) then
          do 10 i=min0(ncol1,ncol2),max(ncol1,ncol2)
  10         call hpdot(i,nrow1,buff,nrtop,nrbot,nrows,nbtr)
          return
        endif
c
c
        if(ncol1.eq.ncol2) then
                do 20 i=min0(nrow1,nrow2),max0(nrow1,nrow2)
  20                    call hpdot(ncol1,i,buff,nrtop,nrbot,nrows,nbtr)
                return
        endif
c
c
        if(ncol1.le.ncol2) then
                nc0=ncol1
                nc1=ncol2
                nr0=nrow1
                nr1=nrow2
        else
                nc0=ncol2
                nc1=ncol1
                nr0=nrow2
                nr1=nrow1
        endif
c
c
        ndc=nc1-nc0
        ndr=nr1-nr0
        nadr=iabs(ndr)
        nrinc=nadr/ndr
        npix=max0(ndc,nadr)+1
c
c                                       ***  SHALLOW  ***
c
        if(ndc.ge.nadr) then
                incd=2*(ndc-nadr)
                incn=-2*nadr
                ia=ndc-2*nadr
                do 30 i=1,npix
                        call hpdot(nc0,nr0,buff,nrtop,nrbot,nrows,nbtr)
                        nc0=nc0+1
                        if(ia.gt.0) then
                                ia=ia+incn
                        else
                                nr0=nr0+nrinc
                                ia=ia+incd
                        endif
  30            continue
                return
c
c                                       ***  STEEP  ***
c
        else
                incd=2*(nadr-ndc)
                incn=-2*ndc
                ia=nadr-2*ndc
                do 40 i=1,npix
                        call hpdot(nc0,nr0,buff,nrtop,nrbot,nrows,nbtr)
                        nr0=nr0+nrinc
                        if(ia.gt.0) then
                                ia=ia+incn
                        else
                                nc0=nc0+1
                                ia=ia+incd
                        endif
  40            continue
                return
        endif
c
c
        end
c&HPDOT
        subroutine hpdot(ncol,nrow,buff,nrtop,nrbot,nrows,nbtr)
c**********************************************************************
c
c
c
c*********************************************************************
c
        character*1 buff(1)
        integer*2 ibtor,nshftr
c
c
        nrow1=nrows-1-nrow
        if(nrow1.lt.nrtop.or.nrow1.gt.nrbot) return
c
c
        nrow1=nrow1-nrtop
        indbuf=nrow1*nbtr+(ncol/8)+1
        ibyte=ichar(buff(indbuf))
c
        buff(indbuf)=char(ibtor(ibyte,nshftr(128,mod(ncol,8))))
c
c
        return
        end       
c&IBTOR
        integer*2 function ibtor(n,m)
c**********************************************************************
c
c
c
c*********************************************************************
c
        ibtor=0
        nn=n
        mm=m
        nck=128
        do 10 i=1,8
                nm=0
                if(nn.ge.nck) then
                        nm=1
                        nn=nn-nck
                endif
                if(mm.ge.nck) then
                        nm=1
                        mm=mm-nck
                endif
                ibtor=ibtor+nm*nck
                nck=nck/2
  10    continue
c
c
        return
        end
c&IBTAND
        integer*2 function ibtand(n,m)
c**********************************************************************
c
c
c
c*********************************************************************
c
        ibtand=0
        nn=n
        mm=m
        nck=128
        do 10 i=1,8
                nm1=0
                nm2=0
                if(nn.ge.nck) then
                        nm1=1
                        nn=nn-nck
                endif
                if(mm.ge.nck) then
                        nm2=1
                        mm=mm-nck
                endif
                ibtand=ibtand+nm1*nm2*nck
                nck=nck/2
  10    continue
c
c
        return
        end
c&NSHFTR
        integer*2 function nshftr(n,ns)
c**********************************************************************
c
c
c
c*********************************************************************
c
        nshftr=n
        if(ns.eq.0) return
c
        do 10 i=1,ns
  10    nshftr=nshftr/2
c
c
        return
        end
