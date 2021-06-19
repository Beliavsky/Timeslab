c&TYPE4
c&TYPE4
        subroutine type4(l)
#include 'tslabc'
        character l(72)*1,vname*15,vn1*15,vn2*15,lck*6,lab40(40)*1
	character vnv*15,vnv15(15)*1
	integer*2 ickse
	equivalence (vnv,vnv15)
        equivalence (lab,lab40)
c
c   find =,<, and > : (we know from typchk theres exactly one of each) :
c
        call findch(l,1,lnb,'=',noceqs,nceqs)
        call findch(l,1,lnb,'<',nocbrl,ncbrl)
        call findch(l,1,lnb,'>',nocbrr,ncbrr)
c
c   form variable name, count commas, and find first comma if is one:
c
        call formnm(l,1,nceqs-1,vname)
        call findch(l,ncbrl,ncbrr,',',noccom,ncloc)
c
c   see if this is to manipulate labels:
c
        call movxy(lck,vname,6)
        call ltoup(lck,6)
        if(lck.ne.'LABEL('.or.l(nceqs-1).ne.')') go to 31
c
c   Yes:
c
        if(l(7).eq.')'.or.noccom.gt.1) then
                write(*,11)
  11            format(' illegal LABEL assignment')
		imflag=1
                go to 99
                endif
        call formnm(l,7,nceqs-2,vname)
        call cklist(vname,xlist,numarr,npv)
        if(npv.eq.-999) then
                write(*,12) vname
  12            format(1x,a15,' doesnt exist')
		imflag=1
                go to 99
                endif
        if(noccom.eq.0) then
                call formnm(l,ncbrl+1,ncbrr-1,vn1)
                call cklist(vn1,xlist,numarr,np)
                if(np.eq.-999) go to 21
                label(npv)=label(np)
                go to 99
  21            write(*,12) vn1
		imflag=1
                go to 99
                endif
        call formnm(l,ncbrl+1,ncloc-1,vn1)
        call formnm(l,ncloc+1,ncbrr-1,vn2)
        call cklist(vn1,xlist,numarr,np1)
        call cklist(vn2,xlist,numarr,np2)
        if(np1.ne.-999.and.np2.ne.-999) then
        call movxy(lab,label(np2),40)
        call movxy(lab40(21),lab40(1),20)
        call movxy(lab,label(np1),20)
        label(npv)=lab
        go to 99
        endif
        write(*,11)
	imflag=1
        go to 99
c
c   Numeric use of TYPE 4:
c
  31    nn=0
        nl=ncbrl+1
        do 4010 i=ncbrl+1,ncbrr
        if(l(i).ne.','.and.l(i).ne.'>') go to 4010
        if(nl.gt.i-1) then
                write(*,4005) vn1
4005    format(1x,a15,' illegal real number in type 4 command')
	imflag=1
                go to 99
                endif
        call formnm(l,nl,i-1,vn1)
	vnv=vn1
	do 201 jj=1,15
 201	if(vnv15(jj).eq.'.') go to 210
	if(ichar(vnv15(1)).ge.48.and.ichar(vnv15(1)).le.57) go to 205
	if(ichar(vnv15(1)).eq.45) go to 205
	go to 210
 205	do 206 jj=15,1,-1
 206	if(vnv15(jj).ne.' ') go to 207
 207	vnv15(jj+1)='.'
	vn1=vnv
 210	continue
        call ckreal(vn1,x)
        if(x.ne.999.) go to 4006
                call cklist(vn1,xlist,numarr,np)
                if(np.eq.-999) then
                write(*,4005) vn1
		imflag=1
                go to 99
                endif
                nleng=nend(np)-nstart(np)+1
		if(ickse(nn+nleng).eq.1) go to 99
                n1=nstart(np)-1
                do 701 j=1,nleng
                nn=nn+1
 701            wk(nn)=array(n1+j)
                nl=i+1
                go to 4010
4006    nn=nn+1
	if(ickse(nn).eq.1) go to 99
        wk(nn)=x
        nl=i+1
4010    continue
  80    continue
        lab='array of real numbers'
c
c   check whether variable name already exists :
c
        call ckadda(vname,nn,lab,1,iref)
c
c   finish :
c
  99    continue
        return
        end
