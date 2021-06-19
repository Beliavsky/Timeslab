c&UTIL
c&ADDARR
        subroutine addarr(n,aname,lab,m)
c*********************************************************************
c
c   Subroutine to add an array to the vector array. The elements to
c   be added are to be in wk(m),...,wk(m+n-1) (passed in common).
c
c**********************************************************************
c
#include 'tslabc'
        character aname*15
        num=num-n
        numarr=numarr+1
        xlist(numarr)=aname
        label(numarr)=lab
        nstart(numarr)=1
        if(numarr.gt.1) nstart(numarr)=nend(numarr-1)+1
        nend(numarr)=nstart(numarr)+n-1
        nn=m - nstart(numarr)
        do 50 i=nstart(numarr),nend(numarr)
  50    array(i)=wk(nn+i)
c
c
  99    continue
        return
        end
c&CLEAN
        subroutine clean(narr)
c********************************************************************
c
c   Subroutine to delete the array having index narr. If narr=0, then
c   the user is given the opportunity to delete any or all of the
c   arrays.
c
c*********************************************************************
c
#include 'tslabc'
        character*1 i10
c
c
        numr=numarr
        do 100 i=numr,1,-1
c
c
        if ( narr .eq. 0 ) goto 5
        if ( narr .eq. i ) goto 40
        goto 100
c
c   we know here that narr = 0:
c
   5    write(*,10) xlist(i),nend(i)-nstart(i)+1
  10    format(' do you want to delete array ',a15,' having ',
     1  i5,' elements?(y/n/q)'$)
        read(*,20) i10
  20    format(a1)
        if(i10.eq.'Q'.or.i10.eq.'q') go to 99
        if(i10.eq.'N'.or.i10.eq.'n') go to 100
        if(i10.eq.'Y'.or.i10.eq.'y') go to 40
        write(*,30)
  30    format(' improper response')
        go to 5
  40    continue
c
c   shift everything (unless we're deleting the last array) :
c
        nn=nend(i)-nstart(i)+1
        numarr=numarr-1
        num=num+nn
        if(narr.ne.0.and.i.gt.numarr) go to 99
        if(i.gt.numarr) go to 100
        do 50 j=i,numarr
        xlist(j)=xlist(j+1)
  50    label(j)=label(j+1)
        do 55 j=i,numarr
        ii=1
        if(j.gt.1) ii=nend(j-1)+1
        nstart(j)=ii
        nend(j)=ii+nend(j+1)-nstart(j+1)
  55    continue
c       do 60 j=i,numarr
c       ns=nstart(j)-1
c       n=nend(j)-nstart(j)+1
c       do 65 k=1,n
c 65    array(ns+k)=array(ns+k+nn)
c 60    continue
        do 60 k = nstart(i),nend(numarr)
  60    array(k) = array(k+nn)
        if(narr.ne.0) go to 99
 100    continue
c
c
  99    continue
        return
        end
c&FINDCH
      subroutine findch(l,n1,n2,char,nocc,ncol)
c********************************************************
c
c   subroutine to find the first occurrence ncol of the 
c   single character char in l(n1),...,l(n2), and the
c   number of occurrences nocc.
c
c*********************************************************
c
      character*1 char,l(n2)
c
      nocc=0
      ncol=0
      do 10 i=n1,n2
      if(l(i).ne.char) go to 10
      nocc=nocc+1
      if(nocc.eq.1) ncol=i
  10  continue
      return
      end
c&FORMNM
      subroutine formnm(l,n1,n2,name)
c**********************************************************
c
c   subroutine to form the character*20 variable name from
c   l(n1),...,l(n2).
c
c***********************************************************
c
      character*1 l(n2)
      character*15 name
        if(n2-n1+1.gt.15.or.n2-n1+1.lt.1) go to 5
      write(name,2) (l(i),i=n1,n2)
   2  format(20a1)
        go to 10
   5    name='***************'
  10    continue
      return
      end
c&LASTNB
        subroutine lastnb(l,n,m)
c********************************************************************
c
c   Subroutine to find the number (m) of the last nonblank character in t
c   the character*1 array l(1),...,l(n).
c
c*********************************************************************
c
        character*1 l(n),apst
	apst=char(39)
c
c   Let m=position of last nonblank character.
c   If last non blank character is apostrophe, rturn m. Otherwise
c   continue on looking for a ;. If there isn't one, return m. If
c   there is, and it's in column one, return m. If there is and it's
c   not in column one, blank out from ; to m and return m=position
c   to left of ;.
c
        do 10 i=n,1,-1
        if(l(i).ne.' ') go to 20
  10    continue
  20    m=i
	if(l(m).eq.apst) go to 99
	do 30 i=m,1,-1
	if(l(i).eq.';') go to 40
  30	continue
	go to 99
  40	if(i.eq.1.or.i.eq.2) go to 99
	do 50 j=i,m
  50	l(j)=' '
	m=i-1
  99	continue
        return
        end
c&LTOUP
        subroutine ltoup(l,n)
c********************************************************************
c
c   Subroutine to convert all characters having ASCII codes between 97
c   and 122 to the corresponding character having code 32 less.
c
c*********************************************************************
c
        character*1 l(n)
        do 10 i=1,n
        ii=ichar(l(i))
  10    if(97.le.ii.and.122.ge.ii) l(i)=char(ii-32)
        return
        end
c&PARSE
      subroutine parse(l,pname,args,nargs,vname,ierr)
c*******************************************************************
c
c   subroutine to find the procedure name (pname), arguments
c   (args), and number of arguments (nargs) for a type 1 entry line
c   l(1),...,l(72) . The error flag ierr is
c   0 if no error is encountered. If the line is type 3 then the
c   variable name is also returned in vname.
c
c*******************************************************************
c
      character*1 l(72),comma,equal
      character*15 pname,args(21),vname
      dimension nleft(21),nright(21)
      data comma,equal/',','='/
      ierr=0
c
c   find left and right parentheses :
c
      call findch(l,1,72,'(',noc1,ncol1)
      call findch(l,1,72,')',noc2,ncol2)
c
c   find where arguments are :
c
      nargs=1
      nleft(1)=ncol1+1
      do 20 i=ncol1,ncol2
      if(l(i).ne.comma) go to 20
      nright(nargs)=i-1
      nargs=nargs+1
      nleft(nargs)=i+1
  20  continue
      nright(nargs)=ncol2-1
c
c   form vector of arguments :
c
      do 30 i=1,nargs
	nrmnl=nright(i)-nleft(i)+1
	if(nrmnl.gt.15.or.nrmnl.lt.1) then
		ierr=1
		go to 99
	endif
  30  call formnm(l,nleft(i),nright(i),args(i))
c
c   find procedure name :
c
      nn=ncol1-1
  40  if(nn.eq.0) go to 50
        if(l(nn).eq.equal) go to 50
      nn=nn-1
      go to 40
  50  continue
	nrmnl=ncol1-1-nn
	if(nrmnl.lt.1.or.nrmnl.gt.15) then
		ierr=1
		go to 99
	endif
      call formnm(l,nn+1,ncol1-1,pname)
c
c   see if type 3 :
c
	if(nn.eq.0) go to 99
        if(nn-1.lt.1.or.nn-1.gt.15) then
		ierr=1
		go to 99
	endif
        call formnm(l,1,nn-1,vname)
  99  return
      end
c&RMBLNK
        subroutine rmblnk(l,n,lnb)
c*********************************************************************
c
c   Subroutine to remove any embedded blanks in the character*1 array
c   l(1),...,l(n). The output integer lnb is the number of the last
c   nonblank element in the new version of the array. The array is
c   assumed to have at least one nonblank character.
c
c   Blanks are not removed from between apostrophes.
c
c*********************************************************************
c
        character*1 l(n),blnk,l1(72)*1,apost*1
        data blnk,apost/' ',''''/
        call lastnb(l,n,lnb)
        j=0
        iap=0
        do 10 i=1,n
        if(l(i).eq.apost) iap=mod(iap+1,2)
        if(l(i).eq.blnk.and.iap.eq.0) go to 10
        j=j+1
        l1(j)=l(i)
  10    continue
        call movct(l,lnb,blnk)
        lnb=j
        call movxy(l,l1,lnb)
c
        return
        end
c&TYPCHK
        subroutine typchk(l,lnb,itype)
c*******************************************************************
c
c   Subroutine to determine the command type of the command line made
c   up of the character*1 array l(1),...,l(lnb).
c   It is assumed that there are no blanks in l and that lnb>1.
c
c   The output integer itype is 1,2,3,4,or 5 according to :
c
c      TYPE 1 : x...x(x....x)        i.e. no equal sign, exactly one (,
c                                    with l(lnb)=) and location of (
c                                    not the first element.
c
c      TYPE 2 : x...x=x...x          i.e. exactly 1 equal sign not at
c                                    beginning or end and no (,),<,>
c
c      TYPE 3 : x...x=x...x(x...x)   i.e. one equal sign not at begin
c                                    or end, exactly one ( not next to
c                                    + and not at end, and exactly one
c                                    which is at end. Also no <.
c
c      TYPE 4 : x...x=<x...>         i.e. exactly one = not at begin or
c                                    end, exactly one < next to = and
c                                    exactly one > at the end.
c
c      TYPE 5 : x......x             i.e. no =,(,),<,>
c
c   The integer itype is 0 if l fits none of these types.
c
c**********************************************************************
c
        character*1 l(lnb),ll
        itype=0
c
c   check for type 5 :
c
        if(lnb.gt.15) go to 11
        do 1 i=1,lnb
        ll=l(i)
        if(ll.eq.'='.or.ll.eq.'('.or.ll.eq.')') go to 11
   1    if(ll.eq.'<'.or.ll.eq.'>') go to 11
        itype=5
        go to 99
c
c   check for '=' (if more than one or at beginning or end then illegal
c   if exactly one and not at beginning or end, could be legal type2-4,
c   if none, could be legal type 1)
c
  11    call findch(l,1,lnb,'=',nocc,neq)
        if(nocc.gt.1.or.neq.eq.1.or.neq.eq.lnb) go to 99
        if(nocc.eq.1) go to 10
c
c   check if legal type 1 :
c
        call findch(l,1,lnb,'(',nocc1,nc1)
        call findch(l,1,lnb,')',nocc2,nc2)
        if(nocc1*nocc2.ne.1.or.lnb.ne.nc2.or.nc1.eq.1) go to 99
        itype=1
        go to 99
c
c   check whether legal type 2,3,4 (first look for <)
c
  10    continue
        call findch(l,1,lnb,'<',nocc1,nc1)
        if(nocc1.eq.0) go to 20
c
c   could be type 4 :
c
        if(nocc1.gt.1.or.nc1.ne.neq+1) go to 99
        call findch(l,1,lnb,'>',nocc2,nc2)
        if(nocc2.ne.1.or.nc2.ne.lnb) go to 99
        itype=4
        go to 99
c
c   could be types 2 or 3 :
c
  20    continue
c
c   type 2 if no >,(,),  type 3 if no >, exactly one ( more than one
c   element to the right of = and exactly one ) which is at the end.
c
        call findch(l,1,lnb,'>',nocc2,nc2)
        if(nocc2.ne.0) go to 99
        call findch(l,1,lnb,'(',nocc3,nc3)
        call findch(l,1,lnb,')',nocc4,nc4)
        if(nocc3+nocc4.ne.0) go to 30
        itype=2
        go to 99
  30    if(nocc3*nocc4.ne.1.or.nc3.le.neq+1.or.nc4.ne.lnb) go to 99
        itype=3
c
c
  99    continue
        return
        end
c&WAITK
        subroutine waitk(n)
c*****************************************************************
c
c
c*****************************************************************
c
        integer*4 status,getc
        character ch

        status=getc(ch)
        n=ichar(ch)
        return
        end
