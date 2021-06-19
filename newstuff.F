	subroutine movxy(y,x,n)
	character*1 x(1),y(1)

	do 10 i=1,n
  10	y(i)=x(i)

	return
	end

	subroutine movct(x,n,c)
	character*1 x(1),c

	do 10 i=1,n
  10	x(i)=c

	return
	end

	subroutine smodea(n)

	return
	end

	double precision function inprod(x,y,n)
	dimension x(1),y(1)

	inprod=0.0
	do 10 i=1,n
 10	inprod=inprod+dble(x(i))*dble(y(i))

	return
	end

	subroutine speak(n,m)
	return
	end

      subroutine unif2(dseed,n,x)
c****************************************************************
c
c   Subroutine to fill the first n elements of the single precision
c   array x with iid U(0,1) random variables.
c 
c   The argument dseed is integer*4 and on input contains
c   a large whole number for a seed. On output it contains a new
c   seed that can be used as input the next time unif2 is called.
c
c****************************************************************
c
         dimension x(n)
         real*8 svnto5,t31m1
         data svnto5,t31m1/16807.d0,2147483647.d0/
	 integer*4 dseed
c
         do 10 i=1,n
            dseed=dmod(svnto5*dseed,t31m1)
 10      x(i)=dseed/t31m1
c
         return
         end

	subroutine colora(n)

	return
	end

	subroutine wrdota(n,m,k)
	return
	end

	subroutine rdcura(i1,i2,i3,i4,i5)
	return
	end

	subroutine curtyp(i,j)
	return
	end

	subroutine scvec(n,x)
	dimension x(1)

	return
	end

	subroutine blnksc(n)
	return
	end

	subroutine stcura(n,m,k)
	return
	end

	subroutine erase
	return
	end

	subroutine scrla(i1,i2,i3,i4,i5,i6,i7)
	return
	end

	subroutine wrtxt1(l,i1,i2,i3,i4)
	character*1 l(1)
	return
	end

	subroutine wrtxt3(l,i1,i2,i3,i4)
	character*1 l(1)
	return
	end

	subroutine wrtxta(l,i1,i2,i3,i4,i5)
	character*1 l(1)
	return
	end

	subroutine wrtxtb(l,i1,i2,i3,i4,i5)
	character*1 l(1)
	return
	end

	integer*2 function rddota(i,j)
	rddota=1
	return
	end

	subroutine echo
	return
	end

	subroutine txtcol(i,j,k)
	return
	end

	subroutine drrect(i1,i2,i3,i4,i5,i6)
	return
	end

	subroutine curdsk(n)
	return
	end

c	subroutine sortl(l,n,m)
c	character*1 l(1)
c	return
c	end

	subroutine rescn(args,nargs)
	return
	end

	subroutine rdscr(chars,wk,i)
	return
	end

	subroutine line(i1,i2,i3,i4,i5)
	return
	end

	integer*2 function inkeya()
	return
	end

	subroutine txtcls
	return
	end

	integer*2 function sswrqq()
	return
	end

	subroutine ckexc(n,m)
	return
	end

	integer*4 function timea()
	return
	end




