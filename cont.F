c&CONT
c&CONT
	subroutine cont(args,nargs)
c*****************************************************************
c
c   Subroutine to handle the command
c
c   CONTOUR(x,y,z,ctr,nx,ny,nctr)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
	integer*2 ickl,icki,ickr
c
c
	if(ickl(args(1),npx,nnx).eq.1) go to 99
	if(ickl(args(2),npy,nny).eq.1) go to 99
	if(ickl(args(3),npz,nnz).eq.1) go to 99
	if(ickl(args(4),npc,nnc).eq.1) go to 99
	if(icki(args(5),5,nx,1,2).eq.1) go to 99
	if(icki(args(6),6,ny,1,2).eq.1) go to 99
	if(icki(args(7),7,nctr,1,2).eq.1) go to 99
	ier=0
	if(nx*ny.gt.nnz) ier=4
	if(nctr.gt.nnc) ier=3
	if(ny.gt.nny) ier=2
	if(nx.gt.nnx) ier=1
	if(ier.ne.0) then
		call error(args,ier,2)
		go to 99
	endif
	if(iptpon.eq.1.and.nplt.gt.0) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973	continue
	call cntour(array(nstart(npx)),array(nstart(npy)),
     1	array(nstart(npz)),nx,ny,nx,npixx,npixy,nx1,ny1,nctr,
     1	array(nstart(npc)),iwk,iwk(1001))
	xmin=array(nstart(npx))
	xmax=array(nstart(npx)+nx-1)
	ymin=array(nstart(npy))
	ymax=array(nstart(npy)+ny-1)
	if(iptpon.eq.1) nplt=nplt+1
c        call wrtxt5(label(npc),40,nx1+2,ny1+npixy+2,1)
c	call pltend(args,nargs,xmin,xmax,ymin,ymax)
c
c
  99    continue
        return
        end
        subroutine cntour(x,y,z,nx,ny,ndim,npx,npy,nx1,ny1,
     1                    nctrs,ctrs,ix,iy)
c******************************************************************
c
c   Subroutine to draw a contour plot on a raster device.
c
c   Input:
c         nx, ny  : Number of grid points in x and y directions
c         x, y    : Vectors containing x and y values of grid points.
c                   Each vector must be in increasing order.
c         z       : Matrix whose (i,j)th element is the value of the
c                   function at grid point (x(i),y(j)).
c         ndim    : Row dimension of z in calling program.
c         npx, npy: Size in pixels of plotting region.
c         nx1, ny1: Pixel column and row numbers of lower left hand
c                   corner of plotting region. Note that the lower left
c                   corner of the screen is column 0, row 0.
c         nctrs   : Number of contours.
c         ctrs    : Vector containing contour values in
c                   increasing order.
c
c   Output:
c         ix, iy  : Vectors defined by: ix(i) = pixel column of x(i).
c                                       iy(i) = pixel row of y(i).
c
c   Notes:
c         1) The outline of the plotting region is drawn, but no grid.
c         2) The contours are all drawn in the same color. To change
c            this, one need only modify the line subroutine.
c         3) The user is responsible for setting the graphics mode
c            before calling cntour and for making sure that the
c            values of npx, npy, nx1, ny1 are selected so that
c            the plotting region fits on the screen.
c         4) The only graphics routine called by cntour is line.
c
c******************************************************************
c
        dimension x(nx),y(ny),z(ndim,1),ctrs(nctrs),ix(nx),iy(ny)
c
c   Draw lines outlining plotting region:
c
c        call line(nx1,ny1,nx1,ny1+npy-1,1)
c        call line(nx1,ny1+npy-1,nx1+npx-1,ny1+npy-1,1)
c        call line(nx1+npx-1,ny1+npy-1,nx1+npx-1,ny1,1)
c        call line(nx1+npx-1,ny1,nx1,ny1,1)
c
        rx=x(nx)-x(1)
        ry=y(ny)-y(1)
        xmin=x(1)
	xmax=x(nx)
        ymin=y(1)
	ymax=y(ny)
	call pgenv(xmin,xmax,ymin,ymax,0,0)

c
c   Find pixel coordinates of grid points:
c
        do 5 i=1,nx
   5    ix(i)=ifpix(x(i),npx,nx1,xmin,rx)
        do 6 i=1,ny
   6    iy(i)=ifpix(y(i),npy,ny1,ymin,ry)
c
c   Loop over i = column of boxes, j = row of boxes, k = which contour
c
        do 20 i=1,nx-1
c
        xleft=x(i)
        xright=x(i+1)
        xdiff=xright-xleft
c
        do 20 j=1,ny-1
c
        ylow=y(j)
        yup=y(j+1)
        ydiff=yup-ylow
c
c       z2 .------------. z3
c           |          |
c           |          |
c           |          |
c           |          |
c       z1 .------------. z4
c
        z1=z(i,j)
        z2=z(i,j+1)
        z3=z(i+1,j+1)
        z4=z(i+1,j)
c
c
        do 10 k=1,nctrs
        cc=ctrs(k)
c
c   nzg = # of z's greater or equal to current contour
c
        nzg=0
        if(z1.ge.cc) nzg=nzg+1
        if(z2.ge.cc) nzg=nzg+1
        if(z3.ge.cc) nzg=nzg+1
        if(z4.ge.cc) nzg=nzg+1
c
c   If nzg = 0, then all later contours are above too so jump out
c   of contour loop:
c
        if(nzg.eq.0) go to 15
        if(nzg.eq.4) go to 10
c
c   ns = sum of vertex numbers of z's that are greater or equal to cc
c
        ns=0
        if(z1.ge.cc) ns=ns+1
        if(z2.ge.cc) ns=ns+2
        if(z3.ge.cc) ns=ns+3
        if(z4.ge.cc) ns=ns+4
c
c   nzg = 3 and ns = j is the same as nzg = 1 and ns = 10 - j
c
        if(nzg.eq.3) then
                ns=10-ns
                nzg=1
        endif
c
        if(nzg.eq.1) then
c
                go to(11,12,13,14) ns
c
c   left to bottom:
c
  11            xx=xpix(xleft,xdiff,cc,z1,z4)
                yy=xpix(ylow,ydiff,cc,z1,z2)
                call mline(xleft,yy,xx,ylow,1)
                go to 10
c
c   left to top:
c
  12            xx=xpix(xleft,xdiff,cc,z2,z3)
                yy=xpix(ylow,ydiff,cc,z1,z2)
                call mline(xleft,yy,xx,yup,1)
                go to 10
c
c   top to right:
c
  13            xx=xpix(xleft,xdiff,cc,z2,z3)
                yy=xpix(ylow,ydiff,cc,z4,z3)
                call mline(xx,yup,xright,yy,1)
                go to 10
c
c   right to bottom:
c
  14            xx=xpix(xleft,xdiff,cc,z1,z4)
                yy=xpix(ylow,ydiff,cc,z4,z3)
                call mline(xx,ylow,xright,yy,1)
                go to 10
c
        endif
c
        if(nzg.eq.2) then
c
c   top to bottom:
c
                if(ns.eq.3.or.ns.eq.7) then
                xx1=xpix(xleft,xdiff,cc,z2,z3)
                xx2=xpix(xleft,xdiff,cc,z1,z4)
                call mline(xx1,yup,xx2,ylow,1)
                go to 10
                endif
c
c   left to right:
c
                if(ns.eq.5) then
                yy1=xpix(ylow,ydiff,cc,z1,z2)
                yy2=xpix(ylow,ydiff,cc,z4,z3)
                call mline(xleft,yy1,xright,yy2,1)
                go to 10
                endif
c
c   two mlines with negative slopes:
c
                if(ns.eq.4) then
                xx=xpix(xleft,xdiff,cc,z1,z4)
                yy=xpix(ylow,ydiff,cc,z1,z2)
                call mline(xleft,yy,xx,ylow,1)
                xx=xpix(xleft,xdiff,cc,z2,z3)
                yy=xpix(ylow,ydiff,cc,z4,z3)
                call mline(xx,yup,xright,yy,1)
                go to 10
                endif
c
c   two lines with positive slopes:
c
                if(ns.eq.6) then
                xx=xpix(xleft,xdiff,cc,z2,z3)
                yy=xpix(ylow,ydiff,cc,z1,z2)
                call mline(xleft,yy,xx,yup,1)
                xx=xpix(xleft,xdiff,cc,z1,z4)
                yy=xpix(ylow,ydiff,cc,z4,z3)
                call mline(xx,ylow,xright,yy,1)
                go to 10
                endif
        endif

c
  10    continue
  15    continue
  20    continue
c
c
        return
        end
	subroutine mline(x1,y1,x2,y2,nc)
	dimension x(2),y(2)
	x(1)=x1
	x(2)=x2
	y(1)=y1
	y(2)=y2
	call pgline(2,x,y)
	return
	end

        real*4 function xpix(x1,xd,c,z1,z2)
c******************************************************************
c
c   Interpolate. x1 is lower x value, xd is difference
c   in x values.
c
c******************************************************************
c
        xpix=x1+xd*(z1-c)/(z1-z2)
        return
        end
