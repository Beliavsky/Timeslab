c&PLOT3D
c&PLT3D
        subroutine plt3d(args,nargs)
c*****************************************************************
c
c   Subroutine to handle the command
c
c   PLOT3D(x,y,z,nx,ny,v)
c
c*****************************************************************
c
#include 'tslabc'
        character args(nargs)*15
        integer*2 ickl,icki,ickr
        dimension vv(3)
c
c
        if(ickl(args(1),npx,nnx).eq.1) go to 99
        if(ickl(args(2),npy,nny).eq.1) go to 99
        if(ickl(args(3),npz,nnz).eq.1) go to 99
        if(icki(args(4),4,nx,1,2).eq.1) go to 99
        if(icki(args(5),5,ny,1,2).eq.1) go to 99
        ier=0
        if(nx*ny.gt.nnz) ier=3
        if(ny.gt.nny) ier=2
        if(nx.gt.nnx) ier=1
        if(ier.ne.0) then
                call error(args,ier,2)
                go to 99
        endif
        if(ickl(args(6),npv,nv).eq.1) go to 99
        if(nv.lt.3) then
                call error(args,6,2)
                go to 99
        endif
        call movxy(vv,array(nstart(npv)),12)
        n1=nx*ny+1
        if(iptpon.eq.1.and.nplt.gt.0) go to 973
        call smodea(6)
        call colora(0)
        call colora(kolor)
 973    continue
c
c
        
        call plot3d(array(nstart(npx)),array(nstart(npy)),
     1  array(nstart(npz)),nx,ny,nx,vv,wk,wk(n1),
     1  nx1,ny1,npixx,npixy,1,0,iwk,iwk(201))
c
c
	nn=nx*ny
	call max(wk,nn,smax,ind)
	call min(wk,nn,smin,ind)
	call max(wk(n1),nn,tmax,ind)
	call min(wk(n1),nn,tmin,ind)
        if(iptpon.eq.1) nplt=nplt+1
c        call wrtxt5(label(npz),40,nx1+2,ny1+npixy+2,1)
c        call pltend(args,nargs,smin,smax,tmin,tmax)
  99    continue
        return
        end
c&PLOT3D
        subroutine plot3d(x,y,z,nx,ny,ndim,v,s,t,nx1,ny1,npx,npy,
     1  kf,kb,nl,nr)
c*******************************************************************
c
c   Subroutine to produce a 3-d plot of values z(i,j) at grid points
c   (x1(i),y1(j)), i=1,...,nx, j=1,...,ny on any raster device.
c   See Note 4 for portability information.
c
c   Input:
c          x,y         : Vectors of length nx and ny with horizontal
c                        and vertical grid point values (both must be
c                        in increasing order).
c          z           : Matrix containing the values to be plotted.
c          nx,ny       : Number of horizontal and vertical grid points.
c          ndim        : Row dimension of z,x, and y in calling program
c                        (greater than or equal to nx).
c          v           : Vector of length 3 containg coordinates of
c                        viewpoint (see Note 1).
c          s,t         : Work matrices (see Note 2).
c          nx1,ny1     : Pixel coordinate of lower left hand corner of
c                        plotting region ( (0,0) is lower left corner
c                        of the screen, see Note 3).
c          npx,npy     : Number of pixels in horizontal and vertical
c                        directions of plotting region (see Note 3).
c          kf,kb       : Integers representing the numbers of the
c                        foreground and background colors.
c          nl,nr       : Integer work vectors whose dimension should
c                        exceed the number of vertical pixels on the
c                        screen being used.
c
c   Notes:
c          1) v(1) and v(2) cannot both be within the range of x and y.
c             If the user asks for a front, back, or side view, and
c             the view is not along a grid line, plot3d will modify
c             v so that it is a view along the next larger grid line.
c
c          2) The row dimension of z, s, and t must all be the same in
c             the calling program. This dimension must be at least nx.
c
c           3) The user is responsible for specifying the location and
c             size of the plotting region so that it fits on the screen.
c             The user is also responsible for setting the graphics mode
c             prior to calling plot3d.
c
c           4) The only non-portable subroutine used is
c
c                subroutine line(ic1,ir1,ic2,ir2,kolor)
c
c              which draws a line from pixel (ic1,ir1) to (ic2,ir2)
c              in color number kolor. This subroutine must consider
c              (0,0) to be the lower left corner of the screen.
c
c   Subprograms Called: prject,formfl,pfill,linep,ifpix,line
c
c*********************************************************************
c
        dimension x(nx),y(ny),s(ndim,1),t(ndim,1),z(ndim,1),v(3)
        dimension nl(1),nr(1)
c
c   Illegal View:
c
	if((x(1).lt.v(1).and.v(1).lt.x(nx)).and.(y(1).lt.v(2).
     1	and.v(2).lt.y(ny))) then
        	call wrtxt1('Illegal View',12,0,184,1)
	endif
c
c   Change view if necessary:
c
	i=1
	j=1
	if(v(1).gt.x(nx)) i=nx
	if(v(2).gt.y(ny)) j=ny
        if(x(1).lt.v(1).and.v(1).lt.x(nx)) then
                do 10 i=1,nx
  10            if(v(1).le.x(i)) go to 11
  11            v(1)=x(i)
        endif
        if(y(1).lt.v(2).and.v(2).lt.y(ny)) then
                do 12 j=1,ny
  12            if(v(2).le.y(j)) go to 13
  13            v(2)=y(j)
        endif
c
c   Do projection:
c
        call prject(x,y,z,nx,ny,ndim,v,s,t,smin,smax,tmin,tmax)
	rs=smax-smin
	rt=tmax-tmin
        call loadpg(smin,smax,tmin,tmax)
        call pgenv(smin,smax,tmin,tmax,0,0)
c        call pgsfs(2)
c
c   Draw Plot:
c
	if((ny-1.ge.j).and.(nx-1.ge.i)) then
		do 40 iy=ny-1,j,-1
		do 40 ix=nx-1,i,-1
  40		call dpatch(s,t,ndim,ix,iy,npx,npy,nx1,ny1,smin,rs,
     1			    tmin,rt,nl,nr,kf,kb)
	endif
c
c
	if((ny-1.ge.j).and.(i-1.ge.1)) then
		do 50 iy=ny-1,j,-1
		do 50 ix=1,i-1
  50		call dpatch(s,t,ndim,ix,iy,npx,npy,nx1,ny1,smin,rs,
     1			    tmin,rt,nl,nr,kf,kb)
	endif
c
c
	if((i-1.ge.1).and.(j-1.ge.1)) then
		do 60 iy=1,j-1
		do 60 ix=1,i-1
  60		call dpatch(s,t,ndim,ix,iy,npx,npy,nx1,ny1,smin,rs,
     1			    tmin,rt,nl,nr,kf,kb)
	endif
c
c
	if((j-1.ge.1).and.(nx-1.ge.i)) then
		do 70 iy=1,j-1
		do 70 ix=nx-1,i,-1
  70		call dpatch(s,t,ndim,ix,iy,npx,npy,nx1,ny1,smin,rs,
     1			    tmin,rt,nl,nr,kf,kb)
	endif
c
c
        return
        end
c&PRJECT
        subroutine prject(x,y,z,nx,ny,ndim,v,s,t,
     1	smin,smax,tmin,tmax)
c********************************************************************
c
c   Subroutine to find real coordinates of perspective projection
c   of the points (x(i),y(j),z(x(i),y(j))), i=1,...,nx,
c   j=1,...,ny when viewed from (v(1),v(2),v(3)). The s and t values
c   are returned in the matrices s and t respectively, while the
c   min's and max's of s and t are returned in smin, smax, tmin, tmax.
c
c********************************************************************
c
        dimension x(nx),y(ny),z(ndim,ny),v(3),s(ndim,ny),t(ndim,ny)
c
c
        v1=v(1)
        v2=v(2)
        v3=v(3)
        d=sqrt(v1**2+v2**2+v3**2)
        d2=d*d
        d1=sqrt(v1**2+v2**2)
c
        do 10 i=1,nx
        do 10 j=1,ny
        r1=x(i)*v1+y(j)*v2+z(i,j)*v3
        r=d2/(d2-r1)
        s(i,j)=r*(y(j)*v1-x(i)*v2)/d1
  10    t(i,j)=(v3+r*(z(i,j)-v3))*d/d1
c
        smin=s(1,1)
        smax=s(1,1)
        tmin=t(1,1)
        tmax=t(1,1)
        do 20 i=1,nx
        do 20 j=1,ny
        if(s(i,j).gt.smax) smax=s(i,j)
        if(s(i,j).lt.smin) smin=s(i,j)
        if(t(i,j).gt.tmax) tmax=t(i,j)
  20    if(t(i,j).lt.tmin) tmin=t(i,j)
c
c
        return
        end
c&DPATCH
        subroutine dpatch(s,t,ndim,i,j,npx,npy,nx1,ny1,smin,rs,
     1  tmin,rt,nl,nr,kf,kb)
c********************************************************************
c
c   Subroutine to draw one of the patches in the process of
c   constructing a 3-D plot.
c
c********************************************************************
c
        dimension s(ndim,1),t(ndim,1),ic(4),ir(4),nl(1),nr(1)
        dimension ss(4),tt(4)
c
        ss(1)=s(i,j)
        tt(1)=t(i,j)
        ss(2)=s(i,j+1)
        tt(2)=t(i,j+1)
        ss(3)=s(i+1,j+1)
        tt(3)=t(i+1,j+1)
        ss(4)=s(i+1,j)
        tt(4)=t(i+1,j)

        call pgsci(0)
        call pgsfs(1)
        call pgpoly(4,ss,tt)

        call pgsci(1)
        call pgsfs(2)
        call pgpoly(4,ss,tt)
        

c
        return
        end
c&IFPIX
        integer*2 function ifpix(x,npx,nx1,xmin,rx)
c********************************************************************
c
c   Convert real world coordinate x to pixel coordinate ifpix.
c
c********************************************************************
c
	ifpix=nx1+min(npx-1,max(0,nint(npx*(x-xmin)/rx)))
c
        return
        end
c&PFILL
        subroutine pfill(ic,ir,nl,nr,kf,kb)
c********************************************************************
c
c   Subroutine to draw the polygon having pixel vertices whose columns
c   are in the vector ic and rows are in the vector ir. The boundary
c   is in color kf, the interior color kb. The boundary is drawn from
c   vertices 1 to 2, 2 to 3, 3 to 4, 4 to 1. The arrays nl and nr must
c   be as long as the largest possible number of rows in the polygon.
c
c********************************************************************
c
        dimension ic(4),ir(4),nl(1),nr(1)
c
c   find top and bottom and left and right of required plotting region:
c
        icmin=ic(1)
        icmax=ic(1)
        irmin=ir(1)
        irmax=ir(1)
        do 10 i=1,4
        if(ic(i).lt.icmin) icmin=ic(i)
        if(ic(i).gt.icmax) icmax=ic(i)
        if(ir(i).lt.irmin) irmin=ir(i)
  10    if(ir(i).gt.irmax) irmax=ir(i)
c
c   polygon is horizontal or vertical line, so no interior:
c
        if(irmin.eq.irmax.or.icmin.eq.icmax) then
                call line(icmin,irmin,icmax,irmax,kf)
                return
        endif
c
c   nl(i), nr(i) to be left and right edge of region for row i-1:
c
c   initialize and find region (including boundary):
c
        do 20 i=irmin+1,irmax+1
        nl(i)=icmax
  20    nr(i)=icmin
        call linep(ic(1),ir(1),ic(2),ir(2),kf,nl,nr)
        call linep(ic(2),ir(2),ic(3),ir(3),kf,nl,nr)
        call linep(ic(3),ir(3),ic(4),ir(4),kf,nl,nr)
        call linep(ic(4),ir(4),ic(1),ir(1),kf,nl,nr)
c
c   have region now, so make all of it color kb:
c
        do 30 i=irmin,irmax
  30    call line(nl(i+1),i,nr(i+1),i,kb)
c
c   draw boundary:
c
        call line(ic(1),ir(1),ic(2),ir(2),kf)
        call line(ic(2),ir(2),ic(3),ir(3),kf)
        call line(ic(3),ir(3),ic(4),ir(4),kf)
        call line(ic(4),ir(4),ic(1),ir(1),kf)
c
c
        return
        end
c&LINEP
        subroutine linep(icol1,irow1,icol2,irow2,kol,nl,nr)
c********************************************************************
c
c   Modification of Bresenham's algorithm for drawing (in color kol)
c   from pixel (icol1,irow1) to (icol2,irow2) where pixel (0,0) is taken
c   to be the lower left corner of the screen. It finds each pixel in
c   the line and instead of turning it on, it updates (in nl and nr)
c   the left and right boundary for that row of pixels. To turn this
c   routine into a standard line drawing routine: 1) ``uncomment'' the
c   4 lines having calls to wrdota, and delete the 2 lines after each of
c   the uncommented lines (make sure to put the 10 and 20 on the
c   first 2 uncommented lines), 2) remove the arguments nl and nr and
c   the dimension statement, and 3) supply a routine called wrdota
c   that writes pixel (ic,ir) in color kol. Note that currently the
c   argument kol is ignored.
c
c********************************************************************
c
        dimension nl(1),nr(1)
c
c   Get left point into (ic0,ir0), right point into (ic1,ir1):
c
        ic0=min0(icol1,icol2)
        ic1=max0(icol1,icol2)
        ir0=irow1
        ir1=irow2
        if(icol1.gt.icol2) then
                ir0=irow2
                ir1=irow1
        endif
c
c   Horizontal line:
c
        if(ir0.eq.ir1) then
                do 10 i=ic0,ic1
c               call wrdota(i,ir0,kol)
                if(i.lt.nl(ir0+1)) nl(ir0+1)=i
  10            if(i.gt.nr(ir0+1)) nr(ir0+1)=i
                return
        endif
c
c   Vertical line:
c
        if(ic0.eq.ic1) then
                do 20 i=min0(ir0,ir1),max0(ir0,ir1)
c               call wrdota(ic0,i,kol)
                if(ic0.lt.nl(i+1)) nl(i+1)=ic0
  20            if(ic0.gt.nr(i+1)) nr(i+1)=ic0
                return
        endif
c
c   General line:
c
        ndc=ic1-ic0
        ndr=ir1-ir0
        indr=iabs(ndr)
        irinc=indr/ndr
        npix=max0(ndc,indr)+1
c
c   Shallow:
c
        if(ndc.ge.indr) then
                incd=2*(ndc-indr)
                incn=-2*indr
                ia=ndc-2*indr
                do 30 i=1,npix
c               call wrdota(ic0,ir0,kol)
                if(ic0.lt.nl(ir0+1)) nl(ir0+1)=ic0
                if(ic0.gt.nr(ir0+1)) nr(ir0+1)=ic0
                        ic0=ic0+1
                        if(ia.gt.0) then
                                ia=ia+incn
                        else
                                ir0=ir0+irinc
                                ia=ia+incd
                        endif
  30            continue
                return
c
c   Steep:
c
        else
                incd=2*(indr-ndc)
                incn=-2*ndc
                ia=indr-2*ndc
                do 40 i=1,npix
c               call wrdota(ic0,ir0,kol)
                if(ic0.lt.nl(ir0+1)) nl(ir0+1)=ic0
                if(ic0.gt.nr(ir0+1)) nr(ir0+1)=ic0
                        ir0=ir0+irinc
                        if(ia.gt.0) then
                                ia=ia+incn
                        else
                                ic0=ic0+1
                                ia=ia+incd
                        endif
  40            continue
        endif
c
        return
        end
