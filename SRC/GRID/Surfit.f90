SUBROUTINE  Read_Val_FktSurfit   !  dasurf.htm

  REAL(8) :: xx(11),yy(11),zz(121) !Bsp:
  REAL(8) :: x_datmin,y_datmin,min_px,min_py
  INTEGER :: i,ier,iopt,is,j,lwrk1,m,mx,my,nc
  INTEGER :: nmax,nxest,nyest
  INTEGER :: u,v,km,ne,bx,by,b1,b2
  REAL(8) :: ai,delta,eps,fp,s,ww,xb,xe,yb,ye
  REAL(8) :: m_half,xmin,xmax,ymin,ymax,breite,laenge
  REAL(8) :: diffmx,diffmy
  INTEGER :: InUnitV2=7
  REAL(8) :: x0H,y0H,dxH,dyH
  CHARACTER*1 :: DummyRow
  INTEGER(16) :: lwrk1_long

  IF(ALLOCATED(c_spl)) THEN
    RETURN
  END IF

  READ(10,*) iopt   ! specify whether a weighted least-square spline variation (-1)
                    !             or  a smoothing spline must be determined (0 or 1)
  READ(10,*) kx     ! degree of the spline x-direction
  READ(10,*) ky     ! degree of the spline y-direction
  READ(10,*) s      ! smoothing factor
  READ(10,*) fst    ! factor for strut of knots

  WRITE(*,*) 'file_namefkt',file_namefkt

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')
  IF (INDEX(file_namefkt,'.utm')==0) THEN
    !open(6,file='surf.out',status='UNKNOWN')

!   m denotes the number of data points
!   m>=(kx+1)*(ky+1) !!!Wichtig!!!
    read(InUnitV2,*) m
    ALLOCATE(x(1:m))
    ALLOCATE(y(1:m))
    ALLOCATE(z(1:m))
    ALLOCATE(w(1:m))

!   fetch the co-ordinate and function values of each data point.
    DO i=1,m
      read(InUnitV2,*) y(i),x(i),z(i)
    END DO
    x_datmin=x(1)
    y_datmin=y(m)
 
    IF (TRIM(file_namefkt).NE. 'dasurf.htm') THEN
      IF (conv_gk=='s'.OR.out_surf=='G') THEN
         DO i=1,m
           Call OGKGEO(y(i)/1000, x(i)/1000, breite, laenge)
           y(i)=(breite*4.0d0*ATAN(1.0d0))/180
           x(i)=(laenge*4.0d0*ATAN(1.0d0))/180
         END DO
      ELSE
        !!Betreff Nullpunkt ermitteln
        !!Wenn kleinste Koordinate aus *.dat beachtet werden soll
        !!z.Zt nur berechnet aber nicht verwendet
        !!gauss[x,y]_min aus grid-Domain nur Bezug,
        !!Nullpunkt wird erst beim OutputGMV verschoben
        IF (gaussx_min<x_datmin) THEN
          min_px=gaussx_min
          diffgx=0.0
        ELSE
          min_px=x_datmin
          diffgx=gaussx_min-x_datmin
        END IF
        IF (gaussy_min<y_datmin) THEN
          min_py=gaussy_min
          diffgy=0.0
        ELSE
          min_py=y_datmin
          diffgy=gaussy_min-y_datmin
        END IF
      END IF
    END IF ! nicht Input dasurf.htm
!   fetch an estimate of the standard deviation of the data values.
    read(InUnitV2,*) delta
  ELSE  
    READ(InUnitV2,*) nxH
    READ(InUnitV2,*) nyH
    READ(InUnitV2,*) x0H
    READ(InUnitV2,*) y0H
    READ(InUnitV2,*) dxH
    READ(InUnitV2,*) dyH
    READ(InUnitV2,*) m
    READ(InUnitV2,'(A1)') DummyRow
    ALLOCATE(x(1:m))
    ALLOCATE(y(1:m))
    ALLOCATE(z(1:m))
    ALLOCATE(w(1:m))
    WRITE(*,*) 'Anzahl m',m
    DO i=1,m
      read(InUnitV2,*) x(i),y(i),z(i)
!     IF (z(i)==0.0d0) THEN
!       z(i)=z(i)-50.0d0
!     END IF  
    END DO
    delta=1.0d0
  END IF

! the weights are set equal to delta**(-1)
      ww = 1./delta
      DO  i=1,m
        w(i) = ww
      END DO

! set up the boundaries of the approximation domain.
      !xb = -2.   !bedingt Bsp.: dasurf.htm
      !xe = 2.    ! ''
      !yb = -2.    ! ''
      !ye = 2.     ! ''
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      DO i=2,m
        xmin=MIN(xmin,x(i))
        xmax=MAX(xmax,x(i))
        ymin=MIN(ymin,y(i))
        ymax=MAX(ymax,y(i))
      END DO
      diffmx=xmax-xmin
      diffmy=ymax-ymin
      WRITE(*,*) 'xmin',xmin
      WRITE(*,*) 'xmax',xmax
      WRITE(*,*) 'ymin',ymin
      WRITE(*,*) 'ymax',ymax
      WRITE(*,*) 'conv_gk',conv_gk
      WRITE(*,*) 'out_surf',out_surf
      IF(conv_gk=='s'.OR.out_surf=='G') THEN
        xb=xmin-0.01*diffmx
        xe=xmax+0.01*diffmx
        yb=ymin-0.01*diffmy
        ye=ymax+0.01*diffmy
      ELSE
        xb=REAL(FLOOR(xmin))   ! nicht fuer Bogenmass
        xe=REAL(CEILING(xmax)) ! nicht fuer Bogenmass
        yb=REAL(FLOOR(ymin))   ! nicht fuer Bogenmass
        ye=REAL(CEILING(ymax)) ! nicht fuer Bogenmass
      END IF
 
! generate a rectangular grid for evaluating the splines.
!      mx = 11
!      my = 11
!      DO i=1,11
!        ai = i-6
!        xx(i) = ai*0.4
!        yy(i) = xx(i)
!      END DO


! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1) or a smoothing spline (iopt=
! 0 or 1) must be determined.
     ! iopt = 0    ! read
     ! kx = 3      ! read
     ! ky = 3      ! read
     ! s = 900000.   !s= 30 !read


! set up the dimension information.
!    -must specify an upper bound for the number of knots required
!     in the x- and y-directions respect.
!    -nxest >= 2*(kx+1), nyest >= 2*(ky+1)
!    -in most practical situation
!     nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient.
      m_half=MIN(m/2,fst)
      nxest = kx+1+CEILING(SQRT(m_half))
      nyest = ky+1+CEILING(SQRT(m_half))
      IF (nxest<(2*kx+2)) THEN
        nxest=2*kx+2
      END IF
      IF (nyest<(2*ky+2)) THEN
        nyest=2*ky+2
      END IF
      nmax = MAX(nxest,nyest)
      ALLOCATE(tx(1:nmax))
      ALLOCATE(ty(1:nmax))
      ALLOCATE(c_spl((nxest-kx-1)*(nyest-ky-1)))
      kwrk = m+(nxest-2*kx-1)*(nyest-2*ky-1)
      ALLOCATE(iwrk(1:kwrk))
      !computation for 'lwrk1'-> must specify the actual dimension of wrk1
      u = nxest-kx-1
      v = nyest-ky-1
      km = max(kx,ky)+1
      ne = nmax
      bx = kx*v+ky+1
      by = ky*u+kx+1
      if(bx.le.by) THEN
           b1 = bx
           b2 = b1+v-ky
      END IF
      if(bx.gt.by) THEN
           b1 = by
           b2 = b1+u-kx
      END IF
      lwrk1_long=16222049172_16
      WRITE(*,*) 'lwrk1 die erste',lwrk1_long
      lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      lwrk1_long=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      WRITE(*,*) 'T0',u*v,(2+b1+b2)
      WRITE(*,*) 'T1',u*v*(2+b1+b2)
      WRITE(*,*) 'T2',2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      WRITE(*,*) 'lwrk1_long berechnet',lwrk1_long
      WRITE(*,*) 'lwrk1 berechnet',lwrk1,u,v,b1,b2,km,m,ne,kx,ky
      ALLOCATE(wrk1(1:lwrk1))
      lwrk2=u*v*(b2+1)+b2
      ALLOCATE(wrk2(1:lwrk2))
      !im Bsp surfit def. : nxest = 15,nyest = 15,nmax = 15
      !                     !kwrk = 300, lwrk1 = 12000, lwrk2 = 6000
 
! choose a value for eps
      eps=0.1e-01

! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1)
 
      IF(iopt==-1) THEN
!if the computation mode iopt=-1 is used, the values tx(kx+2),
!    c          ...tx(nx-kx-1) must be supplied by the user, before entry.
!    c          see also the restrictions (ier=10).
!if iopt=-1: 2*kx+2<=nx<=nxest
!    c                        xb<tx(kx+2)<tx(kx+3)<...<tx(nx-kx-1)<xe
!    c                        2*ky+2<=ny<=nyest
!    c                        yb<ty(ky+2)<ty(ky+3)<...<ty(ny-ky-1)<ye
        !kx = 3
        !ky = 3
        tx_nx=nmax
        ty_ny=nmax
        ai=(xe-xb)/(tx_nx-kx-1-kx-2+2)
        tx(kx+2)=xb+ai
        DO i=kx+3,tx_nx-kx-1
          tx(i)=tx(i-1)+ai
          WRITE(*,*) i,tx(i),xb,xe
        END DO
        ai=(ye-yb)/(ty_ny-ky-1-ky-2+2)
        ty(ky+2)=yb+ai
        DO i=ky+3,ty_ny-ky-1
          ty(i)=ty(i-1)+ai
          WRITE(*,*) i,ty(i),yb,ye
        END DO
      ELSE
        tx_nx=nmax
        ty_ny=nmax
      END IF

! spline approximations of degree k
      WRITE(*,*) 'CALL surfit',iopt
      WRITE(*,*) 'm',m
      WRITE(*,*) 'xb',xb
      WRITE(*,*) 'xe',xe
      WRITE(*,*) 'yb',yb
      WRITE(*,*) 'ye',ye
      WRITE(*,*) 'kx',kx
      WRITE(*,*) 'ky',ky
      WRITE(*,*) 's',s
      WRITE(*,*) 'nxest',nxest
      WRITE(*,*) 'nyest',nyest
      WRITE(*,*) 'nmax',nmax
      WRITE(*,*) 'eps',eps
      WRITE(*,*) 'tx_nx',tx_nx
      WRITE(*,*) 'tx',SIZE(tx)
      WRITE(*,*) 'tx_ny',ty_ny
      WRITE(*,*) 'ty',SIZE(ty)
      WRITE(*,*) 'c_spl',SIZE(c_spl)
      WRITE(*,*) 'fp',fp
      WRITE(*,*) 'wrk1',SIZE(wrk1)
      WRITE(*,*) 'lwrk1',lwrk1
      WRITE(*,*) 'wrk2',SIZE(wrk2)
      WRITE(*,*) 'lwrk2',lwrk2
      WRITE(*,*) 'iwrk',SIZE(iwrk)
      WRITE(*,*) 'kwrk',kwrk
      call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
        nmax,eps,tx_nx,tx,ty_ny,ty,c_spl,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
      WRITE(*,*) 'after surfit',ier

! evaluation of the spline approximation
!      call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xx,mx,yy,my,zz,   &
!         wrk2,lwrk2,iwrk,kwrk,ier)

    CLOSE(UNIT=InUnitV2)
END SUBROUTINE Read_Val_FktSurfit
