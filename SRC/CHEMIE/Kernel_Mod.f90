MODULE Kernel_MOD

  USE Kind_Mod
  USE Parameter_Mod

  IMPLICIT NONE

CONTAINS

FUNCTION Long_cgs(x,y)

  REAL(RealKind) :: Long_cgs

  REAL(RealKind) :: x,y

  REAL(RealKind) :: xd,yd
  REAL(RealKind) :: alpx,betx,alpy,bety
  REAL(RealKind) :: rx,ry,rxt,ryt
  REAL(RealKind) :: vx,vy
  REAL(RealKind) :: ec


! x,y - masses
! xd,yd - diameters
! betx,alpx and bety,alpy for computing the velocities vx,vy (alpx*x**betx)
  xd=(6.d0*x/Pi)**(1.d0/3.d0)
  yd=(6.d0*y/Pi)**(1.d0/3.d0)
  IF (xd.lt.79.37d-4)THEN
    betx=2.d0/3.d0
    alpx=0.45795d+06
  ELSEIF (xd.ge.79.37d-4.and.xd.lt.800.d-4) THEN
    betx=1.d0/3.d0
    alpx=4.962d+03
  ELSEIF (xd.ge.800.d-4.and.xd.lt.4031.74d-4) THEN
    betx=1.d0/6.d0
    alpx=1.732d+03
  ELSEIF (xd.ge.4031.74d-4) THEN
    betx=0.d0
    alpx=917.d0
  ENDIF
  IF (yd.lt.79.37d-4)THEN
     bety=2.d0/3.d0
     alpy=0.45795d+06
  ELSEIF (yd.ge.79.37d-4.and.yd.lt.800.d-4) THEN
     bety=1.d0/3.d0
     alpy=4.962d+03
  ELSEIF (yd.ge.800.d-4.and.yd.lt.4031.74d-4) THEN
     bety=1.d0/6.d0
     alpy=1.732d+03
  ELSEIF (yd.ge.4031.74d-4) THEN
     bety=0.d0
     alpy=917.d0
  ENDIF
  vx=alpx*x**betx
  vy=alpy*y**bety
! rx,ry - radii
  rx=xd/2.d0
  ry=yd/2.d0
  IF (ry.lt.rx) THEN
    rxt=rx
    ryt=ry
    ry=rxt
    rx=ryt
  ENDIF
  IF (ry.ge.50.d-4)THEN
    ec=1.d0
  ELSEIF(rx.le.3.d-4.and.ry.lt.50.d-4)THEN
    ec=0.d0
  ELSE
    ec=4.5d4*ry*ry*(1.d0-3.d-4/rx)
  ENDIF
  IF (ec.lt.0.d0) ec=0.d0

  Long_cgs=Pi*(rx+ry)**2.d0*ec*abs(vx-vy)
END FUNCTION Long_cgs

FUNCTION Long_mks(x,y)

  REAL(RealKind) :: Long_mks

  REAL(RealKind) :: x,y

  REAL(RealKind) :: xd,yd
  REAL(RealKind) :: alpx,betx,alpy,bety
  REAL(RealKind) :: rx,ry,rxt,ryt
  REAL(RealKind) :: vx,vy
  REAL(RealKind) :: ec


! x,y - masses
! xd,yd - diameters
! betx,alpx and bety,alpy for computing the velocities vx,vy (alpx*x**betx)
  xd=(6.d0*x/1.d3/Pi)**(1.d0/3.d0)*1.0d2
  yd=(6.d0*y/1.d3/Pi)**(1.d0/3.d0)*1.0d2
  IF (xd.lt.79.37d-4)THEN
    betx=2.d0/3.d0
    alpx=0.45795d+06
  ELSEIF (xd.ge.79.37d-4.and.xd.lt.800.d-4) THEN
    betx=1.d0/3.d0
    alpx=4.962d+03
  ELSEIF (xd.ge.800.d-4.and.xd.lt.4031.74d-4) THEN
    betx=1.d0/6.d0
    alpx=1.732d+03
  ELSEIF (xd.ge.4031.74d-4) THEN
    betx=0.d0
    alpx=917.d0
  ENDIF
  IF (yd.lt.79.37d-4)THEN
     bety=2.d0/3.d0
     alpy=0.45795d+06
  ELSEIF (yd.ge.79.37d-4.and.yd.lt.800.d-4) THEN
     bety=1.d0/3.d0
     alpy=4.962d+03
  ELSEIF (yd.ge.800.d-4.and.yd.lt.4031.74d-4) THEN
     bety=1.d0/6.d0
     alpy=1.732d+03
  ELSEIF (yd.ge.4031.74d-4) THEN
     bety=0.d0
     alpy=917.d0
  ENDIF
  vx=alpx*x**betx
  vy=alpy*y**bety
! rx,ry - radii
  rx=xd/2.d0
  ry=yd/2.d0
  IF (ry.lt.rx) THEN
    rxt=rx
    ryt=ry
    ry=rxt
    rx=ryt
  ENDIF
  IF (ry.ge.50.d-4)THEN
    ec=1.d0
  ELSEIF(rx.le.3.d-4.and.ry.lt.50.d-4)THEN
    ec=0.d0
  ELSE
    ec=4.5d4*ry*ry*(1.d0-3.d-4/rx)
  ENDIF
  IF (ec.lt.0.d0) ec=0.d0

  Long_mks=Pi*(rx+ry)**2.d0*ec*abs(vx-vy)
  Long_mks=Long_mks/1.0d6
END FUNCTION Long_mks

FUNCTION FUCH(x,y)

  REAL(RealKind):: FUCH
  REAL(RealKind):: rx,ry !radius
  REAL(RealKind):: x,y !masses

  REAL(RealKind) :: x1,y1

  REAL(RealKind):: LAMDA_x, LAMDA_y !mean free path for particles
  REAL(RealKind):: DELTA_x,DELTA_y !variable for fuch interpolation
  REAL(RealKind):: ALPHA, BETA, GAMA, ENT1x,ENT2x,ENT3x,ENT1y,ENT2y,ENT3y
 
  x1=x
  y1=y

  rx=(3.0d0*x1/(4.0d0*PI)/1.d3)**(1.0d0/3.0d0)
  ry=(3.0d0*y1/(4.0d0*PI)/1.d3)**(1.0d0/3.0d0)
 
 

  LAMDA_x=2.0d0*DIFF_coef(rx)/(PI*TERM_vel(x1))
  LAMDA_y=2.0d0*DIFF_coef(ry)/(PI*TERM_vel(y1))

  ENT1x=(2.0d0*rx+LAMDA_x)**3.0d0
  ENT1y=(2.0d0*ry+LAMDA_y)**3.0d0

  ENT2x=(4.0d0*rx**2.0d0+LAMDA_x**2.0d0)**1.5d0
  ENT2y=(4.0d0*ry**2.0d0+LAMDA_y**2.0d0)**1.5d0

  ENT3x=6.0d0*rx*LAMDA_x
  ENT3y=6.0d0*ry*LAMDA_y

  DELTA_x=((ENT1x-ENT2x)/ENT3x)-2.0d0*rx
  DELTA_y=((ENT1y-ENT2y)/ENT3y)-2.0d0*ry


  ALPHA=4.0d0*PI*(rx+ry)*(DIFF_coef(rx)+DIFF_coef(ry))

  BETA=(rx+ry)/(rx+ry+SQRT(DELTA_x**2.0d0+DELTA_y**2.0d0))
  GAMA=4.0d0*(DIFF_coef(rx)+DIFF_coef(ry))/((rx+ry)*SQRT( &
       TERM_vel(x1)**2.0d0+TERM_vel(y1)**2.0d0))

  FUCH=ALPHA/(BETA+GAMA) 

 
END FUNCTION FUCH

END MODULE Kernel_Mod
