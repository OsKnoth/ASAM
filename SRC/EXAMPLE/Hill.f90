MODULE Town_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Thermodynamic_Mod
  USE Physics_Mod
  USE Random_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: NL=1.0d-4 !1.0d-2
  REAL(RealKind) :: NU=2.06d-2 !1.0d-2
  REAL(RealKind) :: H=0.4d0  !1.0d-2
  REAL(RealKind), PARAMETER :: th0=289.2d0 !293.16d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0,SpeedMax=10.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0  
  REAL(RealKind) :: Inflow_lenx=1.0d0  
  REAL(RealKind) :: Inflow_leny=1.0d0  
  REAL(RealKind) :: Inflow_lenz=1.0d0  
  REAL(RealKind) :: offset_x=1.0d0  
  REAL(RealKind) :: offset_x1=1.0d0  
  REAL(RealKind) :: offset_y=1.0d0  
  REAL(RealKind) :: offset_y1=1.0d0  
  REAL(RealKind) :: offset_z=1.0d0  
  REAL(RealKind) :: offset_z1=1.0d0  
  REAL(RealKind) :: intenz=0.1d0
  REAL(RealKind) :: phase(1000)
  
  NAMELIST /Example/    &
                    SpeedMax , &
                    uMax , &
                    vMax , &
                    TkeMax ,&
                    DisMax ,&
                    intenz ,&
                    NL ,&
                    NU ,&
                    N

END MODULE Town_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Town_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells


SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile
SUBROUTINE InputExample(FileName)
  USE Town_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos,i
  REAL(RealKind) :: r
  CHARACTER(300) :: Line
  Inflow_lenx=domain%nx
  Inflow_leny=domain%ny
  Inflow_lenz=domain%nz
  offset_x=domain%x0
  offset_x1=domain%x1
  offset_y=domain%y0
  offset_y1=domain%y1
  offset_z=domain%z0
  offset_z1=domain%z1
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
!  ALLOCATE(TimeEnvi(5))
!  TimeEnvi(1)=StartTime
!  TimeEnvi(2)=10800.0d0
!  TimeEnvi(3)=46800.0d0
!  TimeEnvi(4)=64800.0d0
!  TimeEnvi(5)=EndTime

  DO i=1,1000
    Call Random_NUmber(r)
    phase(i)=r
!WRITE(*,*)'phase',i,phase(i)
  ENDDO  


  ALLOCATE(TimeEnvi(6))
  TimeEnvi(1)=StartTime
  TimeEnvi(2)=10800.0d0
  TimeEnvi(3)=14400.0d0
  TimeEnvi(4)=57600.0d0
  TimeEnvi(5)=64800.0d0
  TimeEnvi(6)=EndTime
!  DO i=2,20
!    TimeEnvi(i)=TimeEnvi(i-1)+20.0d0
!  END DO
END SUBROUTINE InputExample

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  REAL(RealKind) :: SL,SU
  REAL(RealKind) :: t1,t2,pH
  REAL(RealKind) :: s1,s2

!  S=N*N/Grav
!  ThLoc=th0*exp(z*S)
!  IF (N>Zero) THEN
!    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
!  ELSE
!    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
!  END IF
!  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)


  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThLoc=301.0*EXP(SL)
    t1=Grav*Grav*kappa/(Rd*301.0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*301.0*NL*NL)*EXP(-SL)
    PLoc=p0 &
            *(1.0d0+t2-t1)**(1.0d0/kappa)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
!    ThLoc=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
    ThLoc=(301.0)*EXP(SL+SU)
!    t1=Grav*Grav*kappa/(Rd*(290.0 + (z - 800)/100.0)*NL*NL)
!    t2=Grav*Grav*kappa/(Rd*(290.0 + (z - 800)/100.0)*NL*NL)*EXP(-SL)
    t1=Grav*Grav*kappa/(Rd*(301.0)*NL*NL)
    t2=Grav*Grav*kappa/(Rd*(301.0)*NL*NL)*EXP(-SL)
    pH=p0 &
       *(1.0d0+t2-t1)**(1.0d0/kappa)
!    s1=Grav*Grav*kappa/(Rd*(290.0 + (z - 800)/100.0)*NU*NU)*EXP(-SL)
!    s2=Grav*Grav*kappa/(Rd*(290.0 + (z - 800)/100.0)*NU*NU)*EXP(-SL-SU)
    s1=Grav*Grav*kappa/(Rd*(301.0)*NU*NU)*EXP(-SL)
    s2=Grav*Grav*kappa/(Rd*(301.0)*NU*NU)*EXP(-SL-SU)
    PLoc=p0 &
       *(1.0d0+t2-t1+s2-s1)**(1.0d0/kappa)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF
!  IF (z<=H) THEN
!    SL=NL*NL*z/Grav
!    ThLoc=290.0*EXP(SL)
!    pLoc=p0*(One-Grav/(Cpd*290.0*(SL))*(One-EXP(-(SL)*z)))**(Cpd/Rd)
!  ELSE
!    SL=NL*NL*H/Grav
!    SU=NU*NU*(z-H)/Grav
!    ThLoc=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
!    pLoc=p0*(One-Grav/(Cpd*290.0*(SL))*(One-EXP(-(SL)*H)))**(Cpd/Rd) &
!    pLoc=p0*(One-Grav/(Cpd*(290.0 + (z - 800)/100.0)*(SL+SU))*(One-EXP(-(SL+SU)*z)))**(Cpd/Rd)
!  END IF
!  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  REAL(RealKind) :: SL,SU
!  S=N*N/Grav
!  IF(z <= 800) THEN
!  ThLoc=290.0  *exp(z*S) !th0*exp(z*S)
!  ELSE
!  ThLoc=290.0 *exp(z*S) + (z - 800)/100.0 !th0*exp(z*S)
!  ENDIF 
!!  ThLoc=th0*exp(z*S)
!  IF (N>Zero) THEN
!    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
!  ELSE
!    pLoc=p0*(One-kappa*Grav*z/(Rd*ThLoc))**(Cpd/Rd)  !th0))**(Cpd/Rd)
!  END IF
!  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
!  IF (.NOT.Anelastic) THEN
!    RhoProf=Zero
!  END IF  

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThLoc=301.0*EXP(SL)
    pLoc=p0*(One-Grav/(Cpd*301.0*(SL))*(One-EXP(-(SL)*z)))**(Cpd/Rd)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
!    ThLoc=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
    ThLoc=(301.0)*EXP(SL+SU)
!    pLoc=p0*(One-Grav/(Cpd*(290.0 + (z - 800)/100.0)*(SL+SU))*(One-EXP(-(SL+SU)*z)))**(Cpd/Rd)
    pLoc=p0*(One-Grav/(Cpd*(301.0)*(SL+SU))*(One-EXP(-(SL+SU)*z)))**(Cpd/Rd)
  END IF
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S
  REAL(RealKind) :: ThetaLoc
!  REAL(RealKind) :: z
  REAL(RealKind) :: SL,SU
!  S=N*N/Grav
!  ThProfFun=th0*exp(z*S)

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThProfFun=301.0*EXP(SL)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
!    ThProfFun=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
    ThProfFun=(301.0)*EXP(SL+SU)
  END IF

!  IF(z <= 800) THEN
!  ThProfFun=290.0  *exp(z*S) !th0*exp(z*S)
!  ELSE
!  ThProfFun=290.0 *exp(z*S) + (z - 800)/100.0 !th0*exp(z*S)
!  ENDIF 
  IF (.NOT.Anelastic) THEN
    ThProfFun=Zero
  END IF  
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: RhoLoc,ThLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  ThLoc=thProfFun(x,y,z,zHeight,Time)
  QvProfFun=qvr*qvs(RhoLoc,ThLoc,ThLoc)
END FUNCTION QvProfFun

FUNCTION SpeedStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  REAL(REALKIND) :: SpeedStart1
  REAL(RealKind) :: SpeedStart,Speed
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: U0,wx,wy,wz,k,knull,AMP=1.0
  INTEGER        :: ss,s
!  //////turb = ran1()
  CALL Random_Number(r)
  SpeedStart=0.0d0
  Speed = (SQRT(uMax**Two+vMax**Two))
  knull=18.0
  k_end = int(Inflow_leny/4.0)
  DO s = 1,1
     k = 2.0**s !((k_end)*(0.9**s)) !*(1.0 + ((- 0.5d-1 +1.0d-1*r)/((6.0*(0.8**1.0)))))
!   AMP =!   (((k**4.0)*(exp(-2.0*(k/knull)**2.0)))/((knull**4.0)*(exp(-2.0*(knull/knull)**2.0))))!**(1.0/10.0)
      DO ss = 1,2
!     w = ((-1.0)**ss)*k*(SQRT(uMax))/Inflow_leny 
!     wx = ((-1.0)**ss)*k*(vMax)/(offset_x1-offset_x)
     wy = ((-1.0)**ss)*k*(SpeedMax)/(offset_y1-offset_y)
!     wx = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_x1-offset_x)
!     wy = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_y1-offset_y)
!     wz = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_z1-offset_z)
SpeedStart1 =  AMP*SIN(k*((y-offset_y)/(offset_y1-offset_y))*2.0*3.1415 + &
               Time*(wy)*2.0*3.1415)! + k**(2.0))!  &
!             + AMP*SIN(k*((z-offset_z)/(offset_z1-offset_z))*2.0*3.1415 + &
!               Time*(wz)*2.0*3.1415 + k**(2.0))
  SpeedStart=SpeedStart + SpeedStart1
      ENDDO
  ENDDO
  SpeedStart=(intenz*Speed*SpeedStart/((s-1)*(ss-1))) + Speed
!  SpeedStart=(intenz*(SQRT((uMax**2.0)))*SpeedStart/((s-1)*(ss-1))) + uMax
!WRITE(*,*)'offset_x',offset_x,'offset_x1',offset_x1,'Inflow_lenx',Inflow_lenx
!WRITE(*,*)'offset_y',offset_y,'offset_y1',offset_y1,'Inflow_leny',Inflow_leny
! SpeedStart = Speed !+Time/10.0 !+ y/10.0 + !((y/10.0)-2.50)**2.0 !
! -((y/10.0)/2.0) + y/10.0  !Speed

END FUNCTION SpeedStart

FUNCTION SpeedStartU(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  REAL(REALKIND) :: SpeedStartU1
  REAL(RealKind) :: SpeedStartU,Speed
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: U0,wx,wy,wz,k,kmean,knull,AMP=0.20 !0.75
  INTEGER        :: ss,s
!  CALL Random_Number(r)
  SpeedStartU=0.0d0
  Speed = (SQRT(uMax**Two+vMax**Two))
  knull=18.0
  k_end = int(Inflow_leny/10.0)
  DO s = 1,10 !100 !100 !10 !k_end
!  CALL Random_Number(r)
!  IF(s.gt.80) AMP = 1.0
!  IF(s.gt.60.and.s.le.80) AMP = 0.8
!  IF(s.gt.40.and.s.le.60) AMP = 0.6
!  IF(s.gt.20.and.s.le.40) AMP = 0.4
!  IF(s.gt.0.and.s.le.20) AMP = 0.2
     phase(s)=phase(s)
     k = 20.0*((0.8)**s)
!     k = (50.0*(((0.6/50.0)**(1.0/100.0))**s))
!     k = (300.0*(((0.6/300.0)**(1.0/100.0))**(s + r)))
!     k = (Inflow_leny*(((0.1/Inflow_leny)**(1.0/100.0))**s))  !  ((Inflow_leny)*(((kleinsteWellenzahl/Inflow_leny)**(1.0/N))**s))! ohne Einheit, nur k Wellen auf Strecke ! autom. 1/10 Wellen von inflow_length = k_end von 2 - k_end
     kmean = (k/(offset_y1-offset_y)) !Einheit jetzt [1/m]
     k = kmean 
     wy = kmean*(Speed)
     wz = kmean*(Speed)
!     AMP = (0.4**Two) * ((79.0*(wy/(1.0+4.7*(0.4))))/(1.0+263.0*(wy/(1.0+4.7*(0.4)))**(5.0/3.0))) &
!          * ((1.0+2.5*((0.4)**0.6))/((1.0+4.7*(0.4))**(2.0/3.0)))
SpeedStartU1 =  AMP*SIN(kmean*((y-offset_y))*2.0*3.1415 + phase(s)*offset_y1)* &
               (COS(Time*(wy)*2.0*3.1415))
  SpeedStartU=SpeedStartU + SpeedStartU1
  ENDDO
  SpeedStartU=intenz*SpeedStartU/1.0
  SpeedStartU=SpeedStartU + Speed 
!  IF(Time.lt.1.0d-4) SpeedStartU=Speed
!  SpeedStartU=(intenz*Speed*SpeedStartU/((s-1)*(ss-1))) + Speed 
  SpeedStartU=(((log(z/0.1d-1)))/((log(0.66/0.1d-1))))*SpeedStartU ! Zref = 8m
!IF(y.lt.25.0.or.y.gt.325.0) SpeedStartU=(((log(z/0.1d0)))/((log(8.0/0.1d0))))*uMax 
!  SpeedStartU=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(30.0/0.1d0))))*SpeedStartU
 
END FUNCTION SpeedStartU              

FUNCTION SpeedStartV(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  REAL(REALKIND) :: SpeedStartV1
  REAL(RealKind) :: SpeedStartV,Speed
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: U0,wx,wy,wz,k,kmean,knull,AMP=0.20 !0.75
  INTEGER        :: ss,s
!  CALL Random_Number(r)
  SpeedStartV=0.0d0
  Speed = (SQRT(uMax**Two+vMax**Two))
  knull=18.0
  k_end = int(Inflow_leny/10.0)
  DO s = 1,10 !100 !k_end
!  CALL Random_Number(r)
!  IF(s.gt.80) AMP = 1.0
!  IF(s.gt.60.and.s.le.80) AMP = 0.8
!  IF(s.gt.40.and.s.le.60) AMP = 0.6
!  IF(s.gt.20.and.s.le.40) AMP = 0.4
!  IF(s.gt.0.and.s.le.20) AMP = 0.2
     phase(s)=phase(s)
     k = 20.0*((0.8)**s)
!     k = (50.0*(((0.6/50.0)**(1.0/100.0))**s))
!     k = (300.0*(((0.6/300.0)**(1.0/100.0))**(s + r)))
!     k = (Inflow_leny*(((0.1/Inflow_leny)**(1.0/100.0))**s))  !  ((Inflow_leny)*(((kleinsteWellenzahl/Inflow_leny)**(1.0/N))**s))! ohne Einheit, nur k Wellen auf Strecke ! autom. 1/10 Wellen von inflow_length = k_end von 2 - k_end
     kmean = (k/(offset_y1-offset_y)) !Einheit jetzt [1/m]
     k = kmean 
     wy = kmean*(Speed)
     wz = kmean*(Speed)
!     AMP = (0.4**Two) * ((79.0*(wy/(1.0+4.7*(0.4))))/(1.0+263.0*(wy/(1.0+4.7*(0.4)))**(5.0/3.0))) &
!          * ((1.0+2.5*((0.4)**0.6))/((1.0+4.7*(0.4))**(2.0/3.0)))
SpeedStartV1 =  AMP*COS(kmean*((y-offset_y))*2.0*3.1415 + phase(s)*offset_y1)* &
               (SIN(Time*(wy)*2.0*3.1415))
  SpeedStartV=SpeedStartV + SpeedStartV1
  ENDDO
  SpeedStartV=intenz*SpeedStartV/1.0 
  SpeedStartV=SpeedStartV + 0.0 
!  IF(Time.lt.1.0d-4) SpeedStartV=0.0
  SpeedStartV=(((log(z/0.1d-1)))/((log(0.66/0.1d-1))))*SpeedStartV ! Zref = 8m
!IF(y.lt.25.0.or.y.gt.325.0) SpeedStartV=(((log(z/0.1d0)))/((log(8.0/0.1d0))))*vMax 
!  SpeedStartV=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(30.0/0.1d0))))*SpeedStartV

END FUNCTION SpeedStartV

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  COMPLEX(REALKIND) :: UStart1,UStart2
  REAL(REALKIND) :: UStart3,UStart4,UStart5
  REAL(RealKind) :: UStart
  REAL(RealKind) :: SpeedStartU,SpeedStartV
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r,alpha,xs,ys,cent,Times
  INTEGER        :: ss,s

  alpha = atan(vMax/(uMax+Eps))
  cent = offset_x1*sin(alpha)
    ys = y*cos(alpha) + cent -x*sin(alpha)
    xs = x*cos(alpha)+ y*sin(alpha)
  Times = Time - ((xs)/SQRT(uMax**Two+vMax**Two))
IF(cent.eq.offset_x1.or.cent.eq.0.0) Times = Time
  UStart=SpeedStartU(0.0,ys,z,zHeight,Times)*cos(alpha) - SpeedStartV(0.0,ys,z,zHeight,Times)*sin(alpha)
!IF(x.lt.1.0) WRITE(*,*)'u', x,y,Time,Times,cos(Time*(0.2)*2.0*3.1415),SpeedStartU(0.0,ys,z,zHeight,Times)*cos(alpha) &
!           ,SpeedStartV(0.0,ys,z,zHeight,Times)*sin(alpha),UStart

END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!   UStartE=(((log(z/0.1d0)))/((log(8.0/0.1d0))))*uMax
   UStartE=(((log(z/0.1d-1)))/((log(0.66/0.1d-1))))*uMax
! REAL(RealKind) :: x,y,z,zHeight,Time
!  CALL Random_Number(r)
! IF(Time<=10800) THEN
!   UStartE=8.0d0
! ELSE IF(Time>10800.and.TIME<=14400) THEN
!   UStartE=6.0d0
! ELSE IF(Time>14400.and.TIME<=57600) THEN
!   UStartE=6.0d0
! ELSE IF(Time>57600.and.TIME<=64800) THEN
!   UStartE=8.0d0
! ELSE
!   UStartE=8.0d0
! END IF
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  REAL(RealKind) :: VStart
  REAL(RealKind) :: SpeedStartV,SpeedStartU
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r,alpha,cent,ys,xs,Times
  INTEGER        :: ss,s
  alpha = atan(vMax/(uMax+Eps))
  cent = offset_x1*sin(alpha)
    ys = y*cos(alpha) + cent -x*sin(alpha)
    xs = x*cos(alpha)+ y*sin(alpha)
  Times = Time - ((xs)/SQRT(uMax**Two+vMax**Two))
IF(cent.eq.offset_x1.or.cent.eq.0.0) Times = Time

  VStart=SpeedStartV(0.0,ys,z,zHeight,Times)*cos(alpha) + SpeedStartU(0.0,ys,z,zHeight,Times)*sin(alpha)
!   ,SpeedStartU(0.0,ys,z,zHeight,Times)*sin(alpha),VStart

END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Town_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time,r
  CALL Random_Number(r)
!  VStartE=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(1000.0/0.1d0))))*vMax
!  VStartE=VStart(x,y,z,zHeight,Time)
!  VStartE = vMax !- 0.5d0 + 1.0d0*r
!   VStartE=(((log(z/0.1d0)))/((log(8.0/0.1d0))))*vMax
   VStartE=(((log(z/0.1d-1)))/((log(0.66/0.1d-1))))*vMax
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart,turb
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
!  IF(z.LT.900) THEN
!  turb = ran1()
!  WStart=-0.5d0+ 1.0d0*turb
!  ELSE
!  WSTART=0.0d0
!  ENDIF
END FUNCTION WStart

FUNCTION WStartE(x,y,z,zHeight,Time)
  USE Town_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStartE=WStart(x,y,z,zHeight,Time)
END FUNCTION WStartE

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,r
  REAL(RealKind) :: SL,SU
  CALL Random_Number(r)
!  S=N*N/Grav
!!  ThStart=th0*exp(z*S)
!!  ThStart=th0 ! -z/200.0 !th0*exp(z*S)
!  IF(z <= 800) THEN
!  ThStart=290.0  *exp(z*S) !th0*exp(z*S)
!  ELSE
!  ThStart=290.0 *exp(z*S) + (z - 800)/100.0 !th0*exp(z*S)
!  ENDIF 
  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThStart=301.0*EXP(SL)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
!    ThStart=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
    ThStart=(301.0)*EXP(SL+SU)
  END IF

END FUNCTION ThStart

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=0.0d0
END FUNCTION EnStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=Zero
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DHStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DHStart=Zero
END FUNCTION DHStart

FUNCTION DVStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DVStart=Zero
END FUNCTION DVStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  REAL(RealKind) :: SL,SU

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
!    ThStart=290.0*EXP(SL)
    PreStart=p0*(One-Grav/(Cpd*301.0*SL)*(One-EXP(-SL*z)))**(Cpd/Rd)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
!    ThStart=(290.0 + (z - 800)/100.0)*EXP(SL+SU)
    PreStart=p0*(One-Grav/(Cpd*(301.0 + (z - 800)/100.0)*(SL+SU))*(One-EXP(-(SL+SU)*z)))**(Cpd/Rd)
  END IF



!  S=N*N/Grav
!  IF (N>Zero) THEN
!    PreStart=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
!  ELSE
!    PreStart=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
!  END IF
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho

FUNCTION QsStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QsStart = 0.0
END FUNCTION QsStart

FUNCTION NvStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NvStart = 0.0
END FUNCTION NvStart

FUNCTION NcStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NcStart = 0.0
END FUNCTION NcStart

FUNCTION NrStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NrStart = 0.0
END FUNCTION NrStart

FUNCTION NiStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NiStart = 0.0
END FUNCTION NiStart

FUNCTION NsStart(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NsStart = 0.0
END FUNCTION NsStart

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1 = 0.0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2 = 0.0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Town_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3 = 0.0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Town_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
   OmeStart=Zero
END FUNCTION OmeStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  REAL(RealKind) :: Time
!  ThStartSoil=ThSoil
  IF (zSoil<=0.1d0) THEN
!    ThStartSoil=276.50d0
    ThStartSoil=290.50d0
  ELSE IF (zSoil<=0.5d0) THEN
!    ThStartSoil=277.70d0
    ThStartSoil=290.00d0
  ELSE IF (zSoil<=1.0d0) THEN
!    ThStartSoil=280.90d0
    ThStartSoil=289.50d0
  ELSE
!    ThStartSoil=286.50d0
    ThStartSoil=288.00d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  REAL(RealKind) :: Time
  IF (zSoil<=0.1d0) THEN
    QvStartSoil=2.d-1
  ELSE IF (zSoil<=1.0d0) THEN
    QvStartSoil=1.8d-1
  ELSE
    QvStartSoil=1.2d-1
  END IF
END FUNCTION QvStartSoil

FUNCTION QsStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  QsStart=Zero
END FUNCTION QsStart

FUNCTION NvStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NvStart=Zero

END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NcStart=Zero

END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NrStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NrStart=Zero

END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NiStart=Zero

END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NsStart=Zero

END FUNCTION NsStart

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=0.0d0
END FUNCTION OmeStart

