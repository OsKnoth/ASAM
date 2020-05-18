MODULE NVP_Mod

! Owner: Dylia Willink 

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
  REAL(RealKind) :: H=800.0d0  !1.0d-2
  REAL(RealKind), PARAMETER :: th0=289.2d0 !293.16d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: ThInit=301.0d0
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
  CHARACTER*20 :: Problem=''
  LOGICAL :: xDir=.TRUE.
  LOGICAL :: yDir=.FALSE.
! NVP
  REAL(RealKind) :: zRauhNVP=0.75d0
  REAL(RealKind) :: zRauhInFlow=0.01d0
  REAL(RealKind) :: xRauhChange=0.0d0
! Profile
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: ProfInRho=.FALSE.
  LOGICAL :: ProfInThStart=.FALSE.
  LOGICAL :: ProfInQvStart=.FALSE.
! Perturbation of initial profile  
  LOGICAL :: Perturb=.FALSE.
  
  NAMELIST /Example/    &
                    Problem , &
                    uMax , &
                    vMax , &
                    TkeMax ,&
                    DisMax ,&
                    intenz ,&
                    NL ,&
                    NU ,&
                    xDir, &
                    yDir, &
                    ProfIn, & 
                    ProfInRho, & 
                    ProfInThStart, & 
                    ProfInQvStart, & 
                    Perturb, & 
                    zRauhInflow, &
                    xRauhChange, &
                    ThInit, &
                    N

END MODULE NVP_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE NVP_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

  IF (BoundCellLoc%xS<=xRauhChange) THEN
    BoundCellLoc%zRauh=zRauhInFlow
  END IF  

END SUBROUTINE SetBoundCells


SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE NVP_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos,i
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
!  TimeEnvi(4)=64800.0d0 !test marcelk
!  TimeEnvi(5)=EndTime




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
  USE NVP_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfInRho) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoFun=ProfileEqual(cInt,z)
  ELSE
    ThLoc=ThInit
    pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(CpD/Rd)
    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF  

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  ThLoc=ThInit
  pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(CpD/Rd)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE NVP_Mod
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
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE NVP_Mod
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


FUNCTION UStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  COMPLEX(REALKIND) :: UStart1,UStart2
  REAL(REALKIND) :: UStart3,UStart4,UStart5
  complex, parameter    :: i  = (0.0,1.0)
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: uType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: U0,wx,wy,wz,k,knull,AMP=1.0
  INTEGER        :: ss,s
!WRITE(*,*) 'x',offset_x,offset_y,inflow_lenx,inflow_leny
!  //////turb = ran1()
  CALL Random_Number(r)
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,uType,'uProf')
      Load=.FALSE.
    END IF
    UStart=ProfileEqual(cInt,z)
  ELSE  
    UStart=0.0d0
    knull=18.0
    k_end = int(Inflow_leny/4.0)
    DO s = 1,1 !k_end
       k = 2.0 !((10.0)*(0.9**s)) !*(1.0 + ((- 0.5d-1 + 1.0d-1*r)/((6.0*(0.8**1.0)))))
!40.0*(0.97)**s
!       k = 5.0
        DO ss = 1,2
!       w = ((-1.0)**ss)*k*(SQRT(uMax))/Inflow_leny 
       wx = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_x1-offset_x)
       wy = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_y1-offset_y)
       wz = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_z1-offset_z)
!IF(y.GT.offset_y-1.and.x.gt.1)   THEN
!UStart2   = 1.0*AMP*SIN(k*((x-offset_x)/Inflow_lenx)*2.0*3.1415 +
!Time*(w)*2.0*3.1415)   !- 0.5d-1 + 1.0d-1*r  
!UStart   = UStart + UStart2
!ELSE
!       UStart3 = AMP*SIN(k*((x-offset_x)/(offset_x1-offset_x))*2.0*3.1415 + &
!                 Time*(wx)*2.0*3.1415)! + k**(2.0))  &!- 0.5d-1 + 1.0d-1*r  
               UStart3 =  AMP*SIN(k*((y-offset_y)/(offset_y1-offset_y))*2.0*3.1415 + &
                 Time*(wy)*2.0*3.1415) !)*SIN(Time*(wy)*2.0*3.1415)! + k**(2.0))  &
!               + AMP*SIN(k*((z-offset_z)/(offset_z1-offset_z))*2.0*3.1415 + &
!                 Time*(wz)*2.0*3.1415 + k**(2.0))
    UStart=UStart + UStart3
!IF(ss==2)   UStart=UStart/SQRT(2.0)
!ENDIF
        ENDDO
    ENDDO
    UStart=UStart + uMax !(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*UStart/((s-1)*(ss-1))) + uMax
!WRITE(*,*)'offset_z',offset_z,'offset_z1',offset_z1,'Inflow_lenz',Inflow_lenz
!    IF(Time.LT.1.0d-4) UStart = uMax
!    IF(x.GT.(offset_x)) UStart = uMax
!     UStart = uMax - 0.5d0 +1.0d0*r
!IF(z.lt.1000.0)     UStart = (((0.44d0/0.4d0)*(log(z/0.1)))/((0.44d0/0.4d0)*(log(1000.0/0.1))))*UStart 
!     UStart = (((0.44d0/0.4d0)*(log(z/0.02)))/((0.44d0/0.4d0)*(log(2000.0/0.02))))*uMax 
!    0.5d-1 + 1.0d-1*r
!    UStart=
!    (((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(16.0/0.1d0))))*uMax -
!    0.5d-1 + 1.0d-1*r
!     UStart=uMax +0.01*(y-offset_y)
!     UStart = uMax - 0.5d0 +1.0d0*r
!     IF(y.lt.0.51) UStart =uMax -0.5d-1 +1.0d-1*r
!    UStart=UStart - 0.5d0 + 1.0d0*r
!      UStart=uMax !UStart !- 0.5d0 + 1.0d0*r
    IF (xDir) THEN
      UStart=((log((z+1.d-1)/1.d-1))/(log(100.0/1d-1)))*uMax !SpeedStartV
    ELSE  
      UStart=0.0d0
    END IF  
  END IF
  IF (Perturb.AND.Time>1.d-4) THEN
    UStart=UStart - 0.5d0 + 1.0d0*r
  END IF
! IF(Time>100.0)  UStart=uMax/10.0 !UStart !- 0.5d0 + 1.0d0*r
!   IF(z.lt.(3.0/186.0)) UStart = 0.0d0 -0.5d-1 +1.0d-1*r

END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  CALL Random_Number(r)
!  UStartE=Time/100.0   !UStart(x,y,z,zHeight,Time)* 
! IF(Time<=10800) THEN
!   UStartE=8.0d0
! ELSE IF(Time>10800.and.TIME<=46800) THEN
!   UStartE=6.0d0
! ELSE IF(Time>46800.and.TIME<=64800) THEN
!   UStartE=8.0d0
! ELSE
!   UStartE=8.0d0
! END IF
 IF(Time<=10800) THEN
   UStartE=8.0d0
 ELSE IF(Time>10800.and.TIME<=14400) THEN
   UStartE=6.0d0
 ELSE IF(Time>14400.and.TIME<=57600) THEN
   UStartE=6.0d0
 ELSE IF(Time>57600.and.TIME<=64800) THEN
   UStartE=8.0d0
 ELSE
   UStartE=8.0d0
 END IF
END FUNCTION UStartE

!FUNCTION UStart(x,y,z,zHeight,Time)
!  USE NVP_Mod
!  USE Rho_Mod 
!  IMPLICIT NONE
!  REAL(RealKind) :: UStart
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  REAL(RealKind) :: yL,yR
!  REAL(RealKind) :: U0
!!  IF(z.LT.500) THEN
!!  turb = ran1()
!!  UStart=uMax-0.005d0+1.0d-2*turb
!  UStart=uMax
!!  ELSE
!!  USTART=uMax+Four
!!  ENDIF
!END FUNCTION UStart

!FUNCTION UStartE(x,y,z,zHeight,Time)
!  USE NVP_Mod
!  USE UVW_Mod
!  IMPLICIT NONE
!  REAL(RealKind) :: UStartE
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  UStartE= (((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(16.0/0.1d0))))*uMax 
!!  UStartE=UStart(x,y,z,zHeight,Time)
!END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  INTEGER        :: k_end
  COMPLEX(REALKIND) :: VStart1,VStart2
  REAL(REALKIND) :: VStart3,VStart4,VStart5
  complex, parameter    :: i  = (0.0,1.0)
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: V0,wx,wy,wz,k,knull,AMP=1.0
  INTEGER        :: ss,s
!WRITE(*,*) 'x',offset_x,offset_y,inflow_lenx,inflow_leny
!  //////turb = ran1()
  CALL Random_Number(r)
  VStart=0.0d0
  knull=18.0
  k_end = int(Inflow_leny/4.0)
  DO s = 1,1 !k_end
     k = 2.0 !((10.0)*(0.9**s)) !*(1.0 + ((- 0.5d-1 + 1.0d-1*r)/((6.0*(0.8**1.0)))))
!40.0*(0.97)**s
!     k = 5.0
      DO ss = 1,2
!     w = ((-1.0)**ss)*k*(SQRT(uMax))/Inflow_leny 
     wx = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_x1-offset_x)
     wy = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_y1-offset_y)
     wz = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_z1-offset_z)
!IF(y.GT.offset_y-1.and.x.gt.1) THEN
!UStart2 = 1.0*AMP*SIN(k*((x-offset_x)/Inflow_lenx)*2.0*3.1415 +
!Time*(w)*2.0*3.1415) !- 0.5d-1 + 1.0d-1*r  
!UStart = UStart + UStart2
!ELSE
!     VStart3 = AMP*COS(k*((x-offset_x)/(offset_x1-offset_x))*2.0*3.1415 + &
!               Time*(wx)*2.0*3.1415)! + k**(2.0))  &!- 0.5d-1 + 1.0d-1*r  
             VStart3 =  AMP*COS(k*((y-offset_y)/(offset_y1-offset_y))*2.0*3.1415 + &
               Time*(wy)*2.0*3.1415) !*SIN(Time*(wy)*2.0*3.1415-3.1415/2.0)! + k**(2.0))  &
!             + AMP*COS(k*((z-offset_z)/(offset_z1-offset_z))*2.0*3.1415 + &
!               Time*(wz)*2.0*3.1415 + k**(2.0))
  VStart=VStart + VStart3
!IF(ss==2) VStart=VStart/SQRT(2.0)
!ENDIF
      ENDDO
  ENDDO
  VStart= VStart + vMax !(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*VStart/((s-1)*(ss-1))) + vMax
!WRITE(*,*)'offset_z',offset_z,'offset_z1',offset_z1,'Inflow_lenz',Inflow_lenz
!  IF(Time.LT.1.0d-4) VStart = vMax
!  IF(x.GT.(offset_x)) UStart = uMax
!  UStart=
!   VStart = (((0.44d0/0.4d0)*(log(z/0.02d0)))/((0.44d0/0.4d0)*(log(2000.0/0.02d0))))*VStart
!   VStart = (((0.44d0/0.4d0)*(log(z/0.02d0)))/((0.44d0/0.4d0)*(log(2000.0/0.02d0))))*vMax
!  0.5d-1 + 1.0d-1*r
!  UStart=
!  (((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(16.0/0.1d0))))*uMax -
!  0.5d-1 + 1.0d-1*r
!   VStart=0.5d0 !vMax! +0.01*(y-offset_y)
!   VStart = -1.0d0 !vMax !- 0.5d-1 +1.0d-1*r
!   VStart =  vMax - 0.5d-1 +1.0d-1*r
!    VStart =  vMax
  VStart=((log(z/1d-1))/(log(100.0/1d-1)))*vMax !SpeedStartV
  IF (yDir) THEN
    VStart=((log(z/1.d-1))/(log(100.0/1d-1)))*vMax !SpeedStartV
  ELSE  
    VStart=0.0d0
  END IF  
!   IF(z.lt.(3.0/186.0)) VStart =0.0d0 -0.5d-1 +1.0d-1*r
!  VStart=VStart - 0.5d0 + 1.0d0*r

END FUNCTION VStart

!FUNCTION VStart(x,y,z,zHeight,Time)
!  USE NVP_Mod
!  USE Rho_Mod 
!  IMPLICIT NONE
!  REAL(RealKind) :: VStart,turb
!  REAL(RealKind) :: x,y,z,zHeight,Time
!!  IF(z.LT.900) THEN
!!  turb = ran1()
!!  VStart=-0.5d0+ 1.0d0*turb
!!  ELSE
!  VSTART=0.0d0
!!  ENDIF
!END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time,r
  CALL Random_Number(r)
!  VStartE=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(1000.0/0.1d0))))*vMax
!  VStartE=VStart(x,y,z,zHeight,Time)
  VStartE = vMax !- 0.5d0 + 1.0d0*r
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart,turb
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION WStartE(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStartE=WStart(x,y,z,zHeight,Time)
END FUNCTION WStartE

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE ThProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S
  IF (ProfInThStart) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,thType,'ThDensProf')
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    ThStart=ThInit
  END IF  

END FUNCTION ThStart

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=0.0d0
END FUNCTION EnStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=One
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfInQvStart) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'RhoVProf')
      Load=.FALSE.
    END IF
    qvStart=ProfileEqual(cInt,z)
  ELSE
    QvStart=0.0d0
  END IF
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DHStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DHStart=Zero
END FUNCTION DHStart

FUNCTION DVStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DVStart=Zero
END FUNCTION DVStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE NVP_Mod
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
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho

FUNCTION QsStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QsStart = 0.0
END FUNCTION QsStart

FUNCTION NvStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NvStart = 0.0
END FUNCTION NvStart

FUNCTION NcStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NcStart = 0.0
END FUNCTION NcStart

FUNCTION NrStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NrStart = 0.0
END FUNCTION NrStart

FUNCTION NiStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NiStart = 0.0
END FUNCTION NiStart

FUNCTION NsStart(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NsStart = 0.0
END FUNCTION NsStart

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1 = 0.0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2 = 0.0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3 = 0.0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE NVP_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
   OmeStart=Zero
END FUNCTION OmeStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,Time)
  USE NVP_Mod
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
  USE NVP_Mod
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

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE NVP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE NVP_Mod
  USE Kind_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh

FUNCTION DampFun(z,Name)
  USE NVP_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=0.0d0
END FUNCTION DampFun
