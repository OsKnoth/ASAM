PROGRAM Profile
  USE Profile_Mod
  IMPLICIT NONE

  INTEGER :: i
! INTEGER  INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
! REAL(RealKind) :: &
!     T, Y(90), YPRIME(90), TOUT, RWORK(100000), &
!     RPAR(1000), RTOL, ATOL
! REAL(RealKind) :: Theta0,ystart,tout0
! REAL(RealKind) :: Te,ThetaRho,p,e,qt,RhoD,Rho,RhoGes,Const0,Gamm0,KappaM,dThdz,dRHdz
! REAL(RealKind) :: wSave(9,2000)
  REAL(RealKind) :: TeStart     = 288.0d0  
  REAL(RealKind) :: PreStart    = 1013.25d2
  REAL(RealKind) :: RelHumStart = 0.00d0
! REAL(RealKind) :: RelHumInv   = 0.90d0
! REAL(RealKind) :: InvHeight   = 1500.0d0
! REAL(RealKind) :: N_d         = 0.01d0 
! REAL(RealKind) :: DeltaZ      = 5.0d0
  REAL(RealKind), ALLOCATABLE :: c(:,:)
! REAL(RealKind) :: zu,zo,Tu,To
! REAL(RealKind) :: z,zOut
! INTEGER :: iHeight,iz
! LOGICAL :: Inversion=.FALSE.
! LOGICAL :: LinRH=.FALSE.
! LOGICAL :: ConstRH=.TRUE.
! LOGICAL :: Interpol=.TRUE.
! INTEGER :: stat
! INTEGER :: Pos
! INTEGER :: CountSteps=0
! CHARACTER(8) :: ProfType
  CHARACTER(128) :: Filename
! CHARACTER(64) :: par

! INTEGER :: i,nwSave
! INTEGER :: DryPos,VapPos,LiqPos
! REAL(RealKind) :: T1


  WRITE(*,*) "use: ./run outputfilename"
  CALL get_command_argument(1,FileName)
  WRITE(*,*) 'CALL ComputeParameter'
  CALL ComputeParameter
  WRITE(*,*) 'CALL ReadInput'
  CALL ReadInput(FileName)


  WRITE(*,*) 'ProfileType ',ProfileType
  SELECT CASE(ProfileType)
    CASE('MoistRelHum')
      PreStart=Pre(1)
      TeStart=ProfTemp(1)
      RelHumStart=RH(1)
      WRITE(*,*) 'PreStart     ',PreStart
      WRITE(*,*) 'TeStart      ',TeStart
      WRITE(*,*) 'RelHumStart  ',RelHumStart
      CALL ComputeProfile(c,Height,Pre,ProfTemp,RH)
    CASE('MoistQTPotEquiv')
      WRITE(*,*) 'Equiv '
      CALL ComputeProfile(c,Height,Pre,ThE=ThE,rt=rt)
  END SELECT

  WRITE(*,*) 'FileNameOut  ',FileNameOut
  CALL OutputRes(nzProf,zMProf,c,FileNameOut)
      
END PROGRAM Profile
