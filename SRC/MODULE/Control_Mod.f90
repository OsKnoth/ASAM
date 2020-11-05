MODULE Control_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Sp_Mod
  USE Parallel_Mod

  IMPLICIT NONE

! Date
  INTEGER :: StartDay
  INTEGER :: Year
  INTEGER :: Month
  INTEGER :: Day
  REAL(RealKind) :: StartTime
  REAL(RealKind) :: EndTime
  REAL(RealKind) :: lng,lat  
! Restart
  REAL(RealKind) :: RestartTimeIncr
  LOGICAL :: Restart=.FALSE.
! Term Output
  INTEGER :: TermUnit=99
  INTEGER :: OutputUnit=98
  INTEGER :: InputUnit=97
  INTEGER :: OutputUnitS=96

  INTEGER :: ivar = 0        ! Partition variant (0: METIS)
  INTEGER :: ref_glob        ! Minimum refinement level
  CHARACTER(80) :: mi_wgts    ! Name of Weights Input File
  CHARACTER(80) :: ProfileOut ! Name of Profile Output File
  CHARACTER(80) :: ProfileOutGnu ! Name of Profile Output File for GnuPlot
  CHARACTER(40) :: Method
  LOGICAL :: InitPeer
  INTEGER :: CGMaxIterPre
  INTEGER :: CGMaxIterSch
  REAL(RealKind) :: CGTolPre
  REAL(RealKind) :: CGTolSch
  REAL(RealKind) :: BiCGStabTol
  INTEGER :: BiCGStabMaxIter
  REAL(RealKind) :: QMRTol
  REAL(RealKind) :: Fac=2.0d0
  CHARACTER(10) :: MethAdv
  LOGICAL :: JacTransport=.TRUE.
  LOGICAL :: JacTransportX=.TRUE.
  LOGICAL :: JacTransportY=.TRUE.
  LOGICAL :: JacTransportZ=.TRUE.
  LOGICAL :: JacPartial=.FALSE.
  LOGICAL :: JacSound=.TRUE.
  LOGICAL :: MultiEx=.FALSE.
  LOGICAL :: MultiMu=.FALSE.
  LOGICAL :: MultiMuTR=.FALSE.
  LOGICAL :: MultiTriT=.FALSE.
  LOGICAL :: MultiTriTB=.FALSE.
  LOGICAL :: MultiTriTR=.FALSE.
  INTEGER :: GaussIter2=1
  INTEGER :: GaussIter3=1
  INTEGER :: MultIter2=1
  INTEGER :: MultIter3=1
  REAL(RealKind) :: PreFac2=0.5d0
  REAL(RealKind) :: PreFac3=2.0d0
  REAL(RealKind) :: PreFac5=0.5d0
  REAL(RealKind) :: ShiftLoc=0.0d0
  LOGICAL :: Debug=.FALSE.
  CHARACTER*10 :: MethSound='FGMRES'
  LOGICAL :: CoarsePre=.FALSE.
  INTEGER :: StepGMRES=10
  INTEGER :: QMRMaxIter
  INTEGER :: RefLevelG
  LOGICAL :: PD ! positive-definite scheme
! Fast integration  
  INTEGER :: ns
  INTEGER :: nsMin
  LOGICAL :: LinearPressure=.FALSE.
  LOGICAL :: nsCompute=.FALSE.
  REAL(RealKind) :: GammaDiv=0.0d0
  LOGICAL :: GradFull=.TRUE.
  CHARACTER*20 :: IntFast='FB'
  LOGICAL :: Primitive=.TRUE.
  CHARACTER*2 :: AdvPr='V1'
  CHARACTER*8 :: GridType
  LOGICAL :: CloudFast=.TRUE.

  LOGICAL :: VolumeCheck=.FALSE.
  LOGICAL :: PrintNameLists=.FALSE.

  INTEGER :: uPos,uPosL,uPosR
  INTEGER :: vPos,vPosL,vPosR
  INTEGER :: wPos,wPosL,wPosR
! Positionen for turbulence 
! Tke closure 
  INTEGER :: tkePos
  INTEGER :: disPos
  INTEGER :: omePos
! Closure nach Herzog
  INTEGER :: tkeHPos
  INTEGER :: tkeVPos
  INTEGER :: LenPos

  INTEGER :: thPos
  INTEGER :: enPos
  INTEGER :: RhoVPos
  INTEGER :: RhoCPos
  INTEGER :: RhoRPos
  INTEGER :: RhoIPos
  INTEGER :: RhoSPos
  INTEGER :: nvPos
  INTEGER :: ncPos
  INTEGER :: nrPos
  INTEGER :: niPos
  INTEGER :: nsPos
  INTEGER :: prePos
  INTEGER :: rhoPos
  INTEGER :: tracer1Pos
  INTEGER :: tracer2Pos
  CHARACTER*30 :: Environment='Init'
  INTEGER :: uPosEnv
  INTEGER :: vPosEnv
  INTEGER :: wPosEnv
  INTEGER :: tkePosEnv
  INTEGER :: disPosEnv
  INTEGER :: omePosEnv
  INTEGER :: thPosEnv
  INTEGER :: enPosEnv
  INTEGER :: RhoPosEnv
  INTEGER :: RhoVPosEnv
  INTEGER :: RhoCPosEnv
  INTEGER :: RhoRPosEnv
  INTEGER :: RhoIPosEnv
  INTEGER :: RhoSPosEnv
  INTEGER :: nvPosEnv
  INTEGER :: ncPosEnv
  INTEGER :: nrPosEnv
  INTEGER :: niPosEnv
  INTEGER :: nsPosEnv
  INTEGER :: tracer1PosEnv
  INTEGER :: tracer2PosEnv
  INTEGER :: uPosLJac,uPosRJac
  INTEGER :: vPosLJac,vPosRJac
  INTEGER :: wPosLJac,wPosRJac
! Positionen for turbulence 
! Tke closure 
  INTEGER :: tkePosJac
  INTEGER :: disPosJac
  INTEGER :: omePosJac
! Closure nach Herzog
  INTEGER :: tkeHPosJac
  INTEGER :: tkeVPosJac
  INTEGER :: LenPosJac
  INTEGER :: thPosJac
  INTEGER :: enPosJac
  INTEGER :: RhoVPosJac
  INTEGER :: RhoCPosJac
  INTEGER :: RhoRPosJac
  INTEGER :: RhoIPosJac
  INTEGER :: RhoSPosJac
  INTEGER :: nvPosJac
  INTEGER :: ncPosJac
  INTEGER :: nrPosJac
  INTEGER :: niPosJac
  INTEGER :: nsPosJac
  INTEGER :: prePosJac
  INTEGER :: rhoPosJac
  INTEGER :: tracer1PosJac
  INTEGER :: tracer2PosJac
  INTEGER, ALLOCATABLE :: IndexMet(:,:)
  INTEGER, ALLOCATABLE :: IndexTher(:,:)

! Zeitsteuerung
  LOGICAL :: ErrControl
  LOGICAL :: TimeStepControl=.FALSE.
  REAL(RealKind) :: TimeAct
  REAL(RealKind) :: dt
  REAL(RealKind) :: dtP
  REAL(RealKind) :: dtMax
  REAL(RealKind) :: dtStart
  REAL(RealKind) :: dtOld
  REAL(RealKind), ALLOCATABLE :: RTol(:)
  REAL(RealKind), ALLOCATABLE :: ATol(:)
  REAL(RealKind) :: RTolG
  REAL(RealKind) :: ATolG
  REAL(RealKind) :: uRTol
  REAL(RealKind) :: uATol
  REAL(RealKind) :: vRTol
  REAL(RealKind) :: vATol
  REAL(RealKind) :: wRTol
  REAL(RealKind) :: wATol
  REAL(RealKind) :: thRTol
  REAL(RealKind) :: thATol
  REAL(RealKind) :: enRTol
  REAL(RealKind) :: enATol
  REAL(RealKind) :: tkeRTol
  REAL(RealKind) :: tkeATol
  REAL(RealKind) :: disRTol
  REAL(RealKind) :: disATol
  REAL(RealKind) :: tkeHRTol
  REAL(RealKind) :: tkeHATol
  REAL(RealKind) :: tkeVRTol
  REAL(RealKind) :: tkeVATol
  REAL(RealKind) :: LenRTol
  REAL(RealKind) :: LenATol
  REAL(RealKind) :: RhoVRTol
  REAL(RealKind) :: QVRTol
  REAL(RealKind) :: QVATol
  REAL(RealKind) :: QCRTol
  REAL(RealKind) :: QCATol
  REAL(RealKind) :: RhoVATol
  REAL(RealKind) :: QRRTol
  REAL(RealKind) :: QRATol
  REAL(RealKind) :: RhoCRTol
  REAL(RealKind) :: RhoIRTol
  REAL(RealKind) :: RhoCATol
  REAL(RealKind) :: RhoIATol
  REAL(RealKind) :: RhoRRTol
  REAL(RealKind) :: RhoRATol
  REAL(RealKind) :: RhoRTol
  REAL(RealKind) :: RhoATol
  REAL(RealKind) :: PreRTol
  REAL(RealKind) :: PreATol
  INTEGER :: EndIter

! Ausgabe
  LOGICAL :: FillValue=.FALSE.
! Rosenbrock-Method
  REAL(RealKind) :: beta0=1.0e0

! Anzahl Komponenten 
  INTEGER :: VectorComponentsM
  INTEGER :: VectorComponentsME
  INTEGER :: VectorComponentsT
  INTEGER :: VectorComponentsMet
  INTEGER :: VectorComponentsChem
  INTEGER, ALLOCATABLE :: PosE2Pos(:)
  INTEGER :: xOrder,yOrder,zOrder
  INTEGER :: xCoarse,yCoarse,zCoarse
  NAMELIST /ModelControl/ CGMaxIterPre,    &
                          CGTolPre,        &
                          CGMaxIterSch,    &
                          CGTolSch,        &
                          BiCGStabTol,     &
                          BiCGStabMaxIter, &
                          QMRTol,          &
                          QMRMaxIter,      &
                          Fac,             &
                          JacTransport,    &
                          JacTransportX,   &
                          JacTransportY,   &
                          JacTransportZ,   &
                          JacPartial,      &
                          JacSound,        &
                          MultiEx,         &
                          MultiMu,         &
                          MultiMuTR,       &
                          MultiTriT,       &
                          MultiTriTB,      &
                          MultiTriTR,      &
                          GaussIter2,      &
                          GaussIter3,      &
                          MultIter2,       &
                          MultIter3,       &
                          PreFac2,         &
                          PreFac3,         &
                          PreFac5,         &
                          ShiftLoc,           &
                          Debug,           &
                          MethSound,       &
                          CoarsePre,       &
                          StepGMRES,       &
                          RefLevelG,       &
                          PD,       &
                          mi_wgts,         &
                          ProfileOut,      &
                          ProfileOutGnu,   &
                          Method,          &
                          InitPeer,        &
                          MethAdv,          &
                          nsMin,              &
                          nsCompute,              &
                          GammaDiv,        &
                          LinearPressure,  &
                          IntFast,         &
                          CloudFast,       &
                          GradFull,        &
                          Primitive,       &
                          AdvPr,       &
                          Environment, &
                          VolumeCheck,     &
                          PrintNameLists,     &
                          dtMax,           &
                          dtStart,         &
                          dtP,             &
                          TimeStepControl, &
                          ref_glob,        &
                          ivar,            &
                          xOrder,          &
                          yOrder,          &
                          zOrder,          &
                          xCoarse,         &
                          yCoarse,         &
                          zCoarse,         &
                          TermUnit,        &
                          EndTime,         &
                          EndIter,         &
                          StartDay,        &
                          Year,            &
                          Month,           &
                          Day,             &
                          StartTime,       &
                          RestartTimeIncr, &
                          FillValue,       &
                          lng,             &
                          lat
  NAMELIST /ErrorControl/ ErrControl, &
                          RTolG,    &
                          ATolG,    &
                          uRtol,    &
                          uATol,    &
                          vRTol,    &
                          vATol,    &
                          wRTol,    &
                          wATol,    &
                          thRTol,   &
                          thATol,   &
                          enRTol,   &
                          enATol,   &
                          tkeRTol,  &
                          tkeATol,  &
                          disRTol,  &
                          disATol,  &
                          tkeHRTol,  &
                          tkeHATol,  &
                          tkeVRTol,  &
                          tkeVATol,  &
                          LenRTol,  &
                          LenATol,  &                          
                          QVRTol,   &
                          QVATol,   &
                          QCRTol,   &
                          QCATol,   &
                          QRRTol,   &
                          QRATol,   &
                          RhoVRTol,   &
                          RhoVATol,   &
                          RhoCRTol,   &
                          RhoCATol,   &
                          RhoIRTol,   &
                          RhoIATol,   &
                          RhoRRTol,   &
                          RhoRATol,   &
                          RhoRTol,   &
                          RhoATol,   &
                          PreRTol,   &
                          PreATol

  NAMELIST /ModelPos/ uPosL,  &
                      uPosR,  &
                      vPosL,  &
                      vPosR,  &
                      wPosL,  &
                      wPosR,  &
                      thPos,  &
                      enPos,  &
                      tkePos, &
                      disPos, &
                      omePos, &
                      tkeHPos, &
                      tkeVPos, &
                      LenPos, &
                      RhoVPos,  &
                      RhoCPos,  &
                      RhoRPos,  &
                      RhoIPos,  &
                      RhoSPos,  &
                      nvPos,  &
                      ncPos,  &
                      nrPos,  &
                      niPos,  &
                      nsPos,  &
                      tracer1Pos,  &
                      tracer2Pos,  &
                      prePos,   &
                      prePos,   &
                      prePos,   &
                      rhoPos

  NAMELIST /EnvPos/ uPosEnv,   &
                    vPosEnv,   &
                    wPosEnv,   &
                    thPosEnv,  &
!                    enPosEnv,  &
                    tkePosEnv, &
                    disPosEnv, &
                    RhoPosEnv,  &
                    RhoVPosEnv,  &
                    RhoCPosEnv,  &
                    RhoRPosEnv,  &
                    RhoIPosEnv,  &
                    RhoSPosEnv,  &
                    nvPosEnv,  &
                    ncPosEnv,  &
                    nrPosEnv,  &
                    niPosEnv,  &
                    nsPosEnv,  &
                    tracer1PosEnv,  &
                    tracer2PosEnv

                           

CONTAINS

SUBROUTINE InputModelPos(FileName)

  CHARACTER(*) :: FileName
  
  INTEGER :: Pos
  CHARACTER(300) :: Line

  uPosL=0
  uPosR=0
  vPosL=0
  vPosR=0
  wPosL=0
  wPosR=0
  thPos=0
  enPos=0
  tkePos=0
  disPos=0
  omePos=0
  tkeHPos=0
  tkeVPos=0
  LenPos=0
  RhoVPos=0
  RhoCPos=0
  RhoRPos=0
  RhoIPos=0
  RhoSPos=0
  nvPos=0
  ncPos=0
  nrPos=0
  niPos=0
  nsPos=0
  prePos=0
  rhoPos=0
  tracer1Pos=0
  tracer2Pos=0
! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelPos')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelPos)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

  VectorComponentsM=MAX(uPosL,uPosR,vPosL,vPosR,wPosL,wposR,thPos,enPos,tkePos, &
                        disPos,omePos,tkeHPos,tkeVPos,LenPos,                          &
                        RhoVPos,RhoCPos,RhoRPos,RhoIPos,RhoSPos,nvPos,ncPos,nrPos,niPos,nsPos,rhoPos,tracer1Pos,tracer2Pos)

  IF (MyId==0.AND.PrintNameLists) THEN 
    WRITE(TermUnit,*) ' uPosL = ', uPosL
    WRITE(TermUnit,*) ' uPosR = ', uPosR
    WRITE(TermUnit,*) ' vPosL = ', vPosL
    WRITE(TermUnit,*) ' vPosR = ', vPosR
    WRITE(TermUnit,*) ' wPosL = ', wPosL
    WRITE(TermUnit,*) ' wPosR = ', wPosR
    WRITE(TermUnit,*) ' thPos = ', thPos
    WRITE(TermUnit,*) ' enPos = ', enPos
    WRITE(TermUnit,*) 'tracer1Pos = ',tracer1Pos
    WRITE(TermUnit,*) 'tracer2Pos = ',tracer2Pos
    WRITE(TermUnit,*) 'tkePos = ',tkePos
    WRITE(TermUnit,*) 'disPos = ',disPos
    WRITE(TermUnit,*) 'omePos = ',omePos
    WRITE(TermUnit,*) 'tkeHPos = ',tkeHPos
    WRITE(TermUnit,*) 'tkeVPos = ',tkeVPos
    WRITE(TermUnit,*) 'LenPos = ',LenPos  
    WRITE(TermUnit,*) ' RhoVPos = ', RhoVPos
    WRITE(TermUnit,*) ' RhoCPos = ', RhoCPos
    WRITE(TermUnit,*) ' RhoRPos = ', RhoRPos
    WRITE(TermUnit,*) ' RhoIPos = ', RhoIPos
    WRITE(TermUnit,*) ' RhoSPos = ', RhoSPos
    WRITE(TermUnit,*) ' nvPos = ', nvPos
    WRITE(TermUnit,*) ' ncPos = ', ncPos
    WRITE(TermUnit,*) ' nrPos = ', nrPos
    WRITE(TermUnit,*) ' niPos = ', niPos
    WRITE(TermUnit,*) ' nsPos = ', nsPos
    WRITE(TermUnit,*) 'prePos = ',prePos
    WRITE(TermUnit,*) 'rhoPos = ',rhoPos
    WRITE(TermUnit,*)
  END IF

END SUBROUTINE InputModelPos

SUBROUTINE InputEnvPos(FileName)

  CHARACTER(*) :: FileName
  
  INTEGER :: Pos
  CHARACTER(300) :: Line

  uPosEnv=0
  vPosEnv=0
  wPosEnv=0
  thPosEnv=0
!  enPosEnv=0
  tkePosEnv=0
  DisPosEnv=0
  RhoPosEnv=0
  RhoVPosEnv=0
  RhoCPosEnv=0
  RhoRPosEnv=0
  RhoIPosEnv=0
  RhoSPosEnv=0
  nvPosEnv=0
  ncPosEnv=0
  nrPosEnv=0
  niPosEnv=0
  nsPosEnv=0
  tracer1PosEnv=0
  tracer2PosEnv=0
! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&EnvPos')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=EnvPos)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

  VectorComponentsME=MAX(uPosEnv,vPosEnv,wPosEnv,thPosEnv,tkePosEnv, &
                         disPosEnv,omePosEnv,RhoVPosEnv,RhoPosEnv,RhoCPosEnv,RhoRPosEnv,RhoIPosEnv,RhoSPosEnv, &
                         nvPosEnv,ncPosEnv,nrPosEnv,niPosEnv,nsPosEnv,tracer1PosEnv,tracer2PosEnv)
  ALLOCATE(PosE2Pos(VectorComponentsME))
 
  IF (MyID==0.AND.PrintNameLists) THEN
    WRITE(TermUnit,*) '  uPosEnv = ',  uPosEnv
    WRITE(TermUnit,*) '  vPosEnv = ',  vPosEnv
    WRITE(TermUnit,*) '  wPosEnv = ',  wPosEnv
    WRITE(TermUnit,*) ' thPosEnv = ', thPosEnv
!   WRITE(TermUnit,*) ' enPosEnv = ', enPosEnv
    WRITE(TermUnit,*) 'tracer1PosEnv = ',tracer1PosEnv
    WRITE(TermUnit,*) 'tracer2PosEnv = ',tracer2PosEnv
    WRITE(TermUnit,*) 'tkePosEnv = ',tkePosEnv
    WRITE(TermUnit,*) 'disPosEnv = ',disPosEnv
    WRITE(TermUnit,*) 'omePosEnv = ',omePosEnv
    WRITE(TermUnit,*) ' RhoVPosEnv = ', RhoVPosEnv
    WRITE(TermUnit,*) ' RhoPosEnv = ', RhoPosEnv
    WRITE(TermUnit,*) ' RhoCPosEnv = ', RhoCPosEnv
    WRITE(TermUnit,*) ' RhoRPosEnv = ', RhoRPosEnv
    WRITE(TermUnit,*) ' RhoIPosEnv = ', RhoIPosEnv
    WRITE(TermUnit,*) ' RhoSPosEnv = ', RhoSPosEnv
    WRITE(TermUnit,*) ' nvPosEnv = ', nvPosEnv
    WRITE(TermUnit,*) ' ncPosEnv = ', ncPosEnv
    WRITE(TermUnit,*) ' nrPosEnv = ', nrPosEnv
    WRITE(TermUnit,*) ' niPosEnv = ', niPosEnv
    WRITE(TermUnit,*) ' nsPosEnv = ', nsPosEnv
    WRITE(TermUnit,*)
  END IF


END SUBROUTINE InputEnvPos

SUBROUTINE InputModelControl(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: Pos
  CHARACTER(300) :: Line
  
  CGMaxIterPre=20
  CGTolPre=1.d-4
  CGMaxIterSch=20
  CGTolSch=1.d-4
  BiCGStabTol=1.d-6
  BiCGStabMaxIter=30
  QMRTol=1.d-6
  QMRMaxIter=10
  RefLevelG=3
  PD=.TRUE.
  mi_wgts=''
  ProfileOut='ProfileOut'
  ProfileOutGnu='ProfileOutGnu'
  Method='Ros3'
  InitPeer=.TRUE.
  MethAdv='Koren'
  ns=12
  nsMin=12
  nsCompute=.FALSE.
  GammaDiv=Zero
  dtMax=1.0d0
  dtStart=1.0d0
  dtP=1.0d+3
  ref_glob=0
  ivar=0        ! Partition variant (0: METIS)
  xOrder=2
  yOrder=3
  zOrder=1
  xCoarse=2
  yCoarse=2
  zCoarse=2
  EndTime=Zero
  RestartTimeIncr=1.0e20_RealKind
  EndIter=10000
  StartDay=0
  StartTime=0
  lng=12.0d0
  lat=51.0d0

! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelControl')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelControl)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
! WRITE(*,*) 'ivar in Control',ivar,MyId
  IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
    GravComp=Grav
  ELSE
    GravComp=Zero
  END IF
  IF (MultiEx) THEN
    xOrder=1
    yOrder=2
    zOrder=3
  END IF  
  lng=lng*Pi/180.0d0
  lat=lat*Pi/180.0d0

END SUBROUTINE InputModelControl

SUBROUTINE InputErrorControl(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: nc,Pos
  CHARACTER(300) :: Line

  ErrControl=.TRUE.
! RTol = relative change of variables (domain maximum or average) in dt
! ATol = absolute change of variables (domain maximum or average) in dt
  RTolG=1.d-3 ! for uncpecified variables
  ATolG=1.d-5 ! for uncpecified variables
  uRTol=1.d-2
  uATol=1.d-5
  vRTol=1.d-2
  vATol=1.d-5
  wRTol=1.d-2
  wATol=1.d-5
  thRTol=1.d-2
  thATol=1.d-5
  tkeRTol=1.d-2
  tkeATol=1.d-5
  disRTol=1.d-2
  disATol=1.d-5
  tkeHRTol=1.d-2
  tkeHATol=1.d-5
  tkeVRTol=1.d-2
  tkeVATol=1.d-5
  LenRTol=1.d-2
  LenATol=1.d-5
  RhoVRTol=1.d-2
  QVRTol=1.d-2
  QVATol=1.d-2
  QCRTol=1.d-2
  QCATol=1.d-2
  QRRTol=1.d-2
  QRATol=1.d-2
  RhoVATol=1.d-5
  RhoCRTol=1.d-2
  RhoCATol=1.d-5
  RhoIRTol=1.d-2
  RhoIATol=1.d-5
  RhoRRTol=1.d-2
  RhoRATol=1.d-5

! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ErrorControl')>0) THEN
      BACKSPACE(InputUnit)
  WRITE(*,*) 'ErrControl  1',ErrControl
      READ(InputUnit,NML=ErrorControl)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
  WRITE(*,*) 'ErrControl  ',ErrControl
  ALLOCATE(RTol(VectorComponentsT))
  ALLOCATE(ATol(VectorComponentsT))
  RTol=RTolG
  IF (MyID==0) WRITE(*,*) 'ATolG',ATolG
  ATol=ATolG
  IF (uPosL>0) THEN
    RTol(uPosL)=uRTol
    ATol(uPosL)=uATol
  END IF
  IF (uPosR>0) THEN
    RTol(uPosR)=uRTol
    ATol(uPosR)=uATol
  END IF
  IF (vPosL>0) THEN
    RTol(vPosL)=vRTol
    ATol(vPosL)=vATol
  END IF
  IF (vPosR>0) THEN
    RTol(vPosR)=vRTol
    ATol(vPosR)=vATol
  END IF
  IF (wPosL>0) THEN
    RTol(wPosL)=wRTol
    ATol(wPosL)=wATol
  END IF
  IF (wPosR>0) THEN
    RTol(wPosR)=wRTol
    ATol(wPosR)=wATol
  END IF
  IF (thPos>0) THEN
    RTol(thPos)=thRTol
    ATol(thPos)=thATol
  END IF
  IF (enPos>0) THEN
    RTol(enPos)=enRTol
    ATol(enPos)=enATol
  END IF
  IF (tkePos>0) THEN
    RTol(tkePos)=tkeRTol
    ATol(tkePos)=tkeATol
  END IF
  IF (disPos>0) THEN
    RTol(disPos)=disRTol
    ATol(disPos)=disATol
  END IF
  IF (tkeHPos>0) THEN
    RTol(tkeHPos)=tkeHRTol
    ATol(tkeHPos)=tkeHATol
  END IF
  IF (tkeVPos>0) THEN
    RTol(tkeVPos)=tkeVRTol
    ATol(tkeVPos)=tkeVATol
  END IF
  IF (LenPos>0) THEN
    RTol(LenPos)=LenRTol
    ATol(LenPos)=LenATol
  END IF  
  IF (RhoVPos>0) THEN
    RTol(RhoVPos)=RhoVRTol
    ATol(RhoVPos)=RhoVATol
  END IF
  IF (RhoCPos>0) THEN
    RTol(RhoCPos)=RhoCRTol
    ATol(RhoCPos)=RhoCATol
  END IF
  IF (RhoIPos>0) THEN
    RTol(RhoIPos)=RhoIRTol
    ATol(RhoIPos)=RhoIATol
  END IF
  IF (RhoRPos>0) THEN
    RTol(RhoRPos)=RhoRRTol
    ATol(RhoRPos)=RhoRATol
  END IF
  IF (RhoPos>0) THEN
    RTol(RhoPos)=RhoRTol
    ATol(RhoPos)=RhoATol
  END IF
  IF (MyID==0) WRITE(*,*) 'ATOL',ATol(1),RTol(1)

END SUBROUTINE InputErrorControl


END MODULE Control_Mod
