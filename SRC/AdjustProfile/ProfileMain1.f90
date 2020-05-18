PROGRAM Profile
  USE Profile_Mod
  IMPLICIT NONE

  INTEGER  INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
  REAL(RealKind) :: &
      T, Y(90), YPRIME(90), TOUT, RWORK(100000), &
      RPAR(1000), RTOL, ATOL
  REAL(RealKind) :: Theta0,ystart,tout0
  REAL(RealKind) :: Te,ThetaRho,p,e,qv,qt,RhoD,RhoV,Rho,RhoGes,Const0,Gamm0,KappaM,dThdz,dRHdz
  REAL(RealKind) :: wSave(9,2000)
  REAL(RealKind) :: TeStart     = 288.0d0  
  REAL(RealKind) :: PreStart    = 1013.25d2
  REAL(RealKind) :: RelHumStart = 0.00d0
  REAL(RealKind) :: RelHumInv   = 0.90d0
  REAL(RealKind) :: InvHeight   = 1500.0d0
  REAL(RealKind) :: N_d         = 0.01d0 
  REAL(RealKind) :: DeltaZ      = 5.0d0
  REAL(RealKind), ALLOCATABLE :: c(:,:)
  REAL(RealKind) :: zu,zo,Tu,To
  REAL(RealKind) :: z,zOut
  INTEGER :: iHeight,iz
  LOGICAL :: Inversion=.FALSE.
  LOGICAL :: LinRH=.FALSE.
  LOGICAL :: ConstRH=.TRUE.
  LOGICAL :: Interpol=.TRUE.
  INTEGER :: stat
  INTEGER :: Pos
  INTEGER :: CountSteps=0
  CHARACTER(8) :: ProfType
  CHARACTER(128) :: Filename
  CHARACTER(64) :: par

  INTEGER :: i,nwSave
  INTEGER :: DryPos,VapPos,LiqPos
  REAL(RealKind) :: T1


  WRITE(*,*) "use: ./run outputfilename"
  CALL get_command_argument(1,FileName)
  CALL ComputeParameter
  CALL ReadInput(FileName)
  iPar(1)=CaseProfile

  t=0.0
  info(1:15)=0
  lrw=100000
  liw=100000


  PreStart=Pre(1)
  TeStart=Temp(1)
  RelHumStart=RH(1)

  CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=RelHumStart)
  ALLOCATE(c(1:nz,1:neq))
  DO i=1,nz
    WRITE(*,*) 'zM(i)',i,zM(i)
  END DO  

  rtol = 1.0d-4
  atol = 1.0d-8

  tout=1.0d0
  nwSave=0
  idid=0
  yprime=0.0d0
  y=0.0d0
  IF (ipar(1)==2) THEN
    rpar(5)=0.1
    rpar(6)=rpar(5)
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    DO  
      rwork=0.0d0
      info=0.0d0
      INFO(11)=1 ! DDASSL computes consisten initial values
      t=0.0d0
      tout=1.0d0
      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      IF (rpar(5)==RH(1)) THEN
        EXIT
      END IF  
      rpar(5)=MIN(rpar(5)+0.1d0,RH(1))
      rpar(6)=rpar(5)
      CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0.0d0
    info(11)=1 ! DDASSL computes consisten initial values
    t=0.0d0
    tout=1.0d0
    CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF  

  iHeight=1
  iz=1
  ! Forward
  WRITE(*,*) 'Forward integration'
  S1:DO 
    zu=Height(iHeight)
    rpar(1)=zu
    zo=Height(iHeight+1)
    rpar(2)=zo
    z=zu
    rpar(3)=Temp(iHeight)
    rpar(4)=Temp(iHeight+1)
    rpar(5)=RH(iHeight)
    rpar(6)=RH(iHeight+1)
    c(1,:)=y(1:neq)
    S2:DO 
      IF (iz>nz) THEN
        EXIT S1
      END IF  
      IF (zu<zM(iz).AND.zM(iz)<=zO) THEN
        zOut=zM(iz)
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        iz=iz+1
      ELSE IF (zM(iz)<=zu) THEN
        iz=iz+1
      ELSE 
        zOut=zo
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        iHeight=iHeight+1
        EXIT S2
      END IF  
    END DO S2  
  END DO S1 

  ! Rückwärts
  yprime=0.0d0
  y=0.0d0
  IF (ipar(1)==2) THEN
    rpar(5)=0.1
    rpar(6)=rpar(5)
    y(6)=rpar(5)
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    DO  
      rwork=0.0d0
      info=0.0d0
      INFO(11)=1 ! DDASSL computes consisten initial values
      t=0.0d0
      tout=1.0d0
      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      IF (rpar(5)==RH(1)) THEN
        EXIT
      END IF  
      rpar(5)=MIN(rpar(5)+0.1d0,RH(1))
      rpar(6)=rpar(5)
      CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0.0d0
    info(11)=1 ! DDASSL computes consisten initial values
    t=0.0d0
    tout=1.0d0
    CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF  
  rwork=0.0d0
  info=0
  S3:DO 
    zu=Height(1)
    rpar(1)=zu
    zo=Height(2)
    rpar(2)=zo
    z=zu
    rpar(3)=Temp(1)
    rpar(4)=Temp(2)
    rpar(5)=RH(1)
    rpar(6)=RH(2)
    iz=nz
    S4:DO 
      IF (iz<1) THEN
        EXIT S3
      END IF  
      IF (zu>zM(iz)) THEN
        zOut=zM(iz)
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        iz=iz-1
      ELSE IF (zM(iz)>=zu) THEN
        iz=iz-1
      END IF  
    END DO S4  
  END DO S3 

  WRITE(*,*) 'Vor Output',nz
  CALL OutputRes(nz,zM,c,FileNameOut,iPar(1))
  CALL OutputResWind(NumValues,Height,U,V,FileNameOut)

      
END PROGRAM Profile
