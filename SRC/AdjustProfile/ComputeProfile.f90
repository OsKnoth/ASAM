SUBROUTINE ComputeProfile(y,yPrime,c,Height,zM,Temp,RH)

  USE Kind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: y(:)
  REAL(RealKind) :: yPrime(:)
  REAL(RealKind) :: c(:,:)
  REAL(RealKind) :: Height(:)
  REAL(RealKind) :: zM(:)
  REAL(RealKind), OPTIONAL :: Temp(:)
  REAL(RealKind), OPTIONAL :: RH(:)

  INTEGER  INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
  REAL(RealKind) :: T, TOUT, RWORK(100000), &
      RPAR(1000), RTOL, ATOL

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
        IF (zM(iz)==zO) THEN
          iHeight=iHeight+1
          iz=iz+1
          EXIT S2
        END IF 
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

  ! Backward
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
END SUBROUTINE ComputeProfile
