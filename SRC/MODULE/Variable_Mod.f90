MODULE Variable_Mod
  USE Kind_Mod
  USE Control_Mod
  IMPLICIT NONE

  TYPE Variable_T
    CHARACTER*20 :: Name
    CHARACTER*20 :: Type
    CHARACTER*20 :: Unit
    LOGICAL :: ScaleRho=.TRUE.
    REAL(RealKind) :: Fac
  END TYPE Variable_T

  INTEGER :: NumberVariables
  TYPE(Variable_T), ALLOCATABLE :: Variables(:)
  CHARACTER*20, ALLOCATABLE :: SpeciesName(:)

CONTAINS  

SUBROUTINE SetMetVariables
  IF (Position('UCL')>0) THEN
    uPosL=Position('UCL')
    Variables(uPosL)%Name='UCL'
    Variables(uPosL)%Type='Velocity'
    Variables(uPosL)%Unit='m/s'
  END IF
  IF (Position('UCR')>0) THEN
    uPosR=Position('UCR')
    Variables(uPosR)%Name='UCR'
    Variables(uPosR)%Type='Velocity'
    Variables(uPosR)%Unit='m/s'
  END IF
  IF (Position('VCL')>0) THEN
    vPosL=Position('VCL')
    Variables(vPosL)%Name='VCL'
    Variables(vPosL)%Type='Velocity'
    Variables(vPosL)%Unit='m/s'
  END IF
  IF (Position('VCR')>0) THEN
    vPosR=Position('VCR')
    Variables(vPosR)%Name='VCR'
    Variables(vPosR)%Type='Velocity'
    Variables(vPosR)%Unit='m/s'
  END IF
  IF (Position('WCL')>0) THEN
    wPosL=Position('WCL')
  END IF
  IF (Position('WCR')>0) THEN
    wPosR=Position('WCR')
  END IF
  IF (Position('TE')>0) THEN
    thPos=Position('TE')
  END IF
  IF (Position('EN')>0) THEN
    EnPos=Position('EN')
  END IF
  IF (Position('TKE')>0) THEN
    tkePos=Position('TKE')
  END IF
  IF (Position('DIS')>0) THEN
    disPos=Position('DIS')
  END IF
  IF (Position('TKEH')>0) THEN
    tkeHPos=Position('TKEH')
  END IF
  IF (Position('TKEV')>0) THEN
    tkeVPos=Position('TKEV')
  END IF
  IF (Position('LEN')>0) THEN
    LenPos=Position('LEN')
  END IF
  IF (Position('RHO')>0) THEN
    rhoPos=Position('RHO')
  END IF
  IF (Position('PRE')>0) THEN
    prePos=Position('PRE')
  END IF
  IF (Position('RhoV')>0) THEN
    RhoVPos=Position('RhoV')
  END IF
  IF (Position('RhoC')>0) THEN
    RhoCPos=Position('RhoC')
  END IF
  IF (Position('RhoR')>0) THEN
    RhoRPos=Position('RhoR')
  END IF
  IF (Position('NV')>0) THEN
    nvPos=Position('NV')
  END IF
  IF (Position('NC')>0) THEN
    ncPos=Position('NC')
  END IF
  IF (Position('NR')>0) THEN
    nrPos=Position('NR')
  END IF
  IF (Position('TRACER1')>0) THEN
    tracer1Pos=Position('TRACER1')
  END IF
  IF (Position('TRACER2')>0) THEN
    tracer2Pos=Position('TRACER2')
  END IF
END SUBROUTINE SetMetVariables
FUNCTION Position(Species)

  INTEGER :: Position
  CHARACTER(*) :: Species

  INTEGER :: i

  Position=0
  DO i=LBOUND(SpeciesName,1),UBOUND(SpeciesName,1)
    IF (TRIM(Species)==TRIM(SpeciesName(i))) THEN
      Position=i
      EXIT
    END IF
  END DO

END FUNCTION 
END MODULE Variable_Mod
