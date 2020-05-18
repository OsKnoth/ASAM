MODULE Activity_Mod

  USE ActivityAim_Mod
  USE ActivityAimo_Mod
  USE ActivityPitzer_Mod
  USE ActivityPitzer1_Mod

  IMPLICIT NONE
 

CONTAINS

SUBROUTINE ComputeActivity(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

  IF (MethodAct=='Aim') THEN
    CALL ComputeActivityAim(cAqua,ActCoeff,TAbs)
  ELSE IF (MethodAct=='Aimo') THEN
    CALL ComputeActivityAimo(cAqua,ActCoeff,TAbs)
  ELSE IF (MethodAct=='Pit') THEN
    CALL ComputeActivityPitzer(cAqua,ActCoeff)
  ELSE IF (MethodAct=='Pit1') THEN
    CALL ComputeActivityPitzer1(cAqua,ActCoeff)
  ELSE
    ActCoeff=One
  END IF
END SUBROUTINE ComputeActivity  

SUBROUTINE OutputActivity(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

  IF (MethodAct=='Aim') THEN
    CALL OutputActivityAim(cAqua,ActCoeff,TAbs)
  ELSE IF (MethodAct=='Pit') THEN
    CALL OutputActivityPitzer(cAqua,ActCoeff)
  ELSE IF (MethodAct=='Pit1') THEN
    CALL OutputActivityPitzer1(cAqua,ActCoeff)
  ELSE
    ActCoeff=One
  END IF
END SUBROUTINE OutputActivity

SUBROUTINE InitActivity
  CALL AllocateVec4Chemie(Act,VectorComponentsM)
  Act=One
  IF (MethodAct=='Aim') THEN
    CALL InitActivityAim
  ELSE IF (MethodAct=='Aimo') THEN
    CALL InitActivityAimo
  ELSE IF (MethodAct=='Pit') THEN
    CALL InitActivityPitzer
  ELSE IF (MethodAct=='Pit1') THEN
    CALL InitActivityPitzer1
  END IF
END SUBROUTINE InitActivity

END MODULE Activity_Mod
