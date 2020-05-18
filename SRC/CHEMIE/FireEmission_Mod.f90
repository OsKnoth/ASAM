MODULE FireEmission_Mod

  USE Kind_Mod
  USE EmissDeposParameter_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE SetFireEmission(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: i, j
  REAL(RealKind) :: c1
  CHARACTER(40) :: SpeciesName
  CHARACTER(20) :: S1,S2,End
  LOGICAL :: Back

! Input of aerosol emissions from fire
  S1='BEGIN_AEROSOL'
  S2='BEGIN_FIREEMISS'
  End='END_FIREEMISS'
  CALL OpenFile(FileName)
  NumFireAeroEmiss=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      NumFireAeroEmiss=NumFireAeroEmiss+1
    END IF
  END DO
  CALL CloseFile
  ALLOCATE(FireAeroEmiss(NumFireAeroEmiss))
  NumFireAeroEmiss=0

  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=c1)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      NumFireAeroEmiss=NumFireAeroEmiss+1
      FireAeroEmiss(NumFireAeroEmiss)%Pos=Position(SpeciesName)
      FireAeroEmiss(NumFireAeroEmiss)%Konz=c1
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE SetFireEmission

SUBROUTINE FireEmission(U10,FL,ix,iy,iz)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: U10
  REAL(RealKind) :: FL

  INTEGER :: i
  REAL(RealKind) :: EmissionWind
  REAL(RealKind) :: U_thres=6.5d0 ! 10m threshold velocity according Tegen_Fung
  REAL(RealKind) :: C=0.7E-6 ! dimensional factor C

  IF (U10.GT.U_thres) THEN
    EmissionWind=C*(U10-U_thres)*U10**2.d0
    DO i=1,NumFireAeroEmiss
      fVec(FireAeroEmiss(NumFireAeroEmiss)%Pos)%c(ix,iy,iz,:)=fVec(FireAeroEmiss(NumFireAeroEmiss)%Pos)%c(ix,iy,iz,:) &
         +FL/(VolC(ix,iy,iz)+Eps)*EmissionWind
    END DO     
  END IF
END SUBROUTINE FireEmission

END MODULE FireEmission_Mod
