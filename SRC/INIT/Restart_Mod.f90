MODULE Restart_Mod

  USE Parallel_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Output_Mod

  IMPLICIT NONE

  REAL(RealKind) :: RestartTime

CONTAINS

SUBROUTINE WriteRestart(VecMet,VecChem,Vel,Time,FileName)
  TYPE(Vector4Cell_T) :: VecMet(:)
  TYPE(Vector4Cell_T) :: VecChem(:)
  TYPE(VelocityFace_T) :: Vel(:)
  REAL(RealKind) :: Time
  CHARACTER(*) :: FileName

  INTEGER :: ic
  CHARACTER(10) :: iName
  CHARACTER(LEN(FileName)) :: FileNameLoc

  FileNameLoc=FileName(1:INDEX(FileName,'.grid'))
  IF (Time>RestartTime) THEN
    WRITE(iName,'(I8)') MyId
    OPEN(UNIT=OutputUnit,FILE='RESTART/'//TRIM(FileNameLoc)//ADJUSTL(iName),STATUS='REPLACE',FORM='UNFORMATTED')
    WRITE(OutputUnit) Time
    WRITE(OutputUnit) GetGMVStep()
    WRITE(OutputUnit) GetOutputTime()
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ic=1,UBOUND(VecMet(ibLoc)%Vec,1)
        WRITE(OutputUnit) VecMet(ibLoc)%Vec(ic)%c
      END DO  
      DO ic=1,UBOUND(VecChem(ibLoc)%Vec,1)
        WRITE(OutputUnit) VecChem(ibLoc)%Vec(ic)%c
      END DO  
      WRITE(OutputUnit) Vel(ibLoc)%uF,Vel(ibLoc)%vF,Vel(ibLoc)%wF
    END DO  
    CLOSE(OutputUnit)
    RestartTime=RestartTime+RestartTimeIncr
  END IF  

END SUBROUTINE WriteRestart

SUBROUTINE ReadRestart(VecMet,VecChem,Vel,Time,FileName)
  TYPE(Vector4Cell_T) :: VecMet(:)
  TYPE(Vector4Cell_T) :: VecChem(:)
  TYPE(VelocityFace_T) :: Vel(:)
  REAL(RealKind) :: Time
  CHARACTER(*) :: FileName

  INTEGER :: ic
  INTEGER :: GMVStep
  REAL(RealKind) :: OutputTime
  CHARACTER(10) :: iName
  CHARACTER(LEN(FileName)) :: FileNameLoc

  FileNameLoc=FileName(1:INDEX(FileName,'.grid'))

  IF (Restart) THEN
    WRITE(iName,'(I8)') MyId
    OPEN(UNIT=InputUnit,FILE='RESTART/'//TRIM(FileNameLoc)//ADJUSTL(iName),STATUS='OLD',FORM='UNFORMATTED')
    READ(InputUnit) Time
    READ(InputUnit) GMVStep
    CALL SetGMVStep(GMVStep)
    READ(InputUnit) OutputTime
    CALL SetOutputTime(OutputTime)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ic=1,UBOUND(VecMet(ibLoc)%Vec,1)
        READ(InputUnit) VecMet(ibLoc)%Vec(ic)%c
      END DO  
      DO ic=1,UBOUND(VecChem(ibLoc)%Vec,1)
        READ(InputUnit) VecChem(ibLoc)%Vec(ic)%c
      END DO  
      READ(InputUnit) Vel(ibLoc)%uF,Vel(ibLoc)%vF,Vel(ibLoc)%wF
    END DO  
    CLOSE(InputUnit)
  ELSE
    RestartTime=Time+RestartTimeIncr
  END IF
END SUBROUTINE ReadRestart
END MODULE Restart_Mod
