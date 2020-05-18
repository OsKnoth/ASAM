MODULE TestCoeff_Mod

  USE Chemie_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE TestCoeff

  INTEGER :: iPos
  REAL(RealKind) :: c1,c2,c3
  CHARACTER*300 :: Line
  CHARACTER*20 :: SpeciesName
  LOGICAL :: Back
  REAL(RealKind), ALLOCATABLE :: DH(:),DG(:),CP(:)
  REAL(RealKind) :: DHLoc,DGLoc,CPLoc
  REAL(RealKind) :: A,B,C
  REAL(RealKind) :: T1,T2
  REAL(RealKind) :: TRef=298.15d0

  CHARACTER(30) :: FileName

  ALLOCATE(DH(nSpc))
  ALLOCATE(DG(nSpc))
  ALLOCATE(CP(nSpc))
  FileName='Thermo1.dat'
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,'BEGIN_CHEM','BEGIN_DATA','END_DATA', &
                  Name1=SpeciesName,R1=c1,R2=c2,R3=c3)
    IF (Back) THEN
      EXIT
    END IF
    iPos=Position(SpeciesName)
    IF (iPos>0) THEN
      DH(iPos)=c1
      DG(iPos)=c2
      CP(iPos)=c3
    END IF
  END DO
  CALL CloseFile
  DHLoc= &
       -DH(Position('sNACL')) &
       +DH(Position('NAp')) &
       +DH(Position('CLm')) 
  DHLoc=DHLoc*1.d3
  DGLoc= &
       -DG(Position('sNACL')) &
       +DG(Position('NAp')) &
       +DG(Position('CLm')) 
  DGLoc=DGLoc*1.d3
  CPLoc= &
       -CP(Position('sNACL')) &
       +CP(Position('NAp')) &
       +CP(Position('CLm')) 
  A=EXP(-DGLoc/TRef/GASCONST)
  B=CPloc/GASCONST
  C=(-DHLoc+CPLoc*Tref)/GASCONST

  DHLoc= &
       -DH(Position('sNANO3')) &
       +DH(Position('NAp')) &
       +DH(Position('NO3m'))
  DHLoc=DHLoc*1.d3
  DGLoc= &
       -DG(Position('sNANO3')) &
       +DG(Position('NAp')) &
       +DG(Position('NO3m'))
  DGLoc=DGLoc*1.d3
  CPLoc= &
       -CP(Position('sNANO3')) &
       +CP(Position('NAp')) &
       +CP(Position('NO3m'))
  A=EXP(-DGLoc/TRef/GASCONST)
  B=CPloc/GASCONST
  C=(-DHLoc+CPLoc*Tref)/GASCONST


  
END SUBROUTINE TestCoeff
END MODULE TestCoeff_Mod

