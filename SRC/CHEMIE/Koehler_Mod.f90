MODULE Koehler_Mod

  USE MicroPhysics_Mod
  IMPLICIT NONE

CONTAINS

SUBROUTINE InvSaturation(m,S,Mass,tAbs)

  REAL(RealKind) :: m,S,Mass(:),tAbs

  INTEGER :: Iter
  REAL(RealKind) :: alpha,beta
  REAL(RealKind) :: Rad
  REAL(RealKind) :: Fr,FrP
  REAL(RealKind), PARAMETER ::  TolFr=1.d-20 
  INTEGER, PARAMETER ::  MaxIter=20


! F(r)=alpha+beta/m

  Iter=0
  DO
    Iter=Iter+1
    Mass(iWater)=m
    Rad=Radius(Mass) ! Radius(iwater)
    Fr=Kohler1(S+One,Mass,Rad,tAbs) ! ??? Formel dahinter nachvollzogen, ergebnislos
    IF (ABS(Fr)<TolFr) THEN
      Fr=Kohler1(S+One,Mass,Rad,tAbs)
      EXIT
    END IF 
    IF (Iter>MaxIter) THEN
      Fr=Kohler1(S+One,Mass,Rad,tAbs)
      EXIT
    END IF
    FrP=Kohler1P(Mass,Rad,tAbs)
    beta=-FrP*m*m
    alpha=Fr-beta/m
    IF (-beta/alpha<=0.0d0) THEN
      IF (Fr>0.0d0) THEN
        m=2.0d0*m
      ELSE
        m=0.5d0*m
      END IF
    ELSE
      m=-(beta/alpha)
    END IF
  END DO

END SUBROUTINE InvSaturation

SUBROUTINE InitDroplet(nW,mW,cW &
                      ,nD,mD,NumD,MassD,S,tAbs)

  INTEGER :: nW,nD
  REAL(RealKind) :: mW(nW+1),cW(:,:)
  REAL(RealKind) :: mD(nD),NumD(nD),MassD(:,:),S,tAbs

  INTEGER :: iD,iW
  REAL(RealKind) :: Mass(Size(MassD,2))
  REAL(RealKind) :: m,mGes
  
  cW=0.0d0
!  m=0.5d0*(mD(1)+mD(2))
  DO iD=1,nD
    m=mD(iD)
    IF (NumD(iD)>0.0d0) THEN
      MassD(iD,iNC)=NumD(iD)
      Mass(:)=MassD(iD,:)/NumD(iD)
      IF (iWater>0) THEN
        CALL InvSaturation(m,S,Mass,tAbs)
        Mass(iWater)=m ! Wassermasser aus Routine
      END IF
      mGes=TotalMass(Mass)
      DO iW=1,nW
        IF (mW(iW)<=mGes.AND.mGes<mW(iW+1)) THEN
          cW(iW,:)=cW(iW,:)+NumD(iD)*Mass
        END IF
      END DO
    END IF
  END DO
  IF (iWater>0) THEN
    DO iw=1,nW
      IF (cW(iW,iNc)>Zero) THEN
        Mass(:)=cW(iW,:)/cW(iW,iNc) 
        m=Mass(iWater)
        CALL InvSaturation(m,S,Mass,tAbs)
        cW(iW,iWater)=Mass(iWater)*cW(iW,iNc) ! Mass(iWater) durch CALL nicht verändert!!!! In CALL wird m neu berechnet
      END IF
    END DO
  END IF
  IF (iRelax>0) THEN
    cW(:,iRelax)=Zero
  END IF  
END SUBROUTINE InitDroplet

END MODULE Koehler_Mod

