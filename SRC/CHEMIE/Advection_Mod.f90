MODULE Advection_Mod

  USE Kind_Mod
  USE DataType_Mod
  USE Microphysics_Mod
  USE Chemie_Mod 
  USE Ansatz_Mod
  USE InitAerosol_Mod

  IMPLICIT NONE

  REAL(RealKind), PRIVATE, PARAMETER :: EpsQuo=1.d-99

CONTAINS

SUBROUTINE Advection(c,fc,x,HenryTrans,Time)

  TYPE(Vec4_T) :: c(:),fc(:)
  REAL(RealKind) :: x(:)
  LOGICAL :: HenryTrans
  REAL(RealKind) :: Time
  REAL(RealKind) :: Kel,Rao

  INTEGER :: i,is,ig,iReak,l
  REAL(RealKind) :: nC,ImLGes,ImUGes,ImGes
  REAL(RealKind) :: cL,cU,cC
  REAL(RealKind) :: xL,xU,xC
  REAL(RealKind) :: ccL,ccC,ccR,ccF(nAqua),cDiff,ccFSum
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2,r,Flux
  REAL(RealKind) :: RateLoc
  REAL(RealKind) :: Rate(0:nFrac)
  REAL(RealKind) :: RateDiss
  REAL(RealKind) :: mGes(nFrac)
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: dx(0:nFrac+1)
  REAL(RealKind) :: MassAv(nFrac)
  TYPE (Reaction_T), POINTER :: Reaction
  

  Rate=Zero
  DO i=1,nFrac
    Rate(i)=c(iRelax)%c(1,1,1,i)
  END DO
  Rate(0)=Zero
  DO i=1,nFrac-1
    Rate(i)=(Rate(i)*(x(i+2)-x(i+1)) &
            +Rate(i+1)*(x(i+1)-x(i))) &
              /(x(i+2)-x(i))
  END DO
  DO i=1,nFrac
    DO is=1,nAqua
      MM(i,is)=c(is)%c(1,1,1,i)
    END DO
    mGes(i)=TotalMass(MM(i,:))
    IF (mGes(i)>Zero) THEN
      MM(i,:)=MM(i,:)/mGes(i)
    END IF
  END DO
  DO i=2,nFrac
    dx(i)=LOG(x(i+1)/x(i))
!   dx(i)=x(i+1)-x(i)
  END DO
  dx(1)=dx(2)
  dx(0)=dx(1)
  dx(nFrac+1)=dx(nFrac)
  MM(0,:)=MM(1,:)
  MM(nFrac+1,:)=MM(nFrac,:)

  DO i=1,nFrac
    nC=c(iNC)%c(1,1,1,i)
    IF (nc>Numin.AND.mGes(i)>Zero) THEN
      xL=x(i)
      xU=x(i+1)
      CALL Ansatz(nC,mGes(i),xL,xU,cL,cU)
      ImLGes=Rate(i-1)
      ImUGes=Rate(i)
      IF (ImLGes<Zero) THEN
        fc(iNC)%c(1,1,1,i)=fc(iNC)%c(1,1,1,i)+ImLGes*cL
        fc(iNC)%c(1,1,1,i-1)=fc(iNC)%c(1,1,1,i-1)-ImLGes*cL
        ccFSum=Zero
        DO is=3,nAqua
          ccL=MM(i-1,is)
          ccC=MM(i,is)
          ccR=MM(i+1,is)
          hL=dx(i-1)
          hC=dx(i)
          hR=dx(i+1)
          k1=k1F(hR,hC,hL)
          k2=k2F(hR,hC,hL)
          cDiff=(ccC-ccR)
          r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
          ccFSum=ccFSum+ccF(is)
        END DO
        DO is=3,nAqua
          Flux=ImLGes*cL*xL*ccF(is)/ccFSum
          fc(is)%c(1,1,1,i)=fc(is)%c(1,1,1,i)+Flux
          fc(is)%c(1,1,1,i-1)=fc(is)%c(1,1,1,i-1)-Flux
        END DO
      END IF
      IF (ImUGes>Zero) THEN
        fc(iNC)%c(1,1,1,i)=fc(iNC)%c(1,1,1,i)-ImUGes*cU
        IF (i<nFrac) THEN
          fc(iNC)%c(1,1,1,i+1)=fc(iNC)%c(1,1,1,i+1)+ImUGes*cU
        END IF
        ccFSum=Zero
        DO is=3,nAqua
          ccL=MM(i-1,is)
          ccC=MM(i,is)
          ccR=MM(i+1,is)
          hL=dx(i-1)
          hC=dx(i)
          hR=dx(i+1)
          k1=k1F(hL,hC,hR)
          k2=k2F(hL,hC,hR)
          cDiff=(ccC-ccL)
          r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
          ccFSum=ccFSum+ccF(is)
        END DO
        DO is=3,nAqua
          Flux=ImUGes*cU*xU*ccF(is)/ccFSum
          fc(is)%c(1,1,1,i)=fc(is)%c(1,1,1,i)-Flux
          IF (i<nFrac) THEN
            fc(is)%c(1,1,1,i+1)=fc(is)%c(1,1,1,i+1)+Flux
          END IF
        END DO
      END IF
    END IF
  END DO
  fc(iRelax)%c(:,:,:,:)=Zero
CONTAINS
FUNCTION k1F(hL,hC,hR)
  REAL(RealKind) :: k1F,hL,hC,hR
  k1F=(hC/(hC+hR))*((hL+hC)/(hL+hC+hR))
END FUNCTION k1F
FUNCTION k2F(hL,hC,hR)
  REAL(RealKind) :: k2F,hL,hC,hR
  k2F=(hC/(hC+hL))*(hR/(hL+hC+hR))
END FUNCTION k2F
FUNCTION phi(r,k1,k2)
  REAL(RealKind) :: phi,r,k1,k2
  phi=MAX(Zero,Min(r,MIN(k1*r+k2,One)))
END FUNCTION phi
END SUBROUTINE Advection

END MODULE Advection_Mod
