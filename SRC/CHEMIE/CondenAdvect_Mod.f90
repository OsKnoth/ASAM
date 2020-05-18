MODULE CondenAdvect_Mod

  USE Kind_Mod
  USE DataType_Mod
  USE Ansatz_Mod  &
                ,Ansatz=>AnsatzLogLog
!                 ,Ansatz=>AnsatzHyper
  USE Aerosol_Mod
  USE Chemie_Mod 
  USE Microphysics_Mod

  IMPLICIT NONE
  REAL(RealKind), PRIVATE, PARAMETER :: EpsQuo=1.d-99

CONTAINS

SUBROUTINE AdvectionMass(c,fc,HenryTrans,Time)

  TYPE(Vec4_T) :: c(:),fc(:)
  LOGICAL :: HenryTrans
  REAL(RealKind) :: Time
  REAL(RealKind) :: Kel,Rao

  INTEGER :: i,is,ig,iReak,l
  INTEGER :: ix,iy,iz
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
  REAL(RealKind) :: dm(0:nFrac+1)
  REAL(RealKind) :: MassAv(nFrac)
  TYPE (Reaction_T), POINTER :: Reaction
  

  DO i=2,nFrac
    dm(i)=LOG(m(i+1)/m(i))
  END DO
  dm(1)=dm(2)
  dm(0)=dm(1)
  dm(nFrac+1)=dm(nFrac)

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Rate=Zero
        DO i=1,nFrac
          Rate(i)=c(iRelax)%c(ix,iy,iz,i)
        END DO
        Rate(0)=Zero
        DO i=1,nFrac-1
          Rate(i)=(Rate(i)*(m(i+2)-m(i+1)) &
                  +Rate(i+1)*(m(i+1)-m(i))) &
                    /(m(i+2)-m(i))
        END DO
        DO i=1,nFrac
          DO is=1,nAqua
            MM(i,is)=c(is)%c(ix,iy,iz,i)
          END DO
          mGes(i)=TotalMass(MM(i,:))
          IF (mGes(i)>Zero) THEN
            MM(i,:)=MM(i,:)/mGes(i)
          END IF
        END DO
        MM(0,:)=MM(1,:)
        MM(nFrac+1,:)=MM(nFrac,:)

        DO i=1,nFrac
          nC=c(iNC)%c(ix,iy,iz,i)
          IF (nc>Numin.AND.mGes(i)>Zero) THEN
            xL=m(i)
            xU=m(i+1)
            CALL Ansatz(nC,mGes(i),xL,xU,cL,cU)
            ImLGes=Rate(i-1)
            ImUGes=Rate(i)
            IF (ImLGes<Zero) THEN
              ccL=c(iNC)%c(ix,iy,iz,i-1)
              ccC=c(iNC)%c(ix,iy,iz,i)
              ccR=c(iNC)%c(ix,iy,iz,i+1)
              hL=dm(i-1)
              hC=dm(i)
              hR=dm(i+1)
              k1=k1F(hR,hC,hL)
              k2=k2F(hR,hC,hL)
              cDiff=(ccC-ccR)
              r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
              cl=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
              fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)+ImLGes*cL
              fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)+ImLGes*c(iRELAX)%c(ix,iy,iz,i)/(m(i+1)-m(i))
              IF (i>1) THEN
                fc(iNC)%c(ix,iy,iz,i-1)=fc(iNC)%c(ix,iy,iz,i-1)-ImLGes*cL
                fc(iRELAX)%c(ix,iy,iz,i-1)=fc(iRELAX)%c(ix,iy,iz,i-1)-ImLGes*c(iRELAX)%c(ix,iy,iz,i)/(m(i)-m(i-1))
              END IF
              ccFSum=Zero
              DO is=4,nAqua
                ccL=MM(i-1,is)
                ccC=MM(i,is)
                ccR=MM(i+1,is)
                hL=dm(i-1)
                hC=dm(i)
                hR=dm(i+1)
                k1=k1F(hR,hC,hL)
                k2=k2F(hR,hC,hL)
                cDiff=(ccC-ccR)
                r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
                ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
                ccFSum=ccFSum+ccF(is)
              END DO
              ccF(3)=1.0d0-ccFSum
              ccFSum=1.0d0
              DO is=3,nAqua
                Flux=ImLGes*cL*xL*ccF(is)/ccFSum
                fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+Flux
                IF (i>1) THEN
                  fc(is)%c(ix,iy,iz,i-1)=fc(is)%c(ix,iy,iz,i-1)-Flux
                END IF
              END DO
            END IF
            IF (ImUGes>Zero) THEN
              ccL=c(iNC)%c(ix,iy,iz,i-1)
              ccC=c(iNC)%c(ix,iy,iz,i)
              ccR=c(iNC)%c(ix,iy,iz,i+1)
              hL=dm(i-1)
              hC=dm(i)
              hR=dm(i+1)
              k1=k1F(hL,hC,hR)
              k2=k2F(hL,hC,hR)
              cDiff=(ccC-ccL)
              r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
              cU=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
              fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)-ImUGes*cU
              fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)-ImUGes*c(iRELAX)%c(ix,iy,iz,i)/(m(i+1)-m(i))
              IF (i<nFrac) THEN
                fc(iNC)%c(ix,iy,iz,i+1)=fc(iNC)%c(ix,iy,iz,i+1)+ImUGes*cU
                fc(iRELAX)%c(ix,iy,iz,i+1)=fc(iRELAX)%c(ix,iy,iz,i+1)+ImUGes*c(iRELAX)%c(ix,iy,iz,i+1)/(m(i+2)-m(i+1))
              END IF
              ccFSum=Zero
              DO is=4,nAqua
                ccL=MM(i-1,is)
                ccC=MM(i,is)
                ccR=MM(i+1,is)
                hL=dm(i-1)
                hC=dm(i)
                hR=dm(i+1)
                k1=k1F(hL,hC,hR)
                k2=k2F(hL,hC,hR)
                cDiff=(ccC-ccL)
                r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
                ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
!               IF (i==1) THEN
!                 ccF(is)=(ccC*dm(2)+ccR*dm(1))/(dm(1)+dm(2))
!               END IF
                ccFSum=ccFSum+ccF(is)
              END DO
              ccF(3)=1.0d0-ccFSum
              ccFSum=1.0d0
              DO is=3,nAqua
                Flux=ImUGes*cU*xU*ccF(is)/ccFSum
                fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
                IF (i<nFrac) THEN
                  fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
                END IF
              END DO
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
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
END SUBROUTINE AdvectionMass

FUNCTION LatentHeatPot(Rho,RhoV,RhoL,PotT,T)
  REAL(RealKind) :: LatentHeatPot 
  REAL(RealKind) :: Rho,RhoV,RhoL,PotT,T,Lv

  REAL(RealKind) :: RhoD,Cpml,Cvml,Rm,p
    RhoD=Rho-RhoV-RhoL
    Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL
    Cvml=Cvd*RhoD+Cvv*RhoV+Cpl*RhoL
    Rm=Rd*RhoD+Rv*RhoV
    p=Rm*T
    Lv=LatHeat(T)
    LatentHeatPot=-PotT*( &
                    (-Lv/(Cpml*T+Eps) &
                     -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
                     +Rv/Rm                                   &
                    )                                   &
                    -(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
                     )                                   &
                   )
END FUNCTION LatentHeatPot

SUBROUTINE Condensation(c,fc,HenryTrans,Act)

  TYPE(Vec4_T) :: c(0:),fc(0:)
  LOGICAL :: HenryTrans
  TYPE(Vec4_T) :: Act(0:)

  INTEGER :: n
  INTEGER :: i,is,ig,iReak,l
  INTEGER :: ix,iy,iz
  INTEGER :: ip1,im1
  REAL(RealKind) :: nC
  REAL(RealKind) :: RelFac
  REAL(RealKind) :: Rate(0:nFrac),RateLoc
  REAL(RealKind) :: ActLoc
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: mGes(nFrac)
  REAL(RealKind) :: Gas(nGas)
  REAL(RealKind) :: RhoLuft,RhoV,RhoL,T,PotT
  TYPE (Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: RateWater(nFrac)
  REAL(RealKind) :: RateGas(nFrac)
  REAL(RealKind) :: MaxRateGas
  INTEGER :: iRateGas

  REAL(RealKind) :: ttt,Rad

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Rate=Zero
        RateGas=Zero
        MaxRateGas=0.0d0
        DO i=1,nGas
          Gas(i)=c(i+nAqua)%c(ix,iy,iz,1)
        END DO
        Gas(Position('TE')-nAqua)=TAbs%c(ix,iy,iz,1) 
        RhoLuft=c(RhoPos)%c(ix,iy,iz,1)
        RhoL=RhoLC%c(ix,iy,iz,1)
        RhoV=c(RhoVPos)%c(ix,iy,iz,1)
        PotT=c(ThPos)%c(ix,iy,iz,1)
        T=TAbs%c(ix,iy,iz,1)
        DO i=1,nFrac
          nC=c(iNC)%c(ix,iy,iz,i)
          fc(iRelax)%c(ix,iy,iz,i)=-c(iRelax)%c(ix,iy,iz,i)
          DO is=1,nAqua
            MM(i,is)=c(is)%c(ix,iy,iz,i)
          END DO
          mGes(i)=TotalMass(MM(i,:))
          IF (nc>Numin.AND.mGes(i)>Zero) THEN
            ActLoc=Act(iWater)%c(ix,iy,iz,i)
            mFrac=MM(i,:)/nC
            MM(i,:)=MM(i,:)/mGes(i)
            iCell=i
            iF=iz
            RateLoc=CondRateWater(mFrac,Gas,ActLoc)
            Rad=Radius(mFrac)
            fc(iRelax)%c(ix,iy,iz,i)=RelFac1*(RateLoc-c(iRelax)%c(ix,iy,iz,i))
            fc(RhoVPos)%c(ix,iy,iz,1)=fc(RhoVPos)%c(ix,iy,iz,1) &
                              -c(iRelax)%c(ix,iy,iz,i)*nC
            fc(thPos)%c(ix,iy,iz,1)=fc(thPos)%c(ix,iy,iz,1) &
                               +LatentHeatPot(RhoLuft,RhoV,RhoL,PotT,T)*c(iRelax)%c(ix,iy,iz,i)*nC
            fc(iWater)%c(ix,iy,iz,i)=fc(iWater)%c(ix,iy,iz,i)+c(iRelax)%c(ix,iy,iz,i)*nC
            IF (ASSOCIATED(HenryFirst).AND.HenryTrans &
                           .AND.i>=iVec1.AND.i<=iVec2) THEN 
              Reaction=>HenryFirst
              DO iReak=1,nReakHenry
                iEul=iReak
                is=Reaction%SpeciesRight(1)
                ActLoc=Act(is)%c(ix,iy,iz,i)
                RateLoc=CondRateHenry(mFrac,Gas,Reaction,ActLoc,RhoLuft)
                Rate(i)=Rate(i)+RateLoc*MolMass(is)
                ig=Reaction%SpeciesLeft(1)
                IF (Reaction%NumSpeciesLeftAktiv>0) THEN
                  fc(ig)%c(ix,iy,iz,1)=fc(ig)%c(ix,iy,iz,1)-RateLoc*nC
                END IF 
                fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+RateLoc*nC*MolMass(is)
                Reaction=>Reaction%Next
              END DO
            END IF
          END IF
        END DO
        IF (AdvMass) THEN
!         CALL Adv1
!         CALL Adv2
!         CALL Adv3
        END IF  
      END DO
    END DO
  END DO
CONTAINS

SUBROUTINE Adv2

  REAL(RealKind) :: nC,ImLGes,ImUGes,ImGes
  REAL(RealKind) :: cL(1:nFrac),cU(1:nFrac),cC
  REAL(RealKind) :: xL,xU,xC
  REAL(RealKind) :: ccL,ccC,ccR,cDiff,ccFSumL,ccFSumU
  REAL(RealKind) :: ccFL(nFrac,nAqua),ccFU(nFrac,nAqua)
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2,r,Flux
  REAL(RealKind) :: RateLoc
  INTEGER :: i
! REAL(RealKind) :: Rate(0:nFrac)
! REAL(RealKind) :: RateDiss
! REAL(RealKind) :: mGes(nFrac)
! REAL(RealKind) :: mFrac(nAqua)
! REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: dm(0:nFrac+1)
  REAL(RealKind) :: MassAv(nFrac)
  REAL(RealKind) :: RateM(0:nFrac)
  TYPE (Reaction_T), POINTER :: Reaction
  

  IF (.NOT.Dry) THEN 
    DO i=1,nFrac
      Rate(i)=Rate(i)+c(iRelax)%c(ix,iy,iz,i)
    END DO
  END IF
  Rate(0)=Zero
  DO i=1,nFrac-1
    RateM(i)=(Rate(i)*(m(i+2)-m(i+1)) &
            +Rate(i+1)*(m(i+1)-m(i))) &
              /(m(i+2)-m(i))
  END DO
  RateM(0)=Zero

  DO i=2,nFrac
    dm(i)=LOG(m(i+1)/m(i))
    dm(i)=m(i+1)-m(i)
  END DO
  dm(1)=dm(2)
  dm(0)=dm(1)
  dm(nFrac+1)=dm(nFrac)
  MM(0,:)=MM(1,:)
  MM(nFrac+1,:)=MM(nFrac,:)

  cL=Zero
  cu=Zero
  DO i=1,nFrac
    nC=c(iNC)%c(ix,iy,iz,i)
    IF (nc>Numin.AND.mGes(i)>Zero) THEN
      xL=m(i)
      xU=m(i+1)
      CALL Ansatz(nC,mGes(i),xL,xU,cL(i),cU(i))
    END IF  
  END DO  
  DO i=1,nFrac-1
    Flux=Half*(Rate(i)*cU(i)+Rate(i+1)*cL(i+1))+Half*ABS(RateM(i))*(cU(i)-cL(i+1))
    fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)-Flux
    fc(iNC)%c(ix,iy,iz,i+1)=fc(iNC)%c(ix,iy,iz,i+1)+Flux
  END DO    
! Right Side  
  ccFSumU=Zero
  DO i=1,nFrac-1
    hL=dm(i-1)
    hC=dm(i)
    hR=dm(i+1)
    k1=k1F(hR,hC,hL)
    k2=k2F(hR,hC,hL)
    DO is=3,nAqua
      ccL=MM(i-1,is)
      ccC=MM(i,is)
      ccR=MM(i+1,is)
      cDiff=(ccC-ccR)
      r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
      ccFU(i,is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
      ccFSumU=ccFSumU+ccFU(i,is)
    END DO
    IF (ccFSumU>Zero) THEN
      ccFU(i,:)=ccFU(i,:)/ccFSumU
    END IF  
  END DO
! Left Side  
  ccFSumL=Zero
  DO i=2,nFrac
    hL=dm(i-1)
    hC=dm(i)
    hR=dm(i+1)
    k1=k1F(hL,hC,hR)
    k2=k2F(hL,hC,hR)
    DO is=3,nAqua
      ccL=MM(i-1,is)
      ccC=MM(i,is)
      ccR=MM(i+1,is)
      cDiff=(ccC-ccL)
      r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
      ccFL(i,is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
      ccFSumL=ccFSumL+ccFL(i,is)
    END DO
    IF (ccFSumL>Zero) THEN
      ccFL(i,:)=ccFL(i,:)/ccFSumL
    END IF  
  END DO
  DO i=1,nFrac-1
    DO is=3,nAqua
      Flux=Half*m(i+1)*(Rate(i)*cU(i)*ccFU(i,is)+ &
                        Rate(i+1)*cL(i+1)*ccFL(i+1,is)) &
                        +Half*ABS(RateM(i))*m(i+1)*(cU(i)*ccFU(i,is)-cL(i+1)*ccFL(i+1,is))
      fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
      fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
    END DO
  END DO
  fc(iRelax)%c(:,:,:,:)=Zero
END SUBROUTINE Adv2

!SUBROUTINE Adv1
!  REAL(RealKind) :: mFrac(nAqua)
!  REAL(RealKind) :: NN
!  REAL(RealKind) :: ActLoc
!  REAL(RealKind) :: Gas(nGas)
!  REAL(RealKind) :: ImLGes,ImUGes,ImGes
!  REAL(RealKind) :: cL,cU,cC
!  REAL(RealKind) :: xL,xU,xC
!  REAL(RealKind) :: ccL,ccC,ccR,ccF(nAqua),cDiff,ccFSum
!  REAL(RealKind) :: hL,hC,hR
!  REAL(RealKind) :: k1,k2,r,Flux
!  REAL(RealKind) :: RateDiss
!  REAL(RealKind) :: dmLog(0:nFrac+1)
!  REAL(RealKind) :: dm(0:nFrac+1)
!  REAL(RealKind) :: MassAv(nFrac)
!
!
!  DO i=2,nFrac
!    dmLog(i)=LOG(m(i+1)/m(i))
!    dm(i)=m(i+1)-m(i)
!  END DO
!  dm(1)=dm(2)
!  dm(0)=dm(1)
!  dmLog(1)=dmLog(2)
!  dmLog(0)=dmLog(1)
!  dmLog(nFrac+1)=dmLog(nFrac)
!  DO i=1,nFrac
!!   Rate(i)=Rate(i)+c(iRelax)%c(ix,iy,iz,i) OSSI
!  END DO
!  Rate(0)=Zero
!  DO i=1,nFrac-1
!    Rate(i)=(Rate(i)*(m(i+2)-m(i+1)) &
!            +Rate(i+1)*(m(i+1)-m(i))) &
!              /(m(i+2)-m(i))
!  END DO
!  DO i=1,nFrac
!    im1=MAX(i-1,1) 
!    ip1=MIN(i+1,1) 
!    hL=dmLog(i-1)
!    hC=dmLog(i)
!    hR=dmLog(i+1)
!    nC=c(iNC)%c(ix,iy,iz,i)
!    IF (nc>Numin.AND.mGes(i)>Zero) THEN
!      xL=m(i)
!      xU=m(i+1)
!      CALL Ansatz(nC,mGes(i),xL,xU,cL,cU)
!      ImLGes=Rate(i-1)
!      ImUGes=Rate(i)
!      IF (ImLGes<Zero) THEN
!        fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)+ImLGes*cL
!        fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)+ImLGes*c(iRELAX)%c(ix,iy,iz,i)
!        fc(iNC)%c(ix,iy,iz,i-1)=fc(iNC)%c(ix,iy,iz,i-1)-ImLGes*cL
!        fc(iRELAX)%c(ix,iy,iz,i-1)=fc(iRELAX)%c(ix,iy,iz,i-1)-ImLGes*c(iRELAX)%c(ix,iy,iz,i)
!        ccFSum=Zero
!        DO is=4,nAqua
!          ccL=c(is)%c(ix,iy,iz,im1)/(mGes(im1)+1.d-90)
!          ccC=c(is)%c(ix,iy,iz,i)/(mGes(i)+1.d-90)
!          ccR=c(is)%c(ix,iy,iz,ip1)/(mGes(ip1)+1.d-90)
!          k1=k1F(hR,hC,hL)
!          k2=k2F(hR,hC,hL)
!          cDiff=(ccC-ccR)
!          r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
!          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
!          ccFSum=ccFSum+ccF(is)
!        END DO
!        IF (.NOT.Dry) THEN
!          ccF(3)=1.0d0-ccFSum
!          ccFSum=1.0d0
!          DO is=3,nAqua
!            Flux=ImLGes*cL*xL*ccF(is)/ccFSum
!            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+Flux
!            fc(is)%c(ix,iy,iz,i-1)=fc(is)%c(ix,iy,iz,i-1)-Flux
!          END DO
!        ELSE  
!          DO is=4,nAqua
!            Flux=ImLGes*cL*xL*ccF(is)/ccFSum
!            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+Flux
!            fc(is)%c(ix,iy,iz,i-1)=fc(is)%c(ix,iy,iz,i-1)-Flux
!          END DO
!          hL=dmLog(i-1)
!          hC=dmLog(i)
!          hR=dmLog(i+1)
!          ccL=c(3)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
!          ccC=c(3)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
!          ccR=c(3)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
!          k1=k1F(hR,hC,hL)
!          k2=k2F(hR,hC,hL)
!          cDiff=(ccC-ccR)
!          r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
!          ccF(3)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
!          Flux=ImLGes*cL*ccF(3)
!          fc(3)%c(ix,iy,iz,i)=fc(3)%c(ix,iy,iz,i)+Flux
!          fc(3)%c(ix,iy,iz,i-1)=fc(3)%c(ix,iy,iz,i-1)-Flux
!        END IF
!      END IF
!      IF (ImUGes>Zero) THEN
!        fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)-ImUGes*cU
!        fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)-ImUGes*c(iRELAX)%c(ix,iy,iz,i)
!        IF (i<nFrac) THEN
!          fc(iNC)%c(ix,iy,iz,i+1)=fc(iNC)%c(ix,iy,iz,i+1)+ImUGes*cU
!          fc(iRELAX)%c(ix,iy,iz,i+1)=fc(iRELAX)%c(ix,iy,iz,i+1)+ImUGes*c(iRELAX)%c(ix,iy,iz,i+1)
!        END IF
!        ccFSum=Zero
!        DO is=4,nAqua
!          ccL=c(is)%c(ix,iy,iz,im1)/(mGes(im1)+1.d-90)
!          ccC=c(is)%c(ix,iy,iz,i)/(mGes(i)+1.d-90)
!          ccR=c(is)%c(ix,iy,iz,ip1)/(mGes(ip1)+1.d-90)
!          k1=k1F(hL,hC,hR)
!          k2=k2F(hL,hC,hR)
!          cDiff=(ccC-ccL)
!          r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
!          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
!          ccFSum=ccFSum+ccF(is)
!        END DO
!        IF (.NOT.Dry) THEN
!          ccF(3)=1.0d0-ccFSum
!          ccFSum=1.0d0
!          DO is=3,nAqua
!            Flux=ImUGes*cU*xU*ccF(is)/ccFSum
!            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
!            IF (i<nFrac) THEN
!              fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
!            END IF
!          END DO
!        ELSE  
!          DO is=4,nAqua
!            Flux=ImUGes*cU*xU*ccF(is)/ccFSum
!            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
!            IF (i<nFrac) THEN
!              fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
!            END IF
!          END DO
!          hL=dmLog(i-1)
!          hC=dmLog(i)
!          hR=dmLog(i+1)
!          ccL=c(3)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
!          ccC=c(3)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
!          ccR=c(3)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
!          k1=k1F(hL,hC,hR)
!          k2=k2F(hL,hC,hR)
!          cDiff=(ccC-ccL)
!          r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
!          ccF(3)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
!          Flux=ImUGes*cU*ccF(3)
!          fc(3)%c(ix,iy,iz,i)=fc(3)%c(ix,iy,iz,i)-Flux
!          IF (i<nFrac) THEN
!            fc(3)%c(ix,iy,iz,i+1)=fc(3)%c(ix,iy,iz,i+1)+Flux
!          END IF
!        END IF
!      END IF
!    END IF
!  END DO
!END SUBROUTINE Adv1

SUBROUTINE Adv3
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: NN
  REAL(RealKind) :: ActLoc
  REAL(RealKind) :: Gas(nGas)
  REAL(RealKind) :: ImLGes,ImUGes,ImGes
  REAL(RealKind) :: cL,cU,cC
  REAL(RealKind) :: xL,xU,xC
  REAL(RealKind) :: ccL,ccC,ccR,ccF(nAqua),cDiff,ccFSum
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2,r,Flux
  REAL(RealKind) :: RateDiss
  REAL(RealKind) :: dmLog(0:nFrac+1)
  REAL(RealKind) :: dm(0:nFrac+1)
  REAL(RealKind) :: MassAv(nFrac)


  DO i=2,nFrac
    dmLog(i)=LOG(m(i+1)/m(i))
    dm(i)=m(i+1)-m(i)
  END DO
  dm(1)=dm(2)
  dm(0)=dm(1)
  dmLog(1)=dmLog(2)
  dmLog(0)=dmLog(1)
  dmLog(nFrac+1)=dmLog(nFrac)
  DO i=1,nFrac
!   Rate(i)=Rate(i)+c(iRelax)%c(ix,iy,iz,i) OSSI
  END DO
  Rate(0)=Zero
  DO i=1,nFrac-1
!   Rate(i)=(Rate(i)*(m(i+2)-m(i+1)) &
!           +Rate(i+1)*(m(i+1)-m(i))) &
!             /(m(i+2)-m(i))
  END DO
  DO i=1,nFrac
    im1=MAX(i-1,1) 
    ip1=MIN(i+1,1) 
    hL=dmLog(i-1)
    hC=dmLog(i)
    hR=dmLog(i+1)
    nC=c(iNC)%c(ix,iy,iz,i)
    IF (nc>Numin.AND.mGes(i)>Zero) THEN
      xL=m(i)
      xU=m(i+1)
      CALL Ansatz(nC,mGes(i),xL,xU,cL,cU)
      cL=nC/dm(i)
      cU=nC/dm(i)
      ImLGes=Rate(i-1)
      ImLGes=Rate(i)
      ImUGes=Rate(i)
      IF (ImLGes<Zero.AND.i>1) THEN
        fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)+ImLGes*cL
        fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)+ImLGes*c(iRELAX)%c(ix,iy,iz,i)
        fc(iNC)%c(ix,iy,iz,i-1)=fc(iNC)%c(ix,iy,iz,i-1)-ImLGes*cL
        fc(iRELAX)%c(ix,iy,iz,i-1)=fc(iRELAX)%c(ix,iy,iz,i-1)-ImLGes*c(iRELAX)%c(ix,iy,iz,i)
        ccFSum=Zero
        DO is=4,nAqua
          ccL=c(is)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
          ccC=c(is)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
          ccR=c(is)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
          k1=k1F(hR,hC,hL)
          k2=k2F(hR,hC,hL)
          cDiff=(ccC-ccR)
          r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
        END DO
        IF (.NOT.Dry) THEN
          ccF(3)=1.0d0-ccFSum
          ccFSum=1.0d0
          DO is=3,nAqua
            Flux=ImLGes*cL*ccF(is)
            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+Flux
            fc(is)%c(ix,iy,iz,i-1)=fc(is)%c(ix,iy,iz,i-1)-Flux
          END DO
        ELSE  
          DO is=4,nAqua
            Flux=ImLGes*cL*ccF(is)
            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)+Flux
            fc(is)%c(ix,iy,iz,i-1)=fc(is)%c(ix,iy,iz,i-1)-Flux
          END DO
          hL=dmLog(i-1)
          hC=dmLog(i)
          hR=dmLog(i+1)
          ccL=c(3)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
          ccC=c(3)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
          ccR=c(3)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
          k1=k1F(hR,hC,hL)
          k2=k2F(hR,hC,hL)
          cDiff=(ccC-ccR)
          r=(ccL-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(3)=MAX(ccC+phi(r,k1,k2)*(ccC-ccR),Zero)
          Flux=ImLGes*cL*ccF(3)
          fc(3)%c(ix,iy,iz,i)=fc(3)%c(ix,iy,iz,i)+Flux
          fc(3)%c(ix,iy,iz,i-1)=fc(3)%c(ix,iy,iz,i-1)-Flux
        END IF
      END IF
      IF (ImUGes>Zero) THEN
        fc(iNC)%c(ix,iy,iz,i)=fc(iNC)%c(ix,iy,iz,i)-ImUGes*cU
        fc(iRELAX)%c(ix,iy,iz,i)=fc(iRELAX)%c(ix,iy,iz,i)-ImUGes*c(iRELAX)%c(ix,iy,iz,i)
        IF (i<nFrac) THEN
          fc(iNC)%c(ix,iy,iz,i+1)=fc(iNC)%c(ix,iy,iz,i+1)+ImUGes*cU
          fc(iRELAX)%c(ix,iy,iz,i+1)=fc(iRELAX)%c(ix,iy,iz,i+1)+ImUGes*c(iRELAX)%c(ix,iy,iz,i+1)
        END IF
        ccFSum=Zero
        DO is=4,nAqua
          ccL=c(is)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
          ccC=c(is)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
          ccR=c(is)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
          k1=k1F(hL,hC,hR)
          k2=k2F(hL,hC,hR)
          cDiff=(ccC-ccL)
          r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(is)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
          ccFSum=ccFSum+ccF(is)
        END DO
        IF (.NOT.Dry) THEN
          ccF(3)=1.0d0-ccFSum
          ccFSum=1.0d0
          DO is=3,nAqua
            Flux=ImLGes*cU*ccF(is)
            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
            IF (i<nFrac) THEN
              fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
            END IF
          END DO
        ELSE  
          DO is=4,nAqua
            Flux=ImLGes*cU*ccF(is)
            fc(is)%c(ix,iy,iz,i)=fc(is)%c(ix,iy,iz,i)-Flux
            IF (i<nFrac) THEN
              fc(is)%c(ix,iy,iz,i+1)=fc(is)%c(ix,iy,iz,i+1)+Flux
            END IF
          END DO
          hL=dmLog(i-1)
          hC=dmLog(i)
          hR=dmLog(i+1)
          ccL=c(3)%c(ix,iy,iz,im1)/c(1)%c(ix,iy,iz,im1)
          ccC=c(3)%c(ix,iy,iz,i)/c(1)%c(ix,iy,iz,i)
          ccR=c(3)%c(ix,iy,iz,ip1)/c(1)%c(ix,iy,iz,ip1)
          k1=k1F(hL,hC,hR)
          k2=k2F(hL,hC,hR)
          cDiff=(ccC-ccL)
          r=(ccR-ccC+EpsQuo)/(cDiff+EpsQuo)
          ccF(3)=MAX(ccC+phi(r,k1,k2)*(ccC-ccL),Zero)
          Flux=ImUGes*cU*ccF(3)
          fc(3)%c(ix,iy,iz,i)=fc(3)%c(ix,iy,iz,i)-Flux
          IF (i<nFrac) THEN
            fc(3)%c(ix,iy,iz,i+1)=fc(3)%c(ix,iy,iz,i+1)+Flux
          END IF
        END IF
      END IF
    END IF
  END DO
END SUBROUTINE Adv3

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
END SUBROUTINE Condensation 

SUBROUTINE CondensationJac(c,Jac,HenryTrans,Act)

  TYPE(Vec4_T) :: c(0:)
  TYPE(Vec4_T), POINTER :: Jac(:)
  LOGICAL :: HenryTrans
  TYPE(Vec4_T) :: Act(0:)

  INTEGER :: i,l,k,idf,is,istr,iSpecies,ig,iReak
  REAL(RealKind) :: nC,ImLGes,ImUGes
  REAL(RealKind) :: cL,cU,cC
  REAL(RealKind) :: xL,xU,xC
  REAL(RealKind) :: RelFac
  REAL(RealKind) :: MassAv(nFrac)
  REAL(RealKind) :: Rate(0:nFrac,nAqua)
  REAL(RealKind) :: mGes
  REAL(RealKind) :: RhoLuft
  REAL(RealKind) :: ActLoc
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: RateJac(nAqua+nGas)
  REAL(RealKind) :: RateJacHenry(3)
  REAL(RealKind) :: MM(nAqua)
  REAL(RealKind) :: Gas(nGas)
  INTEGER, POINTER :: Struct(:)
  TYPE (Reaction_T), POINTER :: Reaction

  Struct=>WaterFirst%Struct
  FLAG=.FALSE.
  DO izS=iz0+1,iz1
    DO iyS=iy0+1,iy1
      DO ixS=ix0+1,ix1
        Rate=0.0d0
        DO i=1,nGas
          Gas(i)=c(i+nAqua)%c(ixS,iyS,izS,1)
        END DO
        Gas(Position('TE')-nAqua)=TAbs%c(ixS,iyS,izS,1)  !   OSWALD 
        RhoLuft=c(RhoPos)%c(ixS,iyS,izS,1)
        DO i=1,nFrac
          nC=c(iNC)%c(ixS,iyS,izS,i)
          DO is=1,nAqua
            MM(is)=c(is)%c(ixS,iyS,izS,i)
          END DO
          mGes=TotalMass(MM)
          IF (nc>NuMin.AND.mGes>Zero) THEN
            mFrac=MM(:)/nC
            mFrac(iNC)=Zero
            ActLoc=Act(iWater)%c(ixS,iyS,izS,i)
            Rate(i,iWater)=CondRateWater(mFrac,Gas,ActLoc)
            RateJac=CondRateWaterJac(mFrac,Gas,ActLoc)
            RateJac(2)=-nC
            Reaction=>WaterFirst 
            istr=0
            Struct=>Reaction%Struct
            DO k=1,Reaction%NumSpeciesLeftAktiv
              iSpecies=Reaction%SpeciesLeft(k)
              IF (Reaction%SpeciesLeft(k)==Position('aNUMBER')) THEN
                DO l=1,Reaction%NumSpecies
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                  ,-Reaction%Species(l)%Koeff&
                            *SUM(RateJac(1:nAqua)*mFrac)/nC,i)
                END DO
              ELSE IF(iSpecies<=nAqua) THEN
                DO l=1,Reaction%NumSpecies
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,RateJac(iSpecies)*Reaction%Species(l)%Koeff/nC,i) 
                END DO
              ELSE
                DO l=1,Reaction%NumSpecies
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,RateJac(iSpecies)*Reaction%Species(l)%Koeff,i) 
                END DO
              END IF
            END DO
            Reaction=>Reaction%Next
            istr=0
            Struct=>Reaction%Struct
            RateJac(iNC)=c(iRelax)%c(ixS,iyS,izS,i)
            RateJac(iRELAX)=nC
            DO k=1,Reaction%NumSpeciesLeftAktiv
              iSpecies=Reaction%SpeciesLeft(k)
              DO l=1,Reaction%NumSpecies
                istr=istr+1
                idf=Struct(istr)
                Call AddScalar(Jac(idf),Jac(idf) &
                              ,Reaction%Species(l)%Koeff&
                               *RateJac(iSpecies),i)
              END DO
            END DO
            IF (ASSOCIATED(HenryFirst).AND.HenryTrans &
                           .AND.i>=iVec1.AND.i<=iVec2) THEN 
              Reaction=>HenryFirst
              DO iReak=1,nReakHenry
                is=Reaction%SpeciesRight(1)
                ig=Reaction%SpeciesLeft(1)-nAqua
                ActLoc=Act(is)%c(ixS,iyS,izS,i)
                Rate(i,is)=CondRateHenry(mFrac,Gas,Reaction,ActLoc,RhoLuft)
                RateJacHenry=CondRateHenryJac(mFrac,Gas,Reaction,ActLoc,RhoLuft)
                Struct=>Reaction%Struct
                istr=0
!               fcGas(ig)%c(1)=fcGas(ig)%c(1)-Rate(i,is)*nC
!               fc(is+1)%c(i)=fc(is+1)%c(i)+Rate(i,is)*nC*MolMass(is+1)
                IF (Reaction%NumSpeciesLeftAktiv>0) THEN
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,nC*RateJacHenry(1)*Reaction%Species(1)%Koeff,i) 
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,MolMass(is)*nC*RateJacHenry(1)*Reaction%Species(2)%Koeff,i) 
                END IF
                IF (Reaction%NumSpeciesLeftAktiv>0) THEN
!                 right Species
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,RateJacHenry(2)*Reaction%Species(1)%Koeff,i) 
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,MolMass(is)*RateJacHenry(2)*Reaction%Species(2)%Koeff,i) 
!                 aNUMBER
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                  ,Reaction%Species(1)%Koeff&
                                  *(Rate(i,is)-(RateJacHenry(2)*mFrac(is) &
                                         +RateJacHenry(3)*mFrac(iWater))),i)
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                  ,MolMass(is)*Reaction%Species(2)%Koeff&
                                  *(Rate(i,is)-(RateJacHenry(2)*mFrac(is) &
                                         +RateJacHenry(3)*mFrac(iWater))),i)
!                 aH2O
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,RateJacHenry(3)*Reaction%Species(1)%Koeff,i) 
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,MolMass(is)*RateJacHenry(3)*Reaction%Species(2)%Koeff,i) 
                ELSE
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                  ,MolMass(is)*Reaction%Species(1)%Koeff&
                                  *(Rate(i,is)-(RateJacHenry(2)*mFrac(is) &
                                               )),i)
                  istr=istr+1
                  idf=Struct(istr)
                  Call AddScalar(Jac(idf),Jac(idf) &
                                    ,RateJacHenry(2)*Reaction%Species(1)%Koeff*MolMass(is),i) 
                END IF
                Rate(i,is)=Rate(i,is)*MolMass(is)
                Reaction=>Reaction%Next
              END DO
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CondensationJac 

END MODULE CondenAdvect_Mod
