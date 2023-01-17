MODULE Wall_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE Parameter_Mod
  USE ParameterMicrophys_Mod
  USE BoundaryCond_Mod
  USE Control_Mod
  USE Chemie_Mod
  USE Soil_Mod
  USE SoilData_Mod
  USE RevWRF_Mod
  USE Output_Mod
  USE Canopy_Mod
  USE Emission_Mod

IMPLICIT NONE

  INTEGER, PRIVATE :: i,j,k
  INTEGER, PRIVATE :: uPosJac,vPosJac,wPosJac
  REAL(RealKind), PRIVATE, POINTER :: c(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: th(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: T(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThB(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Wg(:,:)
  REAL(RealKind), PRIVATE, POINTER :: tke(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: dis(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoS(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NRain(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NIce(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NSnow(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uFRhs(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vFRhs(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wFRhs(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: D(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoVRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoLRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoRRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoIRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoSRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NRRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NIRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: NSRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoRcBRhs(:,:)
  REAL(RealKind), PRIVATE, POINTER :: RR(:,:)
  REAL(RealKind), PRIVATE, POINTER :: IR(:,:)
  REAL(RealKind), PRIVATE, POINTER :: SR(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Raddir(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Raddif(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Radinf(:,:)
  
  TYPE(Vec4_T), PRIVATE, POINTER :: AS(:)


CONTAINS

SUBROUTINE DragCoeff(VelF,VectorCell,UVec)
 
  TYPE(VelocityFace_T) :: VelF(:)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(Vector4Cell_T), POINTER, OPTIONAL :: UVec(:)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VN,VT

  IF (PRESENT(UVec)) THEN
    uCL=>UVec(ibLoc)%Vec(uPosL)%c
    vCL=>UVec(ibLoc)%Vec(vPosL)%c
    wCL=>UVec(ibLoc)%Vec(wPosL)%c
    uCR=>UVec(ibLoc)%Vec(uPosR)%c
    vCR=>UVec(ibLoc)%Vec(vPosR)%c
    wCR=>UVec(ibLoc)%Vec(wPosR)%c
  ELSE
    uCL=>VectorCell(ibLoc)%Vec(uPosL)%c
    vCL=>VectorCell(ibLoc)%Vec(vPosL)%c
    wCL=>VectorCell(ibLoc)%Vec(wPosL)%c
    uCR=>VectorCell(ibLoc)%Vec(uPosR)%c
    vCR=>VectorCell(ibLoc)%Vec(vPosR)%c
    wCR=>VectorCell(ibLoc)%Vec(wPosR)%c
  END IF
  uF=>VelF(ibLoc)%uF
  vF=>VelF(ibLoc)%vF
  wF=>VelF(ibLoc)%wF
  Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
  Th=>VectorCell(ibLoc)%Vec(ThPos)%c
  T=>TAbsCell(ibLoc)%Vec(1)%c
  ThB=>VectorCell(ibLoc)%Vec(ThPos)%cB
  tke=>VectorCell(ibLoc)%Vec(tkePos)%c
  RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
  IF (RhoCPos>0) THEN
    RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
  ELSE
    RhoL=>RhoLCell(ibLoc)%c
  END IF
  RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
  RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
  RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
! IF (Radiation) THEN
!   Raddir=>RaddirCell(ibLoc)%cB
!   Raddif=>RaddifCell(ibLoc)%cB
!   Radinf=>RadinfCell(ibLoc)%cB
! ELSE
!   ALLOCATE(Raddir(NumBoundCell,Domain%nrSoilLayers))
!   ALLOCATE(Raddif(NumBoundCell,Domain%nrSoilLayers))
!   ALLOCATE(Radinf(NumBoundCell,Domain%nrSoilLayers))
! END IF
  
  DO i=1,NumBoundCell
    ix=BoundCell(i)%ix
    iy=BoundCell(i)%iy
    iz=BoundCell(i)%iz
    n1=BoundCell(i)%n1
    n2=BoundCell(i)%n2
    n3=BoundCell(i)%n3
!   Evaluate velocity tangential and normal  
    uLoc=(FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
          FU(ix,iy,iz)*uF(ix,iy,iz)) / &
         ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(Rho(ix,iy,iz,1)+Eps))
    vLoc=(FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
          FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(Rho(ix,iy,iz,1)+Eps))
    wLoc=(FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
          FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(Rho(ix,iy,iz,1)+Eps))
    VN=uLoc*n1+vLoc*n2+wLoc*n3
    V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
    VT=SQRT(MAX(V-VN*VN,Zero))           
    BoundCell(i)%VT1=VT
    IF (FireEmiss) THEN
      uLoc=(FU(ix-1,iy,iz+1)*uF(ix-1,iy,iz+1)+  &
            FU(ix,iy,iz+1)*uF(ix,iy,iz+1)) / &
           ((FU(ix-1,iy,iz+1)+FU(ix,iy,iz+1)+Eps)*(Rho(ix,iy,iz+1,1)+Eps))
      vLoc=(FV(ix,iy-1,iz+1)*vF(ix,iy-1,iz+1)+  &
            FV(ix,iy,iz+1)*vF(ix,iy,iz+1)) / &
           ((FV(ix,iy-1,iz+1)+FV(ix,iy,iz+1)+Eps)*(Rho(ix,iy,iz+1,1)+Eps))
      wLoc=(FW(ix,iy,iz+1-1)*wF(ix,iy,iz+1-1)+  &
            FW(ix,iy,iz+1)*wF(ix,iy,iz+1)) /  &
           ((FW(ix,iy,iz+1-1)+FW(ix,iy,iz+1)+Eps)*(Rho(ix,iy,iz+1,1)+Eps))
      VN=uLoc*n1+vLoc*n2+wLoc*n3
      V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
      VT=SQRT(MAX(V-VN*VN,Zero))
      BoundCell(i)%VT2=VT
    END IF  

    IF (DynamicSoil) THEN
!     ThB=>VectorCell(ibLoc)%Vec(ThPos)%cB
!     BoundCell(i)%TeS=ThB(i,1)
      IF(BoundCell(i)%LandClass==9) THEN
        CALL DragCoeffSea(i)
      ELSE
!       Wg=>VectorCell(ibLoc)%Vec(RhoVPos)%cB
        CALL DragCoeffSoil(i)
      END IF
    ELSE IF (Canopy) THEN
      ThB=>VectorCell(ibLoc)%Vec(ThPos)%cB
      BoundCell(i)%TeS=ThB(i,1)
      Wg=>VectorCell(ibLoc)%Vec(RhoVPos)%cB
      CALL DragCoeffSoil(i)
    ELSE
      IF (BoundCell(i)%LandClass==9) THEN
        CALL DragCoeffSea(i)
      ELSE
        IF (Wall) THEN
          IF (Radiation) THEN
            CALL DragCoeffWallRad(i)
          ELSE
            CALL DragCoeffWall(i)
          END IF
        ELSE IF (Louis) THEN
          IF (Radiation) THEN
            CALL DragCoeffLouisRad(i)
          ELSE
            CALL DragCoeffLouis(i)
          END IF
        ELSE  
          CALL DragCoeffSoil(i)
        END IF
      END IF
    END IF
  END DO

END SUBROUTINE DragCoeff

SUBROUTINE DragCoeffWallRad(i)
  
  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,dL,zRauh,zRauhT,zPL,logz
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn1,Tgn2
  REAL(RealKind) :: RhoLoc,ThLoc,ThPotSurf,Fpot
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: RiB0,RiB1,RiBT
  REAL(RealKind) :: DragM,DragH
  REAL(RealKind) :: Qg,Qc,Qdir,Qem,Qref,Qdif,RhoInf
  REAL(RealKind) :: P,s,F,Fs,Prt

  REAL(RealKind) ,Parameter :: logmin=0.1d0,logmax=1.d10,zero=0.d0

! Turbulent Prandtl number  
  Prt=0.95
  
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    dL     = BoundCell(i)%dL
    zRauh  = BoundCell(i)%zRauh
    zRauhT = BoundCell(i)%zRauhT
    alb    = BoundCell(i)%alb
    ee     = BoundCell(i)%ee
    shad   = BoundCell(i)%shad
    RhoLoc = Rho(ix,iy,iz,1)
    Tgn1   = BoundCell(i)%TeS
    ThLoc  = Th(ix,iy,iz,1)/RhoLoc
    DragM  = BoundCell(i)%DragM
    DragH  = BoundCell(i)%DragH
    
!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc+Eps))

! OLD!!!
!   Evaluate velocity tangential and normal  
!    uLoc = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoLoc+Eps)
!    vLoc = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoLoc+Eps)
!    wLoc = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoLoc+Eps)
    
    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))           

!   Evaluate height z
    zPL = Half*VolC(ix,iy,iz)/FL
    zPL = dL                      !Oswald
  
    IF (zRauh.GT.zero) THEN
      logz = MAX(LOG(zPL/zRauh+One),logmin)
    ELSE
      logz = logmax
    END IF
    
!   First estimate for potential temperature
!    Fpot=(p0/(RhoLoc*Rd))**kappa
!    ThPotSurf=Tgn1**(Cv/Cp)*Fpot
    Fpot = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
    ThPotSurf = Tgn1*Fpot  
        
    DO it=0,10      
      P = 9.24d0*((PrandtlNumber(thPos)/Prt)**0.75d0-1.0d0)*(1.0d0+0.28d0*EXP(-0.007d0*PrandtlNumber(thPos)/Prt))
      s = logz/Karm
!     Drag coefficient for heat and momentum
      DragH = 1.0d0/((s**2.0d0)*Prt*(1.0d0+P/s))
      DragM = (One/s)**2.0d0

!      ThPotSurf=Tgn1**(Cv/Cp)*Fpot
      ThPotSurf = Tgn1*Fpot

!     Bodenwaermestrom
      Qg   = -BoundCell(i)%D*(Tgn1-BoundCell(i)%T1)/(0.5d0)
!     Waermestrom in die Luft
      Qc   = -RhoLoc*Cpd*DragH*VT*(ThPotSurf-ThLoc)
!     direkte Strahlung
      Qdir = (1-alb)*shad*MAX(Zero,(n1*radn1+n2*radn2+n3*radn3))*BoundCell(i)%raddirekt
!     diffuse Strahlung
      Qdif = (1-alb)*BoundCell(i)%skyviewfactor*BoundCell(i)%raddiffus
      ! downward long wave radiation
      RhoInf = BoundCell(i)%raddiffus
!     Ausstrahlung im Infraroten
      Qem  = -ee*SBsigma*Tgn1**4.0d0

      ! Mehrfachreflexion zwischen Gebaeuden      
      !do j=1,NumBoundCell
      ! Qref=viewfactor(i,j)* &
      !      (BoundCell(j)%SoilCell%alb*BoundCell(j)%shad*BoundCell(i)%raddirekt+ &
      !       BoundCell(j)%SoilCell%alb*BoundCell(i)%raddiffus+ &
      !       BoundCell(j)%SoilCell%ee*5.6704e-8*BoundCell(j)%SoilCell%Tsold**4.)
      !enddo
      Qref = 0.0d0

      F    = Qg+Qc+Qdir+Qdif+RhoInf+Qem+Qref
      Fs   = -BoundCell(i)%D/(0.5d0)-RhoLoc*Cpd*DragH*VT*Fpot-ee*4.0d0*5.6704e-8*Tgn1**3.0d0
      Tgn2 = Tgn1-F/Fs

      IF (ABS(Tgn2-Tgn1)<1.0d-7) THEN
         Tgn1 = Tgn2
         EXIT
      ELSE
        Tgn1 = Tgn2 
      END IF
    END DO
   
   BoundCell(i)%TeS    = Tgn1
   BoundCell(i)%ThetaS = ThPotSurf
   BoundCell(i)%DragM  = DragM
   BoundCell(i)%DragH  = DragH
   BoundCell(i)%DragQ  = DragH
 
END SUBROUTINE DragCoeffWallRad


SUBROUTINE DragCoeffSoil(i)
  
  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it,shad
  INTEGER :: LandClass
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: n1G,n2G,n3G
  REAL(RealKind) :: dL,zRauh,zRauhT,zPL,logz
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn,AbsT
  REAL(RealKind) :: RhoLoc,ThLoc
  REAL(RealKind) :: RiB
  REAL(RealKind) :: DragM,DragH
  REAL(RealKind) :: FCh,FCm
  REAL(RealKind) :: TermM,TermH
  REAL(RealKind) :: Chn,Cmn,fh,fm
  REAL(RealKind) :: Ch,Cm
  REAL(RealKind) :: zL
  REAL(RealKind) :: RichCut
  REAL(RealKind) :: Chi,ChiS
  REAL(RealKind) :: fChi,fChiS,dfdChi
  REAL(RealKind) ,Parameter :: RichCrit=0.25d0
  INTEGER :: Iter


    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    n1G    = BoundCell(i)%n1G
    n2G    = BoundCell(i)%n2G
    n3G    = BoundCell(i)%n3G
    zRauh  = BoundCell(i)%zRauh
    zRauhT = BoundCell(i)%zRauhT
    dL     = BoundCell(i)%dL
    LandClass = BoundCell(i)%LandClass
    RhoLoc = Rho(ix,iy,iz,1)
    IF (DynamicSoil) THEN
      Tgn    = ThB(i,1) 
    ELSE
      Tgn    = BoundCell(i)%TeS
    END IF
    AbsT   = T(ix,iy,iz,1)
    ThLoc  = Th(ix,iy,iz,1)/RhoLoc
    DragM  = BoundCell(i)%DragM
    DragH  = BoundCell(i)%DragH
    
!   Evaluate velocity tangential and normal  
    uLoc = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoLoc+Eps)
    vLoc = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoLoc+Eps)
    wLoc = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoLoc+Eps)

    VN   = uLoc*n1G+vLoc*n2G+wLoc*n3G
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))           
    
    IF (VT>1d-3) THEN

!   Evaluate height z
    zPL = dL                      !Oswald

!   Bulk-Richardson-Number
    RiB  = Grav/Tgn*((AbsT-Tgn)*(dz(iz)-zRauh))/(VT*VT+Eps)

! ****************************************************
    IF (DragScheme=="revLouis") THEN  ! revised Louis scheme (ECMWF)
! ****************************************************

!     Transfer Coefficients for neutral stratification 
      Cmn  = Karm**2.0d0*(LOG(dz(iz)/zRauh))**(-2.0d0)
      Chn  = 0.74d0*Karm**2.0d0*(LOG(dz(iz)/zRauh))**(-1.0d0)*(LOG(dz(iz)/zRauhT))**(-1.0d0)
  
      TermM = 3.0d0*bb*cc*((dz(iz)/zRauh)**(One/Three)-1.0d0)**(Three/Two)*SQRT(ABS(RiB))
      TermH = 3.0d0*bb*cc*((dz(iz)/zRauhT)**(One/Three)-1.0d0)**(Three/Two)*SQRT(ABS(RiB))
!     Stability-Function for
      IF (RiB<Zero) THEN
!       unstable stratification
        fm = 1.0d0+(2.0d0*bb*ABS(RiB))/  &      
                   (1.0d0+Cmn*TermM)   
        fh = 1.0d0+(3.0d0*bb*ABS(RiB))/  &
                   (1.0d0+Chn*TermH)  
      ELSE
!       stable stratification
        fm = 1.0d0/(1.0d0+2.0d0*bb*RiB/(SQRT(1.0d0+dd*RiB))) 
        fh = 1.0d0/(1.0d0+3.0d0*bb*RiB*SQRT(1.0d0+dd*RiB))  
      END IF

!     Transfer Coefficients in case of |v|->0 and Rib->-Inf 
      FCh  = SQRT(Grav/Tgn*ABS(AbsT-Tgn)*(dz(iz)-zRauh))/ &
             (cc*((dz(iz)/zRauhT)**(One/Three)-1.0d0)**(Three/Two))
      FCm  = Two/Three*SQRT(Grav/Tgn*ABS(AbsT-Tgn)*(dz(iz)-zRauh))/ &
             (cc*((dz(iz)/zRauh)**(One/Three)-1.0d0)**(Three/Two))
    
!     Bulk-Aerodynamical Transfer Coefficient for Momentum/Heat
      IF (VT<1.d-3) THEN
        Cm = FCm/(VT+Eps)
        Ch = FCh/(VT+Eps)
      ELSE
        Cm = Cmn*fm
        Ch = Chn*fh
      END IF

! ***************************************************
    ELSE IF (DragScheme=="revWRF") THEN   ! revised WRF scheme
! ***************************************************

      RiB=MAX(RiB,-10.0d0)
      Chi=RiB
      Iter=0
      DO
        fChi=fx(Chi)
        ChiS=Chi*(1.0d0+1.e-8_RealKind)+SIGN(1.e-8_RealKind,Chi)
        fChiS=fx(ChiS)
        dfdChi=(fChiS-fChi)/(ChiS-Chi)
        CALL NewtonMethod(Chi,zL,fChi,dfdChi) 
        Iter=Iter+1
        IF (ABS(Chi-zL)<1.d-6)THEN
          EXIT
        ELSE IF (Iter>20) THEN
          WRITE(*,*) 'ix,iy,iz',ix,iy,iz
          WRITE(*,*) 'RiB',RiB
          WRITE(*,*) 'zRauh',zRauh
          WRITE(*,*) 'zRauhT',zRauhT
          WRITE(*,*) 'dz(iz)',dz(iz)
          EXIT
        ELSE
          Chi=zL
        END IF
      END DO
      TermH = LOG((dz(iz)+zRauh)/zRauhT)      &
             -PsiH(zL*(One+zRauh/dz(iz)),dz(iz),zRauh,zRauhT,RiB,LandClass) &
             +PsiH(zL*zRauhT/dz(iz),dz(iz),zRauh,zRauhT,RiB,LandCLass)
      TermM = LOG((dz(iz)+zRauh)/zRauh)       &
             -PsiM(zL*(One+zRauh/dz(iz)),dz(iz),zRauh,zRauhT,RiB,LandCLass) &
             +PsiM(zL*zRauh/dz(iz),dz(iz),zRauh,zRauhT,RiB,LandCLass)
      Cm    = Karm**Two/TermM/TermM
      Ch    = Karm**Two/TermM/TermH
    END IF
  
    BoundCell(i)%DragM  = Cm
    BoundCell(i)%DragH  = Ch
    BoundCell(i)%DragQ  = Ch
    
    ELSE

    BoundCell(i)%DragM=Zero
    BoundCell(i)%DragH=Zero
    BoundCell(i)%DragQ=Zero

    END IF

    CONTAINS
   
      FUNCTION fx(x)
        REAL(RealKind) :: fx
        REAL(RealKind) :: x
        fx=RiB-RiBulk(dz(iz),zRauh,zRauhT,x,RiB)
      END FUNCTION fx
  
!     FUNCTION dfdx(x)
!       REAL(RealKind) :: dfdx
!       REAL(RealKind) :: x
!       dfdx=-dRiBdChi(dz(iz),zRauh,zRauhT,x,RiB)
!     END FUNCTION dfdx
 
END SUBROUTINE DragCoeffSoil

SUBROUTINE DragCoeffWall(i)
  
  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: FL,dL,zRauh,zPL,logz
  REAL(RealKind) :: DragM
  REAL(RealKind) :: Prt

  REAL(RealKind) ,Parameter :: logmin=0.1d0,logmax=1.d10,zero=0.d0

! Turbulent Prandtl number  
  Prt = 0.95
  
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    FL     = BoundCell(i)%FL
    dL     = BoundCell(i)%dL
    zRauh  = BoundCell(i)%zRauh
    
!   Evaluate height z
    zPL = dL                      !Oswald
    IF (zRauh.GT.zero) THEN
      logz = MAX(LOG(zPL/zRauh+One),logmin)
    ELSE
      logz = logmax
    END IF
!   Set Drag Coefficient
    DragM = (Karm/logz)**Two
    BoundCell(i)%DragM  = DragM
    
END SUBROUTINE DragCoeffWall 

SUBROUTINE DragCoeffLouisRad(i)
  
  INTEGER :: i
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,dL,zRauh,zRauhT,zPL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn1,Tgn2
  REAL(RealKind) :: RhoLoc,ThLoc,ThPotSurf,Fpot
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: Qg,Qc,Qdir,Qem,Qref,Qdif,RhoInf
  REAL(RealKind) :: RiB0,RiB1,RiBT,R,FM,FH,cM,cH
  REAL(RealKind) :: N,dFMdRiB0,dFHdRiB0,dRiBTdTgn,dNdTgn,dNdRiB0,dQgdTgn
  REAL(RealKind) :: F,dDragHdRiB0,dThPotSurfdTgn,dQcdTgn,dQcdRiB0,dQemdTgn,dFdTgn,dFdRiB0
  REAL(RealKind) :: ustar,thstar,DragM,DragH
  REAL(RealKind) :: det
  
  REAL(RealKind) ,Parameter :: zero=0.d0
  
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    dL     = BoundCell(i)%dL
    zRauh  = BoundCell(i)%zRauh
    zRauhT = BoundCell(i)%zRauhT
    alb    = BoundCell(i)%alb
    ee     = BoundCell(i)%ee
    shad   = BoundCell(i)%shad
    RhoLoc = Rho(ix,iy,iz,1)
    Tgn1   = BoundCell(i)%TeS
    ThLoc  = Th(ix,iy,iz,1)/RhoLoc
    DragM  = BoundCell(i)%DragM
    DragH  = BoundCell(i)%DragH

    
!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc+Eps))

    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))         

!   Evaluate height z
    zPL  = Half*VolC(ix,iy,iz)/FL
    zPL  = dL                      !Oswald
    
!   First Estimate for Bulk Richardson Number RiB0 = RiBT
!    Fpot=(p0/(RhoLoc*Rd))**kappa
!    ThPotSurf=Tgn1**(Cv/Cp)*Fpot
    Fpot      = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
    ThPotSurf = Tgn1*Fpot
  
    RiBT      = Grav*zPL*(ThLoc-thPotSurf)/(ThLoc*(VT**2.0d0)+Eps)
    RiB0      = RiBT
    
    DO it=0,5      

!     Louis scheme
!---------------------------------------------------------------
!     Parametric Functions
      IF (RiB0.GT.0.25) THEN  ! critical Richardson number: 0.25
!     Stable Conditions
        FM = One/((One+4.7d0*RiB0)**2.0d0)
        FH = FM
        dFMdRiB0 = -Two*4.7d0/(One+4.7d0*RiB0)**Three
        dFHdRiB0  = dFMdRiB0  
      ELSE 
!     Unstable Conditions
        cM = 7.4d0*(Karm**2.0d0)*9.4d0*SQRT(zPL/(zRauh+Eps))/((LOG(zPL/(zRauh+Eps)))**2.0d0+Eps)
        cH = 5.3d0*cM/7.4d0
        FM = One-9.4d0*RiB0/(One+cM*SQRT(ABS(RiB0)))
        FH = One-9.4d0*RiB0/(One+cH*SQRT(ABS(RiB0)))
        dFMdRiB0 = -9.4d0*(One+cM*SQRT(ABS(RiB0))-0.5d0*RiB0*cM/(SQRT(ABS(RiB0))+Eps))/((One+cM*SQRT(ABS(RiB0)))**Two)
        dFHdRiB0 = -9.4d0*(One+cH*SQRT(ABS(RiB0))-0.5d0*RiB0*cH/(SQRT(ABS(RiB0))+Eps))/((One+cH*SQRT(ABS(RiB0)))**Two)  
      END IF

!     Turbulent Prandtl Number (after Businger)
      R = 7.4d-1
  
!     ThPotSurf=Tgn1**(Cv/Cp)*Fpot  
!      dThPotSurfdTgn=Fpot*(Cv/Cp)*Tgn1**(Cv/Cp-1.0d0)       
      dThPotSurfdTgn = Fpot     
 
!     Bulk Richardson Number  
      RiBT = Grav*zPL*(ThLoc-thPotSurf)/(ThLoc*(VT**Two)+Eps)
      dRiBTdTgn = -Grav*zPL*dThPotSurfdTgn/(ThLoc*(VT**Two)+Eps)

!     Function for Bulk Richardson Number
      N = RiB0*R*LOG(zRauh/(zRauhT+Eps))+(RiB0-RiBT)*R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH
      dNdTgn  = -dRiBTdTgn*R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH
      dNdRiB0 = R*LOG(zRauh/(zRauhT+Eps))+R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH+ &
                (RiB0-RiBT)*R*LOG(zPL/(zRauh+Eps))*(0.5d0*FM**(-0.5d0)*dFMdRiB0*FH-SQRT(FM)*dFHdRiB0)/FH**Two   
  
!     Drag coefficient for heat as initial condition for radiation scheme
      DragH = (Karm/(LOG(zPL/zRauh)+Eps))**Two/(R*(LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*SQRT(FM))+One/FH))
      dDragHdRiB0 = -(Karm/(LOG(zPL/zRauh)+Eps))**Two/(R*(LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*SQRT(FM))+One/FH)**Two) &
                    *(-0.5d0*dFMdRiB0*LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*FM**(Three/Two))-dFHdRiB0/(FH**Two))

!     Radiation scheme to calculate Tgn1 and associated ThPotSurf
!----------------------------------------------------------------
!     Bodenwaermestrom   
      Qg = -BoundCell(i)%D*(Tgn1-BoundCell(i)%T1)/0.5d0 
      dQgdTgn  = -BoundCell(i)%D/0.5d0
!     Waermestrom in die Luft
      Qc = -RhoLoc*Cpd*DragH*VT*(ThPotSurf-ThLoc)
      dQcdTgn  = -RhoLoc*Cpd*DragH*VT*dThPotSurfdTgn
      dQcdRiB0 = -RhoLoc*Cpd*dDragHdRiB0*VT*(ThPotSurf-ThLoc)  

    
!     direkte Strahlung
      Qdir = (One-alb)*shad*MAX(Zero,(n1*radn1+n2*radn2+n3*radn3))*BoundCell(i)%raddirekt  
!     Qdir = 0.0d0
  
!     diffuse Strahlung
      Qdif = (ONE-alb)*BoundCell(i)%raddiffus
!     Qdif = 0.0d0

      ! downward long wave radiation
      RhoInf = BoundCell(i)%radinfred

!     Ausstrahlung im Infraroten   
      Qem  = -ee*SBsigma*Tgn1**4.0d0
      dQemdTgn = -ee*4.0d0*SBsigma*Tgn1**3.0d0
      Qref = 0.0d0     

!     Function for temperature
      F = Qg+Qc+Qdir+Qdif+RhoInf+Qem+Qref           ! Heat balance equation
!     write(*,*) 'Louis, Qg, Qc, Qdir, Qdif, RhoInf, Qem: ', Qg, Qc, Qdir, Qdif, RhoInf, Qem
!     write(*,*) 'Louis, Tgn1: ', Tgn1, it
      dFdTgn  = dQcdTgn+dQemdTgn+dQgdTgn
      dFdRiB0 = dQcdRiB0
    
!     Determinant
      det = dNdRiB0*dFdTgn-dNdTgn*dFdRiB0

!     Numerical Iteration Bulk Richardson Number and temperature at surface 
      RiB1 = RiB0-(dFdTgn*N-dNdTgn*F)/(det+Eps)
      
!     IF(RiB0.LT.0.25d0.AND.RiB1.GE.0.25d0) THEN       !!! CRITICAL RICHARDSON NUMBER PROBLEM
!       RiB0 = 0.25d0
!       GOTO 1
!     ELSE IF(RiB0.GT.0.25d0.AND.RiB1.LE.0.25d0) THEN
!       RiB0 = 0.2500001d0
!       GOTO 1
!     END IF
      
      Tgn2 = Tgn1-(dNdRiB0*F-dFdRiB0*N)/(det+Eps)

!     Friction velocity for momentum and heat        
      ustar  = SQRT(FM*(Karm*VT/(LOG(zPL/zRauh)+Eps))**Two)
!     thstar = Karm**Two*VT*(ThLoc-ThPotSurf)*FH/(LOG(zPL/zRauh)*LOG(zPL/zRauhT)*R*ustar+Eps)

      IF(ABS(RiB1-RiB0)<=10e-6*ABS(RiB1).AND.ABS(Tgn2-Tgn1)<=1.0d-6*Tgn2) THEN  
        RiB0 = RiB1
        Tgn1 = Tgn2
        ThPotSurf = Tgn1*Fpot
        EXIT
      ELSE
        RiB0 = RiB1
        Tgn1 = Tgn2
        ThPotSurf = Tgn1*Fpot
      END IF 
1   CONTINUE
    ENDDO

!   Drag coefficient for momentum and heat    
    DragM=(ustar/(VT+Eps))**Two
    DragH=(Karm/(LOG(zPL/zRauh)+Eps))**Two*FH/(R*(LOG(zRauh/zRauhT)*FH/(LOG(zPL/zRauh)*SQRT(FM))+One))          
    BoundCell(i)%TeS    = Tgn1
    BoundCell(i)%ThetaS = ThPotSurf
    BoundCell(i)%DragM  = DragM
    BoundCell(i)%DragH  = DragH
    BoundCell(i)%DragQ  = DragH

END SUBROUTINE DragCoeffLouisRad

SUBROUTINE DragCoeffLouis(i)
  
  INTEGER :: i
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,dL,zRauh,zRauhT,zPL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn1,Tgn2
  REAL(RealKind) :: RhoLoc,ThLoc,ThPotSurf,Fpot
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: Qg,Qc,Qdir,Qem,Qref,Qdif
  REAL(RealKind) :: RiB0,RiB1,RiBT,R,FM,FH,cM,cH
  REAL(RealKind) :: N,dFMdRiB0,dFHdRiB0,dRiBTdTgn,dNdTgn,dNdRiB0,dQgdTgn
  REAL(RealKind) :: F,dDragHdRiB0,dThPotSurfdTgn,dQcdTgn,dQcdRiB0,dQemdTgn,dFdTgn,dFdRiB0
  REAL(RealKind) :: ustar,thstar,DragM,DragH
  REAL(RealKind) :: det
  
  REAL(RealKind) ,Parameter :: zero=0.d0
  
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    dL     = BoundCell(i)%dL
    zRauh  = BoundCell(i)%zRauh
    zRauhT = BoundCell(i)%zRauhT
    alb    = BoundCell(i)%alb
    ee     = BoundCell(i)%ee
    shad   = BoundCell(i)%shad
    RhoLoc = Rho(ix,iy,iz,1)
    Tgn1   = BoundCell(i)%TeS
    ThLoc  = Th(ix,iy,iz,1)/RhoLoc
    DragM  = BoundCell(i)%DragM
    DragH  = BoundCell(i)%DragH
    
!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    
    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))         

!   Evaluate height z
    zPL  = Half*VolC(ix,iy,iz)/FL
    zPL  = dL                      !Oswald
    
!   First Estimate for Bulk Richardson Number RiB0 = RiBT
!    Fpot=(p0/(RhoLoc*Rd))**kappa
!    ThPotSurf=Tgn1**(Cv/Cp)*Fpot
    Fpot      = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
    ThPotSurf = Tgn1*Fpot
  
    RiBT      = Grav*zPL*(ThLoc-thPotSurf)/(ThLoc*(VT**2.0d0)+Eps)
    RiB0      = RiBT
    
    DO it=0,5      

!     Louis scheme
!---------------------------------------------------------------
!     Parametric Functions
      IF (RiB0.GT.0.25) THEN  ! critical Richardson number: 0.25
!     Stable Conditions
        FM = One/((One+4.7d0*RiB0)**2.0d0)
        FH = FM
        dFMdRiB0 = -Two*4.7d0/(One+4.7d0*RiB0)**Three
        dFHdRiB0  = dFMdRiB0  
      ELSE 
!     Unstable Conditions
        cM = 7.4d0*(Karm**2.0d0)*9.4d0*SQRT(zPL/(zRauh+Eps))/((LOG(zPL/(zRauh+Eps)))**2.0d0+Eps)
        cH = 5.3d0*cM/7.4d0
        FM = One-9.4d0*RiB0/(One+cM*SQRT(ABS(RiB0)))
        FH = One-9.4d0*RiB0/(One+cH*SQRT(ABS(RiB0)))
        dFMdRiB0 = -9.4d0*(One+cM*SQRT(ABS(RiB0))-0.5d0*RiB0*cM/(SQRT(ABS(RiB0))+Eps))/((One+cM*SQRT(ABS(RiB0)))**Two)
        dFHdRiB0 = -9.4d0*(One+cH*SQRT(ABS(RiB0))-0.5d0*RiB0*cH/(SQRT(ABS(RiB0))+Eps))/((One+cH*SQRT(ABS(RiB0)))**Two)  
      END IF

!     Turbulent Prandtl Number (after Businger)
      R = 7.4d-1
  
!     ThPotSurf=Tgn1**(Cv/Cp)*Fpot  
!      dThPotSurfdTgn=Fpot*(Cv/Cp)*Tgn1**(Cv/Cp-1.0d0)       
      dThPotSurfdTgn = Fpot     
 
!     Bulk Richardson Number  
      RiBT = Grav*zPL*(ThLoc-thPotSurf)/(ThLoc*(VT**Two)+Eps)
!!      dRiBTdTgn = -Grav*zPL*dThPotSurfdTgn/(ThLoc*(VT**Two)+Eps)

!     Function for Bulk Richardson Number
      N = RiB0*R*LOG(zRauh/(zRauhT+Eps))+(RiB0-RiBT)*R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH
!!      dNdTgn  = -dRiBTdTgn*R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH
      dNdRiB0 = R*LOG(zRauh/(zRauhT+Eps))+R*LOG(zPL/(zRauh+Eps))*SQRT(FM)/FH+ &
                (RiB0-RiBT)*R*LOG(zPL/(zRauh+Eps))*(0.5d0*FM**(-0.5d0)*dFMdRiB0*FH-SQRT(FM)*dFHdRiB0)/FH**Two   
  
!     Drag coefficient for heat as initial condition for radiation scheme
      DragH = (Karm/(LOG(zPL/zRauh)+Eps))**Two/(R*(LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*SQRT(FM))+One/FH))
      dDragHdRiB0 = -(Karm/(LOG(zPL/zRauh)+Eps))**Two/(R*(LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*SQRT(FM))+One/FH)**Two) &
                    *(-0.5d0*dFMdRiB0*LOG(zRauh/zRauhT)/(LOG(zPL/zRauh)*FM**(Three/Two))-dFHdRiB0/(FH**Two))

!     Radiation scheme to calculate Tgn1 and associated ThPotSurf
!----------------------------------------------------------------
!     Bodenwaermestrom   
!!      Qg = -BoundCell(i)%D*(Tgn1-BoundCell(i)%T1)/0.5d0 
!!      dQgdTgn  = -BoundCell(i)%D/0.5d0
!     Waermestrom in die Luft
!!      Qc = -RhoLoc*Cpd*DragH*VT*(ThPotSurf-ThLoc)
!!      dQcdTgn  = -RhoLoc*Cpd*DragH*VT*dThPotSurfdTgn
!!      dQcdRiB0 = -RhoLoc*Cpd*dDragHdRiB0*VT*(ThPotSurf-ThLoc)  

!     direkte Strahlung
!!      Qdir = (One-alb)*shad*MAX(Zero,(n1*radn1+n2*radn2+n3*radn3))*BoundCell(i)%raddirekt  
!     Qdir = 0.0d0
  
!     diffuse Strahlung
!!      Qdif = (ONE-alb)*BoundCell(i)%raddiffus 
!     Qdif = 0.0d0

!     Ausstrahlung im Infraroten   
!!      Qem  = -ee*5.6704e-8*Tgn1**4.0d0
!!      dQemdTgn = -ee*4.0d0*5.6704e-8*Tgn1**3.0d0
!!      Qref = 0.0d0     

!     Function for temperature
!!      F = Qg+Qc+Qdir+Qdif+Qem+Qref           ! Heat balance equation
!!      dFdTgn  = dQcdTgn+dQemdTgn+dQgdTgn
!!      dFdRiB0 = dQcdRiB0
    
!     Determinant
!!      det = dNdRiB0*dFdTgn-dNdTgn*dFdRiB0
      det = dNdRiB0

!     Numerical Iteration Bulk Richardson Number and temperature at surface 
      RiB1 = RiB0-N/(det+Eps)
      
!     IF(RiB0.LT.0.25d0.AND.RiB1.GE.0.25d0) THEN       !!! CRITICAL RICHARDSON NUMBER PROBLEM
!       RiB0 = 0.25d0
!       GOTO 1
!     ELSE IF(RiB0.GT.0.25d0.AND.RiB1.LE.0.25d0) THEN
!       RiB0 = 0.2500001d0
!       GOTO 1
!     END IF
      
!!      Tgn2 = Tgn1-(dNdRiB0*F-dFdRiB0*N)/(det+Eps)

!     Friction velocity for momentum and heat        
      ustar  = SQRT(FM*(Karm*VT/(LOG(zPL/zRauh)+Eps))**Two)
!     thstar = Karm**Two*VT*(ThLoc-ThPotSurf)*FH/(LOG(zPL/zRauh)*LOG(zPL/zRauhT)*R*ustar+Eps)

!     Drag coefficient for momentum and heat    
      DragM  = (ustar/(VT+Eps))**Two
!!      DragH  = (Karm/(LOG(zPL/zRauh)+Eps))**Two*FH/(R*(LOG(zRauh/zRauhT)*FH/(LOG(zPL/zRauh)*SQRT(FM))+One))          
            
!!      IF(ABS(RiB1-RiB0)<=10e-6*ABS(RiB1).AND.ABS(Tgn2-Tgn1)<=1.0d-6*Tgn2) THEN  
      IF(ABS(RiB1-RiB0)<=10e-6*ABS(RiB1)) THEN  
        RiB0 = RiB1
!!        Tgn1 = Tgn2
!!        ThPotSurf = Tgn1*Fpot
        EXIT
      ELSE
        RiB0 = RiB1
!!        Tgn1 = Tgn2
!!        ThPotSurf = Tgn1*Fpot
      END IF 
1   CONTINUE
    ENDDO
!!    BoundCell(i)%TeS    = Tgn1
!!    BoundCell(i)%ThetaS = ThPotSurf
    BoundCell(i)%DragM  = DragM
    BoundCell(i)%DragH  = DragH
    BoundCell(i)%DragQ  = DragH

END SUBROUTINE DragCoeffLouis

SUBROUTINE DragCoeffSea(i)
  
  INTEGER :: i
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,dL,zRauh1,zRauh2,zRauhT,zPL,zRauhf
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN,U10
  REAL(RealKind) :: Tgn1
  REAL(RealKind) :: RhoLoc,ThLoc,ThPotSurf,Fpot
  REAL(RealKind) :: RiBT,R,FM,FH,cM,cH
  REAL(RealKind) :: ustar,wstar,DragM,DragH,Dragq
  REAL(RealKind) :: Rr
  REAL(RealKind) :: eSurface,RhoDryLoc,pLoc
  
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    dL     = BoundCell(i)%dL
    zRauh1 = BoundCell(i)%zRauh
    zRauhT = BoundCell(i)%zRauhT
    RhoLoc = Rho(ix,iy,iz,1)
    Tgn1   = BoundCell(i)%TeS
    ThLoc  = Th(ix,iy,iz,1)/RhoLoc
    
!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc+Eps))
    
    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))         

!   Turbulent Prandtl Number (after Businger)
    R = 7.4d-1
  
!   Evaluate height z
    zPL  = dL                      !Oswald
    
    Fpot      = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa+Eps))
    ThPotSurf = Tgn1*Fpot
  
!   First Estimate for Bulk Richardson Number RiB0 = RiBT
    RiBT      = Grav*zPL*(ThLoc-ThPotSurf)/(ThLoc*(VT**Two)+Eps)

    IF (RiBT.GT.0.25) THEN  ! critical Richardson number: 0.25
!   Stable Conditions
      FM = One/((One+4.7d0*RiBT)**2.0d0+Eps)
      FH = FM
    ELSE 
!   Unstable Conditions
      cM = 7.4d0*(Karm**2.0d0)*9.4d0*SQRT(zPL/(zRauh1+Eps))/((LOG(zPL/(zRauh1+Eps)))**2.0d0+Eps)
      cH = 5.3d0*cM/7.4d0
      FM = One-9.4d0*RiBT/(One+cM*SQRT(ABS(RiBT))+Eps)
      FH = One-9.4d0*RiBT/(One+cH*SQRT(ABS(RiBT))+Eps)
    END IF

    DO it=0,5      

!     Friction velocity for momentum
      ustar  = SQRT(FM*(Karm*VT/(LOG(zPL/zRauh1)+Eps))**Two)

      IF (RiBT.GT.0.25) THEN
        wstar=Zero
      ELSE
        wstar = SQRT((Karm**Two/(R*LOG(zPl/zRauh1)*LOG(zpl/zRauhT))*VT**3*ABS(RiBT)*FH)**(3/2))
        ustar = MAX(ustar,wstar)
      END IF

!     New roughness lenght for momentum (after Charnock 1955,Nikuradse 1933)
      zRauh2=1.1d-2*ustar**Two/Grav
      zRauh2= MAX(zRauh2,1d-5)

      zRauhT=zRauh1
      zRauhT=MIN(zRauhT,0.1d0)

      IF(ABS(zRauh1-zRauh2)<=10e-6*ABS(zRauh2).AND.zRauh1.GT.1d-5) THEN  
        zRauh1=zRauh2
        EXIT
      ELSE
        zRauh1=zRauh2
      END IF 
1   CONTINUE
    ENDDO

!   Friction velocity for momentum        
    ustar  = SQRT(FM*(Karm*VT/(LOG(zPL/zRauh1)+Eps))**Two)
    U10  = LOG(10d0/zRauh1)*ustar/(SQRT(FM)*Karm)


!   Drag coefficient for momentum, heat and moisture
    DragM  = (ustar/(VT+Eps))**Two
    DragH  = (Karm/(LOG(zPL/zRauh1)+Eps))**Two*FH/(R*(LOG(zRauh1/zRauhT)*FH/(LOG(zPL/zRauh1)*SQRT(FM))+One))          

    BoundCell(i)%U10    = U10
    BoundCell(i)%DragM  = DragM
    BoundCell(i)%DragH  = DragH
    BoundCell(i)%DragQ  = DragH
    BoundCell(i)%zRauh  = zRauh1
    BoundCell(i)%zRauhT = zRauhT

    eSurface=SaturVapor(Tgn1)
    RhoDryLoc=RhoLoc-RhoV(ix,iy,iz,1)-RhoL(ix,iy,iz,1)-RhoR(ix,iy,iz,1)-RhoI(ix,iy,iz,1)-RhoS(ix,iy,iz,1)
    pLoc=PressureTheta(RhoDryLoc,RhoV(ix,iy,iz,1),RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1),RhoI(ix,iy,iz,1)&
    &+RhoS(ix,iy,iz,1),Th(ix,iy,iz,1))+Eps
    BoundCell(i)%qv=0.75*eSurFace/(Rv*Tgn1) ! Barthel_temp

END SUBROUTINE DragCoeffSea

SUBROUTINE Drag(Vector,VelF,Rhs,UVec,Time)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  TYPE(Vector4Cell_T), OPTIONAL :: UVec
  REAL(RealKind) :: Time

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  IF (PRESENT(UVec)) THEN
    uCL=>UVec%Vec(uPosL)%c
    vCL=>UVec%Vec(vPosL)%c
    wCL=>UVec%Vec(wPosL)%c
    uCR=>UVec%Vec(uPosR)%c
    vCR=>UVec%Vec(vPosR)%c
    wCR=>UVec%Vec(wPosR)%c
  ELSE
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
  END IF  
  uRhsL=>Rhs%Vec(uPosL)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsR=>Rhs%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c
  RhoV=>Vector%Vec(RhoVPos)%c
  RhoL=>Vector%Vec(RhoCPos)%c
  RhoR=>Vector%Vec(RhoRPos)%c
  RhoI=>Vector%Vec(RhoIPos)%c
  RhoS=>Vector%Vec(RhoSPos)%c
  Th=>Vector%Vec(ThPos)%c
  T=>TAbsCell(ibLoc)%Vec(1)%c
  ThRhs=>Rhs%Vec(ThPos)%c
  RhoRhs=>Rhs%Vec(RhoPos)%c
  RhoVRhs=>Rhs%Vec(RhoVPos)%c
  IF (RhoCPos>0) THEN
    RhoL=>Vector%Vec(RhoCPos)%c
  ELSE
    RhoL=>RhoLCell(ibLoc)%c
  END IF
  CALL DragCompute(Time)

END SUBROUTINE Drag

SUBROUTINE DragCompute(Time)

  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: DragM,DragH,DragQ
  REAL(RealKind) :: RhoDryLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoILoc
  REAL(RealKind) :: cpml,dPotdt,rrv,eSurface,drvdt,MoistFlux
  REAL(RealKind) :: Rm,PrePi,PotM,PLoc,rrl
  REAL(RealKind) :: TotalFlux,SensFlux,TSurface,RhoVSurface
  REAL(RealKind) :: SensFluxOut,LatFluxOut
  REAL(RealKind) :: SumSensFlux,SumLatFlux
  CHARACTER(4) :: MyIDCHAR

  SumSensFlux=Zero
  SumLatFlux=Zero

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    DragH  = BoundCell(i)%DragH
    DragQ  = BoundCell(i)%DragQ
    DragM  = BoundCell(i)%DragM
    RhoVSurface = BoundCell(i)%qv
    TSurface    = BoundCell(i)%TeS
    RhoLoc      = Rho(ix,iy,iz,1)
    RhoVLoc     = RhoV(ix,iy,iz,1)
    RhoLLoc     = RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
    RhoILoc     = RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
    RhoDryLoc   = RhoLoc-RhoVLoc-RhoLLoc-RhoILoc
    pLoc        = PressureTheta(RhoDryLoc,RhoVLoc,RhoLLoc,RhoILoc,Th(ix,iy,iz,1))+Eps

    uLoc=Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoLoc+Eps)
    vLoc=Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoLoc+Eps)
    wLoc=Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoLoc+Eps)
    VN=uLoc*n1+vLoc*n2+wLoc*n3
    V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
    VT=SQRT(MAX(V-VN*VN,Zero))

    IF (MomFlux) THEN
      uLoc=uCL(ix,iy,iz,1)/(RhoLoc+Eps)
      vLoc=vCL(ix,iy,iz,1)/(RhoLoc+Eps)
      wLoc=wCL(ix,iy,iz,1)/(RhoLoc+Eps)
      VN=uLoc*n1+vLoc*n2+wLoc*n3
      V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
!     uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)- &
!                    RhoLoc*FL*DragM*VT*(uLocL-n1*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(uLoc-n1*VN)
      CALL DistributeFlux(uRhsL(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!     vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)- &
!                  RhoLoc*FL*DragM*VT*(vLocL-n2*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(vLoc-n2*VN)
      CALL DistributeFlux(vRhsL(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!     wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)- &
!                    RhoLoc*FL*DragM*VT*(wLocL-n3*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(wLoc-n3*VN)
      CALL DistributeFlux(wRhsL(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
      uLoc=uCR(ix,iy,iz,1)/(RhoLoc+Eps)
      vLoc=vCR(ix,iy,iz,1)/(RhoLoc+Eps)
      wLoc=wCR(ix,iy,iz,1)/(RhoLoc+Eps)
      VN=uLoc*n1+vLoc*n2+wLoc*n3
      V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
!     uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)- &
!                    RhoLoc*FL*DragM*VT*(uLocR-n1*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(uLoc-n1*VN)
      CALL DistributeFlux(uRhsR(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!     vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)- &
!                    RhoLoc*FL*DragM*VT*(vLocR-n2*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(vLoc-n2*VN)
      CALL DistributeFlux(vRhsR(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!     wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)- &
!                    RhoLoc*FL*DragM*VT*(wLocR-n3*VN)
      TotalFlux=-RhoLoc*FL*DragM*VT*(wLoc-n3*VN)
      CALL DistributeFlux(wRhsR(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
    END IF               
    IF (ScalarFlux) THEN
      rrv=RhoVLoc/RhoDryLoc
      rrl=RhoLLoc/RhoDryLoc
      Rm=Rd+rrv*Rv+Eps
      Cpml=Cpd+Cpv*rrv+Cpl*rrl+Eps
      IF (SensFluxFix) THEN
        SensFlux=SensibleHeatFlux(Time,BoundCell(i)%LandClass,ix,iy,iz)*FL/(RhoLoc*Cpd+Eps)  ! [m3 K s-1]
      ELSE
        SensFlux=-FL*(DragH*VT+1.d-6)*(T(ix,iy,iz,1)-TSurface) ! [m3 K s-1]
      END IF  
      IF (RhoVPos>0) THEN
        IF (LatFluxFix) THEN
          MoistFlux=LatentHeatFlux(Time,BoundCell(i)%LandClass,ix,iy,iz)*FL/L0  ! [kg s-1]
        ELSE
          MoistFlux=-FL*(DragQ*VT+1.d-6)*(RhoV(ix,iy,iz,1)-RhoVSurface) ! [kg s-1]
        END IF  
      ELSE
        MoistFlux=Zero
      END IF  

      PrePi=(pLoc/p0)**(Rm/Cpml)
      TotalFlux=Th(ix,iy,iz,1)*(SensFlux/(T(ix,iy,iz,1)+Eps) + &              ! sens.HeatFlux
                  MoistFlux/(RhoLoc+Eps)*(Rv/Rm-LOG(PrePi)*(Rv/Rm-Cpv/Cpml))) ! lat.HeatFlux
!     ThRhs(ix,iy,iz,1)=ThRhs(ix,iy,iz,1)+TotalFlux
      CALL DistributeFlux(ThRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!     RhoVRhs(ix,iy,iz,1)=RhoVRhs(ix,iy,iz,1)+MoistFlux                ! lat.HeatFlux 
      CALL DistributeFlux(RhoVRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,MoistFlux)
!     RhoRhs(ix,iy,iz,1)=RhoRhs(ix,iy,iz,1)+MoistFlux                 ! lat.HeatFlux
      CALL DistributeFlux(RhoRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                         ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                         ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,MoistFlux)

      SensFluxOut = SensFlux *Cpd*RhoLoc/FL       ! [kg s-3 = W m-2]
      LatFluxOut  = MoistFlux*L0/FL               ! [kg s-3 = W m-2]
      
      SumSensFlux = SumSensFlux+SensFluxOut
      SumLatFlux  = SumLatFlux +LatFluxOut
    END IF  
  END DO

END SUBROUTINE DragCompute

SUBROUTINE JacDrag(Vector,VelF,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vec4_T), POINTER :: Jac(:)

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c
  Th=>Vector%Vec(ThPos)%c
  AS=>Jac

  CALL JacDragCompute

END SUBROUTINE JacDrag

SUBROUTINE BaumDrag(Vector,VelF,Rhs,UVec)

! Analog zu Drag

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  TYPE(Vector4Cell_T), OPTIONAL :: UVec

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  IF (PRESENT(UVec)) THEN
    uCL=>UVec%Vec(uPosL)%c
    vCL=>UVec%Vec(vPosL)%c
    wCL=>UVec%Vec(wPosL)%c
    uCR=>UVec%Vec(uPosR)%c
    vCR=>UVec%Vec(vPosR)%c
    wCR=>UVec%Vec(wPosR)%c
  ELSE
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
  END IF
  uRhsL=>Rhs%Vec(uPosL)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsR=>Rhs%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c

  CALL BaumDragCompute

END SUBROUTINE BaumDrag

SUBROUTINE BaumDragF(Vector,VelF,RhsF)

! Analog zu Drag

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(VelocityFace_T) :: RhsF

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uFRhs=>RhsF%uF
  vFRhs=>RhsF%vF
  wFRhs=>RhsF%wF
  Rho=>Vector%Vec(RhoPos)%c

  CALL BaumDragComputeF

END SUBROUTINE BaumDragF

SUBROUTINE JacBaumDrag(Vector,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c
  Th=>Vector%Vec(ThPos)%c
  AS=>Jac

  CALL JacBaumDragCompute

END SUBROUTINE JacBaumDrag

SUBROUTINE TrafficDrag(Vector,VelF,Rhs,UVec)

! Analog zu Drag

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  TYPE(Vector4Cell_T), OPTIONAL :: UVec

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  IF (PRESENT(UVec)) THEN
    uCL=>UVec%Vec(uPosL)%c
    vCL=>UVec%Vec(vPosL)%c
    wCL=>UVec%Vec(wPosL)%c
    uCR=>UVec%Vec(uPosR)%c
    vCR=>UVec%Vec(vPosR)%c
    wCR=>UVec%Vec(wPosR)%c
  ELSE
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
  END IF
  uRhsL=>Rhs%Vec(uPosL)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsR=>Rhs%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c

  CALL TrafficDragCompute

END SUBROUTINE TrafficDrag

SUBROUTINE JacTrafficDrag(Vector,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c
  AS=>Jac

  CALL JacTrafficDragCompute

END SUBROUTINE JacTrafficDrag
SUBROUTINE NoSlipCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: uLoc,vLoc,wLoc,RhoLoc,FL,DFFV
  REAL(RealKind) ,Parameter :: Two=2.0d0

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    FL     = BoundCell(i)%FL
    RhoLoc = Rho(ix,iy,iz,1)
    
    DFFV   = Two*D(ix,iy,iz,1)*FL*FL/(VolC(ix,iy,iz)*RhoLoc)
    
    uLoc   = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoLoc+Eps)
    vLoc   = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoLoc+Eps)
    wLoc   = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoLoc+Eps)
   
    uRhsL(ix,iy,iz,1) = uRhsL(ix,iy,iz,1)- &
                        DFFV*uLoc
    vRhsL(ix,iy,iz,1) = vRhsL(ix,iy,iz,1)- &
                        DFFV*vLoc
    wRhsL(ix,iy,iz,1) = wRhsL(ix,iy,iz,1)- &
                        DFFV*wLoc
    uRhsR(ix,iy,iz,1) = uRhsR(ix,iy,iz,1)- &
                        DFFV*uLoc
    vRhsR(ix,iy,iz,1) = vRhsR(ix,iy,iz,1)- &
                        DFFV*vLoc
    wRhsR(ix,iy,iz,1) = wRhsR(ix,iy,iz,1)- &
                        DFFV*wLoc
  END DO
  
END SUBROUTINE NoSlipCompute

SUBROUTINE JacNoSlipCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FL,DFV,RhoLoc
  REAL(RealKind) ,Parameter :: Two=2.0d0

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    FL     = BoundCell(i)%FL
    RhoLoc = Rho(ix,iy,iz,1)
    
    DFV = Two*D(ix,iy,iz,1)*FL*FL/(VolC(ix,iy,iz)*VolC(ix,iy,iz)*RhoLoc)
    
    AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
                                                    DFV
    AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1)- &
                                                    DFV
    AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1)- &
                                                    DFV
   
    AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
                                                    DFV
    AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1)- &
                                                    DFV
    AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1)- &
                                                    DFV
  END DO
END SUBROUTINE JacNoSlipCompute

SUBROUTINE JacDragCompute
  
  INTEGER :: i
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN  
  REAL(RealKind) :: FL,DragM,DragH
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: TotalFlux,SensFlux
  

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    DragM  = BoundCell(i)%DragM
    DragH  = BoundCell(i)%DragH
    RhoLoc = Rho(ix,iy,iz,1)

!   Evaluate velocity tangential and normal  
    uLoc=Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoLoc+Eps)
    vLoc=Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoLoc+Eps)
    wLoc=Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoLoc+Eps)
    VN=uLoc*n1+vLoc*n2+wLoc*n3
    V=uLoc*uLoc+vLoc*vLoc+wLoc*wLoc
    VT=SQRT(MAX(V-VN*VN,Zero))
    IF (MomFlux) THEN
      TotalFlux=-RhoLoc*FL*DragM*VT
      CALL JacDistributeMomFlux(AS(IndexMet(uPosLJac,uPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(uPosLJac,vPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(uPosLJac,wPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosLJac,uPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosLJac,vPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosLJac,wPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosLJac,uPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosLJac,vPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosLJac,wPosLJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                               ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1) &
                               ,FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
      CALL JacDistributeMomFlux(AS(IndexMet(uPosRJac,uPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(uPosRJac,vPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(uPosRJac,wPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosRJac,uPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosRJac,vPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(vPosRJac,wPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosRJac,uPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosRJac,vPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,AS(IndexMet(wPosRJac,wPosRJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                               ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                               ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1) &
                               ,FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
    END IF

    IF (ScalarFlux) THEN
      IF (.NOT.SensFluxFix) THEN
!       SensFlux=-FL*(DragH*VT+1.d-6)*(T(ix,iy,iz,1)-TSurface) ! [m3 K s-1]
        SensFlux=-FL*(DragH*VT+1.d-6)
        CALL JacDistributeScalarFlux(AS(IndexMet(ThPosJac,ThPosJac))%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                                    ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                                    ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1) &
                                    ,FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,SensFlux)
      END IF                                                  
    END IF  
  END DO
END SUBROUTINE JacDragCompute

SUBROUTINE DragScalarCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3,FL
  REAL(RealKind) :: DragH,ThetaS
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: TotalFlux
  REAL(RealKind) :: dLoc,Surface,VolLoc,VolDiff
    
  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    DragH  = BoundCell(i)%DragH
    ThetaS = BoundCell(i)%ThetaS
    RhoLoc = Rho(ix,iy,iz,1)
!   ThRhs(ix,iy,iz,1)=ThRhs(ix,iy,iz,1)+ &
!                  FL*DragH*(RhoLoc*ThetaS-Th(ix,iy,iz,1))
    TotalFlux=FL*DragH*(RhoLoc*ThetaS-Th(ix,iy,iz,1))
    CALL DistributeFlux(ThRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                       ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                       ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!   TotalFlux=FL*DragH*(RhoLoc*ThetaS-Th(ix,iy,iz,1))
!   dLoc=ABS(n1)*dx(ix)+ABS(n2)*dy(iy)+ABS(n3)*dz(iz)
!   VolLoc=dx(ix)*dy(iy)*dz(iz)
!   VolLoc=MIN(FL*dLoc,VolLoc)
!   Surface=ABS(FU(ix-1,iy,iz)-FU(ix,iy,iz)) &
!          +ABS(FV(ix,iy-1,iz)-FV(ix,iy,iz)) &
!          +ABS(FW(ix,iy,iz-1)-FW(ix,iy,iz))
!   ThRhs(ix,iy,iz,1)=ThRhs(ix,iy,iz,1)+MIN(VolC(ix,iy,iz)/VolLoc,One)*TotalFlux
!   VolDiff=MAX(VolLoc-VolC(ix,iy,iz),Zero)
!   ThRhs(ix-1,iy,iz,1)=ThRhs(ix-1,iy,iz,1)+MAX(FU(ix-1,iy,iz)-FU(ix,iy,iz),Zero)/Surface*VolDiff/VolLoc*TotalFlux
!   ThRhs(ix+1,iy,iz,1)=ThRhs(ix+1,iy,iz,1)+MAX(FU(ix,iy,iz)-FU(ix-1,iy,iz),Zero)/Surface*VolDiff/VolLoc*TotalFlux
!   ThRhs(ix,iy-1,iz,1)=ThRhs(ix,iy-1,iz,1)+MAX(FV(ix,iy-1,iz)-FV(ix,iy,iz),Zero)/Surface*VolDiff/VolLoc*TotalFlux
!   ThRhs(ix,iy+1,iz,1)=ThRhs(ix,iy+1,iz,1)+MAX(FV(ix,iy,iz)-FV(ix,iy-1,iz),Zero)/Surface*VolDiff/VolLoc*TotalFlux
!   ThRhs(ix,iy,iz-1,1)=ThRhs(ix,iy,iz-1,1)+MAX(FW(ix,iy,iz-1)-FW(ix,iy,iz),Zero)/Surface*VolDiff/VolLoc*TotalFlux
!   ThRhs(ix,iy,iz+1,1)=ThRhs(ix,iy,iz+1,1)+MAX(FW(ix,iy,iz)-FW(ix,iy,iz-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux

  END DO

END SUBROUTINE DragScalarCompute

SUBROUTINE DistributeFlux(Rhs,FU,FV,FW,Vol,FL,dx,dy,dz,n1,n2,n3,TotalFlux)
  REAL(RealKind) :: Rhs(-1:1,-1:1,-1:1)
  REAL(RealKind) :: FU(-1:0),FV(-1:0),FW(-1:0),Vol(-1:1,-1:1,-1:1),FL
  REAL(RealKind) :: dx,dy,dz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: TotalFlux

  REAL(RealKind) :: dLoc,VolLoc,VolDiff,SurFace

  IF (FluxDistribute) THEN
    dLoc=ABS(n1)*dx+ABS(n2)*dy+ABS(n3)*dz
    VolLoc=dx*dy*dz
!   VolLoc=MIN(FL*dLoc,VolLoc)+Eps
    Surface=ABS(FU(-1)-FU(0)) &
           +ABS(FV(-1)-FV(0)) &
           +ABS(FW(-1)-FW(0))+Eps
    Rhs(0,0,0)=Rhs(0,0,0)+MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux
    VolDiff=MAX(VolLoc-Vol(0,0,0),Zero)
    Rhs(-1,0,0)=Rhs(-1,0,0)+MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(1,0,0)=Rhs(1,0,0)+MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,-1,0)=Rhs(0,-1,0)+MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,1,0)=Rhs(0,1,0)+MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,0,-1)=Rhs(0,0,-1)+MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,0,1)=Rhs(0,0,1)+MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
  ELSE  
    Rhs(0,0,0)=Rhs(0,0,0)+TotalFlux
  END IF  
END SUBROUTINE DistributeFlux

SUBROUTINE JacDistributeMomFlux(JacUU,JacUV,JacUW &
                               ,JacVU,JacVV,JacVW &
                               ,JacWU,JacWV,JacWW &
                               ,FU,FV,FW,Vol,FL,dx,dy,dz,n1,n2,n3,TotalFlux)
  REAL(RealKind) :: JacUU(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacUV(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacUW(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacVU(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacVV(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacVW(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacWU(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacWV(-1:1,-1:1,-1:1)
  REAL(RealKind) :: JacWW(-1:1,-1:1,-1:1)
  REAL(RealKind) :: FU(-1:0),FV(-1:0),FW(-1:0),Vol(-1:1,-1:1,-1:1),FL
  REAL(RealKind) :: dx,dy,dz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: TotalFlux

  REAL(RealKind) :: dLoc,VolLoc,VolDiff,SurFace,Fac

  IF (FluxDistribute) THEN
    dLoc=ABS(n1)*dx+ABS(n2)*dy+ABS(n3)*dz
    VolLoc=dx*dy*dz
!   VolLoc=MIN(FL*dLoc,VolLoc)+Eps
    Surface=ABS(FU(-1)-FU(0)) &
           +ABS(FV(-1)-FV(0)) &
           +ABS(FW(-1)-FW(0))+Eps
!   Rhs(0,0,0)=Rhs(0,0,0)+MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux
    Fac=MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux/Vol(0,0,0)
    JacUU(0,0,0)=JacUU(0,0,0)+(One-n1*n1)*Fac
    JacUV(0,0,0)=JacUV(0,0,0)+(    n1*n2)*Fac
    JacUW(0,0,0)=JacUW(0,0,0)+(    n1*n3)*Fac
    JacVU(0,0,0)=JacVU(0,0,0)+(    n2*n1)*Fac
    JacVV(0,0,0)=JacVV(0,0,0)+(One-n2*n2)*Fac
    JacVW(0,0,0)=JacVW(0,0,0)+(    n2*n3)*Fac
    JacWU(0,0,0)=JacWU(0,0,0)+(    n3*n1)*Fac
    JacWV(0,0,0)=JacWV(0,0,0)+(    n3*n2)*Fac
    JacWW(0,0,0)=JacWW(0,0,0)+(One-n3*n3)*Fac

    VolDiff=MAX(VolLoc-Vol(0,0,0),Zero)
!   Rhs(-1,0,0)=Rhs(-1,0,0)+MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(-1,0,0)+Eps)
    JacUU(-1,0,0)=JacUU(-1,0,0)+(One-n1*n1)*Fac
    JacUV(-1,0,0)=JacUV(-1,0,0)+(    n1*n2)*Fac
    JacUW(-1,0,0)=JacUW(-1,0,0)+(    n1*n3)*Fac
    JacVU(-1,0,0)=JacVU(-1,0,0)+(    n2*n1)*Fac
    JacVV(-1,0,0)=JacVV(-1,0,0)+(One-n2*n2)*Fac
    JacVW(-1,0,0)=JacVW(-1,0,0)+(    n2*n3)*Fac
    JacWU(-1,0,0)=JacWU(-1,0,0)+(    n3*n1)*Fac
    JacWV(-1,0,0)=JacWV(-1,0,0)+(    n3*n2)*Fac
    JacWW(-1,0,0)=JacWW(-1,0,0)+(One-n3*n3)*Fac
!   Rhs(1,0,0)=Rhs(1,0,0)+MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(1,0,0)+Eps)
    JacUU(1,0,0)=JacUU(1,0,0)+(One-n1*n1)*Fac
    JacUV(1,0,0)=JacUV(1,0,0)+(    n1*n2)*Fac
    JacUW(1,0,0)=JacUW(1,0,0)+(    n1*n3)*Fac
    JacVU(1,0,0)=JacVU(1,0,0)+(    n2*n1)*Fac
    JacVV(1,0,0)=JacVV(1,0,0)+(One-n2*n2)*Fac
    JacVW(1,0,0)=JacVW(1,0,0)+(    n2*n3)*Fac
    JacWU(1,0,0)=JacWU(1,0,0)+(    n3*n1)*Fac
    JacWV(1,0,0)=JacWV(1,0,0)+(    n3*n2)*Fac
    JacWW(1,0,0)=JacWW(1,0,0)+(One-n3*n3)*Fac
!   Rhs(0,-1,0)=Rhs(0,-1,0)+MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,-1,0)+Eps)
    JacUU(0,-1,0)=JacUU(0,-1,0)+(One-n1*n1)*Fac
    JacUV(0,-1,0)=JacUV(0,-1,0)+(    n1*n2)*Fac
    JacUW(0,-1,0)=JacUW(0,-1,0)+(    n1*n3)*Fac
    JacVU(0,-1,0)=JacVU(0,-1,0)+(    n2*n1)*Fac
    JacVV(0,-1,0)=JacVV(0,-1,0)+(One-n2*n2)*Fac
    JacVW(0,-1,0)=JacVW(0,-1,0)+(    n2*n3)*Fac
    JacWU(0,-1,0)=JacWU(0,-1,0)+(    n3*n1)*Fac
    JacWV(0,-1,0)=JacWV(0,-1,0)+(    n3*n2)*Fac
    JacWW(0,-1,0)=JacWW(0,-1,0)+(One-n3*n3)*Fac
!   Rhs(0,1,0)=Rhs(0,1,0)+MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,1,0)+Eps)
    JacUU(0,1,0)=JacUU(0,1,0)+(One-n1*n1)*Fac
    JacUV(0,1,0)=JacUV(0,1,0)+(    n1*n2)*Fac
    JacUW(0,1,0)=JacUW(0,1,0)+(    n1*n3)*Fac
    JacVU(0,1,0)=JacVU(0,1,0)+(    n2*n1)*Fac
    JacVV(0,1,0)=JacVV(0,1,0)+(One-n2*n2)*Fac
    JacVW(0,1,0)=JacVW(0,1,0)+(    n2*n3)*Fac
    JacWU(0,1,0)=JacWU(0,1,0)+(    n3*n1)*Fac
    JacWV(0,1,0)=JacWV(0,1,0)+(    n3*n2)*Fac
    JacWW(0,1,0)=JacWW(0,1,0)+(One-n3*n3)*Fac
!   Rhs(0,0,-1)=Rhs(0,0,-1)+MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,0,-1)+Eps)
    JacUU(0,0,-1)=JacUU(0,0,-1)+(One-n1*n1)*Fac
    JacUV(0,0,-1)=JacUV(0,0,-1)+(    n1*n2)*Fac
    JacUW(0,0,-1)=JacUW(0,0,-1)+(    n1*n3)*Fac
    JacVU(0,0,-1)=JacVU(0,0,-1)+(    n2*n1)*Fac
    JacVV(0,0,-1)=JacVV(0,0,-1)+(One-n2*n2)*Fac
    JacVW(0,0,-1)=JacVW(0,0,-1)+(    n2*n3)*Fac
    JacWU(0,0,-1)=JacWU(0,0,-1)+(    n3*n1)*Fac
    JacWV(0,0,-1)=JacWV(0,0,-1)+(    n3*n2)*Fac
    JacWW(0,0,-1)=JacWW(0,0,-1)+(One-n3*n3)*Fac
!   Rhs(0,0,1)=Rhs(0,0,1)+MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,0,1)+Eps)
    JacUU(0,0,1)=JacUU(0,0,1)+(One-n1*n1)*Fac
    JacUV(0,0,1)=JacUV(0,0,1)+(    n1*n2)*Fac
    JacUW(0,0,1)=JacUW(0,0,1)+(    n1*n3)*Fac
    JacVU(0,0,1)=JacVU(0,0,1)+(    n2*n1)*Fac
    JacVV(0,0,1)=JacVV(0,0,1)+(One-n2*n2)*Fac
    JacVW(0,0,1)=JacVW(0,0,1)+(    n2*n3)*Fac
    JacWU(0,0,1)=JacWU(0,0,1)+(    n3*n1)*Fac
    JacWV(0,0,1)=JacWV(0,0,1)+(    n3*n2)*Fac
    JacWW(0,0,1)=JacWW(0,0,1)+(One-n3*n3)*Fac
  ELSE  
!   Rhs(0,0,0)=Rhs(0,0,0)+TotalFlux
    Fac=TotalFlux/Vol(0,0,0)
    JacUU(0,0,0)=JacUU(0,0,0)+(One-n1*n1)*Fac
    JacUV(0,0,0)=JacUV(0,0,0)+(    n1*n2)*Fac
    JacUW(0,0,0)=JacUW(0,0,0)+(    n1*n3)*Fac
    JacVU(0,0,0)=JacVU(0,0,0)+(    n2*n1)*Fac
    JacVV(0,0,0)=JacVV(0,0,0)+(One-n2*n2)*Fac
    JacVW(0,0,0)=JacVW(0,0,0)+(    n2*n3)*Fac
    JacWU(0,0,0)=JacWU(0,0,0)+(    n3*n1)*Fac
    JacWV(0,0,0)=JacWV(0,0,0)+(    n3*n2)*Fac
    JacWW(0,0,0)=JacWW(0,0,0)+(One-n3*n3)*Fac
  END IF  
END SUBROUTINE JacDistributeMomFlux

SUBROUTINE JacDistributeScalarFlux(Jac,FU,FV,FW,Vol,FL,dx,dy,dz,n1,n2,n3,TotalFlux)
  REAL(RealKind) :: Jac(-1:1,-1:1,-1:1)
  REAL(RealKind) :: FU(-1:0),FV(-1:0),FW(-1:0),Vol(-1:1,-1:1,-1:1),FL
  REAL(RealKind) :: dx,dy,dz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: TotalFlux

  REAL(RealKind) :: dLoc,VolLoc,VolDiff,SurFace,Fac

  IF (FluxDistribute) THEN
    dLoc=ABS(n1)*dx+ABS(n2)*dy+ABS(n3)*dz
    VolLoc=dx*dy*dz
    VolLoc=MIN(FL*dLoc,VolLoc)
    Surface=ABS(FU(-1)-FU(0)) &
           +ABS(FV(-1)-FV(0)) &
           +ABS(FW(-1)-FW(0))
!   Rhs(0,0,0)=Rhs(0,0,0)+MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux
    Fac=MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux/Vol(0,0,0)
    Jac(0,0,0)=Jac(0,0,0)+Fac

    VolDiff=MAX(VolLoc-Vol(0,0,0),Zero)
!   Rhs(-1,0,0)=Rhs(-1,0,0)+MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(-1,0,0)+Eps)
    Jac(-1,0,0)=Jac(-1,0,0)+Fac
!   Rhs(1,0,0)=Rhs(1,0,0)+MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(1,0,0)+Eps)
    Jac(1,0,0)=Jac(1,0,0)+Fac
!   Rhs(0,-1,0)=Rhs(0,-1,0)+MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,-1,0)+Eps)
    Jac(0,-1,0)=Jac(0,-1,0)+Fac
!   Rhs(0,1,0)=Rhs(0,1,0)+MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,1,0)+Eps)
    Jac(0,1,0)=Jac(0,1,0)+Fac
!   Rhs(0,0,-1)=Rhs(0,0,-1)+MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,0,-1)+Eps)
    Jac(0,0,-1)=Jac(0,0,-1)+Fac
!   Rhs(0,0,1)=Rhs(0,0,1)+MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Fac=MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux/(Vol(0,0,1)+Eps)
    Jac(0,0,1)=Jac(0,0,1)+Fac
  ELSE  
!   Rhs(0,0,0)=Rhs(0,0,0)+TotalFlux
    Fac=TotalFlux/Vol(0,0,0)
    Jac(0,0,0)=Jac(0,0,0)+Fac
  END IF  
END SUBROUTINE JacDistributeScalarFlux

SUBROUTINE BaumDragCompute

! Abbremsung innerhalb der Baumschicht analog zu DragCompute,
! siehe zusaetzlich auch TkeDisBaumCompute,
! lineare Skalierung des Effekts ueber den MASK-Input-Parameter c = 0...1,
! wirksame Groessen innerhalb der Baumschicht:
! zRauhBaum, distBaum, FL=FW (Testwerte siehe in Domain_Mod),
! MASK-Input-Parameter TreePoint(i)%c=1 sollte max. Effekt entsprechen
!
! Keyword = Baum in Namelist ModelPhysics
! Input file = InputFileName.Tree (supplied by GRID model
!              analogously to Weight files)

  INTEGER :: NumPoints
  INTEGER :: ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: c
  REAL(RealKind) :: VAbsL,VAbsR
  REAL(RealKind), PARAMETER :: Cdv=1.0d0 !0.2d0 !(Sectional Drag Coefficient)
  REAL(RealKind), PARAMETER :: LD=1.0d0 !0.2d0  !(Leaf Area Density m^2/m^3)

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
    c=TreePoint(i)%c
!   Evaluate velocity tangential and normal  
    RhoLoc=Rho(ix,iy,iz,1)
    uLocL=uCL(ix,iy,iz,1)/(RhoLoc+Eps)
    uLocR=uCR(ix,iy,iz,1)/(RhoLoc+Eps)
    vLocL=vCL(ix,iy,iz,1)/(RhoLoc+Eps)
    vLocR=vCR(ix,iy,iz,1)/(RhoLoc+Eps)
    wLocL=wCL(ix,iy,iz,1)/(RhoLoc+Eps)
    wLocR=wCR(ix,iy,iz,1)/(RhoLoc+Eps)

    VAbsL=SQRT(uLocL*uLocL+vLocL*vLocL+wLocL*wLocL)
    VAbsR=SQRT(uLocR*uLocR+vLocR*vLocL+wLocR*wLocR)

    uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*uCL(ix,iy,iz,1)*VAbsL
    vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*vCL(ix,iy,iz,1)*VAbsL
    wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*wCL(ix,iy,iz,1)*VAbsL
    uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*uCR(ix,iy,iz,1)*VAbsR
    vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*vCR(ix,iy,iz,1)*VAbsR
    wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)-VolC(ix,iy,iz)*LD*CDv*c*wCR(ix,iy,iz,1)*VAbsR
  END DO

END SUBROUTINE BaumDragCompute
SUBROUTINE BaumDragComputeF

! Abbremsung innerhalb der Baumschicht analog zu DragCompute,
! siehe zusaetzlich auch TkeDisBaumCompute,
! lineare Skalierung des Effekts ueber den MASK-Input-Parameter c = 0...1,
! wirksame Groessen innerhalb der Baumschicht:
! zRauhBaum, distBaum, FL=FW (Testwerte siehe in Domain_Mod),
! MASK-Input-Parameter TreePoint(i)%c=1 sollte max. Effekt entsprechen
!
! Keyword = Baum in Namelist ModelPhysics
! Input file = InputFileName.Tree (supplied by GRID model
!              analogously to Weight files)

  INTEGER :: NumPoints
  INTEGER :: ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: uLoc
  REAL(RealKind) :: vLoc
  REAL(RealKind) :: wLoc
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: c
  REAL(RealKind) :: VAbs
  REAL(RealKind), PARAMETER :: Cdv=1.0d0 !0.2d0 !(Sectional Drag Coefficient)
  REAL(RealKind), PARAMETER :: LD=1.0d0 !0.2d0  !(Leaf Area Density m^2/m^3)

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
    c=TreePoint(i)%c
!   Evaluate velocity tangential and normal  
    RhoLoc=Rho(ix,iy,iz,1)+Eps
    uLoc=(FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+FU(ix,iy,iz)*uF(ix,iy,iz)) &
        /(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)/RhoLoc
    vLoc=(FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+FV(ix,iy,iz)*vF(ix,iy,iz)) &
        /(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)/RhoLoc
    wLoc=(FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+FW(ix,iy,iz)*wF(ix,iy,iz)) &
        /(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)/RhoLoc

    VAbs=SQRT(uLoc*uLoc+vLoc*vLoc+wLoc*wLoc)

    uFRhs(ix,iy,iz)=uFRhs(ix,iy,iz)-FU(ix,iy,iz)*LD*CDv*c*uLoc*VAbs/(FU(ix-1,iy,iz)+FU(ix,iy,iz))
    uFRhs(ix-1,iy,iz)=uFRhs(ix-1,iy,iz)-FU(ix-1,iy,iz)*LD*CDv*c*uLoc*VAbs/(FU(ix-1,iy,iz)+FU(ix,iy,iz))
    vFRhs(ix,iy,iz)=vFRhs(ix,iy,iz)-FV(ix,iy,iz)*LD*CDv*c*vLoc*VAbs/(FV(ix,iy-1,iz)+FV(ix,iy,iz))
    vFRhs(ix,iy-1,iz)=vFRhs(ix,iy-1,iz)-FV(ix,iy-1,iz)*LD*CDv*c*vLoc*VAbs/(FV(ix,iy-1,iz)+FV(ix,iy,iz))
    wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)-FW(ix,iy,iz)*LD*CDv*c*wLoc*VAbs/(FW(ix,iy,iz-1)+FW(ix,iy,iz))
    wFRhs(ix,iy,iz-1)=wFRhs(ix,iy,iz-1)-FW(ix,iy,iz-1)*LD*CDv*c*wLoc*VAbs/(FW(ix,iy,iz-1)+FW(ix,iy,iz))
  END DO

END SUBROUTINE BaumDragComputeF

SUBROUTINE JacBaumDragCompute

  INTEGER :: NumPoints
  INTEGER :: ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,logz
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: DragM
  REAL(RealKind) :: RhoLoc

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
!   AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n1*n1)*VT + (uLoc-n1*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n1*n2)*VT + (uLoc-n1*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(uPosLJac,wPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosLJac,wPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n1*n3)*VT + (uLoc-n1*VN)*(wLoc-n3*VN)/(VT+Eps) )

!   AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n2*n1)*VT + (vLoc-n2*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n2*n2)*VT + (vLoc-n2*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(vPosLJac,wPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosLJac,wPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n2*n3)*VT + (vLoc-n2*VN)*(wLoc-n3*VN)/(VT+Eps) )

!   AS(IndexMet(wPosLJac,uPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosLJac,uPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n3*n1)*VT + (wLoc-n3*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(wPosLJac,vPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosLJac,vPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n3*n2)*VT + (wLoc-n3*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n3*n3)*VT + (wLoc-n3*VN)*(wLoc-n3*VN)/(VT+Eps) )


!   AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n1*n1)*VT + (uLoc-n1*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n1*n2)*VT + (uLoc-n1*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(uPosRJac,wPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(uPosRJac,wPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n1*n3)*VT + (uLoc-n1*VN)*(wLoc-n3*VN)/(VT+Eps) )

!   AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n2*n1)*VT + (vLoc-n2*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n2*n2)*VT + (vLoc-n2*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(vPosRJac,wPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(vPosRJac,wPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n2*n3)*VT + (vLoc-n2*VN)*(wLoc-n3*VN)/(VT+Eps) )

!   AS(IndexMet(wPosRJac,uPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosRJac,uPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n3*n1)*VT + (wLoc-n3*VN)*(uLoc-n1*VN)/(VT+Eps) )

!   AS(IndexMet(wPosRJac,vPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosRJac,vPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (   -n3*n2)*VT + (wLoc-n3*VN)*(vLoc-n2*VN)/(VT+Eps) )

!   AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1) = &
!   AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1) - &
!            RhoLoc*FL*DragM/(VolC(ix,iy,iz)+Eps) * ( (One-n3*n3)*VT + (wLoc-n3*VN)*(wLoc-n3*VN)/(VT+Eps) )

  END DO

END SUBROUTINE JacBaumDragCompute

SUBROUTINE TrafficDragCompute

! Abbremsung innerhalb der Traffic lines

  INTEGER :: i
  INTEGER :: ix,iy,iz

  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: vTr1,vTr2
  REAL(RealKind) :: AbsRelV
  REAL(RealKind) :: Frac,FL
  REAL(RealKind) :: vTraffic=30000.0d0/3600.0d0    ! [m/s]
  REAL(RealKind) :: DensTraffic=600.0d0/3600.0d0  ! [ /s]
  REAL(RealKind) :: cDTraffic=0.45d0    ! []

  DO i=1,NumberOfEmiDomainLoc
    DO iCell=1,SIZE(EmiDomain(i)%BlockCell(ibLoc)%Cell)
      ix=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%ix
      iy=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iy
      iz=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iz
      Frac=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%Frac
      FL=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%FL
      RhoLoc=Rho(ix,iy,iz,1)
      vTr1=vTraffic*EmiDomain(i)%n1
      vTr2=vTraffic*EmiDomain(i)%n2
      uLocL=uCL(ix,iy,iz,1)/(RhoLoc+Eps)
      uLocR=uCR(ix,iy,iz,1)/(RhoLoc+Eps)
      vLocL=vCL(ix,iy,iz,1)/(RhoLoc+Eps)
      vLocR=vCR(ix,iy,iz,1)/(RhoLoc+Eps)
      AbsRelV=SQRT((uLocL-vTr1)**2+(vLocL-vTr2)**2)
      uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)-RhoLoc*VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(uLocL-vTr1)*AbsRelV
      vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-RhoLoc*VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(vLocL-vTr2)*AbsRelV
      AbsRelV=SQRT((uLocR-vTr1)**2+(vLocR-vTr2)**2)
      uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-RhoLoc*VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(uLocR-vTr1)*AbsRelV
      vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-RhoLoc*VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(vLocR-vTr2)*AbsRelV
    END DO                         
  END DO   

END SUBROUTINE TrafficDragCompute

SUBROUTINE JacTrafficDragCompute

! Abbremsung innerhalb der Traffic lines

  INTEGER :: i
  INTEGER :: ix,iy,iz

  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: vTr1,vTr2
  REAL(RealKind) :: AbsRelV
  REAL(RealKind) :: Frac,FL
  REAL(RealKind) :: vTraffic=30000.0d0/3600.0d0    ! [m/s]
  REAL(RealKind) :: DensTraffic=600.0d0/3600.0d0  ! [ /s]
  REAL(RealKind) :: cDTraffic=0.45d0    ! []

  DO i=1,NumberOfEmiDomainLoc
    DO iCell=1,SIZE(EmiDomain(i)%BlockCell(ibLoc)%Cell)
      ix=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%ix
      iy=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iy
      iz=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iz
      Frac=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%Frac
      FL=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%FL
      RhoLoc=Rho(ix,iy,iz,1)
      vTr1=vTraffic*EmiDomain(i)%n1
      vTr2=vTraffic*EmiDomain(i)%n2
      uLocL=uCL(ix,iy,iz,1)/(RhoLoc+Eps)
      uLocR=uCR(ix,iy,iz,1)/(RhoLoc+Eps)
      vLocL=vCL(ix,iy,iz,1)/(RhoLoc+Eps)
      vLocR=vCR(ix,iy,iz,1)/(RhoLoc+Eps)
      AbsRelV=SQRT((uLocL-vTr1)**2+(vLocL-vTr2)**2)
!     uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)-VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(uLocL-vTr1)*AbsRelV
!     vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(vLocL-vTr2)*AbsRelV
      AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) = &
      AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)   &
        -0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*AbsRelV
      AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) = &
      AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1)   &
        -0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*AbsRelV
      AbsRelV=SQRT((uLocR-vTr1)**2+(vLocR-vTr2)**2)
!     uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(uLocR-vTr1)*AbsRelV
!     vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-VolC(ix,iy,iz)*0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*(vLocR-vTr2)*AbsRelV
      AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) = &
      AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)   &
        -0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*AbsRelV
      AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) = &
      AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1)   &
        -0.5d0*Frac*cDTraffic*DensTraffic/vTraffic*AbsRelV
    END DO                         
  END DO   

END SUBROUTINE JacTrafficDragCompute

!
!===========================================================
!----  Rain falling down in Earth
!===========================================================
!

SUBROUTINE RainSurfCompute(Vector,VelF,Rhs,Time)
  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoRLoc,NRLoc,RainFallQr,RainFallNr,RhoLoc
  REAL(RealKind) :: dLoc,VolLoc,VolDiff,Surface

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uC=>uCell(ibLoc)%c
  vC=>vCell(ibLoc)%c
  wC=>wCell(ibLoc)%c
  RhoR=>Vector%Vec(RhoRpos)%c
  Rho=>Vector%Vec(RhoPos)%c
  NRain=>Vector%Vec(nrPos)%c
  RR=>Vector%Vec(RhoRpos)%cB
  RhoRRhs=>Rhs%Vec(RhoRpos)%c
  NRRhs=>Rhs%Vec(nrPos)%c

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoRLoc= RhoR(ix,iy,iz,1)
    NRLoc  = NRain(ix,iy,iz,1)
    RhoLoc = Rho(ix,iy,iz,1)

    RainFallQr=MIN(fall_rain_q(RhoRLoc,NRLoc,RhoLoc)*RhoRLoc*FL*n3,Zero) ! =<0
    RainFallNr=MIN(fall_rain_n(RhoRLoc,NRLoc,RhoLoc)*NRLoc*FL*n3,Zero)   ! =<0

!   Substract rain from lowest cell
    RhoRRhs(ix,iy,iz,1)     = RhoRRhs(ix,iy,iz,1)    +RainFallQr 
    NRRhs(ix,iy,iz,1)       = NRRhs(ix,iy,iz,1)      +RainFallNr

  END DO

END SUBROUTINE RainSurfCompute

SUBROUTINE JacRainSurf(Vector,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoRLoc,NRLoc,RainFall,RhoLoc

  AS=>Jac
  RhoR=>Vector%Vec(RhoRpos)%c
  nRain=>Vector%Vec(nrPos)%c
  Rho=>Vector%Vec(RhoPos)%c

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoRLoc= RhoR(ix,iy,iz,1)
    NRLoc  = NRain(ix,iy,iz,1)
    RhoLoc = Rho(ix,iy,iz,1)  
      ! RhoR
      RainFall=fall_rain_q(RhoRLoc,NRLoc,RhoLoc)*n3*FL ! <0
      AS(IndexMet(RhoRposJac,RhoRposJac))%c(ix,iy,iz,1)=AS(IndexMet(RhoRposJac,RhoRposJac))%c(ix,iy,iz,1)+ &
                                                     RainFall &
                                                     /(VolC(ix,iy,iz) + Eps)
      ! NR
      RainFall=fall_rain_n(RhoRLoc,NRLoc,RhoLoc)*n3*FL ! <0
      AS(IndexMet(nrPosJac,nrPosJac))%c(ix,iy,iz,1)=AS(IndexMet(nrPosJac,nrPosJac))%c(ix,iy,iz,1)+ &
                                                     RainFall &
                                                     /(VolC(ix,iy,iz) + Eps)
  END DO

END SUBROUTINE JacRainSurf

SUBROUTINE IceSurfCompute(Vector,VelF,Rhs,Time)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoILoc,NILoc,IceFallRhoI,IceFallNi,RhoLoc
  REAL(RealKind) :: dLoc,VolLoc,VolDiff,Surface

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uC=>uCell(ibLoc)%c
  vC=>vCell(ibLoc)%c
  wC=>wCell(ibLoc)%c
  Rho=>Vector%Vec(RhoPos)%c
  RhoI=>Vector%Vec(RhoIPos)%c
  NIce=>Vector%Vec(niPos)%c
  IR=>Vector%Vec(RhoIPos)%cB !ice rate 
  RhoIRhs=>Rhs%Vec(RhoIPos)%c
  NIRhs=>Rhs%Vec(niPos)%c

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoILoc= RhoI(ix,iy,iz,1)
    RhoLoc = Rho(ix,iy,iz,1)
    NILoc  = NIce(ix,iy,iz,1)

    IceFallRhoI=MIN(fall_ice_q(RhoILoc,NILoc,RhoLoc)*RhoILoc*FL*n3,Zero) ! =<0
    IceFallNi=MIN(fall_ice_n(RhoILoc,NILoc,RhoLoc)*NILoc*FL*n3,Zero)   ! =<0

    !Substract ice from lowest cell
    RhoIRhs(ix,iy,iz,1)     = RhoIRhs(ix,iy,iz,1)    +IceFallRhoI 
    NIRhs(ix,iy,iz,1)       = NIRhs(ix,iy,iz,1)      +IceFallNi

  END DO

END SUBROUTINE IceSurfCompute

SUBROUTINE JacIceSurf(Vector,Jac)
  
  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)
  
  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoILoc,NILoc,IceFall,RhoLoc
  
  AS=>Jac
  RhoI=>Vector%Vec(RhoIPos)%c
  NIce=>Vector%Vec(niPos)%c
  Rho=>Vector%Vec(RhoPos)%c
  
  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoLoc = Rho(ix,iy,iz,1)
    RhoILoc= RhoI(ix,iy,iz,1)
    NILoc  = NIce(ix,iy,iz,1)
    ! RhoS
    IceFall=fall_ice_q(RhoILoc,NILoc,RhoLoc)*n3*FL ! <0
    AS(IndexMet(RhoIPosJac,RhoIPosJac))%c(ix,iy,iz,1)=AS(IndexMet(RhoIPosJac,RhoIPosJac))%c(ix,iy,iz,1)+ &
                                                       IceFall &
                                                       /(VolC(ix,iy,iz) + Eps)
    ! NI
    IceFall=fall_ice_n(RhoILoc,NILoc,RhoLoc)*n3*FL ! <0
    AS(IndexMet(niPosJac,niPosJac))%c(ix,iy,iz,1)=AS(IndexMet(niPosJac,niPosJac))%c(ix,iy,iz,1)+ &
                                                       IceFall &
                                                       /(VolC(ix,iy,iz) + Eps)
  END DO

END SUBROUTINE JacIceSurf

SUBROUTINE SnowSurfCompute(Vector,VelF,Rhs,Time)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs
  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoSLoc,NSLoc,SnowFallRhoS,SnowFallNs,RhoLoc
  REAL(RealKind) :: dLoc,VolLoc,VolDiff,Surface

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uC=>uCell(ibLoc)%c
  vC=>vCell(ibLoc)%c
  wC=>wCell(ibLoc)%c
  Rho=>Vector%Vec(RhoPos)%c
  RhoS=>Vector%Vec(RhoSPos)%c
  NSnow=>Vector%Vec(nsPos)%c
  SR=>Vector%Vec(RhoSPos)%cB !snow rate 
  RhoSRhs=>Rhs%Vec(RhoSPos)%c
  NSRhs=>Rhs%Vec(nsPos)%c

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoSLoc= RhoS(ix,iy,iz,1)
    RhoLoc = Rho(ix,iy,iz,1)
    NSLoc  = NSnow(ix,iy,iz,1)

    SnowFallRhoS=MIN(fall_snow_q(RhoSLoc,NSLoc,RhoLoc)*RhoSLoc*FL*n3,Zero) ! =<0
    SnowFallNs=MIN(fall_snow_n(RhoSLoc,NSLoc,RhoLoc)*NSLoc*FL*n3,Zero)   ! =<0

    !Substract snow from lowest cell
    RhoSRhs(ix,iy,iz,1)     = RhoSRhs(ix,iy,iz,1)    +SnowFallRhoS 
    NSRhs(ix,iy,iz,1)       = NSRhs(ix,iy,iz,1)      +SnowFallNs

  END DO

END SUBROUTINE SnowSurfCompute

SUBROUTINE JacSnowSurf(Vector,Jac)
  
  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)
  
  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: RhoSLoc,NSLoc,SnowFall,RhoLoc
  
  AS=>Jac
  RhoS=>Vector%Vec(RhoSPos)%c
  NSnow=>Vector%Vec(nsPos)%c
  Rho=>Vector%Vec(RhoPos)%c
  
  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    RhoLoc = Rho(ix,iy,iz,1)
    RhoSLoc= RhoS(ix,iy,iz,1)
    NSLoc  = NSnow(ix,iy,iz,1)
    ! RhoS
    SnowFall=fall_snow_q(RhoSLoc,NSLoc,RhoLoc)*n3*FL ! <0
    AS(IndexMet(RhoSPosJac,RhoSPosJac))%c(ix,iy,iz,1)=AS(IndexMet(RhoSPosJac,RhoSPosJac))%c(ix,iy,iz,1)+ &
                                                       SnowFall &
                                                       /(VolC(ix,iy,iz) + Eps)
    ! NI
    SnowFall=fall_snow_n(RhoSLoc,NSLoc,RhoLoc)*n3*FL ! <0
    AS(IndexMet(nsPosJac,nsPosJac))%c(ix,iy,iz,1)=AS(IndexMet(nsPosJac,nsPosJac))%c(ix,iy,iz,1)+ &
                                                       SnowFall &
                                                       /(VolC(ix,iy,iz) + Eps)
  END DO

END SUBROUTINE JacSnowSurf

FUNCTION SensibleHeatFlux(Time,LandClass,ix,iy,iz)
  REAL(RealKind) :: SensibleHeatFlux
  REAL(RealKind) :: Time
  INTEGER :: LandClass
  INTEGER :: ix,iy,iz

  REAL(RealKind) :: r,VolLocComp,VolLoc

! Input Parameter
  REAL(RealKind) :: FluxIxL,FluxIxR
  REAL(RealKind) :: FluxIyL,FluxIyR
  REAL(RealKind) :: FluxRandomLat

  SELECT CASE(FluxType) ! 'RICO', 'CONST', 'HAT', 'SINUSDRY', 'RANDOMDRY', 'LANDSEA', 'ORO', 'DYCOMS'
    CASE('CONST')
      SensibleHeatFlux=SensFluxConst
    CASE('DYCOMS')
      SensibleHeatFlux=16.0d0   ! 16 W/m**2 sensible heat flu
    CASE('BOMEX')
      SensibleHeatFlux=8.000d0  
!     SensibleHeatFlux=9.356d0  
    CASE('HAT')
      CALL Random_Number(r)
      IF(Time<=3600.0) THEN
        SensibleHeatFlux = (120.0 - 0.5d0 + 1.0 * r)
      ELSEIF(Time.gt.3600.0.and.Time<=9000.0) THEN
        SensibleHeatFlux = ((120.0d0-((Time-3600.0)/45.0)) -0.5d0+1.0d0*r)
      ELSEIF(Time.gt.9000.0.and.Time<=34200.0) THEN
        SensibleHeatFlux = ((0.0d0-((Time-9000.0)/1050.0))-0.5d0+1.0d0*r)
      ELSEIF(Time.gt.34200.0.and.Time<=55800.0) THEN
        SensibleHeatFlux = ((-24.0d0+((Time-34200.0)/900.0))-0.5d0+1.0d0*r)
      ELSEIF(Time.gt.55800.0) THEN
        SensibleHeatFlux = (300.0d0-0.5d0+1.0d0*r)*SIN((Time-55800.0)*Two*Pi/(24.0d0*3600.0d0))
      ELSE
        SensibleHeatFlux = (120.0 - 0.5d0 + 1.0 * r)
      ENDIF
    CASE('SINUSDRY')
      SensibleHeatFlux=200.0d0*SIN(Time*Two*Pi/(24.0d0*3600.0d0))
    CASE('RANDOMDRY')
      CALL Random_Number(r)
      SensibleHeatFlux=(120.0d0 - 0.5d0 + 1.0d0 * r) ! 120 W/m2 = 0.1 k m/s ;  Q(w/m2)  = cp *rho * Q(K m/s)
    CASE('LANDSEA')
      ! ----- Flux parameterization from Engelmann et al. 2011 (Tellus) -----
      IF ((ix>=FluxIxL.AND.ix<=FluxIxR).AND.(iy>=FluxIyL.AND.iy<=FluxIyR)) THEN ! land
        IF (MOD(Time,86400.0d0)>(SunRise+FluxShift).AND.MOD(Time,86400.0d0)<=(SunSet+FluxShift)) THEN 
          SensibleHeatFlux=SensFluxAmp*COS((Time-((SunSet-SunRise)/2.0d0)-(SunRise+FluxShift))/((SunSet-SunRise)/pi))+SensFluxMin
        ELSE 
          SensibleHeatFlux=-77.0d0
        END IF 
      ELSE ! sea
        SensibleHeatFlux=20.0d0 ! sensible
      END IF
    CASE('FIRE')
      IF ((ix.GE.Firex0.AND.ix.LE.Firex1).AND.(iy.GE.Firey0.AND.iy.LE.Firey1).AND.(Time.GT.FireStartTime)) THEN
          SensibleHeatFlux=FireSensFlux
      ELSE
          SensibleHeatFlux=SensFluxConst
      END IF
    CASE('BARBADOS')
      IF (LandClass.NE.9) THEN ! land
        IF (MOD(Time,86400.0d0)>(SunRise+FluxShift).AND.MOD(Time,86400.0d0)<=(SunSet+FluxShift)) THEN ! day
          SensibleHeatFlux=SensFluxAmp*COS((Time-((SunSet-SunRise)/2.0d0)-(SunRise+FluxShift))/((SunSet-SunRise)/pi))+SensFluxMin
        ELSE ! night
          SensibleHeatFlux=SensFluxMin
        END IF 
      ELSE ! sea
        SensibleHeatFlux=6.0d0
      END IF
    CASE('ORO')
      ! ----- Flux parameterization from Engelmann et al. 2011 (Tellus) with orography -----                    
      ! ------------------------------------------------------------------------------------
      ! ----- To differentiate between land and sea, the local cell volume is computed. ---- 
      ! ----- Use this for real topographic data.                                       ----
      ! ------------------------------------------------------------------------------------
      VolLocComp=dx(ix)*dy(iy)*dz(iz) ! max. volume
      VolLoc=VolC(ix,iy,iz)           ! cell volume
      IF (LandClass.NE.9) THEN ! land
        IF (MOD(Time,86400.0d0)>(SunRise+FluxShift).AND.MOD(Time,86400.0d0)<=(SunSet+FluxShift)) THEN 
          SensibleHeatFlux=SensFluxAmp*COS((Time-((SunSet-SunRise)/2.0d0)-(SunRise+FluxShift))/((SunSet-SunRise)/pi))+SensFluxMin
        ELSE 
          SensibleHeatFlux=-77.0d0
        END IF 
      ELSE ! sea
        SensibleHeatFlux=20.0d0
      END IF
    CASE DEFAULT  
      SensibleHeatFlux=Zero
  END SELECT    
END FUNCTION SensibleHeatFlux

FUNCTION LatentHeatFlux(Time,LandClass,ix,iy,iz)
  REAL(RealKind) :: LatentHeatFlux
  REAL(RealKind) :: Time
  INTEGER :: LandClass
  INTEGER :: ix,iy,iz

  REAL(RealKind) :: r,VolLocComp,VolLoc

! Input Parameter
  REAL(RealKind) :: FluxIxL,FluxIxR
  REAL(RealKind) :: FluxIyL,FluxIyR
  REAL(RealKind) :: FluxRandomLat

  SELECT CASE(FluxType) ! 'RICO', 'CONST', 'HAT', 'SINUSDRY', 'RANDOMDRY', 'LANDSEA', 'ORO', 'DYCOMS'
    CASE('CONST')
      LatentHeatFlux=LatFluxConst
    CASE('DYCOMS')
      CALL Random_Number(r)
      LatentHeatFlux=93.0d0+Two*FluxRandomLat*(r-0.5d0)   ! 93 +/-r W/m**2 latent   heat flux 
    CASE('BOMEX')
      LatentHeatFlux=150.00d0
!     LatentHeatFlux=147.62d0
    CASE('LANDSEA')
      ! ----- Flux parameterization from Engelmann et al. 2011 (Tellus) -----
      IF ((ix>=FluxIxL.AND.ix<=FluxIxR).AND.(iy>=FluxIyL.AND.iy<=FluxIyR)) THEN ! land
        LatentHeatFlux=55.0d0 ! constant latent heat flux
      ELSE ! sea
        CALL Random_Number(r)
        LatentHeatFlux=(90.0d0+FluxRandomLat*(r-0.5d0))   ! latent
      END IF  
    CASE('BARBADOS')
      IF (LandClass.NE.9) THEN ! land
        IF (MOD(Time,86400.0d0)>(SunRise+FluxShift).AND.MOD(Time,86400.0d0)<=(SunSet+FluxShift)) THEN ! day
          LatentHeatFlux=LatFluxAmp*COS((Time-((SunSet-SunRise)/2.0d0)-(SunRise+FluxShift))/((SunSet-SunRise)/pi))+LatFluxMin
        ELSE ! night
          LatentHeatFlux=LatFluxMin
        END IF 
      ELSE ! sea
        LatentHeatFlux=56.0d0
      END IF
    CASE('ORO')
      ! ----- Flux parameterization from Engelmann et al. 2011 (Tellus) with orography -----                    
      ! ------------------------------------------------------------------------------------
      ! ----- To differentiate between land and sea, the local cell volume is computed. ---- 
      ! ----- Use this for real topographic data.                                       ----
      ! ------------------------------------------------------------------------------------
      VolLocComp=dx(ix)*dy(iy)*dz(iz) ! max. volume
      VolLoc=VolC(ix,iy,iz)           ! cell volume
      IF (LandClass.NE.9) THEN ! land
        LatentHeatFlux=55.0d0 ! constant latent heat flux
      ELSE ! sea
        CALL Random_Number(r)
        LatentHeatFlux=90.0d0+Two*FluxRandomLat*(r-0.5d0) ! latent
      END IF
    CASE DEFAULT  
      LatentHeatFlux=Zero
  END SELECT    
END FUNCTION LatentHeatFlux


!
!===========================================================
!----  Fire-Routines for constant fluxes in BoundaryCells
!===========================================================
!

SUBROUTINE Fire(Vector,VelF,Rhs)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(VelocityFace_T) :: VelF
  TYPE(Vector4Cell_T) :: Rhs

  uF=>VelF%uF
  vF=>VelF%vF
  wF=>VelF%wF
  uC=>uCell(ibLoc)%c
  vC=>vCell(ibLoc)%c
  wC=>wCell(ibLoc)%c
  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  uRhsL=>Rhs%Vec(uPosL)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsR=>Rhs%Vec(wPosR)%c
  Rho=>Vector%Vec(RhoPos)%c
  Th=>Vector%Vec(ThPos)%c
  ThRhs=>Rhs%Vec(ThPos)%c
  RhoLRhs=>Rhs%Vec(9)%c

  CALL FireCompute

END SUBROUTINE Fire

SUBROUTINE FireCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN,ConstFlux
  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: DragM,DragH,ThetaS
  REAL(RealKind) :: RhoDLoc,RhoVLoc,RhoLLoc,p,T,pVs

  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1
    n2     = BoundCell(i)%n2
    n3     = BoundCell(i)%n3
    FL     = BoundCell(i)%FL
    DragH  = BoundCell(i)%DragH
    DragM  = BoundCell(i)%DragM
    ThetaS = BoundCell(i)%ThetaS
    RhoDLoc= Rho(ix,iy,iz,1)

    IF(ix>=9.and.iy>=29.and.ix<=12..and.iy<=31)THEN
     ConstFlux=100000.d0           ! 25 W/m**2 sensible heat flux
     ThRhs(ix,iy,iz,1) = ThRhs(ix,iy,iz,1)+FL*ConstFlux/(RhoDLoc*Cpd)
     RhoLRhs(ix,iy,iz,1) = 1.0d6
    END IF
   END DO
END SUBROUTINE FireCompute

!
!===========================================================
!---- Sun Angle Compute (cosine of solar zenith angle)
!===========================================================
!

SUBROUTINE SunAngleCompute(Time,phi)

!  REAL(RealKind) :: lng,lat,dekl,H,Time,cosPhi
  REAL(RealKind) :: dekl,H,Time,cosPhi
  REAL(RealKind),INTENT(IN),OPTIONAL :: phi
  REAL(RealKind) :: CurrentTime ! Current time within one day
  INTEGER :: CurrentDay ! Current day considering the passing time
!  lng = 12.0d0*Pi/180.0d0
!  lat = 51.0d0*Pi/180.0d0

  IF(SunMove)THEN
    CurrentTime = Time - INT(Time/86400.0d0)*86400.0d0
    CurrentDay = StartDay + INT(Time/86400.0d0)
    H = Pi*(CurrentTime/3600.0d0-12.0d0)/12.0d0-lng
  ELSE ! If the sun does not move, the time and day are always the start ones.
    CurrentTime = StartTime
    CurrentDay = StartDay
    H = Pi*(-12.0d0+CurrentTime/3600.0d0)/12.0d0-lng
  ENDIF

  IF (Sphere) THEN
    dekl  = COS(2.0d0*PI*(CurrentDay-173.0d0)/365.25d0)
    radn1 = COS(dekl)*COS(H)
    radn2 = -(COS(dekl)*SIN(H))
    radn3 = SIN(dekl)
    cosSun = COS(dekl)*COS(H)*COS(phi)+SIN(dekl)*SIN(phi)
  ELSE
    dekl  = 23.45d0*PI/180.0d0*COS(2.0d0*PI*(CurrentDay-173.0d0)/365.25d0)
    radn1 = -(COS(dekl)*SIN(H))                                ! Vektor der einfallenden Strahlung
    radn2 = -(COS(dekl)*COS(H)*SIN(lat)-SIN(dekl)*COS(lat))
    radn3 = COS(dekl)*COS(H)*COS(lat)+SIN(dekl)*SIN(lat)
    cosSun = radn3
  END IF

END SUBROUTINE SunAngleCompute

END MODULE Wall_Mod
