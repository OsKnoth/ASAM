MODULE Canopy_Mod
  USE Kind_Mod
  USE Parameter_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Chemie_Mod
  USE Megan_version_2
  USE RevWRF_Mod
  USE SoilData_Mod

  IMPLICIT NONE
  REAL(RealKind), PRIVATE, POINTER :: Th(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThB(:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThBRhs(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Wg(:,:)
  REAL(RealKind), PRIVATE, POINTER :: WgRhs(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoS(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoCloud(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoD(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Temp(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Pre(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Raddir(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Raddif(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Radinf(:,:)
  REAL(RealKind), PRIVATE, POINTER :: uF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCR(:,:,:,:)
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
  REAL(RealKind), PRIVATE, POINTER :: RhoDRhs(:,:,:,:)
  TYPE(Vec4_T), PRIVATE, POINTER :: AS(:)

  INTEGER, PRIVATE :: nzS  ! Number of possible soil levels 
  INTEGER, PRIVATE :: nzW  ! Number of possible soil levels 
  REAL(RealKind), ALLOCATABLE, PRIVATE :: dzSoil(:)
  REAL(RealKind), ALLOCATABLE, PRIVATE :: zSoil(:)
  
  INTEGER, PRIVATE, Pointer :: SoilType(:)
  INTEGER, PRIVATE :: LandClass  ! Land Use Class
  INTEGER, PRIVATE :: nrsl ! number of soil layers

  REAL(RealKind), PRIVATE, PARAMETER :: v_depos = 0.025d0 ! deposition velocity
CONTAINS

SUBROUTINE CanopyCompute(Vector,Rhs,Time)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vector4Cell_T) :: Rhs
  REAL(RealKind) :: Time

  !===== Velocity =====!
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
  !===== Density =====!
  Rho=>RhoCell(ibLoc)%c
  RhoV=>Vector%Vec(RhoVPos)%c
  RhoR=>Vector%Vec(RhoRpos)%c
  IF (RhoCPos>0) THEN  ! Barthel
    RhoCloud=>Vector%Vec(RhoCPos)%c
  ELSE
    RhoCloud=>RhoLCell(ibLoc)%c
  END IF
  RhoL=>RhoLCell(ibLoc)%c
  RhoI=>Vector%Vec(RhoIPos)%c
  RhoS=>Vector%Vec(RhoSPos)%c

  RhoRhs=>Rhs%Vec(rhoPos)%c
  RhoVRhs=>Rhs%Vec(RhoVPos)%c
  RhoLRhs=>Rhs%Vec(RhoCPos)%c
  !===== Potential temperature and soil related =====!
  Th=>Vector%Vec(thPos)%c
  ThB=>Vector%Vec(thPos)%cB
  Wg=>Vector%Vec(RhoVPos)%cB
  ! write(*,*) 'CanopyCompute: ThB.', ThB
  ! write(*,*) 'CanopyCompute: Wg.', Wg

  ThRhs=>Rhs%Vec(thPos)%c
  ThBRhs=>Rhs%Vec(thPos)%cB
  WgRhs=>Rhs%Vec(RhoVPos)%cB

  nrsl=Domain%nrsoillayers
  Temp=>TAbsCell(ibLoc)%Vec(1)%c
  Pre=>PreCell(ibLoc)%c

  CALL CanopyFluxCompute(Vector, Rhs, Time)
  CALL CanopySurfCompute_Soil(Time)
END SUBROUTINE CanopyCompute


SUBROUTINE JacCanopyCompute(MatrixMatVectorLU, Vector, Jac)
  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU
  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  Rho=>RhoCell(ibLoc)%c
  Th=>Vector%Vec(ThPos)%c
  AS=>Jac

  CALL JacCanopyFluxCompute(MatrixMatVectorLU)
  CALL JacCanopySurfCompute_Soil
END SUBROUTINE JacCanopyCompute


SUBROUTINE JacCanopyFluxCompute(MatrixMatVectorLU)
  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU
  INTEGER, POINTER :: DiagP(:), Permu(:)
  INTEGER :: diag

  INTEGER, PARAMETER :: EM_SPC_NUM=22
  CHARACTER(LEN=16), PARAMETER :: EM_SPC_NAME(EM_SPC_NUM) = (/ 'ISO             ', &
                                                               'MYRC            ', &
                                                               'SABI            ', &
                                                               'LIMN            ', &
                                                               'CAR3            ', &
                                                               'OCIM            ', &
                                                               'BPI             ', &
                                                               'API             ', &
                                                               'OMTP            ', &
                                                               'FARN            ', &
                                                               'BCAR            ', &
                                                               'OSQT            ', &
                                                               'MBO             ', &
                                                               'MET             ', &
                                                               'ACTO            ', &
                                                               'CH4             ', &
                                                               'NO              ', &
                                                               'ALD             ', &
                                                               'HCHO            ', &
                                                               'CO              ', &
                                                               'CIN             ', &
                                                               'LIN             ' /)

  INTEGER :: i,j,k,nrcl
  INTEGER :: ix,iy,iz
  REAL(RealKind), POINTER :: LAD(:)
  INTEGER, SAVE :: npos(EM_SPC_NUM+1) ! remember the position of each emitted gases in the vector
  LOGICAL, SAVE :: first_round = .TRUE. ! if the code is run the first time

  IF (CanopyEmission) THEN
    IF (first_round) then
      DO k = 1, EM_SPC_NUM
        npos(k) = Position(TRIM(EM_SPC_NAME(k)))
      END DO
      npos(EM_SPC_NUM+1) = Position('O3') ! position of ozone
      first_round = .FALSE.
    END IF

    DiagP => MatrixMatVectorLU%Struct%DiagPtr(:)
    Permu => MatrixMatVectorLU%Struct%Permu(:)

    DO i=1,NumBoundCell
      ix  =  BoundCell(i)%ix
      iy  =  BoundCell(i)%iy
      iz  =  BoundCell(i)%iz

      nrcl=  BoundCell(i)%CanopyCell%NrCanopyLayers
      LAD => BoundCell(i)%CanopyCell%LAD

      DO j=1,nrcl
        !========== Add Jacobian for deposition ==========!
        DO k=1,EM_SPC_NUM+1
          diag = DiagP(Permu(npos(k)))
          AS(diag)%c(ix,iy,iz+j-1,1) = AS(diag)%c(ix,iy,iz+j-1,1) - VolC(ix,iy,iz+j-1)*v_depos*LAD(j)
        END DO ! k
      END DO ! j
    END DO ! i
  END IF
END SUBROUTINE JacCanopyFluxCompute


SUBROUTINE JacCanopySurfCompute_Soil
  INTEGER :: i
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN  
  REAL(RealKind) :: FL,DragH,DragM
  REAL(RealKind) :: RhoAir
  
  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1G
    n2     = BoundCell(i)%n2G
    n3     = BoundCell(i)%n3G
    FL     = BoundCell(i)%FL
    DragH  = BoundCell(i)%DragH
    DragM  = BoundCell(i)%DragM

    RhoAir = Rho(ix,iy,iz,1)
!   Evaluate velocity tangential and normal  
    uLoc = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/RhoAir
    vLoc = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/RhoAir
    wLoc = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/RhoAir 
    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))

    AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n1*n1)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n2*n2)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n3*n3)/(VolC(ix,iy,iz)+Eps)

    AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n1*n1)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n2*n2)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n3*n3)/(VolC(ix,iy,iz)+Eps)
 
    AS(IndexMet(ThPosJac,ThPosJac))%c(ix,iy,iz,1) = AS(IndexMet(ThPosJac,ThPosJac))%c(ix,iy,iz,1)- &
                                                    FL*(DragH*VT+1.d-6)/(VolC(ix,iy,iz)+Eps)
  END DO
END SUBROUTINE JacCanopySurfCompute_Soil

!=======================================!
!----- Flux processes for canopy
!=======================================!
SUBROUTINE CanopyFluxCompute(Vector, Rhs, Time)
  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vector4Cell_T) :: Rhs

  REAL(RealKind)     :: Time
  REAL               :: Time_4

  INTEGER, PARAMETER :: EM_kzMax=20
  REAL,PARAMETER     :: ConvertionFactor=4.766       ! convert solar radiation to PAR, from [W m-2] to [umol m-2 s-1]
  INTEGER            :: EM_Julian                    ! Julian day of the year
  INTEGER            :: EM_Can_Lay                   ! Amount of layers inside the canopy
  REAL               :: EM_Beta                      ! [degree], solar zenith angle
  REAL               :: EM_Lat, EM_Long              ! Latitude and longitude
  INTEGER            :: EM_DATE                      ! Date for new Megan in format YYYDDD scalar
  INTEGER            :: EM_Time_M2                   ! Time for new Megan in hour minutes second format
  REAL               :: EM_PAR                       ! [umol m-2 s-1], Photosynthetically active radiation
  REAL               :: EM_PAR_day                   ! [umol m-2 s-1], Daily average photosynthetically active radiation
  REAL(RealKind)     :: TempK(EM_kzMax)              ! (kz), [K], Array of temperature in RealKind type
  REAL               :: EM_TempK(EM_kzMax)           ! (kz), [K], Array of temperature
  REAL               :: EM_TempK_day                 ! [K], Daily average of temperature
  REAL               :: EM_PRES                      ! [Pa], Pressure in the canopy
  REAL               :: EM_WVM(EM_kzMax)             ! (kz), [g kg-1], Array of water vapour mixing ratio
  REAL               :: EM_WIND(EM_kzMax)            ! (kz), [m s-1], Array of horizontal wind
  REAL               :: EM_SMOIST                    ! [%], Soil moisture
  REAL               :: EM_LAD(EM_kzMax)             ! (kz), [-], Array of vertical distribution for the LAI, values are in the range of [0,1]
  REAL               :: EM_LAI                       ! [m2 m-2], Leave area index
  REAL               :: EM_LAI_past                  ! [m2 m-2], Leave area index of last month
  REAL               :: EM_ER(20,EM_kzMax)           ! (20,kz), Array of output emission buffer
  REAL               :: EM_ER_HB(20,EM_kzMax)        ! (20,kz), Array of output emission buffer
  REAL               :: EM_ER_SB(20,EM_kzMax)        ! (20,kz), Array of output emission buffer
  REAL               :: EM_ER_NT(20,EM_kzMax)        ! (20,kz), Array of output emission buffer
  REAL               :: EM_ER_BT(20,EM_kzMax)        ! (20,kz), Array of output emission buffer
  CHARACTER(LEN=16)  :: EM_VAR(22)                   ! Name of the 22 VOC's
  REAL               :: EM_Ea1pL_M2(4,EM_kzMax)      ! (4,kz), Array of emission activity of light per layer
  REAL               :: EM_Ea1NL_M2(4,EM_kzMax)      ! (4,kz), Array of companied emission activity
  REAL               :: EM_Ea1tL_M2(4,EM_kzMax)      ! (4,kz), Array of emission activity of temperature per layer
  REAL               :: EM_GAM_TMP(20,EM_kzMax)      ! (20,kz), Array of temperature response factor
  REAL               :: EM_GAM_OTHER(20,3)           ! (20,3), Array of other gamma factors
  REAL               :: EM_z(EM_kzMax)               ! (kz), Array of height for each layer (m), e.g., (/1,3,5,7,9/)
  REAL               :: EM_dz(EM_kzMax)              ! (kz), Array of vertical size for each layer (m), e.g., (/1,2,2,2,2/)
  INTEGER            :: EM_kz                        ! Amount of layers in the model run (not only canopy)
  REAL               :: EM_EMI(EM_kzMax,22)          ! (kz,22), [# cm-3 s-1], Matrix of emissions for each layer and each compound
  REAL               :: EM_Sun_Par(EM_kzMax)         ! (kz), Array of sun fraction - 1 above the canopy and decreases inside the canopy
  REAL               :: EM_SunleafTK(EM_kzMax)       ! (kz), Array of temparture for sun leaf in (K) and (C)
  REAL               :: EM_ShadeleafTK(EM_kzMax)     ! (kz), Array of temparture for shade leaf in (K) and (C)
  REAL               :: EM_Sunfrac(EM_kzMax)         ! (kz), Array of the fraction of sun leaves. NB: i = 1 is the top canopy layer
  REAL               :: EM_SH(EM_kzMax)              ! (can_lay), Array of vertical distribution for the sensble heat flux from the canopy
  REAL               :: EM_LH(EM_kzMax)              ! (can_lay), Array of vertical distribution for the latend heat flux from the canopy
  REAL               :: EM_Rbs, EM_Rds, EM_RLs       ! [W m-2], direct, diffuse and longwave radiation on the surface
  INTEGER, PARAMETER :: EM_SPC_NUM=22
  CHARACTER(LEN=16), PARAMETER :: EM_SPC_NAME(EM_SPC_NUM) = (/ 'ISO             ', &
                                                               'MYRC            ', &
                                                               'SABI            ', &
                                                               'LIMN            ', &
                                                               'CAR3            ', &
                                                               'OCIM            ', &
                                                               'BPI             ', &
                                                               'API             ', &
                                                               'OMTP            ', &
                                                               'FARN            ', &
                                                               'BCAR            ', &
                                                               'OSQT            ', &
                                                               'MBO             ', &
                                                               'MET             ', &
                                                               'ACTO            ', &
                                                               'CH4             ', &
                                                               'NO              ', &
                                                               'ALD             ', &
                                                               'HCHO            ', &
                                                               'CO              ', &
                                                               'CIN             ', &
                                                               'LIN             ' /)

  INTEGER :: i,j,k,jj,nrcl
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL
  REAL(RealKind) :: z0
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind) :: LAI
  REAL(RealKind) :: Albc ! Albedo of the canopy
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN,ConstFlux,MoistFlux,TotalFlux
  REAL(RealKind) :: ThLoc,RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc,pLoc,T
  REAL(RealKind) :: VolLoc,SurLoc
  REAL(RealKind) :: raddirekt,raddiffus,radinfred,raddirdif
  REAL(RealKind) :: Rm,Cpml,LvLoc
  REAL(RealKind) :: drvdt,rrv,rrl,dPotdt

  INTEGER, SAVE :: npos(EM_SPC_NUM) ! remember the position of each emitted gases in the vector
  LOGICAL, SAVE :: first_round = .TRUE. ! if the code is run the first time


  EM_Julian=INT(StartDay+Time/86400.0d0)
  EM_Beta=ACOS(cosSun)*180.0/PI
  EM_Lat=lat*180.0/PI
  EM_Long=lng*180.0/PI
  EM_DATE=1000*2010+EM_Julian
  Time_4=Time
  EM_Time_M2=TimeFormatter( Time_4-INT(Time/86400.0d0)*86400.0 )
  EM_PAR_day=500
  EM_TempK_day=278.70 ! Daily temperature in Jan-Dec in Hyytiälä:
                      ! 265.4500, 265.1500, 269.6500, 275.0500, 282.3500,
                      ! 287.4500, 289.4500, 287.3500, 281.8500, 276.8500,
                      ! 271.5500, 267.4500

  IF (first_round) then
    DO k = 1, EM_SPC_NUM
      npos(k) = Position(TRIM(EM_SPC_NAME(k)))
    END DO
    first_round = .FALSE.
  END IF

  DO i=1,NumBoundCell
    ix  =  BoundCell(i)%ix
    iy  =  BoundCell(i)%iy
    iz  =  BoundCell(i)%iz
    n1  =  BoundCell(i)%n1
    n2  =  BoundCell(i)%n2
    n3  =  BoundCell(i)%n3
    FL  =  BoundCell(i)%FL ! surface area of the cut surface on the ground
    nrcl=  BoundCell(i)%CanopyCell%NrCanopyLayers
    LAD => BoundCell(i)%CanopyCell%LAD
    LAI =  BoundCell(i)%CanopyCell%LAI
    Albc = 0.0d0 ! BoundCell(i)%CanopyCell%Alb
    raddirekt = BoundCell(i)%raddirekt
    raddiffus = BoundCell(i)%raddiffus
    radinfred = BoundCell(i)%radinfred
    raddirdif = raddirekt+raddiffus
    z0=zH(ix,iy) ! elevation of orography

    EM_Can_Lay=nrcl ! number of layers inside the canopy
    EM_kz=EM_Can_Lay+1 ! number of layers totally, usually equals to EM_Can_Lay+1

    EM_PAR=raddirdif*(1.0d0-Albc)*ConvertionFactor*0.5 ! Solar*ConvertionFactor/2
    ! write(*,*) 'em_par: ', EM_PAR
    EM_PRES=0.0
    DO j=1,EM_kz
      ThLoc   = Th(ix,iy,iz+j-1,1)
      RhoLoc  = Rho(ix,iy,iz+j-1,1)
      RhoVLoc = RhoV(ix,iy,iz+j-1,1)
      RhoLLoc = RhoCloud(ix,iy,iz+j-1,1) + RhoR(ix,iy,iz+j-1,1)
      RhoDLoc = RhoLoc-RhoVLoc-RhoLLoc
      pLoc    = Pre(ix,iy,iz+j-1,1) ! PressureTheta(RhoDLoc,RhoVLoc,RhoLLoc,ThLoc)+Eps

      TempK(j) = Temp(ix,iy,iz+j-1,1) ! AbsTemp(RhoDLoc,RhoVLoc,pLoc) ! RealKind type temperature
      EM_TempK(j) = TempK(j) ! Real type temperature
      IF (j<=EM_Can_Lay) THEN
        EM_PRES = EM_PRES+pLoc
      END IF
    END DO
    EM_PRES=EM_PRES/EM_Can_Lay ! Use vertical average value of pressure inside the canopy
    ! write(*,*) 'Canopy_Mod: EM_PRES, ', EM_PRES
    EM_WVM(1:EM_kz)=MAX(RhoV(ix,iy,iz:iz+EM_kz-1,1) / &
      (Rho(ix,iy,iz:iz+EM_kz-1,1) - RhoV(ix,iy,iz:iz+EM_kz-1,1) - RhoL(ix,iy,iz:iz+EM_kz-1,1) + Eps)*1000.0, 0.0) ! [g kg-1]
    EM_WIND(1:EM_kz)=SQRT(uC(ix,iy,iz:iz+EM_kz-1,1)*uC(ix,iy,iz:iz+EM_kz-1,1)+ &
                          vC(ix,iy,iz:iz+EM_kz-1,1)*vC(ix,iy,iz:iz+EM_kz-1,1)+ &
                          wC(ix,iy,iz:iz+EM_kz-1,1)*wC(ix,iy,iz:iz+EM_kz-1,1)) ! horizontal or total wind?
    EM_SMOIST=Wg(i,1) ! 0.2
    ! write(*,*) 'EM_SMOIST.', EM_SMOIST
    EM_z(1:EM_kz)=zp(iz:iz+EM_kz-1)-z0
    ! write(*,*) 'EM_z, z0: ', EM_z, z0
    EM_dz(2:EM_kz)=EM_z(2:EM_kz)-EM_z(1:EM_kz-1)
    EM_dz(1)=EM_z(1)
    EM_LAI=LAI
    EM_LAI_past=EM_LAI
    EM_LAD(1:EM_Can_Lay)=LAD*EM_dz(1:EM_Can_Lay)/EM_LAI
    EM_LAD(EM_kz)=0.0
    ! write(*,*) 'EM_LAD', EM_LAD
    ! write(*,*) EM_z, EM_dz

    EM_SH=0.0
    EM_LH=0.0

    ! write(*,*) 'Canopy_Mod: EM_TempK.', EM_TempK
    CALL EMISSION_M2(EM_Julian, EM_Can_Lay, EM_Beta, EM_Lat, EM_Long, EM_DATE, EM_Time_M2, EM_PAR, EM_PAR_day, &
                     EM_TempK(1:EM_kz), EM_TempK_day, EM_PRES, EM_WVM(1:EM_kz), EM_WIND(1:EM_kz), EM_SMOIST, &
                     EM_LAD(1:EM_kz), EM_LAI, EM_LAI_past, &
                     EM_ER(:,1:EM_kz), EM_VAR, EM_ER_HB(:,1:EM_kz), EM_ER_SB(:,1:EM_kz), EM_ER_NT(:,1:EM_kz), EM_ER_BT(:,1:EM_kz), &
                     EM_Ea1pL_M2(:,1:EM_kz), EM_Ea1NL_M2(:,1:EM_kz), EM_Ea1tL_M2(:,1:EM_kz), &
                     EM_GAM_TMP(:,1:EM_kz), EM_GAM_OTHER, EM_z(1:EM_kz), EM_kz, EM_EMI(1:EM_kz,:), &
                     EM_Sun_Par(1:EM_kz), EM_SunleafTK(1:EM_kz), EM_ShadeleafTK(1:EM_kz), EM_Sunfrac(1:EM_kz), &
                     EM_SH(1:EM_kz), EM_LH(1:EM_kz), &
                     EM_Rbs, EM_Rds, EM_RLs)

    BoundCell(i)%CanopyCell%RadPtrRate=EM_Sunfrac(1)
    ! write(*,*) 'Canopy_Mod: EM_Sunfrac: ', EM_Sunfrac(1:EM_Can_Lay+1)
    ! BoundCell(i)%CanopyCell%AvgTem=0.5d0*SUM(EM_SunleafTK(1:EM_Can_Lay)+EM_ShadeleafTK(1:EM_Can_Lay))/EM_Can_Lay ! mean of sunleaf and shadeleaf
    BoundCell(i)%CanopyCell%AvgTem=EM_ShadeleafTK(1) ! SUM(EM_ShadeleafTK(1:EM_Can_Lay))/EM_Can_Lay ! only from shadeleaf
    !===== Set the leaf temperature
    BoundCell(i)%CanopyCell%SunlitLeafT=EM_SunleafTK(1:EM_Can_Lay) ! sunlit leaf
    BoundCell(i)%CanopyCell%ShadedLeafT=EM_ShadeleafTK(1:EM_Can_Lay) ! shaded leaf
    BoundCell(i)%CanopyCell%Hc=EM_SH(1:EM_Can_Lay) ! sunlit leaf
    BoundCell(i)%CanopyCell%LEc=EM_LH(1:EM_Can_Lay) ! shaded leaf
    BoundCell(i)%CanopyCell%Rbs=EM_Rbs
    BoundCell(i)%CanopyCell%Rds=EM_Rds
    BoundCell(i)%CanopyCell%RLs=EM_RLs

    !=====================================================!
    ! Debug
    !=====================================================!
    ! write(*,*) 'Canopy_Mod: EM_Rbs, EM_Rds, EM_RLs.', EM_Rbs, EM_Rds, EM_RLs
    ! write(*,*) 'canopy_mod, EM_TempK: ', EM_TempK(1:EM_kz)
    ! write(*,*) 'canopy_mod, sunleafTK: ', EM_SunleafTK(1:EM_kz-1)
    ! write(*,*) 'canopy_mod, shadeleafTK: ', EM_ShadeleafTK(1:EM_kz-1)
    ! write(*,*) 'canopy_mod, AvgTem: ', BoundCell(i)%CanopyCell%AvgTem
    ! write(*,*) 'canopy_mod, em_sunfrac: ', EM_Sunfrac(1:10)

    ! write(*,*) 'Top: ', raddirekt, raddiffus, radinfred
    ! write(*,*) 'Bot: ', Raddir(i,1), Raddif(i,1), Radinf(i,1)
    ! write(*,*) 'sunleafTK: ', EM_SunleafTK(1:5)
    ! write(*,*) 'shadeleafTK: ', EM_ShadeleafTK(1:5)
    ! write(*,*) 'gapfunc: ', BoundCell(i)%CanopyCell%GapFunc
    ! write(*,*) 'AvgTem: ', BoundCell(i)%CanopyCell%AvgTem
    ! write(*,*) 'EM_SH: ', EM_SH(1:EM_Can_Lay)
    ! write(*,*) 'EM_LH: ', EM_LH(1:EM_Can_Lay)
    ! write(*,*) 'EM_WVM: ', EM_WVM(1:8)

    DO j=1,EM_Can_Lay
      !========== Add heat fluxes to energy equation ==========!
      ThLoc   = Th(ix,iy,iz+j-1,1)
      RhoLoc  = Rho(ix,iy,iz+j-1,1)
      RhoVLoc = RhoV(ix,iy,iz+j-1,1)
      RhoLLoc = RhoCloud(ix,iy,iz+j-1,1) + RhoR(ix,iy,iz+j-1,1)
      RhoDLoc = RhoLoc-RhoVLoc-RhoLLoc
      pLoc    = Pre(ix,iy,iz+j-1,1) ! PressureTheta(RhoDLoc,RhoVLoc,RhoLLoc,ThLoc)+Eps
      T       = Temp(ix,iy,iz+j-1,1) ! AbsTemp(RhoDLoc,RhoVLoc,pLoc)
      Rm      = Rd+Rv*RhoVLoc/RhoDLoc+Eps
      Cpml    = (Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Eps)/RhoDLoc
      LvLoc   = LatHeat(TempK(j)) ! [J kg-1]
      ! write(*,*) 'Canopy_Mod: VolC.', VolC(ix,iy,iz+j-1)
      SurLoc  = VolC(ix,iy,iz+j-1)/dz(iz+j-1) ! EM_SH = SH *LAD*dz/LAI *LAI [W m-2], SH_flux = EM_SH /dz*VolC = EM_SH *SurLoc [W]

      TotalFlux=ThLoc*( EM_SH(j)*SurLoc/(RhoDLoc*Cpml*T) + &
                EM_LH(j)*SurLoc/(RhoDLoc*LvLoc)*(Rv/Rm-Rm/Cpml*log(pLoc/p0)*(Rv/Rm-Cpv/Cpml)) ) ! [kg K s-1]
      MoistFlux=+EM_LH(j)*SurLoc/LvLoc ! [kg s-1]

      ! drvdt=MoistFlux/RhoDLoc
      ! write(*,*) 'TotalFlux: ', TotalFlux
      ! write(*,*) 'MoistFlux: ', MoistFlux
      ThRhs(ix,iy,iz+j-1,1)=ThRhs(ix,iy,iz+j-1,1)+TotalFlux ! [kg K s-1]
      RhoVRhs(ix,iy,iz+j-1,1)=RhoVRhs(ix,iy,iz+j-1,1)+MoistFlux ! [kg s-1]
      RhoRhs(ix,iy,iz+j-1,1)=RhoRhs(ix,iy,iz+j-1,1)+MoistFlux ! [kg s-1]
      !========== Add emissions ==========!
      IF (CanopyEmission) THEN
        DO k=1,EM_SPC_NUM
          Rhs%Vec(npos(k))%c(ix,iy,iz+j-1,1) = Rhs%Vec(npos(k))%c(ix,iy,iz+j-1,1) &
            + VolC(ix,iy,iz+j-1) * ( EM_EMI(iz+j-1,k)*VolC(ix,iy,iz+j-1)*1.0d6 / AVGA & ! Emission, [mol s-1]
            - v_depos * LAD(j) * Vector%Vec(npos(k))%c(ix,iy,iz+j-1,1) ) ! Deposition, [mol s-1]
!         IF (TRIM(EM_SPC_NAME(k)) == 'API') THEN
!           write(*,*) EM_SPC_NAME(k), ' ', j, EM_EMI(j,k)
!           write(*,*) EM_SPC_NAME(k), ' ', j, Vector(ib)%Vec(npos(k))%c(ix,iy,iz+j-1,1)
!         END IF
        END DO ! k
      END IF
    END DO ! j

    !========== Set ground surface temperature ==========!
    ! write(*,*) 'CFC: tes.', BoundCell(i)%TeS ! or other function of time
  END DO ! i: NumBoundCell
END SUBROUTINE CanopyFluxCompute


!=======================================!
!----- Surface processes for canopy
!=======================================!
SUBROUTINE CanopySurfCompute_Soil(Time)

  INTEGER :: i,iw
  INTEGER :: ix,iy,iz,shad
  INTEGER :: nzThB,nzWg
  REAL(RealKind) :: Time
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: n1G,n2G,n3G
  REAL(RealKind) :: FL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: u_shear                ! Shear Velocity
  REAL(RealKind) :: F_SH,F_LH              ! Sensible / Latent Heat Flux
  REAL(RealKind) :: TermW,TermT
  REAL(RealKind) :: zH,zRauh,zRauhT
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: Qdir,Qdif,Qinf,raddir,raddif,radinf
  REAL(RealKind) :: Qlat,Qsens
  REAL(RealKind) :: AlbSoil,AlbVeg
  REAL(RealKind) :: ThAir,TAir
  REAL(RealKind) :: RhoVAir,RhoLAir,RhoIAir,RhoDAir,RhoAir
  REAL(RealKind) :: RhoV_Sat               ! Saturation Humidity at Ground 
  REAL(RealKind) :: RhoV_g                 ! Specific Humdity at Interception Reservoir
  REAL(RealKind) :: DragM,DragH,ThetaS
  REAL(RealKind) :: Wgn                    ! Interception Water [m]
  REAL(RealKind) :: RhsWgn
  REAL(RealKind) :: WSroot                 ! Root Water Content
  REAL(RealKind) :: TS(nzS),WS(nzW)        ! Soil Temperature, Water Content [m**3/m**3]
  REAL(RealKind) :: RhsTS(nzS),RhsWS(nzW)
  REAL(RealKind) :: FSoil                  ! Surface Radiative Flux 
  REAL(RealKind) :: LAI                    ! Leaf Area Index
  REAL(RealKind) :: F_plant                ! Area Covered with Plants
  REAL(RealKind) :: F_i                    ! Area Covered with intercepted water
  REAL(RealKind) :: F_par                  ! Fraction of absorbed photosynthetically active radiation
  REAL(RealKind) :: Cpml,Cvml,Lv,Rm,p
  REAL(RealKind) :: hu                     ! Relative Humdity near Ground
  REAL(RealKind) :: Eva_Pot,Eva_Soil,Eva_Intercept 
  REAL(RealKind) :: Transpiration(nzW),SumTransp
  REAL(RealKind) :: I_perc                 ! Percolation Rate of intercepted Water into Ground
  REAL(RealKind) :: dpotdt,rrl,ploc,potm,rrv,drvdt,ThLoc
  REAL(RealKind) :: TotalFlux,MoistFlux,SensFlux,Tsurface,dThdt
  REAL(RealKind) :: dLoc,Surface,VolLoc,VolDiff
  REAL(RealKind) :: r

  REAL(RealKind) :: Albc


  DO i=1,NumBoundCell
    Wgn        = 0.0d0
    RhsWS      = 0.0d0
    RhsTS      = 0.0d0
    RhsWgn     = 0.0d0
    WSroot     = 0.0d0
    SumTransp  = 0.0d0

    ! write(*,*) 'CSC_Soil: 1'
    ix         = BoundCell(i)%ix
    iy         = BoundCell(i)%iy
    iz         = BoundCell(i)%iz
!    n1         = BoundCell(i)%n1
!    n2         = BoundCell(i)%n2
!    n3         = BoundCell(i)%n3
    n1         = BoundCell(i)%n1G
    n2         = BoundCell(i)%n2G
    n3         = BoundCell(i)%n3G
    n1G        = BoundCell(i)%n1G
    n2G        = BoundCell(i)%n2G
    n3G        = BoundCell(i)%n3G
    FL         = BoundCell(i)%FL
    alb        = BoundCell(i)%alb
    ee         = BoundCell(i)%ee
    shad       = BoundCell(i)%shad
    DragH      = BoundCell(i)%DragH
    DragM      = BoundCell(i)%DragM
    LandClass  = BoundCell(i)%LandClass
    zRauh      = BoundCell(i)%zRauh
    zRauhT     = BoundCell(i)%zRauhT
    Albc       = 0.0d0  ! BoundCell(i)%CanopyCell%Alb ! albedo for the canopy
    RhoAir     = Rho(ix,iy,iz,1)
    RhoVAir    = RhoV(ix,iy,iz,1)
    RhoLAir    = RhoCloud(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
    RhoIAir    = RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
    ThLoc      = Th(ix,iy,iz,1)/(RhoAir+Eps)
    ThetaS     = BoundCell(i)%ThetaS
    SoilType   =>BoundCell(i)%SoilType
    ThAir      = Th(ix,iy,iz,1)
    nzThB      = SIZE(ThB(i,:))
    TS(1:nzThB)= ThB(i,1:)
    nzWg       = SIZE(Wg(i,:))-1
    Wgn        = MIN(MAX(0.0d0,Wg(i,nzWg+1)),5.0d-4*(1.0d0+5.0d0*PlantCover(LandClass)))
    WS(1:nzWg) = MAX(Zero,MIN(Wg(i,1:nzWg),cporv(SoilType(1:nzWg))))
    ! write(*,*) 'CSC_Soil: 2'

!   Evaluate velocity tangential and normal   
    uLocL = uCL(ix,iy,iz,1)/(RhoAir+Eps)
    uLocR = uCR(ix,iy,iz,1)/(RhoAir+Eps)
    uLoc  = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoAir+Eps)
    vLocL = vCL(ix,iy,iz,1)/(RhoAir+Eps)
    vLocR = vCR(ix,iy,iz,1)/(RhoAir+Eps)
    vLoc  = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoAir+Eps)
    wLocL = wCL(ix,iy,iz,1)/(RhoAir+Eps)
    wLocR = wCR(ix,iy,iz,1)/(RhoAir+Eps)
    wLoc  = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoAir+Eps)
    VN    = uLoc*n1+vLoc*n2+wLoc*n3
    V     = uLoc*uLoc+ &
            vLoc*vLoc+ &
            wLoc*wLoc
    VT    = SQRT(MAX(V-VN*VN,Zero))          

    ! write(*,*) 'CSC_Soil: 3'
!   Wind
    uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(uLocL-n1*VN)
    vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)- &
                 RhoAir*FL*DragM*VT*(vLocL-n2*VN)
    wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(wLocL-n3*VN)
    uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(uLocR-n1*VN)
    vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(vLocR-n2*VN)
    wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(wLocR-n3*VN)
    
!   Evaluate height z
    zH = Half*dz(iz) 
    ! write(*,*) 'Canopy_Mod-CanopySurfCompute: TeS.', BoundCell(i)%TeS

    IF (LandClass==9) THEN ! ocean
      TSurface = BoundCell(i)%TeS 
      RhoV_Sat = SaturVapor(TSurface)/(Rv*TSurface)   
      RhoDAir  = RhoAir-RhoVAir-RhoLAir+Eps
      p        = PressureTheta(RhoDAir,RhoVAir,RhoLAir,RhoIAir,ThAir)+Eps
      TAir     = AbsTemp(RhoDAir,RhoVAir,p)
      Cpml     = Cpd*RhoDAir+Cpv*RhoVAir+Cpl*RhoLAir+Eps
      SensFlux = -FL*(DragH*VT+1.d-6)*(TAir-TSurface)                 ! [K m3 s-1]
      MoistFlux= -FL*(DragH*VT+1.d-6)*(RhoVAir-RhoV_Sat)/(RhoAir+Eps) ! [  m3 s-1]
      ! Add random noise on latent heat flux
      CALL Random_Number(r)
      MoistFlux=(MoistFlux*LatHeat(TAir)*RhoAir/FL)-0.5d0+r           ! [kg s-3 = W m-2]
      SensFlux =(SensFlux*Cpml/FL)                                    ! [kg s-3 = W m-2] 
      BoundCell(i)%FluxLat=MoistFlux  ! [W m-2]        
      BoundCell(i)%FluxSens=SensFlux  ! [W m-2]                       
      drvdt    = 1.0d0/RhoDAir*MoistFlux
      rrv      = RhoVAir/RhoDAir
      rrl      = RhoLAir/RhoDAir
      Rm       = Rd*RhoDAir+Rv*RhoVAir+Eps
      Qsens    = SensFlux/Cpml*FL                  ! [K m3 s-1]
      Qlat     = MoistFlux/LatHeat(TAir)/RhoAir*FL ! [  m3 s-1]
      dThdt    = ThLoc*Qlat*(Rv/Rm-LOG((p/p0)**(Rm/Cpml))*(Rv/Rm-Cpv/Cpml)) & ! latent heat flux
                +ThLoc*Qsens/TSurface  ! sensible heat flux 
      TotalFlux= dThdt

      ThRhs(ix,iy,iz,1)=ThRhs(ix,iy,iz,1)+TotalFlux
      RhoVRhs(ix,iy,iz,1)=RhoVRhs(ix,iy,iz,1)+Qlat                ! lat. HeatFlux 
      RhoRhs(ix,iy,iz,1)=RhoRhs(ix,iy,iz,1)+Qlat                  ! lat. HeatFlux

    ELSE  ! land
      ! write(*,*) 'CSC: nrsl.', Domain%nrSoilLayers
      IF (nrsl>1) THEN
!       Soil water flux
        DO iw=1,nzWg-1
          IF (iw.LE.(nzWg-2)) THEN
            TermW       = FlussCond(WS(iw),WS(iw+1),dzSoil(iw),dzSoil(iw+1),SoilType(iw),SoilType(iw+1))
            RhsWS(iw)   = RhsWS(iw)-TermW
            RhsWS(iw+1) = RhsWS(iw+1)+TermW
          ELSE
            TermW       = KappaSoil(WS(iw),SoilType(iw)) ! only gravitational drainage below 5th layer 
            RhsWS(iw)   = RhsWS(iw)-TermW
            RhsWS(iw+1) = RhsWS(iw+1)+TermW
          END IF
        END DO
  
!       Soil temperature flux
        DO iw=1,nzThB-1
          TermT       = FlussDiffT(TS(iw),TS(iw+1),dzSoil(iw),dzSoil(iw+1),WS(iw),WS(iw+1),SoilType(iw),SoilType(iw+1))
          RhsTS(iw)   = RhsTS(iw)-TermT
          RhsTS(iw+1) = RhsTS(iw+1)+TermT
        END DO
!       Lower Boundary Condition
        RhsWS(nzWg) = RhsWS(nzWg)-KappaSoil(WS(nzWg),SoilType(nzWg))
!       Lower Boundary Condition
        TermT         = FlussDiffT(TS(nzThB),297.0_RealKind,Half*dzSoil(nzThB),Half*dzSoil(nzThB), &
                                   WS(nzThB),WS(nzThB),SoilType(nzThB),SoilType(nzThB))
        RhsTS(nzThB)  = RhsTS(nzThB)-TermT
        ! write(*,*) 'CSC_Soil: nrsl>1'
      END IF ! nr_SoilLayers>1
      ! write(*,*) 'CSC_Soil: 5', nrsl

!     Lower Boundary Condition
      RhsWS(nzWg) = RhsWS(nzWg)-KappaSoil(WS(nzWg),SoilType(nzWg))
      ! write(*,*) 'CSC_Soil: 5.1'
!     Lower Boundary Condition
      ! write(*,*) 'TS, TSClim, dzSoil, WS, SoilType'
      ! write(*,*) 'TS', TS(1)
      ! write(*,*) 'TSClim', TSClim
      ! write(*,*) 'WS', WS
      ! write(*,*) 'SoilType', SoilType
      TermT         = FlussDiffT(TS(nzThB),TSClim,Half*dzSoil(nzThB),Half*dzSoil(nzThB), &
                                 WS(nzThB),WS(nzThB),SoilType(nzThB),SoilType(nzThB))
      ! write(*,*) 'CSC_Soil: 5.15'
      RhsTS(nzThB)  = RhsTS(nzThB)-TermT

      ! write(*,*) 'CSC_Soil: 5.2'
!     ==== Init Data =====    
      RhoDAir    = RhoAir-RhoVAir-RhoLAir+Eps
      Cpml       = Cpd*RhoDAir+Cpv*RhoVAir+Cpl*RhoLAir
      Cvml       = Cvd*RhoDAir+Cvv*RhoVAir+Cpv*RhoLAir
      Rm         = Rd*RhoDAir+Rv*RhoVAir
      p          = PressureTheta(RhoDAir,RhoVAir,RhoLAir,RhoIAir,ThAir)
      TAir       = AbsTemp(RhoDAir,RhoVAir,p)
      RhoV_Sat   = SaturVapor(TS(1))/(Rv*TS(1))    ! New one
      RhoV_g     = Wgn*RhoV_Sat+(1.d0-Wgn)*RhoVAir
      Lv         = LatHeat(TAir)
      LAI        = LeafAreaIndex(LandClass)
      F_plant    = PlantCover(LandClass)
      F_i        = (1.0d0-F_plant)*(1.0d0-EXP(-Wgn/1.d-3))
      F_par      = One-EXP(-0.6d0*LAI)
      u_shear    = MAX(0.001d0,((DragM*VT)**Two*(uLoc**Two+vLoc**Two))**(One/Four))
      AlbSoil    = (One-F_plant)*AlbedoSoil(WS(1),SoilType(1))
      AlbVeg     = F_plant*( &
                   (One-F_par)*AlbedoSoil(WS(1),SoilType(1)) &
                   +F_par*AlbedoVeg(F_plant,LandClass))
      BoundCell(i)%CanopyCell%AlbSoil = AlbSoil + AlbVeg
      ! write(*,*) 'AlbSoil, AlbVeg.', AlbSoil, AlbVeg
      ! write(*,*) 'Albc.', Albc
      ! write(*,*) 'CSC_Soil: 5.3'
      ! write(*,*) 'LandClass, PlantCover: ', LandClass, PlantCover(0)
      ! write(*,*) 'LAI, F_plant, F_par, AlbSoil, AlbVeg: ', LAI, F_plant, F_par, AlbSoil, AlbVeg

      ! raddir=BoundCell(i)%raddirekt*(1.0d0-Albc)*BoundCell(i)%CanopyCell%RadPtrRate ! remained direct radiation after penetrating the canopy
      ! raddif=BoundCell(i)%raddiffus*(1.0d0-Albc)*BoundCell(i)%CanopyCell%GapFunc    ! remained diffuse radiation after penetrating the canopy
      ! radinf=BoundCell(i)%radinfred*BoundCell(i)%CanopyCell%GapFunc + &
      !   BoundCell(i)%CanopyCell%emissivity*SBsigma*BoundCell(i)%CanopyCell%AvgTem**4*(1.0d0-BoundCell(i)%CanopyCell%GapFunc)
      raddir=BoundCell(i)%CanopyCell%Rbs
      raddif=BoundCell(i)%CanopyCell%Rds
      radinf=BoundCell(i)%CanopyCell%RLs
      ! write(*,*) 'CanopySurfaceCompute: raddir, raddif, radinf.', raddir, raddif, radinf
!     incoming (direct) solar Radiation
      Qdir       = (One-(AlbSoil+AlbVeg))*raddir*shad &
                   *MAX(Zero,(n1G*radn1+n2G*radn2+n3G*radn3))
!     diffusive Radiation
      Qdif       = (One-(AlbSoil+AlbVeg))*raddif
      ! downward long wave radiation hitting the surface
      Qinf       = radinf

      ! write(*,*) 'CSC_Soil: 6'
      IF (cfcap(SoilType(1))>WS(1)) THEN
        hu=Half*(1.0d0-COS((WS(1)*Pi)/(1.6d0*cfcap(SoilType(1)))))
      ELSE
        hu=1.0d0
      END IF

      Eva_Pot       = -(RhoVAir-RhoV_Sat)*DragH*VT  
      Eva_Soil      = MAX(0.0d0,(One-F_plant-F_i)*(hu*RhoV_Sat-RhoVAir)*DragH*VT) 
      Eva_Intercept = F_i*DragH*VT*(RhoVAir-MIN(RhoV_Sat,RhoV_g)) 

      IF (WS(1)>=cporv(SoilType(1))) THEN
        I_perc = 0.0d0
      ELSE
        I_perc = (One-F_plant)*(Wgn+Eva_Intercept*2.0d0*dtMax/RhoW)*RhoW/tau_perc
      END IF

      DO iw=1,nzWg
        IF (WS(iw).LE.cpwp(SoilType(iw))) THEN
          WSroot = WSroot
        ELSE
          WSroot = WSroot+WS(iw)*MIN(dzSoil(iw), &
                                 MAX(0.d0,(z_root(LandClass)-zSoil(iw-1))))
        END IF
      END DO
      WSroot = WSroot/(z_root(LandClass)+Eps)

      DO iw=1,nzWg
        Transpiration(iw) = TranspirFunction(TAir,DragH*VT,u_shear,WSroot,    &
                                             Qdir,LAI,Eva_Pot,SoilType(iw),   &  
                                             LandClass,WS(iw),iw,F_plant,F_i) 
        SumTransp         = SumTransp+Transpiration(iw)
      END DO
  
!     Upper Boundary Conditions
      F_LH     = +SumTransp+Eva_Soil-Eva_Intercept ! [kg m-2 s-1], define F_LH positive when flux from soil -> air
      F_SH     = -DragH*VT*(TAir-TS(1)) ! [K m s-1], positive when heat flux is soil -> air
      ! write(*,*) 'CM: SumTransp, Eva_Soil, Eva_Intercept.', SumTransp, Eva_Soil, Eva_Intercept

      FSoil    = Qdir+Qdif+Qinf-ee*sigma*TS(1)**4.0d0-Cpair*F_SH*RhoAir-Lv*F_LH ! [W m-2]
      RhsTS(1) = RhsTS(1)+FSoil
      RhsWgn   = RhsWgn-(I_perc-Eva_Intercept)/RhoW   
      RhsWS(1) = RhsWS(1)-(Eva_Soil-I_perc)/RhoW
      RhsTS(1) = RhsTS(1)-(Eva_Soil-I_perc)/RhoW*CapWater*TS(1)
      ! write(*,*) 'ee.', ee

!     Fluxes for output in [W m-2]
      Qlat  = F_LH*Lv 
      Qsens = F_SH*Cpair*RhoAir  ! [kg s-3 = W m-2]

      BoundCell(i)%FluxLat  = Qlat 
      BoundCell(i)%FluxSens = Qsens 
      ! write(*,*) 'Soil_Mod: qlat, qsens.', Qlat, Qsens

      ! F_LH = Qlat/Lv
      ! F_SH = Qsens/(Cpair*RhoAir)

      IF (TS(1)<=T_freeze) RhsWgn = Zero 

      DO iw=1,nzThB
        RhsTS(iw) = (RhsTS(iw)-(RhsWS(iw)+Transpiration(iw))*TS(iw)*CapWater/RhoW) &
                    /dzSoil(iw)/Capacity(WS(iw),SoilType(iw))
      END DO
      DO iw=1,nzWg
        RhsWS(iw) = (RhsWS(iw)-Transpiration(iw)/RhoW)/dzSoil(iw)
        IF (TS(iw)<=T_freeze) RhsWS(iw) = Zero
      END DO
      IF (nrsl==1) THEN
        RhsTS(1)=RhsTS(1)-(TS(1)-TSClim)*Kg
        RhsWS(1)=RhsWS(1)+(0.08d0-WS(1))/tau_wg
      END IF

      TotalFlux = ThAir/TAir*F_SH*FL+F_LH*FL/RhoDAir*ThAir*(Rv/Rm-LOG(p/p0)*Rm/Cpml*(Rv/Rm-Cpv/Cpml)) ! [kg K s-1], ThAir = [kg m-3 K]
      dLoc      = ABS(n1)*dx(ix)+ABS(n2)*dy(iy)+ABS(n3)*dz(iz)
      VolLoc    = dx(ix)*dy(iy)*dz(iz)
      IF (Sphere) THEN
        dLoc    = ABS(n1)*dx(ix)*RadEarth+ABS(n2)*dy(iy)*RadEarth+ABS(n3)*dz(iz)
        VolLoc  = VolCellGlobeCanopy(zP(iz)+RadEarth,   &
                               zP(iz+1)+RadEarth, &
                               yP(iy-1),yP(iy),   &
                               xP(ix-1),xP(ix))
      END IF
      VolLoc    = MIN(FL*dLoc,VolLoc)
      Surface   = ABS(FU(ix-1,iy,iz)-FU(ix,iy,iz)) &
                 +ABS(FV(ix,iy-1,iz)-FV(ix,iy,iz)) &
                 +ABS(FW(ix,iy,iz-1)-FW(ix,iy,iz))
      VolDiff   = MAX(VolLoc-VolC(ix,iy,iz),Zero)
      ThRhs(ix,iy,iz,1)   = ThRhs(ix,iy,iz,1)+MIN(VolC(ix,iy,iz)/VolLoc,One)*TotalFlux
      ThRhs(ix-1,iy,iz,1) = ThRhs(ix-1,iy,iz,1)+MAX(FU(ix-1,iy,iz)-FU(ix,iy,iz),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      ThRhs(ix+1,iy,iz,1) = ThRhs(ix+1,iy,iz,1)+MAX(FU(ix,iy,iz)-FU(ix-1,iy,iz),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      ThRhs(ix,iy-1,iz,1) = ThRhs(ix,iy-1,iz,1)+MAX(FV(ix,iy-1,iz)-FV(ix,iy,iz),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      ThRhs(ix,iy+1,iz,1) = ThRhs(ix,iy+1,iz,1)+MAX(FV(ix,iy,iz)-FV(ix,iy-1,iz),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      ThRhs(ix,iy,iz-1,1) = ThRhs(ix,iy,iz-1,1)+MAX(FW(ix,iy,iz-1)-FW(ix,iy,iz),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      ThRhs(ix,iy,iz+1,1) = ThRhs(ix,iy,iz+1,1)+MAX(FW(ix,iy,iz)-FW(ix,iy,iz-1),Zero) &
                                                /Surface*VolDiff/VolLoc*TotalFlux
      RhoRhs(ix,iy,iz,1)  = RhoRhs(ix,iy,iz,1)+FL*F_LH
      RhoVRhs(ix,iy,iz,1) = RhoVRhs(ix,iy,iz,1)+FL*F_LH
      ThBRhs(i,1:nzThB)   = RhsTS
      WgRhs(i,nzWg+1)     = RhsWgn
      WgRhs(i,1:nzWg)     = RhsWS
    END IF
  END DO

END SUBROUTINE CanopySurfCompute_Soil


!=======================================!
!----- Calculation of drag coefficients
!----- Not used now
!=======================================!
SUBROUTINE DragCoeffCanopy(i)
  INTEGER :: i

  ! CALL DragCoeffCanopy_Wall(i)
  ! CALL DragCoeffCanopy_WallRad(i)
  CALL DragCoeffCanopy_Soil(i)
END SUBROUTINE DragCoeffCanopy


SUBROUTINE DragCoeffCanopy_Wall(i)
  
  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: FL,zRauh,zPL,logz
  REAL(RealKind) :: DragM,DragH
  REAL(RealKind) :: Prt

  REAL(RealKind) ,Parameter :: logmin=0.1d0,logmax=1.d10,zero=0.d0

! Turbulent Prandtl number  
  Prt = 0.95
  
  ix     = BoundCell(i)%ix
  iy     = BoundCell(i)%iy
  iz     = BoundCell(i)%iz
  FL     = BoundCell(i)%FL
  zRauh  = BoundCell(i)%zRauh
  
! Evaluate height z
  zPL = Half*VolC(ix,iy,iz)/FL
  zPL = Half*dz(iz)             !ELMAR
  IF (zRauh.GT.zero) THEN
    logz = MAX(LOG(zPL/zRauh+One),logmin)
  ELSE
    logz = logmax
  END IF
! Set Drag Coefficient, BoundCell(i)%DragH = 0 in default
  DragM = (Karm/logz)**Two
  BoundCell(i)%DragM  = DragM

END SUBROUTINE DragCoeffCanopy_Wall


!===== Calculation of drag coefficients for canopy ground (Copied from DragCoeffWallRad) =====!
SUBROUTINE DragCoeffCanopy_WallRad(i)
  
  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,zRauh,zRauhT,zPL,logz
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn1,Tgn2
  REAL(RealKind) :: RhoLoc,ThLoc,ThPotSurf,Fpot
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: RiB0,RiB1,RiBT
  REAL(RealKind) :: DragM,DragH
  REAL(RealKind) :: Qg,Qc,Qdir,Qem,Qref,Qdif,Qinf
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
  
! Evaluate velocity tangential and normal  
  uLoc = (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
          FU(ix,iy,iz)*uF(ix,iy,iz)) / &
        ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc+Eps))
  vLoc = (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
         FV(ix,iy,iz)*vF(ix,iy,iz)) / &
       ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc+Eps))
  wLoc = (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
         FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
       ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc+Eps))

  VN   = uLoc*n1+vLoc*n2+wLoc*n3
  V    = uLoc*uLoc+ &
         vLoc*vLoc+ &
         wLoc*wLoc
  VT   = SQRT(MAX(V-VN*VN,Zero))           

! Evaluate height z
  zPL = Half*VolC(ix,iy,iz)/FL
  zPL = Half*dz(iz)             !ELMAR
  
  IF (zRauh.GT.zero) THEN
    logz = MAX(LOG(zPL/zRauh+One),logmin)
  ELSE
    logz = logmax
  END IF
  
! First estimate for potential temperature
  Fpot = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
  ThPotSurf = Tgn1*Fpot  

  DO it=0,10      
    P = 9.24d0*((PrandtlNumber(thPos)/Prt)**0.75d0-1.0d0)*(1.0d0+0.28d0*EXP(-0.007d0*PrandtlNumber(thPos)/Prt))
    s = logz/Karm
!   Drag coefficient for heat and momentum
    DragH = 1.0d0/((s**2.0d0)*Prt*(1.0d0+P/s))
    DragM = (One/s)**2.0d0

!   ThPotSurf=Tgn1**(Cv/Cp)*Fpot
    ThPotSurf = Tgn1*Fpot

!   Bodenwaermestrom
    Qg   = -BoundCell(i)%D*(Tgn1-BoundCell(i)%T1)/(0.5d0)
!   Waermestrom in die Luft
    Qc   = -RhoLoc*Cpd*DragH*VT*(ThPotSurf-ThLoc)
!   direkte Strahlung
    Qdir = (1-alb)*shad*MAX(Zero,(n1*radn1+n2*radn2+n3*radn3))*BoundCell(i)%raddirekt
!   diffuse Strahlung
    Qdif = (1-alb)*BoundCell(i)%skyviewfactor*BoundCell(i)%raddiffus
!   downward long wave radiation
    Qinf = BoundCell(i)%raddiffus
!   Ausstrahlung im Infraroten
    Qem  = -ee*SBsigma*Tgn1**4.0d0

    Qref = 0.0d0

    F    = Qg+Qc+Qdir+Qdif+Qinf+Qem+Qref
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

END SUBROUTINE DragCoeffCanopy_WallRad


SUBROUTINE DragCoeffCanopy_Soil(i)

  INTEGER :: i,j
  INTEGER :: ix,iy,iz,it,shad
  REAL(RealKind) :: dt
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: n1G,n2G,n3G
  REAL(RealKind) :: FL,zRauh,zRauhT,zPL,logz
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Tgn
  REAL(RealKind) :: RhoLoc,ThLoc
  REAL(RealKind) :: RiB
  REAL(RealKind) :: DragM,DragH
  REAL(RealKind) :: FCh,FCm
  REAL(RealKind) :: TermM,TermH
  REAL(RealKind) :: Chn,Cmn,fh,fm
  REAL(RealKind) :: Ch,Cm
  REAL(RealKind) :: zL
  REAL(RealKind) :: RichCut
  REAL(RealKind) :: Chi
! REAL(RealKind) :: fx,dfdx
  REAL(RealKind) ,Parameter :: RichCrit=0.25d0


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
    RhoLoc = Rho(ix,iy,iz,1)
    Tgn    = BoundCell(i)%TeS
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

!   Evaluate height z
    zPL  = Half*VolC(ix,iy,iz)/FL
    zPL  = Half*dz(iz)             !ELMAR

!   Bulk-Richardson-Number
    RiB  = Grav/Tgn*((ThLoc-Tgn)*(dz(iz)-zRauh))/(VT*VT+Eps)

  ! write(*,*) 'DCC_Soil: 3'
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
      FCh  = SQRT(Grav/Tgn*ABS(ThLoc-Tgn)*(dz(iz)-zRauh))/ &
             (cc*((dz(iz)/zRauhT)**(One/Three)-1.0d0)**(Three/Two))
      FCm  = Two/Three*SQRT(Grav/Tgn*ABS(ThLoc-Tgn)*(dz(iz)-zRauh))/ &
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

!      CALL FirstGuess(RiB,Chi,dz(iz),zRauh,zRauhT)
      ! write(*,*) 'revWRF'
      Chi=RiB
      DO
        CALL NewtonMethod(Chi,zL,fx(Chi),dfdx(Chi)) 
        IF (ABS(Chi-zL)<1.d-8)THEN
          EXIT
        ELSE
          Chi=zL
        END IF
      END DO
!     TermH = LOG((dz(iz)+zRauh)/zRauhT)      &
!            -Phi('H',zL*(One+zRauh/dz(iz)),dz(iz),zRauh,zRauhT,RiB) &
!            +Phi('H',zL*zRauhT/dz(iz),dz(iz),zRauh,zRauhT,RiB)
!     TermM = LOG((dz(iz)+zRauh)/zRauh)       &
!            -Phi('M',zL*(One+zRauh/dz(iz)),dz(iz),zRauh,zRauhT,RiB) &
!            +Phi('M',zL*zRauh/dz(iz),dz(iz),zRauh,zRauhT,RiB)
      Cm    = Karm**Two/TermM/TermM
      Ch    = Karm**Two/TermM/TermH
    END IF
  
    BoundCell(i)%DragM  = Cm
    BoundCell(i)%DragH  = Ch
    BoundCell(i)%DragQ  = Ch

    CONTAINS
   
      FUNCTION fx(x)
        REAL(RealKind) :: fx
        REAL(RealKind) :: x
        fx=RiB-RiBulk(dz(iz),zRauh,zRauhT,x,RiB)
      END FUNCTION fx
  
      FUNCTION dfdx(x)
        REAL(RealKind) :: dfdx
        REAL(RealKind) :: x
!       dfdx=-dRiBdChi(dz(iz),zRauh,zRauhT,x,RiB)
      END FUNCTION dfdx
 
END SUBROUTINE DragCoeffCanopy_Soil


!=======================================!
!----- Drag for momentum equations
!=======================================!
SUBROUTINE CanopyDrag(Vector, Rhs, UVec)
  TYPE(Vector4Cell_T) :: Vector,Rhs
  TYPE (Vector4Cell_T), OPTIONAL :: UVec

  Rho=>RhoCell(ibLoc)%c
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
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  wRhsR=>Rhs%Vec(wPosR)%c

  CALL CanopyDragCompute
END SUBROUTINE CanopyDrag


SUBROUTINE CanopyDragCompute

  INTEGER :: i,j,ix,iy,iz,nrcl
  REAL(RealKind) :: TempU,TempV,TempW,TempVM ! velocity components and mean velocity
  REAL(RealKind) :: A
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Cd = 0.15d0 ! drag coefficient of the canopy

  DO i=1,NumBoundCell
    ix   =  BoundCell(i)%ix
    iy   =  BoundCell(i)%iy
    iz   =  BoundCell(i)%iz
    nrcl =  BoundCell(i)%CanopyCell%NrCanopyLayers
    LAD  => BoundCell(i)%CanopyCell%LAD

    DO j=iz, iz+nrcl-1
        TempU=(uCR(ix,iy,j,1)*FU(ix,iy,j)+uCL(ix,iy,j,1)*FU(ix-1,iy,j))/(FU(ix-1,iy,j)+FU(ix,iy,j)+Eps)
        TempV=(vCR(ix,iy,j,1)*FV(ix,iy,j)+vCL(ix,iy,j,1)*FV(ix,iy-1,j))/(FV(ix,iy-1,j)+FV(ix,iy,j)+Eps)
        TempW=(wCR(ix,iy,j,1)*FW(ix,iy,j)+wCL(ix,iy,j,1)*FW(ix,iy,j-1))/(FW(ix,iy,j-1)+FW(ix,iy,j)+Eps)
        TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)
        A=0.5*LAD(j-iz+1) ! here LAD is the one-sided value

        uRhsL(ix,iy,j,1)=uRhsL(ix,iy,j,1)-Cd*A*TempVM*TempU/(Rho(ix,iy,j,1)+Eps)
        uRhsR(ix,iy,j,1)=uRhsR(ix,iy,j,1)-Cd*A*TempVM*TempU/(Rho(ix,iy,j,1)+Eps)
        vRhsL(ix,iy,j,1)=vRhsL(ix,iy,j,1)-Cd*A*TempVM*TempV/(Rho(ix,iy,j,1)+Eps)
        vRhsR(ix,iy,j,1)=vRhsR(ix,iy,j,1)-Cd*A*TempVM*TempV/(Rho(ix,iy,j,1)+Eps)
        wRhsL(ix,iy,j,1)=wRhsL(ix,iy,j,1)-Cd*A*TempVM*TempW/(Rho(ix,iy,j,1)+Eps)
        wRhsR(ix,iy,j,1)=wRhsR(ix,iy,j,1)-Cd*A*TempVM*TempW/(Rho(ix,iy,j,1)+Eps)
    END DO
  END DO 

END SUBROUTINE CanopyDragCompute


SUBROUTINE JacCanopyDrag(Vector,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  Rho=>RhoCell(ibLoc)%c
  uCL=>Vector%Vec(uPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCL=>Vector%Vec(vPosL)%c
  vCR=>Vector%Vec(vPosR)%c
  wCL=>Vector%Vec(wPosL)%c
  wCR=>Vector%Vec(wPosR)%c

  AS=>Jac

  CALL JacCanopyDragCompute

END SUBROUTINE JacCanopyDrag


SUBROUTINE JacCanopyDragCompute

  INTEGER :: i,j,ix,iy,iz,nrcl
  REAL(RealKind) :: TempU,TempV,TempW,TempVM,TempVM2
  REAL(RealKind) :: KUL,KUR,KVL,KVR,KWL,KWR
  REAL(RealKind) :: VMRho,UCUCVM2,UCVCVM2,UCWCVM2,VCVCVM2,VCWCVM2,WCWCVM2
  REAL(RealKind) :: ULUL,ULUR,ULVL,ULVR,ULWL,ULWR,VLUL,VLUR,VLVL,VLVR,VLWL,VLWR,WLUL,WLUR,WLVL,WLVR,WLWL,WLWR
  REAL(RealKind) :: A
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Cd = 0.15d0 ! drag coefficient of the canopy

  DO i=1,NumBoundCell
    ix   =  BoundCell(i)%ix
    iy   =  BoundCell(i)%iy
    iz   =  BoundCell(i)%iz
    nrcl =  BoundCell(i)%CanopyCell%NrCanopyLayers
    LAD  => BoundCell(i)%CanopyCell%LAD

    DO j=iz, iz+nrcl-1
      TempU=(uCR(ix,iy,j,1)*FU(ix,iy,j)+uCL(ix,iy,j,1)*FU(ix-1,iy,j))/(FU(ix-1,iy,j)+FU(ix,iy,j)+Eps)
      TempV=(vCR(ix,iy,j,1)*FV(ix,iy,j)+vCL(ix,iy,j,1)*FV(ix,iy-1,j))/(FV(ix,iy-1,j)+FV(ix,iy,j)+Eps)
      TempW=(wCR(ix,iy,j,1)*FW(ix,iy,j)+wCL(ix,iy,j,1)*FW(ix,iy,j-1))/(FW(ix,iy,j-1)+FW(ix,iy,j)+Eps)
      TempVM2=TempU*TempU+TempV*TempV+TempW*TempW
      TempVM=SQRT(TempVM2)
      A=0.5*LAD(j-iz+1)

      KUL=FU(ix-1,iy,j)/(FU(ix-1,iy,j)+FU(ix,iy,j))
      KUR=FU(ix,iy,j)  /(FU(ix-1,iy,j)+FU(ix,iy,j))
      KVL=FV(ix,iy-1,j)/(FV(ix,iy-1,j)+FV(ix,iy,j))
      KVR=FV(ix,iy,j)  /(FV(ix,iy-1,j)+FV(ix,iy,j))
      KWL=FW(ix,iy,j-1)/(FW(ix,iy,j-1)+FW(ix,iy,j))
      KWR=FW(ix,iy,j)  /(FW(ix,iy,j-1)+FW(ix,iy,j))

      VMRho=TempVM/(Rho(ix,iy,j,1)+Eps)
      UCUCVM2=TempU*TempU/(TempVM2+Eps)
      UCVCVM2=TempU*TempV/(TempVM2+Eps)
      UCWCVM2=TempU*TempW/(TempVM2+Eps)
      VCVCVM2=TempV*TempV/(TempVM2+Eps)
      VCWCVM2=TempV*TempW/(TempVM2+Eps)
      WCWCVM2=TempW*TempW/(TempVM2+Eps)

      ULUL=-Cd*A*VMRho*KUL*(UCUCVM2+1)
      ULUR=-Cd*A*VMRho*KUR*(UCUCVM2+1)
      ULVL=-Cd*A*VMRho*KVL*UCVCVM2
      ULVR=-Cd*A*VMRho*KVR*UCVCVM2
      ULWL=-Cd*A*VMRho*KWL*UCWCVM2
      ULWR=-Cd*A*VMRho*KWR*UCWCVM2

      VLUL=-Cd*A*VMRho*KUL*UCVCVM2
      VLUR=-Cd*A*VMRho*KUR*UCVCVM2
      VLVL=-Cd*A*VMRho*KVL*(VCVCVM2+1)
      VLVR=-Cd*A*VMRho*KVR*(VCVCVM2+1)
      VLWL=-Cd*A*VMRho*KWL*VCWCVM2
      VLWR=-Cd*A*VMRho*KWR*VCWCVM2

      WLUL=-Cd*A*VMRho*KUL*UCWCVM2
      WLUR=-Cd*A*VMRho*KUR*UCWCVM2
      WLVL=-Cd*A*VMRho*KVL*VCWCVM2
      WLVR=-Cd*A*VMRho*KVR*VCWCVM2
      WLWL=-Cd*A*VMRho*KWL*(WCWCVM2+1)
      WLWR=-Cd*A*VMRho*KWR*(WCWCVM2+1)

      AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,j,1)+ULUL
      AS(IndexMet(uPosLJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,uPosRJac))%c(ix,iy,j,1)+ULUR
      AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,j,1)+ULVL
      AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,j,1)+ULVR
      AS(IndexMet(uPosLJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,wPosLJac))%c(ix,iy,j,1)+ULWL
      AS(IndexMet(uPosLJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosLJac,wPosRJac))%c(ix,iy,j,1)+ULWR

      AS(IndexMet(uPosRJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,uPosLJac))%c(ix,iy,j,1)+ULUL
      AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,j,1)+ULUR
      AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,j,1)+ULVL
      AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,j,1)+ULVR
      AS(IndexMet(uPosRJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,wPosLJac))%c(ix,iy,j,1)+ULWL
      AS(IndexMet(uPosRJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(uPosRJac,wPosRJac))%c(ix,iy,j,1)+ULWR

      AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,j,1)+VLUL
      AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,j,1)+VLUR
      AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,j,1)+VLVL
      AS(IndexMet(vPosLJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,vPosRJac))%c(ix,iy,j,1)+VLVR
      AS(IndexMet(vPosLJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,wPosLJac))%c(ix,iy,j,1)+VLWL
      AS(IndexMet(vPosLJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosLJac,wPosRJac))%c(ix,iy,j,1)+VLWR
      
      AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,j,1)+VLUL
      AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,j,1)+VLUR
      AS(IndexMet(vPosRJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,vPosLJac))%c(ix,iy,j,1)+VLVL
      AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,j,1)+VLVR
      AS(IndexMet(vPosRJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,wPosLJac))%c(ix,iy,j,1)+VLWL
      AS(IndexMet(vPosRJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(vPosRJac,wPosRJac))%c(ix,iy,j,1)+VLWR

      AS(IndexMet(wPosLJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,uPosLJac))%c(ix,iy,j,1)+WLUL
      AS(IndexMet(wPosLJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,uPosRJac))%c(ix,iy,j,1)+WLUR
      AS(IndexMet(wPosLJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,vPosLJac))%c(ix,iy,j,1)+WLVL
      AS(IndexMet(wPosLJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,vPosRJac))%c(ix,iy,j,1)+WLVR
      AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,j,1)+WLWL
      AS(IndexMet(wPosLJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosLJac,wPosRJac))%c(ix,iy,j,1)+WLWR

      AS(IndexMet(wPosRJac,uPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,uPosLJac))%c(ix,iy,j,1)+WLUL
      AS(IndexMet(wPosRJac,uPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,uPosRJac))%c(ix,iy,j,1)+WLUR
      AS(IndexMet(wPosRJac,vPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,vPosLJac))%c(ix,iy,j,1)+WLVL
      AS(IndexMet(wPosRJac,vPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,vPosRJac))%c(ix,iy,j,1)+WLVR
      AS(IndexMet(wPosRJac,wPosLJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,wPosLJac))%c(ix,iy,j,1)+WLWL
      AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,j,1)=AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,j,1)+WLWR
    END DO 
  END DO 

END SUBROUTINE JacCanopyDragCompute


!=========================================!
!----- Initiation of canopy stuff
!=========================================!
SUBROUTINE ReadCanopy(FileName)
  CHARACTER(*) :: FileName
  CHARACTER(300) :: Line
  INTEGER :: ix,iy,iz,i
  INTEGER :: NrCanopyLayers_temp
  REAL(RealKind) :: FakeBeta(10)
  REAL(RealKind), PARAMETER :: a=1.6d0, b=0.9d0
  REAL(RealKind), ALLOCATABLE :: LAD_temp(:)
  REAL(RealKind), POINTER :: LAD_pt(:)

  ! FakeBeta = (/2.2d0, 2.6d0, 3.35d0, 4.7d0, 6.25d0, 7.25d0, 7.55d0, 7.2d0, 6.05d0, 2.6d0/)
  !===== Read LAD from the grid file =====!
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'#CanopyLayers')>0) THEN
      READ(InputUnit,*) NrCanopyLayers_temp
      ALLOCATE(LAD_temp(NrCanopyLayers_temp))
      READ(InputUnit,*) LAD_temp ! maybe here we need a test for the length of LAD arrays
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

  DO i=1,NumBoundCell
    iz = BoundCell(i)%iz
    ! Set number of canopy layers
    BoundCell(i)%CanopyCell%NrCanopyLayers=NrCanopyLayers_temp
    ! Set LAD values for every layer
    ALLOCATE(BoundCell(i)%CanopyCell%LAD(NrCanopyLayers_temp))
    BoundCell(i)%CanopyCell%LAD = LAD_temp
    ! Calculate LAI for each boundary cell considering the height of every layer
    BoundCell(i)%CanopyCell%LAI = SUM(BoundCell(i)%CanopyCell%LAD*dz(iz:iz+NrCanopyLayers_temp-1))
    BoundCell(i)%CanopyCell%GapFunc = EXP(-0.5d0*a*BoundCell(i)%CanopyCell%LAI**b)
    ! Allocate arrays for leaf temperature
    ALLOCATE(BoundCell(i)%CanopyCell%SunlitLeafT(NrCanopyLayers_temp))
    ALLOCATE(BoundCell(i)%CanopyCell%ShadedLeafT(NrCanopyLayers_temp))
    BoundCell(i)%CanopyCell%SunlitLeafT = BoundCell(i)%TeS
    BoundCell(i)%CanopyCell%ShadedLeafT = BoundCell(i)%TeS
    ! Allocate arrays for H and LE
    ALLOCATE(BoundCell(i)%CanopyCell%Hc(NrCanopyLayers_temp))
    ALLOCATE(BoundCell(i)%CanopyCell%LEc(NrCanopyLayers_temp))
    BoundCell(i)%CanopyCell%Hc = 0.0d0
    BoundCell(i)%CanopyCell%LEc = 0.0d0
  END DO

  !===== Initiate soil stuff =====!
  ! ALLOCATE(dzSoil(1:Domain%nrsoillayers))
  ! ALLOCATE(zSoil(0:Domain%nrsoillayers))

  ! zSoil(0:Domain%nrsoillayers)=Domain%zSDepth(0:Domain%nrsoillayers)
  ! DO i=1,Domain%nrsoillayers
  !   dzSoil(i)=zSoil(i)-zSoil(i-1)
  ! END DO
  ! nzS=Domain%nrsoillayers
  ! nzW=Domain%nrsoillayers

END SUBROUTINE ReadCanopy

SUBROUTINE InitCanopy

  INTEGER :: i

  ALLOCATE(dzSoil(1:Domain%nrsoillayers))
  ALLOCATE(zSoil(0:Domain%nrsoillayers))

  zSoil(0:Domain%nrsoillayers)=Domain%zSDepth(0:Domain%nrsoillayers)
  DO i=1,Domain%nrsoillayers
    dzSoil(i)=zSoil(i)-zSoil(i-1)
  END DO
  nzS=Domain%nrsoillayers
  nzW=Domain%nrsoillayers

END SUBROUTINE InitCanopy

SUBROUTINE VectorInitCanopyFunction(Pos,Vector,Val,Time)
!==================================================
!----  Initialization of Scalar Values
!==================================================
  INTEGER :: Pos
  TYPE (Vector4Cell_T), TARGET :: Vector(:)
  REAL(RealKind) :: Time
  REAL(RealKind) :: Val
  EXTERNAL Val
  INTEGER :: i,iS
  REAL(RealKind) :: zPS
  REAL(RealKind), POINTER :: cB(:,:) 
    
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    cB=>Vector(ibLoc)%Vec(Pos)%cB
    DO i=1,NumBoundCell
      IF (Pos==RhoVPos) THEN
        DO iS=1,SIZE(cB(i,:))-1
          zPS=zSoil(is)-Half*dzSoil(is)
          cB(i,iS) = Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                        ,zH(BoundCell(i)%ix,BoundCell(i)%iy),zPS &
!                        ,zH(BoundCell(i)%xS,BoundCell(i)%yS),zPS &
                        ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(iS),Time)
        END DO
          cB(i,SIZE(cB(i,:)))=Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                                 ,zH(BoundCell(i)%ix,BoundCell(i)%iy),Zero &
!                                 ,zH(BoundCell(i)%xS,BoundCell(i)%yS),Zero &
                                 ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(1),Time)
      ELSE
        DO iS=1,SIZE(cB(i,:))
          zPS=zSoil(is)-Half*dzSoil(is)
          cB(i,iS) = Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                        ,zH(BoundCell(i)%ix,BoundCell(i)%iy),zPS &
!                        ,zH(BoundCell(i)%xS,BoundCell(i)%yS),zPS &
                        ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(iS),Time)
        END DO
      END IF
    END DO
  END DO   ! ibLoc

END SUBROUTINE VectorInitCanopyFunction

!=========================================!
!----- Soil functions
!=========================================!
FUNCTION KappaSoil(Theta_Soil,SoilType)
  REAL(RealKind) :: Theta_Soil
  REAL(RealKind) :: KappaSoil,m,n,WSeffektiv
  INTEGER :: SoilType
  WSEffektiv= (MAX(Eps,(Theta_Soil-cres(SoilType))))/(cporv(SoilType)-cres(SoilType))
  n         = lamb(SoilType)+1
  m         = 1.0d0-1.0d0/n
  SELECT CASE(SoilParam)
    CASE('Campbell')
      KappaSoil = ckw0(SoilType)*(Theta_Soil/cporv(SoilType))**(2.d0*cbedi(SoilType)+3) 
    CASE('CoreyBrooks')
      KappaSoil = ckw0(SoilType)*WSEffektiv**(2.d0*cbedi(SoilType)+3) 
    CASE('VanGenuchten')
      KappaSoil = ckw0(SoilType)*WSEffektiv**(One/Two)*(One-(One-WSEffektiv**(One/m))**m)**Two  
    CASE DEFAULT
      KappaSoil = ckw0(SoilType)*WSEffektiv**(One/Two)*(One-(One-WSEffektiv**(One/m))**m)**Two  
  END SELECT
  KappaSoil = MAX(Eps,KappaSoil)
END FUNCTION KappaSoil

FUNCTION PsiSoil(Theta_Soil,SoilType)
  REAL(RealKind) :: Theta_Soil,WSeffektiv
  REAL(RealKind) :: PsiSoil
  REAL(RealKind) :: m,n
  INTEGER :: SoilType
  n         = lamb(SoilType)+1 
  m         = 1.0d0-1.0d0/n
  WSEffektiv= (MAX(Eps,(Theta_Soil-cres(SoilType))))/(cporv(SoilType)-cres(SoilType))
  SELECT CASE(SoilParam)
    CASE('Campbell')
      PsiSoil = smPot(SoilType)*(MAX(Eps,Theta_Soil)/cporv(SoilType))**(-cbedi(SoilType)) 
    CASE('CoreyBrooks')
      PsiSoil   = smPot(SoilType)*WSEffektiv**(-cbedi(SoilType)) 
    CASE('VanGenuchten')
      PsiSoil   = smPot(SoilType)*(WSeffektiv**(-One/m)-One)**(One/n) 
    CASE DEFAULT
      PsiSoil   = smPot(SoilType)*(WSeffektiv**(-One/m)-One)**(One/n) 
  END SELECT
END FUNCTION PsiSoil

FUNCTION FlussCond(Theta_SoilT,Theta_SoilB,dzT,dzB,SoilTypeT,SoilTypeB) ! Conductive Flux
  REAL(RealKind) :: FlussCond
  REAL(RealKind) :: Theta_SoilT,Theta_SoilB,dzT,dzB
  INTEGER :: SoilTypeT,SoilTypeB
  FlussCond = 2.0d0/(dzT/KappaSoil(Theta_SoilT,SoilTypeT) &
                     +dzB/KappaSoil(Theta_SoilB,SoilTypeB)) &
              *(PsiSoil(Theta_SoilT,SoilTypeT)-PsiSoil(Theta_SoilB,SoilTypeB))  &
              +KappaSoil(Theta_SoilT,SoilTypeT)
END FUNCTION FlussCond

FUNCTION FlussDiffT(TsT,TsB,dzT,dzB,Theta_SoilT,Theta_SoilB,SoilTypeT,SoilTypeB)
  REAL(RealKind) :: FlussDiffT
  REAL(RealKind) :: Theta_SoilT,Theta_SoilB
  REAL(RealKind) :: TsT,TsB,dzT,dzB
  INTEGER :: SoilTypeT,SoilTypeB
  FlussDiffT = 2.0d0/(dzT/ThermCond(Theta_SoilT,SoilTypeT)   &
                     +dzB/ThermCond(Theta_SoilB,SoilTypeB))*(TsT-TsB) &
               +2.0d0/(dzT/KappaSoil(Theta_SoilT,SoilTypeT)  & 
                      +dzB/KappaSoil(Theta_SoilB,SoilTypeB)) &
               *(PsiSoil(Theta_SoilT,SoilTypeT)-PsiSoil(Theta_SoilB,SoilTypeB)) &
               *capwater*(dzT+dzB)/(dzT/TsT+dzB/TsB) &
               +KappaSoil(Theta_SoilT,SoilTypeT)*capwater*TsT
END FUNCTION FlussDiffT

FUNCTION ThermCond(Theta_Soil,SoilType) ! Thermal Conductivity
  REAL(RealKind) :: Theta_Soil,Psi
  REAL(RealKind) :: ThermCond
  INTEGER  :: SoilType
  Psi       = PsiSoil(Theta_Soil,SoilType) 
  IF (SoilType>2) THEN
    ThermCond = MAX(418.0d0*EXP(-LOG10(100.0d0*ABS(Psi))-2.7d0),0.172d0)
  ELSE 
    ThermCond = 0.172d0
  END IF
END FUNCTION ThermCond

FUNCTION Capacity(Theta_Soil,SoilType)
  REAL(RealKind) :: Capacity
  REAL(RealKind) :: Theta_Soil
  INTEGER :: SoilType
! soil capacity = sum of capacity of dry soil, capacity of water and capacity of air (in pores)
  Capacity = (1.0d0-cporv(SoilType))*crhoc(SoilType)+Theta_Soil*capwater+(cporv(SoilType)-Theta_Soil)*capair
END FUNCTION Capacity

FUNCTION  TranspirFunction(TAir,DragCoeff,Ushear,WSroot,RadPAR,   &
                           LAIndex,Epotg,SoilType,LandType,WSoil, &
                           iS,PlantArea,InterceptArea)
  REAL(RealKind) :: TranspirFunction
  REAL(RealKind) :: AtmoResist,FolResist
  REAL(RealKind) :: ReductTransp,Resist_la
  REAL(RealKind) :: StomaResist
  REAL(RealKind) :: DragCoeff
  REAL(RealKind) :: RadPAR
  REAL(RealKind) :: Ushear
  REAL(RealKind) :: Frad,Fwat,Ftem,Fhum
  REAL(RealKind) :: TAir
  REAL(RealKind) :: WSroot,WSoil
  REAL(RealKind) :: wTLP,Epotg
  REAL(RealKind) :: LAIndex
  REAL(RealKind) :: PlantArea,InterceptArea
  INTEGER        :: SoilType,LandType,iS
! turgor loss point of plants:
  wTLP = cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType)) &
!                        *(0.81d0+0.121*ATAN(86400.d0*Epotg-Epotnorm))
                        *(0.81d0+0.121*ATAN(Epotg-Epotnorm/86400.d0))
! describes radiative influence on stomatal resistance,
  Frad = MIN(1.0d0,0.55d0*RadPAR/RadPARcrit) ! RadPAR==0.55*Qdir
! ... of soil water content,
  Fwat = MAX(Zero,MIN(1.0d0,(WSroot-cpwp(SoilType))/(wTLP-cpwp(SoilType))))
! ... of ambient specific humidity,
  Fhum = 1.0d0
! ... and of ambient temperature.
  Ftem = Max(Zero,MIN(1.0d0,4.0d0*((TAir-Temp0)*(TempEnd-TAir))/(TempEnd-Temp0)**2.0d0))
! resistance between the atmosphere and the air between the leafs
  AtmoResist   = DragCoeff**(-1.0d0)
! resistance between leafs and the air between the leafs
  Resist_la    = (MAX(Cdash*Ushear**(Half),1.0d-6))**(-1.0d0)
! stomatal resistance
  StomaResist  = 1.0d0/( ResMax**(-1.0d0)+(ResMin**(-1.0d0)-ResMax**(-1.0d0))      &
                                *Frad*Fwat*Fhum*Ftem )
! reduction of transpiration due to stomatal resistance
  ReductTransp = Resist_la*(Resist_la+StomaResist)**(-1.0d0)
! total resistance = foliage resistance
  FolResist    = 1.0d0/MAX(Eps,(ReductTransp*LAIndex*Resist_la**(-1.0d0)))
  IF (WSoil.LE.cpwp(SoilType)) THEN
    TranspirFunction = Zero
!   WRITE (*,*) 'Tr ZERO'
!   WRITE (*,*) ' WSoil,cpwp(SoilType),SoilType', WSoil,cpwp(SoilType),SoilType
  ELSE
    TranspirFunction = (One-InterceptArea)*PlantArea*            &
                        Epotg*AtmoResist/(AtmoResist+FolResist)  &
                       *MIN(dzSoil(iS),MAX(0.d0,(z_root(LandType)-zSoil(iS-1)))) &
                       /(z_root(LandType)+Eps)  &
                       *WSoil/(WSroot+Eps)
!   WRITE (*,*) 'Tr calc'
  END IF
END FUNCTION TranspirFunction

FUNCTION VolCellGlobeCanopy(r1,r2,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: VolCellGlobeCanopy
  REAL(RealKind) :: r1,r2,phi1,phi2,lam1,lam2
  VolCellGlobeCanopy=MAX(1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1))),Zero)
END FUNCTION VolCellGlobeCanopy


END MODULE Canopy_Mod
