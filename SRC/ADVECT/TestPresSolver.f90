PROGRAM MainProg

  USE Domain_Mod
  USE Init_Mod
  USE JacAccGrav_Mod
  USE Int_Mod
! USE IntPeer_Mod
! USE IntLinPeer_Mod
  USE ReadOutput_Mod
  USE Emission_Mod
  USE Koagulation_Mod
  USE ScalarVectorCellPar_Mod
  USE Microphysics_Mod
  USE Activity_Mod
  USE RhoProf_Mod
  USE IntSub_Mod
  USE Turbulence_Mod
  USE Operator_Mod
  USE Parameter_Mod
  USE Physics_Mod
  USE Tools_Mod
  USE SpecialOutput_Mod
  USE Soil_Mod
  USE Canopy_Mod
  USE ReadWeights_Mod
  USE WindFarm_Mod
  USE Forcing_Mod
  USE Restart_Mod
  USE Grid_Mod
  USE Example_Mod
  USE BoundaryCondition_Mod
  USE OutputMeanProfile_Mod
  USE SingleColumnXProfile_Mod
  USE SingleColumnZProfile_Mod
  USE MeanSurface_Mod
  USE PointSurface_Mod 
  USE TimeStep_Mod
  USE ReadWRF_Mod
  USE ReadWRFnc_Mod
  USE Rhs_Mod
  USE Diagnostic_Mod

!--------------Begin Block-----------------
  USE Subdom_Mod, ONLY: SubDomain
  USE Boundary_Mod
  USE DomDecomp_Mod, ONLY: read_blocks_from_netcdf, vec1d_from_homog, btype, blockranks
  USE Laplace_block_Mod, ONLY: make_operators_static
  USE Compressible_Mod, ONLY: solveM, setup_pressure_solver, update_pressure_solver_dt
  USE Output_Mod, ONLY: output_serial
  USE List_Mod, ONLY: RArray

  USE CoarseGrids_Mod, ONLY: define_coarse_grids, partition_clevels, intialize_coarse_subdomains, &
                             subdomain_lev, isblackijkzero_lev, determine_color_startcell

!--------------End Block-----------------

  IMPLICIT NONE

  TYPE(VelocityFace_T), POINTER :: VelF1(:)
  TYPE(VelocityFace_T), POINTER :: VelF2(:)
  TYPE(VecVelocityFace_T), POINTER :: VecVelF(:)
  TYPE(Vector4Cell_T), POINTER :: VecMet(:)
  TYPE(Vector4Cell_T), POINTER :: VecChem(:)
  TYPE(Vector4Cell_T), POINTER :: VecG(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVecMet(:)
  REAL(RealKind) :: Temp,RhoLoc,rRand
  INTEGER :: Iter,iForce
  INTEGER :: iInt,ix,iy,iz
  INTEGER :: ixVol,iyVol,izVol
  Real(RealKind) :: VolMin
  INTEGER :: i,j,k,iShift
  INTEGER :: iFrac
  INTEGER :: stat
  CHARACTER(80) :: ProfileSTART='ProfileSTART'
  CHARACTER(80) :: ProfileEND='ProfileEND'
  CHARACTER(80) :: InputFileName
  CHARACTER(8) :: ScalarName
  INTEGER :: iFort
  INTEGER :: unr
  INTEGER :: PosDummy
  REAL(RealKind) :: LocTime,MaxTime
  !REAL(RealKind) :: OutputTime !!
  REAL(RealKind) :: dtAct,Time
  REAL(RealKind) :: dtCFL=0.0d0
  INTEGER :: nsCFL

  LOGICAL :: FileExist

  REAL(RealKind) :: TotalStart(1:5),TotalCurrent(1:5)
  INTEGER :: SizeSeed=12
  !INTEGER,Dimension(12) :: seed
  INTEGER, ALLOCATABLE :: seed(:)
  CHARACTER :: ArgRestart

  REAL(RealKind) :: ErrorL1,SumVol
  REAL(RealKind) :: Start,Finish

! Outputs
  LOGICAL :: CheckMean
  LOGICAL :: CheckColumnX 
  LOGICAL :: CheckColumnZ
  LOGICAL :: CheckSurfMean
  LOGICAL :: CheckSurfPoint
 
  LOGICAL :: FinishNudging
  LOGICAL :: FinishTendency
!--------------Begin Block-----------------
! CHARACTER(len=128), PARAMETER :: file_path = "CanyonHof.nc"
  CHARACTER(len=128), PARAMETER :: file_path = "Barn_on_left_high_res3D.nc"
  CHARACTER(len=3) :: file_nr
  TYPE(SubDomain), TARGET :: subdom
  REAL(Realkind), ALLOCATABLE :: rhoM(:), rhothetaMi(:)
  REAL(Realkind), ALLOCATABLE :: rhoM_New(:), rhothetaMi_New(:)
  REAL(Realkind), ALLOCATABLE, TARGET :: rhovelM(:)
  REAL(Realkind), ALLOCATABLE, TARGET :: rhovelM_New(:)
  TYPE(SubDomain), POINTER :: this_subdomain => NULL()
  REAL(Realkind), PARAMETER :: urho0M = 6.0d0, &
                               vrho0M = 6.0d0, &
                               wrho0M = 0.0d0, &
                               rho0M = 1.2d0, &
                               theta0M = 283.0d0

  REAL(RealKind) :: dtM=0.1d0
  REAL(Realkind), PARAMETER :: gamm = 1.0d0
  REAL(Realkind), PARAMETER :: sor_param = 1.5d0
  Real(Realkind) :: resmax_tol = 1.0d-10
  INTEGER, PARAMETER :: NumIter = 100                           
  INTEGER :: n,icellM
  LOGICAL :: MultiGridMicha = .TRUE.
!--------------End Block------------------

  CALL start_MPI


! -- Lesen der Gitterdatei --
  CALL get_command_argument(1,InputFileName)

  CALL MPI_Bcast(InputFileName,80,MPI_CHARACTER,0,MPI_COMM_WORLD,MPIErr)
  IF (MyId==0) THEN
    WRITE(*,*)
    WRITE(*,*) InputFileName
  END IF

! -- Compute Parameter
  CALL ComputeParameter 
! -- Einlesen des Gitters (Block-Struktur) --
  CALL Allocate(Floor)
  CALL inp_part(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL InputExample(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  WRITE(*,*) 'Pres 1'


  CALL Allocate(Floor)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  WRITE(*,*) 'Pres 11'
  CALL ReadWeights(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL MinVol

  WRITE(*,*) 'Pres 111'
  CALL InputModelTransport(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)

  WRITE(*,*) 'Pres 12'
  IF (ChemieFile/='') THEN 
    CALL InputSystem(ChemieFile)
  END IF

  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL SetIndices
  CALL Allocate(JacTrans)
  WRITE(*,*) 'Pres 13'
  CALL AllocateVec4Met(VecMet,VectorComponentsMet)

  ALLOCATE(blockranks(nb))
  DO ib = 1,nb
    blockranks(ib)=blMPI(ib)%proc
  END DO  



  WRITE(*,*) 'Pres 2'
  VectorComponentsM=VectorComponentsT
  CALL InputModelBC(InputFileName)

  IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
      BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
      BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN 
    ALLOCATE(MeanProfile(0:VectorComponentsT,Domain%ix0+1:Domain%ix1))
  END IF

  VecMet=Zero
  cyclic_x = .FALSE.
  cyclic_y = .FALSE.
  cyclic_z = .FALSE.
  IF (BCVel%West=='Period') THEN
    cyclic_x = .TRUE.
  END IF  

  WRITE(*,*) '-- Meteorologie --'
! -- Meteorologie --
  CALL Allocate(VelF1)
  CALL Allocate(VelF2)
  VelF1=Zero
  VelF2=Zero
  VecMet=Zero

  WRITE(*,*) 'Pres 3'
  IF (RhoPos>0) THEN
    ALLOCATE(RhoCell(nbLoc))
    CALL Assign(RhoCell,VecMet,RhoPos)
  END IF
  IF ( ThPos>0) THEN
    CALL Allocate(thProfG)
    CALL ExchangeCell(thProfG)
    CALL Allocate(thProfG)
    CALL Allocate(TAbsCell,1,1)
  END IF  

  IF ( uPosL>0)  THEN
    CALL VectorInit(uPosL,VecMet,UStart,Time)
    CALL Mult(RhoCell,VecMet,uPosL)
    CALL Allocate(uCell)
  END IF
  IF ( uPosR>0)  THEN
    CALL VectorInit(uPosR,VecMet,UStart,Time)
    CALL Mult(RhoCell,VecMet,uPosR)
  END IF
  IF ( vPosL>0) THEN
    CALL VectorInit(vPosL,VecMet,vStart,Time)
    CALL Mult(RhoCell,VecMet,vPosL)
    CALL Allocate(vCell)
  END IF
  IF ( vPosR>0) THEN
    CALL VectorInit(vPosR,VecMet,vStart,Time)
    CALL Mult(RhoCell,VecMet,vPosR)
  END IF
  IF ( wPosL>0) THEN
    CALL VectorInit(wPosL,VecMet,wStart,Time)
    CALL Mult(RhoCell,VecMet,wPosL)
    CALL Allocate(wCell)
  END IF
  IF ( wPosR>0) THEN
    CALL VectorInit(wPosR,VecMet,wStart,Time)
    CALL Mult(RhoCell,VecMet,wPosR)
  END IF
  IF ( ThPos>0) THEN
    CALL VectorInit(ThPos,VecMet,ThStart,Time)
    CALL Mult(RhoCell,VecMet,ThPos)
  END IF
    
  CALL ExchangeCell(VecMet)
  IF (uPosL*uPosR>0) THEN
    CALL VelocityCellToFaceLR(VecMet,VelF1,VelF1,Time)
    CALL BoundaryVelocity(VelF1,Time)
  END IF

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib)) 
    DO i=1,NumBoundCell
      CALL SetBoundCells(BoundCell(i))
    END DO
  END DO
  WRITE(*,*) 'Pres 4'

  ! -- Project VelF1
! CALL PrepareF(VecMet,VelF1,Time)
  dtM = dtP
  CALL IniJacAccGrav


!--------------Begin Block-----------------
! CALL start_MPI(subdom%myrank)
  MultiGridMicha = .TRUE.
  IF (MultiGridMicha) THEN
  subdom%myrank = MyId

  CALL read_blocks_from_netcdf(file_path, subdom)

  CALL define_coarse_grids()
  CALL partition_clevels()
  CALL determine_color_startcell()

  CALL intialize_coarse_subdomains(subdom)
  IF (subdom%myrank .EQ. 0) WRITE(*,*) "Blocks read"

  this_subdomain => subdomain_lev
  CALL this_subdomain%set_cell_type()
  CALL this_subdomain%set_face_type()

  !make the communication interface
  CALL this_subdomain%init_boundcomm_interface(blockranks)

  CALL make_operators_static(this_subdomain, btype, isblackijkzero_lev(1)%data)

  WRITE(*,*) 'Pres 5'
  ALLOCATE(rhoM(this_subdomain%ncells))
  ALLOCATE(rhothetaMi(this_subdomain%ncells))
  ALLOCATE(rhovelM(this_subdomain%nfaces))
  ALLOCATE(rhoM_New(this_subdomain%ncells))
  ALLOCATE(rhothetaMi_New(this_subdomain%ncells))
  ALLOCATE(rhovelM_New(this_subdomain%nfaces))

  icellM = 1
  DO n = 1, this_subdomain%nblocks
    IF (this_subdomain%blockiscomp(n)) THEN
      DO i = 1, this_subdomain%blocks(n)%fld_shape(3)
        DO j = 1, this_subdomain%blocks(n)%fld_shape(2)
          DO k = 1, this_subdomain%blocks(n)%fld_shape(1)
            IF (this_subdomain%blocks(n)%volseff(k, j, i) .LT. 1.001e-20) THEN
              rhoM(icellM) = 0
              rhothetaMi(icellM) = 0
            ELSE
              rhoM(icellM) = rho0M
              rhothetaMi(icellM) = theta0M * rho0M
          END IF
            icellM = icellM + 1
          END DO
        END DO
      END DO
    END IF
  END DO
  END IF

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    VelF1(ibLoc) = 0.0d0
    VelF1(ibLoc)%uF(ix0+1,iy0+1,iz0+1) = 1.0
    VelF1(ibLoc)%uF(ix0+2,iy0+1,iz0+1) = 1.0
    VecMet(ibLoc)%Vec(ThPos)%c(:,:,:,1) = theta0M*VolC/(VolC+1.d-40) 
    VecMet(ibLoc)%Vec(RhoPos)%c(:,:,:,1) = Rho0M*VolC/(VolC+1.d-40) 
  END DO  

  IF (MultiGridMicha) THEN
  CALL CopyVec2VecC(subdom,rhothetaMi,VecMet,ThPos)
  CALL CopyVec2VecC(subdom,rhoM,VecMet,RhoPos)
  CALL CopyVelF2Vec(subdom,VelF1,rhovelM)

  GravComp=0.0d0
  CALL CopyVelF2Vec(subdom,VelF1,rhovelM)
  CALL CopyVecC2Vec(subdom,VecMet,ThPos,rhothetaMi)
  CALL CopyVecC2Vec(subdom,VecMet,RhoPos,rhoM)
  CALL setup_pressure_solver(this_subdomain, dtM, gamm, rhothetaMi, rhoM)
  CALL solveM(this_subdomain, rhothetaMi_new, rhoM_new, rhovelM_new, rhothetaMi, rhoM, rhovelM, dtM, gamm, NumIter, resmax_tol)
  WRITE(*,*) MyId,'rhovelM_New',SUM(ABS(rhovelM_New))
  CALL CopyVec2VelF(subdom,rhovelM_New,VelF2)
  END IF

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
!   VecMet(ibLoc)%Vec(RhoPos)%c(:,:,:,1)=VecMet(ibLoc)%Vec(RhoPos)%c(:,:,:,1)*VolC/(VolC+1.0d-40)
!   VecMet(ibLoc)%Vec(ThPos)%c(:,:,:,1)=VecMet(ibLoc)%Vec(ThPos)%c(:,:,:,1)*VolC/(VolC+1.0d-40)
    Rho=>VecMet(ibLoc)%Vec(RhoPos)%c
    RhoV=>VecMet(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecMet(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecMet(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VecMet(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VecMet(ibLoc)%Vec(RhoSPos)%c
    Th=>VecMet(ibLoc)%Vec(ThPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    CALL AbsTCompute
  END DO  
  CALL ExchangeCell(VecMet)
  CALL JacAccGrav(VecMet)
  WRITE(*,*) 'ProjectVelface ',dtP
  WRITE(*,*) 'VelF ',Dot(VelF1,VelF1)
  CALL ProjectVelface(dtP,VelF1,VecMet,VecMet)
  WRITE(*,*) 'nach ProjectVelface '
! CALL CopyVelF2Vec(subdom,VelF1,rhovelM_New)
! WRITE(*,*) MyId,'rhovelM_New',SUM(ABS(rhovelM_New))

! Output Profile
! IF (ProfOut) THEN
!   CALL VelocityFaceToCellLR(VelF1,VecMet)
!   CALL OutputProfile(VecMet,VecMet,ProfileEND)
! END IF
  CALL MPI_Finalize(MPIErr)

CONTAINS

SUBROUTINE MaxVelF(VelF)
  TYPE(VelocityFace_T), POINTER :: VelF(:)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1-1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (ABS(VelF(ibLoc)%uF(ix,iy,iz))>15.0d0) THEN
            WRITE(*,*) ix,iy,iz,'VelF',VelF(ibLoc)%uF(ix,iy,iz)
            WRITE(*,*) 'Zelle links'
            WRITE(*,*) VolC(ix,iy,iz)
            WRITE(*,*) FU(ix-1,iy,iz),FU(ix,iy,iz)
            WRITE(*,*) FV(ix,iy-1,iz),FV(ix,iy,iz)
            WRITE(*,*) FW(ix,iy,iz-1),FW(ix,iy,iz)
            WRITE(*,*) 'Zelle rechts'
            WRITE(*,*) VolC(ix+1,iy,iz)
            WRITE(*,*) FU(ix+1-1,iy,iz),FU(ix+1,iy,iz)
            WRITE(*,*) FV(ix+1,iy-1,iz),FV(ix+1,iy,iz)
            WRITE(*,*) FW(ix+1,iy,iz-1),FW(ix+1,iy,iz)
          END IF  
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE

FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S,N=1.d-2,th0=300.0d0

  S=N*N/Grav
  ThetaLoc=th0*EXP(S*z)

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S,N=1.0d-2,th0=300.0d0

  S=N*N/Grav
  IF (N>Zero) THEN
    PressLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    PressLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF


END FUNCTION PressLoc

SUBROUTINE Hydrostatic

  INTEGER, PARAMETER :: nzz=10
  REAL(RealKind), PARAMETER :: H=10.0d3
  INTEGER :: i
  REAL(RealKind) :: pLoc 
  REAL(RealKind) :: zz(nzz),ThVert(nzz),RhoVert(nzz)
  REAL(RealKind) :: dzz

  dzz=H/nzz
  zz(1)=0.5d0*dzz
  DO i=2,nzz
    zz(i)=zz(i-1)+dzz
  END DO
  DO i=1,nzz
    ThVert(i)=ThetaLoc(zz(i))
  END DO
  pLoc=PressLoc(zz(1))
  RhoVert(1)=pLoc/((pLoc/p0)**kappa*Rd*ThVert(1))
  DO i=2,nzz
    Call RhoHydro(Rhovert(i),ThVert(i),zz(i) &
                 ,Rhovert(i-1),ThVert(i-1),zz(i-1))
  END DO
  OPEN(UNIT=10,FILE='Profile',STATUS='UNKNOWN')
  WRITE(10,*) 'RhoProf'
  WRITE(10,*) nzz,2 
  WRITE(10,*) 'Equal' 
  DO i=1,nzz
    WRITE(10,*) zz(i),Rhovert(i)
  END DO
  WRITE(10,*) 'thProf'
  WRITE(10,*) nzz,2 
  WRITE(10,*) 'Equal' 
  DO i=1,nzz
    WRITE(10,*) zz(i),Thvert(i)
  END DO
  CLOSE(10)
  
END SUBROUTINE Hydrostatic

SUBROUTINE RhoHydro(Rho,Th,z,Rho1,Th1,z1)

  REAL(RealKind) :: Rho,Th,z,Rho1,Th1,z1

  REAL(RealKind) :: p1,p
  REAL(RealKind) :: RhoNew,RhoL,RhoR
  REAL(RealKind) :: Rhs
  REAL(RealKind) :: F,FL,FR,DF
  REAL(RealKind) :: TolErr

  IF (RealKind==8) THEN
    TolErr=1.e-9_RealKind 
  ELSE  
    TolErr=1.e-4_RealKind 
  END IF  
  p1=PresLoc(Rho1,Th1)
  Rhs=-p1+0.5d0*(z-z1)*Grav*Rho1
  RhoL=0.0d0
  p=PresLoc(RhoL,Th)
  FL=p+0.5d0*(z-z1)*Grav*RhoL+Rhs
  RhoR=Rho1
  p=PresLoc(RhoR,Th)
  FR=p+0.5d0*(z-z1)*Grav*RhoR+Rhs
  F=FR
  Rho=RhoR
  DO 
    DF=DPreDRho(Rho,Th)+0.5d0*(z-z1)*Grav
    RhoNew=Rho-F/DF
    IF (RhoNew<RhoL.OR.RhoNew>RhoR) THEN
      RhoNew=0.5d0*(RhoL+RhoR)
    END IF
    p=PresLoc(RhoNew,Th)
    F=p+0.5d0*(z-z1)*Grav*RhoNew+Rhs
    IF (F<0.0d0) THEN
      FL=F
      RhoL=RhoNew
    ELSE
      FR=F
      RhoR=RhoNew
    END IF
    IF (ABS(F)<=TolErr.OR.ABS(Rho-RhoNew)<=TolErr*Rho) THEN
      EXIT
    ELSE
      Rho=RhoNew
    END IF
  END DO

END SUBROUTINE RhoHydro

FUNCTION PresLoc(RhoLoc,ThLoc)
  REAL(RealKind) :: PresLoc
  REAL(RealKind) :: RhoLoc,ThLoc
  PresLoc=p0*(Rd*RhoLoc*ThLoc/p0)**(One/(One-kappa))
END FUNCTION PresLoc

FUNCTION DPreDRho(RhoLoc,ThLoc)
  REAL(RealKind) :: DPreDrho
  REAL(RealKind) :: RhoLoc,ThLoc
  DPreDrho=p0*(One/(One-kappa))*(Rd*RhoLoc*ThLoc/p0)**(kappa/(One-kappa)) &
           *Rd*ThLoc/p0
END FUNCTION DPreDrho

SUBROUTINE UpdateBoundary(Time)  

  REAL(RealKind) :: Time
  INTEGER, SAVE :: PosTimeEnvi=3
  IF (Time>Time2) THEN
    Time1=Time2
    VecEnv2=Zero
    IF ( uPosEnv>0) THEN
      CALL VectorInit(uPosEnv,VecEnv1,UStartE,Time1)
    END IF
    IF ( vPosEnv>0) THEN
      CALL VectorInit(vPosEnv,VecEnv1,vStartE,Time1)
    END IF
    IF ( wPosEnv>0) THEN
      CALL VectorInit(wPosEnv,VecEnv1,wStart,Time1)
    END IF
    IF (ThPosEnv>0) THEN
      CALL VectorInit(ThPosEnv,VecEnv1,ThStart,Time1)
    END IF
    IF (tkePosEnv>0) THEN
      CALL VectorInit(tkePosEnv,VecEnv1,tkeStart,Time1)
    END IF
    IF (disPosEnv>0) THEN
      CALL VectorInit(disPosEnv,VecEnv1,disStart,Time1)
    END IF
    Time2=TimeEnvi(PosTimeEnvi)
    PosTimeEnvi = PosTimeEnvi + 1
    VecEnv2=Zero
    IF ( uPosEnv>0) THEN
      CALL VectorInit(uPosEnv,VecEnv2,UStartE,Time2)
    END IF
    IF ( vPosEnv>0) THEN
      CALL VectorInit(vPosEnv,VecEnv2,vStartE,Time2)
    END IF
    IF ( wPosEnv>0) THEN
      CALL VectorInit(wPosEnv,VecEnv2,wStart,Time2)
    END IF
    IF (ThPosEnv>0) THEN
      CALL VectorInit(ThPosEnv,VecEnv2,ThStart,Time2)
    END IF
    IF (tkePosEnv>0) THEN
      CALL VectorInit(tkePosEnv,VecEnv2,tkeStart,Time2)
    END IF
    IF (disPosEnv>0) THEN
      CALL VectorInit(disPosEnv,VecEnv2,disStart,Time2)
    END IF
  END IF
END SUBROUTINE UpdateBoundary

SUBROUTINE CopyVelF2Vec(subdom,VelF,x)
  TYPE(SubDomain), INTENT(in) :: subdom
  TYPE(VelocityFace_T) :: VelF(:)
  REAL(Realkind), INTENT(inout) :: x(:)

  INTEGER :: iz, iy, ix
  INTEGER :: iz1, iy1, ix1
  INTEGER :: iz0, iy0, ix0
  INTEGER :: current_ind

  current_ind = 1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ix1=UBOUND(VelF(ibLoc)%uF,1)
    ix0=LBOUND(VelF(ibLoc)%uF,1)
    iy1=UBOUND(VelF(ibLoc)%vF,2)
    iy0=LBOUND(VelF(ibLoc)%vF,2)
    iz1=UBOUND(VelF(ibLoc)%wF,3)
    iz0=LBOUND(VelF(ibLoc)%wF,3)
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          x(current_ind)=VelF(ibLoc)%uF(ix,iy,iz)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          x(current_ind)=VelF(ibLoc)%vF(ix,iy,iz)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          x(current_ind)=VelF(ibLoc)%wF(ix,iy,iz)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CopyVelF2Vec

SUBROUTINE CopyVec2VelF(subdom,x,VelF)
  TYPE(SubDomain), INTENT(in) :: subdom
  REAL(Realkind), INTENT(inout) :: x(:)
  TYPE(VelocityFace_T) :: VelF(:)

  INTEGER :: iz, iy, ix
  INTEGER :: iz1, iy1, ix1
  INTEGER :: iz0, iy0, ix0
  INTEGER :: current_ind

  current_ind = 1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ix1=UBOUND(VelF(ibLoc)%uF,1)
    ix0=LBOUND(VelF(ibLoc)%uF,1)
    iy1=UBOUND(VelF(ibLoc)%vF,2)
    iy0=LBOUND(VelF(ibLoc)%vF,2)
    iz1=UBOUND(VelF(ibLoc)%wF,3)
    iz0=LBOUND(VelF(ibLoc)%wF,3)
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          VelF(ibLoc)%uF(ix,iy,iz)=x(current_ind)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          VelF(ibLoc)%vF(ix,iy,iz)=x(current_ind)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          VelF(ibLoc)%wF(ix,iy,iz)=x(current_ind)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CopyVec2VelF

SUBROUTINE CopyVecC2Vec(subdom,Vec,Pos,x)
  TYPE(SubDomain), INTENT(in) :: subdom
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: Pos
  REAL(Realkind), INTENT(inout) :: x(:)

  INTEGER :: iz, iy, ix
  INTEGER :: current_ind

  current_ind = 1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          x(current_ind)=Vec(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE CopyVecC2Vec

SUBROUTINE CopyVec2VecC(subdom,x,Vec,Pos)
  TYPE(SubDomain), INTENT(in) :: subdom
  REAL(Realkind), INTENT(inout) :: x(:)
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: Pos

  INTEGER :: iz, iy, ix
  INTEGER :: current_ind

  current_ind = 1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          Vec(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)=x(current_ind)
          current_ind=current_ind+1
        END DO
      END DO
    END DO
  END DO
  
END SUBROUTINE CopyVec2VecC

END PROGRAM MainProg
