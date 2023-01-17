MODULE JacAccGrav_Mod

  USE Control_Mod
  USE DataType_Mod
  USE Sp_Mod
  USE TriT_3D_Mod
  USE TriTB_3D_Mod
  USE TriTR_3D_Mod
  USE Reverse_Mod
  USE SparseSolver_Mod
  USE ReadWeights_Mod, ONLY : VolFace
  USE Domain_Mod

  IMPLICIT NONE

  TYPE BlockMatrix
    TYPE(TriTDiag3D_T), POINTER :: MatrixT=>NULL()
    TYPE(TriTBDiag3D_T), POINTER :: MatrixTB=>NULL()
    TYPE(TriTRDiag3D_T), POINTER :: MatrixTR=>NULL()
    TYPE(SpRowCol), POINTER :: MatrixMu(:)=>NULL()
    TYPE(SpRowColMat), POINTER :: MatrixMuMat(:)=>NULL()
    TYPE (DMUMPS_STRUC), POINTER :: MatrixMumps=>NULL()
  END TYPE BlockMatrix

  TYPE(BlockMatrix), POINTER :: LaplNeumann(:)=>NULL()
  TYPE(BlockMatrix), POINTER :: LaplDirichlet(:)=>NULL()
  TYPE(BlockMatrix), POINTER :: MaActual(:)=>NULL()

  TYPE(PressureVelocity), POINTER :: xx(:)=>NULL()
  TYPE(PressureVelocity), POINTER :: bb(:)=>NULL()

  REAL(RealKind), PRIVATE, POINTER :: Theta(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThetaProf(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: p(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Sound(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: T(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: KinEn(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: E(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DTU(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DTV(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DTW(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DUT(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DUR(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DTT(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DUU(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DUV(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DUW(:,:,:)

  TYPE(TriTDiag3D_T),  PRIVATE, POINTER :: AR=>NULL()
  TYPE(TriTBDiag3D_T),  PRIVATE, POINTER :: ATB=>NULL()
  TYPE(TriTRDiag3D_T),  PRIVATE, POINTER :: ATR=>NULL()
  TYPE(SpRowCol),  PRIVATE, POINTER :: A(:)=>NULL()
  TYPE(SpRowColMat),  PRIVATE, POINTER :: AMat(:)=>NULL()
  TYPE (DMUMPS_STRUC), PRIVATE, POINTER :: AMumps
  INTEGER, PRIVATE :: iRef
  REAL(RealKind), PRIVATE, POINTER :: xS(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: bS(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: xSM(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: bSM(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: xS1(:)
  REAL(RealKind), PRIVATE, POINTER :: bS1(:)
  REAL(RealKind), PRIVATE, POINTER :: xSB(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: bSB(:,:,:,:)
  REAL(RealKind), PRIVATE, ALLOCATABLE, TARGET :: bFine(:,:,:)
  REAL(RealKind), PRIVATE, ALLOCATABLE, TARGET :: bFineB(:,:,:,:)
  INTEGER, PRIVATE :: InitJacAccGrav=0

  REAL(RealKind), POINTER :: LaplNeumGlob(:,:)=>NULL()

  INTERFACE MatVecP
    MODULE PROCEDURE MatVecP_TriTB,MatVecP_Tri,MatVecP_Mu,MatVecP_MuTR,MatVecP_TriTR
  END INTERFACE
  INTERFACE LaplCompute
    MODULE PROCEDURE LaplMuCompute,LaplTCompute,LaplTBCompute,LaplTRCompute,LaplMumpsCompute
  END INTERFACE

CONTAINS

SUBROUTINE CoarseNeumann
  
  
  REAL(RealKind) :: Temp
  INTEGER :: in
  INTEGER :: ibMax
  INTEGER :: i,j,info
  TYPE(Nachbar_T), POINTER :: Nachbar
  REAL(RealKind) :: Row(1:nb)
  REAL(RealKind) :: DiagLoc,Diag(1:nb)

  INTEGER :: ix,iy,iz

  
!   Set Boundary
  ibMax=1
  DO ib=1,nb
    CALL Set(Floor(ib))
    IF (igx0==domain%ix0.AND.BCP%West==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE IF (igx1==domain%ix1.AND.BCP%East==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE IF (igy0==domain%iy0.AND.BCP%South==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE IF (igy1==domain%iy1.AND.BCP%North==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE IF (igz0==domain%iz0.AND.BCP%Bottom==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE IF (igz1==domain%iz1.AND.BCP%Top==1) THEN
      Floor(ib)%Boundary=1.0e0
    ELSE
      IF (Anelastic.OR.PseudoIn) THEN
        Floor(ib)%Boundary=1.0e0
      ELSE
        Floor(ib)%Boundary=1.0d0
      END IF  
      ibMax=ib
    END IF
!   Floor(ib)%Boundary=0.0d0 !OSSI
  END DO
  IF (BCP%West/=1.AND.BCP%East/=1.AND. &
      BCP%South/=1.AND.BCP%North/=1.AND. &
      BCP%Bottom/=1.AND.BCP%Top/=1.AND.(Anelastic.OR.PseudoIn)) THEN
      Floor(ibMax)%Boundary=0.0e0
  END IF

  IF (.NOT.ASSOCIATED(LaplNeumGlob)) THEN
    ALLOCATE(LaplNeumGlob(nb,nb))
  END IF
  LaplNeumGlob=Zero
  DO ib=1,nb
    Row=Zero
    CALL Set(Floor(ib))
    IF (blMPI(ib)%Proc==MyId) THEN
      CALL Set(Floor(ib)) 
      ibLoc=blMPI(ib)%ibLoc
      DUU=>DUUG(ibLoc)%uF
      DUV=>DUUG(ibLoc)%vF
      DUW=>DUUG(ibLoc)%wF
      DUT=>DUTG(ibLoc)%p
      DTU=>DTUG(ibLoc)%uF
      DTV=>DTUG(ibLoc)%vF
      DTW=>DTUG(ibLoc)%wF
      DO in=1,AnzahlNachbar
        Nachbar=>Nachbars(in)
        CALL Set(Nachbar)

        IF (Nachbar%nType(2:2)=='w') THEN
          IF (Refine<RefineNachbar .OR. &
             (Refine==RefineNachbar .AND. ibn<ib)) THEN
            jx0=ix0
            Row(ibn)=Row(ibn)- & 
                        SUM(FU(jx0,jy0+1:jy1,jz0+1:jz1)*DTU(jx0,jy0+1:jy1,jz0+1:jz1)* &  
                            FUG(jx0,jy0+1:jy1,jz0+1:jz1)*DUU(jx0,jy0+1:jy1,jz0+1:jz1)/ &
                            (VolFace(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1)+Eps))
          ELSE IF(igx0==domain%ix0.AND.BCP%West==1) THEN
            jx0=ix0
            Row(ibLoc)=Row(ibLoc)+ & 
                        SUM(FU(jx0,jy0+1:jy1,jz0+1:jz1)*DTU(jx0,jy0+1:jy1,jz0+1:jz1)* &  
                            FUG(jx0,jy0+1:jy1,jz0+1:jz1)*DUU(jx0,jy0+1:jy1,jz0+1:jz1)/ &
                            (VolFace(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1)+Eps))
          END IF 
        END IF 

        IF (Nachbar%nType(2:2) == 'e') THEN
          IF (Refine<RefineNachbar .OR. &
             (Refine==RefineNachbar .AND. ibn<ib)) THEN
             jx1=ix1
             Row(ibn)=Row(ibn)- & 
                         SUM(FU(jx1,jy0+1:jy1,jz0+1:jz1)*DTU(jx1,jy0+1:jy1,jz0+1:jz1)* &
                             FUG(jx1,jy0+1:jy1,jz0+1:jz1)*DUU(jx1,jy0+1:jy1,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1)+Eps))
          ELSE IF(igx1==domain%ix1.AND.BCP%East==1) THEN
             jx1=ix1
             Row(ibLoc)=Row(ibLoc)+ & 
                         SUM(FU(jx1,jy0+1:jy1,jz0+1:jz1)*DTU(jx1,jy0+1:jy1,jz0+1:jz1)* &
                             FUG(jx1,jy0+1:jy1,jz0+1:jz1)*DUU(jx1,jy0+1:jy1,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1)+Eps))
          END IF
        END IF

        IF (Nachbar%nType(2:2) == 's') THEN
          IF (Refine<RefineNachbar .OR. &
             (Refine==RefineNachbar .AND. ibn<ib)) THEN
             jy0=iy0
             Row(ibn)=Row(ibn)- & 
                         SUM(FV(jx0+1:jx1,jy0,jz0+1:jz1)*DTV(jx0+1:jx1,jy0,jz0+1:jz1)* &
                             FVG(jx0+1:jx1,jy0,jz0+1:jz1)*DUV(jx0+1:jx1,jy0,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1)+Eps))
          ELSE IF (igy0==domain%iy0.AND.BCP%South==1) THEN
             jy0=iy0
             Row(ibLoc)=Row(ibLoc)+ & 
                         SUM(FV(jx0+1:jx1,jy0,jz0+1:jz1)*DTV(jx0+1:jx1,jy0,jz0+1:jz1)* &
                             FVG(jx0+1:jx1,jy0,jz0+1:jz1)*DUV(jx0+1:jx1,jy0,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1)+Eps))
          END IF
        END IF
       
        IF (Nachbar%nType(2:2) == 'n') THEN
          IF (Refine<RefineNachbar .OR. &
             (Refine==RefineNachbar .AND. ibn<ib)) THEN
             jy1=iy1
             Row(ibn)=Row(ibn)- & 
                         SUM(FV(jx0+1:jx1,jy1,jz0+1:jz1)*DTV(jx0+1:jx1,jy1,jz0+1:jz1)* &
                             FVG(jx0+1:jx1,jy1,jz0+1:jz1)*DUV(jx0+1:jx1,jy1,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1)+Eps))
          ELSE IF (igy1==domain%iy1.AND.BCP%North==1) THEN
             jy1=iy1
             Row(ibLoc)=Row(ibLoc)+ & 
                         SUM(FV(jx0+1:jx1,jy1,jz0+1:jz1)*DTV(jx0+1:jx1,jy1,jz0+1:jz1)* &
                             FVG(jx0+1:jx1,jy1,jz0+1:jz1)*DUV(jx0+1:jx1,jy1,jz0+1:jz1)/ &
                             (VolFace(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1)+Eps))
          END IF
        END IF
        
        IF (Nachbar%nType (2:2) == 'b') THEN
          IF (Refine<RefineNachbar .OR. &
            (Refine==RefineNachbar .AND. ibn<ib)) THEN
            jz0=iz0
            Row(ibn)=Row(ibn)- & 
                         SUM(FW(jx0+1:jx1,jy0+1:jy1,jz0)*DTW(jx0+1:jx1,jy0+1:jy1,jz0)* &
                             FWG(jx0+1:jx1,jy0+1:jy1,jz0)*DUW(jx0+1:jx1,jy0+1:jy1,jz0)/ &
                             (VolFace(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1)+Eps))
          ELSE IF (igz0==domain%iz0.AND.BCP%Bottom==1) THEN
            jz0=iz0
            Row(ibLoc)=Row(ibLoc)+ & 
                         SUM(FW(jx0+1:jx1,jy0+1:jy1,jz0)*DTW(jx0+1:jx1,jy0+1:jy1,jz0)* &
                             FWG(jx0+1:jx1,jy0+1:jy1,jz0)*DUW(jx0+1:jx1,jy0+1:jy1,jz0)/ &
                             (VolFace(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1)+Eps))
          END IF
        END IF
        IF (Nachbar%nType (2:2) == 't') THEN
          IF (Refine<RefineNachbar .OR. &
              (Refine==RefineNachbar .AND. ibn<ib)) THEN
             jz1=iz1
             Row(ibn)=Row(ibn)- & 
                         SUM(FW(jx0+1:jx1,jy0+1:jy1,jz1)*DTW(jx0+1:jx1,jy0+1:jy1,jz1)* &
                             FWG(jx0+1:jx1,jy0+1:jy1,jz1)*DUW(jx0+1:jx1,jy0+1:jy1,jz1)/ &
                             (VolFace(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1)+Eps))
          ELSE IF (igz1==domain%iz1.AND.BCP%Top==1) THEN
             jz1=iz1
             Row(ibLoc)=Row(ibLoc)+ & 
                         SUM(FW(jx0+1:jx1,jy0+1:jy1,jz1)*DTW(jx0+1:jx1,jy0+1:jy1,jz1)* &
                             FWG(jx0+1:jx1,jy0+1:jy1,jz1)*DUW(jx0+1:jx1,jy0+1:jy1,jz1)/ &
                             (VolFace(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1)+Eps))
          END IF
        END IF
      END DO ! in
      DiagLoc=SUM(VolB/(DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)+Eps))

    END IF
    CALL MPI_Bcast(Row(1:nb), nb, MPI_RealKind, blMPI(ib)%Proc,  &
                    MPI_COMM_WORLD, MPIErr) 
    LaplNeumGlob(1:nb,ib)=LaplNeumGlob(1:nb,ib)+Row(1:nb)
    CALL MPI_Bcast(DiagLoc, 1, MPI_RealKind, blMPI(ib)%Proc,  &
                    MPI_COMM_WORLD, MPIErr) 
    Diag(ib)=DiagLoc
  END DO       ! ib

! -- Matrix symmetrisieren --
  DO j=1,nb
    DO i=j+1,nb
      Temp = (LaplNeumGlob(i,j) + LaplNeumGlob(j,i))
      LaplNeumGlob(i,j) = Temp
      LaplNeumGlob(j,i) = Temp
    END DO 
    LaplNeumGlob(j,j)=2.0d0*LaplNeumGlob(j,j)-SUM(LaplNeumGlob(:,j))
  END DO
  DO j=1,nb
    DO i=1,nb
      LaplNeumGlob(i,j)=(beta0*dtP)**2*LaplNeumGlob(i,j)*Floor(i)%Boundary*Floor(j)%Boundary
    END DO
  END DO
! -- Faktorisieren von D(+)*D --
  DO ib=1,nb
    LaplNeumGlob(ib,ib)=LaplNeumGlob(ib,ib)+FacAnela*Diag(ibLoc)
  END DO
  CALL dpofa(LaplNeumGlob,nb,nb,info)
END SUBROUTINE CoarseNeumann

SUBROUTINE IniJacAccGrav

  CALL Allocate(DTUG)
  CALL Allocate(DUTG)
  CALL Allocate(DURG)
  CALL Allocate(DTTG)
  CALL Allocate(DUUG)
  DTUG=One
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DUUG(ibLoc)%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1)=FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)/(FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
    DUUG(ibLoc)%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1)=FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)/(FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)+Eps)
    DUUG(ibLoc)%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1)=FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)/(FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)+Eps)
  END DO  

END SUBROUTINE IniJacAccGrav

SUBROUTINE JacAccGrav(Vec)

  TYPE(Vector4Cell_T), POINTER :: Vec(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Theta=>Vec(ibLoc)%Vec(thPos)%c
    ThetaProf=>ThProfG(ibLoc)%c
    RhoV=>Vec(ibLoc)%Vec(RhoCPos)%c
    RhoL=>Vec(ibLoc)%Vec(RhoCPos)%c
    RhoR=>Vec(ibLoc)%Vec(RhoRPos)%c
    Rho=>RhoCell(ibLoc)%c
    IF (ASSOCIATED(TAbsCell)) THEN
      T=>TAbsCell(ibLoc)%Vec(1)%c
    END IF  
    IF (ASSOCIATED(PreCell)) THEN
      p=>PreCell(ibLoc)%c
    END IF  
    IF (ASSOCIATED(SoundCell)) THEN
      Sound=>SoundCell(ibLoc)%c
    END IF  
    IF (ASSOCIATED(KinEnCell)) THEN
      KinEn=>KinEnCell(ibLoc)%c
    END IF  
    IF (ASSOCIATED(ECell)) THEN
      E=>ECell(ibLoc)%c
    END IF  
    IF (ThetaKind=='Energy') THEN
      Theta=Theta+p
    END IF
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DTT=>DTTG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    ThetaProf=>ThProfG(ibLoc)%c
    CALL WeightDiagCompute 
    IF (ThetaKind=='Energy') THEN
      Theta=Theta-p
    END IF
  END DO
  CALL CoarseNeumann
  IF (.NOT.ASSOCIATED(LaplNeumann)) THEN
    ALLOCATE(LaplNeumann(nbGlob))
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      IF (MultiEx) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixMu)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixMu(1))
        END IF
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixMumps)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixMumps)
          NULLIFY(LaplNeumann(ibLoc)%MatrixMumps%IRN)
          NULLIFY(LaplNeumann(ibLoc)%MatrixMumps%JCN)
        END IF
      ELSE IF (MultiMu) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixMu)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixMu(RefLevel))
        END IF
      ELSE IF (MultiMuTR) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixMuMat)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixMuMat(RefLevel))
        END IF
      ELSE IF (MultiTriT) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixT)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixT)
        END IF
      ELSE IF (MultiTriTB) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixTB)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixTB)
        END IF
      ELSE IF (MultiTriTR) THEN
        IF (.NOT.ASSOCIATED(LaplNeumann(ibLoc)%MatrixTR)) THEN
          ALLOCATE(LaplNeumann(ibLoc)%MatrixTR)
        END IF
      END IF 
    END DO
  END IF 
  IF (.NOT.ASSOCIATED(LaplDirichlet)) THEN
    ALLOCATE(LaplDirichlet(nbGlob))
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      IF (MultiEx) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMu)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixMu(1))
        END IF
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMumps)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixMumps)
          NULLIFY(LaplDirichlet(ibLoc)%MatrixMumps%IRN)
          NULLIFY(LaplDirichlet(ibLoc)%MatrixMumps%JCN)
        END IF
      ELSE IF (MultiMu) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMu)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixMu(RefLevel))
        END IF
      ELSE IF (MultiMuTR) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMuMat)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixMuMat(RefLevel))
        END IF
      ELSE IF (MultiTriT) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixT)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixT)
        END IF
      ELSE IF (MultiTriTB) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixTB)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixTB)
        END IF
      ELSE IF (MultiTriTR) THEN
        IF (.NOT.ASSOCIATED(LaplDirichlet(ibLoc)%MatrixTR)) THEN
          ALLOCATE(LaplDirichlet(ibLoc)%MatrixTR)
        END IF
      END IF 
    END DO
  END IF 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL SetFactorN
    DTU=>DTUG(ibLoc)%uF 
    DTV=>DTUG(ibLoc)%vF 
    DTW=>DTUG(ibLoc)%wF 
    DUT=>DUTG(ibLoc)%p 
    DTT=>DTTG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixMumps)) THEN
      CALL LaplCompute(LaplNeumann(ibLoc)%MatrixMumps,LaplNeumann(ibLoc)%MatrixMu(1))
    ELSE IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixMu)) THEN
      CALL LaplCompute(LaplNeumann(ibLoc)%MatrixMu)
    ELSE IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixMuMat)) THEN
      DUR=>DURG(ibLoc)%p
      CALL LaplComputeMuTR(LaplNeumann(ibLoc)%MatrixMuMat)
    ELSE IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixT)) THEN
      CALL LaplCompute(LaplNeumann(ibLoc)%MatrixT)
    ELSE IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixTB)) THEN
      CALL LaplCompute(LaplNeumann(ibLoc)%MatrixTB)
    ELSE IF (ASSOCIATED(LaplNeumann(ibLoc)%MatrixTR)) THEN
      DUR=>DURG(ibLoc)%p
      CALL LaplCompute(LaplNeumann(ibLoc)%MatrixTR)
    END IF

    CALL SetFactorD
    IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMumps)) THEN
      CALL LaplCompute(LaplDirichlet(ibLoc)%MatrixMumps,LaplDirichlet(ibLoc)%MatrixMu(1))
    ELSE IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMu)) THEN
      CALL LaplCompute(LaplDirichlet(ibLoc)%MatrixMu)
    ELSE IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixMuMat)) THEN
      DUR=>DURG(ibLoc)%p
      CALL LaplComputeMuTR(LaplDirichlet(ibLoc)%MatrixMuMat)
    ELSE IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixT)) THEN
      CALL LaplCompute(LaplDirichlet(ibLoc)%MatrixT)
    ELSE IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixTB)) THEN
      CALL LaplCompute(LaplDirichlet(ibLoc)%MatrixTB)
    ELSE IF (ASSOCIATED(LaplDirichlet(ibLoc)%MatrixTR)) THEN
      DUR=>DURG(ibLoc)%p
      CALL LaplCompute(LaplDirichlet(ibLoc)%MatrixTR)
    END IF
  END DO

END SUBROUTINE JacAccGrav

SUBROUTINE MatVecP_TriTB(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(TriTBDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:,:)
  INTEGER :: nx,ny,nz
  
  nx=SIZE(x,2)
  ny=SIZE(x,3)
  nz=SIZE(x,4)
  CALL ReverseIndices(nx,ny,nz)
  CALL MatVecTBPReverse(y,A,x)
  
END SUBROUTINE MatVecP_TriTB

SUBROUTINE MatVecP_TriTR(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(TriTRDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:,:)
  INTEGER :: nx,ny,nz
  
  nx=SIZE(x,2)
  ny=SIZE(x,3)
  nz=SIZE(x,4)
  CALL ReverseIndices(nx,ny,nz)
  CALL MatVecTRPReverse(y,A,x)
  
END SUBROUTINE MatVecP_TriTR

SUBROUTINE MatVecP_Tri(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:)
  INTEGER :: nx,ny,nz

  nx=SIZE(x,1)
  ny=SIZE(x,2)
  nz=SIZE(x,3)
  CALL ReverseIndices(nx,ny,nz)
  CALL MatVecTPReverse(y,A,x)

END SUBROUTINE MatVecP_Tri

SUBROUTINE MatVecP_Mu(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(SpRowCol), POINTER :: A(:)
  REAL(RealKind) :: x(:,:,:)
  INTEGER :: nx,ny,nz

  nx=SIZE(x,1)
  ny=SIZE(x,2)
  nz=SIZE(x,3)
  CALL ReverseIndices(nx,ny,nz)
  CALL MatVecPMuReverse(y,A,x)

END SUBROUTINE  MatVecP_Mu

SUBROUTINE MatVecPMuReverse(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(SpRowCol), POINTER :: A(:)
  REAL(RealKind) :: x(:,:,:)
  REAL(RealKind) :: xRev(n1,n2,n3)
  REAL(RealKind) :: yRev(n1,n2,n3)
  INTEGER :: nx,ny,nz

  CALL Reverse(xRev,x)
  yRev=0.0d0
  CALL SpAVec(yRev(:,1,1),A(1),xRev(:,1,1))
  CALL ReverseBack(y,yRev)

END SUBROUTINE  MatVecPMuReverse

SUBROUTINE MatVecP_MuTR(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(SpRowColMat), POINTER :: A(:)
  REAL(RealKind) :: x(:,:,:,:)
  INTEGER :: nx,ny,nz

  CALL SpAVec(y(:,:,:,:),A(1),x(:,:,:,:))

END SUBROUTINE  MatVecP_MuTR


SUBROUTINE MatVecTPReverse(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:)
  REAL(RealKind) :: yRev(n1,n2,n3)
  REAL(RealKind) :: xRev(n1,n2,n3)

  CALL Reverse(xRev,x)
  yRev=Zero
  CALL MatVecT(yRev,A,xRev)
  CALL ReverseBack(y,yRev)

END SUBROUTINE MatVecTPReverse

SUBROUTINE MatVecTBPReverse(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(TriTBDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:,:)
  REAL(RealKind) :: yRev(2,n1,n2,n3)
  REAL(RealKind) :: xRev(2,n1,n2,n3)

  CALL ReverseB(xRev,x)
  yRev=Zero
  CALL MatVecTB(yRev,A,xRev)
  CALL ReverseBBack(y,yRev)

END SUBROUTINE MatVecTBPReverse

SUBROUTINE MatVecTRPReverse(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(TriTRDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:,:)
  REAL(RealKind) :: yRev(2,n1,n2,n3)
  REAL(RealKind) :: xRev(2,n1,n2,n3)


  CALL ReverseB(xRev,x)
  yRev=Zero
  CALL MatVecTR(yRev,A,xRev)
  CALL ReverseBBack(y,yRev)

END SUBROUTINE MatVecTRPReverse


SUBROUTINE WeightCoarse(n,m,l,FU,FV,FW,FUT,FVT,FWT,VC,VCP, &
                        nC,mC,lC,FUC,FVC,FWC,FUCT,FVCT,FWCT,VCC,VCCP, &
                        iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: FU(0:n,1:m,1:l),FV(1:n,0:m,1:l),FW(1:n,1:m,0:l)
  REAL(RealKind) :: FUT(0:n,1:m,1:l),FVT(1:n,0:m,1:l),FWT(1:n,1:m,0:l)
  REAL(RealKind) :: VC(n,m,l)
  REAL(RealKind) :: VCP(n,m,l)
  REAL(RealKind) :: FUC(0:nC,1:mC,1:lC),FVC(1:nC,0:mC,1:lC),FWC(1:nC,1:mC,0:lC)
  REAL(RealKind) :: FUCT(0:nC,1:mC,1:lC),FVCT(1:nC,0:mC,1:lC),FWCT(1:nC,1:mC,0:lC)
  REAL(RealKind) :: VCC(nC,mC,lC)
  REAL(RealKind) :: VCCP(nC,mC,lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

! Coarsening U-Faces
  FUC=0.0e0
  FUCT=0.0e0
  DO iC=0,nC
    i=MIN((iShift+1)*iC,n)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        FUC(iC,jC,kC)=FUC(iC,jC,kC)+FU(i,j,k)
        FUCT(iC,jC,kC)=FUCT(iC,jC,kC)+FUT(i,j,k)
      END DO
    END DO 
  END DO

! Coarsening V-Faces
  FVC=0.0e0
  FVCT=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO jC=0,mC
      j=MIN((jShift+1)*jC,m)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        FVC(iC,jC,kC)=FVC(iC,jC,kC)+FV(i,j,k)
        FVCT(iC,jC,kC)=FVCT(iC,jC,kC)+FVT(i,j,k)
      END DO
    END DO 
  END DO
! Coarsening W-Faces
  FWC=0.0e0
  FWCT=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO kC=0,lC
        k=MIN((kShift+1)*kC,l)
        FWC(iC,jC,kC)=FWC(iC,jC,kC)+FW(i,j,k)
        FWCT(iC,jC,kC)=FWCT(iC,jC,kC)+FWT(i,j,k)
      END DO
    END DO
  END DO

  VCC=0.0e0
  VCCP=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        VCC(iC,jC,kC)=VCC(iC,jC,kC)+VC(i,j,k)
        VCCP(iC,jC,kC)=VCCP(iC,jC,kC)+VCP(i,j,k)
      END DO 
    END DO 
  END DO 
END SUBROUTINE WeightCoarse

SUBROUTINE VolCoarse(n,m,l,VC, &
                     nC,mC,lC,VCC, &
                     iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: VC(n,m,l)
  REAL(RealKind) :: VCC(nC,mC,lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

  VCC=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        VCC(iC,jC,kC)=VCC(iC,jC,kC)+VC(i,j,k)
      END DO 
    END DO 
  END DO 
END SUBROUTINE VolCoarse

SUBROUTINE IncrCoarse(n,m,l,d1,d2,d3, &
                     nC,mC,lC,d1C,d2C,d3C, &
                     iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: d1(n),d2(m),d3(l)
  REAL(RealKind) :: d1C(nC),d2C(mC),d3C(lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

  d1C=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    d1C(iC)=d1C(iC)+d1(i)
  END DO 
  d2C=0.0e0
  DO j=1,m
    jc=(j+jShift)/(jShift+1)
    d2C(jC)=d2C(jC)+d2(j)
  END DO 
  d3C=0.0e0
  DO k=1,l
    kc=(k+kShift)/(kShift+1)
    d3C(kC)=d3C(kC)+d3(k)
  END DO 
END SUBROUTINE IncrCoarse

SUBROUTINE CIncrCoarse(n,C,d, &
                     nC,CC,dC,&
                     iShift)
  INTEGER :: n
  INTEGER :: nC
  REAL(RealKind) :: C(n),d(n)
  REAL(RealKind) :: CC(nC),dC(nC)
  INTEGER :: iShift

  INTEGER :: i,j,k,iC,jC,kC

  CC=0.0d0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    CC(iC)=CC(iC)+d(i)*C(i)
  END DO 
  CC=CC/dC
END SUBROUTINE CIncrCoarse

SUBROUTINE CVolCoarse(n,m,l,C,VC, &
                   nC,mC,lC,CC,VCC, &
                   iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: C(n,m,l),VC(n,m,l)
  REAL(RealKind) :: CC(nC,mC,lC),VCC(nC,mC,lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

  CC=0.0d0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        CC(iC,jC,kC)=CC(iC,jC,kC)+C(i,j,k)*VC(i,j,k)
      END DO 
    END DO 
  END DO 
  CC=CC/(VCC+Eps)
END SUBROUTINE CVolCoarse

SUBROUTINE FaceCoarse(n,m,l,FU,FV,FW, &
                        nC,mC,lC,FUC,FVC,FWC, &
                        iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: FU(0:n,1:m,1:l),FV(1:n,0:m,1:l),FW(1:n,1:m,0:l)
  REAL(RealKind) :: FUC(0:nC,1:mC,1:lC),FVC(1:nC,0:mC,1:lC),FWC(1:nC,1:mC,0:lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

! Coarsening U-Faces
  FUC=0.0e0
  DO iC=0,nC
    i=MIN((iShift+1)*iC,n)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        FUC(iC,jC,kC)=FUC(iC,jC,kC)+FU(i,j,k)
      END DO
    END DO 
  END DO

! Coarsening V-Faces
  FVC=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO jC=0,mC
      j=MIN((jShift+1)*jC,m)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        FVC(iC,jC,kC)=FVC(iC,jC,kC)+FV(i,j,k)
      END DO
    END DO 
  END DO
! Coarsening W-Faces
  FWC=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO kC=0,lC
        k=MIN((kShift+1)*kC,l)
        FWC(iC,jC,kC)=FWC(iC,jC,kC)+FW(i,j,k)
      END DO
    END DO
  END DO

END SUBROUTINE FaceCoarse

SUBROUTINE CFaceCoarse(n,m,l,CFU,CFV,CFW,FU,FV,FW, &
                       nC,mC,lC,CFUC,CFVC,CFWC,FUC,FVC,FWC, &
                       iShift,jShift,kShift)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: CFU(0:n,1:m,1:l),CFV(1:n,0:m,1:l),CFW(1:n,1:m,0:l)
  REAL(RealKind) :: FU(0:n,1:m,1:l),FV(1:n,0:m,1:l),FW(1:n,1:m,0:l)
  REAL(RealKind) :: CFUC(0:nC,1:mC,1:lC),CFVC(1:nC,0:mC,1:lC),CFWC(1:nC,1:mC,0:lC)
  REAL(RealKind) :: FUC(0:nC,1:mC,1:lC),FVC(1:nC,0:mC,1:lC),FWC(1:nC,1:mC,0:lC)
  INTEGER :: iShift,jShift,kShift

  INTEGER :: i,j,k,iC,jC,kC

! Coarsening U-Faces
  CFUC=0.0e0
  DO iC=0,nC
    i=MIN((iShift+1)*iC,n)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        CFUC(iC,jC,kC)=CFUC(iC,jC,kC)+CFU(i,j,k)*FU(i,j,k)
      END DO
    END DO 
  END DO
  CFUC=CFUC/(FUC+Eps)

! Coarsening V-Faces
  CFVC=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO jC=0,mC
      j=MIN((jShift+1)*jC,m)
      DO k=1,l
        kc=(k+kShift)/(kShift+1)
        CFVC(iC,jC,kC)=CFVC(iC,jC,kC)+CFV(i,j,k)*FV(i,j,k)
      END DO
    END DO 
  END DO
  CFVC=CFVC/(FVC+Eps)
! Coarsening W-Faces
  CFWC=0.0e0
  DO i=1,n
    ic=(i+iShift)/(iShift+1)
    DO j=1,m
      jc=(j+jShift)/(jShift+1)
      DO kC=0,lC
        k=MIN((kShift+1)*kC,l)
        CFWC(iC,jC,kC)=CFWC(iC,jC,kC)+CFW(i,j,k)*FW(i,j,k)
      END DO
    END DO
  END DO
  CFWC=CFWC/(FWC+Eps)

END SUBROUTINE CFaceCoarse

SUBROUTINE LaplTCompute(A)
  TYPE(TriTDiag3D_T), POINTER :: A
  IF (GradFull) THEN
    CALL LaplTFullCompute(A)
  ELSE
    CALL LaplTDualCompute(A)
  END IF
END SUBROUTINE LaplTCompute

SUBROUTINE LaplTDualCompute(A)

  TYPE(TriTDiag3D_T), POINTER :: A

  INTEGER :: k
  INTEGER :: n1C,n2C,n3C
  INTEGER :: Shift1,Shift2,Shift3
  INTEGER :: Level1,Level2,Level3

  TYPE(TriTDiag3D_T), POINTER :: Neu3
  TYPE(TriTDiag3D_T), POINTER :: Current3
  TYPE(TriTDiag3D_T), POINTER :: Neu2
  TYPE(TriTDiag3D_T), POINTER :: Current2

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)

  CALL ReverseIndices(nx,ny,nz)
  CALL SetCoarse
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  CALL Allocate(A)
  Current3=>A
  DO Level3=0,RefLevel-1
    Shift3=Coarse3**Level3-1
    Shift2=Coarse2**Level3-1
    Shift1=Coarse1**Level3-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL Allocate(Current3,n3C,n2C,n1C)
    CALL Compute(Current3,F1C,F2C,F3C,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
    DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
    Current2=>Current3
    ALLOCATE(Neu2)
    CALL Allocate(Neu2,n3C)
    DO Level2=Level3+1,RefLevel-1
      Shift2=Coarse2**Level2-1
      Shift1=Coarse1**Level2-1
      n1C=(n1+Shift1)/(Shift1+1)
      n2C=(n2+Shift2)/(Shift2+1)
      ALLOCATE(F1C(0:n1C,n2C,n3C))
      ALLOCATE(F2C(n1C,0:n2C,n3C))
      ALLOCATE(F3C(n1C,n2C,0:n3C))
      ALLOCATE(DTU1C(0:n1C,n2C,n3C))
      ALLOCATE(DTU2C(n1C,0:n2C,n3C))
      ALLOCATE(DTU3C(n1C,n2C,0:n3C))
      ALLOCATE(DUU1C(0:n1C,n2C,n3C))
      ALLOCATE(DUU2C(n1C,0:n2C,n3C))
      ALLOCATE(DUU3C(n1C,n2C,0:n3C))
      ALLOCATE(VolGC(n1C,n2C,n3C))
      ALLOCATE(DUTGC(n1C,n2C,n3C))
!     CALL WeightCoarse(n1,n2,n3,F1,F2,F3,F1T,F2T,F3T,VolG,VolGP, &
!                       n1C,n2C,n3C,F1C,F2C,F3C,F1CT,F2CT,F3CT,VolGC,VolGCP, &
!                       Shift1,Shift2,Shift3)
      CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                      n1C,n2C,n3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                      n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                      n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL VolCoarse(n1,n2,n3,VolG, &
                     n1C,n2C,n3C,VolGC, &
                     Shift1,Shift2,Shift3)
      CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                     n1C,n2C,n3C,DUTGC,VolGC, &
                     Shift1,Shift2,Shift3)
      DO k=1,n3C
        CALL Allocate(Current2%d(k)%P%Coarse)
        Current2%d(k)%P%Coarse%Fine=>Current2%d(k)%P
        Neu2%d(k)%P=>Current2%d(k)%P%Coarse
        CALL Allocate(Neu2%d(k)%P,n2C,n1C)
      END DO
      Current2=>Neu2
      CALL Compute(Current2,F1C,F2C,F3C,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
      DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
    END DO
    DEALLOCATE(Neu2%d)
    DEALLOCATE(Neu2)
    IF (Level3<RefLevel-1)  THEN
      CALL Allocate(Current3%Coarse)
      Neu3=>Current3%Coarse
      Neu3%Fine=>Current3
      Current3=>Current3%Coarse
    END IF
  END DO
  DEALLOCATE(F1,F2,F3,DTU1,DTU2,DTU3,DUU1,DUU2,DUU3,VolG,DUTG)

END SUBROUTINE LaplTDualCompute

SUBROUTINE LaplTFullCompute(A)

  TYPE(TriTDiag3D_T), POINTER :: A

  INTEGER :: k
  INTEGER :: n1C,n2C,n3C
  INTEGER :: Shift1,Shift2,Shift3
  INTEGER :: Level1,Level2,Level3

  TYPE(TriTDiag3D_T), POINTER :: Neu3
  TYPE(TriTDiag3D_T), POINTER :: Current3
  TYPE(TriTDiag3D_T), POINTER :: Neu2
  TYPE(TriTDiag3D_T), POINTER :: Current2

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: d1G(:)
  REAL(RealKind), ALLOCATABLE :: d2G(:)
  REAL(RealKind), ALLOCATABLE :: d3G(:)
  REAL(RealKind), ALLOCATABLE :: Metr12G(:)
  REAL(RealKind), ALLOCATABLE :: Metr13G(:)
  REAL(RealKind), ALLOCATABLE :: Metr21G(:)
  REAL(RealKind), ALLOCATABLE :: Metr23G(:)
  REAL(RealKind), ALLOCATABLE :: Metr31G(:)
  REAL(RealKind), ALLOCATABLE :: Metr32G(:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: d1GC(:)
  REAL(RealKind), ALLOCATABLE :: d2GC(:)
  REAL(RealKind), ALLOCATABLE :: d3GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr12GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr13GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr21GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr23GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr31GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr32GC(:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)

  CALL ReverseIndices(nx,ny,nz)
  CALL SetCoarse
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(d1G(1:n1))
  ALLOCATE(d2G(1:n2))
  ALLOCATE(d3G(1:n3))
  ALLOCATE(Metr12G(1:n2))
  ALLOCATE(Metr13G(1:n3))
  ALLOCATE(Metr21G(1:n1))
  ALLOCATE(Metr23G(1:n3))
  ALLOCATE(Metr31G(1:n1))
  ALLOCATE(Metr32G(1:n2))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL ReverseIncr(d1G,d2G,d3G,dx,dy,dz)
  CALL ReverseMetr(Metr12G,Metr13G, &
                   Metr21G,Metr23G, &
                   Metr31G,Metr32G, &
                   MetrXY,MetrXZ, &
                   MetrYX,MetrYZ, &
                   MetrZX,MetrZY)
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  CALL Allocate(A)
  Current3=>A
  DO Level3=0,RefLevel-1
    Shift3=Coarse3**Level3-1
    Shift2=Coarse2**Level3-1
    Shift1=Coarse1**Level3-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(d1GC(n1C))
    ALLOCATE(d2GC(n2C))
    ALLOCATE(d3GC(n3C))
    ALLOCATE(Metr12GC(1:n2C))
    ALLOCATE(Metr13GC(1:n3C))
    ALLOCATE(Metr21GC(1:n1C))
    ALLOCATE(Metr23GC(1:n3C))
    ALLOCATE(Metr31GC(1:n1C))
    ALLOCATE(Metr32GC(1:n2C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL IncrCoarse(n1,n2,n3,d1G,d2G,d3G, &
                    n1C,n2C,n3C,d1GC,d2GC,d3GC, &
                    Shift1,Shift2,Shift3)
    CALL CIncrCoarse(n2,Metr12G,d2G, &
                     n2C,Metr12GC,d2GC, &
                     Shift2)
    CALL CIncrCoarse(n3,Metr13G,d3G, &
                     n3C,Metr13GC,d3GC, &
                     Shift3)
    CALL CIncrCoarse(n1,Metr21G,d1G, &
                     n1C,Metr21GC,d1GC, &
                     Shift1)
    CALL CIncrCoarse(n3,Metr23G,d3G, &
                     n3C,Metr23GC,d3GC, &
                     Shift3)
    CALL CIncrCoarse(n1,Metr31G,d1G, &
                     n1C,Metr31GC,d1GC, &
                     Shift1)
    CALL CIncrCoarse(n2,Metr32G,d2G, &
                     n2C,Metr32GC,d2GC, &
                     Shift2)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL Allocate(Current3,n3C,n2C,n1C)
    CALL Compute(Current3,F1C,F2C,F3C, &
                 d1GC,Metr12GC,Metr13GC, &
                 d2GC,Metr21GC,Metr23GC, &
                 d3GC,Metr31GC,Metr32GC, &
                 VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
    DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,d1GC,d2GC,d3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
    DEALLOCATE(Metr12GC,Metr13GC,Metr21GC,Metr23GC,Metr31GC,Metr32GC)
    Current2=>Current3
    ALLOCATE(Neu2)
    CALL Allocate(Neu2,n3C)
    DO Level2=Level3+1,RefLevel-1
      Shift2=Coarse2**Level2-1
      Shift1=Coarse1**Level2-1
      n1C=(n1+Shift1)/(Shift1+1)
      n2C=(n2+Shift2)/(Shift2+1)
      ALLOCATE(F1C(0:n1C,n2C,n3C))
      ALLOCATE(F2C(n1C,0:n2C,n3C))
      ALLOCATE(F3C(n1C,n2C,0:n3C))
      ALLOCATE(d1GC(n1C))
      ALLOCATE(d2GC(n2C))
      ALLOCATE(d3GC(n3C))
      ALLOCATE(Metr12GC(1:n2C))
      ALLOCATE(Metr13GC(1:n3C))
      ALLOCATE(Metr21GC(1:n1C))
      ALLOCATE(Metr23GC(1:n3C))
      ALLOCATE(Metr31GC(1:n1C))
      ALLOCATE(Metr32GC(1:n2C))
      ALLOCATE(DTU1C(0:n1C,n2C,n3C))
      ALLOCATE(DTU2C(n1C,0:n2C,n3C))
      ALLOCATE(DTU3C(n1C,n2C,0:n3C))
      ALLOCATE(DUU1C(0:n1C,n2C,n3C))
      ALLOCATE(DUU2C(n1C,0:n2C,n3C))
      ALLOCATE(DUU3C(n1C,n2C,0:n3C))
      ALLOCATE(VolGC(n1C,n2C,n3C))
      ALLOCATE(DUTGC(n1C,n2C,n3C))
      CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                      n1C,n2C,n3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                      n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL IncrCoarse(n1,n2,n3,d1G,d2G,d3G, &
                      n1C,n2C,n3C,d1GC,d2GC,d3GC, &
                      Shift1,Shift2,Shift3)
      CALL CIncrCoarse(n2,Metr12G,d2G, &
                       n2C,Metr12GC,d2GC, &
                       Shift2)
      CALL CIncrCoarse(n3,Metr13G,d3G, &
                       n3C,Metr13GC,d3GC, &
                       Shift3)
      CALL CIncrCoarse(n1,Metr21G,d1G, &
                       n1C,Metr21GC,d1GC, &
                       Shift1)
      CALL CIncrCoarse(n3,Metr23G,d3G, &
                       n3C,Metr23GC,d3GC, &
                       Shift3)
      CALL CIncrCoarse(n1,Metr31G,d1G, &
                       n1C,Metr31GC,d1GC, &
                       Shift1)
      CALL CIncrCoarse(n2,Metr32G,d2G, &
                       n2C,Metr32GC,d2GC, &
                       Shift2)
      CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                      n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL VolCoarse(n1,n2,n3,VolG, &
                     n1C,n2C,n3C,VolGC, &
                     Shift1,Shift2,Shift3)
      CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                     n1C,n2C,n3C,DUTGC,VolGC, &
                     Shift1,Shift2,Shift3)
      DO k=1,n3C
        CALL Allocate(Current2%d(k)%P%Coarse)
        Current2%d(k)%P%Coarse%Fine=>Current2%d(k)%P
        Neu2%d(k)%P=>Current2%d(k)%P%Coarse
        CALL Allocate(Neu2%d(k)%P,n2C,n1C)
      END DO
      Current2=>Neu2
      CALL Compute(Current2,F1C,F2C,F3C, &
                   d1GC,Metr12GC,Metr13GC, &
                   d2GC,Metr21GC,Metr23GC, &
                   d3GC,Metr31GC,Metr32GC, &
                   VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
      DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,d1GC,d2GC,d3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
      DEALLOCATE(Metr12GC,Metr13GC,Metr21GC,Metr23GC,Metr31GC,Metr32GC)
    END DO
    DEALLOCATE(Neu2%d)
    DEALLOCATE(Neu2)
    IF (Level3<RefLevel-1)  THEN
      CALL Allocate(Current3%Coarse)
      Neu3=>Current3%Coarse
      Neu3%Fine=>Current3
      Current3=>Current3%Coarse
    END IF
  END DO
  DEALLOCATE(F1,F2,F3,DTU1,DTU2,DTU3,d1G,d2G,d3G,DUU1,DUU2,DUU3,VolG,DUTG)
  DEALLOCATE(Metr12G,Metr13G,Metr21G,Metr23G,Metr31G,Metr32G)

END SUBROUTINE LaplTFullCompute

SUBROUTINE LaplTBCompute(A)

  TYPE(TriTBDiag3D_T), POINTER :: A

  INTEGER :: k
  INTEGER :: n1C,n2C,n3C
  INTEGER :: Shift1,Shift2,Shift3
  INTEGER :: Level1,Level2,Level3

  TYPE(TriTBDiag3D_T), POINTER :: Neu3
  TYPE(TriTBDiag3D_T), POINTER :: Current3
  TYPE(TriTBDiag3D_T), POINTER :: Neu2
  TYPE(TriTBDiag3D_T), POINTER :: Current2

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)

  CALL ReverseIndices(nx,ny,nz)
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(F1G(0:n1,n2,n3))
  ALLOCATE(F2G(n1,0:n2,n3))
  ALLOCATE(F3G(n1,n2,0:n3))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL ReverseWeight(F1G,F2G,F3G,FUG(ix0:ix1,:,:),FVG(:,iy0:iy1,:),FWG(:,:,iz0:iz1))
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  CALL Allocate(A)
  Current3=>A
  DO Level3=0,RefLevel-1
    Shift3=2**Level3-1
    Shift2=2**Level3-1
    Shift1=2**Level3-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(F1GC(0:n1C,n2C,n3C))
    ALLOCATE(F2GC(n1C,0:n2C,n3C))
    ALLOCATE(F3GC(n1C,n2C,0:n3C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL FaceCoarse(n1,n2,n3,F1G,F2G,F3G, &
                    n1C,n2C,n3C,F1GC,F2GC,F3GC, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1G,F2G,F3G, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1GC,F2GC,F3GC, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL Allocate(Current3,n3C,n2C,n1C)
    CALL Compute(Current3,F1C,F2C,F3C,F1GC,F2GC,F3GC,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
    DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,F1GC,F2GC,F3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
    Current2=>Current3
    ALLOCATE(Neu2)
    CALL Allocate(Neu2,n3C)
    DO Level2=Level3+1,RefLevel-1
      Shift2=2**Level2-1
      Shift1=2**Level2-1
      n1C=(n1+Shift1)/(Shift1+1)
      n2C=(n2+Shift2)/(Shift2+1)
      ALLOCATE(F1C(0:n1C,n2C,n3C))
      ALLOCATE(F2C(n1C,0:n2C,n3C))
      ALLOCATE(F3C(n1C,n2C,0:n3C))
      ALLOCATE(DTU1C(0:n1C,n2C,n3C))
      ALLOCATE(DTU2C(n1C,0:n2C,n3C))
      ALLOCATE(DTU3C(n1C,n2C,0:n3C))
      ALLOCATE(F1GC(0:n1C,n2C,n3C))
      ALLOCATE(F2GC(n1C,0:n2C,n3C))
      ALLOCATE(F3GC(n1C,n2C,0:n3C))
      ALLOCATE(DUU1C(0:n1C,n2C,n3C))
      ALLOCATE(DUU2C(n1C,0:n2C,n3C))
      ALLOCATE(DUU3C(n1C,n2C,0:n3C))
      ALLOCATE(VolGC(n1C,n2C,n3C))
      ALLOCATE(DUTGC(n1C,n2C,n3C))
      CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                      n1C,n2C,n3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                      n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL FaceCoarse(n1,n2,n3,F1G,F2G,F3G, &
                      n1C,n2C,n3C,F1GC,F2GC,F3GC, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1G,F2G,F3G, &
                      n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1GC,F2GC,F3GC, &
                      Shift1,Shift2,Shift3)
      CALL VolCoarse(n1,n2,n3,VolG, &
                     n1C,n2C,n3C,VolGC, &
                     Shift1,Shift2,Shift3)
      CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                     n1C,n2C,n3C,DUTGC,VolGC, &
                     Shift1,Shift2,Shift3)
      DO k=1,n3C
        CALL Allocate(Current2%d(k)%P%Coarse)
        Current2%d(k)%P%Coarse%Fine=>Current2%d(k)%P
        Neu2%d(k)%P=>Current2%d(k)%P%Coarse
        CALL Allocate(Neu2%d(k)%P,n2C,n1C)
      END DO
      Current2=>Neu2
      CALL Compute(Current2,F1C,F2C,F3C,F1GC,F2GC,F3GC,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC)
      DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,F1GC,F2GC,F3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC)
    END DO
    DEALLOCATE(Neu2%d)
    DEALLOCATE(Neu2)
    IF (Level3<RefLevel-1)  THEN
      CALL Allocate(Current3%Coarse)
      Neu3=>Current3%Coarse
      Neu3%Fine=>Current3
      Current3=>Current3%Coarse
    END IF
  END DO
  DEALLOCATE(F1,F2,F3,DTU1,DTU2,DTU3,F1G,F2G,F3G,DUU1,DUU2,DUU3,VolG,DUTG)

END SUBROUTINE LaplTBCompute


SUBROUTINE LaplTRCompute(A)

  TYPE(TriTRDiag3D_T), POINTER :: A

  INTEGER :: k
  INTEGER :: n1C,n2C,n3C
  INTEGER :: Shift1,Shift2,Shift3
  INTEGER :: Level1,Level2,Level3

  TYPE(TriTRDiag3D_T), POINTER :: Neu3
  TYPE(TriTRDiag3D_T), POINTER :: Current3
  TYPE(TriTRDiag3D_T), POINTER :: Neu2
  TYPE(TriTRDiag3D_T), POINTER :: Current2

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3G(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DURG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3GC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DURGC(:,:,:)

  CALL ReverseIndices(nx,ny,nz)
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(F1G(0:n1,n2,n3))
  ALLOCATE(F2G(n1,0:n2,n3))
  ALLOCATE(F3G(n1,n2,0:n3))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  ALLOCATE(DURG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL ReverseWeight(F1G,F2G,F3G,FUG(ix0:ix1,:,:),FVG(:,iy0:iy1,:),FWG(:,:,iz0:iz1))
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(DURG,DUR(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  CALL Allocate(A)
  Current3=>A
  DO Level3=0,RefLevel-1
    Shift3=2**Level3-1
    Shift2=2**Level3-1
    Shift1=2**Level3-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(F1GC(0:n1C,n2C,n3C))
    ALLOCATE(F2GC(n1C,0:n2C,n3C))
    ALLOCATE(F3GC(n1C,n2C,0:n3C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    ALLOCATE(DURGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL FaceCoarse(n1,n2,n3,F1G,F2G,F3G, &
                    n1C,n2C,n3C,F1GC,F2GC,F3GC, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1G,F2G,F3G, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1GC,F2GC,F3GC, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DURG,VolG, &
                   n1C,n2C,n3C,DURGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL Allocate(Current3,n3C,n2C,n1C)
    CALL Compute(Current3,F1C,F2C,F3C,F1GC,F2GC,F3GC,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC,DURGC)
    DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,F1GC,F2GC,F3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC,DURGC)
    Current2=>Current3
    ALLOCATE(Neu2)
    CALL Allocate(Neu2,n3C)
    DO Level2=Level3+1,RefLevel-1
      Shift2=2**Level2-1
      Shift1=2**Level2-1
      n1C=(n1+Shift1)/(Shift1+1)
      n2C=(n2+Shift2)/(Shift2+1)
      ALLOCATE(F1C(0:n1C,n2C,n3C))
      ALLOCATE(F2C(n1C,0:n2C,n3C))
      ALLOCATE(F3C(n1C,n2C,0:n3C))
      ALLOCATE(DTU1C(0:n1C,n2C,n3C))
      ALLOCATE(DTU2C(n1C,0:n2C,n3C))
      ALLOCATE(DTU3C(n1C,n2C,0:n3C))
      ALLOCATE(F1GC(0:n1C,n2C,n3C))
      ALLOCATE(F2GC(n1C,0:n2C,n3C))
      ALLOCATE(F3GC(n1C,n2C,0:n3C))
      ALLOCATE(DUU1C(0:n1C,n2C,n3C))
      ALLOCATE(DUU2C(n1C,0:n2C,n3C))
      ALLOCATE(DUU3C(n1C,n2C,0:n3C))
      ALLOCATE(VolGC(n1C,n2C,n3C))
      ALLOCATE(DUTGC(n1C,n2C,n3C))
      ALLOCATE(DURGC(n1C,n2C,n3C))
      CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                      n1C,n2C,n3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                      n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                      Shift1,Shift2,Shift3)
      CALL FaceCoarse(n1,n2,n3,F1G,F2G,F3G, &
                      n1C,n2C,n3C,F1GC,F2GC,F3GC, &
                      Shift1,Shift2,Shift3)
      CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1G,F2G,F3G, &
                      n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1GC,F2GC,F3GC, &
                      Shift1,Shift2,Shift3)
      CALL VolCoarse(n1,n2,n3,VolG, &
                     n1C,n2C,n3C,VolGC, &
                     Shift1,Shift2,Shift3)
      CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                     n1C,n2C,n3C,DUTGC,VolGC, &
                     Shift1,Shift2,Shift3)
      CALL CVolCoarse(n1,n2,n3,DURG,VolG, &
                     n1C,n2C,n3C,DURGC,VolGC, &
                     Shift1,Shift2,Shift3)
      DO k=1,n3C
        CALL Allocate(Current2%d(k)%P%Coarse)
        Current2%d(k)%P%Coarse%Fine=>Current2%d(k)%P
        Neu2%d(k)%P=>Current2%d(k)%P%Coarse
        CALL Allocate(Neu2%d(k)%P,n2C,n1C)
      END DO
      Current2=>Neu2
      CALL Compute(Current2,F1C,F2C,F3C,F1GC,F2GC,F3GC,VolGC,DTU1C,DTU2C,DTU3C,DUU1C,DUU2C,DUU3C,DUTGC,DURGC)
      DEALLOCATE(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,F1GC,F2GC,F3GC,DUU1C,DUU2C,DUU3C,VolGC,DUTGC,DURGC)
    END DO
    DEALLOCATE(Neu2%d)
    DEALLOCATE(Neu2)
    IF (Level3<RefLevel-1)  THEN
      CALL Allocate(Current3%Coarse)
      Neu3=>Current3%Coarse
      Neu3%Fine=>Current3
      Current3=>Current3%Coarse
    END IF
  END DO
  DEALLOCATE(F1,F2,F3,DTU1,DTU2,DTU3,F1G,F2G,F3G,DUU1,DUU2,DUU3,VolG,DUTG,DURG)

END SUBROUTINE LaplTRCompute

FUNCTION Index3(i1,i2,i3,n1,n2,n3)
  INTEGER :: Index3
  INTEGER :: i1,i2,i3,n1,n2,n3
  Index3=i1+((i2-1)+(i3-1)*n2)*n1
END FUNCTION Index3

FUNCTION DivCompute(FU,FV,FW,DTU,DTV,DTW,DTT,VolC)
  

  TYPE(SpRowCol), TARGET :: DivCompute
  REAL(RealKind) :: FU(:,:,:)
  REAL(RealKind) :: FV(:,:,:)
  REAL(RealKind) :: FW(:,:,:)
  REAL(RealKind) :: DTU(:,:,:)
  REAL(RealKind) :: DTV(:,:,:)
  REAL(RealKind) :: DTW(:,:,:)
  REAL(RealKind) :: DTT(:,:,:)
  REAL(RealKind) :: VolC(:,:,:)

  INTEGER :: n,m,l
  INTEGER :: nU,nV,nW
  INTEGER :: ic,iz
  INTEGER :: i,j,k
  INTEGER :: i0,in,j0,jm,k0,kl
  TYPE(SpInd) :: A
  TYPE(SpRowCol), POINTER :: DivComputeP

  n=SIZE(VolC,1)
  m=SIZE(VolC,2)
  l=SIZE(VolC,3)
  nU=(n+1)*m*l
  nV=n*(m+1)*l
  nW=n*m*(l+1)

  A%m=n*m*l
  A%n=nU+nV+nW
  A%NumNonZero=6*m*n*l
  CALL SpNullify(A)
  CALL Allocate(A)
  A%Val=Zero

  ic=1
  iz=1
  DO k=1,l
    DO j=1,m
      DO i=1,n
        ic=Index3(i,j,k,n,m,l)
        A%ColInd(iz)=Index3(i,j,k,n+1,m,l)
        A%RowInd(iz)=ic
        A%Val(iz)=-FU(i,j,k)*DTU(i,j,k)*DTT(i,j,k) 
        iz=iz+1
        A%ColInd(iz)=Index3(i+1,j,k,n+1,m,l)
        A%RowInd(iz)=ic
        A%Val(iz)=FU(i+1,j,k)*DTU(i+1,j,k)*DTT(i,j,k) 
        iz=iz+1

        A%ColInd(iz)=Index3(i,j,k,n,m+1,l)+nU
        A%RowInd(iz)=ic
        A%Val(iz)=-FV(i,j,k)*DTV(i,j,k)*DTT(i,j,k) 
        iz=iz+1
        A%ColInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
        A%RowInd(iz)=ic
        A%Val(iz)=FV(i,j+1,k)*DTV(i,j+1,k)*DTT(i,j,k) 
        iz=iz+1

        A%ColInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
        A%RowInd(iz)=ic
        A%Val(iz)=-FW(i,j,k)*DTW(i,j,k)*DTT(i,j,k) 
        iz=iz+1
        A%ColInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
        A%RowInd(iz)=ic
        A%Val(iz)=FW(i,j,k+1)*DTW(i,j,k+1)*DTT(i,j,k) 
        iz=iz+1
      END DO
    END DO
  END DO
  A%NumNonZero=iz-1 
  A%Val=beta0*dtP*A%Val
  DivComputeP=>DivCompute
  DivComputeP=A
  CALL Deallocate(A)

END FUNCTION DivCompute 

FUNCTION DivTRCompute(FU,FV,FW,DTU,DTV,DTW,VolC)
  

  TYPE(SpRowColMat), TARGET :: DivTRCompute
  REAL(RealKind) :: FU(:,:,:)
  REAL(RealKind) :: FV(:,:,:)
  REAL(RealKind) :: FW(:,:,:)
  REAL(RealKind) :: DTU(:,:,:)
  REAL(RealKind) :: DTV(:,:,:)
  REAL(RealKind) :: DTW(:,:,:)
  REAL(RealKind) :: VolC(:,:,:)

  INTEGER :: n,m,l
  INTEGER :: nU,nV,nW
  INTEGER :: ic,iz
  INTEGER :: i,j,k
  INTEGER :: i0,in,j0,jm,k0,kl
  TYPE(SpIndMat) :: A
  TYPE(SpRowColMat), POINTER :: DivComputeP
  n=SIZE(VolC,1)
  m=SIZE(VolC,2)
  l=SIZE(VolC,3)
  nU=(n+1)*m*l
  nV=n*(m+1)*l
  nW=n*m*(l+1)

  A%m=n*m*l
  A%n=nU+nV+nW
  A%NumNonZero=6*m*n*l
  CALL SpNullify(A)
  CALL Allocate(A)
  A%Val=Zero

  iz=1
  DO k=1,l
    DO j=1,m
      DO i=1,n
        ic=IndexC(i,j,k,n,m,l)
        A%ColInd(iz)=Index3(i,j,k,n+1,m,l)
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=-FU(i,j,k)*DTU(i,j,k) 
        A%Val(iz)%a(2,1)=-FU(i,j,k)
        iz=iz+1
        A%ColInd(iz)=Index3(i+1,j,k,n+1,m,l)
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=FU(i+1,j,k)*DTU(i+1,j,k) 
        A%Val(iz)%a(2,1)=FU(i+1,j,k) 
        iz=iz+1

        A%ColInd(iz)=Index3(i,j,k,n,m+1,l)+nU
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=-FV(i,j,k)*DTV(i,j,k) 
        A%Val(iz)%a(2,1)=-FV(i,j,k) 
        iz=iz+1
        A%ColInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=FV(i,j+1,k)*DTV(i,j+1,k) 
        A%Val(iz)%a(2,1)=FV(i,j+1,k) 
        iz=iz+1

        A%ColInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=-FW(i,j,k)*DTW(i,j,k) 
        A%Val(iz)%a(2,1)=-FW(i,j,k) 
        iz=iz+1
        A%ColInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
        A%RowInd(iz)=ic
        A%Val(iz)%a(1,1)=FW(i,j,k+1)*DTW(i,j,k+1) 
        A%Val(iz)%a(2,1)=FW(i,j,k+1) 
        iz=iz+1
      END DO
    END DO
  END DO
  A%NumNonZero=iz-1 
  A%Val=(beta0*dtP)*A%Val
  DivComputeP=>DivTRCompute
  DivComputeP=A
  CALL Deallocate(A)

END FUNCTION DivTRCompute 

FUNCTION GradDualCompute(FU,FV,FW,DUT,DUU,DUV,DUW,VolC)

  TYPE(SpRowCol), TARGET :: GradDualCompute
  REAL(RealKind) :: FU(:,:,:)
  REAL(RealKind) :: FV(:,:,:)
  REAL(RealKind) :: FW(:,:,:)
  REAL(RealKind) :: DUT(:,:,:)
  REAL(RealKind) :: DUU(:,:,:)
  REAL(RealKind) :: DUV(:,:,:)
  REAL(RealKind) :: DUW(:,:,:)
  REAL(RealKind) :: VolC(:,:,:)

  INTEGER :: n,m,l
  INTEGER :: nU,nV,nW
  INTEGER :: ic,iz
  INTEGER :: i,j,k
  INTEGER :: i0,in,j0,jm,k0,kl
  TYPE(SpInd) :: A
  TYPE(SpRowCol), POINTER :: GradComputeP

  n=SIZE(VolC,1)
  m=SIZE(VolC,2)
  l=SIZE(VolC,3)
  nU=(n+1)*m*l
  nV=n*(m+1)*l
  nW=n*m*(l+1)

  A%m=nU+nV+nW
  A%n=n*m*l
  A%NumNonZero=6*m*n*l
  CALL SpNullify(A)
  CALL Allocate(A)
  A%VAl=Zero

  ic=1
  iz=1
  DO k=1,l
    DO j=1,m
      DO i=1,n
        ic=Index3(i,j,k,n,m,l)
        IF (i>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*FU(i,j,k)*DUT(i,j,k)*DUU(i,j,k) &
                    /(VolC(i-1,j,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE
          A%RowInd(iz)=Index3(i,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac1_L*FU(i,j,k)*DUT(i,j,k)*DUU(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF 
        IF (i<n) THEN
          A%RowInd(iz)=Index3(i+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=Two*FU(i+1,j,k)*DUT(i,j,k)*DUU(i+1,j,k) &
                    /(VolC(i+1,j,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=Fac1_R*FU(i+1,j,k)*DUT(i,j,k)*DUU(i+1,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF

        IF (j>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*FV(i,j,k)*DUT(i,j,k)*DUV(i,j,k) &
                    /(VolC(i,j-1,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac2_L*FV(i,j,k)*DUT(i,j,k)*DUV(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF 
        IF (j<m) THEN
          A%RowInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=Two*FV(i,j+1,k)*DUT(i,j,k)*DUV(i,j+1,k) &
                    /(VolC(i,j+1,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=Fac2_R*FV(i,j+1,k)*DUT(i,j,k)*DUV(i,j+1,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF

        IF (k>1) THEN
          A%RowInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*FW(i,j,k)*DUT(i,j,k)*DUW(i,j,k) &
                    /(VolC(i,j,k-1)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac3_L*FW(i,j,k)*DUT(i,j,k)*DUW(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF
        IF (k<l) THEN
          A%RowInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=Two*FW(i,j,k+1)*DUT(i,j,k)*DUW(i,j,k+1) &
                    /(VolC(i,j,k+1)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=Fac3_R*FW(i,j,k+1)*DUT(i,j,k)*DUW(i,j,k+1) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF
      END DO
    END DO
  END DO
  A%NumNonZero=iz-1 
  A%Val=beta0*dtP*A%Val
  GradComputeP=>GradDualCompute
  GradComputeP=A
  CALL Deallocate(A)

END FUNCTION GradDualCompute 

FUNCTION GradFullCompute(d1,Metr12,Metr13, &
                         d2,Metr21,Metr23, &
                         d3,Metr31,Metr32, &
                         DUT,DUU,DUV,DUW)
  
  TYPE(SpRowCol), TARGET :: GradFullCompute
  REAL(RealKind) :: d1(:),Metr12(:),Metr13(:)
  REAL(RealKind) :: d2(:),Metr21(:),Metr23(:)
  REAL(RealKind) :: d3(:),Metr31(:),Metr32(:)
  REAL(RealKind) :: DUT(:,:,:)
  REAL(RealKind) :: DUU(:,:,:)
  REAL(RealKind) :: DUV(:,:,:)
  REAL(RealKind) :: DUW(:,:,:)

  INTEGER :: n,m,l
  INTEGER :: nU,nV,nW
  INTEGER :: ic,iz
  INTEGER :: i,j,k
  INTEGER :: i0,in,j0,jm,k0,kl
  TYPE(SpInd) :: A
  TYPE(SpRowCol), POINTER :: GradComputeP

  n=SIZE(d1,1)
  m=SIZE(d2,1)
  l=SIZE(d3,1)
  nU=(n+1)*m*l
  nV=n*(m+1)*l
  nW=n*m*(l+1)

  A%m=nU+nV+nW
  A%n=n*m*l
  A%NumNonZero=6*m*n*l
  CALL SpNullify(A)
  CALL Allocate(A)
  A%VAl=Zero

  ic=1
  iz=1
  DO k=1,l
    DO j=1,m
      DO i=1,n
        ic=Index3(i,j,k,n,m,l)
        IF (i>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*DUT(i,j,k)*DUU(i,j,k) &
                    /((d1(i-1)+d1(i))*Metr12(j)*Metr13(k))
          iz=iz+1
        ELSE
          A%RowInd(iz)=Index3(i,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac1_L*DUT(i,j,k)*DUU(i,j,k) &
                    /(d1(i)*Metr12(j)*Metr13(k))
          iz=iz+1
        END IF 
        IF (i<n) THEN
          A%RowInd(iz)=Index3(i+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=Two*DUT(i,j,k)*DUU(i+1,j,k) &
                    /((d1(i+1)+d1(i))*Metr12(j)*Metr13(k))
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)=Fac1_R*DUT(i,j,k)*DUU(i+1,j,k) &
                    /(d1(i)*Metr12(j)*Metr13(k))
          iz=iz+1
        END IF

        IF (j>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*DUT(i,j,k)*DUV(i,j,k) &
                    /((d2(j-1)+d2(j))*Metr21(i)*Metr23(k))
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac2_L*DUT(i,j,k)*DUV(i,j,k) &
                    /(d2(j)*Metr21(i)*Metr23(k))
          iz=iz+1
        END IF 
        IF (j<m) THEN
          A%RowInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=Two*DUT(i,j,k)*DUV(i,j+1,k) &
                    /((d2(j+1)+d2(j))*Metr21(i)*Metr23(k))
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)=Fac2_R*DUT(i,j,k)*DUV(i,j+1,k) &
                    /(d2(j)*Metr21(i)*Metr23(k))
          iz=iz+1
        END IF

        IF (k>1) THEN
          A%RowInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=-Two*DUT(i,j,k)*DUW(i,j,k) &
                    /((d3(k-1)+d3(k))*Metr31(i)*Metr32(j))
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=-Fac3_L*DUT(i,j,k)*DUW(i,j,k) &
                    /(d3(k)*Metr31(i)*Metr32(j))
          iz=iz+1
        END IF
        IF (k<l) THEN
          A%RowInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=Two*DUT(i,j,k)*DUW(i,j,k+1) &
                    /((d3(k+1)+d3(k)*Metr31(i)*Metr32(j)))
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)=Fac3_R*DUT(i,j,k)*DUW(i,j,k+1) &
                    /(d3(k)*Metr31(i)*Metr32(j))
          iz=iz+1
        END IF
      END DO
    END DO
  END DO
  A%NumNonZero=iz-1 
  A%Val=beta0*dtP*A%Val
  GradComputeP=>GradFullCompute
  GradComputeP=A
  CALL Deallocate(A)

END FUNCTION GradFullCompute 

FUNCTION GradTRCompute(FUG,FVG,FWG,DUT,DUR,VolC)
  

  TYPE(SpRowColMat), TARGET :: GradTRCompute
  REAL(RealKind) :: FUG(:,:,:)
  REAL(RealKind) :: FVG(:,:,:)
  REAL(RealKind) :: FWG(:,:,:)
  REAL(RealKind) :: DUT(:,:,:)
  REAL(RealKind) :: DUR(:,:,:)
  REAL(RealKind) :: VolC(:,:,:)

  INTEGER :: n,m,l
  INTEGER :: nU,nV,nW
  INTEGER :: ic,iz
  INTEGER :: i,j,k
  INTEGER :: i0,in,j0,jm,k0,kl
  TYPE(SpIndMat) :: A
  TYPE(SpRowColMat), POINTER :: GradComputeP
  
  n=SIZE(VolC,1)
  m=SIZE(VolC,2)
  l=SIZE(VolC,3)
  nU=(n+1)*m*l
  nV=n*(m+1)*l
  nW=n*m*(l+1)

  A%m=nU+nV+nW
  A%n=n*m*l
  A%NumNonZero=6*m*n*l
  CALL SpNullify(A)
  CALL Allocate(A)

  A%Val=Zero
  
  iz=1
  DO k=1,l
    DO j=1,m
      DO i=1,n
        ic=IndexC(i,j,k,n,m,l)
        IF (i>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-Two*FUG(i,j,k)*DUT(i,j,k) &
                    /(VolC(i-1,j,k)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=-Two*FUG(i,j,k)*DUR(i,j,k) &
                    /(VolC(i-1,j,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE
          A%RowInd(iz)=Index3(1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-FacW*FUG(i,j,k)*DUT(i,j,k) &
                    /(VolC(1,j,k)+Eps)
          A%Val(iz)%a(1,2)=-FacW*FUG(i,j,k)*DUR(i,j,k) &
                    /(VolC(1,j,k)+Eps)
          iz=iz+1
        END IF 
        IF (i<n) THEN
          A%RowInd(iz)=Index3(i+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=Two*FUG(i+1,j,k)*DUT(i,j,k) &
                    /(VolC(i+1,j,k)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=Two*FUG(i+1,j,k)*DUR(i,j,k) &
                    /(VolC(i+1,j,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(n+1,j,k,n+1,m,l)
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=FacE*FUG(i+1,j,k)*DUT(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=FacE*FUG(i+1,j,k)*DUR(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF

        IF (j>1) THEN 
          A%RowInd(iz)=Index3(i,j,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-Two*FVG(i,j,k)*DUT(i,j,k) &
                    /(VolC(i,j-1,k)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=-Two*FVG(i,j,k)*DUR(i,j,k) &
                    /(VolC(i,j-1,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-FacS*FVG(i,j,k)*DUT(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=-FacS*FVG(i,j,k)*DUR(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF 
        IF (j<m) THEN
          A%RowInd(iz)=Index3(i,j+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=Two*FVG(i,j+1,k)*DUT(i,j,k) &
                    /(VolC(i,j+1,k)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=Two*FVG(i,j+1,k)*DUR(i,j,k) &
                    /(VolC(i,j+1,k)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,m+1,k,n,m+1,l)+nU
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=FacN*FVG(i,j+1,k)*DUT(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=FacN*FVG(i,j+1,k)*DUR(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF

        IF (k>1) THEN
          A%RowInd(iz)=Index3(i,j,k,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-Two*FWG(i,j,k)*DUT(i,j,k) &
                    /(VolC(i,j,k-1)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=-Two*FWG(i,j,k)*DUR(i,j,k) &
                    /(VolC(i,j,k-1)+VolC(i,j,k)+Eps)&
                    -GravComp*FWG(i,j,k)/(FWG(i,j,k)+Eps)&
                    *VolC(i,j,k)/(VolC(i,j,k-1)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=-FacB*FWG(i,j,k)*DUT(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=-FacB*FWG(i,j,k)*DUR(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF
        IF (k<l) THEN
          A%RowInd(iz)=Index3(i,j,k+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=Two*FWG(i,j,k+1)*DUT(i,j,k) &
                    /(VolC(i,j,k+1)+VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=Two*FWG(i,j,k+1)*DUR(i,j,k) &
                    /(VolC(i,j,k+1)+VolC(i,j,k)+Eps)&
                    -GravComp*FWG(i,j,k+1)/(FWG(i,j,k+1)+Eps)&
                    *VolC(i,j,k)/(VolC(i,j,k+1)+VolC(i,j,k)+Eps)
          iz=iz+1
        ELSE 
          A%RowInd(iz)=Index3(i,j,l+1,n,m,l+1)+nU+nV
          A%ColInd(iz)=ic
          A%Val(iz)%a(1,1)=FacT*FWG(i,j,k+1)*DUT(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          A%Val(iz)%a(1,2)=FacT*FWG(i,j,k+1)*DUR(i,j,k) &
                    /(VolC(i,j,k)+Eps)
          iz=iz+1
        END IF
      END DO
    END DO
  END DO
  A%NumNonZero=iz-1 
  A%Val=(beta0*dtP)*A%Val
  GradComputeP=>GradTRCompute
  GradComputeP=A
  CALL Deallocate(A)

END FUNCTION GradTRCompute 

SUBROUTINE LaplMuCompute(A)
  TYPE(SpRowCol), POINTER :: A(:)
  IF (GradFull) THEN
    CALL LaplMuFullCompute(A)
  ELSE
    CALL LaplMuDualCompute(A)
  END IF
END SUBROUTINE LaplMuCompute

SUBROUTINE LaplMuDualCompute(A)

  TYPE(SpRowCol), POINTER :: A(:)

  TYPE(SpRowCol), TARGET :: Div
  TYPE(SpRowCol), POINTER :: DivP
  TYPE(SpRowCol), TARGET :: Grad
  TYPE(SpRowCol), POINTER :: GradP

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTTGC(:,:,:)

  TYPE(SpRowCol), POINTER :: AF

  INTEGER :: n1C,n2C,n3C
  INTEGER :: iRef,Shift1,Shift2,Shift3

  CALL ReverseIndices(nx,ny,nz)
  CALL SetCoarse
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DTTG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(DTTG,DTT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  DO iRef=1,RefLevel
    Shift3=Coarse3**(iRef-1)-1
    Shift2=Coarse2**(iRef-1)-1
    Shift1=Coarse1**(iRef-1)-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    CALL SpNullify(Div)
    CALL SpNullify(Grad)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    ALLOCATE(DTTGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DTTG,VolG, &
                   n1C,n2C,n3C,DTTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    AF=>A(iRef)
    DivP=>Div
    DivP=DivCompute(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,DTTGC,VolGC)
    GradP=>Grad
    GradP=GradDualCompute(F1C,F2C,F3C,DUTGC,DUU1C,DUU2C,DUU3C,VolGC)
    CALL SpMm(AF,Div,Grad)
    CALL Axpy(FacAnela,VolGC,AF)
    DEALLOCATE(F1C,F2C,F3C, &
               DTU1C,DTU2C,DTU3C,DTTGC,DUTGC,DUU1C,DUU2C,DUU3C,VolGC)
    CALL Deallocate(Div)
    CALL Deallocate(Grad)
  END DO
  DEALLOCATE(F1,F2,F3, &
             DTU1,DTU2,DTU3,DTTG,DUTG,DUU1,DUU2,DUU3,VolG)

END SUBROUTINE LaplMuDualCompute

SUBROUTINE LaplMumpsCompute(A,AMu)
  TYPE (DMUMPS_STRUC), POINTER :: A
  TYPE(SpRowCol) :: AMu
  IF (GradFull) THEN
    CALL LaplMumpsFullCompute(A,AMu)
  ELSE
!   CALL LaplMuDualCompute(A)
  END IF
END SUBROUTINE LaplMumpsCompute

SUBROUTINE LaplMumpsFullCompute(A,AMu)
  TYPE (DMUMPS_STRUC), POINTER :: A
  TYPE(SpRowCol) :: AMu

  INTEGER :: i,jj
  TYPE(SpRowCol) :: Div
  TYPE(SpRowCol) :: Grad

  Div=DivCompute(FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1) &
                ,DTU,DTV,DTW,DTT &
                ,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  Grad=GradFullCompute(dx,MetrXY,MetrXZ, &
                       dy,MetrYX,MetrYZ, &
                       dz,MetrZX,MetrZY, &
                       DUT,DUU,DUV,DUW)
  CALL SpMm(AMu,Div,Grad)
  CALL Deallocate(Div)
  CALL Deallocate(Grad)
  CALL Axpy(FacAnela,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)+Eps,AMu)

  IF (.NOT.ASSOCIATED(A%IRN)) THEN
    A%COMM=MPI_COMM_SELF%mpi_val
    A%SYM=0
    A%PAR=1
    A%JOB=-1
    A%KEEP=0
    CALL DMUMPS(A)
    ALLOCATE(A%IRN(SIZE(AMu%ColInd)))
    DO i=1,AMu%n
      DO jj=AMu%RowPtr(i),AMu%RowPtr(i+1)-1
        A%IRN(jj)=i
      END DO
    END DO  
  END IF  
  IF (.NOT.ASSOCIATED(A%JCN)) THEN
    A%JCN=>AMu%ColInd
    ! Symbolic phase
    A%JOB=1
    A%N=AMu%m
    A%NZ=SIZE(AMu%ColInd)
    A%ICNTL(1)=0
    A%ICNTL(2)=0
    A%ICNTL(3)=0
    A%ICNTL(5)=0
    A%ICNTL(18)=0
    A%ICNTL(24)=0
    A%ICNTL(28)=1
    CALL DMUMPS(A)
  END IF  
  A%A=>AMu%Val 
  A%JOB=2
  A%ICNTL(14)=20
  CALL DMUMPS(A)
  
END SUBROUTINE LaplMumpsFullCompute


SUBROUTINE LaplMuFullCompute(A)

  TYPE(SpRowCol), POINTER :: A(:)

  TYPE(SpRowCol), TARGET :: Div
  TYPE(SpRowCol), POINTER :: DivP
  TYPE(SpRowCol), TARGET :: Grad
  TYPE(SpRowCol), POINTER :: GradP

  REAL(RealKind), ALLOCATABLE :: F1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: d1G(:)
  REAL(RealKind), ALLOCATABLE :: d2G(:)
  REAL(RealKind), ALLOCATABLE :: d3G(:)
  REAL(RealKind), ALLOCATABLE :: Metr12G(:)
  REAL(RealKind), ALLOCATABLE :: Metr13G(:)
  REAL(RealKind), ALLOCATABLE :: Metr21G(:)
  REAL(RealKind), ALLOCATABLE :: Metr23G(:)
  REAL(RealKind), ALLOCATABLE :: Metr31G(:)
  REAL(RealKind), ALLOCATABLE :: Metr32G(:)
  REAL(RealKind), ALLOCATABLE :: DUU1(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTG(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: F3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: d1GC(:)
  REAL(RealKind), ALLOCATABLE :: d2GC(:)
  REAL(RealKind), ALLOCATABLE :: d3GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr12GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr13GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr21GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr23GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr31GC(:)
  REAL(RealKind), ALLOCATABLE :: Metr32GC(:)
  REAL(RealKind), ALLOCATABLE :: DUU1C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU2C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUU3C(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTTGC(:,:,:)

  TYPE(SpRowCol), POINTER :: AF

  INTEGER :: n1C,n2C,n3C
  INTEGER :: iRef,Shift1,Shift2,Shift3

  CALL ReverseIndices(nx,ny,nz)
  CALL SetCoarse
  ALLOCATE(F1(0:n1,n2,n3))
  ALLOCATE(F2(n1,0:n2,n3))
  ALLOCATE(F3(n1,n2,0:n3))
  ALLOCATE(DTU1(0:n1,n2,n3))
  ALLOCATE(DTU2(n1,0:n2,n3))
  ALLOCATE(DTU3(n1,n2,0:n3))
  ALLOCATE(d1G(1:n1))
  ALLOCATE(d2G(1:n2))
  ALLOCATE(d3G(1:n3))
  ALLOCATE(Metr12G(1:n2))
  ALLOCATE(Metr13G(1:n3))
  ALLOCATE(Metr21G(1:n1))
  ALLOCATE(Metr23G(1:n3))
  ALLOCATE(Metr31G(1:n1))
  ALLOCATE(Metr32G(1:n2))
  ALLOCATE(DUU1(0:n1,n2,n3))
  ALLOCATE(DUU2(n1,0:n2,n3))
  ALLOCATE(DUU3(n1,n2,0:n3))
  ALLOCATE(VolG(n1,n2,n3))
  ALLOCATE(DTTG(n1,n2,n3))
  ALLOCATE(DUTG(n1,n2,n3))
  CALL ReverseWeight(F1,F2,F3,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1))
  CALL ReverseWeight(DTU1,DTU2,DTU3,DTU,DTV,DTW)
  CALL ReverseWeight(DUU1,DUU2,DUU3,DUU,DUV,DUW)
  CALL ReverseIncr(d1G,d2G,d3G,dx,dy,dz)
  CALL ReverseMetr(Metr12G,Metr13G, &
                   Metr21G,Metr23G, &
                   Metr31G,Metr32G, &
                   MetrXY,MetrXZ, &
                   MetrYX,MetrYZ, &
                   MetrZX,MetrZY)
  CALL Reverse(DUTG,DUT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(DTTG,DTT(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  CALL Reverse(VolG,VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

  DO iRef=1,RefLevel
    Shift3=Coarse3**(iRef-1)-1
    Shift2=Coarse2**(iRef-1)-1
    Shift1=Coarse1**(iRef-1)-1
    n1C=(n1+Shift1)/(Shift1+1)
    n2C=(n2+Shift2)/(Shift2+1)
    n3C=(n3+Shift3)/(Shift3+1)
    CALL SpNullify(Div)
    CALL SpNullify(Grad)
    ALLOCATE(F1C(0:n1C,n2C,n3C))
    ALLOCATE(F2C(n1C,0:n2C,n3C))
    ALLOCATE(F3C(n1C,n2C,0:n3C))
    ALLOCATE(DTU1C(0:n1C,n2C,n3C))
    ALLOCATE(DTU2C(n1C,0:n2C,n3C))
    ALLOCATE(DTU3C(n1C,n2C,0:n3C))
    ALLOCATE(d1GC(n1C))
    ALLOCATE(d2GC(n2C))
    ALLOCATE(d3GC(n3C))
    ALLOCATE(Metr12GC(1:n2C))
    ALLOCATE(Metr13GC(1:n3C))
    ALLOCATE(Metr21GC(1:n1C))
    ALLOCATE(Metr23GC(1:n3C))
    ALLOCATE(Metr31GC(1:n1C))
    ALLOCATE(Metr32GC(1:n2C))
    ALLOCATE(DUU1C(0:n1C,n2C,n3C))
    ALLOCATE(DUU2C(n1C,0:n2C,n3C))
    ALLOCATE(DUU3C(n1C,n2C,0:n3C))
    ALLOCATE(VolGC(n1C,n2C,n3C))
    ALLOCATE(DUTGC(n1C,n2C,n3C))
    ALLOCATE(DTTGC(n1C,n2C,n3C))
    CALL FaceCoarse(n1,n2,n3,F1,F2,F3, &
                    n1C,n2C,n3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL CFaceCoarse(n1,n2,n3,DTU1,DTU2,DTU3,F1,F2,F3, &
                    n1C,n2C,n3C,DTU1C,DTU2C,DTU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL IncrCoarse(n1,n2,n3,d1G,d2G,d3G, &
                    n1C,n2C,n3C,d1GC,d2GC,d3GC, &
                    Shift1,Shift2,Shift3)
    CALL CIncrCoarse(n2,Metr12G,d2G, &
                     n2C,Metr12GC,d2GC, &
                     Shift2)
    CALL CIncrCoarse(n3,Metr13G,d3G, &
                     n3C,Metr13GC,d3GC, &
                     Shift3)
    CALL CIncrCoarse(n1,Metr21G,d1G, &
                     n1C,Metr21GC,d1GC, &
                     Shift1)
    CALL CIncrCoarse(n3,Metr23G,d3G, &
                     n3C,Metr23GC,d3GC, &
                     Shift3)
    CALL CIncrCoarse(n1,Metr31G,d1G, &
                     n1C,Metr31GC,d1GC, &
                     Shift1)
    CALL CIncrCoarse(n2,Metr32G,d2G, &
                     n2C,Metr32GC,d2GC, &
                     Shift2)
    CALL CFaceCoarse(n1,n2,n3,DUU1,DUU2,DUU3,F1,F2,F3, &
                    n1C,n2C,n3C,DUU1C,DUU2C,DUU3C,F1C,F2C,F3C, &
                    Shift1,Shift2,Shift3)
    CALL VolCoarse(n1,n2,n3,VolG, &
                   n1C,n2C,n3C,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DUTG,VolG, &
                   n1C,n2C,n3C,DUTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    CALL CVolCoarse(n1,n2,n3,DTTG,VolG, &
                   n1C,n2C,n3C,DTTGC,VolGC, &
                   Shift1,Shift2,Shift3)
    AF=>A(iRef)
    DivP=>Div
    DivP=DivCompute(F1C,F2C,F3C,DTU1C,DTU2C,DTU3C,DTTGC,VolGC)
    GradP=>Grad
    GradP=GradFullCompute(d1GC,Metr12GC,Metr13GC, &
                          d2GC,Metr21GC,Metr23GC, &
                          d3GC,Metr31GC,Metr32GC, &
                          DUTGC,DUU1C,DUU2C,DUU3C)
    CALL SpMm(AF,Div,Grad)
    CALL Axpy(FacAnela,VolGC,AF)
    DEALLOCATE(F1C,F2C,F3C,d1GC,d2GC,d3GC, &
               DTU1C,DTU2C,DTU3C,DTTGC,DUTGC,DUU1C,DUU2C,DUU3C,VolGC)
    DEALLOCATE(Metr12GC,Metr13GC,Metr21GC,Metr23GC,Metr31GC,Metr32GC)
    CALL Deallocate(Div)
    CALL Deallocate(Grad)
  END DO
  DEALLOCATE(F1,F2,F3,d1G,d2G,d3G, &
             DTU1,DTU2,DTU3,DTTG,DUTG,DUU1,DUU2,DUU3,VolG)
  DEALLOCATE(Metr12G,Metr13G,Metr21G,Metr23G,Metr31G,Metr32G)

END SUBROUTINE LaplMuFullCompute

SUBROUTINE LaplComputeMuTR(A)

  TYPE(SpRowColMat), POINTER :: A(:)

  TYPE(SpRowColMat), TARGET :: Div
  TYPE(SpRowColMat), POINTER :: DivP
  TYPE(SpRowColMat), TARGET :: Grad
  TYPE(SpRowColMat), POINTER :: GradP
  TYPE(SpRowColMat) :: DivGrad

  REAL(RealKind), ALLOCATABLE :: FUC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: FVC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: FWC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: FUGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: FVGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: FWGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTUGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTVGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DTWGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DUTGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: DURGC(:,:,:)
  REAL(RealKind), ALLOCATABLE :: VolCC(:,:,:)
  TYPE(SpRowColMat), POINTER :: AF

  INTEGER :: n1C,n2C,n3C
  INTEGER :: iRef,Shift

  n1C=nx
  n2C=ny
  n3C=nz
  DO iRef=1,RefLevel
    Shift=2**(iRef-1)-1
    CALL SpNullify(Div)
    CALL SpNullify(Grad)
    CALL SpNullify(DivGrad)
    ALLOCATE(FUC(0:n1C,1:n2C,1:n3C))
    ALLOCATE(FVC(1:n1C,0:n2C,1:n3C))
    ALLOCATE(FWC(1:n1C,1:n2C,0:n3C))
    ALLOCATE(FUGC(0:n1C,1:n2C,1:n3C))
    ALLOCATE(FVGC(1:n1C,0:n2C,1:n3C))
    ALLOCATE(FWGC(1:n1C,1:n2C,0:n3C))
    ALLOCATE(DTUGC(0:n1C,1:n2C,1:n3C))
    ALLOCATE(DTVGC(1:n1C,0:n2C,1:n3C))
    ALLOCATE(DTWGC(1:n1C,1:n2C,0:n3C))
    ALLOCATE(DUTGC(1:n1C,1:n2C,1:n3C))
    ALLOCATE(DURGC(1:n1C,1:n2C,1:n3C))
    ALLOCATE(VolCC(1:n1C,1:n2C,1:n3C))
    CALL FaceCoarse(nx,ny,nz,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1), &
                    n1C,n2C,n3C,FUC,FVC,FWC, &
                    Shift,Shift,Shift)
    CALL FaceCoarse(nx,ny,nz,FUG(ix0:ix1,:,:),FVG(:,iy0:iy1,:),FWG(:,:,iz0:iz1), &
                    n1C,n2C,n3C,FUGC,FVGC,FWGC, &
                    Shift,Shift,Shift)
    CALL CFaceCoarse(nx,ny,nz,DTU,DTV,DTW,FU(ix0:ix1,:,:),FV(:,iy0:iy1,:),FW(:,:,iz0:iz1), &
                    n1C,n2C,n3C,DTUGC,DTVGC,DTWGC,FUC,FVC,FWC, &
                    Shift,Shift,Shift)
    CALL VolCoarse(nx,ny,nz,VolB, &
                   n1C,n2C,n3C,VolCC, &
                   Shift,Shift,Shift)
    CALL CVolCoarse(nx,ny,nz,DUT,VolB, &
                   n1C,n2C,n3C,DUTGC,VolCC, &
                   Shift,Shift,Shift)
    CALL CVolCoarse(nx,ny,nz,DUR,VolB, &
                   n1C,n2C,n3C,DURGC,VolCC, &
                   Shift,Shift,Shift)
    AF=>A(iRef)
    DivP=>Div
    DivP=DivTRCompute(FUC,FVC,FWC,DTUGC,DTVGC,DTWGC,VolCC)
    GradP=>Grad
    GradP=GradTRCompute(FUGC,FVGC,FWGC,DUTGC,DURGC,VolCC)
    CALL SpMm(AF,Div,Grad)
    CALL Axpy(FacAnela,VolCC,AF)
    DEALLOCATE(FUC,FVC,FWC,FUGC,FVGC,FWGC,DTUGC,DTVGC,DTWGC,DUTGC,DURGC,VolCC)
    CALL Deallocate(Div)
    CALL Deallocate(Grad)
    CALL Deallocate(DivGrad)
    n1C=(n1C+1)/2
    n2C=(n2C+1)/2
    n3C=(n3C+1)/2
  END DO

END SUBROUTINE LaplComputeMuTR

SUBROUTINE ProlongationMu(n,m,l,p,nC,mC,lC,pC)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: p(1:n,1:m,1:l)
  REAL(RealKind) :: pC(1:nC,1:mC,1:lC)

  INTEGER :: i,j,k,iC,jC,kC

  DO i=1,n
    iC=(i+1)/2
    DO j=1,m
      jC=(j+1)/2
      DO k=1,l
        kC=(k+1)/2
        p(i,j,k)=p(i,j,k)+pC(iC,jC,kC)
      END DO
    END DO
  END DO

END SUBROUTINE ProlongationMu

SUBROUTINE ProlongationMuTR(n,m,l,p,nC,mC,lC,pC)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: p(1:2,1:n,1:m,1:l)
  REAL(RealKind) :: pC(1:2,1:nC,1:mC,1:lC)

  INTEGER :: i,j,k,iC,jC,kC

  DO i=1,n
    iC=(i+1)/2
    DO j=1,m
      jC=(j+1)/2
      DO k=1,l
        kC=(k+1)/2
        p(:,i,j,k)=p(:,i,j,k)+pC(:,iC,jC,kC)
      END DO
    END DO
  END DO

END SUBROUTINE ProlongationMuTR

SUBROUTINE RestrictionMu(nC,mC,lC,ResC,n,m,l,Res)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: Res(1:n,1:m,1:l)
  REAL(RealKind) :: ResC(1:nC,1:mC,1:lC)

  INTEGER :: i,j,k,iC,jC,kC

  ResC=0.0e0
  DO i=1,n
    ic=(i+1)/2
    DO j=1,m
      jc=(j+1)/2
      DO k=1,l
        kc=(k+1)/2
        ResC(iC,jC,kC)=ResC(iC,jC,kC)+Res(i,j,k)
      END DO
    END DO
  END DO

END SUBROUTINE RestrictionMu

SUBROUTINE RestrictionMuTR(nC,mC,lC,ResC,n,m,l,Res)

  INTEGER :: n,m,l
  INTEGER :: nC,mC,lC
  REAL(RealKind) :: Res(1:2,1:n,1:m,1:l)
  REAL(RealKind) :: ResC(1:2,1:nC,1:mC,1:lC)

  INTEGER :: i,j,k,iC,jC,kC

  ResC=0.0e0
  DO i=1,n
    ic=(i+1)/2
    DO j=1,m
      jc=(j+1)/2
      DO k=1,l
        kc=(k+1)/2
        ResC(:,iC,jC,kC)=ResC(:,iC,jC,kC)+Res(:,i,j,k)
      END DO
    END DO
  END DO

END SUBROUTINE RestrictionMuTR

SUBROUTINE CellToFaceWeight(Theta,DTU,DTV,DTW,VolC)

  REAL(RealKind) :: Theta(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: DTU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  REAL(RealKind) :: VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1)

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: ThCoarse,ThFine,ThFace

  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(  &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM( &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
             ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
             DTU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix0,jy,jz)= &
                          (Theta(ix0+1,jy,jz,1)*VolC(ix0+1,jy,jz)+ &
                           Theta(ix0,jy,jz,1)*VolC(ix0,jy,jz))/ &
                          (VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM( &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(  &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix1,jy,jz)= &
                          (Theta(ix1+1,jy,jz,1)*VolC(ix1+1,jy,jz)+ &
                           Theta(ix1,jy,jz,1)*VolC(ix1,jy,jz))/ &
                          (VolC(ix1+1,jy,jz)+VolC(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy0,jz)= &
                          (Theta(jx,iy0+1,jz,1)*VolC(jx,iy0+1,jz)+ &
                           Theta(jx,iy0,jz,1)*VolC(jx,iy0,jz))/ &
                          (VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
               ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy1,jz)= &
                          (Theta(jx,iy1+1,jz,1)*VolC(jx,iy1+1,jz)+ &
                           Theta(jx,iy1,jz,1)*VolC(jx,iy1,jz))/ &
                          (VolC(jx,iy1+1,jz)+VolC(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz0)= &
                          (Theta(jx,jy,iz0+1,1)*VolC(jx,jy,iz0+1)+ &
                           Theta(jx,jy,iz0,1)*VolC(jx,jy,iz0))/ &
                          (VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
              (SUM( &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz1)= &
                          (Theta(jx,jy,iz1+1,1)*VolC(jx,jy,iz1+1)+ &
                           Theta(jx,jy,iz1,1)*VolC(jx,jy,iz1))/ &
                          (VolC(jx,jy,iz1+1)+VolC(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO

  IF (TypeW(1:1)=='o') THEN
    IF (BCVel%West/='InFlow') THEN
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Theta(ix+1,iy,iz,1)
        END DO
      END DO
    ELSE
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(Theta(ix,iy,iz,1) &
                             +Theta(ix+1,iy,iz,1)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeE(1:1)=='o') THEN
    IF (BCVel%East/='InFlow') THEN
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Theta(ix,iy,iz,1)
        END DO
      END DO
    ELSE
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(Theta(ix,iy,iz,1) &
                             +Theta(ix+1,iy,iz,1)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    IF (BCVel%South/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Theta(ix,iy+1,iz,1)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(Theta(ix,iy,iz,1) &
                             +Theta(ix,iy+1,iz,1)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1
    IF (BCVel%North/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Theta(ix,iy,iz,1)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(Theta(ix,iy,iz,1) &
                             +Theta(ix,iy+1,iz,1)) 
                       
        END DO
      END DO
    END IF
  END IF
  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    IF (BCVel%Bottom/='InFlow') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=Theta(ix,iy,iz+1,1)
        END DO
      END DO
    ELSE
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF
  IF (TypeT(1:1)=='o') THEN
    IF (BCVel%Top/='InFlow') THEN
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=Theta(ix,iy,iz,1)
        END DO
      END DO
    ELSE 
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF

  DO ix=ix0+1,ix1-1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        DTU(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1) &
                   +VolC(ix+1,iy,iz)*Theta(ix+1,iy,iz,1)) &
                  /(VolC(ix,iy,iz)+VolC(ix+1,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1-1
      DO iz=iz0+1,iz1
        DTV(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1) &
                   +VolC(ix,iy+1,iz)*Theta(ix,iy+1,iz,1)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy+1,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1-1
        DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE CellToFaceWeight


SUBROUTINE CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)

  REAL(RealKind) :: Theta(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: Rho(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: DTU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  REAL(RealKind) :: VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1)

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: ThCoarse,ThFine,ThFace
  REAL(RealKind) :: RhoCoarse,RhoFine,RhoFace

  DTU=0.0d0
  DTV=0.0d0
  DTW=0.0d0
  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
             ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
             DTU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix0,jy,jz)= &
                          (Theta(ix0+1,jy,jz,1)/(Rho(ix0+1,jy,jz,1)+Eps)*VolC(ix0+1,jy,jz)+ &
                           Theta(ix0,jy,jz,1)/(Rho(ix0,jy,jz,1)+Eps)*VolC(ix0,jy,jz))/ &
                          (VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix1,jy,jz)= &
                          (Theta(ix1+1,jy,jz,1)/(Rho(ix1+1,jy,jz,1)+Eps)*VolC(ix1+1,jy,jz)+ &
                           Theta(ix1,jy,jz,1)/(Rho(ix1,jy,jz,1)+Eps)*VolC(ix1,jy,jz))/ &
                          (VolC(ix1+1,jy,jz)+VolC(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy0,jz)= &
                          (Theta(jx,iy0+1,jz,1)/(Rho(jx,iy0+1,jz,1)+Eps)*VolC(jx,iy0+1,jz)+ &
                           Theta(jx,iy0,jz,1)/(Rho(jx,iy0,jz,1)+Eps)*VolC(jx,iy0,jz))/ &
                          (VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
               ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy1,jz)= &
                          (Theta(jx,iy1+1,jz,1)/(Rho(jx,iy1+1,jz,1)+Eps)*VolC(jx,iy1+1,jz)+ &
                           Theta(jx,iy1,jz,1)/(Rho(jx,iy1,jz,1)+Eps)*VolC(jx,iy1,jz))/ &
                          (VolC(jx,iy1+1,jz)+VolC(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz0)= &
                          (Theta(jx,jy,iz0+1,1)*VolC(jx,jy,iz0+1)/(Rho(jx,jy,iz0+1,1)+Eps)+ &
                           Theta(jx,jy,iz0,1)*VolC(jx,jy,iz0)/(Rho(jx,jy,iz0,1)+Eps))/ &
                          (VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            ThFine= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            ThCoarse= &
               SUM(Theta(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz1)= &
                          (Theta(jx,jy,iz1+1,1)*VolC(jx,jy,iz1+1)/(Rho(jx,jy,iz1+1,1)+Eps)+ &
                           Theta(jx,jy,iz1,1)*VolC(jx,jy,iz1)/(Rho(jx,jy,iz1,1)+Eps))/ &
                          (VolC(jx,jy,iz1+1)+VolC(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO

  IF (TypeW(1:1)=='o') THEN
    IF (BCVel%West/='InFlow') THEN
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Theta(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                             +Theta(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeE(1:1)=='o') THEN
    IF (BCVel%East/='InFlow') THEN
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                             +Theta(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    IF (BCVel%South/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Theta(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                             +Theta(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1
    IF (BCVel%North/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                             +Theta(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)) 
                       
        END DO
      END DO
    END IF
  END IF
  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    IF (BCVel%Bottom/='InFlow') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=Theta(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        END DO
      END DO
    ELSE
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF
  IF (TypeT(1:1)=='o') THEN
    IF (BCVel%Top/='InFlow') THEN
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE 
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF
  DO ix=ix0+1,ix1-1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        DTU(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix+1,iy,iz)*Theta(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix+1,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1-1
      DO iz=iz0+1,iz1
        DTV(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy+1,iz)*Theta(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy+1,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1-1
        DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*Theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*Theta(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE CellToFaceRhoWeight
     

SUBROUTINE CellToFaceInvRhoWeight(Rho,DTU,DTV,DTW,VolC)

  REAL(RealKind) :: Rho(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: DTU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  REAL(RealKind) :: DTW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  REAL(RealKind) :: VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1)

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: ThCoarse,ThFine,ThFace
  REAL(RealKind) :: RhoCoarse,RhoFine,RhoFace

  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(One* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
             ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
             DTU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix0,jy,jz)= &
                          (One/(Rho(ix0+1,jy,jz,1)+Eps)*VolC(ix0+1,jy,jz)+ &
                           One/(Rho(ix0,jy,jz,1)+Eps)*VolC(ix0,jy,jz))/ &
                          (VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(One* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/   &
              (SUM(Rho(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)*   &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))    &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            DTU(ix1,jy,jz)= &
                          (One/(Rho(ix1+1,jy,jz,1)+Eps)*VolC(ix1+1,jy,jz)+ &
                           One/(Rho(ix1,jy,jz,1)+Eps)*VolC(ix1,jy,jz))/ &
                          (VolC(ix1+1,jy,jz)+VolC(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy0,jz)= &
                          (One/(Rho(jx,iy0+1,jz,1)+Eps)*VolC(jx,iy0+1,jz)+ &
                           One/(Rho(jx,iy0,jz,1)+Eps)*VolC(jx,iy0,jz))/ &
                          (VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            ThFine= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
              (SUM(Rho(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
               ThFace
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            DTV(jx,iy1,jz)= &
                          (One/(Rho(jx,iy1+1,jz,1)+Eps)*VolC(jx,iy1+1,jz)+ &
                           One/(Rho(jx,iy1,jz,1)+Eps)*VolC(jx,iy1,jz))/ &
                          (VolC(jx,iy1+1,jz)+VolC(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            ThFine= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz0)= &
                          (One/(Rho(jx,jy,iz0+1,1)+Eps)+ &
                           One/(Rho(jx,jy,iz0,1)+Eps))/ &
                          (VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            ThFine= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)) &
               +Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            ThCoarse= &
               SUM(One* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
              (SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1)) &
               +Eps)
            ThFace=(ThCoarse*VolCoarse+ThFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps)
            DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                ThFace
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            DTW(jx,jy,iz1)= &
                          (One/(Rho(jx,jy,iz1+1,1)+Eps)+ &
                           One/(Rho(jx,jy,iz1,1)+Eps))/ &
                          (VolC(jx,jy,iz1+1)+VolC(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO

  IF (TypeW(1:1)=='o') THEN
    IF (BCVel%West/='InFlow') THEN
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=One/(Rho(ix+1,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(One/(Rho(ix,iy,iz,1)+Eps) &
                             +One/(Rho(ix+1,iy,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeE(1:1)=='o') THEN
    IF (BCVel%East/='InFlow') THEN
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=One/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DTU(ix,iy,iz)=Half*(One/(Rho(ix,iy,iz,1)+Eps) &
                             +One/(Rho(ix+1,iy,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    IF (BCVel%South/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=One/(Rho(ix,iy+1,iz,1)+Eps)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(One/(Rho(ix,iy,iz,1)+Eps) &
                             +One/(Rho(ix,iy+1,iz,1)+Eps)) 
        END DO
      END DO
    END IF
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1
    IF (BCVel%North/='InFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=One/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DTV(ix,iy,iz)=Half*(One/(Rho(ix,iy,iz,1)+Eps) &
                             +One/(Rho(ix,iy+1,iz,1)+Eps)) 
                       
        END DO
      END DO
    END IF
  END IF
  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    IF (BCVel%Bottom/='InFlow') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=One/(Rho(ix,iy,iz+1,1)+Eps)
        END DO
      END DO
    ELSE
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*One/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*One/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF
  IF (TypeT(1:1)=='o') THEN
    IF (BCVel%Top/='InFlow') THEN
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)=One/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    ELSE 
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*One/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*One/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        END DO
      END DO
    END IF
  END IF

  DO ix=ix0+1,ix1-1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        DTU(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*One/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix+1,iy,iz)*One/(Rho(ix+1,iy,iz,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix+1,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1-1
      DO iz=iz0+1,iz1
        DTV(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*One/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy+1,iz)*One/(Rho(ix,iy+1,iz,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy+1,iz)+Eps)
      END DO
    END DO
  END DO
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1-1
        DTW(ix,iy,iz)= &
                  (VolC(ix,iy,iz)*One/(Rho(ix,iy,iz,1)+Eps) &
                   +VolC(ix,iy,iz+1)*One/(Rho(ix,iy,iz+1,1)+Eps)) &
                  /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE CellToFaceInvRhoWeight

SUBROUTINE WeightDiagCompute
  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: KappaLoc,Rm,Cpml,Cp_eff
  REAL(RealKind) :: rvLoc,rlLoc,qdLoc,RhoCLoc,qlLoc,Cvml
  REAL(RealKind) :: PotEn,KinEnergy
  REAL(RealKind) :: DpDRho,DpDE,DpDRhoV,DpDRhoL,SoSLoc
  REAL(RealKind) :: ThetaLoc,TLoc,RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc,ELoc,Nume,Deno
  REAL(RealKind) :: Temp
  REAL(RealKind) :: dNumedTheta,dNumedRho
  REAL(RealKind) :: dDenodRho
  REAL(RealKind) :: arg1,darg1dTheta,darg1dRho,ploc,dplocdarg1,cLoc

! IF (RhoCPos>0) THEN
!   DO ix=ix0+1,ix1
!     DO iy=iy0+1,iy1
!       DO iz=iz0+1,iz1
!         RhoVLoc=RhoV(ix,iy,iz,1)
!         RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
!         KappaLoc=(Rd*(Rho(ix,iy,iz,1)-RhoLLoc) &
!                  +(Rv-Rd)*RhoVLoc) &
!                  /(Cpd*Rho(ix,iy,iz,1)+(Cpv-Cpd)*RhoVLoc &
!                  +(Cpl-Cpd)*RhoLLoc+Eps)+Eps
!         DUT(ix,iy,iz)=Rd/(1.0d0-KappaLoc)*(Rd*(theta(ix,iy,iz,1)+Eps)/p0)**(KappaLoc/(1.0d0-KappaLoc))
!         DUR(ix,iy,iz)=0.0d0
!         DTT(ix,iy,iz)=One
!       END DO
!     END DO
!   END DO
!   CALL CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)
!   DUU=One
!   DUV=One
!   DUW=One
  IF (Shallow) THEN
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          DUT(ix,iy,iz)=Grav*theta(ix,iy,iz,1)
          DUR(ix,iy,iz)=0.0d0
          DTT(ix,iy,iz)=One
        END DO
      END DO
    END DO
  ELSE IF (Liquid) THEN
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          RhoLoc=Rho(ix,iy,iz,1)
          ThetaLoc=theta(ix,iy,iz,1)
          IF (RhoLoc>Zero) THEN
            Nume=EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))- alpha_liq/(rho0_liq*kappa_liq*cpl)
            Deno=RhoLoc/Rho0_liq-alpha_liq/(rho0_liq*kappa_liq*cpl)
            dNumedTheta=EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))*(-alpha_liq/RhoLoc)
            dNumedRho=EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))*(alpha_liq*ThetaLoc/RhoLoc/RhoLoc)
            dDenodRho=1.0d0/Rho0_liq
!           pLoc=-1.d0/kappa_liq*LOG(Nume/Deno)+p0_liq
            DUT(ix,iy,iz)=-1.d0/kappa_liq*Deno/Nume*(dNumedTheta*Deno)/Deno/Deno
            DUR(ix,iy,iz)=-1.d0/kappa_liq*Deno/Nume*(dNumedRho*Deno-Nume*dDenodRho)/Deno/Deno
            DTT(ix,iy,iz)=One
          ELSE
            DUT(ix,iy,iz)=Zero
            DUR(ix,iy,iz)=Zero
            DTT(ix,iy,iz)=One
          END IF
!         IF (RhoLoc>Zero) THEN
!           cLoc=alpha_liq/(rho0_liq*kappa_liq*cpl)
!           Nume=(EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))-RhoLoc/Rho0_liq)
!           Deno=RhoLoc/Rho0_liq-cLoc
!           dNumedTheta=EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))*(-alpha_liq/RhoLoc)
!           dNumedRho=EXP(-alpha_liq*(ThetaLoc/RhoLoc-t0_liq))*(alpha_liq*ThetaLoc/RhoLoc/RhoLoc) &
!                  -1.0d0/Rho0_liq
!           dDenodRho=1.0d0/Rho0_liq
!           arg1=Nume/Deno
!           darg1dTheta=dNumedTheta/Deno
!           darg1dRho=(dNumedRho*Deno-Nume*dDenodRho)/Deno/Deno
!           pLoc=arg1*(6.0d0+arg1)/(6.0d0+4.0d0*arg1)
!           dplocdarg1=(6.0d0-2.0d0*arg1)/(6.0d0+4.0d0*arg1)
!           DUT(ix,iy,iz)=dplocdarg1*darg1dTheta
!           DUR(ix,iy,iz)=dplocdarg1*darg1dRho
!         ELSE
!           DUT(ix,iy,iz)=Zero
!           DUR(ix,iy,iz)=Zero
!         END IF
        END DO
      END DO
    END DO
  ELSE IF (LiquidTam) THEN
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
!         p(ix,iy,iz,1)=(Th(ix,iy,iz,1)*(1.0d0-1.0d0/gamma_tam)*Cpl)**(gamma_tam) &
!                       *(1.0d0/(p0+gamma_tam*p0_tam))**(gamma_tam-1.0d0)
          DUT(ix,iy,iz)=gamma_tam*(Theta(ix,iy,iz,1)*(1.0d0-1.0d0/gamma_tam)*Cpl)**(gamma_tam-1.0d0) &
                        *(1.0d0/(p0+gamma_tam*p0_tam))**(gamma_tam-1.0d0)*(1.0d0-1.0d0/gamma_tam)*Cpl
          DUR(ix,iy,iz)=0.0d0
          DTT(ix,iy,iz)=One
        END DO
      END DO
    END DO
  ELSE
    SELECT CASE(ThetaKind)
      CASE('Density') ! Rho*Theta_rho
        CALL CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              ThetaLoc=theta(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
              RhoVLoc=RhoV(ix,iy,iz,1)
              RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
              RhoDLoc=Rho(ix,iy,iz,1)-RhoVLoc-RhoLLoc+Eps
              Rm=Rd+Rv*RhoVLoc/RhoDLoc
              Cpml=Cpd+Cpv*RhoVLoc/RhoDLoc+Cpl*RhoLLoc/RhoDLoc
              KappaLoc=Rm/Cpml
              DUT(ix,iy,iz)=T(ix,iy,iz,1)/(theta(ix,iy,iz,1)+Eps)*(RhoDLoc*Rd+RhoVLoc*Rv)/(One-KappaLoc)
            END DO
          END DO
        END DO
        DTT=One
        DUR=Zero
        DUU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)=FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)/(FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
        DUV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)=FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)/(FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)+Eps)
        DUW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)=FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)/(FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)+Eps)
      CASE('Equiv') ! Rho*Theta_e
        CALL CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              ThetaLoc=theta(ix,iy,iz,1)/Rho(ix,iy,iz,1)
              RhoVLoc=RhoV(ix,iy,iz,1)
              RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
              RhoDLoc=Rho(ix,iy,iz,1)-RhoVLoc-RhoLLoc
              Cp_eff=Cpd+(RhoVLoc+RhoLLoc)/RhoDLoc*Cpl
              Rm=Rd+Rv*RhoVLoc/RhoDLoc
              Cpml=Cpd+Cpv*RhoVLoc/RhoDLoc+Cpl*RhoLLoc/RhoDLoc
              DUT(ix,iy,iz)=T(ix,iy,iz,1)/theta(ix,iy,iz,1) &
                            *(RhoDLoc*Rd+RhoVLoc*Rv)/((Cpml-Rm)/Cp_eff)
            END DO
          END DO
        END DO
        DTT=One
        DUR=Zero
        DUU=One
        DUV=One
        DUW=One
      CASE('Energy') ! Rho*e = E
        SELECT CASE(PreAdv)
        CASE('Outer')       
          CALL CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)
        CASE('Inner')
          CALL CellToFaceRhoWeight(p,Rho,DTU,DTV,DTW,VolC)
        END SELECT
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              RhoVLoc=RhoV(ix,iy,iz,1)
              RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
              RhoDLoc=Rho(ix,iy,iz,1)-RhoVLoc-RhoLLoc
              DUT(ix,iy,iz)=(RhoDLoc*Rd+RhoVLoc*Rv) &
                            /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            END DO
          END DO
        END DO
        DTT=One
        DUR=Zero
        DUU=One
        DUV=One
        DUW=One
      CASE('PreEn') ! p
        SELECT CASE(PreAdv)
        CASE('Outer')
          CALL CellToFaceInvRhoWeight(Rho,DTU,DTV,DTW,VolC)
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              DO iz=iz0+1,iz1
                RhoLoc=Rho(ix,iy,iz,1)
                RhoVLoc=RhoV(ix,iy,iz,1)
                RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
                RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
                TLoc=T(ix,iy,iz,1)
                pLoc=(RhoDLoc*Rd+RhoVLoc*Rv)*TLoc
                Cvml=RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps
                eLoc=(Cvml*TLoc+RhoVLoc*L00)/(RhoLoc+Eps)
                DpDRho=Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                       -Cvd*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 &
                       +eLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
                DpDe=RhoLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
                DpDRhoV=(Rv-Rd)*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                        -L00*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml &
                        -(Cvv-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2
                DpDRhoL=-Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                        -(Cpl-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2
                SoSLoc=SQRT(ABS(DpDRho+pLoc/(RhoLoc+Eps)**2*DpDe+RhoVLoc/(RhoLoc+Eps)*DpDRhoV &
                               +RhoLLoc/(RhoLoc+Eps)*DpDRhoL))
                DTT(ix,iy,iz)=RhoLoc*SoSLoc**2
              END DO
            END DO
          END DO
        CASE('Inner')
          CALL CellToFaceRhoWeight(Sound,Rho,DTU,DTV,DTW,VolC)
          DTT=One
        END SELECT
        DUT=One
        DUR=Zero
        DUU=One
        DUV=One
        DUW=One
      CASE('Anelastic')
        DTU=One 
        DTV=One 
        DTW=One 
        CALL CellToFaceWeight(Rho,DUU,DUV,DUW,VolC)
        DTT=One
        DUT=One  
        DUR=Zero  
      CASE('PseudoIn')
        CALL CellToFaceRhoWeight(ThetaProf,Rho,DTU,DTV,DTW,VolC)
        CALL CellToFaceWeight(ThetaProf,DUU,DUV,DUW,VolC)
        DTT=One
        DUT=One
        DUR=One
     CASE DEFAULT
       DO ix=ix0+1,ix1
         DO iy=iy0+1,iy1
           DO iz=iz0+1,iz1
             DUT(ix,iy,iz)=Rd/(1.0d0-kappa)*(Rd*(theta(ix,iy,iz,1)+Eps)/p0)**(kappa/(1.0d0-kappa))
           END DO
         END DO
       END DO
       CALL CellToFaceRhoWeight(Theta,Rho,DTU,DTV,DTW,VolC)
       DUU=One
       DUV=One
       DUW=One
       DUR=Zero
       DTT=One
     END SELECT ! (ThetaKind)

  END IF

END SUBROUTINE WeightDiagCompute

SUBROUTINE VCycleSet(ibLoc)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibLoc
  INTEGER :: nx,ny,nz

  ib=LocGlob(ibLoc)
  CALL Set(Floor(ib))
  NULLIFY(A)
  NULLIFY(AMat)
  NULLIFY(AMumps)
  NULLIFY(AR)
  NULLIFY(ATB)
  IF (ASSOCIATED(MaActual(ibLoc)%MatrixMumps)) THEN
    AMumps=>MaActual(ibLoc)%MatrixMumps
    xS=>xx(ibLoc)%p(1,:,:,:)
    bS=>bb(ibLoc)%p(1,:,:,:)
    nx=SIZE(xS,1)
    ny=SIZE(xS,2)
    nz=SIZE(xS,3)
  ELSEIF (ASSOCIATED(MaActual(ibLoc)%MatrixMu)) THEN
    iRef=1
    A=>MaActual(ibLoc)%MatrixMu
    xS=>xx(ibLoc)%p(1,:,:,:)
    bS=>bb(ibLoc)%p(1,:,:,:)
    nx=SIZE(xS,1)
    ny=SIZE(xS,2)
    nz=SIZE(xS,3)
    CALL ReverseIndices(nx,ny,nz)
  ELSE IF (ASSOCIATED(MaActual(ibLoc)%MatrixMuMat)) THEN
    iRef=1
    AMat=>MaActual(ibLoc)%MatrixMuMat
    xSB=>xx(ibLoc)%p(:,:,:,:)
    bSB=>bb(ibLoc)%p(:,:,:,:)
    nx=SIZE(xSB,2)
    ny=SIZE(xSB,3)
    nz=SIZE(xSB,4)
  ELSE IF (ASSOCIATED(MaActual(ibLoc)%MatrixT)) THEN
    xS=>xx(ibLoc)%p(1,:,:,:)
    bS=>bb(ibLoc)%p(1,:,:,:)
    nx=SIZE(xS,1)
    ny=SIZE(xS,2)
    nz=SIZE(xS,3)
    AR=>MaActual(ibLoc)%MatrixT
    CALL ReverseIndices(nx,ny,nz)
  ELSE IF (ASSOCIATED(MaActual(ibLoc)%MatrixTB)) THEN
    xSB=>xx(ibLoc)%p(:,:,:,:)
    bSB=>bb(ibLoc)%p(:,:,:,:)
    nx=SIZE(xSB,2)
    ny=SIZE(xSB,3)
    nz=SIZE(xSB,4)
    ATB=>MaActual(ibLoc)%MatrixTB
    CALL ReverseIndices(nx,ny,nz)
  ELSE IF (ASSOCIATED(MaActual(ibLoc)%MatrixTR)) THEN
    xSB=>xx(ibLoc)%p(:,:,:,:)
    bSB=>bb(ibLoc)%p(:,:,:,:)
    nx=SIZE(xSB,2)
    ny=SIZE(xSB,3)
    nz=SIZE(xSB,4)
    ATR=>MaActual(ibLoc)%MatrixTR
    CALL ReverseIndices(nx,ny,nz)
  END IF

END SUBROUTINE VCycleSet

SUBROUTINE VCycleCompute(xD,bD)

  REAL(RealKind) :: xD(:,:,:),bd(:,:,:)

  xD=0.0d0
  CALL Solve(AR,bD,xD)

END SUBROUTINE VCycleCompute

SUBROUTINE VCycleCG(MaxIter,TolRel,TolAbs)

  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: TolRel,TolAbs

  REAL(RealKind), TARGET :: xCG(2,n1,n2,n3)
  REAL(RealKind) :: r(2,n1,n2,n3)
  REAL(RealKind), TARGET :: rhs(nx*ny*nz)

  INTEGER :: Iter,MultIter=1

  IF (ASSOCIATED(AMumps)) THEN
    rhs=RESHAPE(xS,(/nx*ny*nz/))
    AMumps%RHS=>rhs  !OSSI OSSI OSSI
    AMumps%JOB=3
    AMumps%ICNTL(14)=20
    CALL DMUMPS(AMumps)
    xS=RESHAPE(rhs,(/nx,ny,nz/))
  ELSE IF (ASSOCIATED(A)) THEN
    CALL Reverse(xCG(1,:,:,:),xS)
    r(1,:,:,:)=xCG(1,:,:,:)
    CALL VCycleCGMu(xCG(1,:,:,:),MaxIter,TolRel,TolAbs)
    CALL ReverseIndices(nx,ny,nz)
    CALL ReverseBack(xS,xCG(1,:,:,:))
  ELSE IF (ASSOCIATED(AMat)) THEN
    CALL VCycleCGMuTR(MaxIter,TolRel)
  ELSE IF (ASSOCIATED(AR)) THEN
    CALL Reverse(r(1,:,:,:),xS)
    xCG=0.0e0
    CALL Solve(AR,r(1,:,:,:),xCG(1,:,:,:))
!   DO Iter=1,MultIter
!     CALL Multigrid(AR,r(1,:,:,:),xCG(1,:,:,:))
!   END DO
    CALL ReverseBack(xS,xCG(1,:,:,:))
  ELSE IF (ASSOCIATED(ATB)) THEN
    CALL ReverseB(r,xSB)
    xCG=0.0e0
    DO Iter=1,MultIter
      CALL Multigrid(ATB,r,xCG)
    END DO
    CALL ReverseBBack(xSB,xCG)
  ELSE IF (ASSOCIATED(ATR)) THEN
    CALL ReverseB(r,xSB)
    xCG=0.0e0
    DO Iter=1,MultIter
      CALL Multigrid(ATR,r,xCG)
    END DO
    CALL ReverseBBack(xSB,xCG)
  END IF
END SUBROUTINE VCycleCG 

SUBROUTINE VCycleCGM(MaxIter,Tol)


  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: Tol

  REAL(RealKind) ::  xCG(n1,n2,n3),r(n1,n2,n3),p(n1,n2,n3),z(n1,n2,n3),q(n1,n2,n3) 
  REAL(RealKind) ::  VolRev(n1,n2,n3)
  REAL(RealKind) :: alpha,beta,rho,rhoOld
  REAL(RealKind) :: NormRes,NormRes0
  INTEGER :: Iter
  INTEGER :: n

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0

  REAL(RealKind) :: DDOT,DNRM2

! Compute initial Residual r=b-A*x
  CALL Reverse(r,xS)
  CALL Reverse(VolRev,VolB)
  xCG=0.0e0


  NormRes0=SQRT(SUM((r/(VolReV+Eps))*(r/(VolReV+Eps))))

! IF (NormRes0<=1.0e-14) THEN
  IF (NormRes0<=0.0d0) THEN
    MaxIter=0
    Tol=NormRes0
  ELSE
    DO Iter=1,MaxIter
      CALL VCycleCompute(z,r)
      rhoOld=rho
!     rho=z*r
      rho=SUM(z*r)
      IF (Iter==1) THEN
        p=z
      ELSE
!       p=z+beta*p
        beta=rho/rhoOld
        p=z+beta*p
      END IF
!     Compute q=A*p    
      CALL MatVecT(q(:,:,:),AR,p(:,:,:))
!     alpha=rho/(p*q)
      alpha=rho/SUM(p*q)
      xCG=xCG+alpha*p
      r=r-alpha*q
      NormRes=SQRT(SUM((r/(VolRev+Eps))*(r/(VolRev+Eps))))
!     Convergence check
!     WRITE(*,*) MyId,Iter,NormRes,NormRes0,Tol
      IF (NormRes<=Tol*NormRes0.OR.NormRes<=1.e-12) THEN
        MaxIter=Iter
        Tol=NormRes/NormRes0
        EXIT
      END IF
    END DO   
  END IF
  CALL ReverseBack(xS,xCG)
END SUBROUTINE VCycleCGM

SUBROUTINE VCycleCGMu1(MaxIter,TolRel,TolAbs)

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: TolRel,TolAbs


  REAl(RealKind), POINTER :: xCG(:,:,:)
  REAL(RealKind) ::  r(nx,ny,nz),p(nx,ny,nz),z(nx,ny,nz),q(nx,ny,nz) 
  REAL(RealKind) :: alpha,beta,rho,rhoOld
  REAL(RealKind) :: NormRes,NormRes0
  INTEGER :: Iter
  INTEGER :: n

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0

  REAL(RealKind) :: DDOT,DNRM2

  xCG=>xS

! Compute initial Residual r=b-A*x
  r=xCG
  xCG=0.0e0

! NormRes0=SQRT(SUM((r/(VolB+Eps))*(r/(VolB+Eps))))
  NormRes0=SQRT(SUM(r*r))

  IF (NormRes0<=1.0e-10) THEN
    MaxIter=0
    TolRel=NormRes0
  ELSE
    DO Iter=1,MaxIter
      CALL VCycleComputeMu(z,r)
      rhoOld=rho
!     rho=z*r
      rho=SUM(z*r)
      IF (Iter==1) THEN
        p=z
      ELSE
!       p=z+beta*p
        beta=rho/rhoOld
        p=z+beta*p
      END IF
!     Compute q=A*p    
      CALL SpAVec(q(:,:,:),A(1),p(:,:,:))
!     CALL SpAVec_SpRowColL(q(:,1,1),A(1),p(:,1,1))
!     alpha=rho/(p*q)
      alpha=rho/SUM(p*q)
      xCG=xCG+alpha*p
      r=r-alpha*q
!     NormRes=SQRT(SUM((r/(VolB+Eps))*(r/(VolB+Eps))))
      NormRes=SQRT(SUM(r*r))
!     Convergence check
!     WRITE(*,*) Iter,NormRes/NormRes0,NormRes0
      IF (NormRes<=TolRel*NormRes0.OR.NormRes<=1.e-10) THEN
        MaxIter=Iter
        TolRel=NormRes/NormRes0
        EXIT
      END IF
    END DO   
  END IF

END SUBROUTINE VCycleCGMu1

SUBROUTINE VCycleCGMu(x,MaxIter,TolRel,TolAbs)

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: TolRel,TolAbs


  REAl(RealKind), TARGET :: x(:,:,:)
  

  REAL(RealKind) :: r(n1,n2,n3),rTilde(n1,n2,n3),rHat(n1,n2,n3),p(n1,n2,n3) &
                   ,pHat(n1,n2,n3),t(n1,n2,n3),v(n1,n2,n3) 
  REAL(RealKind) :: alpha,beta,omega,omegaOld,rho,rhoOld
  REAL(RealKind) :: rNorm,rNorm0
  INTEGER :: Iter

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0
! x=>xS
! -- Compute initial Residual r=b-A*x --
  r=x
  x=0.0d0

! -- Choose rTilde, for example rTilde=r --
!  CALL Copy(r,rTilde)
  rTilde=r
! CALL Random(rTilde)

  omegaOld=0.0e0
  rNorm0=SQRT(Sum(r*r))
  DO Iter=1,MaxIter
    rhoOld=rho
    rho=SUM(rTilde*r)
    IF (rho==0.0e0) THEN
      MaxIter=Iter
      TolRel=0.0d0
      EXIT
    END IF
    IF (Iter==1) THEN
!      CALL Copy(r,p)
      p=r
    ELSE
      beta=(rho/rhoOld)*(alpha/omega)
      p=r+beta*(p-omega*v) 
    END IF
! -- Solve pHat from preconditioning
    CALL VCycleComputeMu(pHat,p)
    CALL SpAVec(v(:,1,1),A(1),pHat(:,1,1))
    alpha=rho/SUM(rTilde*v)
    r=r-alpha*v
    rNorm=SQRT(SUM(r*r))
    IF (MyId==0) THEN
!     WRITE(TermUnit,*) 'BiCGS T',rNorm,TolRel,rNorm0,SQRT(SUM(x*x))
    END IF
    IF (rNorm<=TolRel*rNorm0+TolAbs) THEN
      x=alpha*pHat+x
      MaxIter=-Iter
      TolRel=rNorm
      EXIT
    END IF
    CALL VCycleComputeMu(rHat,r)
    CALL SpAVec(t(:,1,1),A(1),rHat(:,1,1))
    omegaOld=omega
    omega=SUM(t*r)/SUM(t*t)
    x=alpha*pHat+omega*rHat+x
    r=-omega*t+r
    rNorm=SQRT(SUM(r*r))
    IF (MyId==0) THEN
!     WRITE(TermUnit,*) 'BiCGS T',rNorm,TolRel,rNorm0,SQRT(SUM(x*x))
    END IF
    IF (rNorm<=TolRel*rNorm0+TolAbs) THEN
      MaxIter=Iter
      TolRel=rNorm
      EXIT
    END IF
    IF (Iter==MaxIter) THEN
      MaxIter=2*Iter
      TolRel=rNorm
      EXIT
    END IF
  END DO
END SUBROUTINE VCycleCGMu

SUBROUTINE VCycleCGMuTR(MaxIter,Tol)

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: Tol


  REAl(RealKind), POINTER :: x(:,:,:,:)
  

  REAL(RealKind) :: r(2,nx,ny,nz),rTilde(2,nx,ny,nz),rHat(2,nx,ny,nz),p(2,nx,ny,nz) &
                   ,pHat(2,nx,ny,nz),t(2,nx,ny,nz),v(2,nx,ny,nz) 
  REAL(RealKind) :: alpha,beta,omega,omegaOld,rho,rhoOld
  REAL(RealKind) :: rNorm,rNorm0
  INTEGER :: Iter

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0
  x=>xSB
! -- Compute initial Residual r=b-A*x --
  r=x
  x=0.0d0

! -- Choose rTilde, for example rTilde=r --
!  CALL Copy(r,rTilde)
  rTilde=r
! CALL Random(rTilde)

  omegaOld=0.0e0
  rNorm0=SQRT(Sum(r*r))
  DO Iter=1,MaxIter
    rhoOld=rho
    rho=SUM(rTilde*r)
    IF (rho==0.0e0) THEN
      MaxIter=Iter
      Tol=0.0d0
      EXIT
    END IF
    IF (Iter==1) THEN
!      CALL Copy(r,p)
      p=r
    ELSE
      beta=(rho/rhoOld)*(alpha/omega)
      p=r+beta*(p-omega*v) 
    END IF
! -- Solve pHat from preconditioning
    CALL VCycleComputeMuTR(pHat,p)
    CALL SpAVec(v(:,:,:,:),AMat(1),pHat(:,:,:,:))
    alpha=rho/SUM(rTilde*v)
    r=r-alpha*v
    rNorm=SQRT(SUM(r*r))
    IF (MyId==0) THEN
!     WRITE(TermUnit,*) 'BiCGS TR',rNorm,Tol,rNorm0
    END IF
    IF (rNorm<=Tol*rNorm0) THEN
      x=alpha*pHat+x
      MaxIter=-Iter
      Tol=rNorm
      EXIT
    END IF
    CALL VCycleComputeMuTR(rHat,r)
    CALL SpAVec(t(:,:,:,:),AMat(1),rHat(:,:,:,:))
    omegaOld=omega
    omega=SUM(t*r)/SUM(t*t)
    x=alpha*pHat+omega*rHat+x
    r=-omega*t+r
    rNorm=SQRT(SUM(r*r))
    IF (MyId==0) THEN
!     WRITE(TermUnit,*) 'BiCGS TR',rNorm,Tol,rNorm0
    END IF
    IF (rNorm<=Tol*rNorm0) THEN
      MaxIter=Iter
      Tol=rNorm
      EXIT
    END IF
    IF (Iter==MaxIter) THEN
      MaxIter=2*Iter
      Tol=rNorm
      EXIT
    END IF
  END DO
END SUBROUTINE VCycleCGMuTR
SUBROUTINE VCycleCGMuTR1(MaxIter,Tol)

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: MaxIter
  REAL(RealKind), INTENT(INOUT) :: Tol


  REAl(RealKind), POINTER :: xCG(:,:,:,:)
  REAL(RealKind) ::  r(2,nx,ny,nz),p(2,nx,ny,nz),z(2,nx,ny,nz),q(2,nx,ny,nz) 
  REAL(RealKind) :: alpha,beta,rho,rhoOld
  REAL(RealKind) :: NormRes,NormRes0
  INTEGER :: Iter
  INTEGER :: n

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0

  REAL(RealKind) :: DDOT,DNRM2


  xCG=>xSB

! Compute initial Residual r=b-A*x
  r=xCG
  xCG=0.0e0

! NormRes0=SQRT(SUM((r/(VolB+Eps))*(r/(VolB+Eps))))
  NormRes0=SQRT(SUM(r*r))

  IF (NormRes0<=1.0e-10) THEN
    MaxIter=0
    Tol=NormRes0
  ELSE
    DO Iter=1,MaxIter
      CALL VCycleComputeMuTR(z,r)
      rhoOld=rho
!     rho=z*r
      rho=SUM(z*r)
      IF (Iter==1) THEN
        p=z
      ELSE
!       p=z+beta*p
        beta=rho/rhoOld
        p=z+beta*p
      END IF
!     Compute q=A*p    
      CALL SpAVec(q(:,:,:,:),AMat(1),p(:,:,:,:))
!     CALL SpAVec_SpRowColL(q(:,1,1),A(1),p(:,1,1))
!     alpha=rho/(p*q)
      alpha=rho/SUM(p*q)
      xCG=xCG+alpha*p
      r=r-alpha*q
!     NormRes=SQRT(SUM((r/(VolB+Eps))*(r/(VolB+Eps))))
      NormRes=SQRT(SUM(r*r))
!     Convergence check
!      WRITE(*,*) Iter,NormRes/NormRes0,NormRes0
      IF (NormRes<=Tol*NormRes0.OR.NormRes<=1.e-10) THEN
        MaxIter=Iter
        Tol=NormRes/NormRes0
        EXIT
      END IF
    END DO   
  END IF

END SUBROUTINE VCycleCGMuTR1

RECURSIVE SUBROUTINE VCycleComputeMu(xD,bD)

  REAL(RealKind), TARGET, OPTIONAL :: xD(n1,n2,n3)
  REAL(RealKind), TARGET, OPTIONAL :: bD(n1,n2,n3)

  REAL(RealKind), TARGET :: xC(1:(n1+1)/2,1:(n2+1)/2,1:(n3+1)/2) 
  REAL(RealKind), TARGET :: bC(1:(n1+1)/2,1:(n2+1)/2,1:(n3+1)/2) 
  REAL(RealKind) :: r(1:n1,1:n2,1:n3) 
  REAl(RealKind), POINTER :: xF(:,:,:)
  REAl(RealKind), POINTER :: bF(:,:,:)
  TYPE(SpRowCol), POINTER :: AF

  INTEGER :: i

  INTEGER :: n1C,n2C,n3C
  INTEGER :: n1F,n2F,n3F


  n1F=n1
  n2F=n2
  n3F=n3
  IF (PRESENT(xD)) THEN
    xF=>xD
  ELSE
    xF=>xSM 
  END IF
  IF (PRESENT(bD)) THEN
    bF=>bD
  ELSE IF (iRef==1) THEN
    ALLOCATE(bFine(1:n1,1:n2,1:n3))
    bF=>bFine
    bF=bSM
  ELSE
    bF=>bSM
  END IF
  AF=>A(iRef)

! Relax once
! CALL GaussSeidelF(xF(:,:,:),AF,bF(:,:,:),iRef)
  xF=0.0d0
  CALL IC0Solve(xF(:,1,1),AF,bF(:,1,1),iRef)

  IF (iRef<RefLevel) THEN
!   Compute the Residual
    CALL SpAVec(r(:,1,1),AF,xF(:,1,1))
    r=bF-r  
    iRef=iRef+1
    n1C=(n1+1)/2
    n2C=(n2+1)/2
    n3C=(n3+1)/2

!   Restrict
    CALL RestrictionMu(n1C,n2C,n3C,bC,n1F,n2F,n3F,r)
    xSM=>xC
    bSM=>bC
    n1=n1C
    n2=n2C
    n3=n3C
    CALL VCycleComputeMu
    n1=n1F
    n2=n2F
    n3=n3F
    iRef=iRef-1
!   Prolongate  xF=xF+xC
    CALL ProlongationMu(n1F,n2F,n3F,xF,n1C,n2C,n3C,xC)

  END IF

! Relax once
! CALL GaussSeidelB(xF(:,:,:),AF,bF(:,:,:),iRef)
  CALL IC0Solve(xF(:,1,1),AF,bF(:,1,1),iRef)
  IF (iRef==1.AND..NOT.Present(bD)) THEN
    DEALLOCATE(bFine)
  END IF

END SUBROUTINE VCycleComputeMu

RECURSIVE SUBROUTINE VCycleComputeMuTR(xD,bD)

  REAL(RealKind), TARGET, OPTIONAL :: xD(2,nx,ny,nz)
  REAL(RealKind), TARGET, OPTIONAL :: bD(2,nx,ny,nz)

  REAL(RealKind), TARGET :: xC(2,1:(nx+1)/2,1:(ny+1)/2,1:(nz+1)/2) 
  REAL(RealKind), TARGET :: bC(2,1:(nx+1)/2,1:(ny+1)/2,1:(nz+1)/2) 
  REAL(RealKind) :: r(2,1:nx,1:ny,1:nz) 
  REAl(RealKind), POINTER :: xF(:,:,:,:)
  REAl(RealKind), POINTER :: bF(:,:,:,:)
  TYPE(SpRowColMat), POINTER :: AF

  INTEGER :: i

  INTEGER :: nxC,nyC,nzC
  INTEGER :: nxF,nyF,nzF


  nxF=nx
  nyF=ny
  nzF=nz
  IF (PRESENT(xD)) THEN
    xF=>xD
  ELSE
    xF=>xSB 
  END IF
  IF (PRESENT(bD)) THEN
    bF=>bD
  ELSE IF (iRef==1) THEN
    ALLOCATE(bFineB(2,1:nx,1:ny,1:nz))
    bF=>bFineB
    bF=bSB
  ELSE
    bF=>bSB
  END IF
  AF=>AMat(iRef)

! Relax once
! CALL GaussSeidelF(xF(:,:,:),AF,bF(:,:,:),iRef)
  xF=0.0d0
  CALL IC0Solve(xF(:,:,:,:),AF,bF(:,:,:,:),iRef)

  IF (iRef<RefLevel) THEN
!   Compute the Residual
    CALL SpAVec(r(:,:,:,:),AF,xF(:,:,:,:))
    r=bF-r  
    iRef=iRef+1
    nxC=(nx+1)/2
    nyC=(ny+1)/2
    nzC=(nz+1)/2

!   Restrict
    CALL RestrictionMuTR(nxC,nyC,nzC,bC,nxF,nyF,nzF,r)
    xSB=>xC
    bSB=>bC
    nx=nxC
    ny=nyC
    nz=nzC
    CALL VCycleComputeMuTR
    nx=nxF
    ny=nyF
    nz=nzF
    iRef=iRef-1
!   Prolongate  xF=xF+xC
    CALL ProlongationMuTR(nxF,nyF,nzF,xF,nxC,nyC,nzC,xC)

  END IF

! Relax once
! CALL GaussSeidelB(xF(:,:,:),AF,bF(:,:,:),iRef)
  CALL IC0Solve(xF(:,:,:,:),AF,bF(:,:,:,:),iRef)
  IF (iRef==1.AND..NOT.Present(bD)) THEN
    DEALLOCATE(bFineB)
  END IF

END SUBROUTINE VCycleComputeMuTR



END MODULE JacAccGrav_Mod

