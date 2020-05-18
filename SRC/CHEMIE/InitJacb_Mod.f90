MODULE InitEx_Mod

  USE Floor_Mod
  USE DataType_Mod
  USE Distribution_Mod
  USE Aerosol_Mod
  USE InitAerosol_Mod, ONLY: InitGrid       &
                            ,InitGas        &
                            ,SetIndices      
!                           ,InitAerosol
  USE Kernel_Mod, KernelF=>FUCH

  IMPLICIT NONE

  REAL(RealKind), PARAMETER :: r1=1.0d-3 ! /mu m
 
! Jacobson zum beispiel:mass in mu g ;diameter in mu m
 
  REAL(RealKind) :: M_EC1=67.0d0,     DGV_EC1=10.7d-2,   SIGMAG_EC1=1.66d0
  REAL(RealKind) :: M_OC1=47.6d0,     DGV_OC1=13.6d-2,   SIGMAG_OC1=1.50d0
  REAL(RealKind) :: M_S1E=8.0d-3,     DGV_S1E=4.0d-2,    SIGMAG_S1E=1.8d0
  REAL(RealKind) :: M_S1I=7.2d-2,     DGV_S1I=4.0d-2,    SIGMAG_S1I=1.8d0
  REAL(RealKind) :: M_W1E=1.2d-2,     DGV_W1E=4.0d-2,    SIGMAG_W1E=1.8d0
  REAL(RealKind) :: M_W1I=10.0d-2,    DGV_W1I=4.0d-2,    SIGMAG_W1I=1.8d0
  REAL(RealKind) :: M_N1=3.0d-3,      DGV_N1=4.0d-2,     SIGMAG_N1=1.8d0

  REAL(RealKind) :: M_EC2=0.0d0  ,    DGV_EC2=0.0d0,    SIGMAG_EC2=0.0d0
  REAL(RealKind) :: M_OC2=0.0d0,      DGV_OC2=0.0d0,    SIGMAG_OC2=0.0d0
  REAL(RealKind) :: M_S2E=4.9d-1,     DGV_S2E=3.2d-1,   SIGMAG_S2E=2.16d0
  REAL(RealKind) :: M_S2I=4.43d0,     DGV_S2I=3.2d-1,   SIGMAG_S2I=2.16d0
  REAL(RealKind) :: M_W2E=7.4d-1,     DGV_W2E=3.2d-1,   SIGMAG_W2E=2.16d0
  REAL(RealKind) :: M_W2I=5.0d0  ,    DGV_W2I=3.2d-1,   SIGMAG_W2I=2.16d0
  REAL(RealKind) :: M_N2=2.0d0,       DGV_N2=3.2d-1,    SIGMAG_N2=2.16d0
 
  REAL(RealKind) :: M_EC3=0.0d0,      DGV_EC3=0.0d0,     SIGMAG_EC3=0.0d0
  REAL(RealKind) :: M_OC3=0.0d0,      DGV_OC3=0.0d0,     SIGMAG_OC3=0.0d0
  REAL(RealKind) :: M_S3E=5.0d-2,     DGV_S3E=12.3d0,    SIGMAG_S3E=2.4d0
  REAL(RealKind) :: M_S3I=4.5d-1,     DGV_S3I=12.3d0,    SIGMAG_S3I=2.4d0
  REAL(RealKind) :: M_W3E=7.5d-2,     DGV_W3E=12.3d0,    SIGMAG_W3E=2.4d0
  REAL(RealKind) :: M_W3I=2.0d0,      DGV_W3I=12.3d0,    SIGMAG_W3I=2.4d0
  REAL(RealKind) :: M_N3=2.5d0,       DGV_N3=8.1d0,      SIGMAG_N3=2.4d0

  REAL(RealKind), PRIVATE :: GridFac=0.5d0
  REAL(RealKind) :: TimeOutput=0.0d0,OutputTime=1.0d1
  INTEGER :: iOut=0
  INTEGER, PARAMETER :: nsEx=7
  LOGICAL :: CondProc=.FALSE. 
  LOGICAL :: KoagProc=.TRUE. 

CONTAINS 

SUBROUTINE InitAerosol(VecT,FileName)

  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  CHARACTER(*) :: FileName

  INTEGER :: i,is
  REAL(RealKind) :: mL,mU
  REAL(RealKind) :: rL,rU
  TYPE(Vec4_T), POINTER :: cAero(:)

  REAL(RealKind) :: Num_EC1, DGN_EC1
  REAL(RealKind) :: Num_OC1, DGN_OC1
  REAL(RealKind) :: Num_S1E, DGN_S1E
  REAL(RealKind) :: Num_S1I, DGN_S1I
  REAL(RealKind) :: Num_W1E, DGN_W1E
  REAL(RealKind) :: Num_W1I, DGN_W1I
  REAL(RealKind) :: Num_N1,  DGN_N1
 
  REAL(RealKind) :: Num_EC2, DGN_EC2
  REAL(RealKind) :: Num_OC2, DGN_OC2
  REAL(RealKind) :: Num_S2E, DGN_S2E
  REAL(RealKind) :: Num_S2I, DGN_S2I
  REAL(RealKind) :: Num_W2E, DGN_W2E
  REAL(RealKind) :: Num_W2I, DGN_W2I
  REAL(RealKind) :: Num_N2,  DGN_N2
 
  REAL(RealKind) :: Num_EC3, DGN_EC3
  REAL(RealKind) :: Num_OC3, DGN_OC3
  REAL(RealKind) :: Num_S3E, DGN_S3E
  REAL(RealKind) :: Num_S3I, DGN_S3I
  REAL(RealKind) :: Num_W3E, DGN_W3E
  REAL(RealKind) :: Num_W3I, DGN_W3I
  REAL(RealKind) :: Num_N3,  DGN_N3


  M_EC1=M_EC1*1.d-9
  M_OC1=M_OC1*1.d-9
  M_S1E=M_S1E*1.d-9
  M_S1I=M_S1I*1.d-9
  M_W1E=M_W1E*1.d-9
  M_W1I=M_W1I*1.d-9
  M_N1=M_N1*1.d-9

  M_EC2=M_EC2*1.d-9
  M_OC2=M_OC2*1.d-9
  M_S2E=M_S2E*1.d-9
  M_S2I=M_S2I*1.d-9
  M_W2E=M_W2E*1.d-9
  M_W2I=M_W2I*1.d-9
  M_N2=M_N2*1.d-9

  M_EC3=M_EC3*1.d-9
  M_OC3=M_OC3*1.d-9
  M_S3E=M_S3E*1.d-9
  M_S3I=M_S3I*1.d-9
  M_W3E=M_W3E*1.d-9
  M_W3I=M_W3I*1.d-9
  M_N3=M_N3*1.d-9

  DGV_EC1=DGV_EC1*1.d-6
  DGV_OC1=DGV_OC1*1.d-6
  DGV_S1E=DGV_S1E*1.d-6
  DGV_S1I=DGV_S1I*1.d-6
  DGV_W1E=DGV_W1E*1.d-6
  DGV_W1I=DGV_W1I*1.d-6
  DGV_N1=DGV_N1*1.d-6

  DGV_EC2=DGV_EC2*1.d-6
  DGV_OC2=DGV_OC2*1.d-6
  DGV_S2E=DGV_S2E*1.d-6
  DGV_S2I=DGV_S2I*1.d-6
  DGV_W2E=DGV_W2E*1.d-6
  DGV_W2I=DGV_W2I*1.d-6
  DGV_N2=DGV_N2*1.d-6

  DGV_EC3=DGV_EC3*1.d-6
  DGV_OC3=DGV_OC3*1.d-6
  DGV_S3E=DGV_S3E*1.d-6
  DGV_S3I=DGV_S3I*1.d-6
  DGV_W3E=DGV_W3E*1.d-6
  DGV_W3I=DGV_W3I*1.d-6
  DGV_N3=DGV_N3*1.d-6

  DGN_EC1=DGN(DGV_EC1,SIGMAG_EC1)
  DGN_OC1=DGN(DGV_OC1,SIGMAG_OC1)
  DGN_S1E=DGN(DGV_S1E,SIGMAG_S1E)
  DGN_S1I=DGN(DGV_S1I,SIGMAG_S1I)
  DGN_W1E=DGN(DGV_W1E,SIGMAG_W1E)
  DGN_W1I=DGN(DGV_W1I,SIGMAG_W1I)
  DGN_N1=DGN(DGV_N1,SIGMAG_N1)
 
  DGN_EC2=DGN(DGV_EC2,SIGMAG_EC2)
  DGN_OC2=DGN(DGV_OC2,SIGMAG_OC2)
  DGN_S2E=DGN(DGV_S2E,SIGMAG_S2E)
  DGN_S2I=DGN(DGV_S2I,SIGMAG_S2I)
  DGN_W2E=DGN(DGV_W2E,SIGMAG_W2E)
  DGN_W2I=DGN(DGV_W2I,SIGMAG_W2I)
  DGN_N2=DGN(DGV_N2,SIGMAG_N2)
 
  DGN_EC3=DGN(DGV_EC3,SIGMAG_EC3)
  DGN_OC3=DGN(DGV_OC3,SIGMAG_OC3)
  DGN_S3E=DGN(DGV_S3E,SIGMAG_S3E)
  DGN_S3I=DGN(DGV_S3I,SIGMAG_S3I)
  DGN_W3E=DGN(DGV_W3E,SIGMAG_W3E)
  DGN_W3I=DGN(DGV_W3I,SIGMAG_W3I)
  DGN_N3=DGN(DGV_N3,SIGMAG_N3)
 
  Num_EC1=Mass2Number(M_EC1,DGN_EC1,SIGMAG_EC1,1.0d3)
  Num_OC1=Mass2Number(M_OC1,DGN_OC1,SIGMAG_OC1,1.0d3)
  Num_S1E=Mass2Number(M_S1E,DGN_S1E,SIGMAG_S1E,1.0d3)
  Num_S1I=Mass2Number(M_S1I,DGN_S1I,SIGMAG_S1I,1.0d3)
  Num_W1E=Mass2Number(M_W1E,DGN_W1E,SIGMAG_W1E,1.0d3)
  Num_W1I=Mass2Number(M_W1I,DGN_W1I,SIGMAG_W1I,1.0d3)
  Num_N1=Mass2Number(M_N1,DGN_N1,SIGMAG_N1,1.0d3)

  Num_EC2=Mass2Number(M_EC2,DGN_EC2,SIGMAG_EC2,1.0d3)
  Num_OC2=Mass2Number(M_EC2,DGN_OC2,SIGMAG_OC2,1.0d3)
  Num_S2E=Mass2Number(M_S2E,DGN_S2E,SIGMAG_S2E,1.0d3)
  Num_S2I=Mass2Number(M_S2E,DGN_S2I,SIGMAG_S2I,1.0d3)
  Num_W2E=Mass2Number(M_W2E,DGN_W2E,SIGMAG_W2E,1.0d3)
  Num_W2I=Mass2Number(M_W2I,DGN_W2I,SIGMAG_W2I,1.0d3)
  Num_N2=Mass2Number(M_N2,DGN_N2,SIGMAG_N2,1.0d3)
 
  Num_EC3=Mass2Number(M_EC3,DGN_EC3,SIGMAG_EC3,1.0d3)
  Num_OC3=Mass2Number(M_OC3,DGN_OC3,SIGMAG_OC3,1.0d3)
  Num_S3E=Mass2Number(M_S3E,DGN_S3E,SIGMAG_S3E,1.0d3)
  Num_S3I=Mass2Number(M_S3I,DGN_S3I,SIGMAG_S3I,1.0d3)
  Num_W3E=Mass2Number(M_W3E,DGN_W3E,SIGMAG_W3E,1.0d3)
  Num_W3I=Mass2Number(M_W3I,DGN_W3I,SIGMAG_W3I,1.0d3)
  Num_N3=Mass2Number(M_N3,DGN_N3,SIGMAG_N3,1.0d3)

  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    cAero=>VecT(ibLoc)%Vec
    DO i=1,nFrac
      mL=m(i)
      mU=m(i+1)
      rL=(mL/(4.0d0/3.0d0*Pi*1.0d3))**(1.0d0/3.0d0)
      rU=(mU/(4.0d0/3.0d0*Pi*1.0d3))**(1.0d0/3.0d0)
      cAero(1)%cInt(:,:,:,i)=Moment0(rL,rU,Num_EC1,DGN_EC1,SIGMAG_EC1) &
            +Moment0(rL,rU,Num_OC1,DGN_OC1,SIGMAG_OC1) &
            +Moment0(rL,rU,Num_S1E,DGN_S1E,SIGMAG_S1E) &
            +Moment0(rL,rU,Num_S1I,DGN_S1I,SIGMAG_S1I) &
            +Moment0(rL,rU,Num_W1E,DGN_W1E,SIGMAG_W1E) &
            +Moment0(rL,rU,Num_W1I,DGN_W1I,SIGMAG_W1I) &
            +Moment0(rL,rU,Num_N1,DGN_N1,SIGMAG_N1) &
 
            +Moment0(rL,rU,Num_EC2,DGN_EC2,SIGMAG_EC2) &
            +Moment0(rL,rU,Num_OC2,DGN_OC2,SIGMAG_OC2) &
            +Moment0(rL,rU,Num_S2E,DGN_S2E,SIGMAG_S2E) &
            +Moment0(rL,rU,Num_S2I,DGN_S2I,SIGMAG_S2I) &
            +Moment0(rL,rU,Num_W2E,DGN_W2E,SIGMAG_W2E) &
            +Moment0(rL,rU,Num_W2I,DGN_W2I,SIGMAG_W2I) &
            +Moment0(rL,rU,Num_N2,DGN_N2,SIGMAG_N2) &
 
            +Moment0(rL,rU,Num_EC3,DGN_EC3,SIGMAG_EC3) &
            +Moment0(rL,rU,Num_OC3,DGN_OC3,SIGMAG_OC3) &
            +Moment0(rL,rU,Num_S3E,DGN_S3E,SIGMAG_S3E) &
            +Moment0(rL,rU,Num_S3I,DGN_S3I,SIGMAG_S3I) &
            +Moment0(rL,rU,Num_W3E,DGN_W3E,SIGMAG_W3E) &
            +Moment0(rL,rU,Num_W3I,DGN_W3I,SIGMAG_W3I) &
            +Moment0(rL,rU,Num_N3,DGN_N3,SIGMAG_N3)
      cAero(1)%cInt(:,:,:,i)=cAero(1)%cInt(:,:,:,i)*VolB/(VolB+Eps)

      cAero(2)%cInt(:,:,:,i)=Moment3(rL,rU,Num_EC1,DGN_EC1,SIGMAG_EC1) &
              +Moment3(rL,rU,Num_EC2,DGN_EC2,SIGMAG_EC2) &
              +Moment3(rL,rU,Num_EC3,DGN_EC3,SIGMAG_EC3) &
   
              +Moment3(rL,rU,Num_OC1,DGN_OC1,SIGMAG_OC1) &
              +Moment3(rL,rU,Num_OC2,DGN_OC2,SIGMAG_OC2) &
              +Moment3(rL,rU,Num_OC3,DGN_OC3,SIGMAG_OC3) &
 
                       +Moment3(rL,rU,Num_S1E,DGN_S1E,SIGMAG_S1E) &
               +Moment3(rL,rU,Num_S2E,DGN_S2E,SIGMAG_S2E) &
               +Moment3(rL,rU,Num_S3E,DGN_S3E,SIGMAG_S3E) &
 
                       +Moment3(rL,rU,Num_S1I,DGN_S1I,SIGMAG_S1I) &
               +Moment3(rL,rU,Num_S2I,DGN_S2I,SIGMAG_S2I) &
               +Moment3(rL,rU,Num_S3I,DGN_S3I,SIGMAG_S3I) &
 
                       +Moment3(rL,rU,Num_W1E,DGN_W1E,SIGMAG_W1E) &
               +Moment3(rL,rU,Num_W2E,DGN_W2E,SIGMAG_W2E) &
               +Moment3(rL,rU,Num_W3E,DGN_W3E,SIGMAG_W3E) &
 
                       +Moment3(rL,rU,Num_W1I,DGN_W1I,SIGMAG_W1I) &
               +Moment3(rL,rU,Num_W2I,DGN_W2I,SIGMAG_W2I) &
               +Moment3(rL,rU,Num_W3I,DGN_W3I,SIGMAG_W3I) &
 
                       +Moment3(rL,rU,Num_N1,DGN_N1,SIGMAG_N1) &
               +Moment3(rL,rU,Num_N2,DGN_N2,SIGMAG_N2) &
               +Moment3(rL,rU,Num_N3,DGN_N3,SIGMAG_N3)
      cAero(2)%cInt(:,:,:,i)=cAero(2)%cInt(:,:,:,i)*VolB/(VolB+Eps)
    END DO
  END DO

END SUBROUTINE InitAerosol

END MODULE InitEx_Mod
