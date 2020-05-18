MODULE RungeKuttaMethods_Mod
  USE Parameter_Mod

  IMPLICIT NONE

  TYPE Koeff_T
    REAL(RealKind), POINTER :: Koeff(:)
  END TYPE Koeff_T
  TYPE RKMethod_T
    INTEGER :: nStage
    INTEGER :: nPhi
    REAL(8) :: CFLNumber=0.0d0
    TYPE(Koeff_T), POINTER :: aRunge(:,:)
    REAL(RealKind), POINTER :: cRunge(:)
    REAL(RealKind), POINTER :: dRunge(:,:)
    REAL(RealKind), POINTER :: gRunge(:,:)
    REAL(RealKind), POINTER :: dtRunge(:)
  END TYPE RKMethod_T

  TYPE(RKMethod_T) :: RK1
  TYPE(RKMethod_T) :: RK2a
  TYPE(RKMethod_T) :: RK2Ros
  TYPE(RKMethod_T) :: RK2b
  TYPE(RKMethod_T) :: RK2bT
  TYPE(RKMethod_T) :: RK3
  TYPE(RKMethod_T) :: RK3TVDF
  TYPE(RKMethod_T) :: RK4
  TYPE(RKMethod_T) :: RK43K
  TYPE(RKMethod_T) :: RKEB4
  TYPE(RKMethod_T) :: RKWen
  TYPE(RKMethod_T) :: RKJeb
  TYPE(RKMethod_T) :: RKLan
  TYPE(RKMethod_T) :: RK2Nlopt
  TYPE(RKMethod_T) :: RKNlopt
  TYPE(RKMethod_T) :: RKtvdA
  TYPE(RKMethod_T) :: RKtvdB
  TYPE(RKMethod_T) :: RKJeb1
  TYPE(RKMethod_T) :: RKW23
  TYPE(RKMethod_T) :: RKCel
  TYPE(RKMethod_T) :: RKIMEX
  TYPE(RKMethod_T) :: RKN4EX
  TYPE(RKMethod_T) :: MIS4
  TYPE(RKMethod_T) :: MIS4_4
  TYPE(RKMethod_T) :: MIS2
  TYPE(RKMethod_T) :: MIS3C
  TYPE(RKMethod_T) :: RKN4E2
  TYPE(RKMethod_T) :: RKN3E4


  TYPE(RKMethod_T) :: MIS_SV2_2_3_10_A_fixed
  TYPE(RKMethod_T) :: MIS_SV2_2_4_14_A_fixed
  TYPE(RKMethod_T) :: MIS_SV2_2_5_19_A_fixed
  TYPE(RKMethod_T) :: MIS_SV3_3_4_13_A_fixed
  TYPE(RKMethod_T) :: MIS_SV3_3_5_17_A_fixed

CONTAINS

SUBROUTINE MethodsRK

  INTEGER :: i,j
  INTEGER :: nPhi,nStage
  TYPE(Koeff_T), POINTER :: aRunge(:,:)
  REAL(RealKind), POINTER :: cRunge(:)
  REAL(RealKind), POINTER :: dRunge(:,:)
  REAL(RealKind), POINTER :: gRunge(:,:)
  REAL(RealKind), POINTER :: dtRunge(:)
  REAL(RealKind) :: alpha
  REAL(RealKind) :: ga,de

  nStage=4
  nPhi=2
  RKEB4%nStage=nStage
  RKEB4%nPhi=nPhi
  ALLOCATE(RKEB4%aRunge(nStage+1,nStage))
  ALLOCATE(RKEB4%cRunge(nStage+1))
  ALLOCATE(RKEB4%dRunge(nStage+1,nStage))
  ALLOCATE(RKEB4%gRunge(nStage+1,nStage))
  ALLOCATE(RKEB4%dtRunge(nStage+1))
  aRunge=>RKEB4%aRunge
  cRunge=>RKEB4%cRunge
  dRunge=>RKEB4%dRunge
  gRunge=>RKEB4%gRunge
  dtRunge=>RKEB4%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=0.5d0

  aRunge(3,1)%Koeff(0)=0.5d0
  aRunge(3,1)%Koeff(1)=-1.0d0
  aRunge(3,2)%Koeff(1)=1.0d0

  aRunge(4,1)%Koeff(0)=1.0d0
  aRunge(4,1)%Koeff(1)=-2.0d0
  aRunge(4,3)%Koeff(1)=2.0d0

  aRunge(5,1)%Koeff(0)=1.0d0
  aRunge(5,1)%Koeff(1)=-3.0d0
  aRunge(5,1)%Koeff(2)=2.0d0

  aRunge(5,2)%Koeff(1)=2.0d0
  aRunge(5,2)%Koeff(2)=-2.0d0

  aRunge(5,3)%Koeff(1)=2.0d0
  aRunge(5,3)%Koeff(2)=-2.0d0

  aRunge(5,4)%Koeff(1)=-1.0d0
  aRunge(5,4)%Koeff(2)=2.0d0

  cRunge(2)=0.5d0
  cRunge(3)=0.5d0
  cRunge(4)=1.0d0
  cRunge(5)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=2
  nPhi=0
  RK2Ros%nStage=nStage
  RK2Ros%nPhi=nPhi
  ALLOCATE(RK2Ros%aRunge(nStage+1,nStage))
  ALLOCATE(RK2Ros%cRunge(nStage+1))
  ALLOCATE(RK2Ros%dRunge(nStage+1,nStage))
  ALLOCATE(RK2Ros%gRunge(nStage+1,nStage))
  ALLOCATE(RK2Ros%dtRunge(nStage+1))
  aRunge=>RK2Ros%aRunge
  cRunge=>RK2Ros%cRunge
  dRunge=>RK2Ros%dRunge
  gRunge=>RK2Ros%gRunge
  dtRunge=>RK2Ros%dtRunge
  dRunge=Zero
  gRunge=Zero
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=2.0d0/3.0d0
  aRunge(3,1)%Koeff(0)=1.0d0/4.0d0
  aRunge(3,2)%Koeff(0)=3.0d0/4.0d0
  cRunge(2)=2.0d0/3.0d0
  cRunge(3)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=1
  nPhi=0
  RK1%nStage=nStage
  RK1%nPhi=nPhi
  ALLOCATE(RK1%aRunge(nStage+1,nStage))
  ALLOCATE(RK1%cRunge(nStage+1))
  ALLOCATE(RK1%dRunge(nStage+1,nStage))
  ALLOCATE(RK1%gRunge(nStage+1,nStage))
  ALLOCATE(RK1%dtRunge(nStage+1))
  aRunge=>RK1%aRunge
  cRunge=>RK1%cRunge
  dRunge=>RK1%dRunge
  gRunge=>RK1%gRunge
  dtRunge=>RK1%dtRunge
  dRunge=Zero
  gRunge=Zero
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0
  cRunge(2)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=2
  nPhi=0
  RK2a%nStage=nStage
  RK2a%nPhi=nPhi
  ALLOCATE(RK2a%aRunge(nStage+1,nStage))
  ALLOCATE(RK2a%cRunge(nStage+1))
  ALLOCATE(RK2a%dRunge(nStage+1,nStage))
  ALLOCATE(RK2a%gRunge(nStage+1,nStage))
  ALLOCATE(RK2a%dtRunge(nStage+1))
  aRunge=>RK2a%aRunge
  cRunge=>RK2a%cRunge
  dRunge=>RK2a%dRunge
  gRunge=>RK2a%gRunge
  dtRunge=>RK2a%dtRunge
  dRunge=Zero
  gRunge=Zero
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=0.5d0
  aRunge(3,2)%Koeff(0)=1.0d0
  cRunge(2)=0.5d0
  cRunge(3)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=2
  nPhi=0
  RK2b%nStage=nStage
  RK2b%nPhi=nPhi
  ALLOCATE(RK2b%aRunge(nStage+1,nStage))
  ALLOCATE(RK2b%cRunge(nStage+1))
  ALLOCATE(RK2b%dRunge(nStage+1,nStage))
  ALLOCATE(RK2b%gRunge(nStage+1,nStage))
  ALLOCATE(RK2b%dtRunge(nStage+1))
  aRunge=>RK2b%aRunge
  cRunge=>RK2b%cRunge
  dRunge=>RK2b%dRunge
  gRunge=>RK2b%gRunge
  dtRunge=>RK2b%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0
  aRunge(3,1)%Koeff(0)=0.5d0
  aRunge(3,2)%Koeff(0)=0.5d0
  cRunge(2)=1.0d0
  cRunge(3)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=2
  nPhi=1
  RK2bT%nStage=nStage
  RK2bT%nPhi=nPhi
  ALLOCATE(RK2bT%aRunge(nStage+1,nStage))
  ALLOCATE(RK2bT%cRunge(nStage+1))
  ALLOCATE(RK2bT%dRunge(nStage+1,nStage))
  ALLOCATE(RK2bT%gRunge(nStage+1,nStage))
  ALLOCATE(RK2bT%dtRunge(nStage+1))
  aRunge=>RK2bT%aRunge
  cRunge=>RK2bT%cRunge
  dRunge=>RK2bT%dRunge
  gRunge=>RK2bT%gRunge
  dtRunge=>RK2bT%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0
  aRunge(3,1)%Koeff(0)=1.0d0
  aRunge(3,1)%Koeff(1)=-1.0d0
  aRunge(3,2)%Koeff(1)=1.0d0
  cRunge(2)=1.0d0
  cRunge(3)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO
  
  nStage=3
  nPhi=0
  RK3%nStage=nStage
  RK3%nPhi=nPhi
  RK3%CFLNumber=1.5d0
  ALLOCATE(RK3%aRunge(nStage+1,nStage))
  ALLOCATE(RK3%cRunge(nStage+1))
  ALLOCATE(RK3%dRunge(nStage+1,nStage))
  ALLOCATE(RK3%gRunge(nStage+1,nStage))
  ALLOCATE(RK3%dtRunge(nStage+1))
  aRunge=>RK3%aRunge
  cRunge=>RK3%cRunge
  dRunge=>RK3%dRunge
  gRunge=>RK3%gRunge
  dtRunge=>RK3%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0/3.0d0
  aRunge(3,2)%Koeff(0)=0.5d0
  aRunge(4,3)%Koeff(0)=1.0d0
  cRunge(2)=0.0d0
  cRunge(2)=1.0d0/3.0d0
  cRunge(3)=0.5d0
  cRunge(4)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO
  dRunge=Zero
  gRunge=Zero

  nStage=3
  nPhi=0
  RK3TVDF%nStage=nStage
  RK3TVDF%nPhi=nPhi
  ALLOCATE(RK3TVDF%aRunge(nStage+1,nStage))
  ALLOCATE(RK3TVDF%cRunge(nStage+1))
  ALLOCATE(RK3TVDF%dRunge(nStage+1,nStage))
  ALLOCATE(RK3TVDF%gRunge(nStage+1,nStage))
  ALLOCATE(RK3TVDF%dtRunge(nStage+1))
  aRunge=>RK3TVDF%aRunge
  cRunge=>RK3TVDF%cRunge
  dRunge=>RK3TVDF%dRunge
  gRunge=>RK3TVDF%gRunge
  dtRunge=>RK3TVDF%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0
  aRunge(3,2)%Koeff(0)=1.0d0/4.0d0
  aRunge(4,3)%Koeff(0)=2.0d0/3.0d0
  cRunge(2)=1.0d0
  cRunge(3)=1.0d0/4.0d0
  cRunge(4)=2.0d0/3.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO
  dRunge=Zero
  dRunge(3,2)=1.0d0/4.0d0
  dRunge(4,3)=2.0d0/3.0d0
  gRunge=Zero

  nStage=3
  alpha=2.0d0/3.0d0
  nPhi=0
  RKCel%nStage=nStage
  RKCel%nPhi=nPhi
  ALLOCATE(RKCel%aRunge(nStage+1,nStage))
  ALLOCATE(RKCel%cRunge(nStage+1))
  ALLOCATE(RKCel%dRunge(nStage+1,nStage))
  ALLOCATE(RKCel%gRunge(nStage+1,nStage))
  ALLOCATE(RKCel%dtRunge(nStage+1))
  aRunge=>RKCel%aRunge
  cRunge=>RKCel%cRunge
  dRunge=>RKCel%dRunge
  gRunge=>RKCel%gRunge
  dtRunge=>RKCel%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=1.0d0/3.0d0
  aRunge(3,1)%Koeff(0)=alpha*(2.0d0-3.0d0*alpha)
  aRunge(3,2)%Koeff(0)=alpha*(3.0d0*alpha-1.0d0)
  aRunge(4,1)%Koeff(0)=(3.0d0-5.0d0*alpha)/(6.0d0*alpha)
  aRunge(4,2)%Koeff(0)=3.0d0*(3.0d0*alpha-2.0d0)/(2.0d0*(3.0d0*alpha-1.0d0))
  aRunge(4,3)%Koeff(0)=1.0d0/(2.0d0*alpha*(3.0d0*alpha-1.0d0))
  cRunge(2)=1.0d0/3.0d0
  cRunge(3)=alpha
  cRunge(4)=2.0d0/3.0d0
  dRunge=Zero
  dRunge(4,2)=1.0d0
  gRunge=Zero
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=2
  ga=(2.0d0-SQRT(2.0d0))/2.0d0
  de=1.0d0-1.0d0/(2.0d0*ga)
  nPhi=0
  RKIMEX%nStage=nStage
  RKIMEX%nPhi=nPhi
  ALLOCATE(RKIMEX%aRunge(nStage+1,nStage))
  ALLOCATE(RKIMEX%cRunge(nStage+1))
  ALLOCATE(RKIMEX%dRunge(nStage+1,nStage))
  ALLOCATE(RKIMEX%gRunge(nStage+1,nStage))
  ALLOCATE(RKIMEX%dtRunge(nStage+1))
  aRunge=>RKIMEX%aRunge
  cRunge=>RKIMEX%cRunge
  dRunge=>RKIMEX%dRunge
  gRunge=>RKIMEX%gRunge
  dtRunge=>RKIMEX%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=ga
  aRunge(3,1)%Koeff(0)=-(1.0d0-ga)+de
  aRunge(3,2)%Koeff(0)=(1.0d0-de)
  cRunge=Zero
  cRunge(2)=ga
  cRunge(3)=ga
  dRunge=Zero
  dRunge(3,2)=(1.0d0-ga)/ga
  gRunge=Zero
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=3
  nPhi=0
  RKWen%nStage=nStage
  RKWen%nPhi=nPhi
  ALLOCATE(RKWen%aRunge(nStage+1,nStage))
  ALLOCATE(RKWen%cRunge(nStage+1))
  ALLOCATE(RKWen%dRunge(nStage+1,nStage))
  ALLOCATE(RKWen%gRunge(nStage+1,nStage))
  ALLOCATE(RKWen%dtRunge(nStage+1))
  aRunge=>RKWen%aRunge
  cRunge=>RKWen%cRunge
  dRunge=>RKWen%dRunge
  gRunge=>RKWen%gRunge
  dtRunge=>RKWen%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)= 0.4197138394323076d0
  aRunge(3,1)%Koeff(0)= 0.0494103174250211d0
  aRunge(3,2)%Koeff(0)= 0.5761341869014835d0
  aRunge(4,1)%Koeff(0)=-0.2304896605104612d0
  aRunge(4,2)%Koeff(0)= 0.4282028054551005d0
  aRunge(4,3)%Koeff(0)= 0.6892421263984294d0
  cRunge=Zero
  cRunge(2)=aRunge(2,1)%Koeff(0)
  cRunge(3)=aRunge(3,1)%Koeff(0)+aRunge(3,2)%Koeff(0)
  cRunge(4)=aRunge(4,1)%Koeff(0)+aRunge(4,2)%Koeff(0)+aRunge(4,3)%Koeff(0)
  dRunge=Zero
  dRunge(3,2)=0.0567447894712588d0
  dRunge(4,2)=0.5984205948350817d0
  dRunge(4,3)=-0.0761105269067489d0
  gRunge=Zero
  gRunge(3,2)=0.0784430572133444d0
  gRunge(4,2)=0.6791817901685916d0
  gRunge(4,3)=-0.5441329940536196d0

  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=3
  nPhi=0
  RKJeb%nStage=nStage
  RKJeb%nPhi=nPhi
  RKJeb%CFLNumber=1.5d0
  ALLOCATE(RKJeb%aRunge(nStage+1,nStage))
  ALLOCATE(RKJeb%cRunge(nStage+1))
  ALLOCATE(RKJeb%dRunge(nStage+1,nStage))
  ALLOCATE(RKJeb%gRunge(nStage+1,nStage))
  ALLOCATE(RKJeb%dtRunge(nStage+1))
  aRunge=>RKJeb%aRunge
  cRunge=>RKJeb%cRunge
  dRunge=>RKJeb%dRunge
  gRunge=>RKJeb%gRunge
  dtRunge=>RKJeb%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=2.0492941060709863d-001
  aRunge(3,1)%Koeff(0)=-4.5477553356788974d-001
  aRunge(3,2)%Koeff(0)=9.5613538239378981d-001
  aRunge(4,1)%Koeff(0)=-3.5970281266252929d-002
  aRunge(4,2)%Koeff(0)=-1.5363649484946584d-001
  aRunge(4,3)%Koeff(0)=7.0259062712330234d-001


  cRunge=Zero
  cRunge(2)=aRunge(2,1)%Koeff(0)
  cRunge(3)=aRunge(3,1)%Koeff(0)+aRunge(3,2)%Koeff(0)
  cRunge(4)=aRunge(4,1)%Koeff(0)+aRunge(4,2)%Koeff(0)+aRunge(4,3)%Koeff(0)
  dRunge=Zero
  dRunge(3,2)=7.0302371060435331d-001
  dRunge(4,2)=4.2492220536139252d-001
  dRunge(4,3)=5.4545718243573982d-001
  gRunge=Zero
  gRunge(3,2)=-8.2176071248067006d-001
  gRunge(4,2)=-3.8080670922635063d-001
  gRunge(4,3)=4.5653105107801978d-001

  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=3
  nPhi=0
  RKLan%nStage=nStage
  RKLan%nPhi=nPhi
  ALLOCATE(RKLan%aRunge(nStage+1,nStage))
  ALLOCATE(RKLan%cRunge(nStage+1))
  ALLOCATE(RKLan%dRunge(nStage+1,nStage))
  ALLOCATE(RKLan%gRunge(nStage+1,nStage))
  ALLOCATE(RKLan%dtRunge(nStage+1))
  aRunge=>RKLan%aRunge
  cRunge=>RKLan%cRunge
  dRunge=>RKLan%dRunge
  gRunge=>RKLan%gRunge
  dtRunge=>RKLan%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO

  aRunge(2,1)%Koeff(0)=0.18502450114688243d0
  aRunge(3,1)%Koeff(0)=-0.4479371838488236d0
  aRunge(3,2)%Koeff(0)=0.95552546944728722d0
  aRunge(4,1)%Koeff(0)=-3.07363364865417509d-002
  aRunge(4,2)%Koeff(0)=-0.16954108935313769d0
  aRunge(4,3)%Koeff(0)=0.74204173586290767d0
  cRunge=Zero
  cRunge(2)=aRunge(2,1)%Koeff(0)
  cRunge(3)=aRunge(3,1)%Koeff(0)+aRunge(3,2)%Koeff(0)
  cRunge(4)=aRunge(4,1)%Koeff(0)+aRunge(4,2)%Koeff(0)+aRunge(4,3)%Koeff(0)
  dRunge=Zero
  dRunge(3,2)=0.68540723734731002d0
  dRunge(4,2)=0.37487297544895642d0
  dRunge(4,3)=0.50043909986915092d0
  gRunge=Zero
  gRunge(3,2)=-0.76613049141773570d0
  gRunge(4,2)=-0.43276023666330771d0
  gRunge(4,3)=0.42973260825016535d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  include 'RKNlopt.fort'
  include 'RKtvdA.fortA'
  include 'RKtvdB.fortA'
  include 'RK2Nlopt.fort'
  include 'RKN4EX.fortA'
  include 'MIS4.fortA'
  include 'MIS4_4.fortA'
  include 'MIS2.fortA'
  include 'MIS3C.fortA'
  include 'RKN4E2.fortA'
  include 'RKN3E4.fortA'
  include 'MIS_SV3_3_5_17_A_fixed.fortA'
  include 'MIS_SV3_3_4_13_A_fixed.fortA'
  include 'MIS_SV2_2_3_10_A_fixed.fortA' 
  include 'MIS_SV2_2_4_14_A_fixed.fortA' 
  include 'MIS_SV2_2_5_19_A_fixed.fortA' 
  include 'MIS_SV3_3_5_17_A_fixed.fortA'
  nStage=3
  nPhi=0
  RKJeb1%nStage=nStage
  RKJeb1%nPhi=nPhi
  RKJeb1%CFLNumber=1.5d0
  ALLOCATE(RKJeb1%aRunge(nStage+1,nStage))
  ALLOCATE(RKJeb1%cRunge(nStage+1))
  ALLOCATE(RKJeb1%dRunge(nStage+1,nStage))
  ALLOCATE(RKJeb1%gRunge(nStage+1,nStage))
  ALLOCATE(RKJeb1%dtRunge(nStage+1))
  aRunge=>RKJeb1%aRunge
  cRunge=>RKJeb1%cRunge
  dRunge=>RKJeb1%dRunge
  gRunge=>RKJeb1%gRunge
  dtRunge=>RKJeb1%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO

  aRunge(2,1)%Koeff(0)=0.3391d0
  aRunge(3,1)%Koeff(0)=-0.6126d0
  aRunge(3,2)%Koeff(0)=1.4555d0
  aRunge(4,1)%Koeff(0)=0.4687d0
  aRunge(4,2)%Koeff(0)=0.0642d0
  aRunge(4,3)%Koeff(0)=0.3376d0

  cRunge=Zero
  cRunge(2)=aRunge(2,1)%Koeff(0)
  cRunge(3)=aRunge(3,1)%Koeff(0)+aRunge(3,2)%Koeff(0)
  cRunge(4)=aRunge(4,1)%Koeff(0)+aRunge(4,2)%Koeff(0)+aRunge(4,3)%Koeff(0)
  dRunge=Zero
  dRunge(3,2)=-1.7934d0
  dRunge(4,2)=0.2854d0
  dRunge(4,3)=0.0258d0

  gRunge=Zero
  gRunge(3,2)=1.9215d0
  gRunge(4,2)=-0.8512d0
  gRunge(4,3)=0.3367d0

  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=3
  nPhi=0
  RKW23%nStage=nStage
  RKW23%nPhi=nPhi
  ALLOCATE(RKW23%aRunge(nStage+1,nStage))
  ALLOCATE(RKW23%cRunge(nStage+1))
  ALLOCATE(RKW23%dRunge(nStage+1,nStage))
  ALLOCATE(RKW23%gRunge(nStage+1,nStage))
  ALLOCATE(RKW23%dtRunge(nStage+1))
  aRunge=>RKW23%aRunge
  cRunge=>RKW23%cRunge
  dRunge=>RKW23%dRunge
  gRunge=>RKW23%gRunge
  dtRunge=>RKW23%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO


  aRunge(2,1)%Koeff(0)=2.4094098957512181e-01
  aRunge(3,1)%Koeff(0)=-2.5100027360785432e-01
  aRunge(3,2)%Koeff(0)=1.0013565263154491e+00
  aRunge(4,1)%Koeff(0)=1.8466850240228652e-01
  aRunge(4,2)%Koeff(0)=-1.1401433099445987e-01
  aRunge(4,3)%Koeff(0)=6.9079491510174174e-01

  cRunge=Zero
  cRunge(2)=aRunge(2,1)%Koeff(0)
  cRunge(3)=aRunge(3,1)%Koeff(0)+aRunge(3,2)%Koeff(0)
  cRunge(4)=aRunge(4,1)%Koeff(0)+aRunge(4,2)%Koeff(0)+aRunge(4,3)%Koeff(0)
  dRunge=Zero
  dRunge(3,2)=5.0022548472245065e-01
  dRunge(4,2)=6.3430891601597705e-01
  dRunge(4,3)=-5.1355234244538986e-02

  gRunge=Zero
  gRunge(3,2)=-7.5696295789587087e-01
  gRunge(4,2)=-2.5845191935635359e-01
  gRunge(4,3)=2.6630355867838612e-01

  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=4
  nPhi=0
  RK4%nStage=nStage
  RK4%nPhi=nPhi
  ALLOCATE(RK4%aRunge(nStage+1,nStage))
  ALLOCATE(RK4%cRunge(nStage+1))
  ALLOCATE(RK4%dRunge(nStage+1,nStage))
  ALLOCATE(RK4%gRunge(nStage+1,nStage))
  ALLOCATE(RK4%dtRunge(nStage+1))
  aRunge=>RK4%aRunge
  cRunge=>RK4%cRunge
  dRunge=>RK4%dRunge
  gRunge=>RK4%gRunge
  dtRunge=>RK4%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  cRunge=Zero
  dRunge=Zero
  gRunge=Zero
  aRunge(2,1)%Koeff(0)=0.5d0
  aRunge(3,2)%Koeff(0)=0.5d0
  aRunge(4,3)%Koeff(0)=1.0d0
  aRunge(5,1)%Koeff(0)=1.0d0/6.0d0
  aRunge(5,2)%Koeff(0)=1.0d0/3.0d0
  aRunge(5,3)%Koeff(0)=1.0d0/3.0d0
  aRunge(5,4)%Koeff(0)=1.0d0/6.0d0
  cRunge(2)=0.5d0
  cRunge(3)=0.5d0
  cRunge(4)=1.0d0
  cRunge(5)=1.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    DO j=1,i-1
      aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      gRunge(i,j)=gRunge(i,j)/dtRunge(i)
    END DO
  END DO

  nStage=4
  nPhi=0
  RK43K%nStage=nStage
  RK43K%nPhi=nPhi
  ALLOCATE(RK43K%aRunge(nStage+1,nStage))
  ALLOCATE(RK43K%cRunge(nStage+1))
  ALLOCATE(RK43K%dRunge(nStage+1,nStage))
  ALLOCATE(RK43K%gRunge(nStage+1,nStage))
  ALLOCATE(RK43K%dtRunge(nStage+1))
  aRunge=>RK43K%aRunge
  cRunge=>RK43K%cRunge
  dRunge=>RK43K%dRunge
  gRunge=>RK43K%gRunge
  dtRunge=>RK43K%dtRunge
  DO i=2,nStage+1
    DO j=1,i-1
      ALLOCATE(aRunge(i,j)%Koeff(0:nPhi))
      aRunge(i,j)%Koeff=Zero
    END DO
  END DO
  aRunge(2,1)%Koeff(0)=0.5d0
  aRunge(3,1)%Koeff(0)=-2.0d0/3.0d0
  aRunge(3,2)%Koeff(0)=2.0d0/3.0d0
  aRunge(4,1)%Koeff(0)=0.5d0
  aRunge(4,2)%Koeff(0)=-1.0d0
  aRunge(4,3)%Koeff(0)=1.0d0
  aRunge(5,1)%Koeff(0)=-1.0d0/6.0d0
  aRunge(5,2)%Koeff(0)=2.0d0/3.0d0
  aRunge(5,3)%Koeff(0)=-2.0d0/3.0d0
  aRunge(5,4)%Koeff(0)=1.0d0/6.0d0
  cRunge(2)=0.5d0
  cRunge(3)=0.0d0
  cRunge(4)=0.5d0
  cRunge(5)=0.0d0
  DO i=2,nStage+1
    dtRunge(i)=0.0d0
    DO j=1,i-1
      dtRunge(i)=dtRunge(i)+aRunge(i,j)%Koeff(0)
    END DO  
    IF (dtRunge(i)/=Zero) THEN
      DO j=1,i-1
        aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      END DO
    ELSE
      dtRunge(i)=1.d-5
      DO j=1,i-1
        aRunge(i,j)%Koeff=aRunge(i,j)%Koeff/dtRunge(i)
      END DO
    END IF
  END DO
  dRunge=Zero
  dRunge(2,1)=One
  dRunge(3,2)=One
  dRunge(4,3)=One
  dRunge(5,4)=One
  gRunge=Zero
  
END SUBROUTINE MethodsRK

END MODULE RungeKuttaMethods_Mod

