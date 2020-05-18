MODULE Reverse_Mod

  USE Floor_Mod
  USE BoundaryCond_Mod
  USE Control_Mod

  IMPLICIT NONE
  REAL(RealKind)  :: Fac1_L,Fac1_R
  REAL(RealKind)  :: Fac2_L,Fac2_R
  REAL(RealKind)  :: Fac3_L,Fac3_R
  REAL(RealKind)  :: FacW,FacE
  REAL(RealKind)  :: FacS,FacN
  REAL(RealKind)  :: FacB,FacT

! REAL(RealKind), PRIVATE, PARAMETER :: Fac=1.0e-2
! REAL(RealKind), PRIVATE, PARAMETER :: Fac=2.0e00

  INTEGER :: xIOrder,yIOrder,zIOrder
  INTEGER :: Coarse1,Coarse2,Coarse3
  LOGICAL, PRIVATE :: Input=.TRUE. 
  INTEGER, PRIVATE :: Shape(3)
  INTEGER, PRIVATE, PARAMETER :: First=1
  INTEGER, PRIVATE, PARAMETER :: Second=2
  INTEGER, PRIVATE, PARAMETER :: Third=3
  INTEGER :: n1,n2,n3

CONTAINS

SUBROUTINE SetFactor

  IF (Input) THEN
    CALL InputOrder
    INPUT=.FALSE.
  END IF 
  FacW=Two 
  FacE=Two 
  FacS=Two 
  FacN=Two 
  FacB=Two 
  FacT=Two 
  IF (TypeB=='ob') THEN
    FacB=Fac
  END IF
  IF (TypeT=='ot') THEN
    FacT=Fac
  END IF
  IF (TypeW=='ow') THEN
    FacW=Fac
  END IF
  IF (TypeE=='oe') THEN
    FacE=Fac
  END IF
  IF (TypeS=='os') THEN
    FacS=Fac
  END IF
  IF (TypeN=='on') THEN
    FacN=Fac
  END IF
  SELECT CASE(xOrder)
    CASE(First)
      Fac1_L=FacW 
      Fac1_R=FacE 
    CASE(Second)
      Fac2_L=FacW 
      Fac2_R=FacE 
    CASE(Third)
      Fac3_L=FacW 
      Fac3_R=FacE 
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      Fac1_L=FacS
      Fac1_R=FacN
    CASE(Second)
      Fac2_L=FacS
      Fac2_R=FacN
    CASE(Third)
      Fac3_L=FacS
      Fac3_R=FacN
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      Fac1_L=FacB
      Fac1_R=FacT
    CASE(Second)
      Fac2_L=FacB
      Fac2_R=FacT
    CASE(Third)
      Fac3_L=FacB
      Fac3_R=FacT
  END SELECT

END SUBROUTINE SetFactor

SUBROUTINE SetFactorD

  IF (Input) THEN
    CALL InputOrder
    INPUT=.FALSE.
  END IF 
  FacW=Two 
  FacE=Two 
  FacS=Two 
  FacN=Two 
  FacB=Two 
  FacT=Two 
! IF (TypeB=='ob') THEN
!   FacB=Fac
!   IF (BCVec(wPosL)%Bottom/='OutFlow'.AND.BCVec(wPosL)%Bottom/='MeanValue') THEN ! Hinneburg (vorher BCVel...)
!     FacB=Zero
!   END IF
! END IF
! IF (TypeT=='ot') THEN
!   FacT=Fac
!   IF (BCVec(wPosL)%Top/='OutFlow'.AND.BCVec(wPosL)%Top/='MeanValue') THEN
!     FacT=Zero
!   END IF
! END IF
! IF (TypeW=='ow') THEN
!   FacW=Fac
!   IF (BCVec(uPosL)%West/='OutFlow'.AND.BCVec(uPosL)%West/='MeanValue') THEN
!     FacW=Zero
!   END IF
! END IF
! IF (TypeE=='oe') THEN
!   FacE=Fac
!   IF (BCVec(uPosL)%East/='OutFlow'.AND.BCVec(uPosL)%East/='MeanValue') THEN
!     FacE=Zero
!   END IF
! END IF
! IF (TypeS=='os') THEN
!   FacS=Fac
!   IF (BCVec(vPosL)%South/='OutFlow'.AND.BCVec(vPosL)%South/='MeanValue') THEN
!     FacS=Zero
!   END IF
! END IF
! IF (TypeN=='on') THEN
!   FacN=Fac
!   IF (BCVec(vPosL)%North/='OutFlow'.AND.BCVec(vPosL)%North/='MeanValue') THEN
!     FacN=Zero
!   END IF
! END IF
  IF (TypeB=='ob') THEN
    FacB=Fac
    IF (BCVel%Bottom/='OutFlow') THEN
      FacB=Zero
    END IF
  END IF
  IF (TypeT=='ot') THEN
    FacT=Fac
    IF (BCVel%Top/='OutFlow') THEN
      FacT=Zero
    END IF
  END IF
  IF (TypeW=='ow') THEN
    FacW=Fac
    IF (BCVel%West/='OutFlow') THEN
      FacW=Zero
    END IF
  END IF
  IF (TypeE=='oe') THEN
    FacE=Fac
    IF (BCVel%East/='OutFlow') THEN
      FacE=Zero
    END IF
  END IF
  IF (TypeS=='os') THEN
    FacS=Fac
    IF (BCVel%South/='OutFlow') THEN
      FacS=Zero
    END IF
  END IF
  IF (TypeN=='on') THEN
    FacN=Fac
    IF (BCVel%North/='OutFlow') THEN
      FacN=Zero
    END IF
  END IF
  SELECT CASE(xOrder)
    CASE(First)
      Fac1_L=FacW 
      Fac1_R=FacE 
    CASE(Second)
      Fac2_L=FacW 
      Fac2_R=FacE 
    CASE(Third)
      Fac3_L=FacW 
      Fac3_R=FacE 
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      Fac1_L=FacS
      Fac1_R=FacN
    CASE(Second)
      Fac2_L=FacS
      Fac2_R=FacN
    CASE(Third)
      Fac3_L=FacS
      Fac3_R=FacN
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      Fac1_L=FacB
      Fac1_R=FacT
    CASE(Second)
      Fac2_L=FacB
      Fac2_R=FacT
    CASE(Third)
      Fac3_L=FacB
      Fac3_R=FacT
  END SELECT

END SUBROUTINE SetFactorD

SUBROUTINE SetFactorN

  IF (Input) THEN
    CALL InputOrder
    INPUT=.FALSE.
  END IF 
  FacW=Zero 
  FacE=Zero 
  FacS=Zero 
  FacN=Zero 
  FacB=Zero 
  FacT=Zero 
! IF (TypeB=='ob') THEN
!   FacB=Fac
!   IF (BCVel%Bottom/='OutFlow') THEN
!     FacB=Zero
!   END IF
! END IF
! IF (TypeT=='ot') THEN
!   FacT=Fac
!   IF (BCVel%Top/='OutFlow') THEN
!     FacT=Zero
!   END IF
! END IF
! IF (TypeW=='ow') THEN
!   FacW=Fac
!   IF (BCVel%West/='OutFlow') THEN
!     FacW=Zero
!   END IF
! END IF
! IF (TypeE=='oe') THEN
!   FacE=Fac
!   IF (BCVel%East/='OutFlow') THEN
!     FacE=Zero
!   END IF
! END IF
! IF (TypeS=='os') THEN
!   FacS=Fac
!   IF (BCVel%South/='OutFlow') THEN
!     FacS=Zero
!   END IF
! END IF
! IF (TypeN=='on') THEN
!   FacN=Fac
!   IF (BCVel%North/='OutFlow') THEN
!     FacN=Zero
!   END IF
! END IF
  SELECT CASE(xOrder)
    CASE(First)
      Fac1_L=FacW 
      Fac1_R=FacE 
    CASE(Second)
      Fac2_L=FacW 
      Fac2_R=FacE 
    CASE(Third)
      Fac3_L=FacW 
      Fac3_R=FacE 
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      Fac1_L=FacS
      Fac1_R=FacN
    CASE(Second)
      Fac2_L=FacS
      Fac2_R=FacN
    CASE(Third)
      Fac3_L=FacS
      Fac3_R=FacN
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      Fac1_L=FacB
      Fac1_R=FacT
    CASE(Second)
      Fac2_L=FacB
      Fac2_R=FacT
    CASE(Third)
      Fac3_L=FacB
      Fac3_R=FacT
  END SELECT

END SUBROUTINE SetFactorN

SUBROUTINE SetCoarse

  SELECT CASE(xOrder)
    CASE(First)
      Coarse1=xCoarse
    CASE(Second)
      Coarse2=xCoarse
    CASE(Third)
      Coarse3=xCoarse
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      Coarse1=yCoarse
    CASE(Second)
      Coarse2=yCoarse
    CASE(Third)
      Coarse3=yCoarse
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      Coarse1=zCoarse
    CASE(Second)
      Coarse2=zCoarse
    CASE(Third)
      Coarse3=zCoarse
  END SELECT

END SUBROUTINE SetCoarse


SUBROUTINE InputOrder

  INTEGER :: Temp(3)

  Temp(xOrder)=1
  Temp(yOrder)=2
  Temp(zOrder)=3
  xIOrder=Temp(1)
  yIOrder=Temp(2)
  zIOrder=Temp(3)

END SUBROUTINE InputOrder

SUBROUTINE ReverseIndices(nx,ny,nz)


! INTEGER :: n1,n2,n3
  INTEGER :: nx,ny,nz

  IF (Input) THEN
    CALL InputOrder
    INPUT=.FALSE.
  END IF 

  SELECT CASE(xOrder)
    CASE(First)
      n1=nx
    CASE(Second)
      n2=nx
    CASE(Third)
      n3=nx
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      n1=ny
    CASE(Second)
      n2=ny
    CASE(Third)
      n3=ny
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      n1=nz
    CASE(Second)
      n2=nz
    CASE(Third)
      n3=nz
  END SELECT

END SUBROUTINE ReverseIndices

SUBROUTINE ReverseWeight(F1,F2,F3,FU,FV,FW)

  REAL(RealKind) :: F1(:,:,:),F2(:,:,:),F3(:,:,:)
  REAL(RealKind) :: FU(:,:,:),FV(:,:,:),FW(:,:,:)

  INTEGER :: nx,ny,nz
  INTEGER :: ix,iy,iz

  nx=SIZE(FU,1)-1
  ny=SIZE(FU,2)
  nz=SIZE(FU,3)
  SELECT CASE(xOrder)
    CASE(First)
      F1=RESHAPE(FU,(/n1+1,n2,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Second)
      F2=RESHAPE(FU,(/n1,n2+1,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Third)
      F3=RESHAPE(FU,(/n1,n2,n3+1/),ORDER=(/xOrder,yOrder,zOrder/))
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      F1=RESHAPE(FV,(/n1+1,n2,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Second)
      F2=RESHAPE(FV,(/n1,n2+1,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Third)
      F3=RESHAPE(FV,(/n1,n2,n3+1/),ORDER=(/xOrder,yOrder,zOrder/))
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      F1=RESHAPE(FW,(/n1+1,n2,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Second)
      F2=RESHAPE(FW,(/n1,n2+1,n3/),ORDER=(/xOrder,yOrder,zOrder/))
    CASE(Third)
      F3=RESHAPE(FW,(/n1,n2,n3+1/),ORDER=(/xOrder,yOrder,zOrder/))
  END SELECT

END SUBROUTINE ReverseWeight

SUBROUTINE ReverseIncr(d1,d2,d3,dx,dy,dz)

  REAL(RealKind) :: d1(:),d2(:),d3(:)
  REAL(RealKind) :: dx(:),dy(:),dz(:)

  SELECT CASE(xOrder)
    CASE(First)
      d1=dx
    CASE(Second)
      d2=dx
    CASE(Third)
      d3=dx
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      d1=dy
    CASE(Second)
      d2=dy
    CASE(Third)
      d3=dy
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      d1=dz
    CASE(Second)
      d2=dz
    CASE(Third)
      d3=dz
  END SELECT

END SUBROUTINE ReverseIncr

SUBROUTINE ReverseMetr(Metr12,Metr13  &
                      ,Metr21,Metr23  & 
                      ,Metr31,Metr32  &
                      ,MetrXY,MetrXZ  &
                      ,MetrYX,MetrYZ  &
                      ,MetrZX,MetrZY)

  REAL(RealKind) :: Metr12(:),Metr13(:)
  REAL(RealKind) :: Metr21(:),Metr23(:)
  REAL(RealKind) :: Metr31(:),Metr32(:)
  REAL(RealKind) :: MetrXY(:),MetrXZ(:)
  REAL(RealKind) :: MetrYX(:),MetrYZ(:)
  REAL(RealKind) :: MetrZX(:),MetrZY(:)

  SELECT CASE(xOrder)
    CASE(First)
      IF (yOrder==2) THEN
        Metr12=MetrXY
        Metr13=MetrXZ
      ELSE   
        Metr12=MetrXZ
        Metr13=MetrXY
      END IF  
    CASE(Second)
      IF (yOrder==1) THEN
        Metr21=MetrXY
        Metr23=MetrXZ
      ELSE   
        Metr21=MetrXZ
        Metr23=MetrXY
      END IF  
    CASE(Third)
      IF (yOrder==1) THEN
        Metr31=MetrXY
        Metr32=MetrXZ
      ELSE   
        Metr31=MetrXZ
        Metr32=MetrXY
      END IF  
  END SELECT
  SELECT CASE(yOrder)
    CASE(First)
      IF (xOrder==2) THEN
        Metr12=MetrYX
        Metr13=MetrYZ
      ELSE   
        Metr12=MetrYZ
        Metr13=MetrYX
      END IF  
    CASE(Second)
      IF (xOrder==1) THEN
        Metr21=MetrYX
        Metr23=MetrYZ
      ELSE   
        Metr21=MetrYZ
        Metr23=MetrYX
      END IF  
    CASE(Third)
      IF (xOrder==1) THEN
        Metr31=MetrYX
        Metr32=MetrYZ
      ELSE   
        Metr31=MetrYZ
        Metr32=MetrYX
      END IF  
  END SELECT
  SELECT CASE(zOrder)
    CASE(First)
      IF (xOrder==2) THEN
        Metr12=MetrZX
        Metr13=MetrZY
      ELSE   
        Metr12=MetrZY
        Metr13=MetrZX
      END IF  
    CASE(Second)
      IF (xOrder==1) THEN
        Metr21=MetrZX
        Metr23=MetrZY
      ELSE   
        Metr21=MetrZY
        Metr23=MetrZX
      END IF  
    CASE(Third)
      IF (xOrder==1) THEN
        Metr31=MetrZX
        Metr32=MetrZY
      ELSE   
        Metr31=MetrZY
        Metr32=MetrZX
      END IF  
  END SELECT

END SUBROUTINE ReverseMetr

SUBROUTINE Reverse(pRev,p)

  REAL(RealKind) :: pRev(:,:,:)
  REAL(RealKind) :: p(:,:,:)

  INTEGER :: nx,ny,nz
  INTEGER :: ix,iy,iz

  nx=SIZE(p,1)
  ny=SIZE(p,2)
  nz=SIZE(p,3)
  pRev=RESHAPE(p,(/n1,n2,n3/),ORDER=(/xOrder,yOrder,zOrder/))

END SUBROUTINE Reverse
SUBROUTINE ReverseBack(p,pRev)

  REAL(RealKind) :: p(:,:,:)
  REAL(RealKind) :: pRev(:,:,:)

  INTEGER :: nx,ny,nz
  INTEGER :: ix,iy,iz

  nx=SIZE(p,1)
  ny=SIZE(p,2)
  nz=SIZE(p,3)
  p=RESHAPE(pRev,(/nx,ny,nz/),ORDER=(/xIOrder,yIOrder,zIOrder/))

END SUBROUTINE ReverseBack

SUBROUTINE ReverseB(pRev,p)

  REAL(RealKind) :: pRev(:,:,:,:)
  REAL(RealKind) :: p(:,:,:,:)

  INTEGER :: nx,ny,nz
  INTEGER :: ix,iy,iz

  nx=SIZE(p,2)
  ny=SIZE(p,3)
  nz=SIZE(p,4)
  pRev=RESHAPE(p,(/2,n1,n2,n3/),ORDER=(/1,xOrder+1,yOrder+1,zOrder+1/))

END SUBROUTINE ReverseB

SUBROUTINE ReverseBBack(p,pRev)

  REAL(RealKind) :: p(:,:,:,:)
  REAL(RealKind) :: pRev(:,:,:,:)

  INTEGER :: nx,ny,nz
  INTEGER :: ix,iy,iz

  nx=SIZE(p,2)
  ny=SIZE(p,3)
  nz=SIZE(p,4)
  p=RESHAPE(pRev,(/2,nx,ny,nz/),ORDER=(/1,xIOrder+1,yIOrder+1,zIOrder+1/))

END SUBROUTINE ReverseBBack


END MODULE Reverse_Mod

