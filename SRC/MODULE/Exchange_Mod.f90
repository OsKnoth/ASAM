MODULE Exchange_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Parallel_Mod

  IMPLICIT NONE

  INTEGER, PARAMETER :: FineFine    =6
  INTEGER, PARAMETER :: FineEqual   =5
  INTEGER, PARAMETER :: EqualFine   =4
  INTEGER, PARAMETER :: EqualEqual  =3
  INTEGER, PARAMETER :: EqualCoarse =2
  INTEGER, PARAMETER :: CoarseEqual =1
  INTEGER, PARAMETER :: CoarseCoarse=0

CONTAINS 

SUBROUTINE CopyFace(u1,FU1,u2,FU2) 

  REAL(RealKind) :: u1(:,:),FU1(:,:),u2(:,:),FU2(:,:)

  INTEGER :: i,j
  INTEGER :: iU1,jU1,iU2,jU2
  REAL(RealKind) :: t1,t2,t3,t4,t

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)
  iU2=UBOUND(u2,1)
  jU2=UBOUND(u2,2)

  SELECT CASE(2*(iU1/iU2)+jU1/jU2)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          t=u1(i-1,j-1)*FU1(i-1,j-1) &
           +u1(i  ,j-1)*FU1(i  ,j-1) &
           +u1(i-1,j  )*FU1(i-1,j  ) &
           +u1(i  ,j  )*FU1(i  ,j  )  
          u1(i  ,j  )=u1(i  ,j  )+u2(i/2,j/2)
          u1(i-1,j  )=u1(i-1,j  )+u2(i/2,j/2)
          u1(i  ,j-1)=u1(i  ,j-1)+u2(i/2,j/2)
          u1(i-1,j-1)=u1(i-1,j-1)+u2(i/2,j/2)
          u2(i/2,j/2)=u2(i/2,j/2)+t/(FU2(i/2,j/2)+Eps)
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          t=u1(i-1,j)*FU1(i-1,j) &
           +u1(i  ,j)*FU1(i  ,j)  
          u1(i  ,j)=u1(i  ,j)+u2(i/2,j)
          u1(i-1,j)=u1(i-1,j)+u2(i/2,j)
          u2(i/2,j)=u2(i/2,j)+t/(FU2(i/2,j)+Eps)
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          t=u1(i,j-1)*FU1(i,j-1) &
           +u1(i,j  )*FU1(i,j  )  
          u1(i,j  )=u1(i,j  )+u2(i,j/2)
          u1(i,j-1)=u1(i,j-1)+u2(i,j/2)
          u2(i,j/2)=u2(i,j/2)+t/(FU2(i,j/2)+Eps)
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          t=u1(i,j)
          u1(i,j)=u1(i,j)+u2(i,j)
          u2(i,j)=u2(i,j)+t
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          t=u2(i,2*j-1)*FU2(i,2*j-1) &
           +u2(i,2*j  )*FU2(i,2*j  )  
          u2(i,2*j  )=u2(i,2*j  )+u1(i,j)
          u2(i,2*j-1)=u2(i,2*j-1)+u1(i,j)
          u1(i,j)=u1(i,j)+t/(FU1(i,j)+Eps)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          t=u2(2*i-1,j  )*FU2(2*i-1,j  ) &
           +u2(2*i  ,j  )*FU2(2*i  ,j  )  
          u2(2*i  ,j)=u2(2*i  ,j)+u1(i,j)
          u2(2*i-1,j)=u2(2*i-1,j)+u1(i,j)
          u1(i,j)=u1(i,j)+t/(FU1(i,j)+Eps)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          t=u2(2*i-1,2*j-1)*FU2(2*i-1,2*j-1) &
           +u2(2*i  ,2*j-1)*FU2(2*i  ,2*j-1) &
           +u2(2*i-1,2*j  )*FU2(2*i-1,2*j  ) &
           +u2(2*i  ,2*j  )*FU2(2*i  ,2*j  )  
          u2(2*i  ,2*j  )=u2(2*i  ,2*j  )+u1(i,j)
          u2(2*i-1,2*j  )=u2(2*i-1,2*j  )+u1(i,j)
          u2(2*i  ,2*j-1)=u2(2*i-1,2*j-1)+u1(i,j)
          u2(2*i-1,2*j-1)=u2(2*i  ,2*j-1)+u1(i,j)
          u1(i,j)=u1(i,j)+t/(FU1(i,j)+Eps)
        END DO
      END DO
  END SELECT
END SUBROUTINE CopyFace

SUBROUTINE PackFace(Buff,iBuf,u1,FU1,CopyCase) 

  REAL(RealKind) :: Buff(*)
  INTEGER :: iBuf
  REAL(RealKind) :: u1(:,:),FU1(:,:)
  INTEGER :: CopyCase

  INTEGER :: i,j
  INTEGER :: iU1,jU1

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          iBuf=iBuf+1
          Buff(iBuf)=u1(i-1,j-1)*FU1(i-1,j-1) &
                    +u1(i  ,j-1)*FU1(i  ,j-1) &
                    +u1(i-1,j  )*FU1(i-1,j  ) &
                    +u1(i  ,j  )*FU1(i  ,j  )  
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i-1,j  )*FU1(i-1,j  ) &
                    +u1(i  ,j  )*FU1(i  ,j  )  
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          iBuf=iBuf+1
          Buff(iBuf)=u1(i  ,j-1)*FU1(i  ,j-1) &
                    +u1(i  ,j  )*FU1(i  ,j  )  
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
  END SELECT
END SUBROUTINE PackFace

SUBROUTINE UnpackFace(u1,FU1,Buff,iBuf,CopyCase) 

  REAL(RealKind) :: u1(:,:),FU1(:,:),Buff(*)
  INTEGER :: iBuf
  INTEGER :: CopyCase

  INTEGER :: i,j
  INTEGER :: iU1,jU1
  REAL(RealKind) :: Temp

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          iBuf=iBuf+1
          Temp=Buff(iBuf)
          u1(i  ,j  )=u1(i  ,j  )+Temp
          u1(i-1,j  )=u1(i-1,j  )+Temp
          u1(i  ,j-1)=u1(i  ,j-1)+Temp
          u1(i-1,j-1)=u1(i-1,j-1)+Temp
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          iBuf=iBuf+1
          Temp=Buff(iBuf)
          u1(i  ,j)=u1(i  ,j)+Temp
          u1(i-1,j)=u1(i-1,j)+Temp
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          iBuf=iBuf+1
          Temp=Buff(iBuf)
          u1(i,j  )=u1(i,j  )+Temp
          u1(i,j-1)=u1(i,j-1)+Temp
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)/(FU1(i,j)+Eps)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)/(FU1(i,j)+Eps)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)/(FU1(i,j)+Eps)
        END DO
      END DO
  END SELECT
END SUBROUTINE UnpackFace

SUBROUTINE CopyFlux(c,cN,F) 

  REAL(RealKind) :: c(:,:,:),cN(:,:,:),F(:,:)

  INTEGER :: i,j
  INTEGER :: iU,jU,iUN,jUN
  REAL(RealKind) :: FC(LBOUND(c,3):UBOUND(c,3))

  iU=UBOUND(c,1)
  jU=UBOUND(c,2)
  iUN=UBOUND(cN,1)
  jUN=UBOUND(cN,2)

  SELECT CASE(2*(iU/iUN)+jU/jUN)
    CASE(FineFine)
      DO i=2,iU,2
        DO j=2,jU,2
          FC(:)=cN(i/2,j/2,:)/(F(i,j)+F(i-1,j)+F(i,j-1)+F(i-1,j-1)+Eps)
          c(i  ,j  ,:)=c(i  ,j  ,:)+F(i  ,j  )*FC 
          c(i-1,j  ,:)=c(i-1,j  ,:)+F(i-1,j  )*FC 
          c(i  ,j-1,:)=c(i  ,j-1,:)+F(i  ,j-1)*FC 
          c(i-1,j-1,:)=c(i-1,j-1,:)+F(i-1,j-1)*FC 
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU,2
        DO j=1,jU
          FC(:)=cN(i/2,j,:)/(F(i,j)+F(i-1,j)+Eps)
          c(i  ,j  ,:)=c(i  ,j  ,:)+F(i  ,j  )*FC
          c(i-1,j  ,:)=c(i-1,j  ,:)+F(i-1,j  )*FC
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU
        DO j=2,jU,2
          FC(:)=cN(i,j/2,:)/(F(i,j)+F(i,j-1)+Eps)
          c(i  ,j  ,:)=c(i  ,j  ,:)+F(i  ,j  )*FC
          c(i  ,j-1,:)=c(i  ,j-1,:)+F(i  ,j-1)*FC
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=c(i,j,:)+cN(i,j,:)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=c(i,j,:)+cN(i,2*j-1,:)+cN(i,2*j,:)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=c(i,j,:)+cN(2*i-1,j,:)+cN(2*i,j,:)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=c(i,j,:)+cN(2*i-1,2*j-1,:)+cN(2*i-1,2*j,:) &
                       +cN(2*i,2*j-1,:)+cN(2*i,2*j,:)
        END DO
      END DO
  END SELECT

END SUBROUTINE CopyFlux

SUBROUTINE PackFlux(Buff,iBuf,cN,CopyCase) 

  REAL(RealKind) :: Buff(*)
  INTEGER :: iBuf
  REAL(RealKind) :: cN(:,:,:)
  INTEGER :: CopyCase

  INTEGER :: i,j,k
  INTEGER :: iUN,jUN,it1,it2

  iUN=UBOUND(cN,1)
  jUN=UBOUND(cN,2)
  it1=LBOUND(cN,3)
  it2=UBOUND(cN,3)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,2*iUN,2
        DO j=2,2*jUN,2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i/2,j/2,k)
          END DO
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,2*iUN,2
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i/2,j,k)
          END DO
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iUN
        DO j=2,2*jUN,2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i,j/2,k)
          END DO
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iUN
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i,j,k)
          END DO
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iUN
        DO j=1,jUN/2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i,2*j-1,k)+cN(i,2*j,k)
          END DO
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iUN/2
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(2*i-1,j,k)+cN(2*i,j,k)
          END DO
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iUN/2
        DO j=1,jUN/2
          DO k=it1,it2
           iBuf=iBuf+1
           Buff(iBuf)=cN(2*i-1,2*j-1,k)+cN(2*i-1,2*j,k) &
                     +cN(2*i,2*j-1,k)+cN(2*i,2*j,k)
          END DO
        END DO
      END DO
  END SELECT
END SUBROUTINE PackFlux

SUBROUTINE UnpackFlux(c,Buff,iBuf,F,CopyCase) 

  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: Buff(*)
  INTEGER :: iBuf
  REAL(RealKind) :: F(:,:)
  INTEGER :: CopyCase

  INTEGER :: i,j,k
  INTEGER :: iU,jU,it1,it2
  REAL(RealKind) :: FC

  iU=UBOUND(c,1)
  jU=UBOUND(c,2)
  it1=LBOUND(c,3)
  it2=UBOUND(c,3)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU,2
        DO j=2,jU,2
          DO k=it1,it2
            iBuf=iBuf+1
            FC=Buff(iBuf)/(F(i,j)+F(i-1,j)+F(i,j-1)+F(i-1,j-1)+Eps)
            c(i  ,j  ,k)=c(i  ,j  ,k)+F(i  ,j  )*FC 
            c(i-1,j  ,k)=c(i-1,j  ,k)+F(i-1,j  )*FC 
            c(i  ,j-1,k)=c(i  ,j-1,k)+F(i  ,j-1)*FC 
            c(i-1,j-1,k)=c(i-1,j-1,k)+F(i-1,j-1)*FC 
          END DO
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU,2
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            FC=Buff(iBuf)/(F(i,j)+F(i-1,j)+Eps)
            c(i  ,j  ,k)=c(i  ,j  ,k)+F(i  ,j  )*FC
            c(i-1,j  ,k)=c(i-1,j  ,k)+F(i-1,j  )*FC
          END DO
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU
        DO j=2,jU,2
          DO k=it1,it2
            iBuf=iBuf+1
            FC=Buff(iBuf)/(F(i,j)+F(i,j-1)+Eps)
            c(i  ,j  ,k)=c(i  ,j  ,k)+F(i  ,j  )*FC
            c(i  ,j-1,k)=c(i  ,j-1,k)+F(i  ,j-1)*FC
          END DO
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=c(i,j,k)+Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=c(i,j,k)+Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=c(i,j,k)+Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=c(i,j,k)+Buff(iBuf) 
          END DO
        END DO
      END DO
    END SELECT
END SUBROUTINE UnpackFlux

SUBROUTINE CopyCell(c,cN,Vol) 

REAL(RealKind) :: c(:,:,:),cN(:,:,:),Vol(:,:)

INTEGER :: i,j
INTEGER :: iU,jU,iUN,jUN
REAL(RealKind) :: VC(LBOUND(c,3):UBOUND(c,3))

iU=UBOUND(c,1)
jU=UBOUND(c,2)
iUN=UBOUND(cN,1)
jUN=UBOUND(cN,2)

  SELECT CASE(2*(iU/iUN)+jU/jUN)
    CASE(6)
!     Fine Fine
      DO i=2,iU,2
        DO j=2,jU,2
          c(i  ,j  ,:)=cN(i/2,j/2,:)
          c(i-1,j  ,:)=cN(i/2,j/2,:)
          c(i  ,j-1,:)=cN(i/2,j/2,:)
          c(i-1,j-1,:)=cN(i/2,j/2,:)
        END DO
      END DO
!   END CASE(6)

    CASE(5)
!     Fine Equal
      DO i=2,iU,2
        DO j=1,jU
          c(i  ,j  ,:)=cN(i/2,j,:)
          c(i-1,j  ,:)=cN(i/2,j,:)
        END DO
      END DO
!   END CASE(5)

    CASE(4)
!     Equal Fine
      DO i=1,iU
        DO j=2,jU,2
          c(i  ,j  ,:)=cN(i,j/2,:)
          c(i  ,j-1,:)=cN(i,j/2,:)
        END DO
      END DO
!   END CASE(4)

    CASE(EqualEqual)
!     Equal Equal
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=cN(i,j,:)
        END DO
      END DO
!   END CASE(EqualEqual)

    CASE(EqualCoarse)
!     Equal Coarse
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=(cN(i,2*j-1,:)*Vol(i,2*j-1)+cN(i,2*j,:)*Vol(i,2*j)) &
                /(Vol(i,2*j-1)+Vol(i,2*j)+Eps) 
        END DO
      END DO
!   END CASE(EqualCoarse)

    CASE(CoarseEqual)
!     Coarse Equal
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=(cN(2*i-1,j,:)*Vol(2*i-1,j)+cN(2*i,j,:)*Vol(2*i,j)) &
                /(Vol(2*i-1,j)+Vol(2*i,j)+Eps)
        END DO
      END DO
!   END CASE(CoarseEqual)

    CASE(CoarseCoarse)
!     Coarse Coarse
      DO i=1,iU
        DO j=1,jU
          c(i,j,:)=(cN(2*i-1,2*j-1,:)*Vol(2*i-1,2*j-1) &
                 +cN(2*i-1,2*j,:)*Vol(2*i-1,2*j) &
                 +cN(2*i,2*j-1,:)*Vol(2*i,2*j-1) &
                 +cN(2*i,2*j,:)*Vol(2*i,2*j)) &
                /(Vol(2*i-1,2*j-1)+Vol(2*i-1,2*j) &
                 +Vol(2*i,2*j-1)+Vol(2*i,2*j)+Eps)
        END DO
      END DO
!   END CASE(CoarseEqual)
  END SELECT

END SUBROUTINE CopyCell

SUBROUTINE PackCell(Buff,iBuf,cN,Vol,CopyCase) 

  REAL(RealKind) :: Buff(*)
  INTEGER :: iBuf
  REAL(RealKind) :: cN(:,:,:),Vol(:,:)
  INTEGER :: CopyCase

  INTEGER :: i,j,k
  INTEGER :: iUN,jUN,it1,it2

  iUN=UBOUND(cN,1)
  jUN=UBOUND(cN,2)
  it1=LBOUND(cN,3)
  it2=UBOUND(cN,3)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,2*iUN,2
        DO j=2,2*jUN,2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i/2,j/2,k)
          END DO
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,2*iUN,2
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i/2,j,k)
          END DO
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iUN
        DO j=2,2*jUN,2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i,j/2,k)
          END DO
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iUN
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=cN(i,j,k)
          END DO
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iUN
        DO j=1,jUN/2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=(cN(i,2*j-1,k)*Vol(i,2*j-1)+cN(i,2*j,k)*Vol(i,2*j)) &
                  /(Vol(i,2*j-1)+Vol(i,2*j)+Eps) 
          END DO
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iUN/2
        DO j=1,jUN
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=(cN(2*i-1,j,k)*Vol(2*i-1,j)+cN(2*i,j,k)*Vol(2*i,j)) &
                /(Vol(2*i-1,j)+Vol(2*i,j)+Eps)
          END DO
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iUN/2
        DO j=1,jUN/2
          DO k=it1,it2
            iBuf=iBuf+1
            Buff(iBuf)=(cN(2*i-1,2*j-1,k)*Vol(2*i-1,2*j-1) &
                   +cN(2*i-1,2*j,k)*Vol(2*i-1,2*j) &
                   +cN(2*i,2*j-1,k)*Vol(2*i,2*j-1) &
                   +cN(2*i,2*j,k)*Vol(2*i,2*j)) &
                  /(Vol(2*i-1,2*j-1)+Vol(2*i-1,2*j) &
                   +Vol(2*i,2*j-1)+Vol(2*i,2*j)+Eps)
          END DO
        END DO
      END DO
  END SELECT
END SUBROUTINE PackCell

SUBROUTINE UnpackCell(c,Buff,iBuf,CopyCase) 

  REAL(RealKind) :: c(:,:,:),Buff(*)
  INTEGER :: iBuf
  INTEGER :: CopyCase

  INTEGER :: i,j,k
  INTEGER :: iU,jU,it1,it2

  iU=UBOUND(c,1)
  jU=UBOUND(c,2)
  it1=LBOUND(c,3)
  it2=UBOUND(c,3)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU,2
        DO j=2,jU,2
          DO k=it1,it2
            iBuf=iBuf+1
            c(i  ,j  ,k)=Buff(iBuf)
            c(i-1,j  ,k)=Buff(iBuf)
            c(i  ,j-1,k)=Buff(iBuf)
            c(i-1,j-1,k)=Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU,2
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i  ,j  ,k)=Buff(iBuf)
            c(i-1,j  ,k)=Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU
        DO j=2,jU,2
          DO k=it1,it2
            iBuf=iBuf+1
            c(i  ,j  ,k)=Buff(iBuf)
            c(i  ,j-1,k)=Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=Buff(iBuf)
          END DO
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=Buff(iBuf) 
          END DO
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=Buff(iBuf) 
          END DO
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU
        DO j=1,jU
          DO k=it1,it2
            iBuf=iBuf+1
            c(i,j,k)=Buff(iBuf) 
          END DO
        END DO
      END DO
  END SELECT

END SUBROUTINE UnpackCell

SUBROUTINE CopyFluxVel(u1,u2) 

  REAL(RealKind) :: u1(:,:),u2(:,:)

  INTEGER :: i,j
  INTEGER :: iU1,jU1,iU2,jU2
  REAL(RealKind) :: Temp

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)
  iU2=UBOUND(u2,1)
  jU2=UBOUND(u2,2)

  SELECT CASE(2*(iU1/iU2)+jU1/jU2)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          temp=u1(i-1,j-1) &
              +u1(i  ,j-1) &
              +u1(i-1,j  ) &
              +u1(i  ,j  )  
          u2(i/2,j/2)=u2(i/2,j/2)+temp
          u1(i  ,j  )=u2(i/2,j/2)
          u1(i-1,j  )=u2(i/2,j/2)
          u1(i  ,j-1)=u2(i/2,j/2)
          u1(i-1,j-1)=u2(i/2,j/2)
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          temp=u1(i-1,j) &
              +u1(i  ,j)  
          u2(i/2,j)=u2(i/2,j)+temp
          u1(i  ,j)=u2(i/2,j)
          u1(i-1,j)=u2(i/2,j)
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          temp=u1(i,j-1) &
              +u1(i,j  )  
          u2(i,j/2)=u2(i,j/2)+temp
          u1(i,j  )=u2(i,j/2)
          u1(i,j-1)=u2(i,j/2)
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          u1(i,j)=u1(i,j)+u2(i,j)
          u2(i,j)=u1(i,j)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          temp=u2(i,2*j-1) &
              +u2(i,2*j  )
          u1(i,j)=u1(i,j)+temp
          u2(i,2*j  )=u1(i,j)
          u2(i,2*j-1)=u1(i,j)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          temp=u2(2*i-1,j  ) &
              +u2(2*i  ,j  )
          u1(i,j)=u1(i,j)+temp
          u2(2*i  ,j)=u1(i,j)
          u2(2*i-1,j)=u1(i,j)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          temp=u2(2*i-1,2*j-1) &
              +u2(2*i  ,2*j-1) &
              +u2(2*i-1,2*j  ) &
              +u2(2*i  ,2*j  )  
          u1(i,j)=u1(i,j)+temp
          u2(2*i  ,2*j  )=u1(i,j)
          u2(2*i-1,2*j  )=u1(i,j)
          u2(2*i  ,2*j-1)=u1(i,j)
          u2(2*i-1,2*j-1)=u1(i,j)
        END DO
      END DO
  END SELECT
END SUBROUTINE CopyFluxVel

SUBROUTINE PackFluxVel(Buff,iBuf,u1,CopyCase) 

  REAL(RealKind) :: Buff(*)
  INTEGER :: iBuf
  REAL(RealKind) :: u1(:,:)
  INTEGER :: CopyCase

  INTEGER :: i,j
  INTEGER :: iU1,jU1

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          iBuf=iBuf+1
          Buff(iBuf)=u1(i-1,j-1) &
                    +u1(i  ,j-1) &
                    +u1(i-1,j  ) &
                    +u1(i  ,j  )  
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i-1,j  ) &
                    +u1(i  ,j  )  
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          iBuf=iBuf+1
          Buff(iBuf)=u1(i  ,j-1) &
                    +u1(i  ,j  )  
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          Buff(iBuf)=u1(i,j)
        END DO
      END DO
  END SELECT
END SUBROUTINE PackFluxVel

SUBROUTINE UnpackFluxVel(u1,Buff,iBuf,CopyCase) 

  REAL(RealKind) :: u1(:,:),Buff(*)
  INTEGER :: iBuf
  INTEGER :: CopyCase

  INTEGER :: i,j
  INTEGER :: iU1,jU1
  REAL(RealKind) :: Temp

  iU1=UBOUND(u1,1)
  jU1=UBOUND(u1,2)

  SELECT CASE(CopyCase)
    CASE(FineFine)
      DO i=2,iU1,2
        DO j=2,jU1,2
          iBuf=iBuf+1
          Temp=u1(i-1,j-1) &
              +u1(i  ,j-1) &
              +u1(i-1,j  ) &
              +u1(i  ,j  )  
          Temp=Buff(iBuf)+Temp
          u1(i  ,j  )=Temp
          u1(i-1,j  )=Temp
          u1(i  ,j-1)=Temp
          u1(i-1,j-1)=Temp
        END DO
      END DO
    CASE(FineEqual)
      DO i=2,iU1,2
        DO j=1,jU1
          iBuf=iBuf+1
          Temp=u1(i-1,j) &
              +u1(i  ,j)  
          Temp=Buff(iBuf)+Temp
          u1(i  ,j)=Temp
          u1(i-1,j)=Temp
        END DO
      END DO
    CASE(EqualFine)
      DO i=1,iU1
        DO j=2,jU1,2
          iBuf=iBuf+1
          Temp=u1(i,j-1) &
              +u1(i,j  )  
          Temp=Buff(iBuf)+Temp
          u1(i,j  )=Temp
          u1(i,j-1)=Temp
        END DO
      END DO
    CASE(EqualEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)
        END DO
      END DO
    CASE(EqualCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)
        END DO
      END DO
    CASE(CoarseEqual)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)
        END DO
      END DO
    CASE(CoarseCoarse)
      DO i=1,iU1
        DO j=1,jU1
          iBuf=iBuf+1
          u1(i,j)=u1(i,j)+Buff(iBuf)
        END DO
      END DO
  END SELECT
END SUBROUTINE UnpackFluxVel

END MODULE Exchange_Mod
