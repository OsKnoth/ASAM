MODULE Koagulation_Mod

  USE Kind_Mod
  USE DataType_Mod
  USE Chemie_Mod
  USE Ansatz_Mod
  USE InitAerosol_Mod
  USE Kernel_Mod, KernelF=>FUCH


  IMPLICIT NONE

  REAL(RealKind), PRIVATE, ALLOCATABLE :: Kernel(:,:)

CONTAINS

SUBROUTINE Koagulation(c,fc)

  TYPE (Vec4_T), TARGET :: c(0:)
  TYPE (Vec4_T) :: fc(0:)

  INTEGER :: i,j,k,l,iNzr
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: iFrac,jFrac,kFrac
  REAL(RealKind) :: Temp,Tempk1,Tempk2,Tempk3,Tempk4
  REAL(RealKind) :: mGes(nFrac),mmGes,dmGes

  REAL(RealKind), POINTER :: dNN(:)
  REAL(RealKind), POINTER :: NN(:)
  REAL(RealKind) :: A,b
  REAL(RealKind) :: NNLoc,mLoc,miLoc,mjLoc
  REAL(RealKind) :: mL,mU,mL1,mU1
  REAL(RealKind) :: mLoc1,mLoc2,NNLoc1

  REAL(RealKind), PARAMETER :: Eps=1.d-40

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        NN=>c(iNC)%c(ix,iy,iz,:)
        IF (SUM(NN)>Zero) THEN
          dNN=>fc(iNC)%c(ix,iy,iz,:)
          DO k=1,nFrac
            mGes(k)=Zero
            DO l=2,nAqua
              mGes(k)=mGes(k)+c(l)%c(ix,iy,iz,k)
            END DO
          END DO
          mmGes=SUM(mGes)
 
          DO i=1,nFrac
            DO j=i,nFrac
              NNLoc=Kernel(i,j)*NN(i)*NN(j)
              IF (NNLoc>0.0d0) THEN
                miLoc=Kernel(i,j)*mGes(i)*NN(j)
                mjLoc=Kernel(i,j)*mGes(j)*NN(i)
                mLoc=miLoc+mjLoc
                mL=m(i)+m(j)
                mU=m(i+1)+m(j+1)
                CALL AnsatzLogLogK(NNLoc,mLoc, &
                                   mL,mU,A,b)
                CALL IntLog(NNLoc,mLoc1,A,b,mL,mU,mL)
                mLoc1=mLoc
                mLoc2=Zero
                NNLoc1=Zero
                DO k=j,nFrac
                  mL1=MAX(m(k),mL)
                  mU1=MIN(m(k+1),mU)
                  IF (mL1<mU1) THEN
                    CALL IntLog(NNLoc,mLoc,A,b,mL1,mU1,mL)
                    dNN(k)=dNN(k)+NNloc
                    NNLoc1=NNLoc1+NNLoc
                    DO l=2,nAqua
                      kFrac=(c(l)%c(ix,iy,iz,i)*c(1)%c(ix,iy,iz,j) &
                            +c(l)%c(ix,iy,iz,j)*c(1)%c(ix,iy,iz,i)) &
                            /(mGes(i)*c(1)%c(ix,iy,iz,j) &
                             +mGes(j)*c(1)%c(ix,iy,iz,i)+Eps) 
                      fc(l)%c(ix,iy,iz,k)=fc(l)%c(ix,iy,iz,k)+mLoc*kFrac
                      mLoc2=mLoc2+mLoc*kFrac
                    END DO
                  END IF
                END DO  
                dNN(i)=dNN(i)-NNLoc1
                dNN(j)=dNN(j)-NNLoc1
                DO l=2,nAqua
                  iFrac=c(l)%c(ix,iy,iz,i)/(mges(i)+Eps)
                  jFrac=c(l)%c(ix,iy,iz,j)/(mges(j)+Eps)
                  fc(l)%c(ix,iy,iz,i)=fc(l)%c(ix,iy,iz,i)-miLoc*iFrac*mLoc2/(mLoc1+Eps)
                  fc(l)%c(ix,iy,iz,j)=fc(l)%c(ix,iy,iz,j)-mjLoc*jFrac*mLoc2/(mLoc1+Eps)
                END DO
              END IF
            END DO  
          END DO  
        END IF
      END DO  
    END DO  
  END DO  
      
END SUBROUTINE Koagulation
      
SUBROUTINE InitKoagulation

  CALL KernelCompute(nFrac,m)

END SUBROUTINE InitKoagulation

SUBROUTINE KernelCompute(n,x)

  INTEGER :: n
  REAL(RealKind) :: x(0:n)

  INTEGER :: i,j,ii,jj
  REAL(RealKind) :: xiM,xjM
! REAL(RealKind) :: Kernelf

! Computation of the kernel
  ALLOCATE(Kernel(n,n))
  DO i=1,n
    xiM=(x(i-1)+x(i))/2.0d0
    DO j=1,n
      xjM=(x(j-1)+x(j))/2.0d0
      Kernel(i,j)=0.0d0
      DO ii=1,10
        DO jj=1,10
          xiM=x(i-1)+(ii-0.5d0)*(x(i)-x(i-1))/10
          xjM=x(j-1)+(jj-0.5d0)*(x(j)-x(j-1))/10
          Kernel(i,j)=Kernel(i,j)+KernelF(xiM,xjM)
!         Kernel(i,j)=Kernel(i,j)+Long_cgs(xiM,xjM)/(xiM+xjM)
        END DO
      END DO
      Kernel(i,j)=Kernel(i,j)/100
    END DO
  END DO

END SUBROUTINE KernelCompute

END MODULE Koagulation_Mod
