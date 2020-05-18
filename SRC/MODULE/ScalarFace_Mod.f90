MODULE ScalarFace_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Parameter_Mod
  USE Vector4Cell_Mod
  IMPLICIT NONE

  TYPE ScalarFace_T
    REAL(RealKind), POINTER :: uF(:,:,:,:)
    REAL(RealKind), POINTER :: vF(:,:,:,:)
    REAL(RealKind), POINTER :: wF(:,:,:,:)
  END TYPE ScalarFace_T

  TYPE VectorSFace_T
    TYPE (ScalarFace_T), POINTER :: VecF(:)
  END TYPE VectorSFace_T

  TYPE VecScalarFace_T
    TYPE (ScalarFace_T), POINTER :: VecF(:)
  END TYPE VecScalarFace_T

  TYPE VecVectorSFace_T
    TYPE (VectorSFace_T), POINTER :: VecF(:)
  END TYPE VecVectorSFace_T

  TYPE (ScalarFace_T), POINTER :: ScalarFaceAct(:)

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_ScalarFace &
                    ,Value_VectorSFace
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE ScalarFaceAllocate &
                    ,VectorSFaceAllocate &
                    ,ScalarFaceVectorCellAllocate &
                    ,VectorSFaceVector4CellAllocate
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE ScalarFaceDeallocate
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_ScalarFace
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_ScalarFace &
                    ,Axpy_VectorSFace
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_ScalarFace
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_ScalarFace
  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_ScalarFace
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_ScalarFace &
                    ,Scale_VectorSFace
  END INTERFACE

CONTAINS

FUNCTION Dot_ScalarFace(ScalarFace1,ScalarFace2)

  REAL(RealKind) :: Dot_ScalarFace
  TYPE (ScalarFace_T) :: ScalarFace1,ScalarFace2

  Dot_ScalarFace=SUM(ScalarFace1%uF*ScalarFace2%uF) &
                   +SUM(ScalarFace1%vF*ScalarFace2%vF) &
                   +SUM(ScalarFace1%wF*ScalarFace2%wF) 

END FUNCTION Dot_ScalarFace

SUBROUTINE Copy_ScalarFace(ScalarFace1,ScalarFace2)

  TYPE (ScalarFace_T) :: ScalarFace1,ScalarFace2

  ScalarFace2%uF=ScalarFace1%uF
  ScalarFace2%vF=ScalarFace1%vF
  ScalarFace2%wF=ScalarFace1%wF

END SUBROUTINE Copy_ScalarFace

SUBROUTINE Axpy_ScalarFace(alpha,ScalarFace1,ScalarFace2)

  REAL(RealKind) :: alpha
  TYPE (ScalarFace_T) :: ScalarFace1,ScalarFace2

  IF (alpha/=Zero) THEN
    ScalarFace2%uF=alpha*ScalarFace1%uF+ScalarFace2%uF
    ScalarFace2%vF=alpha*ScalarFace1%vF+ScalarFace2%vF
    ScalarFace2%wF=alpha*ScalarFace1%wF+ScalarFace2%wF
  END IF

END SUBROUTINE Axpy_ScalarFace


SUBROUTINE Axpy_VectorSFace(alpha,VectorSFace1,VectorSFace2)

  REAL(RealKind) :: alpha
  TYPE (VectorSFace_T) :: VectorSFace1,VectorSFace2

  INTEGER :: ic

  IF (alpha/=Zero) THEN
    DO ic=LBOUND(VectorSFace1%VecF,1),UBOUND(VectorSFace1%VecF,1)
      VectorSFace2%VecF(ic)%uF=alpha*VectorSFace1%VecF(ic)%uF+VectorSFace2%VecF(ic)%uF
      VectorSFace2%VecF(ic)%vF=alpha*VectorSFace1%VecF(ic)%vF+VectorSFace2%VecF(ic)%vF
      VectorSFace2%VecF(ic)%wF=alpha*VectorSFace1%VecF(ic)%wF+VectorSFace2%VecF(ic)%wF
    END DO
  END IF

END SUBROUTINE Axpy_VectorSFace

SUBROUTINE Ax1mx2py_ScalarFace(alpha,ScalarFace1,ScalarFace2,ScalarFace3)

  REAL(RealKind) :: alpha
  TYPE (ScalarFace_T) :: ScalarFace1,ScalarFace2,ScalarFace3

  IF (alpha/=Zero) THEN
    ScalarFace3%uF=alpha*(ScalarFace1%uF-ScalarFace2%uF)+ScalarFace3%uF
    ScalarFace3%vF=alpha*(ScalarFace1%vF-ScalarFace2%vF)+ScalarFace3%vF
    ScalarFace3%wF=alpha*(ScalarFace1%wF-ScalarFace2%wF)+ScalarFace3%wF
  END IF

END SUBROUTINE Ax1mx2py_ScalarFace


SUBROUTINE Scale_ScalarFace(alpha,ScalarFace)

  REAL(RealKind) :: alpha
  TYPE (ScalarFace_T) :: ScalarFace

  ScalarFace%uF=alpha*ScalarFace%uF
  ScalarFace%vF=alpha*ScalarFace%vF
  ScalarFace%wF=alpha*ScalarFace%wF

END SUBROUTINE Scale_ScalarFace


SUBROUTINE Scale_VectorSFace(alpha,VectorSFace)

  REAL(RealKind) :: alpha
  TYPE (VectorSFace_T) :: VectorSFace

  INTEGER :: ic

  DO ic=LBOUND(VectorSFace%VecF,1),UBOUND(VectorSFace%VecF,1)
    VectorSFace%VecF(ic)%uF=alpha*VectorSFace%VecF(ic)%uF
    VectorSFace%VecF(ic)%vF=alpha*VectorSFace%VecF(ic)%vF
    VectorSFace%VecF(ic)%wF=alpha*VectorSFace%VecF(ic)%wF
  END DO

END SUBROUTINE Scale_VectorSFace


SUBROUTINE Xpay_ScalarFace(ScalarFace1,alpha,ScalarFace2)

  REAL(RealKind) :: alpha
  TYPE (ScalarFace_T) :: ScalarFace1,ScalarFace2

  ScalarFace2%uF=ScalarFace1%uF+alpha*ScalarFace2%uF
  ScalarFace2%vF=ScalarFace1%vF+alpha*ScalarFace2%vF
  ScalarFace2%wF=ScalarFace1%wF+alpha*ScalarFace2%wF

END SUBROUTINE Xpay_ScalarFace

SUBROUTINE Value_ScalarFace(ScalarFace,Value)

  TYPE (ScalarFace_T), INTENT(OUT) :: ScalarFace
  REAL(RealKind), INTENT(IN) :: Value

  ScalarFace%uF=Value
  ScalarFace%vF=Value
  ScalarFace%wF=Value

END SUBROUTINE Value_ScalarFace

SUBROUTINE Value_VectorSFace(VectorSFace,Value)

  TYPE (VectorSFace_T), INTENT(OUT) :: VectorSFace
  REAL(RealKind), INTENT(IN) :: Value

  INTEGER :: ic

  DO ic=LBOUND(VectorSFace%VecF,1),UBOUND(VectorSFace%VecF,1)
    VectorSFace%VecF(ic)%uF=Value
    VectorSFace%VecF(ic)%vF=Value
    VectorSFace%VecF(ic)%wF=Value
  END DO

END SUBROUTINE Value_VectorSFace

SUBROUTINE ScalarFaceAllocate(ScalarFace)

  TYPE (ScalarFace_T) :: ScalarFace
  ALLOCATE(ScalarFace%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1,1))
  ALLOCATE(ScalarFace%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1,1))
  ALLOCATE(ScalarFace%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1,1))

END SUBROUTINE ScalarFaceAllocate

  
SUBROUTINE ScalarFaceVectorCellAllocate(ScalarFace,Vec4)
  
  TYPE (ScalarFace_T) :: ScalarFace
  TYPE (Vec4_T) :: Vec4

  INTEGER :: it1,it2
  it1=LBOUND(Vec4%c,4)
  it2=UBOUND(Vec4%c,4)
  ALLOCATE(ScalarFace%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1,it1:it2))
  ALLOCATE(ScalarFace%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1,it1:it2))
  ALLOCATE(ScalarFace%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1,it1:it2))

END SUBROUTINE ScalarFaceVectorCellAllocate


SUBROUTINE VectorSFaceAllocate(VectorSFace,VectorComponents1,VectorComponents2)

  TYPE (VectorSFace_T) :: VectorSFace
  INTEGER :: VectorComponents1,VectorComponents2

  INTEGER :: ic

  ALLOCATE(VectorSFace%VecF(VectorComponents1:VectorComponents2))
  DO ic=LBOUND(VectorSFace%VecF,1),UBOUND(VectorSFace%VecF,1)
    ALLOCATE(VectorSFace%VecF(ic)%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1,1))
    ALLOCATE(VectorSFace%VecF(ic)%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1,1))
    ALLOCATE(VectorSFace%VecF(ic)%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1,1))
  END DO

END SUBROUTINE VectorSFaceAllocate

SUBROUTINE VectorSFaceVector4CellAllocate(VectorSFace,Vector4Cell,VectorComponents1,VectorComponents2)

  TYPE (VectorSFace_T) :: VectorSFace
  TYPE(Vector4Cell_T) :: Vector4Cell
  INTEGER, OPTIONAL :: VectorComponents1,VectorComponents2

  INTEGER :: ic,it1,it2,VectorComponentsLoc1,VectorComponentsLoc2

  IF (PRESENT(VectorComponents1)) THEN
    VectorComponentsLoc1=VectorComponents1
    VectorComponentsLoc2=VectorComponents2
  ELSE
    VectorComponentsLoc1=LBOUND(Vector4Cell%Vec,1)
    VectorComponentsLoc2=UBOUND(Vector4Cell%Vec,1)
  END IF

  ALLOCATE(VectorSFace%VecF(VectorComponentsLoc1:VectorComponentsLoc2))
  DO ic=LBOUND(VectorSFace%VecF,1),UBOUND(VectorSFace%VecF,1)
    it1=LBOUND(Vector4Cell%Vec(ic)%c,4)
    it2=UBOUND(Vector4Cell%Vec(ic)%c,4)
    ALLOCATE(VectorSFace%VecF(ic)%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1,it1:it2))
    ALLOCATE(VectorSFace%VecF(ic)%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1,it1:it2))
    ALLOCATE(VectorSFace%VecF(ic)%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1,it1:it2))
  END DO

END SUBROUTINE VectorSFaceVector4CellAllocate


SUBROUTINE ScalarFaceDeallocate(ScalarFace)

  TYPE (ScalarFace_T) :: ScalarFace
  DEALLOCATE(ScalarFace%uF)
  DEALLOCATE(ScalarFace%vF)
  DEALLOCATE(ScalarFace%wF)

END SUBROUTINE ScalarFaceDeallocate

SUBROUTINE ScalarFaceBounds(ScalarFace)

  TYPE (ScalarFace_T) :: ScalarFace

  ix0=LBOUND(ScalarFace%uF,1)
  ix1=UBOUND(ScalarFace%uF,1)
  iy0=LBOUND(ScalarFace%vF,2)
  iy1=UBOUND(ScalarFace%vF,2)
  iz0=LBOUND(ScalarFace%wF,3)
  iz1=UBOUND(ScalarFace%wF,3)

END SUBROUTINE ScalarFaceBounds

FUNCTION NumberVecF(x)

  INTEGER :: NumberVecF
  TYPE(VectorSFace_T) :: x

  INTEGER :: ic

  NumberVecF=0
  DO ic=LBOUND(x%VecF,1),UBOUND(x%VecF,1)
    NumberVecF=NumberVecF+SIZE(x%VecF(ic)%uF,4)
  END DO
END FUNCTION NumberVecF

END MODULE ScalarFace_Mod

MODULE ScalarFacePar_Mod

  USE ScalarFace_Mod
  USE Floor_Mod
  USE Parallel_Mod
  USE Vector4CellPar_Mod
  IMPLICIT NONE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_ScaFace &
                    ,Value_VecSFace
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE ScalarFaceParAllocate &
                    ,VectorSFaceParAllocate &
                    ,VectorSFaceVector4CellParAllocate
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE ScalarFaceParDeallocate 
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_ScaFace
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_ScaFace &
                    ,Axpy_VecFace
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_ScaFace
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_ScaFace &
                    ,Scale_VecFace
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_ScaFace
  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_ScaFace
  END INTERFACE


CONTAINS 

SUBROUTINE ScalarFaceParAllocate(Scalar)

  TYPE (ScalarFace_T), POINTER :: Scalar(:)

  ALLOCATE(Scalar(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(Scalar(ibLoc))
  END DO
END SUBROUTINE ScalarFaceParAllocate

SUBROUTINE VectorSFaceParAllocate(Vector,VectorComponents1,VectorComponents2)

  TYPE (VectorSFace_T), POINTER :: Vector(:)
  INTEGER :: VectorComponents1,VectorComponents2

  ALLOCATE(Vector(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(Vector(ibLoc),VectorComponents1,VectorComponents2)
  END DO
END SUBROUTINE VectorSFaceParAllocate

SUBROUTINE VectorSFaceVector4CellParAllocate(VectorSFace,Vector4Cell &
                                           ,VectorComponents1,VectorComponents2)

  TYPE (VectorSFace_T), POINTER :: VectorSFace(:)
  TYPE(Vector4Cell_T), POINTER :: Vector4Cell(:)
  INTEGER, OPTIONAL :: VectorComponents1,VectorComponents2

  ALLOCATE(VectorSFace(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (PRESENT(VectorComponents1)) THEN
      CALL Allocate(VectorSFace(ibLoc),Vector4Cell(ibLoc),VectorComponents1,VectorComponents2)
    ELSE
      CALL Allocate(VectorSFace(ibLoc),Vector4Cell(ibLoc),VectorComponents1,VectorComponents2)
    END IF
  END DO
END SUBROUTINE VectorSFaceVector4CellParAllocate


SUBROUTINE ScalarFaceParDeallocate(Scalar)

  TYPE (ScalarFace_T), POINTER :: Scalar(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Deallocate(Scalar(ibLoc))
  END DO
  DEALLOCATE(Scalar)
END SUBROUTINE ScalarFaceParDeallocate

FUNCTION Dot_Scaface(Scalar1,Scalar2)

  REAL(RealKind) :: Dot_Scaface
  TYPE(ScalarFace_T), INTENT(IN)  :: Scalar1(:),Scalar2(:)

  REAL(RealKind) :: Temp

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp=Temp+Dot(Scalar1(ibLoc),Scalar2(ibLoc))
  END DO
  Dot_Scaface=Temp

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,Dot_Scaface,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)

END FUNCTION Dot_Scaface


SUBROUTINE  Copy_ScaFace(Scalar1,Scalar2)

  TYPE(ScalarFace_T), INTENT(IN)  :: Scalar1(:)
  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Copy(Scalar1(ibLoc),Scalar2(ibLoc))
  END DO

END SUBROUTINE  Copy_ScaFace

SUBROUTINE  Axpy_ScaFace(alpha,Scalar1,Scalar2)

  REAL(RealKind) :: alpha
  TYPE(ScalarFace_T), INTENT(IN)  :: Scalar1(:)
  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,Scalar1(ibLoc),Scalar2(ibLoc))
  END DO

END SUBROUTINE  Axpy_ScaFace

SUBROUTINE  Axpy_VecFace(alpha,Vector1,Vector2)

  REAL(RealKind) :: alpha
  TYPE(VectorSFace_T), INTENT(IN)  :: Vector1(:)
  TYPE(VectorSFace_T), INTENT(INOUT)  :: Vector2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,Vector1(ibLoc),Vector2(ibLoc))
  END DO

END SUBROUTINE  Axpy_VecFace


SUBROUTINE  Ax1mx2py_ScaFace(alpha,Scalar1,Scalar2,Scalar3)

  REAL(RealKind) :: alpha
  TYPE(ScalarFace_T), INTENT(IN)  :: Scalar1(:),Scalar2(:)
  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar3(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Ax1mx2py(alpha,Scalar1(ibLoc),Scalar2(ibLoc),Scalar3(ibLoc))
  END DO

END SUBROUTINE  Ax1mx2py_ScaFace


SUBROUTINE  Scale_ScaFace(alpha,Scalar)

  REAL(RealKind) :: alpha
  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,Scalar(ibLoc))
  END DO

END SUBROUTINE  Scale_ScaFace

SUBROUTINE  Scale_VecFace(alpha,Vector)

  REAL(RealKind) :: alpha
  TYPE(VectorSFace_T), INTENT(INOUT)  :: Vector(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,Vector(ibLoc))
  END DO

END SUBROUTINE  Scale_VecFace


SUBROUTINE  Xpay_ScaFace(Scalar1,alpha,Scalar2)

  REAL(RealKind) :: alpha
  TYPE(ScalarFace_T), INTENT(IN)  :: Scalar1(:)
  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Xpay(Scalar1(ibLoc),alpha,Scalar2(ibLoc))
  END DO

END SUBROUTINE  Xpay_ScaFace

SUBROUTINE  Value_ScaFace(Scalar,Value)

  TYPE(ScalarFace_T), INTENT(INOUT)  :: Scalar(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Scalar(ibLoc)=Value
  END DO

END SUBROUTINE  Value_ScaFace

SUBROUTINE  Value_VecSFace(Vector,Value)

  TYPE(VectorSFace_T), INTENT(INOUT)  :: Vector(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Vector(ibLoc)=Value
  END DO

END SUBROUTINE  Value_VecSFace


SUBROUTINE RandCopyScalarFace(Vector,NumberOfComponents,Offset)

!---------------------------------------------------------------
!---  Copy Face Values from Neighboring Blocks
!---------------------------------------------------------------

  TYPE(VectorSFace_T), TARGET, INTENT(INOUT)  :: Vector(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: in,ic
  INTEGER :: jx,jy,jz
  TYPE(VectorSFace_T), POINTER  :: subblock
  TYPE (Nachbar_T), POINTER :: Nachbar

  DO ic=1+Offset,NumberOfComponents+Offset
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      subblock=>Vector(ibLoc)
      DO in=1,AnzahlNachbar
        Nachbar=>Nachbars(in)
        CALL Set(Nachbar)
        IF (Nachbar%nType=='iw'.AND.ibLoc<=ibnLoc) THEN           ! west neighbor
          CALL CopyFace(subblock%VecF(ic)%uF(ix0,jy0+1:jy1,jz0+1:jz1,1),FU(ix0,jy0+1:jy1,jz0+1:jz1), &
                        Vector(ibnLoc)%VecF(ic)%uF(jNx1,jNy0+1:jNy1,jNz0+1:jNz1,1),Floor(ibn)%WeiFU(jNx1,jNy0+1:jNy1,jNz0+1:jNz1))
        ELSE IF (Nachbar%nType=='ie'.AND.ibLoc<ibnLoc) THEN       ! east neighbor
          CALL CopyFace(subblock%VecF(ic)%uF(ix1,jy0+1:jy1,jz0+1:jz1,1),FU(ix1,jy0+1:jy1,jz0+1:jz1), &
                        Vector(ibnLoc)%VecF(ic)%uF(jNx0,jNy0+1:jNy1,jNz0+1:jNz1,1),Floor(ibn)%WeiFU(jNx0,jNy0+1:jNy1,jNz0+1:jNz1))
        ELSE IF (Nachbar%nType=='is'.AND.ibLoc<=ibnLoc) THEN      ! south neighbor
          CALL CopyFace(subblock%VecF(ic)%vF(jx0+1:jx1,iy0,jz0+1:jz1,1),FV(jx0+1:jx1,iy0,jz0+1:jz1),  &
                        Vector(ibnLoc)%VecF(ic)%vF(jNx0+1:jNx1,jNy1,jNz0+1:jNz1,1),Floor(ibn)%WeiFV(jNx0+1:jNx1,jNy1,jNz0+1:jNz1))
        ELSE IF (Nachbar%nType=='in'.AND.ibLoc<ibnLoc) THEN       ! north neighbor
          CALL CopyFace(subblock%VecF(ic)%vF(jx0+1:jx1,iy1,jz0+1:jz1,1),FV(jx0+1:jx1,iy1,jz0+1:jz1),  &
                        Vector(ibnLoc)%VecF(ic)%vF(jNx0+1:jNx1,jNy0,jNz0+1:jNz1,1),Floor(ibn)%WeiFV(jNx0+1:jNx1,jNy0,jNz0+1:jNz1))
        ELSE IF (Nachbar%nType=='ib'.AND.ibLoc<=ibnLoc) THEN      ! bottom neighbor
          CALL CopyFace(subblock%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,iz0,1),FW(jx0+1:jx1,jy0+1:jy1,iz0), &
                        Vector(ibnLoc)%VecF(ic)%wF(jNx0+1:jNx1,jNy0+1:jNy1,jNz1,1),Floor(ibn)%WeiFW(jNx0+1:jNx1,jNy0+1:jNy1,jNz1))
        ELSE IF (Nachbar%nType=='it'.AND.ibLoc<ibnLoc) THEN       ! top neighbor
          CALL CopyFace(subblock%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,iz1,1),FW(jx0+1:jx1,jy0+1:jy1,iz1), &
                        Vector(ibnLoc)%VecF(ic)%wF(jNx0+1:jNx1,jNy0+1:jNy1,jNz0,1),Floor(ibn)%WeiFW(jNx0+1:jNx1,jNy0+1:jNy1,jNz0))
        END IF
      END DO ! in
    END DO   ! ibLoc
  END DO
END SUBROUTINE RandCopyScalarFace


SUBROUTINE ExchangeScalarFace(Vector,NumberOfComponents)
   
  TYPE (VectorSFace_T) :: Vector(:) 
  INTEGER, OPTIONAL :: NumberOfComponents

  INTEGER :: i,ip,j   
  INTEGER :: tag
  TYPE(MPI_Request) :: req(2*NumberNeiProc)
  TYPE(MPI_Status) :: StatusArray(2*NumberNeiProc)
  INTEGER :: LocNumberOfComponents
  INTEGER :: LocBufNumberOfComponents
  INTEGER :: Offset

  REAL(RealKind), ALLOCATABLE :: PutBuf(:,:)
  REAL(RealKind), ALLOCATABLE :: GetBuf(:,:)
   
  IF (PRESENT(NumberOfComponents)) THEN
    LocNumberOfComponents=NumberOfComponents
  ELSE
    IF (nbLoc>0) THEN
      LocNumberOfComponents=SIZE(Vector(1)%VecF)
    ELSE
      LocNumberOfComponents=0
    END IF
  END IF
  IF (nbLoc>0) THEN
    LocBufNumberOfComponents=NumberVecF(Vector(1))
    Offset=LBOUND(Vector(1)%VecF,1)-1
  ELSE 
    LocBufNumberOfComponents=0
    Offset=0
  END IF
  
  ALLOCATE(PutBuf(LocBufNumberOfComponents*MaxPutBuf,NumberNeiProc))
  ALLOCATE(GetBuf(LocBufNumberOfComponents*MaxGetBuf,NumberNeiProc))
  CALL RandCopyScalarFace(Vector,LocNumberOfComponents,Offset)
  i=0
  DO ip=0,NumProcs-1
    IF (chain(ip)%nlink > 0) THEN
      i=i+1
      CALL PackBuff(chain(ip)%glied,PutBuf(:,i),Vector,LocNumberOfComponents,Offset)
    END IF
  END DO
                         
  i=0
  DO ip=0,NumProcs-1
    IF (chain(ip)%nlink > 0) THEN
      i=i+1
      tag=MyId+NumProcs*ip
      CALL MPI_IRECV(GetBuf(1,i),LocBufNumberOfComponents*chain(ip)%LenUnPack, &
                     MPI_RealKind,ip,tag,      &
                     MPI_COMM_WORLD,req(i),MPIErr)
    END IF
  END DO
  i=0
  DO ip=0,NumProcs-1
    IF (chain(ip)%nlink > 0) THEN
      i=i+1
      tag=ip+NumProcs*MyId
      CALL MPI_ISEND(PutBuf(1,i),LocBufNumberOfComponents*chain(ip)%LenPack, &
                     MPI_RealKind,ip,tag,      &
                     MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
    END IF
  END DO

  CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

  i=0
  DO ip=0,NumProcs-1
    IF (chain(ip)%nlink > 0) THEN
      i=i+1
      CALL UnPackBuff(chain(ip)%glied,GetBuf(:,i),Vector,LocNumberOfComponents,Offset)
    END IF
  END DO
  DEALLOCATE(PutBuf)
  DEALLOCATE(GetBuf)


CONTAINS 

SUBROUTINE PackBuff(chain,Buffer,Vector,NumberOfComponents,Offset)

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (VectorSFace_T) :: Vector(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: iBuf
  INTEGER :: i,len,ic

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      IF (nType=='pw') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%uF(jxO,jy0+1:jy1,jz0+1:jz1,1),FU(jxO,jy0+1:jy1,jz0+1:jz1),CopyCase)
      ELSE IF (nType=='pe') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%uF(jxI,jy0+1:jy1,jz0+1:jz1,1),FU(jxI,jy0+1:jy1,jz0+1:jz1),CopyCase)
      ELSE IF (nType=='ps') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%vF(jx0+1:jx1,jyO,jz0+1:jz1,1),FV(jx0+1:jx1,jyO,jz0+1:jz1),CopyCase)
      ELSE IF (nType=='pn') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%vF(jx0+1:jx1,jyI,jz0+1:jz1,1),FV(jx0+1:jx1,jyI,jz0+1:jz1),CopyCase)
      ELSE IF (nType=='pb') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,jzO,1),FW(jx0+1:jx1,jy0+1:jy1,jzO),CopyCase)
      ELSE IF (nType=='pt') THEN
        CALL PackFace(Buffer,iBuf,Vector(ibCLoc)%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,jzI,1),FW(jx0+1:jx1,jy0+1:jy1,jzI),CopyCase)
      END IF
    END DO
  END DO
END SUBROUTINE PackBuff

SUBROUTINE UnPackBuff(chain,Buffer,Vector,NumberOfComponents,Offset)

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (VectorSFace_T) :: Vector(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: iBuf
  INTEGER :: i,len,ic

  !WRITE(700+MyId,*) 'Len Buff',SIZE(Buffer,1)
  iBuf=0
  len=SIZE(chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      IF (nType=='pw') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%uF(jxO,jy0+1:jy1,jz0+1:jz1,1),FU(jxO,jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
      ELSE IF (nType=='pe') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%uF(jxI,jy0+1:jy1,jz0+1:jz1,1),FU(jxI,jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
      ELSE IF (nType=='ps') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%vF(jx0+1:jx1,jyO,jz0+1:jz1,1),FV(jx0+1:jx1,jyO,jz0+1:jz1),Buffer,iBuf,CopyCase)
      ELSE IF (nType=='pn') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%vF(jx0+1:jx1,jyI,jz0+1:jz1,1),FV(jx0+1:jx1,jyI,jz0+1:jz1),Buffer,iBuf,CopyCase)
      ELSE IF (nType=='pb') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,jzO,1),FW(jx0+1:jx1,jy0+1:jy1,jzO),Buffer,iBuf,CopyCase)
      ELSE IF (nType=='pt') THEN
        CALL UnpackFace(Vector(ibCLoc)%VecF(ic)%wF(jx0+1:jx1,jy0+1:jy1,jzI,1),FW(jx0+1:jx1,jy0+1:jy1,jzI),Buffer,iBuf,CopyCase)
      END IF
    END DO
  END DO
END SUBROUTINE UnPackBuff

END SUBROUTINE ExchangeScalarFace 

END MODULE ScalarFacePar_Mod

