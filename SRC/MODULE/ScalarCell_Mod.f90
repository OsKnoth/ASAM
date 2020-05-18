MODULE ScalarCell_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Sp_Mod
  USE Parameter_Mod

  IMPLICIT NONE

  TYPE ScalarCell_T
    REAL(RealKind), POINTER :: c(:,:,:,:)
    REAL(RealKind), POINTER :: cB(:,:)
  END TYPE ScalarCell_T
  TYPE VecScalarCell_T
    TYPE(ScalarCell_T), POINTER :: Vec(:)
  END TYPE VecScalarCell_T

  TYPE (ScalarCell_T), POINTER :: ScalarCellAct(:)

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_ScalarCell
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_ScalarCell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_ScalarCell
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE DeallocateScalarCell
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE AllocateScalarCell
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_ScalarCell
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_ScalarCell
  END INTERFACE
  INTERFACE Xpyz
    MODULE PROCEDURE Xpyz_ScalarCell
  END INTERFACE


CONTAINS

SUBROUTINE Xpyz_ScalarCell(ScalarCell1,ScalarCell2,ScalarCell3)

  TYPE (ScalarCell_T) :: ScalarCell1,ScalarCell2,ScalarCell3

  ScalarCell3%c=ScalarCell1%c+ScalarCell2%c

END SUBROUTINE Xpyz_ScalarCell

SUBROUTINE Axpy_ScalarCell(alpha,ScalarCell1,ScalarCell2)

  REAL(RealKind) :: alpha
  TYPE (ScalarCell_T) :: ScalarCell1,ScalarCell2

  ScalarCell2%c=alpha*ScalarCell1%c+ScalarCell2%c

END SUBROUTINE Axpy_ScalarCell

SUBROUTINE Ax1mx2py_ScalarCell(alpha,ScalarCell1,ScalarCell2,ScalarCell3)

  REAL(RealKind) :: alpha
  TYPE (ScalarCell_T) :: ScalarCell1,ScalarCell2,ScalarCell3

  ScalarCell3%c=alpha*(ScalarCell1%c-ScalarCell2%c)+ScalarCell3%c

END SUBROUTINE Ax1mx2py_ScalarCell



SUBROUTINE AllocateScalarCell(ScalarCell,nFrac)

  TYPE (ScalarCell_T) :: ScalarCell
  INTEGER, OPTIONAL :: nFrac
  IF (PRESENT(nFrac)) THEN
    ALLOCATE(ScalarCell%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,nFrac))
  ELSE
    ALLOCATE(ScalarCell%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1))
  END IF
  ScalarCell%c=0.0e0
  ALLOCATE(ScalarCell%cB(1:NumBoundCell,1:7))

END SUBROUTINE AllocateScalarCell

SUBROUTINE DeallocateScalarCell(ScalarCell)

  TYPE (ScalarCell_T) :: ScalarCell
  DEALLOCATE(ScalarCell%c)
  DEALLOCATE(ScalarCell%cB)

END SUBROUTINE DeallocateScalarCell
SUBROUTINE Copy_ScalarCell(ScalarCell1,ScalarCell2)

  TYPE (ScalarCell_T) :: ScalarCell1,ScalarCell2

  ScalarCell2%c=ScalarCell1%c

END SUBROUTINE Copy_ScalarCell
SUBROUTINE Scale_ScalarCell(alpha,ScalarCell1)

  REAL(RealKind) :: alpha
  TYPE (ScalarCell_T) :: ScalarCell1

  ScalarCell1%c=alpha*ScalarCell1%c

END SUBROUTINE Scale_ScalarCell
SUBROUTINE Value_ScalarCell(ScalarCell1,Value)

  IMPLICIT NONE

  TYPE(ScalarCell_T), INTENT(OUT)  :: ScalarCell1
  REAL(RealKind), INTENT(IN)  :: Value

  ScalarCell1%c=Value

END SUBROUTINE Value_ScalarCell
END MODULE ScalarCell_Mod

MODULE ScalarCellPar_Mod
  USE Parallel_Mod
  USE Floor_Mod
  USE Exchange_Mod

  USE ScalarCell_Mod,Mat=>SpDiag

  IMPLICIT NONE

  INTERFACE RandCopyFlux
    MODULE PROCEDURE RandCopyFlux_Cell
  END INTERFACE
  INTERFACE PackBuffFlux
    MODULE PROCEDURE PackBuffFlux_Cell
  END INTERFACE
  INTERFACE UnPackBuffFlux
    MODULE PROCEDURE UnPackBuffFlux_Cell
  END INTERFACE
  INTERFACE ExchangeFlux
    MODULE PROCEDURE ExchangeFlux_Cell
  END INTERFACE
  INTERFACE RandCopyCell
    MODULE PROCEDURE RandCopyCell_Cell
  END INTERFACE
  INTERFACE PackBuffCell
    MODULE PROCEDURE PackBuffCell_Cell
  END INTERFACE
  INTERFACE UnPackBuffCell
    MODULE PROCEDURE UnPackBuffCell_Cell
  END INTERFACE
  INTERFACE ExchangeCell
    MODULE PROCEDURE ExchangeCell_Cell
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_Cell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_Cell
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_Cell
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE Deallocate_Cell
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_Cell
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_Cell
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_Cell
  END INTERFACE
  INTERFACE Xpyz
    MODULE PROCEDURE Xpyz_Cell
  END INTERFACE
  INTERFACE Mult
    MODULE PROCEDURE Mult_Cell
  END INTERFACE
CONTAINS 

SUBROUTINE Mult_Cell(x,y)

  TYPE(ScalarCell_T), INTENT(IN) :: x(:)
  TYPE(ScalarCell_T), INTENT(INOUT)  :: y(:)


  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    y(ibLoc)%c=x(ibLoc)%c*y(ibLoc)%c
  END DO

END SUBROUTINE Mult_Cell

SUBROUTINE  Xpyz_Cell(x,y,z)

  IMPLICIT NONE

  TYPE(ScalarCell_T), INTENT(IN)  :: x(:)
  TYPE(ScalarCell_T), INTENT(IN)  :: y(:)
  TYPE(ScalarCell_T), INTENT(OUT)  :: z(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Xpyz(x(ibLoc),y(ibLoc),z(ibLoc))
  END DO

END SUBROUTINE  Xpyz_Cell

SUBROUTINE  Axpy_Cell(alpha,x,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(ScalarCell_T), INTENT(IN)  :: x(:)
  TYPE(ScalarCell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,x(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Axpy_Cell

SUBROUTINE  Ax1mx2py_Cell(alpha,x1,x2,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(ScalarCell_T), INTENT(IN)  :: x1(:),x2(:)
  TYPE(ScalarCell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Ax1mx2py(alpha,x1(ibLoc),x2(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Ax1mx2py_Cell


SUBROUTINE Scale_Cell(alpha,y)


  REAL(RealKind), INTENT(IN) :: alpha
  TYPE(ScalarCell_T), INTENT(INOUT) :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,y(ibLoc))
  END DO

END SUBROUTINE SCALE_Cell
SUBROUTINE  Copy_Cell(x,y)

  IMPLICIT NONE

  TYPE(ScalarCell_T), INTENT(IN)  :: x(:)
  TYPE(ScalarCell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Copy(x(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Copy_Cell

SUBROUTINE Value_Cell(x,Value)

  IMPLICIT NONE

  TYPE(ScalarCell_T), INTENT(OUT)  :: x(:)
  REAL(RealKind), INTENT(IN)  :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    x(ibLoc)=Value
  END DO

END SUBROUTINE Value_Cell

SUBROUTINE  Allocate_Cell(x,nFrac)

  IMPLICIT NONE

  TYPE(ScalarCell_T), POINTER :: x(:)
  INTEGER, OPTIONAL :: nFrac

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(x(ibLoc),nFrac)
  END DO

END SUBROUTINE  Allocate_Cell

SUBROUTINE  Deallocate_Cell(x)

  IMPLICIT NONE

  TYPE(ScalarCell_T), POINTER :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Deallocate(x(ibLoc))
  END DO
  DEALLOCATE(x)

END SUBROUTINE  Deallocate_Cell

SUBROUTINE RandCopyFlux_Cell(vec)

  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: in
  INTEGER jx,jy,jz

  REAL(RealKind), POINTER :: c(:,:,:,:)
  REAL(RealKind), POINTER :: cN(:,:,:,:)
  REAL(RealKind) :: FUC,FVC,FWC
  TYPE (Nachbar_T), POINTER :: Nachbar

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Vec(ibLoc)%c
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)

      IF (nType=='iw'.OR.nType=='ie') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                      cN(jNxO,jNy0+1:jNy1,jNz0+1:jNz1,:), &
                      FU(jxO,jy0+1:jy1,jz0+1:jz1))
      ELSE IF (nType=='is'.OR.nType=='in') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                      cN(jNx0+1:jNx1,jNyO,jNz0+1:jNz1,:), &
                      FV(jx0+1:jx1,jyO,jz0+1:jz1))
      ELSE IF (nType=='ib'.OR.nType=='it') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                      cN(jNx0+1:jNx1,jNy0+1:jNy1,jNzO,:), &
                      FW(jx0+1:jx1,jy0+1:jy1,jzO))
      END IF
    END DO   ! in
  END DO   ! ibLoc

END SUBROUTINE RandCopyFlux_Cell

SUBROUTINE PackBuffFlux_Cell(Chain,Buffer,vec)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len

  REAL(RealKind), POINTER :: cN(:,:,:,:)
  iBuf=0
  Buffer(:)=Zero
  len=SIZE(Chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    VolC=>Floor(ibC)%VolC
    IF (nType=='pw'.OR.nType=='pe') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackFlux(Buffer,iBuf, &
                    cN(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                    6-CopyCase)
    ELSE IF (nType=='ps'.OR.nType=='pn') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackFlux(Buffer,iBuf, &
                    cN(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                    6-CopyCase)
    ELSE IF (nType=='pb'.OR.nType=='pt') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackFlux(Buffer,iBuf, &
                    cN(jx0+1:jx1,jy0+1:jy1,jzO,:), &
                    6-CopyCase)
    END IF
  END DO ! i

END SUBROUTINE PackBuffFlux_Cell

SUBROUTINE UnPackBuffFlux_Cell(Chain,Buffer,vec)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: i
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len

  REAL(RealKind) :: FUC,FVC,FWC
  REAL(RealKind), POINTER :: c(:,:,:,:)
  INTEGER :: RefineNachbar

  iBuf=0
  len=SIZE(Chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    IF (nType=='pw'.OR.nType=='pe') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                      Buffer,iBuf,                &
                      FU(jxO,jy0+1:jy1,jz0+1:jz1), &
                      CopyCase)
    ELSE IF (nType=='ps'.OR.nType=='pn') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                      Buffer,iBuf,                &
                      FV(jx0+1:jx1,jyO,jz0+1:jz1), &
                      CopyCase)
    ELSE IF (nType=='pb'.OR.nType=='pt') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                      Buffer,iBuf,                &
                      FW(jx0+1:jx1,jy0+1:jy1,jzO), &
                      CopyCase)
    END IF
  END DO ! i

END SUBROUTINE UnPackBuffFlux_Cell

SUBROUTINE ExchangeFlux_Cell(b)

   
  TYPE (ScalarCell_T) :: b(:) 

  INTEGER :: i,ip,j   
  INTEGER :: tag
  TYPE(MPI_Request) :: req(2*NumberNeiProc)
  TYPE(MPI_Status) :: StatusArray(2*NumberNeiProc)

  REAL(RealKind) :: PutBuf(MaxPutBuf,NumberNeiProc)
  REAL(RealKind) :: GetBuf(MaxGetBuf,NumberNeiProc)
  

  CALL RandCopyFlux(b)
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL PackBuffFlux(Chain(ip)%glied,PutBuf(:,i),b)
    END IF
  END DO
                         
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=MyId+NumProcs*ip
       CALL MPI_IRECV(GetBuf(1,i),Chain(ip)%LenUnPack, &
                     MPI_RealKind,ip,tag,      &
                     MPI_COMM_WORLD,req(i),MPIErr)
    END IF
  END DO
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=ip+NumProcs*MyId
       CALL MPI_ISEND(PutBuf(1,i),Chain(ip)%LenPack, &
                    MPI_RealKind,ip,tag,      &
                    MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
    END IF
  END DO

  CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

  i=0
      
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL UnPackBuffFlux(Chain(ip)%glied,GetBuf(:,i),b)
    END IF
  END DO

END SUBROUTINE ExchangeFlux_Cell


SUBROUTINE RandCopyCell_Cell(vec)

!---------------------------------------------------------------
!---  Copy Concentrations from Neighboring Blocks
!---------------------------------------------------------------


  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: in
  INTEGER :: jx,jy,jz

  REAL(RealKind), POINTER :: c(:,:,:,:)
  REAL(RealKind), POINTER :: cN(:,:,:,:)
  TYPE (Nachbar_T), POINTER :: Nachbar
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Vec(ibLoc)%c
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      VolC=>Floor(ibN)%VolC

      IF (nType=='iw'.OR.nType=='ie') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyCell(c(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                      cN(jNxI,jNy0+1:jNy1,jNz0+1:jNz1,:), &
                      VolC(jNxI,jNy0+1:jNy1,jNz0+1:jNz1))
      ELSE IF (nType=='is'.OR.nType=='in') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyCell(c(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                      cN(jNx0+1:jNx1,jNyI,jNz0+1:jNz1,:), &
                      VolC(jNx0+1:jNx1,jNyI,jNz0+1:jNz1))
      ELSE IF (nType=='ib'.OR.nType=='it') THEN ! west or east neighbor
        cN=>Vec(ibNLoc)%c
        CALL CopyCell(c(jx0+1:jx1,jy0+1:jy1,jzO,:), &
                      cN(jNx0+1:jNx1,jNy0+1:jNy1,jNzI,:), &
                      VolC(jNx0+1:jNx1,jNy0+1:jNy1,jNzI))
      END IF

    END DO ! in
  END DO   ! ibLoc

END SUBROUTINE RandCopyCell_Cell


SUBROUTINE PackBuffCell_Cell(Chain,Buffer,vec)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len

  REAL(RealKind), POINTER :: cN(:,:,:,:)

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(Chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    IF (nType=='pw'.OR.nType=='pe') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackCell(Buffer,iBuf, &
                    cN(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                    VolC(jxI,jy0+1:jy1,jz0+1:jz1), &
                    6-CopyCase)
    ELSE IF (nType=='ps'.OR.nType=='pn') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackCell(Buffer,iBuf, &
                    cN(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                    VolC(jx0+1:jx1,jyI,jz0+1:jz1), &
                    6-CopyCase)
    ELSE IF (nType=='pb'.OR.nType=='pt') THEN
      cN=>Vec(ibCLoc)%c
      CALL PackCell(Buffer,iBuf, &
                    cN(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                    VolC(jx0+1:jx1,jy0+1:jy1,jzI), &
                    6-CopyCase)
    END IF
  END DO ! i

END SUBROUTINE PackBuffCell_Cell

SUBROUTINE UnPackBuffCell_Cell(Chain,Buffer,vec)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (ScalarCell_T), TARGET :: vec(:)

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len

  REAL(RealKind), POINTER :: c(:,:,:,:)

  iBuf=0
  len=SIZE(Chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    IF (nType=='pw'.OR.nType=='pe') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackCell(c(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                      Buffer,iBuf,                &
                      CopyCase)
    ELSE IF (nType=='ps'.OR.nType=='pn') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackCell(c(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                      Buffer,iBuf,                &
                      CopyCase)
    ELSE IF (nType=='pb'.OR.nType=='pt') THEN
      c=>Vec(ibCLoc)%c
      CALL UnpackCell(c(jx0+1:jx1,jy0+1:jy1,jzO,:), &
                      Buffer,iBuf,                &
                      CopyCase)
    END IF
  END DO ! i

END SUBROUTINE UnPackBuffCell_Cell

SUBROUTINE ExchangeCell_Cell(b)

  TYPE (ScalarCell_T), POINTER :: b(:) 

  INTEGER :: i,ip,j   
  INTEGER :: tag
  TYPE(MPI_Request) :: req(2*NumberNeiProc)
  TYPE(MPI_Status) :: StatusArray(2*NumberNeiProc)

  REAL(RealKind) :: PutBuf(MaxPutBuf,NumberNeiProc)
  REAL(RealKind) :: GetBuf(MaxGetBuf,NumberNeiProc)
  
  IF (ASSOCIATED(b)) THEN
    CALL RandCopyCell(b)
    i=0
    DO ip=0,NumProcs-1
      IF (Chain(ip)%nlink > 0) THEN
        i=i+1
        CALL PackBuffCell(Chain(ip)%glied,PutBuf(:,i),b)
      END IF
    END DO
                         
    i=0
    DO ip=0,NumProcs-1
      IF (Chain(ip)%nlink > 0) THEN
        i=i+1
        tag=MyId+NumProcs*ip
         CALL MPI_IRECV(GetBuf(1,i),Chain(ip)%LenUnPack, &
                       MPI_RealKind,ip,tag,      &
                       MPI_COMM_WORLD,req(i),MPIErr)
      END IF
    END DO
    i=0
    DO ip=0,NumProcs-1
      IF (Chain(ip)%nlink > 0) THEN
        i=i+1
        tag=ip+NumProcs*MyId
         CALL MPI_ISEND(PutBuf(1,i),Chain(ip)%LenPack, &
                      MPI_RealKind,ip,tag,      &
                      MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
      END IF
    END DO
  
    CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

    i=0
      
    DO ip=0,NumProcs-1
      IF (Chain(ip)%nlink > 0) THEN
        i=i+1
        CALL UnPackBuffCell(Chain(ip)%glied,GetBuf(:,i),b)
      END IF
    END DO
  END IF

END SUBROUTINE ExchangeCell_Cell

END MODULE ScalarCellPar_Mod

