MODULE Vector4CellPar_Mod
  USE Parallel_Mod
  USE Exchange_Mod
  USE Floor_Mod

  USE Vector4Cell_Mod

  IMPLICIT NONE

  INTERFACE RandCopyCell
    MODULE PROCEDURE RandCopyCell_VecCell
  END INTERFACE
  INTERFACE PackBuffCell
    MODULE PROCEDURE PackBuffCell_VecCell
  END INTERFACE
  INTERFACE UnPackBuffCell
    MODULE PROCEDURE UnPackBuffCell_VecCell
  END INTERFACE
  INTERFACE ExchangeCell
    MODULE PROCEDURE ExchangeCell_VecCell
  END INTERFACE
  INTERFACE RandCopyFlux
    MODULE PROCEDURE RandCopyFlux_VecCell
  END INTERFACE
  INTERFACE PackBuffFlux
    MODULE PROCEDURE PackBuffFlux_VecCell
  END INTERFACE
  INTERFACE UnPackBuffFlux
    MODULE PROCEDURE UnPackBuffFlux_VecCell
  END INTERFACE
  INTERFACE ExchangeFlux
    MODULE PROCEDURE ExchangeFlux_VecCell
  END INTERFACE
  
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_VecCell
  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_VecCell
  END INTERFACE
  INTERFACE Sum
    MODULE PROCEDURE Sum_VecCell
  END INTERFACE
  INTERFACE Error
    MODULE PROCEDURE Error_VecCell
  END INTERFACE
  INTERFACE ErrorMax
    MODULE PROCEDURE ErrorMax_VecCell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_VecCell
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_VecCell
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_VecCell
  END INTERFACE
  INTERFACE Ax1px2py
    MODULE PROCEDURE Ax1px2py_VecCell
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_VecCell
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_VecCell
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_VecCell
  END INTERFACE
  INTERFACE Random
    MODULE PROCEDURE Random_VecCell
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE Deallocate_VecCell
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_VecCell1 &
                    ,Allocate_VecCell2 &
                    ,Allocate_VecCellVecCell
  END INTERFACE
CONTAINS 

FUNCTION Dot_VecCell(x,y)

  REAL(RealKind) :: Dot_VecCell
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:),y(:)

  REAL(RealKind) :: Temp

  Temp=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp=Temp+Dot(x(ibLoc),y(ibLoc))
  END DO
  Dot_VecCell=Temp

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,Dot_VecCell,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)

END FUNCTION Dot_VecCell

FUNCTION Sum_VecCell(x,Pos)

  REAL(RealKind) :: Sum_VecCell
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:)
  INTEGER, INTENT(IN)  :: Pos

  REAL(RealKind) :: Temp

  Temp=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Temp=Temp+Sum(x(ibLoc),Pos)
  END DO
  Sum_VecCell=Temp

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,Sum_VecCell,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)

END FUNCTION Sum_VecCell

FUNCTION Error_VecCell(x,y,z,ATol,RTol,t)

  REAL(RealKind) :: Error_VecCell
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:),y(:),z(:),t(:)
  REAL(RealKind) :: ATol(:),RTol(:)

  REAL(RealKind) :: LenVec,LenVecLoc
  REAL(RealKind) :: Temp

  Temp=0.0d0
  LenVecLoc=0.0d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Temp=Temp+Error(x(ibLoc),y(ibLoc),z(ibLoc),ATol,RTol,LenVecLoc,t(ibLoc))
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(LenVecLoc,LenVec,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(Temp,Error_VecCell,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Error_VecCell=SQRT(Error_VecCell/LenVec)

END FUNCTION Error_VecCell

FUNCTION ErrorMax_VecCell(x,y,z,ATol,RTol,t)

  REAL(RealKind) :: ErrorMax_VecCell
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:),y(:),z(:),t(:)
  REAL(RealKind) :: ATol(:),RTol(:)

  REAL(RealKind) :: Temp

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Temp=MAX(Temp,ErrorMax(x(ibLoc),y(ibLoc),z(ibLoc),ATol,RTol,t(ibLoc)))
  END DO

! -- Global max over all processors --
  CALL MPI_Allreduce(Temp,ErrorMax_VecCell,1,MPI_RealKind, &
&                    MPI_MAX,MPI_Comm_World,MPIErr)

END FUNCTION ErrorMax_VecCell

SUBROUTINE Scale_VecCell(alpha,y)


  REAL(RealKind), INTENT(IN) :: alpha
  TYPE(Vector4Cell_T), INTENT(INOUT) :: y(:)


  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,y(ibLoc))
  END DO

END SUBROUTINE SCALE_VecCell
SUBROUTINE  Axpy_VecCell(alpha,x,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,x(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Axpy_VecCell
SUBROUTINE  Ax1mx2py_VecCell(alpha,x1,x2,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(Vector4Cell_T), INTENT(IN)  :: x1(:),x2(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Ax1mx2py(alpha,x1(ibLoc),x2(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Ax1mx2py_VecCell
SUBROUTINE  Ax1px2py_VecCell(alpha,x1,x2,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(Vector4Cell_T), INTENT(IN)  :: x1(:),x2(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Ax1px2py(alpha,x1(ibLoc),x2(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Ax1px2py_VecCell
SUBROUTINE  Axpby_VecCell(alpha,x,beta,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha,beta
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpby(alpha,x(ibLoc),beta,y(ibLoc))
  END DO

END SUBROUTINE  Axpby_VecCell
SUBROUTINE  Xpay_VecCell(x,alpha,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha
  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Xpay(x(ibLoc),alpha,y(ibLoc))
  END DO

END SUBROUTINE  Xpay_VecCell
SUBROUTINE  Copy_VecCell(x,y)

  IMPLICIT NONE

  TYPE(Vector4Cell_T), INTENT(IN)  :: x(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Copy(x(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Copy_VecCell

SUBROUTINE Value_VecCell(x,Value)

  IMPLICIT NONE

  TYPE(Vector4Cell_T), INTENT(OUT) :: x(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    x(ibLoc)=Value
  END DO

END SUBROUTINE Value_VecCell

SUBROUTINE Random_VecCell(x)

  IMPLICIT NONE

  TYPE(Vector4Cell_T), INTENT(OUT)  :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CAll Random(x(ibLoc))
  END DO

END SUBROUTINE Random_VecCell

SUBROUTINE  Allocate_VecCell1(x,VectorComponents)

  TYPE(Vector4Cell_T), POINTER :: x(:)
  INTEGER :: VectorComponents

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(x(ibLoc),VectorComponents)
  END DO

END SUBROUTINE  Allocate_VecCell1

SUBROUTINE  Allocate_VecCell2(x,VectorComponents1,VectorComponents2)

  TYPE(Vector4Cell_T), POINTER :: x(:)
  INTEGER :: VectorComponents1,VectorComponents2

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(x(ibLoc),VectorComponents1,VectorComponents2)
  END DO

END SUBROUTINE  Allocate_VecCell2


SUBROUTINE  Allocate_VecCellVecCell(x,xO,VectorComponents1,VectorComponents2)

  TYPE(Vector4Cell_T), POINTER :: x(:)
  TYPE(Vector4Cell_T), POINTER :: xO(:)
  INTEGER, OPTIONAL :: VectorComponents1,VectorComponents2

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (PRESENT(VectorComponents1)) THEN
      CALL Allocate(x(ibLoc),xO(ibLoc),VectorComponents1,VectorComponents2)
    ELSE
      CALL Allocate(x(ibLoc),xO(ibLoc))
    END IF
  END DO

END SUBROUTINE  Allocate_VecCellVecCell

SUBROUTINE  Deallocate_VecCell(x)

  IMPLICIT NONE

  TYPE(Vector4Cell_T), POINTER :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Deallocate(x(ibLoc))
  END DO
  DEALLOCATE(x)

END SUBROUTINE  Deallocate_VecCell

SUBROUTINE RandCopyFlux_VecCell(VecP,NumberOfComponents,Offset)

  TYPE (Vector4Cell_T), TARGET :: vecP(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: in
  INTEGER :: jx,jy,jz
  INTEGER :: ic

  REAL(RealKind), POINTER :: c(:,:,:,:)
  REAL(RealKind), POINTER :: cN(:,:,:,:)
  REAL(RealKind) :: FUC,FVC,FWC
  TYPE (Nachbar_T), POINTER :: Nachbar

  DO ic=1+Offset,NumberOfComponents+Offset
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      c=>VecP(ibLoc)%Vec(ic)%c
      DO in=1,AnzahlNachbar
        Nachbar=>Nachbars(in)
        CALL Set(Nachbar)

        IF (nType=='iw') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                        cN(jNxO,jNy0+1:jNy1,jNz0+1:jNz1,:), &
                        FU(jxO,jy0+1:jy1,jz0+1:jz1))
        ELSE IF (nType=='ie') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                        cN(jNxO,jNy0+1:jNy1,jNz0+1:jNz1,:), &
                        FU(jxI,jy0+1:jy1,jz0+1:jz1))
        ELSE IF (nType=='is') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                        cN(jNx0+1:jNx1,jNyO,jNz0+1:jNz1,:), &
                        FV(jx0+1:jx1,jyO,jz0+1:jz1))
        ELSE IF (nType=='in') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                        cN(jNx0+1:jNx1,jNyO,jNz0+1:jNz1,:), &
                        FV(jx0+1:jx1,jyI,jz0+1:jz1))
        ELSE IF (nType=='ib') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                        cN(jNx0+1:jNx1,jNy0+1:jNy1,jNzO,:), &
                        FW(jx0+1:jx1,jy0+1:jy1,jzO))
        ELSE IF (nType=='it') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                        cN(jNx0+1:jNx1,jNy0+1:jNy1,jNzO,:), &
                        FW(jx0+1:jx1,jy0+1:jy1,jzI))
        END IF 
      END DO   ! in
    END DO   ! ibLoc
  END DO   ! ic

END SUBROUTINE RandCopyFlux_VecCell


SUBROUTINE PackBuffFlux_VecCell(Chain,Buffer,Vec,NumberOfComponents,Offset)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (Vector4Cell_T), TARGET :: vec(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len
  INTEGER :: ic 

  REAL(RealKind), POINTER :: cN(:,:,:,:)

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(Chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      VolC=>Floor(ibC)%VolC
      IF (nType=='pw'.OR.nType=='pe') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackFlux(Buffer,iBuf, &
                      cN(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                      6-CopyCase) 
      ELSE IF (nType=='ps'.OR.nType=='pn') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackFlux(Buffer,iBuf, &
                      cN(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                      6-CopyCase) 
      ELSE IF (nType=='pb'.OR.nType=='pt') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackFlux(Buffer,iBuf, &
                      cN(jx0+1:jx1,jy0+1:jy1,jzO,:), & 
                      6-CopyCase) 
      END IF
    END DO ! i
  END DO ! ic

END SUBROUTINE PackBuffFlux_VecCell

SUBROUTINE UnPackBuffFlux_VecCell(Chain,Buffer,Vec,NumberOfComponents,Offset)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (Vector4Cell_T), TARGET :: vec(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len
  INTEGER :: ic

  REAL(RealKind) :: FUC,FVC,FWC
  REAL(RealKind), POINTER :: c(:,:,:,:)

  iBuf=0
  len=SIZE(Chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      IF (nType=='pw') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        FU(jxO,jy0+1:jy1,jz0+1:jz1), &
                        CopyCase)
      ELSE IF (nType=='pe') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        FU(jxI,jy0+1:jy1,jz0+1:jz1), &
                        CopyCase)
      ELSE IF (nType=='ps') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        FV(jx0+1:jx1,jyO,jz0+1:jz1), &
                        CopyCase)
      ELSE IF (nType=='pn') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        FV(jx0+1:jx1,jyI,jz0+1:jz1), &
                        CopyCase)
      ELSE IF (nType=='pb') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                        Buffer,iBuf,                &
                        FW(jx0+1:jx1,jy0+1:jy1,jzO), &
                        CopyCase)
      ELSE IF (nType=='pt') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackFlux(c(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                        Buffer,iBuf,                &
                        FW(jx0+1:jx1,jy0+1:jy1,jzI), &
                        CopyCase)
      END IF
    END DO ! i
  END DO ! ic

END SUBROUTINE UnPackBuffFlux_VecCell

SUBROUTINE ExchangeFlux_VecCell(b,NumberOfComponents)

   
  TYPE (Vector4Cell_T) :: b(:) 
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
      LocNumberOfComponents=SIZE(b(1)%Vec)
    ELSE
      LocNumberOfComponents=0
    END IF
  END IF
  IF (nbLoc>0) THEN
    LocBufNumberOfComponents=NumberVec(b(1))
    Offset=LBOUND(b(1)%Vec,1)-1
  ELSE 
    LocBufNumberOfComponents=0
    Offset=0
  END IF
 
  ALLOCATE(PutBuf(LocBufNumberOfComponents*MaxPutBuf,NumberNeiProc))
  ALLOCATE(GetBuf(LocBufNumberOfComponents*MaxGetBuf,NumberNeiProc))

  CALL RandCopyFlux(b,LocNumberOfComponents,Offset)
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL PackBuffFlux(Chain(ip)%glied,PutBuf(:,i),b,LocNumberOfComponents,Offset)
    END IF
  END DO
                         
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=MyId+NumProcs*ip
       CALL MPI_IRECV(GetBuf(1,i),LocBufNumberOfComponents*Chain(ip)%LenUnPack, &
                     MPI_RealKind,ip,tag,      &
                     MPI_COMM_WORLD,req(i),MPIErr)
    END IF
  END DO
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=ip+NumProcs*MyId
       CALL MPI_ISEND(PutBuf(1,i),LocBufNumberOfComponents*Chain(ip)%LenPack, &
                    MPI_RealKind,ip,tag,      &
                    MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
    END IF
  END DO

  CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

  i=0
      
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL UnPackBuffFlux(Chain(ip)%glied,GetBuf(:,i),b,LocNumberOfComponents,Offset)
    END IF
  END DO

  DEALLOCATE(PutBuf)
  DEALLOCATE(GetBuf)

END SUBROUTINE ExchangeFlux_VecCell

SUBROUTINE RandCopyCell_VecCell(VecP,NumberOfComponents,Offset)

!---------------------------------------------------------------
!---  Copy Concentrations from Neighboring Blocks
!---------------------------------------------------------------

  TYPE (Vector4Cell_T), TARGET :: VecP(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: in
  INTEGER :: jx,jy,jz
  INTEGER :: ic

  REAL(RealKind), POINTER :: c(:,:,:,:)
  REAL(RealKind), POINTER :: cN(:,:,:,:)
  TYPE (Nachbar_T), POINTER :: Nachbar
  
  DO ic=1+Offset,NumberOfComponents+Offset
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      c=>VecP(ibLoc)%Vec(ic)%c
      DO in=1,AnzahlNachbar
        Nachbar=>Nachbars(in)
        CALL Set(Nachbar)
        VolC=>Floor(ibN)%VolC

        IF (nType=='iw'.OR.nType=='ie') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyCell(c(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                        cN(jNxI,jNy0+1:jNy1,jNz0+1:jNz1,:), &
                        VolC(jNxI,jNy0+1:jNy1,jNz0+1:jNz1))
        ELSE IF (nType=='is'.OR.nType=='in') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyCell(c(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                        cN(jNx0+1:jNx1,jNyI,jNz0+1:jNz1,:), &
                        VolC(jNx0+1:jNx1,jNyI,jNz0+1:jNz1))
        ELSE IF (nType=='ib'.OR.nType=='it') THEN ! west or east neighbor
          cN=>VecP(ibNLoc)%Vec(ic)%c
          CALL CopyCell(c(jx0+1:jx1,jy0+1:jy1,jzO,:), &
                        cN(jNx0+1:jNx1,jNy0+1:jNy1,jNzI,:), &
                        VolC(jNx0+1:jNx1,jNy0+1:jNy1,jNzI))
        END IF

      END DO ! in
    END DO   ! ibLoc
  END DO   ! ic

END SUBROUTINE RandCopyCell_VecCell

SUBROUTINE PackBuffCell_VecCell(Chain,Buffer,vec,NumberOfComponents,Offset)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (Vector4Cell_T), TARGET :: vec(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len

  REAL(RealKind), POINTER :: cN(:,:,:,:)
  INTEGER :: ic

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(Chain)
  iBuf=0
  Buffer(:)=Zero
  len=SIZE(Chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      IF (nType=='pw'.OR.nType=='pe') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackCell(Buffer,iBuf, &
                      cN(jxI,jy0+1:jy1,jz0+1:jz1,:), &
                      VolC(jxI,jy0+1:jy1,jz0+1:jz1), &
                      6-CopyCase)
      ELSE IF (nType=='ps'.OR.nType=='pn') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackCell(Buffer,iBuf, &
                      cN(jx0+1:jx1,jyI,jz0+1:jz1,:), &
                      VolC(jx0+1:jx1,jyI,jz0+1:jz1), &
                      6-CopyCase)
      ELSE IF (nType=='pb'.OR.nType=='pt') THEN
        cN=>Vec(ibCLoc)%Vec(ic)%c
        CALL PackCell(Buffer,iBuf, &
                      cN(jx0+1:jx1,jy0+1:jy1,jzI,:), &
                      VolC(jx0+1:jx1,jy0+1:jy1,jzI), &
                      6-CopyCase)
      END IF

    END DO ! i
  END DO ! ic

END SUBROUTINE PackBuffCell_VecCell

SUBROUTINE UnPackBuffCell_VecCell(Chain,Buffer,vec,NumberOfComponents,Offset)

  TYPE (gliedT) :: Chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (Vector4Cell_T), TARGET :: vec(:)
  INTEGER :: NumberOfComponents,Offset

  INTEGER :: i,ib
  INTEGER :: jx,jy,jz
  INTEGER :: iBuf,len
  INTEGER :: ic

  REAL(RealKind), POINTER :: c(:,:,:,:)

  iBuf=0
  len=SIZE(Chain)
  DO ic=1+Offset,NumberOfComponents+Offset
    DO i=1,len
      CALL Set(Chain(i))
      CALL Set(Floor(ibC))
      CALL Set(Nachbar)
      IF (nType=='pw'.OR.nType=='pe') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackCell(c(jxO,jy0+1:jy1,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        CopyCase)
      ELSE IF (nType=='ps'.OR.nType=='pn') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackCell(c(jx0+1:jx1,jyO,jz0+1:jz1,:), &
                        Buffer,iBuf,                &
                        CopyCase)
      ELSE IF (nType=='pb'.OR.nType=='pt') THEN
        c=>Vec(ibCLoc)%Vec(ic)%c
        CALL UnpackCell(c(jx0+1:jx1,jy0+1:jy1,jzO,:), &
                        Buffer,iBuf,                &
                        CopyCase)
      END IF
    END DO ! i
  END DO ! ic
END SUBROUTINE UnPackBuffCell_VecCell

SUBROUTINE ExchangeCell_VecCell(b,NumberOfComponents)

  TYPE (Vector4Cell_T) :: b(:) 
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
      LocNumberOfComponents=SIZE(b(1)%Vec)
    ELSE
      LocNumberOfComponents=0
    END IF
  END IF
  IF (nbLoc>0) THEN
    LocBufNumberOfComponents=NumberVec(b(1))
    Offset=LBOUND(b(1)%Vec,1)-1
  ELSE 
    LocBufNumberOfComponents=0
    Offset=0
  END IF
  
  ALLOCATE(PutBuf(LocBufNumberOfComponents*MaxPutBuf,NumberNeiProc))
  ALLOCATE(GetBuf(LocBufNumberOfComponents*MaxGetBuf,NumberNeiProc))
  CALL RandCopyCell(b,LocNumberOfComponents,Offset)
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL PackBuffCell(Chain(ip)%glied,PutBuf(:,i),b,LocNumberOfComponents,Offset)
    END IF
  END DO
                         
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=MyId+NumProcs*ip
       CALL MPI_IRECV(GetBuf(1,i),LocBufNumberOfComponents*Chain(ip)%LenUnPack, &
                     MPI_RealKind,ip,tag,      &
                     MPI_COMM_WORLD,req(i),MPIErr)
    END IF
  END DO
  i=0
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      tag=ip+NumProcs*MyId
       CALL MPI_ISEND(PutBuf(1,i),LocBufNumberOfComponents*Chain(ip)%LenPack, &
                    MPI_RealKind,ip,tag,      &
                    MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
    END IF
  END DO

  CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

  i=0
      
  DO ip=0,NumProcs-1
    IF (Chain(ip)%nlink > 0) THEN
      i=i+1
      CALL UnPackBuffCell(Chain(ip)%glied,GetBuf(:,i),b,LocNumberOfComponents,Offset)
    END IF
  END DO
  DEALLOCATE(PutBuf)
  DEALLOCATE(GetBuf)

END SUBROUTINE ExchangeCell_VecCell
END MODULE Vector4CellPar_Mod

