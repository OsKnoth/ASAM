MODULE VelocityFace_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Parameter_Mod
  USE Exchange_Mod
  IMPLICIT NONE

  TYPE VelocityFace_T
    REAL(RealKind), POINTER :: uF(:,:,:)
    REAL(RealKind), POINTER :: vF(:,:,:)
    REAL(RealKind), POINTER :: wF(:,:,:)
  END TYPE VelocityFace_T

  TYPE VectorFace_T
    TYPE (VelocityFace_T), POINTER :: VecF(:)
  END TYPE VectorFace_T

  TYPE VecVelocityFace_T
    TYPE (VelocityFace_T), POINTER :: VecF(:)
  END TYPE VecVelocityFace_T

  TYPE VecVectorFace_T
    TYPE (VectorFace_T), POINTER :: VecF(:)
  END TYPE VecVectorFace_T

  TYPE (VelocityFace_T), POINTER :: VelocityFaceAct(:)

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_VelocityFace &
                    ,Value_VectorFace
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE VelocityFaceAllocate &
                    ,VectorFaceAllocate
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE VelocityFaceDeallocate
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_VelocityFace
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_VelocityFace &
                    ,Axpy_VectorFace
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_VelocityFace
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_VelocityFace
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_VelocityFace
  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_VelocityFace
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_VelocityFace &
                    ,Scale_VectorFace
  END INTERFACE

CONTAINS

FUNCTION Dot_VelocityFace(VelocityFace1,VelocityFace2)

  REAL(RealKind) :: Dot_VelocityFace
  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2

  Dot_VelocityFace=SUM(VelocityFace1%uF*VelocityFace2%uF) &
                   +SUM(VelocityFace1%vF*VelocityFace2%vF) &
                   +SUM(VelocityFace1%wF*VelocityFace2%wF) 

END FUNCTION Dot_VelocityFace

SUBROUTINE Copy_VelocityFace(VelocityFace1,VelocityFace2)

  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2

  VelocityFace2%uF=VelocityFace1%uF
  VelocityFace2%vF=VelocityFace1%vF
  VelocityFace2%wF=VelocityFace1%wF

END SUBROUTINE Copy_VelocityFace

SUBROUTINE Axpy_VelocityFace(alpha,VelocityFace1,VelocityFace2)

  REAL(RealKind) :: alpha
  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2

  IF (alpha/=Zero) THEN
    VelocityFace2%uF=alpha*VelocityFace1%uF+VelocityFace2%uF
    VelocityFace2%vF=alpha*VelocityFace1%vF+VelocityFace2%vF
    VelocityFace2%wF=alpha*VelocityFace1%wF+VelocityFace2%wF
  END IF

END SUBROUTINE Axpy_VelocityFace


SUBROUTINE Axpy_VectorFace(alpha,VectorFace1,VectorFace2)

  REAL(RealKind) :: alpha
  TYPE (VectorFace_T) :: VectorFace1,VectorFace2

  INTEGER :: ic

  IF (alpha/=Zero) THEN
    DO ic=LBOUND(VectorFace1%VecF,1),UBOUND(VectorFace1%VecF,1)
      VectorFace2%VecF(ic)%uF=alpha*VectorFace1%VecF(ic)%uF+VectorFace2%VecF(ic)%uF
      VectorFace2%VecF(ic)%vF=alpha*VectorFace1%VecF(ic)%vF+VectorFace2%VecF(ic)%vF
      VectorFace2%VecF(ic)%wF=alpha*VectorFace1%VecF(ic)%wF+VectorFace2%VecF(ic)%wF
    END DO
  END IF

END SUBROUTINE Axpy_VectorFace

SUBROUTINE Ax1mx2py_VelocityFace(alpha,VelocityFace1,VelocityFace2,VelocityFace3)

  REAL(RealKind) :: alpha
  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2,VelocityFace3

  IF (alpha/=Zero) THEN
    VelocityFace3%uF=alpha*(VelocityFace1%uF-VelocityFace2%uF)+VelocityFace3%uF
    VelocityFace3%vF=alpha*(VelocityFace1%vF-VelocityFace2%vF)+VelocityFace3%vF
    VelocityFace3%wF=alpha*(VelocityFace1%wF-VelocityFace2%wF)+VelocityFace3%wF
  END IF

END SUBROUTINE Ax1mx2py_VelocityFace


SUBROUTINE Scale_VelocityFace(alpha,VelocityFace)

  REAL(RealKind) :: alpha
  TYPE (VelocityFace_T) :: VelocityFace

  VelocityFace%uF=alpha*VelocityFace%uF
  VelocityFace%vF=alpha*VelocityFace%vF
  VelocityFace%wF=alpha*VelocityFace%wF

END SUBROUTINE Scale_VelocityFace


SUBROUTINE Scale_VectorFace(alpha,VectorFace)

  REAL(RealKind) :: alpha
  TYPE (VectorFace_T) :: VectorFace

  INTEGER :: ic

  DO ic=LBOUND(VectorFace%VecF,1),UBOUND(VectorFace%VecF,1)
    VectorFace%VecF(ic)%uF=alpha*VectorFace%VecF(ic)%uF
    VectorFace%VecF(ic)%vF=alpha*VectorFace%VecF(ic)%vF
    VectorFace%VecF(ic)%wF=alpha*VectorFace%VecF(ic)%wF
  END DO

END SUBROUTINE Scale_VectorFace


SUBROUTINE Xpay_VelocityFace(VelocityFace1,alpha,VelocityFace2)

  REAL(RealKind) :: alpha
  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2

  VelocityFace2%uF=VelocityFace1%uF+alpha*VelocityFace2%uF
  VelocityFace2%vF=VelocityFace1%vF+alpha*VelocityFace2%vF
  VelocityFace2%wF=VelocityFace1%wF+alpha*VelocityFace2%wF

END SUBROUTINE Xpay_VelocityFace

SUBROUTINE Axpby_VelocityFace(alpha,VelocityFace1,beta,VelocityFace2)

  REAL(RealKind) :: alpha,beta
  TYPE (VelocityFace_T) :: VelocityFace1,VelocityFace2

  VelocityFace2%uF=alpha*VelocityFace1%uF+beta*VelocityFace2%uF
  VelocityFace2%vF=alpha*VelocityFace1%vF+beta*VelocityFace2%vF
  VelocityFace2%wF=alpha*VelocityFace1%wF+beta*VelocityFace2%wF

END SUBROUTINE Axpby_VelocityFace

SUBROUTINE Value_VelocityFace(VelocityFace,Value)

  TYPE (VelocityFace_T), INTENT(OUT) :: VelocityFace
  REAL(RealKind), INTENT(IN) :: Value

  VelocityFace%uF=Value
  VelocityFace%vF=Value
  VelocityFace%wF=Value

END SUBROUTINE Value_VelocityFace

SUBROUTINE Value_VectorFace(VectorFace,Value)

  TYPE (VectorFace_T), INTENT(OUT) :: VectorFace
  REAL(RealKind), INTENT(IN) :: Value

  INTEGER :: ic

  DO ic=LBOUND(VectorFace%VecF,1),UBOUND(VectorFace%VecF,1)
    VectorFace%VecF(ic)%uF=Value
    VectorFace%VecF(ic)%vF=Value
    VectorFace%VecF(ic)%wF=Value
  END DO

END SUBROUTINE Value_VectorFace

SUBROUTINE VelocityFaceAllocate(VelocityFace)

  TYPE (VelocityFace_T) :: VelocityFace
  ALLOCATE(VelocityFace%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1))
  ALLOCATE(VelocityFace%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1))
  ALLOCATE(VelocityFace%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1))

END SUBROUTINE VelocityFaceAllocate

SUBROUTINE VectorFaceAllocate(VectorFace,VectorComponents1,VectorComponents2)

  TYPE (VectorFace_T) :: VectorFace
  INTEGER :: VectorComponents1,VectorComponents2

  INTEGER :: ic

  ALLOCATE(VectorFace%VecF(VectorComponents1:VectorComponents2))
  DO ic=LBOUND(VectorFace%VecF,1),UBOUND(VectorFace%VecF,1)
    ALLOCATE(VectorFace%VecF(ic)%uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1))
    ALLOCATE(VectorFace%VecF(ic)%vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1))
    ALLOCATE(VectorFace%VecF(ic)%wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1))
  END DO

END SUBROUTINE VectorFaceAllocate

SUBROUTINE VelocityFaceDeallocate(VelocityFace)

  TYPE (VelocityFace_T) :: VelocityFace
  DEALLOCATE(VelocityFace%uF)
  DEALLOCATE(VelocityFace%vF)
  DEALLOCATE(VelocityFace%wF)

END SUBROUTINE VelocityFaceDeallocate

SUBROUTINE VelocityFaceBounds(VelocityFace)

  TYPE (VelocityFace_T) :: VelocityFace

  ix0=LBOUND(VelocityFace%uF,1)
  ix1=UBOUND(VelocityFace%uF,1)
  iy0=LBOUND(VelocityFace%vF,2)
  iy1=UBOUND(VelocityFace%vF,2)
  iz0=LBOUND(VelocityFace%wF,3)
  iz1=UBOUND(VelocityFace%wF,3)

END SUBROUTINE VelocityFaceBounds

END MODULE VelocityFace_Mod

MODULE VelocityFacePar_Mod

  USE VelocityFace_Mod
  USE Floor_Mod
  USE Parallel_Mod

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_VelFace &
                    ,Value_VecFace
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE VelocityFaceParAllocate &
                    ,VectorFaceParAllocate
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE VelocityFaceParDeallocate
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_VelFace
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_VelFace &
                    ,Axpy_VecFace
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_VelFace 
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_VelFace
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_VelFace &
                    ,Scale_VecFace
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_VelFace
  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_VelFace
  END INTERFACE


CONTAINS 

SUBROUTINE VelocityFaceParAllocate(Velocity)

  TYPE (VelocityFace_T), POINTER :: Velocity(:)

  ALLOCATE(Velocity(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(Velocity(ibLoc))
  END DO
END SUBROUTINE VelocityFaceParAllocate

SUBROUTINE VectorFaceParAllocate(Vector,VectorComponents1,VectorComponents2)

  TYPE (VectorFace_T), POINTER :: Vector(:)
  INTEGER :: VectorComponents1,VectorComponents2

  ALLOCATE(Vector(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(Vector(ibLoc),VectorComponents1,VectorComponents2)
  END DO
END SUBROUTINE VectorFaceParAllocate

SUBROUTINE VelocityFaceParDeallocate(Velocity)

  TYPE (VelocityFace_T), POINTER :: Velocity(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Deallocate(Velocity(ibLoc))
  END DO
  DEALLOCATE(Velocity)
END SUBROUTINE VelocityFaceParDeallocate

FUNCTION Dot_Velface(Velocity1,Velocity2)

  REAL(RealKind) :: Dot_Velface
  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:),Velocity2(:)

  REAL(RealKind) :: Temp

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp=Temp+Dot(Velocity1(ibLoc),Velocity2(ibLoc))
  END DO
  Dot_Velface=Temp

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,Dot_Velface,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)

END FUNCTION Dot_Velface


SUBROUTINE  Copy_VelFace(Velocity1,Velocity2)

  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:)
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Copy(Velocity1(ibLoc),Velocity2(ibLoc))
  END DO

END SUBROUTINE  Copy_VelFace

SUBROUTINE  Axpy_VelFace(alpha,Velocity1,Velocity2)

  REAL(RealKind) :: alpha
  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:)
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,Velocity1(ibLoc),Velocity2(ibLoc))
  END DO

END SUBROUTINE  Axpy_VelFace

SUBROUTINE  Axpby_VelFace(alpha,Velocity1,beta,Velocity2)

  REAL(RealKind) :: alpha,beta
  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:)
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpby(alpha,Velocity1(ibLoc),beta,Velocity2(ibLoc))
  END DO

END SUBROUTINE  Axpby_VelFace


SUBROUTINE  Axpy_VecFace(alpha,Vector1,Vector2)

  REAL(RealKind) :: alpha
  TYPE(VectorFace_T), INTENT(IN)  :: Vector1(:)
  TYPE(VectorFace_T), INTENT(INOUT)  :: Vector2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpy(alpha,Vector1(ibLoc),Vector2(ibLoc))
  END DO

END SUBROUTINE  Axpy_VecFace


SUBROUTINE  Ax1mx2py_VelFace(alpha,Velocity1,Velocity2,Velocity3)

  REAL(RealKind) :: alpha
  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:),Velocity2(:)
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity3(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Ax1mx2py(alpha,Velocity1(ibLoc),Velocity2(ibLoc),Velocity3(ibLoc))
  END DO

END SUBROUTINE  Ax1mx2py_VelFace


SUBROUTINE  Scale_VelFace(alpha,Velocity)

  REAL(RealKind) :: alpha
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,Velocity(ibLoc))
  END DO

END SUBROUTINE  Scale_VelFace

SUBROUTINE  Scale_VecFace(alpha,Vector)

  REAL(RealKind) :: alpha
  TYPE(VectorFace_T), INTENT(INOUT)  :: Vector(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,Vector(ibLoc))
  END DO

END SUBROUTINE  Scale_VecFace


SUBROUTINE  Xpay_VelFace(Velocity1,alpha,Velocity2)

  REAL(RealKind) :: alpha
  TYPE(VelocityFace_T), INTENT(IN)  :: Velocity1(:)
  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity2(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Xpay(Velocity1(ibLoc),alpha,Velocity2(ibLoc))
  END DO

END SUBROUTINE  Xpay_VelFace

SUBROUTINE  Value_VelFace(Velocity,Value)

  TYPE(VelocityFace_T), INTENT(INOUT)  :: Velocity(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Velocity(ibLoc)=Value
  END DO

END SUBROUTINE  Value_VelFace

SUBROUTINE  Value_VecFace(Vector,Value)

  TYPE(VectorFace_T), INTENT(INOUT)  :: Vector(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Vector(ibLoc)=Value
  END DO

END SUBROUTINE  Value_VecFace

SUBROUTINE RandCopyFace(Velocity)

!---------------------------------------------------------------
!---  Copy Face Values from Neighboring Blocks
!---------------------------------------------------------------

  TYPE(VelocityFace_T), TARGET, INTENT(INOUT)  :: Velocity(:)

  INTEGER :: in
  INTEGER :: jx,jy,jz
  TYPE(VelocityFace_T), POINTER  :: subblock
  TYPE (Nachbar_T), POINTER :: Nachbar

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    subblock=>Velocity(ibLoc)
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      IF (Nachbar%nType=='iw'.AND.ibLoc<=ibnLoc) THEN           ! west neighbor
        CALL CopyFace(subblock%uF(ix0,jy0+1:jy1,jz0+1:jz1),FU(ix0,jy0+1:jy1,jz0+1:jz1), &
                      Velocity(ibnLoc)%uF(jNx1,jNy0+1:jNy1,jNz0+1:jNz1),Floor(ibn)%WeiFU(jNx1,jNy0+1:jNy1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='ie'.AND.ibLoc<ibnLoc) THEN       ! east neighbor
        CALL CopyFace(subblock%uF(ix1,jy0+1:jy1,jz0+1:jz1),FU(ix1,jy0+1:jy1,jz0+1:jz1), &
                      Velocity(ibnLoc)%uF(jNx0,jNy0+1:jNy1,jNz0+1:jNz1),Floor(ibn)%WeiFU(jNx0,jNy0+1:jNy1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='is'.AND.ibLoc<=ibnLoc) THEN      ! south neighbor
        CALL CopyFace(subblock%vF(jx0+1:jx1,iy0,jz0+1:jz1),FV(jx0+1:jx1,iy0,jz0+1:jz1),  &
                      Velocity(ibnLoc)%vF(jNx0+1:jNx1,jNy1,jNz0+1:jNz1),Floor(ibn)%WeiFV(jNx0+1:jNx1,jNy1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='in'.AND.ibLoc<ibnLoc) THEN       ! north neighbor
        CALL CopyFace(subblock%vF(jx0+1:jx1,iy1,jz0+1:jz1),FV(jx0+1:jx1,iy1,jz0+1:jz1),  &
                      Velocity(ibnLoc)%vF(jNx0+1:jNx1,jNy0,jNz0+1:jNz1),Floor(ibn)%WeiFV(jNx0+1:jNx1,jNy0,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='ib'.AND.ibLoc<=ibnLoc) THEN      ! bottom neighbor
        CALL CopyFace(subblock%wF(jx0+1:jx1,jy0+1:jy1,iz0),FW(jx0+1:jx1,jy0+1:jy1,iz0), &
                      Velocity(ibnLoc)%wF(jNx0+1:jNx1,jNy0+1:jNy1,jNz1),Floor(ibn)%WeiFW(jNx0+1:jNx1,jNy0+1:jNy1,jNz1))
      ELSE IF (Nachbar%nType=='it'.AND.ibLoc<ibnLoc) THEN       ! top neighbor
        CALL CopyFace(subblock%wF(jx0+1:jx1,jy0+1:jy1,iz1),FW(jx0+1:jx1,jy0+1:jy1,iz1), &
                      Velocity(ibnLoc)%wF(jNx0+1:jNx1,jNy0+1:jNy1,jNz0),Floor(ibn)%WeiFW(jNx0+1:jNx1,jNy0+1:jNy1,jNz0))
      END IF
    END DO ! in
  END DO   ! ibLoc
END SUBROUTINE RandCopyFace


SUBROUTINE ExchangeFace(Velocity)
   
  TYPE (VelocityFace_T) :: Velocity(:) 

   INTEGER :: i,ip,j   
   INTEGER :: tag
   TYPE(MPI_Request) :: req(2*NumberNeiProc)
   TYPE(MPI_Status) :: StatusArray(2*NumberNeiProc)

   REAL(RealKind) :: PutBuf(MaxPutBuf,NumberNeiProc)
   REAL(RealKind) :: GetBuf(MaxGetBuf,NumberNeiProc)
   
   CALL RandCopyFace(Velocity)
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         CALL PackBuff(chain(ip)%glied,PutBuf(:,i),Velocity)
      END IF
   END DO
                         
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         tag=MyId+NumProcs*ip
         CALL MPI_IRECV(GetBuf(1,i),chain(ip)%LenUnPack, &
                       MPI_RealKind,ip,tag,      &
                       MPI_COMM_WORLD,req(i),MPIErr)
      END IF
   END DO
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         tag=ip+NumProcs*MyId
         CALL MPI_ISEND(PutBuf(1,i),chain(ip)%LenPack, &
                       MPI_RealKind,ip,tag,      &
                       MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
      END IF
   END DO

   CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

   i=0
      
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         CALL UnPackBuff(chain(ip)%glied,GetBuf(:,i),Velocity)
      END IF
   END DO

CONTAINS 

SUBROUTINE PackBuff(chain,Buffer,Velocity)

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (VelocityFace_T) :: Velocity(:)

  INTEGER :: iBuf
  INTEGER :: i, len

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    IF (nType=='pw') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%uF(jx0,jy0+1:jy1,jz0+1:jz1),FU(jx0,jy0+1:jy1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pe') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%uF(jx1,jy0+1:jy1,jz0+1:jz1),FU(jx1,jy0+1:jy1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='ps') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%vF(jx0+1:jx1,jy0,jz0+1:jz1),FV(jx0+1:jx1,jy0,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pn') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%vF(jx0+1:jx1,jy1,jz0+1:jz1),FV(jx0+1:jx1,jy1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pb') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%wF(jx0+1:jx1,jy0+1:jy1,jz0),FW(jx0+1:jx1,jy0+1:jy1,jz0),CopyCase)
    ELSE IF (nType=='pt') THEN
      CALL PackFace(Buffer,iBuf,Velocity(ibCLoc)%wF(jx0+1:jx1,jy0+1:jy1,jz1),FW(jx0+1:jx1,jy0+1:jy1,jz1),CopyCase)
    END IF
  END DO
END SUBROUTINE PackBuff

SUBROUTINE UnPackBuff(chain,Buffer,Velocity)

  IMPLICIT NONE

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (VelocityFace_T) :: Velocity(:)

  INTEGER :: iBuf
  INTEGER :: i,len

  iBuf=0
  len=SIZE(chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)

    IF (nType=='pw') THEN
      CALL UnpackFace(Velocity(ibCLoc)%uF(jx0,jy0+1:jy1,jz0+1:jz1),FU(jx0,jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pe') THEN
      CALL UnpackFace(Velocity(ibCLoc)%uF(jx1,jy0+1:jy1,jz0+1:jz1),FU(jx1,jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='ps') THEN
      CALL UnpackFace(Velocity(ibCLoc)%vF(jx0+1:jx1,jy0,jz0+1:jz1),FV(jx0+1:jx1,jy0,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pn') THEN
      CALL UnpackFace(Velocity(ibCLoc)%vF(jx0+1:jx1,jy1,jz0+1:jz1),FV(jx0+1:jx1,jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pb') THEN
      CALL UnpackFace(Velocity(ibCLoc)%wF(jx0+1:jx1,jy0+1:jy1,jz0),FW(jx0+1:jx1,jy0+1:jy1,jz0),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pt') THEN
      CALL UnpackFace(Velocity(ibCLoc)%wF(jx0+1:jx1,jy0+1:jy1,jz1),FW(jx0+1:jx1,jy0+1:jy1,jz1),Buffer,iBuf,CopyCase)
    END IF
  END DO
END SUBROUTINE UnPackBuff

END SUBROUTINE ExchangeFace 

END MODULE VelocityFacePar_Mod

