MODULE IntSub_Mod

  USE Parallel_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
! USE Iter_Mod,CG=>FGMRES
! USE Iter_Mod,CG=>GMRES
  USE Iter_Mod
  USE JacAccGrav_Mod
  USE Function_Mod
  USE Rhs_Mod
  USE Output_Mod

  IMPLICIT NONE

  TYPE(JacSpMatrix4_T), POINTER :: JacMet(:)
  TYPE(JacSpMatrix4_T), POINTER :: JacChem(:)

  TYPE(PressureVelocity), PRIVATE, SAVE, POINTER :: b(:),x(:)

  INTEGER, PRIVATE, SAVE :: Init=0

CONTAINS

SUBROUTINE ProjectVelFace(dt,Vel,VecC,VelC,Tol,VecG)

  REAL(RealKind) :: dt
  TYPE (VelocityFace_T) :: Vel(:)
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VelC(:)
  REAL(RealKind), OPTIONAL :: Tol
  TYPE(Vector4Cell_T), OPTIONAL :: VecG(:)

  INTEGER :: MaxIter
  REAL(RealKind) :: TolAct
  REAL(RealKind) :: Temp

  IF (Init==0) THEN
    Init=1
    CALL Allocate(x)
    CALL Allocate(b)
  END IF
!   Projection
  b=Zero
  x=Zero
  CALL rhs(b,Vel,VecC,VecG)
  Temp=DOT2(b,b)
  MaxIter=QMRMaxIter
  IF (PRESENT(Tol)) THEN
    TolAct=Tol
  ELSE 
    TolAct=QMRTol
  END IF
  CALL SolveSound(x,b,MaxIter,TolAct)
  IF (PRESENT(VecG)) THEN
    CALL UpdateVelocity(x,b,Vel,VecC,VelC,VecG)
  ELSE  
    CALL UpdateVelocity(x,b,Vel,VecC,VelC)
  END IF

END SUBROUTINE ProjectVelFace

END MODULE IntSub_Mod
