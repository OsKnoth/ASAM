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

  TYPE(JacSpMatrix4_T), POINTER :: JacTrans(:)

  TYPE(PressureVelocity), PRIVATE, SAVE, POINTER :: b(:),x(:)

  INTEGER, PRIVATE, SAVE :: Init=0

CONTAINS

SUBROUTINE ProjectVelFace(dt,Vel,IncrVecMet,VecMet,IncrVecChem,VecChem,Tol,VecG)

  REAL(RealKind) :: dt
  TYPE (VelocityFace_T) :: Vel(:)
  TYPE(Vector4Cell_T) :: IncrVecMet(:)
  TYPE(Vector4Cell_T) :: VecMet(:)
  TYPE(Vector4Cell_T), OPTIONAL :: IncrVecChem(:)
  TYPE(Vector4Cell_T), OPTIONAL :: VecChem(:)
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
  WRITE(*,*) 'Vel',DOT(Vel,Vel)
  CALL rhs(b,Vel,IncrVecMet,VecG)
  Temp=DOT2(b,b)
  WRITE(*,*) 'Temp',Temp
  MaxIter=QMRMaxIter
  IF (PRESENT(Tol)) THEN
    TolAct=Tol
  ELSE 
    TolAct=QMRTol
  END IF
  CALL SolveSound(x,b,MaxIter,TolAct)
  IF (PRESENT(VecG)) THEN
    IF (PRESENT(IncrVecChem)) THEN 
      CALL UpdateVelocity(x,b,Vel,IncrVecMet,VecMet,IncrVecChem,VecChem,VecG)
    ELSE
      CALL UpdateVelocity(x,b,Vel,IncrVecMet,VecMet,VecG)
    END IF  
  ELSE  
    IF (PRESENT(IncrVecChem)) THEN 
      CALL UpdateVelocity(x,b,Vel,IncrVecMet,VecMet,IncrVecChem,VecChem)
    ELSE
      CALL UpdateVelocity(x,b,Vel,IncrVecMet,VecMet)
    END IF  
  END IF

END SUBROUTINE ProjectVelFace

END MODULE IntSub_Mod
