MODULE Vortex_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Parallel_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE ReadProfileRad_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: alpha=2.5d0
  REAL(RealKind) :: Lz=6.0d+3  
  REAL(RealKind) :: Xi0=2.34d-3
  REAL(RealKind) :: b=53.5d+3
  REAL(RealKind) :: f=10.0d-4  
  REAL(RealKind) :: A=0.5d0
  REAL(RealKind) :: rb=50.0d+3
  REAL(RealKind) :: zb=6.0d+3
  REAL(RealKind) :: sigmar=15.0d+3
  REAL(RealKind) :: sigmaz=3.0d+3
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: th0=285.0d0
  LOGICAL :: ProfIn=.FALSE.
  INTEGER :: nr=2400, nnz=88, deltar=250, deltaz=250, maxiter=50
  
  NAMELIST /Example/ alpha &
                    ,Lz &
                    ,Xi0 &
                    ,b &
                    ,f &
                    ,ProfIn &
                    ,A &
                    ,rb &
                    ,zb &
                    ,sigmar &
                    ,sigmaz &
                    ,nr    &
                    ,nnz
CONTAINS

SUBROUTINE BalancedThermodynamics

  INTEGER :: i, k, iter, nProf
  REAL(8) , ALLOCATABLE :: VV(:,:), P(:,:), rho(:,:), T(:,:)
  REAL(8) , ALLOCATABLE ::  PRight(:), rhoRight(:), TRight(:)
  REAL(8) , ALLOCATABLE ::  rP(:), zP(:)
  REAL(8) :: z, r, Vrz, Vr
  REAL(8) :: Phi,fCorLoc,S,ThLoc,pLoc,RhoLoc,VLoc
  REAL(8) :: zLoc, Temp, PotTemp
  REAL(8) , ALLOCATABLE :: PProf(:), rhoProf(:), TProf(:), zProf(:)
  CHARACTER(33) :: line
  REAL(8) :: Res

  ALLOCATE(rP(1:nr+1))
  ALLOCATE(zP(1:nnz+1))
  ALLOCATE(VV(1:nr+1,1:nnz+1))
  ALLOCATE(P(1:nr+1,1:nnz+1))
  ALLOCATE(rho(1:nr,1:nnz+1))
  ALLOCATE(T(1:nr,1:nnz+1))
  ALLOCATE(PRight(1:nnz+1))
  ALLOCATE(RhoRight(1:nnz+1))
  ALLOCATE(TRight(1:nnz+1))

  rP(1)=deltar/2.0d0
  DO i=2,nr+1
    rP(i)=rP(i-1)+deltar
  END DO

  zP(1)=deltaz/2.0d0
  DO i=2,nnz+1
    zP(i)=zP(i-1)+deltaz
  END DO

  OPEN(UNIT=10,FILE="profiljordan.csv")   
  DO i=1,1                        
    READ(10,*) line                     
  END DO                          
  READ(10,*) nProf
  ALLOCATE(zProf(1:nProf))
  ALLOCATE(PProf(1:nProf))
  ALLOCATE(rhoProf(1:nProf))
  ALLOCATE(TProf(1:nProf))
  
  
  DO i=nProf,1,-1                       
    READ(10,*) zProf(i), rhoProf(i), PProf(i)                    
  END DO                          
 
   S=N*N/Grav
   DO k=1,nnz+1
    zLoc=zP(k)
    zLoc=MIN(zLoc,zProf(nProf))
    IF (zLoc<=zProf(1)) THEN
      pRight(k)=pProf(1)
      rhoRight(k)=rhoProf(1)
    ELSE
      DO i=2,nProf
        IF (zLoc<=zProf(i)) THEN
          pRight(k)=Int(zLoc,zProf(i-1),pProf(i-1),zProf(i),pProf(i))
          rhoRight(k)=Int(zLoc,zProf(i-1),rhoProf(i-1),zProf(i),rhoProf(i))
          EXIT
        END IF
      END DO
    END IF 
    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*zP(k))))**(Cpd/Rd)
    ThLoc=th0*EXP(zP(k)*S)
    pRight(k)=pLoc
    rhoRight(k)=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    TRight(k)=ThLoc*(pLoc/p0)**kappa
   END DO

  fCorLoc=fCor(LBOUND(fCor,1))
  DO i=1,nr+1
     Vr=(Xi0*(b*b)*(1.0d0-EXP(-(rP(i)/b)*(rP(i)/b))))/(2.0d0*rP(i))
     DO k=1,nnz+1
        Vrz=Vr*EXP(-zP(k)**alpha/(alpha*Lz**alpha))
        VV(i,k)=(fCorLoc*Vrz+(Vrz*Vrz)/rP(i))
     END DO
  END DO


  DO k=1,nnz+1
    rho(:,k)=rhoRight(k)
    P(:,k)=PRight(k)
    T(:,k)=TRight(k)
  END DO 

  DO iter=1,1
    DO i=nr,1,-1
      DO k=nnz,1,-1
        VLoc=0.25d0*(VV(i,k)+VV(i+1,k)+VV(i,k+1)+VV(i+1,k+1))
!       Half*(P(i+1,k+1)-P(i,k+1)+P(i+1,k)-P(i,k))/deltar=-Half*(P(i,k+1)-P(i,k)+P(i+1,k+1)-P(i+1,k))/deltaz*VV/Grav
        P(i,k)=(-Half*(P(i+1,k+1)-P(i,k+1)+P(i+1,k))/deltar &
               -Half*(P(i,k+1)+P(i+1,k+1)-P(i+1,k))/deltaz*VLoc/Grav) &
               /(-Half/deltar-Half/deltaz*VLoc/Grav)
      END DO
    END DO
  END DO
  DO i=nr,1,-1
    DO k=nnz,1,-1
      Rho(i,k)=-Two*(p(i,k+1)-p(i,k))/deltaz/Grav-Rho(i,k+1)
      T(i,k)=p(i,k)/(Rd*Rho(i,k))
    END DO
  END DO

! Check


  OPEN(UNIT=20,FILE=Profile)
  OPEN(UNIT=30,FILE='PlotPotTemp')
  WRITE(20,*) 'RhoProf' 
  WRITE(20,*)  nnz,nr
  DO i=1,nnz
   WRITE(20,*) zP(i)
  END DO
  DO i=1,nr
   WRITE(20,*) rP(i)
  END DO
  DO i=1,nr
    DO k=1,nnz
      WRITE(20,*) p(i,k)/(Rd*T(i,k)) 
    END DO
  END DO

  WRITE(20,*) 'PreProf' 
  WRITE(20,*)  nnz,nr
  DO i=1,nnz
   WRITE(20,*) zP(i)
  END DO
  DO i=1,nr
   WRITE(20,*) rP(i)
  END DO
  DO i=1,nr
    DO k=1,nnz
      WRITE(20,*) P(i,k) 
    END DO
  END DO

  WRITE(20,*) 'ThProf' 
  WRITE(20,*)  nnz,nr
  DO i=1,nnz
   WRITE(20,*) zP(i)
  END DO
  DO i=1,nr
   WRITE(20,*) rP(i)
  END DO
  DO i=1,nr
    DO k=1,nnz
      PotTemp=T(i,k)*(p0/p(i,k))**kappa
      WRITE(20,*) PotTemp
      WRITE(30,*) rP(i),zP(k),PotTemp
    END DO
  END DO
  CLOSE(10)
  CLOSE(20)
  CLOSE(30)

  DEALLOCATE(rP)
  DEALLOCATE(zP)
  DEALLOCATE(VV)
  DEALLOCATE(P)
  DEALLOCATE(rho)
  DEALLOCATE(T)
  DEALLOCATE(PRight)
  DEALLOCATE(RhoRight)
  DEALLOCATE(TRight)
  DEALLOCATE(zProf)
  DEALLOCATE(PProf)
  DEALLOCATE(rhoProf)
  DEALLOCATE(TProf)


END SUBROUTINE BalancedThermodynamics

FUNCTION Int(zM,zL,cL,zR,cR)
  REAL(8) :: Int,zM,zL,zR,cL,cR
  Int=((zM-zL)*cR+(zR-zM)*cL)/(zR-zL)
END FUNCTION Int

FUNCTION PreLoc(x,y,z,zHeight,Time)

  REAL(RealKind) :: PreLoc

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:),rP(:),zzP(:)
  CHARACTER*10, SAVE :: PreType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: r

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfileRad(zzP,rP,cInt,PreType,'PreProf')
      Load=.FALSE.
    END IF
    r=SQRT(x*x+y*y)
    PreLoc=ProfileRad(r,z,cInt,rP,zzP)
  ELSE
   PreLoc=0.0d0
  END IF

END FUNCTION PreLoc                    
FUNCTION ThLoc(x,y,z,zHeight,Time)

  REAL(RealKind) :: ThLoc
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun
  REAL(RealKind), POINTER, SAVE :: cInt(:,:),rP(:),zzP(:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: r, fak

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfileRad(zzP,rP,cInt,thType,'ThProf')
      Load=.FALSE.
    END IF
    r=SQRT(x*x+y*y)
    ThLoc=ProfileRad(r,z,cInt,rP,zzP)
  ELSE
    ThLoc=0.0d0
  END IF

!  potential temperature perturbation:
  fak=(x/(SQRT(x*x+y*y)))*((4*x*x)/(x*x+y*y)-3)  ! COS(3 \lambda)
! fak=Two*x*y/(x*x+y*y)                          ! COS(2 \lambda)
  ThLoc=ThLoc+A*fak*exp(-(((r-rb)*(r-rb))/(sigmar*sigmar) &
    +((z-zb)*(z-zb)*(z-zb)*(z-zb))/(sigmaz*sigmaz*sigmaz*sigmaz)))

END FUNCTION ThLoc

END MODULE Vortex_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Vortex_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells



SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Vortex_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)
END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE Vortex_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
  IF (MyId==0) THEN
    CALL BalancedThermodynamics
  END IF
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  WRITE(*,*) 'MyID in InputExample',MyID,MPIerr
END SUBROUTINE InputExample


FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE ReadProfileRad_Mod
  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,pLoc,ThLoc1

  S=N*N/Grav
  ThLoc1=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc1)
  RhoProf=Zero


END FUNCTION RhoProf


FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE ReadProfileRad_Mod
  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*exp(z*S)
  ThProfFun=Zero

END FUNCTION ThProfFun


FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Vortex_Mod
  USE ReadProfileRad_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,Th

  pLoc=PreLoc(x,y,z,zHeight,Time)
  Th=ThLoc(x,y,z,zHeight,Time)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*Th)

END FUNCTION RhoFun


FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  USE ReadProfileRad_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThStart=ThLoc(x,y,z,zHeight,Time)
END FUNCTION ThStart


FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  qvStart=0.0d0

END FUNCTION QvStart


FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  qcStart=0.0d0

END FUNCTION QcStart


FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=0.0d0

END FUNCTION QiStart


FUNCTION UStart(x,y,z)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z

  REAL(RealKind) :: r,Vr,Vrz
  
  r=SQRT(x*x+y*y)
 
  Vr=(Xi0*(b*b)*(1.0d0-EXP(-(r/b)*(r/b))))/(2.0d0*r)
  Vrz=Vr*EXP(-z**alpha/(alpha*Lz**alpha))

  UStart=-Vrz*y/r

END FUNCTION UStart


FUNCTION VStart(x,y,z)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z

  REAL(RealKind) :: r,Vr,Vrz

  r=SQRT(x*x+y*y)
  
  Vr=(Xi0*(b*b)*(1.0d0-EXP(-(r/b)*(r/b))))/(2.0d0*r)
  Vrz=Vr*EXP(-z**alpha/(alpha*Lz**alpha))
    
  VStart=Vrz*x/r

END FUNCTION VStart


FUNCTION DStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=0.0d0
END FUNCTION DStart


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE


FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE


FUNCTION WStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=0.0d0
END FUNCTION TkeStart


FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=0.0d0
END FUNCTION DisStart


FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=0.0d0
END FUNCTION TkeHStart


FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=0.0d0
END FUNCTION TkeVStart


FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=0.0d0
END FUNCTION LenStart


FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart


FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart


FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Rho1,Rho2
  Rho1=RhoFun(x,y,z,zHeight,Time)
  Rho2=RhoProf(x,y,z,zHeight,Time)
! RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
  RhoStart=Rho1-Rho2
END FUNCTION RhoStart


FUNCTION TStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun =0.0d0
END FUNCTION HeightFun
