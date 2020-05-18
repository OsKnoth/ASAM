MODULE TimeStep_Mod
  USE Parallel_Mod
  USE Domain_Mod
  USE Names_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Init_Mod

  IMPLICIT NONE

  REAL(RealKind) :: CFLNumber

CONTAINS 
SUBROUTINE TimeStepCFL(VelF,dtAct,nsCFL)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dtAct
  INTEGER :: nsCFL

  REAL(RealKind) :: dtCFL
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: uMax,px,pxMax,uMaxC
  REAL(RealKind) :: uSMax
  REAL(RealKind) :: uMaxG,pxMaxG
  REAL(RealKind) :: uSMaxG
  REAL(RealKind) :: uLoc,vLoc,wLoc
  REAL(RealKind) :: uLocC,vLocC,wLocC
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc
  REAL(RealKind) :: uSLoc,vSLoc,wSLoc
  REAL(RealKind) :: dtS
  REAL(RealKind) :: Mach,uC,sC
  REAL(RealKind), SAVE :: MachGlob=0.0d0
  INTEGER :: nsLoc
  INTEGER :: iStage
  REAL(RealKind) :: nxSFac
  REAL(RealKind) :: nySFac
  REAL(RealKind) :: nzSFac

  nxSFac=1.0d0
  IF (domain%nx<=1) nxSFac=0.0d0
  nySFac=1.0d0
  IF (domain%ny<=1) nySFac=0.0d0
  nzSFac=1.0d0
  IF (domain%nz<=1) nzSFac=0.0d0
  dtCFL=Zero
  uMax=Zero
  uMaxC=Zero
  uSMax=Zero
  Mach=0.0d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    pxMax=Zero
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          uLoc=Half*(FU(ix,iy,iz)*VelF(ibLoc)%uF(ix,iy,iz)+FU(ix-1,iy,iz)*VelF(ibLoc)%uF(ix-1,iy,iz)) &
              /(VolC(ix,iy,iz)*RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
          uLocC=Half*(FU(ix,iy,iz)*VelF(ibLoc)%uF(ix,iy,iz)+FU(ix-1,iy,iz)*VelF(ibLoc)%uF(ix-1,iy,iz)) &
               /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps) &
               /(Half*(FU(ix,iy,iz)+FU(ix-1,iy,iz))+Eps) 
          vLoc=Half*(FV(ix,iy,iz)*VelF(ibLoc)%vF(ix,iy,iz)+FV(ix,iy-1,iz)*VelF(ibLoc)%vF(ix,iy-1,iz)) &
              /(VolC(ix,iy,iz)*RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
          vLocC=Half*(FV(ix,iy,iz)*VelF(ibLoc)%vF(ix,iy,iz)+FV(ix,iy-1,iz)*VelF(ibLoc)%vF(ix,iy-1,iz)) &
               /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps) &
               /(Half*(FV(ix,iy,iz)+FV(ix,iy-1,iz))+Eps) 
          wLoc=Half*(FW(ix,iy,iz)*VelF(ibLoc)%wF(ix,iy,iz)+FW(ix,iy,iz-1)*VelF(ibLoc)%wF(ix,iy,iz-1)) &
              /(VolC(ix,iy,iz)*RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
          wLocC=Half*(FW(ix,iy,iz)*VelF(ibLoc)%wF(ix,iy,iz)+FW(ix,iy,iz-1)*VelF(ibLoc)%wF(ix,iy,iz-1)) &
               /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps) &
               /(Half*(FW(ix,iy,iz)+FW(ix,iy,iz-1))+Eps) 
          dxLoc=dx(ix)
          dyLoc=dy(iy)
          dzLoc=dz(iz)
          dxLoc=VolC(ix,iy,iz)/(Half*(FU(ix,iy,iz)+FU(ix-1,iy,iz))+Eps)
          dyLoc=VolC(ix,iy,iz)/(Half*(FV(ix,iy,iz)+FV(ix,iy-1,iz))+Eps)
          dzLoc=VolC(ix,iy,iz)/(Half*(FW(ix,iy,iz)+FW(ix,iy,iz-1))+Eps)
          dxLoc=VolC(ix,iy,iz)/(MAX(FU(ix,iy,iz),FU(ix-1,iy,iz))+Eps)
          dyLoc=VolC(ix,iy,iz)/(MAX(FV(ix,iy,iz),FV(ix,iy-1,iz))+Eps)
          dzLoc=VolC(ix,iy,iz)/(MAX(FW(ix,iy,iz),FW(ix,iy,iz-1))+Eps)
          uSLoc=cS0*SQRT(nxSFac/(dxLoc*dxLoc+Eps)+nySFac/(dyLoc*dyLoc+Eps)+nzSFac/(dzLoc*dzLoc+Eps))*&
          &VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)    
          uMax=MAX(uMax,ABS(uLoc)+ABS(vLoc)+ABS(wLoc))
          uMaxC=MAX(uMaxC,SQRT(uLocC**2+vLocC**2+wLocC**2))
          uC=SQRT(uLocC**2+vLocC**2+wLocC**2)
          sC=SQRT(1.402*PreCell(ibLoc)%c(ix,iy,iz,1)/RhoCell(ibLoc)%c(ix,iy,iz,1))
          Mach=MAX(Mach,uC/sC)
          uSMax=MAX(uSMax,uSLoc)
        END DO  
      END DO  
    END DO  
  END DO
  CALL MPI_Allreduce(uMax,uMaxG,1,MPI_RealKind, &
&                    MPI_MAX,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(uSMax,uSMaxG,1,MPI_RealKind, &
&                    MPI_MAX,MPI_Comm_World,MPIErr)
  dtCFL=MIN(CFLNumber/(uMaxG+Eps),dtMax)
  IF (MyID==0) WRITE(*,*) 'uSMaxG',uSMaxG
  dtS=0.9d0/uSMaxG
  MachGlob=MAX(Mach,MachGlob)
  IF (TimeStepControl) THEN
    dtAct=MIN(dtCFL,dtMax)
  END IF
  nsCFL=MAX(CEILING(dtAct/dtS),nsMin)
  IF (MyID==0) THEN
    WRITE(*,*) 'uMaxG     ',uMaxG
    WRITE(*,*) 'dtcFL     ',dtcFL
    WRITE(*,*) 'dtS       ',dtS
    WRITE(*,*) 'dtAct     ',dtAct
    WRITE(*,*) 'nsCFL     ',nsCFL
    WRITE(*,*) 'CFLNumber ',CFLNumber
  END IF  
END SUBROUTINE TimeStepCFL
END MODULE TimeStep_Mod
