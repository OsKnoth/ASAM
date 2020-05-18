MODULE Radiation_Mod
  USE Chemie_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Parameter_Mod
  USE Physics_Mod
  USE RadiationFuLiou_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: i,j,k
  INTEGER, PRIVATE :: ic,icE

  REAL(RealKind), PRIVATE, POINTER :: c(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qc(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoL(:,:,:,:)

 CONTAINS

  SUBROUTINE RadiationCompute(Vector,VelF,Rhs)
    TYPE(Vector4Cell_T) :: Vector
    TYPE(VelocityFace_T) :: VelF
    TYPE(Vector4Cell_T) :: Rhs
    INTEGER        :: ix,iy,iz,OutputUnit2
    REAL(RealKind) :: z,F0,F1,D,alphaZ,kappaC,zinv,Q0,Q1
    REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoInv
    REAL(RealKind) :: HeatingRate
    REAL(RealKind),POINTER :: MedianFlux(:)
    CHARACTER(80) :: ProfName

    Rho=>RhoCell(ibLoc)%c
    RhoV=>Vector%Vec(RhoVPos)%c
     IF (RhoCPos>0) THEN  ! Barthel
       RhoL=>Vector%Vec(RhoCPos)%c
     ELSE
       RhoL=>RhoLCell(ibLoc)%c
     END IF

    ALLOCATE(MedianFlux(iz0+1:iz1))

    MedianFlux=Zero
    F0=70
    F1=22
    D=3.75e-6       !  subsidence
    alphaZ=1
    kappaC=85        !  m**2/kg
    zinv=800        !  height of inversion
    RhoInv=1.13     !

    write(ProfName,'(I8)') ib
    OPEN(UNIT=OutputUnit,FILE='heatflux/'//TRIM(ADJUSTL(ProfName)),STATUS='UNKNOWN')
    DO iy=iy0,iy1
      DO ix=ix0,ix1
        Q1=0            !  cloudwater z->unlimited
        Q0=0            !  cloudwater 0->z
        DO iz=iz0+1,iz1
          Q1=Q1+kappaC*RhoL(ix,iy,iz,1)*dz(iz)
        END DO
        DO iz=iz0+1,iz1
          MedianFlux(iz)=MedianFlux(iz)+F0*exp(-Q1)+F1*exp(-Q0)+&
          &RhoInv*Cpd*D*alphaZ*((zP(iz)-zinv)**(4/3)/4+zinv*(zP(iz)-zinv)**(1/3))
          Q0=Q0+kappaC*RhoL(ix,iy,iz,1)*dz(iz)
          Q1=Q1-kappaC*RhoL(ix,iy,iz,1)*dz(iz)
        END DO
      END DO
    END DO
    DO iz=iz0+1,iz1
      write(OutputUnit,*) iz,zP(iz),MedianFlux(iz)/((iy1-iy0+1)*(ix1-ix0+1))
    END DO
    CLOSE(OutputUnit)

    DEALLOCATE(MedianFlux)
  END SUBROUTINE RadiationCompute

END MODULE Radiation_mod
