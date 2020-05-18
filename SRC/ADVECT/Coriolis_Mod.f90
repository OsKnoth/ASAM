MODULE Coriolis_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Physics_Mod
  USE Buoyancy_Mod

CONTAINS

SUBROUTINE CoriolisFreeCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        uRhs(ix,iy,iz,1)=uRhs(ix,iy,iz,1)+Temp*vC(ix,iy,iz,1) 
        vRhs(ix,iy,iz,1)=vRhs(ix,iy,iz,1)-Temp*uC(ix,iy,iz,1) 
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisFreeCompute

SUBROUTINE CoriolisFreeComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        TempV=Temp*(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        TempU=Temp*(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+TempV*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+TempV*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-TempU*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-TempU*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisFreeComputeLR

SUBROUTINE CoriolisCylFreeComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV,Temp1

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        Temp1=BousinesqF1(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))/Rho(ix,iy,iz,1)
        TempV=Temp1*Temp*(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+TempV*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+TempV*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        TempU=Temp1*Temp*(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-TempU*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-TempU*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisCylFreeComputeLR

SUBROUTINE CoriolisCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        uRhs(ix,iy,iz,1)=uRhs(ix,iy,iz,1) &
                        +Temp*(vC(ix,iy,iz,1)-vE(ix,iy,iz,1))
        vRhs(ix,iy,iz,1)=vRhs(ix,iy,iz,1) &
                        -Temp*(uC(ix,iy,iz,1)-uE(ix,iy,iz,1))
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisCompute

SUBROUTINE CoriolisComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+Temp*(TempV-Rho(ix,iy,iz,1)*vE(ix,iy,iz,1))
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+Temp*(TempV-Rho(ix,iy,iz,1)*vE(ix,iy,iz,1))
        TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-Temp*(TempU-Rho(ix,iy,iz,1)*uE(ix,iy,iz,1))
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Temp*(TempU-Rho(ix,iy,iz,1)*uE(ix,iy,iz,1))
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisComputeLR

SUBROUTINE CoriolisProfileComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+Temp*(TempV-Rho(ix,iy,iz,1)*vEProf(iz))
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+Temp*(TempV-Rho(ix,iy,iz,1)*vEProf(iz))
        TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-Temp*(TempU-Rho(ix,iy,iz,1)*uEProf(iz))
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Temp*(TempU-Rho(ix,iy,iz,1)*uEProf(iz))
      END DO
    END DO
  END DO

END SUBROUTINE CoriolisProfileComputeLR

SUBROUTINE JacCoriolisCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        AS(IndexMet(uPosJac,vPosJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosJac,vPosJac))%c(ix,iy,iz,1)+ &
          Temp
        AS(IndexMet(vPosJac,uPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,uPosJac))%c(ix,iy,iz,1)- &
          Temp
      END DO
    END DO
  END DO

END SUBROUTINE JacCoriolisCompute

SUBROUTINE JacCoriolisComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=fCor(ix,iy)
        AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE JacCoriolisComputeLR

SUBROUTINE JacCoriolisCylCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=Omega
        AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
          Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE JacCoriolisCylCompute

SUBROUTINE CurvatureCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: phi,Temp
  REAL(RealKind) :: phiL,phiR,TanPhi

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      phiL=yP(iy-1)
      phiR=yP(iy)
      phi=Half*(phiL+phiR)
      TanPhi=-Curv*(COS(phiR)-COS(phiL))/(SIN(phiR)-SIN(phiL))
      DO ix=ix0+1,ix1
        Temp=TAN(phi)*uC(ix,iy,iz,1)/RadEarth/(Rho(ix,iy,iz,1)+Eps)
        uRhs(ix,iy,iz,1)=uRhs(ix,iy,iz,1)+Temp*vC(ix,iy,iz,1)
        vRhs(ix,iy,iz,1)=vRhs(ix,iy,iz,1)-Temp*uC(ix,iy,iz,1)
      END DO
    END DO
  END DO

END SUBROUTINE CurvatureCompute

SUBROUTINE CurvatureComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: phi,Temp,TempU,TempV
  REAL(RealKind) :: TempE,TempUE,TempVE
  REAL(RealKind) :: phiL,phiR,TanPhi

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      phiL=yP(iy-1)
      phiR=yP(iy)
      phi=Half*(phiL+phiR)
      TanPhi=-Curv*(COS(phiR)-COS(phiL))/(SIN(phiR)-SIN(phiL))
      DO ix=ix0+1,ix1
        TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        Temp=TanPhi/RadEarth/(Rho(ix,iy,iz,1)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+Temp*uCL(ix,iy,iz,1)*TempV
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+Temp*uCR(ix,iy,iz,1)*TempV
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-Temp*TempU*TempU
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Temp*TempU*TempU
      END DO
    END DO
  END DO

END SUBROUTINE CurvatureComputeLR

SUBROUTINE JacCurvatureCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: phi,Temp,TempU,TempV
  REAL(RealKind) :: phiL,phiR,TanPhi

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      phi=yP(iy-1)+0.5e0*dy(iy)
      phiL=yP(iy-1)
      phiR=yP(iy)
      TanPhi=-Curv*(COS(phiR)-COS(phiL))/(SIN(phiR)-SIN(phiL))
      DO ix=ix0+1,ix1
        TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        Temp=TanPhi/RadEarth/(Rho(ix,iy,iz,1)+Eps)

        AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)+ &
          Temp*TempV
        AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)+ &
          Temp*TempV

        AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*uCL(ix,iy,iz,1)*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosLJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*uCL(ix,iy,iz,1)*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosLJac))%c(ix,iy,iz,1)+ &
          Temp*uCR(ix,iy,iz,1)*FV(ix,iy-1,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosRJac,vPosRJac))%c(ix,iy,iz,1)+ &
          Temp*uCR(ix,iy,iz,1)*FV(ix,iy,iz)/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
          
        AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
          Two*Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*TempU
        AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosLJac,uPosRJac))%c(ix,iy,iz,1)- &
          Two*Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*TempU
        AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosLJac))%c(ix,iy,iz,1)- &
          Two*Temp*FU(ix-1,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*TempU
        AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
          Two*Temp*FU(ix,iy,iz)/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*TempU
      END DO
    END DO
  END DO

END SUBROUTINE JacCurvatureCompute

SUBROUTINE CurvatureCylComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV
  REAL(RealKind) :: Rad

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Rad=Half*(yP(iy-1)+yP(iy))
      DO ix=ix0+1,ix1
        TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        Temp=TempU/Rad/(Rho(ix,iy,iz,1)+Eps)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+Temp*TempV
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+Temp*TempV
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-Temp*TempU
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Temp*TempU
      END DO
    END DO
  END DO

END SUBROUTINE CurvatureCylComputeLR

SUBROUTINE JacCurvatureCylCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp
  REAL(RealKind) :: Rad

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Rad=Half*(yP(iy-1)+yP(iy))
      DO ix=ix0+1,ix1
        Temp=uC(ix,iy,iz,1)/Rad/(Rho(ix,iy,iz,1)+Eps)
        AS(IndexMet(uPosJac,vPosJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosJac,vPosJac))%c(ix,iy,iz,1)+ &
          Temp
        AS(IndexMet(vPosJac,uPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,uPosJac))%c(ix,iy,iz,1)- &
              2.0d0*Temp
        Temp=vC(ix,iy,iz,1)/Rad/(Rho(ix,iy,iz,1)+Eps)
        AS(IndexMet(uPosJac,uPosJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosJac,uPosJac))%c(ix,iy,iz,1)+ &
          Temp
      END DO
    END DO
  END DO

END SUBROUTINE JacCurvatureCylCompute

SUBROUTINE CentrifugalComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV,RhoLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Temp=Half*fCor(ix,iy)
      DO ix=ix0+1,ix1
        RhoLoc=BousinesqF2(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
        TempU=RhoLoc*TemP*Temp*(xP(ix-1)+Half*dx(ix))
        TempV=RhoLoc*TemP*Temp*(yP(iy-1)+Half*dy(iy))
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+TempU
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+TempU
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)+TempV
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)+TempV
      END DO
    END DO
  END DO

END SUBROUTINE CentrifugalComputeLR

SUBROUTINE JacCentrifugalComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV,Temp1

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Temp=Half*fCor(ix,iy)
      DO ix=ix0+1,ix1
        Temp1=BousinesqJac(ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
        TempU=Temp1*TemP*Temp*(xP(ix-1)+Half*dx(ix))
        TempV=Temp1*TemP*Temp*(yP(iy-1)+Half*dy(iy))
        IF (RhoPosJac>0) THEN
          AS(IndexMet(uPosJac,RhoPosJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosJac,RhoPosJac))%c(ix,iy,iz,1)+ &
            TempU
          AS(IndexMet(vPosJac,RhoPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,RhoPosJac))%c(ix,iy,iz,1)+ &
            TempV
        ELSE IF (ThPosJac>0) THEN
          AS(IndexMet(uPosJac,ThPosJac))%c(ix,iy,iz,1)=AS(IndexMet(uPosJac,ThPosJac))%c(ix,iy,iz,1)+ &
            TempU
          AS(IndexMet(vPosJac,ThPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,ThPosJac))%c(ix,iy,iz,1)+ &
            TempV
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE JacCentrifugalComputeLR

SUBROUTINE CentrifugalCylComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV,RhoLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Temp=Half*fCor(ix,iy)
      DO ix=ix0+1,ix1
        RhoLoc=BousinesqF2(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
        TempV=RhoLoc*TemP*Temp*(yP(iy-1)+Half*dy(iy))
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)+TempV
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)+TempV
      END DO
    END DO
  END DO

END SUBROUTINE CentrifugalCylComputeLR

SUBROUTINE JacCentrifugalCylComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,TempU,TempV,Temp1

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=Half*fCor(ix,iy)
        Temp1=BousinesqJac2(ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
        TempV=Temp1*Temp*Temp*(yP(iy-1)+Half*dy(iy))
        IF (RhoPosJac>0) THEN
          AS(IndexMet(vPosJac,RhoPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,RhoPosJac))%c(ix,iy,iz,1)+ &
            TempV
        ELSE IF (ThPosJac>0) THEN
          AS(IndexMet(vPosJac,ThPosJac))%c(ix,iy,iz,1)=AS(IndexMet(vPosJac,ThPosJac))%c(ix,iy,iz,1)+ &
            TempV
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE JacCentrifugalCylComputeLR

END MODULE Coriolis_Mod
