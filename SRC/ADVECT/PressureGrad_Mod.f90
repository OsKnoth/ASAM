MODULE PressureGrad_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Physics_Mod

CONTAINS

SUBROUTINE PGradComputeC
  IF (GradFull) THEN
    CALL PGradComputeCFull
  ELSE
    CALL PGradComputeC1
  END IF
END SUBROUTINE PGradComputeC

SUBROUTINE PGradComputeC1

  INTEGER :: ix,iy,iz,in,jx,jy,jz
  REAL(RealKind) :: Grad
  REAL(RealKind) :: F,VolFine,VolCoarse,pFine,pCoarse

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Grad=FU(ix,iy,iz)*(p(ix+1,iy,iz,1)-p(ix,iy,iz,1))/(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)+Eps)
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Two*VolC(ix+1,iy,iz)*Grad
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Grad=FV(ix,iy,iz)*(p(ix,iy+1,iz,1)-p(ix,iy,iz,1))/(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)+Eps)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Two*VolC(ix,iy+1,iz)*Grad
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Grad=FW(ix,iy,iz)*(p(ix,iy,iz+1,1)-p(ix,iy,iz,1))/(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Two*VolC(ix,iy,iz+1)*Grad
        wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO


  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            uRhsL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)=uRhsL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1) &
                -Two*VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix0,jy,jz)*(p(ix0+1,jy,jz,1)-p(ix0,jy,jz,1))/(VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps) 
            uRhsL(ix0+1,jy,jz,1)=uRhsL(ix0+1,jy,jz,1)-Two*VolC(ix0+1,jy,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            uRhsR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)=uRhsR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1) &
                -Two*VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix1,jy,jz)*(p(ix1+1,jy,jz,1)-p(ix1,jy,jz,1))/(VolC(ix1,jy,jz)+VolC(ix1+1,jy,jz)+Eps) 
            uRhsR(ix1,jy,jz,1)=uRhsR(ix1,jy,jz,1)-Two*VolC(ix1,jy,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* & 
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            vRhsL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)=vRhsL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1) &
                -Two*VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy0,jz)*(p(jx,iy0+1,jz,1)-p(jx,iy0,jz,1))/(VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
            vRhsL(jx,iy0+1,jz,1)=vRhsL(jx,iy0+1,jz,1)-Two*VolC(jx,iy0+1,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            pFine= & 
               SUM(p(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            vRhsR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)=vRhsR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1) &
                -Two*VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy1,jz)*(p(jx,iy1+1,jz,1)-p(jx,iy1,jz,1))/(VolC(jx,iy1,jz)+VolC(jx,iy1+1,jz)+Eps)
            vRhsR(jx,iy1,jz,1)=vRhsR(jx,iy1,jz,1)-Two*VolC(jx,iy1,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (Volcoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            wRhsL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)=wRhsL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1) &
               -Two*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz0)*(p(jx,jy,iz0+1,1)-p(jx,jy,iz0,1))/(VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
            wRhsL(jx,jy,iz0+1,1)=wRhsL(jx,jy,iz0+1,1)-Two*VolC(jx,jy,iz0+1)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            wRhsR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)=wRhsR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1) &
               -Two*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz1)*(p(jx,jy,iz1+1,1)-p(jx,jy,iz1,1))/(VolC(jx,jy,iz1)+VolC(jx,jy,iz1+1)+Eps)
            wRhsR(jx,jy,iz1,1)=wRhsR(jx,jy,iz1,1)-Two*VolC(jx,jy,iz1)*Grad
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE PGradComputeC1

SUBROUTINE PGradComputeCFull

  INTEGER :: ix,iy,iz,in,jx,jy,jz
  REAL(RealKind) :: Grad
  REAL(RealKind) :: F,VolFine,VolCoarse,pFine,pCoarse
  REAL(RealKind) :: MetrLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Grad=DUU(ix,iy,iz)*FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)*(p(ix+1,iy,iz,1)-p(ix,iy,iz,1))/(MetrXY(iy)*(dx(ix+1)+dx(ix)+Eps)) 
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Two*VolC(ix+1,iy,iz)*Grad
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Grad=DUV(ix,iy,iz)*FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)*(p(ix,iy+1,iz,1)-p(ix,iy,iz,1))/(MetrYX(ix)*(dy(iy+1)+dy(iy)+Eps)) 
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Two*VolC(ix,iy+1,iz)*Grad
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Grad=DUW(ix,iy,iz)*FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)*(p(ix,iy,iz+1,1)-p(ix,iy,iz,1))/(dz(iz+1)+dz(iz)+Eps) 
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Two*VolC(ix,iy,iz+1)*Grad
        wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)-Two*VolC(ix,iy,iz)*Grad
      END DO
    END DO
  END DO


  DO in=1,AnzahlNachbar
    CALL Set(Nachbars(in))
    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrXY(jy:jy+IncrY-1)*dy(jy:jy+IncrY-1)) &
                   /SUM(dy(jy:jy+IncrY-1))
            pFine= &
               SUM(p(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(MetrLoc*(dLoc+dx(ix0+1)+Eps))
            uRhsL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)=uRhsL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1) &
                -Two*VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)*(p(ix0+1,jy,jz,1)-p(ix0,jy,jz,1))/(MetrXY(jy)*(dx(ix0+1)+dLoc+Eps)) 
            uRhsL(ix0+1,jy,jz,1)=uRhsL(ix0+1,jy,jz,1)-Two*VolC(ix0+1,jy,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrXY(jy:jy+IncrY-1)*dy(jy:jy+IncrY-1)) &
                   /SUM(dy(jy:jy+IncrY-1))
            pFine= &
               SUM(p(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(MetrLoc*(dx(ix1)+dLoc+Eps))
            uRhsR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)=uRhsR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1) &
                -Two*VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)*(p(ix1+1,jy,jz,1)-p(ix1,jy,jz,1))/(MetrXY(jy)*(dx(ix1)+dLoc+Eps)) 
            uRhsR(ix1,jy,jz,1)=uRhsR(ix1,jy,jz,1)-Two*VolC(ix1,jy,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrYX(jx:jx+IncrX-1)*dx(jx:jx+IncrX-1)) &
                   /SUM(dx(jx:jx+IncrX-1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* & 
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(MetrLoc*(dLoc+dy(iy0+1)+Eps))
            vRhsL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)=vRhsL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1) &
                -Two*VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)*(p(jx,iy0+1,jz,1)-p(jx,iy0,jz,1))/(MetrYX(jx)*(dy(iy0+1)+dLoc+Eps))
            vRhsL(jx,iy0+1,jz,1)=vRhsL(jx,iy0+1,jz,1)-Two*VolC(jx,iy0+1,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrYX(jx:jx+IncrX-1)*dx(jx:jx+IncrX-1)) &
                   /SUM(dx(jx:jx+IncrX-1))
            pFine= & 
               SUM(p(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(MetrLoc*(dy(iy1)+dLoc+Eps))
            vRhsR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)=vRhsR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1) &
                -Two*VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)*(p(jx,iy1+1,jz,1)-p(jx,iy1,jz,1))/(MetrYX(jx)*(dLoc+dy(iy1)+Eps))
            vRhsR(jx,iy1,jz,1)=vRhsR(jx,iy1,jz,1)-Two*VolC(jx,iy1,jz)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (Volcoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(dz(iz0+1)+dLoc+Eps)
            wRhsL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)=wRhsL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1) &
               -Two*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)*(p(jx,jy,iz0+1,1)-p(jx,jy,iz0,1))/(dz(iz0+1)+dLoc+Eps)
            wRhsL(jx,jy,iz0+1,1)=wRhsL(jx,jy,iz0+1,1)-Two*VolC(jx,jy,iz0+1)*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(dLoc+dz(iz1)+Eps)
            wRhsR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)=wRhsR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1) &
               -Two*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)*(p(jx,jy,iz1+1,1)-p(jx,jy,iz1,1))/(dz(iz1)+dLoc+Eps)
            wRhsR(jx,jy,iz1,1)=wRhsR(jx,jy,iz1,1)-Two*VolC(jx,jy,iz1)*Grad
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE PGradComputeCFull

SUBROUTINE PGradComputeFLin

  INTEGER :: ix,iy,iz,in,jx,jy,jz
  REAL(RealKind) :: Grad
  REAL(RealKind) :: F,VolFine,VolCoarse,pFine,pCoarse
  REAL(RealKind) :: VolLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VolLoc=VolC(ix+1,iy,iz)+VolC(ix,iy,iz)
        Grad=pFU(ix,iy,iz)*FU(ix,iy,iz)*(th(ix+1,iy,iz,1)-th(ix,iy,iz,1))/(VolLoc+Eps) 
        uFRhs(ix,iy,iz)=uFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VolLoc=VolC(ix,iy+1,iz)+VolC(ix,iy,iz)
        Grad=pFV(ix,iy,iz)*FV(ix,iy,iz)*(th(ix,iy+1,iz,1)-th(ix,iy,iz,1))/(VolLoc+Eps)
        vFRhs(ix,iy,iz)=vFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VolLoc=VolC(ix,iy,iz+1)+VolC(ix,iy,iz)
        Grad=pFW(ix,iy,iz)*FW(ix,iy,iz)*(th(ix,iy,iz+1,1)-th(ix,iy,iz,1))/(VolLoc+Eps) 
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(pFU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(th(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(th(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=pFU(ix0,jy,jz)*FU(ix0,jy,jz)*(th(ix0+1,jy,jz,1)-th(ix0,jy,jz,1))/(VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps) 
            uFRhs(ix0,jy,jz)=uFRhs(ix0,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(pFU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(th(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(th(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=pFU(ix1,jy,jz)*FU(ix1,jy,jz)*(th(ix1+1,jy,jz,1)-th(ix1,jy,jz,1))/(VolC(ix1,jy,jz)+VolC(ix1+1,jy,jz)+Eps) 
            uFRhs(ix1,jy,jz)=uFRhs(ix1,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(pFV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            pFine= &
               SUM(th(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* & 
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(th(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=pFV(jx,iy0,jz)*FV(jx,iy0,jz)*(th(jx,iy0+1,jz,1)-th(jx,iy0,jz,1))/(VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
            vFRhs(jx,iy0,jz)=vFRhs(jx,iy0,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(pFV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            pFine= & 
               SUM(th(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(th(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=pFV(jx,iy1,jz)*FV(jx,iy1,jz)*(th(jx,iy1+1,jz,1)-th(jx,iy1,jz,1))/(VolC(jx,iy1,jz)+VolC(jx,iy1+1,jz)+Eps)
            vFRhs(jx,iy1,jz)=vFRhs(jx,iy1,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(pFW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            pFine= &
               SUM(th(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            pCoarse= &
               SUM(th(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (Volcoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=pFW(jx,jy,iz0)*FW(jx,jy,iz0)*(th(jx,jy,iz0+1,1)-th(jx,jy,iz0,1))/(VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
            wFRhs(jx,jy,iz0)=wFRhs(jx,jy,iz0)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(pFW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            pFine= &
               SUM(th(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            pCoarse= &
               SUM(th(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=pFW(jx,jy,iz1)*FW(jx,jy,iz1)*(th(jx,jy,iz1+1,1)-th(jx,jy,iz1,1))/(VolC(jx,jy,iz1)+VolC(jx,jy,iz1+1)+Eps)
            wFRhs(jx,jy,iz1)=wFRhs(jx,jy,iz1)-Two*Grad
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE PGradComputeFLin

SUBROUTINE PGradComputeF

  IF (GradFull) THEN
    CALL PGradComputeFFull
  ELSE
    CALL PGradComputeF1
  END IF

END SUBROUTINE PGradComputeF

SUBROUTINE PGradComputeFFull

  INTEGER :: ix,iy,iz,in,jx,jy,jz
  REAL(RealKind) :: Grad
  REAL(RealKind) :: F,VolFine,VolCoarse,pFine,pCoarse
  REAL(RealKind) :: MetrLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Grad=FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)*(p(ix+1,iy,iz,1)-p(ix,iy,iz,1))/(MetrXY(iy)*(dx(ix+1)+dx(ix)+Eps))
        uFRhs(ix,iy,iz)=uFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Grad=FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)*(p(ix,iy+1,iz,1)-p(ix,iy,iz,1))/(MetrYX(ix)*(dy(iy+1)+dy(iy)+Eps))
        vFRhs(ix,iy,iz)=vFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Grad=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)*(p(ix,iy,iz+1,1)-p(ix,iy,iz,1))/(dz(iz+1)+dz(iz)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO


  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrXY(jy:jy+IncrY-1)*dy(jy:jy+IncrY-1)) &
                   /SUM(dy(jy:jy+IncrY-1))
            pFine= &
               SUM(p(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(MetrLoc*(dLoc+dx(ix0+1)+Eps))
            uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)*(p(ix0+1,jy,jz,1)-p(ix0,jy,jz,1))/(MetrXY(jy)*(dx(ix0+1)+dLoc+Eps)) 
            uFRhs(ix0,jy,jz)=uFRhs(ix0,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrXY(jy:jy+IncrY-1)*dy(jy:jy+IncrY-1)) &
                   /SUM(dy(jy:jy+IncrY-1))
            pFine= &
               SUM(p(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(MetrLoc*(dx(ix1)+dLoc+Eps))
            uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)*(p(ix1+1,jy,jz,1)-p(ix1,jy,jz,1))/(MetrXY(jy)*(dx(ix1)+dLoc+Eps)) 
            uFRhs(ix1,jy,jz)=uFRhs(ix1,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrYX(jx:jx+IncrX-1)*dx(jx:jx+IncrX-1)) &
                   /SUM(dx(jx:jx+IncrX-1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* & 
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(MetrLoc*(dLoc+dy(iy0+1)+Eps))
            vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)*(p(jx,iy0+1,jz,1)-p(jx,iy0,jz,1))/(MetrYX(jx)*(dy(iy0+1)+dLoc+Eps))
            vFRhs(jx,iy0,jz)=vFRhs(jx,iy0,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            MetrLoc=SUM(MetrYX(jx:jx+IncrX-1)*dx(jx:jx+IncrX-1)) &
                   /SUM(dx(jx:jx+IncrX-1))
            pFine= & 
               SUM(p(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(MetrLoc*(dy(iy1)+dLoc+Eps))
            vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)*(p(jx,iy1+1,jz,1)-p(jx,iy1,jz,1))/(MetrYX(jx)*(dLoc+dy(iy1)+Eps))
            vFRhs(jx,iy1,jz)=vFRhs(jx,iy1,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (Volcoarse+Eps)
            Grad=F/(F+Eps)*(pFine-pCoarse)/(dz(iz0+1)+dLoc+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)*(p(jx,jy,iz0+1,1)-p(jx,jy,iz0,1))/(dz(iz0+1)+dLoc+Eps)
            wFRhs(jx,jy,iz0)=wFRhs(jx,jy,iz0)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            Grad=F/(F+Eps)*(pCoarse-pFine)/(dLoc+dz(iz1)+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            Grad=FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)*(p(jx,jy,iz1+1,1)-p(jx,jy,iz1,1))/(dz(iz1)+dLoc+Eps)
            wFRhs(jx,jy,iz1)=wFRhs(jx,jy,iz1)-Two*Grad
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE PGradComputeFFull

SUBROUTINE PGradComputeF1

  INTEGER :: ix,iy,iz,in,jx,jy,jz
  REAL(RealKind) :: Grad
  REAL(RealKind) :: F,VolFine,VolCoarse,pFine,pCoarse

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Grad=FU(ix,iy,iz)*(p(ix+1,iy,iz,1)-p(ix,iy,iz,1))/(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)+Eps) 
        uFRhs(ix,iy,iz)=uFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Grad=FV(ix,iy,iz)*(p(ix,iy+1,iz,1)-p(ix,iy,iz,1))/(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)+Eps)
        vFRhs(ix,iy,iz)=vFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Grad=FW(ix,iy,iz)*(p(ix,iy,iz+1,1)-p(ix,iy,iz,1))/(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps) 
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)-Two*Grad
      END DO
    END DO
  END DO
  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            Grad=FU(ix0,jy,jz)*(p(ix0+1,jy,jz,1)-p(ix0,jy,jz,1))/(VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps) 
            uFRhs(ix0,jy,jz)=uFRhs(ix0,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)=uFRhs(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            Grad=FU(ix1,jy,jz)*(p(ix1+1,jy,jz,1)-p(ix1,jy,jz,1))/(VolC(ix1,jy,jz)+VolC(ix1+1,jy,jz)+Eps) 
            uFRhs(ix1,jy,jz)=uFRhs(ix1,jy,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* & 
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy0,jz)*(p(jx,iy0+1,jz,1)-p(jx,iy0,jz,1))/(VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)
            vFRhs(jx,iy0,jz)=vFRhs(jx,iy0,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            pFine= & 
               SUM(p(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)=vFRhs(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1) &
                -Two*Grad
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            Grad=FV(jx,iy1,jz)*(p(jx,iy1+1,jz,1)-p(jx,iy1,jz,1))/(VolC(jx,iy1,jz)+VolC(jx,iy1+1,jz)+Eps)
            vFRhs(jx,iy1,jz)=vFRhs(jx,iy1,jz)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jx=jx0+1,jx1,IncrX
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (Volcoarse+Eps)
            Grad=F*(pFine-pCoarse)/(VolFine+VolCoarse+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jy=jy0+1,jy1
            Grad=FW(jx,jy,iz0)*(p(jx,jy,iz0+1,1)-p(jx,jy,iz0,1))/(VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)
            wFRhs(jx,jy,iz0)=wFRhs(jx,jy,iz0)-Two*Grad
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jx=jx0+1,jx1,IncrX
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            pFine= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            pCoarse= &
               SUM(p(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            Grad=F*(pCoarse-pFine)/(VolFine+VolCoarse+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1) &
               -Two*Grad
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jy=jy0+1,jy1
            Grad=FW(jx,jy,iz1)*(p(jx,jy,iz1+1,1)-p(jx,jy,iz1,1))/(VolC(jx,jy,iz1)+VolC(jx,jy,iz1+1)+Eps)
            wFRhs(jx,jy,iz1)=wFRhs(jx,jy,iz1)-Two*Grad
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE PGradComputeF1


END MODULE PressureGrad_Mod
