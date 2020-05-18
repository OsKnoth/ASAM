MODULE Operator_Mod

  USE Parameter_Mod
  IMPLICIT NONE

CONTAINS
FUNCTION GradCentrN(cL,cC,cR,FL,FR,VL,VCe,VR)
  REAL(RealKind) :: GradCentrN
  REAL(RealKind) :: cL,cC,cR,FL,FR,VL,VCe,VR
  GradCentrN=Two*(FR*FR*(cR-cC)/(VCe+VCe+Eps)+FL*FL*(cC-cL)/(VCe+VCe+Eps))&
            /(FR+FL+Eps)
END FUNCTION GradCentrN

FUNCTION Gradient(cL,cR,F,VL,VR)
  REAL(RealKind) :: Gradient
  REAL(RealKind) :: cL,cR,F,VL,VR
  Gradient=Two*F*(cR-cL)/(VL+VR+Eps)
END FUNCTION Gradient

FUNCTION GradCentrFull(cL,cC,cR,dzL,dzC,dzR)
  REAL(RealKind) :: GradCentrFull
  REAL(RealKind) :: cL,cC,cR,dzL,dzC,dzR
  GradCentrFull=(cR-cC)/(dzR+dzC+Eps)+(cC-cL)/(dzC+dzL+Eps) 
END FUNCTION GradCentrFull

FUNCTION GradCentr(cL,cC,cR,FL,FR,VL,VCe,VR)
  REAL(RealKind) :: GradCentr
  REAL(RealKind) :: cL,cC,cR,FL,FR,VL,VCe,VR
  GradCentr=Two*(FR*FR*(cR-cC)/(VR+VCe+Eps)+FL*FL*(cC-cL)/(VCe+VL+Eps))&
            /(FR+FL+Eps)
END FUNCTION GradCentr

FUNCTION GradCentrT(cL,cC,cR,DL,DC,DR,FL,FR,VL,VCe,VR)
  REAL(RealKind) :: GradCentrT,cL,cC,cR,DL,DC,DR,FL,FR,VL,VCe,VR
  GradCentrT=(FR*FR*(DR+DC)*(cR-cC)/(VR+VCe+Eps)+FL*FL*(DL+DC)*(cC-cL)/(VCe+VL+Eps))&
            /(FR+FL+Eps)/(DC+Eps)
END FUNCTION GradCentrT

FUNCTION IntCellToFace(cL,cR,VL,VR)
  REAL(RealKind) :: IntCellToFace,cL,cR,VL,VR
  IntCellToFace=(VL*cR+VR*cL)/(VL+VR+Eps)
END FUNCTION IntCellToFace

FUNCTION CellToFaceVol(cL,cR,VL,VR)
  REAL(RealKind) :: CellToFaceVol,cL,cR,VL,VR
  CellToFaceVol=(VL*cL+VR*cR)/(VL+VR+Eps)
END FUNCTION CellToFaceVol

FUNCTION DF(DL,VL,DR,VR)
  REAL(RealKind) :: DF,DL,VL,DR,VR
! DF=DL*DR*(VL+VR)/(DL*VL+DR*VR+Eps)! harmonisches Mittel
  DF=(DL*VR+DR*VL)/(VL+VR+Eps)       ! arithmetsiches Mittel
END FUNCTION DF

FUNCTION FaceToCell(cL,cR,FL,FR)
  REAL(RealKind) :: FaceToCell,cL,cR,FL,FR
  FaceToCell=(FL*cL+FR*cR)/(FL+FR+Eps)
END FUNCTION FaceToCell

SUBROUTINE GradxP(p,VolC,FU,FV,FW,Grad)
  
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(0:1,-1:1,-1:1), INTENT(IN) :: p,VolC 
  REAL(RealKind),DIMENSION(-1:0), INTENT(IN) :: FU,FV,FW
  REAL(RealKind), INTENT(OUT) :: Grad 

!      Zentrale Groessen der betrachteten Box
!       VCe  = VolC(0 , 0 , 0  )
!       VW   = VolC(-1 , 0 , 0  )
!       VE   = VolC( 1 , 0 , 0 )
!       VS   = VolC( 0 ,-1 , 0 )
!       VN   = VolC( 0 , 1 , 0 )
!       VB   = VolC( 0 , 0 ,-1 )
!       VT   = VolC( 0 , 0 , 1 )
!       FWe  = FU(-1 )
!       FE   = FU( 0 )
!       FS   = FV(-1 )
!       FN   = FV( 0 )
!       FB   = FW(-1 )
!       FT   = FW( 0 )
!       u    = Half*(uCL(0,0,0)+uCR(0,0,0))/(Rho(0,0,0)+Eps) 
!       v    = Half*(vCL(0,0,0)+vCR(0,0,0))/(Rho(0,0,0)+Eps) 
!       w    = Half*(wCL(0,0,0)+wCR(0,0,0))/(Rho(0,0,0)+Eps) 
!       dudx = GradCentr(Half*(uCL(-1,0,0)+uCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
!                       ,u &
!                       ,Half*(uCL(1,0,0)+uCR(1,0,0))/(Rho(1,0,0)+Eps)    &
!                       ,FWe,FE,VW,VCe,VE)
!       dudy = GradCentr(Half*(uCL(0,-1,0)+uCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
!                       ,u &
!                       ,Half*(uCL(0,1,0)+uCR(0,1,0))/(Rho(0,1,0)+Eps)    &
!                       ,FS,FN,VS,VCe,VN)
!       dudz = GradCentr(Half*(uCL(0,0,-1)+uCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
!                       ,u &
!                       ,Half*(uCL(0,0,1)+uCR(0,0,1))/(Rho(0,0,1)+Eps)    &
!                       ,FB,FT,VB,VCe,VT)
!       dvdx = GradCentr(Half*(vCL(-1,0,0)+vCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
!                       ,v &
!                       ,Half*(vCL(1,0,0)+vCR(1,0,0))/(Rho(1,0,0)+Eps)    &
!                       ,FWe,FE,VW,VCe,VE)
!       dvdy = GradCentr(Half*(vCL(0,-1,0)+vCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
!                       ,v &
!                       ,Half*(vCL(0,1,0)+vCR(0,1,0))/(Rho(0,1,0)+Eps)    &
!                       ,FS,FN,VS,VCe,VN)
!       dvdz = GradCentr(Half*(vCL(0,0,-1)+vCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
!                       ,v &
!                       ,Half*(vCL(0,0,1)+vCR(0,0,1))/(Rho(0,0,1)+Eps)    &
!                       ,FB,FT,VB,VCe,VT)
!       dwdx = GradCentr(Half*(wCL(-1,0,0)+wCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
!                       ,w &
!                       ,Half*(wCL(1,0,0)+wCR(1,0,0))/(Rho(1,0,0)+Eps)    &
!                       ,FWe,FE,VW,VCe,VE)
!       dwdy = GradCentr(Half*(wCL(0,-1,0)+wCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
!                       ,w &
!                       ,Half*(wCL(0,1,0)+wCR(0,1,0))/(Rho(0,1,0)+Eps)    &
!                       ,FS,FN,VS,VCe,VN)
!       dwdz = GradCentr(Half*(wCL(0,0,-1)+wCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
!                       ,w &
!                       ,Half*(wCL(0,0,1)+wCR(0,0,1))/(Rho(0,0,1)+Eps)    &
!                       ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
!       S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
!                    (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
!                    (dvdz+dwdy)*(dvdz+dwdy))
!       O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
!                (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE GradxP

SUBROUTINE DeformVort(uCL,vCL,wCL,uCR,vCR,wCR,Rho,VolC,FU,FV,FW,S,O)
  
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uCL,vCL,wCL,uCR,vCR,wCR
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: Rho,VolC 
  REAL(RealKind),DIMENSION(-1:0), INTENT(IN) :: FU,FV,FW
  REAL(RealKind), INTENT(OUT) :: S,O 

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU(-1 )
        FE   = FU( 0 )
        FS   = FV(-1 )
        FN   = FV( 0 )
        FB   = FW(-1 )
        FT   = FW( 0 )
        u    = Half*(uCL(0,0,0)+uCR(0,0,0))/(Rho(0,0,0)+Eps) 
        v    = Half*(vCL(0,0,0)+vCR(0,0,0))/(Rho(0,0,0)+Eps) 
        w    = Half*(wCL(0,0,0)+wCR(0,0,0))/(Rho(0,0,0)+Eps) 
        dudx = GradCentr(Half*(uCL(-1,0,0)+uCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
                        ,u &
                        ,Half*(uCL(1,0,0)+uCR(1,0,0))/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentr(Half*(uCL(0,-1,0)+uCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
                        ,u &
                        ,Half*(uCL(0,1,0)+uCR(0,1,0))/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentr(Half*(uCL(0,0,-1)+uCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
                        ,u &
                        ,Half*(uCL(0,0,1)+uCR(0,0,1))/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentr(Half*(vCL(-1,0,0)+vCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
                        ,v &
                        ,Half*(vCL(1,0,0)+vCR(1,0,0))/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentr(Half*(vCL(0,-1,0)+vCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
                        ,v &
                        ,Half*(vCL(0,1,0)+vCR(0,1,0))/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentr(Half*(vCL(0,0,-1)+vCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
                        ,v &
                        ,Half*(vCL(0,0,1)+vCR(0,0,1))/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentr(Half*(wCL(-1,0,0)+wCR(-1,0,0))/(Rho(-1,0,0)+Eps) &
                        ,w &
                        ,Half*(wCL(1,0,0)+wCR(1,0,0))/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentr(Half*(wCL(0,-1,0)+wCR(0,-1,0))/(Rho(0,-1,0)+Eps) &
                        ,w &
                        ,Half*(wCL(0,1,0)+wCR(0,1,0))/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentr(Half*(wCL(0,0,-1)+wCR(0,0,-1))/(Rho(0,0,-1)+Eps) &
                        ,w &
                        ,Half*(wCL(0,0,1)+wCR(0,0,1))/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
                     (dvdz+dwdy)*(dvdz+dwdy))
        O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
                 (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE DeformVort

SUBROUTINE DeformVortC(uC,vC,wC,Rho,VolC,FU1,FV1,FW1,S,O)
  
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uC,vC,wC
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: Rho,VolC 
  REAL(RealKind),DIMENSION(-1:0,0:0,0:0), INTENT(IN) :: FU1
  REAL(RealKind),DIMENSION(0:0,-1:0,0:0), INTENT(IN) :: FV1
  REAL(RealKind),DIMENSION(0:0,0:0,-1:0), INTENT(IN) :: FW1
  REAL(RealKind), INTENT(OUT) :: S,O 

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU1(-1,0,0 )
        FE   = FU1( 0,0,0 )
        FS   = FV1(0,-1,0 )
        FN   = FV1(0, 0,0 )
        FB   = FW1(0,0,-1 )
        FT   = FW1(0,0, 0 )
        u    = uC(0,0,0)/(Rho(0,0,0)+Eps) 
        v    = vC(0,0,0)/(Rho(0,0,0)+Eps) 
        w    = wC(0,0,0)/(Rho(0,0,0)+Eps) 
        dudx = GradCentr(uC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,u &
                        ,uC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentr(uC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,u &
                        ,uC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentr(uC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,u &
                        ,uC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentr(vC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,v &
                        ,vC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentr(vC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,v &
                        ,vC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentr(vC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,v &
                        ,vC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentr(wC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,w &
                        ,wC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentr(wC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,w &
                        ,wC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentr(wC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,w &
                        ,wC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
                     (dvdz+dwdy)*(dvdz+dwdy))
        O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
                 (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE DeformVortC

SUBROUTINE StrainRateRho(uC,vC,wC,Rho,VolC,FU1,FV1,FW1,S,Sij)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uC,vC,wC
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: Rho,VolC
  REAL(RealKind),DIMENSION(-1:0,0:0,0:0), INTENT(IN) :: FU1
  REAL(RealKind),DIMENSION(0:0,-1:0,0:0), INTENT(IN) :: FV1
  REAL(RealKind),DIMENSION(0:0,0:0,-1:0), INTENT(IN) :: FW1
  REAL(RealKind),DIMENSION(1:6)  :: Sij
  REAL(RealKind), INTENT(OUT) :: S

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU1(-1,0,0 )
        FE   = FU1( 0,0,0 )
        FS   = FV1(0,-1,0 )
        FN   = FV1(0, 0,0 )
        FB   = FW1(0,0,-1 )
        FT   = FW1(0,0, 0 )
        u    = uC(0,0,0)/(Rho(0,0,0)+Eps)
        v    = vC(0,0,0)/(Rho(0,0,0)+Eps)
        w    = wC(0,0,0)/(Rho(0,0,0)+Eps)
        dudx = GradCentr(uC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,u &
                        ,uC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentr(uC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,u &
                        ,uC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentr(uC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,u &
                        ,uC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentr(vC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,v &
                        ,vC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentr(vC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,v &
                        ,vC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentr(vC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,v &
                        ,vC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentr(wC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,w &
                        ,wC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentr(wC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,w &
                        ,wC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentr(wC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,w &
                        ,wC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
!        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
!                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
!                     (dvdz+dwdy)*(dvdz+dwdy))
        Sij(1) = Half*(dudx+dudx)
        Sij(2) = Half*(dvdy+dvdy)
        Sij(3) = Half*(dwdz+dwdz)
        Sij(4) = Half*(dudy+dvdx)
        Sij(5) = Half*(dudz+dwdx)
        Sij(6) = Half*(dvdz+dwdy)
        S = SQRT(Two*((Sij(1)*Sij(1))+(Sij(2)*Sij(2))+(Sij(3)*Sij(3))+ &
                      Two*(Sij(4)*Sij(4))+Two*(Sij(5)*Sij(5))+Two*(Sij(6)*Sij(6))))
END SUBROUTINE StrainRateRho

SUBROUTINE StrainRate(uC,vC,wC,VolC,FU1,FV1,FW1,S,Sij)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uC,vC,wC
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: VolC
  REAL(RealKind),DIMENSION(-1:0,0:0,0:0), INTENT(IN) :: FU1
  REAL(RealKind),DIMENSION(0:0,-1:0,0:0), INTENT(IN) :: FV1
  REAL(RealKind),DIMENSION(0:0,0:0,-1:0), INTENT(IN) :: FW1
  REAL(RealKind),DIMENSION(1:6)  :: Sij
  REAL(RealKind), INTENT(OUT) :: S

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU1(-1,0,0 )
        FE   = FU1( 0,0,0 )
        FS   = FV1(0,-1,0 )
        FN   = FV1(0, 0,0 )
        FB   = FW1(0,0,-1 )
        FT   = FW1(0,0, 0 )
        u    = uC(0,0,0)
        v    = vC(0,0,0)
        w    = wC(0,0,0)
        dudx = GradCentr(uC(-1,0,0) &
                        ,u &
                        ,uC(1,0,0)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentr(uC(0,-1,0) &
                        ,u &
                        ,uC(0,1,0)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentr(uC(0,0,-1) &
                        ,u &
                        ,uC(0,0,1)    &
                        ,FB,FT,VB,VCe,VT)
      dvdx = GradCentr(vC(-1,0,0) &
                        ,v &
                        ,vC(1,0,0)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentr(vC(0,-1,0) &
                        ,v &
                        ,vC(0,1,0)   &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentr(vC(0,0,-1) &
                        ,v &
                        ,vC(0,0,1)   &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentr(wC(-1,0,0) &
                        ,w &
                        ,wC(1,0,0)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentr(wC(0,-1,0) &
                        ,w &
                        ,wC(0,1,0)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentr(wC(0,0,-1) &
                        ,w &
                        ,wC(0,0,1)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
!        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
!                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
!                     (dvdz+dwdy)*(dvdz+dwdy))
        Sij(1) = Half*(dudx+dudx)
        Sij(2) = Half*(dvdy+dvdy)
        Sij(3) = Half*(dwdz+dwdz)
        Sij(4) = Half*(dudy+dvdx)
        Sij(5) = Half*(dudz+dwdx)
        Sij(6) = Half*(dvdz+dwdy)
        S = SQRT(Two*((Sij(1)*Sij(1))+(Sij(2)*Sij(2))+(Sij(3)*Sij(3))+ &                                                                                                                                                  
                      Two*(Sij(4)*Sij(4))+Two*(Sij(5)*Sij(5))+Two*(Sij(6)*Sij(6))))                                                                                                                                       
END SUBROUTINE StrainRate

SUBROUTINE DeformVortN(uC,vC,wC,Rho,VolC,FU1,FV1,FW1,S,O,dudz)
  
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uC,vC,wC
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: Rho,VolC 
  REAL(RealKind),DIMENSION(-1:0,0:0,0:0), INTENT(IN) :: FU1
  REAL(RealKind),DIMENSION(0:0,-1:0,0:0), INTENT(IN) :: FV1
  REAL(RealKind),DIMENSION(0:0,0:0,-1:0), INTENT(IN) :: FW1
  REAL(RealKind), INTENT(OUT) :: S,O 

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU1(-1,0,0 )
        FE   = FU1( 0,0,0 )
        FS   = FV1(0,-1,0 )
        FN   = FV1(0, 0,0 )
        FB   = FW1(0,0,-1 )
        FT   = FW1(0,0, 0 )
        u    = uC(0,0,0)/(Rho(0,0,0)+Eps) 
        v    = vC(0,0,0)/(Rho(0,0,0)+Eps) 
        w    = wC(0,0,0)/(Rho(0,0,0)+Eps) 
        dudx = GradCentrN(uC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,u &
                        ,uC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentrN(uC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,u &
                        ,uC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentrN(uC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,u &
                        ,uC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentrN(vC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,v &
                        ,vC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentrN(vC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,v &
                        ,vC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentrN(vC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,v &
                        ,vC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentrN(wC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,w &
                        ,wC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentrN(wC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,w &
                        ,wC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentrN(wC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,w &
                        ,wC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
                     (dvdz+dwdy)*(dvdz+dwdy))
        O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
                 (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE DeformVortN


SUBROUTINE DeformVortCW(uC,vC,wC,Rho,VolC,FU1,FV1,FW1,S,O,dudz)
  
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: uC,vC,wC
  REAL(RealKind),DIMENSION(-1:1,-1:1,-1:1), INTENT(IN) :: Rho,VolC 
  REAL(RealKind),DIMENSION(-1:0,0:0,0:0), INTENT(IN) :: FU1
  REAL(RealKind),DIMENSION(0:0,-1:0,0:0), INTENT(IN) :: FV1
  REAL(RealKind),DIMENSION(0:0,0:0,-1:0), INTENT(IN) :: FW1
  REAL(RealKind), INTENT(OUT) :: S,O 

!      Zentrale Groessen der betrachteten Box
        VCe  = VolC(0 , 0 , 0  )
        VW   = VolC(-1 , 0 , 0  )
        VE   = VolC( 1 , 0 , 0 )
        VS   = VolC( 0 ,-1 , 0 )
        VN   = VolC( 0 , 1 , 0 )
        VB   = VolC( 0 , 0 ,-1 )
        VT   = VolC( 0 , 0 , 1 )
        FWe  = FU1(-1,0,0 )
        FE   = FU1( 0,0,0 )
        FS   = FV1(0,-1,0 )
        FN   = FV1(0, 0,0 )
        FB   = FW1(0,0,-1 )
        FT   = FW1(0,0, 0 )
        u    = uC(0,0,0)/(Rho(0,0,0)+Eps) 
        v    = vC(0,0,0)/(Rho(0,0,0)+Eps) 
        w    = wC(0,0,0)/(Rho(0,0,0)+Eps) 
        dudx = GradCentr(uC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,u &
                        ,uC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dudy = GradCentr(uC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,u &
                        ,uC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentr(uC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,u &
                        ,uC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentr(vC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,v &
                        ,vC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentr(vC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,v &
                        ,vC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentr(vC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,v &
                        ,vC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentr(wC(-1,0,0)/(Rho(-1,0,0)+Eps) &
                        ,w &
                        ,wC(1,0,0)/(Rho(1,0,0)+Eps)    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentr(wC(0,-1,0)/(Rho(0,-1,0)+Eps) &
                        ,w &
                        ,wC(0,1,0)/(Rho(0,1,0)+Eps)    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentr(wC(0,0,-1)/(Rho(0,0,-1)+Eps) &
                        ,w &
                        ,wC(0,0,1)/(Rho(0,0,1)+Eps)    &
                        ,FB,FT,VB,VCe,VT)

!       Deformations(S)- und Vorticity(O)-Invarianten:
        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
                     (dvdz+dwdy)*(dvdz+dwdy))
        O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
                 (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE DeformVortCW

FUNCTION uFace(uF,RhoL,RhoR,VolL,VolR)
  REAL(RealKind) :: uFace
  REAL(RealKind) :: uF,RhoL,RhoR,VolL,VolR

  uFace=uF

END FUNCTION uFace

FUNCTION PreFace(pL,pR,RhoL,RhoR,VolL,VolR)

  REAL(RealKind) :: PreFace
  REAL(RealKind) :: pL,pR,RhoL,RhoR,VolL,VolR

  PreFace=(pL/(RhoL+Eps)*RhoR+pR/(RhoR+Eps)*RhoL)/(RhoL+RhoR+Eps)

END FUNCTION PreFace  

FUNCTION WFace(WeiL,WeiR,RhoL,RhoR,VolL,VolR)

  REAL(RealKind) :: WFace
  REAL(RealKind) :: WeiL,WeiR,RhoL,RhoR,VolL,VolR

  WFace=(WeiL/(RhoL+Eps)*VolL+WeiR/(RhoR+Eps)*VolR)/(VolL+VolR+Eps)

END FUNCTION WFace  
END MODULE Operator_Mod


