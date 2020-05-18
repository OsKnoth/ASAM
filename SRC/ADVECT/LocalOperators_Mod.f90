MODULE LocalOperators_Mod

 USE Kind_Mod

 IMPLICIT NONE 

CONTAINS 
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

SUBROUTINE DeformVortC(uC,vC,wC,Rho,VolC,FU1,FV1,FW1,S,O,uvw)
  
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
  REAL(RealKind), INTENT(OUT) :: S,O,uvw 

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
        uvw = SQRT(u**two+v**two+w**two)

!       Deformations(S)- und Vorticity(O)-Invarianten:
        S = SQRT(two*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+                     &
                     (dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+     &
                     (dvdz+dwdy)*(dvdz+dwdy))
        O = SQRT((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+         &
                 (dvdz-dwdy)*(dvdz-dwdy))

END SUBROUTINE DeformVortC


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

MODULE LocalOperators_Mod


