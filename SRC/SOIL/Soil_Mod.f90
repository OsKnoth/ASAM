MODULE Soil_Mod
  USE Kind_Mod
  USE Parameter_Mod
  USE SoilData_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod

  IMPLICIT NONE
  REAL(RealKind), PRIVATE, POINTER :: Th(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThB(:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThBRhs(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Wg(:,:)
  REAL(RealKind), PRIVATE, POINTER :: WgRhs(:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoS(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoD(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoVRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoLRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoDRhs(:,:,:,:)
  TYPE(Vec4_T), PRIVATE, POINTER :: AS(:)
  INTEGER :: nzS  ! Number of possible soil levels 
  INTEGER :: nzW  ! Number of possible soil levels 
  REAL(RealKind), ALLOCATABLE :: dzSoil(:)
  REAL(RealKind), ALLOCATABLE :: zSoil(:)
  
  INTEGER, Pointer :: SoilType(:)
  INTEGER :: LandClass  ! Land Use Class
  INTEGER :: nrsoillay

  CHARACTER(80) :: OutputNameSoil
  CHARACTER(80) :: OutputNameCut


CONTAINS

SUBROUTINE InitSoil

  INTEGER :: i

  ALLOCATE(dzSoil(1:Domain%nrsoillayers))
  ALLOCATE(zSoil(0:Domain%nrsoillayers))

  zSoil(0:Domain%nrsoillayers)=Domain%zSDepth(0:Domain%nrsoillayers)
  DO i=1,Domain%nrsoillayers
    dzSoil(i)=zSoil(i)-zSoil(i-1)
  END DO
  nzS=Domain%nrsoillayers
  nzW=Domain%nrsoillayers

END SUBROUTINE InitSoil

SUBROUTINE Soil(Vector,Rhs,UVec,Time)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vector4Cell_T) :: Rhs
  TYPE(Vector4Cell_T), OPTIONAL :: UVec
  REAL(RealKind) :: Time

  IF (PRESENT(UVec)) THEN
    uCL=>UVec%Vec(uPosL)%c
    vCL=>UVec%Vec(vPosL)%c
    wCL=>UVec%Vec(wPosL)%c
    uCR=>UVec%Vec(uPosR)%c
    vCR=>UVec%Vec(vPosR)%c
    wCR=>UVec%Vec(wPosR)%c
  ELSE
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
  END IF  
  uRhsL=>Rhs%Vec(uPosL)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsR=>Rhs%Vec(wPosR)%c
  Rho=>RhoCell(ibLoc)%c
  RhoV=>Vector%Vec(RhoVPos)%c
  RhoR=>Vector%Vec(RhoRPos)%c
  IF(RhoCPos>0) THEN ! Barthel
    RhoC=>Vector%Vec(RhoCPos)%c
  ELSE
    RhoC=>RhoLCell(ibLoc)%c
  END IF
  RhoI=>Vector%Vec(RhoIPos)%c
  RhoS=>Vector%Vec(RhoSpos)%c

  nrsoillay=Domain%nrsoillayers
  RhoRhs=>Rhs%Vec(RhoPos)%c
  RhoVRhs=>Rhs%Vec(RhoVPos)%c
  RhoLRhs=>Rhs%Vec(RhoCpos)%c
  Th=>Vector%Vec(ThPos)%c
  ThB=>Vector%Vec(ThPos)%cB
  ThBRhs=>Rhs%Vec(ThPos)%cB
  ThRhs=>Rhs%Vec(ThPos)%c
  Wg=>Vector%Vec(RhoVPos)%cB
! write(*,*) 'Soil_Mod: Wg, ', Wg
  WgRhs=>Rhs%Vec(RhoVPos)%cB
  CALL SoilCompute(Time)
END SUBROUTINE Soil

SUBROUTINE SoilCompute(Time)

  REAL(RealKind) :: Time

  INTEGER :: i,iw
  INTEGER :: ix,iy,iz,shad
  INTEGER :: nzThB,nzWg
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: n1G,n2G,n3G
  REAL(RealKind) :: FL,dL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind) :: u_shear                ! Shear Velocity
  REAL(RealKind) :: F_SH,F_LH              ! Sensible / Latent Heat Flux
  REAL(RealKind) :: TermW,TermT
  REAL(RealKind) :: zH,zRauh,zRauhT
  REAL(RealKind) :: alb,ee
  REAL(RealKind) :: Qdir,Qdif,Qinf
  REAL(RealKind) :: Qlat,Qsens
  REAL(RealKind) :: AlbSoil,AlbVeg
  REAL(RealKind) :: ThAir,TAir
  REAL(RealKind) :: RhoVAir,RhoLAir,RhoIAir,RhoDAir,RhoAir
  REAL(RealKind) :: RhoV_Sat               ! Saturation Humidity at Ground 
  REAL(RealKind) :: RhoV_g                 ! Specific Humdity at Interception Reservoir
  REAL(RealKind) :: DragM,DragH,DragQ,ThetaS
  REAL(RealKind) :: Wgn                    ! Interception Water [m]
  REAL(RealKind) :: RhsWgn
  REAL(RealKind) :: WSroot                 ! Root Water Content
  REAL(RealKind) :: TS(nzS),WS(nzW)        ! Soil Temperature, Water Content [m**3/m**3]
  REAL(RealKind) :: RhsTS(nzS),RhsWS(nzW)
  REAL(RealKind) :: FSoil                  ! Surface Radiative Flux 
  REAL(RealKind) :: LAI                    ! Leaf Area Index
  REAL(RealKind) :: F_plant                ! Area Covered with Plants
  REAL(RealKind) :: F_i                    ! Area Covered with intercepted water
  REAL(RealKind) :: F_par                  ! Fraction of absorbed photosynthetically active radiation
  REAL(RealKind) :: Cpml,Cvml,Lv,Rm,p
  REAL(RealKind) :: hu                     ! Relative Humdity near Ground
  REAL(RealKind) :: Eva_Pot,Eva_Soil,Eva_Intercept 
  REAL(RealKind) :: Transpiration(nzW),SumTransp
  REAL(RealKind) :: I_perc                 ! Percolation Rate of intercepted Water into Ground
  REAL(RealKind) :: dpotdt,rrl,ploc,potm,rrv,ThLoc
  REAL(RealKind) :: TotalFlux,MoistFlux,SensFlux,Tsurface,dThdt
  REAL(RealKind) :: dLoc,Surface,VolLoc,VolDiff
  REAL(RealKind) :: r

  DO i=1,NumBoundCell
    Wgn        = 0.0d0
    RhsWS      = 0.0d0
    RhsTS      = 0.0d0
    RhsWgn     = 0.0d0
    WSroot     = 0.0d0
    SumTransp  = 0.0d0

    ix         = BoundCell(i)%ix
    iy         = BoundCell(i)%iy
    iz         = BoundCell(i)%iz
!    n1         = BoundCell(i)%n1
!    n2         = BoundCell(i)%n2
!    n3         = BoundCell(i)%n3
    n1         = BoundCell(i)%n1G
    n2         = BoundCell(i)%n2G
    n3         = BoundCell(i)%n3G
    n1G        = BoundCell(i)%n1G
    n2G        = BoundCell(i)%n2G
    n3G        = BoundCell(i)%n3G
    FL         = BoundCell(i)%FL+Eps
    dL         = BoundCell(i)%dL+Eps
    alb        = BoundCell(i)%alb
    ee         = BoundCell(i)%ee
    shad       = BoundCell(i)%shad
    DragH      = BoundCell(i)%DragH
    DragQ      = BoundCell(i)%DragQ
    DragM      = BoundCell(i)%DragM
    LandClass  = BoundCell(i)%LandClass
    zRauh      = BoundCell(i)%zRauh
    zRauhT     = BoundCell(i)%zRauhT
    RhoAir     = Rho(ix,iy,iz,1)
    RhoVAir    = RhoV(ix,iy,iz,1)
    RhoLAir    = RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
    RhoIAir    = RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
    RhoDAir    = RhoAir-RhoVAir-RhoLAir+Eps
    ThLoc      = Th(ix,iy,iz,1)/(RhoAir+Eps)
    ThetaS     = BoundCell(i)%ThetaS
    SoilType   =>BoundCell(i)%SoilType
    ThAir      = Th(ix,iy,iz,1)
    nzThB      = SIZE(ThB(i,:))
    TS(1:nzThB)= ThB(i,1:)
    nzWg       = SIZE(Wg(i,:))-1
    Wgn        = MIN(MAX(0.0d0,Wg(i,nzWg+1)),5.0d-4*(1.0d0+5.0d0*PlantCover(LandClass)))
    WS(1:nzWg) = MAX(Zero,MIN(Wg(i,1:nzWg),cporv(SoilType(1:nzWg))))
    IF  (nzS>=1) THEN
      TSurface = TS(1)
    ELSE
      TSurface = BoundCell(i)%TeS 
    END IF
    RhoV_Sat = SaturVapor(TSurface)/(Rv*TSurface)   
    p        = PressureTheta(RhoDAir,RhoVAir,RhoLAir,RhoIAir,ThAir)+Eps
    TAir     = AbsTemp(RhoDAir,RhoVAir,p)
    Cpml     = Cpd*RhoDAir+Cpv*RhoVAir+Cpl*RhoLAir+Eps
    rrv      = RhoVAir/RhoDAir
    rrl      = RhoLAir/RhoDAir
    Rm       = Rd*RhoDAir+Rv*RhoVAir+Eps
    Lv       = LatHeat(TAir)

!   Evaluate velocity tangential and normal   
    uLocL = uCL(ix,iy,iz,1)/(RhoAir+Eps)
    uLocR = uCR(ix,iy,iz,1)/(RhoAir+Eps)
    uLoc  = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(RhoAir+Eps)
    vLocL = vCL(ix,iy,iz,1)/(RhoAir+Eps)
    vLocR = vCR(ix,iy,iz,1)/(RhoAir+Eps)
    vLoc  = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(RhoAir+Eps)
    wLocL = wCL(ix,iy,iz,1)/(RhoAir+Eps)
    wLocR = wCR(ix,iy,iz,1)/(RhoAir+Eps)
    wLoc  = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(RhoAir+Eps)
    VN    = uLoc*n1+vLoc*n2+wLoc*n3
    V     = uLoc*uLoc+ &
            vLoc*vLoc+ &
            wLoc*wLoc
    VT    = SQRT(MAX(V-VN*VN,Zero))          

    WRITE(*,*) 'DragM ',DragM
!   Wind
    uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(uLocL-n1*VN)
    vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)- &
                 RhoAir*FL*DragM*VT*(vLocL-n2*VN)
    wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(wLocL-n3*VN)
    uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(uLocR-n1*VN)
    vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(vLocR-n2*VN)
    wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)- &
                   RhoAir*FL*DragM*VT*(wLocR-n3*VN)
    
!   Evaluate height z
    zH = dL
     
    IF (LandClass==9) THEN ! ocean
      Qsens    = -FL*(DragH*VT+1.d-6)*(TAir-TSurface)                 ! [K m3 s-1]
      Qlat     = -FL*(DragQ*VT+1.d-6)*(RhoVAir-RhoV_Sat)/(RhoAir+Eps) ! [  m3 s-1]
    ELSE  ! land
      IF (nzS>1) THEN
!       Soil water flux
        DO iw=1,nzWg-1
          IF (iw.LE.(nzWg-2)) THEN
            TermW       = FlussCond(WS(iw),WS(iw+1),dzSoil(iw),dzSoil(iw+1),SoilType(iw),SoilType(iw+1))
            RhsWS(iw)   = RhsWS(iw)-TermW
            RhsWS(iw+1) = RhsWS(iw+1)+TermW
          ELSE
            TermW       = KappaSoil(WS(iw),SoilType(iw)) ! only gravitational drainage below 5th layer 
            RhsWS(iw)   = RhsWS(iw)-TermW
            RhsWS(iw+1) = RhsWS(iw+1)+TermW
          END IF
        END DO
  
!       Soil temperature flux
        DO iw=1,nzThB-1
          TermT       = FlussDiffT(TS(iw),TS(iw+1),dzSoil(iw),dzSoil(iw+1),WS(iw),WS(iw+1),SoilType(iw),SoilType(iw+1))
          RhsTS(iw)   = RhsTS(iw)-TermT
          RhsTS(iw+1) = RhsTS(iw+1)+TermT
        END DO
!       Lower Boundary Condition
        RhsWS(nzWg) = RhsWS(nzWg)-KappaSoil(WS(nzWg),SoilType(nzWg))
!       Lower Boundary Condition
        TermT         = FlussDiffT(TS(nzThB),TS(nzThB),Half*dzSoil(nzThB),Half*dzSoil(nzThB), &
                                   WS(nzThB),WS(nzThB),SoilType(nzThB),SoilType(nzThB))
        RhsTS(nzThB)  = RhsTS(nzThB)-TermT
      ELSE ! nr_SoilLayers=1
!       Lower Boundary Condition
        RhsWS(nzWg) = RhsWS(nzWg)-KappaSoil(WS(nzWg),SoilType(nzWg))
!       Lower Boundary Condition
        TermT         = FlussDiffT(TS(nzThB),TSClim,Half*dzSoil(nzThB),Half*dzSoil(nzThB), &
                                   WS(nzThB),WS(nzThB),SoilType(nzThB),SoilType(nzThB))
        RhsTS(nzThB)  = RhsTS(nzThB)-TermT
      END IF

!     ==== Init Data =====    
      RhoDAir    = RhoAir-RhoVAir-RhoLAir+Eps
      Cpml       = Cpd*RhoDAir+Cpv*RhoVAir+Cpl*RhoLAir
      Cvml       = Cvd*RhoDAir+Cvv*RhoVAir+Cpv*RhoLAir
      Rm         = Rd*RhoDAir+Rv*RhoVAir
      p          = PressureTheta(RhoDAir,RhoVAir,RhoLAir,RhoIAir,ThAir)
      TAir       = AbsTemp(RhoDAir,RhoVAir,p)
      TSurface   = TS(1)
      RhoV_Sat   = SaturVapor(TS(1))/(Rv*TS(1))    ! New one
      RhoV_g     = Wgn*RhoV_Sat+(1.d0-Wgn)*RhoVAir
      LAI        = LeafAreaIndex(LandClass)
      F_plant    = PlantCover(LandClass)
      F_i        = (1.0d0-F_plant)*(1.0d0-EXP(-Wgn/1.d-3))
      F_par      = One-EXP(-0.6d0*LAI)
      u_shear    = MAX(0.001d0,((DragM*VT)**Two*(uLoc**Two+vLoc**Two))**(One/Four))
      AlbSoil    = (One-F_plant)*AlbedoSoil(WS(1),SoilType(1))
      AlbVeg     = F_plant*( &
                   (One-F_par)*AlbedoSoil(WS(1),SoilType(1)) &
                   +F_par*AlbedoVeg(F_plant,LandClass))
      BoundCell(i)%alb = AlbSoil+AlbVeg
      
  ! WRITE (*,*) 'After flux: i,LandClass,TS(1),TSurface', i, LandClass, TS(1),TSurface

!     incoming (direct) solar Radiation
      Qdir       = (One-(AlbSoil+AlbVeg))*BoundCell(i)%raddirekt*shad &
                   *MAX(Zero,(n1G*radn1+n2G*radn2+n3G*radn3))
!     diffusive Radiation
      Qdif       = (One-(AlbSoil+AlbVeg))*BoundCell(i)%raddiffus
      ! downward long wave radiation hitting the surface
      Qinf       = BoundCell(i)%radinfred

!     WRITE (*,*) 'Soil: AlbSoil,AlbVeg,1-Sum,raddir,raddif,radinf',AlbSoil,AlbVeg,(One-(AlbSoil+AlbVeg)),BoundCell(i)%raddirekt,BoundCell(i)%raddiffus,BoundCell(i)%radinfred
!     WRITE (*,*) 'Soil: Qdir,Qdif,Qinf',Qdir,Qdif,Qinf

      IF (cfcap(SoilType(1))>WS(1)) THEN
        hu=Half*(1.0d0-COS((WS(1)*Pi)/(1.6d0*cfcap(SoilType(1)))))
      ELSE
        hu=1.0d0
      END IF

      Eva_Pot       = -(RhoVAir-RhoV_Sat)*DragH*VT  
      Eva_Soil      = MAX(0.0d0,(One-F_plant-F_i)*(hu*RhoV_Sat-RhoVAir)*DragH*VT) 
      Eva_Intercept = F_i*DragH*VT*(RhoVAir-MIN(RhoV_Sat,RhoV_g)) 

!     WRITE(*,*) 'Eva_Pot, Eva_Soil, Eva_Intercept',Eva_Pot, Eva_Soil, Eva_Intercept,F_plant,F_i
!     WRITE(*,*) 'i,SoilType',i,SoilType

!     WRITE (*,*) 'LC,WS(1),cporv',LandClass,WS(1),cporv(SoilType(1))
      IF (WS(1)>=cporv(SoilType(1))) THEN
        I_perc = 0.0d0
      ELSE
        I_perc = (One-F_plant)*(Wgn+Eva_Intercept*2.0d0*dtMax/RhoW)*RhoW/tau_perc
      END IF

      DO iw=1,nzWg
        IF (WS(iw).LE.cpwp(SoilType(iw))) THEN
          WSroot = WSroot
        ELSE
          WSroot = WSroot+WS(iw)*MIN(dzSoil(iw), &
                                 MAX(0.d0,(z_root(LandClass)-zSoil(iw-1))))
        END IF
      END DO
      WSroot = WSroot/(z_root(LandClass)+Eps)

      DO iw=1,nzWg
        Transpiration(iw) = TranspirFunction(TAir,DragH*VT,u_shear,WSroot,    &
                                             Qdir+Qdif,LAI,Eva_Pot,SoilType(iw),   &  
                                             LandClass,WS(iw),iw,F_plant,F_i) 
        SumTransp         = SumTransp+Transpiration(iw)
      END DO
  
!     === OLD ===
!     Upper Boundary Conditions (define fluxes F_SH and F_LF positive when from soil -> air)
!     F_SH     = -DragH*VT*(TAir-TS(1))            ! [K m s-1]
!     F_LH     = +SumTransp+Eva_Soil-Eva_Intercept ! [kg m-2 s-1]
!     FSoil    = Qdir+Qdif+Qinf-ee*sigma*TS(1)**4.0d0-Cpd*F_SH*RhoAir-Lv*F_LH ! [W m-2]
!     RhsTS(1) = RhsTS(1)+FSoil
!     RhsWgn   = RhsWgn-(I_perc-Eva_Intercept)/RhoW   
!     RhsWS(1) = RhsWS(1)-(Eva_Soil-I_perc)/RhoW
!     RhsTS(1) = RhsTS(1)-(Eva_Soil-I_perc)/RhoW*CapWater*TS(1)

!     === NEW ===
!     Upper Boundary Conditions (define fluxes F_SH and F_LF positive when from soil -> air)
      F_SH      = -DragH*VT*(TAir-TS(1))            ! [K m s-1]
      F_LH      = +SumTransp+Eva_Soil-Eva_Intercept ! [kg m-2 s-1]
      Qsens     = F_SH*FL        ! [K m3 s-1]
      Qlat      = F_LH/RhoAir*FL ! [  m3 s-1]
      dLoc      = ABS(n1)*dx(ix)+ABS(n2)*dy(iy)+ABS(n3)*dz(iz)
      VolLoc    = dx(ix)*dy(iy)*dz(iz)
      VolLoc    = MIN(FL*dLoc,VolLoc)
      MoistFlux = Qlat *Lv *RhoAir/FL*MIN(VolC(ix,iy,iz)/VolLoc,One) ! [kg s-3 = W m-2]
      SensFlux  = Qsens*Cpd*RhoAir/FL*MIN(VolC(ix,iy,iz)/VolLoc,One) ! [kg s-3 = W m-2]
      FSoil     = Qdir+Qdif+Qinf-ee*sigma*TS(1)**4.0d0-SensFlux-MoistFlux ! [W m-2]
      RhsTS(1)  = RhsTS(1)+FSoil ! [W m-2]
!     RhsTS(1)  = RhsTS(1)-(Eva_Soil-I_perc)/RhoW*CapWater*TS(1)
      RhsWgn    = RhsWgn-(I_perc-Eva_Intercept)/RhoW   
!     RhsWS(1)  = RhsWS(1)-(Eva_Soil-I_perc)/RhoW
      RhsWS(1)  = RhsWS(1)-(Eva_Soil)/RhoW

!     WRITE(*,*) 'Upper BC: ',i,LandClass,RhsTS(1),RhsWS(1),Eva_Soil,I_perc,F_LH

      IF (TS(1)<=T_freeze) RhsWgn = Zero 

!     Other Layers
      DO iw=1,nzThB
        RhsTS(iw) = (RhsTS(iw)-(RhsWS(iw)+Transpiration(iw))*TS(iw)*CapWater/RhoW) &
                    /dzSoil(iw)/Capacity(WS(iw),SoilType(iw))   ! [K s-1]
      END DO
      DO iw=1,nzWg
        RhsWS(iw) = (RhsWS(iw)-Transpiration(iw)/RhoW)/dzSoil(iw)
        IF (TS(iw)<=T_freeze) RhsWS(iw) = Zero
      END DO
      IF (nzS==1) THEN
        RhsTS(1)=RhsTS(1)-(TS(1)-TSClim)*Kg
        RhsWS(1)=RhsWS(1)+(0.08d0-WS(1))/tau_wg
      END IF
      ThBRhs(i,1:nzThB)   = RhsTS
      WgRhs(i,nzWg+1)     = RhsWgn
      WgRhs(i,1:nzWg)     = RhsWS
  

    END IF ! landclass 0-9
!   Compute Total Flux for virt. pot. temperature equation
    TotalFlux= Qlat*ThLoc*(Rv/Rm-LOG((p/p0)**(Rm/Cpml))*(Rv/Rm-Cpv/Cpml)) & ! latent heat flux
              +Qsens*ThLoc/Tair  ! sensible heat flux 

!   Compute RHS
!   ThRhs(ix,iy,iz,1)=ThRhs(ix,iy,iz,1)+TotalFlux
    CALL DistrFlux(ThRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                       ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                       ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,TotalFlux)
!   RhoVRhs(ix,iy,iz,1)=RhoVRhs(ix,iy,iz,1)+Qlat               ! lat.HeatFlux 
    CALL DistrFlux(RhoVRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                       ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                       ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,Qlat)
!   RhoRhs(ix,iy,iz,1)=RhoRhs(ix,iy,iz,1)+Qlat                 ! lat.HeatFlux
    CALL DistrFlux(RhoRhs(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1) &
                       ,FU(ix-1:ix,iy,iz),FV(ix,iy-1:iy,iz),FW(ix,iy,iz-1:iz) &
                       ,VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),FL,dx(ix),dy(iy),dz(iz),n1,n2,n3,Qlat)

!   Diagnostics for output
    dLoc=ABS(n1)*dx(ix)+ABS(n2)*dy(iy)+ABS(n3)*dz(iz)
    VolLoc=dx(ix)*dy(iy)*dz(iz)
    VolLoc=MIN(FL*dLoc,VolLoc)
    MoistFlux = Qlat *Lv *RhoAir/FL*MIN(VolC(ix,iy,iz)/VolLoc,One) ! [kg s-3 = W m-2]
    SensFlux  = Qsens*Cpd*RhoAir/FL*MIN(VolC(ix,iy,iz)/VolLoc,One) ! [kg s-3 = W m-2]
    BoundCell(i)%FluxLat  = MoistFlux 
    BoundCell(i)%FluxSens = SensFlux 
    BoundCell(i)%TSoil = TS(1) 
    BoundCell(i)%QVSoil = WS(1) 
    SensFluxCell(ibLoc)%cB(i,1)=BoundCell(i)%FluxSens
    LatFluxCell(ibLoc)%cB(i,1)=BoundCell(i)%FluxLat
    LandClassCell(ibLoc)%cB(i,1)=BoundCell(i)%LandClass
    RoughnessLengthCell(ibLoc)%cB(i,1)=BoundCell(i)%zRauh
    AlbedoCell(ibLoc)%cB(i,1)=BoundCell(i)%alb
    BulkCoeffDragCell(ibLoc)%cB(i,1)=BoundCell(i)%DragM
    BulkCoeffHeatCell(ibLoc)%cB(i,1)=BoundCell(i)%DragH
    BulkCoeffMoistCell(ibLoc)%cB(i,1)=BoundCell(i)%DragQ
  ! WRITE (*,*) 'i,MoistFlux,Qlat,Lv',i,MoistFlux,Qlat,Lv,F_LH,+SumTransp,+Eva_Soil,-Eva_Intercept,RhoVAir,RhoV_Sat,DragH,VT
  END DO
! Test output
! IF (MOD(INT(Time),60)==0) THEN
!   OPEN(123,FILE='SoilData.txt',STATUS='UNKNOWN',ACTION='WRITE',POSITION='APPEND')
!   WRITE(123,*) Time,TS(1),TAir,BoundCell(1)%FluxSens,RhoV_Sat,RhoVAir,BoundCell(1)%FluxLat,BoundCell(1)%DragM,BoundCell(1)%DragH
!   CLOSE(123)
! END IF

END SUBROUTINE SoilCompute

SUBROUTINE JacSoil(Vector,Jac)
  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  uCL=>Vector%Vec(uPosL)%c
  vCL=>Vector%Vec(vPosL)%c
  wCL=>Vector%Vec(wPosL)%c
  uCR=>Vector%Vec(uPosR)%c
  vCR=>Vector%Vec(vPosR)%c
  wCR=>Vector%Vec(wPosR)%c
  Rho=>RhoCell(ibLoc)%c
  Th=>Vector%Vec(ThPos)%c
  AS=>Jac
  CALL JacSoilCompute
END SUBROUTINE JacSoil

SUBROUTINE JacSoilCompute
  INTEGER :: i
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN  
  REAL(RealKind) :: FL,DragH,DragM
  REAL(RealKind) :: RhoAir
  
  DO i=1,NumBoundCell
    ix     = BoundCell(i)%ix
    iy     = BoundCell(i)%iy
    iz     = BoundCell(i)%iz
    n1     = BoundCell(i)%n1G
    n2     = BoundCell(i)%n2G
    n3     = BoundCell(i)%n3G
    FL     = BoundCell(i)%FL+Eps
    DragH  = BoundCell(i)%DragH
    DragM  = BoundCell(i)%DragM

    RhoAir = Rho(ix,iy,iz,1)
!   Evaluate velocity tangential and normal  
    uLoc = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/RhoAir
    vLoc = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/RhoAir
    wLoc = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/RhoAir 
    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero))

    AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosLJac,uPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n1*n1)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosLJac,vPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n2*n2)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosLJac,wPosLJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n3*n3)/(VolC(ix,iy,iz)+Eps)

    AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(uPosRJac,uPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n1*n1)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(vPosRJac,vPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n2*n2)/(VolC(ix,iy,iz)+Eps)
    AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(wPosRJac,wPosRJac))%c(ix,iy,iz,1)- &
                                                    RhoAir*FL*DragM*VT* &
                                                    (One-n3*n3)/(VolC(ix,iy,iz)+Eps)
 
!   AS(IndexMet(ThPosJac,ThPosJac))%c(ix,iy,iz,1) = AS(IndexMet(ThPosJac,ThPosJac))%c(ix,iy,iz,1)- &
!                                                   FL*(DragH*VT+1.d-6)/(VolC(ix,iy,iz)+Eps)
  END DO
END SUBROUTINE JacSoilCompute

FUNCTION KappaSoil(Theta_Soil,SoilType)
  REAL(RealKind) :: Theta_Soil
  REAL(RealKind) :: KappaSoil,m,n,WSeffektiv
  INTEGER :: SoilType
  WSEffektiv= (MAX(Eps,(Theta_Soil-cres(SoilType))))/(cporv(SoilType)-cres(SoilType))
  n         = lamb(SoilType)+1
  m         = 1.0d0-1.0d0/n
  SELECT CASE(SoilParam)
    CASE('Campbell')
      KappaSoil = ckw0(SoilType)*(Theta_Soil/cporv(SoilType))**(2.d0*cbedi(SoilType)+3) 
    CASE('CoreyBrooks')
      KappaSoil = ckw0(SoilType)*WSEffektiv**(2.d0*cbedi(SoilType)+3) 
    CASE('VanGenuchten')
      KappaSoil = ckw0(SoilType)*WSEffektiv**(One/Two)*(One-(One-WSEffektiv**(One/m))**m)**Two  
    CASE DEFAULT
      KappaSoil = ckw0(SoilType)*WSEffektiv**(One/Two)*(One-(One-WSEffektiv**(One/m))**m)**Two  
  END SELECT
  KappaSoil = MAX(Eps,KappaSoil)
END FUNCTION KappaSoil

FUNCTION PsiSoil(Theta_Soil,SoilType)
  REAL(RealKind) :: Theta_Soil,WSeffektiv
  REAL(RealKind) :: PsiSoil
  REAL(RealKind) :: m,n
  INTEGER :: SoilType
  n         = lamb(SoilType)+1 
  m         = 1.0d0-1.0d0/n
  WSEffektiv= (MAX(Eps,(Theta_Soil-cres(SoilType))))/(cporv(SoilType)-cres(SoilType))
  SELECT CASE(SoilParam)
    CASE('Campbell')
      PsiSoil = smPot(SoilType)*(MAX(Eps,Theta_Soil)/cporv(SoilType))**(-cbedi(SoilType)) 
    CASE('CoreyBrooks')
      PsiSoil   = smPot(SoilType)*WSEffektiv**(-cbedi(SoilType)) 
    CASE('VanGenuchten')
      PsiSoil   = smPot(SoilType)*(WSeffektiv**(-One/m)-One)**(One/n) 
    CASE DEFAULT
      PsiSoil   = smPot(SoilType)*(WSeffektiv**(-One/m)-One)**(One/n) 
  END SELECT
END FUNCTION PsiSoil

FUNCTION FlussCond(Theta_SoilT,Theta_SoilB,dzT,dzB,SoilTypeT,SoilTypeB) ! Conductive Flux
  REAL(RealKind) :: FlussCond
  REAL(RealKind) :: Theta_SoilT,Theta_SoilB,dzT,dzB
  INTEGER :: SoilTypeT,SoilTypeB
  FlussCond = 2.0d0/(dzT/KappaSoil(Theta_SoilT,SoilTypeT) &
                     +dzB/KappaSoil(Theta_SoilB,SoilTypeB)) &
              *(PsiSoil(Theta_SoilT,SoilTypeT)-PsiSoil(Theta_SoilB,SoilTypeB))  &
              +KappaSoil(Theta_SoilT,SoilTypeT)
END FUNCTION FlussCond

FUNCTION FlussDiffT(TsT,TsB,dzT,dzB,Theta_SoilT,Theta_SoilB,SoilTypeT,SoilTypeB)
  REAL(RealKind) :: FlussDiffT
  REAL(RealKind) :: Theta_SoilT,Theta_SoilB
  REAL(RealKind) :: TsT,TsB,dzT,dzB
  INTEGER :: SoilTypeT,SoilTypeB
  FlussDiffT = 2.0d0/(dzT/ThermCond(Theta_SoilT,SoilTypeT)   &
                     +dzB/ThermCond(Theta_SoilB,SoilTypeB))*(TsT-TsB) &
               +2.0d0/(dzT/KappaSoil(Theta_SoilT,SoilTypeT)  & 
                      +dzB/KappaSoil(Theta_SoilB,SoilTypeB)) &
               *(PsiSoil(Theta_SoilT,SoilTypeT)-PsiSoil(Theta_SoilB,SoilTypeB)) &
               *capwater*(dzT+dzB)/(dzT/TsT+dzB/TsB) &
               +KappaSoil(Theta_SoilT,SoilTypeT)*capwater*TsT
END FUNCTION FlussDiffT

FUNCTION ThermCond(Theta_Soil,SoilType) ! Thermal Conductivity
  REAL(RealKind) :: Theta_Soil,Psi
  REAL(RealKind) :: ThermCond
  INTEGER  :: SoilType
  Psi       = PsiSoil(Theta_Soil,SoilType) 
  IF (SoilType>2) THEN
    ThermCond = MAX(418.0d0*EXP(-LOG10(100.0d0*ABS(Psi))-2.7d0),0.172d0)
  ELSE 
    ThermCond = 0.172d0
  END IF
END FUNCTION ThermCond

FUNCTION Capacity(Theta_Soil,SoilType)
  REAL(RealKind) :: Capacity
  REAL(RealKind) :: Theta_Soil
  INTEGER :: SoilType
! soil capacity = sum of capacity of dry soil, capacity of water and capacity of air (in pores)
  Capacity = (1.0d0-cporv(SoilType))*crhoc(SoilType)+Theta_Soil*capwater+(cporv(SoilType)-Theta_Soil)*capair
END FUNCTION Capacity

FUNCTION  TranspirFunction(TAir,DragCoeff,Ushear,WSroot,RadPAR,   &
                           LAIndex,Epotg,SoilType,LandType,WSoil, &
                           iS,PlantArea,InterceptArea)
  REAL(RealKind) :: TranspirFunction
  REAL(RealKind) :: AtmoResist,FolResist
  REAL(RealKind) :: ReductTransp,Resist_la
  REAL(RealKind) :: StomaResist
  REAL(RealKind) :: DragCoeff
  REAL(RealKind) :: RadPAR
  REAL(RealKind) :: Ushear
  REAL(RealKind) :: Frad,Fwat,Ftem,Fhum
  REAL(RealKind) :: TAir
  REAL(RealKind) :: WSroot,WSoil
  REAL(RealKind) :: wTLP,Epotg
  REAL(RealKind) :: LAIndex
  REAL(RealKind) :: PlantArea,InterceptArea
  INTEGER        :: SoilType,LandType,iS
! turgor loss point of plants:
  wTLP = cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType)) &
!                        *(0.81d0+0.121*ATAN(86400.d0*Epotg-Epotnorm))
                        *(0.81d0+0.121*ATAN(Epotg-Epotnorm/86400.d0))
! describes radiative influence on stomatal resistance,
  Frad = MIN(1.0d0,0.55d0*RadPAR/RadPARcrit) ! RadPAR==0.55*Qdir ! MJ: 55% d. Globalstrahlung -> Qdir+Qdif ?
! ... of soil water content,
  Fwat = MAX(Zero,MIN(1.0d0,(WSroot-cpwp(SoilType))/(wTLP-cpwp(SoilType))))
! ... of ambient specific humidity,
  Fhum = 1.0d0
! ... and of ambient temperature.
  Ftem = Max(Zero,MIN(1.0d0,4.0d0*((TAir-Temp0)*(TempEnd-TAir))/(TempEnd-Temp0)**2.0d0))
! resistance between the atmosphere and the air between the leafs
  AtmoResist   = DragCoeff**(-1.0d0)
! resistance between leafs and the air between the leafs
  Resist_la    = (MAX(Cdash*Ushear**(Half),1.0d-6))**(-1.0d0)
! stomatal resistance
  StomaResist  = 1.0d0/( ResMax**(-1.0d0)+(ResMin**(-1.0d0)-ResMax**(-1.0d0))      &
                                *Frad*Fwat*Fhum*Ftem )
! reduction of transpiration due to stomatal resistance
  ReductTransp = Resist_la*(Resist_la+StomaResist)**(-1.0d0)
! total resistance = foliage resistance
  FolResist    = 1.0d0/MAX(Eps,(ReductTransp*LAIndex*Resist_la**(-1.0d0)))
  IF (WSoil.LE.cpwp(SoilType)) THEN
    TranspirFunction = Zero
!   WRITE (*,*) 'Tr ZERO'
  ELSE
    TranspirFunction = (One-InterceptArea)*PlantArea*            &
                        Epotg*AtmoResist/(AtmoResist+FolResist)  &
                       *MIN(dzSoil(iS),MAX(0.d0,(z_root(LandType)-zSoil(iS-1)))) &
                       /(z_root(LandType)+Eps)  &
                       *WSoil/(WSroot+Eps)
!   WRITE (*,*) ' WSoil,cpwp(SoilType),SoilType,Epotg,AtmoResist', WSoil,cpwp(SoilType),SoilType,Epotg,AtmoResist,TranspirFunction
!   WRITE (*,*) 'Tr calc'
  END IF
END FUNCTION TranspirFunction

FUNCTION VolCellGlobe(r1,r2,phi1,phi2,lam1,lam2)
  REAL(RealKind) :: VolCellGlobe
  REAL(RealKind) :: r1,r2,phi1,phi2,lam1,lam2
  VolCellGlobe=MAX(1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1))),Zero)
END FUNCTION VolCellGlobe

SUBROUTINE VectorInitSoilFunction(Pos,Vector,Val,Time)
!
!==================================================
!----  Initialization of Scalar Values
!==================================================
!
  INTEGER :: Pos
  TYPE (Vector4Cell_T), TARGET :: Vector(:)
  REAL(RealKind) :: Time
  REAL(RealKind) :: Val
  EXTERNAL Val
  INTEGER :: i,iS
  REAL(RealKind) :: zPS
  REAL(RealKind), POINTER :: cB(:,:) 
    
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    cB=>Vector(ibLoc)%Vec(Pos)%cB
    DO i=1,NumBoundCell
      IF (Pos==RhoVPos) THEN
        DO iS=1,SIZE(cB(i,:))-1
          zPS=zSoil(is)-Half*dzSoil(is)
          cB(i,iS) = Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                        ,zH(BoundCell(i)%ix,BoundCell(i)%iy),zPS &
!                        ,zH(BoundCell(i)%xS,BoundCell(i)%yS),zPS &
                        ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(iS),Time)
        END DO
          cB(i,SIZE(cB(i,:)))=Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                                 ,zH(BoundCell(i)%ix,BoundCell(i)%iy),Zero &
!                                 ,zH(BoundCell(i)%xS,BoundCell(i)%yS),Zero &
                                 ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(1),Time)
      ELSE
        DO iS=1,SIZE(cB(i,:))
          zPS=zSoil(is)-Half*dzSoil(is)
          cB(i,iS) = Val(BoundCell(i)%xS,BoundCell(i)%yS,BoundCell(i)%zS &
                        ,zH(BoundCell(i)%ix,BoundCell(i)%iy),zPS &
!                        ,zH(BoundCell(i)%xS,BoundCell(i)%yS),zPS &
                        ,BoundCell(i)%LandClass,BoundCell(i)%SoilType(iS),Time)
        END DO
      END IF
    END DO
  END DO   ! ibLoc

END SUBROUTINE VectorInitSoilFunction


SUBROUTINE DistrFlux(Rhs,FU,FV,FW,Vol,FL,dx,dy,dz,n1,n2,n3,TotalFlux)
  REAL(RealKind) :: Rhs(-1:1,-1:1,-1:1)
  REAL(RealKind) :: FU(-1:0),FV(-1:0),FW(-1:0),Vol(-1:1,-1:1,-1:1),FL
  REAL(RealKind) :: dx,dy,dz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: TotalFlux

  REAL(RealKind) :: dLoc,VolLoc,VolDiff,SurFace

  IF (FluxDistribute) THEN
    dLoc=ABS(n1)*dx+ABS(n2)*dy+ABS(n3)*dz
    VolLoc=dx*dy*dz
    VolLoc=MIN(FL*dLoc,VolLoc)
    Surface=ABS(FU(-1)-FU(0)) &
           +ABS(FV(-1)-FV(0)) &
           +ABS(FW(-1)-FW(0))
    Rhs(0,0,0)=Rhs(0,0,0)+MIN(Vol(0,0,0)/VolLoc,One)*TotalFlux
    VolDiff=MAX(VolLoc-Vol(0,0,0),Zero)
    Rhs(-1,0,0)=Rhs(-1,0,0)+MAX(FU(-1)-FU(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(1,0,0)=Rhs(1,0,0)+MAX(FU(0)-FU(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,-1,0)=Rhs(0,-1,0)+MAX(FV(-1)-FV(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,1,0)=Rhs(0,1,0)+MAX(FV(0)-FV(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,0,-1)=Rhs(0,0,-1)+MAX(FW(-1)-FW(0),Zero)/Surface*VolDiff/VolLoc*TotalFlux
    Rhs(0,0,1)=Rhs(0,0,1)+MAX(FW(0)-FW(-1),Zero)/Surface*VolDiff/VolLoc*TotalFlux
  ELSE
    Rhs(0,0,0)=Rhs(0,0,0)+TotalFlux
  END IF
END SUBROUTINE DistrFlux

END MODULE Soil_Mod
