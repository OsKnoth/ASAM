MODULE MarineEmission_Mod

  USE EmissDeposParameter_Mod

  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: Species_typ(:) ! 1:Salt 2:OM
  REAL(RealKind), ALLOCATABLE :: GesFrac(:)
  REAL(RealKind) :: dr80drform, drformdrdry, Omegah
  REAL(RealKind) :: Rho_SS=2.165d3 ! salt density [kg/m^3] ! more exakt for seasalt in general; 2.2 is only NaCl
  REAL(RealKind) :: Rho_W=1.d3 ! water density [kg/m^3]
  REAL(RealKind) :: Rho_OM=1.3d3 ! density of the organic material [kg/m^3] ! not sure, value could be between 1.0 and 1.3
  REAL(RealKind) :: Chl_a=0.d0 ! Chlorophyllconcentration in sea water [mg/L]
  REAL(RealKind) :: SedExp ! Exponent for Sedimentationkorrektion
  LOGICAL :: PMOM=.FALSE.

CONTAINS

SUBROUTINE SetMarineEmission(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: i, j
  REAL(RealKind) :: c1,c2,c3,c4,c5,c6,c7,c8
  REAL(RealKind) :: Temp=288.15d0-273.15d0
  CHARACTER(20) :: S1,S2,End,SpeciesName
  LOGICAL :: Back
  REAL(RealKind) :: Anteil(2,nAqua)
  REAL(RealKind) :: SaltChargeGes,Work, Sort, Anteilsumme(2)
  REAL(RealKind) :: salinity=35.d0 ! salinity [promille]

! Inpur of gas emissions from the oceans
  S1='BEGIN_GAS'
  S2='BEGIN_SEATRANS'
  End='END_SEATRANS'
  CALL OpenFile(FileName)
  NumSeaTrans=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      NumSeaTrans=NumSeaTrans+1
    END IF
  END DO
  CALL CloseFile
  ALLOCATE(SeaGasEmi(NumSeaTrans))
  NumSeaTrans=0

  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=c1,R2=c2,R3=c3,R4=c4,R5=c5,R6=c6,R7=c7,R8=c8)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      NumSeaTrans=NumSeaTrans+1
      SeaGasEmi(NumSeaTrans)%Pos=Position(SpeciesName)
      SeaGasEmi(NumSeaTrans)%Par=c1/AVGA*1.d6
      SeaGasEmi(NumSeaTrans)%Konz=c2*c7
      SeaGasEmi(NumSeaTrans)%Schmidt=SchmidtGas(Temp,c3,c4,c5,c6)
      SeaGasEmi(NumSeaTrans)%Henry=c7
      SeaGasEmi(NumSeaTrans)%Depos=c8
    END IF
  END DO
  CALL CloseFile

! Initializing of freshly emitted sea salt composition
  ALLOCATE(Species_typ(nAqua+1))
  Species_typ=0
  Anteil=0.d0
  S1='BEGIN_AERO'
  S2='BEGIN_SALT'
  End='END_SALT'
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=Work,R2=Sort)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      Species_typ(Position(SpeciesName))=NINT(Sort)
      Anteil(Species_typ(Position(SpeciesName)),Position(SpeciesName))=Work
    END IF
  END DO
  CALL CloseFile

  SaltChargeGes=Zero
  DO i=1,nAqua
    DO j=1,2
      SaltChargeGes=SaltChargeGes+Charge(i)*Anteil(j,i)/MolMass(i)
    END DO
  END DO
  IF (Neutral) THEN
    IF (SaltChargeGes>Zero) THEN
      Anteil(1,Position('OHm'))=SaltChargeGes*MolMass(Position('OHm'))+Anteil(1,Position('OHm'))
    ELSE
      Anteil(1,Position('Hp'))=-SaltChargeGes*MolMass(Position('Hp'))+Anteil(1,Position('Hp'))
    END IF
  END IF

  ALLOCATE(GesFrac(nAqua+1))
  GesFrac=0.0d0
  Anteilsumme(1)=MAX(SUM(Anteil(1,:)),1.d-9)
  Anteilsumme(2)=MAX(SUM(Anteil(2,:)),1.d-9)
  DO i=1,nAqua
    GesFrac(i)=MAX(Anteil(1,i)/Anteilsumme(1),Anteil(2,i)/Anteilsumme(2))
  END DO

! Calculation of not changed variables
! Calculate conversion factor of dry sea salt to humidity of 80%
  dr80drform=1.d0/3.7d0*((2.d0-8.d-1)/(1.d0-8.d-1))**(1.d0/3.d0) ! Glg. 3 from Lewis and Schwartz 2006 
  drformdrdry=(Rho_SS/Rho_W*1.d3/salinity)**(1.d0/3.d0) ! Lewis and Schwartz 2004 page 54 left top
  Omegah=dr80drform*drformdrdry

! set Chlorophyllparameter to zero if orgmaterial tracer is not present in inputfile
  IF (POSITION('sWISOC')==0) THEN
    Chl_a=0.d0
    IF (MAXVAL(Species_typ(:)).GE.2) PMOM=.TRUE.
  END IF

  CALL CloseFile

END SUBROUTINE SetMarineEmission

SUBROUTINE SeaSaltEmission(U10,SST,ustern,zPL_zRauh,FL,ix,iy,iz)

  INTEGER :: ix,iy,iz
  INTEGER :: nD
  REAL(RealKind), INTENT (IN) :: U10, SST ! 10m windspeed, watertemperature 
  REAL(RealKind), INTENT (IN) :: ustern,zPL_zRauh
  REAL(RealKind) :: radtL, radtR, dr, interv ! lower, upper radius of the dry particle
  REAL(RealKind) :: rL, rR, rC ! radius at left, right side and in the center of the bin of dry particle
  REAL(RealKind) :: Dpd_L, Dpd_R, Dpd_C ! diameter of dry particle
  REAL(RealKind) :: Dp80_L, Dp80_R, Dp80_C ! particle diameter at RH of 80%
  REAL(RealKind) :: VSS_L, VSS_R, VSS_C ! volume of sea salt in the particle
  REAL(RealKind) :: VOM_L, VOM_R, VOM_C ! volume of organic material in the particle
  REAL(RealKind) :: dFNdDpdrya, dFNdDpdryb
  REAL(RealKind) :: FL

  INTEGER :: i, is
  REAL(RealKind), ALLOCATABLE :: mD(:)
  REAL(RealKind), ALLOCATABLE :: NumD(:),MassD(:,:)
  REAL(RealKind), ALLOCATABLE :: c(:,:)
  REAL(RealKind) :: MassGes
  REAL(RealKind), POINTER :: RhoLoc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVLoc(:,:,:,:)
  REAL(RealKind), POINTER :: Temp(:,:,:,:)
  REAL(RealKind) :: SattW, eSatt
  REAL(RealKind) :: RhoPart=1.2d3 ! kg/m^3
  REAL(RealKind) :: Vsed

  RhoLoc=>RhoC_T%c
  RhoVLoc=>cVec(RhoVPos)%c
  Temp=>TAbs%c

! set parameter for the flux calculation
  interv=4.d4
  nD=NINT(interv)
  radtL=5.d-9
  radtR=5.d-5
  dr=(LOG(radtR)-LOG(radtL))/interv
  rL=radtL


  ALLOCATE(mD(nD+1))
  ALLOCATE(NumD(nD))
  NumD=0.d0
  ALLOCATE(MassD(nD,nAqua))
  ALLOCATE(c(nFrac,nAqua))

! Flux at the left side of the first bin
  Dpd_L=rL*2.d6
  CALL SSA_Volume(Dpd_L,U10,Dp80_L,VOM_L,VSS_L) ! Volume of sea salt and organic material
  Vsed=VFinalF(RhoPart,RhoLoc(ix,iy,iz,1),Dp80_L/2.d6)
  mD(1)=Rho_SS*VSS_L+Rho_OM*VOM_L ! kg
  dFNdDpdrya=SSA_EmissFlux(Dp80_L,Dpd_L,U10,SST)*Dp80_L/Dpd_L ! numberflux at the lower limit
  dFNdDpdrya = dFNdDpdrya * (zPL_zRauh)**(-Vsed/(Karm*ustern)) ! F_int => F_eff (half height of lowest grid cell)

  DO i=1,nD

    rR=EXP(LOG(rL)+dr)
    rC=5.d-1*(rL+rR)

    Dpd_R=rR*2.d6
    CALL SSA_Volume(Dpd_R,U10,Dp80_R,VOM_R,VSS_R)
    Vsed=VFinalF(RhoPart,RhoLoc(ix,iy,iz,1),Dp80_R/2.d6)
    dFNdDpdryb=SSA_EmissFlux(Dp80_R,Dpd_R,U10,SST)*Dp80_R/Dpd_R ! numberflux at the upper limit
    dFNdDpdryb = dFNdDpdryb * (zPL_zRauh)**(-Vsed/(Karm*ustern)) ! F_int => F_eff (half height of lowest grid cell)

    NumD(i)=(dFNdDpdrya+dFNdDpdryb)*(Dpd_R-Dpd_L)/2.d0 ! total numberflux [m^-2*s^-2]

    Dpd_C=rC*2.d6
    CALL SSA_Volume(Dpd_C,U10,Dp80_C,VOM_C,VSS_C)

    MassD(i,1)=Rho_SS*VSS_C+Rho_OM*VOM_C*NumD(i)  ! kg/m^2/s
    MassGes=MassD(i,1)  ! mass of all particles kg/m^2/s 

    DO is=1,nAqua
      IF (Species_typ(is).EQ.1) THEN
        MassD(i,is)=Rho_SS*VSS_C*NumD(i)*GesFrac(is) ! Mass of ions included in sea salt
      ELSE IF (is==POSITION('sWISOC') .OR. Species_typ(is).EQ.2) THEN
        MassD(i,is)=Rho_OM*VOM_C*NumD(i)*GesFrac(is) ! mass of organic enrichment
      END IF
    END DO

! right side is the new left side
    dFNdDpdrya=dFNdDpdryb
    rL=rR
    Dpd_L = Dpd_R
    VOM_L = VOM_R
    VSS_L = VSS_R
    mD(i+1)=Rho_SS*VSS_R+Rho_OM*VOM_R ! kg

  ENDDO

! add water to the emitted particles
  eSatt=SaturVapor(Temp(ix,iy,iz,1))
  SattW=(Rv*Temp(ix,iy,iz,1)*RhoVLoc(ix,iy,iz,1))/eSatt-1.d0
  CALL InitDroplet(nFrac,m,c &
                  ,nD,mD,NumD,MassD,SattW,Temp(ix,iy,iz,1))
! c: kg/m²/s
  DO is=1,nAqua
    fVec(is)%c(ix,iy,iz,:)=fVec(is)%c(ix,iy,iz,:)+c(:,is)*FL/(VolC(ix,iy,iz)+Eps)*RhoLoc(ix,iy,iz,1)
  END DO

! fVec: kg/m^3/s

  DEALLOCATE(c)
  DEALLOCATE(mD)
  DEALLOCATE(NumD)
  DEALLOCATE(MassD)

END SUBROUTINE SeaSaltEmission


FUNCTION SSA_EmissFlux(Dp80,Dpdry,U10,Tw)

  REAL(RealKind) :: Dp80, Dpdry, U10, Tw ! Units: µm, µm, m/s, K
  REAL(RealKind) :: SSA_EmissFlux
  REAL(RealKind) :: Wind_part, Size_part, SST_part
  REAL(RealKind) :: Sofiev_1, Sofiev_2
  REAL(RealKind) :: P_Long

  Wind_part = 2.d-8*U10**3.74d0

  IF (Tw.LE.278.15d0) THEN
    Sofiev_1 = 0.092d0*Dpdry**(-0.96d0) ! valid for Tw=271.15K ; -2°C
    Sofiev_2 = 0.15d0*Dpdry**(-0.88d0) ! valid for Tw=278.15K ; 5°C
    SST_Part = (Sofiev_1*(278.15d0-Tw)+Sofiev_2*(Tw-271.15d0))/7.d0
  ELSE IF (Tw.GT.278.15d0 .AND. Tw.LE.288.15d0) THEN
    Sofiev_1 = 0.15d0*Dpdry**(-0.88d0) ! valid for Tw=278.15K ; 5°C
    Sofiev_2 = 0.48d0*Dpdry**(-0.36d0) ! valid for Tw=288.15K ; 15°C
    SST_Part = (Sofiev_1*(288.15d0-Tw)+Sofiev_2*(Tw-278.15d0))/1.d1
  ELSE IF (Tw.GT.288.15d0) THEN
    Sofiev_1 = 0.48d0*Dpdry**(-0.36d0) ! valid for Tw=288.15K ; 15°C
    Sofiev_2 = 1.d0 ! valid for Tw=298.15K ; 25°C
    SST_Part = (Sofiev_1*(298.15d0-Tw)+Sofiev_2*(Tw-288.15d0))/1.d1
  END IF

  IF (Dp80<=One) THEN
    P_Long=1.46d0*(Log10(Dp80))**3+1.33d0*(Log10(Dp80))**2-1.82d0*(Log10(Dp80))+8.83d0
  ELSE
    P_Long=-1.53d0*(Log10(Dp80))**3-8.1d-2*(Log10(Dp80))**2-4.26d-1*(Log10(Dp80))+8.84d0
  ENDIF
  Size_part=1.d1**(P_Long) ! dF/dlogDp80
  Size_part=Size_part/(Dp80*Log(1.d1)) ! dF/dDp80

  SSA_EmissFlux = Wind_part*Size_part*SST_Part

END FUNCTION SSA_EmissFlux

!--------------------------------------------------------------------------------------------
SUBROUTINE SSA_Volume(Dpd,U10,Dp80,VOM,VSS) ! mass relation with OM_SSA=M_OM/(M_OM+M_SS)

  REAL(8), INTENT (IN) ::  Dpd,U10
  REAL(8), INTENT (OUT) :: Dp80,VOM,VSS

  REAL(8) :: Vp80, VSS80, Vpd, Vp80_alt
  REAL(8) :: Teil_1, Teil_2, Teil_3
  REAL(8) :: OM_SSA
  REAL(8) :: Rho_OM = 1.3E3
  REAL(8) :: Rho_SS = 2.165E3
  INTEGER :: i

  Pi=Four*atan(One)

  Dp80=Dpd*Omegah            ! µm
  IF (PMOM) THEN
    Vpd=Pi/6.d0*Dpd**3     ! µm^3
    Vp80=Pi/6.d0*Dp80**3   ! µm^3

    Teil_1 = One/(One+exp(-2.63E0*Chl_a+0.18E0*U10))
    Teil_3 = 0.03E0/(One+exp(-2.63E0*Chl_a+0.18E0*U10))

    Vp80_alt=Zero
    DO i=1,20
      Teil_2 = One+0.03E0*exp(6.81E0*Dp80)
      OM_SSA = Teil_1/Teil_2 + Teil_3
      VOM = Vpd*Rho_SS*OM_SSA/((Rho_SS-Rho_OM)*OM_SSA+Rho_OM)
      VSS = Vpd - VOM
      VSS80 = VSS*Omegah**3
      Vp80=VSS80+VOM
      Dp80=(6.d0/Pi*Vp80)**(One/Three)
      IF (Vp80.EQ.Vp80_alt) EXIT
      Vp80_alt=Vp80
    END DO
  ELSE
    VSS=Pi/6.d0*Dpd**3  ! µm^3
    VOM=Zero              ! µm^3
  END IF

  VSS = VSS*1.d1**(-1.8E1) ! µm^3 > m^3
  VOM = VOM*1.d1**(-1.8E1) ! µm^3 > m^3

END SUBROUTINE SSA_Volume

SUBROUTINE SeaGasEmission(U10,FL,ix,iy,iz)

  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: U10
  REAL(RealKind) :: FL
  REAL(RealKind), Pointer :: Ca
  REAL(RealKind) :: Kw,Fg

  DO i=1,NumSeaTrans

    IF (SeaGasEmi(NumSeaTrans)%Pos>0) THEN
      IF (SeaGasEmi(NumSeaTrans)%Schmidt.NE.0.d0) THEN
        Kw=TransCoeff(U10,SeaGasEmi(NumSeaTrans)%Schmidt)
        Ca=>cVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,1)

        Fg=Kw*(SeaGasEmi(NumSeaTrans)%Konz-Ca)*RhoC_T%c(ix,iy,iz,1) ! mol/m²/s
!    Komment: Henrykonstante: dimensionslose Form, daher kann Ca ursprüngliche Größe bleiben
        fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)=fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:) &
                                +Fg*FL/(VolC(ix,iy,iz)+Eps) ! fVec mol/m^3/s
      ELSE
        Fg=SeaGasEmi(NumSeaTrans)%Par*RhoC_T%c(ix,iy,iz,1)
        fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)=fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)+Fg ! fVec mol/m^3/s
      END IF
    END IF

  END DO

END SUBROUTINE SeaGasEmission


END MODULE MarineEmission_Mod
