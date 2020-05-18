MODULE MarineEmission_Mod

  USE EmissDeposParameter_Mod
  USE Control_Mod,      ONLY: RhoVpos

  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: Species_typ(:) ! 1:Salt 2:OM
  REAL(RealKind), ALLOCATABLE :: GesFrac(:)
  REAL(RealKind) :: dr80drform, drformdrdry, Omegah
  REAL(RealKind) :: Rho_SS=2.165d3 ! salt density [kg/m^3] ! more exakt for seasalt in general; 2.2 is only NaCl
  REAL(RealKind) :: Rho_W=1.d3 ! water density [kg/m^3]
  REAL(RealKind) :: Rho_SW=1.024d3 ! sea water density [kg/m^3]
  REAL(RealKind) :: Rho_OM=1.3d3 ! density of the organic material [kg/m^3] ! not sure, value could be between 1.0 and 1.3
  REAL(RealKind) :: Chl_a=0.d0 ! Chlorophyllconcentration in sea water [mg/L]
  REAL(RealKind) :: salinity=3.5d-2 ! salinity [g/g]
  REAL(RealKind) :: SedExp ! Exponent for Sedimentationkorrektion
  LOGICAL :: PMOM=.FALSE.
  REAL, PARAMETER :: Avog = 6.0221e23

  REAL(RealKind), ALLOCATABLE :: Ads_koeff(:) ! mol/m^3
  REAL(RealKind), ALLOCATABLE :: Conc_SML(:), Conc_bulk(:) ! mol/m^3
  REAL(RealKind), ALLOCATABLE :: Spec_area(:) ! kg/mol
  REAL(RealKind), ALLOCATABLE :: Mol_area(:) ! kg/mol
  REAL(RealKind), ALLOCATABLE :: Dens(:) ! kg/m^3
  REAL(RealKind), ALLOCATABLE :: MolVol(:) ! m^3/mol
  REAL(RealKind), ALLOCATABLE :: Rel_Vol(:,:)
  REAL(RealKind), ALLOCATABLE :: BulkFrac_Vol(:), SMLFrac_Vol(:)

CONTAINS

SUBROUTINE SetMarineEmission(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: i, j
  REAL(RealKind) :: c1,c2,c3,c4,c5,c6,c7,c8
  REAL(RealKind) :: Temp=288.15d0-273.15d0
  CHARACTER(20) :: S1,S2,End,SpeciesName
  LOGICAL :: Back
  REAL(RealKind) :: SaltChargeGes,Work, Sort, Anteilsumme(2)
  REAL(RealKind), ALLOCATABLE :: Conc_temp_bulk(:,:)
  REAL(RealKind), ALLOCATABLE :: Conc_temp_SML(:,:)
  REAL(RealKind) :: Salt_conc_bulk, Salt_conc_SML ! Total salt concentration in input file

! Input of gas emissions from the oceans
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

! Initializing composition of the emitted sea spray aerosol
  ALLOCATE(Species_typ(nAqua))
  Species_typ=0
  ALLOCATE(Conc_bulk(nAqua))
  Conc_bulk=0.d0
  ALLOCATE(Conc_SML(nAqua))
  Conc_SML=0.d0
  ALLOCATE(Spec_area(nAqua))
  Spec_area=0.d0
  ALLOCATE(Mol_area(nAqua))
  Mol_area=0.d0
  ALLOCATE(Ads_koeff(nAqua))
  Ads_koeff=0.d0
  ALLOCATE(Dens(nAqua))
  Dens=0.d0
  Salt_conc_bulk=0.d0
  Salt_conc_SML=0.d0
  S1='BEGIN_AERO'
  S2='BEGIN_SALT'
  End='END_SALT'
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=Sort,R2=c1,R3=c2,R4=c3,R5=c4,R6=c5)
    IF (Back) THEN
      EXIT
    END IF
    IF (Position(SpeciesName)>0) THEN
      Species_typ(Position(SpeciesName))=NINT(Sort)
      Conc_bulk(Position(SpeciesName))=c1
      Conc_SML(Position(SpeciesName))=c2
      Spec_area(Position(SpeciesName))=c3
      Ads_koeff(Position(SpeciesName))=c4
      Dens(Position(SpeciesName))=c5
    END IF
    IF (NINT(Sort).EQ.1) THEN ! Calculate total salt concentration
      Salt_conc_bulk=Salt_conc_bulk+c1*1.d-3
      Salt_conc_SML=Salt_conc_SML+c2*1.d-3
    END IF
  END DO
  CALL CloseFile

! Convert variables from input units to units needed for modelrun
  Spec_area = Spec_area*1.d-20       ! m^2/molec. = A^2/molec. * 10^-20 m^2/A^2
  Mol_area  = Spec_area*Avog      ! m^2/mol = m^2/molec. * molec./mol
  Conc_bulk = Conc_bulk*1.d-3          ! kg/m^3 = mg/L * 10^3 L/m^3 * 10^-6 kg/mg
  Conc_SML  = Conc_SML*1.d-3            ! kg/m^3 = mg/L * 10^3 L/m^3 * 10^-6 kg/mg
  Ads_koeff(:) = Ads_koeff(:)/MolMass(1:nAqua)       ! m^3/kg = m^3/mol * mol/kg
  DEALLOCATE(Spec_area)

! Calculate molare Volume
  ALLOCATE(MolVol(nAqua))
  WHERE (Dens.LE.0.d0) Dens=1.d3
  MolVol(:)=MolMass(1:nAqua)/Dens(:)
  
! Calculate Charge and neutralise aerosols
  SaltChargeGes=SUM(Charge(1:nAqua)*Conc_bulk(1:nAqua)/MolMass(1:nAqua))
  IF (Neutral) THEN
    IF (SaltChargeGes>Zero) THEN
      Conc_bulk(Position('OHm')) = SaltChargeGes*MolMass(Position('OHm'))+Conc_bulk(Position('OHm'))
    ELSE
      Conc_bulk(Position('Hp'))  = -SaltChargeGes*MolMass(Position('Hp'))+Conc_bulk(Position('Hp'))
    END IF
  END IF

  SaltChargeGes=SUM(Charge(1:nAqua)*Conc_SML(1:nAqua)/MolMass(1:nAqua))
  IF (Neutral) THEN
    IF (SaltChargeGes>Zero) THEN
      Conc_SML(Position('OHm')) = SaltChargeGes*MolMass(Position('OHm'))+Conc_SML(Position('OHm'))
    ELSE
      Conc_SML(Position('Hp'))  = -SaltChargeGes*MolMass(Position('Hp'))+Conc_SML(Position('Hp'))
    END IF
  END IF

! Calculate relative Concentrations
  ALLOCATE (Conc_temp_bulk(2,nAqua))
  Conc_temp_bulk=0.0d0
  ALLOCATE (Conc_temp_SML(2,nAqua))
  Conc_temp_SML=0.0d0
  WHERE (Species_typ(:).EQ.1)
    Conc_temp_bulk(1,:)=Conc_bulk(:)
    Conc_temp_SML(1,:)=Conc_SML(:)
  ELSE WHERE (Species_typ(:).EQ.2)
    Conc_temp_bulk(2,:)=Conc_bulk(:)
    Conc_temp_SML(2,:)=Conc_SML(:)
  END WHERE
  WHERE (Species_typ(:).EQ.1)
    Conc_bulk(:)=Conc_temp_bulk(1,:)*Salt_conc_bulk/SUM(Conc_temp_bulk(1,:))
    Conc_SML(:)=Conc_temp_SML(1,:)*Salt_conc_SML/SUM(Conc_temp_SML(1,:))
  END WHERE

  salinity=Salt_conc_bulk*1.0d-3

  !Anteilsumme(1)=MAX(SUM(Conc_temp_bulk(1,:)*MolVol(:)),1.d-9)
  !Anteilsumme(2)=MAX(SUM(Conc_temp_bulk(2,:)*MolVol(:)),1.d-9)
  !ALLOCATE(GesFrac(nAqua))
  !GesFrac(:)=MAX((Conc_temp_bulk(1,:)*MolVol(:))/Anteilsumme(1), & ! Relative amount of salts to total salt concentration
!&                (Conc_temp_bulk(2,:)*MolVol(:))/Anteilsumme(2)) ! Relative amount of organic species to total concentration of organics
  DEALLOCATE (Conc_temp_bulk)

  ALLOCATE (BulkFrac_Vol(nAqua))
  BulkFrac_Vol(:)=Conc_bulk(:)/Dens(1:nAqua)  ! m^3/m^3 = kg/m^3 / kg/m^3
  IF (SUM(BulkFrac_Vol(:)).GT.0.d0) BulkFrac_Vol(:)=BulkFrac_Vol(:)/SUM(BulkFrac_Vol(:))

  ALLOCATE (SMLFrac_Vol(nAqua))
  SMLFrac_Vol(:)=Conc_SML(:)/Dens(1:nAqua)  ! m^3/m^3 = kg/m^3 / kg/m^3
  IF (SUM(SMLFrac_Vol(:)).GT.0.d0) SMLFrac_Vol(:)=SMLFrac_Vol(:)/SUM(SMLFrac_Vol(:))

! Calculate organic amount emitted with film droplets (here as long as regional and temporal constant)
  ALLOCATE (Rel_Vol(nAqua,3))
  Rel_Vol=0.d0
  CALL Calc_Volume_ratio_film
  DEALLOCATE (Mol_area)
  DEALLOCATE (Ads_koeff)
  DEALLOCATE (MolVol)

! Calculation of not changed variables
! Calculate conversion factor of dry sea salt to humidity of 80%
  dr80drform=1.d0/3.7d0*((2.d0-8.d-1)/(1.d0-8.d-1))**(1.d0/3.d0) ! Glg. 3 from Lewis and Schwartz 2006 
  drformdrdry=(Rho_SS/Rho_W*1.d3/(salinity*1.d3))**(1.d0/3.d0) ! Lewis and Schwartz 2004 page 54 left top
  Omegah=dr80drform*drformdrdry

! set Chlorophyllparameter to zero if orgmaterial tracer is not present in inputfile
  IF (POSITION('sWISOC')==0) THEN
    Chl_a=0.d0
    IF (MAXVAL(Species_typ(:)).GE.2) PMOM=.TRUE.
  END IF

  CALL CloseFile

END SUBROUTINE SetMarineEmission

SUBROUTINE Calc_Volume_ratio_film

  IMPLICIT NONE

  INTEGER :: i

  REAL(RealKind) :: Sum_Ads_SML, Sum_Ads_bub ! einheitenlos, einheitenlos
  REAL(RealKind) :: SML_Flaeche, SML_bub_faktor, Vol_all_Subst, Vol_Film_SS ! m^2/mol, einheitenlos, m^3/mol , m^3/mol
  REAL(RealKind) :: Ads_eq_SML(nAqua), Ads_eq_bub(nAqua) ! einheitenlos, einheitenlos
  REAL(RealKind) :: Langmuir_rel_SML(nAqua), Langmuir_rel_bub(nAqua) ! einheitenlos, einheitenlos
  REAL(RealKind) :: Vol_ML_subst(nAqua)
  REAL(RealKind) :: Vol_Film_subst(nAqua)

  REAL(RealKind) :: Sum_Ads_area
  REAL(RealKind) :: Sum_Ads_Vol

  REAL(RealKind) :: Film_dicke(2) ! m
  DATA  Film_dicke & ! m
&     / 5.d-8, 2.d-7 /

    Ads_eq_SML(:) = Ads_koeff(:)*Conc_SML(:)                                             ! einheitenlos = m^3/kg * kg/m^3
    Ads_eq_bub(:) = Ads_koeff(:)*Conc_bulk(:)                                            ! einheitenlos = m^3/kg * kg/m^3
    Sum_Ads_SML = SUM(Ads_eq_SML(:))                                                     ! einheitenlos
    Sum_Ads_bub = SUM(Ads_eq_bub(:))                                                     ! einheitenlos
    Langmuir_rel_SML(:) = Ads_eq_SML(:)/(1.d0+Sum_Ads_SML)                               ! einheitenlos
    Langmuir_rel_bub(:) = Ads_eq_bub(:)/(1.d0+Sum_Ads_bub)                               ! einheitenlos
    SML_Flaeche = SUM(Langmuir_rel_SML(:)*Mol_area(:))                                   ! m^2/mol = einheitenlos * m^2/mol
    SML_bub_faktor = SML_Flaeche/SUM(Langmuir_rel_bub(:)*Mol_area(:))                    ! einheitenlos = m^2/mol / m^2/mol
    Vol_ML_subst(:)=(Langmuir_rel_SML(:)+Langmuir_rel_bub(:)*SML_bub_faktor)*MolVol(:)  ! m^3/mol = einheitenlos * m^3/mol

  DO i=1,2
    Vol_Film_subst(:)=Conc_SML(:)*Film_dicke(i)*SML_Flaeche/Dens(:)              ! m^3/mol = kg/m^3 * m * m^2/mol / kg/m^3
    Vol_all_Subst = SUM(Vol_ML_subst(:)+Vol_Film_subst(:))                               ! m^3/mol = m^3/mol + m^3/mol ...
    Rel_Vol(:,i)=(Vol_ML_subst(:)+Vol_Film_subst(:))/Vol_all_Subst         ! einheitenlos = m^3/mol + m^3/mol / (m^3/mol + m^3/mol)
  END DO


END SUBROUTINE Calc_Volume_ratio_film

SUBROUTINE Calc_Volume_ratio_jet(Dp_dry)

  IMPLICIT NONE

  REAL(RealKind), INTENT (IN) :: Dp_dry

  INTEGER :: i

  REAL(RealKind) :: SML_depth=1.d2 ! µm
  REAL(RealKind) :: R_bub
  REAL(RealKind) :: PreFak_A, PreFak_B

  R_bub=1.d1*Dp_dry*4.d0*5.d-1

  PreFak_A=MIN(MAX(((SML_depth-R_bub)/SML_depth+5.d-1)*2.d0,0.d0),2.d0)
  PreFak_B=MIN(MAX(((R_bub-SML_depth)/SML_depth+5.d-1)*2.d0,0.d0),2.d0)
  Rel_Vol(:,3)=(PreFak_A*SMLFrac_Vol(:)+PreFak_B*BulkFrac_Vol(:))/2.d0

END SUBROUTINE Calc_Volume_ratio_jet


FUNCTION Salter15(Dp_dry,U10,SST,Mode)  RESULT(SS_num_flux)

  IMPLICIT NONE

  INTEGER, INTENT (IN)  :: Mode
  REAL(RealKind), INTENT (IN) :: Dp_dry, U10, SST

  REAL(RealKind) :: A(3), B(3), C(3), D(3), sigma(3), Dpmod(3)
  REAL(RealKind) :: Temp_Wind_para, SS_num_flux

  A =     (/  -5.2168d5,      0.d0, 0.d0    /)
  B =     (/  3.31725d7,   7.374d5, 1.421d4 /)
  C =     (/ -6.95275d8, -2.4803d7, 1.4662d7/)
  D =     (/  1.0684d10,  7.7373d8, 1.7075d8/)
  Dpmod = (/     9.5d-2,     6.d-1, 1.5d0   /)
  sigma = (/      2.1d0,    1.72d0, 1.6d0   /)

  SS_num_flux = 0.d0
  Temp_Wind_para = 2.d-8*U10**3.74*(A(Mode)*SST**3+B(Mode)*SST**2+C(Mode)*SST+D(Mode))
  SS_num_flux = SS_num_flux + Temp_Wind_para/(SQRT(2.d0*Pi)*log(sigma(Mode)))* & ! dFNdlogDpdry
&               exp(-1.d0/2.d0*((log(Dp_dry)-log(Dpmod(Mode)))/log(sigma(Mode)))**2)

END FUNCTION Salter15


SUBROUTINE SeaSaltEmission(U10,SST,ustern,zPL,zRauh,FL,ix,iy,iz)

  INTEGER :: ix,iy,iz
  INTEGER :: nD
  REAL(RealKind), INTENT (IN) :: U10, SST ! 10m windspeed, watertemperature 
  REAL(RealKind), INTENT (IN) :: ustern,zPL,zRauh
  REAL(RealKind) :: radtL, radtR, dr, interv ! lower, upper radius of the dry particle
  REAL(RealKind) :: rL, rR, rC ! radius at left, right side and in the center of the bin of dry particle
  REAL(RealKind) :: Dpd_L, Dpd_R, Dpd_C ! diameter of dry particle
  REAL(RealKind) :: Dp80_L, Dp80_R, Dp80_C ! particle diameter at RH of 80%
  REAL(RealKind) :: VSS_L, VSS_R, VSS_C ! volume of sea salt in the particle
  REAL(RealKind) :: VOM_L, VOM_R, VOM_C ! volume of organic material in the particle
  REAL(RealKind) :: dFN_dDpdrya(4), dFN_dDpdryb(4)
  REAL(RealKind) :: dFVol_dDpdrya(3), dFVol_dDpdryb(3)
  REAL(RealKind) :: FL

  INTEGER :: i, is
  REAL(RealKind), ALLOCATABLE :: mD(:)
  REAL(RealKind), ALLOCATABLE :: NumD(:),MassD(:,:)
  REAL(RealKind), ALLOCATABLE :: c(:,:)
  REAL(RealKind), ALLOCATABLE :: dFMass_dDpdry(:,:)
  REAL(RealKind) :: MassGes
  REAL(RealKind), POINTER :: RhoVLoc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoLoc(:,:,:,:)
  REAL(RealKind), POINTER :: Temp(:,:,:,:)
  REAL(RealKind) :: SattW, eSatt
  REAL(RealKind) :: RhoPart=1.2d3 ! kg/m^3
  REAL(RealKind) :: Vsed

  RhoLoc=>cVec(Rhopos)%c
  RhoVLoc=>cVec(RhoVpos)%c
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
  ALLOCATE(dFMass_dDpdry(2,nAqua))

! Flux at the left side of the first bin
  Dpd_L=rL*2.d6
  CALL Calc_Volume_ratio_jet(Dpd_L)

  Vsed=VFinalF(RhoPart,RhoLoc(ix,iy,iz,1),rL*2.d0)
  dFN_dDpdrya(1) = Salter15(Dpd_L,U10,SST,1) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
  dFN_dDpdrya(2) = Salter15(Dpd_L,U10,SST,2) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
  dFN_dDpdrya(3) = Salter15(Dpd_L,U10,SST,3) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
  WHERE (dFN_dDpdrya.LE.0.d0) dFN_dDpdrya = 0.d0
  dFN_dDpdrya(4)=SUM(dFN_dDpdrya(1:3))

  dFVol_dDpdrya(:) = dFN_dDpdrya(1:3)*4.d0*Pi/3.d0*rL**3
  dFMass_dDpdry(1,:) = (dFVol_dDpdrya(1)*Rel_Vol(:,1)+dFVol_dDpdrya(2)*Rel_Vol(:,2)+dFVol_dDpdrya(3)*Rel_Vol(:,3))*Dens(:)

  mD(1)=SUM(dFMass_dDpdry(1,:))/dFN_dDpdrya(4)

  DO i=1,nD

    rR=EXP(LOG(rL)+dr)
    rC=5.d-1*(rL+rR)

    Dpd_R=rR*2.d6
    CALL Calc_Volume_ratio_jet(Dpd_L)
    !CALL SSA_Volume(Dpd_R,U10,Dp80_R,VOM_R,VSS_R)
    Vsed=VFinalF(RhoPart,RhoLoc(ix,iy,iz,1),rR*2.d0)
    dFN_dDpdryb(1) = Salter15(Dpd_L,U10,SST,1) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
    dFN_dDpdryb(2) = Salter15(Dpd_L,U10,SST,2) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
    dFN_dDpdryb(3) = Salter15(Dpd_L,U10,SST,3) * (zPL/1.d1)**(-Vsed/(Karm*ustern))
    WHERE (dFN_dDpdryb.LE.0.d0) dFN_dDpdryb = 0.d0
    dFN_dDpdryb(4)=SUM(dFN_dDpdryb(1:3))

    dFVol_dDpdryb(:) = dFN_dDpdryb(1:3)*4.d0*Pi/3.d0*rL**3
    dFMass_dDpdry(2,:) = (dFVol_dDpdryb(1)*Rel_Vol(:,1)+dFVol_dDpdryb(2)*Rel_Vol(:,2)+dFVol_dDpdryb(3)*Rel_Vol(:,3))*Dens(:)

    NumD(i)=(dFN_dDpdrya(4)+dFN_dDpdryb(4))*(Dpd_R-Dpd_L)/2.d0 ! total numberflux [m^-2*s^-2]
    MassD(i,:)=(dFMass_dDpdry(1,:)+dFMass_dDpdry(2,:))*(Dpd_R-Dpd_L)/2.d0
    MassGes=SUM(MassD(i,:))  ! mass of all particles kg/m^2/s 

! right side is the new left side
    dFN_dDpdrya=dFN_dDpdryb
    dFMass_dDpdry(1,:)=dFMass_dDpdry(2,:)
    rL=rR
    Dpd_L = Dpd_R
    VOM_L = VOM_R
    VSS_L = VSS_R
    mD(i+1)=SUM(dFMass_dDpdry(2,:))/dFN_dDpdryb(4)

  ENDDO

! add water to the emitted particles
  eSatt=SaturVapor(Temp(ix,iy,iz,1))
  SattW=(Rv*Temp(ix,iy,iz,1)*RhoVLoc(ix,iy,iz,1))/eSatt-1.d0
  CALL InitDroplet(nFrac,m,c,nD,mD,NumD,MassD,SattW,Temp(ix,iy,iz,1))
! c: kg/m²/s
  DO is=1,nAqua
    fVec(is)%c(ix,iy,iz,:)=fVec(is)%c(ix,iy,iz,:)+c(:,is)*FL/(VolC(ix,iy,iz)+Eps)*rhoLoc(ix,iy,iz,1)
  END DO

! fVec: kg/m^3/s

  DEALLOCATE(dFMass_dDpdry)
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

        Fg=Kw*(SeaGasEmi(NumSeaTrans)%Konz-Ca)*cVec(RhoPos)%c(ix,iy,iz,1) ! mol/m²/s
!    Komment: Henrykonstante: dimensionslose Form, daher kann Ca ursprüngliche Größe bleiben
        fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)=fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:) &
                                +Fg*FL/(VolC(ix,iy,iz)+Eps) ! fVec mol/m^3/s
      ELSE
        Fg=SeaGasEmi(NumSeaTrans)%Par*cVec(RhoPos)%c(ix,iy,iz,1)
        fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)=fVec(SeaGasEmi(NumSeaTrans)%Pos)%c(ix,iy,iz,:)+Fg ! fVec mol/m^3/s
      END IF
    END IF

  END DO

END SUBROUTINE SeaGasEmission


END MODULE MarineEmission_Mod
