Module M2_Parameter
!=======================================================================!
! This file includes:
!   M2_LD_FCT.f90
!   M2_REL_EM_ACT.f90
!   M2_RHO_SPC.f90
!   M2_SPC_MGN.f90
!   M2_TEMPD_PRM.f90
!   canopy.txt
!=======================================================================!

!=======================================================================!
!  LD_FCT.EXT
!  This include file contains "light dependent" factors.
!
!  This is an input file for MEGAN v2.0
!  This is version v210 of this file
!  Created by A. Guenther 8/11/07
!=======================================================================!
INTEGER, PARAMETER ::   N_LDF_SPC=20
CHARACTER*16       ::   LDF_SPC(20)
REAL               ::   LDF_FCT(20)
INTEGER            ::   LDF_MAP(20)

DATA LDF_SPC(  1), LDF_FCT(  1), LDF_MAP(  1) / 'ISOP            ', 0.9999, 1 /
DATA LDF_SPC(  2), LDF_FCT(  2), LDF_MAP(  2) / 'MYRC            ', 0.05  , 2 /
DATA LDF_SPC(  3), LDF_FCT(  3), LDF_MAP(  3) / 'SABI            ', 0.1   , 3 /
DATA LDF_SPC(  4), LDF_FCT(  4), LDF_MAP(  4) / 'LIMO            ', 0.05  , 4 /
DATA LDF_SPC(  5), LDF_FCT(  5), LDF_MAP(  5) / '3CAR            ', 0.05  , 5 /
DATA LDF_SPC(  6), LDF_FCT(  6), LDF_MAP(  6) / 'OCIM            ', 0.8   , 6 /
DATA LDF_SPC(  7), LDF_FCT(  7), LDF_MAP(  7) / 'BPIN            ', 0.1   , 7 /
DATA LDF_SPC(  8), LDF_FCT(  8), LDF_MAP(  8) / 'APIN            ', 0.1   , 8 /
DATA LDF_SPC(  9), LDF_FCT(  9), LDF_MAP(  9) / 'OMTP            ', 0.1   , 9 /
DATA LDF_SPC( 10), LDF_FCT( 10), LDF_MAP( 10) / 'FARN            ', 0.5   , 10/
DATA LDF_SPC( 11), LDF_FCT( 11), LDF_MAP( 11) / 'BCAR            ', 0.5   , 11/
DATA LDF_SPC( 12), LDF_FCT( 12), LDF_MAP( 12) / 'OSQT            ', 0.5   , 12/
DATA LDF_SPC( 13), LDF_FCT( 13), LDF_MAP( 13) / 'MBO             ', 0.9999, 13/
DATA LDF_SPC( 14), LDF_FCT( 14), LDF_MAP( 14) / 'MEOH            ', 0.75  , 14/
DATA LDF_SPC( 15), LDF_FCT( 15), LDF_MAP( 15) / 'ACTO            ', 0.25  , 15/
DATA LDF_SPC( 16), LDF_FCT( 16), LDF_MAP( 16) / 'CH4             ', 0.75  , 16/
DATA LDF_SPC( 17), LDF_FCT( 17), LDF_MAP( 17) / 'NO              ', 0.0   , 17/
DATA LDF_SPC( 18), LDF_FCT( 18), LDF_MAP( 18) / 'ACTA            ', 0.5   , 18/
DATA LDF_SPC( 19), LDF_FCT( 19), LDF_MAP( 19) / 'FORM            ', 0.5   , 19/
DATA LDF_SPC( 20), LDF_FCT( 20), LDF_MAP( 20) / 'CO              ', 0.5   , 20/

!=======================================================================!
!  REL_EM_ACT.EXT
!  This include file contains "Emission factors for trees at different periods"
!  factors.
!
!  MEGAN v2.0
!
!  Created by Tan 11/30/06
!=======================================================================!
INTEGER, PARAMETER :: N_CAT=5
REAL :: Anew(5)
REAL :: Agro(5)
REAL :: Amat(5)
REAL :: Aold(5)

DATA Anew(1), Agro(1), Amat(1), Aold(1) /1.0 , 1.0, 1.0  , 1.0/
DATA Anew(2), Agro(2), Amat(2), Aold(2) /2.0 , 1.8, 0.95 , 1.0/
DATA Anew(3), Agro(3), Amat(3), Aold(3) /0.4 , 0.6, 1.075, 1.0/
DATA Anew(4), Agro(4), Amat(4), Aold(4) /3.0 , 2.6, 0.85 , 1.0/
DATA Anew(5), Agro(5), Amat(5), Aold(5) /0.05, 0.6, 1.125, 1.0/

!=======================================================================!
!  PDT_LOS_CP.EXT
!  This include file contains "production and loss within canopy"
!  factors.
!
!  MEGAN v2.0
!
!  Created by Tan 11/30/06
!=======================================================================!
INTEGER, PARAMETER :: N_RHO_SPC=20
CHARACTER*16 ::   RHO_SPC(20)
REAL        ::   RHO_FCT(20)
INTEGER     ::   RHO_MAP(20)

DATA RHO_SPC(  1), RHO_FCT(  1), RHO_MAP(  1) / 'ISOP            ', 1.0, 1 /
DATA RHO_SPC(  2), RHO_FCT(  2), RHO_MAP(  2) / 'MYRC            ', 1.0, 2 /
DATA RHO_SPC(  3), RHO_FCT(  3), RHO_MAP(  3) / 'SABI            ', 1.0, 3 /
DATA RHO_SPC(  4), RHO_FCT(  4), RHO_MAP(  4) / 'LIMO            ', 1.0, 4 /
DATA RHO_SPC(  5), RHO_FCT(  5), RHO_MAP(  5) / '3CAR            ', 1.0, 5 /
DATA RHO_SPC(  6), RHO_FCT(  6), RHO_MAP(  6) / 'OCIM            ', 1.0, 6 /
DATA RHO_SPC(  7), RHO_FCT(  7), RHO_MAP(  7) / 'BPIN            ', 1.0, 7 /
DATA RHO_SPC(  8), RHO_FCT(  8), RHO_MAP(  8) / 'APIN            ', 1.0, 8 /
DATA RHO_SPC(  9), RHO_FCT(  9), RHO_MAP(  9) / 'OMTP            ', 1.0, 9 /
DATA RHO_SPC( 10), RHO_FCT( 10), RHO_MAP( 10) / 'FARN            ', 1.0, 10/
DATA RHO_SPC( 11), RHO_FCT( 11), RHO_MAP( 11) / 'BCAR            ', 1.0, 11/
DATA RHO_SPC( 12), RHO_FCT( 12), RHO_MAP( 12) / 'OSQT            ', 1.0, 12/
DATA RHO_SPC( 13), RHO_FCT( 13), RHO_MAP( 13) / 'MBO             ', 1.0, 13/
DATA RHO_SPC( 14), RHO_FCT( 14), RHO_MAP( 14) / 'MEOH            ', 1.0, 14/
DATA RHO_SPC( 15), RHO_FCT( 15), RHO_MAP( 15) / 'ACTO            ', 1.0, 15/
DATA RHO_SPC( 16), RHO_FCT( 16), RHO_MAP( 16) / 'CH4             ', 1.0, 16/
DATA RHO_SPC( 17), RHO_FCT( 17), RHO_MAP( 17) / 'NO              ', 1.0, 17/
DATA RHO_SPC( 18), RHO_FCT( 18), RHO_MAP( 18) / 'ACTA            ', 1.0, 18/
DATA RHO_SPC( 19), RHO_FCT( 19), RHO_MAP( 19) / 'FORM            ', 1.0, 19/
DATA RHO_SPC( 20), RHO_FCT( 20), RHO_MAP( 20) / 'CO              ', 1.0, 20/

!=======================================================================!
!  SPC_MGN.EXT
!  This include file contains MEGAN species
!
!  MEGAN v2.0
!
!  Created by Tan 12/02/06
!=======================================================================!
INTEGER, PARAMETER :: N_MGN_SPC=20
CHARACTER*16 :: MGN_SPC(20)
REAL :: MGN_MWT(20)

DATA MGN_SPC(  1), MGN_MWT(  1) / 'ISOP            ', 1.0/
DATA MGN_SPC(  2), MGN_MWT(  2) / 'MYRC            ', 1.0/
DATA MGN_SPC(  3), MGN_MWT(  3) / 'SABI            ', 1.0/
DATA MGN_SPC(  4), MGN_MWT(  4) / 'LIMO            ', 1.0/
DATA MGN_SPC(  5), MGN_MWT(  5) / '3CAR            ', 1.0/
DATA MGN_SPC(  6), MGN_MWT(  6) / 'OCIM            ', 1.0/
DATA MGN_SPC(  7), MGN_MWT(  7) / 'BPIN            ', 1.0/
DATA MGN_SPC(  8), MGN_MWT(  8) / 'APIN            ', 1.0/
DATA MGN_SPC(  9), MGN_MWT(  9) / 'OMTP            ', 1.0/
DATA MGN_SPC( 10), MGN_MWT( 10) / 'FARN            ', 1.0/
DATA MGN_SPC( 11), MGN_MWT( 11) / 'BCAR            ', 1.0/
DATA MGN_SPC( 12), MGN_MWT( 12) / 'OSQT            ', 1.0/
DATA MGN_SPC( 13), MGN_MWT( 13) / 'MBO             ', 1.0/
DATA MGN_SPC( 14), MGN_MWT( 14) / 'MEOH            ', 1.0/
DATA MGN_SPC( 15), MGN_MWT( 15) / 'ACTO            ', 1.0/
DATA MGN_SPC( 16), MGN_MWT( 16) / 'CH4             ', 1.0/
DATA MGN_SPC( 17), MGN_MWT( 17) / 'NO              ', 1.0/
DATA MGN_SPC( 18), MGN_MWT( 18) / 'ACTA            ', 1.0/
DATA MGN_SPC( 19), MGN_MWT( 19) / 'FORM            ', 1.0/
DATA MGN_SPC( 20), MGN_MWT( 20) / 'CO              ', 1.0/

!=======================================================================!
!  TEMPD_PRM.EXT
!  This include file contains "temperature dependent" parameter for 
!  light-independent emissions 
!
!  This is an input file for MEGAN v2.0
!  This is version v210 of this file
!  Created by A. Guenther 8/11/07
!=======================================================================!
INTEGER, PARAMETER     ::   N_TDF_SPC=20
CHARACTER*16 ::  TDF_SPC(20)
REAL        ::   TDF_PRM(20)
INTEGER     ::   TDF_MAP(20)

DATA TDF_SPC(  1), TDF_PRM(  1), TDF_MAP(  1) / 'ISOP            ', 0.09, 1 /
DATA TDF_SPC(  2), TDF_PRM(  2), TDF_MAP(  2) / 'MYRC            ', 0.1 , 2 /
DATA TDF_SPC(  3), TDF_PRM(  3), TDF_MAP(  3) / 'SABI            ', 0.1 , 3 /
DATA TDF_SPC(  4), TDF_PRM(  4), TDF_MAP(  4) / 'LIMO            ', 0.1 , 4 /
DATA TDF_SPC(  5), TDF_PRM(  5), TDF_MAP(  5) / '3CAR            ', 0.1 , 5 /
DATA TDF_SPC(  6), TDF_PRM(  6), TDF_MAP(  6) / 'OCIM            ', 0.1 , 6 /
DATA TDF_SPC(  7), TDF_PRM(  7), TDF_MAP(  7) / 'BPIN            ', 0.1 , 7 /
DATA TDF_SPC(  8), TDF_PRM(  8), TDF_MAP(  8) / 'APIN            ', 0.1 , 8 /
DATA TDF_SPC(  9), TDF_PRM(  9), TDF_MAP(  9) / 'OMTP            ', 0.1 , 9 /
DATA TDF_SPC( 10), TDF_PRM( 10), TDF_MAP( 10) / 'FARN            ', 0.17, 10/
DATA TDF_SPC( 11), TDF_PRM( 11), TDF_MAP( 11) / 'BCAR            ', 0.17, 11/
DATA TDF_SPC( 12), TDF_PRM( 12), TDF_MAP( 12) / 'OSQT            ', 0.17, 12/
DATA TDF_SPC( 13), TDF_PRM( 13), TDF_MAP( 13) / 'MBO             ', 0.09, 13/
DATA TDF_SPC( 14), TDF_PRM( 14), TDF_MAP( 14) / 'MEOH            ', 0.08, 14/
DATA TDF_SPC( 15), TDF_PRM( 15), TDF_MAP( 15) / 'ACTO            ', 0.11, 15/
DATA TDF_SPC( 16), TDF_PRM( 16), TDF_MAP( 16) / 'CH4             ', 0.05, 16/
DATA TDF_SPC( 17), TDF_PRM( 17), TDF_MAP( 17) / 'NO              ', 0.11, 17/
DATA TDF_SPC( 18), TDF_PRM( 18), TDF_MAP( 18) / 'ACTA            ', 0.13, 18/
DATA TDF_SPC( 19), TDF_PRM( 19), TDF_MAP( 19) / 'FORM            ', 0.09, 19/
DATA TDF_SPC( 20), TDF_PRM( 20), TDF_MAP( 20) / 'CO              ', 0.09, 20/

!=======================================================================!
! Canopychar
!=======================================================================!
INTEGER, PARAMETER :: NrCha=16, NrTyp=9
REAL,DIMENSION(NrCha,NrTyp) :: Canopychar

!===== Original canopy.txt =====!
! DATA Canopychar(1 ,:) /12.  , 16.  , 16.  , 6.   /    ! Canopy depth
! DATA Canopychar(2 ,:) /0.005, 0.015, 0.05 , 0.015/    ! Leaf width
! DATA Canopychar(3 ,:) /0.08 , 0.1  , 0.1  , 0.1  /    ! Leaf length
! DATA Canopychar(4 ,:) /13.87, 24.  , 24.  , 12.  /    ! Canopy height
! DATA Canopychar(5 ,:) /0.1  , 0.2  , 0.2  , 0.2  /    ! Scattering coefficients for visibld light
! DATA Canopychar(6 ,:) /0.8  , 0.8  , 0.8  , 0.8  /    ! Scattering coefficients for near-IR light
! DATA Canopychar(7 ,:) /0.01 , 0.057, 0.057, 0.057/    ! reflection coefficients for diffused visible light
! DATA Canopychar(8 ,:) /0.2  , 0.389, 0.389, 0.389/    ! reflection coefficients for diffused near-IR light
! DATA Canopychar(9 ,:) /0.95 , 0.87 , 0.9  , 0.85 /    ! Cluster factor (clumping factor)
! DATA Canopychar(10,:) /0.5  , 0.95 , 0.95 , 0.95 /    ! Eps, short for emissivity
! DATA Canopychar(11,:) /1.25 , 1.25 , 1.25 , 1.25 /    ! Transpire type
! DATA Canopychar(12,:) /0.08 , 0.06 , 0.06 , 0.06 /    ! Stability (Daytime temperature lapse rate)
! DATA Canopychar(13,:) /-0.06, -0.06, -0.06, -0.06/    ! Stability (Nighttime temperature lapse rate)
! DATA Canopychar(14,:) /800. , 700. , 700. , 700. /    !
! DATA Canopychar(15,:) /150. , 150. , 150. , 150. /    !
! DATA Canopychar(16,:) /0.7  , 0.7  , 0.7  , 0.7  /    !

!===== Ditte's canopy.txt =====!
DATA Canopychar(1 ,:) /8.76  , 16.  , 16.  , 6.   , 1.   , 0.5  , 1.   , 16.  , 9.    /  ! Canopy depth: 9.0 - it is really 8.76
DATA Canopychar(2 ,:) /0.0012, 0.015, 0.03 , 0.015, 0.015, 0.01 , 0.02 , 0.05 , 0.0005/  ! [m], Leaf width
DATA Canopychar(3 ,:) /0.065 , 0.1  , 0.055, 0.1  , 0.1  , 0.15 , 0.15 , 0.1  , 0.065 /  ! [m], Leaf length
DATA Canopychar(4 ,:) /18.53 , 24.  , 24.  , 12.  , 2.   , 0.5  , 1.   , 24.  , 18.0  /  ! Canopy height: 18.0
DATA Canopychar(5 ,:) /0.1   , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.2  , 0.1   /  !  Scattering coefficient for PPFD
DATA Canopychar(6 ,:) /0.8   , 0.8  , 0.8  , 0.8  , 0.8  , 0.8  , 0.8  , 0.8  , 0.8   /  !  Scattering coefficient for near IR
DATA Canopychar(7 ,:) /0.01  , 0.057, 0.057, 0.057, 0.057, 0.057, 0.057, 0.057, 0.01  /  !  Reflection coefficient for diffuse PPFD
DATA Canopychar(8 ,:) /0.2   , 0.389, 0.389, 0.389, 0.389, 0.389, 0.389, 0.389, 0.2   /  !  Reflection coefficient for diffuse near IR
DATA Canopychar(9 ,:) /0.95  , 0.87 , 0.9  , 0.85 , 0.85 , 0.7  , 0.7  , 0.6  , 0.95  /  ! Clustering coefficient (accounts for leaf clumping
                                                                                         ! influence on mean projected leaf area
                                                                                         ! in the direction of the suns beam)
                                                                                         ! use 0.85 for default, corn=0.4-0.9; Pine=0.6-1.0;
                                                                                         ! oak=0.53-0.67; tropical rainforest=1.1
DATA Canopychar(10,:) /0.95  , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.5   /  ! Leaf IR emissivity
DATA Canopychar(11,:) /1.25  , 1.25 , 1.25 , 1.25 , 1    , 1.25 , 1.25 , 1.25 , 1.25  /  ! Leaf stomata and cuticle factor:
                                                                                         ! 1=hypostomatous, 2=amphistomatous, 1.25=hypostomatous
                                                                                         ! but with some transpiration through cuticle
DATA Canopychar(12,:) /0.08  , 0.06 , 0.06 , 0.06 , 0.06 , 0.06 , 0.06 , 0.06 , 0.08  /  ! Daytime temperature lapse rate (K m-1)
DATA Canopychar(13,:) /-0.06 , -0.06, -0.06, -0.06, -0.06, -0.06, -0.06, -0.06, -0.06 /  ! Nighttime temperature lapse rate (K m-1)
DATA Canopychar(14,:) /800.  , 700. , 700. , 700. , 700. , 700  , 700  , 700  , 800   /  ! Warm (>283K) canopy total humidity change (Pa)
DATA Canopychar(15,:) /150.  , 150. , 150. , 150. , 150. , 150  , 150  , 150  , 150   /  ! Cool (>= 283K) canopy total humidity change (Pa)
DATA Canopychar(16,:) /0.7   , 0.7  , 0.7  , 0.7  , 0.7  , 0.7  , 0.7  , 0.7  , 0.7   /  ! Normalized canopy depth where wind is negligible


!=======================================================================!
! Other parameters
!=======================================================================!
REAL, PARAMETER :: ug2tonne = 1e-12    ! convert microgram to metric tonne
REAL, PARAMETER :: hr2sec = 3600       ! convert hr to second
REAL, PARAMETER :: ug2g = 1e-6         ! convert microgram to gram
REAL, PARAMETER :: Avog = 6.0221e23    ! Avogadro-number [molecules mol-1]
INTEGER, PARAMETER :: NEMIS = N_MGN_SPC    ! Number of output emission variable and number of MEGAN species
INTEGER, PARAMETER :: NPFT = NrTyp    ! Number of plant functional types
REAL, PARAMETER :: leaf_biomass_hyy = 0.0538    ! [g cm-2] = 0.538 [kg m-2]
REAL, PARAMETER :: nggh2gcm2s = 1e-9 / 3600. * leaf_biomass_hyy    ! [ng g-1 h-1] to [g cm-2 s-1]
REAL, PARAMETER :: Sb = 5.67e-8        ! Stefan-Boltzman constant (W/m2/K4)
REAL, PARAMETER :: ConvertWm2toUmolm2s = 4.766, &    ! Convert radiation from [W/m2] to [umol/m2/s1]
                   SolarConstant       = 1367., &    ! Solar constant [W/m2]
                   Td                  = 301.,  &    ! Temperature used in Ea1t99
                   Cce                 = 0.57
REAL, PARAMETER :: WaterAirRatio = 18.016 / 28.97    ! Ratio between water and air molecules
REAL, PARAMETER :: ConvertPPFD = 4.766    ! Solar raidation to PPFD

END MODULE M2_Parameter
