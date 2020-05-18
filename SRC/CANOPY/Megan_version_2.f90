MODULE Megan_version_2
!***********************************************************************
!
!   CALL:
!      MODULE GAMMA_ETC
!         GAMMA_LAI
!         GAMMA_P
!         GAMMA_TISOP
!         GAMMA_TNISP
!         GAMMA_A
!         GAMMA_S
!
!   Created by Jack Chen 11/04
!   Modified by Tan 11/21/06 for MEGAN v2.0
!
!   History:
!
!       Jun, 2010  modifications to I/O and general structure - Anton Rusanen
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Scientific algorithm
!
!             Emission = [EF] [GAMMA] [RHO]
!                        where
!                        [EF]    = emission factor (ug/m2h)
!                        [GAMMA] = emission activity factor (non-dimension)
!                        [RHO]   = production and loss within plant canopies (non-dimensino)
!                        Assumption: [RHO] = 1  (11/27/06)   (See PDT_LOT_CP.EXT)
!
!             GAMMA    = [GAMMA_CE] [GAMMA_age] [GAMMA_SM]
!                        where
!                        [GAMMA_CE]  = canopy correction factor
!                        [GAMMA_age] = leaf age correction factor
!                        [GAMMA_SM]  = soil moisture correction factor
!                        Assumption: [GAMMA_SM]  = 1  (11/27/06)
!             
!             GAMMA_CE = [GAMMA_LAI] [GAMMA_P] [GAMMA_T]    # this is the PCEEA algorithm
!                        where
!                        [GAMMA_LAI] = leaf area index factor
!                        [GAMMA_P]   = PPFD emission activity factor
!                        [GAMMA_T]   = temperature response factor
!
!             Emission = [EF] [GAMMA_LAI] [GAMMA_P] [GAMMA_T] [GAMMA_age] [GAMMA_SM] [RHO]
!                        Derivation:
!                        Emission = [EF] [GAMMA_etc] (1-LDF) + [EF] [GAMMA_etc] [LDF] [GAMMA_P]    # [GAMMA_SM] = 1, [RHO] = 1
!                        Emission = [EF] [GAMMA_etc] {(1-LDF) + [LDF] [GAMMA_P]}
!                        where LDF = light dependent function (non-dimension)
!                                    (See LD_FCT.EXT)
!
!             Final Equation
!             
!             Emission = [EF] [GAMMA_LAI] [GAMMA_T] [GAMMA_age] [GAMMA_SM] [RHO] * { (1-LDF) + [LDF] [GAMMA_P] }
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE M2_Parameter
  USE M2_GAMMA_ETC           ! Module containing gamma functions
  USE M2_Canopy
  
  !How does this thing use the functions from canopy?

  IMPLICIT NONE

  !..  Constants
  ! Length of the time step (days) 
  ! if the model produces odd values the reason is probably this and/or the removal of the time loop + the layer things I did ( will have done, in the future) in canopy.F

  INTEGER, PARAMETER   :: TSTLEN = 1

  ! parameter for unit conversion

  LOGICAL, PARAMETER :: TONPHR = .FALSE.  ! output in tons/hr flag

  !..  Program I/O files

  CONTAINS 

  SUBROUTINE EMISSION_M2(Day, Layersin, Beta, LATin, LONGin, DATEin, TIMEin, PPFDin, D_PPFDin, &
                         TEMPin, DTEMPin, PRESin, HUMin, WINDin, SMOISTin, LADpin, LAIcin, LAIpin, &
                         ERout, VARout, ER_HB_out, ER_SB_out, ER_NT_out, ER_BT_out, GAM_PHO_out, GAM_CE_out, GAM_T_out, &
                         GAM_TMP_out, GAM_OTHER_out, Zin, kzin, EMI, sun_par, Sunleaftk, Shadeleaftk, sunfrac, &
                         SH_out, LH_out, &
                         Rbs, Rds, RLs)

    INTEGER, INTENT(IN) ::  DAY                   ! Day from YYYYDDD
    REAL, INTENT(IN)    ::  Beta                  ! Solar elevation angle
    INTEGER, INTENT(IN) ::  Layersin ,kzin        ! Number of Layers with vegetation in them, presumed from ground 0 to Layersin.
                                                  ! kzin is usually equal to (Layersin+1)
    INTEGER, INTENT(IN) ::  DATEin, TIMEin        ! Date (yyyyddd) and time (hhmmss)


    REAL, SAVE          ::  EFin(366,22)          ! emission factors, 1-dimensional array, 20 slots by default
    REAL, INTENT(IN)    ::  LATin, LONGin         ! LATitutde and LONGitude
    REAL, INTENT(IN)    ::  DTEMPin               ! TEMPerature (K) and Daily average TEMPerature (K)
    REAL, INTENT(IN)    ::  PPFDin, D_PPFDin      ! PAR and Daily average of PAR in (umol/m2/s)
    REAL, INTENT(IN)    ::  PRESin                ! Pressure (Pa)
    REAL, INTENT(IN)    ::  LAIcin , LAIpin       ! LAIc is LAI for current month, LAIp for previous
    REAL, INTENT(IN)    ::  SMOISTin              ! [m3 m-3], Soil moisture, (a number belonging to [0,1])

    REAL, INTENT(IN), DIMENSION(kzin) ::  LADpin  ! Leaf area density as [0,1] in the canopy
    REAL, INTENT(IN), DIMENSION(kzin) ::  TEMPin  ! TEMPerature [K]
    REAL, INTENT(IN), DIMENSION(kzin) ::  HUMin   ! Mixing ratio of water vapor
    REAL, INTENT(IN), DIMENSION(kzin) ::  Zin     ! Layer top array

    REAL, INTENT(IN), DIMENSION(kzin) :: WINDin   ! Wind speed (m/s)

    REAL, INTENT(INOUT) :: ERout(N_MGN_SPC, kzin) ! Emission rates output (normally g/s, can be set with a boolean to ton/hr)

    CHARACTER*16       VARout(N_MGN_SPC)          !VOC names in the same order as ERout

    REAL, DIMENSION(N_MGN_SPC, kzin), INTENT(OUT) ::  ER_HB_out, ER_SB_out, ER_NT_out, ER_BT_out  ! emissions by plant functional type
    REAL, DIMENSION(4, kzin),         INTENT(OUT) ::  GAM_PHO_out, GAM_CE_out, GAM_T_out          ! Gamma factors for photons, CE, and temperature
    REAL, DIMENSION(N_MGN_SPC, kzin), INTENT(OUT) ::  GAM_TMP_out                                 ! Gamma factor for what?
    REAL, DIMENSION(N_MGN_SPC, 3),    INTENT(OUT) ::  GAM_OTHER_out                               ! Gamma factors for soil moisture, leaf age and lai correction
    REAL, DIMENSION(kzin),            INTENT(OUT) ::  SH_out, LH_out                              ! Heat fluxes from the canopy

    REAL,DIMENSION(kzin) :: sun_par, Sunleaftk, Shadeleaftk, sunfrac

    REAL, DIMENSION(Layersin) :: midpoints ! Layer centerpoints for GAMMA_CE

    REAL, DIMENSION(N_MGN_SPC, kzin) :: ERtemporary

    REAL :: ISO_EMI, MYR_EMI, SAB_EMI, LIM_EMI, CAR_EMI, &
            OCI_EMI, BPI_EMI, API_EMI, OMT_EMI, FAR_EMI, &
            BCA_EMI, OSQ_EMI, MBO_EMI, MET_EMI, ACT_EMI, &
            CH4_EMI, NO_EMI , ACA_EMI, FOR_EMI, CO_EMI , &
            CIN_EMI, LIN_EMI, &
            emi_factor, time

    REAL, DIMENSION(kzin,22) :: EMI

    ! Program name
    CHARACTER*16  :: PROGNAME = 'MEGAN'
     
    !..  Internal parameters
    ! internal paramters (status and buffer)
    INTEGER       IOS                              ! i/o status
    CHARACTER*256 MESG                             ! message buffer
    
    ! parameter for output species

    ! local variables and their descriptions:
    REAL          GAREA                            ! Area in one grid (meter^2)

    REAL          LDF                              ! Light dependent factor
    REAL          RHO                              ! Production and loss within canopy

    REAL :: EF              ! input annual emission factor

    REAL, ALLOCATABLE    :: GAM_PHO(:)       ! light correction factor
    REAL, ALLOCATABLE    :: GAM_TMP(:)             ! temperature correction factor
    REAL :: GAM_LHT         ! LAI correction factor
    REAL :: GAM_AGE         ! leaf age correction factor
    REAL :: GAM_SMT         ! Soil moilture correction factor
    REAL, ALLOCATABLE    :: GAM_T(:)
    REAL, ALLOCATABLE    :: GAM_CE (:)
    REAL                 :: DI
    REAL, ALLOCATABLE    :: GAMMAfactors(:,:,:)  ! gamma factors (  GAMMA_T, GAMMA_PHO, GAMMA_CE), In layers, for different plant types 

    INTEGER :: ISLTYP

    INTEGER          ILINE                         ! current line
    CHARACTER(LEN=1000) LINE                       ! input line buffer
    INTEGER       CID, INX, INY                    ! Input grid x and y
    INTEGER, PARAMETER :: MXTCOL = 9               ! Columns in an input line
    CHARACTER*30     SEGMENT( MXTCOL )             ! Input line fields

    INTEGER, PARAMETER :: NVARS = MXTCOL - 3       ! Number of LATin, LONGin, and
                                                   ! PFT factor variables
    CHARACTER*16 VNAME( NVARS )                    ! Variable names
    REAL :: PFTF(NPFT)                   ! PFT factor array

    CHARACTER*16   VARIABLENAMES(30)

!   loop indices
    INTEGER       T, S, I, J ,K, VAR, VAR2, lay     ! Counters
    INTEGER       NMAP            ! Index

    ! times
    INTEGER       MON             ! Month from YYYYDDD

    INTEGER :: test1,test2,test3

    ! Added by Zhou Putian
    LOGICAL, SAVE :: first_time = .TRUE. ! reading 'EF_day.txt' and 'canopy.txt' only once
    REAL :: temp_parameter ! used for simplify code
    REAL :: dz(Layersin) ! Zin(lay) - Zin(lay-1) for lay>1, Zin(1) for lay==1

    REAL :: Rbs, Rds, RLs

    ! months
    CHARACTER*3   MONTHS( 12 )
    DATA          MONTHS &
    &  / 'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , &
    &    'JUL' , 'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'   /
    CHARACTER*2   MONNUM( 12 )
    DATA          MONNUM &
    &  /  '1 ' ,  '2 ' ,  '3 ' ,  '4 ' ,  '5 ' ,  '6 ' , &
    &     '7 ' ,  '8 ' ,  '9 ' ,  '10' ,  '11' ,  '12'   /

    ! #ifdef PARALLEL
    !     CHARACTER(LEN=*), PARAMETER :: &
    !          filename1 = 'sosa_in', &
    !          filename2 = 'sosa_out'
    ! #elif LINUX
    !    CHARACTER(LEN=*), PARAMETER :: &
    !         filename1 = 'sosa_in', &
    !         filename2 = 'sosa_out'
    !     CHARACTER(LEN=*), PARAMETER :: &
    !          filename1 = './sosa_in', &
    !          filename2 = './sosa_out'
    ! #else
    !     CHARACTER(LEN=*), PARAMETER :: &
    !          filename1 = './sosa_in', &
    !          filename2 = './sosa_out'
    ! #endif

!**********************************************************************

!======================================================================
!..  Begin program
!======================================================================
!----------------------------------------------------------------------
!....1) File set up and assign I/O parameters
!----------------------------------------------------------------------
    GAREA = 1. ! Calculate everthing for 1 [m2], then multiply the area in the explicit models

    !..  Get input parameters

    DO S = 1, NEMIS
      VARIABLENAMES(S) = TRIM(MGN_SPC(S))
    END DO
    VARIABLENAMES(NEMIS + 1) = 'D_TEMP'
    VARIABLENAMES(NEMIS + 2) = 'D_PPFD'
!----------------------------------------------------------------------
!....2) Process emission rates
!----------------------------------------------------------------------
!..  Allocate memory
    ALLOCATE ( GAM_PHO(NPFT), STAT = IOS )
    ALLOCATE ( GAM_TMP(Layersin), STAT = IOS )
    ALLOCATE ( GAM_T(NPFT), STAT = IOS )      
    ALLOCATE ( GAM_CE (NPFT), STAT = IOS )

    ALLOCATE (GAMMAfactors(3, NPFT, Layersin))

    !.. Partial fractions of different plant types 
    !.. (doesn't need to sum to 100%, because parts of the grid can be barren)
    !.. You should hardcode appropriate values for your model here, or this can be added to input parameters

    PFTF = 0.
    PFTF(1) = 100.  ! All are Scots pines
    PFTF(2) = 100.  ! All are Scots pines

    !.. Calculate daily PPDF
    !  PPFD: SRAD - short wave from sun (W/m2)
    !  assuming 4.766 (umol m-2 s-1) per (W m-2)
    !  assume 1/2 of SRAD is in 400-700nm band
    !  PPFD = SRAD * 4.766 * 0.5

!------INPUT OTHER VARIABLES----------------
      
    DI     = 0.0       ! something related to DIstomata, value 0 came with model
    ISLTYP = 6         ! Soiltype (I don't know if 6 is right, in fact I don't know what type of soil it is, value came with model)
    
    ! write(*,*) 'Megan_version_2: Zin, Layersin', Zin, Layersin
    midpoints= centralizer(Zin, Layersin)
    !rescaling to [0,1]
    midpoints= midpoints/(Zin(Layersin))
    ! calculate dz for each layer, added by Zhou Putian
    dz(2:Layersin) = Zin(2:Layersin) - Zin(1:Layersin-1)
    dz(1) = Zin(1) ! suppose the ground is 0m
    ! write(*,*) 'Megan_version_2: midpoints, dz, ', midpoints, dz

    !--------INPUT standard emission potential
    ! 1 =  Isoprene
    ! 2 =  MBO (2methyl-3buten-2ol)
    ! 3 =  MYRC
    ! 4 =  Sabinene
    ! 5 =  Limonen
    ! 6 =  3-Carene
    ! 7 =  Ocimene
    ! 8 =  Beta-pinene
    ! 9 =  Alpha-pinene
    ! 10 = FARN
    ! 11 = Betacarophylene
    ! 12 = Methanol
    ! 13 = Aceton
    ! 14 = Acetaldehyde
    ! 15 = Formaldehyde
    ! 16 = Methan
    ! 17 = NO
    ! 18 = Other monoterpene
    ! 19 = Other sesquiterpenes
    ! 20 = CO
    ! 21 = Cineole
    ! 22 = Linalool
    IF (first_time) THEN
      OPEN(unit=5,file='EF_day.txt')
      DO J = 1,366
        READ(5,*) (EFin(J,I), I = 1,22)
      END DO
      CLOSE(5)

      first_time = .FALSE.
    END IF

    !..Go over all the chemical species
    DO S = 1, NEMIS
      !..Initialize variables
      EF = 0.
      GAM_PHO = 0.
      GAM_TMP = 0.
      GAM_LHT = 0.
      GAM_AGE = 0.
      GAM_SMT = 0.
      GAM_T = 0.
      GAM_CE = 0.

      !..Get EF
      EF = 1. ! EFin(S)

      GAMMAfactors = -1.0

      !..Select algorithms for differet chemical species
      !  Due to the differences in the calculation for gamma,
      !  this logical condition will check for species and choose
      !  (call) the right algorithms to calculate the gamma and
      !  return back the main.
      IF (S == 1) THEN
        CALL GAMMA_CE(DATEin, TIMEin, Beta, LATin, LONGin, &    ! input
                      PPFDin, TEMPin, LAIcin, WINDin, PRESin, HUMin, &    ! input
                      DI, PFTF, LADpin, midpoints, Layersin, &    ! input
                      GAM_T, GAM_PHO, GAM_CE, GAMMAfactors, sun_par, &    ! output
                      Sunleaftk, Shadeleaftk, sunfrac, SH_out, LH_out, &    ! output
                      Rbs, Rds, RLs)    ! output
      END IF
      SELECT CASE ( VARIABLENAMES(S) )
        CASE ( 'ISOP' )
        ! NOTE isoprene is not getting calculated. I don't know enough about isoprene emissions to fix this. I could just uncomment
        ! those calulations, but I can't tell if the results will be scientifically valid. The calculation was commented out when I got this model.
        !   CALL GAMMA_TISOP( TEMPin(1:Layersin), DTEMPin, Layersin, GAM_T )
        CASE ('MYRC', 'SABI', 'LIMO', '3CAR', 'OCIM', 'BPIN', 'APIN', 'OMTP',                     &
              'FARN', 'BCAR', 'OSQT', 'MBO', 'MEOH', 'ACTO', 'CH4', 'NO', 'ACTA', 'FORM', 'CO')
          CALL GAMMA_TNISP(VARIABLENAMES(S), TEMPin(1:Layersin), Layersin, GAM_TMP)
        CASE DEFAULT
          MESG = 'Error: Chemical species, invalid variable: ' // TRIM(VARIABLENAMES(S))
          WRITE(*,*) MESG
          STOP
      END SELECT

      CALL GAMMA_LAI(LAIcin, GAM_LHT)
      ! CALL GAMMA_P(DATEin, TIMEin, LATin, LONGin, PPFDin, D_PPFDin, GAM_PHO)
      CALL GAMMA_A(DATEin, TIMEin, VARIABLENAMES(S), LAIpin, LAIcin, TSTLEN, DTEMPin, GAM_AGE)
      ! CALL GAMMA_S( GAM_SMT )
      CALL GAMMA_S(SMOISTin, ISLTYP, GAM_SMT) ! ISLTYP, soil type, 6 here

      DO VAR = 1, NEMIS
        IF (TRIM(VARIABLENAMES(S)) .EQ. TRIM(LDF_SPC(VAR))) THEN
          LDF = LDF_FCT(VAR)
        END IF

        IF (TRIM(VARIABLENAMES(S)) .EQ. TRIM(RHO_SPC(VAR))) THEN
          RHO = RHO_FCT(VAR)
        END IF
      END DO

      !Same thing but separately for every layer.
      SELECT CASE ( VARIABLENAMES(S) )
      CASE ('ISOP')
        DO lay=1, Layersin
          IF (TONPHR) THEN
            ER_BT_out(S,lay) = (GAMMAfactors(3,2,lay) * GAM_AGE * GAM_SMT * RHO) * GAREA * ug2tonne
            ER_NT_out(S,lay) = (GAMMAfactors(3,1,lay) * GAM_AGE * GAM_SMT * RHO) * GAREA * ug2tonne
            ER_SB_out(S,lay) = (GAMMAfactors(3,3,lay) * GAM_AGE * GAM_SMT * RHO) * GAREA * ug2tonne
            ER_HB_out(S,lay) = (GAMMAfactors(3,4,lay) * GAM_AGE * GAM_SMT * RHO) * GAREA * ug2tonne
          ELSE
            ER_BT_out(S,lay) = (GAMMAfactors(3,2,lay) * GAM_AGE * GAM_SMT * RHO)
            ER_NT_out(S,lay) = (GAMMAfactors(3,1,lay) * GAM_AGE * GAM_SMT * RHO)
            ER_SB_out(S,lay) = (GAMMAfactors(3,3,lay) * GAM_AGE * GAM_SMT * RHO)
            ER_HB_out(S,lay) = (GAMMAfactors(3,4,lay) * GAM_AGE * GAM_SMT * RHO)
          END IF
        END DO
      CASE ('MYRC', 'SABI', 'LIMO', '3CAR', 'OCIM', 'BPIN', 'APIN', 'OMTP',                     &
            'FARN', 'BCAR', 'OSQT', 'MBO', 'MEOH', 'ACTO', 'CH4', 'NO', 'ACTA', 'FORM', 'CO')
        DO lay=1, Layersin
          IF ( TONPHR ) THEN
            ER_BT_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,2,lay)*LDF)) * GAREA * ug2tonne
            ER_NT_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,1,lay)*LDF)) * GAREA * ug2tonne
            ER_SB_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,3,lay)*LDF)) * GAREA * ug2tonne
            ER_HB_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,4,lay)*LDF)) * GAREA * ug2tonne
          ELSE
            ER_BT_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,2,lay)*LDF))
            ER_NT_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,1,lay)*LDF))
            ER_SB_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,3,lay)*LDF))
            ER_HB_out(S,lay) = (EF * GAM_TMP(lay) * GAM_AGE * GAM_LHT * GAM_SMT * RHO) * &
                               ((1-LDF) + (GAMMAfactors(2,4,lay)*LDF))
          ENDIF
        END DO
      CASE DEFAULT
        WRITE(*,*) 'Error: Chemical species, invalid variable: '
        WRITE(*,*) TRIM(VARIABLENAMES(S))
        STOP
      END SELECT

!----------------------------------------------------------------------
!....3) Write out the calculated ER and met data
!----------------------------------------------------------------------
      ! Write emission to file

      IF (S <= N_MGN_SPC) THEN
        VARout(S) = VARIABLENAMES( S)
      END IF

      IF (S .EQ. 9) THEN ! OMTP
        DO lay = 1, 1 !looping over plant functional types
          DO J = 1,layersin
            ! GAM_T_out(lay, J)     = GAMMAfactors(1, lay, (layersin-J+1)) ! Note GAM_T is not actually used anywhere
            ! GAM_PHO_out(lay, J)   = GAMMAfactors(2, lay, (layersin-J+1))
            ! GAM_CE_out(lay, J)    = GAMMAfactors(3, lay, (layersin-J+1))
            GAM_T_out(lay, J)     = GAMMAfactors(1, lay, J) ! Note GAM_T is not actually used anywhere
            GAM_PHO_out(lay, J)   = GAMMAfactors(2, lay, J)
            GAM_CE_out(lay, J)    = GAMMAfactors(3, lay, J)
          ENDDO
        ENDDO
      ENDIF

      GAM_OTHER_out(S, 1) = GAM_AGE
      GAM_OTHER_out(S, 2) = GAM_LHT
      GAM_OTHER_out(S, 3) = GAM_SMT
      
      GAM_TMP_out(S, 1:Layersin) = GAM_TMP(1:Layersin)
    END DO ! End loop for emission species (S)

    !ERtemporary = ER_HB_out + ER_NT_out + ER_SB_out + ER_BT_out
    ERtemporary = ER_NT_out
    
    !Erout= Ertemporary
    ERout(:, 1:Layersin) = Ertemporary(:, 1:Layersin)
    !Rout = Ertemporary

    ! Standard emission potential [ng/g(needledryweight)/h into g/cm2/s] for Hyytiala spring from Hakola et al. Biogeosc., 3, 2006
    ! (with g(needledryweight into cm2 by CarboEuroflux data: Standing leaf biomass Hyytiala 0.538 kg/m2 or 0.0538 g/cm2)

    ISO_EMI = EFin(Day,1)  * nggh2gcm2s
    MYR_EMI = EFin(Day,2)  * nggh2gcm2s
    SAB_EMI = EFin(Day,3)  * nggh2gcm2s
    LIM_EMI = EFin(Day,4)  * nggh2gcm2s
    CAR_EMI = EFin(Day,5)  * nggh2gcm2s
    OCI_EMI = EFin(Day,6)  * nggh2gcm2s
    BPI_EMI = EFin(Day,7)  * nggh2gcm2s
    API_EMI = EFin(Day,8)  * nggh2gcm2s
    OMT_EMI = EFin(Day,9)  * nggh2gcm2s
    FAR_EMI = EFin(Day,10) * nggh2gcm2s
    BCA_EMI = EFin(Day,11) * nggh2gcm2s
    OSQ_EMI = EFin(Day,12) * nggh2gcm2s
    MBO_EMI = EFin(Day,13) * nggh2gcm2s
    MET_EMI = EFin(Day,14) * nggh2gcm2s
    ACT_EMI = EFin(Day,15) * nggh2gcm2s
    CH4_EMI = EFin(Day,16) * nggh2gcm2s
    NO_EMI  = EFin(Day,17) * nggh2gcm2s
    ACA_EMI = EFin(Day,18) * nggh2gcm2s
    FOR_EMI = EFin(Day,19) * nggh2gcm2s
    CO_EMI  = EFin(Day,20) * nggh2gcm2s
    CIN_EMI = EFin(Day,21) * nggh2gcm2s
    LIN_EMI = EFin(Day,22) * nggh2gcm2s

    ! Emission with EF in g/cm2/s with ER into molecules/cm3/s
    DO lay = 1,Layersin
      temp_parameter = Avog / dz(lay) / 100. * LADpin(lay)

      EMI(lay,1)  = ERout(1,lay)  * ISO_EMI /  68 * temp_parameter
      EMI(lay,2)  = ERout(2,lay)  * MYR_EMI / 136 * temp_parameter
      EMI(lay,3)  = ERout(3,lay)  * SAB_EMI / 136 * temp_parameter
      EMI(lay,4)  = ERout(4,lay)  * LIM_EMI / 136 * temp_parameter
      EMI(lay,5)  = ERout(5,lay)  * CAR_EMI / 136 * temp_parameter
      EMI(lay,6)  = ERout(6,lay)  * OCI_EMI / 136 * temp_parameter
      EMI(lay,7)  = ERout(7,lay)  * BPI_EMI / 136 * temp_parameter
      EMI(lay,8)  = ERout(8,lay)  * API_EMI / 136 * temp_parameter
      EMI(lay,9)  = ERout(9,lay)  * OMT_EMI / 136 * temp_parameter
      EMI(lay,10) = ERout(10,lay) * FAR_EMI / 204 * temp_parameter
      EMI(lay,11) = ERout(11,lay) * BCA_EMI / 204 * temp_parameter
      EMI(lay,12) = ERout(12,lay) * OSQ_EMI / 204 * temp_parameter
      EMI(lay,13) = ERout(13,lay) * MBO_EMI /  86 * temp_parameter
      EMI(lay,14) = ERout(14,lay) * MET_EMI /  32 * temp_parameter
      EMI(lay,15) = ERout(15,lay) * ACT_EMI /  58 * temp_parameter
      EMI(lay,16) = ERout(16,lay) * CH4_EMI /  16 * temp_parameter
      EMI(lay,17) = ERout(17,lay) * NO_EMI  /  30 * temp_parameter
      EMI(lay,18) = ERout(18,lay) * ACA_EMI /  44 * temp_parameter
      EMI(lay,19) = ERout(19,lay) * FOR_EMI /  30 * temp_parameter
      EMI(lay,20) = ERout(20,lay) * CO_EMI  /  28 * temp_parameter
      EMI(lay,21) = 0.            * CIN_EMI / 154 * temp_parameter
      EMI(lay,22) = 0.            * LIN_EMI / 154 * temp_parameter
    END DO

    sun_par(Layersin+1) = 1.
    EMI(Layersin+1, :) = 0.

    DEALLOCATE ( GAM_PHO )   ! light correction factor
    DEALLOCATE ( GAM_TMP )   ! temperature correction factor
    DEALLOCATE ( GAM_T )

    DEALLOCATE ( GAM_CE )
    DEALLOCATE ( GAMMAfactors)
!======================================================================
!..  End subroutine
!======================================================================
  END SUBROUTINE EMISSION_M2


!======================================================================
!..  Some added functionality, can be moved to a new module
!======================================================================
Pure function Timeformatter(input) result(output)

  IMPLICIT NONE

  REAL, INTENT(IN) :: input
  INTEGER :: output

  REAL :: tmp
  INTEGER :: h, m ,s

  tmp= input
  output=0

  h   = INT(tmp/3600) ! truncation of values should work as intented
  tmp = tmp - (h*3600)
  m   = INT(tmp/60)
  tmp = tmp - (m*60)
  s   = INT(tmp)

  output = (10000*h) + (100*m) + s

END function timeformatter


Pure function centralizer(Inputheightarray, Inputlayers) result(output)

  IMPLICIT NONE

  REAL, DIMENSION(:), INTENT(IN) :: Inputheightarray
  INTEGER, INTENT(IN) :: Inputlayers
  REAL, DIMENSION(Inputlayers) :: output

  INTEGER :: ind
  REAL :: temp

  output= Inputheightarray(1:Inputlayers)


  DO ind=1, Inputlayers
     temp=0

    IF(ind == 1) THEN
      temp= Inputheightarray(ind)
    ELSE
      temp = Inputheightarray(ind) - Inputheightarray(ind-1)
    END IF

    output(ind) = output(ind) - (temp*0.5)

  END DO

END function centralizer

END MODULE Megan_version_2
