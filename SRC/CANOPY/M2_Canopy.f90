!   By Alex Guenther --April 2005
!   Modified by Xuemei -July 2007  
!
! Modified for different IO
!================================================================================================
!
!   Input and output files have to be selected before starting the program
!
!================================================================================================
!!
!   Input varibles
!
!   Day                                         Julian day
!   Lat                                         Latitude
!   Hour                                        Hour of the day
!   Tc                                          Temperature [C]
!   PPFD                                        Incoming photosynthetic active radiation [umol/m2/s1]
!   Wind                                        Wind speed [m s-1]
!   Humidity                                    Mixing ratio of water vapor, [g kg-1]
!   Cantypye                                    Defines set of canopy characteristics
!   LAI                                         Leaf area index [m2 per m2 ground area]
!   DI                                          ???
!   Pres                                        Pressure [Pa]
!
!
!   Used variables:
!
!   Solar                                       Solar radiation [W/m2]
!   Maxsolar                                    Maximum of solar radiation
!   Beta                                        Solar angle above horizon
!   Sinbeta                                     Sine of Solar angle above horizon
!   TairK0                                      Above canopy air temperature [K]
!   TairK                                       Array of canopy air temperature [K]
!   Ws0                                         Above canopy wind speed [m/s]
!   Ws                                          Array of canopy wind speed [m/s]
!   HumidairPa0                                 Above canopy ambient humidity [Pa]
!   HumidairPa                                  Array of canopy ambient humidity in [Pa]
!   StomataDI                                   An index indicating water status of leaves. used to modify stomatal conductance
!   Transmis                                    Transmission of PPFD that is diffuse
!   Difffrac                                    Fraction of PPFD that is diffuse
!   PPFDfrac                                    Fraction of solar rad that is PPFD
!   Trate                                       Stability of boundary ???
!   SH                                          Sensible heat flux ???
!   VPgausWt                                    Array of gaussian weighting factors
!   VPgausDis                                   Array of gaussian weighting factors
!   VPslwWT                                     Array of gaussian weighting factors, SLW is the specific leaf weight?
!   SunFrac                                     Array of the fraction of sun leaves. i = 1 is the top canopy layer, 2 is the next layer, etc.
!   SunPPFD                                     Array of incoming (NOT absorbed) PPFD on a sun leaf [umol/m2/s]
!   ShadePPFD                                   Array of incoming (NOT absorbed) PPFD on a shade leaf [umol/m2/s]
!   SunQv                                       Array of visible radiation (in and out) fluxes on sun leaves
!   ShadeQv                                     Array of absorbed visible radiation (in and out) fluxes on shade leaves
!   SunQn                                       Array of absorbed near IR radiation (in and out) fluxes on sun leaves
!   ShadeQn                                     Array of absorbed near IR radiation (in and out) fluxes on shade leaves
!   SunleafTK                                   Array of leaf temperature for sun leaves [K]
!   SunleafSH                                   Array of sensible heat flux for sun leaves [W/m2]
!   SunleafLH                                   Array of latent heat flux for sun leaves [W/m2]
!   SunleafIR                                   Array of infrared flux for sun leaves [W/m2]
!   ShadeleafTK                                 Array of leaf temperature for shade leaves [K]
!   ShadeleafSH                                 Array of sensible heat flux for shade leaves [W/m2]
!   ShadeleafLH                                 Array of latent heat flux for shade leaves [W/m2]
!   ShadeleafIR                                 Array of infrared flux for shade leaves [W/m2]
!   QbAbsV, QbAbsN                              Absorbed direct beam light for visible and near infra red
!   QdAbsV, QdAbsN                              Array of absorbed diffuse light for visible and near infra red
!   QsAbsV, QsAbsN                              Array of absorbed scattered light for visible and near infra red
!   QBeamV, QBeamN                              Above canopy beam (direct) light for visible and near infra red
!   QDiffV, QDiffN                              Above canopy diffuse light for visible and near infra red
!   Ea1pLayer                                   Array of emission activity of light per layer
!   Ea1tLayer                                   Array of emission activity of temperature per layer
!   Ea1Layer                                    Array of companied emission activity
!   Ea1pCanopy                                  Total emission activity of light
!   Ea1tCanopy                                  Total emission activity of temperature
!   Ea1Canopy                                   Total companied emission activity
!
!   Calcbeta                                    Function: Calculation of solar zenith angle
!   WaterVapPres                                Function: Convert water mixing ratio (kg/kg) to water vapor pressure
!   Stability                                   Function: Temperature lapse rate
!   Ea1t99                                      Function: Temperature dependence activity factor for emission type 1
!   Ea1p99                                      Function:
!   DIstomata                                   Function:
!   CalcEccentricity                            Function: Calculate the eccentricity of sun for current day since orbit is not circle
!
!======================================================================
MODULE M2_Canopy

  USE M2_Parameter

CONTAINS

SUBROUTINE GAMMA_CE(JDATE, JTIME, Beta, LAT, LONG, PPFD0, Tc, LAI, Wind, Pres, Humidity, DI, &    ! input
                    PFTF, LADp, midpoint, Layers,            &    ! input 
                    Ealt, Ealp, Ealpt, gammasout, sun_par, Sunleaftk, Shadeleaftk, Sunfrac, SH, LH, &    ! output
                    Rbs, Rds, RLs)    ! output 

  IMPLICIT NONE

  INTEGER :: Layers   ! Number of layers inside the canopy


  INTEGER :: II, JJ, KK,I,J,JDATE,JTIME,K, counter
  INTEGER  Day
  REAL :: Pres, LAI, DI
  REAL :: PPFD0
  REAL :: LONG, Hour

  REAL :: Beta
  REAL :: LAT, Sinbeta
  REAL, DIMENSION(Layers) :: VPslwWT, Sunfrac, QdAbsV, QsAbsV, QdAbsn,           &
                            QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD,              &
                            HumidairPa, Sunleaftk, SunleafSH, SunleafLH, SunleafIR, Shadeleaftk, &
                            ShadeleafSH, ShadeleafLH, ShadeleafIR, Ea1pLayer, Ea1tLayer, Ea1Layer,   &
                            SH, LH

  REAL :: Solar, Maxsolar, Transmis, Difffrac, PPFDfrac, QbAbsn, Trate, StomataDI, QBeamV, &
                     QDiffV, Qbeamn, Qdiffn, QbAbsV, Ea1tCanopy, Ea1pCanopy, Ea1Canopy, TairK0,       &
                     HumidairPa0, Ws0

  REAL, DIMENSION(NPFT)::  Ealpt, Ealt, Ealp, PFTF
                  
  INTEGER :: Cantype

  REAL, DIMENSION(Layers)     :: LADp, midpoint, sun_par
  REAL, DIMENSION(Layers+1)   :: Tc, Wind, Humidity ! tempreature, wind, mixing ratio [g kg-1]
  REAL, DIMENSION(3,NrTyp,Layers) :: gammasout ! gamma factors bundled

  REAL, INTENT(OUT) :: Rbs, Rds, RLs    ! Rds: diffuse radiation, Rbs: direct radiation, RLs: longwave radiation. All are on the surface
  REAL :: sunfrac_s ! Sunfrac at the surface, equals to EXP(-Kd*LAI)

  DO counter=1, Layers
     HumidairPa(counter) = WaterVapPres(Humidity(counter), Pres)
  END DO

  Day =  MOD(JDATE,1000)
  ! Convert from XXXXXX format to XX.XX (solar hour)
  ! HOUR = 0 -> 23.xx
  Hour = JTIME/10000 + LONG/15   ! Solar hour
      IF ( Hour .LT. 0.0 ) THEN
        Hour = Hour + 24.0
        Day = Day - 1
      ELSEIF ( Hour .GE. 24.0 ) THEN
        write(*,*) 'Invalid hour: HOUR(I,J) is ', Hour
      ENDIF

      ! Sinbeta now with Beta values from Scadis
      Sinbeta = SIN((90.-Beta) / 57.29578) ! equal to cosSun in ASAM
      ! Why set it to 0???
      IF (Sinbeta .EQ. 1) THEN
        Sinbeta = 0.
      END IF

      TairK0    = Tc(Layers+1)                           !Air temperature above canopy
      Ws0       = Wind(Layers+1)                         ! Wind speed above canopy
      StomataDI = DIstomata(DI)
      Solar     = PPFD0/ConvertWm2toUmolm2s*2.
      Maxsolar  = Sinbeta * SolarConstant * CalcEccentricity(Day)
  ! Call GaussianIntegration(LADp, midpoint, Layers,NCLOS,NROWS)

  ! NOTE you might want to see if this( LADp and midpoint) is correct. I think it is.

  ! write(*,*) 'GAMMA_CE: LADp, ', LADp
  ! write(*,*) 'GAMMA_CE: midpoint, ', midpoint

  Call SolarFractions(0, Solar, Maxsolar, Transmis, Difffrac, PPFDfrac)

  !=================================================!
  ! How the solar radiation is divided into visible
  ! and near IR, and also diffuse and direct
  !=================================================!
  QDiffV = PPFDfrac * Solar * Difffrac
  QBeamV = PPFDfrac * Solar * (1 - Difffrac)
  QDiffN = (1 - PPFDfrac) * Solar * Difffrac
  QBeamN = (1 - PPFDfrac) * Solar * (1 - Difffrac)

  Call WeightSLW(midpoint, LADp, LAI, Layers, VPslwWT)

  DO K = 1, 1
    Cantype = K
    Call CanopyRad(midpoint, LADp, Layers, LAI, Sinbeta, QBeamV, QDiffV, Qbeamn, Qdiffn, Cantype, Sunfrac, &
                   QbAbsV, QdAbsV, QsAbsV, QbAbsn, QdAbsn, QsAbsn, SunQv, ShadeQv, SunQn, ShadeQn, SunPPFD, &
                   ShadePPFD, PFTF(K), Rbs, Rds, sunfrac_s)

    HumidairPa0 =  WaterVapPres(Humidity(Layers+1), Pres)
    Trate       =  Stability(Cantype, Solar)

    Call CanopyEB(Trate, Layers, midpoint, Cantype, StomataDI, Tc, HumidairPa, Wind, SunPPFD, ShadePPFD, & 
                  SunQv, ShadeQv, SunQn, ShadeQn, Sunleaftk, SunleafSH, SunleafLH,SunleafIR, Shadeleaftk, ShadeleafSH, &
                  ShadeleafLH, ShadeleafIR, Ws0, TairK0, HumidairPa0, RLs, sunfrac_s)

    Ea1tCanopy        = 0.
    Ea1pCanopy        = 0.
    Ea1Canopy         = 0.
    !SH                = 0.
    !LH                = 0.

    !Layers are what we want out. Ea1player = emission activity of light per layer,
    !                             Ea1tlayer = emission activity of temperature per layer,
    !                             Ea1layer  = companied emission activity per layer

    DO II = 1, Layers
      Ea1tLayer(II) = Ea1t99(Sunleaftk(II), Td)   * Sunfrac(II) +     &
                          Ea1t99(Shadeleaftk(II), Td) * (1 - Sunfrac(II))

      Ea1pLayer(II) = Ea1p99(midpoint(II) * LAI, SunPPFD(II))   * Sunfrac(II) +     &
                          Ea1p99(midpoint(II) * LAI, ShadePPFD(II)) * (1 - Sunfrac(II))

      SH(II)        = (SunleafSH(II) * Sunfrac(II) + ShadeleafSH(II) * (1 - Sunfrac(II))) * &
                          LAI * LADp(II)

      LH(II)        = (SunleafLH(II) * Sunfrac(II) + ShadeleafLH(II) * (1 - Sunfrac(II))) * &
                          LAI * LADp(II)

      Ea1Layer(II)  = Ea1t99(Sunleaftk(II), Td)   * Ea1p99(midpoint(II) * LAI, SunPPFD(II))   * Sunfrac(II) + & 
                          Ea1t99(Shadeleaftk(II), Td) * Ea1p99(midpoint(II) * LAI, ShadePPFD(II)) * (1 - Sunfrac(II))

      IF (K .eq. 1) then
        IF (SunPPFD(1) .GT. 0.) THEN
          sun_par(II) = SunPPFD(II) / SunPPFD(Layers)
        ELSE
          sun_par(II) = 1.
        ENDIF
      ENDIF
      !What is VPslwWT and should it be used?

      Ea1pCanopy   = Ea1pCanopy + Ea1pLayer(II) * VPslwWT(II) * LADp(II)
      Ea1tCanopy   = Ea1tCanopy + Ea1tLayer(II) * VPslwWT(II) * LADp(II)
      Ea1Canopy    = Ea1Canopy  + Ea1Layer(II)  * VPslwWT(II) * LADp(II)

      gammasout(1, K ,II) = Ea1tLayer(II)                    ! * VPslwWT(II)  ! * LADp(II)
      gammasout(2, K ,II) = Ea1pLayer(II)                    ! * VPslwWT(II)  ! * LADp(II)
      gammasout(3, K ,II) = Ea1Layer(II) !* Cce * LAI    ! *LADp(II) ! * LADp(II)
    ENDDO
    Ea1Canopy  = Ea1Canopy !* Cce * LAI ! and Cce is a parameter, about what?
    Ealpt(K)   = Ea1Canopy
    Ealt(K)    = Ea1tCanopy
    Ealp(K)    = Ea1pCanopy 

  ENDDO
  RETURN
END SUBROUTINE GAMMA_CE
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Calculates the solar zenith angle
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION Calcbeta(Day, Lat, Hour)
  IMPLICIT NONE

  REAL :: Pi, Rpi, Rpi180, Hour, Lat, SinDelta, &
          CosDelta, A, B, Sinbeta, Calcbeta
  INTEGER :: Day
  
  Pi       = 3.14159
  Rpi180   = 57.29578
  
  SinDelta = -SIN(0.40907) * COS(6.28 * (Day + 10) / (365))
  CosDelta = (1 - SinDelta**2.)**0.5
  
  A = SIN(Lat / Rpi180) * SinDelta
  B = COS(Lat / Rpi180) * Cosdelta
  Sinbeta = A + B * COS(2 * Pi * (Hour - 12) / 24)
  Calcbeta = ASIN(Sinbeta) * 57.29578
END FUNCTION Calcbeta

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION DIstomata
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

       FUNCTION DIstomata(DI)

         IMPLICIT NONE
         REAL :: DIhigh, DIlow, DI, DIstomata

! > -.5 incipient,  mild or no drought; < -4 extreme drought
         DIhigh = -0.5
         DIlow  = -5

         IF (DI > DIhigh) THEN
         DIstomata = 1.0  ! no drought
         ELSEIF (DI > DIlow) THEN
         DIstomata = 1.0 - (0.9 * ((DI - DIhigh) / (DIlow - DIhigh))) ! interpolate
         ELSE
         DIstomata = 0  ! Maximum drought, maximum stomatal resistance
         ENDIF

       END FUNCTION DIstomata
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

       FUNCTION CalcEccentricity(Day)
         IMPLICIT NONE

         INTEGER :: Day
         REAL :: CalcEccentricity

         CalcEccentricity = 1. + 0.033 * COS(2. * 3.14 * (Day - 10.) / 365.)
       END FUNCTION CalcEccentricity

!--This subroutine (Gaussianintegration) is replaced with values for the layers coming from Sosa--
!--------------------------------------------------------------------------------
!        SUBROUTINE GaussianIntegration(Weightgauss, Distgauss, 
!     &                                   Layers,NCLOS,NROWS)
!
!        IMPLICIT NONE
!        INTEGER :: I,J, II, Layers,NCLOS,NROWS
!
!        REAL, DIMENSION(NCLOS,NROWS,Layers) :: 
!     &                     Weightgauss, Distgauss
!
!        DO I=1, NCLOS
!        DO J=1, NROWS
!
!       IF (Layers .EQ. 1) THEN
!       Weightgauss(I,J,1) = 1
!       Distgauss(I,J,1)   = 0.5
!       ELSEIF (Layers .EQ. 3) THEN
!       Weightgauss(I,J,1) = 0.277778
!       Weightgauss(I,J,2) = 0.444444
!       Weightgauss(I,J,3) = 0.277778
!       Distgauss(I,J,1)   = 0.112702
!       Distgauss(I,J,2)   = 0.5
!       Distgauss(I,J,3)   = 0.887298
!        ELSEIF (Layers .EQ. 5) THEN
!        Weightgauss(I,J,1) = 0.1184635
!        Weightgauss(I,J,2) = 0.2393144
!        Weightgauss(I,J,3) = 0.284444444
!        Weightgauss(I,J,4) = 0.2393144
!        Weightgauss(I,J,5) = 0.1184635
!        Distgauss(I,J,1)   = 0.0469101
!        Distgauss(I,J,2)   = 0.2307534
!        Distgauss(I,J,3)   = 0.5
!        Distgauss(I,J,4)   = 0.7692465
!        Distgauss(I,J,5)   = 0.9530899
!        ELSE
!        DO II = 1, Layers
!        Weightgauss(I,J,II) = 1. / Layers
!        Distgauss(I,J,II)   = (II - 0.5) / Layers
!         ENDDO
!         ENDIF
!         ENDDO
!         ENDDO
!         RETURN
!         END SUBROUTINE GaussianIntegration


!-------------------------------------------------------------------------------
!   SUBROUTINE WeightSLW
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        !This should be checked, Anton


        SUBROUTINE WeightSLW(Distgauss, Weightgauss, LAI, Layers, SLW)
          IMPLICIT NONE
  
          INTEGER :: II, Layers,I,J,NCLOS,NROWS
          REAL ::  SLWsum, RealII
          REAL, DIMENSION(Layers) :: &
                          Distgauss, Weightgauss, SLW
          REAL :: LAI

          SLWsum = 0
          DO II = 1, Layers
             RealII = II
             SLW(II) = 0.63 + 0.37 * EXP(-(LAI * Distgauss(II)) * ((RealII - 1) / Layers))
             SLWsum = SLWsum + SLW(II) * Weightgauss(II)
          ENDDO
          DO II = 1, Layers
             SLW(II) = SLW(II) / SLWsum
          ENDDO
          ! write(*,*) 'WeightSLW: SLW, ', SLW
        END SUBROUTINE WeightSLW

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   SUBROUTINE SolarFractions
!
!   Transmission, fraction of PPFD that is diffuse, fraction of solar rad that is PPFD
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        SUBROUTINE SolarFractions(Timeperiod, Solar, Maxsolar, Transmis, FracDiff, PPFDfrac)
          IMPLICIT NONE

          INTEGER :: Timeperiod, I, J
          REAL :: Solar, Maxsolar, Transmis, FracDiff, PPFDfrac
          REAL ::  TransMin, TransSlope

          IF (Timeperiod .EQ. 1) THEN     ! Daily transmission
            TransMin   = 0.26
            TransSlope= 1.655
          ELSE                       ! Hourly transmission
            TransMin   = 0.26
            TransSlope = 1.655
          ENDIF

          IF (Maxsolar <= 0) THEN
            Transmis = 0.5
          ELSE
            Transmis = Solar / Maxsolar
          ENDIF

          ! Estimate diffuse fraction based on daily transmission (Roderick 1999, Goudriann and Van Laar 1994- P.33)
          IF (Transmis > 0.81) THEN
            FracDiff = 0.05
          ELSEIF (Transmis > TransMin) THEN
            FracDiff = 0.96-TransSlope * (Transmis - TransMin)
          ELSE
            FracDiff = 0.96
          ENDIF

          !The fraction of total solar radiation that is PPFD (43% to 55%) G. and L. 84
          PPFDfrac = 0.43 + FracDiff * 0.12
        END SUBROUTINE SolarFractions
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   Subroutine CanopyRad
!
!   Canopy light environment model
!   Code developed by Alex Guenther, based on Spitters et al. (1986), Goudrian and Laar (1994), Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

SUBROUTINE CanopyRad(midpoint, LADp, Layers, LAI, Sinbeta, &
                     QBeamV, QDiffV, Qbeamn, Qdiffn, Cantype, &
                     Sunfrac, QbAbsV, QdAbsV, QsAbsV, & 
                     QbAbsN, QdAbsN, QsAbsN, SunQv, &  
                     ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, & 
                     PFTP, &
                     Rbs, Rds, sunfrac_s)

  IMPLICIT NONE

  INTEGER :: II, JJ, Layers

  REAL :: ScatV, ScatN, RefldV, RefldN, ReflbV, &
          ReflbN,Kb, Kd, KbpV, KbpN, KdpV, KdpN, LAIdepth, &
          Cluster, QdAbsVL, QsAbsVL, QdAbsNL, &
          QsAbsNL, QqV, Qqn

  REAL,DIMENSION(Layers) :: midpoint, LADp, &
                Sunfrac, QdAbsV, QsAbsV, QdAbsn, QsAbsn, SunQv, &
                ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD

  REAL  :: Solar, Maxsolar, PFTP, &
           Difffrac,PPFDfrac, QbAbsn,QBeamV,Sinbeta,QDiffV, &
           Qbeamn, Qdiffn, QbAbsV

  INTEGER :: Cantype

  REAL :: LAI

  REAL, INTENT(OUT) :: Rbs, Rds, sunfrac_s

  ! Scattering coefficients (scatV,scatN), diffuse and beam reflection coefficients (ref..) for visible or near IR
  ScatV   = Canopychar(5,Cantype)
  ScatN   = Canopychar(6,Cantype)
  RefldV  = Canopychar(7,Cantype)
  RefldN  = Canopychar(8,Cantype)
  Cluster = Canopychar(9,Cantype)
  ! Extinction coefficients for black leaves for beam (kb) or diffuse (kd)
  Kb = Cluster * 0.5 / Sinbeta    ! (0.5 assumes a spherical leaf angle distribution (0.5 = cos (60 deg))
  Kd = 0.8 * Cluster    ! (0.8 assumes a spherical leaf angle distribution)
  !===============================!
  ! Set Rbs and Rds as 0 initially
  !===============================!
  Rbs = 0.
  Rds = 0.

  IF ( PFTP.NE.0.0) THEN
    IF (((QBeamV + QDiffV) > 0.001) .AND. (Sinbeta > 0.00002) .AND. (LAI > 0.001)) THEN    ! Daytime

      CALL CalcExtCoeff(QBeamV, ScatV, Kb, Kd, ReflbV, KbpV, KdpV, QbAbsV)
      CALL CalcExtCoeff(QBeamN, ScatN, Kb, Kd, ReflbN, KbpN, KdpN, QbAbsN)

      DO II = 1, Layers
        !===== LAI depth at this layer, Distgauss is the midpoint of each layer
        ! LAIdepth = LAI * midpoint(II)    ! Original equation, but Tian thinks it is not correct, and this is fixed in MEGAN 2.1
        LAIdepth = LAI * (SUM(LADp(II+1:Layers))+0.0*LADp(II))    ! Added by Tian, considering the midpoint can not represent the LAI depth

        Sunfrac(II) = EXP(-Kb * LAIdepth)    ! fraction of leaves that are sunlit

        CALL CalcRadComponents(QDiffV, QBeamV, KdpV, & 
                               KbpV, Kb, ScatV, RefldV, &
                               reflbV, LAIdepth, QbAbsV, QdAbsVL, QsAbsVL) 

        CALL CalcRadComponents(Qdiffn, Qbeamn, KdpN, &
                               KbpN, Kb, ScatN, RefldN, &
                               reflbN, LAIdepth, QbAbsn, QdAbsNL, QsAbsNL) 

        ShadePPFD(II) = (QdAbsVL + QsAbsVL) * &
                          ConvertPPFD / (1 - ScatV)
        SunPPFD(II) = ShadePPFD(II) + (QbAbsV * & 
                      ConvertPPFD / (1 - ScatV))
        QdAbsV(II)    = QdAbsVL
        QsAbsV(II)    = QsAbsVL
        QdAbsn(II)    = QdAbsNL
        QsAbsn(II)    = QsAbsNL
        ShadeQv(II)   = QdAbsVL + QsAbsVL
        SunQv(II)     = ShadeQv(II) + QbAbsV
        ShadeQn(II)   = QdAbsNL + QsAbsNL
        SunQn(II)     = ShadeQn(II) + QbAbsn
      ENDDO
      !=============================!
      ! Calculate Rbs, Rds and sunfrac_s
      !=============================!
      Rbs = (1.-ReflbV)*QBeamV*EXP(-KbpV*LAI) + (1.-ReflbN)*QBeamN*EXP(-KbpN*LAI)
      Rds = (1.-RefldV)*QDiffV*EXP(-Kd*LAI) + (1.-RefldN)*QDiffN*EXP(-Kd*LAI)
      sunfrac_s = EXP(-Kd*LAI)    ! Use Kd since this fraction does not depend on the sun angle
    ELSE    ! Night time
      QbAbsV = 0
      QbAbsn = 0

      DO II = 1, Layers
        Sunfrac(II)   = 0.2 ! 0.2
        SunQn(II)     = 0
        ShadeQn(II)   = 0
        SunQv(II)     = 0
        ShadeQv(II)   = 0
        SunPPFD(II)   = 0
        ShadePPFD(II) = 0
        QdAbsV(II)    = 0
        QsAbsV(II)    = 0
        QdAbsn(II)    = 0
        QsAbsn(II)    = 0
      END DO
      !=============================!
      ! Calculate Rbs, Rds and sunfrac_s
      !=============================!
      Rbs = 0.
      Rds = 0.
      sunfrac_s = EXP(-Kd*LAI)    ! Use Kd since this fraction does not depend on the sun angle
    END IF    ! if daytime or night time
  ELSE    ! PFTP == 0
    DO II = 1, Layers
      Sunfrac(II)   = 0.
      SunQn(II)     = 0
      ShadeQn(II)   = 0
      SunQv(II)     = 0
      ShadeQv(II)   = 0
      SunPPFD(II)   = 0
      ShadePPFD(II) = 0
      QdAbsV(II)    = 0
      QsAbsV(II)    = 0
      QdAbsn(II)    = 0
      QsAbsn(II)    = 0
    ENDDO
    QbAbsV = 0
    QbAbsn  = 0
  ENDIF
END SUBROUTINE CanopyRad

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   Subroutine CalcExtCoeff
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        SUBROUTINE CalcExtCoeff(Qbeam, scat, Kb, Kd, Reflb, Kbp, Kdp, QbeamAbsorb)

          IMPLICIT NONE

          REAL, INTENT(IN) :: Qbeam, scat, Kb, Kd
          REAL, INTENT(OUT) :: Reflb, Kbp, Kdp, QbeamAbsorb
          REAL :: P

          P     = (1. - scat)**0.5
          Reflb = 1. - Exp(-(1. - P) / (1. + P) * 2. * kb / (1. + kb))

          ! Extinction coefficients
          Kbp   = Kb * P
          Kdp   = Kd * P
          QbeamAbsorb = kb * Qbeam * (1 - scat) !calculate Absorbed beam radiation
        END SUBROUTINE CalcExtCoeff

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   Subroutine CalcRadComponents
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        SUBROUTINE CalcRadComponents(Qdiff, Qbeam, kdp, kbp, kb, & 
             scat, Refld, Reflb, LAIdepth, QbAbs, QdAbs, QsAbs)

          IMPLICIT NONE

          REAL :: Qdiff, Qbeam, kdp, kbp, kb, scat, refld, reflb, &
                         LAIdepth, QbAbs, QdAbs, QsAbs

          QdAbs = Qdiff * Kdp * (1 - Refld) * Exp(-Kdp * LAIdepth)
          QsAbs = Qbeam * ( Kbp * (1 - Reflb) * Exp(-Kbp * LAIdepth) - & 
                            Kb * (1 - Scat) * Exp(-Kb * LAIdepth) )
        END SUBROUTINE CalcRadComponents


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Convert water mixing ratio (g/kg) to water vapor pressure (Pa or Kpa depending on units of input )
!   Mixing ratio (g/kg), temp (C), pressure (KPa)!
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION WaterVapPres(Dens, Pres)
  IMPLICIT NONE
  
  REAL :: Dens, Dens001, Pres, WaterVapPres
  
  Dens001 = Dens * 0.001 !uncommented rusan, [g kg-1] -> [kg kg-1]
  WaterVapPres = (Dens001 / (Dens001 + WaterAirRatio)) * Pres
  ! WaterVapPres = Dens001 / (1+Dens001) / WaterAirRatio * Pres

END FUNCTION WaterVapPres

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!   FUNCTION
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION Stability(Cantype, Solar)
  IMPLICIT NONE

  INTEGER :: Cantype
  REAL :: Solar, Trateboundary, Stability

  Trateboundary = 500

  IF (Solar > Trateboundary) THEN
    ! Daytime temperature lapse rate
    Stability = Canopychar(12, Cantype)
  ELSEIF (Solar > 0) THEN
    Stability = Canopychar(12, Cantype) - &
                ((Trateboundary - Solar) / Trateboundary) * &
                (Canopychar(12, Cantype) - Canopychar(13, Cantype))
  ELSE
    ! Nightime temperature lapse rate
    Stability = Canopychar(13, Cantype)
  ENDIF

END FUNCTION Stability

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   Subroutine CanopyEB
!
!   Canopy energy balance model for estimating leaf temperature
!   Code developed by Alex Guenther, based on Goudrian and Laar (1994), Leuning (1997)
!   Initial code 8-99, modified 7-2000 and 12-2001
!
!   Note: i denotes an array containing a vertical profile through the canopy with 0 (above canopy conditions) plus 1 to number of canopy layers
!
!
! NOTE winnd and temperature are now directly give to this subroutine asa layers, It doesn't calculate them anymore,
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
SUBROUTINE CanopyEB(Trate, Layers, Distgauss, &
                    Cantype, StomataDI, TairK, HumidairPa, Ws, &
                    SunPPFD, ShadePPFD, SunQv, ShadeQv, SunQn, ShadeQn, &
                    Sunleaftk, SunleafSH, SunleafLH, &
                    SunleafIR, Shadeleaftk, ShadeleafSH, &
                    ShadeleafLH, ShadeleafIR, Ws0, &
                    TairK0, HumidairPa0, &
                    RLs, sunfrac_s)


  IMPLICIT NONE

  INTEGER :: II, Layers,TAGI1,TAGI2

  REAL :: Cdepth, Lwidth, Llength, Cheight, Eps, &
          TranspireType,Deltah, Ldepth, Wsh, &
          IRin,IRout

  REAL,DIMENSION(Layers) ::  Distgauss, SunQv,ShadeQv, SunQn, ShadeQn, SunPPFD, ShadePPFD, &
                             TairK, HumidairPa, Ws, Sunleaftk, SunleafSH, SunleafLH, &
                             SunleafIR, Shadeleaftk, ShadeleafSH, ShadeleafLH, ShadeleafIR

  REAL ::  Trate, StomataDI, TairK0, HumidairPa0, Ws0
  INTEGER :: Cantype

  REAL :: RLs, sunfrac_s

  Cdepth        = Canopychar(1, Cantype)    ! Never used
  Lwidth        = Canopychar(2, Cantype)    ! Never used
  Llength       = Canopychar(3, Cantype)
  Cheight       = Canopychar(4, Cantype)    ! Never used
  Eps           = Canopychar(10,Cantype)
  TranspireType = Canopychar(11,Cantype)

  DO II = 1, Layers
    IRin            = UnexposedLeafIRin(TairK(II), Eps)
    ShadeleafIR(II) = 2 * IRin
    SunleafIR(II)   = 0.5 * ExposedLeafIRin(HumidairPa0, TairK0) + 1.5 * IRin

    ! Sun
    ! write(*,*) '===== Sunleaf, layer, SunQv, SunQn, TairK:', II, SunPPFD(II), SunQv(II), SunQn(II), TairK(II)
    CALL LeafEB(SunPPFD(II), SunQv(II) + SunQn(II), &
                SunleafIR(II), Eps, TranspireType, Lwidth, Llength, &
                TairK(II), HumidairPa(II), Ws(II), &
                Sunleaftk(II), &
                SunleafSH(II),SunleafLH(II), IRout, &
                StomataDI)

    SunleafIR(II) = SunleafIR(II) - IRout

    ! Shade
    ! write(*,*) '----- Shadeleaf, layer, ShadeQv, ShadeQn, TairK:', II, ShadePPFD(II), ShadeQv(II), ShadeQn(II), TairK(II)
    CALL LeafEB(ShadePPFD(II), ShadeQv(II) + ShadeQn(II), &
                ShadeleafIR(II), Eps, TranspireType, Lwidth, Llength, &
                TairK(II), HumidairPa(II), Ws(II), &
                Shadeleaftk(II), &
                ShadeleafSH(II),ShadeleafLH(II), &
                IRout, StomataDI)

    ShadeleafIR(II) = ShadeleafIR(II) - IRout
  ENDDO
  ! write(*,*) 'CanopyEB: SunleafIR(1), ShadeleafIR(1), IRout, sunfrac_s.', SunleafIR(1), ShadeleafIR(1), IRout, sunfrac_s
  !===== Calculate RLs =====!
  RLs = sunfrac_s * ExposedLeafIRin(HumidairPa0, TairK0) + (1.-sunfrac_s)*UnexposedLeafIRin(TairK(1), Eps)
END SUBROUTINE CanopyEB

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION UnexposedLeafIRin
!
!   Calculate IR into leaf that is not exposed to the sky
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

FUNCTION UnexposedLeafIRin(Tk, Eps)
  IMPLICIT NONE

  REAL :: Eps, Tk, UnexposedLeafIRin

  UnexposedLeafIRin = Eps * Sb * (Tk**4.) !

END FUNCTION UnexposedLeafIRin


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION ExposedLeafIRin
!
!   Calculate IR into leaf that is exposed to the sky
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

FUNCTION ExposedLeafIRin(HumidPa, Tk)
  IMPLICIT NONE

  REAL :: Tk, HumidPa, EmissAtm, ExposedLeafIRin

  ! Apparent atmospheric emissivity for clear skies: function of water vapor pressure (Pa) and ambient
  ! Temperature (K) based on Brutsaert(1975) referenced in Leuning 1997

  EmissAtm        = 0.642 * (HumidPa / Tk)**(1./7.)
  ExposedLeafIRin = EmissAtm * Sb * (Tk**4.)

END FUNCTION ExposedLeafIRin


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   Subroutine LeafEB
!
!   Leaf energy balance
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
SUBROUTINE LeafEB(PPFD, Q, IRin, Eps, TranspireType, &
                  Lwidth, Llength, TairK, HumidairPa, Ws, Tleaf, &
                  SH, LH, IRout, StomataDI)

  IMPLICIT NONE

  INTEGER :: II
  REAL :: PPFD, Q, IRin, Eps, TranspireType, Lwidth, &
          Llength, TairK, HumidairPa, &
          Ws, Tleaf, SH, LH, IRout, StomataDI, HumidAirKgm3, & 
          GHforced, StomRes, IRoutairT, &
          LatHv, LHairT, &  
          Tdelt, Balance, E1, GH1, SH1, &
          LH1, IRout1, GH
  REAL :: Ws1

  IF (Ws <= 0) Then
    Ws1 = 0.001
  ELSE
    Ws1 = Ws
  END IF

  ! Air vapor density kg m-3
  HumidAirKgm3 = ConvertHumidityPa2kgm3(HumidairPa, TairK)

  ! Heat convection coefficient (W m-2 K-1) for forced convection. Nobel page 366
  GHforced = 0.0259 / (0.004 * ((Llength / Ws1)**0.5))

  ! Stomatal resistence s m-1
  StomRes  = ResSC(PPFD, StomataDI)

  IRoutairT = LeafIROut(TairK, Eps)
  ! Latent heat of vaporization (J Kg-1)
  LatHv = LHV(TairK)

  ! Latent heat flux
  LHairT = LeafLE(TairK, HumidAirKgm3, LatHv, &
                  GHforced, StomRes, TranspireType)

  E1 = (Q + IRin - IRoutairT - LHairT)
  IF (E1 .EQ. 0.) THEN
    E1 = -1.
  END IF

  Tdelt = 1 ! Tleaf-Tair
  Balance = 10
  ! write(*,*) 'Balance: Q, IRin, IRout, SH1, LH, Tdelt, Balance'
  DO II = 1, 10
    IF (ABS(Balance) > 2) THEN
      GH1 = LeafBLC(GHforced, Tdelt, Llength)       ! Boundary layer conductance
      SH1 = LeafH(Tdelt, GH1)                       ! Convective heat flux
      LatHv = LHV(TairK + Tdelt)                      ! Latent heat of vaporization (J Kg-1)
      LH = LeafLE(TairK + Tdelt, HumidAirKgm3, &
                       LatHv, GH1, StomRes, TranspireType)
      LH1 = LH - LHairT

      IRout = LeafIROut(TairK + Tdelt, Eps)
      IRout1  = IRout - IRoutairT

      Tdelt = E1 / ((SH1 + LH1 + IRout1) / Tdelt)
      Balance = Q + IRin - IRout - SH1 - LH
      ! write(*,*) 'Balance: ', Q, IRin, IRout, SH1, LH, Tdelt, Balance
    ELSE
      EXIT
    ENDIF
  ENDDO

  If (Tdelt > 10)  Tdelt = 10
  If (Tdelt < -10) Tdelt = -10

  Tleaf = TairK + Tdelt

  GH    = LeafBLC(GHforced, Tleaf - TairK, Llength)
  SH    = LeafH(Tleaf - TairK, GH)
  LH    = LeafLE(Tleaf, HumidAirKgm3, LatHv, &
                 GH, StomRes, TranspireType)
  IRout = LeafIROut(Tleaf, Eps)
END SUBROUTINE LeafEB

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!   FUNCTION
!
!   Saturation vapor density  (kg/m3)
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION ConvertHumidityPa2kgm3(Pa, Tk)
  IMPLICIT NONE

  REAL :: ConvertHumidityPa2kgm3, Pa, Tk

  ConvertHumidityPa2kgm3 = 0.002165 * Pa / Tk ! P=rho*R*T, Rv=R0/Mv, R0 = 8.3144621 [J K-1 mol-1]
END FUNCTION ConvertHumidityPa2kgm3
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Leaf stomatal cond. resistence s m-1
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

          FUNCTION ResSC(Par, StomataDI)
          IMPLICIT NONE

           REAL :: Par, StomataDI, SCadj, ResSC


         SCadj = StomataDI * ((0.0027 * 1.066 * Par) / &
                  ((1 + 0.0027 * 0.0027 * Par**2.)**0.5))

         IF (SCadj < 0.1) THEN
         ResSC = 2000
         ELSE
         ResSC = 200 / SCadj
         ENDIF

        END FUNCTION ResSC


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   IR thermal radiation energy output by leaf
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        FUNCTION LeafIROut(Tleaf, Eps)
        IMPLICIT NONE
        REAL :: Tleaf, Eps, LeafIROut

        LeafIROut = Eps * Sb * (2 * (Tleaf**4.))

        END FUNCTION LeafIROut


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Latent Heat of vaporization(J Kg-1)from Stull p641
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

        FUNCTION LHV(Tk)
        IMPLICIT NONE

        REAL :: Tk, LHV
        LHV = 2.501e6 - (2370. * (Tk - 273.15))

        END FUNCTION LHV


!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Latent energy term in Energy balance
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION LeafLE(Tleaf, Ambvap, LatHv, GH, StomRes, TranspireType)
  IMPLICIT NONE

  REAL :: Tleaf, Ambvap, LatHv, GH, StomRes, &
          TranspireType, LeafRes, Vapdeficit, LeafLE, LE

  LeafRes    = (1 / (1.075 * (GH / 1231.))) + StomRes
  Vapdeficit = (SvdTk(Tleaf) - Ambvap)
  ! Latent heat of vap (J Kg-1) * vap deficit(Kg m-3) / leaf resistence (s m-1)
  LE = TranspireType * (1 / LeafRes) * LatHv * Vapdeficit

  LeafLE = MAX(LE, 0.0)

END FUNCTION  LeafLE
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Boundary layer conductance
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

FUNCTION LeafBLC(GHforced, Tdelta, Llength)
  IMPLICIT NONE

  REAL :: GHforced, Tdelta, Llength, Ghfree, LeafBLC

  ! This is based on Leuning 1995 p.1198 except using molecular conductivity (.00253 W m-1 K-1 Stull p 640) instead of molecular
  ! diffusivity so that you end up with a heat convection coefficient (W m-2 K-1) instead of a conductance for free convection

  IF (Tdelta >= 0) THEN
    GhFree = 0.5 * 0.00253 * ((1.6e8 * Tdelta / & 
             (Llength**3.))**0.25) / Llength
  ELSE
    GhFree = 0
  ENDIF

  LeafBLC = GHforced + GhFree
END FUNCTION LeafBLC

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Convective energy term in Energy balance (W m-2 heat flux from both sides of leaf)
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION LeafH(Tdelta, GH)
  IMPLICIT NONE

  REAL :: Tdelta, GH, LeafH
  ! 2 sides X conductance X Temperature gradient
  LeafH = 2 * GH * Tdelta
END FUNCTION LeafH

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Saturation vapor density  (kg/m3)
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION SvdTk(Tk)
  IMPLICIT NONE

  REAL :: Tk, Svp, SvdTk

  ! Saturation vapor pressure (millibars)
  Svp = 10**((-2937.4 / Tk) - (4.9283 * LOG10(Tk)) + 23.5518)  
  SvdTk = 0.2165 * Svp / Tk
END FUNCTION  SvdTk

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!   Temperature dependence activity factor for emission type 1 (e.g. isoprene, MBO)
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!

         FUNCTION Ea1t99(Tm, Td)

           IMPLICIT NONE

           REAL :: Tm, Td, Ea1t99, Ctm1, Ctm2, Topt, X, Eopt

           IF (Tm < 260) THEN
             Ea1t99 = 0
           ELSE
             ! Energy of activation and deactivation
             Ctm1 = 95
             Ctm2 = 230

             ! Temperature at which maximum emission occurs
             Topt = 312.5 + 0.5 * (Td - 301)
             X    = ((1 / Topt) - (1 / Tm)) / 0.00831
    
             ! Maximum emission (relative to emission at 30 C)
             Eopt   = 1.9 * EXP(0.125 * (Td - 301))
             Ea1t99 = Eopt * Ctm2 * Exp(Ctm1 * X) / &
                      (Ctm2 - Ctm1 * (1 - EXP(Ctm2 * X)))
           END IF
         END FUNCTION  Ea1t99

!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!
!   FUNCTION
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
FUNCTION Ea1p99(LAI, PPFD)
  IMPLICIT NONE

  REAL :: LAI, PPFD, Alpha, C1, Ea1p99

  Alpha  = 0.001 + 0.00085 * LAI
  C1     = 1.42 * EXP(-0.3 * LAI)
  Ea1p99 = (Alpha * C1 * PPFD) / &
           ((1 + Alpha**2. * PPFD**2.)**0.5)
END FUNCTION  Ea1p99

END MODULE M2_Canopy
