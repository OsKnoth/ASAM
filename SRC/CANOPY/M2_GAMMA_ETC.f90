! =======================================================================
!     MODULE GAMMA_ETC
! 
!     This module contain the function to calculate
!     GAMMA_P, GAMMA_T, GAMMA_L, GAMMA_A for BVOCs.
! 
!     CONTAINS: 1)GAMMA_LAI
!               2)GAMMA_P
!               3)GAMMA_TISOP
!               4)GAMMA_TNISP
!               5)GAMMA_A
!               6)GAMMA_S
! 
!     Note:
! 
!     Requirement:
! 
!     CALL: SOLARANGLE
! 
!     Created by Tan 11/21/06 for MEGAN v2.0
! 
!     History:
! 
!     Made some changes in gamma_TNSIP to make it take temprature in layers
! =======================================================================

MODULE M2_GAMMA_ETC
  USE M2_Parameter

  IMPLICIT NONE

! ...  Program I/O parameters

! ...  External parameters

CONTAINS
! ***********************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1 (11/27/06)
!
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }
!           where LDF = light dependent function (non-dimension)
!
!     For ISOPRENE
!                 Assumption: LDF = 1 for isoprene            (11/27/06)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!
!     For NON-ISOPRENE
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! =======================================================================
! ...  Begin module
! =======================================================================


!-----------------------------------------------------------------------
! .....1) Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
! 
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------
SUBROUTINE GAMMA_LAI(LAI, GAM_L)

  IMPLICIT NONE

  REAL, INTENT(IN)  :: LAI
  REAL, INTENT(OUT) :: GAM_L

  GAM_L = 0.49*LAI / SQRT(1+0.2*(LAI**2))

END SUBROUTINE GAMMA_LAI

!-----------------------------------------------------------------------
! .....3) Calculate GAM_T (GAMMA_T) for isoprene (Guenther et al, 2006, ACP)
!-----------------------------------------------------------------------
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
! 
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB
! 
!     SUBROUTINE GAMMA_TISOP returns the GAMMA_T value for isoprene
!-----------------------------------------------------------------------
SUBROUTINE GAMMA_TISOP( TEMP, D_TEMP, LAYERS, GAM_T )
  IMPLICIT NONE

  INTEGER LAYERS
  REAL    TEMP(LAYERS)                   ! hourly surface temperature
  REAL    D_TEMP                         ! daily surface temperature
  REAL    GAM_T(LAYERS)                  ! GAMMA_T

!  Local parameters
  REAL    Eopt, Topt, X(LAYERS)
  REAL    AAA(LAYERS), BBB(LAYERS)
  REAL, PARAMETER :: CT1 = 80.0
  REAL, PARAMETER :: CT2 = 200.0

  Eopt = 1.75 * exp(0.08*(D_TEMP-297.0))
  Topt = 313.0 + 0.6*(D_TEMP-297.0)
  X = ( 1./Topt-1./TEMP ) / 0.00831

  AAA = Eopt*CT2*exp(CT1*X)
  BBB = ( CT2-CT1*( 1-exp(CT2*X) ) )

  GAM_T = AAA/BBB
END SUBROUTINE GAMMA_TISOP

!-----------------------------------------------------------------------
! .....4) Calculate GAM_T (GAMMA_T) for non-isoprene
!-----------------------------------------------------------------------
! 
!             GAMMA_T =  exp[BETA*(T-Ts)]
!           where BETA   = temperature dependent parameter
!                 Ts     = standard temperature (normally 303K, 30C )
! 
!     SUBROUITINE GAMMA_TNISP returns the GAMMA_T value for non-isoprene
!-----------------------------------------------------------------------
SUBROUTINE GAMMA_TNISP(SPCNAM, TEMP, LAYERS, GAM_T)

  IMPLICIT NONE

  INTEGER       INDEX1
  EXTERNAL      INDEX1
  CHARACTER*16  SPCNAM
  INTEGER       VAR, LAYERS
  INTEGER       SPCNUM                             ! Species number
  REAL          TEMP(LAYERS), GAM_T(LAYERS)
  REAL, PARAMETER :: Ts = 303.0

  DO VAR = 1, N_TDF_SPC
    IF (TRIM(SPCNAM) .EQ. TRIM(TDF_SPC(VAR))) THEN
      SPCNUM = VAR
    END IF
  END DO

  ! SPCNUM = INDEX1(SPCNAM,N_TDF_SPC,TDF_SPC)

  GAM_T = EXP( TDF_PRM(SPCNUM)*(TEMP-Ts) )
END SUBROUTINE GAMMA_TNISP
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! .....5) Calculate GAM_A (GAMMA_age)
!-----------------------------------------------------------------------
! 
!            GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!           where Fnew = new foliage fraction
!                 Fgro = growing foliage fraction
!                 Fmat = mature foliage fraction
!                 Fold = old foliage fraction
!                 Anew = relative emission activity for new foliage
!                 Agro = relative emission activity for growing foliage
!                 Amat = relative emission activity for mature foliage
!                 Aold = relative emission activity for old foliage
! 
! 
!             For foliage fraction
!             Case 1) LAIc = LAIp
!               Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
! 
!             Case 2) LAIp > LAIc
!               Fnew = 0.0  , Fgro = 0.0
!               Fmat = 1-Fold
!               Fold = (LAIp-LAIc)/LAIp
! 
!             Case 3) LAIp < LAIc
!               Fnew = 1-(LAIp/LAIc)                       t <= ti
!                    = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
! 
!               Fmat = LAIp/LAIc                           t <= tm
!                    = (LAIp/LAIc) +
!                        ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
! 
!               Fgro = 1 - Fnew - Fmat
!               Fold = 0.0
! 
!           where
!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!                = 2.9                                  Tt >  303
!             tm = 2.3*ti
! 
!             t  = length of the time step (days)
!             ti = number of days between budbreak and the induction of
!                  emission
!             tm = number of days between budbreak and the initiation of
!                  peak emissions rates
!             Tt = average temperature (K) near top of the canopy during
!                  current time period (daily ave temp for this case)
! 
! 
!             For relative emission activity
!             Case 1) Constant
!               Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
! 
!             Case 2) Monoterpenes
!               Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
! 
!             Case 3) Sesquiterpenes
!               Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
! 
!             Case 4) Methanol
!               Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
! 
!             Case 5) Isoprene
!               Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
! 
!     SUBROUTINE GAMMA_A returns GAMMA_A
!-----------------------------------------------------------------------
SUBROUTINE GAMMA_A(JDATE, JTIME, SPC_NAME, LAIARp, LAIARc, TSTLEN, D_TEMP, GAM_A )

  IMPLICIT NONE

  CHARACTER*16  :: FUNCNAME = 'GAMMA_A'

! ...  External parameters
  INTEGER     JDATE, JTIME
  CHARACTER*8 SPC_NAME
  REAL        D_TEMP
  REAL        LAIARp, LAIARc
  INTEGER     TSTLEN
  REAL        GAM_A

! ...  Local parameters
  REAL  Fnew, Fgro, Fmat, Fold

  INTEGER        AINDX          ! relative emission acitivity index
  CHARACTER*256  MESG           ! message buffer

  REAL     LAIp                 ! LAI at previous time step
  REAL     LAIc                 ! LAI at current time step
  INTEGER  t                    ! time step
  REAL     ti                   ! number of days between budbreak and the induction of emission
  REAL     tm                   ! number of days between budbreak and the initiation of peak emissions rates
  REAL     Tt                   ! average temperature (K), daily ave temp

! ...  Choose relative emission activity
  SELECT CASE ( TRIM(SPC_NAME) )
  CASE ('ACTO','ACTA','FORM','CH4','NO','CO')
    AINDX = 1
  CASE ('MYRC','SABI','LIMO','3CAR','OCIM','BPIN','APIN','OMTP')
    AINDX = 2
  CASE ('FARN','BCAR','OSQT')
    AINDX = 3
  CASE ('MEOH')
    AINDX = 4
  CASE ( 'ISOP','MBO' )
    AINDX = 5
  CASE DEFAULT
    WRITE(*,*) 'Error: Chemical species, invalid variable: ', TRIM(SPC_NAME)
    !CALL M3EXIT(FUNCNAME,JDATE,JTIME,MESG,2)
    STOP
  END SELECT

  t = TSTLEN
  LAIc = LAIARc
  LAIp = LAIARp
  Tt   = D_TEMP
! ...  Calculate foliage fraction
  IF (LAIp .EQ. LAIc) THEN
    Fnew = 0.0
    Fgro = 0.1
    Fmat = 0.8
    Fold = 0.1
  ELSEIF (LAIp .GT. LAIc) THEN
    Fnew = 0.0
    Fgro = 0.0
    Fold = ( LAIp-LAIc ) / LAIp
    Fmat = 1. - Fold
  ELSEIF (LAIp .LT. LAIc) THEN
!    Calculate ti and tm
     IF (Tt .LE. 303.0) THEN
       ti = 5.0 + 0.7*(300-Tt)
     ELSEIF (Tt .GT. 303.0) THEN
       ti = 2.9
     END IF
     tm = 2.3*ti

!    Calculate Fnew and Fmat, then Fgro and Fold
!    Fnew
     IF (t .LE. ti) THEN
       Fnew = 1.0 - (LAIp/LAIc)
     ELSEIF (t .GT. ti) THEN
       Fnew = (ti/t) * ( 1-(LAIp/LAIc) )
     END IF

!    Fmat
     IF (t .LE. tm) THEN
       Fmat = LAIp/LAIc
     ELSEIF (t .GT. tm) THEN
       Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )
     END IF

     Fgro = 1.0 - Fnew - Fmat
     Fold = 0.0
  END IF

! ...  Calculate GAMMA_A
  GAM_A = Fnew*Anew(AINDX) + Fgro*Agro(AINDX) + &
          Fmat*Amat(AINDX) + Fold*Aold(AINDX)

END SUBROUTINE GAMMA_A

!-----------------------------------------------------------------------
! .....6) Calculate GAM_SMT (GAMMA_SM)
!-----------------------------------------------------------------------
! 
!             GAMMA_SM =     1.0   (non-dimension)
! 
!    SOILW: [m3 m-3], wilting point (the soil moisture level below which
!           plants can not extract water from the soil [Chen and Dudhia,
!           2001, MWR]
!    DSOIL: [m3 m-3], empirical paramter, =0.06 [Pegoraro etal, 2004, FPB]
!    ISLTYP: soil type
!
!    SUBROUTINE GAMMA_S returns the GAMMA_SM values
!---------------------------------------modify the GAM_SMT_----------
SUBROUTINE GAMMA_S( SMOIS, ISLTYP, GAM_S )

  IMPLICIT NONE

  REAL  GAM_S, SMOIS, SOILW
  REAL, PARAMETER :: DSOIL = 0.06
  REAL :: SOIL1
  INTEGER  ISLTYP

  IF (ISLTYP .EQ. 1) THEN
    SOILW = 0.01
  ELSEIF (ISLTYP .EQ. 14) THEN
    SOILW = 1.0
  ELSEIF ((ISLTYP .EQ. 15) .OR. (ISLTYP .EQ. 16)) THEN
    SOILW = 0.05
  ELSE
    SOILW = 0.138
  END IF

  SOIL1 =  SOILW + DSOIL
  IF (SMOIS >= SOIL1) THEN
    GAM_S = 1.0
  ELSEIF (SMOIS < SOIL1 .AND. SMOIS >= SOILW) THEN
    GAM_S = (SMOIS - SOILW) / DSOIL
  ELSE
    GAM_S =0.0
  END IF

END SUBROUTINE GAMMA_S

END MODULE M2_GAMMA_ETC
