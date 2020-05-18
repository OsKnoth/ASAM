!===========================================================================================
!===  Modul, das saemtliche sea-salt-emission-Parameterisierungen enthaelt, 
!===  und wenn noetig um weitere ergaenzt werden kann
!===  Autor:  Stefan Barthel, Maerz 2010
!============================================================================================
!=== Inhalt: 5 Teile
!=== 1. Parameter und Bemerkungen zu den einzelnen Parameterisierungen
!=== 2. Parameterisierungen
!=== 3. Größenabhängige Funktionen für die teilbaren Parameterisierungen
!=== 4. Windabhängige Funktionen für die teilbaren Parameterisierungen
!=== 5. Funktionen für die OM-abhängige Größeneberechnung
!============================================================================================

MODULE  Salt_Mod

!============================================================================================
!=== 1. Parameter und Bemerkungen zu den einzelnen Parameterisierungen
!============================================================================================

  IMPLICIT NONE

! Parametrisierung
! verwendete Quelle
! urspruengliche Quelle, falls nicht identisch
! Alle Angaben zum Gueltigkeitsbereich entstammen verschiedenen Arbeiten, 
! und entsprechen meist dem Bereich der Messung, aus derer die Parameterisierung 
! gewonnen wurde, daher kann sie auch mehr oder weniger ausserhalb dieses Bereichs 
! Anwendung finden

  REAL(8), PARAMETER :: Zero   = 0.E0
  REAL(8), PARAMETER :: One    = 1.E0
  REAL(8), PARAMETER :: Two    = 2.E0
  REAL(8), PARAMETER :: Three  = 3.E0
  REAL(8), PARAMETER :: Four   = 4.E0
  REAL(8), PARAMETER :: Five   = 5.E0
  REAL(8), PARAMETER :: Six    = 6.E0
  REAL(8), PARAMETER :: Seven  = 7.E0
  REAL(8), PARAMETER :: Eight  = 8.E0
  REAL(8), PARAMETER :: Nine   = 9.E0
  REAL(8), PARAMETER :: Ten    = 1.E1
  REAL(8), PARAMETER :: Eleven = 11.E0
  REAL(8), PARAMETER :: Half   = 5.E-1
  REAL(8) :: ustern, DynViscAir, MeanFreePathAir, z_rauh, Rho_particle, Rho_air
  REAL(8) :: Sed_correct_exp
  REAL(8) :: Karm = 4.E-1
  REAL(8) :: H_Half
  REAL(8) :: dustern, dhwind
  REAL(8), PARAMETER :: ustern_min = 1.E-2
  REAL(8), PARAMETER :: hwind_min = 3.E0
  INTEGER, PARAMETER :: hwind_num = 10000
  REAL(8), ALLOCATABLE :: ustern_tab(:)
  REAL(8) :: ustar_tab(hwind_num)
  CHARACTER (99) :: Func_to_use='Long11_Long11_Sofiev11'

  INTERFACE
    FUNCTION Emiss_flux(Dp80,U10,SST) RESULT(SS_num_flux)
      REAL(8), INTENT(in) :: Dp80,U10,SST
      REAL(8) :: SS_num_flux
    END FUNCTION Emiss_flux

    FUNCTION Emiss_windpart(U10) RESULT(Windpart)
      REAL(8), INTENT(in) :: U10
      REAL(8) :: Windpart
    END FUNCTION Emiss_windpart

    FUNCTION Emiss_sizepart(Dp80,SST) RESULT(SizePart)
      REAL(8), INTENT(in) :: Dp80,SST
      REAL(8) :: SizePart
    END FUNCTION Emiss_sizepart

    FUNCTION SST_correction(Dp80,SST) RESULT(SSTPart)
      REAL(8), INTENT(in) :: Dp80,SST
      REAL(8) :: SSTPart
    END FUNCTION SST_correction

    SUBROUTINE Volume_calculation(Dpd,Omega,ChlAlpha,U10,Dp80,VOM,VSS)
      REAL(8), INTENT (IN) ::  Dpd,Omega,ChlAlpha,U10
      REAL(8), INTENT (OUT) :: Dp80,VOM,VSS
    END SUBROUTINE Volume_calculation
  END INTERFACE
  PROCEDURE(Emiss_flux), POINTER :: EmissFunction
  PROCEDURE(SST_correction), POINTER :: SSTFunction
  PROCEDURE(Emiss_sizepart), POINTER :: SizeFunction
  PROCEDURE(Emiss_windpart), POINTER :: WindFunction
  PROCEDURE(Volume_calculation), POINTER :: VolFunction

!============================================================================================
CONTAINS

!============================================================================================
!=== 2. Parameterisierungen
!============================================================================================
!=== Parameterisierungen nach Zeit und dann nach Alphabet geordnet (geteilt):
!=== Monahan86 (ja), Andreas92 (nein), Smith93 (nein), Fairall94 (ja), Andreas98
!=== (nein), Smith98 (nein), PattBel99 (ja), DeLeeuw00 (ja), Guelle01 (nein),
!=== Vigantti01 (nein), DeLeeuw03 (nein), Gong03 (ja), Martensson03 (ja),
!=== Clarke06 (ja), Long11 (ja), Sofiev11 (ja)
!============================================================================================

FUNCTION Monahan86(Dp80,U10,Tw) RESULT(SS_num_flux)   ! Funktion von Monahan 86

! Monahan, E.C., Spiel, D.E., Davidson, K.L., A model of marine aerosol generation 
!    via whitecapsand wave distribution, in Oceanic Whitecaps and Their Role in Air-Sea Exchange, 
!    E.C. Monhan and G. Mac Niocaill. Eds., D.Reidel, 167-174, 1986
! Gueltigkeitsbereich: Funktion gilt fuer den Bereich des bubble burstings, also fuer 
! alle Seesalzpartikel unterhalb eines Trockenradiuses von etwa 10-15 Mikrometer 
! (Laut Monahan et al. 86 bis 5 Mikrometer Raddry), sowie oberhalb 0,5 Mikrometer (Raddry) (Schulz et al. 04)
! Laut Arbeit von Monahan et al. 86 existiert noch ein zweiter Teil fuer die Spumeregion, 
! der allerdings, wie bspw. in Andreas 92 und 98 gezeigt, die Produktion ueberschaetzt, 
! daher wurde hier darauf verzichtet
! laut Schulz et al. 04 gilt Funktion im Zusammenhang mit dem Bildungsradius, allerdings 
! verwenden verschiedene Modelle (Bsp: HAM, Geos), sowie verschiedene Autoren 
! (Bsp: Gong 03, Zhang et al. 05) die Funktion bei einem Tropfenradius bei 80% Luftfeuchte
! Monahan et al. 86 macht keine Aussage zu diesem Sachverhalt
! daher wurde die Funktion hier bei rel. Feuchte von 80% verwendet

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Monahan86Wind
  SizeFunction => Monahan86_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)
 
END FUNCTION Monahan86


!--------------------------------------------------------------------------------------------
FUNCTION Andreas92(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Andreas 92

! Funktion arbeitet mit Radius bei 80%
! Andreas, E.L. Sea spray and the turbulent air-sea heat fluxes. 
! Journal of Geophysical Research, 97(C7), pp. 11429-11441, 1992

  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: C1,C2,C3,B0,B1,B2,B3,B4,Radius
  REAL(8), PARAMETER :: r1=10.0E0
  REAL(8), PARAMETER :: r2=37.5E0
  REAL(8), PARAMETER :: r3=100.0E0
  REAL(8), PARAMETER :: r4=250.0E0

  Radius = Dp80/Two

! Andreas 92 hat zu diesen Parameter nur einzelne Werte zu best. Windgeschwindigkeiten angegeben. 
! Daher angegebene Gleichungen nur Naeherungen
  B0 = 1.5841E0+1.1759E0*LOG(U10)
  B1 = -1.25E-3*U10**3+0.0384E0*U10**2-0.1766E0*U10-3.7443E0
  B2 = 8.3816E-3*U10**3-0.2612E0*U10**2+1.2477E0*U10+4.0963E0
  B3 = -0.0145E0*U10**3+0.4723E0*U10**2-2.9635E0*U10+2.7509E0
  B4 = 7.5208E-3*U10**3-0.2501E0*U10**2+1.8143E0*U10-3.794E0

  IF (Radius.LE.r1) THEN
    SS_num_flux = Ten**(B0+B1*LOG10(Radius)+B2*(LOG10(Radius))**2+B3*(LOG10(Radius))**3+B4*(LOG10(Radius))**4)
  ELSEIF (Radius.GT.r1.AND.Radius.LE.r2) THEN
    SS_num_flux = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*Radius**(-1)
  ELSEIF (Radius.GT.r2.AND.Radius.LE.r3) THEN
    SS_num_flux = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*r2**(-1)
    C2        = SS_num_flux*r2**2.8
    SS_num_flux = C2*Radius**(-2.8)
  ELSEIF (Radius.GT.r3.AND.Radius.LE.r4) THEN  
    SS_num_flux = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*r2**(-1)
    C2        = SS_num_flux*r2**2.8
    SS_num_flux = C2*r3**(-2.8)
    C3        = SS_num_flux*r3**8
    SS_num_flux = C3*Radius**(-8)
  ELSE
    SS_num_flux = Zero
  ENDIF

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two ! dF/dRp80 -> dF/dDp80

END FUNCTION Andreas92


!--------------------------------------------------------------------------------------------
FUNCTION  Smith93(Dp80,U10,Tw) RESULT(SS_num_flux)  ! Funktion von Smith 93

! Funktion arbeitet mit Radius bei 80%
! Da Smith nur log(A) angibt, und weiterfuehrende Literatur sich widerspricht hier log als log10 
! nach Hoppel 2002 angesehen, und entsprechend bei der Berechnung von A 10^ und nicht e^ gerechnet
! Smith, M.H.,Park, P.M., Consterdine, I.E., Marine aerosol concentrations and estimated 
! fluxes over the sea, Q.J.R. Meteorol. Soc., 119, 809-824, 1993
! Gueltigkeitsbereich: Die Messung, aus der die Parametrisierung gewonnen wurde erfolgte 
! zwischen etwa 1 und 25 Mikrometer Rad80 (Smith und Harrison 98) demnach im Breich von 0,5 
! und etwa 13 Mikrometer Trockenradius
! Die Windgeschwindigkeit ist auf unter 30 m/s begrenzt (Smith et al. 93)

  REAL(8), PARAMETER :: asf1  = 3.1E0
  REAL(8), PARAMETER :: asf2  = 3.3E0
  REAL(8), PARAMETER :: asr01 = 2.1E0  ! Mikrometer
  REAL(8), PARAMETER :: asr02 = 9.2E0  ! Mikrometer
  REAL(8), PARAMETER :: asA11 = 6.76E-2
  REAL(8), PARAMETER :: asA12 = 2.43E0
  REAL(8), PARAMETER :: asA21 = 9.59E-1
  REAL(8), PARAMETER :: asA22 = 1.476E0
  
  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: asA1,asA2,Radius

  Radius=Dp80/Two

  asA1 = Ten**(asA11*U10+asA12)
  asA2 = Ten**(asA21*U10**0.5-asA22)

  SS_num_flux = asA1*exp(-asf1*(log(Radius/asr01))**2)+asA2*exp(-asf2*(log(Radius/asr02))**2) ! dF/dRp80

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two ! dF/dDp80

END FUNCTION Smith93


!--------------------------------------------------------------------------------------------
FUNCTION Fairall94(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Fairall 1994

!  Funktion arbeitet eigentlich mit Bildungsradius, setzt aber bei Andreas 92 an, 
!  daher einige Umformungen noetig, die sich teilweise wieder aufheben
!  vgl. Andreas 2002
! Andreas, E.L., A review of the sea spray generartion function for the open ocean 
! Fairall, C.W., Kepert, J.D., Holland, G.J., The effect of sea spray  on surfacs energy 
! transports over the ocean, The Global Atmosphere and Ocean System, 2(2-3), pp. 121-142, 1994

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Fairall94Wind
  SizeFunction => Fairall94_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Fairall94


!--------------------------------------------------------------------------------------------
FUNCTION Andreas98(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Andreas 98

! Funktion arbeitet mit Radius bei 80%
! Andreas, E.L., A new sea spray generation function for wind speeds up to 32 m*s^-1, J. Phys. Ocenogr. 28, 2175-2184, 1998
! Gueltigkeitsbereich: zwischen 1 und 250 Mikrometer Rad80 (0,5 und 125 Raddry) und bis zu 32 m/s (Andreas 98)
! Funktion ist oberhalb dieses Bereichs indentisch Null

  REAL(8), PARAMETER :: r1=10.0E0
  REAL(8), PARAMETER :: r2=37.5E0
  REAL(8), PARAMETER :: r3=100.0E0
  REAL(8), PARAMETER :: r4=250.0E0

  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: C1,C2,C3,Radius

  Radius=Dp80/Two

  IF (Radius.LE.r1) THEN
    SS_num_flux = 3.5E0*Smith93(Dp80,U10,Tw)*Two
  ELSEIF (Radius.GT.r1.AND.Radius.LE.r2) THEN
    SS_num_flux = 3.5E0*Smith93(r1*Two,U10,Tw)*Two
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*Radius**(-1)
  ELSEIF (Radius.GT.r2.AND.Radius.LE.r3) THEN
    SS_num_flux = 3.5E0*Smith93(r1*Two,U10,Tw)*Two
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*r2**(-1)
    C2        = SS_num_flux*r2**2.8
    SS_num_flux = C2*Radius**(-2.8)
  ELSEIF (Radius.GT.r3.AND.Radius.LE.r4) THEN  
    SS_num_flux = 3.5E0*Smith93(r1*Two,U10,Tw)*Two
    C1        = SS_num_flux*r1**1
    SS_num_flux = C1*r2**(-1)
    C2        = SS_num_flux*r2**2.8
    SS_num_flux = C2*r3**(-2.8)
    C3        = SS_num_flux*r3**8
    SS_num_flux = C3*Radius**(-8)
  ELSE
    SS_num_flux = Zero
  ENDIF

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two ! dF/dRp80 -> dF/dDp80

END FUNCTION Andreas98


!--------------------------------------------------------------------------------------------
FUNCTION  Smith98(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Smith 98

! Funktion arbeitet mit Radius bei 80%

! Funktion nach Smith, M.H., Harrison, N.M., The Seespray Generation Function, J. Aerossol Sci., 29, 189-190, 1998
! Gueltigkeitsbereich: 1 bis 150 Mikrometer Rad 80 (0,5 bis 75 Raddry) und weniger als 20 m/s Wind (Smith und Harrison 98)

  REAL(8), PARAMETER :: ashf1=-1.5E0
  REAL(8), PARAMETER :: ashf2=-1.0E0
  REAL(8), PARAMETER :: ashr01=3.0E0  ! Mikrometer
  REAL(8), PARAMETER :: ashr02=30.0E0  ! Mikrometer
  REAL(8), PARAMETER :: ashA11=0.2E0
  REAL(8), PARAMETER :: ashA21=6.8E-3

  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: ashA1,ashA2,Radius

  Radius=Dp80/Two

  ashA1 = ashA11*U10**3
  ashA2 = ashA21*U10**3.5

  SS_num_flux = ashA1*exp(ashf1*(log(Radius/ashr01))**2)+ashA2*exp(ashf2*(log(Radius/ashr02))**2) ! dF/dRp80

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two ! dF/dDp80

END FUNCTION Smith98


!--------------------------------------------------------------------------------------------
FUNCTION PattBel99(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Pattison und Belcher 1999

! Funktion gilt mit dem Durchmesser bei 80%

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => PattBel99Wind
  SizeFunction => PattBel99_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION PattBel99


!--------------------------------------------------------------------------------------------
FUNCTION DeLeeuw00(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von De Leeuw 2000
! Funktion gilt mit dem Bildungsdurchmesser (Originalpaper)
! De Leeuw, G., Neele, F.P., Hill, M., Smith, M.H., Vignati, Sea spray arosol production by 
! owaves breaking in the surf zone, J. Geophys. Res., 105, 29397-29409, 2000
! Gueltigkeitsbereich: zwischen 0,2 und 3 Mikrometer Trockenradius 
! (entspr. 1,6 bis 20 Mikrometer Bildungsdurchmesser) (de Leeuw et al. 2000) und unterhalb 9 m/s (Schulz et al. 2004)

  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => DeLeeuw00Wind
  SizeFunction => DeLeeuw00_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION DeLeeuw00


!--------------------------------------------------------------------------------------------
FUNCTION Guelle01(Dp80,U10,Tw) RESULT(SS_num_flux) ! Guelle verwendete einen Mix aus Monahan86 und Smith98

! Bis zu einem Trockenradius von 4 Mikrometer gilt Monahan, Schulz et al. 2004
! Darueber Smith
! Schulz, M., De Leeuw, G., Balkanski, Y., Emissions of Atmospheric Trace Compounds, chap. 
! Sea-salt aerosol source funtions and emissions, Ed. Kluwer, 333-359, 2004
! Guelle, W., Schulz, M., Balkanski, Y., Dentener, F., Influence of the source formulation 
! on meddling the atmospheric global distribution of sea salt aerosol, J. Geophys. Res., 106, 27509-27524, 2001
! Gueltigkeitsbereich: entsprechend der beiden verwendten Funktionen zwischen 0,5 und 75 Raddry sowie unter 20 m/s  


  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw

  IF (Dp80/Four .LT. Four) THEN
    SS_num_flux = Monahan86(Dp80,U10,Tw)
  ELSE
    SS_num_flux = Smith98(Dp80,U10,Tw)
  ENDIF

END FUNCTION Guelle01


!--------------------------------------------------------------------------------------------
FUNCTION Vignati01(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Vignati 2001

! Funktion arbeitet mit Radius bei 80%
! Schulz, M., De Leeuw, G., Balkanski, Y., Emissions of Atmospheric Trace Compounds, chap. 
! Sea-salt aerosol source funtions and emissions, Ed. Kluwer, 333-359, 2004
! Vignati, E., de Leeuw, G., Berkowicz, R., Modelling coastal aerosol transport and 
! effects of surf-produced aerosols on processes in the marine boundary layer, 
! J. Geophys. Res., 106, 20225-20238, 2001
! Gueltigkeitsbereich: 6-17 m/s (Schulz et al. 2004)

  REAL(8), PARAMETER :: VN11=0.095E0
  REAL(8), PARAMETER :: VN12=0.283E0
  REAL(8), PARAMETER :: VN21=0.0422E0
  REAL(8), PARAMETER :: VN22=0.288E0
  REAL(8), PARAMETER :: VN31=0.069E0
  REAL(8), PARAMETER :: VN32=3.5E0
  REAL(8), PARAMETER :: VR1=0.2E0
  REAL(8), PARAMETER :: VR2=2.0E0
  REAL(8), PARAMETER :: VR3=12.0E0
  REAL(8), PARAMETER :: VSigma1=1.9E0
  REAL(8), PARAMETER :: VSigma2=2.0E0
  REAL(8), PARAMETER :: VSigma3=3.0E0

  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: VN1,VN2,VN3,pi,Radius

  pi = Four*atan(One)

  Radius = Dp80/Two

  VN1 = Ten**(VN11*U10+VN12)
  VN2 = Ten**(VN21*U10+VN22)
  VN3 = Ten**(VN31*U10-VN32)

  SS_num_flux = (VN1/((Two*pi)**0.5*log10(VSigma1)))*exp(-(log10(Radius)-log10(VR1))**2/(Two*(log10(VSigma1))**2))+&
&             (VN2/((Two*pi)**0.5*log10(VSigma2)))*exp(-(log10(Radius)-log10(VR2))**2/(Two*(log10(VSigma2))**2))+&
&             (VN3/((Two*pi)**0.5*log10(VSigma3)))*exp(-(log10(Radius)-log10(VR3))**2/(Two*(log10(VSigma3))**2)) 
! Vignati01: dF(logR80)/d(logR80)

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/(Dp80*log(Ten)) ! dF(logR80)/d(logR80) -> dF(logR80)/dDp80
! Umrechnung nicht vollständig -> nicht richtig
! dF(logR80)/dDp80 -> dF/dDp80 fehlt noch

END FUNCTION Vignati01


!--------------------------------------------------------------------------------------------
FUNCTION DeLeeuw03(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von De Leeuw 2003
! Funktion gilt mit dem Bildungsdurchmesser (Originalpaper)

  REAL(8), PARAMETER :: DeLeeuw03c1=2.4E-1
  REAL(8), PARAMETER :: DeLeeuw03c2=4.E-1
  REAL(8), PARAMETER :: DeLeeuw03c3=1.E4
  REAL(8), PARAMETER :: DeLeeuw03A11=1.41E0
  REAL(8), PARAMETER :: DeLeeuw03A12=9.8E-1
  REAL(8), PARAMETER :: DeLeeuw03A21=5.1E-1
  REAL(8), PARAMETER :: DeLeeuw03A22=1.82E0
  REAL(8), PARAMETER :: DeLeeuw03B11=-1.E-1
  REAL(8), PARAMETER :: DeLeeuw03B12=1.69E0
  REAL(8), PARAMETER :: DeLeeuw03B2=1.09E0
  REAL(8), INTENT (IN) :: U10, Dp80, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Radius
  REAL(8) :: c,A1,A2,B1
  REAL(8) :: R1=One, R2=One


  Radius = Dp80/Two

  c  = (DeLeeuw03c1*U10+DeLeeuw03c2)*DeLeeuw03c3
  A1 = DeLeeuw03A11*U10+DeLeeuw03A12
  A2 = DeLeeuw03A21*U10-DeLeeuw03A22
  B1 = DeLeeuw03A11*U10+DeLeeuw03B12

  SS_num_flux = c*(A1*exp(-B1*(log(Radius/R1))**2)+A2*exp(-DeLeeuw03B2*(log(Radius/R2))**2))
  SS_num_flux = SS_num_flux*1.E4

  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two

END FUNCTION DeLeeuw03


!--------------------------------------------------------------------------------------------
FUNCTION Gong03(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Gong 03

! da Monahan Ausgangspunkt ist gilt auch die hiesige Funktion bei rel. Feuchte von 80%
! Gong S.L., A parameterization of sea-salt aerosol source function for sub- and super-micron particles, 
! Global Biochemical Cycles, 17, 1097-1103, 2003
 
  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Gong03Wind
  SizeFunction => Gong03_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Gong03


!--------------------------------------------------------------------------------------------
FUNCTION Martensson03(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Martensson 2003

! Funktion benoetigt den Durchmesser des trockenen Partikels allerdings in Meter und nicht wie sonst in Mikrometer
! Martensson, E.M., Nilsson, E.D., de Leeuw, G., Cohen, L.H., Hansson, H.-C., Laboratory simulations 
! and parameterization of primary marine aerosol production, 
! J. Geophys. Res., 108 (D9), 4297, doi: 10.1029/2002JD002263, 2003
! Gueltigkeitsbereich: 0,01 bis 1,4 Mikrometer Trockenradius (Martenson et al. 2003)

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Martensson03Wind
  SizeFunction => Martensson03_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Martensson03


!--------------------------------------------------------------------------------------------
FUNCTION Clarke06(Dp80,U10,Tw) RESULT(SS_num_flux)  ! Funktion von Clarke 06

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Clarke06Wind
  SizeFunction => Clarke06_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Clarke06


!--------------------------------------------------------------------------------------------
FUNCTION Long11(Dp80,U10,Tw)   RESULT(SS_num_flux) ! Funktion von Long 11
! Long M.S., Keene W.D., Kieber D.J., Erickson D.J., Maring H., A sea-state
! based source function for size- and composition-resolved marine aerosol
! production, Atmospheric Chemistry and Physics, 11, 1203-1216, 2011

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Long11Wind
  SizeFunction => Long11_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Long11


!--------------------------------------------------------------------------------------------
FUNCTION Sofiev11(Dp80,U10,Tw) RESULT(SS_num_flux) ! Funktion von Sofiev 2011
! Der Windabhängige Teil ist der gleiche wie bei Monahan86

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

  WindFunction => Sofiev11Wind
  SizeFunction => Sofiev11_sizepart
  SS_num_flux = WindFunction(U10)*SizeFunction(Dp80,Tw)

END FUNCTION Sofiev11


!--------------------------------------------------------------------------------------------
FUNCTION LMS12(Dp80,U10,Tw)  RESULT(SS_num_flux) ! Funktion von Fan 2011 modifiziert mit Long 11

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Wind_part, Size_part

 IF (Dp80.LT.4.E-1) THEN  ! entspr. Dpd=0.8µm
   SS_num_flux = Long11(Dp80,U10,Tw)
 ELSE IF (Dp80.GE.4.E-1.AND.Dp80.LE.3.E0) THEN
   SS_num_flux = MAX(Long11(Dp80,U10,Tw),Monahan86(Dp80,U10,Tw))
 ELSE
   IF (U10.GE.9.E0) THEN
     IF (Dp80.LT.1.2E1) THEN  ! entspr. Dpd=0.8µm
       SS_num_flux = Monahan86(Dp80,U10,Tw)
     ELSE IF (Dp80.GE.1.44E1) THEN
       SS_num_flux = Smith93(Dp80,U10,Tw)
     ELSE
       SS_num_flux = MAX(Monahan86(Dp80,U10,Tw),Smith93(Dp80,U10,Tw))
     END IF
   ELSE
     SS_num_flux = Monahan86(Dp80,U10,Tw)
   END IF
 END IF

END FUNCTION LMS12


!--------------------------------------------------------------------------------------------
FUNCTION ModGuelle(Dp80,U10,Tw)  RESULT(SS_num_flux)

! Funktion arbeitet wie Guelle, nur mit dem Unterschied, dass hier der Uebergang durch 
! das Maximum der beiden Funktionen flexibel gestaltet wird.
! Als wesentlicher Unterschied zu Guelle ist der Verwendung von Gong03 statt Monahan86. 
! Ursache hierfuer ist, das Monahan fuer sehr kleine Tropfen zu große Werte ergibt. Das wurde von Gong korrigiert.

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux

  SS_num_flux = MAX(Gong03(Dp80,U10,Tw),Smith98(Dp80,U10,Tw))

END FUNCTION ModGuelle


!--------------------------------------------------------------------------------------------
FUNCTION MMS(Dp80,U10,Tw) RESULT(SS_num_flux)

! Parameterisation how it is used in COSMO-ART
! For details look in the dissertation of Kristina Lundgren

  REAL(8), INTENT (IN) :: Dp80, U10, Tw
  REAL(8) :: SS_num_flux
  REAL(8) :: Radius

  Radius = Dp80/Two

  IF(Radius<One) THEN
    SS_num_flux = Martensson03(Dp80,U10,Tw)
  ELSE IF(Radius>Ten) THEN
    SS_num_flux = Smith93(Dp80,U10,Tw)
  ELSE
    SS_num_flux = Monahan86(Dp80,U10,Tw)
  END IF

END FUNCTION MMS


!--------------------------------------------------------------------------------------------
! Parameterisierung Grythe 2013
FUNCTION Grythe14(Dp80,U10,Tw) RESULT(SS_num_flux)

  REAL(8), PARAMETER :: GyF1 = 2.35d2
  REAL(8), PARAMETER :: GyF2 = 2.d-2
  REAL(8), PARAMETER :: GyF3 = 6.8d0
  REAL(8), PARAMETER :: GyEF1 = -5.5d-1
  REAL(8), PARAMETER :: GyEF2 = -1.5d0
  REAL(8), PARAMETER :: GyEF3 = -1.d0
  REAL(8), PARAMETER :: GyEQ1 = 1.d-1
  REAL(8), PARAMETER :: GyEQ2 = 3.d0
  REAL(8) :: Radius
  REAL(8), PARAMETER :: GyEQ3 = 3.d1
  REAL(8) :: SS_num_flux
  REAL(8), INTENT (IN) :: U10, Dp80, Tw

  Radius = Dp80/Two

  SS_num_flux = GyF1*U10**3.5*exp(GyEF1*(LOG(Radius/GyEQ1)**2))+ &
&            GyF2*U10**3.5*exp(GyEF2*(LOG(Radius/GyEQ2)**2))+ &
&            GyF3*U10**3*exp(GyEF3*(LOG(Radius/GyEQ3)**2))

  SS_num_flux = SS_num_flux * (Ten/z_rauh)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SS_num_flux = SS_num_flux/Two ! dF/dDp80 = dF/dDpdry * dDpdry/dDp80 = dF/dDpdry * 1/2

END FUNCTION Grythe14


!============================================================================================
!=== 3. Größenabhängige Funktionen für die teilbaren Parameterisierungen
!============================================================================================
!=== Parameterisierungen nach Zeit und dann nach Alphabet geordnet:
!=== Monahan86, Fairall94, PattBel99, DeLeeuw00, Gong03, Martensson03, Clarke06,
!=== Long11, Sofiev11
!============================================================================================

! Functions for the Calculation of the Emissionparameters bevor starting
! Timeintegration

!--------------------------------------------------------------------------------------------
FUNCTION Monahan86_sizepart(Dp80,Tw) RESULT(SizePart)  ! Funktion von Monahan 86

! laut Schulz et al. 04 gilt Funktion im Zusammenhang mit dem Bildungsradius, allerdings 
! verwenden verschiedene Modelle (Bsp: HAM, Geos), sowie verschiedene Autoren 
! (Bsp: Gong 03, Zhang et al. 05) die Funktion bei einem Tropfenradius bei 80% Luftfeuchte
! Monahan et al. 86 macht keine Aussage zu diesem Sachverhalt
! daher wurde die Funktion hier bei rel. Feuchte von 80% verwendet
 
  REAL(8), PARAMETER :: am3 = 5.7E-2
  REAL(8), PARAMETER :: am5 = 1.19E0
  REAL(8), PARAMETER :: am6 = 3.80E-1
  REAL(8), PARAMETER :: am7 = 6.50E-1

  REAL(8) :: SizePart
  REAL(8), INTENT (IN) :: Dp80, Tw
  REAL(8) :: B, Radius

  Radius = Dp80/Two

  B = (am6-log10(Radius))/am7

  SizePart = Radius**(-3)*(One+am3*Radius**1.05)*Ten**(am5*exp(-(B**2))) ! dF/dRp80
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/Two ! dF/dDp80 = dF/dr80 * dr80/dDp80 = dF/dr80 * 1/2

END FUNCTION Monahan86_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Fairall94_sizepart(Dp80,Tw) RESULT(SizePart) ! Funktion von Fairall 1994

!  Funktion arbeitet eigentlich mit Bildungsradius, setzt aber bei Andreas 92 an, 
!  daher einige Umformungen noetig, die sich teilweise wieder aufheben
!  vgl. Andreas 2002

  REAL(8), PARAMETER :: FB0= 4.405E0
  REAL(8), PARAMETER :: FB1=-2.646E0
  REAL(8), PARAMETER :: FB2=-3.156E0
  REAL(8), PARAMETER :: FB3= 8.902E0
  REAL(8), PARAMETER :: FB4=-4.482E0
  REAL(8), PARAMETER :: FC1= 1.02E4
  REAL(8), PARAMETER :: FC2= 6.95E6
  REAL(8), PARAMETER :: FC3= 1.75E17
  REAL(8), PARAMETER :: r1=10.0E0
  REAL(8), PARAMETER :: r2=37.5E0
  REAL(8), PARAMETER :: r3=100.0E0
  REAL(8), PARAMETER :: r4=250.0E0

  REAL(8) :: Andreas92
  REAL(8) :: SizePart
  REAL(8), INTENT (IN) :: Dp80, Tw
  REAL(8) :: C1,C2,C3,B0,B1,B2,B3,B4,Radius

  Radius=Dp80/Two

  IF (Radius.LE.r1) THEN
    Andreas92 = Ten**(B0+B1*LOG10(Radius)+B2*(LOG10(Radius))**2+B3*(LOG10(Radius))**3+B4*(LOG10(Radius))**4)
  ELSEIF (Radius.GT.r1.AND.Radius.LE.r2) THEN
    Andreas92 = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = Andreas92*r1**1
    Andreas92 = C1*Radius**(-1)
  ELSEIF (Radius.GT.r2.AND.Radius.LE.r3) THEN
    Andreas92 = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = Andreas92*r1**1
    Andreas92 = C1*r2**(-1)
    C2        = Andreas92*r2**2.8
    Andreas92 = C2*Radius**(-2.8)
  ELSEIF (Radius.GT.r3.AND.Radius.LE.r4) THEN  
    Andreas92 = Ten**(B0+B1*LOG10(r1)+B2*(LOG10(r1))**2+B3*(LOG10(r1))**3+B4*(LOG10(r1))**4)
    C1        = Andreas92*r1**1
    Andreas92 = C1*r2**(-1)
    C2        = Andreas92*r2**2.8
    Andreas92 = C2*r3**(-2.8)
    C3        = Andreas92*r3**8
    Andreas92 = C3*Radius**(-8)
  ELSE
    Andreas92 = Zero
  ENDIF

  SizePart=Andreas92/Two
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)

END FUNCTION Fairall94_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION PattBel99_sizepart(Dp80,Tw) RESULT(SizePart) ! Funktion von Pattison und Belcher 1999

! Funktion gilt mit dem Durchmesser bei 80%
! Schulz, M., De Leeuw, G., Balkanski, Y., Emissions of Atmospheric Trace Compounds, chap. 
! Sea-salt aerosol source funtions and emissions, Ed. Kluwer, 333-359, 2004
! Pattison, M.J., Belcher, S.E., Production rates of sea-spray droplets, 
! J. Geophys. Res., 104, 18397-18407, 1999
 
  REAL(8), INTENT (IN) :: Dp80, Tw
  REAL(8) :: Radius
  REAL(8) :: SizePart

  Radius = Dp80/Two

  SizePart = (Radius*Two)**(-2.74)

  SizePart = SizePart/Two
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)

END FUNCTION PattBel99_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION DeLeeuw00_sizepart(Dp80,Tw) RESULT(SizePart)! Funktion von De Leeuw 2000
! Funktion gilt mit dem Bildungsdurchmesser

  REAL(8), INTENT (IN) :: Dp80,Tw
  REAL(8) :: Dpform
  REAL(8) :: SizePart

  Dpform=Dp80*Two

  SizePart = Dpform**(-1.65)  ! dF/dDpform ! [#/cm²/s]
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart*1.E4*Two ! dF/dD80 ! [#/m²/s]

END FUNCTION DeLeeuw00_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Gong03_sizepart(Dp80,Tw) RESULT(SizePart)  ! Funktion von Gong 03

! da Monahan Ausgangspunkt ist gilt auch die hiesige Funktion bei rel. Feuchte von 80%

  REAL(8), PARAMETER :: Gong1=1.373E0
  REAL(8), PARAMETER :: Gong2=3.41E0
  REAL(8), PARAMETER :: Gong3=0.057E0
  REAL(8), PARAMETER :: Gong4=3.45E0
  REAL(8), PARAMETER :: Gong5=1.607E0
  REAL(8), PARAMETER :: AGong1=4.7E0
  REAL(8), PARAMETER :: AGong2=30E0
  REAL(8), PARAMETER :: AGong3=-0.017E0
  REAL(8), PARAMETER :: AGong4=-1.44E0
  REAL(8), PARAMETER :: BGong=0.433E0
 
  REAL(8) :: SizePart
  REAL(8),INTENT (IN) :: Dp80,Tw
  REAL(8) :: A,B,Radius

  Radius=Dp80/Two

  A = AGong1*(One+AGong2*Radius)**(AGong3*Radius**AGong4)*(-One)
  B = (BGong-log10(Radius))/BGong

  SizePart = Radius**A*(One+Gong3*Radius**Gong4)*Ten**(Gong5*exp(-B**2))
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/Two

END FUNCTION Gong03_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Martensson03_sizepart(Dp80,Tw) RESULT(SizePart) ! Funktion von Martensson 2003

! Funktion benoetigt den Durchmesser des trockenen Partikels allerdings in Meter und nicht wie sonst in Mikrometer

  REAL(8), PARAMETER :: d10= 7.609E8
  REAL(8), PARAMETER :: d11= 1.829E16
  REAL(8), PARAMETER :: d12= 6.791E23
  REAL(8), PARAMETER :: d13=-1.616E31
  REAL(8), PARAMETER :: d14= 7.188E37
  REAL(8), PARAMETER :: c10=-2.881E6
  REAL(8), PARAMETER :: c11=-3.003E13
  REAL(8), PARAMETER :: c12=-2.867E21
  REAL(8), PARAMETER :: c13= 5.932E28
  REAL(8), PARAMETER :: c14=-2.576E35
  REAL(8), PARAMETER :: d20= 2.279E9
  REAL(8), PARAMETER :: d21=-3.787E16
  REAL(8), PARAMETER :: d22= 2.528E23
  REAL(8), PARAMETER :: d23=-7.310E29
  REAL(8), PARAMETER :: d24= 7.368E35
  REAL(8), PARAMETER :: c20=-6.743E6
  REAL(8), PARAMETER :: c21= 1.183E14
  REAL(8), PARAMETER :: c22=-8.148E20
  REAL(8), PARAMETER :: c23= 2.404E27
  REAL(8), PARAMETER :: c24=-2.452E33
  REAL(8), PARAMETER :: d30=-5.800E8
  REAL(8), PARAMETER :: d31= 1.105E15
  REAL(8), PARAMETER :: d32=-8.297E20
  REAL(8), PARAMETER :: d33= 2.601E26
  REAL(8), PARAMETER :: d34=-2.859E31
  REAL(8), PARAMETER :: c30= 2.181E6
  REAL(8), PARAMETER :: c31=-4.165E12
  REAL(8), PARAMETER :: c32= 3.132E18
  REAL(8), PARAMETER :: c33=-9.841E23
  REAL(8), PARAMETER :: c34= 1.085E29

  REAL(8) :: SizePart
  REAL(8) :: Durchdry
  REAL(8), INTENT (IN) :: Dp80,Tw
  REAL(8) :: A,B

  Durchdry=Dp80/2.E6

  IF (Durchdry*1.E6.GE.0.020E0 .AND. Durchdry*1.E6.LT.0.145E0) THEN
    A = c14*Durchdry**4+c13*Durchdry**3+c12*Durchdry**2+c11*Durchdry+c10
    B = d14*Durchdry**4+d13*Durchdry**3+d12*Durchdry**2+d11*Durchdry+d10
  ELSEIF (Durchdry*1.E6.GE.0.145E0 .AND. Durchdry*1.E6.LT.0.419E0) THEN
    A = c24*Durchdry**4+c23*Durchdry**3+c22*Durchdry**2+c21*Durchdry+c20
    B = d24*Durchdry**4+d23*Durchdry**3+d22*Durchdry**2+d21*Durchdry+d20
  ELSEIF (Durchdry*1.E6.GE.0.419E0 .AND. Durchdry*1.E6.LE.2.8E0) THEN
    A = c34*Durchdry**4+c33*Durchdry**3+c32*Durchdry**2+c31*Durchdry+c30
    B = d34*Durchdry**4+d33*Durchdry**3+d32*Durchdry**2+d31*Durchdry+d30
  ELSE
    A = Zero
    B = Zero
  ENDIF

  SizePart = A*Tw+B
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/(Dp80*log(Ten))

END FUNCTION Martensson03_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Clarke06_sizepart(Dp80,Tw) RESULT(SizePart) ! Funktion von Clarke 2006

! Funktion benoetigt den Durchmesser des trockenen Partikels
! Clarke, A. D., S. R. Owens, and J. Zhou (2006) An ultrafine sea.salt flux
! from breaking waves: Implications for cloud condensation nuclei in the remote
! marine atmosphere, J. Geophys. Res., 111, D06202, doi:10.1029/2005JD006565 

  REAL(8), PARAMETER :: CBeta10=-5.001E3
  REAL(8), PARAMETER :: CBeta11=0.808E6
  REAL(8), PARAMETER :: CBeta12=-1.980E7
  REAL(8), PARAMETER :: CBeta13=2.188E8
  REAL(8), PARAMETER :: CBeta14=-1.144E9
  REAL(8), PARAMETER :: CBeta15=2.290E9
  REAL(8), PARAMETER :: CBeta20=3.854E3
  REAL(8), PARAMETER :: CBeta21=1.168E4
  REAL(8), PARAMETER :: CBeta22=-6.572E4
  REAL(8), PARAMETER :: CBeta23=1.003E5
  REAL(8), PARAMETER :: CBeta24=-6.407E4
  REAL(8), PARAMETER :: CBeta25=1.493E4
  REAL(8), PARAMETER :: CBeta30=4.498E2
  REAL(8), PARAMETER :: CBeta31=0.839E3
  REAL(8), PARAMETER :: CBeta32=-5.394E2
  REAL(8), PARAMETER :: CBeta33=1.218E2
  REAL(8), PARAMETER :: CBeta34=-1.213E1
  REAL(8), PARAMETER :: CBeta35=4.514E-1

  REAL(8) :: SizePart
  REAL(8) :: Durchdry
  REAL(8), INTENT (IN) :: Dp80,Tw

  Durchdry=Dp80/Two

  IF (Durchdry.GE.0.010E0 .AND. Durchdry.LT.0.132E0) THEN
    SizePart = CBeta15*Durchdry**5 + CBeta14*Durchdry**4 + &
&                       CBeta13*Durchdry**3 + CBeta12*Durchdry**2 + CBeta11*Durchdry + CBeta10
  ELSEIF (Durchdry.GE.0.132E0 .AND. Durchdry.LT.1.2E0) THEN
    SizePart = CBeta25*Durchdry**5 + CBeta24*Durchdry**4 + &
&                       CBeta23*Durchdry**3 + CBeta22*Durchdry**2 + CBeta21*Durchdry + CBeta20
  ELSEIF (Durchdry.GE.1.2E0 .AND. Durchdry.LE.8.E0) THEN
    SizePart = CBeta35*Durchdry**5 + CBeta34*Durchdry**4 + &
&                       CBeta33*Durchdry**3 + CBeta32*Durchdry**2 + CBeta31*Durchdry + CBeta30
  ELSE
    SizePart = Zero
  ENDIF

  SizePart = SizePart*1.E4
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/(Dp80*log(Ten)) ! dF/dDp80

END FUNCTION Clarke06_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Long11_sizepart(Dp80,Tw) RESULT (SizePart)  ! Funktion von Long 10

  REAL(8), PARAMETER :: LP11=1.46E0
  REAL(8), PARAMETER :: LP12=1.33E0
  REAL(8), PARAMETER :: LP13=-1.82E0
  REAL(8), PARAMETER :: LP14=8.83E0
  REAL(8), PARAMETER :: LP21=-1.53E0
  REAL(8), PARAMETER :: LP22=-8.1E-2
  REAL(8), PARAMETER :: LP23=-4.26E-1
  REAL(8), PARAMETER :: LP24=8.84E0

  REAL(8) :: SizePart
  REAL(8), INTENT (IN) :: Dp80,Tw
  REAL(8) :: P_Long

  IF(Dp80<=One) THEN
    P_Long = LP11*(Log10(Dp80))**3+LP12*(Log10(Dp80))**2+LP13*(Log10(Dp80))+LP14
  ELSE
    P_Long = LP21*(Log10(Dp80))**3+LP22*(Log10(Dp80))**2+LP23*(Log10(Dp80))+LP24
  ENDIF
  SizePart = Ten**(P_Long) ! dF/dlogDp80
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/(Dp80*Log(Ten)) ! dF/dDp80

END FUNCTION Long11_sizepart


!--------------------------------------------------------------------------------------------
FUNCTION Sofiev11_sizepart(Dp80,Tw) RESULT(SizePart)! Funktion von Sofiev 2011
! Funktion arbeitet mit Trockenradius (Originalpaper)
! Sofiev, M., J. Sopares, M.Prank, G. de Leeuw, and J. Kokkonen (2011), A
! regional-to-global model of emission and transport of sea salt particles in
! the atmosphere, J. Geophys. Res. 116, D21302, doi:10.1029/2010JD014713

  REAL(8), PARAMETER :: SofA=1.E6
  REAL(8), PARAMETER :: SofB11=9.E-2
  REAL(8), PARAMETER :: SofB12=3.E-3
  REAL(8), PARAMETER :: SofC11=5.E-2
  REAL(8), PARAMETER :: SofC12D11=1.05E0
  REAL(8), PARAMETER :: SofD12=0.27E0
  REAL(8), PARAMETER :: SofD13=1.1E0

  REAL(8), INTENT (IN) :: Dp80, Tw
  REAL(8) :: SizePart
  REAL(8) :: Dp_dry
  REAL(8) :: SofB,SofC,SofD

  Dp_dry = Dp80/Two

  SofB = (exp(-SofB11/(Dp_dry+SofB12)))/(Two+exp(-Five/Dp_dry))
  SofC = (One+SofC11*Dp_dry**1.05)/Dp_dry**3
  SofD = Ten**(SofC12D11*exp(-((SofD12-log10(Dp_dry))/SofD13)**2))

  SizePart = SofA*SofB*SofC*SofD
  SizePart = SizePart * (Ten/z_rauh)**Sed_correct_exp ! F_int => F_eff (10m)
  SizePart = SizePart * (H_Half/Ten)**Sed_correct_exp ! F_eff (10m) => F_eff (halbe Zellhöhe)
  SizePart = SizePart/Two

END FUNCTION Sofiev11_sizepart



!============================================================================================
!=== 4. Windabhängige Funktionen für die teilbaren Parameterisierungen
!============================================================================================

! Functions for the Calculation of the Winddependent Parameter for the
! Emissionfluxes

!--------------------------------------------------------------------------------------------
FUNCTION Monahan86Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10
  REAL(8), PARAMETER :: am1 = 1.373E0 
  REAL(8), PARAMETER :: am2 = 3.41E0 

  Windpart = am1*U10**am2

END FUNCTION Monahan86Wind


!--------------------------------------------------------------------------------------------
FUNCTION Fairall94Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = (U10/Eleven)**3.4E0

END FUNCTION Fairall94Wind


!--------------------------------------------------------------------------------------------
FUNCTION PattBel99Wind(U10) RESULT(Windpart) ! Funktion von Pattison und Belcher 1999

  REAL(8), PARAMETER :: PB1=69.2E0
  REAL(8), PARAMETER :: PB2=3.25E0
  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = PB1*U10**PB2

END FUNCTION PattBel99Wind


!--------------------------------------------------------------------------------------------
FUNCTION DeLeeuw00Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10
  REAL(8), PARAMETER :: Leeuw1=1.1E0
  REAL(8), PARAMETER :: Leeuw2=0.23E0

  Windpart = Leeuw1*exp(Leeuw2*U10)

END FUNCTION DeLeeuw00Wind


!--------------------------------------------------------------------------------------------
FUNCTION Martensson03Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = 3.84E-6*U10**3.41

END FUNCTION Martensson03Wind


!--------------------------------------------------------------------------------------------
FUNCTION Gong03Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = Monahan86Wind(U10)

END FUNCTION Gong03Wind


!--------------------------------------------------------------------------------------------
FUNCTION Clarke06Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = Martensson03Wind(U10)

END FUNCTION Clarke06Wind


!--------------------------------------------------------------------------------------------
FUNCTION Long11Wind(U10)  RESULT(Windpart)

  REAL(8), PARAMETER :: LW1=2.E-8
  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = LW1*U10**3.74

END FUNCTION Long11Wind


!--------------------------------------------------------------------------------------------
FUNCTION Sofiev11Wind(U10) RESULT(Windpart)

  REAL(8) :: Windpart
  REAL(8), INTENT (IN) :: U10

  Windpart = Martensson03Wind(U10)

END FUNCTION Sofiev11Wind



!============================================================================================
!=== 5. Funktionen für die OM-abhängige Größeneberechnung
!============================================================================================

! Routines for the Calculation of the total Dropletsize, the SS- and the 
! OM-Volume

!--------------------------------------------------------------------------------------------
SUBROUTINE VolAfterLong11(Dpd,Omega,ChlAlpha,U10,Dp80,VOM,VSS) ! volume relation with OM_SSA=V_OM/V_SS

  REAL(8), INTENT (IN) ::  Dpd,Omega,ChlAlpha,U10
  REAL(8), INTENT (OUT) :: Dp80,VOM,VSS

  REAL(8), PARAMETER :: LGam1=-2.01E0
  REAL(8), PARAMETER :: LGam2=4.E1
  REAL(8), PARAMETER :: LDel11=0.306E0
  REAL(8), PARAMETER :: LDel21=5.6E-2
  REAL(8), PARAMETER :: LDel22=2.08E1
  REAL(8), PARAMETER :: LDpStr=8.E0
  REAL(8) :: Vp80, VSS80, Vpd, Vp80_alt
  REAL(8) :: Pi
  REAL(8) :: GammaN,DeltaN
  INTEGER :: i

  Pi=Four*atan(One)

  Dp80=Dpd*Omega            ! µm
  IF (ChlAlpha>Zero) THEN

    Vpd=Pi/Six*Dpd**3     ! µm^3
    Vp80=Pi/Six*Dp80**3   ! µm^3

    IF (Dp80>One) THEN
      DeltaN=LDel21*LDel22*ChlAlpha/(One+LDel22*ChlAlpha)

      VSS=Vpd/(One+DeltaN)  ! µm^3
      VOM=Vpd-VSS            ! µm^3
      Vp80=VSS*Omega**3+VOM  ! µm^3
      Dp80=(Six/Pi*Vp80)**(One/Three) ! µm

    ELSE
      GammaN=LGam1*LGam2*ChlAlpha/(One+LGam2*ChlAlpha)
      Vp80_alt=Zero
      DO i=1,30
        DeltaN=LDel11*Dp80**GammaN
        VSS80=Omega**3/(Omega**3+DeltaN)*Vp80
        VSS=VSS80/Omega**3
        VOM=Vpd-VSS
        Vp80=VSS80+VOM
        Dp80=(Six/Pi*Vp80)**(One/Three)

        IF (Vp80.EQ.Vp80_alt) EXIT
        Vp80_alt=Vp80

      ENDDO
      DeltaN=LDel11*Dp80**GammaN
      VOM=Pi/Six*Dp80**3*(DeltaN/(Omega**3+DeltaN)) ! VOM=Vp80-VSS80  ! µm^3
      VSS=Pi/Six*Dp80**3/(Omega**3+DeltaN)                            ! µm^3
    END IF

  ELSE
    VSS=Pi/Six*Dpd**3  ! µm^3
    VOM=Zero              ! µm^3
  END IF

  VSS = VSS*Ten**(-1.8E1) ! µm^3 > m^3
  VOM = VOM*Ten**(-1.8E1) ! µm^3 > m^3

END SUBROUTINE VolAfterLong11


!--------------------------------------------------------------------------------------------
SUBROUTINE VolAfterGantt11(Dpd,Omega,ChlAlpha,U10,Dp80,VOM,VSS) ! mass relation with OM_SSA=M_OM/(M_OM+M_SS)

  REAL(8), INTENT (IN) ::  Dpd,Omega,ChlAlpha,U10
  REAL(8), INTENT (OUT) :: Dp80,VOM,VSS

  REAL(8) :: Vp80, VSS80, Vpd, Vp80_alt
  REAL(8) :: Pi
  REAL(8) :: Teil_1, Teil_2, Teil_3
  REAL(8) :: OM_SSA
  REAL(8) :: Rho_OM = 1.3E3
  REAL(8) :: Rho_SS = 2.165E3
  INTEGER :: i

  Pi=Four*atan(One)

  Dp80=Dpd*Omega            ! µm
  IF (ChlAlpha>Zero) THEN
    Vpd=Pi/Six*Dpd**3     ! µm^3
    Vp80=Pi/Six*Dp80**3   ! µm^3

    Teil_1 = One/(One+exp(-2.63E0*ChlAlpha+0.18E0*U10))
    Teil_3 = 0.03E0/(One+exp(-2.63E0*ChlAlpha+0.18E0*U10))

    Vp80_alt=Zero
    DO i=1,20
      Teil_2 = One+0.03E0*exp(6.81E0*Dp80)
      OM_SSA = Teil_1/Teil_2 + Teil_3
      VOM = Vpd*Rho_SS*OM_SSA/((Rho_SS-Rho_OM)*OM_SSA+Rho_OM)
      VSS = Vpd - VOM
      VSS80 = VSS*Omega**3
      Vp80=VSS80+VOM
      Dp80=(Six/Pi*Vp80)**(One/Three)
      IF (Vp80.EQ.Vp80_alt) EXIT
      Vp80_alt=Vp80
    END DO
  ELSE
    VSS=Pi/Six*Dpd**3  ! µm^3
    VOM=Zero              ! µm^3
  END IF

  VSS = VSS*Ten**(-1.8E1) ! µm^3 > m^3
  VOM = VOM*Ten**(-1.8E1) ! µm^3 > m^3

END SUBROUTINE VolAfterGantt11


!--------------------------------------------------------------------------------------------
SUBROUTINE VolAfterGantt12(Dpd,Omega,ChlAlpha,U10,Dp80,VOM,VSS) ! mass relation with OM_SSA=M_OM/(M_OM+M_SS)

  REAL(8), INTENT (IN) ::  Dpd,Omega,ChlAlpha,U10
  REAL(8), INTENT (OUT) :: Dp80,VOM,VSS

  REAL(8) :: Vp80, VSS80, Vpd, Vp80_alt
  REAL(8) :: Pi
  REAL(8) :: Teil_1, Teil_2, Teil_3
  REAL(8) :: OM_SSA
  REAL(8) :: Rho_OM = 1.3E3
  REAL(8) :: Rho_SS = 2.165E3
  REAL(8) :: A=6.d0
  REAL(8) :: B=3.d0
  INTEGER :: i

  Pi=Four*atan(One)

  Dp80=Dpd*Omega            ! µm
  IF (ChlAlpha>Zero) THEN
    Vpd=Pi/Six*Dpd**3     ! µm^3
    Vp80=Pi/Six*Dp80**3   ! µm^3

    Teil_1 = One/(One+exp(B*(-2.63E0*ChlAlpha+0.18E0*U10)))
    Teil_3 = 0.03E0/(One+exp(B*(-2.63E0*ChlAlpha+0.18E0*U10)))

    Vp80_alt=Zero
    DO i=1,20
      Teil_2 = One+0.03E0*exp(6.81E0*Dp80)
      OM_SSA = A*(Teil_1/Teil_2 + Teil_3)
      VOM = Vpd*Rho_SS*OM_SSA/((Rho_SS-Rho_OM)*OM_SSA+Rho_OM)
      VSS = Vpd - VOM
      VSS80 = VSS*Omega**3
      Vp80=VSS80+VOM
      Dp80=(Six/Pi*Vp80)**(One/Three)
      IF (Vp80.EQ.Vp80_alt) EXIT
      Vp80_alt=Vp80
    END DO
  ELSE
    VSS=Pi/Six*Dpd**3  ! µm^3
    VOM=Zero              ! µm^3
  END IF

  VSS = VSS*Ten**(-1.8E1) ! µm^3 > m^3
  VOM = VOM*Ten**(-1.8E1) ! µm^3 > m^3

END SUBROUTINE VolAfterGantt12


SUBROUTINE Sed_correct_exp_sub(Dp80) ! Routine to calculate the exponent for the height correction

  REAL(8) :: Dp80
  REAL(8) :: R_80
  REAL(8) :: Grav = 9.80665e0
  REAL(8) :: v_term, Knud, Slip

  R_80 = Dp80/Two
  Knud = MeanFreePathAir / (Dp80*1.E-6)
  Slip = 1.e0 + Knud * (1.257e0 + 0.4e0*exp(-1.1e0/Knud))
  v_term = Two*(Rho_particle-Rho_air)*Grav*(R_80*1.E-6)**2*Slip/(Nine*Rho_air*DynViscAir)

  Sed_correct_exp = -v_term/(Karm*ustern) ! entspr. Glg. 2.9-4 Lewis Schwartz 2004 

END SUBROUTINE Sed_correct_exp_sub


FUNCTION SST_correction_Jaeg(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction after Jaegle et al. 2011

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart
  REAL(8) :: Temp_surf

  Temp_surf = T_surf - 273.15 ! °C
  SSTPart = 0.3e0+1.e-1*Temp_surf-7.6e-3*Temp_surf**2+2.1e-4*Temp_surf**3

END FUNCTION SST_correction_Jaeg


FUNCTION SST_correction_Sof(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction after Sofiev et al. 2011

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart
  REAL(8) :: SST_correction_sof_1, SST_correction_sof_2
  REAL(8) :: Vorfak, expfak

  IF (T_surf .LE. 278.15) THEN
    SST_correction_sof_1 = 0.092e0*Dp_dry**(-0.96e0) ! valid for Tw=271.15K ; -2°C
    SST_correction_sof_2 = 0.15e0*Dp_dry**(-0.88e0) ! valid for Tw=278.15K ; 5°C
    SSTPart = (SST_correction_sof_1*(278.15e0-T_Surf)+SST_correction_sof_2*(T_Surf-271.15e0))/7.e0
  ELSE IF (T_surf .GT. 278.15 .AND. T_surf .LE. 288.15) THEN
    SST_correction_sof_1 = 0.15e0*Dp_dry**(-0.88e0) ! valid for Tw=278.15K ; 5°C
    SST_correction_sof_2 = 0.48e0*Dp_dry**(-0.36e0) ! valid for Tw=288.15K ; 15°C
    SSTPart = (SST_correction_sof_1*(288.15e0-T_Surf)+SST_correction_sof_2*(T_Surf-278.15e0))/1.e1
  ELSE IF (T_surf .GT. 288.15) THEN
    SST_correction_sof_1 = 0.48e0*Dp_dry**(-0.36e0) ! valid for Tw=288.15K ; 15°C
    SST_correction_sof_2 = 1.e0 ! valid for Tw=298.15K ; 25°C
    SSTPart = (SST_correction_sof_1*(298.15e0-T_Surf)+SST_correction_sof_2*(T_Surf-288.15e0))/1.e1
  END IF

END FUNCTION SST_correction_Sof


FUNCTION SST_correction_Mart(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction extrapolated from Martensson 2003 SSA flux parameterisation

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart
  REAL(8) :: A_Tsurf, y0_Tsurf, Tsurf

  Tsurf = MIN(MAX(T_surf,2.7e2),3.08e2)
  A_Tsurf = -1.658511738e-1*Tsurf + 4.94485513217e1
  y0_Tsurf = 3.09110104e-2*Tsurf - 8.2161126574e0

  SSTPart = y0_Tsurf + A_Tsurf * exp(-Dp_dry/9.02736222e-2)

END FUNCTION SST_correction_Mart


FUNCTION SST_correction_Zab_1(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction extrapolated Zabori 2012 with exponential Interpolation

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart
  !REAL(8) :: Dg=1.3574314425842e-1 ! dN/dlogDp
  !REAL(8) :: log10sigma=4.7081975948974e-1 ! dN/dlogDp
  !REAL(8) :: Num=1.8385887618471e+1 ! dN/dlogDp
  REAL(8) :: Dg=0.439686316235d0 ! dN/dDp
  REAL(8) :: log10sigma=0.4708239830537d0 ! dN/dDp
  REAL(8) :: Num=4.491750890136d0 ! dN/dDp
  REAL(8) :: SST_correction_0
  REAL(8) :: Pi

  Pi=Four*atan(One)

  IF (T_surf .LE. 286.15) THEN
    !SST_correction_0 = Num/((2*Pi)**0.5*log10sigma)*exp(-(log10(Dp_dry/Dg)/log10sigma)**2/2) ! dN/dlogDp
    SST_correction_0 = Num/((2*Pi)**0.5*Dp_dry*log10sigma)*exp(-(log10(Dp_dry/Dg)/log10sigma)**2/2)     ! dN/dDp
    SSTPart = MAX(exp(LOG(SST_correction_0)*(286.15-T_surf)/13.d0),1.d-2)
  ELSE
    SSTPart = 1.e0
  END IF

END FUNCTION SST_correction_Zab_1


FUNCTION SST_correction_Zab_2(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction extrapolated from Zabori 2012 with linear Interpolation

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart
  !REAL(8) :: Dg=1.3574314425842e-1 ! dN/dlogDp
  !REAL(8) :: log10sigma=4.7081975948974e-1 ! dN/dlogDp
  !REAL(8) :: Num=1.8385887618471e+1 ! dN/dlogDp
  REAL(8) :: Dg=0.439686316235d0 ! dN/dDp
  REAL(8) :: log10sigma=0.4708239830537d0 ! dN/dDp
  REAL(8) :: Num=4.491750890136d0 ! dN/dDp
  REAL(8) :: SST_correction_0
  REAL(8) :: SST_correction_13 = 1.d0
  REAL(8) :: Pi

  Pi=Four*atan(One)

  IF (T_surf .LE. 286.15) THEN
    !SST_correction_0 = Num/((2*Pi)**0.5*log10sigma)*exp(-(log10(Dp_dry/Dg)/log10sigma)**2/2) ! dN/dlogDp
    SST_correction_0 = Num/((2*Pi)**0.5*Dp_dry*log10sigma)*exp(-(log10(Dp_dry/Dg)/log10sigma)**2/2)     ! dN/dDp
    SSTPart = MAX((SST_correction_13*(T_surf-273.15)+SST_correction_0*(286.15-T_surf))/13.d0,1.d-2)
  ELSE
    SSTPart = 1.e0
  END IF

END FUNCTION SST_correction_Zab_2

FUNCTION SST_correction_None(Dp_dry,T_surf) RESULT(SSTPart) ! Calculation of the sst-correction extrapolated from Zabori 2012

  REAL(8), INTENT (IN) :: Dp_dry,T_surf
  REAL(8) :: SSTPart

  SSTPart = 1.e0

END FUNCTION SST_correction_None

SUBROUTINE SizeAfterLong11_rev(Dp80,Omega,ChlAlpha,U10,Dpd) ! volume relation with OM_SSA=V_OM/V_SS

  REAL(8), INTENT (IN) ::  Dp80,Omega,ChlAlpha,U10
  REAL(8), INTENT (OUT) :: Dpd

  REAL(8), PARAMETER :: LGam1=-2.01E0
  REAL(8), PARAMETER :: LGam2=4.E1
  REAL(8), PARAMETER :: LDel11=0.306E0
  REAL(8), PARAMETER :: LDel21=5.6E-2
  REAL(8), PARAMETER :: LDel22=2.08E1
  REAL(8), PARAMETER :: LDpStr=8.E0
  REAL(8) :: Vp80, VSS80, Vpd, Vp80_alt, VOM,VSS
  REAL(8) :: Pi
  REAL(8) :: GammaN,DeltaN
  INTEGER :: i

  Pi=Four*atan(One)

  Vp80=Pi/Six*Dp80**3   ! µm^3
  IF (Dp80>One) THEN
    DeltaN=LDel21*LDel22*ChlAlpha/(One+LDel22*ChlAlpha)
  ELSE
    GammaN=LGam1*LGam2*ChlAlpha/(One+LGam2*ChlAlpha)
    DeltaN=LDel11*Dp80**GammaN
  END IF
  Vpd=Vp80*(1.E0+DeltaN)/(Omega**3+DeltaN)
  Dpd=(Six/Pi*Vpd)**(One/Three) ! µm


END SUBROUTINE SizeAfterLong11_rev


!--------------------------------------------------------------------------------------------
END MODULE Salt_Mod
