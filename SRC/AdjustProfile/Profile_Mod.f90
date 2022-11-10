MODULE Profile_Mod

  USE Kind_Mod
  USE Parameter_Mod
  IMPLICIT NONE

  INTEGER :: neq
  INTEGER, PARAMETER :: Dry=1
  INTEGER, PARAMETER :: Moist=2

  CHARACTER(30) :: ProfileType

  INTEGER, PRIVATE :: RhoPos=0
  INTEGER, PRIVATE :: ThPos=0
  INTEGER, PRIVATE :: ThDPos=0
  INTEGER, PRIVATE :: RhoVPos=0
  INTEGER, PRIVATE :: PrePos=0
  INTEGER, PRIVATE :: TempPos=0
  INTEGER, PRIVATE :: TPos=0
  INTEGER, PRIVATE :: QVPos=0
  INTEGER, PRIVATE :: rvPos=0
  INTEGER, PRIVATE :: rtPos=0
  INTEGER, PRIVATE :: ThEPos=0
  INTEGER, PRIVATE :: RhPos=0
  INTEGER, PRIVATE :: uWindPos=0
  INTEGER, PRIVATE :: vWindPos=0
  INTEGER, PRIVATE :: wWindPos=0

  INTEGER :: nzProf
  REAL(RealKind), ALLOCATABLE :: zPProf(:)
  REAL(RealKind), ALLOCATABLE :: zMProf(:)
  REAL(RealKind), ALLOCATABLE :: dzProf(:)

  INTEGER, PRIVATE :: InputUnit=10
  INTEGER :: OutputUnit=11
  CHARACTER(200) :: NameSounding
  CHARACTER(80) :: FileNameOut
  LOGICAL :: Sounding=.FALSE.
  LOGICAL :: WRFData=.FALSE.
  REAL(RealKind) :: Tropopause  = 3000.0d0
  REAL(RealKind) :: PotTemp0=-9.d9, Temp0=-9.d9, Pres0=-9.d9, RelHum0=-9.d9, QV0=-9.d9
  REAL(RealKind) :: dTemp_BL=0.d0, dRH_BL=0.d0, dQV_BL=0.d0
  REAL(RealKind) :: InvLayerHeight=-9.d9, dzProf_ManProf=1.d1, dzProf_Inv=-9.d9
  REAL(RealKind) :: dTemp_Inv=0.d0, dRH_Inv=0.d0, dQV_Inv=0.d0
  REAL(RealKind) :: dTemp_FT=0.0d0, dRH_FT=-9.d9, dQV_FT=-9.d9
  REAL(RealKind) :: PotTempE0=-9.d9, rt0=-9.d9

  NAMELIST /AdjustCtrl/ NameSounding   &
                       ,Sounding       &
                       ,WRFData        &
                       ,FileNameOut    &
                       ,Tropopause     &
                       ,ProfileType    &
                       ,PotTempE0      & ! EquivalentPotential temperature at surface in °C or K
                       ,PotTemp0       & ! Potential temperature at surface in °C or K
                       ,Temp0          & ! Absolute air temperature at surface in °C or K
                       ,Pres0          & ! Air pressure at surface in hPa or Pa
                       ,RelHum0        & ! Relative humidity at surface in %
                       ,QV0            & ! Specific humidity at surface in kg/kg
                       ,rt0            & ! Relative total water content at the surface
                       ,dTemp_BL       & ! Gradient of temperature (potential or absolute dependend on give input) in boundary layer in K/m
                       ,dRH_BL         & ! Gradient of relative humidity in boundary layer in %/m
                       ,dQV_BL         & ! Gradient of specific humidity in boundary layer in kg/kg/m
                       ,InvLayerHeight & ! Altitude of Inversionlayer at top of boundary layer
                       ,dTemp_Inv      & ! Temperature jump (potential or absolute dependend on give input) between the two layers below and above inversion in K/m
                       ,dRH_Inv        & ! Gradient of relative humidity at inversion in %/m
                       ,dQV_Inv        & ! Gradient of specific humidity at inversion in kg/kg/m
                       ,dTemp_FT       & ! Gradient of temperature (potential or absolute dependend on give input) in free troposphere in K/m
                       ,dRH_FT         & ! Gradient of relative humidity in free troposphere
                       ,dQV_FT         & ! Gradient of specific humidity in free troposphere
                       ,dzProf_ManProf     & ! Gridspacing of manual set profile
                       ,dzProf_Inv           ! Thickness of inversion layer

  INTEGER :: NumSounding
  REAL(RealKind), ALLOCATABLE :: Height(:)
  REAL(RealKind), ALLOCATABLE :: Pre(:)
  REAL(RealKind), ALLOCATABLE :: ProfTemp(:)
  REAL(RealKind), ALLOCATABLE :: RH(:)
  REAL(RealKind), ALLOCATABLE :: RhoV(:)
  REAL(RealKind), ALLOCATABLE :: QV(:)
  REAL(RealKind), ALLOCATABLE :: rt(:)
  REAL(RealKind), ALLOCATABLE :: rva(:)
  REAL(RealKind), ALLOCATABLE :: ThE(:)
  INTEGER :: NumValues

CONTAINS


SUBROUTINE ReadInput(FileName)

  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line
  REAL(RealKind) :: z0,z1,zw0,zw1
  INTEGER :: iz,k,nzProfz,nLines
  CHARACTER(5) :: indata_type,DummyChar
  INTEGER :: stat


  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&AdjustCtrl')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=AdjustCtrl)
      EXIT
    END IF
  END DO
1 CONTINUE
  REWIND(InputUnit)
  WRITE(*,*) 'ProfileType',ProfileType
  DO
    READ(InputUnit,*,END=2) Line
    IF (INDEX(Line,'#Gitter')>0) THEN
    !-------------------------------------
       READ(InputUnit,*) indata_type
       READ(InputUnit,*) z0,z1,nzProf
       READ(InputUnit,*) z0,z1,nzProf
       READ(InputUnit,*) z0,z1,nzProf
       EXIT
    END IF   
  END DO     
2 CONTINUE  
  REWIND(InputUnit)
  ALLOCATE(zPProf(0:nzProf))
  ALLOCATE(dzProf(nzProf))
  ALLOCATE(zMProf(nzProf))
  dzProf=(z1-z0)/nzProf
  zPProf(0)=z0
  DO iz=1,nzProf
    zPProf(iz)=zPProf(iz-1)+dzProf(iz)
  END DO  
  DO 
    READ(InputUnit,*,END=3) Line
    !........................................................................
    IF (INDEX(Line,'#zGrid')>0) THEN   ! set different distance z-direction
      READ(InputUnit,*) nlines
      nzProf=0
      zPProf(0)=z0
      DO k=1,nlines
        READ(InputUnit,*) zw0,zw1,nzProfz
        DO iz=nzProf+1,nzProf+nzProfz
          dzProf(iz)=(zw1-zw0)/nzProfz
          zPProf(iz)=zPProf(iz-1)+dzProf(iz)
        END DO
        nzProf=nzProf+nzProfz
      END DO
      EXIT
    END IF   
  END DO     
3 CONTINUE  
  DO iz=1,nzProf
    zMProf(iz)=0.5d0*(zPProf(iz-1)+zPProf(iz))
  END DO  
  CLOSE(UNIT=InputUnit)
  WRITE(*,*) 'Profildatei einlesen'
! Profildatei einlesen
  IF (Sounding) THEN
    OPEN(InputUnit,FILE=TRIM(NameSounding)//'.txt',ACTION='read',IOSTAT=stat)
    IF (stat==0) THEN
      DO k=1,19
        READ(InputUnit,'(A)') DummyChar 
      END DO
      NumSounding=10000
      ALLOCATE(Height(NumSounding))
      ALLOCATE(Pre(NumSounding))
      ALLOCATE(ProfTemp(NumSounding))
      ALLOCATE(RH(NumSounding))
      ALLOCATE(RhoV(NumSounding))
      ALLOCATE(QV(NumSounding))
      DO iz=1,NumSounding
        READ(InputUnit,*) Height(iz),Pre(iz),ProfTemp(iz),RH(iz),DummyChar,DummyChar
        Pre(iz)=Pre(iz)*100.0d0 ! Umrechnung hPa Pa
        ProfTemp(iz)=ProfTemp(iz)+273.15d0
        RH(iz)=RH(iz)/100.0d0
        IF (Height(iz)>Tropopause+10.0d0) THEN
          NumValues=iz
          WRITE (*,*) NumValues,'data rows read from file ',Filename
          EXIT
        END IF
      END DO
    ELSE 
      WRITE (*,*) 'Fehler beim Lesen der Datei ',filename
      STOP
    END IF
    CLOSE(InputUnit)
  ELSE IF (WRFData) THEN
    ALLOCATE(Height(NumSounding))
    ALLOCATE(Pre(NumSounding))
    ALLOCATE(ProfTemp(NumSounding))
    ALLOCATE(RH(NumSounding))
    ALLOCATE(RhoV(NumSounding))
    ALLOCATE(QV(NumSounding))
  ELSE
    WRITE(*,*) 'Create Profile'
    NumSounding=NINT((z1-z0)/dzProf_ManProf+1.d0)
    ALLOCATE(Height(NumSounding))
    DO iz=1,NumSounding
      Height(iz)=z0+dzProf_ManProf*(iz-1)
    END DO
    ALLOCATE(Pre(NumSounding))
    ALLOCATE(ProfTemp(NumSounding))
    ALLOCATE(RH(NumSounding))
    ALLOCATE(rt(NumSounding))
    ALLOCATE(rva(NumSounding))
    ALLOCATE(ThE(NumSounding))
    CALL CreateProfile() ! Values defined on borders of the gridcells
  END IF
END SUBROUTINE ReadInput


FUNCTION  LLv(T)

  REAL(RealKind) :: LLv,T
  LLv=L0-(Cpl-Cpv)*(T-Tes0)
END FUNCTION  LLv

FUNCTION SatVapor(T)
  REAL(RealKind) :: SatVapor,T
  REAL(RealKind) :: T_C
  ! Guide to Meteorological Instruments and Methods of Observation (CIMO Guide) (WMO, 2008)
  T_C=T-273.15d0
  SatVapor=6.112d2*EXP(17.62d0*T_C/(243.12d0+T_C))
END FUNCTION  SatVapor


SUBROUTINE Res(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  SELECT CASE(ProfileType)
    CASE('Dry')
      CALL ResDry(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistRelHum')
      CALL ResMoistRelHum(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistRelHum2')
      CALL ResMoistRelHum2(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistRhoV')
      CALL ResMoistRhoV(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistQV')
      CALL ResMoistQV(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistQVPot')
      CALL ResMoistQVPot(z,y,yprime,delta,ires,rpar,ipar)
    CASE('MoistQTPotEquiv')
      CALL ResMoistQTPotEquiv(z,y,yprime,delta,ires,rpar,ipar)
  END SELECT  

END SUBROUTINE Res

SUBROUTINE JAC(z,y,yprime,PD,CJ,RPAR,IPAR)
  INTEGER :: ipar(:)
  REAL(RealKind) :: z,y(1:neq),yprime(1:neq),PD(1:neq,1:neq),CJ,rpar(:)
  SELECT CASE(ProfileType)
    CASE('MoistQVPot')
      CALL JacResMoistQVPot(z,y,yprime,pd,cj,rpar,ipar)
  END SELECT  
END SUBROUTINE Jac

SUBROUTINE Res1(t,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: t,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Theta0
  REAL(RealKind) :: p,RhoD,RhoV,RhoVL,qv,qt,Te,Theta,ThetaM,ThetaE,pPrime
  REAL(RealKind) :: e,qt0,RelHum,pvs,fac

  p        =y(1)
  pPrime   =yPrime(1)
  RhoD     =y(2)
  Te       =y(3)
  Theta    =y(4)
  RhoV     =y(5)
  RhoVL    =y(6)

  RelHum=Rpar(2)
  Theta0=Rpar(3)
  fac=Rpar(4)

  Delta(1)=pPrime+(RhoD+RhoVL)*Grav
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)
  Delta(3)=(289.0*exp(2.47e6/1004./Te*(RhoVL-RhoV)/(RhoD+RhoVL)))-Theta
  Delta(4)=Te-Theta*(p/p0)**(Rd/Cpd)
  pvs=SatVapor(Te) 
  Delta(5)=RhoV-MIN(RhoVL,pvs/(Rv*Te))
  Delta(6)=RhoVL-fac*9.0e-3*(RhoD+RhoVL)
END SUBROUTINE RES1


SUBROUTINE Res2(t,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: t,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Theta0
  REAL(RealKind) :: p,RhoD,RhoV,RhoVL,qv,qt,Te,Theta,ThetaM,ThetaE,pPrime
  REAL(RealKind) :: e,qt0,RelHum,pvs,fac

  p        =y(1)
  pPrime   =yPrime(1)
  RhoD     =y(2)
  Te       =y(3)
  Theta    =y(4)
  RhoV     =y(5)
  RhoVL    =y(6)

  RelHum=Rpar(2)
  Theta0=Rpar(3)
  fac=Rpar(4)

  Delta(1)=pPrime+(RhoD+RhoVL)*Grav
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)
  Delta(3)=((Theta0+(t-840.)**(1./3.))*exp(LLv(Te)/Cpd/Te*(RhoVL-RhoV)/(RhoD+RhoVL)))-Theta
  Delta(4)=Te-Theta*(p/p0)**(Rd/Cpd) 
  pvs=SatVapor(Te)
  Delta(5)=RhoV-RhoVL !MIN(RhoVL,pvs/(Rv*Te))
  Delta(6)=RhoVL-fac*9.0e-3*(RhoD+RhoVL)

END SUBROUTINE RES2

SUBROUTINE Res3(t,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: t,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Theta0
  REAL(RealKind) :: p,RhoD,RhoV,RhoVL,qv,qt,Te,Theta,ThetaM,ThetaE,pPrime
  REAL(RealKind) :: e,qt0,RelHum,pvs,fac
  REAL(RealKind) :: Const0,Gamm0 

  p        =y(1)
  pPrime   =yPrime(1)
  RhoD     =y(2)
  Te       =y(3)
  RhoV     =y(4)

  RelHum=Rpar(2)
  Theta0=Rpar(3)
  fac=Rpar(4)
 

  Delta(1)=pPrime+(RhoD+RhoV)*Grav        ! dp/dzProf = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)         ! p = rho*R*T
  Delta(3)=Theta0-Te*(p0/p)**(Rd/Cpd)     ! Potential temperature
  Delta(4)=RhoV-RelHum*SatVapor(Te)/(Rv*Te) ! RhoV

END SUBROUTINE Res3

SUBROUTINE ResDry(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
 
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*Rd*Rho                     ! p = rho*R*T
  Delta(3)=Theta-Te*(p0/p)**(Rd/Cpd)        ! Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)

END SUBROUTINE ResDry

SUBROUTINE ResMoistRelHum(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: RelHu,RelHo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho,RhoD,RhoV,RelH
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  RelH     =y(6)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  RelHu=Rpar(5)
  RelHo=Rpar(6)
  WRITE(*,*) 'y     ',y(1:6)
  WRITE(*,*) 'zu,zo',zu,zo
  WRITE(*,*) 'tu,to',tu,to
  WRITE(*,*) 'RelHu,RelHo',RelHu,RelHo
 
  RhoD=Rho-RhoV
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
  Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
  Delta(5)=RhoV-RelH*SatVapor(Te)/(Rv*Te) ! RhoV
  Delta(6)=RelH-((z-zu)*RelHo+(zo-z)*RelHu)/(zo-zu)
  WRITE(*,*) 'p     ',p,Te*Rd*RhoD,Te*Rv*RhoV
  WRITE(*,*) 'Delta ',delta(1:6)

END SUBROUTINE ResMoistRelHum

SUBROUTINE ResMoistRelHum2(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: RelHu,RelHo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho,RhoD,RhoV,RelH
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Rho      =y(3)
  RelH     =y(4)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  RelHu=Rpar(5)
  RelHo=Rpar(6)
  WRITE(*,*) 'HUM2 y     ',y(1:4)
  WRITE(*,*) 'zu,zo',zu,zo
  WRITE(*,*) 'tu,to',tu,to
  WRITE(*,*) 'RelHu,RelHo',RelHu,RelHo
 
  RhoD=Rho-RelH*SatVapor(Te)/(Rv*Te)
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*Rd*RhoD-RelH*SatVapor(Te)          ! p = rho*R*T
  Delta(3)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
  Delta(4)=RelH-((z-zu)*RelHo+(zo-z)*RelHu)/(zo-zu)
  WRITE(*,*) 'p     ',p,Te*Rd*RhoD,Te*Rv*RhoV
  WRITE(*,*) 'HUM2 Delta ',delta(1:4)

END SUBROUTINE ResMoistRelHum2

SUBROUTINE ResMoistRhoV(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: RelHu,RelHo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho,RhoD,RhoV,RelH
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  RelH     =y(6)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  RelHu=Rpar(5)
  RelHo=Rpar(6)
 
  RhoD=Rho-RhoV
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
  Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
  Delta(5)=RhoV-RelH*SatVapor(Te)/(Rv*Te) ! RhoV
  Delta(6)=RelH-((z-zu)*RelHo+(zo-z)*RelHu)/(zo-zu)

END SUBROUTINE ResMoistRhoV

SUBROUTINE JacResMoistRhoV(z,y,yprime,pd,cj,rpar,ipar)

  INTEGER :: ires,ipar(:)
  REAL(RealKind) :: z,y(:),yprime(:),cj,rpar(:)
  REAL(RealKind) :: pd(:,:)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: RelHu,RelHo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho,RhoD,RhoV,RelH
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
! PD=DG/DY+CJ*DG/DYPRIME
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  RelH     =y(6)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  RelHu=Rpar(5)
  RelHo=Rpar(6)
 
  RhoD=Rho-RhoV
! Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
! g(1)=yp(1)+y(4)*Grav  
  pd(1,1)=cj
  pd(1,4)=Grav
! Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
! g(2)=y(1)-y(2)*(Rd*(y(4)-y(5))+Rv*y(5))          ! p = rho*R*T
  pd(2,1)=1.0d0
  pd(2,2)=-(Rd*RhoD+Rv*RhoV)
  pd(2,4)=-Te*Rd
  pd(2,5)=-Te*(Rv-Rd)*Rhov
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
! Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
! g(3)=y(3)-y(2)*(p0/y(1))**(Rd*(y(4)-y(5)+Rv*y(5))/Cpd*(y(4)-y(5))+Cpv*y(5))*((y(4)-y(5))+y(5)*Rd/Rv)/y(4)
  pd(3,1)=-Te*(RhoD+RhoV*Rd/Rv)/Rho*(p0/p)**((Rm/Cpml)-1.0d0)*(Rm/Cpml)*(-p0/(p*p))
  pd(3,2)=-(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho
  pd(3,3)=1.0d0
! Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
! Delta(5)=RhoV-RelH*SatVapor(Te)/(Rv*Te) ! RhoV
! Delta(6)=RelH-((z-zu)*RelHo+(zo-z)*RelHu)/(zo-zu)

END SUBROUTINE JacResMoistRhoV

SUBROUTINE ResMoistQV(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: Tu,To,zu,zo
  REAL(RealKind) :: QVu,QVo
  REAL(RealKind) :: p,pPrime,Te,Theta,Rho,RhoD,RhoV,QV
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  QV       =y(6)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  QVu=Rpar(5)
  QVo=Rpar(6)
 
  RhoD=Rho-RhoV
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
  Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
  Delta(5)=RhoV-Rho*QV
  Delta(6)=QV-((z-zu)*QVo+(zo-z)*QVu)/(zo-zu)

END SUBROUTINE ResMoistQV

SUBROUTINE ResMoistQVPot(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: ThetaDu,ThetaDo,zu,zo
  REAL(RealKind) :: QVu,QVo
  REAL(RealKind) :: uWindo,uWindu,vWindo,vWindu,wWindo,wWindu
  REAL(RealKind) :: uWind,vWind,wWind
  REAL(RealKind) :: p,pPrime,Te,Theta,ThetaD,Rho,RhoD,RhoV,QV
  REAL(RealKind) :: Rm,Cpml

! p      y(1)
! T      y(2)
! \Theta y(3)
! \rho   y(4)
  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  QV       =y(6)
  ThetaD   =y(7)
  uWind    =y(8)
  vWind    =y(9)
  wWind    =y(10)

  zu=Rpar(1)
  zo=Rpar(2)
  ThetaDu=Rpar(3)
  ThetaDo=Rpar(4)
  QVu=Rpar(5)
  QVo=Rpar(6)
  uWindu=Rpar(7)
  uWindo=Rpar(8)
  vWindu=Rpar(9)
  vWindo=Rpar(10)
  wWindu=Rpar(11)
  wWindo=Rpar(12)
 
  RhoD=Rho-RhoV
  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
  Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
  Delta(4)=ThetaD-((z-zu)*ThetaDo+(zo-z)*ThetaDu)/(zo-zu)
  Delta(5)=RhoV-Rho*QV
  Delta(6)=QV-((z-zu)*QVo+(zo-z)*QVu)/(zo-zu)
  Delta(7)=ThetaD-Te*(p0/p)**(Rd/Cpd)                              ! Dry Potential Temperature
  Delta(8)=uWind-((z-zu)*uWindo+(zo-z)*uWindu)/(zo-zu)
  Delta(9)=vWind-((z-zu)*vWindo+(zo-z)*vWindu)/(zo-zu)
  Delta(10)=wWind-((z-zu)*wWindo+(zo-z)*wWindu)/(zo-zu)

END SUBROUTINE ResMoistQVPot

SUBROUTINE JacResMoistQVPot(z,y,yprime,pd,cj,rpar,ipar)

  INTEGER :: ires,ipar(:)
  REAL(RealKind) :: z,y(:),yprime(:),pd(:,:),cj,rpar(:)

  REAL(RealKind) :: ThetaDu,ThetaDo,zu,zo
  REAL(RealKind) :: QVu,QVo
  REAL(RealKind) :: p,pPrime,Te,Theta,ThetaD,Rho,RhoD,RhoV,QV
  REAL(RealKind) :: Rm,Cpml
  REAL(RealKind) :: DRmDRho,DRmDRhoV,DCpmlDRho,DCpmlDRhoV
  REAL(RealKind) :: DRhoDDRho,DRhoDDRhoV
  REAL(RealKind) :: t1,Dt1Dp
  REAL(RealKind) :: t2,Dt2Dp,Dt2DRho,Dt2DRhoV
  REAL(RealKind) :: t3,Dt3Dp,Dt3DRho,Dt3DRhoV
  REAL(RealKind) :: t4,Dt4Dp,Dt4DRho,Dt4DRhoV
  REAL(RealKind) :: Res2,DRes2Dp,DRes2DRho,DRes2DRhoV,DRes2DTe
  REAL(RealKind) :: Res3,DRes3Dp,DRes3DRho,DRes3DRhoV,DRes3DTe,DRes3DTheta
  REAL(RealKind) :: Res7,DRes7Dp,DRes7DTe,DRes7DThetaD

  p        =y(1)
  pPrime   =yPrime(1)
  Te       =y(2)
  Theta    =y(3)
  Rho      =y(4)
  RhoV     =y(5)
  QV       =y(6)
  ThetaD   =y(7)

  zu=Rpar(1)
  zo=Rpar(2)
  ThetaDu=Rpar(3)
  ThetaDo=Rpar(4)
  QVu=Rpar(5)
  QVo=Rpar(6)

! PD=DG/DY+CJ*DG/DYPRIME
  RhoD=Rho-RhoV
  DRhoDDRho=1.0d0
  DRhoDDRhov=-1.0d0
!  Delta(1)=pPrime+Rho*Grav                 ! dp/dzProf = - rho*g
!  g(1)=yp(1)+y(4)*Grav  
  pd(1,1)=cj
  pd(1,4)=Grav

! Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  DRes2Dp=1.0d0
  DRes2DTe=-(Rd*RhoD+Rv*RhoV)
  DRes2DRho=-Te*Rd*DRhodDRho
  DRes2DRhoV=-Te*Rd*DRhoDDRhoV-Te*Rv
  pd(2,1)=DRes2Dp
  pd(2,2)=DRes2DTe
  pd(2,4)=DRes2DRho
  pd(2,5)=DRes2DRhoV

  Rm=Rd*RhoD+Rv*RhoV
  DRmDRho=Rd*DRhoDDRho
  DRmDRhoV=Rv
  Cpml=Cpd*RhoD+Cpv*RhoV
  DCpmlDRho=Cpd*DRhoDDRho
  DCpmlDRhoV=Cpv
  t1=p0/p
  Dt1Dp=-p0/(p*p)
  t2=Rm/Cpml
  Dt2DRho=(DRmDRho*Cpml-Rm*DCpmlDRho)/(Cpml*Cpml)
  Dt2DRhoV=(DRmDRhoV*Cpml-Rm*DCpmlDRhoV)/(Cpml*Cpml)
  t3=t1**t2
  Dt3Dp=t3*(Dt1Dp/t1)*t2
  Dt3DRho=t3*log(t1)*Dt2DRho
  Dt3DRhoV=t3*log(t1)*Dt2DRhoV
  t4=(RhoD+RhoV*Rv/Rd)/Rho
  Dt4DRho=(DRhoDDRho*Rho-(RhoD+RhoV*Rv/Rd))/(Rho*Rho)
  Dt4DRhoV=(DRhoDDRhoV+Rv/Rd)/Rho
! Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho        ! Density Potential temperature
  Res3=Theta-Te*t3*t4
  DRes3Dp=-Te*Dt3Dp*t4
  DRes3DTe=-t3*t4
  DRes3DTheta=1.0d0
  DRes3DRho=-Te*(Dt3Drho*t4+t3*Dt4DRho)
  DRes3DRhoV=-Te*(Dt3DrhoV*t4+t3*Dt4DRhoV)
  pd(3,1)=DRes3Dp
  pd(3,2)=DRes3DTe
  pd(3,3)=DRes3DTheta
  pd(3,4)=DRes3DRho
  pd(3,5)=DRes3DRhoV

! Delta(4)=ThetaD-((z-zu)*ThetaDo+(zo-z)*ThetaDu)/(zo-zu)
  pd(4,7)=1.0d0

! Delta(5)=RhoV-Rho*QV
  pd(5,4)=-QV
  pd(5,5)=1.0d0
  pd(5,6)=-Rho

! Delta(6)=QV-((z-zu)*QVo+(zo-z)*QVu)/(zo-zu)
  pd(6,6)=1.0d0

! Delta(7)=ThetaD-Te*(p0/p)**(Rd/Cpd)                              ! Dry Potential Temperature
  t1=p0/p
  Dt1Dp=-p0/(p*p)
  t2=t1**(Rd/Cpd)
  Dt2Dp=t2*Dt1Dp/t1*(Rd/Cpd)
  Res7=ThetaD-Te*t2
  DRes7Dp=-Te*Dt2Dp
  DRes7DTe=-t2
  DRes7DThetaD=1.0d0
  pd(7,1)=DRes7Dp
  pd(7,2)=DRes7DTe
  pd(7,7)=DRes7DThetaD

  pd(8,8)=1.0d0
  pd(9,9)=1.0d0
  pd(10,10)=1.0d0
END SUBROUTINE JacResMoistQVPot

SUBROUTINE ResMoistQTPotEquiv(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  REAL(RealKind) :: p,rho,T,r_t,r_v,rho_qv,Theta_e,pPrime
  REAL(RealKind) :: Rho_d,p_d,T_C,p_vs,L,r_t0,Theta_e0
  REAL(RealKind) :: a,b

  ! p      y(1)
  ! T      y(2)
  ! \Theta y(3)
  ! \rho   y(4)
  p = y(1)
  rho = y(2)
  T = y(3)
  r_t = y(4)
  r_v = y(5)
  rho_qv = y(6)
  Theta_e = y(7)
  pPrime = yPrime(1)
  WRITE(*,*) 'y      ',y(1:7)

  Theta_e0=Rpar(3)
  r_t0=Rpar(5)

  rho_d = rho / (1 + r_t)
  p_d = Rd * rho_d * T
  T_C = T - 273.15d0
  p_vs = 611.2d0 * EXP(17.62d0 * T_C / (243.12d0 + T_C))
  L = L00 - (cpl - cpv) * T
  Delta(1) = pPrime + grav * rho
  Delta(2) = p - (Rd * rho_d + Rv * rho_qv) * T
  Delta(3) = Theta_e - T * (p_d / p0)**(-Rd/  &
       (cpd + cpl * r_t)) * EXP(L * r_v / ((cpd + cpl * r_t) * T))
  Delta(4) = r_t - r_t0
  Delta(5) = rho_qv-rho_d * r_v
  Delta(6) = Theta_e - Theta_e0
  a = p_vs / (Rv * T) - rho_qv
  b = rho - rho_qv - rho_d
  Delta(7)=a+b-SQRT(a*a+b*b)
  WRITE(*,*) 'Delta  ',Delta(1:7)

END SUBROUTINE ResMoistQTPotEquiv

SUBROUTINE InputRes(y,yprime,rpar,ipar,PreStart,RhStart,TeStart,ThDStart,QVStart,&
                    ThEStart,rtStart,uWindStart,vWindStart,wWindStart)

  REAL(RealKind) :: y(:),yprime(:),rpar(:)
  INTEGER :: ipar(:)
  REAL(RealKind), OPTIONAL :: PreStart
  REAL(RealKind), OPTIONAL :: RHStart
  REAL(RealKind), OPTIONAL :: TeStart
  REAL(RealKind), OPTIONAL :: ThDStart
  REAL(RealKind), OPTIONAL :: QVStart
  REAL(RealKind), OPTIONAL :: ThEStart
  REAL(RealKind), OPTIONAL :: rtStart
  REAL(RealKind), OPTIONAL :: uWindStart
  REAL(RealKind), OPTIONAL :: vWindStart
  REAL(RealKind), OPTIONAL :: wWindStart

  REAL(RealKind) :: Theta0,p,Te,e
  CHARACTER*20:: par
  SELECT CASE(ProfileType)
    CASE('Dry')
      neq=4
      PrePos=1
      TempPos=2
      ThPos=3
      RhoPos=4
      IF (PRESENT(preStart)) THEN
        y(1)=preStart
      ELSE
        y(1)=1.0d5 
      END IF  
      IF (PRESENT(teStart)) THEN
        y(2)=teStart
        rpar(3)=teStart
        rpar(4)=teStart
      ELSE
        y(2)=293.15
        rpar(3)=293.15
        rpar(4)=293.15
      END IF  
      yPrime(1)=-y(RhoPos)*Grav
      y(3)=y(2)*(p0/y(1))**(Rd/Cpd)
      y(4)=y(1)/(Rd*y(2))
      rpar(1)=0.0d0
      rpar(2)=1.0d0
    CASE('MoistRelHum','MoistQV')
      neq=6
      PrePos=1
      TempPos=2
      ThPos=3
      RhoPos=4
      RhoVPos=5
      WRITE(*,*) 'PRESENT',PRESENT(preStart),PRESENT(teStart),PRESENT(RHStart),PRESENT(QVStart)
      WRITE(*,*) 'PRESENT',preStart,teStart,RHStart
      IF (PRESENT(preStart)) THEN
        y(PrePos)=preStart
      ELSE
        IF (y(PrePos)==0.0d0) THEN
          y(PrePos)=1.0d5 
        END IF  
      END IF  
      IF (PRESENT(teStart)) THEN
        y(2)=teStart
        rpar(3)=teStart
        rpar(4)=teStart
      ELSE
        y(2)=293.15
        rpar(3)=293.15
        rpar(4)=293.15
      END IF  
      IF (PRESENT(RHStart)) THEN
        y(6)=RHStart
        rpar(5)=RHStart
        rpar(6)=RHStart
      END IF  
      IF (PRESENT(QVStart)) THEN
        y(6)=QVStart
        rpar(5)=QVStart
        rpar(6)=QVStart
      END IF  
      IF (y(3)==0.0d0) THEN
        y(3)=y(2)*(p0/y(1))**(Rd/Cpd)
      END IF  
      IF (y(4)==0.0d0) THEN
        y(4)=y(1)/(Rd*y(2))
        WRITE(*,*) 'RhoD in Start',y(1),y(2),y(4),y(4)*Rd*y(2),Rd
      END IF  
      IF (PRESENT(RHStart)) THEN
        y(5)=y(6)*SatVapor(y(2))/(Rv*y(2))
      ELSE IF (PRESENT(QVStart).AND.y(5)==0.0d0) THEN
        y(5)=y(6)
      END IF  
      y(4)=(y(1)-Rv*y(2)*y(5))/(Rd*y(2))
      y(4)=y(4)+y(5)
      yPrime(1)=-y(RhoPos)*Grav
      rpar(1)=0.0d0
      rpar(2)=1.0d0
      WRITE(*,*) 'y',y(1:6)
    CASE('MoistRelHum2')
      neq=4
      PrePos=1
      TempPos=2
      RhoPos=3
      RhoVPos=4
      WRITE(*,*) 'PRESENT',PRESENT(preStart),PRESENT(teStart),PRESENT(RHStart),PRESENT(QVStart)
      WRITE(*,*) 'PRESENT',preStart,teStart,RHStart
      IF (PRESENT(preStart)) THEN
        y(PrePos)=preStart
      ELSE
        IF (y(PrePos)==0.0d0) THEN
          y(PrePos)=1.0d5 
        END IF  
      END IF  
      IF (PRESENT(teStart)) THEN
        y(2)=teStart
        rpar(3)=teStart
        rpar(4)=teStart
      ELSE
        y(2)=293.15
        rpar(3)=293.15
        rpar(4)=293.15
      END IF  
      IF (PRESENT(RHStart)) THEN
        y(4)=RHStart
        rpar(5)=RHStart
        rpar(6)=RHStart
      END IF  
      IF (y(3)==0.0d0) THEN
        y(3)=y(1)/(Rd*y(2))
        WRITE(*,*) 'RhoD in Start',y(1),y(2),y(4),y(4)*Rd*y(2),Rd
      END IF  
      y(3)=(y(1)-Rv*y(2)*y(5))/(Rd*y(2))
      y(3)=y(3)+y(4)
      yPrime(1)=-y(RhoPos)*Grav
      rpar(1)=0.0d0
      rpar(2)=1.0d0
      WRITE(*,*) 'y',y(1:4)
    CASE('MoistQVPot')
      neq=10
      PrePos=1
      TempPos=2
      ThPos=3
      RhoPos=4
      RhoVPos=5
      QvPos=6
      ThDPos=7
      uWindPos=8
      vWindPos=9
      wWindPos=10
      IF (PRESENT(preStart)) THEN
        y(PrePos)=preStart
      ELSE
        IF (y(PrePos)==0.0d0) THEN
          y(PrePos)=1.0d5 
        END IF  
      END IF  
      IF (PRESENT(teStart)) THEN
        y(TempPos)=teStart
        rpar(3)=teStart
        rpar(4)=teStart
        y(ThDPos)=teStart
      ELSE
        y(TempPos)=293.15
        y(ThDPos)=293.15
        rpar(3)=293.15
        rpar(4)=293.15
      END IF  
      IF (PRESENT(ThDStart)) THEN
        y(ThDPos)=ThDStart
        y(ThPos)=ThDStart
        rpar(3)=ThDStart
        rpar(4)=ThDStart
        y(TempPos)=y(ThDPos)*(y(PrePos)/p0)**(Rd/Cpd)
        y(RhoPos)=y(PrePos)/(Rd*y(TempPos))
      END IF
      IF (PRESENT(RHStart)) THEN
        y(6)=RHStart
        rpar(5)=RHStart
        rpar(6)=RHStart
      END IF  
      IF (PRESENT(QVStart)) THEN
        y(QvPos)=QVStart
        rpar(5)=QVStart
        rpar(6)=QVStart
      END IF  
      IF (PRESENT(uWindStart)) THEN
        y(uWindPos)=uWindStart
        rpar(7)=uWindStart
        rpar(8)=uWindStart
      END IF  
      IF (PRESENT(vWindStart)) THEN
        y(vWindPos)=vWindStart
        rpar(9)=vWindStart
        rpar(10)=vWindStart
      END IF  
      IF (PRESENT(wWindStart)) THEN
        y(wWindPos)=wWindStart
        rpar(11)=wWindStart
        rpar(12)=wWindStart
      END IF  
      IF (y(ThPos)==0.0d0) THEN
        y(ThPos)=y(TempPos)*(p0/y(PrePos))**(Rd/Cpd)
      END IF  
      IF (y(RhoPos)==0.0d0) THEN
        y(RhoPos)=y(PrePos)/(Rd*y(2))

      END IF  
      IF (PRESENT(RHStart)) THEN
        y(5)=y(6)*SatVapor(y(2))/(Rv*y(2))
      ELSE IF (PRESENT(QVStart).AND.y(RhoVPos)==0.0d0) THEN
        y(RhoVPos)=y(QvPos)*y(RhoPos)
      END IF  
      yPrime(1)=-y(RhoPos)*Grav
      rpar(1)=0.0d0
      rpar(2)=1.0d0
    CASE('MoistQTPotEquiv')
      y=0.0d0
      yPrime=0.0d0
      neq=7
      PrePos=1
      RhoPos=2
      TPos=3
      rtPos=4
      rvPos=5
      RhoVPos=6
      ThEPos=7
      IF (PRESENT(preStart)) THEN
        y(PrePos)=preStart
      ELSE
        y(PrePos)=1.0d5 
      END IF  
      y(TPos)=TheStart
      y(RhoPos)=y(prePos)/(Rd*y(TPos))
      y(rtPos)=rtStart
      y(rvPos)=rtStart
      y(RhoVPos)=y(RhoPos)*y(rvPos)
      y(ThEPos)=TheStart
      rpar(3)=TheStart
      rpar(5)=rtStart
    CASE('MoistRhoV')
    CASE DEFAULT  
      y=0.0d0
      yPrime=0.0d0
  END SELECT  

END SUBROUTINE InputRes

SUBROUTINE OutputRes(nzProf,zMProf,c,FileNameOut)


  INTEGER :: nzProf
  REAL(RealKind) ::zMProf(:),c(:,:)
  CHARACTER(*) :: FileNameOut

  INTEGER :: i
  REAL(8) :: Rho,RhoD,RhoV,RhoL,Rm,Cpml,L,p,pD,T,Te,Theta,r_t,r_v,ThE

  OPEN(unit=OutputUnit,FILE=TRIM(FileNameOut)//'.prof')
  SELECT CASE(ProfileType)
    CASE('Dry')
      WRITE(OutputUnit,*) 'RhoProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'ThDensProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,ThPos)
      END DO
      WRITE(OutputUnit,*) 'TempProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,TempPos)
      END DO
      WRITE(OutputUnit,*) 'PreProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,PrePos)
      END DO
    CASE('MoistRelHum','MoistRhoV')
      WRITE(OutputUnit,*) 'RhoProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'ThDensProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,ThPos)
      END DO
      WRITE(OutputUnit,*) 'TempProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,TempPos)
      END DO
      WRITE(OutputUnit,*) 'RhoVProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,RhoVPos)/c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'PreProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,PrePos)
      END DO
    CASE('MoistQTPotEquiv')
      WRITE(OutputUnit,*) 'RhoProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'PreProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,prePos)
      END DO
      WRITE(OutputUnit,*) 'TempProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,TPos)
      END DO
      WRITE(OutputUnit,*) 'ThProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        Rho=c(i,RhoPos)
        RhoV=c(i,RhoVPos)
        RhoD =Rho/(1.0d0+c(i,rtPos))
        RhoL=Rho-RhoD-RhoV
        Rm=Rd*RhoD+Rv*RhoV
        Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL
        p=c(i,PrePos)
        Te=c(i,TPos)
        Theta=Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rv/Rd)/Rho
        WRITE(OutputUnit,*) zMProf(i),Theta
      END DO
      WRITE(OutputUnit,*) 'QvProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        WRITE(OutputUnit,*) zMProf(i),c(i,RhoVPos)/c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'QcProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        Rho=c(i,RhoPos)
        RhoV=c(i,RhoVPos)
        RhoD =Rho/(1.0d0+c(i,rtPos))
        RhoL=Rho-RhoD-RhoV
        WRITE(OutputUnit,*) zMProf(i),RhoL/Rho
      END DO  
      WRITE(OutputUnit,*) 'ThEProf'
      WRITE(OutputUnit,*) nzProf,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nzProf
        Rho=c(i,RhoPos)
        RhoV=c(i,RhoVPos)
        RhoD =Rho/(1.0d0+c(i,rtPos))
        RhoL=Rho-RhoD-RhoV
        T=c(i,TPos)
        pD=RhoD*T*Rd
        r_t=c(i,rtPos)
        r_v=c(i,rvPos)
        L = L00 - (cpl - cpv) * T
        ThE=T*(pD/p0)**(-Rd/  &
            (cpd + cpl * r_t)) * EXP(L * r_v / ((cpd + cpl * r_t) * T))
        WRITE(OutputUnit,*) zMProf(i),ThE
      END DO  
  END SELECT    
  CLOSE(OutputUnit)

END SUBROUTINE OutputRes


SUBROUTINE CreateProfile()

  INTEGER :: i, iz
  REAL(RealKind) :: PotTemp(NumSounding)
  REAL(RealKind) :: QV(NumSounding)
  REAL(RealKind) :: RhoL, RhoV, VaporPres
  REAL(RealKind) :: dTemp, dRelHum, dQV

! Check Inputdata for unambiguousness
  IF (PotTemp0.GT.0.d0 .AND. Temp0.GT.0.d0) &
&   STOP 'ERROR: Startvalue for potential and absolut temperature at surface given; only one of them possible'

  IF (dTemp_FT.LT.0.d0) &
&   STOP 'ERROR: No temperature gradient defined for free troposphere'

  IF (RelHum0.GT.0.d0 .AND. QV0.GT.0.d0) &
&   STOP 'ERROR: Startvalue for relative and specific humidity at surface given; only one of them possible'

  IF (dRH_FT.GT.-9.d8 .AND. dQV_FT.GT.-9.d8) &
&   STOP 'ERROR: Gradient for relative and specific humidity in free troposphere given; only one of them possible'

  IF (Pres0.LT.1.d4) Pres0=Pres0*1.d2
  IF (InvLayerHeight.GT.0.d0) dzProf_Inv=MAX(dzProf_Inv,dzProf_ManProf)

  Pre(1) = Pres0
  IF (PotTemp0.GT.0.d0) THEN ! Profile calculated with given potential temperature
    WRITE(*,*) 'Case 1'
    IF (PotTemp0.LT.1.d2) PotTemp0=PotTemp0+273.15
    PotTemp(1) = PotTemp0
    ProfTemp(1)=(Pres0/1000.d2)**(Rd/Cpd)*PotTemp(1)
    DO iz=2,NumSounding
      IF (Height(iz).LE.InvLayerHeight) THEN
        dTemp=dTemp_BL
      ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dzProf_Inv) THEN
        dTemp=dTemp_Inv
      ELSE IF (Height(iz).GT.InvLayerHeight+dzProf_Inv) THEN
        dTemp=dTemp_FT
      END IF
      PotTemp(iz)= PotTemp(iz-1)+dTemp*dzProf_ManProf
      ProfTemp(iz)=ProfTemp(iz-1)
      DO i=1,20 ! Itertation, because of two unknown variables
        Pre(iz) = Pre(iz-1)*EXP((-1.d0*grav*dzProf_ManProf)/(Rd*(ProfTemp(iz-1)+ProfTemp(iz))/2.d0))
        ProfTemp(iz) = PotTemp(iz)*(Pre(iz)/1000.d2)**(Rd/Cpd)
      END DO
    END DO
  ELSE IF(Temp0>0) THEN ! Profile calculated with given absolute temperature
    WRITE(*,*) 'Case 2'
    IF (Temp0.LT.1.d2) Temp0=Temp0+273.15
    ProfTemp(1)=Temp0
    DO iz=2,NumSounding
      IF (Height(iz).LE.InvLayerHeight) THEN
        dTemp=dTemp_BL
      ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dzProf_Inv) THEN
        dTemp=dTemp_Inv
      ELSE IF (Height(iz).GT.InvLayerHeight+dzProf_Inv) THEN
        dTemp=dTemp_FT
      END IF
      ProfTemp(iz)= ProfTemp(iz-1)+dTemp*dzProf_ManProf
      Pre(iz) = Pre(iz-1)*EXP((-1.d0*grav*dzProf_ManProf)/(Rd*(ProfTemp(iz-1)+ProfTemp(iz))/2.d0))
    END DO
  END IF
  IF (rt0>0.0d0.AND.PotTempE0>0.0d0) THEN
    WRITE(*,*) 'Case 3'
    DO iz=1,NumSounding
      rt(iz)=rt0
      ThE(iz)=PotTempE0
    END DO
    ProfileType='MoistQTPotEquiv'
  ELSE IF (RelHum0>0.d0 .OR. QV0>0.d0) THEN ! Moist case
    WRITE(*,*) 'Case 3'
!   CaseProfile=Moist
    RhoL = Pre(1)/(Rd*ProfTemp(1))
    IF (RelHum0.GT.0.d0) THEN
      RH(1)=RelHum0
      VaporPres = RH(1)/1.d2*SatVapor(ProfTemp(1))
      RhoV = VaporPres/(Rv*ProfTemp(1))
      QV(1)=RhoV/(RhoV+RhoL)
    ELSE IF (QV0>0.d0) THEN
      QV(1)=QV0
      RhoV = RhoL*(QV(1)/(1.d0-QV(1)))
      VaporPres = RhoV*Rv*ProfTemp(1)
      RH(1) = VaporPres/SatVapor(ProfTemp(1))*1.d2
    END IF


    IF (dRH_FT.GE.0.d0) THEN
      DO iz=2,NumSounding
        IF (Height(iz)<=InvLayerHeight) THEN
          dRelHum=dRH_BL
        ELSE IF (Height(iz)>InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dzProf_Inv) THEN
          dRelHum=dRH_Inv
        ELSE IF (Height(iz)>InvLayerHeight+dzProf_Inv) THEN
          dRelHum=dRH_FT
        END IF
        RH(iz)= RH(iz-1)+dRelHum*dzProf_ManProf
      END DO
    ELSE IF (dQV_FT>=0.d0) THEN
      DO iz=2,NumSounding
        RhoL = Pre(iz)/(Rd*ProfTemp(iz))
        IF (Height(iz).LE.InvLayerHeight) THEN
          dQV=dQV_BL
        ELSE IF (Height(iz)>InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dzProf_Inv) THEN
          dQV=dQV_Inv
        ELSE IF (Height(iz)>InvLayerHeight+dzProf_Inv) THEN
          dQV=dQV_FT
        END IF
        QV(iz)= QV(iz-1)+dQV*dzProf_ManProf
        RhoV = RhoL*(QV(iz)/(1.d0-QV(iz)))
        VaporPres = RhoV*Rv*ProfTemp(iz)
        RH(iz) = VaporPres/SatVapor(ProfTemp(iz))*1.d2
      END DO
    ELSE
      STOP 'ERROR: No gradient for humidity defined for free troposphere'
    END IF

    RH=RH/1.d2
    ProfileType='MoistRelHum2'
    WRITE(*,*) '===== Wet profile calculated'
  ELSE ! Dry case
    WRITE(*,*) '===== Dry profile calculated'
    ProfileType='Dry'
  END IF

END SUBROUTINE CreateProfile

SUBROUTINE ComputeProfile(c,Height,Pre,Temp,ThD,RH,QV,ThE,rt,uWind,vWind,wWind)

  USE Kind_Mod

  REAL(RealKind), ALLOCATABLE :: c(:,:)
  REAL(RealKind) :: Height(:)
  REAL(RealKind), OPTIONAL :: Pre(:)
  REAL(RealKind), OPTIONAL :: Temp(:)
  REAL(RealKind), OPTIONAL :: ThD(:)
  REAL(RealKind), OPTIONAL :: RH(:)
  REAL(RealKind), OPTIONAL :: QV(:)
  REAL(RealKind), OPTIONAL :: ThE(:)
  REAL(RealKind), OPTIONAL :: rt(:)
  REAL(RealKind), OPTIONAL :: uWind(:)
  REAL(RealKind), OPTIONAL :: vWind(:)
  REAL(RealKind), OPTIONAL :: wWind(:)

  INTEGER :: iHeight,iz
  INTEGER, PARAMETER :: lrw=1000
  INTEGER, PARAMETER :: liw=1000
  INTEGER  INFO(15), IDID, IWORK(liw), IPAR(100)
  REAL(RealKind) :: T, TOUT, RWORK(lrw), &
      RPAR(100), RTOL, ATOL
  REAL(RealKind) :: z,zu,zo,zOut    
  REAL(RealKind) :: TeStart     = 288.0d0  
  REAL(RealKind) :: ThDStart     = 288.0d0  
  REAL(RealKind) :: PreStart    = 1013.25d2
  REAL(RealKind) :: RelHumStart = 0.00d0
  REAL(RealKind) :: QVStart=1.d-4
  REAL(RealKind) :: ThEStart=320.0d0
  REAL(RealKind) :: rtStart=1.d-2
  REAL(RealKind) :: uWindStart=0.00d0
  REAL(RealKind) :: vWindStart=0.00d0
  REAL(RealKind) :: wWindStart=0.00d0
  REAL(RealKind) :: y(20),yPrime(20)


  t=0.0
  info(1:15)=0


  IF (PRESENT(Pre)) THEN
    PreStart=Pre(1)
  END IF
  IF (PRESENT(Temp)) THEN
    TeStart=Temp(1)
  END IF
  IF (PRESENT(ThD)) THEN
    ThDStart=ThD(1)
  END IF
  IF (PRESENT(RH)) THEN
    RelHumStart=RH(1)
  END IF
  IF (PRESENT(QV)) THEN
    QVStart=QV(1)
  END IF
  IF (PRESENT(ThE)) THEN
    ThEStart=ThE(1)
  END IF
  IF (PRESENT(rt)) THEN
    rtStart=rt(1)
  END IF
  IF (PRESENT(uWind)) THEN
    uWindStart=uWind(1)
  END IF
  IF (PRESENT(vWind)) THEN
    vWindStart=vWind(1)
  END IF
  IF (PRESENT(wWind)) THEN
    wWindStart=wWind(1)
  END IF


  rtol = 1.0d-5
  atol = 1.0d-6

  tout=1.0d0
  yprime=0.0d0
  y=0.0d0
  WRITE(*,*) 'ProfileType ',ProfileType,PRESENT(RH),PRESENT(QV),PRESENT(rt)
  IF (ProfileType/='Dry') THEN
    IF (PRESENT(RH)) THEN
      RH(1)=0.01
      RH(2)=0.01
      RelHumStart=RH(1)
      WRITE(*,*) 'RH Case'
      CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=RelHumStart)
      WRITE(*,*) 'RH Case',y
    ELSE IF (PRESENT(QV)) THEN
      QvStart=QV(1) !9.d-4 !QV(1) !2.d-4
      CALL InputRes(y,yprime,rpar,ipar,PreStart,ThDStart=ThDStart,QVStart=QvStart,&
                    uWindStart=uWindStart,vWindStart=vWindStart,wWindStart=wWindStart)
    ELSE IF (PRESENT(rt)) THEN
      rtStart=1.d-8*rt(1) !9.d-4 !QV(1) !2.d-4
      CALL InputRes(y,yprime,rpar,ipar,ThEStart=ThEStart,rtStart=rtStart,&
                    uWindStart=uWindStart,vWindStart=vWindStart,wWindStart=wWindStart)
    END IF  
    DO  
      rwork=0.0d0
      info=0.0d0
      INFO(11)=1 ! DDASSL computes consisten initial values
      INFO(5)=0
      z=Height(1)
      zOut=Height(1)+1.d-6
      WRITE(*,*) 'DDASSL 1',NEQ,Y
      CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOut, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      WRITE(*,*) 'nach DDASSL 1',NEQ
      IF (PRESENT(RH)) THEN
        IF (rpar(5)==RH(1)) THEN
          EXIT
        END IF  
        RelHumStart=MIN(rpar(5)+0.01d0,RH(1))
        CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=RelHumStart)
      ELSE IF (PRESENT(QV)) THEN
        IF (rpar(5)==QV(1)) THEN
          EXIT
        END IF  
        QVStart=MIN(rpar(5)+9.d-5,QV(1))
        rpar(5)=QVStart
        rpar(6)=QVStart
        y(QvPos)=QVStart
!       CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,QVStart=QVStart)
      ELSE IF (PRESENT(rt)) THEN
        IF (rpar(5)==rt(1)) THEN
          EXIT
        END IF  
        rtStart=MIN(rpar(5)+9.d-7,rt(1))
        y(rtPos)=rtStart
        rpar(5)=rtStart
        rpar(6)=rtStart
      END IF  
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0
    info(5)=0
    info(11)=1 ! DDASSL computes consisten initial values
    z=Height(1)
    zOut=Height(1)+1.d-6
    WRITE(*,*) 'DDASSL 2'
    CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOut, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF 
  IF (ALLOCATED(c)) THEN
    DEALLOCATE(c)
    ALLOCATE(c(SIZE(zMProf),neq))
  ELSE  
    ALLOCATE(c(SIZE(zMProf),neq))
  END IF  
  c=0.0d0
  iHeight=1
  iz=1
  ! Forward
  S1:DO 
    zu=Height(iHeight)
    rpar(1)=zu
    zo=Height(iHeight+1)
    rpar(2)=zo
    z=zu
    IF (PRESENT(Temp)) THEN
      rpar(3)=Temp(iHeight)
      rpar(4)=Temp(iHeight+1)
    END IF  
    IF (PRESENT(ThD)) THEN
      rpar(3)=ThD(iHeight)
      rpar(4)=ThD(iHeight+1)
    END IF  
    IF (PRESENT(RH)) THEN
      rpar(5)=RH(iHeight)
      rpar(6)=RH(iHeight+1)
    END IF  
    IF (PRESENT(QV)) THEN
      rpar(5)=QV(iHeight)
      rpar(6)=QV(iHeight+1)
    END IF  
    IF (PRESENT(rt)) THEN
      rpar(5)=rt(iHeight)
      rpar(6)=rt(iHeight+1)
    END IF  
    IF (PRESENT(uWind)) THEN
      rpar(7)=uWind(iHeight)
      rpar(8)=uWind(iHeight+1)
    END IF  
    IF (PRESENT(vWind)) THEN
      rpar(9)=vWind(iHeight)
      rpar(10)=vWind(iHeight+1)
    END IF  
    IF (PRESENT(wWind)) THEN
      rpar(11)=wWind(iHeight)
      rpar(12)=wWind(iHeight+1)
    END IF  
    S2:DO 
      IF (iz>nzProf) THEN
        EXIT S1
      END IF  
      IF (zu<zMProf(iz).AND.zMProf(iz)<=zO) THEN
        zOut=zMProf(iz)
        WRITE(*,*) 'DDASSL 3'
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        IF (zMProf(iz)==zO) THEN
          iHeight=iHeight+1
          iz=iz+1
          EXIT S2
        END IF 
        iz=iz+1
      ELSE IF (zMProf(iz)<=zu) THEN
        iz=iz+1
      ELSE 
        zOut=zo
        WRITE(*,*) 'DDASSL 4'
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        iHeight=iHeight+1
        EXIT S2
      END IF  
    END DO S2  
  END DO S1 

  ! Backward
  yprime=0.0d0
  y=0.0d0
  IF (ProfileType/='Dry') THEN
    IF (PRESENT(RH)) THEN
      RelHumStart=0.01d0
      CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=RelHumStart)
    ELSE IF (PRESENT(QV)) THEN
      QvStart=QV(1)
      CALL InputRes(y,yprime,rpar,ipar,PreStart,ThDStart=ThDStart,QVStart=QvStart,&
                    uWindStart=uWindStart,vWindStart=vWindStart,wWindStart=wWindStart)
    ELSE IF (PRESENT(rt)) THEN
      rtStart=1.d-8*rt(1) !9.d-4 !QV(1) !2.d-4
      CALL InputRes(y,yprime,rpar,ipar,ThEStart=ThEStart,rtStart=rtStart,&
                    uWindStart=uWindStart,vWindStart=vWindStart,wWindStart=wWindStart)
    END IF  
    DO  
      rwork=0.0d0
      info=0
      info(5)=0
      INFO(11)=1 ! DDASSL computes consisten initial values
      z=Height(1)
      zOut=Height(1)-1.d-6
      CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOut, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      IF (PRESENT(RH)) THEN
        IF (rpar(5)==RH(1)) THEN
          EXIT
        END IF  
        RelHumStart=MIN(rpar(5)+0.01d0,RH(1))
        CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=RelHumStart)
      ELSE IF (PRESENT(QV)) THEN
        IF (rpar(5)==QV(1)) THEN
          EXIT
        END IF  
        QVStart=MIN(rpar(5)+1.d-6,QV(1))
        rpar(5)=QVStart
        rpar(6)=QVStart
        y(QvPos)=QVStart
!       CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,QVStart=QVStart)
      ELSE IF (PRESENT(rt)) THEN
        IF (rpar(5)==rt(1)) THEN
          EXIT
        END IF  
        rtStart=MIN(rpar(5)+9.d-7,rt(1))
        y(rtPos)=rtStart
        rpar(5)=rtStart
        rpar(6)=rtStart
      END IF  
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0
    info(5)=0
    info(11)=1 ! DDASSL computes consisten initial values
    z=Height(1)
    zOut=Height(1)-1.d-6
    CALL DDASSL (RES, NEQ, T, z, YPRIME, zOut, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF  
  rwork=0.0d0
  info=0
  S3:DO 
    zu=Height(1)
    rpar(1)=zu
    zo=Height(2)
    rpar(2)=zo
    z=zu
    IF (PRESENT(Temp)) THEN
      rpar(3)=Temp(1)
      rpar(4)=Temp(2)
    END IF  
    IF (PRESENT(ThD)) THEN
      rpar(3)=ThD(1)
      rpar(4)=ThD(2)
    END IF  
    IF (PRESENT(RH)) THEN
      rpar(5)=RH(1)
      rpar(6)=RH(2)
    END IF
    IF (PRESENT(QV)) THEN
      rpar(5)=QV(1)
      rpar(6)=QV(2)
    END IF
    IF (PRESENT(rt)) THEN
      rpar(5)=rt(iHeight)
      rpar(6)=rt(iHeight+1)
    END IF  
    IF (PRESENT(uWind)) THEN
      rpar(7)=uWind(1)
      rpar(8)=uWind(2)
    END IF  
    IF (PRESENT(vWind)) THEN
      rpar(9)=vWind(1)
      rpar(10)=vWind(2)
    END IF  
    IF (PRESENT(wWind)) THEN
      rpar(11)=wWind(1)
      rpar(12)=wWind(2)
    END IF  
    iz=nzProf
    S4:DO 
      IF (iz<1) THEN
        EXIT S3
      END IF  
      IF (zu>zMProf(iz)) THEN
        zOut=zMProf(iz)
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        C(iz,uWindPos)=uWind(1)
        C(iz,vWindPos)=vWind(1)
        C(iz,wWindPos)=wWind(1)
        iz=iz-1
      ELSE IF (zMProf(iz)>=zu) THEN
        iz=iz-1
      END IF  
    END DO S4  
  END DO S3 
END SUBROUTINE ComputeProfile

END MODULE Profile_Mod

