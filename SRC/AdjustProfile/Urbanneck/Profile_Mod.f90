MODULE Profile_Mod

  USE Kind_Mod
  USE Parameter_Mod
  IMPLICIT NONE

  INTEGER :: neq
  INTEGER, PARAMETER :: Dry=1
  INTEGER, PARAMETER :: Moist=2

  INTEGER, PRIVATE :: RhoPos=0
  INTEGER, PRIVATE :: ThPos=0
  INTEGER, PRIVATE :: RhoVPos=0
  INTEGER, PRIVATE :: PrePos=0
  INTEGER, PRIVATE :: TempPos=0

  INTEGER :: nz
  REAL(RealKind), ALLOCATABLE :: zP(:)
  REAL(RealKind), ALLOCATABLE :: zM(:)
  REAL(RealKind), ALLOCATABLE :: dz(:)

  INTEGER :: InputUnit=10
  INTEGER :: OutputUnit=11
  CHARACTER(200) :: NameSounding
  CHARACTER(80) :: FileNameOut
  INTEGER :: CaseProfile
  LOGICAL :: Sounding=.FALSE.
  REAL(RealKind) :: Tropopause  = 3000.0d0
  REAL(RealKind) :: PotTemp0=-9.d9, Temp0=-9.d9, Pres0=-9.d9, RelHum0=-9.d9, QV0=-9.d9
  REAL(RealKind) :: dTemp_BL=0.d0, dRH_BL=0.d0, dQV_BL=0.d0
  REAL(RealKind) :: InvLayerHeight=-9.d9, dz_ManProf=1.d1, dz_Inv=-9.d9
  REAL(RealKind) :: dTemp_Inv=0.d0, dRH_Inv=0.d0, dQV_Inv=0.d0
  REAL(RealKind) :: dTemp_FT=0.0d0, dRH_FT=-9.d9, dQV_FT=-9.d9

  NAMELIST /AdjustCtrl/ NameSounding   &
                       ,Sounding       &
                       ,FileNameOut    &
                       ,Tropopause     &
                       ,CaseProfile    &
                       ,PotTemp0       & ! Potential temperature at surface in °C or K
                       ,Temp0          & ! Absolute air temperature at surface in °C or K
                       ,Pres0          & ! Air pressure at surface in hPa or Pa
                       ,RelHum0        & ! Relative humidity at surface in %
                       ,QV0            & ! Specific humidity at surface in kg/kg
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
                       ,dz_ManProf     & ! Gridspacing of manual set profile
                       ,dz_Inv           ! Thickness of inversion layer

  INTEGER :: NumSounding
  REAL(RealKind), ALLOCATABLE :: Height(:)
  REAL(RealKind), ALLOCATABLE :: Pre(:)
  REAL(RealKind), ALLOCATABLE :: Temp(:)
  REAL(RealKind), ALLOCATABLE :: RH(:)
  INTEGER :: NumValues

CONTAINS


SUBROUTINE ReadInput(FileName)

  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line
  REAL(RealKind) :: z0,z1,zw0,zw1
  INTEGER :: iz,k,nzz,nLines
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
  DO
    READ(InputUnit,*,END=2) Line
    IF (INDEX(Line,'#Gitter')>0) THEN
    !-------------------------------------
       READ(InputUnit,*) indata_type
       READ(InputUnit,*) z0,z1,nz
       READ(InputUnit,*) z0,z1,nz
       READ(InputUnit,*) z0,z1,nz
       EXIT
    END IF   
  END DO     
2 CONTINUE  
  REWIND(InputUnit)
  ALLOCATE(zP(0:nz))
  ALLOCATE(dz(nz))
  ALLOCATE(zM(nz))
  dz=(z1-z0)/nz
  zP(0)=z0
  DO iz=1,nz
    zP(iz)=zP(iz-1)+dz(iz)
  END DO  
  DO 
    READ(InputUnit,*,END=3) Line
    !........................................................................
    IF (INDEX(Line,'#zGrid')>0) THEN   ! set different distance z-direction
      READ(InputUnit,*) nlines
      nz=0
      zP(0)=z0
      DO k=1,nlines
        READ(InputUnit,*) zw0,zw1,nzz
        DO iz=nz+1,nz+nzz
          dz(iz)=(zw1-zw0)/nzz
          zP(iz)=zP(iz-1)+dz(iz)
        END DO
        nz=nz+nzz
      END DO
      EXIT
    END IF   
  END DO     
3 CONTINUE  
  DO iz=1,nz
    zM(iz)=0.5d0*(zP(iz-1)+zP(iz))
  END DO  
  CLOSE(UNIT=InputUnit)

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
      ALLOCATE(Temp(NumSounding))
      ALLOCATE(RH(NumSounding))
      DO iz=1,NumSounding
        READ(InputUnit,*) Height(iz),Pre(iz),Temp(iz),RH(iz),DummyChar,DummyChar
        Pre(iz)=Pre(iz)*100.0d0 ! Umrechnung hPa Pa
        Temp(iz)=Temp(iz)+273.15d0
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
  ELSE
    NumSounding=NINT((z1-z0)/dz_ManProf+1.d0)
    ALLOCATE(Height(NumSounding))
    DO iz=1,NumSounding
      Height(iz)=z0+dz_ManProf*(iz-1)
    END DO
    ALLOCATE(Pre(NumSounding))
    ALLOCATE(Temp(NumSounding))
    ALLOCATE(RH(NumSounding))
    CALL CreateProfile() ! Values defined on borders of the gridcells
  END IF
END SUBROUTINE ReadInput


FUNCTION  LLv(T)

  REAL(RealKind) :: LLv,T
  LLv=L0-(Cpl-Cpv)*(T-Tes0)
END FUNCTION  LLv

FUNCTION  SatVap(T)
  REAL(RealKind) :: SatVap,T
  SatVap=p0Star*(T/Tes0)**((Cpv-Cpl)/Rv) &
        *EXP((L00/Rv)*(One/Tes0-One/T))
END FUNCTION  SatVap


SUBROUTINE Res(z,y,yprime,delta,ires,rpar,ipar)

  INTEGER :: ires,ipar(*)
  REAL(RealKind) :: z,y(*),yprime(*),delta(*),rpar(*)

  INTEGER :: Case

  Case=ipar(1)
  SELECT CASE(Case)
    CASE(Dry)
      CALL ResDry(z,y,yprime,delta,ires,rpar,ipar)
    CASE(Moist)
      WRITE(*,*) 'ResMoist'
      CALL ResMoist(z,y,yprime,delta,ires,rpar,ipar)
  END SELECT  

END SUBROUTINE Res

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
  pvs=SatVap(Te) 
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
  pvs=SatVap(Te)
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
 

  Delta(1)=pPrime+(RhoD+RhoV)*Grav        ! dp/dz = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)         ! p = rho*R*T
  Delta(3)=Theta0-Te*(p0/p)**(Rd/Cpd)     ! Potential temperature
  Delta(4)=RhoV-RelHum*SatVap(Te)/(Rv*Te) ! RhoV

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
 
  Delta(1)=pPrime+Rho*Grav                 ! dp/dz = - rho*g
  Delta(2)=p-Te*Rd*Rho                     ! p = rho*R*T
  Delta(3)=Theta-Te*(p0/p)**(Rd/Cpd)        ! Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)

END SUBROUTINE ResDry

SUBROUTINE ResMoist(z,y,yprime,delta,ires,rpar,ipar)

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
  WRITE(*,*) 'y  ',y(1:6)
  WRITE(*,*) 'yPr',yPrime(1:6)

  zu=Rpar(1)
  zo=Rpar(2)
  Tu=Rpar(3)
  To=Rpar(4)
  RelHu=Rpar(5)
  RelHo=Rpar(6)
 
  RhoD=Rho-RhoV
  Delta(1)=pPrime+Rho*Grav                 ! dp/dz = - rho*g
  Delta(2)=p-Te*(Rd*RhoD+Rv*RhoV)          ! p = rho*R*T
  Rm=Rd*RhoD+Rv*RhoV
  Cpml=Cpd*RhoD+Cpv*RhoV
  Delta(3)=Theta-Te*(p0/p)**(Rm/Cpml)*(RhoD+RhoV*Rd/Rv)/Rho        ! Density Potential temperature
  Delta(4)=Te-((z-zu)*To+(zo-z)*Tu)/(zo-zu)
  Delta(5)=RhoV-RelH*SatVap(Te)/(Rv*Te) ! RhoV
  Delta(6)=RelH-((z-zu)*RelHo+(zo-z)*RelHu)/(zo-zu)

END SUBROUTINE ResMoist

SUBROUTINE InputRes(y,yprime,rpar,ipar,prestart,rhstart,teStart)

  REAL(RealKind) :: y(:),yprime(:),rpar(:)
  INTEGER :: ipar(:)
  REAL(RealKind), OPTIONAL :: preStart
  REAL(RealKind), OPTIONAL :: rhStart
  REAL(RealKind), OPTIONAL :: teStart

  REAL(RealKind) :: Theta0,p,Te,e
  CHARACTER*20:: par
  INTEGER :: Case

  Case=ipar(1)
  SELECT CASE(Case)
    CASE(Dry)
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
    CASE(Moist)
      neq=6
      PrePos=1
      TempPos=2
      ThPos=3
      RhoPos=4
      RhoVPos=5
      IF (PRESENT(preStart)) THEN
        y(1)=preStart
      ELSE
        IF (y(1)==0.0d0) THEN
          y(1)=1.0d5 
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
      IF (PRESENT(rhStart)) THEN
        y(6)=rhStart
        rpar(5)=rhStart
        rpar(6)=rhStart
      ELSE
        y(6)=0.7d0
        rpar(5)=0.7d0
        rpar(6)=0.7d0
      END IF  
      IF (y(3)==0.0d0) THEN
        y(3)=y(2)*(p0/y(1))**(Rd/Cpd)
      END IF  
      IF (y(4)==0.0d0) THEN
        y(4)=y(1)/(Rd*y(2))
      END IF  
      y(5)=y(6)*SatVap(y(2))/(Rv*y(2))
      yPrime(1)=-y(RhoPos)*Grav
      rpar(1)=0.0d0
      rpar(2)=1.0d0
    CASE DEFAULT  
      y=0.0d0
      yPrime=0.0d0
  END SELECT  

END SUBROUTINE InputRes

SUBROUTINE OutputRes(nz,zM,c,FileNameOut,Case)


  INTEGER :: nz
  REAL(RealKind) ::zM(:),c(:,:)
  CHARACTER(*) :: FileNameOut
  INTEGER :: Case

  INTEGER :: i

  OPEN(unit=OutputUnit,FILE=TRIM(FileNameOut)//'.prof')
  SELECT CASE(Case)
    CASE(Dry)
      WRITE(OutputUnit,*) 'RhoProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'ThDensProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,ThPos)
      END DO
      WRITE(OutputUnit,*) 'TempProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,TempPos)
      END DO
      WRITE(OutputUnit,*) 'PreProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,PrePos)
      END DO
    CASE(Moist)
      WRITE(OutputUnit,*) 'RhoProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,RhoPos)
      END DO
      WRITE(OutputUnit,*) 'ThDensProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,ThPos)
      END DO
      WRITE(OutputUnit,*) 'TempProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,TempPos)
      END DO
      WRITE(OutputUnit,*) 'RhoVProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,RhoVPos)
      END DO
      WRITE(OutputUnit,*) 'PreProf'
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
      DO i=1,nz
        WRITE(OutputUnit,*) zM(i),c(i,PrePos)
      END DO
  END SELECT    
  CLOSE(OutputUnit)

END SUBROUTINE OutputRes


SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)

  INTEGER :: ipar(*)
  REAL(RealKind) :: t,y(*),yprime(*),PD(*),CJ,rpar(*)

END SUBROUTINE JAC


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
  IF (InvLayerHeight.GT.0.d0) dz_Inv=MAX(dz_Inv,dz_ManProf)

  Pre(1) = Pres0
  IF (PotTemp0.GT.0.d0) THEN ! Profile calculated with given potential temperature
    IF (PotTemp0.LT.1.d2) PotTemp0=PotTemp0+273.15
    PotTemp(1) = PotTemp0
    Temp(1)=(Pres0/1000.d2)**(Rd/Cpd)*PotTemp(1)
    DO iz=2,NumSounding
      IF (Height(iz).LE.InvLayerHeight) THEN
        dTemp=dTemp_BL
      ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dz_Inv) THEN
        dTemp=dTemp_Inv
      ELSE IF (Height(iz).GT.InvLayerHeight+dz_Inv) THEN
        dTemp=dTemp_FT
      END IF
      PotTemp(iz)= PotTemp(iz-1)+dTemp*dz_ManProf
      Temp(iz)=Temp(iz-1)
      DO i=1,20 ! Itertation, because of two unknown variables
        Pre(iz) = Pre(iz-1)*EXP((-1.d0*grav*dz_ManProf)/(Rd*(Temp(iz-1)+Temp(iz))/2.d0))
        Temp(iz) = PotTemp(iz)*(Pre(iz)/1000.d2)**(Rd/Cpd)
      END DO
    END DO
  ELSE ! Profile calculated with given absolute temperature
    IF (Temp0.LT.1.d2) Temp0=Temp0+273.15
    Temp(1)=Temp0
    DO iz=2,NumSounding
      IF (Height(iz).LE.InvLayerHeight) THEN
        dTemp=dTemp_BL
      ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dz_Inv) THEN
        dTemp=dTemp_Inv
      ELSE IF (Height(iz).GT.InvLayerHeight+dz_Inv) THEN
        dTemp=dTemp_FT
      END IF
      Temp(iz)= Temp(iz-1)+dTemp*dz_ManProf
      Pre(iz) = Pre(iz-1)*EXP((-1.d0*grav*dz_ManProf)/(Rd*(Temp(iz-1)+Temp(iz))/2.d0))
    END DO
  END IF

  IF (RelHum0.GT.0.d0 .OR. QV0.GT.0.d0) THEN ! Moist case
    CaseProfile=Moist
    RhoL = Pre(1)/(Rd*Temp(1))
    IF (RelHum0.GT.0.d0) THEN
      RH(1)=RelHum0
      VaporPres = RH(1)/1.d2*SatVap(Temp(1))
      RhoV = VaporPres/(Rv*Temp(1))
      QV(1)=RhoV/(RhoV+RhoL)
    ELSE IF (QV0.GT.0.d0) THEN
      QV(1)=QV0
      RhoV = RhoL*(QV(1)/(1.d0-QV(1)))
      VaporPres = RhoV*Rv*Temp(1)
      RH(1) = VaporPres/SatVap(Temp(1))*1.d2
    END IF


    IF (dRH_FT.GE.0.d0) THEN
      DO iz=2,NumSounding
        IF (Height(iz).LE.InvLayerHeight) THEN
          dRelHum=dRH_BL
        ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dz_Inv) THEN
          dRelHum=dRH_Inv
        ELSE IF (Height(iz).GT.InvLayerHeight+dz_Inv) THEN
          dRelHum=dRH_FT
        END IF
        RH(iz)= RH(iz-1)+dRelHum*dz_ManProf
      END DO
    ELSE IF (dQV_FT.GE.0.d0) THEN
      DO iz=2,NumSounding
        RhoL = Pre(iz)/(Rd*Temp(iz))
        IF (Height(iz).LE.InvLayerHeight) THEN
          dQV=dQV_BL
        ELSE IF (Height(iz).GT.InvLayerHeight .AND. Height(iz-1).LE.InvLayerHeight+dz_Inv) THEN
          dQV=dQV_Inv
        ELSE IF (Height(iz).GT.InvLayerHeight+dz_Inv) THEN
          dQV=dQV_FT
        END IF
        QV(iz)= QV(iz-1)+dQV*dz_ManProf
        RhoV = RhoL*(QV(iz)/(1.d0-QV(iz)))
        VaporPres = RhoV*Rv*Temp(iz)
        RH(iz) = VaporPres/SatVap(Temp(iz))*1.d2
      END DO
    ELSE
      STOP 'ERROR: No gradient for humidity defined for free troposphere'
    END IF

    RH=RH/1.d2
    WRITE(*,*) '===== Wet profile calculated'
  ELSE ! Dry case
    WRITE(*,*) '===== Dry profile calculated'
    CaseProfile=Dry
  END IF

END SUBROUTINE CreateProfile

END MODULE Profile_Mod

