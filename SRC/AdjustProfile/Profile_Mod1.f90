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
  LOGICAL :: ProfIn=.TRUE.
  LOGICAL :: HorizShear=.FALSE.
  REAL(RealKind) :: Tropopause  = 3000.0d0
  REAL(RealKind) :: dd,ff

  NAMELIST /AdjustCtrl/ NameSounding  &
                       ,FileNameOut   &
                       ,Tropopause    &
                       ,HorizShear    &
                       ,CaseProfile

  INTEGER,PARAMETER :: NumSounding=10000
  REAL(RealKind) :: Height(NumSounding)
  REAL(RealKind) :: Pre(NumSounding)
  REAL(RealKind) :: Temp(NumSounding)
  REAL(RealKind) :: RH(NumSounding)
  REAL(RealKind) :: U(NumSounding)
  REAL(RealKind) :: V(NumSounding)
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
  IF (ProfIn) THEN
    OPEN(InputUnit,FILE=TRIM(NameSounding)//'.txt',ACTION='read',IOSTAT=stat)
    IF (stat==0) THEN
      DO k=1,19
        READ(InputUnit,'(A)') DummyChar 
      END DO
      DO iz=1,NumSounding
        READ(InputUnit,*) Height(iz),Pre(iz),Temp(iz),RH(iz),dd,ff
        Pre(iz)=Pre(iz)*100.0d0 ! Umrechnung hPa Pa
        Temp(iz)=Temp(iz)+273.15d0
        RH(iz)=RH(iz)/100.0d0
        IF (HorizShear) THEN
          U(iz)=-SIN(dd/180.0d0*Pi)*ff
          V(iz)=-COS(dd/180.0d0*Pi)*ff
        ELSE
          U(iz)=ff
          V(iz)=0.0d0
        END IF
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
! WRITE(*,*) 'y       ',y(1:neq)
! WRITE(*,*) 'yPrime  ',yPrime(1:neq)
! WRITE(*,*) 'Delta   ',Delta(1:neq)

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
      WRITE(*,*) 'Init Case Dry'
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
  END SELECT  

END SUBROUTINE InputRes

SUBROUTINE OutputResWind(nz,zM,U,V,FileNameOut)

  INTEGER :: nz
  REAL(RealKind) :: U(:),V(:),zM(:)
  CHARACTER(*) :: FileNameOut
 
  INTEGER :: i

  OPEN(unit=OutputUnit,FILE=TRIM(FileNameOut)//'.prof',POSITION='APPEND')
  WRITE(OutputUnit,*) 'UProf'
  WRITE(OutputUnit,*) nz+1,2
  WRITE(OutputUnit,*) 'Equal'
  WRITE(OutputUnit,*) 0.0d0,0.0d0
  DO i=1,nz
    WRITE(OutputUnit,*) zM(i),U(i)
  END DO
  WRITE(OutputUnit,*) 'VProf'
  WRITE(OutputUnit,*) nz+1,2
  WRITE(OutputUnit,*) 'Equal'
  WRITE(OutputUnit,*) 0.0d0,0.0d0
  DO i=1,nz
    WRITE(OutputUnit,*) zM(i),V(i)
  END DO
  CLOSE(OutputUnit)
END SUBROUTINE OutputResWind

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
END MODULE Profile_Mod

