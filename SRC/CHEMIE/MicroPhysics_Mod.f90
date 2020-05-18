MODULE MicroPhysics_mod

  USE Kind_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Chemie_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Aerosol_Mod

  IMPLICIT NONE
  LOGICAL :: Flag
  REAL(RealKind) :: RhoFStart
  REAL(RealKind), PARAMETER :: DeltaV=1.096d-7 ! Freie Wegl�nge
  REAL(RealKind), PARAMETER :: DeltaT=2.16d-7  ! Freie Wegl�nge
  REAL(RealKind), PARAMETER :: alphaC=1.0d0    ! Accomodation coefficient for 
! REAL(RealKind), PARAMETER :: alphaC=0.042   ! Accomodation coefficient for 
                                        ! the condensation of water vapour
  REAL(RealKind), PARAMETER :: alphaT=0.70     ! Thermal accomodation coefficient  
  REAL(RealKind), PARAMETER :: Press=0.9d5     ! Total pressure
  REAL(RealKind), PARAMETER :: RUniv=0.08206   ! allgemeine Gaskonstante [L atm/(mol K)]
! REAL(RealKind), PARAMETER :: Rv=461.5d0      ! Gaskonstante f�r Wasserdampf
! REAL(RealKind), PARAMETER :: Rd=287.04d0     ! Gaskonstante f�r trockene Luft
! REAL(RealKind), PARAMETER :: Cpd=1.005e3
  REAL(RealKind), PARAMETER :: kth=2.43d-2     ! Thermische Leitf�higkeit
  REAL(RealKind), PARAMETER :: RhoW=1.0d3      ! Dichte von Wasser
  REAL(RealKind), PARAMETER :: RhoWL=1.0d0      ! Dichte von Wasser
  REAL(RealKind), PARAMETER :: RhoSalt=1500.d0    ! Dichte des Salzes
! REAL(RealKind), PARAMETER :: RhoSalt=1800.d0    ! Dichte des Salzes
  REAL(RealKind), PARAMETER :: MolW=18.0d0     ! Molare Masse von Wasser 
  REAL(RealKind), PARAMETER :: MolS=132.13d0   ! Molare Masse des Salzes
                                        ! Hier Ammoniumsulfat
  REAL(RealKind) :: MassSalt
  REAL(RealKind), PARAMETER :: VantS=3.0d0     ! van't Hoff Zahl des Salzes 
  INTEGER :: iCell,iF,iEul
  
  REAL(RealKind) :: Rad1,quot,rate


  INTERFACE CondRateJac
    MODULE PROCEDURE CondRateWaterJac
  END INTERFACE  CondRateJac


CONTAINS 
FUNCTION Kohler(Satt,Mass,Rad,TAbs,Act)
  REAL(RealKind) :: Kohler
  REAL(RealKind) :: Satt,Mass(:),Rad,TAbs,Act
  REAL(RealKind) :: Kel,Rao

  Kel=KelvinFak(tAbs,Mass)
  Rao=RaoultFak(Mass)
  Kohler=Satt-Act*Rao*EXP(Kel/Rad)
END FUNCTION Kohler

FUNCTION KohlerP(Mass,Rad,TAbs,Act)
  REAL(RealKind) :: KohlerP
  REAL(RealKind) :: Mass(:),Rad,TAbs,Act
  REAL(RealKind) :: Kel,Rao
  REAL(RealKind) :: KelP,RaoP,RadP

  Kel=KelvinFak(tAbs,Mass)
  Rao=RaoultFak(Mass)
  KelP=KelvinFakP(tAbs,Mass)
  RaoP=RaoultFakP(Mass)
  RadP=RadiusP(Mass)
  KohlerP=-Act*EXP(Kel/Rad) &
          *(RaoP+Rao*(KelP*Rad-Kel*RadP)/(Rad*Rad))

END FUNCTION KohlerP

FUNCTION Kohler1(Satt,Mass,Rad,TAbs)
  REAL(RealKind) :: Kohler1
  REAL(RealKind) :: Satt,Mass(:),Rad,TAbs
  REAL(RealKind) :: Kel,Rao

  Kel=KelvinFak(tAbs,Mass)
  Rao=RaoultFak(Mass)
  Kohler1=Satt-Rao*EXP(Kel/Rad)
END FUNCTION Kohler1

FUNCTION Kohler1P(Mass,Rad,TAbs)
  REAL(RealKind) :: Kohler1P
  REAL(RealKind) :: Mass(:),Rad,TAbs
  REAL(RealKind) :: Kel,Rao
  REAL(RealKind) :: KelP,RaoP,RadP

  Kel=KelvinFak(tAbs,Mass)
  Rao=RaoultFak(Mass)
  KelP=KelvinFakP(tAbs,Mass)
  RaoP=RaoultFakP(Mass)
  RadP=RadiusP(Mass)
  Kohler1P=-EXP(Kel/Rad) &
          *(RaoP+Rao*(KelP*Rad-Kel*RadP)/(Rad*Rad))

END FUNCTION Kohler1P


FUNCTION RaoultFak(Mass) 

  REAL(RealKind) :: RaoultFak
  REAL(RealKind) :: Mass(:)
  REAL(RealKind) :: MolGes,MolWater
  INTEGER :: i

  MolWater=Mass(iWater)/MolMass(iWater)
  MolGes=MolWater
  DO i=iWater+1,nSoluble
    MolGes=MolGes &
           +ABS(Mass(i)/MolMass(i))
  END DO
  RaoultFak=MolWater/MolGes

END FUNCTION RaoultFak

FUNCTION RaoultFakP(Mass)

  REAL(RealKind) :: RaoultFakP
  REAL(RealKind) :: Mass(:)
  REAL(RealKind) :: MolSalz,MolWater
  INTEGER :: i

  MolSalz=Zero
  DO i=iWater+1,nSoluble
    MolSalz=MolSalz &
           +Mass(i)/MolMass(i)
  END DO
  MolWater=Mass(iWater)/MolMass(iWater)

  RaoultFakP=MolSalz/((MolWater+MolSalz)*(MolWater+MolSalz)*MolMass(iWater))

END FUNCTION RaoultFakP

FUNCTION KelvinFak(tAbs,Mass) ! unter Bedingung, dass Sigma passt, entspricht dieser Term dem Kelvinterm i.O.

  REAL(RealKind) :: KelvinFak,Mass(:)
  REAL(RealKind) :: tAbs

  REAL(RealKind) :: tau,Sigma
  
  tau=1.0d0-tAbs/647.069d0
  Sigma=235.8d-3*tau**1.256d0*(1.0d0-0.625d0*tau) 
  KelvinFak=2.d0*Sigma/(tAbs*Rw*RhoW)

END FUNCTION KelvinFak

FUNCTION KelvinFakP(tAbs,Mass)

  REAL(RealKind) :: KelvinFakP,Mass(:)
  REAL(RealKind) :: tAbs


  KelvinFakP=Zero

END FUNCTION KelvinFakP

FUNCTION TotalMass(Mass)

  REAL(RealKind) :: TotalMass
  REAL(RealKind) :: Mass(:)

  INTEGER :: i

  TotalMass=Zero
  IF (iWater>0.AND..NOT.Dry) THEN
    TotalMass=Mass(iWater)
  END IF
  DO i=1+nSoluble,nAqua
    TotalMass=TotalMass+Mass(i)
  END DO
  DO i=MAX(iWater,iNC)+1,nSoluble
    TotalMass=TotalMass+Mass(i)
  END DO

END FUNCTION TotalMass

FUNCTION Radius(Mass)

  REAL(RealKind) :: Radius
  REAL(RealKind) :: Mass(:)

  INTEGER :: i,n
  REAL(RealKind) :: Vol

  n=SIZE(Mass)
  Vol=Zero
  IF (iWater>0) THEN
    Vol=Mass(iWater)/RhoSpecies(iWater)
  END IF
  DO i=nSoluble+1,nAqua
    Vol=Vol+Mass(i)/RhoSpecies(i)
  END DO
  IF (Vol<Zero) THEN
    WRITE(*,*) 'Mass(iWater)',Mass(iWater)
    WRITE(*,*) 'Bin: ', iCell, 'z-Cell:', iF
    WRITE(*,*) 'nSoluble',nSoluble,nAqua
    DO i=nSoluble+1,nAqua
      WRITE(*,*) 'Mass(i)',Mass(i)
    END DO
    STOP 'Vol kleiner Null'
  END IF
  Radius=(3.0d0*Vol/(4.0d0*Pi))**(1.0d0/3.0d0)

END FUNCTION Radius

FUNCTION RadiusP(Mass)

  REAL(RealKind) :: RadiusP
  REAL(RealKind) :: Mass(:)

  INTEGER :: i,n
  REAL(RealKind) :: Vol

  n=SIZE(Mass)
  IF (iWater>0) THEN
    Vol=Mass(iWater)/RhoSpecies(iWater)
    DO i=nSoluble+1,nAqua
      Vol=Vol+Mass(i)/RhoSpecies(i)
    END DO
    RadiusP=(3.0d0*Vol/(4.0d0*Pi))**(-2.0d0/3.0d0)*3.0d0/(4.0d0*Pi)/RhoSpecies(iWater)/3.0d0
  ELSE
    RadiusP=Zero
  END IF

END FUNCTION RadiusP

FUNCTION RadiusWater(Mass)

  REAL(RealKind) :: RadiusWater
  REAL(RealKind) :: Mass(:)

  REAL(RealKind) :: Vol

  IF (iWater>0) THEN
    Vol=Mass(iWater)/RhoSpecies(iWater)
    RadiusWater=(3.0d0*Vol/(4.0d0*Pi))**(1.0d0/3.0d0)
  ELSE
    RadiusWater=Zero
  END IF

END FUNCTION RadiusWater

FUNCTION CondRateHenry(Mass,Gas,Reaction,Act,RhoLuft)

  REAL(RealKind) :: CondRateHenry
  REAL(RealKind) :: Mass(:) ! kg
  REAL(RealKind) :: Gas(:)
  TYPE (Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: Act,RhoLuft

  INTEGER :: i,is,iGas
  REAL(RealKind) :: Rad,Temp,kt,va 
  TYPE(Vec4_T) :: w
  TYPE(Vec4_T) :: TVec
  REAL(RealKind), POINTER :: Constants(:)

  CALL Allocate(w,1)
  CALL Allocate(TVec,1)
  Temp=Gas(thPos-nAqua)
  TVec%c=Temp
  Constants=>Reaction%Constants(1,1,1,:)
  SELECT CASE (Reaction%TypeR)
    CASE ('CONST')
      CALL ConstCompute(w,Constants)
    CASE ('TEMP3')
      CALL Temp3Compute(w,Constants,TVec)
  END SELECT
! Computation of the Schwartz-Coefficient
  Rad=Radius(Mass)
  iGas=Reaction%SpeciesLeft(1)-nAqua
  is=Reaction%SpeciesRight(1)
! Thermal Velocity
  va=SQRT((8.0d0*BolzMol*Temp)/(Pi*MolMass(is))) ! m/s
  kt=One/(Rad*Rad/(Three*Dg(iGas)) & ! 1/s
         +4.0d0*Rad/(Three*va*Accomod(iGas)))
  CondRateHenry=kt*Mass(iWater)*(Gas(iGas)/RhoLuft/(RhoSpecies(iWater)) &
                   -(Act*Mass(is)/MolMass(is)/Mass(iWater)) & 
                    /(w%c(1,1,1,1)*Temp*RUniv))
  CALL Deallocate(w)
  CALL Deallocate(TVec)
  
END FUNCTION CondRateHenry


FUNCTION CondRateHenryJac(Mass,Gas,Reaction,Act,RhoLuft)
 
  REAL(RealKind) :: CondRateHenryJac(3)
  REAL(RealKind) :: Mass(:)
  REAL(RealKind) :: Gas(:)
  TYPE (Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: Act,RhoLuft
 
  INTEGER :: i,is,iGas
  REAL(RealKind) :: Rad,Temp,kt,va
  REAL(RealKind) :: dRaddWater,A,B,dktdWater
  TYPE(Vec4_T) :: w
  TYPE(Vec4_T) :: TVec
  REAL(RealKind), POINTER :: Constants(:)
 
  CALL Allocate(w,1)
  CALL Allocate(TVec,1)
  Temp=Gas(thPos-nAqua)
  TVec%c=Temp
  Constants=>Reaction%Constants(1,1,1,:)
  SELECT CASE (Reaction%TypeR)
    CASE ('CONST')
      CALL ConstCompute(w,Constants)
    CASE ('TEMP3')
      CALL Temp3Compute(w,Constants,TVec)
  END SELECT
! Computation of the Schwartz-Coefficient
  Rad=Radius(Mass)
  dRaddWater=RadiusP(Mass)
  iGas=Reaction%SpeciesLeft(1)-nAqua
  is=Reaction%SpeciesRight(1)
! Thermal Velocity
  va=SQRT((8.0d0*BolzMol*Temp)/(Pi*MolMass(is)))
  A=One/(Three*Dg(iGas))
  B=4.0d0/(Three*va*Accomod(iGas))
  kt=One/(A*Rad*Rad+B*Rad)
  dktdWater=-kt*kt*(2.0d0*A*Rad+B)*dRaddWater
! CondRateHenry=kt*(Gas(iGas)*Mass(iWater)/(RhoSpecies(iWater)) &
!                  -(Mass(is)/MolMass(is)) &
!                   /(w%c(1)*Temp*RUniv))
  CondRateHenryJac(1)=kt*Mass(iWater)/(RhoSpecies(iWater))
  CondRateHenryJac(2)=-kt*Act/MolMass(is) &
                    /(w%c(1,1,1,1)*Temp*RUniv)
  CondRateHenryJac(3)=kt*Gas(iGas)/RhoLuft/(RhoSpecies(iWater)) &
                     +dktdWater                &
                     *(Gas(iGas)/RhoLuft*Mass(iWater)/(RhoSpecies(iWater)) &
                   -(Act*Mass(is)/MolMass(is)) &
                    /(w%c(1,1,1,1)*Temp*RUniv))
  CALL Deallocate(w)
  CALL Deallocate(TVec)
 
END FUNCTION CondRateHenryJac

SUBROUTINE Adjust(dqv,qv,p,tAbs)

  REAL(RealKind) :: dqv,qv,p,tAbs

  REAL(RealKind) :: dqvOld,es,qvStar,qvNeu,tAbsNeu,RhoF
  INTEGER :: Iter

  dqv=Zero
  dqvOld=dqv
  Iter=0
  DO 
    Iter=Iter+1
    qvNeu=qv+dqv
    tAbsNeu=tAbs-lv/Cpd*dqv
    RhoF=RhoFF(tAbsNeu,p,qvNeu)
    es=611.2d0*EXP(17.67d0*(tAbsNeu-273.15d0)/(tAbsNeu-29.65d0))
    qvStar=es/(rhoF*rv*tAbsNeu)
    dqvOld=dqv
    dqv=qvStar-qv 
    IF (ABS(dqv-dqvOld)<=qvStar*1.d-5.OR.Iter>1000) THEN
      EXIT
    END IF
  END DO

END SUBROUTINE Adjust


FUNCTION CondRateWater(Mass,Gas,Act)

  REAL(RealKind) :: CondRateWater
  REAL(RealKind) :: Mass(:)
  REAL(RealKind) :: Gas(:)
  REAL(RealKind) :: Act

  REAL(RealKind) :: RhoV,Rho,RhoD,RhoL,TAbs,p,pV,pVs,PotM
  REAL(RealKind) :: Diff,kthDyn,dmFak,Satt,Rad
  REAL(RealKind) :: KoLoc
  INTEGER :: iLoc

  RhoV=Gas(RhoVPos-nAqua)
  TAbs=Gas(thPos-nAqua)
  Rho=Gas(rhoPos-nAqua)
! p=Gas(prePos-nAqua)
! RhoF=RhoFF(tAbs,p,qv)
  RhoD=Rho-RhoV  !  -RhoL+Eps  OSSI
  RhoL=0.0d0     !             OSSI  
  pVs=SaturVapor(TAbs)
  pV=Rv*TAbs*RhoV
  Satt=pV/pVs
  p=pv+Rd*TAbs*RhoD
  
! Berechnung der Kondensation/Deposition
! ohne Zusatzeffekte (Kruemmung, Loesung, Ventilation, Kinetisch)
! Hilfsgroessen/zusammengesetzte Ausdruecke
  Rad=Radius(Mass)
! Diff=8.79d-5*(tAbs**1.81D0)/p
  Diff=4.0122d-5*(tAbs**1.94d0)/p
  kthDyn=(5.69d0+0.017d0*(tAbs-273.15d0))*418.5d-5
  Diff=Diff/(Rad/(Rad+DeltaV)+Diff/Rad/alphaC*SQRT(2.0d0*Pi/tAbs/rv))
  kthDyn=kthDyn/(Rad/(Rad+DeltaT) &
        +kthDyn/Rad/alphaT/Cpd/RhoD*SQRT(2.0d0*Pi/tAbs/rd))
  dmFak=4.d0*Pi/(rv*tAbs/(Diff*pVs)+(lv/(rv*tAbs)-1.D0)*lv/(kthDyn*tAbs))

  KoLoc=Kohler(Satt,Mass,Rad,TAbs,Act)
  
  CondRateWater=dmFak*Rad*KoLoc

END FUNCTION CondRateWater

FUNCTION CondRateWaterJac(Mass,Gas,Act)

  REAL(RealKind) :: Mass(:)
  REAL(RealKind) :: Gas(:)
  REAL(RealKind) :: CondRateWaterJac(SIZE(Mass)+SIZE(Gas))
  REAL(RealKind) :: Act

  REAL(RealKind) :: f,f0,Temp
  INTEGER :: iMass,LenMass
  INTEGER :: iGas,LenGas
  REAL(RealKind), PARAMETER :: h=1.d-7

  LenMass=SIZE(Mass)
  LenGas=SIZE(Gas)
  f0=CondRateWater(Mass,Gas,Act)
  CondRateWaterJac=0.0d0
  DO iMass=1,LenMass
    Temp=Mass(iMass)
    Mass(iMass)=Mass(iMass)*(1.0d0+h)+1.d-50
    f=CondRateWater(Mass,Gas,Act)
    CondRateWaterJac(iMass)=(f-f0)/(Mass(iMass)-Temp)
    Mass(iMass)=Temp
  END DO

  DO iGas=1,LenGas
    Temp=Gas(iGas)
    Gas(iGas)=Gas(iGas)*(1.0d0+h)+1.d-50
    f=CondRateWater(Mass,Gas,Act)
    CondRateWaterJac(iGas+SIZE(Mass))=(f-f0)/(Gas(iGas)-Temp)
    Gas(iGas)=Temp
  END DO

END FUNCTION CondRateWaterJac

FUNCTION RhoAF(tAbs,p)

   REAL(RealKind) :: RhoAF
   
   REAL(RealKind) :: tAbs,p

   RhoAF=p/(Rd*tAbs)

END FUNCTION RhoAF

FUNCTION RhoFF(tAbs,p,qv)

   REAL(RealKind) :: RhoFF

   REAL(RealKind) :: tAbs,p,qv

   RhoFF=p/(Rd*tAbs*(1.0d0+(Rv-Rd)/Rd*qv))

END FUNCTION RhoFF

FUNCTION SattF(qv,tAbs,p)

  REAL(RealKind) :: SattF

  REAL(RealKind) :: qv,tAbs,p

  REAL(RealKind) :: RhoF,ee,es

  RhoF=RhoFF(tAbs,p,qv)
  es=qvs(RhoF,tAbs,tAbs)
  ee=qv
  SattF=ee/es
END FUNCTION SattF


END MODULE MicroPhysics_mod
