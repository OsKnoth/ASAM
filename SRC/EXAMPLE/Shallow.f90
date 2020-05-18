MODULE Shallow_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  LOGICAL :: ProfIn=.FALSE.
 
  REAL(RealKind), PARAMETER :: GravH02=2.94d4 ! m^2/s^2
  REAL(RealKind), PARAMETER :: GravH03=2.94d4 ! m^2/s^2
  REAL(RealKind), PARAMETER :: GravH04=1.0d5 ! m^2/s^2
  REAL(RealKind), PARAMETER :: H05=5960.0d0   ! m
  REAL(RealKind), PARAMETER :: H06=8000.0d0   ! m
  REAL(RealKind), PARAMETER :: H0G=10000.0d0   ! m
  REAL(RealKind), PARAMETER :: GravH0MB=5.6768d04   ! m^2/s^2
  REAL(RealKind), PARAMETER :: VMB=20                ! m/s
  REAL(RealKind), PARAMETER :: sigma4=12.74244d0*12.74244d0
  REAL(RealKind), PARAMETER :: Omega4=7.2920d-5
  REAL(RealKind), PARAMETER :: omega6=7.8480d-6  ! s^-1 
  REAL(RealKind), PARAMETER :: K6=7.8480d-6   ! s^-1
  REAL(RealKind), PARAMETER :: R6=4.0d0   
  REAL(RealKind), PARAMETER :: xe3=0.3d0
  REAL(RealKind) :: phib3
  REAL(RealKind) :: phie3
  REAL(RealKind) :: U03
  REAL(RealKind) :: U04=20.0d0
  REAL(RealKind) :: lam0,phi0
  REAL(RealKind) :: RadiusC
  REAL(RealKind) :: alpha1=0.0d0
  REAL(Realkind) :: alpha2=0.0d0
  REAL(Realkind) :: alpha3=0.0d0
  REAL(Realkind) :: alphaG=1.0d0/3.0d0
  REAL(Realkind) :: betaG=1.0d0/15.0d0
  REAL(Realkind) :: hH=120.0d0
  REAL(RealKind) :: lam1=270.0d0
  REAL(RealKind) :: phi1=0.0d0
  REAL(RealKind) :: U1=38.610d0
  REAL(RealKind) :: R1=10000.0d0
  REAL(RealKind) :: H01=1000.0d0
  REAL(RealKind) :: AA=10.0d0
  REAL(RealKind) :: uMax
  REAL(RealKind) :: hS=2000.0d0   ! m
  CHARACTER*40 :: Problem
  REAL(RealKind) :: rotation_angle
  !Jet
  REAL(8) :: aJet=6.0d6   ! m
  REAL(8) :: yCJet=6.0d6
  REAL(8) :: WidthJet=1.50d6
  REAL(8) :: betaJet=2.0d-11 ! 1/(ms) 
  REAL(8) :: h0Jet=3000.0d0  ! m
  REAL(8) :: GravJet=10.0d0  ! m/s^2
  REAL(8) :: u0Jet=80.0d0    ! m/s
  REAL(8) :: wXJet=1.5d6
  REAL(8) :: wYJet=2.0d6
  REAL(8) :: xCPJet=24.0d6
  REAL(8) :: yCPJet=6.0d6
  REAL(8) :: hHatJet=120.0d0
  
  NAMELIST /Example/ uMax &
                    ,hS   &
                    ,hHatJet &
                    ,u0Jet &
                    ,alpha1 &
                    ,U1 &
                    ,Problem &
                    ,alpha2 &
                    ,alpha3 

END MODULE Shallow_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Shallow_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Shallow_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Shallow_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line


  R1=RadEarth/Three
  U1=Two*Pi*RadEarth/(12.0d0*24.0d0*3600.0d0)
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
  lam0=1.5d0*Pi
  phi0=Pi/6.0d0
  RadiusC=Pi/9.0d0
  alpha1=Two*Pi*alpha1/360.0d0
  alpha3=Two*Pi*alpha3/360.0d0
  lam1=Two*Pi*lam1/360.0d0
  phi1=Two*Pi*phi1/360.0d0
  phib3=-1.0d0/6.0d0*Pi
  phie3=Half*Pi
  U03=Two*Pi*RadEarth/(12.0d0*24.0d0*3600.0d0)
  rotation_angle = RotAngle
  SELECT CASE (Problem)
    CASE ('Jet')
      GravJet=Grav
  END SELECT  
END SUBROUTINE InputExample

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Shallow_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time

  RhoProf=0.0d0

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Shallow_Mod
  USE Start_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time


  ThProfFun=0.0d0

END FUNCTION ThProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Shallow_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  INTERFACE
    FUNCTION Simpson(x0,x1,dx,Integrand)
      USE Kind_Mod
      REAL(RealKind) :: Simpson
      REAL(RealKind) :: x0,x1,dx
      INTERFACE
        FUNCTION Integrand(tau)
          USE Kind_Mod
          REAL(RealKind) :: Integrand,tau
        END FUNCTION Integrand
      END INTERFACE
    END FUNCTION Simpson
  END INTERFACE
  INTERFACE
    FUNCTION Integrand3(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand3,tau
    END FUNCTION Integrand3
  END INTERFACE
  INTERFACE
    FUNCTION Integrand4(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand4,tau
    END FUNCTION Integrand4
  END INTERFACE
  INTERFACE
    FUNCTION IntegrandG(tau)
      USE Kind_Mod
      REAL(RealKind) :: IntegrandG,tau
    END FUNCTION IntegrandG
  END INTERFACE

  REAL(RealKind) :: lam,phi,phiPr,r,HeightLoc
  REAL(RealKind) :: A,B,C,SurPot,tetaQ,f,f04,teta04,uQ,hQ,C4,phi04,phi2G
  REAL(RealKind) :: lam04,TimeLoc
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: rot_lam0,rot_phi0
  REAL(RealKind) :: yHat,Pert

  SELECT CASE (Problem)
    CASE ('Damm')
      IF (x<=Zero) THEN
        RhoFun=10.0d0
      ELSE
        RhoFun=9.5d0
      END IF
    CASE ('William1')
      RhoFun=One
    CASE ('William2')
      phi=y
      RhoFun=(GravH02-(RadEarth*Omega*UMax+Half*UMax*Umax)*SIN(phi)*SIN(phi))/Grav
    CASE ('William3')
      lam=x
      phi=y
      phiPr=ASIN(SIN(phi)*COS(alpha3)-COS(phi)*COS(lam)*SIN(alpha3))
      RhoFun=(GravH03-RadEarth*Simpson(-Half*pi,phiPr,Pi/100.0d0,Integrand3))/Grav
    CASE ('William4')
      lam=x
      phi=y
      phi04=Pi/4.0d0
      f=Two*Omega4*SIN(phi)
      f04=Two*Omega4*SIN(Pi/4)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*(SIN(2.0d0*phi))**14.0d0
      hQ=(GravH04-Simpson(-Half*pi,phi,Pi/100.0d0,Integrand4))/Grav
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1-C4)/(1+C4)))
      RhoFun=(Grav*hQ+f*tetaQ)/Grav
    CASE ('William5')
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      CALL Rotate(lam0,phi0,rot_lam0,rot_phi0)
      r=SQRT(MIN(RadiusC*RadiusC, &
                 (rot_lam-rot_lam0)*(rot_lam-rot_lam0)+(rot_phi-rot_phi0)*(rot_phi-rot_phi0)))
      HeightLoc=hS*(One-r/RadiusC)
      RhoFun=(Grav*H05-(RadEarth*Omega*UMax+Half*UMax*Umax)*SIN(rot_phi)*SIN(rot_phi))/Grav-HeightLoc
    CASE ('William6')
      lam=x
      phi=y
      A=Half*omega6*(Two*Omega+omega6)*COS(phi)*COS(phi) &
       +0.25d0*K6*K6*(COS(phi))**(Two*R6)*((R6+One)*COS(phi)*COS(phi) &
       +(Two*R6*R6-R6-Two)-Two*R6*R6*(COS(phi))**(-Two))
      B=Two*(Omega+omega6)*K6/(R6+One)/(R6+Two) &
       *(COS(phi))**R6*((R6*R6+Two*R6+Two) &
       -(R6+One)*(R6+One)*COS(phi)*COS(phi))
      C=0.25d0*K6*K6*(COS(phi))**(Two*R6)*((R6+One)*COS(phi)*COS(phi) &
       -(R6+Two))
      RhoFun=(Grav*H06+RadEarth*RadEarth &
            *(A+B*COS(R6*lam)+C*COS(Two*R6*lam)))/Grav
    CASE ('Jet')        
      yHat=(y-yCJet)/WidthJet
      IF (yHat<=-1.0d0) THEN
        RhoFun=h0Jet
      ELSE IF (yHat<=1.0d0) THEN
        RhoFun=h0Jet-WidthJet*betaJet*u0Jet/GravJet &
            *(yCJet*(16.0d0/35.0d0+yHat-yHat**3.0d0+3.0d0/5.0d0*yHat**5.0d0 &
                     -1.0d0/7.0d0*yHat**7.0d0) & 
             +WidthJet*(-1.0d0/8.0d0+0.5d0*yHat**2.0d0-3.0d0/4.0d0*yHat**4.0d0 &
                       +0.5d0*yHat**6.0d0-1.0d0/8.0d0*yHat**8.0d0)            &
             )
        Pert=hHatJet*COS(0.5d0*Pi*yHat)**2*EXP(-((x-xCPJet)/wXJet)**2)*EXP(-((y-yCPJet)/wYJet)**2)
        RhoFun=RhoFun+Pert
      ELSE
        RhoFun=h0Jet-32.0d0/35.0d0*WidthJet*betaJet*u0Jet*yCJet/GravJet
      END IF  
    CASE ('Galewsky')
      lam=x
      phi=y
      phi2G=Pi/4.0d0
      RhoFun=(Grav*H0G-(Simpson(-Half*Pi,phi,Pi/100.0d0,IntegrandG)))/Grav &
            +hH*COS(phi)*EXP(-((lam-Pi)/alphaG)**2.0d0)*EXP(-((phi2G-phi)/betaG)**2.0d0)
    CASE ('Tsunami')
      phi=y
      SurPot=-40000.d0+20000.d0*EXP(-1000.d0*(phi+Pi/6.d0)**2.d0)  
      RhoFun=(-SurPot+AA*EXP(-1000.d0*(phi-Pi/6.d0)**2.d0))/Grav
    CASE ('Dancing')
      lam=x
      phi=y
      RhoFun=(GravH0MB+Two*Omega*RadEarth*vMB*(SIN(phi))**Three*COS(phi)*SIN(lam))/Grav
    CASE ('Twin')
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      RhoFun=(Grav*H05+8.0d0*UMax*UMax*COS(rot_phi)*COS(rot_phi)*(1.0d0-4.0d0/3.0d0*COS(rot_phi) &
            +Half*COS(rot_phi)*COS(rot_phi))+4.0d0*RadEarth*Omega*UMax*COS(rot_phi)*COS(rot_phi) &
            *(1.0d0-2.0d0/3.0d0*COS(rot_phi)))/Grav
    CASE DEFAULT
      RhoFun=1.0d0
  END SELECT 

END FUNCTION RhoFun

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: lam,phi,r
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: rot_lam0,rot_phi0

  SELECT CASE (Problem)
    CASE ('Damm')
      IF (x<=Zero) THEN
        HeightFun=0.0d0
      ELSE
        HeightFun=(20.0d0*x+(700.0d0-x)*0.0d0)/700.0d0
      END IF
    CASE ('William2')
      HeightFun=0.0d0
    CASE ('William5')
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      CALL Rotate(lam0,phi0,rot_lam0,rot_phi0)
      r=SQRT(MIN(RadiusC*RadiusC, &
                 (rot_lam-rot_lam0)*(rot_lam-rot_lam0)+(rot_phi-rot_phi0)*(rot_phi-rot_phi0)))
      HeightFun=hS*(One-r/RadiusC)
    CASE ('Tsunami')
      phi=y
      HeightFun=(-40000.d0+20000.d0*EXP(-1000.d0*(phi+Pi/6.d0)**2.d0))/Grav  
    CASE DEFAULT
      HeightFun=0.0d0
  END SELECT
END FUNCTION HeightFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun

  ThStart=One

END FUNCTION ThStart

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  qvStart=0.0d0

END FUNCTION QvStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Shallow_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  qvProfFun=0.0d0

END FUNCTION QvProfFun


FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

 qcStart=0.0d0

END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKInd) :: lam,phi

  INTERFACE
    FUNCTION Simpson(x0,x1,dx,Integrand)
      USE Kind_Mod
      REAL(RealKind) :: Simpson
      REAL(RealKind) :: x0,x1,dx
      INTERFACE
        FUNCTION Integrand(tau)
          USE Kind_Mod
          REAL(RealKind) :: Integrand,tau
        END FUNCTION Integrand
      END INTERFACE
    END FUNCTION Simpson
  END INTERFACE
  INTERFACE
    FUNCTION Integrand4(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand4,tau
    END FUNCTION Integrand4
  END INTERFACE

  REAL(RealKind) :: phi04,f04,lam04,teta04,uQ,TimeLoc,C4,C4lam,C4phi,C4lamt
  REAL(RealKind) :: tetaQ,tetaQphi,tetaQlam,tetaQt,deruSt,C4phit,derhSlam,vS,C4t,f,uS,hS4,hQ
  REAL(RealKind) :: deruSlam,deruSphi,duSt,C4phiphi,C4lamphi
  REAL(RealKind) :: derhSt,derhSphi,dhSt


  SELECT CASE (Problem)
    CASE ('William4')
      lam=x
      phi=y
      lam04=0.0d0
      phi04=Pi/4.0d0
      f=Two*Omega4*SIN(phi)
      f04=2.0d0*Omega4*SIN(Pi/4.0d0)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*(SIN(2.0d0*phi))**14.0d0
      TimeLoc=Time
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1.0d0-C4)/(1.0d0+C4)))
      hQ=(GravH04-Simpson(-Half*pi,phi,Pi/100.0d0,Integrand4))/Grav
      hS4=(Grav*hQ+f*tetaQ)/Grav
      !derivatives of C(C4)
      C4lam=-COS(phi04)*COS(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phi=SIN(phi04)*COS(phi)-COS(phi04)*SIN(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4t=COS(phi04)*COS(phi)*SIN(lam-U04/RadEarth*TimeLoc-lam04)*U04/RadEarth
      C4phiphi=-SIN(phi04)*SIN(phi)-COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4lamphi=COS(phi04)*SIN(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phit=-COS(phi04)*SIN(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)*(U04/RadEarth)
      !derivatives of tetaQ-->uS,vS
      tetaQt=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4t
      tetaQlam=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4lam
      tetaQphi=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4phi
      uS=uQ-tetaQphi/RadEarth
      vS=tetaQlam/(RadEarth*COS(phi))
      !partial derivatives of u
      deruSphi=U04*28.0d0*(SIN(2.0d0*phi))**13.0d0*COS(2.0d0*phi) &
              -tetaQ*(4.0d0*sigma4*C4phi**2.0d0*(sigma4-1.0d0-C4)+(1.0d0+C4)**2.0d0*C4phiphi*2.0d0*sigma4) &
              /((1.0d0+C4)**4.0d0*RadEarth)
      deruSt=-tetaQ*(4.0d0*sigma4*C4t*C4phi*(sigma4-1.0d0-C4)+(1.0d0+C4)**2.0d0*C4phit*2.0d0*sigma4) &
            /((1.0d0+C4)**4.0d0*RadEarth)
      deruSlam=-tetaQ*(4.0d0*sigma4*C4lam*C4phi*(sigma4-1.0d0-C4)+(1.0d0*C4)**2.0d0*C4lamphi*2.0d0*sigma4) &
              /((1.0d0+C4)**4.0d0*RadEarth)
      !substantial derivative of u
      duSt=deruSt+uS/(RadEarth*COS(phi))*deruSlam+vS/RadEarth*deruSphi-uS*vS*TAN(phi)/RadEarth
      !partial derivatives of h
      derhSt=f*tetaQt/Grav
      derhSphi=(f*tetaQphi+2.0d0*Omega4*COS(phi)*tetaQ)/Grav-uQ*(RadEarth*f+uQ*TAN(phi))/Grav
      derhSlam=f*tetaQlam/Grav
      !substantial derivative of h
      dhSt=derhSt+uS/(RadEarth*COS(phi))*derhSlam+vS/RadEarth*derhSphi
      !forcing term 
      ForceU=hS4*duSt+uS*dhSt
  END SELECT



END FUNCTION ForceU


FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time

  INTERFACE
    FUNCTION Simpson(x0,x1,dx,Integrand)
      USE Kind_Mod
      REAL(RealKind) :: Simpson
      REAL(RealKind) :: x0,x1,dx
      INTERFACE
        FUNCTION Integrand(tau)
          USE Kind_Mod
          REAL(RealKind) :: Integrand,tau
        END FUNCTION Integrand
      END INTERFACE
    END FUNCTION Simpson
  END INTERFACE
  INTERFACE
    FUNCTION Integrand4(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand4,tau
    END FUNCTION Integrand4
  END INTERFACE

  REAL(RealKInd) :: lam,phi
  REAL(RealKind) :: TimeLoc,tetaQphi,tetaQlam,tetaQlamphi,C4lamt,C4t,tetaQ,tetaQt,dervSt,derhSphi,uS,f04
  REAL(RealKind) :: f,teta04,phi04,lam04,C4lam,C4phi,uQ,C4,dvSt,vS,dervSphi
  REAL(RealKind) :: C4lamlam,C4lamphi,dervSlam,hS4,hQ
  REAL(RealKind) :: derhSt,derhSlam,dhSt

 SELECT CASE (Problem)
    CASE ('William4')
      lam=x
      phi=y
      phi04=Pi/4.0d0
      lam04=0.0d0
      f04=2.0d0*Omega4*SIN(Pi/4.0d0)
      f=2.0d0*Omega4*SIN(phi)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*(SIN(2.0d0*phi))**14.0d0
      TimeLoc=Time
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1.0d0-C4)/(1.0d0+C4)))
      hQ=(GravH04-Simpson(-Half*pi,phi,Pi/100.0d0,Integrand4))/Grav
      hS4=(Grav*hQ+f*tetaQ)/Grav
      !derivatives of C(=C4)
      C4lam=-COS(phi04)*COS(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phi=SIN(phi04)*COS(phi)-COS(phi04)*SIN(phi) &
            *COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4t=COS(phi04)*COS(phi)*SIN(lam-U04/RadEarth*TimeLoc-lam04)*U04/RadEarth
      C4lamt=COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)*(U04/RadEarth)
      C4lamlam=-COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4lamphi=COS(phi04)*SIN(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      !derivatives of tetaQ-->uS,vS
      tetaQlam=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4lam
      tetaQphi=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4phi
      tetaQt=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4t
      tetaQlamphi=tetaQ*(2.0d0*C4lam*C4phi*(2.0d0*sigma4**2.0d0-1.0d0-C4)+(1.0d0+C4)**2.0d0*C4lamphi) &
                 /((1.0d0+C4)**4.0d0)
      uS=uQ-tetaQphi/RadEarth
      vS=tetaQlam/(RadEarth*COS(phi))
      !partial derivatives of v
      dervSlam=tetaQ*(4.0d0*sigma4*C4lam**2.0d0*(sigma4-1.0d0-C4)+(1.0d0+C4)**2.0d0*C4lamlam*2.0d0*sigma4) &
              /(RadEarth*COS(phi)*(1.0d0+C4)**4.0d0)
      dervSphi=(COS(phi)*tetaQlamphi+tetaQlam*SIN(phi))/(RadEarth*(COS(phi))**2.0d0)
      dervSt=tetaQ*(4.0d0*sigma4*C4t*C4lam*(sigma4-1.0d0-C4)+(1.0d0+C4)**2.0d0*C4lamt*2.0d0*sigma4) &
            /((1.0d0+C4)**4.0d0*(RadEarth*COS(phi)))
      !partial derivatives of h
      derhSt=f*tetaQt/Grav
      derhSphi=(f*tetaQphi+2.0d0*Omega4*COS(phi)*tetaQ)/Grav-uQ*(RadEarth*f+uQ*TAN(phi))/Grav
      derhSlam=f*tetaQlam/Grav
      !substantial derivative of v
      dvSt=dervSt+uS/(RadEarth*COS(phi))*dervSlam+vS/RadEarth*dervSphi &
          +uS**2*TAN(phi)/RadEarth+Grav/RadEarth*derhSphi+f*uS
      !substantial derivative of h
      dhSt=derhSt+uS/(RadEarth*COS(phi))*derhSlam+vS/RadEarth*derhSphi
      !forcing term
      ForceV=hS4*dvSt+vS*dhSt
 END SELECT



END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time

  ForceW=Zero

END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho,lam,phi
  REAL(RealKind) :: x,y,z,zHeight,Time

  INTERFACE
    FUNCTION Simpson(x0,x1,dx,Integrand)
      USE Kind_Mod
      REAL(RealKind) :: Simpson
      REAL(RealKind) :: x0,x1,dx
      INTERFACE
        FUNCTION Integrand(tau)
          USE Kind_Mod
          REAL(RealKind) :: Integrand,tau
        END FUNCTION Integrand
      END INTERFACE
    END FUNCTION Simpson
  END INTERFACE
  INTERFACE
    FUNCTION Integrand4(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand4,tau
    END FUNCTION Integrand4
  END INTERFACE

  REAL(RealKind) :: TimeLoc,hS4,C4lamphi,dervSphi,deruSlam
  REAL(RealKind) :: derhSt,tetaQ,tetaQt,hQ,f,teta04,uQ,f04,phi04,C4phi,C4t,C4lam,C4,lam04
  REAL(RealKind) :: dhSt,derhSlam,tetaQphi,tetaQlam,uS,vS,derhSphi
  SELECT CASE (Problem)
    CASE ('William4')
      lam=x
      phi=y
      phi04=Pi/4.0d0
      lam04=0.0d0
      f04=2.0d0*Omega4*SIN(Pi/4.0d0)
      f=2.0d0*Omega4*SIN(phi)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*(SIN(2.0d0*phi))**14.0d0
      TimeLoc=Time
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1.0d0-C4)/(1.0d0+C4)))
      hQ=(GravH04-Simpson(-Half*pi,phi,Pi/100.0d0,Integrand4))/Grav
      hS4=(Grav*hQ+f*tetaQ)/Grav
      !derivatives of C(=C4)
      C4lam=-COS(phi04)*COS(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phi=SIN(phi04)*COS(phi)-COS(phi04)*SIN(phi) &
            *COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4t=COS(phi04)*COS(phi)*SIN(lam-U04/RadEarth*TimeLoc-lam04)*U04/RadEarth
      C4lamphi=COS(phi04)*SIN(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      !derivatives of tetaQ-->uS,vS
      tetaQt=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4t
      tetaQlam=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4lam
      tetaQphi=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4phi
      uS=uQ-tetaQphi/RadEarth
      vS=tetaQlam/(RadEarth*COS(phi))
      !partial derivatives of h
      derhSt=f*tetaQt/Grav
      derhSphi=(f*tetaQphi+tetaQ*2.0d0*Omega4*COS(phi))/Grav-uQ*(RadEarth*f+uQ*TAN(phi))/Grav
      derhSlam=f*tetaQlam/Grav
      !substantial derivative of h
      dhSt=derhSt+uS/(RadEarth*COS(phi))*derhSlam+vS/RadEarth*derhSphi
      !forcing term
      ForceRho=dhSt
  END SELECT

END FUNCTION ForceRho



FUNCTION UStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: lam,phi,uPr,b0,b1,x0,xe,lamPr,phiPr,phib,phie,lampPr,f
  REAL(RealKind) :: VStart,hQ,t4,uS,vS,hS4
! REAL(RealKind) :: C4lam,C4phi,C4t,tetaQlam,tetaQphi,tetaQt
  REAL(RealKind) :: dervSt,deruSt,derhSt
  REAL(RealKind) :: deruSlam,dervSphi,derhSlam,derhSphi,FvR,FuR,FhR,alpha4
  REAL(RealKind) :: uM,eN,phi0G,phi1G
  REAL(RealKind) :: phi04,f04,lam04,teta04,uQ,TimeLoc,C4,C4lam,C4phi
  REAL(RealKind) :: tetaQ,tetaQPhi
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: u_lat,v_lat,v_tmp
  REAL(RealKind) :: yHat

  SELECT CASE (Problem)
    CASE ('Damm')
      UStart=uMax
    CASE ('William1')
      lam=x
      phi=y
      UStart=U1*(COS(alpha1)*COS(phi)+SIN(alpha1)*COS(lam)*SIN(phi))
    CASE ('William2') 
      phi=y
      UStart=uMax*COS(phi)
    CASE ('William3')
      lam=x
      phi=y
      phiPr=ASIN(SIN(phi)*COS(alpha3)-COS(phi)*COS(lam)*SIN(alpha3))
      IF (COS(alpha3)*COS(lam)*COS(phi)+SIN(alpha3)*SIN(phi)>=0) THEN
        lamPr=ASIN(SIN(lam)*COS(phi)/COS(phiPr))
      ELSE
        lamPr=Pi-ASIN(SIN(lam)*COS(phi)/COS(phiPr))
      END IF 
      x0=xe3*(phiPr-phib3)/(phie3-phib3)
      IF (x0<=Zero) THEN
         b0=Zero
      ELSE
         b0=EXP(-One/x0)
      END IF
      IF ((xe3-x0)<=Zero) THEN
         b1=Zero
      ELSE
         b1=EXP(-One/(xe3-x0))
      END IF
      uPr=U03*b0*b1*EXP(4.0d0/xe3)
      VStart=(-uPr*SIN(alpha3)*SIN(lamPr))/COS(phi)
      IF (ABS(COS(lam))<=1-10.0d0) THEN
      END IF
      UStart=(VStart*SIN(phi)*SIN(lam)+uPr*COS(lamPr))/COS(lam)
    CASE ('William4')
      lam=x
      phi=y
      phi04=Pi/4.0d0
      lam04=0.0d0
      f04=2.0d0*Omega4*SIN(Pi/4.0d0)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*((SIN(2.0d0*phi))**14.0d0)
      TimeLoc=0.0d0 
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4lam=-COS(phi04)*COS(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phi=SIN(phi04)*COS(phi)-COS(phi04)*SIN(phi) &
            *COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1.0d0-C4)/(1.0d0+C4)))
      tetaQphi=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4phi
      UStart=uQ-tetaQphi/RadEarth 
    CASE ('William5') 
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      u_lat=uMax*COS(rot_phi)
      v_lat=0.0d0
      IF (ABS(rotation_angle)<1.0E-8) THEN
        UStart = u_lat
      ELSE
        ! rotate wind components
        CALL turnwi(u_lat,v_lat,UStart,v_tmp,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        IF (ABS(UStart)<1.0E-10) UStart=0.0d0
      ENDIF
    CASE ('William6')
      lam=x
      phi=y
      UStart=RadEarth*omega6*COS(phi) &
            +RadEarth*K6*(COS(phi))**(R6-One) &
            *(R6*SIN(phi)*SIN(phi)-COS(phi)*COS(phi))*COS(R6*lam)
    CASE ('Jet')
      yHat=(y-yCJet)/WidthJet
      IF (ABS(yHat)<=1.0d0) THEN
        UStart=u0Jet*(1.0d0-3.0d0*yHat**2.0d0+3.0d0*yHat**4.0d0-yHat**6.0d0)
      ELSE  
        UStart=0.0d0
      END IF  
    CASE ('Galewsky')
      phi=y
      uM=80.0d0
      phi0G=Pi/7.0d0
      phi1G=Pi/2.0d0-phi0G
      eN=EXP(-4.0d0/(phi1G-phi0G)**2.0d0)
      IF ((phi<=phi0G).OR.(phi>=phi1G)) THEN
         UStart=0
      ELSE
         UStart=uM/eN*EXP(1.0d0/((phi-phi0G)*(phi-phi1G)))
      END IF
    CASE ('Dancing')
      lam=x
      phi=y
      UStart=-vMB*(Three*SIN(phi)*COS(phi)*COS(phi)-(SIN(phi))**Three)*SIN(lam)
    CASE('Tsunami')
      UStart=Zero
    CASE ('Twin') 
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      u_lat=4.0d0*uMax*COS(rot_phi)*(1.0d0-COS(rot_phi))
      v_lat=0.0d0
      IF (ABS(rotation_angle)<1.0E-8) THEN
        UStart = u_lat
      ELSE
        ! rotate wind components
        CALL turnwi(u_lat,v_lat,UStart,v_tmp,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        IF (ABS(UStart)<1.0E-10) UStart=0.0d0
      ENDIF
    CASE DEFAULT
      UStart=Zero
  END SELECT 

END FUNCTION UStart

 

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: uPr,b0,b1,x0,xe,lamPr,phiPr,phib,phie,lampPr
  REAL(RealKind) :: SurPot,GeoHei 
  REAL(RealKind) :: phi04,f04,lam04,teta04,uQ,TimeLoc,C4,C4lam,C4phi
  REAL(RealKind) :: tetaQ,tetaQlam,tetaQPhi
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: u_lat,v_lat,u_tmp


  SELECT CASE (Problem)
    CASE ('William1')
      lam=x
      phi=y
      VStart=-U1*SIN(alpha1)*SIN(lam)
    CASE ('William3')
      lam=x
      phi=y
      phiPr=ASIN(SIN(phi)*COS(alpha3)-COS(phi)*COS(lam)*SIN(alpha3))
      IF (COS(alpha3)*COS(lam)*COS(phi)+SIN(alpha3)*SIN(phi)>=0) THEN
        lamPr=ASIN(SIN(lam)*COS(phi)/COS(phiPr))
      ELSE
        lamPr=Pi-ASIN(SIN(lam)*COS(phi)/COS(phiPr))
      END IF 
      x0=xe3*(phiPr-phib3)/(phie3-phib3)
      IF (x0<=Zero) THEN
         b0=Zero
      ELSE
         b0=EXP(-One/x0)
      END IF
      IF ((xe3-x0)<=Zero) THEN
         b1=Zero
      ELSE
         b1=EXP(-One/(xe3-x0))
      END IF
      uPr=U03*b0*b1*EXP(4.0d0/xe3)
      VStart=(-uPr*SIN(alpha3)*SIN(lamPr))/COS(phi)
    CASE ('William4')
      lam=x
      phi=y
      phi04=Pi/4.0d0
      lam04=0.0d0
      f04=2.0d0*Omega4*SIN(Pi/4.0d0)
      teta04=-0.03d0*(GravH04/f04)
      uQ=U04*(SIN(2.0d0*phi))**14.0d0
      TimeLoc=0.0d0
      C4=SIN(phi04)*SIN(phi)+COS(phi04)*COS(phi)*COS(lam-U04*TimeLoc/RadEarth-lam04)
      C4lam=-COS(phi04)*COS(phi)*SIN(lam-U04*TimeLoc/RadEarth-lam04)
      C4phi=SIN(phi04)*COS(phi)-COS(phi04)*SIN(phi) &
            *COS(lam-U04*TimeLoc/RadEarth-lam04)
      tetaQ=teta04*EXP(-sigma4*((1.0d0-C4)/(1.0d0+C4)))
      tetaQlam=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4lam
      tetaQphi=2.0d0*sigma4/((1.0d0+C4)*(1.0d0+C4))*tetaQ*C4phi
      VStart=tetaQlam/(RadEarth*COS(phi))
    CASE ('William5')
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      u_lat=uMax*COS(rot_phi)
      IF (ABS(rotation_angle)<1.0E-8) THEN
        VStart = 0.0d0
      ELSE
        v_lat=0.0d0
        ! pole point velocities are not well-defined
        IF (ABS(Pi*0.5d0-phi)<1.0E-8.OR.ABS(Pi*0.5d0+phi)<1.0E-8) THEN
          VStart = 0.0d0
        ELSE
          ! rotate wind components
          CALL turnwi(u_lat,v_lat,u_tmp,VStart,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        ENDIF
      ENDIF
    CASE ('William6')
      lam=x
      phi=y
      VStart=-RadEarth*K6*R6*(COS(phi))**(R6-One)*SIN(phi)*SIN(R6*lam)
    CASE ('Jet')
      VStart=Zero
    CASE ('Galewsky')
      VStart=Zero
    CASE DEFAULT
      VStart=Zero
    CASE ('Dancing')
      lam=x
      phi=y
      VStart=vMB*SIN(phi)*SIN(phi)*COS(lam)
    CASE ('Tsunami')
      phi=y
      SurPot=-40000.d0+20000.d0*EXP(-1000.d0*(phi+Pi/6.d0)**2.d0)  
      GeoHei=(-SurPot+AA*EXP(-1000.d0*(phi-Pi/6.d0)**2.d0))
      VStart=-200.0d0*AA*EXP(-1000.0d0*(phi-Pi/6.0d0)*(phi-Pi/6.0d0))/GeoHei
    CASE ('Twin')
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      u_lat=4.0d0*uMax*COS(rot_phi)*(1.0d0-COS(rot_phi))
      IF (ABS(rotation_angle)<1.0E-8) THEN
        VStart = 0.0d0
      ELSE
        v_lat=0.0d0
        ! pole point velocities are not well-defined
        IF (ABS(Pi*0.5d0-phi)<1.0E-8.OR.ABS(Pi*0.5d0+phi)<1.0E-8) THEN
          VStart = 0.0d0
        ELSE
          ! rotate wind components
          CALL turnwi(u_lat,v_lat,u_tmp,VStart,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        ENDIF
      ENDIF
  END SELECT 

END FUNCTION VStart

FUNCTION DStart(x,y,z,zHeight,Time)

  USE Shallow_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Shallow_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: lam,phi,rLoc

  SELECT CASE (Problem)
    CASE ('William1')
      lam=x
      phi=y
      rLoc=DistGreatCircle(lam1,lam,phi1,phi)
      IF (rLoc<=R1) THEN
        DummyStart=Half*H01*(One+COS(Pi*rLoc/R1))
      ELSE
        DummyStart=0.0d0
      END IF
    CASE DEFAULT
      DummyStart=0.0d0
  END SELECT   
CONTAINS
FUNCTION DistGreatCircle(lambda0,lambda1,phi0,phi1)
    REAL(8) :: DistGreatCircle
    REAL(8) :: lambda0,lambda1,phi0,phi1

    DistGreatCircle = RadEarth*ACOS(SIN(phi0)*SIN(phi1)+COS(phi0)*COS(phi1)*COS(lambda1-lambda0))
END FUNCTION DistGreatCircle
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=Zero
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.0d0
END FUNCTION TStart

FUNCTION Simpson(x0,x1,dx,Integrand)
  USE Kind_Mod
  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: Simpson

  REAL(RealKind) :: x0,x1,dx
  INTERFACE
    FUNCTION Integrand(tau)
      USE Kind_Mod
      REAL(RealKind) :: Integrand,tau
    END FUNCTION Integrand
  END INTERFACE

  INTEGER :: i,n
  REAL(RealKind) :: xi,h

  n=(x1-x0)/dx+1
  h=(x1-x0)/n
  Simpson=Half*(Integrand(x0)+Integrand(x1))
  xi=x0
  DO i=1,n-1
    xi=xi+h
    Simpson=Simpson+Integrand(xi)
  END DO
  xi=x0-Half*h
  DO i=1,n
    xi=xi+h
    Simpson=Simpson+Two*Integrand(xi)
  END DO
  Simpson=h/3.0d0*Simpson

END FUNCTION Simpson

FUNCTION Integrand3(tau)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Integrand3

  REAL(RealKind) :: tau

  REAL(RealKind) :: x0,b0,b1,uPr

  x0=xe3*(tau-phib3)/(phie3-phib3)
  IF (x0<=Zero) THEN
     b0=Zero
  ELSE
     b0=EXP(-One/x0)
  END IF
  IF ((xe3-x0)<=Zero) THEN
     b1=Zero
  ELSE
     b1=EXP(-One/(xe3-x0))
  END IF
  uPr=U03*b0*b1*EXP(4.0d0/xe3)
  IF (ABS(tau)<Half*Pi) THEN
    Integrand3=(Two*Omega*SIN(tau)+uPr*TAN(tau)/RadEarth)*uPr
  ELSE
    Integrand3=Zero
  END IF 

END FUNCTION Integrand3

FUNCTION Integrand4(tau)

  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Integrand4

  REAL(RealKind) :: tau

  REAL(RealKind) :: uQ,f

  uQ=U04*(SIN(2.0d0*tau))**14.0d0
  f=2.0d0*Omega4*SIN(tau)
  IF (ABS(tau)<Half*Pi) THEN
    Integrand4=(RadEarth*f+uQ*TAN(tau))*uQ
  ELSE
    Integrand4=Zero
  END IF
  
END FUNCTION Integrand4

FUNCTION IntegrandG(tau)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: IntegrandG

  REAL(RealKind) :: tau

  REAL(RealKind) :: UStart,f,phi0G,phi1G,uM,eN

  uM=80.0d0
  phi0G=Pi/7.0d0
  phi1G=Pi/2.0d0-phi0G
  eN=EXP(-4.0d0/(phi1G-phi0G)**2.0d0)
  f=2*Omega4*SIN(tau)
  IF ((tau<=phi0G).OR.(tau>=phi1G)) THEN
    UStart=0.0d0
  ELSE
    UStart=uM/eN*EXP(1.0d0/((tau-phi0G)*(tau-phi1G)))
  END IF
  IF (ABS(tau)<Half*Pi) THEN
    IntegrandG=(RadEarth*f+UStart*TAN(tau))*UStart
  ELSE
    IntegrandG=Zero
  END IF

END FUNCTION IntegrandG

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ThStartSoil=0.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4


FUNCTION NvStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NvStart = 0.0
END FUNCTION NvStart

FUNCTION NcStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NcStart = 0.0
END FUNCTION NcStart

FUNCTION NrStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NrStart = 0.0
END FUNCTION NrStart

FUNCTION NiStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NiStart = 0.0
END FUNCTION NiStart

FUNCTION NsStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NsStart = 0.0
END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart

FUNCTION QsStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QsStart = 0.0
END FUNCTION QsStart

FUNCTION EnStart(lam,phi,z,zHeight,Time)
  USE Shallow_Mod
  USE ThProf_Mod
  REAL(RealKind) :: EnStart
  EnStart=Zero
END FUNCTION EnStart


SUBROUTINE Rotate(lam,phi,rot_lam,rot_phi)

  USE Shallow_Mod
  IMPLICIT NONE!

  REAL(RealKind) :: lam,phi,rot_lam,rot_phi

  IF (ABS(rotation_angle)<1.0E-8) THEN
    rot_lam = lam
    rot_phi = phi
  ELSE
     CALL regrot(lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*pi+rotation_angle,1)
  ENDIF
END SUBROUTINE Rotate

SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
  USE Shallow_Mod
  IMPLICIT NONE!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
  INTEGER :: kxdim,kydim,kx,ky,kcall
  REAL(RealKind) :: pxreg,pyreg,&
                    pxrot,pyrot,&
                    pxcen,pycen
  REAL(RealKind) :: zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
                    zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
  INTEGER  :: jy,jx

  zpih = Pi*0.5d0
  zsycen = SIN((pycen+zpih))
  zcycen = COS((pycen+zpih))
!
  IF (kcall.EQ.1) THEN
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
     zsyrot = max(zsyrot,-1.d0)
     zsyrot = min(zsyrot,+1.d0)

     pyrot = ASIN(zsyrot)

     zcyrot = COS(pyrot)
     zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
     zcxrot = max(zcxrot,-1.d0)
     zcxrot = min(zcxrot,+1.d0)
     zsxrot = zcyreg*zsxmxc/zcyrot

     pxrot = ACOS(zcxrot)

     IF (zsxrot<0.0) pxrot = -pxrot
  ELSE IF (kcall.EQ.-1) THEN
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
     zsyreg = max(zsyreg,-1.d0)
     zsyreg = min(zsyreg,+1.d0)

     pyreg = ASIN(zsyreg)

     zcyreg = COS(pyreg)
     zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
     zcxmxc = max(zcxmxc,-1.d0)
     zcxmxc = min(zcxmxc,+1.d0)
     zsxmxc = zcyrot*zsxrot/zcyreg
     zxmxc  = ACOS(zcxmxc)
     IF (zsxmxc<0.0) zxmxc = -zxmxc

     pxreg = zxmxc + pxcen

  ELSE
     WRITE(6,'(1x,''invalid kcall in regrot'')')
     STOP
  ENDIF
END SUBROUTINE regrot

SUBROUTINE turnwi(puarg,pvarg,pures,pvres,         &
                  pxreg,pyreg,pxrot,pyrot,   &
                  pxcen,pycen,kcall)
!
! USE Shallow_Mod
  IMPLICIT NONE
  INTEGER, PARAMETER :: RealKind=8
  REAL(RealKind) :: Pi
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    pa    : transformation coefficients
!*    pb    :    -"-
!*    pc    :    -"-
!*    pd    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

  INTEGER kxdim,kydim,kx,ky,kcall
  REAL(RealKind) puarg,pvarg,    &
                 pures,pvres,    &
                 pa,   pb,       &
                 pc,   pd,       &
                 pxreg,pyreg,    &
                 pxrot,pyrot
  REAL(RealKind) pxcen,pycen
!-----------------------------------------------------------------------
  INTEGER jy,jx
  REAL(RealKind) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
                 zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!-----------------------------------------------------------------------
  Pi=4.0d0*ATAN(1.0d0)
  IF (kcall.EQ.1) THEN
     zpih = Pi*0.5d0
     zsyc = SIN(pycen+zpih)
     zcyc = COS(pycen+zpih)
     !
     zsxreg = SIN(pxreg)
     zcxreg = COS(pxreg)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     !
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     !
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     !
     pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
     pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
          zsxmxc*zsyreg*zcxrot
     pc = zsyc*zsxmxc/zcyrot
     pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
     !
     pures = pa*puarg + pb*pvarg
     pvres = pc*puarg + pd*pvarg
  ELSE IF (kcall.EQ.-1) THEN
     zpih = Pi*0.5d0
     zsyc = SIN(pycen+zpih)
     zcyc = COS(pycen+zpih)
     !
     zsxreg = SIN(pxreg)
     zcxreg = COS(pxreg)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     !
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     !
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     !
     pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
     pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
          zcxmxc*zsxrot*zsyrot
     pc =-zsyc*zsxrot/zcyreg
     pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
     !
     pures = pa*puarg + pb*pvarg
     pvres = pc*puarg + pd*pvarg
  ELSE
     WRITE(6,'(1x,''invalid kcall in turnwi'')')
     STOP
  ENDIF
END SUBROUTINE turnwi

