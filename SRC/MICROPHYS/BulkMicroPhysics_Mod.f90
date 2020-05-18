MODULE BulkMicroPhysics_Mod

  USE Control_Mod
  USE Parameter_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Chemie_Mod

  IMPLICIT NONE

CONTAINS 

SUBROUTINE MicroProcess(c,f,Process,Task,lprint) 

  REAL(RealKind) :: c(:),f(:)
  CHARACTER*8 :: Process,Task
  LOGICAL, OPTIONAL :: lprint

  REAL(RealKind) :: RhoD,p,qD,qV,qL,qR,T
  REAL(RealKind) :: pVs,qVs,Rm,Cpml,Cvml,Lv,Cp_eff,RDens
  REAL(RealKind) :: Qcond,Qauto,Qacc,Qdep
  REAL(RealKind) :: Rho,RhoV,RhoL,RhoR,RhoI,RhoS,PotM,Kinetic,zh,eLoc
  REAL(RealKind) :: RhsRhoV,RhsRhoL,RhsPotM
  REAL(RealKind) :: DpDqv,DpDql,DpDRhoV,DpDRhoL 

  REAL(RealKind) :: Nd,Dd,Psi,Dia,Mass,SedVel,Rey
  REAL(RealKind) :: Lambda,Nmean,Mmean,VentFactor,ThermoG
  REAL(RealKind) :: DmDtAcc,DmDtDep 
  REAL(RealKind) :: temptest 

  REAL(RealKind), PARAMETER :: ar=523.5987756d0 
  REAL(RealKind), PARAMETER :: br=3.0d0 
  REAL(RealKind), PARAMETER :: cr=130.0d0 
  REAL(RealKind), PARAMETER :: dr=0.5d0 
  REAL(RealKind), PARAMETER :: as=2.5d-2 
  REAL(RealKind), PARAMETER :: bs=2.0d0 
  REAL(RealKind), PARAMETER :: cs=4.0d0 
  REAL(RealKind), PARAMETER :: ds=0.25d0 
  REAL(RealKind), PARAMETER :: N0=1.0d7 
  REAL(RealKind), PARAMETER :: Er=0.8d0 
  REAL(RealKind), PARAMETER :: Alphar=One 
  REAL(RealKind), PARAMETER :: Betar=Two
  REAL(RealKind), PARAMETER :: Es=0.2d0 
  REAL(RealKind), PARAMETER :: Alphas=0.3d0 
  REAL(RealKind), PARAMETER :: Betas=Three



  SELECT CASE(Process)
    CASE('COND') ! condensation/evaporation process (with present RhoV,RhoC,RhoR)
      Rho  = c(1)
      RhoV = c(2)
      RhoL = c(3)
      PotM = c(4)
      RhoR = c(5)
      RhoL = RhoL + RhoR + Eps
      RhoI = Zero
      RhoD = Rho - RhoV - RhoL - RhoI + Eps
      qD   = RhoD / (Rho+Eps)
      qV   = RhoV / (Rho+Eps)
      qL   = RhoL / (Rho+Eps)
      Cpml = Cpd*RhoD + Cpv*RhoV + Cpl*RhoL + RhoI*Cpi+Eps
      Cvml = Cvd*RhoD + Cvv*RhoV + Cpl*RhoL + RhoI*Cpi+Eps
      Rm   = Rd*RhoD + Rv*RhoV + Eps

      IF ( Task == 'Function' ) THEN
        T = c(6)
        p = c(7)
      ELSE
        SELECT CASE(ThetaKind)
          CASE('Density')
            p = PressureTheta(RhoD,RhoV,RhoL,RhoI,PotM) + Eps
            T = AbsTemp(RhoD,RhoV,p) + Eps
          CASE('Equiv')
            T =c (6)
            CALL AbsTNewton(T,PotM/(Rho+Eps),RhoD,RhoV,RhoL)
          CASE('Energy')
            Kinetic = c(8)
            zh = c(9)
            T  = AbsTEnergy(PotM+Eps,Rho+Eps,RhoV,RhoL,Kinetic,zh)
          CASE('PreEn')
            p = PotM 
            T = p / ((RhoD*Rd + RhoV*Rv) + Eps)
        END SELECT
      END IF
  
      pVs   = SaturVapor(T)
      RhoL  = c(3) 
      Qcond = RelCloud * ( (pVs/(Rv*T+Eps)-RhoV) + RhoL - &
            SQRT((pVs/(Rv*T+Eps)-RhoV) * (pVs/(Rv*T+Eps)-RhoV)+RhoL*RhoL) )
      RhoL=RhoL+RhoR+Eps
  
      SELECT CASE(ThetaKind)
        CASE('Density')
          Lv      = LatHeat(T)
          RhsPotM = PotM*( (-Lv/(Cpml*T+Eps) &
                        -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
                        +Rv/Rm                                   &
                        )*Qcond                                   &
                        +(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
                        )*(-Qcond)                                   &
                      )
        CASE('Equiv')
          Cp_eff  = Cpd + (RhoV+RhoL)/RhoD*Cpl
          RhsPotM = -PotM/(Cp_eff*RhoD)*Rv*LOG(RelHumidity(T,RhoV))*Qcond
        CASE('Energy')
          RhsPotM = Zero
        CASE('PreEn')
          eLoc    = Cvml*T+RhoV*L00
          DpDRhoV = (Rv-Rd)*(eLoc-RhoV*L00)/Cvml &
                    -L00*(RhoD*Rd+RhoV*Rv)/Cvml &
                    -(Cvv-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
          DpDRhoL = -Rd*(eLoc-RhoV*L00)/Cvml &
                    -(Cpl-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
          RhsPotM = (DpDRhoV-DpDRhoL)*Qcond
      END SELECT

      RhsRhoV = Qcond
      RhsRhoL = -Qcond
      f(1) = RhsRhoV
      f(2) = RhsRhoL
      f(3) = RhsPotM
    CASE('AUTO') ! autoconversion process
      RhoL=c(1) 
      Psi=1.0d3*RhoL                             ! RhoR in g m^(-3)
      Nd=200.0d0                                 ! 50=maritime clouds, 2000=continental clouds
      Dd=0.146d0-5.964d0*1.0d-2*LOG(Nd/2000.0d0) ! relative dispersion
      IF (RhoL<Zero) THEN
        Qauto=Zero
      ELSE
        Qauto=1.67d0*1.0d-5*Psi**3/(5.0d0*Psi+0.036d0*Nd/Dd)
      END IF
      f(1)=-Qauto*1.0d0
      f(2)=Qauto*1.0d0
      IF (PRESENT(lprint)) THEN
        WRITE(*,*) 'Qauto,RhoL ',Qauto,RhoL
      END IF
    CASE('ACC') ! accretion process
      RhoL=c(1)
      RhoR=c(2)
      IF (RhoL<Zero) THEN
        RhoL=Zero
      END IF
      IF (RhoR<Zero) THEN
        RhoR=Zero
      END IF
      Lambda=((ar*N0*EXP(gammln(br+One)))/(RhoR+Eps))**(One/(br+One))
      Nmean=N0/Lambda   ! mean concentration of precipitating particles
      Mmean=RhoR/Nmean  ! mean mass of a precipitating particle
      Dia=(Mmean/ar)**(One/br) 
      SedVel=cr*Dia**dr
      DmDtAcc=Pi/Four*Dia**Two*SedVel*Er*Alphar*RhoL
      IF (c(1)<Zero) THEN
        Qacc=Zero
      ELSE
        Qacc=Nmean*DmDtAcc
      END IF
      f(1)=-Qacc
      f(2)=Qacc
      IF (PRESENT(lprint)) THEN
        WRITE(*,*) 'Qacc:',Qacc
      END IF
    CASE('DEP') ! deposition process
      Rho=c(1)
      RhoV=c(2)
      RhoL=c(3)
      PotM=c(4)
      RhoR=c(5)
      IF (RhoL<Zero) THEN
        RhoL=Zero
      END IF
      IF (RhoR<Zero) THEN
        RhoR=Zero
      END IF
      RhoL=RhoL+RhoR+Eps
      RhoI=Zero
      RhoD=Rho-RhoV-RhoL+Eps
      qD=RhoD/(Rho+Eps)
      qV=RhoV/(Rho+Eps)
      qL=RhoL/(Rho+Eps)
      Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL+Eps
      Cvml=Cvd*RhoD+Cvv*RhoV+Cpl*RhoL+Eps
      Rm=Rd*RhoD+Rv*RhoV+Eps
      IF (Task=='Function') THEN
        T=c(6)
        p=c(7)
      ELSE
        SELECT CASE(ThetaKind)
        CASE('Density')
          p=PressureTheta(RhoD,RhoV,RhoL,RhoI,PotM)+Eps
          T=AbsTemp(RhoD,RhoV,p)+Eps
        CASE('Equiv')
          T=c(6)
          CALL AbsTNewton(T,PotM/(Rho+Eps),RhoD,RhoV,RhoL)
        CASE('Energy')
          Kinetic=c(8)
          zh=c(9)
          T=AbsTEnergy(PotM+Eps,Rho+Eps,RhoV,RhoL,Kinetic,zh)
        CASE('PreEn')
          p=PotM 
          T=p/((RhoD*Rd+RhoV*Rv)+Eps)
        END SELECT
      END IF
      pVs=SaturVapor(T)
      qVs=Rd/Rv*pVs/(p-pVs+Eps+Eps)
      IF (Rho>Zero) THEN
        Lambda=((ar*N0*EXP(gammln(br+One)))/(RhoR+Eps))**(One/(br+One))
      ELSE
        Lambda=Eps
      END IF 
      Nmean=N0/Lambda   ! mean concentration of precipitating particles
      Mmean=RhoR/Nmean  ! mean mass of a precipitating particle
      Dia=(Mmean/ar)**(One/br) 
      SedVel=cr*Dia**dr
      Rey=Dia*SedVel/2.0d-5
      ThermoG=1.0d-7*(2.2d0*T/(pVs+Eps)+2.2d2/(T+Eps))**(-One)
      DmDtDep=Four*Pi*Dia/Betar*(qV/(qVs+Eps)-One)*(0.78d0+0.27d0*Rey**(Half))*ThermoG
      Qdep=Nmean*DmDtDep
      IF (c(5)<Zero) Qdep=Zero 
      SELECT CASE(ThetaKind)
        CASE('Density')
          Lv=LatHeat(T)
          RhsPotM=PotM*( &
                        (-Lv/(Cpml*T+Eps) &
                        -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
                        +Rv/Rm                                   &
                        )*Qdep                                   &
                        +(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
                        )*(-Qdep)                                   &
                      )
        CASE('Equiv')
          Cp_eff=Cpd+(RhoV+RhoL)/RhoD*Cpl
          RhsPotM=-PotM/(Cp_eff*RhoD)*Rv*LOG(RelHumidity(T,RhoV))*Qdep
        CASE('Energy')
          RhsPotM=Zero
        CASE('PreEn') 
          eLoc=Cvml*T+RhoV*L00
          DpDRhoV=(Rv-Rd)*(eLoc-RhoV*L00)/Cvml &
                  -L00*(RhoD*Rd+RhoV*Rv)/Cvml &
                  -(Cvv-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
          DpDRhoL=-Rd*(eLoc-RhoV*L00)/Cvml &
                  -(Cpl-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
          RhsPotM=(DpDRhoV-DpDRhoL)*Qdep
      END SELECT
      f(1)=-Qdep ! RhoV
      f(2)=Qdep  ! RhoR
      f(3)=-RhsPotM
      IF (PRESENT(lprint)) THEN
        WRITE(*,*) 'Qdep:',Qdep
        WRITE(*,*) 'RhsPotM:',RhsPotM
      END IF
    CASE DEFAULT
      IF (MyId==0) THEN
        WRITE(*,*) '  Error in .chem file - wrong process: Process :: ', TRIM(Process)
        WRITE(*,*) '  Stopped in BulkMicroPhysics_Mod.'
        STOP  
      END IF
  END SELECT

END SUBROUTINE MicroProcess

SUBROUTINE BulkMicroJac(Gas,Jac) 
  TYPE(Vec4_T) :: Gas(:)
  TYPE(Vec4_T) :: Jac(:)

  INTEGER :: ix,iy,iz
  INTEGER :: istr,k,l,idf,iSpecies,iReak
  REAL(RealKind) :: Temp,Df
  REAL(RealKind) :: c(20),f(20),fC(20)
  TYPE(Reaction_T), POINTER :: Reaction

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Reaction=>MicroFirst
        DO iReak=1,nReakMicro 
          istr=0
          DO k=1,Reaction%NumSpeciesLeftAktiv
            iSpecies=Reaction%SpeciesLeft(k)
            c(k)=Gas(iSpecies)%c(ix,iy,iz,1)
            IF (iSpecies==RhoPos) THEN
              c(k)=c(k)+RhoProfG(ibLoc)%c(ix,iy,iz,1) 
            END IF
          END DO
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2)=PreCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+3)=KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4)=zP(iz-1)+zP(iz)
          CALL MicroProcess(c,fC,Reaction%TypeR,'Function')
          DO k=1,Reaction%NumSpeciesLeftAktiv
            Temp=c(k)
            c(k)=c(k)*(1.0d0+1.e-8_RealKind)+SIGN(1.e-8_RealKind,c(k))
            CALL MicroProcess(c,f,Reaction%TypeR,'Jacobi  ')
            DO l=1,Reaction%NumSpecies
              istr=istr+1
              idf=Reaction%Struct(istr)
              Df=(f(l)-fC(l))/(c(k)-Temp+Eps)
              Jac(idf)%c(ix,iy,iz,1)=Jac(idf)%c(ix,iy,iz,1)+Df
            END DO
            c(k)=Temp 
          END DO
          Reaction=>Reaction%Next
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE BulkMicroJac

SUBROUTINE BulkMicro(Gas,fGas) 
  ! Eps = 1.0d-40 from Parameter_Mod 

  TYPE(Vec4_T) :: Gas(:),fGas(:)

  INTEGER :: ix,iy,iz
  INTEGER :: k,iSpecies,iReak
  TYPE(Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: c(20),f(20),Rho
  LOGICAL :: printf

  ! LOOP over all grid cells
  DO ix = ix0+1 , ix1
    DO iy = iy0+1 , iy1
      DO iz = iz0+1 , iz1
        
        Reaction => MicroFirst
        
        DO iReak = 1 , nReakMicro 
      
          DO k = 1 , Reaction%NumSpeciesLeftAktiv
            iSpecies = Reaction%SpeciesLeft(k)
            c(k)     = Gas(iSpecies)%c(ix,iy,iz,1)
            
            IF (iSpecies==RhoPos) c(k) = c(k) + RhoProfG(ibLoc)%c(ix,iy,iz,1) 
          END DO
      
          c(Reaction%NumSpeciesLeftAktiv+1) = TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2) = PreCell(ibLoc)%c(ix,iy,iz,1)
! OSSI    c(Reaction%NumSpeciesLeftAktiv+3) = KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4) = zP(iz-1) + zP(iz)
      
          CALL MicroProcess(c,f,Reaction%TypeR,'Function')
      
          DO k = 1 , Reaction%NumSpecies
            iSpecies = Reaction%Species(k)%Species
            fGas(iSpecies)%c(ix,iy,iz,1) = fGas(iSpecies)%c(ix,iy,iz,1) + f(k)
          END DO
      
          Reaction => Reaction%Next  ! iReac = iReac + 1
      
        END DO
      
      END DO
    END DO
  END DO

END SUBROUTINE BulkMicro

END MODULE BulkMicroPhysics_Mod
