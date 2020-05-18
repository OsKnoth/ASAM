MODULE LSCMicroPhysics_Mod

  USE Control_Mod
  USE Parameter_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Chemie_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: ShiftInd

CONTAINS 

SUBROUTINE MicroProcess(c,f,Process,Task,lprint) 

  REAL(RealKind) :: c(:),f(:)
  CHARACTER*8 :: Process,Task
  LOGICAL, OPTIONAL :: lprint

  REAL(RealKind) :: RhoD,p,qD,qV,qL,qR,T
  REAL(RealKind) :: pVs,qVs,Rm,Cpml,Cvml,Lv,Cp_eff,RDens
  REAL(RealKind) :: Qcond,Qauto,Qacc,Qdep
  REAL(RealKind) :: Rho,RhoV,RhoL,RhoR,RhoI,PotM,Kinetic,zh,eLoc
  REAL(RealKind) :: RhsRhoV,RhsRhoL,RhsPotM
  REAL(RealKind) :: DpDqv,DpDql,DpDRhoV,DpDRhoL 

  REAL(RealKind) :: Nd,Dd,Psi,Dia,Mass,SedVel,Rey
  REAL(RealKind) :: Lambda,Nmean,Mmean,VentFactor,ThermoG
  REAL(RealKind) :: DmDtAcc,DmDtDep 

  SELECT CASE(Process)
  CASE('COND') 
    Rho=c(1)
    RhoV=c(2)
    PotM=c(3)
    RhoL=Zero
    RhoI=Zero
    RhoD=Rho-RhoV+Eps
    Cpml=Cpd*RhoD+Cpv*RhoV+Eps
    Cvml=Cvd*RhoD+Cvv*RhoV+Eps
    Rm=Rd*RhoD+Rv*RhoV+Eps
    IF (Task=='Function') THEN
      T=c(ShiftInd+1)
      p=c(ShiftInd+2)
    ELSE
      SELECT CASE(ThetaKind)
      CASE('Density')
        p=PressureTheta(RhoD,RhoV,RhoL,RhoI,PotM)+Eps
        T=AbsTemp(RhoD,RhoV,p)+Eps
      CASE('Equiv')
        T=c(ShiftInd+1)
        CALL AbsTNewton(T,PotM/(Rho+Eps),RhoD,RhoV,RhoL)
      CASE('Energy')
        Kinetic=c(ShiftInd+3)
        zh=c(ShiftInd+4)
        T=AbsTEnergy(PotM+Eps,Rho+Eps,RhoV,RhoL,Kinetic,zh)
      CASE('PreEn')
        p=PotM 
        T=p/((RhoD*Rd+RhoV*Rv)+Eps)
      END SELECT
    END IF
    pVs=SaturVapor(T)
    Qcond=RelCloud*((pVs/(Rv*T+Eps)-RhoV)+RhoL- &
          SQRT((pVs/(Rv*T+Eps)-RhoV)*(pVs/(Rv*T+Eps)-RhoV)+RhoL*RhoL))
!   WRITE (*,*) 'Qcond,p,T,pVs/(Rv*T+Eps),RhoV',Qcond,p,T,pVs/(Rv*T+Eps),RhoV
    SELECT CASE(ThetaKind)
    CASE('Density')
      Lv=LatHeat(T)
      RhsPotM=PotM*(Rv/Rm-LOG((p+Eps)/p0)*(Rv/Rm-Cpv/Cpml)-Lv/(Cpml*T+Eps))*Qcond
!     RhsPotM=PotM*( &
!                   (-Lv/(Cpml*T+Eps) &
!                    -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
!                    +Rv/Rm                                   &
!                   )*Qcond                                   &
!                   +(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
!                    )*(-Qcond)                                   &
!                  )
    CASE('Equiv')
      Cp_eff=Cpd+(RhoV+RhoL)/RhoD*Cpl
      RhsPotM=-PotM/(Cp_eff*RhoD)*Rv*LOG(RelHumidity(T,RhoV))*Qcond
    CASE('Energy')
      RhsPotM=Zero
    CASE('PreEn')
      eLoc=Cvml*T+RhoV*L00
      DpDRhoV=(Rv-Rd)*(eLoc-RhoV*L00)/Cvml &
              -L00*(RhoD*Rd+RhoV*Rv)/Cvml &
              -(Cvv-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
      DpDRhoL=-Rd*(eLoc-RhoV*L00)/Cvml &
              -(Cpl-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
      RhsPotM=(DpDRhoV-DpDRhoL)*Qcond
    END SELECT
    f(1)=Qcond   ! <0, QV
    f(2)=Qcond   ! <0, RHO
    f(3)=RhsPotM ! TH
    f(4)=-QCond  ! >0, QC
  CASE DEFAULT
  IF (MyId==0) THEN
    WRITE(*,*) 'Error in .chem file - wrong process. Stopped.'
    STOP  
  END IF
  END SELECT

END SUBROUTINE MicroProcess

SUBROUTINE LSCMicroJac(Gas,Jac) 
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
!         IF (ix==100.AND.iz==10) THEN
!           WRITE (*,*) 'ix=100, iz=10'
!           WRITE (*,*) 'nReakMicro:',iReak,'von',nReakMicro
!           WRITE (*,*) 'Reaction%NumSpeciesLeftAktiv:',Reaction%NumSpeciesLeftAktiv
!         END IF
          istr=0
          DO k=1,Reaction%NumSpeciesLeftAktiv
!           WRITE (*,*) '(1) Reaction%NumSpeciesLeftAktiv:',k,'von',Reaction%NumSpeciesLeftAktiv
            iSpecies=Reaction%SpeciesLeft(k)
            c(k)=Gas(iSpecies)%c(ix,iy,iz,1)
            IF (iSpecies==RhoPos) THEN
              c(k)=c(k)+RhoProfG(ibLoc)%c(ix,iy,iz,1) 
            END IF
          END DO
          ShiftInd=Reaction%NumSpeciesLeftAktiv
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2)=PreCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+3)=KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4)=zP(iz-1)+zP(iz)
          DO k=1,Reaction%NumSpeciesLeftAktiv
            CALL MicroProcess(c,fC,Reaction%TypeR,'Function')
            Temp=c(k)
            c(k)=c(k)*(1.0d0+1.e-8_RealKind)+SIGN(1.e-8_RealKind,c(k))
            CALL MicroProcess(c,f,Reaction%TypeR,'Jacobi  ')
            DO l=1,Reaction%NumSpecies
              istr=istr+1
              idf=Reaction%Struct(istr)
              Df=(f(l)-fC(l))/(c(k)-Temp+Eps)
              Jac(idf)%c(ix,iy,iz,1)=Jac(idf)%c(ix,iy,iz,1)+Df
!             IF (ix==100.AND.iz==20.AND.ABS(f(l)>0.0d) THEN
!             IF (ABS(f(l))>0.0d0) THEN
!               WRITE(*,*) 'ix iz',ix,iz
!               WRITE(*,*) 'k',k,SpeciesName(Reaction%SpeciesLeft(k))
!               WRITE(*,*) 'l',l,SpeciesName(Reaction%Species(l)%Species)
!               WRITE (*,*) 'k,l,idf,Jac,Df:',k,l,idf,Jac(idf)%c(ix,iy,iz,1),Df
!               WRITE (*,*) 'VolC',VolC(ix,iy,iz),f(l),fc(l)
!             END IF
            END DO
            c(k)=Temp 
          END DO
          Reaction=>Reaction%Next
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE LSCMicroJac

SUBROUTINE LSCMicro(Gas,fGas) 

  TYPE(Vec4_T) :: Gas(:),fGas(:)

  INTEGER :: ix,iy,iz
  INTEGER :: k,iSpecies,iReak
  TYPE(Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: c(20),f(20),Rho
  LOGICAL :: printf

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        Reaction=>MicroFirst
        DO iReak=1,nReakMicro 
          DO k=1,Reaction%NumSpeciesLeftAktiv
            iSpecies=Reaction%SpeciesLeft(k)
            c(k)=Gas(iSpecies)%c(ix,iy,iz,1)
            IF (iSpecies==RhoPos) THEN
              c(k)=c(k)+RhoProfG(ibLoc)%c(ix,iy,iz,1) 
            END IF
          END DO
          ShiftInd=Reaction%NumSpeciesLeftAktiv
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2)=PreCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+3)=KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4)=zP(iz-1)+zP(iz)
!         IF (ix>=97.AND.ix<=103.AND.iz==5) THEN
!           printf=.TRUE.
!         ELSE
!           printf=.FALSE.
!         END IF
          CALL MicroProcess(c,f,Reaction%TypeR,'Function')
          DO k=1,Reaction%NumSpecies
            iSpecies=Reaction%Species(k)%Species
!           IF (Reaction%TypeR/='AUTO') THEN
              fGas(iSpecies)%c(ix,iy,iz,1)=fGas(iSpecies)%c(ix,iy,iz,1)+f(k)
!           END IF
          END DO
          Reaction=>Reaction%Next 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE LSCMicro

END MODULE LSCMicroPhysics_Mod
