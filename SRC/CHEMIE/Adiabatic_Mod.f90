MODULE Adiabatic_Mod

  USE Control_Mod
  USE Parameter_Mod
  USE InputTool_Mod
  USE Rates_Mod
  USE Chemie_Mod
  USE Physics_Mod
  USE Microphysics_Mod
  USE Transport_Mod

  IMPLICIT NONE

  TYPE Trajektory_T
    REAL(RealKind) :: Time
    REAL(RealKind) :: lat
    REAL(RealKind) :: lon
    REAL(RealKind) :: Altitude
    REAL(RealKind) :: Pressure
    REAL(RealKind) :: u,v,w
    REAL(RealKind) :: Te
    REAL(RealKind) :: RhoV
    REAL(RealKind) :: RelHum
    REAL(RealKind) :: Dust
  END TYPE Trajektory_T
  TYPE(Trajektory_T), ALLOCATABLE :: Trajektory(:)
  INTEGER :: nInt,nBerg,nV 
  REAL(RealKind) :: PStart,PEnd
  REAL(RealKind) :: tAbsStart
  REAL(RealKind) :: RelHum
  REAL(RealKind), ALLOCATABLE :: hP(:),sP(:),tp(:),pp(:)
  REAL(RealKind), ALLOCATABLE :: V(:),sV(:)
  REAL(RealKind), ALLOCATABLE :: hSP(:),hEP(:)
  REAL(RealKind), ALLOCATABLE :: vSP(:),vEP(:)
  REAL(RealKind) :: LatLoc,LonLoc,tSimul,Hour
  REAL(RealKind) :: JNO2
  REAL(RealKind), PARAMETER :: Cdil=0.0d0 !3.9d-4 ! Parameter for dilution
  REAL(RealKind), PARAMETER :: RelaxPar=1.0d-1 
  REAL(RealKind) :: p_akt,h_akt,l_akt,i_akt
  REAL(RealKind) :: RhoTot,TVirt,TAbsLoc,RhoV_Loc
  REAL(RealKind) :: TimePrev

  INTEGER :: Date
  NAMELIST /MetData/ LatLoc  &
                  ,LonLoc   &
                  ,Hour &
                  ,Date &
                  ,pStart &
                  ,tAbsStart &
                  ,RelHum &
                  ,pEnd &
                  ,JNO2 
  


CONTAINS 

SUBROUTINE InputAdiabatic(FileName,VecT)

  CHARACTER(*) :: FileName
  TYPE(Vector4Cell_T), POINTER :: VecT(:)

  INTEGER :: i,i1,j,iV
  INTEGER :: InputUnit
  INTEGER :: nTraj
  REAL(RealKind) :: hS,hE
  CHARACTER*300 :: Line
  CHARACTER*40 :: DummyName
  CHARACTER*40 :: ModelName
  REAL(RealKind) :: RhoA,RhoV,esw,ew
  REAL(RealKind) :: KappaLoc
  REAL(RealKind) :: c1,Work(20)
  LOGICAL :: Back
  REAL(RealKind) :: pLoc,TeLoc,RhoVLoc
  REAL(RealKind) :: QLoc,RhoDLoc,QLLoc,RhoThetaLoc
  REAL(RealKind) :: Rm,Cpml
  TYPE(Trajektory_T) :: TempTra
  REAL(RealKind) :: EndTimeTraj

  Date     = 010621       ! Date: yymmdd  (21.June 2001)
  LatLoc      = 4.5e+01      ! latitude  [grad] (Schmuecke)
  LonLoc      = 0.0d+01
  Hour=9.0d0
  JNO2=Zero
  TimePrev = 0.0d0

  InputUnit=1
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&MetData')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=MetData)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(InputUnit)
  tSimul=Hour*3600.0d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    esw=SaturVapor(TAbsStart)
    ew=esw*RelHum
    RhoA=RhoAF(tAbsStart,pStart-ew)
    RhoV=ew/(Rv*tAbsStart)
    RhoV_Loc=RhoV
    RhoTot=RhoA+RhoV
    VecT(ibLoc)%Vec(RhoVPos)%c=RhoV
    KappaLoc=(RhoA*Rd+RhoV*Rv)/(RhoA*Cpd+RhoV*Cpv)
    TVirt=RhoTot*tAbsStart*(p0/pStart)**KappaLoc*(RhoA+Rv/Rd*RhoV)/(RhoA+RhoV)
    VecT(ibLoc)%Vec(ThPos)%c=TVirt
    VecT(ibLoc)%Vec(RhoPos)%c=RhoTot
    TAbsLoc=tAbsStart
    TAbsCell(ibLoc)%Vec(1)%c=tAbsStart
  END DO

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=2) Line
    IF (INDEX(Line,'TOPOGRAPHY')>0) THEN
      EXIT
    END IF
  END DO
  READ(InputUnit,*) nInt
  READ(InputUnit,*) nBerg
  ALLOCATE(hP(nBerg))
  READ(InputUnit,*) hP
  ALLOCATE(sP(nBerg))
  READ(InputUnit,*) sP
  ALLOCATE(pp(nBerg))
  READ(InputUnit,*) pp
  READ(InputUnit,*) nV
  ALLOCATE(sV(nV))
  ALLOCATE(V(nV))
  READ(InputUnit,*) V 
  READ(InputUnit,*) sV 
2 CONTINUE
  CLOSE(InputUnit)

  nInt=MIN(nInt,nBerg)
  ALLOCATE(tp(nInt))
  ALLOCATE(hSP(nBerg-1))
  ALLOCATE(hEP(nBerg-1))
  ALLOCATE(vSP(nBerg-1))
  ALLOCATE(vEP(nBerg-1))
  iV=2
  i1=1
  DO i=1,nBerg 
    IF (sP(i)>=sV(iV)) THEN
      DO j=i1,i-1
        hSP(j)=hP(i1) 
        hEP(j)=hP(i)
        vSP(j)=V(iV-1) 
        vEP(j)=V(iV)
      END DO
      i1=i
      iV=iV+1
    END IF
  END DO
  tp(1)=0.0d0
  DO i=2,nInt
    tp(i)=t1(sP(i))
  END DO
  hS=hP(1)
  hE=hP(1)
  DO i=2,nBerg
    hE=MAX(hE,hP(i))
  END DO
    
  IF (JNO2>Zero) THEN
    chi=Zenith(LatLoc,LonLoc,Date,tSimul)
    Dust=DustFac()
  END IF
  p_akt = pp(1)
  h_akt = hP(1)
  l_akt = sP(1)
  i_akt = 1

!  IF (TrajFile/='') THEN
  IF (TrajIn) THEN
    CALL OpenFile(TrajFile) !OpenFile('Sopran_09.entr')
    CALL LineFile(Back,'BEGIN_MET',End='END_MET' &
                 ,Name1=DummyName,R1=c1)
    nTraj=c1
    i=nTraj
    i=1
    ModelName=TRIM(DummyName)
    ALLOCATE(Trajektory(nTraj))
    DO
      CALL LineFile(Back,'BEGIN_MET',End='END_MET' &
                   ,Name1=DummyName,R=Work)
      IF (Back) THEN
        EXIT
      END IF
      IF (ModelName=='COSMO_MUSCAT') THEN
        Trajektory(i)%Time=Work(1)
        Trajektory(i)%lat=Work(2)
        Trajektory(i)%lon=Work(3)
        Trajektory(i)%Altitude=Work(4)
        Trajektory(i)%Pressure=Work(5)*100.0d0
        Trajektory(i)%u=Work(6)
        Trajektory(i)%v=Work(7)
        Trajektory(i)%w=Work(8)
        Trajektory(i)%Te=Work(9)
        Trajektory(i)%RhoV=Work(10)
        Trajektory(i)%RelHum=Work(11)
        Trajektory(i)%Dust=Work(12)
      ELSE IF (ModelName=='HYSPLIT') THEN
        Trajektory(i)%Time=Work(8)*(-3600.d0)
        Trajektory(i)%lat=Work(9)
        Trajektory(i)%lon=Work(10)
        Trajektory(i)%Altitude=Work(11)
        Trajektory(i)%Pressure=Work(12)*100.0d0
        Trajektory(i)%Te=Work(13)*(Trajektory(i)%Pressure/p0)**(Rd/Cpd)
        Trajektory(i)%RelHum=Work(14)
        Trajektory(i)%RhoV=Work(15)*1.d-3
      END IF
      i=i+1
    END DO
    CALL CloseFile
    EndTimeTraj=Trajektory(nTraj)%Time
    DO i=1,nTraj/2
      TempTra=Trajektory(i)
      Trajektory(i)=Trajektory(nTraj-i+1)
      Trajektory(nTraj-i+1)=TempTra
      Trajektory(i)%Time=EndTimeTraj-Trajektory(i)%Time
      Trajektory(nTraj-i+1)%Time=EndTimeTraj-Trajektory(i)%Time
    END DO  
    pLoc=Trajektory(1)%Pressure
    TeLoc=Trajektory(1)%Te
    RhoVLoc=Trajektory(1)%RhoV
    QLoc=pLoc/(Rd+(Rv-Rd)*RhoVLoc)/TeLoc
    RhoVLoc=QLoc*RhoVLoc
    RhoDLoc=QLoc-RhoVLoc
    QLLoc=0.0d0
    Rm=Rd*RhoDLoc+Rv*RhoVLoc+Eps
    Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*QLLoc+Eps
    RhoThetaLoc=TeLoc*(p0/pLoc)**(Rm/Cpml)*(RhoDLoc+Rv/Rd*RhoVLoc)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      VecT(ibLoc)%Vec(RhoVPos)%c=RhoVLoc
      VecT(ibLoc)%Vec(ThPos)%c=RhoThetaLoc
      VecT(ibLoc)%Vec(RhoPos)%c=QLoc
      TAbsCell(ibLoc)%Vec(1)%c=TeLoc
    END DO
    EndTime=MIN(EndTime,Trajektory(nTraj)%Time)
  END IF
END SUBROUTINE InputAdiabatic

SUBROUTINE UpdraftVelocity(dpdt,Time)

  REAL(RealKind) :: dpdt,Time

  INTEGER :: i,it,nTraj
  REAL(RealKind) :: tL,tU,pL,pU
  REAL(RealKind) :: TimeLoc
  REAL(RealKind) :: h,hL,hU,sL,sU
  REAL(RealKind) :: sFind
  REAL(RealKind) :: t
  REAL(RealKind) :: sNeu
  REAL(RealKind) :: tNeu
  REAL(RealKind) :: v,w,dpdh
  REAL(RealKind), PARAMETER :: TolRel=1.d-4
  REAL(RealKind), PARAMETER :: TolAbs=1.d-8
 
  IF (TrajIn) THEN
    dpdt=0.0d0 ! OSSI
  ELSE
    TimeLoc=MIN(Time,tP(nInt))
    DO i=1,nInt
      IF (TimeLoc<=tP(i+1)) THEN
        it=i
        EXIT
      END IF
    END DO
    tL=tP(it)
    sL=sP(it)
    hL=hP(it)
    pL=pp(it)
    tU=tP(it+1)
    sU=sP(it+1)
    hU=hP(it+1)
    pU=pp(it+1)
    IF (hL==hU) THEN
      dpdh=0.0d0
    ELSE
      dpdh=(pU-pL)/(hU-hL)
    END IF
    IF (TimeLoc-tL<tU-TimeLoc) THEN
      sFind=sL
    ELSE
      sFind=sU
    END IF
    DO
      IF (MIN(TimeLoc-tL,tU-TimeLoc)<=TolRel*TimeLoc+TolAbs) THEN
        EXIT
      ELSE
        sNeu=sL+(sU-sL)/(tU-tL)*(TimeLoc-tL)
        tNeu=t1(sNeu)
        sFind=sNeu
        IF (tNeu<TimeLoc) THEN
          sL=sNeu
          tL=tNeu
        ELSE
          sU=sNeu
          tU=tNeu
        END IF
      END IF
    END DO
    IF (hU==hL) THEN
      w=Zero
    ELSE
      sL=sP(it)
      sU=sP(it+1)
      h=hL+(sFind-sL)/(sU-sL)*(hU-hL)
      v=((h-hSP(it))/(hEP(it)-hSP(it)))**1.0d0*(vEP(it)-vSP(it))+vSP(it)
      w=(hU-hL)/(sU-sL)*v
    END IF
    dpdt=w*dpdh
  END IF

END SUBROUTINE UpdraftVelocity

SUBROUTINE AdiabaticJac(Gas,JacInter,Time)
  TYPE(Vec4_T) :: Gas(0:)
  TYPE(Vec4_T) :: JacInter(:)
  REAL(RealKind) :: Time

  INTEGER :: istr,k,l,idf
  REAL(RealKind) :: Rate,RhoV,tabs,p,RateJac(3)
  REAL(RealKind) :: w,dpdh,dpdt
  REAL(RealKind) :: pLoc,Th
  REAL(RealKind) :: QLoc,QL,RhoD
  REAL(RealKind) :: kappaLoc
  REAL(RealKind) :: RateRho,RateRhoV,RateQL
  INTEGER :: ic,ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Th=Gas(thPos)%c(ix,iy,iz,1)
        tAbs=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
        QLoc=Gas(RhoPos)%c(ix,iy,iz,1)
        RhoV=Gas(RhoVPos)%c(ix,iy,iz,1)
        IF (RhoCPos>0) THEN 
          QL=Gas(RhoCPos)%c(ix,iy,iz,1)
        ELSE
          QL=RhoLCell(ibLoc)%c(ix,iy,iz,1)
        END IF
        RhoD=QLoc-RhoV-QL
        pLoc=PreCell(ibLoc)%c(ix,iy,iz,1)
        kappaLoc=(RhoD*Rd+RhoV*Rv)/(RhoD*Cpd+RhoV*Cpv)
        CALL UpdraftVelocity(dpdt,Time)
!       RateRho=(1.d0-kappaLoc)*QLoc/pLoc*dpdt 
!       RateRhoV=(1.d0-kappaLoc)*RhoV/pLoc*dpdt 
!       fGas(thPos)%c(ix,iy,iz,1)=fGas(thPos)%c(ix,iy,iz,1)+Rate
! Derivative with respect to RhoV
!       RateJac(1)=-0.2d0*rd*tAbs*dpdt/(Cpd*p*(1.0d0+0.8d0*RhoV)*(1.0d0+0.8d0*RhoV))
! Derivative with respect to tAbs
!       RateJac(2)=rd*(1.0d0+0.6d0*RhoV)/(Cpd*(1.0d0+0.8d0*RhoV)*p)*dpdt
! Derivative with respect to tAbs
!       RateJac(3)=-rd*tAbs*(1.0d0+0.6d0*RhoV)/(Cpd*(1.0d0+0.8d0*RhoV)*p*p)*dpdt
!       istr=0
!       DO k=1,AdiabaticFirst%NumSpeciesLeftAktiv
!         DO l=1,AdiabaticFirst%NumSpecies
!           istr=istr+1
!           idf=AdiabaticFirst%Struct(istr)
!           Call AddScalar(JacInter(idf),JacInter(idf) &
!                         ,RateJac(istr)*AdiabaticFirst%Species(l)%Koeff,1)
!         END DO
!       END DO
      END DO
    END DO
  END DO
END SUBROUTINE AdiabaticJac

SUBROUTINE Adiabatic(Gas,fGas,Time)

  TYPE(Vec4_T) :: Gas(0:),fGas(0:)
  REAL(RealKind) :: Time

  REAL(RealKind) :: tabs, pLoc, Th
  REAL(RealKind) :: QLoc, RhoV, QL, RhoD
  REAL(RealKind) :: kappaLoc
  REAL(RealKind) :: Rate, RateRho, RateRhoV, RateQL
  REAL(RealKind) :: w,dpdh,dpdt,RhoCs2
  INTEGER :: ic,ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Th=Gas(thPos)%c(ix,iy,iz,1)
        tAbs=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
        QLoc=Gas(RhoPos)%c(ix,iy,iz,1)
        RhoV=Gas(RhoVPos)%c(ix,iy,iz,1)
        IF (RhoCPos>0) THEN 
          QL=Gas(RhoCPos)%c(ix,iy,iz,1)
        ELSE
          QL=RhoLCell(ibLoc)%c(ix,iy,iz,1)
        END IF
        RhoD=QLoc-RhoV-QL
        pLoc=PreCell(ibLoc)%c(ix,iy,iz,1)
        kappaLoc=(RhoD*Rd+RhoV*Rv)/(RhoD*Cpd+RhoV*Cpv)
        CALL UpdraftVelocity(dpdt,Time)
        RhoCs2=SoundCell(ibLoc)%c(ix,iy,iz,1)
        DO ic=1,UBOUND(fGas,1)
          IF (ic/=iRelax) THEN
            fGas(ic)%c(ix,iy,iz,1)=fGas(ic)%c(ix,iy,iz,1)+Gas(ic)%c(ix,iy,iz,1)/RhoCs2*dpdt
          END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE Adiabatic

SUBROUTINE MetAmbientComput(c,Time)
  TYPE(Vec4_T) :: c(0:)
  REAL(RealKind) :: Time

  REAL(RealKind) :: pLoc,TeLoc
  REAL(RealKind) :: TimeLoc
  REAL(RealKind) :: RhoLoc,RhoDLoc,RhoLLoc,RhoVLoc,RhoThetaLoc
  REAL(RealKind) :: Rm,Cpml
  REAL(RealKind) :: tL,tU
  REAL(RealKind) :: pL,pU
  REAL(RealKind) :: TeL,TeU
  REAL(RealKind) :: RhoVL,RhoVU
  INTEGER :: ix,iy,iz
  INTEGER :: nTraj,i,it
  REAL(RealKind) :: FVIRT,Cpm,KAPPA_m
  REAL(RealKind) :: DELTAT,DELTABS
  REAL(RealKind) :: DRhoV,DRHO,DELT,RhoA
  REAL(RealKind) :: v_akt,w_akt
  REAL(RealKind) :: CurP,CurRhoV
  REAL(RealKind) :: esw,ew
  REAL(RealKind) :: DP,dpdh,dpdt

  nTraj=SIZE(Trajektory)
  TimeLoc=MIN(Time,Trajektory(nTraj)%Time)
  it=nTraj-1
  DO i=1,nTraj-1
    IF (TimeLoc<=Trajektory(i+1)%Time) THEN
      it=i
      EXIT
    END IF
  END DO
  tL=Trajektory(it)%Time
  tU=Trajektory(it+1)%Time
  pL=Trajektory(it)%Pressure
  pU=Trajektory(it+1)%Pressure
  pLoc=pL+(pU-pL)*(Time-tL)/(tU-tL)
  TeL=Trajektory(it)%Te
  TeU=Trajektory(it+1)%Te
  TeLoc=TeL+(TeU-TeL)*(Time-tL)/(tU-tL)
  RhoVL=Trajektory(it)%RhoV
  RhoVU=Trajektory(it+1)%RhoV
  RhoVLoc=RhoVL+(RhoVU-RhoVL)*(Time-tL)/(tU-tL)
  RhoLoc=pLoc/(Rd+(Rv-Rd)*RhoVLoc)/TeLoc
  RhoVLoc=RhoLoc*RhoVLoc
  RhoDLoc=RhoLoc-RhoVLoc
  RhoLLoc=0.0d0
  Rm=Rd*RhoDLoc+Rv*RhoVLoc+Eps
  Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Eps
  RhoThetaLoc=TeLoc*(p0/pLoc)**(Rm/Cpml)*(RhoDLoc+Rv/Rd*RhoVLoc)

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
          VecAmb(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1)=RhoLoc
          VecAmb(ibLoc)%Vec(RhoVPos)%c(ix,iy,iz,1)=RhoVLoc
          VecAmb(ibLoc)%Vec(ThPos)%c(ix,iy,iz,1)=RhoThetaLoc
      END DO
    END DO
  END DO

END SUBROUTINE MetAmbientComput

SUBROUTINE DilutionParcel(Gas,fGas,Time)
  REAL(RealKind) :: Dilut
  TYPE(Vec4_T) :: Gas(0:),fGas(0:)
  REAL(RealKind) :: Time

  INTEGER :: Pos,k
  INTEGER :: ic,ix,iy,iz

  IF (TIME>0.0d0) THEN
    Dilut=Cdil/SQRT(Time)
  ELSE
    Dilut=0.0d0
  END IF

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        fGas(1)%c(ix,iy,iz,1)=fGas(1)%c(ix,iy,iz,1) &
                              -Dilut*(Gas(1)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(1)%c(ix,iy,iz,1)) 
        DO k=3,nAqua    ! k=1 aNUMBER, k=2 aRELAX, k>2 aqu species
          Pos=k
          fGas(Pos)%c(ix,iy,iz,1)=fGas(Pos)%c(ix,iy,iz,1) &
                                  -Dilut*(Gas(Pos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)) 
        END DO
        DO k=nAqua+NumMet+1,nAqua+nGas
          Pos=k
          fGas(Pos)%c(ix,iy,iz,1)=fGas(Pos)%c(ix,iy,iz,1) &
                                  -Dilut*(Gas(Pos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)) 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE DilutionParcel

SUBROUTINE DilutionParcelJac(Gas,Jac,Jac1,Time)
  REAL(RealKind) :: Dilut,Dilutdt
  TYPE(Vec4_T) :: Gas(0:)
  TYPE(Vec4_T) :: Jac(:)
  TYPE(SpMatrix4Cell_T) :: Jac1
  REAL(RealKind) :: Time

  INTEGER :: diag
  INTEGER :: Pos,k
  INTEGER :: ic,ix,iy,iz
  INTEGER, POINTER :: DiagP(:),Permu(:)

  DiagP=>Jac1%Struct%DiagPtr(:)
  Permu=>Jac1%Struct%Permu(:)
  IF (TIME>0.0d0) THEN
    Dilut=Cdil/SQRT(Time+Eps)
  ELSE
    Dilut=0.0d0
  END IF

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        diag=DiagP(Permu(1))
        Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                               -Dilut 
        DO k=3,nAqua
          Pos=k
          diag=DiagP(Permu(Pos))
          Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                                 -Dilut 
        END DO
        DO k=nAqua+NumMet+1,nAqua+nGas
          Pos=k
          diag=DiagP(Permu(Pos))
          Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                                 -Dilut 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE DilutionParcelJac

SUBROUTINE RelaxationParcel(Gas,fGas,Time)
  REAL(RealKind) :: Dilut
  TYPE(Vec4_T) :: Gas(0:),fGas(0:)
  REAL(RealKind) :: Time

  INTEGER :: Pos,k
  INTEGER :: ic,ix,iy,iz

  IF (TIME>0.0d0) THEN
    Dilut=Cdil/SQRT(Time)
  ELSE
    Dilut=0.0d0
  END IF

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        fGas(1)%c(ix,iy,iz,1)=fGas(1)%c(ix,iy,iz,1) &
                              -RelaxPar*Gas(1)%c(ix,iy,iz,1)/Gas(RhoPos)%c(ix,iy,iz,1)  &
                              *(Gas(RhoPos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1))
        DO k=3,nAqua    ! k=1 aNUMBER, k=2 aRELAX, k>2 aqu species
          Pos=k
          fGas(Pos)%c(ix,iy,iz,1)=fGas(Pos)%c(ix,iy,iz,1) &
                                  -RelaxPar*Gas(Pos)%c(ix,iy,iz,1)/Gas(RhoPos)%c(ix,iy,iz,1)  &
                                  *(Gas(RhoPos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1))
        END DO
        DO k=nAqua+1,nAqua+NumMet
          Pos=k
          fGas(Pos)%c(ix,iy,iz,1)=fGas(Pos)%c(ix,iy,iz,1) &
                                  -RelaxPar*(Gas(Pos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(Pos)%c(ix,iy,iz,1))
        END DO
        DO k=nAqua+NumMet+1,nAqua+nGas
          Pos=k
          fGas(Pos)%c(ix,iy,iz,1)=fGas(Pos)%c(ix,iy,iz,1) &
                                  -RelaxPar*Gas(Pos)%c(ix,iy,iz,1)/Gas(RhoPos)%c(ix,iy,iz,1)  &
                                  *(Gas(RhoPos)%c(ix,iy,iz,1)-VecAmb(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1))
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE RelaxationParcel

SUBROUTINE RelaxationParcelJac(Gas,Jac,Jac1,Time)
  REAL(RealKind) :: Dilut,Dilutdt
  TYPE(Vec4_T) :: Gas(0:)
  TYPE(Vec4_T) :: Jac(:)
  TYPE(SpMatrix4Cell_T) :: Jac1
  REAL(RealKind) :: Time

  INTEGER :: diag
  INTEGER :: Pos,k
  INTEGER :: ic,ix,iy,iz
  INTEGER, POINTER :: DiagP(:),Permu(:)

  DiagP=>Jac1%Struct%DiagPtr(:)
  Permu=>Jac1%Struct%Permu(:)
  IF (TIME>0.0d0) THEN
    Dilut=Cdil/SQRT(Time+Eps)
  ELSE
    Dilut=0.0d0
  END IF

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        diag=DiagP(Permu(1))
        Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                               -Dilut 
        DO k=3,nAqua
          Pos=k
          diag=DiagP(Permu(Pos))
          Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                                 -Dilut 
        END DO
        DO k=nAqua+NumMet+1,nAqua+nGas
          Pos=k
          diag=DiagP(Permu(Pos))
          Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1) &
                                 -Dilut 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE RelaxationParcelJac

FUNCTION t1(s)
  REAL(RealKind) :: t1
  REAL(RealKind) :: s

  INTEGER :: i,ih
  REAL(RealKind) :: a,h,lu,lo

  DO i=1,nInt
    IF (s<=sP(i+1)) THEN
      ih=i
      EXIT
    END IF
  END DO
  IF (vEP(ih)/=vSP(ih)) THEN
    h=(s-sP(ih))/(sP(ih+1)-sP(ih))*(hP(ih+1)-hP(ih))+hP(ih)
    lu=(hP(ih)-hSP(ih))/(hEP(ih)-hSP(ih))
    lo=(h-hSP(ih))/(hEP(ih)-hSP(ih))
    a=vSP(ih)/(vEP(ih)-vSP(ih))
    t1=((hEP(ih)-hSP(ih))/(vEP(ih)-vSP(ih)))*((sP(ih+1)-sP(ih))/(hP(ih+1)-hP(ih)))            &
      *LOG((lo+a)/(lu+a))+tP(ih)
  ELSE
    IF (hP(ih+1)==hP(ih)) THEN
      t1=(s-sP(ih))/vSP(ih)+tP(ih)
    ELSE
      h=(s-sP(ih))/(sP(ih+1)-sP(ih))*(hP(ih+1)-hP(ih))+hP(ih)
      lu=(hP(ih)-hSP(ih))/(hEP(ih)-hSP(ih))
      lo=(h-hSP(ih))/(hEP(ih)-hSP(ih))
      t1=(hEP(ih)-hSP(ih))*((sP(ih+1)-sP(ih))/(hP(ih+1)-hP(ih)))*(lo-lu)/vSP(ih)+tP(ih)
    END IF
  END IF
    
END FUNCTION t1

FUNCTION DustFac()

!=========================================================
!===  Scaling of NO2 photolysis rate JNO2 by measurement         
!===  PHOTABC:  A: 7.67e-03   B: 1.773179e-00   C:  0.77233e-00
!=========================================================
!

  REAL(RealKind) :: DustFac
  REAL(RealKind) ::  chiz,ychiz,eychiz
  REAL(RealKind) ::  JNO2_Rate

  chiz=chi*0.77233d-00
  IF (chiz<PiHalf) THEN
    ychiz=1.773179d-00*(One-(One/COS(chiz)))
    IF (ychiz>-30.0d0) THEN
      eychiz=EXP(ychiz)
    ELSE
      eychiz=9.357d-14
    END IF
  ELSE
    eychiz=9.357d-14
  END IF
  JNO2_Rate=7.67d-03*eychiz

  DustFac=JNO2/JNO2_Rate
END FUNCTION DustFac


        REAL(RealKind) FUNCTION  Zenith(LAT,LONG,IDAT,Time)
!  this subroutine calculates solar zenith and azimuth angles for a particular
!  time and location.  Must specify:
!  INPUT:
!       LAT  - latitude in decimal degrees
!       LONG - longitude in decimal degrees
!       IDAT - Date at Greenwich - specify year (19yy), month (mm), day (dd)
!              format is six-digit integer:  yymmdd
!       Time   c$simulation time in seconds
!       GMT  - Greenwich mean time - decimal military eg.
!               22.75 = 45 min after ten pm gmt
!  OUTPUT
!       Zenith
!       Azimuth
!*******************************************************************************
!
! subroutine received from G. MAUERSBERGER 1/98 REPLACES sonne() subroutine in 
! surf program
!
!*******************************************************************************
!       USE mo_control

        IMPLICIT REAL(RealKind) (A-H,O-Z)
        IMPLICIT INTEGER (I-N)
       
        INTEGER :: IDAT
        REAL(RealKind) :: LAT,LONG
        REAL(RealKind) :: LBGMT,LZGMT
        REAL(RealKind) :: ML
!
        INTEGER :: IMN(12)
        DATA IMN/31,28,31,30,31,30,31,31,30,31,30,31/

!--------------------------------------------------------
!  set GMT
        GMT = Time / 3600.e0
!  convert to radians
        DR=Pi/180.0d0
        RLT = LAT*DR
        RPHI = LONG*DR
!  parse date
        IIYEAR = IDAT/10000
        IYEAR = 19*100 + IIYEAR
        IF (IIYEAR <= 50) IYEAR = IYEAR + 100 
        IMTH = (IDAT - IIYEAR*10000)/100
        IDAY = IDAT - IIYEAR*10000 - IMTH*100
!  identify and correct leap years
        IIY = (IIYEAR/4)*4
        IF(IIY.EQ.IIYEAR) IMN(2) = 29
!  count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5
        YREF =  2442047.5
        NYEARS = IYEAR - 1974
        LEAP = (NYEARS+1)/4
        IF(NYEARS.LE.-1) LEAP = (NYEARS-2)/4
        NOLEAP = NYEARS - LEAP
        YR = YREF + 365.*NOLEAP + 366.*LEAP
!
        IJD = 0
        IN = IMTH - 1
        IF(IN.EQ.0) GO TO 40
        DO 30 I=1,IN
        IJD = IJD + IMN(I)
   30   CONTINUE
        IJD = IJD + IDAY
        GO TO 50
   40   IJD = IDAY
   50   IJ = IYEAR - 1973
!      print julian days current "ijd"
        JD = IJD + (YR - YREF)
        D = JD + GMT/24.0
!      calc geom mean longitude
        ML = 279.2801988 + .9856473354*D + 2.267E-13*D*D
        RML = ML*DR
!
!      calc equation of time in sec
!      w = mean long of perigee
!      e = eccentricity
!      epsi = mean obliquity of ecliptic
        W = 282.4932328 + 4.70684E-5*D + 3.39E-13*D*D
        WR = W*DR
        EC = 1.6720041E-2 - 1.1444E-9*D - 9.4E-17*D*D
        EPSI = 23.44266511 - 3.5626E-7*D - 1.23E-15*D*D
        PEPSI = EPSI*DR
        YT = (TAN(PEPSI/2.0))**2
        CW = COS(WR)
        SW = SIN(WR)
        SSW = SIN(2.0*WR)
        EYT = 2.*EC*YT
        FEQT1 = SIN(RML)*(-EYT*CW - 2.*EC*CW)
        FEQT2 = COS(RML)*(2.*EC*SW - EYT*SW)
        FEQT3 = SIN(2.*RML)*(YT - (5.*EC**2/4.)*(CW**2-SW**2))
        FEQT4 = COS(2.*RML)*(5.*EC**2*SSW/4.)
        FEQT5 = SIN(3.*RML)*(EYT*CW)
        FEQT6 = COS(3.*RML)*(-EYT*SW)
        FEQT7 = -SIN(4.*RML)*(.5*YT**2)
        FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7
        EQT = FEQT*13751.0
!
!   convert eq of time from sec to deg
        REQT = EQT/240.
!   calc right ascension in rads
        RA = ML - REQT
        RRA = RA*DR
!   calc declination in rads, deg
        TAB = 0.43360*SIN(RRA)
        RDECL = ATAN(TAB)
        DECL = RDECL/DR
!   calc local hour angle
        LBGMT = 12.0 - EQT/3600. + LONG*24./360.
        LZGMT = 15.0*(GMT - LBGMT)
        ZPT = LZGMT*DR
        CSZ = SIN(RLT)*SIN(RDECL) + COS(RLT)*COS(RDECL)*COS(ZPT)
        ZR = ACOS(CSZ)
!   calc local solar azimuth
        CAZ = (SIN(RDECL) - SIN(RLT)*COS(ZR))/(COS(RLT)*SIN(ZR))
        RAZ = ACOS(CAZ)
        AZIMUTH = RAZ/DR

!--- set Zenith Angle
        Zenith =  1.745329252D-02 * ZR/DR
!       WRITE(*,*) 'Zenith',Zenith*180.d0/Pi
        
!---------------------------------------------------------------
200     FORMAT(' ',F7.2,2(12X,F7.2))

        END FUNCTION  Zenith
END MODULE Adiabatic_Mod
