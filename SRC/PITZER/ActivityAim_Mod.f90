MODULE ActivityAim_Mod

  USE Control_Mod
  USE Chemie_Mod
  USE Aerosol_Mod
  USE Microphysics_Mod

  IMPLICIT NONE

  INTEGER, ALLOCATABLE :: AimInd(:)
  INTEGER :: nN,nA,nC
  INTEGER, ALLOCATABLE :: ChargeCat(:),ChargeAn(:) 
  REAL(RealKInd), ALLOCATABLE :: Vca(:,:) &
                                ,Vac(:,:)
  CHARACTER(20), ALLOCATABLE :: nameA(:),nameC(:),nameN(:)
  TYPE IntVec_T
    INTEGER, POINTER :: Ind(:)=>NULL()
  END TYPE IntVec_T
  TYPE(IntVec_T), ALLOCATABLE :: IndCAT(:)
  TYPE(IntVec_T), ALLOCATABLE :: IndAN(:)
  INTEGER :: zMax

!     Arrays holding stoichiometry ratios
!     (Vca = za/(zc+za), Vac = zc/(zc+za)) for the ions used
!     in the mole fraction model.
  REAL(RealKind), PARAMETER :: RhoAim=13.0d0
  REAL(RealKind), ALLOCATABLE :: A1ca(:,:) &
                                 ,A2ca(:,:) &
                                 ,B1ca(:,:) &
                                 ,B2ca(:,:) &
                                 ,Wnca(:,:,:) &
                                 ,Unca(:,:,:) &
                                 ,Vnca(:,:,:) &
                                 ,Wcca(:,:,:) &
                                 ,Waac(:,:,:) &
                                 ,Qncca(:,:,:,:) &
                                 ,Qnaac(:,:,:,:) &
                                 ,Ucca(:,:,:) &
                                 ,Uaac(:,:,:) &
                                 ,Wnn(:,:) &
                                 ,Unn(:,:) &
                                 ,Ynnca(:,:,:,:) &
                                 ,Xaaac(:,:,:,:) &
                                 ,Xccca(:,:,:,:) &
                                 ,Zccaa(:,:,:,:) &
                                 ,Cnnn(:,:,:)



CONTAINS

SUBROUTINE InitSystem

  INTEGER :: iA,iC,IndSpc

!
!  Ref: Cations           Ref: Anions              Ref: Neutrals
!    1  H                   1  HSO4                  1  H2O
!    2  NH4                 2  SO4
!    3  Na                  3  NO3
!                           4  Cl
!                           5  Br

  nA=5
  ALLOCATE(nameA(nA))
  nameA(1)='HSO4m'
  nameA(2)='SO4mm'
  nameA(3)='NO3m'
  nameA(4)='CLm'
  nameA(5)='BRm'
  ALLOCATE(ChargeAn(nA))
  ChargeAn(1)=1
  ChargeAn(2)=2
  ChargeAn(3)=1
  ChargeAn(4)=1
  ChargeAn(5)=1

  nC=3
  ALLOCATE(nameC(nC))
  nameC(1)='Hp'
  nameC(2)='NH4p'
  nameC(3)='NAp'
  ALLOCATE(ChargeCat(nC))
  ChargeCat(1)=1
  ChargeCat(2)=1
  ChargeCat(3)=1

  nN=1
  ALLOCATE(nameN(nN))
  nameN(1)='aH2O'

  DO iC=1,nC
    IndSpc=Position(nameC(iC))
    IF (IndSpc>0) THEN
      AimInd(IndSpc)=iC
    END IF
  END DO
  DO iA=1,nA
    IndSpc=Position(nameA(iA))
    IF (IndSpc>0) THEN
      AimInd(IndSpc)=-iA
    END IF
  END DO

  CALL SetVCA
  CALL SetUnsymm

END SUBROUTINE InitSystem

SUBROUTINE SetVCA

  INTEGER :: iA,iC

  ALLOCATE(Vca(nC,nA))
  ALLOCATE(Vac(nA,nC))
  DO iA=1,nA
    DO iC=1,nC
      Vca(iC,iA)=(ChargeCat(iC)+ChargeAn(iA))/FLOAT(ChargeAn(iA))
      Vac(iA,iC)=(ChargeCat(iC)+ChargeAn(iA))/FLOAT(ChargeCat(iC))
    END DO
  END DO

END SUBROUTINE SetVCA

SUBROUTINE SetUnsymm

  INTEGER :: iA,iC,iz,izNum
  INTEGER :: zMaxLoc

  zMaxLoc=0
  DO iA=1,nA
    zMaxLoc=MAX(zMaxLoc,ChargeAn(iA))
  END DO
  ALLOCATE(IndAN(zMaxLoc))
  DO iz=1,zMaxLoc
    izNum=0
    DO iA=1,nA
      IF (ChargeAn(iA)==iz) THEN
        izNum=izNum+1
      END IF
    END DO
    ALLOCATE(IndAN(iz)%Ind(izNum))
    izNum=0
    DO iA=1,nA
      IF (ChargeAn(iA)==iz) THEN
        izNum=izNum+1
        IndAN(iz)%Ind(izNum)=iA
      END IF
    END DO
  END DO
  zMax=zMaxLoc

  zMaxLoc=0
  DO iC=1,nC
    zMaxLoc=MAX(zMaxLoc,ChargeCat(iC))
  END DO
  ALLOCATE(IndCAT(zMaxLoc))
  DO iz=1,zMaxLoc
    izNum=0
    DO iC=1,nC
      IF (ChargeCat(iC)==iz) THEN
        izNum=izNum+1
      END IF
    END DO
    ALLOCATE(IndCAT(iz)%Ind(izNum))
    izNum=0
    DO iC=1,nC
      IF (ChargeCat(iC)==iz) THEN
        izNum=izNum+1
        IndCAT(iz)%Ind(izNum)=iC
      END IF
    END DO
  END DO
  zMax=MAX(zMax,zMaxLoc)
END SUBROUTINE SetUnsymm

SUBROUTINE InitActivityAim
  
  ALLOCATE(AimInd(nAqua))
  AimInd=0
  CALL InitSystem
  CALL AllocateCoefficient 

END SUBROUTINE InitActivityAim

FUNCTION FUNCI(xCAT,xAN)
  REAL(RealKind) :: FUNCI
  REAL(RealKind) :: xCAT(:),xAN(:)

  INTEGER :: i

  FUNCI=0.D0
  DO i=1,nC
    FUNCI=FUNCI+xCAT(i)*ChargeCat(i)**2
  END DO
  DO i=1,nA
    FUNCI=FUNCI+xAN(i)*ChargeAn(i)**2
  END DO
  FUNCI = FUNCI/2.D0

END FUNCTION FUNCI


SUBROUTINE AllocateCoefficient

  ALLOCATE(A1ca(nC,nA))
  ALLOCATE(A2ca(nC,nA))
  ALLOCATE(B1ca(nC,nA))
  ALLOCATE(B2ca(nC,nA))
  ALLOCATE(Wnca(nN,nC,nA))
  ALLOCATE(Unca(nN,nC,nA))
  ALLOCATE(Vnca(nN,nC,nA))
  ALLOCATE(Wcca(nC,nC,nA))
  ALLOCATE(Waac(nA,nA,nC))
  ALLOCATE(Qncca(nN,nC,nC,nA))
  ALLOCATE(Qnaac(nN,nA,nA,nC))
  ALLOCATE(Ucca(nC,nC,nA))
  ALLOCATE(Uaac(nA,nA,nC))
  ALLOCATE(Wnn(nN,nN))
  ALLOCATE(Unn(nN,nN))
  ALLOCATE(Ynnca(nN,nN,nC,nA))
  ALLOCATE(Xaaac(nA,nA,nA,nC))
  ALLOCATE(Xccca(nC,nC,nC,nA))
  ALLOCATE(Zccaa(nC,nC,nA,nA))
  ALLOCATE(Cnnn(nN,nN,nN))

  A1ca=Zero
  A2ca=Zero
  B1ca=Zero
  B2ca=Zero
  Wnca=Zero
  Unca=Zero
  Vnca=Zero
  Wcca=Zero
  Waac=Zero
  Qncca=Zero
  Qnaac=Zero
  Ucca=Zero
  Uaac=Zero
  Wnn=Zero
  Unn=Zero
  Ynnca=Zero
  Xaaac=Zero
  Xccca=Zero
  Zccaa=Zero
  Cnnn=Zero

END SUBROUTINE AllocateCoefficient


FUNCTION AFT(T)

  REAL(RealKInd) :: AFT
  REAL(RealKInd) :: T
!
!    ============================================================
!   |  Function yields mf based DH coefficient for absolute      |
!   |  temperature T (K). Values below 273.15 K are extrapolated |
!   |                                                            |
!   |          Tr must be >= 234.15K,                            |
!   |          Aphi   = DH constant at Tr,                       |
!   |          AlRTr  = AL / (RTr) at Tr,                        |
!   |          AjR    = AJ / R at Tr,                            |
!   |          dAjRdT = d (AJ/R) / dT at Tr.                     |
!    ============================================================
!
  REAL(RealKind), PARAMETER :: R=8.3144D0  &
                              ,Tr=273.15D0  &
                              ,AphiTr=0.376421485D0  &
                              ,AlRTr=0.600305325D0 &
                              ,AjR=1.677818299D0  &
                              ,dAjRdT=0.1753875763D0 &
                              ,XMIN=234.15D0  &
                              ,XMAX=373.15D0
  INTEGER, PARAMETER :: N=17

  REAL(RealKind) :: AA(0:18)= &
      (/0.797256081240D+00,0.573389669896D-01, &
       0.977632177788D-03, 0.489973732417D-02,-0.313151784342D-02, &
       0.179145971002D-02,-0.920584241844D-03, 0.443862726879D-03, &
      -0.203661129991D-03, 0.900924147948D-04,-0.388189392385D-04, &
       0.164245088592D-04,-0.686031972567D-05, 0.283455806377D-05, &
      -0.115641433004D-05, 0.461489672579D-06,-0.177069754948D-06, &
       0.612464488231D-07,-0.175689013085D-07/)  

  INTEGER :: i
  REAL(RealKind) :: a,x
  REAL(RealKind) :: AlTr,AjTr,dAJdT
  REAL(RealKind) :: TX(0:18)
!
!
!   ..for T>Tr, polynomial reproduces Archer's values:
  IF (T>=Tr) THEN
    X = (2.D0*T - XMAX - XMIN)/(XMAX - XMIN)
    TX(0) = 1.D0
    TX(1) = X
    AFT   = AA(0)*TX(0)/2.D0 + AA(1)*TX(1)
    DO I=1,N
      TX(I+1) = 2*X*TX(I) - TX(I-1)
      AFT     = AFT + AA(I+1)*TX(I+1)
    END DO
!
!   ..for T<Tr, we have an empirical extrapolation:
  ELSE
    AlTr  = AlRTr * R * Tr
    AjTr  = AjR * R
    dAJdT = dAjRdT * R 
    a     = 1.45824246467D0
    x     = (dAJdT - a)/(2.D0*Tr)
!
    AFT = AphiTr + AlTr/(4.D0*R)*(1.D0/Tr-1.D0/T)          &
        + AjTr/(4.D0*R)*(LOG(T/Tr)+Tr/T-1.D0)              &
        + a/(8.D0*R)*(T-Tr**2/T-2.D0*Tr*LOG(T/Tr))         &
        + x/(4.D0*R)*(T**2/6.D0+Tr**2/2.D0-Tr**2*LOG(T/Tr) &
                      -2.D0/3.D0*Tr**3/T)
!
  ENDIF
!
!   ..convert to mole fraction basis:
  AFT = AFT * 7.4504148D0
!
END FUNCTION AFT

FUNCTION gFunc(x)

  REAL(RealKind) :: gFunc
  REAL(RealKind) :: x

  gFunc=2.0d0*(1.0d0-(1.0d0+x)*EXP(-x))/(x*x)

END FUNCTION gFunc

SUBROUTINE JFunc(x,J,JPrime)

  REAL(RealKind) :: x
  REAL(RealKind) :: J,JPrime

  REAL(RealKind) :: T1,T2,T3,Z

  REAL(RealKind), PARAMETER :: C1=4.581D0
  REAL(RealKind), PARAMETER :: C2=-7.238D-1
  REAL(RealKind), PARAMETER :: C3=-1.2D-2
  REAL(RealKind), PARAMETER :: C4=5.28D-1
 
  T1=C1*x**C2
  T2=C3*x**C4
  T3=EXP(T2)
  Z=4.0d0+T1*T3
  J=x/Z
  JPrime=(Z-T3*(C2*T1+C4*T1*T2))/(Z*Z)

END SUBROUTINE JFunc


SUBROUTINE ExcessGibbs(xCAT,xAN,xNEUT,T &
                      ,ActCat,ActAn,ActNeut &
                      ,Excess)

  REAL(RealKind) :: xCAT(:),xAN(:),xNEUT(:)
  REAL(RealKind) :: T
  REAL(RealKind) :: ActCat(:),ActAn(:),ActNeut(:)
  REAL(RealKind) :: Excess

  INTEGER :: iA,i1A,i2A,i3A
  INTEGER :: iC,i1C,i2C,i3C
  INTEGER :: iN,i1N,i2N,i3N
  INTEGER :: zi,zj
  REAL(RealKind) :: AxCalc,IxCalc,SQRTIxCalc
  REAL(RealKind) :: Sum,DSumDIx,SumCatZ,SumAnZ
  REAL(RealKind) :: xij
  REAL(RealKind) :: JCalc(zMax,zMax),JCalcDx(zMax,zMax)
  REAL(RealKind) :: ThetaCalc(zMax,zMax),ThetaCalcDIx(zMax,zMax)
  REAL(RealKind) :: SumAC(zMax)
  REAL(RealKind) :: Temp1,Temp2
  REAL(RealKind) :: Temp3,Temp4
  REAL(RealKind) :: DH1,DH2,DDHDIx
  REAL(RealKind) :: HOE,DHOEDIx
  REAL(RealKind) :: xBar,uBar,uPrimeBar
  REAL(RealKind) :: gCalc1,gCalc2
  REAL(RealKind) :: DgCalc1DIx,DgCalc2DIx
  REAL(RealKind) :: ActCatLoc(nC),ActAnLoc(nA),ActNeutLoc(nN)
  REAL(RealKind) :: ActCatLoc1(nC),ActAnLoc1(nA),ActNeutLoc1(nN)
  REAL(RealKind) :: ActCatLoc2(nC),ActAnLoc2(nA),ActNeutLoc2(nN)
  REAL(RealKind) :: ActCatLoc3(nC),ActAnLoc3(nA),ActNeutLoc3(nN)
  REAL(RealKind) :: SumWcca,SumWaac
  REAL(RealKind) :: SumUcca,SumUaac
  REAL(RealKind) :: SumWnca,SumUnca,SumVnca
  REAL(RealKind) :: SumQncca,SumQnaac
  REAL(RealKind) :: SumYnnca
  REAL(RealKind) :: SumXccca,SumXaaac
  REAL(RealKind) :: SumZccaa
  REAL(RealKind) :: SumWnn,SumUnn
  REAL(RealKind) :: Sumnn,SumCnnn
  REAL(RealKind) :: CorrAct


  AxCalc=AFT(T)
  IxCalc=FUNCI(xCAT,xAN)
  SQRTIxCalc = SQRT(IxCalc)

!     --------------------------------------
!    | (1) Charge * mole fraction sums that |
!    | we will need later.                  |
!     --------------------------------------
  SumCatZ=0.D0
  DO iC=1,nC
    SumCatZ=SumCatZ+xCat(iC)*ChargeCat(iC)
  END DO
  SumAnZ=0.D0
  DO iA=1,nA
    SumAnZ=SumAnZ+xAn(iA)*ChargeAn(iA)
  END DO
  ActCat=0.0d0
  ActAn=0.0d0
  ActNeut=0.0d0
  

!
!   #(01)# Calculate Debye-Huckel (long-range) terms.
!          ** Act coeff contribution is +SUMDH ** 
!
  Temp1=1.0d0+RhoAim*SQRTIxCalc
  Temp2=LOG(Temp1)
  DH1=-4.0d0*(AxCalc/RhoAim)*IxCalc*Temp2   
  DDHDIx=-4.0d0*(AxCalc/RhoAim)*Temp2 &
           -2.0d0*AxCalc*SQRTIxCalc/Temp1

  Temp1=-DDHDIx*(IxCalc-0.5d0)+DH1
  ActCatLoc=0.0d0
  DH2=0.0d0
  DO iA=1,nA
    Sum=0.0d0
    DSumDIx=0.0d0
    DO iC=1,nC
      IF (A1ca(iC,iA)>0.0d0) THEN
        xBar=A1ca(iC,iA)*SQRTIxCalc
        Temp1=EXP(-xBar)
        uBar=2.0d0*B1ca(iC,iA)*(1.0d0-(1.0d0+xBar)*Temp1)
        uPrimeBar=2.0d0*B1ca(iC,iA)*xBar*Temp1
        gCalc1=uBar/(xBar*xBar)
        DgCalc1DIx=(uPrimeBar-2.0d0*uBar/xBar)*(0.5d0*A1ca(iC,iA)/SQRTIxCalc)/(xBar*xBar)
      ELSE 
        gCalc1=B1ca(iC,iA)
        DgCalc1DIx=0.0d0
      END IF
      IF (A2ca(iC,iA)>0.0d0) THEN
        xBar=A2ca(iC,iA)*SQRTIxCalc
        Temp1=EXP(-xBar)
        uBar=2.0d0*B2ca(iC,iA)*(1.0d0-(1.0d0+xBar)*Temp1)
        uPrimeBar=2.0d0*B2ca(iC,iA)*xBar*Temp1
        gCalc2=uBar/(xBar*xBar)
        DgCalc2DIx=(uPrimeBar-2.0d0*uBar/xBar)*(0.5d0*A2ca(iC,iA)/SQRTIxCalc)/(xBar*xBar)
      ELSE
        gCalc2=B2ca(iC,iA)
        DgCalc2DIx=0.0d0
      END IF
      Sum=Sum &
         +xCat(iC)*(gCalc1+gCalc2)
      ActCatLoc(iC)=ActCatLoc(iC) &
         +xAn(iA)*(gCalc1+gCalc2)
      DSumDIx=DSumDIx &
         +xCat(iC)*(DgCalc1DIx+DgCalc2DIx)
    END DO
    ActAnLoc(iA)=Sum
    DH2=DH2+xAn(iA)*Sum
    DDHDIx=DDHDIx+xAn(iA)*DSumDIx
  END DO  

  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-DDHDIx*IxCalc+DH1-DH2
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-DDHDIx*(IxCalc-0.5d0*ChargeCat(iC)**2)+ActCatLoc(iC)+DH1-DH2
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-DDHDIx*(IxCalc-0.5d0*ChargeAn(iA)**2)+ActAnLoc(iA)+DH1-DH2
  END DO

!
!     
!   #(02)# Calculate unsymmetrical mixing terms. 
!          ** Act coeff contributions are +UNSYMc, -UNSYMcc and -UNSYMaa ** 
!
  IF (zMax>1) THEN
    DO zi=1,zMax
      DO zj=zi,zMax
        xij=6.0d0*zi*zj*AxCalc*SQRTIxCalc
        CALL JFunc(xij,JCalc(zi,zj),JCalcDx(zi,zj))
      END DO
    END DO
    DO zi=1,zMax
      DO zj=zi+1,zMax
        xij=6.0d0*zi*zj*AxCalc*SQRTIxCalc
        ThetaCalc(zi,zj)=((zi*zj)/(4.0d0*IxCalc))*(JCalc(zi,zj)-0.5d0*(JCalc(zi,zi)+JCalc(zj,zj)))
        ThetaCalcDIx(zi,zj)= &
             -ThetaCalc(zi,zj)/IxCalc &
             +((zi*zj)/(8.0d0*IxCalc*IxCalc)) &
             *(zi*zj*JCalcDx(zi,zj)-0.5d0*(zi*zi*JCalcDx(zi,zi)+zj*zj*JCalcDx(zj,zj))) &
             *6.0d0*AxCalc*SQRTIxCalc
      END DO
    END DO

    ActCatLoc=0.0d0
    ActAnLoc=0.0d0
    HOE=0.0d0
    DHOEDIx=0.0d0
    IF (SIZE(IndCat)>1) THEN
      DO zi=1,SIZE(IndCat)
        SumAC(zi)=0.0d0
        DO iC=1,SIZE(IndCat(zi)%Ind)
          SumAC(zi)=SumAC(zi)+xCat(IndCat(zi)%Ind(iC))
        END DO
      END DO
      DO zi=1,SIZE(IndCat)
        DO zj=zi+1,SIZE(IndCat)
          HOE=HOE+2.0d0*SumAC(zi)*SumAC(zj)*ThetaCalc(zi,zj)
          DHOEDIx=DHOEDIx+2.0d0*SumAC(zi)*SumAC(zj)*ThetaCalcDIx(zi,zj)
          DO iC=1,SIZE(IndCat(zi)%Ind)
            ActCatLoc(IndCat(zi)%Ind(iC))=ActCatLoc(IndCat(zi)%Ind(iC))+2.0d0*SumAC(zj)*ThetaCalc(zi,zj)
          END DO
          DO iC=1,SIZE(IndCat(zj)%Ind)
            ActCatLoc(IndCat(zj)%Ind(iC))=ActCatLoc(IndCat(zj)%Ind(iC))+2.0d0*SumAC(zi)*ThetaCalc(zi,zj)
          END DO
        END DO
      END DO
    END IF

    IF (SIZE(IndAn)>1) THEN
      DO zi=1,SIZE(IndAn)
        SumAC(zi)=0.0d0
        DO iA=1,SIZE(IndAn(zi)%Ind)
          SumAC(zi)=SumAC(zi)+xAn(IndAn(zi)%Ind(iA))
        END DO
      END DO
      DO zi=1,SIZE(IndAn)
        DO zj=zi+1,SIZE(IndAn)
          HOE=HOE+2.0d0*SumAC(zi)*SumAC(zj)*ThetaCalc(zi,zj)
          DHOEDIx=DHOEDIx+2.0d0*SumAC(zi)*SumAC(zj)*ThetaCalcDIx(zi,zj)
          DO iA=1,SIZE(IndAn(zi)%Ind)
            ActAnLoc(IndAn(zi)%Ind(iA))=ActAnLoc(IndAn(zi)%Ind(iA))+2.0d0*SumAC(zj)*ThetaCalc(zi,zj)
          END DO
          DO iA=1,SIZE(IndAn(zj)%Ind)
            ActAnLoc(IndAn(zj)%Ind(iA))=ActAnLoc(IndAn(zj)%Ind(iA))+2.0d0*SumAC(zi)*ThetaCalc(zi,zj)
          END DO
        END DO
      END DO
    END IF
    DO iN=1,nN
      ActNeut(iN)=ActNeut(iN)-DHOEDIx*IxCalc-HOE
    END DO
    DO iC=1,nC
      ActCat(iC)=ActCat(iC)-DHOEDIx*(IxCalc-0.5d0*ChargeCat(iC)**2)+ActCatLoc(iC)-HOE
    END DO
    DO iA=1,nA
      ActAn(iA)=ActAn(iA)-DHOEDIx*(IxCalc-0.5d0*ChargeAn(iA)**2)+ActAnLoc(iA)-HOE
    END DO
  END IF
!
!
!   #(3)# Calculate summation for Wcca parameters: 
!         ** Act coeff contribution is +SUMWCCA **
!
  SumWcca=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iA=1,nA
    Temp2=0.0d0
    ActCatLoc1=0.0d0
    DO i2C=1,nC
      Temp1=0.0d0
      DO i1C=i2C+1,nC
        Temp1=Temp1+xCat(i1C)*Wcca(i1C,i2C,iA)
        ActCatLoc1(i1C)=ActCatLoc1(i1C)+xCat(i2C)*Wcca(i1C,i2C,iA)
      END DO
      Temp2=Temp2+xCat(i2C)*Temp1
      ActCatLoc1(i2C)=ActCatLoc1(i2C)+Temp1
    END DO
    SumWcca=SumWcca+ChargeAn(iA)*xAn(iA)*Temp2
    ActCatLoc=ActCatLoc+ActCatLoc1*ChargeAn(iA)*xAn(iA)
    ActAnLoc(iA)=ChargeAn(iA)*Temp2
  END DO
  SumWcca=2.0d0*SumWcca/SumAnZ
  ActCatLoc=2.0d0*ActCatLoc/SumAnZ
  ActAnLoc=2.0d0*ActAnLoc/SumAnZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-SumWcca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-SumWcca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(1.0d0+ChargeAn(iA)/SumAnZ)*SumWcca
  END DO
!
!
!   #(4)# Calculate summation for Waac parameters: 
!         ** Act coeff contribution is -SUMWAAC **
!
  SumWaac=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iC=1,nC
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO i2A=1,nA
      Temp1=0.0d0
      DO i1A=i2A+1,nA
        Temp1=Temp1+xAn(i1A)*Waac(i1A,i2A,iC)
        ActAnLoc1(i1A)=ActAnLoc1(i1A)+xAn(i2A)*Waac(i1A,i2A,iC)
      END DO
      Temp2=Temp2+xAn(i2A)*Temp1
      ActAnLoc1(i2A)=ActAnLoc1(i2A)+Temp1
    END DO
    SumWaac=SumWaac+ChargeCat(iC)*xCat(iC)*Temp2
    ActAnLoc=ActAnLoc+ActAnLoc1*ChargeCat(iC)*xCat(iC)
    ActCatLoc(iC)=ChargeCat(iC)*Temp2
  END DO
  SumWaac=2.0d0*SumWaac/SumCatZ
  ActCatLoc=2.0d0*ActCatLoc/SumCatZ
  ActAnLoc=2.0d0*ActAnLoc/SumCatZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-SumWaac
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(1.0d0+ChargeCat(iC)/SumCatZ)*SumWaac
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-SumWaac
  END DO
!
!   #(5)# Calculate summation for Ucca parameters: 
!         ** Act coeff contribution is +SUMUCCA **
!
  SumUcca=0.0d0  
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iA=1,nA
    Temp2=0.0d0
    ActCatLoc1=0.0d0
    DO i2C=1,nC
      Temp1=0.0d0
      DO i1C=i2C+1,nC
        Temp1=Temp1+xCat(i1C)*(xCat(i1C)*Vca(i1C,iA)-xCat(i2C)*Vca(i2C,iA))*Ucca(i1C,i2C,iA)
        ActCatLoc1(i1C)=ActCatLoc1(i1C)-xCat(i2C)*(xCat(i2C)*Vca(i2C,iA)-2.0d0*xCat(i1C)*Vca(i1C,iA))*Ucca(i1C,i2C,iA)
        ActCatLoc1(i2C)=ActCatLoc1(i2C)-xCat(i2C)*xCat(i1C)*Vca(i2C,iA)*Ucca(i1C,i2C,iA)
      END DO
      Temp2=Temp2+xCat(i2C)*Temp1
      ActCatLoc1(i2C)=ActCatLoc1(i2C)+Temp1
    END DO
    SumUcca=SumUcca+ChargeAn(iA)*xAn(iA)*Temp2
    ActCatLoc=ActCatLoc+ActCatLoc1*ChargeAn(iA)*xAn(iA)
    ActAnLoc(iA)=ChargeAn(iA)*Temp2
  END DO
  SumUcca=2.0d0*SumUcca/SumAnZ
  ActCatLoc=2.0d0*ActCatLoc/SumAnZ
  ActAnLoc=2.0d0*ActAnLoc/SumAnZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-2.0d0*SumUcca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-2.0d0*SumUcca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(2.0d0+ChargeAn(iA)/SumAnZ)*SumUcca
  END DO
!
!
!   #(6)# Calculate summation for Uaac parameters: 
!         ** Act coeff contribution is -SUMUAAC **
!
  SumUaac=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iC=1,nC
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO i2A=1,nA
      Temp1=0.0d0
      DO i1A=i2A+1,nA
        Temp1=Temp1+xAn(i1A)*(xAn(i1A)*Vac(i1A,iC)-xAn(i2A)*Vac(i2A,iC))*Uaac(i1A,i2A,iC)
        ActAnLoc1(i1A)=ActAnLoc1(i1A)-xAn(i2A)*(xAn(i2A)*Vac(i2A,iC)-2.0d0*xAn(i1A)*Vac(i1A,iC))*Uaac(i1A,i2A,iC)
        ActAnLoc1(i2A)=ActAnLoc1(i2A)-xAn(i2A)*xAn(i1A)*Vac(i2A,iC)*Uaac(i1A,i2A,iC)
      END DO
      Temp2=Temp2+xAn(i2A)*Temp1
      ActAnLoc1(i2A)=ActAnLoc1(i2A)+Temp1
    END DO
    SumUaac=SumUaac+ChargeCat(iC)*xCat(iC)*Temp2
    ActAnLoc=ActAnLoc+ActAnLoc1*ChargeCat(iC)*xCat(iC)
    ActCatLoc(iC)=ChargeCat(iC)*Temp2
  END DO
  SumUaac=2.0d0*SumUaac/SumCatZ
  ActCatLoc=2.0d0*ActCatLoc/SumCatZ
  ActAnLoc=2.0d0*ActAnLoc/SumCatZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-2.0d0*SumUaac
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(2.0d0+ChargeCat(iC)/SumCatZ)*SumUaac
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-2.0d0*SumUaac
  END DO
!
!
!   #(7)# Calculate summation for Wnca parameters: 
!         ** Act coeff contribution is +SUMWnca **
!
  SumWnca=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iN=1,nN
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO iC=1,nC
      Temp1=0.0d0
      DO iA=1,nA
        Temp1=Temp1+xAn(iA)*(ChargeCat(iC)+ChargeAn(iA))*Wnca(iN,iC,iA)
        ActAnLoc2(iA)=(ChargeCat(iC)+ChargeAn(iA))*Wnca(iN,iC,iA)
      END DO
      Temp2=Temp2+xCat(iC)*Temp1
      ActAnLoc1=ActAnLoc1+xCat(iC)*ActAnLoc2
      ActCatLoc1(ic)=Temp1
    END DO
    SumWnca=SumWnca+xNeut(iN)*Temp2
    ActNeutLoc(iN)=Temp2
    ActAnLoc=ActAnLoc+xNeut(iN)*ActAnLoc1
    ActCatLoc=ActCatLoc+xNeut(iN)*ActCatLoc1
  END DO
  Temp1=0.5d0*(1.0d0/SumAnZ+1.0d0/SumCatZ)
  SumWnca=SumWnca*Temp1
  ActCatLoc=ActCatLoc*Temp1
  ActAnLoc=ActAnLoc*Temp1
  ActNeutLoc=ActNeutLoc*Temp1
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ActNeutLoc(iN)-SumWnca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(1.0d0+0.5d0*ChargeCat(iC)/SumCatZ)*SumWnca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(1.0d0+0.5d0*ChargeAn(iA)/SumAnZ)*SumWnca
  END DO
!
!   #(8)# Calculate summation for Unca parameters:
!         ** Act coeff contribution is +SUMUnca **
!
  SumUnca=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iN=1,nN
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO iC=1,nC
      Temp1=0.0d0
      DO iA=1,nA
        Temp1=Temp1+xAn(iA)*(ChargeCat(iC)+ChargeAn(iA))**2 &
             /(FLOAT(ChargeCat(iC))*FLOAT(ChargeAn(iA)))*Unca(iN,iC,iA)
        ActAnLoc2(iA)=(ChargeCat(iC)+ChargeAn(iA))**2 &
                     /(FLOAT(ChargeCat(iC))*FLOAT(ChargeAn(iA)))*Unca(iN,iC,iA)
      END DO
      Temp2=Temp2+xCat(iC)*Temp1
      ActAnLoc1=ActAnLoc1+xCat(iC)*ActAnLoc2
      ActCatLoc1(ic)=Temp1
    END DO
    SumUnca=SumUnca+xNeut(iN)*Temp2
    ActNeutLoc(iN)=Temp2
    ActAnLoc=ActAnLoc+xNeut(iN)*ActAnLoc1
    ActCatLoc=ActCatLoc+xNeut(iN)*ActCatLoc1
  END DO
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ActNeutLoc(iN)-2.0d0*SumUnca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-2.0d0*SumUnca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-2.0d0*SumUnca
  END DO
!
!   #(9)# Calculate summation for Vnca parameters: 
!         ** Act coeff contribution is +SUMVnca **
!
  SumVnca=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iN=1,nN
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO iC=1,nC
      Temp1=0.0d0
      DO iA=1,nA
        Temp1=Temp1+xAn(iA)*Vnca(iN,iC,iA)
        ActAnLoc2(iA)=Vnca(iN,iC,iA)
      END DO
      Temp2=Temp2+xCat(iC)*Temp1
      ActAnLoc1=ActAnLoc1+xCat(iC)*ActAnLoc2
      ActCatLoc1(ic)=Temp1
    END DO
    SumVnca=SumVnca+xNeut(iN)**2*Temp2
    ActNeutLoc(iN)=2.0d0*xNeut(iN)*Temp2
    ActAnLoc=ActAnLoc+xNeut(iN)**2*ActAnLoc1
    ActCatLoc=ActCatLoc+xNeut(iN)**2*ActCatLoc1
  END DO
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+4.0d0*(ActNeutLoc(iN)-3.0d0*SumVnca)
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+4.0d0*(ActCatLoc(iC)-3.0d0*SumVnca)
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+4.0d0*(ActAnLoc(iA)-3.0d0*SumVnca)
  END DO
!
!   #(10)# Calculate summation for Qncca parameters: 
!         ** Act coeff contribution is +SUMQncca **
!
  SumQncca=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iN=1,nN
    Temp3=0.0d0
    ActCatLoc1=0.0d0
    DO iA=1,nA
      Temp2=0.0d0
      ActCatLoc2=0.0d0
      DO i2C=1,nC
        Temp1=0.0d0
        DO i1C=i2C+1,nC
          Temp1=Temp1+xCat(i1C)*Qncca(iN,i1C,i2C,iA)
          ActCatLoc2(i1C)=ActCatLoc2(i1C)+xCat(i2C)*Qncca(iN,i1C,i2C,iA)
        END DO
        Temp2=Temp2+xCat(i2C)*Temp1
        ActCatLoc2(i2C)=ActCatLoc2(i2C)+Temp1
      END DO
      Temp3=Temp3+ChargeAn(iA)*xAn(iA)*Temp2
      ActCatLoc1=ActCatLoc1+ChargeAn(iA)*xAn(iA)*ActCatLoc2
      ActAnLoc1(iA)=ChargeAn(iA)*Temp2
    END DO
    SumQncca=SumQncca+xNeut(iN)*Temp3
    ActCatLoc=ActCatLoc+xNeut(iN)*ActCatLoc1
    ActAnLoc=ActAnLoc+xNeut(iN)*ActAnLoc1
    ActNeutLoc(iN)=Temp3
  END DO
  SumQncca=4.0d0*SumQncca/SumCatZ
  ActCatLoc=4.0d0*ActCatLoc/SumCatZ
  ActAnLoc=4.0d0*ActAnLoc/SumCatZ
  ActNeutLoc=4.0d0*ActNeutLoc/SumCatZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ActNeutLoc(iN)-2.0d0*SumQncca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-2.0d0*SumQncca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(2.0d0+ChargeAn(iA)/SumAnZ)*SumQncca
  END DO
!
!   #(11)# Calculate summation for Qnaac parameters: 
!          ** Act coeff contribution is -SUMQnaac **
!
  SumQnaac=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iN=1,nN
    Temp3=0.0d0
    ActAnLoc1=0.0d0
    DO iC=1,nC
      Temp2=0.0d0
      ActAnLoc2=0.0d0
      DO i2A=1,nA
        Temp1=0.0d0
        DO i1A=i2A+1,nA
          Temp1=Temp1+xAn(i1A)*Qnaac(iN,i1A,i2A,iC)
          ActAnLoc2(i1A)=ActAnLoc2(i1A)+xAn(i2A)*Qnaac(iN,i1A,i2A,iC)
        END DO
        Temp2=Temp2+xAn(i2A)*Temp1
        ActAnLoc2(i2A)=ActAnLoc2(i2A)+Temp1
      END DO
      Temp3=Temp3+ChargeCat(iC)*xCat(iC)*Temp2
      ActAnLoc1=ActAnLoc1+ChargeCat(iC)*xCat(iC)*ActAnLoc2
      ActCatLoc1(iC)=ChargeCat(iC)*Temp2
    END DO
    SumQnaac=SumQnaac+xNeut(iN)*Temp3
    ActCatLoc=ActCatLoc+xNeut(iN)*ActCatLoc1
    ActAnLoc=ActAnLoc+xNeut(iN)*ActAnLoc1
    ActNeutLoc(iN)=Temp3
  END DO
  SumQnaac=4.0d0*SumQnaac/SumAnZ
  ActCatLoc=4.0d0*ActCatLoc/SumAnZ
  ActAnLoc=4.0d0*ActAnLoc/SumAnZ
  ActNeutLoc=4.0d0*ActNeutLoc/SumAnZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ActNeutLoc(iN)-2.0d0*SumQnaac
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(2.0d0+ChargeCat(iC)/SumCatZ)*SumQnaac
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-2.0d0*SumQnaac
  END DO
!
!   #(12)# Calculate summation for Ynnca parameters: 
!         ** Act coeff contribution is +SUMYnnca **
!

! noch nicht getestet nN>=2 !
  SumYnnca=0.D0
  ActAnLoc=0.0d0
  ActCatLoc=0.0d0
  ActNeutLoc=0.0d0
  DO  i2N=1,nN
    Temp3=0.0d0
    ActAnLoc1=0.0d0
    ActCatLoc1=0.0d0
    DO i1N=i2N+1,nN
      Temp2=0.0d0
      ActAnLoc2=0.0d0
      DO iC=1,nC
        Temp1=0.0d0
        DO iA=1,nA
          Temp1=Temp1+xAn(iA)*(ChargeAn(iA)+ChargeCat(iC))*Ynnca(i2N,i1N,iC,iA)
          ActAnLoc3(iA)=(ChargeAn(iA)+ChargeCat(iC))*Ynnca(i2N,i1N,iC,iA)
        END DO
        Temp2=Temp2+xCat(iC)*Temp1
        ActAnLoc2=ActAnLoc2+xCat(iC)*ActAnLoc3
        ActCatLoc2=Temp1
      END DO
      Temp3=Temp3+xNeut(i1N)*Temp2
      ActAnLoc1=ActAnLoc1+xNeut(i1N)*ActAnLoc2
      ActCatLoc1=ActCatLoc1+xNeut(i1N)*ActCatLoc2
      ActNeutLoc(i1N)=ActNeutLoc(i1N)+xNeut(i2N)*Temp2
    END DO
    SumYnnca=SumYnnca+xNeut(i2N)*Temp3
    ActAnLoc=ActAnLoc+xNeut(i2N)*ActAnLoc1
    ActCatLoc=ActCatLoc+xNeut(i2N)*ActCatLoc1
    ActNeutLoc(i2N)=ActNeutLoc(i2N)+Temp3
  END DO
  Temp1=0.5d0*(1.0d0/SumAnZ+1.0d0/SumCatZ)
  SumYnnca=SumYnnca*Temp1
  ActCatLoc=ActCatLoc*Temp1
  ActAnLoc=ActAnLoc*Temp1
  ActNeutLoc=ActNeutLoc*Temp1
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ActNeutLoc(iN)-2.0d0*SumYnnca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(2.0d0+0.5d0*ChargeCat(iC)/SumCatZ)*SumYnnca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(2.0d0+0.5d0*ChargeAn(iA)/SumAnZ)*SumYnnca
  END DO
!
!   #(13)# Calculate summation for Xccca parameters: 
!          ** Act coeff contribution is +SUMXccca **
!
  SumXccca=0.D0
  ActAnLoc=0.0d0
  ActCatLoc=0.0d0
  DO iA=1,nA
    ActCatLoc1=0.0d0
    Temp3=0.0d0 
    DO i3C=1,nC
      Temp2=0.0d0 
      DO i2C=i3C+1,nC
        Temp1=0.0d0 
        DO i1C=i2C+1,nC
          Temp1=Temp1+xCat(i1C)*Xccca(i1C,i2C,i3C,iA)
          ActCatLoc1(i1C)=ActCatLoc1(i1C)+xCat(i2C)*xCat(i3C)*Xccca(i1C,i2C,i3C,iA)
        END DO
        Temp2=Temp2+xCat(i2C)*Temp1
        ActCatLoc1(i2C)=ActCatLoc1(i2C)+xCat(i3C)*Temp1
      END DO
      Temp3=Temp3+xCat(i3C)*Temp2
      ActCatLoc1(i3C)=ActCatLoc1(i3C)+Temp2
    END DO
    SumXccca=SumXccca+ChargeAn(iA)*xAn(iA)*Temp3
    ActAnLoc(iA)=ChargeAn(iA)*Temp3
    ActCatLoc=ActCatLoc+ChargeAn(iA)*xAn(iA)*ActCatLoc1
  END DO
  SumXccca=SumXccca/SumAnZ
  ActCatLoc=ActCatLoc/SumAnZ
  ActAnLoc=ActAnLoc/SumAnZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-2.0d0*SumXccca
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-2.0d0*SumXccca
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(2.0d0+ChargeAn(iA)/SumAnZ)*SumXccca
  END DO
!
!
!   #(14)# Calculate summation for Xaaac parameters: 
!         ** Act coeff contribution is -SUMXaaac **
!
  SumXaaac=0.D0
  ActAnLoc=0.0d0
  ActCatLoc=0.0d0
  DO iC=1,nC
    ActAnLoc1=0.0d0
    Temp3=0.0d0 
    DO i3A=1,nA
      Temp2=0.0d0 
      DO i2A=i3A+1,nA
        Temp1=0.0d0 
        DO i1A=i2A+1,nA
          Temp1=Temp1+xAn(i1A)*Xaaac(i1A,i2A,i3A,iC)
          ActAnLoc1(i1A)=ActAnLoc1(i1A)+xAn(i2A)*xAn(i3A)*Xaaac(i1A,i2A,i3A,iC)
        END DO
        Temp2=Temp2+xAn(i2A)*Temp1
        ActAnLoc1(i2A)=ActAnLoc1(i2A)+xAn(i3A)*Temp1
      END DO
      Temp3=Temp3+xAn(i3A)*Temp2
      ActAnLoc1(i3A)=ActAnLoc1(i3A)+Temp2
    END DO
    SumXaaac=SumXaaac+ChargeCat(iC)*xCat(iC)*Temp3
    ActCatLoc(iC)=ChargeCat(iC)*Temp3
    ActAnLoc=ActAnLoc+ChargeCat(iC)*xCat(iC)*ActAnLoc1
  END DO
  SumXaaac=SumXaaac/SumCatZ
  ActCatLoc=ActCatLoc/SumCatZ
  ActAnLoc=ActAnLoc/SumCatZ
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-2.0d0*SumXaaac
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(2.0d0+ChargeCat(iC)/SumCatZ)*SumXaaac
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-2.0d0*SumXaaac
  END DO
!
!
!   #(15)# Calculate summation for Zccaa parameters: 
!         ** Act coeff contribution is +SUMZccaa **
!
  SumZccaa=0.D0
  ActAnLoc=0.0d0
  ActCatLoc=0.0d0
  DO i2C=1,nC
    Temp3=0.0d0
    ActAnLoc2=0.0d0
    DO i1C=i2C+1,nC
      ActAnLoc1=0.0d0
      Temp2=0.0d0
      DO i2A=1,nA
        Temp1=0.0d0
        DO i1A=i2A+1,nA
          Temp1=Temp1+xAn(i1A)*Zccaa(i1C,i2C,i1A,i2A)
          ActAnLoc1(i1A)=ActAnLoc1(i1A)+xAn(i2A)*Zccaa(i1C,i2C,i1A,i2A)
        END DO
        Temp2=Temp2+xAn(i2A)*Temp1
        ActAnLoc1(i2A)=ActAnLoc1(i2A)+Temp1
      END DO
      Temp3=Temp3+xCat(i1C)*Temp2
      ActAnLoc2=ActAnLoc2+xCat(i1C)*ActAnLoc1
      ActCatLoc(i1C)=ActCatLoc(i1C)+xCat(i2C)*Temp2
    END DO
    SumZccaa=SumZccaa+xCat(i2C)*Temp3
    ActCatLoc(i2C)=ActCatLoc(i2C)+Temp3
    ActAnLoc=ActAnLoc+xCat(i2C)*ActAnLoc2
  END DO
  SumZccaa=SumZccaa/(0.5d0*(SumCatZ+SumAnZ))
  ActCatLoc=ActCatLoc/(0.5d0*(SumCatZ+SumAnZ))
  ActAnLoc=ActAnLoc/(0.5d0*(SumCatZ+SumAnZ))
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)-2.0d0*SumZccaa
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)-(2.0d0+ChargeCat(iC)/(SumCatZ+SumAnZ))*SumZccaa
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)-(2.0d0+ChargeAn(iA)/(SumCatZ+SumAnZ))*SumZccaa
  END DO
!
!
!   #(16)# Calculate summation for Wnn and Unn parameters: 
!         ** Act coeff contribution is -SUMnn **
!
  Sumnn = 0.D0
  ActNeutLoc=0.0d0
  DO i2N=1,nN
    Temp1=0.0d0
    Temp2=0.0d0
    DO i1N=i2N+1,nN
      Temp1=Temp1+xNeut(i1N)*(Wnn(i1N,i2N)+2.0d0*Unn(i1N,i2N)*xNeut(i1N))
      Temp2=Temp2-2.0d0*xNeut(i1N)*Unn(i1N,i2N)
      ActNeutLoc(i1N)=ActNeutLoc(i1N) &
                     +xNeut(i2N)*(Wnn(i1N,i2N)+4.0d0*Unn(i1N,i2N)*xNeut(i1N) &
                     -2.0d0*xNeut(i2N)*Unn(i1N,i2N))
    END DO
    Sumnn=Sumnn+xNeut(i2N)*(Temp1+xNeut(i2N)*Temp2)
    ActNeutLoc(i2N)=ActNeutLoc(i2N)+Temp1+2.0d0*xNeut(i2N)*Temp2
  END DO
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ ActNeutLoc(iN)-2.0d0*Sumnn
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-2.0d0*Sumnn
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-2.0d0*Sumnn
  END DO
!
!
!   #(17)# Calculate summation for Cnnn parameters: 
!         ** Act coeff contribution is -SUMCnnn **
!
  SumCnnn=0.0d0
  ActNeutLoc=0.0d0
  DO i3N=1,nN
    Temp2=0.0d0
    DO i2N=i3N+1,nN
      DO i1N=i2N+1,nN
        Temp1=Temp1-2.0d0*xNeut(i1N)*Cnnn(i1N,i2N,i3N)
        ActNeutLoc(i1N)=ActNeutLoc(i1N)-2.0d0*xNeut(i2N)*xNeut(i3N)*Cnnn(i1N,i2N,i3N)
      END DO
      Temp2=Temp2+xNeut(i2N)*Temp1
      ActNeutLoc(i2N)=ActNeutLoc(i2N)+xNeut(i3N)*Temp1
    END DO
    SumCnnn=SumCnnn+xNeut(i3N)*Temp2
    ActNeutLoc(i3N)=ActNeutLoc(i3N)+Temp1
  END DO
  DO iN=1,nN
    ActNeut(iN)=ActNeut(iN)+ ActNeutLoc(iN)-2.0d0*Sumnn
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-2.0d0*Sumnn
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-2.0d0*Sumnn
  END DO
!
!
!   #(18)# Calculate correction to infinite dilution reference state:
!          ** Act coeff contribution is -CORRact **
!          ** Expression below is formulated assuming ref. state of
!             infinite dilution with respect to a pure (single) solvent
!             which *must* be neutral '1' **

  CorrAct=0.D0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iC=1,nC
    Temp1=0.0d0
    DO iA=1,nA
      Temp1=Temp1+xAn(iA)*(ChargeCat(iC)+ChargeAn(iA))*Wnca(1,iC,iA)
      ActAnLoc1(iA)=(ChargeCat(iC)+ChargeAn(iA))*Wnca(1,iC,iA)
    END DO
    CorrAct=CorrAct+xCat(iC)*Temp1
    ActAnLoc=ActAnLoc+xCat(iC)*ActAnLoc1
    ActCatLoc(ic)=Temp1
  END DO

  Temp1=0.5d0*(1.0d0/SumAnZ+1.0d0/SumCatZ)
  CorrAct=CorrAct*Temp1
  ActCatLoc=ActCatLoc*Temp1
  ActAnLoc=ActAnLoc*Temp1
  ActNeut(1)=EXP(ActNeut(1))
  DO iN=2,nN
    ActNeut(iN)=ActNeut(iN)-(Wnn(iN,1)-uNn(iN,1))
    ActNeut(iN)=EXP(ActNeut(iN))
  END DO
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-ActCatLoc(iC)+0.5d0*ChargeCat(iC)/SumCatZ*CorrAct
    ActCat(iC)=EXP(ActCat(iC))
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-ActAnLoc(iA)+0.5d0*ChargeAn(iA)/SumAnZ*CorrAct
    ActAn(iA)=EXP(ActAn(iA))
  END DO

END SUBROUTINE ExcessGibbs

SUBROUTINE SetCoefficient(T)

  REAL(RealKind) :: T ! Temperature 

  REAL(RealKind) :: D1,D4
  REAL(RealKind) :: fn1,fn2,fn3,fn4,fn5,fn6,fn7
  REAL(RealKind) :: TR,BL1,BL2,WL,UL,VL
  REAL(RealKind) :: BJ1,BJ2,WJ,UJ,VJ
  REAL(RealKind) :: d2B1,d2B2,d2W,d2U,d2V
!
!
!
!  ----------------------------------------------------------------
! |  SECTION (1): H2SO4, HNO3, HCl and HBr. Parameters from the    |
! |  Carslaw et al 1995 paper, and more recent work by me (HBr).   |
!  ----------------------------------------------------------------
      D1 = T-328.15D0
 
!  ** H - HSO4 **
 
      A1ca(1,1) = 17.D0
      B1ca(1,1) =            0.178334467D2 + &
                        D1*(-0.625268629D1*1.D-1 + & 
                   D1*(0.5D0*0.591429323D0*1.D-2 + &
                         D1*(0.223751841D0*1.D-3/6.D0 + & 
                         D1*(0.D0         *1.D-3/12.D0 + &
                          D1*0.D0         *1.D-3/20.D0)))) 
      Wnca(1,1,1) =         -0.998416390D1 + & 
                         D1*(0.348821776D0  *1.D-1 +  &
                 D1*(0.5D0*(-0.119526172D-1)*1.D-2 + &
                         D1*(0.909425662D-2 *1.D-3/6.D0 + & 
                         D1*(0.149166944D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0)))) 
      Unca(1,1,1) =         -0.143238371D1 +  &
                        D1*(-0.201636224D0  *1.D-1 + & 
                 D1*(0.5D0*(-0.443804232D-1)*1.D-2 + &
                         D1*(0.641847819D-2 *1.D-3/6.D0 + & 
                         D1*(0.296327801D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0)))) 
      Vnca(1,1,1) =         -0.207474566D1 + & 
                         D1*(0.594737744D0  *1.D-1 + & 
                   D1*(0.5D0*0.674052221D-1 *1.D-2 + &
                         D1*(0.D0           *1.D-3/6.D0 + & 
                        D1*(-0.394845016D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
!
!
!  ** H - SO4 **
!
      A1ca(1,2) = 9.5D0
      B1ca(1,2) =           -0.982408701D2 + &
                        D1*(-0.205401806D2  *1.D-1 + & 
                 D1*(0.5D0*(-0.207137292D1) *1.D-2 + &
                        D1*(-0.376521937D-1 *1.D-3/6.D0 + & 
                        D1*(-0.139689758D-1 *1.D-3/12.D0 + &
                         D1* 0.D0           *1.D-3/20.D0))))
      Wnca(1,1,2) =         -0.107752155D2   + & 
                        D1*(-0.879298257D0  *1.D-1 + & 
                 D1*(0.5D0*(-0.440528485D0) *1.D-2 + &
                        D1*(-0.544913927D-1 *1.D-3/6.D0 + & 
                        D1*(-0.173541364D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
      Unca(1,1,2) =         -0.133603464D2 + & 
                        D1*(-0.459479578D1  *1.D-1 + & 
                 D1*(0.5D0*(-0.146220346D1) *1.D-2 + &
                        D1*(-0.157872023D0  *1.D-3/6.D0 + & 
                        D1*(-0.162230945D-3 *1.D-3/12.D0 + &
                         D1* 0.D0           *1.D-3/20.D0))))
      Vnca(1,1,2) =          0.310121997D1    + & 
                         D1*(0.446189009D1  *1.D-1 + & 
                   D1*(0.5D0*0.975254718D0  *1.D-2 + &
                         D1*(0.588748231D-2 *1.D-3/6.D0 + & 
                        D1*(-0.901983372D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
!
!
!  ** H - NO3 **
!
      D4  = T-330.D0
      fn1 = 1.D-1
      fn2 = 1.D-2/2.D0
      fn3 = 1.D-3/6.D0
      fn4 = 1.D-4/12.D0
      fn5 = 1.D-5/20.D0
      fn6 = 1.D-6/30.D0
      fn7 = 1.D-7/42.D0
!      
      A1ca(1,3) = 29.03169775D0
      B1ca(1,3) =             0.2470883216D2 + & 
                         d4*( 0.8537056072D1 *fn1 + &
                         d4*( 0.D0           *fn2 + &
                         d4*(-0.5686459898D1 *fn3 + &
                         d4*(-0.2282434425D1 *fn4 + &
                         d4*(-0.2606309141D0 *fn5 + &
                         d4*(-0.9318482814D-2*fn6 + &
                         d4*  0.D0           *fn7))))))
      Wnca(1,1,3) =          -0.3301004522D1 + & 
                         d4*( 0.1812617746D0 *fn1 + &
                         d4*( 0.1138792155D-1*fn2 + &
                         d4*( 0.3112011367D-2*fn3 + &
                         d4*( 0.D0           *fn4 + &
                         d4*( 0.D0           *fn5 + &
                         d4*( 0.5981401888D-5*fn6 + &
                          d4* 0.D0           *fn7))))))
      Unca(1,1,3) =          -0.1009606757D1  + & 
                         d4*(-0.9191466033D-1*fn1 + &
                         d4*( 0.3513876226D-1*fn2 + &
                         d4*( 0.D0           *fn3 + &
                         d4*( 0.D0           *fn4 + &
                         d4*( 0.4419544700D-3*fn5 + &
                         d4*( 0.4469965190D-4*fn6 + &
                          d4* 0.D0           *fn7))))))
      Vnca(1,1,3) =           0.1807488767D1 + & 
                         d4*( 0.D0           *fn1 + &
                         d4*(-0.2615013798D0 *fn2 + &
                         d4*(-0.5744350823D-1*fn3 + &
                         d4*( 0.D0           *fn4 + &
                         d4*( 0.5151608617D-3*fn5 + &
                         d4*( 0.D0           *fn6 + &
                          d4* 0.D0           *fn7))))))
!
!
!  ** H-Cl **
!
      D4  = T-330.D0
      fn1 = 1.D-1
      fn2 = 1.D-2/2.D0
      fn3 = 1.D-3/6.D0
      fn4 = 1.D-4/12.D0
      fn5 = 1.D-5/20.D0
      fn6 = 1.D-6/30.D0
      fn7 = 1.D-7/42.D0
!
      A1ca(1,4) = 2.43288D0
      B1ca(1,4) = 0.637130453D2 + d4*(-0.127961571D1*fn1 + &
                      d4*(-0.855184587D0*fn2 + &
                       d4*(-0.483462259D0*fn3 + &
                        d4*(-0.769606918D-1*fn4 + &
                         d4*(0.227428690D-3*fn5 + &
                          d4*(0.D0*fn6 + &
                           d4*0.D0*fn7))))))
      Wnca(1,1,4) = 0.D0 + d4*(0.712414208D0*fn1 + &
                      d4*(0.D0*fn2 + &
                       d4*(-0.253763077D-1*fn3 + &
                        d4*(-0.725301184D-2*fn4 + &
                         d4*(0.D0*fn5 + &
                          d4*(0.D0*fn6 + &
                           d4*0.D0*fn7))))))
      Unca(1,1,4) = 0.757587583D1 + d4*(0.121792108D1*fn1 + &
                      d4*(0.230535499D0*fn2 + &
                       d4*(0.576304748D-1*fn3 + &
                        d4*(0.D0*fn4 + &
                         d4*(0.D0*fn5 + &
                          d4*(0.D0*fn6 + &
                           d4*0.D0*fn7))))))
      Vnca(1,1,4) = -0.153591245D2 + d4*(-0.389164072D0*fn1 + &
                      d4*(-0.891045013D-1*fn2 + &
                       d4*(0.D0*fn3 + &
                        d4*(0.654055055D-2*fn4 + &
                         d4*(0.D0*fn5 + &
                          d4*(0.D0*fn6 + &
                           d4*0.D0*fn7))))))
!
!   ** H - Br ** (added 9/5/98)
      A1ca(1,5) = 6.19176D0
      B1ca(1,5) =     0.327847880D2  + & 
                     d4*(0.919647588D0    *fn1 + &
                     d4*(0.742200662D0    *fn2 + &
                     d4*(0.263972123D0    *fn3 + &
                     d4*(0.D0             *fn4 + &
                     d4*(0.616206313D-3   *fn5 + &
                     d4*(0.D0             *fn6 + &
                     d4* 0.D0             *fn7))))))
      Wnca(1,1,5) =  -0.131504151D2  + & 
                     d4*(0.116060928D1    *fn1 + &
                     d4*(0.312135679D0    *fn2 + &
                     d4*(0.146033030D0    *fn3 + &
                     d4*(0.708445887D-2   *fn4 + &
                     d4*(0.D0             *fn5 + &
                     d4*(0.D0             *fn6 + &
                     d4* 0.D0             *fn7))))))
      Unca(1,1,5) =  -0.699721143D1  + &
                     d4*(0.215281904D1    *fn1 + &
                     d4*(0.902856022D0    *fn2 + &
                     d4*(0.410269215D0    *fn3 + &
                     d4*(0.214487024D-1   *fn4 + &
                     d4*(0.D0             *fn5 + &
                     d4*(0.D0             *fn6 + &
                     d4* 0.D0             *fn7))))))
      Vnca(1,1,5) =  -0.414565596D1  + & 
                     d4*(-0.124431743D1   *fn1 + &
                     d4*(-0.748579644D0   *fn2 + &
                     d4*(-0.319716995D0   *fn3 + &
                     d4*(-0.161436912D-1  *fn4 + &
                     d4*(0.D0             *fn5 + &
                     d4*(0.D0             *fn6 + &
                     d4* 0.D0             *fn7))))))
!
!
!  ** HSO4 - NO3 - H **
      Waac(1,3,1) = -4.280D0
      Uaac(1,3,1) =  0.201362D0 + 0.084830D0*(T-273.15D0)
      Waac(3,1,1) =  Waac(1,3,1)
      Uaac(3,1,1) = -Uaac(1,3,1)
!
!  ** SO4 - NO3 - H **
      Waac(2,3,1) = -0.033291D0*(T-273.15D0)
      Waac(3,2,1) =  Waac(2,3,1)
!
!
!.....HSO4-Br-H interactions: (14 Feb 10, 2001)
      Waac(1,5,1) = -17.1133D0
      Waac(5,1,1) = Waac(1,5,1)
      Qnaac(1,1,5,1) = 15.4803D0 - 0.558433D-2*(T - 298.15D0)
      Qnaac(1,5,1,1) = Qnaac(1,1,5,1)
!
!.....SO4-Br-H interactions: (14 Feb 10, 2001)
      Waac(2,5,1) = -5.31288D0 
      IF(T .LE. 298.15D0) THEN
        Waac(2,5,1) = Waac(2,5,1) - 0.488750D0*(T - 298.15D0)
      ELSE
!     ..halving the slope above 25oC is a fix to roughly improve
!       predictions. It may be better to have an entirely separate
!       parameterisation above this temperature.
        Waac(2,5,1) = Waac(2,5,1) - 0.488750D0*(T - 298.15D0)/2.D0
      ENDIF
      Waac(5,2,1) = Waac(2,5,1)
!
!
!  -------------------------------------------------------------
! |  SECTION (2): OTHER ELEMENTS OF THE H-NH4-NO3-HSO4-SO4-H2O  |
! |  MODEL.                                                     |
!  -------------------------------------------------------------
!
!  ** NH4 - SO4 ** (Values from 1995 JCED paper)
!
      TR        = 298.15D0
      A1ca(2,2) = 13.D0
      B1ca(2,2) = 0.1399385D+02 + (T-TR)*(0.D0 - TR*(-0.02103172D0)) &
                    + 0.5D0*(T**2 - TR**2)*(-0.02103172D0)
!
      A2ca(2,2) = 1.5D0
      B2ca(2,2)=-0.1713243D2 + (T-TR)*(0.8461758D0-TR*(-0.9880561D-3)) &
                   + 0.5D0*(T**2 - TR**2)*(-0.9880561D-3)
!
      Wnca(1,2,2)=-0.1904921D1 + (T-TR)*(0.07911813D0-TR*(0.1173141D-3)) &
                    + 0.5D0*(T**2 - TR**2)*(0.1173141D-3) 
!
      Unca(1,2,2)= 0.2125957D1 + (T-TR)*(-7.160173D-3-TR*(0.3064704D-3)) &
                    + 0.5D0*(T**2 - TR**2)*(0.3064704D-3) 
!
      Vnca(1,2,2)=-0.2291087D1+(T-TR)*(-0.03772486D0-TR*(-1.865885D-4)) &
                    + 0.5D0*(T**2 - TR**2)*(-1.865885D-4)
!
!
!  ** NH4 - HSO4 ** (June '97 revisions)
!
      D1 = T-298.15D0
!
      A1ca(2,1)   =  13.D0
      B1ca(2,1)   =  0.495871301D+01 + D1*(0.4D0*1.D-1) 
      Wnca(1,2,1) = -0.163376858D+01 + D1*(-0.539327161D-1*1.D-1) 
      Unca(1,2,1) =  0.362701710D+00
      Vnca(1,2,1) =  0.367166597D+00 + D1*(-0.74241954D-1*1.D-1) 
!
!
!  ** NH4 - NO3 ** (all parameters revised May '97)
!
      TR  =  298.15D0
      BL1 =  0.283192D0
      BL2 = -0.0716627D0
      WL  = -0.00696723D0
      UL  =  0.171736D-2 
      VL  =  0.221706D-2 
!
      BJ1 = -0.352077D-02
      BJ2 =  0.D0
      WJ  =  0.489933D-5 
      UJ  =  0.D0
      VJ  = -0.359406D-4 
!
      d2B1 = BJ1 - (2/TR)*BL1
      d2B2 = BJ2 - (2/TR)*BL2
      d2W  = WJ - (2/TR)*WL
      d2U  = UJ - (2/TR)*UL
      d2V  = VJ - (2/TR)*VL
!
      A1ca(2,3)=7.D0
      A2ca(2,3)=13.D0
      B1ca(2,3)=0.130466D2+(T-TR)*(BL1-TR*d2B1)+0.5D0*(T**2-TR**2)*d2B1 
      B2ca(2,3)=-0.162254D2+(T-TR)*(BL2-TR*d2B2)+0.5D0*(T**2-TR**2)*d2B2 
      Wnca(1,2,3)=0.616136D0+(T-TR)*(WL-TR*d2W)+0.5D0*(T**2-TR**2)*d2W
      Unca(1,2,3)=-0.403564D-1+(T-TR)*(UL-TR*d2U)+0.5D0*(T**2-TR**2)*d2U
      Vnca(1,2,3)=-0.680507D0+(T-TR)*(VL-TR*d2V)+0.5D0*(T**2-TR**2)*d2V 
!
!
!  ** HSO4 - SO4 - NH4 ** (revised June '97)
      D1 = T - 298.15D0
      Waac(1,2,2) = -8.17099584D0
      Uaac(1,2,2) = -11.3871516D0 + D1*(0.633515474D0*1.D-1) 
      Waac(2,1,2) = Waac(1,2,2)
      Uaac(2,1,2) = -Uaac(1,2,2)
!
!  ** H - NH4 - HSO4 ** (revised June '97)
      D1 = T - 298.15D0
      Wcca(1,2,1)    = -16.4807039D0 + D1*(0.494791553D0*1.D-1) 
      Qncca(1,1,2,1) = 7.084464D0 
      Wcca(2,1,1)    = Wcca(1,2,1)
      Qncca(1,2,1,1) = Qncca(1,1,2,1)
!
!  ** H - NH4 - SO4 ** (revised June '97)
      Wcca(1,2,2)    = -8.33078088D0
      Qncca(1,1,2,2) = 2.60624775D0
      Wcca(2,1,2)    = Wcca(1,2,2)
      Qncca(1,2,1,2) = Qncca(1,1,2,2)
!
!  ** H - NH4 - NO3 ** (revised May and June '97):
      D1 = T - 298.15D0
      Wcca(1,2,3)    = -0.40728D1 + 0.42031D-01*D1 & 
                       - 0.32245D-03*D1**2
      Qncca(1,1,2,3) = 0.75481D0
!   ..NH4 first, to correspond to order in fitting program:
      Ucca(2,1,3) = 0.10055D1
      Wcca(2,1,3) = Wcca(1,2,3)
      Qncca(1,2,1,3) =  Qncca(1,1,2,3)
      Ucca(1,2,3)    = -Ucca(2,1,3)
!
!  ** SO4 - NO3 - NH4 ** (revised 2 June '97)
      Waac(2,3,2)    = 0.32359D0+0.16567D0*(T-298.15D0)
      Qnaac(1,2,3,2) = 2.6800D0-0.13322D0*(T-298.15D0)
!   ..sign reversal below, as NO3 before SO4 in fitting program:
      Uaac(2,3,2) = -0.31582D0
      Waac(3,2,2)    =  Waac(2,3,2)
      Qnaac(1,3,2,2) =  Qnaac(1,2,3,2)
      Uaac(3,2,2)    = -Uaac(2,3,2)
!
!  ** HSO4 - NO3 - NH4 ** (June '97)
      Waac(1,3,2) = -3.07081D0
      Waac(3,1,2) =  Waac(1,3,2)
!
END SUBROUTINE SetCoefficient

SUBROUTINE SetCoefficient_25(T)

  REAL(RealKind) :: T ! Temperature 

  REAL(RealKind) :: D1

!
!
!  *********************************************************************
!  *                                                                   *
!  * This routine last updated: 10/3/97                                *
!  * for 25oC limited H-Na-NH4-Cl-NO3-SO4-H2O fit                      *
!  *                                                                   *
!  *********************************************************************
!
!
!
! ====================================================================
! ====================================================================
!    SECTION (1) FOR H2SO4, HNO3 AND HCL. 
!

!
!  ** H - HSO4 ** (Carslaw et al. model)
!
      D1=T-328.15D0
!
      A1ca(1,1)=17.D0
      B1ca(1,1)=             0.178334467D2 + & 
                        D1*(-0.625268629D1*1.D-1 + & 
                   D1*(0.5D0*0.591429323D0*1.D-2 + &
                         D1*(0.223751841D0*1.D-3/6.D0 + & 
                         D1*(0.D0         *1.D-3/12.D0 + &
                          D1*0.D0         *1.D-3/20.D0))))
      Wnca(1,1,1)=          -0.998416390D1 + & 
                         D1*(0.348821776D0  *1.D-1 + & 
                 D1*(0.5D0*(-0.119526172D-1)*1.D-2 + &
                         D1*(0.909425662D-2 *1.D-3/6.D0 + & 
                         D1*(0.149166944D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
      Unca(1,1,1)=          -0.143238371D1 + & 
                        D1*(-0.201636224D0  *1.D-1 + & 
                 D1*(0.5D0*(-0.443804232D-1)*1.D-2 + &
                         D1*(0.641847819D-2 *1.D-3/6.D0 + & 
                         D1*(0.296327801D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
      Vnca(1,1,1)=          -0.207474566D1 +  &
                         D1*(0.594737744D0  *1.D-1 +  &
                   D1*(0.5D0*0.674052221D-1 *1.D-2 + &
                         D1*(0.D0           *1.D-3/6.D0 + & 
                        D1*(-0.394845016D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
!
!
!  ** H - SO4 **
!
      A1ca(1,2)=9.5D0
      B1ca(1,2)=            -0.982408701D2 + &
                        D1*(-0.205401806D2  *1.D-1 +  &
                 D1*(0.5D0*(-0.207137292D1) *1.D-2 + &
                        D1*(-0.376521937D-1 *1.D-3/6.D0 +  &
                        D1*(-0.139689758D-1 *1.D-3/12.D0 + &
                         D1* 0.D0           *1.D-3/20.D0))))
      Wnca(1,1,2)=          -0.107752155D2   + & 
                        D1*(-0.879298257D0  *1.D-1 + & 
                 D1*(0.5D0*(-0.440528485D0) *1.D-2 + &
                        D1*(-0.544913927D-1 *1.D-3/6.D0 + & 
                        D1*(-0.173541364D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
      Unca(1,1,2)=          -0.133603464D2 + & 
                        D1*(-0.459479578D1  *1.D-1 + & 
                 D1*(0.5D0*(-0.146220346D1) *1.D-2 + &
                        D1*(-0.157872023D0  *1.D-3/6.D0 + & 
                        D1*(-0.162230945D-3 *1.D-3/12.D0 + &
                         D1* 0.D0           *1.D-3/20.D0))))
      Vnca(1,1,2)=           0.310121997D1    + & 
                         D1*(0.446189009D1  *1.D-1 + & 
                   D1*(0.5D0*0.975254718D0  *1.D-2 + &
                         D1*(0.588748231D-2 *1.D-3/6.D0 + & 
                        D1*(-0.901983372D-3 *1.D-3/12.D0 + &
                          D1*0.D0           *1.D-3/20.D0))))
!
!
!.....HNO3...
!     Parameters from limited 25oC only fit:
!
      A1ca(1,3)=17.D0
      B1ca(1,3)=13.53417796D0
      Wnca(1,1,3)=-3.071864721D0
      Unca(1,1,3)=1.965818001D0 
      Vnca(1,1,3)=-1.411912043D0
!
!.....H-Cl...
!     Parameters from limited 25oC only fit:
!
      A1ca(1,4)=13.D0
      B1ca(1,4)=17.5347093D0
      Wnca(1,1,4)=-14.9654933D0
      Unca(1,1,4)=-13.7294155D0
      Vnca(1,1,4)=3.20778857D0
!
!.....HSO4-NO3-H interactions:
!     I assume these are the same as determined in the Carslaw et al.
!     paper:
!
      Waac(1,3,1)=-4.280D0
      Uaac(1,3,1)=0.201362D0 + 0.084830D0*(T-273.15D0)
      Waac(3,1,1)=Waac(1,3,1)
      Uaac(3,1,1)=-Uaac(1,3,1)
!
!.....SO4-NO3-H interactions:
      Waac(2,3,1)=-0.033291D0*(T-273.15D0)
      Waac(3,2,1)=Waac(2,3,1)
!
!
! ====================================================================
! ====================================================================
!    SECTION (2), H-NH4-HSO4-SO4.
!
!..NH4-HSO4..(22 Feb [5] result):
!
      A1ca(2,1)=19.D0
      B1ca(2,1)=14.2261681D0
      Wnca(1,2,1)=-2.56359462D0
      Unca(1,2,1)=-0.796273529D0
      Vnca(1,2,1)=0.663584552D0
!
!
!.....NH4-SO4..(7 Feb [5] result):
!
      A1ca(2,2)=13.D0
      B1ca(2,2)= -2.858988D0
      A2ca(2,2)=1.5D0
      B2ca(2,2)= 0.D0
      Wnca(1,2,2)=-0.7401490D0
      Unca(1,2,2)=0.9408600D0
      Vnca(1,2,2)=-2.587430D0
!
!
!.....NH4-NO3...(14 Feb [3] result):
! 
      A1ca(2,3)=7.D0
      B1ca(2,3)=24.7529D0
      A2ca(2,3)=13.D0
      B2ca(2,3)=-29.9961D0
      Wnca(1,2,3)=0.900729D0
      Unca(1,2,3)=0.379736D0
      Vnca(1,2,3)=-1.42646D0
!
!
!   ..Now the mixture parameters:
!
!     ..HSO4-SO4-NH4 interactions (22 Feb [5]):
        Waac(1,2,2)=-14.753D0
        Qnaac(1,1,2,2)=4.7204D0
        Uaac(1,2,2)=-16.317D0
        Waac(2,1,2)=Waac(1,2,2)
        Qnaac(1,2,1,2)=Qnaac(1,1,2,2)
        Uaac(2,1,2)=-Uaac(1,2,2)
!
!     ..H-NH4-HSO4 interactions (22 Feb [5]):
        Wcca(1,2,1)=-19.494D0 
        Qncca(1,1,2,1)=8.7607D0 
        Wcca(2,1,1)=Wcca(1,2,1)
        Qncca(1,2,1,1)=Qncca(1,1,2,1)

!     ..H-NH4-SO4 interactions (22 Feb [5]):
        Wcca(1,2,2)=-4.3507D0
        Ucca(1,2,2)=6.5216D0
        Wcca(2,1,2)=Wcca(1,2,2)
        Ucca(2,1,2)=-Ucca(1,2,2)
!
!     ..H-NH4-NO3 interactions (18 Feb [13]):
        Wcca(1,2,3)=-3.0708D0
        Qncca(1,1,2,3)=0.28491D0
        Ucca(1,2,3)=-0.46338D0
        Wcca(2,1,3)=Wcca(1,2,3)
        Qncca(1,2,1,3)=Qncca(1,1,2,3)
        Ucca(2,1,3)=-Ucca(1,2,3)
!
!     ..SO4-NO3-NH4 interactions (17 Feb [4]):
        Qnaac(1,2,3,2)=2.9795D0
        Uaac(2,3,2)=-1.2163D0
        Qnaac(1,3,2,2)=Qnaac(1,2,3,2)
        Uaac(3,2,2)=-Uaac(2,3,2)
!
!     ..HSO4-NO3-NH4 interactions (17 Feb [4]):
!       **updated 16 Aug **
        Waac(1,3,2)=-2.9369D0
        Waac(3,1,2)=Waac(1,3,2)
!
!
! ===================================================================
! ===================================================================
!    SECTION (3), parameters for Na-Cl-NO3-SO4 Mixtures (25oC)
!
!.....Na-Cl (4 Feb [5], "limit1"):
!
      A1ca(3,4)=5.D0
      B1ca(3,4)=19.93376D0
      Wnca(1,3,4)=-5.646077D0
      Unca(1,3,4)=-3.609246D0
      Vnca(1,3,4)=-2.459821D0
!
!.....Na-SO4 (5 Feb [2], "limit"):
!
      A1ca(3,2)=8.D0
      B1ca(3,2)=34.46602D0
      Wnca(1,3,2)=-3.725962D0
      Unca(1,3,2)=-1.959160D0
      Vnca(1,3,2)=-4.860570D0
!
!.....Na-NO3 (10 Feb [4]):
!
      A1ca(3,3)=5.D0
      B1ca(3,3)=26.99939D0
      A2ca(3,3)=13.D0
      B2ca(3,3)=-21.60500D0
      Wnca(1,3,3)=0.5269081D0
      Unca(1,3,3)=0.2666436D0
      Vnca(1,3,3)=-2.302876D0
!
!  ...SO4-Cl-Na (10 Feb [2]):
      Waac(2,4,3)=4.827D0
      Waac(4,2,3)=Waac(2,4,3)
      Qnaac(1,2,4,3)=0.05163D0
      Qnaac(1,4,2,3)=Qnaac(1,2,4,3)
!
!  ...NO3-Cl-Na (11 Feb [1]):
      Waac(3,4,3)=-6.923D0
      Waac(4,3,3)=Waac(3,4,3)
      Qnaac(1,3,4,3)=4.181D0
      Qnaac(1,4,3,3)=Qnaac(1,3,4,3)
!
!  ...SO4-NO3-Na (11 Feb [3]):
      Waac(2,3,3)=-9.498D0
      Waac(3,2,3)=Waac(2,3,3)
      Qnaac(1,2,3,3)=8.528D0
      Qnaac(1,3,2,3)=Qnaac(1,2,3,3)
!
! ===================================================================
! ===================================================================
!     SECTION (4), other mixtures (work in progress):
!
!  ...Na-HSO4 (20 Feb [9]):
      A1ca(3,1)=19.D0
      B1ca(3,1)= 62.27961D0
      Wnca(1,3,1)= -2.932425D0
      Unca(1,3,1)= -4.857197D0
      Vnca(1,3,1)=  4.888311D0
!
!  ...NH4-Cl (4 Feb [1]):
      A1ca(2,4)=15.D0
      B1ca(2,4)= 4.659688D0
      Wnca(1,2,4)= -0.5682911D0
      Unca(1,2,4)= 2.072437D0
      Vnca(1,2,4)= -1.25D0
!
!  ...H-Na-HSO4 (20 Feb [9]):
      Wcca(1,3,1)=-8.96894D0
      Qncca(1,1,3,1)=4.16202D0
      Ucca(1,3,1)=-2.92819D0
      Wcca(3,1,1)=Wcca(1,3,1)
      Qncca(1,3,1,1)=Qncca(1,1,3,1)
      Ucca(3,1,1)=-Ucca(1,3,1)
!
!  ...H-Na-SO4 (20 Feb [9]):
      Wcca(1,3,2)=15.9075D0
      Wcca(3,1,2)=Wcca(1,3,2)
      Qncca(1,1,3,2)=-8.82425D0
      Qncca(1,3,1,2)=Qncca(1,1,3,2)
!
!  ...HSO4-SO4-Na (20 Feb [9]):
      Qnaac(1,1,2,3)=-4.68641D0
      Qnaac(1,2,1,3)=Qnaac(1,1,2,3)
!
!  ...H-Na-Cl (12 Feb [6]):
      Wcca(1,3,4)=2.2490D0
      Wcca(3,1,4)=Wcca(1,3,4)
      Qncca(1,1,3,4)=-0.25080D0
      Qncca(1,3,1,4)=Qncca(1,1,3,4)
!
!  ...H-Na-NO3 (18 Feb [3]):
      Wcca(1,3,3)=0.46039D0
      Wcca(3,1,3)=Wcca(1,3,3)
      Ucca(1,3,3)=1.1749D0
      Ucca(3,1,3)=-Ucca(1,3,3)
!
!  ...H-NH4-Cl (12 Feb [8]):
      Wcca(1,2,4)=-19.977D0
      Wcca(2,1,4)=Wcca(1,2,4)
      Qncca(1,1,2,4)=10.233D0
      Qncca(1,2,1,4)=Qncca(1,1,2,4)
!
!  ...NH4-Cl-NO3 (14 Feb [6]):
      Waac(4,3,2)=-0.2207D0
      Waac(3,4,2)=Waac(4,3,2)
      Qnaac(1,4,3,2)=-0.1173D0
      Qnaac(1,3,4,2)=Qnaac(1,4,3,2)
!
!  ...NH4-Na-NO3 (15 Feb [2], .RS2 result):
      Wcca(2,3,3)=-0.35411D0
      Wcca(3,2,3)=Wcca(2,3,3)
      Qncca(1,2,3,3)=0.046254D0
      Qncca(1,3,2,3)=Qncca(1,2,3,3)
      Ucca(2,3,3)=0.2130D0
      Ucca(3,2,3)=-Ucca(2,3,3)
!
!  ...NH4-Na-Cl (13 Feb [2]):
      Wcca(2,3,4)=-5.6414D0
      Wcca(3,2,4)=Wcca(2,3,4)
      Qncca(1,2,3,4)=3.2919D0
      Qncca(1,3,2,4)=Qncca(1,2,3,4)
!
!  ...NH4-Na-SO4 (13 Feb [3]):
      Wcca(2,3,2)=-1.4832D0
      Wcca(3,2,2)=Wcca(2,3,2)
      Qncca(1,2,3,2)=0.76211D0
      Qncca(1,3,2,2)=Qncca(1,2,3,2)
!
!  ...NH4-Cl-SO4 (13 Feb [5]):
      Uaac(2,4,2)=-1.0709D0
      Uaac(4,2,2)=-Uaac(2,4,2)
      Qnaac(1,4,2,2)=1.0869D0 
      Qnaac(1,2,4,2)=Qnaac(1,4,2,2)
!
      RETURN
END SUBROUTINE SetCoefficient_25

SUBROUTINE ComputeActivityAim(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

  INTEGER :: i,iFrac,ix,iy,iz,is
  REAL(RealKind) :: xCat(nC),xAn(nA),xNeut(nN)
  REAL(RealKind) :: ActNeut(nN),ActCat(nC),ActAn(nA)

  REAL(RealKind) :: TempLoc,Excess
  REAL(RealKind) :: MolGes,MolWater,Fak,SumX

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        TempLoc=TAbs%c(ix,iy,iz,1) 
        CALL SetCoefficient_25(TempLoc)
        DO iFrac=1,nFrac
          IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
            IF (cAqua(iWater)%c(ix,iy,iz,iFrac)<Zero) THEN
               WRITE(*,*) 'Negativ',cAqua(iWater)%c(ix,iy,iz,iFrac),iFrac,ix,iy,iz
               STOP
            END IF
            SumX=Zero
            xCat=Zero
            xAn=Zero
            xNeut=Zero
            DO i=1,nAqua
              IF (AimInd(i)>0) THEN
                xCat(AimInd(i))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
                SumX=SumX+xCat(AimInd(i))
              ELSE IF (AimInd(i)<0) THEN
                xAn(ABS(AimInd(i)))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
                SumX=SumX+xAn(ABS(AimInd(i)))
              END IF
            END DO
            xNeut(1)=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
            SumX=SumX+xNeut(1)
            xCat=xCat/SumX
            xAn=xAn/SumX
            xNeut=xNeut/SumX
            CALL ExcessGibbs(xCAT,xAN,xNEUT,TempLoc &
                            ,ActCat,ActAn,ActNeut &
                            ,Excess)
            MolWater=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
            MolGes=0.0d0
            DO i=iWater+1,nSoluble
              MolGes=MolGes &
                    +ABS(cAqua(i)%c(ix,iy,iz,iFrac))/MolMass(i)
            END DO
            Fak=One/(MolGes/MolWater+One)
            DO i=1,nAqua
              IF (AimInd(i)>0) THEN
                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActCat(AimInd(i))*Fak
              ELSE IF (AimInd(i)<0) THEN
                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActAn(ABS(AimInd(i)))*Fak
              END IF
            END DO
            ActCoeff(iWater)%c(ix,iy,iz,iFrac)=ActNeut(1)
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeActivityAim

SUBROUTINE OutputActivityAim(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

  INTEGER :: i,iFrac,ix,iy,iz,is
  REAL(RealKind) :: xCat(nC),xAn(nA),xNeut(nN)
  REAL(RealKind) :: ActNeut(nN),ActCat(nC),ActAn(nA)

  REAL(RealKind) :: TempLoc,Excess
  REAL(RealKind) :: Fak,SumX

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        TempLoc=TAbs%c(ix,iy,iz,1) ! OSSI Umrechnen
        CALL SetCoefficient_25(TempLoc)
        DO iFrac=1,nFrac
          IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
            SumX=Zero
            xCat=Zero
            xAn=Zero
            xNeut=Zero
            DO i=1,nAqua
              IF (AimInd(i)>0) THEN
                xCat(AimInd(i))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
                SumX=SumX+xCat(AimInd(i))
              ELSE IF (AimInd(i)<0) THEN
                xAn(ABS(AimInd(i)))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
                SumX=SumX+xAn(ABS(AimInd(i)))
              END IF
            END DO
            xNeut(1)=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
            SumX=SumX+xNeut(1)
            xCat=xCat/SumX
            xAn=xAn/SumX
            xNeut=xNeut/SumX
            CALL ExcessGibbs(xCAT,xAN,xNEUT,TempLoc &
                            ,ActCat,ActAn,ActNeut &
                            ,Excess)
            DO i=3,nAqua
              IF (AimInd(i)>0) THEN
                WRITE(*,*) SpeciesName(i),xCat(AimInd(i)),ActCat(AimInd(i))
              ELSE IF (AimInd(i)<0) THEN
                WRITE(*,*) SpeciesName(i),xAn(-AimInd(i)),ActAn(-AimInd(i))
              END IF
            END DO
            WRITE(*,*) xNeut,ActNeut
            DO i=3,nAqua
              WRITE(*,*) SpeciesName(i),cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)/cAqua(iWater)%c(ix,iy,iz,iFrac), &
                         ActCoeff(i)%c(ix,iy,iz,iFrac)
            END DO
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE OutputActivityAim

END MODULE ActivityAim_Mod


