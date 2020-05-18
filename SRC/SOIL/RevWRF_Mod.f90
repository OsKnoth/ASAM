MODULE RevWRF_Mod
  USE Kind_Mod
  USE Control_Mod

  IMPLICIT NONE

  CONTAINS

    FUNCTION PsiH(Chi,z,z0,z0h,RiB,LandClass)
      REAL(RealKind) :: PsiH
      REAL(RealKind) :: Chi
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) :: RiB
      INTEGER :: LandClass

      REAL(RealKind) :: a,b

      IF (LandClass<=8) THEN
        a = -5.3d0
        b =  1.1d0
          IF (RiB>=Zero) THEN
            PsiH  = a*LOG(Chi+(One+Chi**b)**(One/b))
          ELSE
            PsiH  = (PsiHK(Chi)+Chi**Two*PsiH_C(Chi))/(One+Chi**Two)
          END IF
        PsiH = MIN(0.9d0*LOG((z+z0h)/z0),PsiH)
      ELSE IF (LandClass==10) THEN
        PsiH = PsiH_SHEBA(Chi)
      ELSE
        STOP 'Falsche Landklasse'
      END IF  
    END FUNCTION PsiH

    FUNCTION PsiM(Chi,z,z0,z0h,RiB,LandClass)
      REAL(RealKind) :: PsiM
      REAL(RealKind) :: RiB
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) :: Chi
      INTEGER :: LandClass

      REAL(RealKind) :: a,b

      IF (LandClass<=8) THEN
        a = -6.1d0
        b =  2.5d0
          IF (RiB>=Zero) THEN
            PsiM  = a*LOG(Chi+(One+Chi**b)**(One/b))
          ELSE
            PsiM  = (PsiMK(Chi)+Chi**Two*PsiM_C(Chi))/(One+Chi**Two)
        END IF
        PsiM = MIN(0.9d0*LOG((z+z0h)/z0),PsiM)
      ELSE IF (LandClass==10) THEN
        PsiM = PsiM_SHEBA(Chi)
      ELSE
        STOP 'Falsche Landklasse'
      END IF  
    END FUNCTION PsiM

    FUNCTION PsiHK(Chi)
      REAL(RealKind) :: PsiHK
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var1
      REAL(RealKind) :: alpha
      
      alpha = 16.d0
      Var1  = (One-alpha*Chi)**(One/Four)
      PsiHK = Two*LOG((One+Var1**Two)/Two)
    END FUNCTION PsiHK

    FUNCTION PsiMK(Chi)
      REAL(RealKind) :: PsiMK
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var1
      REAL(RealKind) :: alpha

      alpha = 16.d0
      Var1  = (One-alpha*Chi)**(One/Four)
      PsiMK = Two*LOG((One+Var1)/Two)   &
              +LOG((One+Var1**Two)/Two)  &
              -Two*ATAN(Var1)+PiHalf
    END FUNCTION PsiMK

    FUNCTION PsiH_C(Chi)
      REAL(RealKind) :: PsiH_C
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var2
      REAL(RealKind) :: beta
      
      beta = 34.d0
      Var2 = (One-beta*Chi)**(One/Three)
      PsiH_C = Three/Two*LOG((Var2**Two+Var2+One)/Three)    &
             -SQRT(Three)*ATAN((Two*Var2+One)/SQRT(Three)) &
             +Pi/SQRT(Three)
    END FUNCTION PsiH_C

    FUNCTION PsiM_C(Chi)
      REAL(RealKind) :: PsiM_C
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var2
      REAL(RealKind) :: beta
       
      beta = 10.d0
      Var2 = (One-beta*Chi)**(One/Three)
      PsiM_C = Three/Two*LOG((Var2**Two+Var2+One)/Three)    &
             -SQRT(Three)*ATAN((Two*Var2+One)/SQRT(Three)) &
             +Pi/SQRT(Three)
    END FUNCTION PsiM_C


    ! formula from SHEBA flux–profile relationships in the stable
    ! atmospheric boundary layer by Gravech et al. 2007
    FUNCTION PsiH_SHEBA(Chi)
      REAL(RealKind) :: PsiH_SHEBA
      REAL(RealKind) :: Chi
      REAL(RealKind) :: grbh
      REAL(RealKind) :: ah
      REAL(RealKind) :: bh
      REAL(RealKind) :: ch

      ah = 5.d0
      bh = 5.d0
      ch = 3.d0
      grbh = SQRT(ch**Two-Four)

      PsiH_SHEBA = - bh/Two * LOG(One+ch*Chi+Chi)  &
                   + (-ah/grbh + bh*ch/(Two*grbh))  &
                   * (LOG((Two*Chi+ch-grbh)/(Two*Chi+ch+grbh))-  &
                      LOG((ch-grbh)/(ch+grbh)))

    END FUNCTION PsiH_SHEBA


    FUNCTION PsiM_SHEBA(Chi)
      REAL(RealKind) :: PsiM_SHEBA
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var3
      REAL(RealKind) :: grbm
      REAL(RealKind) :: bm
      REAL(RealKind) :: am

      am = 5.d0
      bm = am/6.5d0
      grbm = ((One-bm)/bm)**(One/Three)
      Var3 = (One+Chi)**(One/Three)

      PsiM_SHEBA = -Three*am/bm*(Var3-One)+(am*grbm)/(Two*bm)   &
                   *(Two*LOG((Var3+grbm)/(One+grbm))  &
                   -LOG(Var3**Two-Var3*grbm+grbm**Two)/(One-grbm+grbm**Two)  &
                   +Two*SQRT(Three)*(ATAN((Two*Var3-grbm)/(SQRT(Three)*grbm))  &
                   -ATAN((Two-grbm)/(SQRT(Three)*grbm))))

    END FUNCTION PsiM_SHEBA  
     


!   FUNCTION dPsidChi(Typ,Chi,z0,z0h,z,RiB)
!     REAL(RealKind) :: dPsidChi
!     REAL(RealKind) :: RiB
!     REAL(RealKind) :: z,z0,z0h
!     REAL(RealKind) :: Chi,Chi1,Chi2
!     CHARACTER(1)   :: Typ
!     Chi1   = Chi+Eps*z
!     Chi2   = Chi-Eps*z
!     IF (RiB<Zero) THEN
!       dPsidChi=( ((PsiK(Typ,Chi1)-PsiK(Typ,Chi2))/(Two*Eps*z) &
!                   +(Two*Chi*Psi_C(Typ,Chi)+Chi**Two*(Psi_C(Typ,Chi1)-Psi_C(Typ,Chi2))/(Two*Eps*z))) &
!                 *(One+Chi**Two) &
!                 -(PsiK(Typ,Chi)+Chi**Two*Psi_C(Typ,Chi))*Two*Chi) &
!                /((One+Chi**Two)**Two)
!     ELSE
!       dPsidChi=(Psi(Typ,Chi1,z,z0,z0h,RiB)-Psi(Typ,Chi2,z,z0,z0h,RiB))/(Two*Eps*z)
!     END IF
!   END FUNCTION dPsidChi

    FUNCTION RiBulk(z,z0,z0h,Chi,RiB)
      REAL(RealKind) :: RiBulk
      REAL(RealKind) :: z,z0,z0h,Chi
      REAL(RealKind) :: RiB
      INTEGER :: LandClass=1

      REAL(RealKind) :: PsiHeat,PsiMom
      PsiHeat= LOG((z+z0)/z0h)                   &
              -PsiH(Chi*(One+z0/z),z,z0,z0h,RiB,LandClass) &
              +PsiH(Chi*z0h/z,z,z0,z0h,RiB,LandClass)
      PsiMom = LOG((z+z0)/z0)                    &
              -PsiM(Chi*(One+z0/z),z,z0,z0h,RiB,LandClass) &
              +PsiM(Chi*z0/z,z,z0,z0h,RiB,LandClass)
      RiBulk = Chi*PsiHeat/PsiMom**Two
    END FUNCTION

!   FUNCTION dRiBdChi(z,z0,z0h,Chi,RiB)
!     REAL(RealKind) :: dRiBdChi
!     REAL(RealKind) :: RiB
!     REAL(RealKind) :: z,z0,z0h
!     REAL(RealKind) :: Chi
!     REAL(RealKind) :: PsiHeat,PsiMom
!     PsiHeat= LOG((z+z0)/z0h)                   &
!             -Psi('H',Chi*(One+z0/z),z,z0,z0h,RiB) &
!             +Psi('H',Chi*z0h/z,z,z0,z0h,RiB)
!     PsiMom = LOG((z+z0)/z0)                    &
!             -Psi('M',Chi*(One+z0/z),z,z0,z0h,RiB) &
!             +Psi('M',Chi*z0/z,z,z0,z0h,RiB)
!     dRiBdChi = RiBulk(z,z0,z0h,Chi,RiB)/Chi &
!               +Chi*( (-dPsidChi('H',Chi*(One+z0/z),z0,z0h,z,RiB) &
!                        +dPsidChi('H',Chi*(z0h/z),z0,z0h,z,RiB)) &
!                      *PsiMom**Two &
!                     -Two*PsiMom &
!                         *(-dPsidChi('M',Chi*(One+z0/z),z0,z0h,z,RiB) &
!                           +dPsidChi('M',Chi*(z0/z),z0,z0h,z,RiB)) &
!                      *PsiHeat ) &
!                    /PsiMom**Four
!   END FUNCTION

    SUBROUTINE NewtonMethod(x_n,x_n1,fx,dfdx)
      REAL(RealKind) :: x_n
      REAL(RealKind) :: x_n1
      REAL(RealKind) :: fx,dfdx
      x_n1=x_n-fx/dfdx
      !WRITE(*,*) 'x_n1',x_n1,'x_n',x_n,fx,dfdx
    END SUBROUTINE

    SUBROUTINE FirstGuess(RiB,Chi,z,z0,z0h)
      REAL(RealKind) :: RiB,Chi
      REAL(RealKind) :: RichCut
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) ,Parameter :: RichCrit=0.25d0

      IF (Zero<RiB.AND.RiB<=0.2d0) THEN
        Chi=LOG(z/z0)*RiB/(One-RiB/RichCrit)
      ELSE IF (RiB>0.2d0) THEN
        RichCut=(LOG(z/z0)+One/RichCrit)**(-One)
        Chi=LOG(z/z0)*RiB/(One-RichCut/RichCrit)
      ELSE
        Chi=(LOG(z/z0)**Two/LOG(z/z0h)-0.55d0)*RiB
      END IF
    END SUBROUTINE FirstGuess


END MODULE RevWRF_Mod
