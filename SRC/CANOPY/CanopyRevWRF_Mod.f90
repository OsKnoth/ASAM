MODULE CanopyRevWRF_Mod
  USE Kind_Mod
  USE Control_Mod

  IMPLICIT NONE

  CONTAINS

    FUNCTION Phi(Typ,Chi,z,z0,z0h,RiB)
      REAL(RealKind) :: Phi
      REAL(RealKind) :: RiB
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) :: Chi
      REAL(RealKind) :: a,b
      CHARACTER(1) :: Typ
      IF (Typ.EQ.'M') THEN
        a = -6.1d0
        b =  2.5d0
      ELSE IF (Typ.EQ.'H') THEN
        a = -5.3d0
        b =  1.1d0
      END IF
      IF (RiB>=Zero) THEN
        Phi  = a*LOG(Chi+(One+Chi**b)**(One/b))
      ELSE
        Phi  = (PhiK(Typ,Chi)+Chi**Two*Phi_C(Typ,Chi))/(One+Chi**Two)
      END IF
      Phi = MIN(0.9d0*LOG((z+z0h)/z0),Phi)
    END FUNCTION Phi

    FUNCTION PhiK(Typ,Chi)
      REAL(RealKind) :: PhiK
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var1
      REAL(RealKind) :: alpha
      CHARACTER(1) :: Typ
      alpha= 16.d0
      Var1  = (One-alpha*Chi)**(One/Four)
      IF (Typ.EQ.'H') THEN
        PhiK = Two*LOG((One+Var1**Two)/Two)
      ELSE IF (Typ.EQ.'M') THEN
        PhiK = Two*LOG((One+Var1)/Two)   &
              +LOG((One+Var1**Two)/Two)  &
              -Two*ATAN(Var1)+PiHalf
      END IF
    END FUNCTION PhiK

    FUNCTION Phi_C(Typ,Chi)
      REAL(RealKind) :: Phi_C
      REAL(RealKind) :: Chi
      REAL(RealKind) :: Var2
      REAL(RealKind) :: beta
      CHARACTER(1)   :: Typ
      IF (Typ.EQ.'M') THEN
        beta = 10.d0
      ELSE IF (Typ.EQ.'H') THEN
        beta = 34.d0
      END IF
      Var2 = (One-beta*Chi)**(One/Three)
      Phi_C = Three/Two*LOG((Var2**Two+Var2+One)/Three)    &
             -SQRT(Three)*ATAN((Two*Var2+One)/SQRT(Three)) &
             +Pi/SQRT(Three)
    END FUNCTION Phi_C

    FUNCTION dPhidChi(Typ,Chi,z0,z0h,z,RiB)
      REAL(RealKind) :: dPhidChi
      REAL(RealKind) :: RiB
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) :: Chi,Chi1,Chi2
      CHARACTER(1)   :: Typ
      Chi1   = Chi+Eps*z
      Chi2   = Chi-Eps*z
      IF (RiB<Zero) THEN
        dPhidChi=( ((PhiK(Typ,Chi1)-PhiK(Typ,Chi2))/(Two*Eps*z) &
                    +(Two*Chi*Phi_C(Typ,Chi)+Chi**Two*(Phi_C(Typ,Chi1)-Phi_C(Typ,Chi2))/(Two*Eps*z))) &
                  *(One+Chi**Two) &
                  -(PhiK(Typ,Chi)+Chi**Two*Phi_C(Typ,Chi))*Two*Chi) &
                 /((One+Chi**Two)**Two)
      ELSE
        dPhidChi=(Phi(Typ,Chi1,z,z0,z0h,RiB)-Phi(Typ,Chi2,z,z0,z0h,RiB))/(Two*Eps*z)
      END IF
    END FUNCTION dPhidChi

    FUNCTION RiBulk(z,z0,z0h,Chi,RiB)
      REAL(RealKind) :: RiBulk
      REAL(RealKind) :: RiB
      REAL(RealKind) :: z,z0,z0h,Chi
      REAL(RealKind) :: PhiHeat,PhiMom
      PhiHeat= LOG((z+z0)/z0h)                   &
              -Phi('H',Chi*(One+z0/z),z,z0,z0h,RiB) &
              +Phi('H',Chi*z0h/z,z,z0,z0h,RiB)
      PhiMom = LOG((z+z0)/z0)                    &
              -Phi('M',Chi*(One+z0/z),z,z0,z0h,RiB) &
              +Phi('M',Chi*z0/z,z,z0,z0h,RiB)
      RiBulk = Chi*PhiHeat/PhiMom**Two
    END FUNCTION

    FUNCTION dRiBdChi(z,z0,z0h,Chi,RiB)
      REAL(RealKind) :: dRiBdChi
      REAL(RealKind) :: RiB
      REAL(RealKind) :: z,z0,z0h
      REAL(RealKind) :: Chi
      REAL(RealKind) :: PhiHeat,PhiMom
      PhiHeat= LOG((z+z0)/z0h)                   &
              -Phi('H',Chi*(One+z0/z),z,z0,z0h,RiB) &
              +Phi('H',Chi*z0h/z,z,z0,z0h,RiB)
      PhiMom = LOG((z+z0)/z0)                    &
              -Phi('M',Chi*(One+z0/z),z,z0,z0h,RiB) &
              +Phi('M',Chi*z0/z,z,z0,z0h,RiB)
      dRiBdChi = RiBulk(z,z0,z0h,Chi,RiB)/Chi &
                +Chi*( (-dPhidChi('H',Chi*(One+z0/z),z0,z0h,z,RiB) &
                         +dPhidChi('H',Chi*(z0h/z),z0,z0h,z,RiB)) &
                       *PhiMom**Two &
                      -Two*PhiMom &
                          *(-dPhidChi('M',Chi*(One+z0/z),z0,z0h,z,RiB) &
                            +dPhidChi('M',Chi*(z0/z),z0,z0h,z,RiB)) &
                       *PhiHeat ) &
                     /PhiMom**Four
    END FUNCTION

    SUBROUTINE NewtonMethod(x_n,x_n1,fx,dfdx)
      REAL(RealKind) :: x_n
      REAL(RealKind) :: x_n1
      REAL(RealKind) :: fx,dfdx
      x_n1=x_n-fx/dfdx
!      WRITE(*,*) 'x_n1',x_n1,'x_n',x_n,fx,dfdx
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


END MODULE CanopyRevWRF_Mod
