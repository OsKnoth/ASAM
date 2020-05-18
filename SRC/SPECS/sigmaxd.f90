! Oberflaechenspannung sigma von Wasser in Luft als Funktion der
! Temperatur nach Pruppacher und Klett (1997) S.130.

DOUBLE PRECISION FUNCTION sigmaxd(TT)

USE data_mp_ice, ONLY : T0

IMPLICIT NONE
DOUBLE PRECISION :: TT,Tx

DOUBLE PRECISION :: a0,a1,a2,a3,a4,a5,a6
PARAMETER(a0 = 75.93D0, a1 = 0.115D0, a2 = 6.818D-2)
PARAMETER(a3 = 6.511D-3, a4 = 2.933D-4, a5 = 6.283D-6)
PARAMETER(a6 = 5.285D-8)


Tx = TT - T0
! sigma in [erg/cm**2] = [dyn/cm]
sigmaxd = a0 + a1*Tx + a2*Tx**2.0D0 + a3*Tx**3.0D0 + a4*Tx**4.0D0       &
      + a5*Tx**5.0D0 + a6*Tx**6.0D0
! sigma in [N/m]
sigmaxd = sigmaxd*1.0D-3


RETURN
END
