! Waermekapazitaet von Wasser nach Pruppacher und Klett (1997) S.93
! ACHTUNG! Gilt nur fuer TT <= T0.

DOUBLE PRECISION FUNCTION cwasserxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT
  
  
  cwasserxd = 4186.8D0*(1.000938D0 - 2.7052D-3*(TT - T0)                    &
              - 2.3235D-5*(TT - T0)**2.0D0 + 4.3778D-6*(TT - T0)**3.0D0     &
              + 2.7136D-7*(TT - T0)**4.0D0)
  
  RETURN
END FUNCTION cwasserxd
