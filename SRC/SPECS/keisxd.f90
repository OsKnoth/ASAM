! Waermeleitfaehigkeit von Eis nach Pruppacher und Klett (1997) S.676

DOUBLE PRECISION FUNCTION keisxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT
  
  keisxd = 418.68D0*(537.4D0 - 1.48D0*(TT - T0) + 0.028*(TT - T0)**2.0D0)   &
                          *1.0D-05
  
  RETURN
END FUNCTION keisxd
