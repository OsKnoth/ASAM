! Waermekapazitaet von Eis nach Pruppacher und Klett (1997) S.86

DOUBLE PRECISION FUNCTION ceisxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT
  
  ceisxd = 4186.8D0*(0.503D0 + 0.00175D0*(TT - T0))
  
  RETURN
END FUNCTION ceisxd
