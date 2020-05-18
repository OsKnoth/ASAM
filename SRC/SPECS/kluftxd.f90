! Waermeleitfaehigkeit von Luft nach Pruppacher und Klett (1997) S.508

DOUBLE PRECISION FUNCTION kluftxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT,Tx
  
  kluftxd = 4.1868D+02*(5.69D0 + 0.017D0*(TT - T0))*1.0D-05
  
  RETURN
END FUNCTION kluftxd
