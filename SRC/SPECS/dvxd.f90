! Berechnung des Diffusionskoeffizienten Dv von Wasserdampf 
!(siehe PK 97, S.503)

DOUBLE PRECISION FUNCTION dvxd(TT,pp)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT,pp
  
  DOUBLE PRECISION :: p0
  PARAMETER(p0 = 101325.0D0)
  
  
  dvxd = 1.0D-4*0.211D0*(p0/pp)*(TT/T0)**1.94D0
  
  RETURN
END FUNCTION dvxd
