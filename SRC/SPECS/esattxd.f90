! Berechnung des Saettigungsdampfdruckes ueber Wasser

DOUBLE PRECISION FUNCTION esattxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT
  
  
  ! esattxd in [hPa]
  !!esattxd = 6.1D0*DEXP(17.27D0*(TT - T0)/(TT - 35.86D0))
  ! esattxd in [Pa]
  !!esattxd = esattxd*100.0D0
  
  ! CMS consistent version
  esattxd = 6.112D2 * EXP(17.67D0*(TT - T0)/(TT - 29.65D0))
  
  RETURN
END FUNCTION esattxd
