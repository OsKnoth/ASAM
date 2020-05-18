! dynamischen Viskositaet nach Pruppacher und Klett (1997) S.417
! in cgs units (poise) instead of SI units (Pa s): 1 poise = 0.1 Pa s
! function is used by the following program parts
! - ventixd.f90: 
! - vstokesxd.f90:

DOUBLE PRECISION FUNCTION etaxd(TT)
  
  USE data_mp_ice, ONLY : T0
  
  IMPLICIT NONE
  DOUBLE PRECISION :: TT,Tx
  
  Tx = TT - T0
  
  ! etaxd in poise
  IF(Tx.GE.0.0D0) etaxd = (1.718 + 0.0049*Tx)*1.0D-4
  IF(Tx.LT.0.0D0) etaxd = (1.718 + 0.0049*Tx - 1.2D-5*Tx**2.0D0)*1.0D-4
  
  ! etaxd in Pa s 
  etaxd = etaxd * 0.1D0
  
  RETURN
END FUNCTION etaxd
