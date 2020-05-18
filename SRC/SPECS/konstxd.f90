subroutine konstxd   

! Subroutine for the definition of all globally used constants
! handoff is in Common-Blocks, that are merged in konst.com 
!
! called by wmain.f
!
! calling none

INCLUDE 'HEADER3'

implicit none

INTEGER :: i


! physical constants
T0    = 273.16D0
CP    = 1005.D0
GRAV  = 9.81D0
KTH   = 2.43D-2
LS    = 2.834D6
LV    = 2.50078D6
LV_ICE = LS - LV 
mol_w = 18.D-3
PI    = 2.D0*ASIN(1.D0)
RDRY  = 287.D0
RW    = 461.5D0
RHOW  = 1000.D0
RHOI  =  900.D0
SIGMA0= 7.28D-2
E0    = 6.112D2
C1    = 17.67D0
C2    =-29.65D0
! AP chemical parameters
inacl = 0
! NaCl
if(inacl.eq.1) then
  bb = 1.64D-3        ! NaCl
  mol_AS = 58.44D-3
  rho_AS = 2165.D0
  vantS = 2.D0
  mol_g(1) = 0.0D0
  phi_g(1) = 1.D0
  mol_g(2) = 0.1D0
  phi_g(2) = 0.9324D0
  mol_g(3) = 0.5D0
  phi_g(3) = 0.9209D0
  mol_g(4) = 1.0D0
  phi_g(4) = 0.9355D0
  mol_g(5) = 2.0D0
  phi_g(5) = 0.9833D0
  mol_g(6) = 5.0D0
  phi_g(6) = 1.1916D0
  mol_g(7) = 12.0D0
  phi_g(7) = 1.55D0
endif
! (NH_4)_2 SO_4
if(inacl.eq.0) then
  bb = 0.67D-3        ! AmmS
  mol_AS = 132.13D-3
  rho_AS = 1770.D0
  vantS = 3.D0
  mol_g(1) = 0.0D0
  phi_g(1) = 1.D0
  mol_g(2) = 0.1D0
  phi_g(2) = 0.767D0
  mol_g(3) = 0.5D0
  phi_g(3) = 0.677D0
  mol_g(4) = 1.0D0
  phi_g(4) = 0.640D0
  mol_g(5) = 2.0D0
  phi_g(5) = 0.623D0
  mol_g(6) = 5.0D0
  phi_g(6) = 0.672D0
  mol_g(7) = 12.0D0
  phi_g(7) = 1.D0
endif
! average biomass burning aerosol with soluble part   
rho_AI = 1500.D0
rho_AI = rho_AS
mol_AI = 100.D-3

CAPPA = RDRY / CP
FAC34 = 3.D0/(4.D0*PI)
FACT  = RHOW/FAC34
FACTI = RHOI/FAC34
fact_AS = rho_AS/FAC34
fact_AI = rho_AI/FAC34
! entrainment
GAMMA=0.5D0
RAD_BUB0=2.D3
! for terminal velocity of drops
k1_vel=1.19D8
k2_vel=2.20D2
rho_0=1.20D0
! empirical parameters for Breakup
PC1  = 2.94D-7
PC2  = 3400.D0
QC1  = 145.37D0
QC2  = 7.D0
! smallest size
R0=1.D-9
!      R0 = R0 * 2.D0**(.25D0/((NMAX-1.D0)*3.D0))
! parameters for AP distributions in apvert.f


RAmin=1.D-9
RAmax=1.D-5
! VG 
!PR=1.001D0
PR=1.005D0
! parameters for surface tension
TC=647.069D0
mu=1.256D0
b1=235.8D-3
b2=-0.625D0
! parameters for diffusion and heat conduction
ALPHA_C=0.036D0
ALPHA_C=0.042D0
ALPHA_C=1.0D0
ALPHA_T=0.7D0
DEL_V=1.04D-7
DEL_T=2.16D-7
! weighting factor for LDM in cond_shift.f
! WEI_Q = 1 Verschiebung der Aerosolmasse proportional zur Wassermasse
! WEI_Q = 0 Verschiebung der Aerosolmasse proportional zur Anzahl
! Fuer 1d-Version
! WEI_Q = 0.5D0 
! Fuer 2d-Version
! WEI_Q = 0.0D0
! VG
!WEI_Q = 0.5D0
WEI_Q = 1.0D0
WEI_N = 1.D0 - WEI_Q
! stop criterion for iteration loops
SC_cond=1.D-6
SC_ggw =1.D-6
SATTabs=1.D-10
size_fact=1.00001D0
accu_fact=1.0001D0     ! immer deutlich groesser als size_fact
step_fact=2.D0
ITERMAX=500
! Numbers
!SMALL1= 1.D-20
SMALL1= 1.D0
QU1D3=1.D0/3.D0
QU4D3=4.D0/3.D0
QU5D3=5.D0/3.D0

! Calculations of some parameters
! piecewise linear function for osmotic coefficient
DO i=1,6
  a_phi(i) = (phi_g(i+1) - phi_g(i))/(mol_g(i+1) - mol_g(i))
  b_phi(i) = phi_g(i) - a_phi(i)*mol_g(i)
END DO
a_phi(7) = 0.D0
b_phi(7) = phi_g(7)
      

return
end
