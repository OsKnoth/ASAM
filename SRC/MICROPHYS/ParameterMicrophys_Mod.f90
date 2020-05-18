MODULE ParameterMicrophys_Mod

  USE Thermodynamic_Mod
  USE Parameter_Mod

  IMPLICIT NONE

!-------------------------------------------------------------------------------- 
!Summary of all parameter needed in the module 'TwoMomentBulkMicroPhysics_Mod.f90'
!-------------------------------------------------------------------------------- 
  REAL(RealKind), PARAMETER :: Lvi=3.330d5
  REAL(RealKind), PARAMETER :: k_rr=7.12d0
  REAL(RealKind), PARAMETER :: K_T  = 2.500d-2      ! Conductivity of heat [J m^-1 s^-1 K^-1]
  REAL(RealKind), PARAMETER :: nu_l = 1.460d-5      ! Kinem. Visc. von Luft [m^2 s^-1]
  REAL(RealKind), PARAMETER :: N_sc = 0.710         ! Schmidt-Zahl (PK, S.541)
  REAL(RealKind), PARAMETER :: D_v  = 3.0d-5        ! diffusivity of water vapour [m^2 s^-1]
  REAL(RealKind), PARAMETER :: n_f  = 0.333         ! Exponent von N_sc im Vent-koeff. (PK, S.541)
  REAL(RealKind), PARAMETER :: m_f  = 0.500         ! Exponent von N_re im Vent-koeff. (PK, S.541)
  REAL(RealKind), PARAMETER :: q_krit = 1.000d-9
  REAL(RealKind), PARAMETER :: q_krit_fr = 1.000d-6 ! critical value for freezing of raindrops
  REAL(RealKind), PARAMETER :: k_cr   = 5.25d0      ! Parameter in Accretion function 
  REAL(RealKind), PARAMETER :: q_krit_ii = 1.000d-6 ! q-Schwellenwert fuer ice_selfcollection 
  REAL(RealKind), PARAMETER :: D_krit_ii = 100.0d-6 ! D-Schwellenwert fuer ice_selfcollection 
  REAL(RealKind), PARAMETER :: D_krit_i  = 1.0d-5 ! D-Schwellenwert fuer ice_selfcollection

  REAL(RealKind), PARAMETER :: a_q = 0.429250753972163723
  REAL(RealKind), PARAMETER :: b_q = 0.180944125457962318
  REAL(RealKind), PARAMETER :: N_M92=1.0d3          !factor for the deposition-condensation formula after Meyers et al. (1992) m^-3
  REAL(RealKind), PARAMETER :: a_M92=-0.639         !factor a for the deposition-condensation formula after Meyers et al. (1992)
  REAL(RealKind), PARAMETER :: b_M92=12.96          !factor b for the deposition-condensation formula after Meyers et al. (1992)
  REAL(RealKind), PARAMETER :: c_M92=-2.80          !factor c for the deposition-condensation formula after Meyers et al. (1992)
  REAL(RealKind), PARAMETER :: d_M92=0.262          !factor d for the deposition-condensation formula after Meyers et al. (1992)
  REAL(RealKind), PARAMETER :: c_r = 1.0            ! p_rain.cap 
  REAL(RealKind), PARAMETER :: A_het = 2.0d-1       !kg^-1 s^-1 parameter for heterogeneous freezing, first 2.0d-1  
  REAL(RealKind), PARAMETER :: B_het = 6.5d-1       ! K^-1 parameter for heterogeneous freezing
  REAL(RealKind), PARAMETER :: T_f = 233.16d0       ! freezing point -40°C
  REAL(RealKind), PARAMETER :: T_3 = 273.16d0       ! tripel point

! Limiting raindrop size distributions (Seifert and Beheng, 2006, Appendix D)
  REAL(RealKind), PARAMETER :: xr_min     = 2.6d-10 ! kg 
  REAL(RealKind), PARAMETER :: N0_min     = 2.5d5   ! m^-4
  REAL(RealKind), PARAMETER :: N0_max     = 2.0d7   ! m^-4
  REAL(RealKind), PARAMETER :: lambda_min = 1.0d3   ! m^-1
  REAL(RealKind), PARAMETER :: lambda_max = 1.0d4   ! m^-1

! Parameters for different distribution functions and particle properties for the classes cloud water, rain water, ice and snow 
  REAL(RealKind), PARAMETER :: p_cloud_nu    = 0.333333
  REAL(RealKind), PARAMETER :: p_cloud_mu    = 0.666666
  REAL(RealKind), PARAMETER :: p_cloud_x_max = 2.60d-10
  REAL(RealKind), PARAMETER :: p_cloud_x_min = 4.20d-15
  REAL(RealKind), PARAMETER :: p_cloud_a_geo = 1.24d-01
  REAL(RealKind), PARAMETER :: p_cloud_b_geo = 0.333333
  REAL(RealKind), PARAMETER :: p_cloud_a_vel = 3.75d+05
  REAL(RealKind), PARAMETER :: p_cloud_b_vel = 0.666667
  REAL(RealKind), PARAMETER :: p_cloud_a_ven = 0.780000
  REAL(RealKind), PARAMETER :: p_cloud_b_ven = 0.308000
  REAL(RealKind), PARAMETER :: p_cloud_cap   = 2.0

  REAL(RealKind), PARAMETER :: p_rain_nu    = -0.666666
  REAL(RealKind), PARAMETER :: p_rain_mu    =  0.333333
  REAL(RealKind), PARAMETER :: p_rain_x_max = 3.00d-06
  REAL(RealKind), PARAMETER :: p_rain_x_min = 2.60d-10
  REAL(RealKind), PARAMETER :: p_rain_a_geo = 1.24d-01
  REAL(RealKind), PARAMETER :: p_rain_b_geo = 0.333333
  REAL(RealKind), PARAMETER :: p_rain_a_vel = 1.59d+02
  REAL(RealKind), PARAMETER :: p_rain_b_vel = 0.266667
  REAL(RealKind), PARAMETER :: p_rain_a_ven = 0.780000
  REAL(RealKind), PARAMETER :: p_rain_b_ven = 0.308000
  REAL(RealKind), PARAMETER :: p_rain_cap   = 2.0
  REAL(RealKind), PARAMETER :: p_rain_k_au  = 0.0
  REAL(RealKind), PARAMETER :: p_rain_k_sc  = 0.0
  REAL(RealKind), PARAMETER :: p_rain_k_c   = 0.0
  REAL(RealKind), PARAMETER :: p_rain_k_1   = 0.0
  REAL(RealKind), PARAMETER :: p_rain_k_2   = 2.0

! HK87 'hex plates'
  REAL(RealKind), PARAMETER :: p_ice_nu    = One     !after SB2006 -0.333333 coefficient for Gamma distribution
  REAL(RealKind), PARAMETER :: p_ice_mu    =  0.333333 !coefficient for Gamma distribution
  REAL(RealKind), PARAMETER :: p_ice_x_max = 1.00d-07  !maximum ice mass [kg]
  REAL(RealKind), PARAMETER :: p_ice_x_min = 1.00d-12  !minimum ice mass [kg]
  REAL(RealKind), PARAMETER :: p_ice_a_geo = 2.17d-01  !power law variable
  REAL(RealKind), PARAMETER :: p_ice_b_geo = 0.302115  !power law variable
  REAL(RealKind), PARAMETER :: p_ice_a_vel = 4.19d+01  !power law variable for fall speed
  REAL(RealKind), PARAMETER :: p_ice_b_vel = 0.260000  !power law variable for fall speed
  REAL(RealKind), PARAMETER :: p_ice_a_vel_SB = 317.0   !power law variable for fall speed after Seifert and Beheng (2005)
  REAL(RealKind), PARAMETER :: p_ice_b_vel_SB = 0.363  !power law variable for fall speed after Seifert and Beheng (2005)
  REAL(RealKind), PARAMETER :: p_ice_a_ven = 0.860000  !ventilation index a for the ventilation parameter for thin plates
  REAL(RealKind), PARAMETER :: p_ice_b_ven = 0.280000  !ventilation index b for the ventilation parameter for thin plates
  REAL(RealKind), PARAMETER :: p_ice_cap   = 3.141593
  REAL(RealKind), PARAMETER :: p_ice_k_au  = 0.0
  REAL(RealKind), PARAMETER :: p_ice_k_sc  = 0.0
  REAL(RealKind), PARAMETER :: p_ice_k_c   = 0.0
  REAL(RealKind), PARAMETER :: p_ice_k_1   = 0.0
  REAL(RealKind), PARAMETER :: p_ice_k_2   = 2.0

! Locatelli and Hobbs
  REAL(RealKind), PARAMETER :: p_snow_nu    = 0.500000
  REAL(RealKind), PARAMETER :: p_snow_mu    = 0.50000
  REAL(RealKind), PARAMETER :: p_snow_x_max = 2.00d-06
  REAL(RealKind), PARAMETER :: p_snow_x_min = 1.73d-09
  REAL(RealKind), PARAMETER :: p_snow_a_geo = 8.16d-00
  REAL(RealKind), PARAMETER :: p_snow_b_geo = 0.526316
  REAL(RealKind), PARAMETER :: p_snow_a_vel = 2.77d+01
  REAL(RealKind), PARAMETER :: p_snow_b_vel = 0.215790
  REAL(RealKind), PARAMETER :: p_snow_a_ven = 0.780000 !ventilation index a for the ventilation parameter
  REAL(RealKind), PARAMETER :: p_snow_b_ven = 0.308000 !ventilation index b for the ventilation parameter
  REAL(RealKind), PARAMETER :: p_snow_cap   = 2.0
  REAL(RealKind), PARAMETER :: p_snow_k_au  = 0.0
  REAL(RealKind), PARAMETER :: p_snow_k_sc  = 0.0
  REAL(RealKind), PARAMETER :: p_snow_k_c   = 0.0
  REAL(RealKind), PARAMETER :: p_snow_k_1   = 0.0
  REAL(RealKind), PARAMETER :: p_snow_k_2   = 2.0

! further parameters needed in autoconversion and selfcollection kernels 
  REAL(RealKind), PARAMETER :: p_cloud_k_c  = 4.44d+9     ! long-Kernel
  REAL(RealKind), PARAMETER :: p_cloud_k_1  = 4.00d+2     ! parameter for Phi-function.
  REAL(RealKind), PARAMETER :: p_cloud_k_2  = 0.7d+0      ! parameter for Phi-function.
  REAL(RealKind), PARAMETER :: p_cloud_k_au = 1.99047d+19 ! parameter for autoconversion.
  REAL(RealKind), PARAMETER :: p_cloud_k_sc = 2.05131d+10 ! parameter for selfcollection.
  REAL(RealKind), PARAMETER :: ice_s_vel    = 0.25
  REAL(RealKind), PARAMETER :: snow_s_vel   = 0.25

! parameters needed for fall speeds of rain and snow
  REAL(RealKind), PARAMETER :: l_min = 1.0d03 ! m^-1
  REAL(RealKind), PARAMETER :: l_max = 1.0d04 ! m^-1
  REAL(RealKind), PARAMETER :: N_min = 2.5d05 ! m^-4
  REAL(RealKind), PARAMETER :: N_max = 1.0d07 ! m^-4
  REAL(RealKind), PARAMETER :: a_r = 9.65d0   ! m s^-1
  REAL(RealKind), PARAMETER :: b_r = 10.3d0   ! m s^-1
  REAL(RealKind), PARAMETER :: c_rain = 600.0d0  ! m^-1

  TYPE MicroPhys_Variables_T
    REAL(RealKind) :: Rho
    REAL(RealKind) :: RhoV
    REAL(RealKind) :: RhoC
    REAL(RealKind) :: RhoR
    REAL(RealKind) :: RhoI
    REAL(RealKind) :: RhoS 
    REAL(RealKind) :: NC
    REAL(RealKind) :: NR
    REAL(RealKind) :: NI
    REAL(RealKind) :: NV
    REAL(RealKind) :: NS
    REAL(RealKind) :: TE 
  END TYPE MicroPhys_Variables_T

  TYPE(MicroPhys_Variables_T), ALLOCATABLE :: MPhys(:)

CONTAINS

FUNCTION delta_i(n)
  
  REAL(RealKind) :: delta_i 
  REAL(RealKind) :: n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  
  delta_i=(gfct((2.0*bi+nui+1.0+n)/mui)/gfct((nui+1.0)/mui)) &
         & *(gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(2.0*bi+n)

END FUNCTION delta_i

FUNCTION delta_ii(n)
  
  REAL(RealKind) :: delta_ii,n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
            
  delta_ii=2.0*(gfct((bi+nui+1.0+n)/mui)/gfct((nui+1.0)/mui))   &
          & *(gfct((bi+nui+1.0)/mui)/gfct((nui+1.0)/mui))         &
          & *(gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(bi+n)  &
          & *(gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(bi) 
END FUNCTION delta_ii

FUNCTION theta_i(n)
  
  REAL(RealKind) :: theta_i,n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  REAL(RealKind) :: betai=p_ice_b_vel_SB

  theta_i=gfct((2.0*betai+2.0*bi+nui+1+n)/mui)/gfct((2.0*bi+nui+1+n)/mui)* &
          (gfct((nui+1.0)/mui)/gfct((nui+2)/mui))**(2.0*betai)

END FUNCTION theta_i

FUNCTION theta_ii(n)

  REAL(RealKind) :: theta_ii,n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  REAL(RealKind) :: betai=p_ice_b_vel_SB
  
  theta_ii=2.0*(gfct((betai+bi+nui+1.0+n)/mui)/gfct((bi+nui+1.0+n)/mui))* &
               (gfct((betai+bi+nui+1.0)/mui)/gfct((bi+nui+1.0)/mui))* &
               (gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(2.0*betai)


END FUNCTION theta_ii

FUNCTION a_vent(n)

  REAL(RealKind) :: a_vent,n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu

  a_vent=p_ice_a_ven*(gfct((nui+bi+n)/mui)/(gfct((nui+1.0)/mui))) &
        & *(gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(bi+n-1.0)
END FUNCTION a_vent

FUNCTION b_vent(n)

  REAL(RealKind) :: b_vent,n
  REAL(RealKind) :: bi = p_ice_b_geo
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  REAL(RealKind) :: betai=p_ice_b_vel_SB

  b_vent=p_ice_b_ven*(gfct((nui+1.5*bi+n+0.5*betai)/mui))/(gfct((nui+1.0)/mui)) &
                    & *(gfct((nui+1.0)/mui)/gfct((nui+2.0)/mui))**(1.5*bi+0.5*betai+n-1.0)

END FUNCTION b_vent

FUNCTION fall_snow_n(RhoS,NS,Rho)
!Fall speed of NS
  REAL(RealKind) :: fall_snow_n,RhoS,NS,alf_n,c_lam,lam,Rho
  REAL(RealKind) :: nus=p_snow_nu
  REAL(RealKind) :: mus=p_snow_mu
  REAL(RealKind) :: alphas=p_snow_a_vel
  REAL(RealKind) :: betas=p_snow_b_vel
  REAL(RealKind) :: x_s
  c_lam=(gfct((nus+One)/mus)/gfct((nus+Two)/mus))
  alf_n=alphas*(gfct((nus+betas+One)/mus)/gfct((nus+One)/mus))

  IF (RhoS>1.0d-10.AND.PrecipSnow)THEN
    x_s=RhoS/(NS+Eps)
    x_s=MIN(MAX(x_s,p_snow_x_min),p_snow_x_max) 
    
    lam=(c_lam*x_s)**(betas)*SQRT(Rho0/(Rho+Eps))

    fall_snow_n=alf_n*lam
    fall_snow_n=MAX(fall_snow_n,0.1d0)
    fall_snow_n=MIN(fall_snow_n,3.0d0)
    fall_snow_n=-fall_snow_n
  ELSE
    fall_snow_n=Eps
  END IF

END FUNCTION fall_snow_n

FUNCTION fall_snow_q(RhoS,NS,Rho)
!Fall speed of RhoS
  REAL(RealKind) :: fall_snow_q,RhoS,NS,alf_q,c_lam,lam,Rho
  REAL(RealKind) :: nus=p_snow_nu
  REAL(RealKind) :: mus=p_snow_mu
  REAL(RealKind) :: alphas=p_snow_a_vel
  REAL(RealKind) :: betas=p_snow_b_vel
  REAL(RealKind) :: x_s
  c_lam= (gfct((nus+One)/mus)/gfct((nus+Two)/mus))
  alf_q=alphas*(gfct((nus+betas+Two)/mus)/gfct((nus+Two)/mus))

  IF (RhoS>1.0d-10.AND.PrecipSnow) THEN
    x_s=RhoS/(NS+Eps)
    x_s=MIN(MAX(x_s,p_snow_x_min),p_snow_x_max)   

    lam=(c_lam*x_s)**(betas)*SQRT(Rho0/(Rho+Eps))

    fall_snow_q=alf_q*lam
    fall_snow_q=MAX(fall_snow_q,0.1d0)
    fall_snow_q=MIN(fall_snow_q,3.0d0)
    fall_snow_q=-fall_snow_q
  ELSE
    fall_snow_q=Eps
  END IF

END FUNCTION fall_snow_q

FUNCTION ice_sedimentation(RhoI,NI)
!Fall speed of ice crystals just for the ISDAC case
  REAL(RealKind) :: ice_sedimentation,x_i,D_i,RhoI,NI
  REAL(RealKind) :: a_v=12.0d0
  REAL(RealKind) :: b_v=0.5d0
  REAL(RealKind) :: a_m=44.2d0
  REAL(RealKind) :: b_m=3.0d0

  x_i=RhoI/(NI+Eps)
  x_i=MIN(MAX(x_i,p_ice_x_min),p_ice_x_max)
  D_i=a_m*x_i**b_m

  IF(RhoI>1.0d-12.AND.D_i.GE.D_krit_i.AND.PrecipIce) THEN

    ice_sedimentation=-a_v*D_i**b_v
  ELSE
    ice_sedimentation=Eps
  END IF
END FUNCTION ice_sedimentation

FUNCTION fall_ice_n(RhoI,NI,Rho)
!Fall speed of ice crystals just for the ISDAC case
  REAL(RealKind) :: ice_sedimentation,D_i
  REAL(RealKind) :: a_v=12.0d0
  REAL(RealKind) :: b_v=0.5d0
  REAL(RealKind) :: a_m=44.2d0
  REAL(RealKind) :: b_m=3.0d0
!Fall speed of NI
  REAL(RealKind) :: fall_ice_n,RhoI,NI,alf_n,c_lam,lam,Rho
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  REAL(RealKind) :: alphai=p_ice_a_vel
  REAL(RealKind) :: betai=p_ice_b_vel
  REAL(RealKind) :: x_i
  c_lam=(gfct((nui+One)/mui)/gfct((nui+Two)/mui))
  alf_n=alphai*(gfct((nui+betai+One)/mui)/gfct((nui+One)/mui))

  IF (MicroScheme=='ISDAC') THEN
    x_i=RhoI/(NI+Eps)
    x_i=MIN(MAX(x_i,p_ice_x_min),p_ice_x_max)
!   D_i=a_m*x_i**b_m
!   IF(RhoI>1.0d-12.AND.D_i.GE.D_krit_i.AND.PrecipIce) THEN
!     fall_ice_n=-a_v*D_i**b_v
!   ELSE
!     fall_ice_n=Eps
!   END IF
    IF (RhoI>1.0d-12.AND.NI>0.0d0.AND.PrecipIce) THEN
      fall_ice_n=-6.387d0*x_i**(0.1666d0)
    ELSE
      fall_ice_n=Zero
    END IF
  ELSE
    IF (RhoI>1.0d-12.AND.NI>0.0d0.AND.PrecipIce) THEN
      x_i=RhoI/(NI+Eps)
      x_i=MIN(MAX(x_i,p_ice_x_min),p_ice_x_max)
      
      lam=(c_lam*x_i)**(betai)*SQRT(Rho0/(Rho+Eps))
  
      fall_ice_n=alf_n*lam
      fall_ice_n=MAX(fall_ice_n,0.1d0)
      fall_ice_n=MIN(fall_ice_n,3.0d0)
      fall_ice_n=-fall_ice_n
    ELSE
      fall_ice_n=Eps
    END IF
  END IF

END FUNCTION fall_ice_n

FUNCTION fall_ice_q(RhoI,NI,Rho)
!Fall speed of ice crystals just for the ISDAC case
  REAL(RealKind) :: ice_sedimentation,D_i
  REAL(RealKind) :: a_v=12.0d0
  REAL(RealKind) :: b_v=0.5d0
  REAL(RealKind) :: a_m=44.2d0
  REAL(RealKind) :: b_m=3.0d0
!Fall speed of RhoI
  REAL(RealKind) :: fall_ice_q,RhoI,NI,alf_q,c_lam,lam,Rho
  REAL(RealKind) :: nui=p_ice_nu
  REAL(RealKind) :: mui=p_ice_mu
  REAL(RealKind) :: alphai=p_ice_a_vel
  REAL(RealKind) :: betai=p_ice_b_vel
  REAL(RealKind) :: x_i

  c_lam=(gfct((nui+One)/mui)/gfct((nui+Two)/mui))
  alf_q=alphai*(gfct((nui+betai+Two)/mui)/gfct((nui+Two)/mui))

  IF (MicroScheme=='ISDAC') THEN
    x_i=RhoI/(NI+Eps)
    x_i=MIN(MAX(x_i,p_ice_x_min),p_ice_x_max)
!   D_i=a_m*x_i**b_m
!   IF(RhoI>1.0d-12.AND.D_i.GE.D_krit_i.AND.PrecipIce) THEN
!     fall_ice_q=-a_v*D_i**b_v
!   ELSE
!     fall_ice_q=Eps
!   END IF
    IF (RhoI>1.0d-12.AND.NI>0.0d0.AND.PrecipIce) THEN
      fall_ice_q=-6.387d0*x_i**(0.1666d0)
    ELSE
      fall_ice_q=Zero
    END IF
  ELSE
    IF (RhoI>1.0d-12.AND.NI>0.0d0.AND.PrecipIce) THEN
      x_i=RhoI/(NI+Eps)
      x_i=MIN(MAX(x_i,p_ice_x_min),p_ice_x_max)
      
      lam=(c_lam*x_i)**(betai)*SQRT(Rho0/(Rho+Eps))
  
      fall_ice_q=alf_q*lam
      fall_ice_q=MAX(fall_ice_q,0.1d0)
      fall_ice_q=MIN(fall_ice_q,3.0d0)
      fall_ice_q=-fall_ice_q
    ELSE
      fall_ice_q=Eps
    END IF
  END IF

END FUNCTION fall_ice_q

FUNCTION fall_rain_n(RhoR,NR,Rho)
!Fall speed of NR
  REAL(RealKind) :: fall_rain_n,RhoR,NR,Rho
  REAL(RealKind) :: N_0,x_r,lambda_r
  
  IF (RhoR>1.0d-10.AND.PrecipRain) THEN
    x_r=RhoR/(NR+Eps)
    x_r=MIN(MAX(x_r,p_rain_x_min),p_rain_x_max)
    N_0=MAX(N_min,MIN(N_max,NR*((Pi*rho0_liq)/x_r)**(One/3.0d0)))
    lambda_r=MAX(l_min,MIN(l_max,((Pi*rho0_liq*N_0)/RhoR)**(One/4.0d0)))
    
    fall_rain_n=a_r-b_r/(One+c_rain/lambda_r)
    fall_rain_n=SQRT(Rho0/Rho)*fall_rain_n
    fall_rain_n=MAX(fall_rain_n,1.0d-1)
    fall_rain_n=MIN(fall_rain_n,2.0d+1)
    fall_rain_n=-fall_rain_n
  ELSE
    fall_rain_n=Eps
  END IF

END FUNCTION fall_rain_n

FUNCTION fall_rain_q(RhoR,NR,Rho)
!Fall speed of RhoR
  REAL(RealKind) :: fall_rain_q,RhoR,NR,Rho
  REAL(RealKind) :: N_0,x_r,lambda_r
  
  IF (RhoR>1.0d-10.AND.PrecipRain) THEN
    x_r=RhoR/(NR+Eps)
    x_r=MIN(MAX(x_r,p_rain_x_min),p_rain_x_max)
    N_0=MAX(N_min,MIN(N_max,NR*((Pi*rho0_liq)/x_r)**(One/3.0d0)))
    lambda_r=MAX(l_min,MIN(l_max,((Pi*rho0_liq*N_0)/RhoR)**(One/4.0d0)))
            
    fall_rain_q=a_r-b_r*(One+c_rain/lambda_r)**(-4.0d0)
    fall_rain_q=SQRT(Rho0/Rho)*fall_rain_q
    fall_rain_q=MAX(fall_rain_q,1.0d-1)
    fall_rain_q=MIN(fall_rain_q,2.0d+1)
    fall_rain_q=-fall_rain_q
  ELSE
    fall_rain_q=Eps
  END IF

END FUNCTION fall_rain_q

END MODULE ParameterMicrophys_Mod
