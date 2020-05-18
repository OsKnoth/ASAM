subroutine ggw_newxd(m_ggw,r_ggw,SATT,TABS,r_qsa,r_sol,r_dry,m_sol,   &
                     a_sig,r_krit,s_krit)
  
  ! Msol is the solved mass of the particle, SATT the supersaturation
  ! iterative calculation of the euilibrium size of an AP for
  ! a given saturation ratio
  ! called by cond.f
  !           apvert.f
  ! calling molal.f
  
  INCLUDE 'HEADER2'
  
  implicit none
  
  integer :: i_w,j_w,i_S,j_S,ITER
  double precision :: SATT,m_ggw,r_ggw,Msol,m_ggw_i,m_ggw_j,fact_i0,fact_i1
  double precision :: fact_j0,fact_j1,fact_S0,fact_S1,sat_i,sat_j,a0,n_S,b_S
  double precision :: r_krit,s_krit,r_qsa,r_sol,r_dry,r_wet,m_sol
  double precision :: a_sig,molality,phi_S
  double precision :: r_min,r_max,s_ggw,TABS,SATTtop,SATTbot
  
  
  !print *,'ggw_newxd:SATT,s_krit',SATT,s_krit
  ! no if as satt smaller, s_krit is zero
  if(SATT.gt.s_krit) then
     r_ggw = r_krit
     m_ggw = FACT*r_ggw**3
     ! exits subroutine!
     return
  endif

  ! SATT boundaries for iteration
  SATTtop = MAX(SATT*(1.D0+SC_ggw),SATT*(1.D0-SC_ggw))
  SATTbot = MIN(SATT*(1.D0+SC_ggw),SATT*(1.D0-SC_ggw))
  if(ABS(SATT).lt.SATTabs) then
     SATTtop =  SATTabs
     SATTbot = -SATTabs
  endif

  ! bisectional method for the calculation of the AP equilibrium size
  ! within a size interval [r_min,r_max] with r_max=r_krit und r_min small
  ! Start with interval center, calculation of the corresponding SATT
  ! iteration until given accuracy in SATT is reached
  
  ! Estimation of the boundaries
  r_max = r_krit
  r_min = r_max
  ITER=0

111 continue
  ITER=ITER+1

!!!      if(ITER.gt.ITERMAX)
!!!     &   write(*,*) "ggw_new.f a: ITER zu groﬂ",ITER,s_ggw
!!!     &  ,r_min,r_wet,molality

  if(ITER.gt.500) goto 222
  
  r_min = r_min*0.5D0
  ! new equation for AP equilibrium size (cf. drop growth equation:
  ! Kelvin/Raoult)
  r_wet = (r_dry**3 + r_min**3)**QU1D3
  !      r_wet=MAX(r_qsa,r_wet)
  m_ggw = FACT*r_min**3
  
  molality=m_sol/(mol_AS*m_ggw)
  call molalxd(molality,Phi_S)
  if(icond.eq.1) s_ggw = exp(2.D0*(a_sig+ bb*molality)/(RHOW*RW*TABS*r_wet)  &
                                         - molality*vantS*mol_w*Phi_S)-1.D0
  if(icond.eq.2) s_ggw = 2.D0*(a_sig+ bb*molality)/(RHOW*RW*TABS*r_wet)      &
                                         - molality*vantS*mol_w*Phi_S
  if(icond.eq.3) s_ggw = 2.D0*a_sig/(RHOW*RW*TABS*r_min) - molality*vantS*mol_w
!!!      if(s_ggw.gt.10.D0) then
!!!        write(*,*) ITER,"ggw_new.f: s_ggw 1",s_ggw,1.D9*r_wet
!!!     &  ,1.D9*r_dry,1.D9*r_min,1.D9*r_krit
!!!      endif
  if(s_ggw.ge.SATT) then
     r_max = r_min
     goto 111
  endif
  ! equilibrium
  ITER=0
222 continue

  ITER=ITER+1
!!!      if(ITER.gt.500)
!!!     &   write(*,*) "ggw_new.f b: ITER zu groﬂ",ITER,s_ggw
!!!     &  ,r_ggw,r_wet,molality
  
  if(ITER.gt.500) goto 333
  ! new equation for AP equilibrium size (cf. drop growth equation:
  ! Kelvin/Raoult)
  r_ggw = (r_max + r_min)*0.5D0
  r_wet = (r_dry**3 + r_ggw**3)**QU1D3
  !      r_wet=MAX(r_qsa,r_wet)
  m_ggw = FACT*r_ggw**3
  !print *,'ggw_newxd 2:m_ggw',m_ggw
  molality = m_sol/(mol_AS*m_ggw)
  call molalxd(molality,Phi_S)
  if(icond.eq.1) s_ggw = exp(2.D0*(a_sig+ bb*molality)/(RHOW*RW*TABS*r_wet)  &
                                         - molality*vantS*mol_w*Phi_S)-1.D0
  if(icond.eq.2) s_ggw=2.D0*(a_sig+ bb*molality)/(RHOW*RW*TABS*r_wet)      &
                                         - molality*vantS*mol_w*Phi_S
  if(icond.eq.3) s_ggw=2.D0*a_sig/(RHOW*RW*TABS*r_ggw) - molality*vantS*mol_w
  if(s_ggw.gt.SATTtop) then
     r_max = r_ggw
     goto 222
  endif
  if(s_ggw.lt.SATTbot) then
     r_min = r_ggw
     goto 222
  endif
  
333 continue
  
  return
end subroutine ggw_newxd
