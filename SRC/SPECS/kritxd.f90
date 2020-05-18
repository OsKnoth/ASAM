SUBROUTINE kritxd(r_qsa,r_sol,r_dry,m_sol,S_krit,r_krit,m_krit,a_sig,TABS,J,I,i_lm,j_lm,k_lm)
  
  ! calculation of the critical radius and supersaturation using 
  ! new Kelvin and Raoult terms
  ! for a given AP
  !
  ! called by cond_ini.f
  !           apvert.f
  !          (monito.f)
  !
  ! calling molal.f
  

  INCLUDE 'HEADER2'
  
  IMPLICIT NONE
  
  integer :: ITER,J,I,i_lm,j_lm,k_lm,l_lm,k_lm2
  double precision :: TABS,tau,a_sig,molality,phi_S,r_dry,r_qsa,r_sol,m_sol
  double precision :: r_krit,S_krit,m_krit,r_ini1,r_ini2,r_old,r_min,r_max
  double precision :: r_wet1,r_wet2,m_ggw1,m_ggw2,s_ggw1,s_ggw2,ds_ggw
  double precision :: ds_ggw_old,molality1,molality2,Phi_S1,Phi_S2
  
  ! temperature dependence of surface tension
  tau=1.D0-TABS/TC
  a_sig=b1*tau**mu*(1.D0+b2*tau)

  ! calculation of S_krit and r_krit/m_krit only in the 1st step of the
  ! Saturation iteration
  ! new version: iterative method for S_krit, r_krit:
  ! r_ini is chosen. r_ini1=r_ini * size_fact with size_fact = 1 + a small number.
  ! The difference between the saturation ratios calculated for r_ini and r_ini1
  ! is calculated as a measure for the slope of the function at r_ini.
  ! If slope < 0, then r_ini_new=r_ini / step_fact, 
  ! if slope > 0, then r_ini_new=r_ini*step_fact.
  ! Iteration steps are undertaken until the sign of the slope changes. 
  ! Then, the interval with r_krit is found. 
  ! Now the interval has to shrink in every iteration step until the
  ! required accuracy in r or s is reached.
  ! The geometric mean r_ini1 of the interval and the slope for this
  ! value are calculated. If the slope < 0, then r_max_new=r_ini2, if 
  ! slope > 0 then r_min_new=r_ini1 and the iteration goes on 
  ! until the required accuracy in r or s is reached.
  
  ! initialization of iteration
  ds_ggw_old=0.D0
  r_ini1=r_sol

  ! first: finding the correct interval

  IF (icond == 1) THEN 
     
     DO ITER=1,ITERMAX
        ! smaller particle
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        ! Exponent is done out of loop
        s_ggw1 = (2.D0*(a_sig+ bb*molality1)/(RHOW*RW*TABS*r_wet1)  &
             - molality1*vantS*mol_w*Phi_S1)
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)
        
        CALL molalxd(molality2,Phi_S2)
        
        ! Exponent is done out of loop
        s_ggw2 = (2.D0*(a_sig+ bb*molality2)/(RHOW*RW*TABS*r_wet2)  &
             - molality2*vantS*mol_w*Phi_S2)
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1
        
        ! if sign of slope has changed then correct interval is found
        IF (.NOT. (ds_ggw * ds_ggw_old < 0.D0 .OR. ds_ggw == 0.D0) ) EXIT
        
        ! for next iteration step
        IF(ds_ggw > 0.D0) THEN 
           r_old=r_ini1
           r_ini1=r_ini1*step_fact
        ELSE
           r_old=r_ini2
           r_ini1=r_ini2/step_fact
        ENDIF
        
        ds_ggw_old=ds_ggw
     ENDDO

     ! Just getting exp out of loop
     s_ggw1 = exp(s_ggw1) - 1.D0
     s_ggw2 = exp(s_ggw2) - 1.D0

     
  ELSE IF (icond == 2) THEN 
     
     DO ITER=1,ITERMAX
        
        ! smaller particle
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        s_ggw1=2.D0*(a_sig+ bb*molality1)/(RHOW*RW*TABS*r_wet1)      &
             - molality1*vantS*mol_w*Phi_S1
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)
        
        CALL molalxd(molality2,Phi_S2)
        
        s_ggw2=2.D0*(a_sig+ bb*molality2)/(RHOW*RW*TABS*r_wet2)      &
             - molality2*vantS*mol_w*Phi_S2
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1
        
        ! if sign of slope has changed then correct interval is found
        IF (.NOT. (ds_ggw*ds_ggw_old < 0.D0.or.ds_ggw == 0.D0) ) EXIT
        
        ! for next iteration step
        IF (ds_ggw > 0.D0) THEN
           r_old=r_ini1
           r_ini1=r_ini1*step_fact
        ELSE
           r_old=r_ini2
           r_ini1=r_ini2/step_fact
        ENDIF
        
        ds_ggw_old=ds_ggw
        
     ENDDO
     
  ELSE IF (icond == 3) THEN 

     DO ITER=1,ITERMAX
        
        ! smaller particle
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        s_ggw1=2.D0*a_sig/(RHOW*RW*TABS*r_ini1)                      &
             - molality1*vantS*mol_w
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)
        
        CALL molalxd(molality2,Phi_S2)
     
        s_ggw2=2.D0*a_sig/(RHOW*RW*TABS*r_ini2)                      &
             - molality2*vantS*mol_w
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1
        
        ! if sign of slope has changed then correct interval is found
        IF (.NOT. (ds_ggw*ds_ggw_old < 0.D0.or.ds_ggw == 0.D0) ) EXIT
        
        ! for next iteration step
        IF (ds_ggw > 0.D0) THEN
           r_old=r_ini1
           r_ini1=r_ini1*step_fact
        ELSE
           r_old=r_ini2
           r_ini1=r_ini2/step_fact
        ENDIF
        
        ds_ggw_old=ds_ggw
        
     ENDDO
     
  ENDIF

  ! If iteration goes wrong, print some values
  IF (ITER >= ITERMAX) THEN
     write(*,*) "krit.f a1: ",s_ggw1,ds_ggw,r_ini1*1.D9,r_wet1*1.D9,molality1
     write(*,*) "krit.f a2: ",s_ggw2,ds_ggw,r_ini2*1.D9,r_wet2*1.D9,molality2
     write(*,*) J,I," krit.f a3: ",ITER,"r_qsa",r_qsa*1.D9,"r_sol",     &
          r_sol*1.D9,"r_dry",r_dry*1.D9,"m_sol",m_sol,"a_sig",a_sig
  ENDIF
  
  
  r_min = MIN(r_old,r_ini1)
  r_max = MAX(r_old,r_ini2)
  
  IF (icond == 1) THEN
     
     DO ITER=1,500
        r_ini1=sqrt(r_min*r_max)
        r_ini2=r_ini1*size_fact
        
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        ! Exponent is done out of loop
        s_ggw1 = (2.D0*(a_sig+ bb*molality1)/(RHOW*RW*TABS*r_wet1)  &
             - molality1*vantS*mol_w*Phi_S1)
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)
        
        CALL molalxd(molality2,Phi_S2)
        
        ! Exponent is done out of loop
        s_ggw2 = (2.D0*(a_sig+ bb*molality2)/(RHOW*RW*TABS*r_wet2)  &
             - molality2*vantS*mol_w*Phi_S2)
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1

        if(ds_ggw < 0.D0) then
           r_max=r_ini2
        else
           r_min=r_ini1
        endif
        if((r_max/r_min) <= accu_fact) exit
     ENDDO
     
     ! Just getting exp out of loop
     s_ggw1 = exp(s_ggw1) - 1.D0
     s_ggw2 = exp(s_ggw2) - 1.D0

  ELSE IF (icond == 2) THEN 
     
     DO ITER=1,500
        r_ini1=sqrt(r_min*r_max)
        r_ini2=r_ini1*size_fact
        
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        s_ggw1=2.D0*(a_sig+ bb*molality1)/(RHOW*RW*TABS*r_wet1)      &
             - molality1*vantS*mol_w*Phi_S1
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)

        CALL molalxd(molality2,Phi_S2)
        
        s_ggw2=2.D0*(a_sig+ bb*molality2)/(RHOW*RW*TABS*r_wet2)      &
             - molality2*vantS*mol_w*Phi_S2
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1
        IF (ds_ggw < 0.D0) THEN
           r_max=r_ini2
        ELSE
           r_min=r_ini1
        ENDIF
        
        IF ((r_max/r_min) <= accu_fact) EXIT
        
     ENDDO
     
  ELSE IF (icond == 3) THEN 
     
     DO ITER=1,500
        r_ini1=sqrt(r_min*r_max)
        r_ini2=r_ini1*size_fact
        
        r_wet1=(r_dry**3 + r_ini1**3)**QU1D3
        m_ggw1=FACT*r_ini1**3
        molality1=m_sol/(mol_AS*m_ggw1)
        
        CALL molalxd(molality1,Phi_S1)
        
        s_ggw1=2.D0*a_sig/(RHOW*RW*TABS*r_ini1)                      &
             - molality1*vantS*mol_w
        
        ! larger particle
        r_ini2=r_ini1*size_fact
        r_wet2=(r_dry**3 + r_ini2**3)**QU1D3
        m_ggw2=FACT*r_ini2**3
        molality2=m_sol/(mol_AS*m_ggw2)

        CALL molalxd(molality2,Phi_S2)
        
        s_ggw2=2.D0*a_sig/(RHOW*RW*TABS*r_ini2)                      &
             - molality2*vantS*mol_w
        
        ! slope
        ds_ggw=s_ggw2-s_ggw1
        if(ds_ggw < 0.D0) then
           r_max=r_ini2
        else
           r_min=r_ini1
        endif
        
        if((r_max/r_min) <= accu_fact) exit
        
     ENDDO
     
  ENDIF
  

  ! final choice of critical values
  S_krit=MAX(s_ggw2,s_ggw1)
  r_krit=MIN(r_max,r_min)
  m_krit=r_krit**3*FACT


  
  
  RETURN
END subroutine kritxd








