! subroutine to calculate which fraction of drops freeze after collision
! with aps; temperature dependant
! for different insoluble aerosol particles
! factor should not be negative and not larger than 1 => 1 means 100% of 
! colliding drops freeze
! temperature in C!

subroutine contact_tempxd(TABS,fac_contact,ifreeze,IT)
  
  INCLUDE 'HEADER2'   
  
  implicit none
  
  integer :: I,IT,IP,ifreeze(ITMAX)
  DOUBLE PRECISION :: TABS, fac_contact, T_C, T_C0
  DOUBLE PRECISION :: a,b
  parameter (T_C0 = 273.d0)
  
  ! to exclude contact freezing !!!!!!!!!!!!!!!
  IP=1
  fac_contact = 0.d0
  if(ikofr.eq.0) return
  if(ifreeze(IT).eq.0) return
  
  ! bacteria
  !      if(ikofr.eq.1) then
  if(ifreeze(IT).eq.1) then
     a = -26.4
     b = -74.2
  endif
  
  ! leaf litter
  !      if(ikofr.eq.2) then
  if(ifreeze(IT).eq.2) then
     a = -8.91
     b = -46.72
  endif
  
  ! pollen
  !      if(ikofr.eq.3) then
  if(ifreeze(IT).eq.3) then
     a = -8.64
     b = -56.3
  endif
  
  ! montmorillonite
  !      if(ikofr.eq.5) then
  if(ifreeze(IT).eq.5) then
     a = -10.14
     b = -32.77
  endif
  
  ! kaolinite
  !      if(ikofr.eq.6) then
  if(ifreeze(IT).eq.6) then
     a = -10.07
     b = -69.35
  endif
  
  ! soot
  !      if(ikofr.eq.7) then
  if(ifreeze(IT).eq.7) then
     a = -6.2
     b = -57.3
  endif
  
  T_C = TABS- T_C0
  fac_contact = 1.d-2 * (a*T_C + b)
  fac_contact = max(0.d0,fac_contact)
  fac_contact = min(1.d0,fac_contact)
  
  !      print *,'temp_contact.f: TABS,fac_contact'
  !     &       ,TABS,fac_contact
  
  if(iice.eq.0) fac_contact = 0.d0    ! no ice
  
  
  return
end subroutine contact_tempxd
