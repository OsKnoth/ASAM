! calculate here if particle is larger than the drop
! if so, contact freezing is possible

subroutine contact_freezing_matrixxd
  
  INCLUDE 'HEADER3'
  
  implicit none
  
  INTEGER :: I,J
  double precision :: r_fact
  
  r_fact=0.91D0    ! m_drop < 0.75 m_ap
  do I = 1,SMAX
     do J = 1,JMAX
        if(RMITTE(J).LT.r_fact*RSMITTE(I)) then
           CONTACT(J,I)=1
        else 
           CONTACT(J,I)=0
        endif
     enddo
  enddo
  ! 111  FORMAT(2E14.6,3I14.2)
  
  return
end subroutine contact_freezing_matrixxd
