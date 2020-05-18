subroutine molalxd(molality,Phi_S)

! Calculation of the osmotic coefficient for the Raoult term

INCLUDE 'HEADER2'  

implicit none

integer :: i
double precision :: Phi_S,molality


! ideal solution
if(iideal.eq.1) then
  Phi_S=1.D0
  return
endif
! linear function
do i=1,6
  if(molality.ge.mol_g(i).and.molality.lt.mol_g(i+1)) then
    Phi_S=a_phi(i)*molality + b_phi(i)
    goto 100
  endif
enddo
Phi_S=b_phi(7)
100 continue


return
end
