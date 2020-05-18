subroutine term_velxd(vel_part,r_part,RHOdyn)

! calculation of the terminal drop velocity

INCLUDE 'HEADER2' 

implicit none

double precision :: vel_part,r_part,RHOdyn ,vel_stokes,vel_large

! Roger & Yau 
! Stokes terminal velocity 
vel_stokes=k1_vel * r_part**2
! For large drops
vel_large=k2_vel * SQRT(r_part)
vel_large=MIN(vel_large,10.D0)
! density (height) dependence
vel_large=vel_large * SQRT(rho_0/RHOdyn )
! choice of regime
vel_part=MIN(vel_stokes,vel_large)

!      write(77,*) vel_part,height,vel_large,vel_stokes
!      write(78,*) r_part*1.D6,vel_part,vel_large,vel_stokes

      
return
end
