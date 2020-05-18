subroutine vent_coeffxd(f_vent,DIFFSTAR,vel_part,r_part,TABS,RHOdyn)

! calculation of the ventilation coefficient

INCLUDE 'HEADER2'  

implicit none

double precision :: f_vent,X_vent,N_Sc,N_Re,eta_env,DIFFSTAR,vel_part
DOUBLE PRECISION :: r_part,TABS ,RHOdyn 

! Dynamic viscosity of air
eta_env=1.72D-5*(393.D0/(TABS +120.D0))*(TABS /273.D0)**1.5D0
! Reynolds number
N_Re=2.D0*vel_part*r_part*RHOdyn /eta_env
! Schmidt number
N_Sc=eta_env/(DIFFSTAR*RHOdyn )
X_vent=N_Sc**QU1D3 * N_Re **0.5D0
if(X_vent.lt.1.4D0) then
  f_vent=1.D0 + 0.108D0 * X_vent**2
else
  X_vent = MIN(X_vent,51.4D0)
  f_vent=0.78D0 + 0.308 * X_vent
endif

!      write(30,*) X_vent,f_vent,r_part*1.D6,N_Sc,N_Re
      

return
end
