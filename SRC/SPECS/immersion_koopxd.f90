! calculate here if drop freezes by immersion freezing
! depends on temperature, cooling rate, drop volume and particle type
! Immersion freezing occurs still instantaneously

SUBROUTINE immersion_koopxd(TABS,TABSold,NW,QW,QS,QA,IMMERQ,IMMERN,   &
                            IMMERA,IMMERS,ifreeze,ip)
  
  INCLUDE 'HEADER3'
  
  IMPLICIT NONE
  
  INTEGER :: I,J,IP,ifreeze(ITMAX),JSTART
  DOUBLE PRECISION :: IMMERQ(JMAX,SMAX),IMMERN(JMAX,SMAX)
  DOUBLE PRECISION :: IMMERA(JMAX,SMAX),IMMERS(JMAX,SMAX)
  DOUBLE PRECISION :: NW(JMAX,SMAX,IPMAX),QW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: TABS, T_C, T_C0,TABSold
  DOUBLE PRECISION :: b
  DOUBLE PRECISION :: DROPVOL(JMAX),DROPMASS(JMAX)
  DOUBLE PRECISION :: APMASS(JMAX,SMAX),MOLALITY(JMAX,SMAX)
  DOUBLE PRECISION :: PREFAC, SMALL
  DOUBLE PRECISION :: a_w, T_freeze,T_freeze_water,temp_dep_koop
  DOUBLE PRECISION :: molal
  PARAMETER (T_C0 = 273.16d0)
  
  IF(TABS.GE.T_C0) RETURN
  IF(TABSold.LE.TABS) RETURN
  !JSTART = 1     ! Bin where immersion freezing can start
  JSTART = JMAX/2
  ! number of liquid drops = NW should not be too small!!
! VG
!  SMALL = 1.d-5
  SMALL = 1.d0
  
  ! units in cm!!
  ! for bacteria
  IF(iimfr.EQ.1) b = 6.19
  
  ! for leaf litter
  IF(iimfr.EQ.2) b = 0.438
  
  ! for pollen
  IF(iimfr.EQ.3) b = 1.01d-2
  
  ! for illite
  IF(iimfr.EQ.4) b = 6.19d-5
  
  ! for montmorillonite
  IF(iimfr.EQ.5) b = 3.23d-5 
  
  ! for kaolinite
  IF(iimfr.EQ.6) b = 6.15d-8 
  
  ! for soot
  IF(iimfr.EQ.7) b = 2.91d-9
  
  ! -dnf/dt = nu * b *vd * dexp(-T) * dT/dt   ! freezing rate
  DO J=JSTART,JMAX
     dropvol(J) = 4.0D0/3.0D0*PI*(RMITTE(J)*1.0D+2)**3.0D0   ! in cm3
     !  dropmass(J) = dropvol(J)*1.0D-6*RHOW                    ! in kg
     !  dropvol(J) = RMITTE(J)**3.0D0/FAC34                     ! in m3
     dropmass(J) = MMITTE(J)
  END DO
  
  ! NW unit is m-3!!
  ! eigentlich kg-1 ...
  DO I=1,SMAX
     DO J=JSTART,JMAX
        IF(NW(J,I,IP).GT.SMALL) THEN
           apmass(J,I) = QS(J,I,IP)/NW(J,I,IP)               ! in kg
           molality(J,I) = (QS(J,I,IP)/mol_AS)/QW(J,I,IP) 
           !        molality(J,I) = (apmass(J,I)/mol_AS)/dropmass(J)
           ! calculate freezing point depression! First the freezing temperature of 
           ! pure water
           ! to calculate freezing temperature of pure water, call koop2 with a_w =1
           ! later it will be called with other values for a_w!!
           !        a_w = 1       ! activity of pure water
           !        CALL koop2(a_w,T_freeze_water)
           T_freeze_water = 240.739D0
           molal = molality(J,I)
           
           CALL koop_freeze(molal,a_w,T_freeze)
           
           temp_dep_koop = T_freeze - T_freeze_water
           !        write(99,*) molal, a_w,temp_dep_koop
           ! to consider freezing point depression in the equation of the freezing rate:
           ! add  the freezing point depression to the actual temperature! 
           ! then the value of the exponent will be less and less drops will freeze
           ! CAUTION! if value of freezing point depression larger than value of 
           ! temperature exp will be positive -> wrong values! equation works only 
           ! as long as exp is negative or zero
           ! if freezing point depression is larger than supercooling, 
           ! drop should not freeze!
           IF(molal.GE.5.0D0) THEN
              PREFAC = 0.0D0
           ELSE
              IF(DABS(temp_dep_koop).GT.DABS(TABS - T_C0)) THEN
                 PREFAC = b*DEXP(0.0D0)*(TABS - TABSold)/DELTAT
                 PREFAC = 0.0D0
              ELSE
                 PREFAC = b*DEXP(-(TABS - T_C0) + temp_dep_koop)          &
                      *(TABS - TABSold)/DELTAT
                 !               IMMERN(J,I)= - PREFAC * DROPVOL(J) * (NW(J,I,IP) * 1.d-6)
                 !               IMMERN(J,I)= IMMERN(J,I) * 1.d6 ! in m-3
                 ! DROPVOL muss in cm^3 eingesetzt werden!
                 ! TEST: Verringerung um Faktor 10
                 PREFAC = 0.001D0 * PREFAC 
                 IMMERN(J,I) = -PREFAC*DROPVOL(J)*NW(J,I,IP) ! 1%: *0.01
                 IMMERN(J,I) = MIN(IMMERN(J,I),NW(J,I,IP))
                 IMMERN(J,I) = MAX(IMMERN(J,I),0.D0)
                 IMMERQ(J,I) = IMMERN(J,I)/NW(J,I,IP)*QW(J,I,IP)
                 IMMERQ(J,I) = MIN(IMMERQ(J,I),QW(J,I,IP))
                 IMMERA(J,I) = IMMERN(J,I)/NW(J,I,IP)*QA(J,I,IP)
                 IMMERA(J,I) = MIN(IMMERA(J,I),QA(J,I,IP))
                 IMMERS(J,I) = IMMERN(J,I)/NW(J,I,IP)*QS(J,I,IP)
                 IMMERS(J,I) = MIN(IMMERS(J,I),QS(J,I,IP))
                 !          write(*,*) J,I,'imme1',PREFAC,molal
                 !          write(*,*) temp_dep_koop,TABS - TABSold,TABS
                 !          write(*,*) IMMERN(J,I)/NW(J,I,IP),NW(J,I,IP),IMMERN(J,I)
              END IF ! freezing point depression
           END IF ! molality
        END IF ! number
     END DO
  END DO
  
  RETURN
END SUBROUTINE immersion_koopxd


!------------------------------------------------------------------------
! subroutine koop2 
!
! this subroutine calculate the ice nucleation in solution drops
! follwing the approach of 
! Koop, Luo, Tsias, and Peter, Nature, vol 406, 2000, p. 611 ff
!-------------------------------------------------------------------------

subroutine koop2(a_w,T_freeze)
  
  implicit none
  
  double precision :: temp_step, temp_step2
  parameter (temp_step=0.25D0)   !Temperature step first loop  (T decreases until j>1)
  parameter (temp_step2=0.001D0) !Temperature step second loop (T increases until j<1)
  double precision :: T_start
  parameter (T_start=241.00D0) !Largest possible T for a_w=1 (240.739)
  double precision :: T            ! temperature in K
  double precision :: a_wi         ! activity coefficient of water/ice mix
  double precision :: a_w          ! activity of pure water, 
  ! or different solutions with molal e.g. 0.1 
  double precision :: delta_a_w_i  ! difference of the potentials of water 
  !  in pure ice and pure liquid
  double precision :: delta_mu
  double precision :: gas_const 
  parameter (gas_const = 8.3143D0)    ! j/(mol *K), from Mortimer
  double precision :: J, log_J        ! ice nucelation rate coefficient
  double precision :: J_threshold
  parameter (J_threshold = 1.D0)
  double precision :: T_freeze        ! freezing temperature
  
  T = T_start
  J=0.D0
  do while(J.lt.J_threshold)
     delta_mu = 210368 + 131.438 * T - 3.32373e6 * 1./T - 41729.1 * dlog(T)
     a_wi = dexp(delta_mu/(gas_const * T))
     ! instead of a_w put here the activity of the solution 
     delta_a_w_i = a_w - a_wi
     log_J = -906.7 + 8502 * delta_a_w_i - 26924 * (delta_a_w_i)**2     &
             + 29180 * (delta_a_w_i )**3 
     J = dexp(log_J)
     if(J.ge.1.) then
        T_freeze= T
     endif
     T = T - temp_step
  enddo ! end of while loop
  
  ! set temperature for next iteration
  T = T_freeze
  
  do while(J.gt.J_threshold)
     delta_mu = 210368 + 131.438 * T - 3.32373e6 * 1./T - 41729.1 * dlog(T)
     a_wi = dexp(delta_mu/(gas_const * T))
     ! instead of a_w put here the activity of the solution 
     delta_a_w_i = a_w - a_wi
     log_J = -906.7 + 8502 * delta_a_w_i - 26924 * (delta_a_w_i)**2    &
          + 29180 * (delta_a_w_i )**3 
     J = dexp(log_J)
     if(J.ge.1.) then
        T_freeze = T
     endif
     T = T_freeze + temp_step2
  enddo ! while loop
  
  return
end subroutine koop2

!------------------------------------------------------------------------
! subroutine koop_freeze
!
! this subroutine contains data fields for NaCl and (NH4)2SO4 for molality
! and activity coefficients
! values taken from Pruppacher and Klett 1997, p.112-113
! with n+2 values, value 1 = pure water, value n+2 that of saturated solution 
!  (p.113) inbetween from table on p. 112
!
! 1) look where the drop molality fits in
! 2) interpolate activity coefficient for drop molality
! 3) call koop2, which is the routine for the freezing temperature
!
!------------------------------------------------------------------------

subroutine koop_freeze(molal,a_w,T_freeze)
  
  IMPLICIT NONE
  double precision :: molal
  double precision :: a_w
  double precision :: T_freeze
  integer :: nval
  parameter(nval=7) ! n+2
  double precision :: molal_tab(nval)
  double precision :: activity_tab(nval)
  double precision :: TINY ! to avoid unnecessary interpolation
  double precision :: checkit
  integer :: iflag
  integer :: i, istart,iend
  integer :: j1,j2,j3 ! for interpolation
  
  
  !      data molal_tab/0.,0.1,0.5,1.0,2.0,5.0,6.2/ ! NaCl
  !      data activity_tab/1.,0.997,0.984,0.967,0.932,0.807,0.753/
  
  data molal_tab/0.,0.1,0.5,1.0,2.0,5.0,5.7/ ! (NH4)2SO4
  data activity_tab/1.,0.996,0.982,0.966,0.935,0.831,0.78/
  
  iflag = 0
  istart = 1
  iend = 1
  
  ! avoid that molality larger than saturation molality
  
  if(molal.ge.molal_tab(nval)) then !if1
     a_w = activity_tab(nval)
     !         print *,'molality larger than saturation molality!!!'
  else ! if1
     do i =1,nval-1
        ! avoid unnecessary interpolation
        if(molal.lt.0.1) then
           TINY=1.e-3
        else
           TINY = 1.e-2
        endif
        checkit = abs(molal - molal_tab(i))
        if(checkit.le.TINY) then ! if2
           a_w=activity_tab(i)
           !               print *,'take activity value directly from table'
           iflag = 1 ! do not interpolate
        else ! if2
           if(molal.ge.molal_tab(i).and.molal.lt.molal_tab(i+1)) then !if3
              istart=i
              iend=i+1
           endif            !if3
        endif               ! if2
     enddo
     if(istart.lt.iend)then
        if(iflag.eq.0)then
           ! 3 point interpolation
           ! 
           ! avoid here to get out of range of table
           j1 = istart - 1
           j1 = max(1,j1)
           j2 = istart
           j3 = iend
           j3 = min(iend,nval)
           call interpol(molal,molal_tab(j1),molal_tab(j2),molal_tab(j3),   &
                         a_w,activity_tab(j1),activity_tab(j2),activity_tab(j3))  
        endif                  ! iflag
     endif
     !         print *,'molality,istart,iend,molality(i,i+1)'
     !         print *, molal, istart, iend
     !     &        ,molal_tab(istart)
     !     &        ,molal_tab(iend)
     if(a_w.gt.1.) print *, 'a_w.gt.1.!!!'
  endif                     ! if1
  
  !      print *,'a_w',a_w
  !      print *,'molality',molal
  
  call koop2(a_w,T_freeze)
    
  return
end subroutine koop_freeze


!************************************************************
!*                                                          *
!*    subroutine interpol
!*    
!*    version from Nov. 2002  Sabine Wurzler                *
!*                                                          *
!************************************************************

subroutine interpol(x,x1,x2,x3,y,y1,y2,y3)
  
  !     interpolates from a table
  !     1 dimensional linear Lagrange Interpolation
  !     this version interpolates between the values directly
  !     includes a check for over- and undershooting 
  !
  ! want to know y at value x
  ! Stuetzstellen rlx1, rlx2, rlx3
  ! Differenzen x-xi: tex1, tex2, tex3
  ! Produkte: dex1, dex2, dex3 
  
  IMPLICIT NONE
  double precision :: x,x1,x2,x3
  double precision :: y,y1,y2,y3
  double precision :: rln1,rlx1,rlx2,rlx3
  double precision :: tex1,tex2,tex3
  double precision :: dex1,dex2,dex3
  integer :: nstuetz  ! number of points for interpolation
  parameter (nstuetz = 3)
  double precision :: psix(nstuetz)
  double precision :: xy(nstuetz)
  double precision :: TINY
  parameter (TINY=1.e-20)
  
  rln1=x
  rlx1=x1
  rlx2=x2
  rlx3=x3
  
  tex1=rln1-rlx1
  tex2=rln1-rlx2
  tex3=rln1-rlx3
  
  dex1=(rlx1-rlx2)*(rlx1-rlx3)
  dex2=(rlx2-rlx1)*(rlx2-rlx3)
  dex3=(rlx3-rlx1)*(rlx3-rlx2)
  
  ! hier muss nulldivision abgefangen werden!!!!
  if(dex1.ge.TINY) then
     psix(1)=abs(tex2*tex3)/abs(dex1)
  endif
  
  if(dex2.ge.TINY) then
     psix(2)=abs(tex1*tex3)/abs(dex2)
  endif
  
  if(dex3.ge.TINY) then
     psix(3)=abs(tex1*tex2)/abs(dex3)
  endif
  
  if(abs(rlx1-rlx2).le.TINY.and.abs(rlx1-rlx3).ge.TINY) then
     psix(3)=rln1/abs(rlx1-rlx3)
     psix(2)=0.
     psix(1)=1.-psix(3)
     !         print *,'version1'
  endif
  
  if(abs(rlx1-rlx2).ge.TINY.and.abs(rlx1-rlx3).le.TINY) then
     psix(1)=rln1/abs(rlx1-rlx2)
     psix(2)=1.-psix(1)
     psix(3)=0.
     !         print *,'version2'
  endif
  
  if(tex1.le.TINY.and.abs(rlx2-rlx3).gt.TINY) then
     psix(1)=rln1/abs(rlx2-rlx3)
     psix(2)=1.-psix(1)
     psix(3)=0.
     !         print *,'version 6'
  endif
  
  if(tex2.le.TINY.and.abs(rlx1-rlx3).gt.TINY) then
     psix(1)=0.
     psix(2)=rln1/abs(rlx1-rlx3)
     psix(3)=1.-psix(2)
     !         print *,'version 7'
  endif
  
  if(tex3.le.TINY.and.abs(rlx1-rlx2).gt.TINY) then
     psix(3)=rln1/abs(rlx1-rlx2)
     if(psix(3).gt.1.) then
        psix(3)=1./psix(3)
        psix(3)=1.-psix(3)
        !            print *,'turn psix3'
     endif
     psix(2)=1-psix(3)
     psix(1)=0.
     !         print *,'version 8'
     !         print *, 'version 8, psix(1),psix(2),psix(3)', 
     !     &             psix(1),psix(2),psix(3)
  endif
  
  ! maybe check necessary for rlx1-rlx3 > 0., rlx2-rlx3 > 0.rlx1 -rlx2 > 0.
  
  goto 11
  
  if(psix(1).lt.TINY) then
     psix(1)=0.
     psix(2)=rln1/abs(rlx1-rlx3)
     psix(3)=1.-psix(2)
     !         print *,'version3'
  endif
  
  if(psix(2).lt.TINY) then
     psix(1)=rln1/abs(rlx2-rlx3)
     psix(2)=0.
     psix(3)=1.-psix(1)
     !         print *,'version4'
  endif
  
  if(psix(3).le.TINY) then
     psix(1)=rln1/abs(rlx1-rlx2)
     psix(2)=1.-psix(1)
     psix(3)=0.
     !         print *,'version5'
  endif
  
11 continue
  
  psix(1)=max(psix(1),0.0D0)
  psix(2)=max(psix(2),0.0D0)
  psix(3)=max(psix(3),0.0D0)
  
  psix(1)=min(psix(1),1.0D0)
  psix(2)=min(psix(2),1.0D0)
  psix(3)=min(psix(3),1.0D0)
  
  xy(1)=psix(1)*y1
  xy(2)=psix(2)*y2
  xy(3)=psix(3)*y3
  
  y=xy(1)+xy(2)+xy(3)
  
  return
end subroutine interpol
