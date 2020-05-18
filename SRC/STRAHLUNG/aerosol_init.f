
      subroutine aerosol_init
c
c                        8/14/95, 4/1/97 , 2/10/2000
c
c  **********************************************************************
c  Subroutine to create aerosol optical properties.  There are several
c  inputs and 6 outputs.
c
c    INPUTS FROM COMMON BLOCKS OR HEADER FILE:
c
c    a_tau(nwi) :  The input column aerosol optical depth
c    (real)           (common block "aer_tau" - see header file).
c
c    a_wli(nwi) :  Wavelength in microns corresponding to aerosol tau in "a_tau"
c
c    aprof(# layers): The input aerosol optical depth profile - LAYERS
c    (real)           (common block "aer_prof").
c
c    itp:       Aerosol type, given in header file rad_0598.h.
c
c    ifg:       The table will compute vertical distributions based on
c    (integer)  relative humidity (see explanation below).  If ifg is
c               set to 0, each layer will have properties calculated
c               based on the relative humidity of that layer.  If ifg
c               is set equal to another integer (1 through the number of
c               relative humidities given in the block data "aerosol")
c               the routine will calculate a vertical profile of optical
c               properties based on the relative humidity corresponding
c               to the index given.  The indices are: 1: 0%; 2: 50%;
c               3: 70%; 4: 80%; 5:90%; 6: 95%; 7: 98%; and 8: 99%.
c               If the number of relative humidities changes, these
c               numbers will have to be modified.
c
c    ivd:       Vertical tau distribution flag.  If set to zero, the
c               distribution is based on Jim Spinhirne's marine
c               distribution formulation, and no user input is required.
c               If set to one, the user's own vertical distribution is
c               used, and must be present in the array aprof(nlayers).
c               NOTE: This vertical distribution is used as a weighting
c               factor ONLY, to distribute input column optical depths!
c
c----------------------------------------------------------------------------
c    a_ssa, a_ext, a_asy:  Input single-scattering albedos, extinction
c           coefficients, and asymmetry parameters.  These variables
c           are dimensioned (# of bands, # of relative humidities,
c           # of aerosol types). An x or y is appended on these
c           variable names: if x, the numbers correspond to the 18
c           original bands.  If y, the numbers are for the 10
c           sub-intervals in the first shortwave band (.2-.7 microns).
c           All of these variables come from the block data statements
c           aerosol# (# corresponds to an integer, eg. aerosol1) and
c           are in common blocks aer_optx and aer_opty.
c
c    nv,mb,pp,pt,ph,dz: number of layers, number of bands, and the
c           pressure, temperature, humidity and thickness profiles.
c           These are shared by several subroutines.
c
c    OUTPUTS:
c
c    a_tau1,a_ssa1,a_asy1:  The optical depth, single-scattering albedo,
c       and asymmetry parameter vertical profiles for 18 bands.  These
c       are dimensioned (nvx, 18)  These are in the common block
c       aer_initx, which is shared by the subroutine "aerosolx".
c
c    a_tau2,a_ssa2,a_asy2:  Properties for SW band 1's 10 subintervals.
c       These are dimensioned (nvx, 10)  These are in the common block
c       aer_inity, which is shared by the subroutine "aerosoly".
c
c  **********************************************************************
      USE RadParams
c##      include 'rad_0698.h'
	implicit none
      integer iq,mtop,n,m,ict,ix,iy,irh,krh,iac,itp
      real, dimension(mbx,nrh,naer) :: a_ssax,a_extx,a_asyx
      real, dimension(mby,nrh,naer) :: a_ssay,a_exty,a_asyy

      real, dimension(nvx) :: tauxxx
      real, dimension(nvx,mbx,mxac) :: a_tau1,a_ext1,a_ssa1,a_asy1
      real, dimension(nvx,mby,mxac) :: a_tau2,a_ext2,a_ssa2,a_asy2

      real ,dimension(nvx)  :: taux1,taux2,rh,ht,rhp
      real sumxxx

      real,dimension(mxat) :: a_wli,a_tau
      real,dimension(nvx)  :: aprof
	real,dimension(nvx,mbx) :: wvd_x
	real,dimension(nvx,mby) :: wvd_y

      real pp,pt,ph,po,dz,p1,h1,z,sig,tp
      real rhx(nrh)
      real wts(4),tau3(2),tau3y(4)
      real aotf,wlf,sump,rirh
      real spinhirne_sig, spinhirne_tau

      common /aer_optx/ a_ssax,a_extx,a_asyx
      common /aer_opty/ a_ssay,a_exty,a_asyy
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2
      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /thick/ dz(nvx)

      common /tau_spline_aot/ aotf(15),wlf(15)

      data rhx /0.,50.,70.,80.,90.,95.,98.,99./
      data wts /.23015,.28274,.25172,.23539/


c  Initialize.

	 rh     = -9999.
	 a_ssa1 = 0. ; a_ext1 = 0. ; a_asy1 = 0. ; a_tau1 = 0.
	 a_ssa2 = 0. ; a_ext2 = 0. ; a_asy2 = 0. ; a_tau2 = 0.

	if(nac<0.or.nac>mxac) stop 'nac:# Aerosol Constituents'
	if(n_atau<0.or.n_atau>mxat) stop 'n_atau:# AerosolTau/Wavelengths'
	if(ifg<0.or.ifg>8) stop 'ifg: Aerosol RH% Flag'
	AEROSOL_CONSTITUENTS : do iac = 1,nac

	a_wli(1:n_atau) = a_wlis(1:n_atau,iac)
	a_tau(1:n_atau) = a_taus(1:n_atau,iac)
	aprof(1:nvx)  = aprofs(1:nvx,iac)
	itp	      = itps(iac)
	if ( itp < 1 .or. itp > naer ) stop ' itp : Bad Aerosol Type'
!	print*,'CONSTITUENTS',iac,itp

! FOR Aerosol Optical Properties types that are constant with RH
	if (itp==1  .or. itp==2 .or. itp==3 .or.
     &      itp==10 .or. itp==12 .or.itp==13 .or. itp==18 ) then
!!       Has already been filled in Block data
	else
	do krh=2,8
	 a_extx(1:mbx,krh,itp)= a_extx(1:mbx,1,itp)
	 a_ssax(1:mbx,krh,itp)= a_ssax(1:mbx,1,itp)
	 a_asyx(1:mbx,krh,itp)= a_asyx(1:mbx,1,itp)

	 a_exty(1:mby,krh,itp)= a_exty(1:mby,1,itp)
	 a_ssay(1:mby,krh,itp)= a_ssay(1:mby,1,itp)
	 a_asyy(1:mby,krh,itp)= a_asyy(1:mby,1,itp)

	enddo

	endif
!	if ( ifg .ne.0) print*,'CHECK',ifg,itp,a_ssax(1:mbx,ifg,itp)

c  ******************************************************************
c  Calculate heights at center of layer - find highest layer to place
c  aerosols (15 km) - calculate relative humidities of each layer as
c  needed.  Values of RH > 99% will be set equal to 99% to make table
c  lookup easier. "mtop" is the highest aerosol layer.
c  ******************************************************************
      z=0.
      m=nv
      iq=0
      do while (iq.eq.0)
       ht(m)=(z*2.+dz(m))/2.
       z=z+dz(m)
       if (z.gt.15.) then
         iq=1
         mtop=m
       endif
       p1=(pp(m)+pp(m+1))/2.
       tp=(pt(m)+pt(m+1))/2.
       h1=(ph(m)+ph(m+1))/2.
       call ql_rh(rh(m),tp,p1,h1)
       if (rh(m).gt.98.9) rh(m)=98.9
       if ((rh(m).lt..01).and.(rh(m).gt.-999.)) rh(m)=0.
       m=m-1
       end do

c  *************************************************************
c  Calculate vertical distribution of asymmetry, ss albedo and
c  extinction, based on aerosol type and relative humidity.
c  If ifg is not equal to 0, parameters  will corresponds to a
c  single RH, as described in header file. Loop 31 deals with
c  the 18 original bands, loop 32 with the 10 band 1 subintervals.
c  *************************************************************
      do 30 m=mtop,nv
       do 31 n=1,mbx
        if (rh(m).eq.-9999.) then
          a_ext1(m,n,iac)=-9999.
          a_ssa1(m,n,iac)=-9999.
          a_asy1(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext1(m,n,iac)=a_extx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_extx(n,ict,itp)-a_extx(n,ict-1,itp))
            a_ssa1(m,n,iac)=a_ssax(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_ssax(n,ict,itp)-a_ssax(n,ict-1,itp))
            a_asy1(m,n,iac)=a_asyx(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_asyx(n,ict,itp)-a_asyx(n,ict-1,itp))
	  rhp(m) = rh(m)
          else                        ! Dependence on prescribed RH.
            a_ext1(m,n,iac)=a_extx(n,ifg,itp)
            a_ssa1(m,n,iac)=a_ssax(n,ifg,itp)
            a_asy1(m,n,iac)=a_asyx(n,ifg,itp)
          endif
        endif
 31     continue
!-------------------------------------------
       do 32 n=1,mby
        if (rh(m).eq.-9999.) then
          a_ext2(m,n,iac)=-9999.
          a_ssa2(m,n,iac)=-9999.
          a_asy2(m,n,iac)=-9999.
        else
          if (ifg.eq.0) then          ! Dependence on layer RH.
            ict=2
            do while (rh(m).ge.rhx(ict))
             ict=ict+1
             end do
            a_ext2(m,n,iac)=a_exty(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_exty(n,ict,itp)-a_exty(n,ict-1,itp))
            a_ssa2(m,n,iac)=a_ssay(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_ssay(n,ict,itp)-a_ssay(n,ict-1,itp))
            a_asy2(m,n,iac)=a_asyy(n,ict-1,itp)+(rh(m)-rhx(ict-1))/
     1     (rhx(ict)-rhx(ict-1))*(a_asyy(n,ict,itp)-a_asyy(n,ict-1,itp))
          else                        ! Dependence on prescribed RH.
            a_ext2(m,n,iac)=a_exty(n,ifg,itp)
            a_ssa2(m,n,iac)=a_ssay(n,ifg,itp)
            a_asy2(m,n,iac)=a_asyy(n,ifg,itp)
          endif
        endif
 32     continue

 30    continue

c  ******************************************************************
c  Vertical distribution of aerosol optical depths - CAGEX and CERES.
c       --------------------------------------------------------------
c       Use Spinhirne's vertical distribution of scattering properties
c       to calculate vertical distribution of optical depths.  The
c       distribution gives a scattering coefficient ("sig"). Use this,
c       along with the single-scattering albedo, to produce an
c       RH-dependent extinction coefficient (extx, exty, etc.), from
c       which optical depth is calculated (taux, tauy, etc.).  This
c       optical depth is summed (sum1, sumy2, sum, etc.) to give
c       column tau for weighting purposes.
c       --------------------------------------------------------------

     	select case (ivd)
	case default
	 stop ' ivd : Aerosol Profile flag'
	case (0)  !! DEFAULT VERTICAL DISTRIBUTION Spinhirne

	sumxxx=0.0

        do  m=mtop,nv

	 sig = spinhirne_sig( ht(m))
	 tauxxx(m) = spinhirne_tau(sig,a_ssa2(m,9,iac),dz(m))
	 sumxxx   = sumxxx + tauxxx(m)
!		print*,m,sig,a_ssa2(m,9,iac)
        enddo

	do m=mtop,nv
	 tauxxx(m) = tauxxx(m)  / sumxxx
!!!	 aprofs(m,iac) = tauxxx(m) !! See what the Sphinhirne profiles look like
	enddo

! ----------------------------------------------------------------
	case (1)   ! USER'S OWN VERTICAL DISTRIBUTION IVD=1

         sump =   sum( aprof(mtop:nv) )
	 tauxxx(mtop:nv)= aprof(mtop:nv) / sump

         if(sump.eq.0.)stop 'No VERTICAL Profile OF AEROSOL TAU '

	end select

c  ********************************************************************
c  IAFORM=2
c
c  Distribute optical depth spectrally into the first 2 Fu-Liou bands.
c  Band 1 will consist of the first 4 MFRSR bands, weighted with
c  respect to energy.  Band two will be the fifth MFRSR band.
c
c  Also, distribute optical depths into 4 of the 10 band 1 subintervals.
c  Subinterval 7 is directly inserted, since there is one MFRSR
c  measurement within the range of this band.  Subintervals 7 and 8
c  straddle the .497 micron MFRSR measurement, so interpolated values
c  are inserted into these, using .409 and .497 measurements for 7, and
c  .497 and .606 for 8.  Subinterval 10 contains two MFRSR measurements,
c  so it is filled using an energy-weighted average.  This is all
c  hardwired, so we need all of the MFRSR bands (.409, .497, .606, and
c  .661) for it to work. (The .855 micron band is also needed, but not
c  for this interval distribution.
c  ********************************************************************

	select case ( iaform )
	 case default
	 stop ' iaform : Bad value of iaform '
	case(1)        ! CERES
!! No operations necessary

	case(2)        ! For CAGEX

        tau3(1)=a_tau(1)*wts(1)+a_tau(2)*wts(2)+
     1          a_tau(3)*wts(3)+a_tau(4)*wts(4)
        tau3(2)=a_tau(5)
        tau3y(1)=a_tau(1)      ! For subinterval 7 of 1st band (.409)
        tau3y(2)=a_tau(1)+.6705*(a_tau(2)-a_tau(1)) ! Subi 8 of band 1
        tau3y(3)=a_tau(2)+.4541*(a_tau(3)-a_tau(2)) ! Subi 9 of band 1
        tau3y(4)=a_tau(3)*.5175+a_tau(4)*.4825      ! Subi 10 of band 1

	case(3)        ! For AOT_SPLINEFIT

	if ( ifg == 0 ) then ! Find Aerosol weighted collumn mean RH index
	 rirh=0
	 do m =mtop,nv
	 rirh = rirh + rhp(m)* tauxxx(m)  !! Aerosol Profile weighted mean RH
!	 print*,m,rhp(m),tauxxx(m)
	 enddo

	   irh =1
	   do ix= 1,7
	   if( rirh >= rhx(ix) .and. rirh < rhx(ix+1) ) irh=ix
	   enddo
	   if( rirh >= rhx(8) )irh =8

	else  ! Use assigned RH index
	 irh = ifg
	endif

! Can't handle ZERO in Log interpolation
	where ( a_tau .lt. 1.0E-20) a_tau = 1.0E-20

       call atau_spline_iaform3(mxat,n_atau,a_wli,a_tau,itp,irh)

!	write(22,'(a20,15f8.3)') 'AOT in Fu Bands',aotf(1:15)

!!! ACCOUNT FOR VERTICAL EXTINCTION VARIABILITY WITH HUMIDITY ABOUT THE MEAN RH "irh"
!!! ( IAFORM==3) only
	do iy = 1,mby
	wvd_y(mtop:nv,iy)=tauxxx(mtop:nv)
     &             *a_ext2(mtop:nv,iy,iac)/a_exty(iy,irh,itp)
	sump =   sum( wvd_y(mtop:nv,iy) )
	wvd_y(mtop:nv,iy) =  wvd_y(mtop:nv,iy) /sump
	enddo

	do ix = 1,mbx
	wvd_x(mtop:nv,ix)=tauxxx(mtop:nv)
     &             *a_ext1(mtop:nv,ix,iac)/a_extx(ix,irh,itp)
	sump =   sum( wvd_x(mtop:nv,ix) )
	wvd_x(mtop:nv,ix) =  wvd_x(mtop:nv,ix) /sump
	enddo

	end select


! ----------------------------------------------------------------
c       Use weighted optical depths  to distribute our input
c       column optical depths vertically and spectrally where needed.
c       For bands with "measured" input, we simply do the weighting.
c       For the remaining bands, we weight according to our vertically
c       distributed extinction coefficients (calculated in loop 30),
c       which carry all the spectral resolution we need.  a_tau1 is for
c       the 18 original bands, a_tau2 is for the 10 band 1 subintervals.
! ----------------------------------------------------------------
        VERTICAL : do  m=mtop,nv

	select case ( iaform )

	case(1)       ! For CERES

           a_tau1(m,1,iac)   = a_tau(1) * tauxxx(m)
           a_tau1(m,2:18,iac)= a_tau1(m,1,iac)*
     &                 a_ext1(m,2:18,iac)/a_ext1(m,1,iac)

           a_tau2(m,9,iac)  = a_tau(1) * tauxxx(m)

           a_tau2(m,1:10,iac)=a_tau2(m,9,iac)*
     &                 a_ext2(m,1:10,iac)/a_ext2(m,9,iac)

        case(2)        ! For CAGEX

	    a_tau1(m,1:2,iac) = tau3(1:2) * tauxxx(m)
            a_tau1(m,3:18,iac)=a_tau1(m,2,iac)*
     &                a_ext1(m,3:18,iac)/a_ext1(m,2,iac)

	    a_tau2(m,7:10,iac) = tau3y(1:4) * tauxxx(m)
            a_tau2(m,1:6,iac)  = a_tau2(m,7,iac)*
     &                a_ext2(m,1:6,iac)/a_ext2(m,7,iac)

	 case(3)       ! For AOT_SPLINEFIT



!	a_tau2(m,1:10,iac) = aotf(1:10)  * tauxxx(m)
	a_tau2(m,1:10,iac) = aotf(1:10)  * wvd_y(m,1:10)
!	a_tau1(m,1,iac)    = aotf(9)     * tauxxx(m)
	a_tau1(m,1,iac)    = aotf(9)     * wvd_x(m,1)
!	a_tau1(m,2:6,iac)  = aotf(11:15) * tauxxx(m)
	a_tau1(m,2:6,iac)  = aotf(11:15) * wvd_x(m,2:6)
            a_tau1(m,7:18,iac) =a_tau1(m,2,iac)*
     &                 a_ext1(m,7:18,iac)/a_ext1(m,2,iac)

         end select

!	print'(3I4,2f8.2,16f7.3)', m,iac,itp,dz(m),rh(m),
!     & (wvd_y(m,iy),iy=1,10),(wvd_x(m,ix),ix=1,6)

        enddo VERTICAL

!------------------------------------------------------------------------------
!!!--- Diagnostic Output of Atau
!	do ii=1,10
!	xxx=0
!	 do jj=1,nv
!	 xxx =xxx+ a_tau2(jj,ii,iac)
!	 enddo
!	aotf(ii)=xxx
!	enddo

!	do ii=2,6
!	xxx=0
!	 do jj=1,nv
!	 xxx =xxx+ a_tau1(jj,ii,iac)
!	 enddo
!	aotf(9+ii)=xxx
!	enddo

!	write(22,'(a20,15f8.3)') 'AOT in Fu Bands',aotf(1:15)

	enddo AEROSOL_CONSTITUENTS

      return
      end

!===========================================================================
      subroutine aerosolxy ( ib,cmode )
c *********************************************************************
c                      Modified 2/14/00
c
c tae, wae, and wwae are the optical depth, single scattering albedo,
c and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
c due to the Mie scattering of aerosols for a given layer.
c
c  This subroutine is called for bands 2 - 18 (ib)
c  or vis subbands 1-10 (ig)
c *********************************************************************
      USE RadParams
      implicit none
	character*1 cmode
      integer i,ib,iac
      real x1,x2,x3,x4,y1,y2,y3,y4,tae,wae,wwae
      real ,dimension(nvx,18,mxac) :: a_tau1,a_ssa1,a_asy1
      real ,dimension(nvx,10,mxac) :: a_tau2,a_ssa2,a_asy2
      common /aer_initx/ a_tau1,a_ssa1,a_asy1
      common /aer_inity/ a_tau2,a_ssa2,a_asy2

      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)

      AEROSOL_CONSTITUENTS  : do iac=1,nac

      LEVELS : do  i = 1, nv
       select case (cmode)
	case ('x')
         tae(i,iac) = a_tau1(i,ib,iac)
         wae(i,iac) = a_ssa1(i,ib,iac)
         x1     = a_asy1(i,ib,iac)
	case ('y')
         tae(i,iac) = a_tau2(i,ib,iac)
         wae(i,iac) = a_ssa2(i,ib,iac)
         x1     = a_asy2(i,ib,iac)
       end select

       x2 = x1 * x1
       x3 = x2 * x1
       x4 = x3 * x1
       y1 = 3.0 * x1
       y2 = 5.0 * x2
       y3 = 7.0 * x3
       y4 = 9.0 * x4

       wwae(i,1,iac) = y1
       wwae(i,2,iac) = y2
       wwae(i,3,iac) = y3
       wwae(i,4,iac) = y4

	enddo LEVELS
	enddo AEROSOL_CONSTITUENTS

      return
      end
!----------------------------------------------------------------
	real function spinhirne_sig(ht)

	data sig0,a,ap,b,bp,f/0.025,0.4,2981.0,1.6,2.5,1.5e-7/

         t1=  sig0*(1+a)**2
         t4 = f*(1+ap)**2

         t2 = exp(ht/b)
         t3 = (a+exp(ht/b))**2
         t5 = exp(ht/bp)
         t6 = (a+exp(ht/bp))**2
         spinhirne_sig=t1*t2/t3+t4*t5/t6   ! scattering coefficient

	return
	end
!---------------------------------------------
	real function spinhirne_tau(sig,ssa,dz)
	ext = sig / ssa
	spinhirne_tau = ext / dz
	return
	end
