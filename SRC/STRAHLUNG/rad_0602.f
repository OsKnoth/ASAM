C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -c 6-24-98 (1)
      subroutine rad ( as,as1, u0, ss, pts, ee, ur, ITriErr )
c +++ PAUSE FIX      subroutine rad ( as, u0, ss, pts, ee , ur )
c 6-24-98 (1)
c **********************************************************************
c MODIFIED 0598 Tho compute Window RADIANCES ( ONLY WITH 2-Stream OPTION)
c MODIFIED 0598 A) Outputs Window (800-1250cm-1) Flux Profiles., (F. Rose) 
c		B) Includes new Window K's from Dave Kratz,  (F. Rose) 
c 		C) H2O continuum absorbtion extended to (5-2200cm-1) 
c		   according to CKD_2.1, (F. Rose) 
c MODIFIED 4/2/97  to include CKD continuum for LW and associated 
c                  parameterizations, (F. Rose) 
c MODIFIED 4/1/97 to apportion 1st SW band into 10 sub-intervals 
c                  for aerosols. Removed commented out sections. 
c                  from older versions (T. Alberta)
c MODIFIED 10/26/96 to include aerosols and direct/diffuse. (T. Alberta)
c
c MODIFIED 10/25/96 - changed logicals for 2 and 4 stream
c                     options (T. Alberta)
c MODIFIED 10/23/96 - many cosmetic changes. (T. Alberta)
c MODIFIED 10/23/96 - many cosmetic changes. (T. Alberta)
c
c MODIFIED 9/96 to include variable levels after compilation. (F. Rose)
c
c In this radiation scheme,  six  and  12 bands are selected for solar 
c and thermal IR regions, respectively. The spectral division is below: 
c 0.2 - 0.7 um, 0.7 - 1.3 um, 1.3 - 1.9 um, 1.9 - 2.5 um, 2.5 -3.5 um,
c 3.5 - 4.0 um, and 2200 - 1900 cm**-1, 1900 - 1700 cm**-1, 1700 -1400
c cm**-1,  1400 - 1250 cm**-1,  1250 - 1100 cm**-1, 1100 - 980 cm**-1,
c 980 - 800 cm**-1,  800 - 670 cm**-1,  670 - 540 cm**-1, 540 - 400 cm
c **-1,  400 - 280 cm**-1,  280 - 0 cm**-1,  where  the index  for the
c spectral band ( ib = 1, 2, ..., 18 ) is defined.
c
c                                                 ***********
c                                                 * Common  *
c             **********************              * block   *  ******
c             *  INPUT PARAMETERS  *              * or call *  *TYPE*
c             **********************              ***********  ******
c   as(mbs)   Surface albedo, mbs = 6                CALL       REAL
c   as1(10)   Surface albedo in 10 vis_subbands      CALL       REAL
c             (Spectral reflectances for SW bands)
c   u0        cosine of solar zenith angle           CALL       REAL
c   ss        solar constant (W/m**2)                CALL       REAL
c   pts       surface temperature (K)                CALL       REAL
c   ee(mbir)  IR surface emissivity, mbir=12         CALL       REAL
c   pp(nv1x)  atmospheric pressure (millibars)      atmos       REAL
c   pt(nv1x)  atmospheric temperature (K)           atmos       REAL
c   ph(nv1x)  water vapor mixing ratio (kg/kg)      atmos       REAL
c   po(nv1x)  ozone mixing ratio (kg/kg)            atmos       REAL
c   pre(nvx)  effective radius of water cloud (um)  clouds      REAL
c   plwc(nvx) liquid water content (g/m**3)         clouds      REAL
c   pde(nvx)  effective diameter of ice cloud (um)  clouds      REAL
c   piwc(nvx) ice water content (g/m**3)            clouds      REAL
c   prwc(nvx) rain water content (g/m**3)           rains       REAL
c   pgwc(nvx) graupel water content (g/m**3)        graups      REAL
c   umco2     concentration of CO2 (ppmv)           umcon       REAL
c   umch4     concentration of CH4 (ppmv)           umcon       REAL
c   umn2o     concentration of N2O (ppmv)           umcon       REAL
c   fourssl   if true, run four stream for solar    tsfslog    LOGICAL
c             if false, two stream is run
c   foursir   if true, run four stream for IR       tsfslog    LOGICAL
c             if false, run two/four stream for IR  
c   nv        number of LAYERS                      levels     INTEGER
c   itp       aerosol type (1, 2, 3 - see header)   aer_tau    INTEGER  
c   ivd       Aerosol optical depth vertical        aer_tau    INTEGER
c             distribution (see header)
c   ipr       Aerosol properties flag (see header)  aer_tau    INTEGER  
c   irobckd   LW continuum option                   irobckd    INTEGER 
c             1: Roberts continuum
c             2: Exact CKD_2.1 continuum
c             3: No continuum
c             4: Parameterized CKD_2.1 continuum
!             5: Parameterized CKD_2.4 continuum 
!  nhb	      #of "hidden bands for the thermal Longwave > 2200cm-1
!	      0: none - operates as older versions
!	      1: 0-2500cm-1 
!	      2: 0-2850cm-1   
c   a_tau(5)  MFRSR-derived aerosol optical depths  aer_tau     REAL
c             (For CAGEX only)
c   a_taux    Column aerosol optical depth          aer_tau     REAL
c             (For CERES only)
c
c NEW FOR 0598  NEW K's For WINDOW band 11,12,13 
c
c   idkfr     Window K's Option			   /dkfrwn/     INTEGER  
c	      0: orig Fu_liou 
c	      1:Kratz Window K's :: SLOW but ACCURATE
c             2: Hybred          :: FAST   SUGGEST idkfr=2
c
c H20 Continuum absorbtion now for entire LW (5-2200cm-1)
c    iwtas    1:Linear				    /cont_tas/   INTEGER
c             2:log ,
c             3:linear w/Plank wgt   :: SUGGEST iwtas=3
c             4: Log w/plank wgt 
!     OBSOLETE for irobckd =5  CKD_2.4 continuum linear w/Plank wgt is used
c Fu suggests use of linear averaging of continuum spectral Taus weighted by
c plank function 
c
c cfc_conc(3)  Window K's now include CFC's         /cfcs/      REAL
c	      ( F11,F12,F22 ) in ppv
c
c Suggested data cfc_conc/0.268e-09 , 0.503e-09 ,0.105e-09/
c  
c NEW FOR 0698 RADIANCE
c Input:
c	ur - Cosine of View Zenith Angle 
c Output:: 
c	fiurt(nv1x)  - TOTAL LW (0-2200cm-1) radiance [Wm-2sr-1]
c	fiurw(nv1x)  - WINDOW (800-1250cm-1) radiance [Wm-2sr-1]
c 
c   Check these common blocks in this program and in the header file
c   for proper format.
c   
c
c Note:  (1)  as(mbs) and ee(mbir) consider the substantial wavelength
c             dependence of surface albedos and emissivities.
c        (2)  For CO2, CH4 and N2O, uniform mixing is assumed  through
c             the atmosphere. The  concentrations  can be changed
c             through 'common /umcon/ umco2, umch4, umn2o '.
c        (3)  nvx, nv1x, ndfsx, mdfsx, ndfs4x, mbx, mbsx, mbirx,  
c             and  ncx  are given through 'rad_0598.h'.  These
c             variables are used for dimensioning only!
c        (4)  nv1 and 1 are the surface and top levels, respectively.
c
c                       **********************
c                       *  OUTPUT PARAMETERS  *
c                       **********************
c              fds(nv1x)   downward solar flux ( W / m ** 2 )
c              fus(nv1x)   upward solar flux ( W / m **2 )
c              dts(nvx)    solar heating rate ( K / day )
c              fdir(nv1x)  downward IR flux ( W / m ** 2 )
c              fuir(nv1x)  upward IR flux ( W / m **2 )
c              dtir(nvx)   IR heating rate ( K / day )
c              fd(nv1x)    downward net flux ( W / m ** 2 )
c              fu(nv1x)    upward net flux ( W / m **2 )
c              fdsdr(nv1x) downward direct SW flux ( W / m **2 )
c              fdsdf(nv1x) downward diffuse SW flux ( W / m **2 )
c              dt(nvx)     net heating rate ( K / day )
c              fdwn(nv1x)  downward WINDOW flux ( W / m ** 2 )
c              fuwn(nv1x)  upward WINDOW flux ( W / m **2 )
c
c Note:  Solar, IR, and net represent 0.2 - 0.4 um, 2200 - 0 cm**-1,
c        and  entire spectral regions, respectively.
c 	 Window represents (800-1250cm-1)
c 6-24-98 (2)	
c                  ********************************* 
c                  * OUTPUT FOR IR WINDOW RADIANCE *
c                  *********************************
c              fiurw(nv1) upward radiance at ur in IR window (W/m**2/Sr)
c              fiurt(nv1) upward radiance at ur in TOTAL LW (W/m**2/Sr)
c
c Note: The IR window is defined as the spectral interval between
c       800 to 1250 cm**-1.
c 6-24-98 (2)
c
c Fu 07-08-98
c The improved parameterization of cirrus radiative properties in
c both solar and IR spectra (Fu 1996; Fu et al. 1998) has been 
c incorporated into the radiation model.  Note that the definition
c of the generalized effective size used in the new parameterization
c (Eq. 3.10 or Eq. 2.3 in Fu 1996) is different from the mean
c effective size defined in Eq. 2.1 of Fu and Liou (1993).  Now 
c you can make choice between the two versions of cirrus para-
c meterization through the logical variable "fl93i".  Use appropriate
c effective sizes of ice clouds through input "pde" for the two
c different versions of parameterization.
c Fu 07-08-98
c
c *********************************************************************
	use TRMM_WINDOW
c **********************************************************************
      USE RadParams
      implicit none
C##      include 'rad_0698.h'
c Fu 07-08-98
        logical fl93i
c Fu 07-08-98
      integer i,ib,ig,mbn,kg(mbx+2),iac
      real pp,pt,ph,po,pre,plwc,pde,piwc,prwc,pgwc,umco2,umch4,umn2o
      real fds,fus,dts,fdir,fuir,dtir,fd,fu,dt,fu1,fd1,bf,bs
      real as(mbsx),as1(12), ee(mbirx),u0,ss,pts,f0,fuq1,fuq2,xx,asx
      real hk

!BEGIN  Fu1001
! Add "hidden bands" 19 & 20 in SOLAR/IR Overlap region
! hidden band 19 = 2850:2500cm-1
! hidden band 20 = 2500:2200cm-1
! Use 1-Albedo of band 6 as sfc emissivity in hidden band 19.
! emissivity of band 7 used in hidden band 20.

! Cloud and Aerosol Properties taken from band 7 1800:2200cm-1
! Handles extended  Plank Source Function 
! Gas absorption H20  2200:2500cm-1(Kratz Correlated-K).
! Gas absorption H20/CO2/N2  2200:2500cm-1(Kratz Correlated-K).
! Continuum absorption for 2500:2850cm-1 set to zero.

	integer ibt
	real eeband
!END  Fu1001



! Chou arrays :5/99
	integer npc,ntop,nbot
	real dn1(nv1x),up1(nv1x),u0d
	data u0d / 0.50 /
! Chou arrays :5/99

c---------- 10/29/96 (1)
c---------- 4/1/97 (1)
      real ctau, dtau, tae, wae, wwae 
	real aprop
      common /aer_proc/ ctau(18), dtau(10)
c---------- 4/1/97 (1)
      common /aer/ tae(nvx,mxac), wae(nvx,mxac), wwae(nvx,4,mxac)
      common /aerpout/aprop(nvx,mxac,3)
c---------- 10/29/96 (1)

c---------- 10/28/96 (1)
      real fdsdr,fdsdf,ffdr,ffdf
      common /radiat/ fds(nv1x), fus(nv1x), dts(nvx),
     1                fdir(nv1x), fuir(nv1x), dtir(nvx),
     1                fd(nv1x), fu(nv1x), dt(nvx),  
     1                fdsdr(nv1x), fdsdf(nv1x)
      common /dirdiff/ffdr(nv1x),ffdf(nv1x)

c---------- 10/28/96 (1)

      common /atmos/ pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
      common /clouds/ pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
      common /rains/ prwc(nvx)
      common /graups/ pgwc(nvx)
      common /umcon/ umco2, umch4, umn2o 
      common /dfsout/ fu1(nv1x), fd1(nv1x)
      common /planci/ bf(nv1x), bs

	real otau,tw,ti
      	common /outtau/ otau(mbx) !! 7-21-96
	common /wat/ tw(nvx)
	common /ic/  ti(nvx)


!!! NEW FOR 0598 ----------------------------------------------------------


  	real cfc_conc(3)
	real fuwn,fdwn
	integer idkfr,iwtas
c +++ PAUSE FIX 
        INTEGER ITriErr
c
        common /wndow/ fuwn(nv1x),fdwn(nv1x)
	common /dkfrwn/ idkfr  
	common /cfcs/ cfc_conc
        common /cont_tas/ iwtas
c 6-24-98 (4a)	
	real fiurt,fiurw,fiur,ur
  	common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x) 
	real trwn_flt_r, trwn_unf_r, trwn_f
	common /TRMMWNOUT/ trwn_flt_r, trwn_unf_r, trwn_f 
c 6-24-98 (4a)
	integer isolar_spectrum
	real hk1, fk1o3,sol_spect,fk1h2o

        common /band1/ hk1(10), fk1o3(10),sol_spect(0:7),fk1h2o(10)
 	common /select_solar_spectra/ isolar_spectrum

	logical lchou,lband6a,lpar,lray
	common/CHOU/lchou,lband6a,lpar,lray

c   uvfu(nv1,10) upward flux within 10 UV-VIS subintervals (W/m**2)
c   uvfd(nv1,10) downward flux within 10 UV-VIS subintervals (W/m**2)
      real uvfu,uvfd,uvdir,uvdif
      common /uvflux/ uvfu(nv1x,10),uvfd(nv1x,10),
     &               uvdir(nv1x,10),uvdif(nv1x,10)
	real swfu,swfd,swdir,swdif
	common /swflux/ swfu(nv1x,6),swfd(nv1x,6),
     &               swdir(nv1x,6),swdif(nv1x,6)
	real rlwfu,rlwfd,sbf,sbs
	common /lwflux/ rlwfu(nv1x,7:20),rlwfd(nv1x,7:20),
     & sbf(nv1x,7:20),sbs(7:20)
	integer isksw
	
	common /seijik/ isksw
c  ******************************************************************
c  kg(mb+2) is the number of intervals to perform the g-quadrature in
c  each band to consider the nongray gaseous absorption.  In total,
c  we need to perform "sum(kg)" spectral calculations in  the  scattering
c  problem for each atmospheric profile.
c  ******************************************************************
      data kg / 10, 8, 12, 7, 12, 5, 
     1            2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8,
     1            5,5 /


!! Depending on K's option # of Solver loops needed changes...
	if ( idkfr == 0) kg(11:13) = (/3,5,2/)
	if ( idkfr == 1) kg(11:13) = (/80,40,10/)
	if ( idkfr == 2) kg(11:13) = (/5,5,2/)

!!	isksw=0
	if ( isksw > 0) then
	 kg(2:6) = (/7,8,7,8,7/) !ht02 sav
  	 lchou = .false.
	endif
c  *********************************************************************
c  The following variables are declared here, instead of in the calling
c  program, since any change in these would require modifications
c  the code anyway.  A check is inserted here, to make sure the 
c  number of layers given by the user, nv, is less than or equal
c  to the number of levels nvx, given in the header file rad_0598.h.
c  nvx is used for dimensioning arrays.  Uncomment to use.
c  *********************************************************************
      nv1=nv+1        ! Number of levels
      mb=18           ! Number of bands
      mbs=6           ! number of shortwave bands
      mbir=12         ! number of longwave bands
      nc=8            ! number of cloud types
      ndfs=nv         ! number of layers
      mdfs=nv1        ! number of levels
      ndfs4=4*ndfs    ! number of layers * 4
      ndfs2=2*ndfs    ! number of layers * 2



cc      if (nv.gt.nvx) then
cc        print *,'ARRAY ERROR: number of levels specified (nv) '
cc        print *,'is greater than number allowed (nvx)!!!'
cc        stop
cc      endif
cc      ITriErr = 0
c  ******************************************************************
c  ******************** 10/23/96 (end) ******************************
c  ******************************************************************

c Fu 07-08-98

 	fl93i = .false. !! False = New 98Ice cld 
!	fl93i = .true.  ! True = Old93 ice
	if ( fl93i ) then
! Nothing
	else
! Modify Ice particle size Input for FU98ICE is Generalized
! Effective diameter (Dge) while input for FU93ICE is Mean
! Effective diamter (De).
! This Parmeterization allows input of (De) in this case
	do i=1,nv	
	 if ( pde(i) > 0)
     &	 pde(i) = -2.4+0.7566*pde(i)+9.22E-04*pde(i)*pde(i)
	enddo
	endif

c Fu 07-08-98
      f0 = 1.0 / 3.14159
      do 10 i = 1, nv1
       fds(i) = 0.0
       fus(i) = 0.0
       fdir(i) = 0.0
       fuir(i) = 0.0
	fuwn(i)=0.0
	fdwn(i)=0.0
c---------- 10/28/96 (2)
       fdsdr(i) = 0.0        
       fdsdf(i) = 0.0  
       if (i.lt.nv1) then
         dts(i)=0.      
         dtir(i)=0.      
         dt(i)=0.  
       endif    
c---------- 10/28/96 (2)
c 6-24-98 (5)	
           fiurw(i) = 0.0
	   fiurt(i) = 0.0
c 6-24-98 (5)	
10     continue

	swfu=0.0
	swfd=0.0
	swdir=0.0
	swdif=0.0
	rlwfu=0.0
	rlwfd=0.0
	sbf=0.0 ;sbs=0
      call thicks
      call rayle2

      if ( u0 .le. 1.0e-4 ) then
        mbn = mbs + 1    ! Longwave only if sun is too low
      else
        mbn = 1          ! Shortwave and longwave
      endif


c---------- 10/29/96 (2)
      call aerosol_init
c---------- 10/29/96 (2)

!TRMM (begin)
	fu_sr=0.0
	fu_sf=0.0
!TRMM (end)

	call ckd1_init( isolar_spectrum,lray )


 	do 20 ibt = mbn,mb+nhb
!! WATCH CAREFULLY the use of ib and ibt
! ibt (1-20) band index over added hidden bands 
! ib (1-18) band index over reported bands
	ib=ibt
	if ( ibt == 19  ) ib =7 !!!! (19)2850:2500cm-1 
	if ( ibt == 20  ) ib =7 !!!! (20)2500:2200cm-1

           if ( fl93i ) then
	      call ice ( ib )
	   else
              call icenew ( ib )
           endif


! Use Yong Hu's Water cloud optical properties for SW bands
	if ( ib <=6 ) then
	  call water_hu ( ib ) ! CALL Yong Hu's WATER CLOUD OPTICS For SW bands ib=1:6
	else
 	  call water ( ib )
	endif

       call rain ( ib )
       call graup ( ib )

c---------- 4/1/97 (3)
! No more ipr option
       if (ib.ne.1) then
       call aerosolxy (ib,'x')

	if ( ibt <= 18 ) then
         ctau(ib)=0.
         do i=1,nv
	 do iac=1,nac
         ctau(ib)=ctau(ib)+tae(i,iac)
	 end do
         end do
	endif

       endif
c---------- 4/1/97 (3)

c --- 7-21-97 -- Vertically Integrated Cloud Water/Ice optical Depth
	if ( ibt <= 18 ) then
	otau(ib)=0.0
	do i=1,nv
        otau(ib)=otau(ib)+(tw(i)+ti(i))
	enddo
	endif

c --- 7-21-97
	if (ib ==1 .and. lchou ) then
! Chou scheme3 
	ntop=0
	nbot=0
	npc=0
	 do i=1,nv
	  if  ( (tw(i)+ti(i)) .gt. 0.0 )then
	   ntop =i
	   npc= i
	   exit
	  endif
	 enddo
	 do i=nv,1,-1
	  if  ( (tw(i)+ti(i)) .gt. 0.0 )then
	   nbot =i+1
	   npc= i
	   exit
	  endif
	 enddo
! Chou scheme3
	endif

! MOVED i       call rayle ( ib, u0  )
c---------- 4/2/97 (2)
c      call gascon ( ib )  ! REPLACED WITH CONDITIONAL BELOW.
c
      if (irobckd .eq.1) then
        call gascon ( ib )  !! OLD ROBERTS CONTINUUM
      elseif (irobckd .eq.2) then
        call gascon_ckd ( ib ) !! EXACT CKD_2.1 (SLOW)
      elseif (irobckd .eq.3) then
        call gascon_off      !! Option to Turn off Continuum Absorbtion
      elseif (irobckd .gt. 3 ) then
        call gascon_ckd_parm(ib) !! Parameterized CKD_2.1 Cont. Abs.
      else
        stop 
     1  'Set Continuum Model 1=Roberts 2=CKD_2.1 3=NONE, 4=PARM CKD_2.1'
      endif
			
c---------- 4/2/97 (2)

       if ( ibt .gt. mbs ) then
	call planck ( ibt, pts )
!	print*,'Planck',ibt,ib,bs*3.1415
      endif
		
	if ( ibt > 18 ) call gascon_off ! No H20 Continuum Absorption > 2200cm-1

		
       do 30 ig = 1, kg(ibt)
c---------- 4/1/97 (4)
! No more ipr option
        if (ib.eq.1) then
          call aerosolxy (ig,'y')

	if(ig ==9) then
	 do i=1,nv
	 do iac=1,nac
	 aprop(i,iac,1)= tae(i,iac)
	 aprop(i,iac,2)= wae(i,iac)
	 aprop(i,iac,3)= wwae(i,1,iac) *0.3333
	 enddo;enddo
	endif

          dtau(ig)=0.
          do i=1,nv
	   do iac=1,nac
           dtau(ig)=dtau(ig)+tae(i,iac)
	   end do
           end do
        endif
c---------- 4/1/97 (4)
	call rayle ( ib, u0, ig , lray)
        

	if ( isksw > 0 .and. ( (ibt>=2 .and. ibt<=6) 
     & .or. (ibt==1 .and. ig >=9) )  ) then
	 call seijifu_ht02a_sav( ibt, ig, hk )
	else
 	 call gases ( ibt, ig, hk )
	endif

	
        call comscp 

c  ---------------------------
c  10/25/96 -- 11/4/95 (begin)
c  ---------------------------
        if ( ib .le. mbs ) then 

	asx = as(ib)
    !   write(*,*) 'strahlung 0602: ig',ig,as1(ig)
!! If as1 array is filled with good albedo for band 1 use it !
	if ( ib == 1        .and. 
     &     as1(ig) .ge. 0.0 .and.
     &     as1(ig) .le. 1.0 ) asx = as1(ig)

          if ( fourssl ) then
		
            call qfts ( ib, asx, u0, f0 )  !! FOUR STREAM SOLIR
		
          else
            quadra = .false.
            hemisp = .false.
            edding = .true.
c +++ PAUSE FIX           call qftsts ( ib, asx, u0, f0 ) ! TWO STREAM SOLIR
            call qftsts ( ib, asx, u0, f0, ITriErr ) ! TWO STREAM SOLIR
	            IF (ITriErr .NE. 0) RETURN
          endif
	
          do 40 i = 1, nv1
           fds(i) = fds(i) + fd1(i) * hk
           fus(i) = fus(i) + fu1(i) * hk
c---------- 10/28/96 (3)
           fdsdr(i) = fdsdr(i) + ffdr(i) * hk
           fdsdf(i) = fdsdf(i) + ffdf(i) * hk
c---------- 10/28/96 (3)
	if ( ib == 1 ) then
	 uvfu(i,ig)  = fu1(i)  * hk
	 uvfd(i,ig)  = fd1(i)  * hk
	 uvdir(i,ig) = ffdr(i) * hk
	 uvdif(i,ig) = ffdf(i) * hk
	endif
	if ( ib <= 6 ) then
	 swfu(i,ib)  = swfu(i,ib)  +fu1(i)  * hk
	 swfd(i,ib)  = swfd(i,ib)  +fd1(i)  * hk
	 swdir(i,ib) = swdir(i,ib) +ffdr(i) * hk
	 swdif(i,ib) = swdif(i,ib) +ffdf(i) * hk
	endif
 40        continue

        else !!!THERMAL
!-------------------------------------------------------------------
	eeband =   ee(ib-mbs)
	if (ibt == 19 ) eeband = 1.0-as(6)
	
        if ( foursir ) then
            call qfti ( ib, eeband ) !FOUR STREAM IR
          else
            quadra = .false.
            edding = .false.
            hemisp = .true.
            mquadr = .false.
c 6-24-98 (6)	
c                  call qftisf ( ib, eeband )!TWO-FOUR STREAM COMB IR
c +++ PAUSE FIX                   call qftisf ( ib, eeband, ur )
                   call qftisf ( ib, eeband, ur, ITriErr )
c 6-24-98 (6)

c +++ PAUSE FIXc               call qftits ( ib, eeband ) !TWO-STREAM IR
c               call qftits ( ib, eeband, ITriErr ) !TWO-STREAM IR
          endif
c  -------------------------
c  10/25/96 -- 11/4/95 (end)
c  -------------------------

          do i = 1, nv1
           fdir(i) = fdir(i) + fd1(i) * hk
           fuir(i) = fuir(i) + fu1(i) * hk
	   fiurt(i) = fiurt(i) + fiur(i) * hk * 0.1591549
           end do

	if ( ibt >= 7 .and. ibt <=20) then
   	 do i = 1, nv1
           rlwfu(i,ibt) =  rlwfu(i,ibt) + fu1(i) * hk
           rlwfd(i,ibt) =  rlwfd(i,ibt) + fd1(i) * hk
	   sbf(i,ibt) = bf(i)
	   sbs(ibt) = bs
  	enddo
	endif

	if ( ib.eq.11 .or. ib.eq.12 .or. ib.eq.13) then
		do i=1,nv1
	        fdwn(i) = fdwn(i) + fd1(i) * hk
                fuwn(i) = fuwn(i) + fu1(i) * hk
c 6-24-98 (7)	
           fiurw(i) = fiurw(i) + fiur(i) * hk * 0.1591549
c 6-24-98 (7)
		enddo
!TRMM (begin)
! GET Spectral Flux and Radiance for Filtering to TRMM WINDOW
          fu_sf(ib) = fu_sf(ib) + fu1(1)  * hk
          fu_sr(ib) = fu_sr(ib) + fiur(1) * hk * 0.1591549
!TRMM (end)      
	endif 


        endif

30      continue  
20     continue

c  ------------------------------------------------------------------
c  In this model, we used the solar spectral irradiance determined by
c  Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
c  in the spectral region 0.2 - 4.0 um.
c
c  fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
c  hk of 0.03.
c  ------------------------------------------------------------------
      if ( lband6a ) then
!	fuq1 = ss / ( 1340.0 +11.0 )
	fuq1 = ss /  sol_spect(0)  
	else
	fuq1 = ss / ( sol_spect(0) - sol_spect(7) )
       endif


c       fuq2 = bs * 0.03 * 3.14159 * ee(12)
c fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
c hk of 0.03.
        fuq2 = 0.0

      do 60 i = 1, nv1
       fds(i) = fds(i) * fuq1
       fus(i) = fus(i) * fuq1
       fuir(i) = fuir(i) + fuq2
       fd(i) = fds(i) + fdir(i)
       fu(i) = fus(i) + fuir(i)
c---------- 10/28/96 (4)
       fdsdr(i) = fdsdr(i) * fuq1
       fdsdf(i) = fdsdf(i) * fuq1
c---------- 10/28/96 (4)
60     continue

!! Correct Spectral Output
	 uvfu  = uvfu  * fuq1
	 uvfd  = uvfd  * fuq1
	 uvdir = uvdir * fuq1
	 uvdif = uvdif * fuq1
	
	 swfu  = swfu  * fuq1
	 swfd  = swfd  * fuq1
	 swdir = swdir * fuq1
	 swdif = swdif * fuq1

!---------------------------------------------------------------------
! Account for SW absorption by CO2 & O2 after CHOU : Journal Climate Feb90 & CLIRAD-SW NASA TM
	if (lchou ) then

	do i =1,nv1
	 if ( fdsdr(i)/fds(i)  < 0.50 ) then
	  npc=i 
	if ( ntop==0 .or. nbot ==0 ) then
	  ntop = npc 
	  nbot = npc+1
	endif
	  exit
	 endif
	enddo

	if ( npc .ne.0 .and. npc < ntop) npc = ntop
	if ( npc .ne.0 .and. npc > nbot) npc = nbot

	call add_chou3(nv,umco2,u0,pp,pt,ph,fds,fus,dn1,up1,npc,ntop,nbot)

! Here we assume ALL of the CO2 & O2 absorption comes from the Direct Beam
!	fdsdr(1:nv1)=fdsdr(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )

! Here we assume the CO2 & O2 absorption is split evenly between Direct & diffuse.
	fdsdr(1:nv1)=fdsdr(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )* 
     &  fdsdr(1:nv1)/fds(1:nv1)
	fdsdf(1:nv1)=fdsdf(1:nv1)- (fds(1:nv1)-dn1(1:nv1) )* 
     &  fdsdf(1:nv1)/fds(1:nv1)

! ADJUST SPECTRAL DIRECT  FlUX by CHOU ( VERY COARSE APROXIMATION)
	swdir(1:nv1,2)=swdir(1:nv1,2)-(fds(1:nv1)-dn1(1:nv1))*0.550* 
     &  fdsdr(1:nv1)/fds(1:nv1)
	swdir(1:nv1,4)=swdir(1:nv1,4)-(fds(1:nv1)-dn1(1:nv1))*0.450* 
     &  fdsdr(1:nv1)/fds(1:nv1)
! ADJUST SPECTRAL DIFFUSE  FlUX by CHOU ( VERY COARSE APROXIMATION)
	swdif(1:nv1,2)=swdif(1:nv1,2)-(fds(1:nv1)-dn1(1:nv1))*0.550* 
     &  fdsdf(1:nv1)/fds(1:nv1)
	swdif(1:nv1,4)=swdif(1:nv1,4)-(fds(1:nv1)-dn1(1:nv1))*0.450*
     &  fdsdf(1:nv1)/fds(1:nv1)


! Spectral SW UP
	swfu(1:nv1,2) = swfu(1:nv1,2) - (fus(1:nv1)-up1(1:nv1) ) *0.550 
	swfu(1:nv1,4) = swfu(1:nv1,4) - (fus(1:nv1)-up1(1:nv1) ) *0.450 

! Spectral SW DN
	swfd(1:nv1,2) = swfd(1:nv1,2) - (fds(1:nv1)-dn1(1:nv1) ) *0.550 
	swfd(1:nv1,4) = swfd(1:nv1,4) - (fds(1:nv1)-dn1(1:nv1) ) *0.450 


! DEBUG PRINT
!	do i=1,nv+1	
!	print'(A4,2i3,f8.1,3x,2(3f8.2,3x),4(2f8.2,3x))','CPRF',i,npc,pp(i),
!     &      fds(i),dn1(i),dn1(i)-fds(i),
!     &      fus(i),up1(i),up1(i)-fus(i)
!     ,      ,fdsdr(i),fdsdf(i)
!	enddo

!	print'(a30,4f8.4,i8)',
!     & 'CPRF: Toa ALBEDO (w/o, with)',fus(1)/fds(1),up1(1)/dn1(1)
!     &,u0,fus(nv+1)/fds(nv+1),npc

	fds(1:nv1)=dn1(1:nv1) 
	fus(1:nv1)=up1(1:nv1)

!	do i=1,nv1
!	if ( abs (fds(i)-fdsdr(i)-fdsdf(i) ) > 1.0) stop ' CHOU: DIR/DIFFUSE'
!	enddo

!End Chou absorption
	endif !lchou
!---------------------------------------------------------------------


      do 70 i = 1, nv
       xx = fds(i) -fus(i) - fds(i+1) + fus(i+1)
       dts(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
       xx = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
       dtir(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
       dt(i) = dts(i) + dtir(i)
70     continue

!TRMM (begin)
	iwncld(11:13) = 0
	where ( otau(11:13) > 1.0 ) iwncld(11:13) = 1

	call trmm_wnflt(ur)

	 trwn_flt_r = sat_flt_r
	 trwn_unf_r = sat_unf_r
	 trwn_f     = sat_f

!TRMM (end)

      return
      end

