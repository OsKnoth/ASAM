!================================================================
	module pckd24a
	parameter (ncoef=7,nreg=2,nband=12)

	real  h2obnd(nband)
	real ck24_3(ncoef,nreg,nband)! lin Plank wgt

	data h2obnd /-5,-3.5,-2.0,-2,-1.,-4,-4,-4,-3,-3.5,-3,-2/

       data ck24_3/  !ckd24fu.fuliou.lin.plnk.out
! band        1
     x 1.667e+00, 9.421e-01,-7.358e-03, 1.355e+00,
     x 2.557e+03, 5.798e+01,-4.570e-01,
! band        1
     x 6.417e+00, 1.002e+00,-6.991e-03, 1.010e+00,
     x 1.203e+01, 4.501e-02,-2.428e-02,
! band        2
     x 2.390e+00, 9.528e-01,-6.058e-03, 1.071e+00,
     x 2.676e+02, 9.848e+00,-1.459e-01,
! band        2
     x 4.849e+00, 1.002e+00,-6.910e-03, 8.961e-01,
     x 1.635e+01, 2.115e-02, 7.243e-02,
! band        3
     x 2.326e+00, 9.720e-01,-6.551e-03, 8.739e-01,
     x 6.984e+01, 8.346e-01, 4.824e-02,
! band        3
     x 5.002e+00, 1.005e+00,-9.286e-03, 6.222e-01,
     x 1.168e+01, 3.611e-03, 3.148e-01,
! band        4
     x-4.865e+00, 8.455e-01,-6.911e-03, 1.475e+00,
     x 2.905e+02, 7.078e+00,-6.846e-01,
! band        4
     x 4.596e+00, 1.012e+00,-1.152e-02, 5.713e-01,
     x 1.270e+01,-1.395e-03, 3.447e-01,
! band        5
     x-5.396e+00, 8.596e-01,-8.479e-03, 1.619e+00,
     x 1.664e+02, 3.236e+00,-7.782e-01,
! band        5
     x 7.478e+00, 1.007e+00,-1.963e-02, 2.771e-01,
     x 6.021e+00,-4.489e-03, 6.709e-01,
! band        6
     x 1.262e+00, 2.347e-01,-2.360e-02, 1.655e-01,
     x 5.068e+02, 2.462e+01, 3.920e-01,
! band        6
     x 9.334e+00, 1.002e+00,-2.429e-02, 3.575e-02,
     x 2.751e-01,-1.189e-03, 9.593e-01,
! band        7
     x-1.222e+00, 5.423e-01,-2.327e-02, 5.197e-01,
     x 6.423e+02, 5.038e+01, 1.502e-01,
! band        7
     x 8.506e+00, 1.000e+00,-2.339e-02, 8.891e-03,
     x-6.805e-01,-1.639e-04, 9.917e-01,
! band        8
     x-3.638e+00, 8.534e-01,-1.344e-02, 6.816e-01,
     x 5.385e+02, 4.428e+01,-6.366e-03,
! band        8
     x 6.921e+00, 1.002e+00,-1.974e-02, 6.350e-02,
     x 6.838e-01,-1.121e-03, 9.237e-01,
! band        9
     x-2.329e+00, 7.893e-01,-2.588e-03, 1.017e+00,
     x 1.525e+02, 1.029e+01,-1.486e-01,
! band        9
     x 6.742e-01, 1.008e+00,-3.376e-03, 9.105e-01,
     x 1.074e+01,-3.307e-03, 5.741e-02,
! band       10
     x-1.677e+00, 9.173e-01,-5.780e-03, 1.504e+00,
     x 7.886e+02, 2.288e+01,-5.999e-01,
! band       10
     x 3.396e+00, 1.005e+00,-3.433e-03, 1.012e+00,
     x 7.635e+00, 3.010e-03,-2.418e-02,
! band       11
     x 7.943e-01, 9.260e-01,-5.050e-03, 1.141e+00,
     x 2.221e+02, 1.021e+01,-2.246e-01,
! band       11
     x 3.356e+00, 1.002e+00,-4.719e-03, 9.578e-01,
     x 6.164e+00, 1.186e-03, 2.264e-02,
! band       12
     x-5.874e+00, 7.060e-01,-1.532e-03, 1.141e+00,
     x 1.463e+02, 6.534e+00,-4.308e-01,
! band       12
     x 4.709e-01, 1.010e+00,-6.067e-03, 8.513e-01,
     x 1.161e+01,-6.629e-03, 8.885e-02
     &/

	contains
c=====================================================================
        function parm_ckd24(iband,amnt,patm,temp,dz)
c Parameterization of CKD_2.4 continuum over Fu-Liou Bands
c Input:
c iband  =  integer (1-12) where
c	  Band 1 ='  5:280cm-1'
c         Band 2 ='280:400cm-1'
c         Band 3 ='400:540cm-1'
c         Band 4 ='540:670cm-1'
c         Band 5 ='670:800cm-1'
c         Band 6 ='800:980cm-1'
c         Band 7 ='980:1100cm-1'
c         Band 8 ='1100:1250cm-1'
c         Band 9 ='1250:1400cm-1'
c         Band10 ='1400:1700cm-1'
c         Band11 ='1700:1900cm-1'
c         Band12 ='1900:2200cm-1'
c amnt = h2O ammount (g/cm**2)
c patm = pressure (atm)
c temp = temperature (k)
c dz   = pathlength (Km)
c Output:
c parm_ckd24 = parameterized CKD_2.4optical depth for band
c234567890123456789012345678901234567890123456789012345678901234567890

	

! These Regressions are more sensitive to pathlength
! So accomodations for very Thin or Thick layers are made.
	 dz1  =dz
	 factor=1.000
	   if ( dz < 0.25  ) then
		factor = 0.25/dz
		 dz1   = 0.25
	   elseif (dz > 1.50) then
	        factor = 1.50/dz
	        dz1    = 1.50
	   endif

	amnt1=amnt*factor

! Regression is now broken up into TWO parts one for small
! one for large	water vapor ammounts.

	ireg=1
	if ( log(amnt1) > h2obnd(iband) ) ireg=2

       ph2o = amnt1 *(8.314d+07 *temp )/
     &              (dz1*1.0d+05*18.01534 *1.01325d+06)

	patmx = log(patm)

!	print'(8f8.3)',log(amnt1),temp,patmx,(ph2o),amnt1,log(ph2o),dz1
!	print*, ireg,iband,aa(1:7,ireg,iband)

 	tau_log	    = ck24_3(1,ireg,iband)	      +
     $		      ck24_3(2,ireg,iband)* log(amnt1)  +
     $		      ck24_3(3,ireg,iband)* temp       +
     $                ck24_3(4,ireg,iband)* patmx       +
     $		      ck24_3(5,ireg,iband)* (ph2o)     +
     $		      ck24_3(6,ireg,iband)* amnt1       +
     $		      ck24_3(7,ireg,iband)* log(ph2o) 
!	print*,tau_log
	parm_ckd24 = exp ( tau_log )
	parm_ckd24 = parm_ckd24/factor
	return
	end function parm_ckd24
!-------------------------------------------------------------------------------
	end module pckd24a
