MODULE RadiationFuLiou_Mod
  USE DataType_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Physics_Mod

  INTEGER, PRIVATE :: i,j,k
  INTEGER, PRIVATE :: ic,icE

  REAL(RealKind), PRIVATE, POINTER :: c(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qc(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)

 CONTAINS

  SUBROUTINE RadiationFuLiouCompute(NumBound)
   USE RadParams
   USE TRMM_WINDOW
!   USE Physics_Mod, ONLY: raddiffus,raddirekt,cosSun
   USE Physics_Mod, ONLY: cosSun

   INTEGER,INTENT(IN),OPTIONAL :: NumBound

   common /atmos/    pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
   common /clouds/   pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
   common /rains/    prwc(nvx)
   common /graups/   pgwc(nvx)
   common /umcon/    umco2, umch4, umn2o
   common /radiat/   fds(nv1x), fus(nv1x), dts(nvx),&
                     fdir(nv1x), fuir(nv1x), dtir(nvx),&
                     fd(nv1x), fu(nv1x), dt(nvx),&
                     fdsdr(nv1x), fdsdf(nv1x)
   common /dfsout/   fu1(nv1x), fd1(nv1x)
   common /planci/   bf(nv1x), bs
   common /wndow/    fuwn(nv1x),fdwn(nv1x)!! WINDOW FLUX
   common /dkfrwn/   idkfr
   common /cfcs/     cfc_conc(3)
   common /cont_tas/ iwtas
   common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x)
   common /TRMMWNOUT/ trwn_flt_r, trwn_unf_r, trwn_f
   common /select_solar_spectra/ isolar_spectrum
   common /uvflux/   uvfu(nv1x,10),uvfd(nv1x,10), &
                     uvdir(nv1x,10),uvdif(nv1x,10)
   common /swflux/   swfu(nv1x,6),swfd(nv1x,6),    &
                    swdir(nv1x,6),swdif(nv1x,6)
   common /lwflux/   rlwfu(nv1x,7:20),rlwfd(nv1x,7:20),&
                    sbf(nv1x,7:20),sbs(7:20)
   common /aerpout/  aprop(nvx,mxac,3)
   common /seijik/   isksw

   real*4  as(mbsx),as1(12), ee(mbirx)
   real rh(nv1x)
   real sh_aer(3)

101	format(f7.5)
102	format(I2)

   logical,save :: first=.TRUE.
   real z,g,Rho,t0,Gamma,cld

   if(first) then
!   open(1,file='../testatms/barbados.lay',status='old')
    open(1,file='../testatms/kmls.lay',status='old')
    read (1,*) skint,nlay
    nlev=nlay+1
    do i=1,nlev
     read(1,*) pp(i),pt(i),ph(i),po(i),cld
    enddo
    close(1)
   endif

   g=9.81
   Rho=1.23
   t0=skint
   z=z0

   nlev=nlay+1
   nv=nlay-1
   pts=290
   isolar_spectrum=4         !
   nhb=2                     !
   isksw=0
   icld=72
   tau_vis=0
   tau_aer=0
   itp=3
   icont=5
   no_rayle=0
   umco2=350
   idkfr=2
   iwtas=3
   ur=0.9
   iaform=3
   ifg=0
   nv1=nv+1
   ndfs = nv !! set up rad24a
   mdfs = nv1
   ndfs4 = 4*ndfs
   ndfs2=  2*ndfs
   mb   = 18
   mbs  = 6
   mbir = 12
   nc   = 8
   ss=1368
   umco2 = 330.0
   umch4 = 1.75
   umn2o = 0.31
   cfc_conc= (/0.268e-09 , 0.503e-09 ,0.105e-09/)
   as(1:mbsx)=0.0
   as1(1:12) =0.0
!  as(1:mbsx)=0.15 ! computed in soil_mod
!  as1(1:10) =0.15 ! computed in soil_mod
   u0=cosSun
   ee(1:mbirx)=1.0
   fourssl = .false. ! two stream
   foursir = .false.
   ipha=1
   clwp=0.0
   fraca=0.10
   iclrcld=1
   irobckd=icont
   a_wlis =-9999.
   a_taus =-9999.
   itps = -9999
   nac=0
   itps(1)=itp ! PRIMARY
   itps(2)=11  ! soot
   plwc = 0
   pre  = 0

   call rad (as ,as1, u0, ss, pts, ee , ur ,ITriErr)

         !do i=1,nv
	 !call ql_rh(rh(i),pt(i),pp(i),ph(i))
 	 !print '(i3,f8.2,f8.1,f8.1,5f8.2,1x,3(2x,3f7.2))',i,pp(i),pt(i),rh(i) &
         !     ,fds(i),fus(i),fdsdr(i),fdsdf(i),dts(i),dt(i),fdir(i),fuir(i),fdwn(i),fuwn(i)
 	 !write(*,'(i3,15f9.3)') i,pp(i),rh(i),fds(i),fus(i),fdir(i),fuir(i),fdsdr(i) &
         !     ,fdsdf(i),fdwn(i),fuwn(i),fd(i),fu(i),dt(i),dts(i),dtir(i)
	 !xxs(i) = fds(i) -fus(i) - fds(i+1) + fus(i+1)
         !dts(i) = 8.4392 * xxs(i) / ( pp(i+1) - pp(i) )
         !xxir(i) = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
         !dtir(i) = 8.4392 * xxir(i) / ( pp(i+1) - pp(i) )
         !enddo

!    raddiffus=fdsdf(nv)+fdir(nv) !// diffuse sw + infrared from atmosphere
!    raddirekt=fdsdr(nv)          !// direkte sw
    BoundCell(NumBound)%raddirekt=fdsdr(nv)          !// direkte sw
    BoundCell(NumBound)%raddiffus=fdsdf(nv)          !// diffuse sw
    BoundCell(NumBound)%radinfred=fdir(nv)           !// infrared from atmosphere

!   WRITE (*,*) 'FuLiou: raddirekt, raddiffus, radinfred'
!   WRITE (*,*) BoundCell(NumBound)%raddirekt, BoundCell(NumBound)%raddiffus,  BoundCell(NumBound)%radinfred

  END SUBROUTINE

END MODULE RadiationFuLiou_Mod
