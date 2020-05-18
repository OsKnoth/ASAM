        SUBROUTINE driver
	USE RadParams
	use TRMM_WINDOW

	common /atmos/    pp(nv1x), pt(nv1x), ph(nv1x), po(nv1x)
	common /clouds/   pre(nvx), plwc(nvx), pde(nvx), piwc(nvx)
	common /rains/    prwc(nvx)
	common /graups/   pgwc(nvx)
        common /umcon/    umco2, umch4, umn2o 
	common /radiat/   fds(nv1x), fus(nv1x), dts(nvx),
     &                    fdir(nv1x), fuir(nv1x), dtir(nvx),
     &	                  fd(nv1x), fu(nv1x), dt(nvx),
     &                    fdsdr(nv1x), fdsdf(nv1x)  
	common /dfsout/   fu1(nv1x), fd1(nv1x)
	common /planci/   bf(nv1x), bs
        common /wndow/    fuwn(nv1x),fdwn(nv1x)!! WINDOW FLUX
	common /dkfrwn/   idkfr
        common /cfcs/     cfc_conc(3)
	common /cont_tas/ iwtas
  	common /radiance/ fiurt(nv1x),fiurw(nv1x),fiur(nv1x) 
	common /TRMMWNOUT/ trwn_flt_r, trwn_unf_r, trwn_f 
        common /select_solar_spectra/ isolar_spectrum
        common /uvflux/   uvfu(nv1x,10),uvfd(nv1x,10),   
     &                    uvdir(nv1x,10),uvdif(nv1x,10)
        common /swflux/   swfu(nv1x,6),swfd(nv1x,6),     
     &                    swdir(nv1x,6),swdif(nv1x,6)
        common /lwflux/   rlwfu(nv1x,7:20),rlwfd(nv1x,7:20),
     &                    sbf(nv1x,7:20),sbs(7:20)
	common /aerpout/  aprop(nvx,mxac,3)
	common /seijik/   isksw
	
        real*4  as(mbsx),as1(10), ee(mbirx)	
	real rh(nv1x)
        real sh_aer(3)
        
101	format(f7.5)
102	format(I2)

        !beispielmodell mit 20 schichten in z-richtung und Profile für p und t im Modellbereich
        integer nz
        real ,dimension(100) :: dz,iz0,iz1,pmod,tmod,dtmod,hmod
        real z,p0,g,Rho,t0,z0,Gamma,cld
                
        nz=100
        z0=0
        do i=1,nz
         dz(i)=20
        enddo
        
        p0=1013.23
        g=9.81
        Rho=1.23
        t0=294
        z=z0
        Gamma=0.00582  ! passt grad an kmls.lay
         
        do i=1,nz+1      
         pmod(i)=p0*(1-2./7.*9.81/287./Gamma*log(1+Gamma*z/t0))**(7./2.)
         tmod(i)=(t0+Gamma*z)*(1-2./7.*9.81/287./Gamma*log(1+Gamma*z/t0))
         hmod(i)=1.23597E-02-0.000004*z
         z=z+dz(i);
        enddo
        
        write(*,*) "reading atmosphere Mid-latitude summer..."
        !restliche Athmosphäre einlesen
        open(1,file='./testatms/kmls.lay',status='old')
        read (1,*) skint,nlay
        nlev=nlay+1
        do i=1,nlev
         read(1,*) pp(i),pt(i),ph(i),po(i),cld
         !write(*,*) pp(i),pt(i),ph(i),po(i),cld
        enddo
        close(1)
    
        i=0
        do while (pp(i+1)<pmod(nz+1))    ! Index ab wann modellatmosphäre losgeht
          i=i+1
        enddo
                
        do j=1,nz+1                    ! Modellatmosphäre anfügen
         pp(i+j)=pmod(nz+1-j)
         pt(i+j)=tmod(nz+1-j)
         ph(i+j)=hmod(nz+1-j)
         po(i+j)=po(i)
        enddo
        
        nlay=i+nz
        nlev=i+nz+1
                
        !call utrhpwx(nlay,pp,pt,ph,po,pts)
        !write(*,*) "resulting atmosphere------------------"
        nv=nlay-1
        do i=1,nlay
        ! write(*,*) i,pp(i),pt(i),ph(i),po(i)
        enddo
        
        write(*,*) "calculating fluxes...."
        
        pts=294
	isolar_spectrum=4         !
	nhb=2                     !
	isksw=0
        u0=0.866025
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
        !umco2 = 330.0
	umch4 = 1.75
	umn2o = 0.31
        cfc_conc= (/0.268e-09 , 0.503e-09 ,0.105e-09/)
	as(1:mbsx)=0.15
	as1(1:10) =0.15
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
	     
               
        do i=1,nv+1
	 call ql_rh(rh(i),pt(i),pp(i),ph(i))
	 !print '(i3,f8.2,f8.1,f8.1,5f8.2,1x,3(2x,3f7.2))',i,pp(i),pt(i),rh(i)
     &   !       ,fds(i),fus(i),fdsdr(i),fdsdf(i),dts(i),fdir(i),fuir(i),fdwn(i),fuwn(i)
         !write(*,'(i3,15f9.3)') i,pp(i),rh(i),fds(i),fus(i),fdir(i),fuir(i),fdsdr(i)
     &   !       ,fdsdf(i),fdwn(i),fuwn(i),fd(i),fu(i),dt(i),dts(i),dtir(i)
         !write(1,'(f8.2,f8.1,f8.1,5f8.2,1x,3(2x,3f7.2))') pp(i),pt(i),rh(i)
     &   !       ,fds(i),fus(i),fdsdr(i),fdsdf(i),dts(i),fdir(i),fuir(i),fdwn(i),fuwn(i)
        enddo
        
        write(*,*) "iz   heatingrate in K/day"
        do i=1,nz
         dtmod(i)=dt(nv+1-i)      ! Heatingrate in K/day für die z Schichten des Modells
         write(*,'(i3,f8.2)') i,dtmod(i)
        enddo
        
        !print'(a9,14f6.2,2x,2f8.1)','LWSPECTOA',rlwfu( 1,7:20),sum(rlwfu(1,7:20)),fiurt(1)
	!print'(a9,14f6.2,2x,2f8.1)','LWSPECSFC',rlwfu(nv1,7:20),sum(rlwfu(nv1,7:20)),fiurt(nv1)
        
        write(*,*) "ready."
        stop	
        end
        
        
