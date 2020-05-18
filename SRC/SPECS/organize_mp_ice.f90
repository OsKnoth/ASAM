MODULE SpectralMicro_Mod

  USE Kind_Mod
  USE Control_Mod
  USE Parallel_Mod
  USE Domain_Mod
  USE Transport_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE Chemie_Mod

  USE data_mp_ice,  ONLY:                                              &
      tmeteors

  USE data_mp_ice,  ONLY: &
      anz,   &
      nmax,  &
      tmax,  &
      epsi,  &
      ipmax,  &            
      itmax,  &      
      smax,  &      
      simax,  &
      icond, &
      ikoll, & 
      ibrea, &
      iice,  &  ! Eis an/aus
      idepo, &  ! Depositionsgefrieren
      iimfr, &  ! Immersionsgefrieren
      ikofr, &  ! Eis - Tropfen Kollision
      iinsol,&  ! Kontaktgefrieren 
      ikeis, &  ! Eis - Eis Kollision
      imelt, &  ! Schmelzen
      ifrier,&  ! Gefrieren des Wassermantels
      isolutens, &  ! Oberflaechenspannung im Kelvinterm: 
      iideal,    &  ! Ideale Loesung (Tropfen)
      DELTAT,                     &
      mquer,                      &   
      rquer,                      &
      squer,                      &
      mfquer,                     &
      rfquer,                     &
      sfquer,                     &
      ifreeze,                    &
      miv

  USE utilities_mp_ice,       ONLY:    &
      saturation,                 &
      mp_allocate

  IMPLICIT NONE

  INTEGER,PRIVATE             ::                               &
     l,           & ! loop indices
     istep_lm, itime,   & ! LM and MP time step, input from LM
     k, it, si, ip,     & ! loop indices
     nsat,              & ! number of saturated gridpoints
     indx,              & !  fills in for i and j as a counter for filled gridpoints
     idummy
     
  INTEGER,PRIVATE    :: initialized=0

  REAL (RealKind),PRIVATE      ::       &
! MAKROSKOPISCH
     Talt,                       &
     TT,                         &
     usatt,                      &
     qq,                         &
     dTTdyn,                     &
     drho,                       &
     dpp,                        &
     dTT,                        &
     ddqq,                       &
     ddqqdyn

! MIKROPHYSIK
! dwasser
  REAL (RealKind),POINTER     ::       &
     dqq(:,:),              &
     dqqdyn(:,:),           & 
     dni(:,:,:,:),   &
     dqi(:,:,:,:),   &
     qq2(:,:),               &
     ddni(:,:),          &
     ddqi(:,:),          &
     nni(:,:),           &
     qqi(:,:),           & 
     dnw(:,:,:,:),     &
     dqw(:,:,:,:),     &
     dqws(:,:,:,:),    &
     dqwa(:,:,:,:),    &
     ddnw(:,:),            &
     ddqw(:,:),            &
     ddqws(:,:),           &
     ddqwa(:,:),           &
! deis
     dnf(:,:,:,:),     &
     dqf(:,:,:,:),     &
     dqfs(:,:,:,:),    &
     dqfa(:,:,:,:),    &
     dqfw(:,:,:,:),    &
     ddnf(:,:),            &
     ddqf(:,:),            &
     ddqfs(:,:),           &
     ddqfa(:,:),           &
     ddqfw(:,:),           &
! dinsoluble
! quergroessen
     mmquer(:,:),          &
     rrquer(:,:),          &
     ssquer(:,:),          &
     mmfquer(:,:),         &
     rrfquer(:,:),         &
     ssfquer(:,:),         &
! input vars
     nnw(:,:),             &
     qqw(:,:),             &
     qqws(:,:),            &
     qqwa(:,:),            &
     nnf(:,:),             &
     qqf(:,:),             &
     qqfs(:,:),            &
     qqfa(:,:),            &
     qqfw(:,:),            &
     nw(:,:,:,:),      &
     qw(:,:,:,:),      &
     qwa(:,:,:,:),     &
     nf(:,:,:,:),      &
     qf(:,:,:,:),      &     
     qfa(:,:,:,:)

  REAL(RealKind), PRIVATE, POINTER::  &
     vvtw(:,:),                 &
     vvtf(:,:)
     
  REAL(RealKind), PRIVATE, POINTER :: Th(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ThRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoVRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: RhoLRhs(:,:,:,:)
  
CONTAINS

SUBROUTINE InitSpectralMicro(p,Vector,ix,iy)
   TYPE(Vec4_T) :: Vector(:)
   REAL(RealKind), POINTER :: p(:,:,:,:)
   REAL(RealKind) :: qv(nz),rhoz(nz)
   INTEGER :: ix,iy,iz
   ! set microphysical timestep
   ! attention! data_microphysics has still to be merged with 
   ! data_mp_ice. At the moment, deltat is set in both!
   ! And initialized in lmorg as well due to this..
   ! set number of bins
   ipmax=1
   itmax=1
   smax=1
   simax=1

   Th=>Vector(ThPos)%c
   Rho=>Vector(RhoPos)%c
   RhoV=>Vector(RhoVPos)%c
   RhoL=>Vector(RhoCPos)%c

   IF (MyId == 0) THEN
      PRINT *,'ALLOCATING MICROPHYSICS QUANTITIES'
   ENDIF
   CALL mp_allocate (nz,nFrac)
      
   IF (MyId == 0) THEN
      PRINT *,'INITIALIZING MICROPHYSICS'
   ENDIF
   
   ALLOCATE(vvtw(1:nz,nFrac))
   ALLOCATE(vvtf(1:nz,nFrac))

   ALLOCATE(dqq(1,ipmax))
   ALLOCATE(dqqdyn(1,ipmax))
   ALLOCATE(dni(1,simax,itmax,ipmax))
   ALLOCATE(dqi(1,simax,itmax,ipmax))
   ALLOCATE(qq2(1,ipmax))
   ALLOCATE(ddni(simax,itmax))
   ALLOCATE(ddqi(simax,itmax))
   ALLOCATE(nni(simax,itmax))
   ALLOCATE(qqi(simax,itmax))

   ALLOCATE(dnw(1,nFrac,smax,ipmax))
   ALLOCATE(dqw(1,nFrac,smax,ipmax))
   ALLOCATE(dqws(1,nFrac,smax,ipmax))
   ALLOCATE(dqwa(1,nFrac,smax,ipmax))
   ALLOCATE(ddnw(nFrac,smax))
   ALLOCATE(ddqw(nFrac,smax))
   ALLOCATE(ddqws(nFrac,smax))
   ALLOCATE(ddqwa(nFrac,smax))
! deis
   ALLOCATE(dnf(1,nFrac,smax,ipmax))
   ALLOCATE(dqf(1,nFrac,smax,ipmax))
   ALLOCATE(dqfs(1,nFrac,smax,ipmax))
   ALLOCATE(dqfa(1,nFrac,smax,ipmax))
   ALLOCATE(dqfw(1,nFrac,smax,ipmax))
   ALLOCATE(ddnf(nFrac,smax))
   ALLOCATE(ddqf(nFrac,smax))
   ALLOCATE(ddqfs(nFrac,smax))
   ALLOCATE(ddqfa(nFrac,smax))
   ALLOCATE(ddqfw(nFrac,smax))
! dinsoluble
! quergroessen
   ALLOCATE(mmquer(nFrac,smax))
   ALLOCATE(rrquer(nFrac,smax))
   ALLOCATE(ssquer(nFrac,smax))
   ALLOCATE(mmfquer(nFrac,smax))
   ALLOCATE(rrfquer(nFrac,smax))
   ALLOCATE(ssfquer(nFrac,smax))
! input vars
   ALLOCATE(nnw(nFrac,smax))
   ALLOCATE(qqw(nFrac,smax))
   ALLOCATE(qqws(nFrac,smax))
   ALLOCATE(qqwa(nFrac,smax))
   ALLOCATE(nnf(nFrac,smax))
   ALLOCATE(qqf(nFrac,smax))
   ALLOCATE(qqfs(nFrac,smax))
   ALLOCATE(qqfa(nFrac,smax))
   ALLOCATE(qqfw(nFrac,smax))
   ALLOCATE(nw(1,nFrac,smax,ipmax))
   ALLOCATE(qw(1,nFrac,smax,ipmax))
   ALLOCATE(qwa(1,nFrac,smax,ipmax))
   ALLOCATE(nf(1,nFrac,smax,ipmax))
   ALLOCATE(qf(1,nFrac,smax,ipmax))
   ALLOCATE(qfa(1,nFrac,smax,ipmax))

   DO iz=iz0+1,iz1
    qv(iz)=RhoV(ix,iy,iz,1)/(Rho(ix,iy,iz,1))
    rhoz(iz)=Rho(ix,iy,iz,1)
   END DO
   
   CALL init_mp(                  &
              Vector(Position('a_LP_N'))%c(ix,iy,:,1:nFrac), & ! nw, & ! LM   
              Vector(Position('a_LP_Rho'))%c(ix,iy,:,1:nFrac), & ! qw, & ! LM   
              Vector(Position('a_LP_RhoAS'))%c(ix,iy,:,1:nFrac), & ! qws, & ! LM   
              Vector(Position('a_LP_RhoAI'))%c(ix,iy,:,1:nFrac), & ! qwa, & ! LM   
              vvtw(1:nz,:)                                   , & ! vtw,& ! LM   
              Vector(Position('a_MPP_N'))%c(ix,iy,:,1:nFrac), & ! nf, & ! LM   
              Vector(Position('a_MPP_Rho'))%c(ix,iy,:,1:nFrac), & ! qf, & ! LM   
              Vector(Position('a_MPP_RhoAS'))%c(ix,iy,:,1:nFrac), & ! qfs, & ! LM   
              Vector(Position('a_MPP_RhoAI'))%c(ix,iy,:,1:nFrac), & ! qfa, & ! LM   
              Vector(Position('a_MPP_RhoAF'))%c(ix,iy,:,1:nFrac), & ! qfw, & ! LM   
              vvtf(1:nz,:)                                   , & ! vtf,& ! LM   
              rhoz(:)                                        , & ! rho,& ! LM   
              Th(ix,iy,:,1)                                  , & ! TT, & ! LM   
              p(ix,iz,:,1)                                   , & ! pp, & ! LM   
              Vector(Position('a_IP_N'))%c(ix,iy,:,1:nFrac), & ! ni, & ! LM   
              Vector(Position('a_IP_Rho'))%c(ix,iy,:,1:nFrac), & ! qia, & ! LM   
              nz                                             , & ! kz, & ! LM 
              zP(:)                                          , & ! dz, & ! LM
              nFrac                                          , & ! jmax = nFrac
              qv(:)                                            & ! for satt, ! LM
              )
   
   ! Set top and bottom velocities
!   tmeteors(1)%v (1,:) = 0.0d0
!   tmeteors(2)%v (1,:) = 0.0d0
!   tmeteors(1)%v (nz+1,:) = tmeteors(1)%v (nz,:)
!   tmeteors(2)%v (nz+1,:) = tmeteors(2)%v (nz,:)
   
   ! set velocities vor insolubles equal to velocity of mixed phase
!   tmeteors(3)%v (:,:) = tmeteors(2)%v (:,:)


   ! set density at ground level (used in koll)
!   rhotot0 = 1.014d0
   
   IF (MyId == 0) THEN 
      OPEN (12,FILE="mp_par.dat")
      WRITE (12,*) 
      WRITE (12,*) 'itmax     = ',itmax
      WRITE (12,*) 'SMAX      = ',SMAX
      WRITE (12,*) 'SIMAX     = ',SIMAX
      WRITE (12,*) 'anz       = ',anz
      WRITE (12,*) 'ipmax     = ',ipmax
      WRITE (12,*) 'nmax      = ',nmax
      WRITE (12,*) 'tmax      = ',tmax
      WRITE (12,*) 'epsi      = ',epsi
      WRITE (12,*) 'icond     = ',icond
      WRITE (12,*) 'ikoll     = ',ikoll
      WRITE (12,*) 'ibrea     = ',ibrea
      WRITE (12,*) 'iice      = ',iice
      WRITE (12,*) 'idepo     = ',idepo
      WRITE (12,*) 'iimfr     = ',iimfr
      WRITE (12,*) 'ikofr     = ',ikofr
      WRITE (12,*) 'iinsol    = ',iinsol
      WRITE (12,*) 'ikeis     = ',ikeis
      WRITE (12,*) 'imelt     = ',imelt
      WRITE (12,*) 'ifrier    = ',ifrier
      WRITE (12,*) 'isolutens = ',isolutens
      WRITE (12,*) 'iideal    = ',iideal
      CLOSE (12)
   ENDIF
END SUBROUTINE InitSpectralMicro

SUBROUTINE CalcSpectralMicro(p,Vector,Rhs,ix,iy,iz)
  TYPE(Vec4_T) :: Vector(:),Rhs(:)
  REAL(RealKind), POINTER :: p(:,:,:,:)
  INTEGER :: ix,iy,iz
   !-----------------------------------------
   ! Test for non empty k planes (water only)
   !-----------------------------------------
  ! IF (SUM(tmeteors(1)%qc (:,:,height, 1,:)) .NE. 0.0) THEN 
  !    nsat = 0
  !    DO j = jstartpar , jendpar
  !       DO i = istartpar , iendpar           
  !          qc_dummy = SUM(tmeteors(1)%qc(i,j,height,1,:))
  !          IF (qc_dummy .NE. 0.0) THEN
  !             nsat = nsat+1
  !             iwrk (nsat) = i
  !             jwrk (nsat) = j
  !          ENDIF
  !       ENDDO
  !    ENDDO
  ! ENDIF
   
   !-----------------------------------------
   ! loop over gridpoint starts
   !-----------------------------------------

     
      
      ! Calculte averaged values for mass and radius
      
      DO si=1,smax
         DO ip=1,ipmax
 !           nw(1,:,si,ip)  = tmeteors(1)%qc(i,j,height,1,1          :  nFrac)
 !           qw(1,:,si,ip)  = tmeteors(1)%qc(i,j,height,1,1+  nFrac:2*nFrac)
 !           qwa(1,:,si,ip) = tmeteors(1)%qc(i,j,height,1,3*nFrac+1:4*nFrac)
 !           nf(1,:,si,ip)  = tmeteors(2)%qc(i,j,height,1,1          :  nFrac)
 !           qf(1,:,si,ip)  = tmeteors(2)%qc(i,j,height,1,1+  nFrac:2*nFrac)
 !           qfa(1,:,si,ip) = tmeteors(2)%qc(i,j,height,1,3*nFrac+1:4*nFrac)
         ENDDO
      ENDDO
      
      CALL xquerxd(mquer,rquer,squer,mfquer,rfquer,sfquer,nw,qw,qwa,nf,qf,qfa,1)         
      
      ! Just to be sure that dimensions are right
      IF (ipmax .NE. 1) THEN
         print *,'FATAL ERROR IN ORGANIZE_MP_ICE: IPMAX .NE. 1'
         print *,'MODEL WILL BE ABORTED  '
         STOP
      ENDIF
      IF (itmax .NE. 1) THEN
         print *,'FATAL ERROR IN ORGANIZE_MP_ICE: ITMAX .NE. 1'
         print *,'MODEL WILL BE ABORTED  '
         STOP
      ENDIF
      IF (smax .NE. 1) THEN
         print *,'FATAL ERROR IN ORGANIZE_MP_ICE: SMAX .NE. 1'
         print *,'MODEL WILL BE ABORTED  '
         STOP
      ENDIF
      IF (ipmax .NE. 1) THEN
         print *,'FATAL ERROR IN ORGANIZE_MP_ICE: ip .NE. 1'
         print *,'MODEL WILL BE ABORTED  '
         STOP
      ENDIF
      
      ! calculate density and pressure for microphysics due to dynamics
!      rho_mp(i,j)        = rho_old(i,j,height)                 +         &
!                                  ((itime*deltat) * &
!                                             (rho(i,j,height     )-rho_old(i,j,height     )) / (itimemax))
!      p_mp               = p0(i,j,height)+pp(i,j,height,nnow)  +         &
!                                  ((itime*deltat) * &
!                                             (pp (i,j,height,nnew)-pp     (i,j,height,nnow)) / (itimemax))
      
      ! Define input values for microphysics
      
      ! temperature and water vapor mixing ratio
!      qq       = qv(i,j,height,1)   ! QV      
!      TT       = t (i,j,height,nnew  )   ! TABS   
      
      ! calculate saturation, t already updated dynamically
!      CALL saturation (usatt,TT, p_mp, qq)  ! SATT scalar
      
      ! macroscopic quantities
      ddqq    = 0.0d0
      ddqqdyn = 0.0d0
      dTT     = 0.0d0
      dTTdyn  = 0.0d0
      drho    = 0.0d0
      dpp     = 0.0d0
      
      DO si=1,smax
         ! water
         nnw(:,si)  = Vector(Position('a_LP_N'))%c(ix,iy,iz,1:nFrac)
         qqw(:,si)  = Vector(Position('a_LP_Rho'))%c(ix,iy,iz,1:nFrac)
         qqws(:,si) = Vector(Position('a_LP_RhoAS'))%c(ix,iy,iz,1:nFrac)
         qqwa(:,si) = Vector(Position('a_LP_RhoAI'))%c(ix,iy,iz,1:nFrac)
         
         ! ice
         nnf(:,si)  = Vector(Position('a_MPP_N'))%c(ix,iy,iz,1:nFrac)
         qqf(:,si)  = Vector(Position('a_MPP_Rho'))%c(ix,iy,iz,1:nFrac)
         qqfs(:,si) = Vector(Position('a_MPP_RhoAS'))%c(ix,iy,iz,1:nFrac)
         qqfa(:,si) = Vector(Position('a_MPP_RhoAI'))%c(ix,iy,iz,1:nFrac)
         qqfw(:,si) = Vector(Position('a_MPP_RhoAF'))%c(ix,iy,iz,1:nFrac)
         
         ! dwater
         ddnw(:,si)    = 0.0d0
         ddqw(:,si)    = 0.0d0
         ddqws(:,si)   = 0.0d0
         ddqwa(:,si)   = 0.0d0
         
         ! dice
         ddnf(:,si)    = 0.0d0
         ddqf(:,si)    = 0.0d0
         ddqfs(:,si)   = 0.0d0
         ddqfa(:,si)   = 0.0d0
         ddqfw(:,si)   = 0.0d0
         
         ! avaraged
         mmquer(:,si)  = mquer(1,:,si)
         rrquer(:,si)  = rquer(1,:,si)
         ssquer(:,si)  = squer(1,:,si)
         mmfquer(:,si) = mfquer(1,:,si)
         rrfquer(:,si) = rfquer(1,:,si)
         ssfquer(:,si) = sfquer(1,:,si)
         
         DO it=1,itmax
            ! insoluble aerosols
            nni(:,it)  = Vector(Position('a_IP_N'))%c(ix,iy,iz,1:nFrac)
            qqi(:,it)  = Vector(Position('a_IP_Rho'))%c(ix,iy,iz,1:nFrac)
            
            ! dinsoluble aerosols
            ddni(:,it) = 0.0d0
            ddqi(:,it) = 0.0d0
         END DO
      END DO
      
      DELTAT=1.0d0
      !-------------------------------------   
      ! Call the microphysics
      CALL cloudxd(    &
           ! makroskopic
           dTT,     & ! DELT
           dTTdyn,  & ! deltdyn
           dpp,     & ! DP
           drho,    & ! DRHO
           ! dwater
           ddnw,    & ! DNW
           ddqw,    & ! DQW
           ddqws,   & ! DQS
           ddqwa,   & ! DQA
           ! makroskopic
           ddqq,    & ! DQV
           ddqqdyn, & !dqvdyn
           usatt,   & ! SATT
           ! water
           nnw,     & ! NW
           qqw,     & ! QW
           qqws,    & ! QS
           qqwa,    & ! QA
           ! makroskopisch
           qq,& ! 
           p(ix,iy,iz,1),  &  ! PTOT
           rho(ix,iy,iz,1),&  !rho(height)
           TT,&               ! TABS
           ! averaged
           mmquer,  & ! MQUER
           rrquer,  & ! RQUER
           ssquer,  & ! SQUER
           mmfquer, & ! MIQUER
           rrfquer, & ! RIQUER
           ssfquer, & ! SIQUER
           ! dice
           ddnf,    & ! DNFROD
           ddqf,    & ! DQFROD
           ddqfs,   & ! DQSFROD
           ddqfa,   & ! DQAFROD
           ddqfw,   & ! dqwfrod
           ! ice
           nnf,     & ! NFROD 
           qqf,     & ! QFROD 
           qqfs,    & ! QSFROD
           qqfa,    & ! QAFROD
           qqfw,    & ! qwfrod
           ! dinsoluble
           ddni,    & ! DNWINS
           ddqi,    & ! DQAINS
           ! insoluble
           nni,     & ! NWINS
           qqi,     & ! QAINS
           ! makroskopic
           TT,      & ! TABSold
           ! microphysics
           ifreeze, & ! ifreeze
           miv,     & ! miv
           1,       & ! ip = 1
           ! index                            
           iz,      & ! kk
           ! microphysics
           vvtw,    & ! vvtw
           vvtf,    &
           ix,iy,iz,0.0d0 &
           )      ! vvtf
      
!------------------------------------------------------------------------------
! End of module procedure src_microphysics
!------------------------------------------------------------------------------
END SUBROUTINE CalcSpectralMicro

SUBROUTINE SpectralMicro(p,Vector,Rhs)
  TYPE(Vec4_T) :: Vector(:), Rhs(:)
  REAL(RealKind), POINTER :: p(:,:,:,:)
  INTEGER :: ix,iy,iz

  IF(initialized.eq.0)THEN
   DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
     CALL InitSpectralMicro(p,Vector,ix,iy)   
    END DO
   END DO
   initialized=1
  END IF
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        CALL CalcSpectralMicro(p,Vector,Rhs,ix,iy,iz) 
      END DO
    END DO
  END DO
END SUBROUTINE SpectralMicro

END MODULE SpectralMicro_Mod
