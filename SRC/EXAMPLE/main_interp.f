MODULE parameters
!==============================================================================
IMPLICIT NONE

  INTEGER, PARAMETER       ::                                                 &
       ireals    = SELECTED_REAL_KIND (12,200),                               &
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

       iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers

  INTEGER (KIND=iintegers), POINTER, SAVE ::                                  &
   i_index(:,:),   & ! i-index of coarse mesh grid point which is lower left
                     ! to a given ASAM (fine mesh) grid point
   j_index(:,:)      ! j-index of coarse mesh grid point which is lower left
                     ! to a given ASAM (fine mesh) grid point
                                                                                                         
  INTEGER (KIND=iintegers)    ::                                              &
   ie2lm,           & ! ie for ASAM, local domain
   je2lm,           & ! je for ASAM, local domain
   ie_in,           & ! ie for LM, total domain
   je_in              ! je for LM, total domain

  REAL (KIND=ireals), POINTER, SAVE       ::                                  &
   x_wght (:,:),   & ! relative distance between x- (i-) coarse mesh and
                     ! fine mesh grid points
   y_wght (:,:)      ! relative distance between y- (j-) coarse mesh and
                     ! fine mesh grid points
  
  REAL (KIND=ireals)          ::                                              &
   pollat,         & ! latitude of the rotated north pole (in degrees, N>0)
   pollon,         & ! longitude of the rotated north pole (in degrees, E>0)
   dlat,           & ! grid point distance in zonal direction (in degrees)
   dlon,           & ! grid point distance in meridional direction (in degrees)
   startlat,       & ! transformed latitude of the lower left grid point
                     ! of the total domain (in degrees, N>0)
   startlon          ! transformed longitude of the lower left grid point
                     ! of the total domain (in degrees, E>0)
                                                                                                         
  REAL (KIND=ireals)          ::                                              &
   pollat_in,      & ! latitude of the rotated north pole (in degrees, N>0)
   pollon_in,      & ! longitude of the rotated north pole (in degrees, E>0)
   dlat_in,        & ! grid point distance in zonal direction (in degrees)
   dlon_in,        & ! grid point distance in meridional direction
   startlat_in_tot,& ! transformed latitude of the lower left grid point
                     ! of the total domain (in degrees, N>0)
   startlon_in_tot,& ! transformed longitude of the lower left grid point
                     ! of the total domain (in degrees, E>0)
   endlat_in_tot,  & ! transformed latitude of the upper right grid point
                     ! of the total domain (in degrees, N>0)
   endlon_in_tot,  & ! transformed longitude of the upper right grid point
                     ! of the total domain (in degrees, E>0)
   startlat_in,    & ! transformed latitude of the lower left grid point
                     ! of the local domain (in degrees, N>0)
   startlon_in,    & ! transformed longitude of the lower left grid point
                     ! of the local domain (in degrees, E>0)
   endlat_in,      & ! transformed latitude of the upper right grid point
                     ! of the local domain (in degrees, N>0)
   endlon_in         ! transformed longitude of the upper right grid point
                     ! of the local domain (in degrees, E>0)
 
  INTEGER (KIND=iintegers)    ::                                              &
   ie_in_tot,      & ! ie for input grid, total domain
   je_in_tot,      & ! je for input grid, total domain
   ie_in_max,      & ! Max. of ie_in (local) on all processors
   je_in_max         ! Max. of je_in (local) on all processors
 
  REAL (KIND=ireals), POINTER       ::                                        &
   latlm_m   (:,:), & ! latitudes of the LM grid points
   lonlm_m   (:,:), & ! longitude of the LM grid points
   zlat_out  (:,:), &
   zlon_out  (:,:)
                                                                                                       
  REAL (KIND=ireals), POINTER, SAVE    ::                                     &
   u_lm (:,:,:,:),v_lm (:,:,:,:),D_lm (:,:,:,:),T_lm (:,:,:,:),               &
   rho_lm (:,:,:,:),qv_lm (:,:,:,:),qc_lm (:,:,:,:),z(:)
    
  REAL (KIND=ireals), POINTER, SAVE  ::                                       &
   u_zw (:,:,:,:),v_zw (:,:,:,:),D_zw (:,:,:,:),T_zw (:,:,:,:),               &
   rho_zw (:,:,:,:),qv_zw (:,:,:,:),qc_zw (:,:,:,:),                          &
   u (:,:,:,:),v (:,:,:,:),D (:,:,:,:),TETA (:,:,:,:),                        &
   rho (:,:,:,:),qv (:,:,:,:),qc (:,:,:,:),hp (:,:),                          &
   h (:,:), hstr (:,:,:), h2str(:,:,:), hzw(:)
    
  REAL    (KIND=ireals)    :: dTETA
    
  REAL (KIND=ireals), POINTER, SAVE :: value (:,:)
    
  CHARACTER(128) :: Line

END MODULE parameters


PROGRAM main_interp
!=============================================================================
USE parameters

!Parameterlist:
INTEGER :: I,J,K,TI,X,Y,WAHL,ielm,jelm,ke,te

CHARACTER*180 datname

!--------------------------------------------------------
!  1.)  Stretching/Interpolation 
!--------------------------------------------------------
  
CALL StretchInterp (u,v,D,T,rho,qv,qc,hp,ielm,jelm,ke,te)

!------------------------------------------------------
!  2.)  Output in file asam.in
!------------------------------------------------------
    DO TI=1,te
     WRITE(datname,'(I1,A)')TI,'_asam.in'
     OPEN (TI,file=TRIM(datname))
     DO I=1,ielm
      DO J=1,jelm
       DO WAHL=1,8
        IF (WAHL.EQ.1) THEN
         WRITE(9,*)'uProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),u(TI,I,J,K)
         END DO
         WRITE(9,*) '  '

        ELSEIF (WAHL.EQ.2) THEN
         WRITE(9,*)'vProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),v(TI,I,J,K)
         END DO
         WRITE(9,*) '  '
        
        ELSEIF (WAHL.EQ.3) THEN
       ! WRITE(9,*)'wProf'
       ! WRITE(9,222)ke+1
         DO K=0,ke
       !  WRITE(9,333)hp(I,J),w(TI,I,J,K)
         END DO
       ! WRITE(9,*) '  '

        ELSEIF (WAHL.EQ.4) THEN
         WRITE(9,*)'ThProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),TETA(TI,I,J,K)
         END DO
         WRITE(9,*) '  '
 
        ELSEIF (WAHL.EQ.5) THEN
         WRITE(9,*)'RhoProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),rho(TI,I,J,K)
         END DO
         WRITE(9,*) '  '

        ELSEIF (WAHL.EQ.6) THEN
         WRITE(9,*)'DProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),D(TI,I,J,K)
         END DO
         WRITE(9,*) '  '

        ELSEIF (WAHL.EQ.7) THEN
         WRITE(9,*)'qvProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),qv(TI,I,J,K)
         END DO
         WRITE(9,*) '  '

        ELSEIF (WAHL.EQ.8) THEN
         WRITE(9,*)'qcProf'
         WRITE(9,222)ke
         DO K=1,ke
          WRITE(9,333)hp(I,J),qc(TI,I,J,K)
         END DO
         WRITE(9,*) '  '
        END IF
       END DO
      END DO
     END DO
    CLOSE(10)
    END DO

222 FORMAT (I2)
333 FORMAT (F7.2,F9.3)

END


SUBROUTINE StretchInterp (u,v,D,TETA,rho,qv,qc,hp,ie2lm,je2lm,ke2lm,te)
!=============================================================================
USE parameters

!Parameterlist:
INTEGER, PARAMETER          :: ke_in=22,htop=2000                              

INTEGER (KIND=iintegers)    :: l,l1,l2,j1,j2,di,dj,k,kn,ke2lm,t,te

INTEGER (KIND=iintegers)    ::                                                &
   horizInt,        &
   vertINT
!-----------------------------------------------------------------------------
! End of header
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Begin Subroutine StretchInterp
!-----------------------------------------------------------------------------
  ke2lm=81
  te=2
  dTETA=0.01

CALL ReadNamelist(value,'Namelist IN')
  ie_in=value(3,1)
  je_in=value(3,2)
  
CALL ReadNamelist(value,'Namelist OUT')
  ie2lm=value(3,1)
  je2lm=value(3,2)
  
  

ALLOCATE (h(ie_in,je_in),hstr(ie_in,je_in,ke_in),h2str(ie_in,je_in,ke_in))

ALLOCATE (u_lm(te,ie_in,je_in,ke_in),v_lm(te,ie_in,je_in,ke_in),    &
          D_lm(te,ie_in,je_in,ke_in),T_lm(te, ie_in, je_in, ke_in), & 
          rho_lm(te,ie_in,je_in,ke_in),qv_lm(te,ie_in,je_in,ke_in), &
          qc_lm(te,ie_in,je_in,ke_in)) 

ALLOCATE (u_zw(te,ie_in,je_in,ke2lm),v_zw(te,ie_in,je_in,ke2lm),    &
          D_zw(te,ie_in,je_in,ke2lm),T_zw(te,ie_in,je_in,ke2lm),    &
          rho_zw(te,ie_in,je_in,ke2lm),qv_zw(te,ie_in,je_in,ke2lm), &
          qc_zw (te,ie_in,je_in,ke2lm))
   

ALLOCATE (u(te,ie2lm,je2lm,ke2lm),v(te,ie2lm,je2lm,ke2lm),          &
          D(te,ie2lm,je2lm,ke2lm),TETA(te,ie2lm,je2lm,ke2lm),          &
          rho(te,ie2lm,je2lm,ke2lm),qv(te,ie2lm,je2lm,ke2lm),       &
          qc(te,ie2lm,je2lm,ke2lm),hp(ie2lm,je2lm),hzw(ke2lm)) 

CALL ReadFile (h,hp,z,u_lm,v_lm,D_lm,T_lm,rho_lm,qv_lm,qc_lm)  


CALL organize_grids (x_wght,y_wght,i_index,j_index)

DO t=1,te
DO l2=1,je2lm
 DO l1=1,ie2lm
  j1=i_index(l1,l2)
  j2=j_index(l1,l2)
                                                                                                                                                             
  DO dj=0,1
   DO di=0,1
    DO k=1,ke_in
                                                                                                                                                             
!-----------------------------------------------------------------------------
! Stretch: u,v,D
!-----------------------------------------------------------------------------
     hstr(j1+di,j2+dj,k)=h(j1+di,j2+dj)+z(k)
     h2str(j1+di,j2+dj,k)=((htop-hp(l1,l2))/(htop-h(j1+di,j2+dj)))*hstr(j1+di,j2+dj,k)+ &
                          htop*(1-((htop-hp(l1,l2))/(htop-h(j1+di,j2+dj))))
                                                                                                                                                             
!----------------------------------------------------------------------------
! Vertical Interpolation: u,v,D
!----------------------------------------------------------------------------
     DO kn=1,ke2lm
      hzw(kn)=(kn-1)*20+400
     
      IF (hzw(kn).GE.h2str(j1+di,j2+dj,k).AND.hzw(kn).LE.h2str(j1+di,j2+dj,k+1)) THEN
                                                                                                                                                             
       u_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),h2str(j1+di,j2+dj,k),h2str(j1+di,j2+dj,k+1), &
                                  u_lm(t,j1+di,j2+dj,k),u_lm(t,j1+di,j2+dj,k+1))           
                                                                                                                                                             
       v_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),h2str(j1+di,j2+dj,k),h2str(j1+di,j2+dj,k+1), &
                                v_lm(t,j1+di,j2+dj,k),v_lm(t,j1+di,j2+dj,k+1))
                                                                                                                                                         
       D_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),h2str(j1+di,j2+dj,k),h2str(j1+di,j2+dj,k+1), &
                                D_lm(t,j1+di,j2+dj,k),D_lm(t,j1+di,j2+dj,k+1))
                                                                                                                                                             
      ELSE
       EXIT
      END IF
     END DO
    END DO
   END DO
  END DO

!---------------------------------------------------------------------------
! Horizontal Interpolation: u,v,D
!---------------------------------------------------------------------------
  DO kn=1,ke2lm
 
   u(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),u_zw(t,j1,j2,kn),u_zw(t,j1+1,j2,kn), &
                        u_zw(t,j1,j2+1,kn),u_zw(t,j1+1,j2+1,kn))
 
   v(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),v_zw(t,j1,j2,kn),v_zw(t,j1+1,j2,kn), &
                        v_zw(t,j1,j2+1,kn),v_zw(t,j1+1,j2+1,kn))
 
   D(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),D_zw(t,j1,j2,kn),D_zw(t,j1+1,j2,kn), &
                        D_zw(t,j1,j2+1,kn),D_zw(t,j1+1,j2+1,kn))
  END DO
 
!--------------------------------------------------------------------------
! Vertical Interpolation: T,rho,qv,qc
!--------------------------------------------------------------------------
  DO dj=0,1
   DO di=0,1
    DO k=1,ke_in
     DO kn=1,ke2lm
      IF (hzw(kn).GE.hstr(j1+di,j2+dj,k).AND.hzw(kn).LE.hstr(j1+di,j2+dj,k+1)) THEN
        
       T_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),hstr(j1+di,j2+dj,k),hstr(j1+di,j2+dj,k+1), &
                            T_lm(t,j1+di,j2+dj,k),T_lm(t,j1+di,j2+dj,k+1))
 
       rho_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),hstr(j1+di,j2+dj,k),hstr(j1+di,j2+dj,k+1), &
                              rho_lm(t,j1+di,j2+dj,k),rho_lm(t,j1+di,j2+dj,k+1))
 
       qv_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),hstr(j1+di,j2+dj,k),hstr(j1+di,j2+dj,k+1), &
                             qv_lm(t,j1+di,j2+dj,k),qv_lm(t,j1+di,j2+dj,k+1))
        
       qc_zw(t,j1+di,j2+dj,kn)=vertInt(hzw(kn),hstr(j1+di,j2+dj,k),hstr(j1+di,j2+dj,k+1), &
                             qc_lm(t,j1+di,j2+dj,k),qc_lm(t,j1+di,j2+dj,k+1))
 
      ELSE
       EXIT
      END IF
     END DO
    END DO
   END DO
  END DO

!-------------------------------------------------------------------------
! Horizontal Interpolation: T,rho,qv,qc
!-------------------------------------------------------------------------
  DO kn=1,ke2lm
          
   TETA(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),T_zw(t,j1,j2,kn),T_zw(t,j1+1,j2,kn), &
                        T_zw(t,j1,j2+1,kn),T_zw(t,j1+1,j2+1,kn))
          
   rho(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),rho_zw(t,j1,j2,kn),rho_zw(t,j1+1,j2,kn), &
                        rho_zw(t,j1,j2+1,kn),rho_zw(t,j1+1,j2+1,kn))
          
   qv(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),qv_zw(t,j1,j2,kn),qv_zw(t,j1+1,j2,kn), &
                        qv_zw(t,j1,j2+1,kn),qv_zw(t,j1+1,j2+1,kn))
          
   qc(t,l1,l2,kn)=horizInt(x_wght(l1,l2),y_wght(l1,l2),qc_zw(t,j1,j2,kn),qc_zw(t,j1+1,j2,kn), &
                        qc_zw(t,j1,j2+1,kn),qc_zw(t,j1+1,j2+1,kn))        
  END DO

!------------------------------------------------------------------------
! Fill up: T,rho,qv,qc
!------------------------------------------------------------------------
  DO kn=1,ke2lm
   IF (TETA(t,l1,l2,kn).GT.0.) THEN    
    IF (hzw(kn).EQ.hp(l1,l2)) THEN
     EXIT
    ELSE
     l=kn
     DO WHILE (hzw(l).GT.hp(l1,l2))
      l=l-1      
      TETA(t,l1,l2,l)=TETA(t,l1,l2,kn)+dTETA*(kn-l)*20
     END DO
!     T(t,l1,l2,l)=T(t,l1,l2,kn)+dT*(hzw(l)+20-hp(l1,l2))
    END IF
   ELSE
    EXIT
   END IF
  END DO

  DO kn=1,ke2lm
   IF (rho(t,l1,l2,kn).GT.0.) THEN    
    IF (hzw(kn).EQ.hp(l1,l2)) THEN
     EXIT
    ELSE
     l=kn
     DO WHILE (hzw(l).GT.hp(l1,l2))
      l=l-1      
      rho(t,l1,l2,l)=rho(t,l1,l2,kn)
     END DO
!     rho(t,l1,l2,l)=rho(t,l1,l2,kn)
    END IF
   ELSE
    EXIT
   END IF
  END DO

  DO kn=1,ke2lm
   IF (qv(t,l1,l2,kn).GT.0.) THEN    
    IF (hzw(kn).EQ.hp(l1,l2)) THEN
     EXIT
    ELSE
     l=kn
     DO WHILE (hzw(l).GT.hp(l1,l2))
      l=l-1      
      qv(t,l1,l2,l)=qv(t,l1,l2,kn)
     END DO
!     qv(t,l1,l2,l)=qv(t,l1,l2,kn)
    END IF
   ELSE
    EXIT
   END IF
  END DO
 
  DO kn=1,ke2lm       
   IF (qc(t,l1,l2,kn).LT.1111.) THEN   
    IF (hzw(kn).EQ.hp(l1,l2)) THEN
     EXIT
    ELSE
     l=kn
     DO WHILE (hzw(l).GT.hp(l1,l2))
      l=l-1
      qv(t,l1,l2,l)=qv(t,l1,l2,kn)
!      qv_hyp(l1,l2,l)=e/(rho(l1,l2,l)*T(l1,l2,l))*2.1668e-3
!      dqc(l1,l2,l)=qv_hyp(l1,l2,l)-qv(l1,l2,kn)
!      qc(l1,l2,l)=qc(l1,l2,kn)-dqc(l1,l2,l)
!      IF (qv(l1,l2,l).LE.0.) qv(l1,l2,l)=0
     END DO
!     qv(l1,l2,l)=qv(l1,l2,kn)
    END IF
   ELSE
    EXIT
   END IF    
  END DO

 END DO
END DO
END DO

END SUBROUTINE StretchInterp


SUBROUTINE ReadNamelist(value,Name)
!=======================================================================
  USE parameters

  INTEGER :: i,n,InputUnit
  CHARACTER(*) :: Name

  InputUnit=1
  OPEN(UNIT=InputUnit,FILE='Namelist',STATUS='OLD')
  DO 
    READ (InputUnit,*,END=1) Line
    IF (INDEX(Line,TRIM(Name))>0) THEN
      EXIT
    END IF
  END DO
  READ (InputUnit,*) n
  ALLOCATE (value(n,2))
  DO i=1,n
    READ (InputUnit,*) value (i,:)
  END DO
  CLOSE (InputUnit)
1 CONTINUE

END SUBROUTINE ReadNamelist


SUBROUTINE ReadFile (hp,z,u_lm,v_lm,TKVM_lm,T_lm,rho_lm,qv_lm,qc_lm)
!=======================================================================
IMPLICIT NONE
INTEGER :: I,J,K,TI,X,Y,WAHL
REAL ::                        &        
        a ,                    & ! Dummy
        p_lm(24,10,10,250),    & ! Luftdruck [hPa]
        u_lm(24,10,10,250),    & ! zonaler Wind [m/s]
        v_lm(24,10,10,250),    & ! meridionaler Wind [m/s]
        w_lm(24,10,10,250),    & ! vertikale Windgeschw. [m/s]
        qc_lm(24,10,10,250),   & ! Wolkenfluessigwasser [g/kg]
        qv_lm(24,10,10,250),   & ! Wasserdampfgehalt [g/kg]
        T_lm(24,10,10,250),    & ! Temperatur [Grad Celsius/K]
        TETA_lm(24,10,10,250), & ! Potent. Temperatur [K] 
        TKVM_lm(24,10,10,250), & ! Turb. Diff.-koeff. fuer Impuls [m**2/s]
        TKVH_lm(24,10,10,250), & ! Turb. Diff.-koeff. fuer Waerme/Feuchte [m**2/s]
        rh_lm(24,10,10,250),   & ! Rel. Feuchte [%] 
        rho_lm(24,10,10,250),  & ! Luftdichte [kg/m**3]
        hp(10,10),             & ! Hoehe LM-Punkt ueber NN [m]
        hoehe(24,10,10,250),   & ! Hoehe LM-Punkt ueber NN [m]
        lat_lm(10,10),         & ! geogr. Breite LM-Punkt [Grad]
        lon_lm(10,10),         & ! geogr. Laenge LM-Punkt [Grad]
        z(51)                    ! Hoehe der LM-z-Level [m]
INTEGER,PARAMETER :: S=29,E=52,PANZ=6,hours0=1,hours=1,SSP=438,ESP=460      !474,496
INTEGER,PARAMETER:: ANZW=22
CHARACTER*180 test

!-------------------------------------------------
!  1.)  Open files
!-------------------------------------------------
OPEN (10,file='01_10_02utc_YUPRGRPT')
OPEN (11,file='YUSPECIF')
 
!---------------------------------------------------------
!  2.) Read 'YUSPECIF','YUPRGRPT'
!---------------------------------------------------------
13 FORMAT (A180)
      DO I=1,ESP
       READ(11,13) test
         IF (I.GE.SSP) THEN
         READ(test,*)a,a,z(ANZW-(I-SSP))
         PRINT*,z(ANZW-(I-SSP)),(ANZW-(I-SSP))
         END IF
      END DO
 
      DO I=1,PANZ
       DO J=1,PANZ
        DO TI=1,hours
         DO K=1,E
          READ(10,13) test
          IF (K.GE.S .AND. K.LT.E) THEN
           READ(test,*)hp(I,J),lat_lm(I,J),lon_lm(I,J),p_lm(TI,I,J,(ANZW-(K-S))), &
                       T_lm(TI,I,J,(ANZW-(K-S))),rho_lm(TI,I,J,(ANZW-(K-S))),     &
                       qv_lm(TI,I,J,(ANZW-(K-S))),qc_lm(TI,I,J,(ANZW-(K-S))),     &
                       rh_lm(TI,I,J,(ANZW-(K-S))),a,a,u_lm(TI,I,J,(ANZW-(K-S))),  &
                       v_lm(TI,I,J,(ANZW-(K-S))),w_lm(TI,I,J,(ANZW-(K-S))),a,     & 
                       TKVM_lm(TI,I,J,(ANZW-(K-S))),TKVH_lm(TI,I,J,(ANZW-(K-S)))

!------------------------------------------------------
!  2.1) Conversion of temperature
!------------------------------------------------------
           T_lm(TI,I,J,(ANZW-(K-S)))=T_lm(TI,I,J,(ANZW-(K-S)))+273.15

!------------------------------------------------------
!  2.2) Conversion into Potentielle Temperatur
!------------------------------------------------------
           TETA_lm(TI,I,J,(ANZW-(K-S)))=T_lm(TI,I,J,(ANZW-(K-S)))*(1000/          &
                                        p_lm(TI,I,J,(ANZW-(K-S))))**0.286
!------------------------------------------------------
!  2.3) Calculation of air density 
!------------------------------------------------------
         ! rho_lm(TI,I,J,(ANZW-(K-S)))=287.*(1+0.608*rh_lm(TI,I,J,(ANZW-(K-S)))/100.)* &
         ! T_lm(TI,I,J,(ANZW-(K-S)))/(p_lm(TI,I,J,(ANZW-(K-S)))*100.)
          END IF
         END DO
        END DO
       END DO
      END DO

CLOSE(10)
CLOSE(11)

END SUBROUTINE ReadFile


SUBROUTINE organize_grids (x_wght,y_wght,i_index,j_index)
!==============================================================================
!------------------------------------------------------------------------------
! Description:
!   This external subroutine organizes the generation of the different grids
!   and computes the grid correspondences.
!------------------------------------------------------------------------------
USE parameters
!Parameterlist:
                      
INTEGER  :: nzbounds,izstart_in,izend_in,jzstart_in,jzend_in,i,j,n,nzstat,izerror

REAL    (KIND=ireals)    :: &
  zlats, zlons,             &
  zminlat_m, zmaxlat_m, zminlon_m, zmaxlon_m,  &
  zminlat,   zmaxlat,   zminlon,   zmaxlon,    &
  phirot2phi, rlarot2rla, phi2phirot, rla2rlarot
    

! Errors
CHARACTER (LEN= 80)       ::  &
  yzerror       ! error message

CHARACTER (LEN= 80)       ::  &
  yerrmsg      ! error message

INTEGER (KIND=iintegers)  ::  &
  ierror       ! error status

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine organize_grids
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Section 1: Initializations
!------------------------------------------------------------------------------

  ierror  = 0
  yerrmsg = '         '
  izerror = 0
  yzerror = '         '

!------------------------------------------------------------------------------
! Section 2: Allocate and compute geographical coordinates (latitude and
!            longitude) for ASAM output grid
!------------------------------------------------------------------------------
CALL ReadNamelist(value,'Namelist IN')
  pollat_in=value(1,1)
  pollon_in=value(1,2)
  startlat_in_tot=value(2,1)
  startlon_in_tot=value(2,2)
  ie_in_tot=value(3,1)
  je_in_tot=value(3,2)
  dlat_in=value(4,1)
  dlon_in=value(4,2)

CALL ReadNamelist(value,'Namelist OUT')
  pollat=value(1,1)
  pollon=value(1,2)
  startlat=value(2,1)
  startlon=value(2,2)
  ie2lm=value(3,1)
  je2lm=value(3,2)
  dlat=value(4,1)
  dlon=value(4,2)

  !  Compute the latitude and longitude of the LM-grid
  !  (Arakawa C-grid) at the locations of:
  !   - mass points (arrays with "subnames" '_m')
  ALLOCATE (latlm_m(ie2lm,je2lm),  &
            lonlm_m(ie2lm,je2lm),  &
            zlat_out(ie2lm,je2lm), &
            zlon_out(ie2lm,je2lm), &
            STAT=nzstat)
  IF (nzstat /= 0) THEN
    ierror  = 1
    yerrmsg = 'cannot allocate lat/lon in organize_grids'
    RETURN
  ENDIF

  ! geographical coordinates for the mass grid points
  DO j = 1,je2lm
    zlats = startlat + (j-1)*dlat
    DO i = 1,ie2lm
      zlons = startlon + (i-1)*dlon
      IF(zlons > 180.0_ireals) THEN
        zlons = zlons - 360.0_ireals
      ENDIF
      latlm_m(i,j)   = phirot2phi(zlats, zlons, pollat, pollon)
      lonlm_m(i,j)   = rlarot2rla(zlats, zlons, pollat, pollon)
    ENDDO
  ENDDO

!------------------------------------------------------------------------------
! Section 3: Generation and decomposition of coarse grid (not GME)
!------------------------------------------------------------------------------
  ! This subroutine computes the start- and end-indices of the grid points 
  ! from the coarse grid that are necessary to interpolate to the fine grid.
  ! We add 2 additional rows and columns for the interpolation.
  ! The data from the coarse grid then is decomposed accordingly.
  
  nzbounds = 2     ! this number can be changed if it turns out that
                   ! the subdomains of the coarse grid are not big enough

  !----------------------------------------------------------------------------
  !- Section 3.1: Compute latitudes and longitudes of the fine grid points 
  !               in the coordinates of the coarse system
  !----------------------------------------------------------------------------

  IF ((pollat_in == 90.0_ireals) .AND. (pollon_in == 180.0_ireals)) THEN
    ! coarse grid is unrotated
    ! zlat_out = latlm_* and zlon_out = lonlm_*
    zlat_out (:,:)= latlm_m(:,:)
    zlon_out (:,:) = lonlm_m(:,:)
  ELSEIF ((pollat_in == pollat) .AND. (pollon_in == pollon)) THEN
    ! coarse grid has the same rotated pole like the fine grid
    DO j = 1, je2lm
      zlat_out (:,j) = startlat + (j-1) * dlat      
    ENDDO
    DO i = 1, ie2lm
      zlon_out (i,:) = startlon + (i-1) * dlon      
    ENDDO
  ELSE
    ! coarse grid has a different rotated pole from the fine grid
    ! rotate the geographical coordinates latlm_* and lonlm_* to the 
    ! rotation of the coarse system
    DO j = 1, je2lm
      DO i = 1, ie2lm
        zlat_out (i,j) = phi2phirot (latlm_m(i,j), lonlm_m(i,j), pollat_in, pollon_in)        
        zlon_out (i,j) = rla2rlarot (latlm_m(i,j), lonlm_m(i,j), pollat_in, pollon_in)        
      ENDDO
    ENDDO
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 3.2: Determine MIN and MAX of latitude and longitude
  !               and check whether coarse grid is big enough
  !----------------------------------------------------------------------------

  ! for the mass grid points
  zminlat = MINVAL (zlat_out(1:ie2lm,1:je2lm))
  zmaxlat = MAXVAL (zlat_out(1:ie2lm,1:je2lm))
  zminlon = MINVAL (zlon_out(1:ie2lm,1:je2lm))
  zmaxlon = MAXVAL (zlon_out(1:ie2lm,1:je2lm))
  
  ! Check whether coarse grid is big enough
  IF (zminlat-nzbounds*dlat_in < startlat_in_tot) THEN
    ierror  = 2
    yerrmsg = 'Coarse domain is not big enough'
  ENDIF
  IF (zmaxlat+nzbounds*dlat_in > (startlat_in_tot+(je_in_tot-1)*dlat_in)) THEN
    ierror  = 3
    yerrmsg = 'Coarse domain is not big enough'
  ENDIF
  IF (zminlon-nzbounds*dlon_in < startlon_in_tot) THEN
    ierror  = 4
    yerrmsg = 'Coarse domain is not big enough'
  ENDIF
  IF (zmaxlon+nzbounds*dlon_in > (startlon_in_tot+(ie_in_tot-1)*dlon_in)) THEN
    ierror  = 5
    yerrmsg = 'Coarse domain is not big enough'
  ENDIF

  !----------------------------------------------------------------------------
  !- Section 3.3: Indices of coarse grid covering this part of fine grid
  !----------------------------------------------------------------------------

  ! looking for southern grid point
  DO j = 1, je_in_tot
    zlats  = startlat_in_tot + (j-1) * dlat_in
    IF (zlats >= zminlat) THEN
      jzstart_in = j-nzbounds-1
      EXIT
    ENDIF
  ENDDO

  ! looking for northern grid point
  DO j = 1, je_in_tot
    zlats  = startlat_in_tot + (j-1) * dlat_in
    IF (zlats >  zmaxlat) THEN
      jzend_in   = j+nzbounds
      EXIT
    ENDIF
  ENDDO

  ! looking for western grid point
  DO i = 1, ie_in_tot
    zlons  = startlon_in_tot + (i-1) * dlon_in
    IF (zlons >= zminlon) THEN
      izstart_in = i-nzbounds-1
      EXIT
    ENDIF
  ENDDO

  ! looking for eastern grid point
  DO i = 1, ie_in_tot
    zlons  = startlon_in_tot + (i-1) * dlon_in
    IF (zlons >  zmaxlon) THEN
      izend_in   = i+nzbounds
      EXIT
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
! Section 3.4: Determine size of this part of coarse grid
!------------------------------------------------------------------------------

  ! start- and end-latitude and -longitude

  startlat_in = startlat_in_tot + (jzstart_in-1) * dlat_in
  endlat_in   = startlat_in_tot + (jzend_in  -1) * dlat_in
  startlon_in = startlon_in_tot + (izstart_in-1) * dlon_in
  endlon_in   = startlon_in_tot + (izend_in  -1) * dlon_in

!------------------------------------------------------------------------------
! Section 4: Generation of coarse and fine grid correspondence 
!------------------------------------------------------------------------------

  ! Allocate memory
  ALLOCATE (i_index(ie2lm,je2lm), j_index(ie2lm,je2lm),   STAT=nzstat)
  ALLOCATE (x_wght (ie2lm,je2lm), y_wght (ie2lm,je2lm),   STAT=nzstat)

  ! For every LM (fine mesh) grid point the grid point of the coarse grid
  ! is determined that is just to the lower left of the LM grid point.
  ! This is done also for the staggered positions, which gives 5 possibilities:

  !  1:  mass grid points:
  DO j = 1, je2lm
    DO i = 1, ie2lm
      i_index(i,j) = INT ((zlon_out(i,j) - startlon_in) / dlon_in) + 1
      j_index(i,j) = INT ((zlat_out(i,j) - startlat_in) / dlat_in) + 1
      x_wght (i,j) = (zlon_out(i,j) - startlon_in) / dlon_in + 1.0_ireals - &
                         REAL (i_index(i,j), ireals)
      y_wght (i,j) = (zlat_out(i,j) - startlat_in) / dlat_in + 1.0_ireals - &
                         REAL (j_index(i,j), ireals)
    ENDDO
  ENDDO
     
END SUBROUTINE organize_grids

!******************************************************************************
!******************************************************************************

FUNCTION  phirot2phi ( phirot, rlarot, polphi, pollam, polgam )
                                                                                                                                                             
!-----------------------------------------------------------------------------
! Description:
!   This function converts phi from one rotated system to phi in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!-----------------------------------------------------------------------------                                                                                                                                                         
                                                                                                                                                          
!Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system
                                                                                                                                                             
REAL, INTENT (IN), OPTIONAL      ::        &
  polgam      ! angle between the north poles of the systems
                                                                                                                                                             
REAL                   ::        &
  phirot2phi  ! latitude in the geographical system
                                                                                                                                                             
! Local variables
REAL                   ::        &
  zsinpol, zcospol, zphis, zrlas, zarg, zgam
                                                                                                                                                             
REAL                   ::        &
  zrpi18 = 57.2957795,                  &
  zpir18 = 0.0174532925
                                                                                                                                                               
!------------------------------------------------------------------------------           
! Begin function phirot2phi
!------------------------------------------------------------------------------
              
  zsinpol     = SIN (zpir18 * polphi)
  zcospol     = COS (zpir18 * polphi)
              
  zphis       = zpir18 * phirot
  IF (rlarot > 180.0) THEN
    zrlas = rlarot - 360.0
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas       = zpir18 * zrlas
              
  IF ( PRESENT (polgam) ) THEN
    zgam  = zpir18 * polgam
    zarg  = zsinpol*SIN (zphis) +                                           &
        zcospol*COS(zphis) * ( COS(zrlas)*COS(zgam) - SIN(zgam)*SIN(zrlas) )
  ELSE
    zarg  = zcospol * COS (zphis) * COS (zrlas) + zsinpol * SIN (zphis)
  ENDIF
              
  phirot2phi  = zrpi18 * ASIN (zarg)
              
END FUNCTION phirot2phi

!******************************************************************************
!******************************************************************************

FUNCTION  phi2phirot ( phi, rla, polphi, pollam )
              
!------------------------------------------------------------------------------
! Description:
!   This routine converts phi from the real geographical system to phi
!   in the rotated system.
!------------------------------------------------------------------------------

! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the rotated system
  rla        ! longitude in the rotated system
              
REAL                   ::        &
  phi2phirot ! longitude in the rotated system
              
! Local variables
REAL                       ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
              
REAL                       ::    &
  zrpi18 = 57.2957795,       & !
  zpir18 = 0.0174532925
              
!------------------------------------------------------------------------------            
! Begin function phi2phirot
!------------------------------------------------------------------------------
              
  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0) THEN
    zrla1  = rla - 360.0
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1
              
  zarg1    = SIN (zphi) * zsinpol
  zarg2    = COS (zphi) * zcospol * COS (zrla - zlampol)
              
  phi2phirot = zrpi18 * ASIN (zarg1 + zarg2)
              
END FUNCTION phi2phirot
              
!******************************************************************************
!******************************************************************************

FUNCTION  rlarot2rla (phirot, rlarot, polphi, pollam, polgam)
              
!------------------------------------------------------------------------------
! Description:
!   This function converts lambda from one rotated system to lambda in another
!   system. If the optional argument polgam is present, the other system
!   can also be a rotated one, where polgam is the angle between the two
!   north poles.
!   If polgam is not present, the other system is the real geographical
!   system.
!------------------------------------------------------------------------------

! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,   & ! latitude of the rotated north pole
  pollam,   & ! longitude of the rotated north pole
  phirot,   & ! latitude in the rotated system
  rlarot      ! longitude in the rotated system
              
REAL, INTENT (IN), OPTIONAL      ::        &
  polgam      ! latitude of the rotated north pole
              
REAL                   ::        &
  rlarot2rla  ! latitude in the geographical system
              
! Local variables
REAL                   ::        &
  zsinpol, zcospol, zlampol, zphis, zrlas, zarg1, zarg2, zgam
              
REAL                   ::        &
  zrpi18 = 57.2957795,   & !
  zpir18 = 0.0174532925

!------------------------------------------------------------------------------             
! Begin function rlarot2rla
!------------------------------------------------------------------------------
              
  zsinpol = SIN (zpir18 * polphi)
  zcospol = COS (zpir18 * polphi)
              
  zlampol = zpir18 * pollam
  zphis   = zpir18 * phirot
  IF (rlarot > 180.0) THEN
    zrlas = rlarot - 360.0
  ELSE
    zrlas = rlarot
  ENDIF
  zrlas   = zpir18 * zrlas
              
  IF ( PRESENT(polgam) ) THEN
    zgam    = zpir18 * polgam
    zarg1   = SIN (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    - COS (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
              
    zarg2   = COS (zlampol) *                                                &
      (- zsinpol*COS(zphis) * (COS(zrlas)*COS(zgam) - SIN(zrlas)*SIN(zgam))  &
       + zcospol * SIN(zphis))                                               &
    + SIN (zlampol)*COS(zphis) * (SIN(zrlas)*COS(zgam) + COS(zrlas)*SIN(zgam))
  ELSE
    zarg1   = SIN (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) -    &
              COS (zlampol) *             SIN(zrlas) * COS(zphis)
    zarg2   = COS (zlampol) * (-zsinpol * COS(zrlas) * COS(zphis)  +    &
                                zcospol *              SIN(zphis)) +   &
              SIN (zlampol) *             SIN(zrlas) * COS(zphis)
  ENDIF
              
  IF (zarg2 == 0.0) zarg2 = 1.0E-20
              
  rlarot2rla = zrpi18 * ATAN2(zarg1,zarg2)
              
END FUNCTION rlarot2rla

!******************************************************************************
!******************************************************************************
                                                                             
FUNCTION  rla2rlarot ( phi, rla, polphi, pollam )
                                                                                                                                                             
!------------------------------------------------------------------------------
! Description:
!   This routine converts lambda from the real geographical system to lambda
!   in the rotated system.
!
! Method:
!   Transformation formulas for converting between these two systems.
!------------------------------------------------------------------------------
!
! Parameter list:
REAL, INTENT (IN)      ::        &
  polphi,  & ! latitude of the rotated north pole
  pollam,  & ! longitude of the rotated north pole
  phi,     & ! latitude in the rotated system
  rla        ! longitude in the rotated system
                                                                                                                                                             
REAL                   ::        &
  rla2rlarot ! latitude in the geographical system
                                                                                                                                                             
! Local variables
REAL                        ::    &
  zsinpol, zcospol, zlampol, zphi, zrla, zarg1, zarg2, zrla1
           
REAL                        ::    &
  zrpi18 = 57.2957795,       & !
  zpir18 = 0.0174532925
           
!------------------------------------------------------------------------------
! Begin function rla2rlarot
!------------------------------------------------------------------------------
           
  zsinpol  = SIN (zpir18 * polphi)
  zcospol  = COS (zpir18 * polphi)
  zlampol  =      zpir18 * pollam
  zphi     =      zpir18 * phi
  IF (rla > 180.0) THEN
    zrla1  = rla - 360.0
  ELSE
    zrla1  = rla
  ENDIF
  zrla     = zpir18 * zrla1
           
  zarg1    = - SIN (zrla-zlampol) * COS(zphi)
  zarg2    = - zsinpol * COS(zphi) * COS(zrla-zlampol) + zcospol * SIN(zphi)
           
  IF (zarg2 == 0.0) zarg2 = 1.0E-20
           
  rla2rlarot = zrpi18 * ATAN2 (zarg1,zarg2)
           
END FUNCTION rla2rlarot
           
!******************************************************************************
!******************************************************************************

FUNCTION vertInt(zM,zL,zR,cL,cR)
  REAL :: vertInt,zM,zL,zR,cL,cR
  vertInt=((zM-zL)*cR+(zR-zM)*cL)/(zR-zL)
END FUNCTION vertInt

!******************************************************************************
!******************************************************************************

FUNCTION horizInt(x_w,y_w,var,varR,varO,varOR)
  REAL :: horizInt,x_w,y_w,var,varR,varO,varOR
  horizInt=(1-x_wght)*(1-y_wght)*var+x_wght*(1-y_wght)*varR+ &
           (1-x_wght)*y_wght*varO+x_wght*y_wght*varOR
END FUNCTION horizInt
