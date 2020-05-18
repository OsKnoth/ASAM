  MODULE Control_Mod

  USE Kind_Mod

  IMPLICIT NONE
  !Input GridDistsCtrl
  REAL(8) :: dist_fscv=1.0d-12     ! for fine scaling dist to in_out-def, CheckVertex 
  REAL(8) :: distx_coeff=0.01      ! Distance Point(x) coefficient
  REAL(8) :: disty_coeff=0.01      ! Distance Point(y) coefficient
  REAL(8) :: distz_coeff=0.01      ! Distance Point(z) coefficient
  REAL(8) :: dxViewLoc=1.0d-8      ! Distance Point(x) coefficient of model border
  REAL(8) :: dyViewLoc=1.0d-8      ! Distance Point(y) coefficient of model border
  REAL(8) :: dzViewLoc=1.0d-8      ! Distance Point(z) coefficient of model border
  INTEGER :: IncrVol=10            ! counter for cuts of volume-analysis (1 für Einheitswürfel 
  REAL(8) :: dist_scMaxCell=1.0d-12 ! adjustment value to the filter(screen) of cells 
                                    ! with roughly max. Vol
  REAL(8) :: RelVol=1.d-6           ! Adjust volume to zero
  REAL(8) :: ShiftPoint=1.0d0
  INTEGER :: Depth=8
  INTEGER :: NumCoarse=0
  INTEGER :: NumSmooth=0

  !Input OutGMVControl
  CHARACTER :: out_wahlgrid          ! C/G -> Cartesien,Globe
  CHARACTER :: out_type              ! a/b -> ascii,binary
  CHARACTER(9) :: conv_ingrid        ! convert input grid -> spherical,
  CHARACTER :: invalue_to_out        ! convert 'rad' to input value (z.Zt.'geo') 
                                     ! for gmv-Output  
                                     !                    -> nonspherical (cartesian)
  REAL(8) :: RadOutput               !Erdradius
                                     !RadOutput, value to numerator for parametrization
  REAL(8) :: ScaleRad                !ScaleRad, value to denominator for parametrization
  REAL(8) :: ScaleSoil               !ScaleSoil, value scaling dzi_soil /part, std. 1%
  REAL(8) :: zRauh=0.1d0             !Roughness length
  INTEGER :: MoveXYGrid              !MoveXYGrid, move grid to scale-x-y, 0 .or. 1  as factor
  REAL(8) :: x0DomainHaus
  REAL(8) :: y0DomainHaus
  REAL(8) :: x1DomainHaus
  REAL(8) :: y1DomainHaus
  REAL(8) :: x0Tree
  REAL(8) :: y0Tree
  REAL(8) :: z0Tree
  REAL(8) :: x1Tree
  REAL(8) :: y1Tree
  REAL(8) :: z1Tree
  REAL(8) :: xOffset=0.0d0
  REAL(8) :: yOffset=0.0d0

  !Input GridFileOut
  LOGICAL :: WNull         ! -> *.WNull
  LOGICAL :: GCut          ! -> *.Cut.out.gmvG
  LOGICAL :: GCut2         ! -> *.Cut2.out.gmvG
  LOGICAL :: GSoil         ! -> *.Soil.out.gmvG
  LOGICAL :: GONull        ! -> *.ONull.out.gmvG
  LOGICAL :: GOro          ! -> *.Oro.out.gmvG
  LOGICAL :: Bound         ! -> *.bound
  LOGICAL :: Pbound        ! -> *.pva.bound
  LOGICAL :: Pgall         ! -> *.pva.gall
  LOGICAL :: Ptropo        ! -> *.pva.tropo   
  LOGICAL :: Emission      ! -> *.Emission   
  LOGICAL :: gmv           ! -> *.gmv   
  LOGICAL :: vtk           ! -> *.vtk   


  NAMELIST /GridDistsCtrl/ dist_fscv, distx_coeff, disty_coeff, distz_coeff, &
                           dxViewLoc, dyViewLoc, dzViewLoc, IncrVol, &
                           dist_scMaxCell,RelVol,ShiftPoint,Depth,NumCoarse,NumSmooth

  NAMELIST /OutGMVControl/ out_wahlgrid   &
                          ,out_type       &
                          ,conv_ingrid    &
                          ,invalue_to_out &
                          ,RadOutput      &
                          ,ScaleRad       &
                          ,ScaleSoil      &
                          ,zRauh          &
                          ,MoveXYGrid     &
                          ,x0DomainHaus   &
                          ,y0DomainHaus   &
                          ,x1DomainHaus   &
                          ,y1DomainHaus   &
                          ,x0Tree         &
                          ,y0Tree         &
                          ,z0Tree         &
                          ,x1Tree         &
                          ,y1Tree         &
                          ,z1Tree         &
                          ,xOffset        &
                          ,yOffset

  NAMELIST /GridFileOut/   WNull,GCut,GCut2,GSoil,GONull,GOro, &
                           Bound,Pbound,Pgall,Ptropo,Emission, &
                           gmv,vtk

  INTEGER, PRIVATE :: InputUnit=10                         

CONTAINS

SUBROUTINE InputControl(FileName)
  CHARACTER(*) :: FileName
  CALL input_GMVControl(FileName)
  CALL input_dists_ctrl(FileName)
  CALL input_GridFileOut(FileName)
END SUBROUTINE InputControl

SUBROUTINE input_GMVControl(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  out_wahlgrid='C'
  out_type='b'
  conv_ingrid='n'
  invalue_to_out='n'
  RadOutput=0.0d0
  ScaleRad=1.0d0 
  ScaleSoil=1.0d0
  MoveXYGrid=0
  x0DomainHaus=-1.d20
  y0DomainHaus=-1.d20   
  x1DomainHaus=1.d20   
  y1DomainHaus=1.d20
  x0Tree=1.d20
  y0Tree=1.d20   
  z0Tree=1.d20   
  x1Tree=-1.d20   
  y1Tree=-1.d20
  z1Tree=-1.d20
  WRITE(*,*) 'x0Tree 1',x0Tree
  WRITE(*,*) 'y0Tree 1',y0Tree
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&OutGMVControl')>0) THEN
      BACKSPACE(InputUnit)
      WRITE(*,*) 'Line ',Line(1:50)
      READ(InputUnit,NML=OutGMVControl)
  WRITE(*,*) 'x0Tree 2',x0Tree
  WRITE(*,*) 'y0Tree 2',y0Tree
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

END SUBROUTINE input_GMVControl

SUBROUTINE input_dists_ctrl(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&GridDistsCtrl')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=GridDistsCtrl)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
END SUBROUTINE input_dists_ctrl

SUBROUTINE input_GridFileOut(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  WNull    =.TRUE.    ! -> *.WNull
  GCut     =.TRUE.    ! -> *.Cut.out.gmvG
  GCut2    =.TRUE.    ! -> *.Cut2.out.gmvG
  GSoil    =.TRUE.    ! -> *.Soil.out.gmvG
  GONull   =.TRUE.    ! -> *.ONull.out.gmvG
  GOro     =.TRUE.    ! -> *.Oro.out.gmvG
  Bound    =.TRUE.    ! -> *.bound
  Pbound   =.TRUE.    ! -> *.pva.bound
  Pgall    =.TRUE.    ! -> *.pva.gall
  Ptropo   =.TRUE.    ! -> *.pva.tropo   
  Emission =.FALSE.   ! -> *.pva.tropo   
  gmv      =.TRUE.    ! -> *.gmv
  vtk      =.TRUE.    ! -> *.vtk
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&GridFileOut')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=GridFileOut)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
END SUBROUTINE input_GridFileOut

END MODULE ControL_Mod

