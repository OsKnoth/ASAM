!###############################################################################
! OutputWeightBlkJoins.f90
!       Output aller Bloecke zusammengefasst als globales Gitter, entsprechend
!       der Angabe Vergröberung/Verfeinerung
!
!       Berechnung nx_ges,ny_ges,nz_ges:  
!              Durchsuchen aller Nachbarbloecke ueber  'ib%nachbars' 
!              in 'east,west,south,north,top,bottom', wenn vorhanden und addieren 
!       Ausgabe dx,dy,dz:
!             Durchsuchen aller Nachbarbloecke in 'east,west,south,north,top,bottom'
!       Ausgabe Faces und Cells:
!             Generell nach aufsteigender Nummerierung Bloecke
!             
!
!################################################################################

SUBROUTINE WriteWeightBlkJoins(FileName)

!------------------------------------------------------------------------------
! Description:
!------------------------------------------------------------------------------

  USE Floor_Mod

  USE F_Mod ,  ONLY:   &
               conv_gk

  USE Parametric_Mod ,  ONLY:  &
               indata_type,OutGrid

  USE Grid_Mod, ONLY: &
               VolFace_XY,VolFace_YZ,VolFace_ZX
!==============================================================================
  IMPLICIT NONE
!==============================================================================
!
! Subroutine arguments:
! --------------------
  CHARACTER*50 :: FileName
!
! Local :
! -----------------------------
  INTEGER :: k,j,i
  INTEGER :: ib,ibs
  INTEGER :: nx_ges,ny_ges,nz_ges
  INTEGER :: in_out
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine WriteWeightBlkJoins
!------------------------------------------------------------------------------

  zRauh=0.1d0
  Albedo=0.5
  Emissivity=0.8

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.WeightBlkJoins',STATUS='unknown')
  !WRITE(*,*) "      ",TRIM(FileName),'.WeightBlkJoins'
  WRITE(*,*) "              * ",TRIM(FileName),'.WeightBlkJoins'
  !IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. &
  !   (indata_type=='geo'.AND.out_wahlgrid=='C').OR. &
  !   (indata_type=='rad'.AND.out_wahlgrid=='C')) THEN

  IF((indata_type=='gk'.AND.conv_gk=='s').OR. &
     (indata_type=='geo'.AND.OutGrid=="Cart").OR. &
     (indata_type=='rad'.AND.OutGrid=="Cart")) THEN
    WRITE(10,*) 'G'
  ELSE
    WRITE(10,*) OutGrid(1:1)
  END IF
  !................................................................................
  WRITE(10,*) Domain%x0,Domain%y0,Domain%z0  !0.0d0
!  nx_ges=0
!  ny_ges=0
!  nz_ges=0
!  DO ib=1,nb
!    CALL Set(Floor(ib))
!    nx_ges=nx_ges+nx
!    ny_ges=ny_ges+ny
!    nz_ges=nz_ges+nz
!  END DO  !ib

  !................................nx_ges.........................................
  ibs=1
  nx_ges=Floor(ibs)%nx
  IF (Floor(ibs)%typee=="ie") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(2)%ntype=='ie') THEN
         nx_ges=nx_ges+Floor(Floor(ibs)%nachbars(2)%ib)%nx
         ibs=Floor(ibs)%nachbars(2)%ib
      ELSE
         ibs=0
      END IF
    END DO
  END IF
  ibs=1
  IF (Floor(ibs)%typew=="iw") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(1)%ntype=="iw") THEN
         nx_ges=nx_ges+Floor(Floor(ibs)%nachbars(1)%ib)%nx
         ibs=Floor(ibs)%nachbars(1)%ib
      ELSE
          ibs=0
      END IF
    END DO
  END IF
  !...............................ny_ges..........................................
  ibs=1
  ny_ges=Floor(ibs)%ny
  IF (Floor(ibs)%typen=="in") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(4)%ntype=="in") THEN
         ny_ges=ny_ges+Floor(Floor(ibs)%nachbars(4)%ib)%ny
         ibs=Floor(ibs)%nachbars(4)%ib
      ELSE
         ibs=0
      END IF
    END DO
  END IF
  ibs=1
  IF (Floor(ibs)%types=="is") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(3)%ntype=="is") THEN
         ny_ges=ny_ges+Floor(Floor(ibs)%nachbars(3)%ib)%ny
         ibs=Floor(ibs)%nachbars(3)%ib
      ELSE
         ibs=0
      END IF
    END DO
  END IF
  !..............................nz_ges............................................. 
  ibs=1
  nz_ges=Floor(ibs)%nz
  IF (Floor(ibs)%typet=="it") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(6)%ntype=="it") THEN
         nz_ges=nz_ges+Floor(Floor(ibs)%nachbars(6)%ib)%nz
         ibs=Floor(ibs)%nachbars(6)%ib
      ELSE
         ibs=0
      END IF
    END DO
  END IF
  ibs=1
  IF (Floor(ibs)%typeb=="ib") THEN
    DO WHILE (ibs>0)
      IF(Floor(ibs)%nachbars(5)%ntype=="ib") THEN
         nz_ges=nz_ges+Floor(Floor(ibs)%nachbars(5)%ib)%nz
         ibs=Floor(ibs)%nachbars(5)%ib
      ELSE
         ibs=0
      END IF
    END DO
  END IF
  !...............................dx()...............................................
  WRITE(10,*) nx_ges

  !DO ib=1,nb
  !  CALL Set(Floor(ib))
  !  DO i=ix0+1,ix1
  !     WRITE(10,*) dx(i)
  !  END DO
  !END DO  !ib

  ib=1
  IF (Floor(ib)%typew=="iw") THEN
    DO WHILE (Floor(ib)%nachbars(1)%ntype=="iw")
         ib=Floor(ib)%nachbars(1)%ib
    END DO
  END IF
  !IF (Floor(ib)%typee=="ie") THEN
    DO WHILE (ib>0)
      CALL Set(Floor(ib))
      DO i=ix0+1,ix1
        WRITE(10,*) dx(i)
      END DO
      IF(Floor(ib)%nachbars(2)%ntype=="ie") THEN
         ib=Floor(ib)%nachbars(2)%ib
      ELSE
         ib=0
      END IF
    END DO
  !END IF
  !...............................dy()...............................................
  WRITE(10,*) ny_ges
  !DO ib=1,nb
  !  CALL Set(Floor(ib))
  !  DO j=iy0+1,iy1
  !    WRITE(10,*) dy(j)
  !  END DO
  !END DO  !ib

  ib=1
  IF (Floor(ib)%types=="is") THEN
    DO WHILE (Floor(ib)%nachbars(3)%ntype=="is")
       ib=Floor(ib)%nachbars(3)%ib
    END DO
  END IF
  !IF (Floor(ib)%typen=="in") THEN
    DO WHILE (ib>0)
      CALL Set(Floor(ib))
      DO j=iy0+1,iy1
        WRITE(10,*) dy(j)
      END DO
      IF(Floor(ib)%nachbars(4)%ntype=="in") THEN
         ib=Floor(ib)%nachbars(4)%ib
      ELSE
         ib=0
      END IF
    END DO
  !END IF
  !...............................dz()...............................................
  WRITE(10,*) nz_ges
  !!WRITE(10,*) nzLang     !!! Filaus
  !DO ib=1,nb
  !  CALL Set(Floor(ib))
  !  !DO k=iz0+1,iz0+nzLang        !!! Filaus
  !  !  WRITE(10,*) domain%dz(k)   !!! Filaus
  !  !END DO                       !!! Filaus
  !  DO k=iz0+1,iz1
  !    WRITE(10,*) dz(k)
  !  END DO
  !END DO  !ib
  ib=1
  IF (Floor(ib)%typeb=="ib") THEN
    DO WHILE (Floor(ib)%nachbars(5)%ntype=="ib")
      ib=Floor(ib)%nachbars(5)%ib
    END DO
  END IF
  !IF (Floor(ib)%typet=="it") THEN
    DO WHILE (ib>0)
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        WRITE(10,*) dz(k)
      END DO
      IF(Floor(ib)%nachbars(6)%ntype=="it") THEN
         ib=Floor(ib)%nachbars(6)%ib
      ELSE
         ib=0
      END IF
    END DO
  !END IF
  !............................FacesYZ...............................................
  WRITE(10,*) NrW_All_FYZ, "   NrW_All_FYZ"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0,ix1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            IF ((Faces_YZ(i,j,k)%Face%in_out<4.AND.Faces_YZ(i,j,k)%Face%ec>0).OR. &
                (Faces_YZ(i,j,k)%Face%in_out<0.AND.Faces_YZ(i,j,k)%Face%ec==-1).OR. &
                 Faces_YZ(i,j,k)%Face%Vol==0.0d0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_YZ(i,j,k)%Face%Vol
            END IF
          ELSE
            in_out=Vertices(i,j-1,k-1)%in_out &
                  +Vertices(i,j,k-1)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) 0.0d0
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  WRITE(10,*) NrMP_All_FYZ, "   NrMP_All_FYZ"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0,ix1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            IF (Faces_YZ(i,j,k)%Face%mp>0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_YZ(i,j,k)%Face%MidPoint
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  !............................FacesZX...............................................
  WRITE(10,*) NrW_All_FZX, "   NrW_All_FZX"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            IF ((Faces_ZX(i,j,k)%Face%in_out<4.AND.Faces_ZX(i,j,k)%Face%ec>0).OR. &
                (Faces_ZX(i,j,k)%Face%in_out<0.AND.Faces_ZX(i,j,k)%Face%ec==-1).OR. &
                 Faces_ZX(i,j,k)%Face%Vol==0.0d0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_ZX(i,j,k)%Face%Vol
            END IF
          ELSE
            in_out=Vertices(i-1,j,k-1)%in_out &
                  +Vertices(i,j,k-1)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) 0.0d0
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  WRITE(10,*) NrMP_All_FZX, "   NrMP_All_FZX"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            IF (Faces_ZX(i,j,k)%Face%mp>0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_ZX(i,j,k)%Face%MidPoint
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  !............................FacesXY...............................................
  WRITE(10,*) NrW_All_FXY, "   NrW_All_FXY"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            IF ((Faces_XY(i,j,k)%Face%in_out<4.AND.Faces_XY(i,j,k)%Face%ec>0).OR. &
                (Faces_XY(i,j,k)%Face%in_out<0.AND.Faces_XY(i,j,k)%Face%ec==-1).OR. &
                 Faces_XY(i,j,k)%Face%Vol==0.0d0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_XY(i,j,k)%Face%Vol
            END IF
          ELSE
            in_out=Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) 0.0d0
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  WRITE(10,*) NrMP_All_FXY, "   NrMP_All_FXY"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            IF (Faces_XY(i,j,k)%Face%mp>0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Faces_XY(i,j,k)%Face%MidPoint
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  !..................................Cells...............................................
  WRITE(10,*) NrW_All_Cells, "    NrW_All_Cells"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            IF (Cell(i,j,k)%Cell%in_out<8) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Cell(i,j,k)%Cell%Vol
              !WRITE(10,*) Cell(i,j,k)%Cell%Vol,Cell(i,j,k)%Cell%MidPoint
            END IF
          ELSE
            in_out=Vertices(i-1,j-1,k-1)%in_out &
                  +Vertices(i,j-1,k-1)%in_out &
                  +Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i-1,j,k-1)%in_out &
                  +Vertices(i-1,j,k)%in_out   &
                  +Vertices(i,j-1,k)%in_out   &
                  +Vertices(i,j,k-1)%in_out   &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) 0.0d0
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  WRITE(10,*) NrMP_All_Cells, "    NrMP_All_Cells"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            IF (Cell(i,j,k)%Cell%in_out<8.AND.Cell(i,j,k)%Cell%Vol>0.0d0.AND. &
                  (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6)<1.0d0-1.d-2)) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Cell(i,j,k)%Cell%MidPoint
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  WRITE(10,*) NrB_All_Cells, "    NrB_All_Cells"
  DO ib=1,nb
    CALL Set(Floor(ib))
    !WRITE(10,*) ib, " Block"
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            !IF (OutGrid=="Globe")
            !ELSE
            !END IF
            IF (Cell(i,j,k)%Cell%in_out<8.AND.Cell(i,j,k)%Cell%Vol>0.0d0.AND. &
                  (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5)<1.0d0-1.d-2 &
               .OR.Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6)<1.0d0-1.d-2)) THEN
              WRITE(10,*) i,j,k
              IF (xP(i-1)+0.5d0*dx(i)<=-2.0d4) THEN
!               WRITE(10,'(6d15.7)') 1.0d-99,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
                WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              ELSE IF (xP(i-1)+0.5d0*dx(i)<=-1.0d4) THEN
!               WRITE(10,'(6d15.7)') (xP(i-1)+0.5d0*dx(i)+2.0d4)/1.0d4*zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
                WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              ELSE
                WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  CLOSE(10)

END SUBROUTINE WriteWeightBlkJoins
!##########################################################################################
