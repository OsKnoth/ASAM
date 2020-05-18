!###############################################################################
! OutputWeightBlk_Mod.f90
!
!
!------------------------------------------------------------------------------
! Description:
!------------------------------------------------------------------------------
!################################################################################

MODULE OutputWeightBlk_Mod

  USE Floor_Mod

  USE F_Mod ,  ONLY:   &
               conv_gk

  USE Parametric_Mod ,  ONLY:  &
               indata_type,OutGrid

  USE Grid_Mod, ONLY: &
               VolFace_XY,VolFace_YZ,VolFace_ZX

  USE IOControl_Mod

  IMPLICIT NONE

!================================================================================
CONTAINS
!================================================================================

SUBROUTINE CountWeightFaceYZ(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Floor(ib)%Faces_YZ(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
        NrRW_FacesYZ=NrRW_FacesYZ+1
     END IF
   ELSE
        in_out=Vertices(i,j-1,k-1)%in_out &
              +Vertices(i,j,k-1)%in_out &
              +Vertices(i,j-1,k)%in_out &
              +Vertices(i,j,k)%in_out
        IF (in_out<=0) THEN
           NrRW_FacesYZ=NrRW_FacesYZ+1
        END IF
   END IF
END SUBROUTINE CountWeightFaceYZ


SUBROUTINE OutWeightFaceYZ(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Floor(ib)%Faces_YZ(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
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
END SUBROUTINE OutWeightFaceYZ

SUBROUTINE CountWeightFaceZX(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
        NrRW_FacesZX=NrRW_FacesZX+1
     END IF
   ELSE
       in_out=Vertices(i-1,j,k-1)%in_out &
             +Vertices(i,j,k-1)%in_out &
             +Vertices(i-1,j,k)%in_out &
             +Vertices(i,j,k)%in_out
       IF (in_out<=0) THEN
          NrRW_FacesZX=NrRW_FacesZX+1
       END IF
   END IF
END SUBROUTINE CountWeightFaceZX

SUBROUTINE OutWeightFaceZX(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out

   IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. &
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
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
END SUBROUTINE OutWeightFaceZX

SUBROUTINE CountWeightFaceXY(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       NrRW_FacesXY=NrRW_FacesXY+1
     END IF
   ELSE
      in_out=Vertices(i-1,j-1,k)%in_out &
            +Vertices(i,j-1,k)%in_out &
            +Vertices(i-1,j,k)%in_out &
            +Vertices(i,j,k)%in_out
      IF (in_out<=0) THEN
        NrRW_FacesXY=NrRW_FacesXY+1
      END IF
   END IF
END SUBROUTINE CountWeightFaceXY

SUBROUTINE OutWeightFaceXY(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
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
END SUBROUTINE OutWeightFaceXY

SUBROUTINE CountWeightFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       NrRW_FacesXY=NrRW_FacesXY+1
     END IF
   !ELSE
   !   in_out=Vertices(i-1,j-1,k)%in_out &
   !         +Vertices(i,j-1,k)%in_out &
   !         +Vertices(i-1,j,k)%in_out &
   !         +Vertices(i,j,k)%in_out
   !   IF (in_out<=0) THEN
   !     NrRW_FacesXY=NrRW_FacesXY+1
   !   END IF
   END IF
END SUBROUTINE CountWeightFaceXYBotton

SUBROUTINE OutWeightFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<4.AND.Face%ec>0).OR. &
         (Face%in_out<0.AND.Face%ec==-1).OR. & 
          Face%Vol==0.0d0) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Face%Vol
     END IF
   !ELSE
   !  in_out=Vertices(i-1,j-1,k)%in_out &
   !        +Vertices(i,j-1,k)%in_out &
   !        +Vertices(i-1,j,k)%in_out &
   !        +Vertices(i,j,k)%in_out
   !  IF (in_out<=0) THEN
   !    WRITE(10,*) i,j,k
   !    WRITE(10,*) 0.0d0
   !  END IF
   END IF
END SUBROUTINE OutWeightFaceXYBotton

SUBROUTINE CountWeightCell(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       NrRW_Cells=NrRW_Cells+1
     ELSE
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
       NrRW_Cells=NrRW_Cells+1
     END IF
   ELSE
        in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%in_out &
              +Floor(ib)%Vertices(i,j-1,k-1)%in_out &
              +Floor(ib)%Vertices(i-1,j-1,k)%in_out &
              +Floor(ib)%Vertices(i-1,j,k-1)%in_out &
              +Floor(ib)%Vertices(i-1,j,k)%in_out   &
              +Floor(ib)%Vertices(i,j-1,k)%in_out   &
              +Floor(ib)%Vertices(i,j,k-1)%in_out   &
              +Floor(ib)%Vertices(i,j,k)%in_out
        IF (in_out<=0) THEN
          NrRW_Cells=NrRW_Cells+1
        END IF
   END IF

END SUBROUTINE CountWeightCell

SUBROUTINE CountWeightCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       NrRW_Cells=NrRW_Cells+1
     ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
      ! NrRW_Cells=NrRW_Cells+1
     END IF
   ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
      !ELSE auch weg lassen
      !  in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k-1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k-1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k-1)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%in_out
      !  IF (in_out<=0) THEN
      !    !checke Cell z0+1, da Grenzfall 'ob'
      !    in_out=Floor(ib)%Vertices(i-1,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k+1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k+1)%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k+1)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k+1)%in_out
      !    IF (in_out<=0) THEN
      !      NrRW_Cells=NrRW_Cells+1
      !    END IF
      !  END IF
   END IF

END SUBROUTINE CountWeightCellBotton

SUBROUTINE OutWeightCell(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
         VolC(i,j,k)=dx(i)*dy(j)*dz(k) ! only for Emission and only Cartesian
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%Vol
       !WRITE(10,*) yes_v,Cell%Vol
       !WRITE(10,*) Cell(i,j,k)%Cell%Vol,Cell(i,j,k)%Cell%MidPoint
       VolC(i,j,k)=MAX(Cell%Vol,0.d0) ! only for Emission
     ELSE
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
       WRITE(10,*) i,j,k
       WRITE(10,*) 0.0d0
       !WRITE(10,*) no_v,0.0d0
       VolC(i,j,k)=0.0d0 ! only for Emission
     END IF
   ELSE
        in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%in_out &
              +Floor(ib)%Vertices(i,j-1,k-1)%in_out &
              +Floor(ib)%Vertices(i-1,j-1,k)%in_out &
              +Floor(ib)%Vertices(i-1,j,k-1)%in_out &
              +Floor(ib)%Vertices(i-1,j,k)%in_out   &
              +Floor(ib)%Vertices(i,j-1,k)%in_out   &
              +Floor(ib)%Vertices(i,j,k-1)%in_out   &
              +Floor(ib)%Vertices(i,j,k)%in_out
        IF (in_out<=0) THEN
          WRITE(10,*) i,j,k
          !WRITE(10,*) no_v,0.0d0
          WRITE(10,*) 0.0d0
          VolC(i,j,k)=0.0d0 ! only for Emission
        ELSE
          VolC(i,j,k)=dx(i)*dy(j)*dz(k) ! only for Emission and only Cartesian
        END IF
   END IF

END SUBROUTINE OutWeightCell

SUBROUTINE OutWeightCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out
INTEGER :: yes_v,no_v

   yes_v=1    ! yes_Volume
   no_v=0     ! no_volume
   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
         !Cellen mit Edge od. Point an Berg
         !Cellen mit MaxVol, nicht unter NrW_ ausgeben
         !maxVol je Face für RandCellen zu checken noch nicht möglich
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       WRITE(10,*) i,j,k-1
       WRITE(10,*) Cell%Vol
       !WRITE(10,*) yes_v,Cell%Vol
       !WRITE(10,*) Cell(i,j,k)%Cell%Vol,Cell(i,j,k)%Cell%MidPoint
     ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
        !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
        !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
        ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
      ! WRITE(10,*) i,j,k-1
      ! WRITE(10,*) 0.0d0
     END IF
   ELSE
      !wenn celle Vol=0.d0 ist, ist für Botton nicht von Interesse im ASAM 
      !ELSE auch weg lassen
      !  in_out=Floor(ib)%Vertices(i-1,j-1,k-1)%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k-1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k-1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k-1)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%in_out
      !  IF (in_out<=0) THEN
      !    !WRITE(10,*) i,j,k-1
      !    !WRITE(10,*) 0.0d0
      !    !checke Cell z0+1, da Grenzfall 'ob'
      !    in_out=Floor(ib)%Vertices(i-1,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i,j-1,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j-1,k+1)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k)%in_out &
      !        +Floor(ib)%Vertices(i-1,j,k+1)%in_out   &
      !        +Floor(ib)%Vertices(i,j-1,k+1)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k)%in_out   &
      !        +Floor(ib)%Vertices(i,j,k+1)%in_out
      !    IF (in_out<=0) THEN
      !       WRITE(10,*) i,j,k-1
      !       WRITE(10,*) Cell%Vol
      !    END IF
      !  END IF
   END IF

END SUBROUTINE OutWeightCellBotton

SUBROUTINE OutMPCellBorder(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out<7.AND.Cell%Vol>0.0d0.AND.Cell%vc>0 ) THEN
       !ausschalten Celle mit VolMax,Edge-Grenze,1Point Grenze
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%MidPoint
     ELSE IF (Cell%in_out==6.AND.Cell%vc>0) THEN
       !Cell%in_out==6 mit Schnittfläche
       WRITE(10,*) i,j,k
       WRITE(10,*) Cell%MidPoint
     END IF
   END IF

END SUBROUTINE OutMPCellBorder


SUBROUTINE WriteWeightBlk(FileName)
  CHARACTER*50 :: FileName
! -----------------------------
  INTEGER :: nl,k,j,i,izv,widthx,widthy,widthz
  INTEGER :: ib,in,posb
  INTEGER :: nx_ges,ny_ges,nz_ges
  INTEGER :: in_out
  INTEGER :: yes_v,no_v
!------------------------------------------------------------------------------
  zRauh=0.1d0
  !Albedo=0.5
  Albedo=0.08
  Emissivity=0.8
  yes_v=1
  no_v=0

  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Weight2',STATUS='unknown')
  CALL WriteAuswWeightToProt 

  IF((indata_type=='gk'.AND.conv_gk=='s').OR. &
     (indata_type=='gk'.AND.OutGrid=="Globe") .OR. &
      indata_type=='geo' .OR. indata_type=='rad' .OR. &
      indata_type=='wgs84') THEN
    WRITE(10,*) 'Globe'
  ELSE IF (indata_type=='cyl') THEN
    WRITE(10,*) 'Cyl' 
    !WRITE(10,*) OutGrid(1:3)
  ELSE
    If(OutGrid(1:4)/='Cart') THEN
       WRITE(*,*) "Output Grid stimmt nicht mit definition Cart überein!"
    END IF
    WRITE(10,*) 'Cart'
  END IF

  WRITE(10,*) nb,"  ! Number of blocks"  
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      IF (Nachbars(in)%nType(2:2)=="b") THEN
        posb=in
      END IF
    END DO
    Write(*,*) "                                     Block : ",ib,"\/",nb
    WRITE(10,*)  "#.... Block",ib,"....................................................."
    WRITE(10,*)  xP(ix0),yP(iy0),zP(iz0)
    CALL WriteBlkNrToProt(ib)

    WRITE(10,*) nx+2       !+rand 'w,e'
    WRITE(10,*) dx(ix0)
    DO i=ix0+1,ix1
       WRITE(10,*) dx(i)
    END DO
    WRITE(10,*) dx(ix1+1)

    WRITE(10,*) ny+2       !+rand 'n,s'
    WRITE(10,*) dy(iy0) 
    DO j=iy0+1,iy1
       WRITE(10,*) dy(j)
    END DO
    WRITE(10,*) dy(iy1+1)

    WRITE(10,*) nz+2       !+rand 't,b'
    WRITE(10,*) dz(iz0)
    DO k=iz0+1,iz1
       WRITE(10,*) dz(k)
    END DO
    WRITE(10,*) dz(iz1+1)

    !WRITE(10,*) nzLang     !!! Filaus
    !DO k=iz0+1,iz0+nzLang        !!! Filaus
    !  WRITE(10,*) domain%dz(k)   !!! Filaus
    !END DO
    !............................FacesYZ...............................................
    !Count akt. NrRW_FacesYZ
    !.......................
    ! count east,west-border
    widthx=((ix1+1)-(ix0-1))
    DO i=ix0-1,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL CountWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! count inside block
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL CountWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    IF(NrW_FacesYZ/=NrRW_FacesYZ) THEN
       Write(*,*) leerzei24,"(Warning  NrW_FacesYZ)"
       CALL WriteWarnNrFace(NrW_FacesYZ,NrRW_FacesYZ,"FacesYZ")
       NrW_FacesYZ=NrRW_FacesYZ
    ENDIF
    !Output Nr_FacesYZ
    !...................
    WRITE(10,*) NrW_FacesYZ, "   NrW_FacesYZ"
    ! east,west-border
    widthx=((ix1+1)-(ix0-1))
    DO i=ix0-1,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    !!1 ! north, south-Rand
    !!1 widthy=((iy1+1)-(iy0))
    !!1 DO k=iz0+1,iz1
    !!1    DO j=iy0,iy1+1,widthy  !border north,south
    !!1      DO i=ix0,ix1
    !!1       CALL OutWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    !!1 ! top, botton Rand
    !!1 widthz=((iz1+1)-(iz0))
    !!1 DO k=iz0,iz1+1,widthz    !border top,botton 
    !!1    DO j=iy0+1,iy1
    !!1     DO i=ix0,ix1
    !!1       CALL OutWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    !block
    !WRITE(10,*) NrW_FacesYZ, "   NrW_FacesYZ"   !WRITE(10,*) NrAll_FacesYZ
    ! block
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
          !IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
          !  IF ((Faces_YZ(i,j,k)%Face%in_out<4.AND.Faces_YZ(i,j,k)%Face%ec>0).OR. &
          !      (Faces_YZ(i,j,k)%Face%in_out<0.AND.Faces_YZ(i,j,k)%Face%ec==-1).OR. &
          !       Faces_YZ(i,j,k)%Face%Vol==0.0d0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) Faces_YZ(i,j,k)%Face%Vol
          !  END IF
          !ELSE
          !  in_out=Vertices(i,j-1,k-1)%in_out &
          !        +Vertices(i,j,k-1)%in_out &
          !        +Vertices(i,j-1,k)%in_out &
          !        +Vertices(i,j,k)%in_out
          !  IF (in_out<=0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) 0.0d0
          !  END IF
          !END IF
        END DO
      END DO
    END DO
    WRITE(10,*) NrMP_FacesYZ, "   NrMP_FacesYZ"
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
    !............................FacesZX...............................................
    !Count akt. NrRW_FacesZX
    !........................
    ! count north,south-border
    widthy=((iy1+1)-(iy0-1))
    DO j=iy0-1,iy1+1,widthy
      DO i=ix0+1,ix1
         DO k=iz0+1,iz1
          CALL CountWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! count inside block
    DO j=iy0,iy1
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL CountWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    IF(NrW_FacesZX/=NrRW_FacesZX) THEN
       Write(*,*) leerzei24,"(Warning  NrW_FacesZX)"
       CALL WriteWarnNrFace(NrW_FacesZX,NrRW_FacesZX,"FacesZX")
       NrW_FacesZX=NrRW_FacesZX
    ENDIF
    !Output all to  FacesZX
    !......................
    WRITE(10,*) NrW_FacesZX, "   NrW_FacesZX"
    !!1 ! east,west-border
    !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2 ok.
    !!1 widthx=((ix1+1)-(ix0))
    !!1 DO k=iz0+1,iz1
    !!1    DO j=iy0,iy1
    !!1      DO i=ix0,ix1+1,widthx   !border east,west
    !!1       CALL OutWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    ! north, south-Rand
    widthy=((iy1+1)-(iy0-1))
    DO j=iy0-1,iy1+1,widthy  !border north,south
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL OutWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    !!1 ! top, botton Rand
    !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2 ok.
    !!1 widthz=((iz1+1)-(iz0))
    !!1 DO k=iz0,iz1+1,widthz    !border top,botton 
    !!1    DO j=iy0,iy1
    !!1     DO i=ix0+1,ix1
    !!1       CALL OutWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    !block
    !WRITE(10,*) NrW_FacesZX, "   NrW_FacesZX"   !WRITE(10,*) NrAll_FacesZX
    ! block
     DO j=iy0,iy1
       DO i=ix0+1,ix1
         DO k=iz0+1,iz1
          CALL OutWeightFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
          !IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
          !  IF ((Faces_ZX(i,j,k)%Face%in_out<4.AND.Faces_ZX(i,j,k)%Face%ec>0).OR. &
          !      (Faces_ZX(i,j,k)%Face%in_out<0.AND.Faces_ZX(i,j,k)%Face%ec==-1).OR. &
          !       Faces_ZX(i,j,k)%Face%Vol==0.0d0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) Faces_ZX(i,j,k)%Face%Vol
          !  END IF
          !ELSE
          !  in_out=Vertices(i-1,j,k-1)%in_out &
          !        +Vertices(i,j,k-1)%in_out &
          !        +Vertices(i-1,j,k)%in_out &
          !        +Vertices(i,j,k)%in_out
          !  IF (in_out<=0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) 0.0d0
          !  END IF
          !END IF
        END DO
      END DO
    END DO
    WRITE(10,*) NrMP_FacesZX, "   NrMP_FacesZX"
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
    !............................FacesXY...............................................
    !Count akt. NrRW_FacesXY
    !.......................
    ! count border top,botton
    widthz=((iz1+1)-(iz0-1))
    !DO k=iz0-1,iz1+1,widthz
    !   DO j=iy0+1,iy1
    !    DO i=ix0+1,ix1
    !      CALL CountWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
    !    END DO
    !  END DO
    !END DO
    IF (Nachbars(posb)%nType=='ob') THEN
      ! ob- Berg-Grenze unter zg0 oder Standard
      k=iz0
      DO j=iy0+1,iy1
       DO i=ix0+1,ix1
         in_out=Vertices(i-1,j-1,k)%in_out &
               +Vertices(i  ,j-1,k)%in_out &
               +Vertices(i-1,j  ,k)%in_out &
               +Vertices(i  ,j  ,k)%in_out
           !Berg-Grenze  unterhalb z0
           IF (in_out==4) THEN
           !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN
           !IF (in_out>0) THEN !prot_lnok2
!27.10.11              NrRW_FacesXY=NrRW_FacesXY+1   ! iz0  zählend
!27.10.11              NrRW_FacesXY=NrRW_FacesXY+1   ! iz0-1 zählend
           ELSE 
             CALL CountWeightFaceXYBotton(Faces_XY(i,j,k)%Face,i,j,k-1,ib)
           END IF
        END DO
      END DO
     ELSE  ! 'ib'
       k=iz0-1
         DO j=iy0+1,iy1
           DO i=ix0+1,ix1
             CALL CountWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
           END DO
         END DO
    END IF 
    ! it/ot 
    k=iz1+1
      DO j=iy0+1,iy1
       DO i=ix0+1,ix1
          CALL CountWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    ! count inside block
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    IF(NrW_FacesXY/=NrRW_FacesXY) THEN
       Write(*,*) leerzei24,"(Warning  NrW_FaceXY)"
       CALL WriteWarnNrFace(NrW_FacesXY,NrRW_FacesXY,"FacesXY")
       NrW_FacesXY=NrRW_FacesXY
    ELSE 
       NrW_FacesXY=NrRW_FacesXY ! NrRN_FXY wird im ob Rand mitgezählt aber nicht ausgegeben!
    ENDIF
    !Output all to FacesXY
    !......................
    WRITE(10,*) NrW_FacesXY, "   NrW_FacesXY"
    !!1 ! east,west-border
    !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2 ok.
    !!1 widthx=((ix1+1)-(ix0))
    !!1 DO k=iz0,iz1
    !!1    DO j=iy0+1,iy1
    !!1      DO i=ix0,ix1+1,widthx   !border east,west
    !!1       CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    !!1 ! north, south-Rand
    !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2 ok.
    !!1 widthy=((iy1+1)-(iy0))
    !!1 DO k=iz0,iz1
    !!1    DO j=iy0,iy1+1,widthy  !border north,south
    !!1      DO i=ix0+1,ix1
    !!1       CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
    !!1     END DO
    !!1   END DO
    !!1 END DO
    ! Output: top, botton Rand
    widthz=((iz1+1)-(iz0-1))
    !DO k=iz0-1,iz1+1,widthz
    !   DO j=iy0+1,iy1
    !    DO i=ix0+1,ix1
    !      CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
    !    END DO
    !  END DO
    !END DO
    ! it/ot
    k=iz1+1
    DO j=iy0+1,iy1
      DO i=ix0+1,ix1
         CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
      END DO
    END DO
    ! Output: 'ob'-Rand 
    ! Speziell: Berg-Grenze unter z0
    ! Faces_XY: z0-1
    IF (Nachbars(posb)%nType=='ob') THEN
       k=iz0
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
           in_out=Vertices(i-1,j-1,k)%in_out &
                 +Vertices(i  ,j-1,k)%in_out &
                 +Vertices(i-1,j  ,k)%in_out &
                 +Vertices(i  ,j  ,k)%in_out
           IF (in_out==4) THEN    !Berg-Grenze unter zg0
           !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN    !Berg-Grenze unter zg0
           !IF (in_out>0) THEN    !Berg-Grenze unter zg0
!27.10.11             WRITE(10,*) i,j,k-1  !Rand iz0-1
!27.10.11             WRITE(10,*) 0.0d0
           ELSE
             CALL OutWeightFaceXYBotton(Faces_XY(i,j,k)%Face,i,j,k-1,ib)
           END IF
         END DO
       END DO
    ELSE  ! 'ib'
       k=iz0-1
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
            CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
         END DO
       END DO
    END IF
    !Output block: z0-Faces_XY
    !Speziell !Berg-Grenze unter zg0
    IF (Nachbars(posb)%nType=='ob') THEN
       k=iz0
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
           in_out=Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j,k)%in_out
           IF (in_out==4) THEN     !Berg-Grenze unter zg0
           !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN     !Berg-Grenze unter zg0
           !IF (in_out>0) THEN     !Berg-Grenze unter zg0
!27.10.11             WRITE(10,*) i,j,k     !Block-z0-faces_xy
!27.10.11             WRITE(10,*) 0.0d0
           END IF
         END DO
       END DO
    END IF

    ! Output: block
    !WRITE(10,*) NrW_FacesXY, "   NrW_FacesXY"    !WRITE(10,*) NrAll_FacesXY
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL OutWeightFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
          !IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
          !  IF ((Faces_XY(i,j,k)%Face%in_out<4.AND.Faces_XY(i,j,k)%Face%ec>0).OR. &
          !      (Faces_XY(i,j,k)%Face%in_out<0.AND.Faces_XY(i,j,k)%Face%ec==-1).OR. & 
          !       Faces_XY(i,j,k)%Face%Vol==0.0d0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) Faces_XY(i,j,k)%Face%Vol
          !  END IF
          !ELSE
          !  in_out=Vertices(i-1,j-1,k)%in_out &
          !        +Vertices(i,j-1,k)%in_out &
          !        +Vertices(i-1,j,k)%in_out &
          !        +Vertices(i,j,k)%in_out
          !  IF (in_out<=0) THEN
          !    WRITE(10,*) i,j,k
          !    WRITE(10,*) 0.0d0
          !  END IF
          !END IF
        END DO
      END DO
    END DO
    WRITE(10,*) NrMP_FacesXY, "   NrMP_FacesXY"
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
    !.............................Cells...............................................
    widthx=((ix1+1)-(ix0))
    widthy=((iy1+1)-(iy0))
    widthz=((iz1+1)-(iz0))
    !Count akt. NrRW_Cells
    !----------------------
    ! east,west-border
    DO i=ix0,ix1+1,widthx
       DO j=iy0+1,iy1
         DO k=iz0+1,iz1
          CALL CountWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! north, south-border
    DO j=iy0,iy1+1,widthy
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL CountWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! top border
    k=iz1+1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    ! botton border
    IF (Nachbars(posb)%nType=='ob') THEN
       k=iz0
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
           !in_out=Vertices(i-1,j-1,k)%in_out &
           !       +Vertices(i,j-1,k)%in_out &
           !       +Vertices(i-1,j,k)%in_out &
           !       +Vertices(i,j,k)%in_out
           !Berg-Grenze unterhalb z0
           !IF (in_out==4) THEN
           ! Cells: z0
            CALL CountWeightCellBotton(Cell(i,j,k)%Cell,i,j,k,ib)
!27.10.11             NrRW_Cells=NrRW_Cells+1
           !END IF
         END DO
       END DO
    ELSE !'ib' 
       k=iz0
       DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL CountWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
          END DO
       END DO
    END IF
    ! block 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    IF(NrW_Cells/=NrRW_Cells) THEN
       Write(*,*) leerzei24,"(Warning  NrW_Cells)"
       CALL WriteWarnNrCells(NrW_Cells,NrRW_Cells)
       NrW_Cells=NrRW_Cells
    ENDIF
    ! Output Cells
    !-------------
    WRITE(10,*) NrW_Cells, "    NrW_Cells"
    ! east,west-border
    DO i=ix0,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! north, south border
    DO j=iy0,iy1+1,widthy
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL OutWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! top border 
    !DO k=iz0,iz1+1,widthz
    k=iz1+1
    DO j=iy0+1,iy1
      DO i=ix0+1,ix1
        CALL OutWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
      END DO
    END DO
    ! botton border
    IF (Nachbars(posb)%nType=='ob') THEN
       k=iz0
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
          CALL OutWeightCellBotton(Cell(i,j,k)%Cell,i,j,k,ib)
!27.10.11           WRITE(10,*) i,j,0
!27.10.11           WRITE(10,*) 0.0d0
         END DO
       END DO
     ELSE !'ib' 
       k=iz0
       DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL OutWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
          END DO
       END DO
    END IF 
    ! block 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL OutWeightCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO

!    WRITE(10,*) NrMP_Cells, "    NrMP_Cells"
    !!! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
    !!! east,west-border
    !!DO k=iz0+1,iz1
    !!  DO j=iy0+1,iy1
    !!    DO i=ix0,ix1+1,widthx   !border east,west
    !!      CALL OutMPCellBorder(Cell(i,j,k)%Cell,i,j,k,ib)
    !!    END DO
    !!  END DO
    !!END DO
    !!! north,south-border
    !!DO k=iz0+1,iz1
    !!  DO j=iy0,iy1+1,widthy  !border north,south
    !!    DO i=ix0+1,ix1
    !!      CALL OutMPCellBorder(Cell(i,j,k)%Cell,i,j,k,ib)
    !!    END DO
    !!  END DO
    !!END DO
    !!! top,botton border 
    !!DO k=iz0,iz1+1,widthz    !border top,botton 
    !!  DO j=iy0+1,iy1
    !!    DO i=ix0+1,ix1
    !!      CALL OutMPCellBorder(Cell(i,j,k)%Cell,i,j,k,ib)
    !!    END DO
    !!  END DO
    !!END DO
    ! Block 
!    DO k=iz0+1,iz1
!      DO j=iy0+1,iy1
!        DO i=ix0+1,ix1
!          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!            IF (Cell(i,j,k)%Cell%in_out<8.AND.Cell(i,j,k)%Cell%Vol>0.0d0.AND. &
!                  (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1)<1.d0-dist_scMaxCell &
!               .OR.Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2)<1.d0-dist_scMaxCell  &
!               .OR.Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3)<1.d0-dist_scMaxCell  &
!               .OR.Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4)<1.d0-dist_scMaxCell  &
!               .OR.Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5)<1.d0-dist_scMaxCell  &
!               .OR.Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6)<1.d0-dist_scMaxCell)) THEN
!              WRITE(10,*) i,j,k
!              WRITE(10,*) Cell(i,j,k)%Cell%MidPoint
!            END IF
!          END IF
!        END DO
!      END DO
!    END DO
    IF( NrMP_Cells/=NrB_Cells) THEN
       Write(*,*)  "ib=",ib,"  NrMP_Cells=", NrMP_Cells,"  NrB_Cells=",NrB_Cells 
    END IF
    IF( NrMP_Cells/=NrB_Cells) STOP  'NrMP_Cells != NrB_Cells'
    WRITE(10,*) NrB_Cells, "      NrB_Cells", struct_bound, "   Bound_Struct"
    ! Ausgabe: Reihenfolge Schleifen gleich Bound-Schleife Output-ASAM anpassen
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
        IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            !Orig IF (Cell(i,j,k)%Cell%in_out<8.AND.Cell(i,j,k)%Cell%Vol>0.0d0.AND. &
            !Special hinzu !!
            IF (((Cell(i,j,k)%Cell%in_out<5.AND.Cell(i,j,k)%Cell%Vol>0.0d0) .OR. &
                 (Cell(i,j,k)%Cell%in_out==6.AND.Cell(i,j,k)%Cell%vc>0) .OR. &
                 (Cell(i,j,k)%Cell%in_out==5.AND.Cell(i,j,k)%Cell%vc>0)) &
                .AND. &
                (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1) &
                   & <1.d0-dist_scMaxCell     .OR. &
                 Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2) &
                   & <1.d0-dist_scMaxCell     .OR. &
                 Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3) &
                   & <1.d0-dist_scMaxCell     .OR. &
                 Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4) &
                   & <1.d0-dist_scMaxCell     .OR. &
                 Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5) &
                   & <1.d0-dist_scMaxCell     .OR. &
                 Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6) &
                   & <1.d0-dist_scMaxCell) &
               ) THEN
              WRITE(10,*) i,j,k
              WRITE(10,*) Cell(i,j,k)%Cell%MidPoint
              !IF (xP(i-1)+0.5d0*dx(i)<=-2.0d4) THEN
!             !  WRITE(10,'(6d15.7)') 1.0d-99,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              !  WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              !ELSE IF (xP(i-1)+0.5d0*dx(i)<=-1.0d4) THEN
!             !  WRITE(10,'(6d15.7)') (xP(i-1)+0.5d0*dx(i)+2.0d4)/1.0d4*zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              !  WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              !ELSE
              !  WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,0.5d0,0.5d0,0.0d0
              !END IF
              WRITE(10,'(6d15.7)') zRauh,Albedo,Emissivity,Cell(i,j,k)%Cell%CutF_MidP
              IF(nr_soildef>0) THEN
                DO nl=1,nr_soildef
                   IF(i>(s_ixa(nl)*2.e0**RefineX).AND.i<=(s_ixe(nl)*2.e0**RefineX) .AND. &
                      j>(s_iya(nl)*2.e0**RefineY).AND.j<=(s_iye(nl)*2.e0**RefineY)) THEN
                      WRITE(10,'(I3,I3,I3,I3,I3,I3,I3,I3)') Cell(i,j,k)%Cell%LandClass,soil_type(nl,1:7)
                      !WRITE(10,'(I3,I3,I3,I3,I3,I3,I3,I3)') soil_type(nl,1:7)
                      !Write(10,*) "ib=", ib, "s_ixa(nl)=", s_ixa(nl), " s_ixe(nl)=",s_ixe(nl)
                      !WRITE(10,*) "RefineX=", RefineX, " RefineY=",RefineY, " RefineZ=",RefineZ 
                   END IF
                END DO
              ELSE IF (nr_landdef>0 .AND. nr_soildef==0) THEN
                WRITE(10,'(I3)') Cell(i,j,k)%Cell%LandClass
              END IF
            END IF
        END IF    ! IF (ASSOCIATED....
        END DO
      END DO
    END DO

  END DO  !...nb 
  CALL WriteEndWeightToProt
  CLOSE(10)

END SUBROUTINE WriteWeightBlk


END MODULE OutputWeightBlk_Mod
!##########################################################################################
