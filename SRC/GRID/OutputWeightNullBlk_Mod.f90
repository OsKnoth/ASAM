!###############################################################################
! File: OutputWeightNullBlk_Mod.f90
!  - Ausgabe Weight-Null Face Cellen zum Input ASAM-Model 
!
!------------------------------------------------------------------------------
! Description:  -Ausgabe dx,dy,dz;
!                        Faces,Cellen unterhalb Berg mit Inhalt Volumen=0
!                     -> Block/Blöcke
!                     -> In-/Out Border east,-west,-north,-south,-top 
!                     -> Border 'ob', 'out bottom' nicht
!               -Ausgabe Numbers Weight:
!                     -> NrRW_FacesYZ, NrRW_FacesZX, NrRW_FacesXY
!                     -> NrRW_Cells 
!               -Ausgabe Numbers MittelPoint, Bound, Bound_Struct: auf '0' gesetzt
!                     -> NrMP_FacesYZ, NrMP_FacesZX, NrMP_FacesXY
!                     -> NrB_Cells, struct_bound
!------------------------------------------------------------------------------
!################################################################################

MODULE OutputWeightNullBlk_Mod

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

SUBROUTINE CountWeightNullFaceYZ(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Floor(ib)%Faces_YZ(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE CountWeightNullFaceYZ


SUBROUTINE OutWeightNullFaceYZ(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Floor(ib)%Faces_YZ(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE OutWeightNullFaceYZ

SUBROUTINE CountWeightNullFaceZX(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE CountWeightNullFaceZX

SUBROUTINE OutWeightNullFaceZX(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out

   IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE OutWeightNullFaceZX

SUBROUTINE CountWeightNullFaceXY(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE CountWeightNullFaceXY

SUBROUTINE OutWeightNullFaceXY(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE OutWeightNullFaceXY

SUBROUTINE CountWeightNullFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
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
END SUBROUTINE CountWeightNullFaceXYBotton

SUBROUTINE OutWeightNullFaceXYBotton(Face,i,j,k,ib)
TYPE (Face_T) :: Face
INTEGER       :: i,j,k,ib

INTEGER :: in_out
   IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
     IF ((Face%in_out<0.AND.Face%ec==-1).OR. Face%Vol==0.0d0) THEN
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
END SUBROUTINE OutWeightNullFaceXYBotton

SUBROUTINE CountWeightNullCell(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out

   ! RandCellen,
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
       ! max raus
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       !NrRW_Cells=NrRW_Cells+1
     ELSE
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

END SUBROUTINE CountWeightNullCell

SUBROUTINE OutWeightNullCell(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out

   ! RandCellen,
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
       ! max raus
     ELSE IF (Cell%in_out<6.AND.Cell%Vol>0.0d0.OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       !WRITE(10,*) i,j,k
       !WRITE(10,*) Cell%Vol
     ELSE
       WRITE(10,*) i,j,k
       WRITE(10,*) 0.0d0
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
          WRITE(10,*) 0.0d0
        END IF
   END IF

END SUBROUTINE OutWeightNullCell

SUBROUTINE CountWeightNullCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out

   ! RandCellen,
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
       ! max raus
     ELSE IF ((Cell%in_out<6.AND.Cell%Vol>0.0d0.AND.Cell%vc>0).OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       !NrRW_Cells=NrRW_Cells+1
     ELSE
      !wenn celle Vol=0.d0 ist, o-Botton ist, für ASAM nicht von Interesse
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

END SUBROUTINE CountWeightNullCellBotton

SUBROUTINE OutWeightNullCellBotton(Cell,i,j,k,ib)
TYPE(Cell_T) :: Cell
INTEGER :: i,j,k,ib

INTEGER :: in_out

   ! RandCellen, alternative da Rand-Faces noch nicht bearbeitet
   IF (ASSOCIATED(Floor(ib)%Cell(i,j,k)%Cell)) THEN
     IF (Cell%in_out>6.OR.(Cell%in_out==6.AND.Cell%vc==0)) THEN
       ! max raus
     ELSE IF ((Cell%in_out<6.AND.Cell%Vol>0.0d0.AND.Cell%vc>0).OR. &
              (Cell%in_out==6.AND.Cell%vc>0) ) THEN
       !WRITE(10,*) i,j,k-1
       !WRITE(10,*) Cell%Vol
     ELSE
      !wenn celle Vol=0.d0 ist, o-Botton ist, für ASAM nicht von Interesse
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

END SUBROUTINE OutWeightNullCellBotton




SUBROUTINE WriteWeightNullBlk(FileName)
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

  OPEN(UNIT=10,FILE=TRIM(FileName)//'.WNull',STATUS='unknown')
  CALL WriteAuswWeightNullToProt 

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

    !............................FacesYZ...............................................
    !Count akt. NrRW_FacesYZ
    !.......................
    ! count east,west-border
    widthx=((ix1+1)-(ix0-1))
    DO i=ix0-1,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL CountWeightNullFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! count inside block
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL CountWeightNullFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    !Output Nr_FacesYZ
    !...................
    WRITE(10,*) NrRW_FacesYZ, "   NrW_FacesYZ"
    ! east,west-border
    widthx=((ix1+1)-(ix0-1))
    DO i=ix0-1,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightNullFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! block
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightNullFaceYZ(Faces_YZ(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    WRITE(10,*) 0 , "   NrMP_FacesYZ"
    !............................FacesZX...............................................
    !Count akt. NrRW_FacesZX
    !........................
    ! count north,south-border
    widthy=((iy1+1)-(iy0-1))
    DO j=iy0-1,iy1+1,widthy
      DO i=ix0+1,ix1
         DO k=iz0+1,iz1
          CALL CountWeightNullFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! count inside block
    DO j=iy0,iy1
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL CountWeightNullFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    !Output all to  FacesZX
    !......................
    WRITE(10,*) NrRW_FacesZX, "   NrW_FacesZX"
    ! north, south-Rand
    widthy=((iy1+1)-(iy0-1))
    DO j=iy0-1,iy1+1,widthy  !border north,south
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL OutWeightNullFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    ! block
    DO j=iy0,iy1
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
         CALL OutWeightNullFaceZX(Faces_ZX(i,j,k)%Face,i,j,k,ib)
       END DO
     END DO
    END DO
    WRITE(10,*) 0, "   NrMP_FacesZX"
    !............................FacesXY...............................................
    !Count akt. NrRW_FacesXY
    !.......................
    ! count border top,botton
    widthz=((iz1+1)-(iz0-1))
    IF (Nachbars(posb)%nType=='ob') THEN
      ! ob- Berg-Grenze unter zg0 oder Standard
     !k=iz0
     !DO j=iy0+1,iy1
     ! DO i=ix0+1,ix1
     !   in_out=Vertices(i-1,j-1,k)%in_out &
     !        +Vertices(i  ,j-1,k)%in_out &
     !         +Vertices(i-1,j  ,k)%in_out &
     !         +Vertices(i  ,j  ,k)%in_out
     !     !Berg-Grenze  unterhalb z0
     !     IF (in_out==4) THEN
     !     !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN
     !     !IF (in_out>0) THEN
     !         !NrRW_FacesXY=NrRW_FacesXY+1   ! iz0  zählend
     !         !NrRW_FacesXY=NrRW_FacesXY+1   ! iz0-1 zählend
     !     ELSE 
     !       !CALL CountWeightNullFaceXYBotton(Faces_XY(i,j,k)%Face,i,j,k-1,ib)
     !     END IF
     !  END DO
     !END DO
     ELSE  ! 'ib'
       k=iz0-1
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
            CALL CountWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
         END DO
       END DO
    END IF 
    ! it/ot 
    k=iz1+1
    DO j=iy0+1,iy1
      DO i=ix0+1,ix1
         CALL CountWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
      END DO
    END DO
    ! count inside block
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    !Output all to FacesXY
    !......................
    WRITE(10,*) NrRW_FacesXY, "   NrW_FacesXY"
    ! Output: top, botton Rand
    widthz=((iz1+1)-(iz0-1))
    ! it/ot
    k=iz1+1
    DO j=iy0+1,iy1
      DO i=ix0+1,ix1
         CALL OutWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
      END DO
    END DO
    ! Output: 'ob'-Rand 
    ! Speziell: Berg-Grenze unter z0
    ! Faces_XY: z0-1
    IF (Nachbars(posb)%nType=='ob') THEN
      !k=iz0    ! 'iz0' zum Check für Zuweisung 'iz0-1' 
      !DO j=iy0+1,iy1
      !  DO i=ix0+1,ix1
      !    in_out=Vertices(i-1,j-1,k)%in_out &
      !          +Vertices(i  ,j-1,k)%in_out &
      !          +Vertices(i-1,j  ,k)%in_out &
      !          +Vertices(i  ,j  ,k)%in_out
      !    IF (in_out==4) THEN
      !    !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN
      !    !IF (in_out>0) THEN    !Berg-Grenze unter zg0
      !       !WRITE(10,*) i,j,k-1  !Rand iz0-1
      !       !WRITE(10,*) 0.0d0
      !    ELSE
      !      !CALL OutWeightNullFaceXYBotton(Faces_XY(i,j,k)%Face,i,j,k-1,ib)
      !    END IF
      !  END DO
      !END DO
    ELSE  ! 'ib'
       k=iz0-1
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
            CALL OutWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
         END DO
       END DO
    END IF
    !Output block: z0-Faces_XY
    !Speziell: Berg-Grenze unter zg0
    IF (Nachbars(posb)%nType=='ob') THEN
       k=iz0
       DO j=iy0+1,iy1
         DO i=ix0+1,ix1
           in_out=Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j,k)%in_out
           IF (in_out==4) THEN
           !IF (in_out>0.AND.(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face))) THEN     
           !IF (in_out>0) THEN 
             !WRITE(10,*) i,j,k     !Block-z0-faces_xy
             !WRITE(10,*) 0.0d0
           END IF
         END DO
       END DO
    END IF

    ! Output: block
    DO k=iz0,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL OutWeightNullFaceXY(Faces_XY(i,j,k)%Face,i,j,k,ib)
        END DO
      END DO
    END DO
    WRITE(10,*) 0 , "   NrMP_FacesXY"

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
          CALL CountWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! north, south-border
    DO j=iy0,iy1+1,widthy
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL CountWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! top border
    k=iz1+1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    ! botton border
    IF (Nachbars(posb)%nType=='ob') THEN
      !k=iz0
      !DO j=iy0+1,iy1
      !  DO i=ix0+1,ix1
      !    ! CALL CountWeightNullCellBotton(Cell(i,j,k)%Cell,i,j,k,ib)
      !  END DO
      !END DO
    ELSE !'ib' 
       k=iz0
       DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL CountWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
          END DO
       END DO
    END IF
    ! block 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL CountWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! Output Cells
    !-------------
    WRITE(10,*) NrRW_Cells, "    NrW_Cells"
    ! east,west-border
    DO i=ix0,ix1+1,widthx
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          CALL OutWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! north, south border
    DO j=iy0,iy1+1,widthy
      DO i=ix0+1,ix1
        DO k=iz0+1,iz1
          CALL OutWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO
    ! top border 
    k=iz1+1
    DO j=iy0+1,iy1
      DO i=ix0+1,ix1
        CALL OutWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
      END DO
    END DO
    ! botton border
    IF (Nachbars(posb)%nType=='ob') THEN
      !k=iz0
      !DO j=iy0+1,iy1
      !  DO i=ix0+1,ix1
      !   !CALL OutWeightNullCellBotton(Cell(i,j,k)%Cell,i,j,k,ib)
      !  END DO
      !END DO
     ELSE !'ib' 
       k=iz0
       DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL OutWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
          END DO
       END DO
    END IF 
    ! block 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          CALL OutWeightNullCell(Cell(i,j,k)%Cell,i,j,k,ib)
        END DO
      END DO
    END DO

    WRITE(10,*) 0 , "      NrB_Cells", 0, "   Bound_Struct"

  END DO  !...nb 
  CALL WriteEndWeightNullToProt
  CLOSE(10)

END SUBROUTINE WriteWeightNullBlk


END MODULE OutputWeightNullBlk_Mod
!##########################################################################################
