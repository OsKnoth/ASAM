SUBROUTINE AnalyzeAllVertices1

  INTEGER :: ix,iy,iz,jx,jy,jz
  INTEGER :: ib,in
  INTEGER :: withx,withy,withz
  INTEGER :: ListP(3,5000),NumList,iNum,in_out
  ! ListP(3,1000) 1000 -> Allocatable:: (((ix1+1-ix0)*(iy1+1-iy0)*(iz1+1-iz0))*2)*3)
  ! DresdenNeustadt  300*200*5 Gitter als Bsp:

  WRITE(OutUnitChL,*) "---AnalyzeAllVertices-------------------------------------"
  dxLoc=(Domain%x1-Domain%x0)/Domain%nx*distx_coeff
  dyLoc=(Domain%y1-Domain%y0)/Domain%ny*disty_coeff
  dzLoc=(Domain%z1-Domain%z0)/Domain%nz*distz_coeff
  !! dist(x,y,z)_coeff=0.05d0  ! bei Surfit zu gross
  !! dist(x,y,z)_coeff=0.01d0   ! bei Surfit okey

  calc_z=0.0
  ncz=0
  nr_out=0
  nr_cutplane=0
  nr_inside=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    !Write(*,*) "Block :  ", ib
    DO ix=ix0,ix1
      DO iy=iy0,iy1
        DO iz=iz0,iz1
           IF (Vertices(ix,iy,iz)%nrP==-1) THEN
              CALL CheckVertex(Vertices(ix,iy,iz))
!           ELSE
!             CALL CheckVertexInOut(Vertices(ix,iy,iz))
           END IF
        END DO
      END DO
    END DO
!OSSI   CALL CheckVertRandBlock 

    !WRITE(*,*) nr_out
    NumList=0
!   Reduce vertices
!   XY-Faces
    Write (OutUnitChL,*)"ib, Reduce vertices XY-Faces", ib
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
           IF (Vertices(ix,iy,iz)%in_out==-1 .AND. &
               Vertices(ix-1,iy-1,iz)%in_out==-1 .AND. &
               Vertices(ix-1,iy,iz)%in_out>=0 .AND. &
               Vertices(ix,iy-1,iz)%in_out>=0.) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix-1
               ListP(2,NumList)=iy-1
               ListP(3,NumList)=iz
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz
             ELSE IF (Vertices(ix-1,iy,iz)%in_out==-1 .AND. &
                 Vertices(ix,iy-1,iz)%in_out==-1 .AND. &
               Vertices(ix,iy,iz)%in_out>=0 .AND. &
               Vertices(ix-1,iy-1,iz)%in_out>=0.) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix-1
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy-1
               ListP(3,NumList)=iz
           END IF
        END DO
      END DO
    END DO
    WRITE(OutUnitChL,*) "NumList=",NumList
!   XZ-Faces
    Write (OutUnitChL,*)"ib ,Reduce vertices XZ-Faces", ib
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
             IF (Vertices(ix,iy,iz)%in_out==-1 .AND. &
                 Vertices(ix-1,iy,iz-1)%in_out==-1 .AND. &
                 Vertices(ix-1,iy,iz)%in_out>=0 .AND. &
                 Vertices(ix,iy,iz-1)%in_out>=0 ) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix-1
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz-1
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz
             ELSE IF (Vertices(ix-1,iy,iz)%in_out==-1.AND. &
                 Vertices(ix,iy,iz-1)%in_out==-1 .AND. &
                 Vertices(ix-1,iy,iz-1)%in_out>=0 .AND. &
                 Vertices(ix,iy,iz)%in_out>=0 ) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix-1
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz-1
             END IF
        END DO
      END DO
    END DO
    WRITE(OutUnitChL,*) "NumList=",NumList
!   YZ-Faces
    Write (OutUnitChL,*)"ib, Reduce vertices YZ-Faces", ib
    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
             IF (Vertices(ix,iy,iz)%in_out==-1 .AND. &
                 Vertices(ix,iy-1,iz-1)%in_out==-1 .AND. &
                 Vertices(ix,iy,iz-1)%in_out>=0 .AND. &
                 Vertices(ix,iy-1,iz)%in_out>=0) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy-1
               ListP(3,NumList)=iz-1
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz
             ELSE IF (Vertices(ix,iy-1,iz)%in_out==-1 .AND. &
                 Vertices(ix,iy,iz-1)%in_out==-1 .AND. &
                 Vertices(ix,iy-1,iz-1)%in_out>=0 .AND. &
                 Vertices(ix,iy,iz)%in_out>=0) THEN
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy-1
               ListP(3,NumList)=iz
               NumList=NumList+1
               ListP(1,NumList)=ix
               ListP(2,NumList)=iy
               ListP(3,NumList)=iz-1
             END IF
        END DO
      END DO
    END DO
    WRITE(OutUnitChL,*) "NumList=",NumList
    DO iNum=1,NumList
      ix=ListP(1,iNum)
      iy=ListP(2,iNum)
      iz=ListP(3,iNum)
      Vertices(ix,iy,iz)%in_out=0
      !Vertices(ix,iy,iz)%in_out=1
      Vertices(ix,iy,iz)%Shift=0
      IF (Vertices(ix,iy,iz)%nrP<=0) THEN
        nr_out=nr_out+1
        Vertices(ix,iy,iz)%nrP=nr_out
      END IF
      IF (Vertices(ix,iy,iz)%nrCutP<=0) THEN
        nr_cutplane=nr_cutplane+1
        Vertices(ix,iy,iz)%nrCutP=nr_cutplane
      END IF
      IF (Vertices(ix,iy,iz)%nrInP<=0) THEN
        nr_inside=nr_inside+1
        Vertices(ix,iy,iz)%nrInP=nr_inside
      END IF
    END DO

 calc_z=0.0
 ncz=0
 !OSSI CALL SearchMaxCalcZ()
 Write(OutUnitChL,*) "ib, Vergebe Nummern an RandBloecke", ib
!   Vergebe Nummern in Block ib

    DO in=1,AnzahlNachbar

      CALL Set(Nachbars(in))
      Write(OutUnitChL,*) "  ib -> ", ib,",",in, " <- nachbar"

!     Kopiere Nummern in gleich Vertices (Nodes)
      !IF (Nachbars(in)%nType(1:1)=='i'.OR.Nachbars(in)%nType(1:1)=='p') THEN

      SELECT CASE (Nachbars(in)%nType(2:2))
      CASE ("w","e")
        !IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw'.OR. &
        !    Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
        !   IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
        IF ((Nachbars(in)%nType=='iw'.AND.TypeW=='iw').OR. &
            (Nachbars(in)%nType=='ie'.AND.TypeE=='ie')) THEN
           IF (Nachbars(in)%nType=='iw') THEN
             ix=Floor(ibn)%ix1
             jx=jx1      ! 'j' Coordinates of border current block
           ELSE
             ix=Floor(ibn)%ix0
             jx=jx0
           END IF
   !       -------------------
   !                |
   !          D     |    N
   !                |
   !       -------------------
           IF (RefineY>RefineNachbarY) THEN
   !         -------------        -------------
   !         |     |     |  'iw'  |     |     |  'ie'
   !         |     |-----|        |-----|     |
   !         |     |     |        |     |     |
   !         -------------        -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jy=jy0,jy1,IncrY
                 iy=jy/IncrY
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineZ<=RefineNachbarZ
               DO jy=jy0,jy1,IncrY
                 iy=jy/IncrY
                 DO jz=jz0,jz1
                   iz=jz*Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineY==RefineNachbarY) THEN
   !         -------------
   !         |     |     |  'iw,ie'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jy=jy0,jy1
                 iy=jy
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineZ<=RefineNachbarZ
               DO jy=jy0,jy1
                 iy=jy
                 DO jz=jz0,jz1
                   iz=jz*Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineY<RefineNachbarY) THEN
   !         -------------
   !         |  |  |     |   'iw,ie'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jy=jy0,jy1
                 iy=IncrY*jy
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/IncrZ
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE    ! RefineZ<=RefineNachbarZ
               DO jy=jy0,jy1
                 iy=IncrY*jy
                 DO jz=jz0,jz1
                   iz=jz*IncrZ
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF  ! RefineY
   
           IF (RefineZ>RefineNachbarZ) THEN
   !         -------------
   !         |     |  |  |  'iw,ie'
   !         |     |-----|
   !         |     |  |  |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jy=jy0,jy1,IncrY
                 iy=jy/IncrY
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineY<=RefineNachbarY
               DO jy=jy0,jy1
                 iy=jy*IncrY
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineZ==RefineNachbarZ) THEN
   !         -------------
   !         |     |     |   'iw,ie'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jy=jy0,jy1,IncrY
                 iy=jy/IncrY
                 DO jz=jz0,jz1
                   iz=jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineY<=RefineNachbarY
               DO jy=jy0,jy1
                 iy=jy*IncrY
                 DO jz=jz0,jz1
                   iz=jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineZ<RefineNachbarZ) THEN
   !         -------------
   !         |  |  |     |  'iw,ie'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jy=jy0,jy1,IncrY
                 iy=jy/IncrY
                 DO jz=jz0,jz1
                   iz=IncrZ*jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineY<=RefineNachbarY
               DO jy=jy0,jy1
                 iy=jy*IncrY
                 DO jz=jz0,jz1
                   iz=IncrZ*jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF ! RefineZ
        ELSE  !'ow' or. 'oe'
           !IF (Nachbars(in)%nType=='ow'.OR.Nachbars(in)%nType=='pw') THEN
           IF (Nachbars(in)%nType=='ow') THEN
             ix=ix0-1
             DO iy=iy0,iy1
               DO iz=iz0,iz1
                 Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix0,iy,iz)%in_out
               END DO
             END DO
           ELSE IF(Nachbars(in)%nType=='oe') THEN
             ix=ix1+1
             DO iy=iy0,iy1
               DO iz=iz0,iz1
                 Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix1,iy,iz)%in_out
               END DO
             END DO
           END IF
        END IF  ! 'iw, ie, ow, oe'
!......................................................................
      CASE ("s","n")
        !IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps'.OR. &
        !    Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
        !   IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
        IF ((Nachbars(in)%nType=='is'.AND.TypeS=='is').OR. &
            (Nachbars(in)%nType=='in'.AND.TypeN=='in')) THEN
           IF (Nachbars(in)%nType=='is') THEN
             iy=Floor(ibn)%iy1
             jy=jy1
           ELSE
             iy=Floor(ibn)%iy0
             jy=jy0
           END IF
   !       -------------------
   !                |
   !          N     |    D
   !                |
   !       -------------------
           IF (RefineX>RefineNachbarX) THEN
   !         -------------
   !         |     |  |  |  'is,in'
   !         |     |-----|
   !         |     |  |  |
   !         -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineZ<=RefineNachbarZ
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jz=jz0,jz1
                   iz=jz*Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineX==RefineNachbarX) THEN
   !         -------------
   !         |     |     |  'is,in'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jx=jx0,jx1
                 ix=jx
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineZ<=RefineNachbarZ
               DO jx=jx0,jx1
                 ix=jx
                 DO jz=jz0,jz1
                   iz=jz*Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineX<RefineNachbarX) THEN
   !         -------------
   !         |  |  |     |   'is,in'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineZ>RefineNachbarZ) THEN
               DO jx=jx0,jx1
                 ix=IncrX*jx
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/IncrZ
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineZ<=RefineNachbarZ
               DO jx=jx0,jx1
                 ix=IncrX*jx
                 DO jz=jz0,jz1
                   iz=jz*IncrZ
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF  ! RefineX
   
           IF (RefineZ>RefineNachbarZ) THEN
   !         -------------
   !         |     |  |  |   'is,in'
   !         |     |-----|
   !         |     |  |  |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/IncrZ
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jz=jz0,jz1,IncrZ
                   iz=jz/Incrz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineZ==RefineNachbarZ) THEN
   !         -------------
   !         |     |     |   'is,in'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jz=jz0,jz1
                   iz=jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jz=jz0,jz1
                   iz=jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineZ<RefineNachbarZ) THEN
   !         -------------
   !         |  |  |     |   'is,in'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jz=jz0,jz1
                   iz=IncrZ*jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jz=jz0,jz1
                   iz=IncrZ*jz
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF ! RefineZ
        ELSE
           !IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
           IF (Nachbars(in)%nType=='os') THEN
             iy=iy0-1
             DO ix=ix0,ix1
               DO iz=iz0,iz1
                   Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix,iy0,iz)%in_out 
               END DO
             END DO
           ELSE IF(Nachbars(in)%nType=='on') THEN
             iy=iy1+1
             DO ix=ix0,ix1
               DO iz=iz0,iz1
                  Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix,iy1,iz)%in_out 
               END DO
             END DO
           END IF
        END IF  ! 'is, in, os, on'
!...........................................................................
      CASE ("b","t")
        !IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb'.OR. &
        !    Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
        !   IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
        IF (Nachbars(in)%nType=='ib'.AND.TypeB=='ib'.OR. &
            Nachbars(in)%nType=='it'.AND.TypeT=='it') THEN
           IF (Nachbars(in)%nType=='ib') THEN
             iz=Floor(ibn)%iz1
             jz=jz1
           ELSE
             iz=Floor(ibn)%iz0
             jz=jz0
           END IF
   !       -------------------
   !                |
   !          N     |    D
   !                |
   !       -------------------
           IF (RefineX>RefineNachbarX) THEN
   !         -------------
   !         |     |  |  |   'ib,it'
   !         |     |-----|
   !         |     |  |  |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jy=jy0,jy1,IncrY
                   iy=jy/Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineY<=RefineNachbarY
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jy=jy0,jy1
                   iy=jy*Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineX==RefineNachbarX) THEN
   !         -------------
   !         |     |     |  'ib,it'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jx=jx0,jx1
                 ix=jx
                 DO jy=jy0,jy1,IncrY
                   iy=jy/Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineY<=RefineNachbarY
               DO jx=jx0,jx1
                 ix=jx
                 DO jy=jy0,jy1
                   iy=jy*Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineX<RefineNachbarX) THEN
   !         -------------
   !         |  |  |     |  'ib,it'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineY>RefineNachbarY) THEN
               DO jx=jx0,jx1
                 ix=IncrX*jx
                 DO jy=jy0,jy1,IncrY
                   iy=jy/IncrY
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineY<=RefineNachbarY
               DO jx=jx0,jx1
                 ix=IncrX*jx
                 DO jy=jy0,jy1
                   iy=jy*IncrY
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF  ! RefineX
   
           IF (RefineY>RefineNachbarY) THEN
   !         -------------
   !         |     |  |  |  'ib,it'
   !         |     |-----|
   !         |     |  |  |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jy=jy0,jy1,IncrY
                   iy=jy/Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jy=jy0,jy1,IncrY
                   iy=jy/Incry
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineY==RefineNachbarY) THEN
   !         -------------
   !         |     |     |  'ib,it'
   !         |---------- |
   !         |     |     |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jy=jy0,jy1
                   iy=jy
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE     ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jy=jy0,jy1
                   iy=jy
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           ELSE IF (RefineY<RefineNachbarY) THEN
   !         -------------
   !         |  |  |     |  'ib,it'
   !         |---- |     |
   !         |  |  |     |
   !         -------------
             IF (RefineX>RefineNachbarX) THEN
               DO jx=jx0,jx1,IncrX
                 ix=jx/IncrX
                 DO jy=jy0,jy1
                   iy=IncrY*jy
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             ELSE      ! RefineX<=RefineNachbarX
               DO jx=jx0,jx1
                 ix=jx*IncrX
                 DO jy=jy0+1,jy1
                   iy=IncrY*jy
                   IF ((Floor(ibn)%Vertices(ix,iy,iz)%nrP==-1) .OR. &
                       (Floor(ibn)%Vertices(ix,iy,iz)%nrP==0) ) THEN
                      Floor(ibn)%Vertices(ix,iy,iz)%nrP=Vertices(jx,jy,jz)%nrP
                      Floor(ibn)%Vertices(ix,iy,iz)%in_out=Vertices(jx,jy,jz)%in_out
                      Floor(ibn)%Vertices(ix,iy,iz)%nrInP=Vertices(jx,jy,jz)%nrInP
                      Floor(ibn)%Vertices(ix,iy,iz)%nrCutP=Vertices(jx,jy,jz)%nrCutP
                   END IF
                   Floor(ibn)%Vertices(ix,iy,iz)%Shift=Vertices(jx,jy,jz)%Shift
                 END DO
               END DO
             END IF
           END IF ! RefineY
        ELSE
           !IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
           IF (Nachbars(in)%nType=='ob') THEN
             iz=iz0-1
             DO ix=ix0,ix1
               DO iy=iy0,iy1
                 Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix,iy,iz0)%in_out
               END DO
             END DO
           ELSE IF(Nachbars(in)%nType=='ot') THEN
             iz=iz1+1
             DO ix=ix0,ix1
               DO iy=iy0,iy1
                 Floor(ib)%Vertices(ix,iy,iz)%in_out=Floor(ib)%Vertices(ix,iy,iz1)%in_out
               END DO
             END DO
           END IF
        END IF  ! 'ib, it, ob, ot'
    END SELECT  ! [i,o]- w,e,s,n,b,t 
    !END IF
!........................................................................................
    END DO  ! AnzahlNachbar
  END DO   ! ib

END SUBROUTINE AnalyzeAllVertices1
