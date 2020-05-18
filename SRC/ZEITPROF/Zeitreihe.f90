PROGRAM Container_Built
!  Programm, use for time average Gnu-output-files
  INTEGER ix,cells,nFrac
  INTEGER, PARAMETER :: RealKind=8 
  REAL(REALKIND) ::  a
  CHARACTER(80) :: post

! Count the files
  CALL FileCount(ix)

  CALL getarg(3,post)

    IF (Post=='Met'.OR.Post=='Special'.OR.Post=='Besonderes') THEN
! Read Time from Gridfile -> a
! Read number of vertical cells from the Gridfile -> cells
      CALL ReadGridGas(a,cells)
! Create time series in every vertical cell
      CALL ZeitreiheMet(ix,a,cells)
    ELSE
      DO k=1,80
        IF (post(k:k).EQ.'p'.OR.post(k:k)=='m'.OR.post(k:k)=='a'.OR.post(k:k)=='s') THEN ! Aerosol
! Read Time from Gridfile -> a
! Read number of vertical cells from the Gridfile -> cells
          CALL ReadGridFl(a,cells,nFrac)

! Create time series in every vertical cell
          CALL ZeitreiheFl(ix,a,cells,nFrac)
          EXIT
        ENDIF
        IF (post(k:k).EQ.' ') THEN ! Gas
! Read Time from Gridfile -> a
! Read number of vertical cells from the Gridfile -> cells
          CALL ReadGridGas(a,cells)

! Create time series in every vertical cell
          CALL ZeitreiheGas(ix,a,cells)
          EXIT
        ENDIF
      ENDDO
    ENDIF

END


SUBROUTINE FileCount(ix)

  CHARACTER filename*100
  INTEGER ix

  ix=0

  OPEN(UNIT=77, FILE='INPUT_ZEITMITTELUNG')
  DO 
    READ(77,*,END=100) filename
    ix=ix+1
  ENDDO
100    CONTINUE
  CLOSE(77)

END SUBROUTINE FileCount


SUBROUTINE ReadGridGas(a,cells)

  INTEGER, PARAMETER :: RealKind=8 
  CHARACTER dummy*1000
  CHARACTER InputFileName*1000
  CHARACTER Timeinput*3
  INTEGER cells
  REAL(REALKIND) a

  CALL getarg(1,InputFileName)

  OPEN(UNIT=66, FILE=InputFileName)
  DO
    READ(66,*) dummy
    IF (dummy(1:15).EQ.'OutputTimeStep=') THEN
      DO k=1,1000
        IF (dummy(k:k).EQ.'=') THEN
          dummy(k:k) =' '
        ELSEIF (dummy(k:k).EQ.',') THEN
          dummy(k:k) =' '
        ENDIF
      ENDDO
      READ (dummy,*) dummy(1:14),a
      EXIT
     ENDIF

  ENDDO

REWIND(66)

  DO
    READ(66,*) dummy
    IF (dummy.EQ.'#OutputDomain') THEN
    READ(66,*) dummy
    READ(66,*) dummy
    READ(66,*) dummy
    READ(66,*) dummy,cells
      EXIT
    ENDIF
  ENDDO
  CLOSE(66)

  CALL getarg(2,Timeinput)

  IF (Timeinput.EQ.'sek') THEN
    a=a
  ELSEIF (Timeinput.EQ.'min') THEN
    a=a/60  ! Minuten
  ELSEIF (Timeinput.EQ.'hou') THEN
    a=a/3600 ! Stunden
  ELSEIF (Timeinput.EQ.'day') THEN
    a=a/86400 ! Tage
  ENDIF

END SUBROUTINE ReadGridGas

SUBROUTINE ReadGridFl(a,cells,nFrac)

  INTEGER, PARAMETER :: RealKind=8 
  CHARACTER dummy*1000
  CHARACTER InputFileName*1000
  CHARACTER Timeinput*3
  INTEGER cells,nFrac
  REAL(REALKIND) a

  CALL getarg(1,InputFileName)

  OPEN(UNIT=66, FILE=InputFileName)
  DO
    READ(66,*) dummy
    IF (dummy(1:6).EQ.'nFrac=') THEN
      DO k=1,1000
        IF (dummy(k:k).EQ.'=') THEN
          dummy(k:k) =' '
        ELSEIF (dummy(k:k).EQ.',') THEN
          dummy(k:k) =' '
        ENDIF
      ENDDO
      READ (dummy,*) dummy(1:5),nFrac
      EXIT
     ENDIF

  ENDDO

REWIND(66)

  DO
    READ(66,*) dummy
    IF (dummy(1:15).EQ.'OutputTimeStep=') THEN
      DO k=1,1000
        IF (dummy(k:k).EQ.'=') THEN
          dummy(k:k) =' '
        ELSEIF (dummy(k:k).EQ.',') THEN
          dummy(k:k) =' '
        ENDIF
      ENDDO
      READ (dummy,*) dummy(1:14),a
      EXIT
     ENDIF

  ENDDO

REWIND(66)

  DO
    READ(66,*) dummy
    IF (dummy.EQ.'#OutputDomain') THEN
    READ(66,*) dummy
    READ(66,*) dummy
    READ(66,*) dummy
    READ(66,*) dummy,cells
      EXIT
    ENDIF
  ENDDO
  CLOSE(66)

  CALL getarg(2,Timeinput)

  IF (Timeinput.EQ.'sek') THEN
    a=a
  ELSEIF (Timeinput.EQ.'min') THEN
    a=a/60  ! Minuten
  ELSEIF (Timeinput.EQ.'hou') THEN
    a=a/3600 ! Stunden
  ELSEIF (Timeinput.EQ.'day') THEN
    a=a/86400 ! Tage
  ENDIF

END SUBROUTINE ReadGridFl

SUBROUTINE ReadMode(r,nFrac)

  INTEGER :: i,nFrac
  INTEGER, PARAMETER :: RealKind=8
  CHARACTER dummy*1000
  CHARACTER InputFileName*1000
  REAL(RealKind) :: mC,rStart,pFac
  REAL(REALKIND) :: r(nFrac),m(nFrac+1)
  REAL(RealKind) :: RhoW=1.0d3
  REAL(RealKind) :: Pi=2.0d0*ASIN(1.0d0)

  CALL getarg(1,InputFileName)

  OPEN(UNIT=66, FILE=InputFileName)

  DO
    READ(66,*) dummy
    IF (dummy(1:7).EQ.'rStart=') THEN
      DO k=1,1000
        IF (dummy(k:k).EQ.'=') THEN
          dummy(k:k) =' '
        ELSEIF (dummy(k:k).EQ.',') THEN
          dummy(k:k) =' '
        ENDIF
      ENDDO
      READ (dummy,*) dummy(1:6),rStart
      EXIT
     ENDIF

  ENDDO

REWIND(66)

  DO
    READ(66,*) dummy
    IF (dummy(1:5).EQ.'pFac=') THEN
      DO k=1,1000
        IF (dummy(k:k).EQ.'=') THEN
          dummy(k:k) =' '
        ELSEIF (dummy(k:k).EQ.',') THEN
          dummy(k:k) =' '
        ENDIF
      ENDDO
      READ (dummy,*) dummy(1:4),pFac
      EXIT
     ENDIF

  ENDDO

  CLOSE(66)

  m(1)=4.0d0/3.0d0*Pi*RhoW*rStart**3
  DO i=1,nFrac
    m(i+1)=m(i)*pFac
    mC=0.5d0*(m(i)+m(i+1))
    r(i)=(3.0d0*mC/4.0d0/Pi/RhoW)**(1.0d0/3.0d0) ! in m
  END DO

END SUBROUTINE ReadMode

SUBROUTINE ZeitreiheFl(ix,a,cells,Werte)

  INTEGER i0,ix,i,x,r,z,cells
  INTEGER, PARAMETER :: RealKind=8 
  INTEGER :: Werte
  REAL(REALKIND) ::  GnuplotOut(0:Werte)
  REAL(REALKIND) ::  Zeitausgabe(cells,Werte,ix)
  REAL(REALKIND) ::  a
  REAL(REALKIND) ::  Rad(Werte)
  CHARACTER filename*100
  CHARACTER dummy*1000
  CHARACTER hoehe*10

  CALL ReadMode(Rad,Werte)

  i0=1

  DO z=1,cells
    WRITE(hoehe,'(I8)') z ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_timeaverage_cell_'//TRIM(ADJUSTL(hoehe)), status='replace')
    OPEN(UNIT=77, FILE='INPUT_ZEITMITTELUNG')

! Erzeugung der Zeitreihe in der angewählten Höhe
    DO x=i0,ix

      READ(77,*) filename
      OPEN(UNIT=88, FILE=filename)

! Weglesen der unteren Zeilen wenn höhere Ausgelesen werden
        DO k=1,z
          READ(88,*) dummy
        ENDDO

! Zuweisen der Werte zur entsprechenden Zeitreihe
      READ(88,*) GnuplotOut(:)
      DO i=1,Werte
        Zeitausgabe(z,i,x)=GnuplotOut(i)
      ENDDO

      CLOSE(88)
      WRITE(*,*)filename
      WRITE(*,*)z,x,GnuplotOut(0)
    ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
    DO r=1,ix
      WRITE(99,*)(r-1)*a,Zeitausgabe(z,:,r)
    ENDDO
    CLOSE(77)
    CLOSE(99)
  ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
  DO z=1,cells
    WRITE(hoehe,'(I8)') z ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_cellaverage_trC_'//TRIM(ADJUSTL(hoehe)), status='replace')
    DO r=1,ix ! Zeit
    DO i=1,Werte ! Klassen
      WRITE(99,*) (r-1)*a,Rad(i),Zeitausgabe(z,i,r)
    ENDDO
      WRITE(99,*) '  '
    ENDDO
    CLOSE(99)
  ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
  DO r=1,ix ! Zeit
    WRITE(hoehe,'(I8)') r ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_timeaverage_rzC_'//TRIM(ADJUSTL(hoehe)), status='replace')
    DO i=1,Werte ! Klassen
    DO z=1,cells
      WRITE(99,*) Rad(i),z*10,Zeitausgabe(z,i,r)
    ENDDO
      WRITE(99,*) '  '
    ENDDO
    CLOSE(99)
  ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
  DO i=1,Werte ! Klassen
    WRITE(hoehe,'(I8)') i ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_modeaverage_tzC_'//TRIM(ADJUSTL(hoehe)), status='replace')
    DO r=1,ix ! Zeit
    DO z=1,cells
      WRITE(99,*) (r-1)*a,z*10,Zeitausgabe(z,i,r)
    ENDDO
      WRITE(99,*) '  '
    ENDDO
    CLOSE(99)
  ENDDO

END SUBROUTINE ZeitreiheFl

SUBROUTINE ZeitreiheGas(ix,a,cells)

  INTEGER i0,ix,i,x,r,z,cells
  INTEGER, PARAMETER :: RealKind=8 
  INTEGER, PARAMETER :: Werte=1
  REAL(REALKIND) ::  GnuplotOut(0:Werte)
  REAL(REALKIND) ::  Zeitausgabe(cells,Werte,ix)
  REAL(REALKIND) ::  a
  CHARACTER filename*100
  CHARACTER dummy*1000
  CHARACTER hoehe*10

  i0=1

  DO z=1,cells
    WRITE(hoehe,'(I8)') z ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_timeaverage_cell_'//TRIM(ADJUSTL(hoehe)), status='replace')
    OPEN(UNIT=77, FILE='INPUT_ZEITMITTELUNG')

! Erzeugung der Zeitreihe in der angewählten Höhe
    DO x=i0,ix

      READ(77,*) filename
      OPEN(UNIT=88, FILE=filename)

! Weglesen der unteren Zeilen wenn höhere Ausgelesen werden
        DO k=1,z
          READ(88,*) dummy
        ENDDO

! Zuweisen der Werte zur entsprechenden Zeitreihe
      READ(88,*) GnuplotOut(:)
      DO i=1,Werte
        Zeitausgabe(z,i,x)=GnuplotOut(i)
      ENDDO

      CLOSE(88)
      WRITE(*,*)filename
      WRITE(*,*)z,x,GnuplotOut(0)
    ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
    DO r=1,ix
      WRITE(99,*)(r-1)*a,Zeitausgabe(z,1,r)
    ENDDO
    CLOSE(77)
    CLOSE(99)
  ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
  OPEN(UNIT=99, FILE='Gas_3d_Darstellung', status='replace')
  DO r=1,ix ! Zeit
  DO z=1,cells
    WRITE(99,*) (r-1)*a,z,Zeitausgabe(z,1,r)
  ENDDO
    WRITE(99,*) '  '
  ENDDO
  CLOSE(99)

END SUBROUTINE ZeitreiheGas

SUBROUTINE ZeitreiheMet(ix,a,cells)

  INTEGER i0,ix,i,x,r,z,cells
  INTEGER, PARAMETER :: RealKind=8 
  INTEGER :: Werte=-1
  REAL(REALKIND), ALLOCATABLE ::  GnuplotOut(:)
  REAL(REALKIND), ALLOCATABLE ::  Zeitausgabe(:,:,:)
  REAL(REALKIND) ::  a
  CHARACTER filename*100
  CHARACTER dummy*10000
  CHARACTER hoehe*10

  i0=1

! Bestimmen der Anzahl der ausgegebenen Spalten
  OPEN(UNIT=77, FILE='INPUT_ZEITMITTELUNG')
  READ(77,*) filename
  OPEN(UNIT=88, FILE=filename)
  READ(88,*) dummy
  READ(88,"(A10000)") dummy
  DO k=1,LEN_TRIM(dummy)
    IF (dummy(k:k)=='.') THEN
      Werte=Werte+1
    END IF
  END DO
  CLOSE(88)
  CLOSE(77)
  ALLOCATE (GnuplotOut(0:Werte))
  ALLOCATE (Zeitausgabe(cells,Werte,ix))

  DO z=1,cells
    WRITE(hoehe,'(I8)') z ! Durchnummerieren der Ausgabedatei
    OPEN(UNIT=99, FILE='Gnu_timeaverage_cell_'//TRIM(ADJUSTL(hoehe)), status='replace')
    OPEN(UNIT=77, FILE='INPUT_ZEITMITTELUNG')

! Erzeugung der Zeitreihe in der angewählten Höhe
    DO x=i0,ix

      READ(77,*) filename
      OPEN(UNIT=88, FILE=filename)

! Weglesen der unteren Zeilen wenn höhere Ausgelesen werden
        DO k=1,z
          READ(88,*) dummy
        ENDDO

! Zuweisen der Werte zur entsprechenden Zeitreihe
      READ(88,*) GnuplotOut(:)
      DO i=1,Werte
        Zeitausgabe(z,i,x)=GnuplotOut(i)
      ENDDO

      CLOSE(88)
      WRITE(*,*)filename
      WRITE(*,*)z,x,GnuplotOut(0)
    ENDDO

! Schreiben der Zeitreihe für alle Variablen in entsprechend nummerierte Datei
    DO r=1,ix
      WRITE(99,*)(r-1)*a,Zeitausgabe(z,:,r)
    ENDDO
    CLOSE(77)
    CLOSE(99)
  ENDDO
  DEALLOCATE (GnuplotOut)
  DEALLOCATE (Zeitausgabe)

END SUBROUTINE ZeitreiheMet
