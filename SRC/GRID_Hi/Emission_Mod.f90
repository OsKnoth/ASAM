MODULE Emission_Mod

  USE F_Mod
  USE Floor_Mod

  IMPLICIT NONE

  INTEGER :: ns,nlm
  INTEGER, ALLOCATABLE :: nl(:),nlength(:,:),nwidth(:,:),nthick(:,:)
  REAL(8), ALLOCATABLE :: LPoint(:,:,:,:),Width(:),Thick(:)
  REAL(8), ALLOCATABLE :: Estreet(:,:,:),Estreet0(:,:,:)
  REAL(8), ALLOCATABLE :: EPoint0(:,:,:),lstep(:,:,:),wstep(:,:,:),hstep(:,:)
  REAL(8) :: WidthM,ThickM ! Richtwerte fuer Umverteilungen in SR WriteEmission
  REAL(8), ALLOCATABLE :: EmissM(:)
  CHARACTER(5), ALLOCATABLE :: SpeciesName(:)

CONTAINS

SUBROUTINE ReadEmission

  INTEGER :: nls, is,il, ip,Number
  INTEGER :: InUnitEmi=9
  REAL(8) :: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax, hilf
  CHARACTER*5 :: ch

  Xmin= 1.d20
  Xmax=-1.d20
  Ymin= 1.d20
  Ymax=-1.d20
  Zmin= 1.d20
  Zmax=-1.d20

  OPEN(UNIT=InUnitEmi,FILE=TRIM(file_nameemi),STATUS='OLD')

! Ermittlung MAX(Linienzahl) fuer Dimensionalisierung
  nlm=0
  READ(InUnitEmi,*) SpeciesNum
  IF (SpeciesNum<=0) THEN
    WRITE(*,*) '*** STOP: SpeciesNum<=0 as read in from file',TRIM(file_nameemi),' ***'
    STOP
  END IF
  READ(InUnitEmi,*) (ch,il=1,SpeciesNum)
  READ(InUnitEmi,*) ns
  DO is=1,ns
    READ(InUnitEmi,*) Number,nls,ch
    READ(InUnitEmi,*) (ip,il=1,nls),ch
    READ(InUnitEmi,*) (ch,il=1,nls)
    DO il=1,nls
      READ(InUnitEmi,*) hilf,hilf,hilf,hilf,hilf,hilf,ch
    END DO
    READ(InUnitEmi,*) hilf,hilf,ch
    READ(InUnitEmi,*) (hilf,il=1,SpeciesNum),ch
    nlm=MAX(nlm,nls)
  END DO
  ALLOCATE(nl(ns))
  ALLOCATE(LPoint(ns,nlm,2,3),Thick(ns),Width(ns))
  ALLOCATE(SpeciesName(SpeciesNum),EmissM(SpeciesNum))
  ALLOCATE(Estreet(ns,nlm,SpeciesNum))

  REWIND(InUnitEmi)
  WidthM=0.d0
  ThickM=0.d0
  EmissM=0.d0

! ReadEmission
  READ(InUnitEmi,*) SpeciesNum
  READ(InUnitEmi,*) SpeciesName
  READ(InUnitEmi,*) ns
  DO is=1,ns
    READ(InUnitEmi,*) Number,nl(is),ch
    READ(InUnitEmi,*) (ip,il=1,nl(is)),ch
    READ(InUnitEmi,*) (ch,il=1,nl(is))
    DO il=1,nl(is)
      READ(InUnitEmi,*) LPoint(is,il,1,1),LPoint(is,il,2,1) &
                       ,LPoint(is,il,1,2),LPoint(is,il,2,2) &
                       ,LPoint(is,il,1,3),LPoint(is,il,2,3),ch
      Xmin=MIN(LPoint(is,il,1,1),LPoint(is,il,2,1),Xmin)
      Xmax=MAX(LPoint(is,il,1,1),LPoint(is,il,2,1),Xmax)
      Ymin=MIN(LPoint(is,il,1,2),LPoint(is,il,2,2),Ymin)
      Ymax=MAX(LPoint(is,il,1,2),LPoint(is,il,2,2),Ymax)
      Zmin=MIN(LPoint(is,il,1,3),LPoint(is,il,2,3),Zmin)
      Zmax=MAX(LPoint(is,il,1,3),LPoint(is,il,2,3),Zmax)
    END DO
    READ(InUnitEmi,*) Thick(is),Width(is),ch
    READ(InUnitEmi,*) (Estreet(is,1,il),il=1,SpeciesNum),ch ! µg/ms
    IF (Number<0) THEN
      Estreet(is,1,:) = 0.d0
    ELSE IF (Thick(is)<=0.d0.OR.Width(is)<=0.d0) THEN
      Estreet(is,1,:) = 0.d0
    ELSE
      Estreet(is,1,:) = MAX(Estreet(is,1,:),0.d0)
    END IF
    WidthM=MAX(  Width(is)  ,WidthM)
    ThickM=MAX(  Thick(is)  ,ThickM)
    DO il=1,SpeciesNum
      EmissM(il)=MAX(Estreet(is,1,il),EmissM(il))
    END DO
    DO il=1,nl(is)
      Estreet(is,il,:)=Estreet(is,1,:) ! Attention: emission value is for the sum of lines il
    END DO
  END DO
  CLOSE(UNIT=InUnitEmi)

  WRITE(*,*) '    Streets:   ',ns
  WRITE(*,8) '    Expansion x:',Xmin,'...',Xmax
  WRITE(*,8) '              y:',Ymin,'...',Ymax
  WRITE(*,8) '              z:',Zmin,'...',Zmax
    IF (nr_wahlemi==nr_emi_street_oro) THEN
  WRITE(*,*) '      ORO considered (but not in this scope)'
    END IF
8 FORMAT(1x,a16,f11.2,8x,a3,7x,f11.2)

  IF (     Xmin>=Domain%x0.AND.Xmax<=Domain%x1.AND. &
           Ymin>=Domain%y0.AND.Ymax<=Domain%y1.AND. &
           Zmin>=Domain%z0.AND.Zmax<=Domain%z1) THEN ! totally inside
    WRITE(*,*) '      inside of domain'
  ELSE IF (Xmin>=Domain%x1.OR. Xmax<=Domain%x0.OR. &
           Ymin>=Domain%y1.OR. Ymax<=Domain%y0) THEN ! beyond xy-domain
    STOP '***** Stop: beyond xy-domain *****'
  ELSE IF (Zmin>=Domain%z1.OR. Zmax< Domain%z0) THEN ! beyond z-domain
    WRITE(*,*) '      beyond z-domain'
  ELSE IF (Zmin>=Domain%z0.AND.Zmax<=Domain%z1) THEN ! partially inside z-domain
    WRITE(*,*) '      partially outside of xy-domain'
  ELSE
    WRITE(*,*) '      partially outside of z-domain'
  END IF

  WRITE(*,9) '    Maxima         '
  WRITE(*,9) '    Emission thick:',ThickM
  WRITE(*,9) '    Emission width:',WidthM
  WRITE(*,7) '    Species number:',SpeciesNum
  DO il=1,SpeciesNum
    WRITE(*,6) SpeciesName(il),' µg/ms',EmissM(il)
  END DO
6 FORMAT(5x,a5,1x,a6,f12.2)
9 FORMAT(1x,a19,f8.2)
7 FORMAT(1x,a19,i8)

END SUBROUTINE ReadEmission

SUBROUTINE AnalyzeEmission

!==============================================================================
!  IMPLICIT NONE
!==============================================================================
!
! Subroutine arguments:
! --------------------
  CHARACTER*50 :: FileName
!
! Local :
! -----------------------------
  INTEGER :: k,j,i, ib,ib1,ib1v,is,il
  INTEGER :: iwidth,ilength,ithick,nwidth1,nlength1,nthick1
  INTEGER :: iwidthact,ilengthact,ithickact,iwidthanf,ilengthanf,ithickanf
  INTEGER :: idir,jdir,idirw,jdirw,idirl,jdirl
  INTEGER :: iend,jend,iendw,jendw,iendl,jendl
  INTEGER :: iact,jact,iactw,jactw,iactl,jactl
  INTEGER :: ianf,janf
  INTEGER :: kdir,kend,kanf,kact,kactw,kactl
  REAL(8) :: eps=1.d-40,lmin,hmin
  REAL(8) :: disp
  REAL(8) :: lvectorx,lvectory,lvectorz,wvectorx,wvectory,hvectorz
  REAL(8) :: epoint0x,epoint0y,epoint0z
  REAL(8) :: epointw(3),epointl(3),epointh(3),epointwv(3),epointlv(3)
  REAL(8) :: length1,width1,thick1,dlength,dwidth,dthick !,lengthg
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine AnalyzeEmission 
!------------------------------------------------------------------------------

  Emiss =0.d0 ! Emission rate per cell (µg/s)
  Emiss0=0.d0 ! Emission rate density (µg/m³s) maximal possible in the cell
              ! (= Estreet0, accumulative only for different streets)
  Nemiss=0    ! Street number of foregoing street covering the cell
              ! (to detect a different street for inclusion in Emiss0)

  ALLOCATE(Estreet0(ns,nlm,SpeciesNum))
  ALLOCATE(EPoint0(ns,nlm,3),lstep(ns,nlm,3),wstep(ns,nlm,2),hstep(ns,nlm))
  ALLOCATE(nlength(ns,nlm),nwidth(ns,nlm),nthick(ns,nlm))

! Kleinste horizontale und vertikale Laengeneinheiten
  lmin=1.d20
  hmin=1.d20
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO i=ix0+1,ix1
      lmin=MIN(lmin,dx(i))
    END DO
    DO j=iy0+1,iy1
      lmin=MIN(lmin,dy(j))
    END DO
    DO k=iz0+1,iz1
      hmin=MIN(hmin,dz(k))
    END DO
  END DO

! 3D-Spreading of streets into emission points
! (epoint = corner epoint0 + width step + length step + thick step)

  DO is=1,ns
    !lengthg=0.d0
    !DO il=1,nl(is)
    !  lvectorx=LPoint(is,il,2,1)-LPoint(is,il,1,1)
    !  lvectory=LPoint(is,il,2,2)-LPoint(is,il,1,2)
    !  lvectorz=LPoint(is,il,2,3)-LPoint(is,il,1,3)
    !  length1=SQRT(lvectorx*lvectorx+lvectory*lvectory) ! linear extension
    !  lengthg=lengthg+length1 ! Street length relevant for emission
    !END DO

    DO il=1,nl(is)
      lvectorx=LPoint(is,il,2,1)-LPoint(is,il,1,1)
      lvectory=LPoint(is,il,2,2)-LPoint(is,il,1,2)
      lvectorz=LPoint(is,il,2,3)-LPoint(is,il,1,3)

      length1=SQRT(lvectorx*lvectorx+lvectory*lvectory) ! linear extension
      width1 =Width(is)
      thick1 =Thick(is)

      nlength1=MAX(10.d0*length1/lmin,1.1d0) ! numbers of partitions
      nwidth1 =MAX(10.d0*width1 /lmin,1.1d0) ! 10.d0 = 10 points per cell
      nthick1 =MAX(10.d0*thick1 /hmin,1.1d0)

      dlength=length1/nlength1 ! step sizes
      dwidth =width1 /nwidth1
      dthick =thick1 /nthick1

      lvectorx=lvectorx/(length1+eps) ! unit length vector
      lvectory=lvectory/(length1+eps)
      lvectorz=lvectorz/(length1+eps)
      wvectorx=-lvectory ! unit width vector (perpendicular to length vector)
      wvectory= lvectorx

      epoint0x=LPoint(is,il,1,1)-wvectorx*width1*0.5d0 ! street corner point
      epoint0y=LPoint(is,il,1,2)-wvectory*width1*0.5d0
      epoint0z=LPoint(is,il,1,3)

!      nlength1=length1+width1+1.d0            ! Case of very small cells:
!      dlength=(length1+width1)/nlength1       ! length extension by width/2 at each end
!      epoint0x=epoint0x-lvectorx*width1*0.5d0
!      epoint0y=epoint0y-lvectory*width1*0.5d0

      lvectorx=lvectorx*dlength ! step vectors (incl. thick)
      lvectory=lvectory*dlength
      lvectorz=lvectorz*dlength
      wvectorx=wvectorx*dwidth
      wvectory=wvectory*dwidth
      hvectorz=dthick

      epoint0x=epoint0x+lvectorx*0.5d0+wvectorx*0.5d0 ! first emission point
      epoint0y=epoint0y+lvectory*0.5d0+wvectory*0.5d0 ! (= epoint0 shifted by half steps)
      epoint0z=epoint0z+lvectorz*0.5d0+hvectorz*0.5d0

      EPoint0(is,il,1)=epoint0x
      EPoint0(is,il,2)=epoint0y
      EPoint0(is,il,3)=epoint0z
      lstep  (is,il,1)=lvectorx
      lstep  (is,il,2)=lvectory
      lstep  (is,il,3)=lvectorz
      wstep  (is,il,1)=wvectorx
      wstep  (is,il,2)=wvectory
      hstep  (is,il)  =hvectorz
      nlength(is,il)  =nlength1
      nwidth (is,il)  =nwidth1
      nthick (is,il)  =nthick1

      !   Up to here: Estreet  = Emission rate of street per length (µg/ms)
      !   Henceforth: Estreet0 = Emission rate of street per volume unit (µg/m³s)
      !               Estreet  = Emission rate of individual street points (µg/s)

      Estreet0(is,il,:)=Estreet (is,il,:)/(width1*thick1) ! (lengthg*width1*thick1)
      Estreet (is,il,:)=Estreet0(is,il,:)*(dlength*dwidth*dthick)

!     Estreet=0 for street parts (is,il) definitely outside of domain:
!     (to shorten the following main loop)
      IF ((LPoint(is,il,1,1)<=domain%x0.AND.LPoint(is,il,2,1)<=domain%x0).OR. &
          (LPoint(is,il,1,1)>=domain%x1.AND.LPoint(is,il,2,1)>=domain%x1).OR. &
          (LPoint(is,il,1,2)<=domain%y0.AND.LPoint(is,il,2,2)<=domain%y0).OR. &
          (LPoint(is,il,1,2)>=domain%y1.AND.LPoint(is,il,2,2)>=domain%y1)) THEN
        Estreet (is,il,:)=0.d0
        Estreet0(is,il,:)=0.d0
      END IF
    END DO
  END DO

  DEALLOCATE(LPoint)

! Filling the emission points in the grid cells

  DO is=1,ns
    DO il=1,nl(is)

      disp=0.d0
      DO i=1,SpeciesNum
        disp=disp+Estreet(is,il,i)
      END DO
      IF (disp==0.d0) GOTO 12

      iwidthact =1
      ilengthact=1
      ithickact =1
      nlength1  =nlength(is,il)
      nwidth1   =nwidth (is,il)
      nthick1   =nthick (is,il)
      epointw(1)=epoint0(is,il,1)
      epointw(2)=epoint0(is,il,2)
      epointw(3)=epoint0(is,il,3)
      epointl(1)=epoint0(is,il,1)
      epointl(2)=epoint0(is,il,2)
      epointl(3)=epoint0(is,il,3)
      epointh(1)=epoint0(is,il,1)
      epointh(2)=epoint0(is,il,2)
      epointh(3)=epoint0(is,il,3)
      IF (nr_wahlemi==nr_emi_street_oro) THEN ! + Orographie
        epointh(3)=epointh(3)-foro(epointh(1),epointh(2),0.d0)
      END IF

      ib1=-1

1     CONTINUE ! search for block

!     Block-Bestimmung ib1 (primaer in x-y-Richtung)
      ib1v=ib1
      ib1=0 ! darf nur bei xy-offside so bleiben
            ! (bei z-offside muss xy-Block gefunden werden)
      DO ib=1,nb
        CALL Set(Floor(ib))
        IF   (epointh(1)>=xP(ix0).AND.epointh(1)<xP(ix1)) THEN
          IF (epointh(2)>=yP(iy0).AND.epointh(2)<yP(iy1)) THEN
            ib1=ib
            ! sichert Erreichen des thick-Loops ohne staendiges Springen
            ! zur Blocksuche im Falle von nur in z-Richtung vorliegender Abseitslage
            IF (epointh(3)>=zP(iz0).AND.epointh(3)<zP(iz1)) THEN
              ib1=ib ! Block gefunden
              EXIT
            END IF
          END IF
        END IF
      END DO ! Falls ausserhalb der 3d-domain --> Block nb, aber ib1=0
      IF (ib1v/=-1.AND.(ib1==0.OR.ib1==ib1v)) STOP ' *** Error: Blocks do not cover whole domain ***'
      IF (ib1>0) CALL Set(Floor(ib1))
!     --> ausserhalb xy-domain: kein zutreffender xy-Block
!         innerhalb  xy-domain:      zutreffender xy-Block
!                               (evt. aber ausserhalb z-domain)

      iwidthanf =iwidthact
      ilengthanf=ilengthact
      ithickanf =ithickact

!     Anfang,Ende,Richtung der ij-Loops bzgl. width steps (max. Bereich)
      IF (wstep(is,il,1)>=0.d0) THEN
        idirw=1
        iendw=ix1
        iactw=ix0+1
      ELSE
        idirw=-1
        iendw=ix0+1
        iactw=ix1
      END IF
      IF (wstep(is,il,2)>=0.d0) THEN
        jdirw=1
        jendw=iy1
        jactw=iy0+1
      ELSE
        jdirw=-1
        jendw=iy0+1
        jactw=iy1
      END IF
!     Anfang,Ende,Richtung der ij-Loops bzgl. length steps (max. Bereich)
      IF (lstep(is,il,1)>=0.d0) THEN
        idirl=1
        iendl=ix1
        iactl=ix0+1
      ELSE
        idirl=-1
        iendl=ix0+1
        iactl=ix1
      END IF
      IF (lstep(is,il,2)>=0.d0) THEN
        jdirl=1
        jendl=iy1
        jactl=iy0+1
      ELSE
        jdirl=-1
        jendl=iy0+1
        jactl=iy1
      END IF
!     Anfang des k-Loops
      kactw=iz0+1
      kactl=iz0+1

      DO iwidth=iwidthanf,nwidth1 ! width-Loop ====================
        iwidthact=iwidth

        epointwv(1)=epointw(1)
        epointwv(2)=epointw(2)
        epointwv(3)=epointw(3)
        epointw (1)=epoint0(is,il,1)+(iwidth-1)*wstep(is,il,1)
        epointw (2)=epoint0(is,il,2)+(iwidth-1)*wstep(is,il,2)
        epointw (3)=epoint0(is,il,3)
        IF (nr_wahlemi==nr_emi_street_oro) THEN ! + Orographie
          epointw(3)=epointw(3)-foro(epointw(1),epointw(2),0.d0)
        END IF

        DO ilength=ilengthanf,nlength1 ! length-Loop ====================
          ilengthact=ilength

          epointlv(1)=epointl(1)
          epointlv(2)=epointl(2)
          epointlv(3)=epointl(3)
          epointl (1)=epointw(1)+(ilength-1)*lstep(is,il,1)
          epointl (2)=epointw(2)+(ilength-1)*lstep(is,il,2)
          epointl (3)=epoint0(is,il,3)+(ilength-1)*lstep(is,il,3)
          IF (nr_wahlemi==nr_emi_street_oro) THEN ! + Orographie
            epointl(3)=epointl(3)-foro(epointl(1),epointl(2),0.d0)
          END IF

          epointh(1)=epointl(1)
          epointh(2)=epointl(2)
          epointh(3)=epointl(3) ! unterster Punkt des thick-Loops

          IF (epointh(1)<domain%x0.OR.epointh(1)>=domain%x1.OR. &
              epointh(2)<domain%y0.OR.epointh(2)>=domain%y1) THEN
            ib1=0
            GOTO 11 ! outside xy-domain (definitely outside 3d-domain)
          ELSE IF (ib1==0) THEN
            ithickact=ithickanf
            GOTO 1 ! afore outside xy-domain (now inside)
          ELSE IF (epointh(1)<xP(ix0).OR.epointh(1)>=xP(ix1).OR. &
                   epointh(2)<yP(iy0).OR.epointh(2)>=yP(iy1)) THEN
            ithickact=ithickanf
            GOTO 1 ! outside actual xy-block (selbst wenn ausserhalb z-domain)
          END IF
          ! --> xy-Block zutreffend (aber uneindeutig oder evt. ausserhalb z-domain)

          IF (ilength==1) THEN
!           Anfang,Ende,Richtung der ij-Loops (fuer width-Loop)
            ianf=iactw
            janf=jactw
            iend=iendw
            jend=jendw
            idir=idirw
            jdir=jdirw
!           Anfang,Ende,Richtung des k-Loops (fuer width-Loop)
            kact=kactw
            IF (epointw(3)>=epointwv(3)) THEN ! (fuer ithick=1)
              kdir=1
              kend=iz1
            ELSE
              kdir=-1
              kend=iz0+1
            END IF
          ELSE
!           Anfang,Ende,Richtung der ij-Loops (fuer length-Loop)
            ianf=iactl
            janf=jactl
            iend=iendl
            jend=jendl
            idir=idirl
            jdir=jdirl
!           Anfang,Ende,Richtung des k-Loops (fuer length-Loop)
            kact=kactl
            IF (epointl(3)>=epointlv(3)) THEN ! (fuer ithick=1)
              kdir=1
              kend=iz1
            ELSE
              kdir=-1
              kend=iz0+1
            END IF
          END IF

          DO i=ianf,iend,idir ! i-Zellen-Findung --------------------
            iact=i
            IF (epointh(1)>=xP(i-1).AND.epointh(1)<xP(i)) THEN
              EXIT ! i found, goto j-Loop
            END IF
          END DO ! i
          DO j=janf,jend,jdir ! j-Zellen-Findung --------------------
            jact=j
            IF (epointh(2)>=yP(j-1).AND.epointh(2)<yP(j)) THEN
              EXIT ! j found, goto k-Loop
            END IF
          END DO ! j

          DO ithick=ithickanf,nthick1 ! thick-Loop ====================
            ithickact=ithick

            epointh(3)=epointl(3)+(ithick-1)*hstep(is,il)

            IF (epointh(3)<domain%z0.OR.epointh(3)>=domain%z1) THEN
              ib1=0
              GOTO 10 ! outside z-domain (definetely outside 3d-domain)
            ELSE IF (ib1==0) THEN
              GOTO 1 ! afore outside z-domain (now inside)
            ELSE IF (epointh(3)<zP(iz0).OR.epointh(3)>=zP(iz1)) THEN
              GOTO 1 ! outside actual z-block (definitely inside 3d-domain)
            END IF
            ! --> xyz-Block zutreffend

            kanf=kact
            IF (ithick>1) THEN
              kdir=1
              kend=iz1
            END IF

            DO k=kanf,kend,kdir ! k-Zellen-Findung --------------------
              kact=k
              IF (epointh(3)>=zP(k-1).AND.epointh(3)<zP(k)) THEN
                EXIT ! k found, goto emission
              END IF
            END DO ! k

            ! Emiss-Rate (nicht Dichte)
            Emiss(iact,jact,kact,:)=Emiss(iact,jact,kact,:)+Estreet(is,il,:)

            ! Maximal moegliche Emiss-Dichte (in einer Zelle) speichern (Kreuzungen)
            IF (is>Nemiss(iact,jact,kact)) THEN
              !Nemiss(iact,jact,kact)  =is
              !Emiss0(iact,jact,kact,:)=Emiss0(iact,jact,kact,:)+Estreet0(is,il,:) &
              !                        *MIN(Thick(is),dz(kact))/dz(kact)
              Nemiss(iact,jact,kact)  =0
              DO k=1,SpeciesNum
                Emiss0(iact,jact,kact,k)=MAX(Emiss0(iact,jact,kact,k),Estreet0(is,il,k) &
                                             *MIN(Thick(is),dz(kact))/dz(kact))
              END DO
            END IF

            IF (ithick==1) THEN
              iactl=iact
              jactl=jact
              kactl=kact
              IF (ilength==1) THEN
                iactw=iact
                jactw=jact
                kactw=kact
              END IF
            END IF

10          CONTINUE

          END DO ! thick

11        CONTINUE
          ithickanf=1

        END DO ! length
        ilengthanf=1

      END DO ! width

12    CONTINUE

    END DO ! il
  END DO ! is

  DEALLOCATE(EPoint0,lstep,wstep,hstep,nl,estreet,estreet0)
  DEALLOCATE(nlength,nwidth,nthick,Width,Thick)

END SUBROUTINE AnalyzeEmission

SUBROUTINE WriteEmission(FileName)

!==============================================================================
!  IMPLICIT NONE
!==============================================================================
!
! Subroutine arguments:
! --------------------
  CHARACTER*50 :: FileName
!
! Local :
! -----------------------------
  INTEGER :: k,j,i, ib,is, inum
  INTEGER :: ii,jj,kk,i1,i2,j1,j2,k1,k2,iid,jjd,kkd
  INTEGER :: OutUnitEmi=10
  REAL(8) :: de0,de1 ,deTol=0.01d0 ! relative Toleranz
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine WriteEmission 
!------------------------------------------------------------------------------

  WRITE(*,*) '  Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=OutUnitEmi,FILE=TRIM(FileName)//'.Emission',STATUS='unknown')
  WRITE(*,*) '  > ',TRIM(FileName),'.Emission'
! Emiss0 here similar to Emiss: max. Emission rate (not density)
! Nemiss here marking cells with Emiss>0

  Nemiss=0
  WRITE(OutUnitEmi,*) SpeciesNum
  WRITE(OutUnitEmi,*) (SpeciesName(i),' ',i=1,SpeciesNum)

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          Emiss0(i,j,k,:)=(1.d0+deTol)*Emiss0(i,j,k,:)*VolC(i,j,k) ! max. moegl. Emission rate (not density)
        END DO
      END DO
    END DO
    is=0
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          de0=0.d0
          DO inum=1,SpeciesNum
            de0=MAX(Emiss(i,j,k,inum),de0)
          END DO
          IF (de0>0.d0) THEN
            Nemiss(i,j,k)=1 ! Emission>0
          END IF
          IF (Nemiss(i,j,k)>0) THEN
            Nemiss(i,j,k)=0
            DO inum=1,SpeciesNum
              de0=Emiss(i,j,k,inum)-Emiss0(i,j,k,inum) ! Umverteilungsmasse
              IF (de0>0.d0) THEN
                iid=WidthM/dx(i)+2.d0
                jjd=WidthM/dy(j)+2.d0
                kkd=ThickM/dz(k)+1.d0
                k1=MAX(k-kkd,iz0+1)
                j1=MAX(j-jjd,iy0+1)
                i1=MAX(i-iid,ix0+1)
                k2=MIN(k+kkd,iz1)
                j2=MIN(j+jjd,iy1)
                i2=MIN(i+iid,ix1)
                de1=0.d0
                DO kk=k1,k2
                  DO jj=j1,j2
                    DO ii=i1,i2
                      de1=de1+MAX(Emiss0(ii,jj,kk,inum)-Emiss(ii,jj,kk,inum),0.d0) ! Aufnahmevermoegen
                    END DO
                  END DO
                END DO
                IF (VolC(i,j,k)==0.d0.AND.inum==1) THEN
                  WRITE(*,11) '  Emission value in closed cell ',ib,i,j,k
                END IF
                IF (de0>de1) WRITE(*,10) '  Emission loss:',ib,i,j,k,inum,VolC(i,j,k),de0-de1,Emiss(i,j,k,inum)
11              FORMAT(1x,a32,4i5)
10              FORMAT(1x,a16,i4,3i5,4x,i5,f11.2,5x,f8.3,'  of',f8.3)
                IF (de1>0.d0) THEN
                  de0=MIN(de0,de1)/de1 ! Maximale Umverteilungsmasse/Aufnahmevermoegen
                  DO kk=k1,k2
                    DO jj=j1,j2
                      DO ii=i1,i2
                        de1=Emiss0(ii,jj,kk,inum)-Emiss(ii,jj,kk,inum)
                        IF (de1>0.d0) THEN
                          Emiss(ii,jj,kk,inum)=Emiss(ii,jj,kk,inum)+de0*de1
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
                Emiss(i,j,k,inum)=Emiss0(i,j,k,inum)
              END IF
              IF (Emiss(i,j,k,inum)>0.d0.AND.Nemiss(i,j,k)==0) THEN ! numbering
                Nemiss(i,j,k)=1
                is=is+1
              END IF
            END DO
          END IF
        END DO
      END DO
    END DO

    WRITE(OutUnitEmi,*) ib,is,'!block,cells (EMISSION)' ! point sources, blockwise
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (Nemiss(i,j,k)>0) THEN
            WRITE(OutUnitEmi,*) i,j,k
            WRITE(OutUnitEmi,*) Emiss(i,j,k,:) ! µg/s
          END IF
        END DO
      END DO
    END DO
  END DO ! ib
  DEALLOCATE(SpeciesName)
  CLOSE(UNIT=OutUnitEmi)

END SUBROUTINE WriteEmission

END MODULE Emission_Mod
