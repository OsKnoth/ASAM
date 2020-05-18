MODULE Masking_Mod

  USE F_Mod
  USE Floor_Mod

  IMPLICIT NONE

! Analogie Emission_Mod(Street) <--> Masking_Mod(Baumbestand)
  INTEGER :: nsMask,nlmMask
  INTEGER, ALLOCATABLE :: nlMask(:),nlengthMask(:,:),nwidthMask(:,:),nthickMask(:,:)
  REAL(8), ALLOCATABLE :: LPointMask(:,:,:,:),WidthMask(:),ThickMask(:)
  REAL(8), ALLOCATABLE :: Ebaum(:,:),Ebaum0(:,:)
  REAL(8), ALLOCATABLE :: EPoint0Mask(:,:,:),lstepMask(:,:,:),wstepMask(:,:,:),hstepMask(:,:)
  REAL(8) :: WidthMaskM,ThickMaskM,MaskM

CONTAINS

SUBROUTINE ReadMasking

  INTEGER :: nls, is,il, ip,Number
  INTEGER :: InUnitMask=9
  REAL(8) :: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax, hilf
  CHARACTER*5 :: ch

  Xmin= 1.d20
  Xmax=-1.d20
  Ymin= 1.d20
  Ymax=-1.d20
  Zmin= 1.d20
  Zmax=-1.d20

  OPEN(UNIT=InUnitMask,FILE=TRIM(file_namemask),STATUS='OLD')

! Ermittlung MAX(Linienzahl) fuer Dimensionalisierung
  nlmMask=0
  READ(InUnitMask,*) nsMask
  DO is=1,nsMask
    READ(InUnitMask,*) Number,nls,ch
    READ(InUnitMask,*) (ip,il=1,nls),ch
    READ(InUnitMask,*) (ch,il=1,nls)
    DO il=1,nls
      READ(InUnitMask,*) hilf,hilf,hilf,hilf,hilf,hilf,ch
    END DO
    READ(InUnitMask,*) hilf,hilf,ch
    READ(InUnitMask,*) hilf,ch
    nlmMask=MAX(nlmMask,nls)
  END DO
  ALLOCATE(nlMask(nsMask))
  ALLOCATE(LPointMask(nsMask,nlmMask,2,3),ThickMask(nsMask),WidthMask(nsMask))
  ALLOCATE(Ebaum(nsMask,nlmMask))

  REWIND(InUnitMask)
  WidthMaskM=0.d0
  ThickMaskM=0.d0
  MaskM=0.d0

! ReadMask
  READ(InUnitMask,*) nsMask
  DO is=1,nsMask
    READ(InUnitMask,*) Number,nlMask(is),ch
    READ(InUnitMask,*) (ip,il=1,nlMask(is)),ch
    READ(InUnitMask,*) (ch,il=1,nlMask(is))
    DO il=1,nlMask(is)
      READ(InUnitMask,*) LPointMask(is,il,1,1),LPointMask(is,il,2,1) &
                        ,LPointMask(is,il,1,2),LPointMask(is,il,2,2) &
                        ,LPointMask(is,il,1,3),LPointMask(is,il,2,3),ch
      Xmin=MIN(LPointMask(is,il,1,1),LPointMask(is,il,2,1),Xmin)
      Xmax=MAX(LPointMask(is,il,1,1),LPointMask(is,il,2,1),Xmax)
      Ymin=MIN(LPointMask(is,il,1,2),LPointMask(is,il,2,2),Ymin)
      Ymax=MAX(LPointMask(is,il,1,2),LPointMask(is,il,2,2),Ymax)
      Zmin=MIN(LPointMask(is,il,1,3),LPointMask(is,il,2,3),Zmin)
      Zmax=MAX(LPointMask(is,il,1,3),LPointMask(is,il,2,3),Zmax)
    END DO
    READ(InUnitMask,*) ThickMask(is),WidthMask(is),ch
    READ(InUnitMask,*) Ebaum(is,1),ch ! 0...1
    IF (Number<0) THEN
      Ebaum(is,1) = 0.d0
    ELSE IF (ThickMask(is)<=0.d0.OR.WidthMask(is)<=0.d0) THEN
      Ebaum(is,1) = 0.d0
    ELSE
      Ebaum(is,1) = MAX(Ebaum(is,1),0.d0)
    END IF
    WidthMaskM=MAX(  WidthMask(is)  ,WidthMaskM)
    ThickMaskM=MAX(  ThickMask(is)  ,ThickMaskM)
    MaskM=MAX(Ebaum(is,1),MaskM)
    DO il=1,nlMask(is)
      Ebaum(is,il)=Ebaum(is,1)
    END DO
  END DO
  CLOSE(UNIT=InUnitMask)

  WRITE(*,*) '    Trees:     ',nsMask
  WRITE(*,8) '    Expansion x:',Xmin,'...',Xmax
  WRITE(*,8) '              y:',Ymin,'...',Ymax
  WRITE(*,8) '              z:',Zmin,'...',Zmax
    IF (nr_wahlmask==nr_mask_baum_oro) THEN
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
  WRITE(*,9) '    Tree thick:    ',ThickMaskM
  WRITE(*,9) '    Tree width:    ',WidthMaskM
  WRITE(*,9) '    Tree density:  ',MaskM
9 FORMAT(1x,a19,f8.2)

END SUBROUTINE ReadMasking

SUBROUTINE AnalyzeMasking

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
  REAL(8) :: lvectorx,lvectory,lvectorz,wvectorx,wvectory,hvectorz
  REAL(8) :: epoint0x,epoint0y,epoint0z
  REAL(8) :: epointw(3),epointl(3),epointh(3),epointwv(3),epointlv(3)
  REAL(8) :: length1,width1,thick1,dlength,dwidth,dthick
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine AnalyzeMasking 
!------------------------------------------------------------------------------

  Mask =0.d0 ! Tree portion (density*dVol) per cell
  Mask0=0.d0 ! Tree density maximal possible in the cell
             ! (= Ebaum0, maximum of different trees in the cell)

  ALLOCATE(Ebaum0(nsMask,nlmMask))
  ALLOCATE(EPoint0Mask(nsMask,nlmMask,3),lstepMask(nsMask,nlmMask,3) &
            ,wstepMask(nsMask,nlmMask,2),hstepMask(nsMask,nlmMask))
  ALLOCATE(nlengthMask(nsMask,nlmMask),nwidthMask(nsMask,nlmMask),nthickMask(nsMask,nlmMask))

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

! 3D-Spreading of tree area to individual points
! (epoint = corner epoint0 + width step + length step + thick step)

  DO is=1,nsMask
    DO il=1,nlMask(is)
      lvectorx=LPointMask(is,il,2,1)-LPointMask(is,il,1,1)
      lvectory=LPointMask(is,il,2,2)-LPointMask(is,il,1,2)
      lvectorz=LPointMask(is,il,2,3)-LPointMask(is,il,1,3)

      length1=SQRT(lvectorx*lvectorx+lvectory*lvectory) ! linear extension
      width1 =WidthMask(is)
      thick1 =ThickMask(is)

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

      epoint0x=LPointMask(is,il,1,1)-wvectorx*width1*0.5d0 ! corner point
      epoint0y=LPointMask(is,il,1,2)-wvectory*width1*0.5d0
      epoint0z=LPointMask(is,il,1,3)

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

      epoint0x=epoint0x+lvectorx*0.5d0+wvectorx*0.5d0 ! first point
      epoint0y=epoint0y+lvectory*0.5d0+wvectory*0.5d0 ! (= epoint0 shifted by half steps)
      epoint0z=epoint0z+lvectorz*0.5d0+hvectorz*0.5d0

      EPoint0Mask(is,il,1)=epoint0x
      EPoint0Mask(is,il,2)=epoint0y
      EPoint0Mask(is,il,3)=epoint0z
      lstepMask  (is,il,1)=lvectorx
      lstepMask  (is,il,2)=lvectory
      lstepMask  (is,il,3)=lvectorz
      wstepMask  (is,il,1)=wvectorx
      wstepMask  (is,il,2)=wvectory
      hstepMask  (is,il)  =hvectorz
      nlengthMask(is,il)  =nlength1
      nwidthMask (is,il)  =nwidth1
      nthickMask (is,il)  =nthick1

      ! Up to here: Ebaum  = Density of a tree (0...1)
      ! Henceforth: Ebaum0 = Ebaum (maximal density possible in a cell)
      !             Ebaum  = Tree volume distributed on the street points

      Ebaum0(is,il)=Ebaum (is,il)
      Ebaum (is,il)=Ebaum0(is,il)*dlength*dwidth*dthick

!     Ebaum=0 for tree parts (is,il) definitely outside of domain:
!     (to shorten the following main loop)
      IF ((LPointMask(is,il,1,1)<=domain%x0.AND.LPointMask(is,il,2,1)<=domain%x0).OR. &
          (LPointMask(is,il,1,1)>=domain%x1.AND.LPointMask(is,il,2,1)>=domain%x1).OR. &
          (LPointMask(is,il,1,2)<=domain%y0.AND.LPointMask(is,il,2,2)<=domain%y0).OR. &
          (LPointMask(is,il,1,2)>=domain%y1.AND.LPointMask(is,il,2,2)>=domain%y1)) THEN
        Ebaum (is,il)=0.d0
        Ebaum0(is,il)=0.d0
      END IF
    END DO
  END DO

  DEALLOCATE(LPointMask)

! Filling the tree points in the grid cells

  DO is=1,nsMask
    DO il=1,nlMask(is)

      IF (Ebaum(is,il)==0.d0) GOTO 12

      iwidthact =1
      ilengthact=1
      ithickact =1
      nlength1  =nlengthMask(is,il)
      nwidth1   =nwidthMask (is,il)
      nthick1   =nthickMask (is,il)
      epointw(1)=EPoint0Mask(is,il,1)
      epointw(2)=EPoint0Mask(is,il,2)
      epointw(3)=EPoint0Mask(is,il,3)
      epointl(1)=EPoint0Mask(is,il,1)
      epointl(2)=EPoint0Mask(is,il,2)
      epointl(3)=EPoint0Mask(is,il,3)
      epointh(1)=EPoint0Mask(is,il,1)
      epointh(2)=EPoint0Mask(is,il,2)
      epointh(3)=EPoint0Mask(is,il,3)
      IF (nr_wahlmask==nr_mask_baum_oro) THEN ! + Orographie
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
      IF (wstepMask(is,il,1)>=0.d0) THEN
        idirw=1
        iendw=ix1
        iactw=ix0+1
      ELSE
        idirw=-1
        iendw=ix0+1
        iactw=ix1
      END IF
      IF (wstepMask(is,il,2)>=0.d0) THEN
        jdirw=1
        jendw=iy1
        jactw=iy0+1
      ELSE
        jdirw=-1
        jendw=iy0+1
        jactw=iy1
      END IF
!     Anfang,Ende,Richtung der ij-Loops bzgl. length steps (max. Bereich)
      IF (lstepMask(is,il,1)>=0.d0) THEN
        idirl=1
        iendl=ix1
        iactl=ix0+1
      ELSE
        idirl=-1
        iendl=ix0+1
        iactl=ix1
      END IF
      IF (lstepMask(is,il,2)>=0.d0) THEN
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
        epointw (1)=EPoint0Mask(is,il,1)+(iwidth-1)*wstepMask(is,il,1)
        epointw (2)=EPoint0Mask(is,il,2)+(iwidth-1)*wstepMask(is,il,2)
        epointw (3)=EPoint0Mask(is,il,3)
        IF (nr_wahlmask==nr_mask_baum_oro) THEN ! + Orographie
          epointw(3)=epointw(3)-foro(epointw(1),epointw(2),0.d0)
        END IF

        DO ilength=ilengthanf,nlength1 ! length-Loop ====================
          ilengthact=ilength

          epointlv(1)=epointl(1)
          epointlv(2)=epointl(2)
          epointlv(3)=epointl(3)
          epointl (1)=epointw(1)+(ilength-1)*lstepMask(is,il,1)
          epointl (2)=epointw(2)+(ilength-1)*lstepMask(is,il,2)
          epointl (3)=EPoint0Mask(is,il,3)+(ilength-1)*lstepMask(is,il,3)
          IF (nr_wahlmask==nr_mask_baum_oro) THEN ! + Orographie
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

            epointh(3)=epointl(3)+(ithick-1)*hstepMask(is,il)

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

            ! Tree volume
            Mask(iact,jact,kact)=Mask(iact,jact,kact)+Ebaum(is,il)

            ! Maximal moegliche Dichte (in einer Zelle) speichern
            Mask0(iact,jact,kact)=MAX(Mask0(iact,jact,kact),Ebaum0(is,il) &
                                      *MIN(ThickMask(is),dz(kact))/dz(kact))

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

  DEALLOCATE(EPoint0Mask,lstepMask,wstepMask,hstepMask,nlMask,Ebaum,Ebaum0)
  DEALLOCATE(nlengthMask,nwidthMask,nthickMask,WidthMask,ThickMask)

END SUBROUTINE AnalyzeMasking

SUBROUTINE WriteMasking(FileName)

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
  INTEGER :: OutUnitMask=10
  REAL(8) :: de0,de1 ,deTol=0.01d0,Eps=1.d-40 ! absolute Toleranz
!
! End of header
!==============================================================================

!------------------------------------------------------------------------------
! Begin Subroutine WriteMasking 
!------------------------------------------------------------------------------

  WRITE(*,*) '  Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=OutUnitMask,FILE=TRIM(FileName)//'.Mask',STATUS='unknown')
  WRITE(*,*) '  > ',TRIM(FileName),'.Mask'
! Mask here similar to Mask0: tree density (not volume)

  DO ib=1,nb
    CALL Set(Floor(ib))
    is=0
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          Mask0(i,j,k)=Mask0(i,j,k)*VolC(i,j,k)/(VolC(i,j,k)+Eps) ! max. possible density
          Mask (i,j,k)=Mask (i,j,k)*VolC(i,j,k)/(VolC(i,j,k)*VolC(i,j,k)+Eps) ! actual density
          IF (Mask(i,j,k)>deTol) THEN
            de0=Mask(i,j,k)-Mask0(i,j,k) ! Ueberschuss
            IF (de0>deTol) THEN
              WRITE(*,1) '  Baum loss:    ',ib,i,j,k,VolC(i,j,k),de0,Mask(i,j,k)
1             FORMAT(1x,a16,i4,3i5,4x,5x,f11.2,5x,f8.3,'  of',f8.3)
              Mask(i,j,k)=Mask0(i,j,k)
            END IF
            IF (Mask(i,j,k)>deTol) THEN ! numbering
              is=is+1
            END IF
          END IF
        END DO
      END DO
    END DO
    WRITE(OutUnitMask,*) ib,is,'!block,cells (MASK)' ! densities, blockwise
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        DO i=ix0+1,ix1
          IF (Mask(i,j,k)>deTol) THEN
            WRITE(OutUnitMask,*) i,j,k
            WRITE(OutUnitMask,*) MIN(Mask(i,j,k),1.d0) ! Tree density (0...1)
          END IF
        END DO
      END DO
    END DO
  END DO ! ib
  CLOSE(UNIT=OutUnitMask)

END SUBROUTINE WriteMasking

END MODULE Masking_Mod
