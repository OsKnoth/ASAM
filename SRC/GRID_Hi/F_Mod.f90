MODULE f_mod

  USE Geometry_Mod
  USE Haus_Mod
  USE Tree_Mod
  USE Domain_Mod

  IMPLICIT NONE

  REAL(8), PRIVATE :: H
  REAL(8), PRIVATE :: a
  REAL(8), PRIVATE :: aUP
  REAL(8), PRIVATE :: aDOWN
  REAL(8), PRIVATE :: Lambda
  REAL(8), PRIVATE :: xStart
  REAL(8), PRIVATE :: xEnd
  REAL(8), PRIVATE :: xL,xC,xR,hMin
  REAL(8), PRIVATE :: r1,r2
! Parameter Valley
  REAL(8), PRIVATE :: hP,Vx,Sx,Px,Sy,Py
! Parameter Inlet
  REAL(8), PRIVATE :: h1Inlet,h2Inlet,h3Inlet,h4Inlet,h5Inlet
  REAL(8), PRIVATE :: x1Inlet,x2Inlet,x3Inlet,x4Inlet,x5Inlet
  REAL(8), PRIVATE :: xOut1,xOut2,xOut3
  REAL(8), PRIVATE :: hOut1,hOut2
! Parameter BaroIn
  REAL(8), PRIVATE :: RotAngle


  REAL(8), PRIVATE :: b
  REAL(8), PRIVATE :: c
  REAL(8), PRIVATE :: d
  REAL(8), PRIVATE :: alpha
  REAL(8), PRIVATE :: beta
  REAL(8), PRIVATE :: nn
  REAL(8), PRIVATE :: mu
  REAL(8), PRIVATE :: phi

  INTEGER :: NumberHaus
  TYPE(Haus_T),POINTER :: Haus(:)
  TYPE(Box_T), POINTER :: Root
  INTEGER :: Depth=8

  CHARACTER(20) :: name_fkt
  CHARACTER(50) :: read_namefkt
  CHARACTER(50) :: file_namefkt
  INTEGER, PARAMETER :: anz_fkt=18
  INTEGER, PARAMETER :: nr_fkt_haus=10
  INTEGER, PARAMETER :: nr_fkt_gml=11
  INTEGER, PARAMETER :: nr_fkt_oro=12
  INTEGER, PARAMETER :: nr_fkt_kugel=13
  INTEGER, PARAMETER :: nr_fkt_sierra=14
  INTEGER, PARAMETER :: nr_fkt_surfit=15
  INTEGER, PARAMETER :: nr_fkt_wgs_oro=17
  INTEGER, PARAMETER :: nr_fkt_wgs_surf=18
  CHARACTER(20), DIMENSION(anz_fkt) :: w_namefkt=(/"Agnesi       ",&
                                                   "Agnesi3D     ",&
                                                   "Agnesi_2sided",&
                                                   "AgnesiDouble ",&
                                                   "Schaer       ",&
                                                   "Concave      ",&
                                                   "Bannon2D     ",&
                                                   "Bannon3D     ",&
                                                   "EigeneFunc   ",&
                                                   "Haus         ",&
                                                   "Haus_GML     ",&
                                                   "Oro          ",&
                                                   "Kugel        ",&
                                                   "Sierra       ",&
                                                   "SurfitProx   ",&
                                                   "HillFroehlich",&
                                                   "WGS84_Oro    ",&
                                                   "WGS84_Surf   "/)
  INTEGER, PARAMETER :: nr_fkt_gml_oro=anz_fkt+1 ! used externally and automatically

! #Emission:
  CHARACTER(20) :: name_emi
  CHARACTER(50) :: read_nameemi
  CHARACTER(50) :: file_nameemi
  INTEGER, PARAMETER :: anz_emi=1
  INTEGER, PARAMETER :: nr_emi_street=1
  CHARACTER(20), DIMENSION(anz_emi) :: w_nameemi=(/"Emi_Street   "/)
  INTEGER, PARAMETER :: nr_emi_street_oro=anz_emi+1 ! used externally and automatically

! #Baum:
  CHARACTER(20) :: name_mask
  CHARACTER(50) :: read_namemask
  CHARACTER(50) :: file_namemask
  INTEGER, PARAMETER :: anz_mask=1
  INTEGER, PARAMETER :: nr_mask_baum=1
  CHARACTER(20), DIMENSION(anz_mask) :: w_namemask=(/"Mask_Baum    "/)
  INTEGER, PARAMETER :: nr_mask_baum_oro=anz_mask+1 ! used externally and automatically

  !INTEGER       :: anz_parfkt(anz_fkt)=(/2,0,3,7/)
  INTEGER :: OutUnitCheck24=24
  !.......................................................................
  ! def 'Oro'
  REAL(8), ALLOCATABLE :: xPH(:),yPH(:),Height(:,:)
  REAL(8) :: x0H,y0H,x1H,y1H,dxH,dyH,maxHOro
  INTEGER :: nxH,nyH
  !........................................................................
  ! def 'WGS84_Oro', 'WGS84_Surf'
  REAL(8) :: grdy1,miny1,seky1
  REAL(8) :: grdy0,miny0,seky0
  REAL(8) :: grdx1,minx1,sekx1
  REAL(8) :: grdx0,minx0,sekx0
  INTEGER :: nxWgs,nyWgs
  CHARACTER    :: fout_wahl      ! Schnittstelle: out_wahlgrid -> C/G -> Cartesien,Globe 
  CHARACTER(9) :: fconv_ingrid   ! Schnittstelle: conv_ingrid
  CHARACTER(5) :: findata_type   ! Schnittstelle: indata_type 
  ! def 'Corine'
  !GRID_CODE,CLC_CODE,LABEL1,LABEL2,LABEL3,RGB
  CHARACTER(8)  :: str_corine
  CHARACTER(8)  :: logical_corine
  CHARACTER(50) :: file_ncorine
  CHARACTER(25)  :: file_clc_list="CLC_IfT.csv"    !"CLC.csv"  ! Def.: CLC1990
  LOGICAL :: ifcorine=.FALSE.
  INTEGER, PARAMETER :: Nr_Land_Cl=44 ! Anz.def.LandNutzungsKlassen, CLC1990
  INTEGER, ALLOCATABLE :: LandCover(:,:)
  TYPE rgb_T
    INTEGER :: nr_r
    INTEGER :: nr_g
    INTEGER :: nr_b
    INTEGER :: mat_rgb
  END TYPE
  TYPE Corine_T
    INTEGER :: grid_c
    INTEGER :: clc_c
    INTEGER :: mat_c          ! Erweitert für File=CLC_IfT.csv
    CHARACTER(32) :: lab1
    CHARACTER(50) :: lab2
    CHARACTER(90) :: lab3
    CHARACTER(11) :: lrgb
  END TYPE
  TYPE (Corine_T), POINTER :: v_clc(:)    ! (1:Nr_Land_Cl)
  TYPE (rgb_T),    POINTER :: v_rgb(:)    ! (1:Nr_Land_Cl)
  !........................................................................
  !def 'Kugel'
  REAL(8) :: phi0,lam0
  CHARACTER*20 :: NameFun
  REAL(8) :: RIn,ROut
  !........................................................................
  ! for sierra (approximation with curfit,splev)
  REAL(8), DIMENSION(:), ALLOCATABLE :: x, y, w      ! with 'm' allocate
  REAL(8), DIMENSION(:), ALLOCATABLE :: wrk          ! with 'lwrk' allocate
  !REAL(8), DIMENSION(:), ALLOCATABLE :: sp           ! with 'm' allocate
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwrk         ! with 'nest' allocate
  REAL(8), DIMENSION (:), ALLOCATABLE :: t, b_coef   ! with 'nest' allocate
  INTEGER :: n  ! total number of knots
  INTEGER :: k  ! degree of s(x) smoothing
  REAL(8) :: s  ! specify the smoothing factor
  INTEGER :: iopt  ! must specify whether a weighted
  !.........................................................................
  ! for        (approximation with surfit,bispev)a
  REAL(8), DIMENSION(:), ALLOCATABLE :: z    !,x, y, w
  REAL(8), DIMENSION(:), ALLOCATABLE :: tx,ty,c_spl,wrk1,wrk2  !,iwrk
  INTEGER :: kx,ky,tx_nx,ty_ny
  INTEGER :: kwrk,lwrk2
  CHARACTER :: out_surf,conv_gk
  REAL(8) :: gaussx_min=0.0,gaussy_min=0.0
  REAL(8) :: diffgx=0.0,diffgy=0.0  ! Differenz: Gauss-grid-min(x,y) and *.dat-min(x,y) 
             ! !!diffgx,diffgy also summand at cart-parametric-calculations !!! 
  !.........................................................................

CONTAINS


SUBROUTINE Read_Val_Fkt    ! für  Agnesi
    READ(10,*) H
    READ(10,*) aUP   ! Halbwertsbreite vor Berg
    READ(10,*) aDOWN ! Halbwertsbreite nach Berg
END SUBROUTINE Read_Val_Fkt

SUBROUTINE Read_Val_Fkt2    ! für  Agnesi 3D
    READ(10,*) H
    READ(10,*) a
    READ(10,*) b
END SUBROUTINE Read_Val_Fkt2

SUBROUTINE Read_Val_Fkt3    ! für  Agnesi 2sided 2D
    READ(10,*) H
    READ(10,*) aUP
    READ(10,*) aDOWN
END SUBROUTINE Read_Val_Fkt3

SUBROUTINE Read_Val_Fkt4    ! für  Agnesi Double
    READ(10,*) H
    READ(10,*) a
    READ(10,*) b
    READ(10,*) d  !!!Abstand Bergmaxima
END SUBROUTINE Read_Val_Fkt4

SUBROUTINE Read_Val_Fkt5   ! für  Schaer
  READ(10,*) H
  READ(10,*) a
  READ(10,*) Lambda
END SUBROUTINE Read_Val_Fkt5

SUBROUTINE Read_Val_Fkt6   ! für ! Concave
  READ(10,*) H
  READ(10,*) a
  READ(10,*) b
  READ(10,*) c
  READ(10,*) d
  READ(10,*) alpha
  READ(10,*) beta
END SUBROUTINE Read_Val_Fkt6

SUBROUTINE Read_Val_Fkt7   ! für Bannon 2D
  READ(10,*) H
  READ(10,*) a
  READ(10,*) nn
  READ(10,*) mu
  READ(10,*) phi
END SUBROUTINE Read_Val_Fkt7

SUBROUTINE Read_Val_Fkt8   ! für Bannon 3D
  READ(10,*) H
  READ(10,*) a
  READ(10,*) b
  READ(10,*) nn
  READ(10,*) mu
  READ(10,*) phi
END SUBROUTINE Read_Val_Fkt8

SUBROUTINE Read_Val_Fkt9   ! linearer Anstieg
  REAL(8) :: A,B,C,D
  REAL(8) :: x4Inlet1
  REAL(8) :: x4Inlet2

  READ(10,*) NameFun
  SELECT CASE (NameFun)
    CASE ('Annulus')
      READ(10,*) ROut
      READ(10,*) RIn
      READ(10,*) H
    CASE ('Kegel')
      READ(10,*) H
      READ(10,*) xL
      READ(10,*) xC
      READ(10,*) xR
      READ(10,*) hMin
    CASE ('Leer')
      READ(10,*) r1
      READ(10,*) r2
    CASE ('ValleyTwo')
      WRITE(*,*) 'ValleyTwo'
      READ(10,*) hP
      READ(10,*) Vx
      READ(10,*) Sx
      READ(10,*) Px
      READ(10,*) Sy
      READ(10,*) Py
    CASE ('Valley3D')
      WRITE(*,*) 'Valley3D'
      READ(10,*) hP
      READ(10,*) Vx
      READ(10,*) Sx
      READ(10,*) Sy
    CASE ('Valley2D')
      WRITE(*,*) 'Valley2D'
      READ(10,*) hP
      READ(10,*) Sy
    CASE ('Inlet')
      WRITE(*,*) 'Inlet'
      READ(10,*) h1Inlet
      READ(10,*) h2Inlet
      READ(10,*) h3Inlet
      READ(10,*) h4Inlet
      READ(10,*) h5Inlet
      READ(10,*) x1Inlet
      READ(10,*) x2Inlet
      READ(10,*) x3Inlet
      READ(10,*) x4Inlet
      READ(10,*) x5Inlet
      READ(10,*) xOut1
      READ(10,*) xOut2
      READ(10,*) xOut3
      READ(10,*) hOut1
      READ(10,*) hOut2
!       left derivative
!       f9Prime=-2.0d0*(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x4Inlet-x3Inlet)
!       right derivative
!       f9Prime=-h4Inlet/(x5Inlet-x4Inlet)
!       A*(x4-B)=C/(D-x4)
!       (x4-B)*(D-x4)=C/A
!       x4^2-(B+D)*x4+BD+C/A=0
!       x4=(B+D)/2+-SQRT(1/4*(B*D)^2-BD-C/A)
      A=-2.0d0*(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0
      B=x3Inlet
      C=-h4Inlet
      D=x5Inlet
      x4Inlet=0.5d0*(B+D)-SQRT(0.25d0*(B+D)**2.0d0-B*D-C/A)
    CASE ('BaroIn')
      READ(10,*) RotAngle
      RotAngle=RotAngle/180.0d0*Pi
  END SELECT
END SUBROUTINE Read_Val_Fkt9


SUBROUTINE  Read_Val_Haus
  INTEGER:: h,i,k,nr_p,nr_h,Number
  TYPE(Point_T) ::Pdis
  REAL(8):: disp
  REAL(8):: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
  INTEGER:: InUnitV2=7

  Xmin= 1.d20
  Xmax=-1.d20
  Ymin= 1.d20
  Ymax=-1.d20
  Zmin= 1.d20
  Zmax=-1.d20

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')

  READ(InUnitV2,*) NumberHaus
  ALLOCATE(Haus(1:NumberHaus))
  nr_h=0 ! number of erroneous objects (not analyzed)
  WRITE(*,*) 'NumberHaus',NumberHaus
  DO h=1,NumberHaus
    WRITE(*,*) 'h ',h
    READ(InUnitV2,*) Number,Haus(h)%NumberOfPoints
    ALLOCATE(Haus(h)%Points(Haus(h)%NumberOfPoints))
    DO i=1,Haus(h)%NumberOfPoints
      READ(InUnitV2,*) Haus(h)%Points(i)%x,Haus(h)%Points(i)%y,Haus(h)%Points(i)%z
      Xmin=MIN(Haus(h)%Points(i)%x,Xmin)
      Xmax=MAX(Haus(h)%Points(i)%x,Xmax)
      Ymin=MIN(Haus(h)%Points(i)%y,Ymin)
      Ymax=MAX(Haus(h)%Points(i)%y,Ymax)
      Zmin=MIN(Haus(h)%Points(i)%z,Zmin)
      Zmax=MAX(Haus(h)%Points(i)%z,Zmax)
    END DO
    READ(InUnitV2,*) Haus(h)%NumberOfFaces
    ALLOCATE(Haus(h)%Faces(Haus(h)%NumberOfFaces))
    READ(InUnitV2,*)(Haus(h)%Faces(i)%NumberOfPoints, i=1,Haus(h)%NumberOfFaces)
    DO i=1,Haus(h)%NumberOfFaces
       nr_p=Haus(h)%Faces(i)%NumberOfPoints
       ALLOCATE(Haus(h)%Faces(i)%ListOfPoints(nr_p))
       READ(InUnitV2,*) Haus(h)%Faces(i)%type &
                       ,(Haus(h)%Faces(i)%ListOfPoints(k), k=1,nr_p)
    END DO   
    !Anpassung neue Syntax(Read_Val_Haus_GML)
    DO i=1,Haus(h)%NumberOfFaces
      nr_p=Haus(h)%Faces(i)%NumberOfPoints
      ALLOCATE(Haus(h)%Faces(i)%Points(nr_p))
      DO k=1,nr_p
        Haus(h)%Faces(i)%Points(k)=Haus(h)%Points(Haus(h)%Faces(i)%ListOfPoints(k))
      END DO
    END DO
    !Ende Anpassung
    CALL NormalForm(Haus(h))
    IF (Number<0) Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z=-1.d0
    IF (Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z<=0.d0) nr_h=nr_h+1
    ! Hier noch erforderlich: Erstellung der Haus-Liste für DistAll (see Read_Val_Haus_GML)
  END DO
  CLOSE(UNIT=InUnitV2)

  WRITE(*,*) '    Objects:   ',NumberHaus
  WRITE(*,*) '      faulty:  ',nr_h ,' (not analyzed)'
  WRITE(*,9) '    Expansion x:',Xmin,'...',Xmax
  WRITE(*,9) '              y:',Ymin,'...',Ymax
  WRITE(*,9) '              z:',Zmin,'...',Zmax
9 FORMAT(1x,a16,f11.2,8x,a3,7x,f11.2)

  IF (     Xmin>=Domain%x0.AND.Xmax<=Domain%x1.AND. &
           Ymin>=Domain%y0.AND.Ymax<=Domain%y1.AND. &
           Zmin>=Domain%z0.AND.Zmax<=Domain%z1) THEN ! totally inside
    WRITE(*,*) '      inside of domain'
  ELSE IF (Xmin>=Domain%x1.OR. Xmax<=Domain%x0.OR. &
           Ymin>=Domain%y1.OR. Ymax<=Domain%y0) THEN ! beyond xy-domain
    STOP '***** Stop: beyond xy-domain *****'
  ELSE IF (Zmin>=Domain%z1.OR. Zmax<=Domain%z0) THEN ! beyond z-domain
    STOP '***** Stop: beyond z-domain *****'
  ELSE IF (Zmin>=Domain%z0.AND.Zmax<=Domain%z1) THEN ! partially inside z-domain
    WRITE(*,*) '      partially outside of xy-domain'
  ELSE
    WRITE(*,*) '      partially outside of z-domain'
  END IF

END SUBROUTINE Read_Val_Haus


SUBROUTINE  Read_Val_Haus_GML
  INTEGER:: h,i,j,k,l,nr_p,nr_h,nr_h1,Number
  TYPE(Point_T) ::Pdis
  REAL(8):: disp,disp1
  REAL(8):: Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,ZminOro,ZmaxOro
  REAL(8):: Xmin1,Xmax1,Ymin1,Ymax1,Zmin1,Zmax1 ! Individual houses
  INTEGER:: InUnitV2=7
  !Read-Syntax:
  !numberHaus
  !numberFlächen
  !numberEckenJeFläche ..
  !array's x- y- z-KoordinatenFläche

  Xmin= 1.d20
  Xmax=-1.d20
  Ymin= 1.d20
  Ymax=-1.d20
  Zmin= 1.d20
  Zmax=-1.d20
  ZminOro= 1.d20
  ZmaxOro=-1.d20

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')

  READ(InUnitV2,*) NumberHaus
  ALLOCATE(Haus(1:NumberHaus))
  nr_h =0 ! number of faulty  objects (not analyzed)
  nr_h1=0 ! number of outside objects (not analyzed)
  DO h=1,NumberHaus
    READ(InUnitV2,*) Number,Haus(h)%NumberOfFaces
    ALLOCATE(Haus(h)%Faces(Haus(h)%NumberOfFaces))
    READ(InUnitV2,*) (Haus(h)%Faces(i)%NumberOfPoints, i=1,Haus(h)%NumberOfFaces)
    READ(InUnitV2,*) (Haus(h)%Faces(i)%type, i=1,Haus(h)%NumberOfFaces)
    Xmin1= 1.d20
    Xmax1=-1.d20
    Ymin1= 1.d20
    Ymax1=-1.d20
    Zmin1= 1.d20
    Zmax1=-1.d20
    DO i=1,Haus(h)%NumberOfFaces
      ALLOCATE(Haus(h)%Faces(i)%Points(Haus(h)%Faces(i)%NumberOfPoints))
      READ(InUnitV2,*)  (Haus(h)%Faces(i)%Points(j)%x, j=1,Haus(h)%Faces(i)%NumberOfPoints) &
                       ,(Haus(h)%Faces(i)%Points(k)%y, k=1,Haus(h)%Faces(i)%NumberOfPoints) & 
                       ,(Haus(h)%Faces(i)%Points(l)%z, l=1,Haus(h)%Faces(i)%NumberOfPoints)
      DO j=1,Haus(h)%Faces(i)%NumberOfPoints
        Xmin1=MIN(Haus(h)%Faces(i)%Points(j)%x,Xmin1)
        Xmax1=MAX(Haus(h)%Faces(i)%Points(j)%x,Xmax1)
        Ymin1=MIN(Haus(h)%Faces(i)%Points(j)%y,Ymin1)
        Ymax1=MAX(Haus(h)%Faces(i)%Points(j)%y,Ymax1)
        Zmin1=MIN(Haus(h)%Faces(i)%Points(j)%z,Zmin1)
        Zmax1=MAX(Haus(h)%Faces(i)%Points(j)%z,Zmax1)
      END DO
    END DO

    IF (nr_wahlfkt==nr_fkt_gml_oro) THEN
!     Addition der Orografiewerte (Mittelwert)
      disp=0.d0
      disp1=0.d0
      DO i=1,Haus(h)%NumberOfFaces
        DO j=1,Haus(h)%Faces(i)%NumberOfPoints
          disp=disp-fOro(Haus(h)%Faces(i)%Points(j)%x &
                        ,Haus(h)%Faces(i)%Points(j)%y &
                        ,0.d0)
          disp1=disp1+1.d0
        END DO
      END DO
      disp=disp/disp1
      DO i=1,Haus(h)%NumberOfFaces
        DO j=1,Haus(h)%Faces(i)%NumberOfPoints
          Haus(h)%Faces(i)%Points(j)%z=Haus(h)%Faces(i)%Points(j)%z + disp
          ZminOro=MIN(Haus(h)%Faces(i)%Points(j)%z,ZminOro)
          ZmaxOro=MAX(Haus(h)%Faces(i)%Points(j)%z,ZmaxOro)
        END DO
      END DO
    ELSE
      ZminOro=Zmin
      ZmaxOro=Zmax
    END IF

    CALL NormalForm(Haus(h))
    IF (Number<0) Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z=-1.d0
    IF (Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z==-1.d0) nr_h =nr_h +1 ! faulty

    CALL BoundingBox(Haus(h)) 
    IF (Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z==-2.d0) nr_h1=nr_h1+1 ! outside

    IF (Haus(h)%Faces(Haus(h)%NumberOfFaces)%Normal%z>0.d0) THEN
      Xmin=MIN(Xmin1,Xmin)
      Xmax=MAX(Xmax1,Xmax)
      Ymin=MIN(Ymin1,Ymin)
      Ymax=MAX(Ymax1,Ymax)
      Zmin=MIN(Zmin1,Zmin)
      Zmax=MAX(Zmax1,Zmax)
    END IF
  END DO
  CLOSE(UNIT=InUnitV2)
  ALLOCATE(Root)
  NULLIFY(Root%List)
  Root%P0%x=Domain%x0
  Root%P0%y=Domain%y0
  Root%P0%z=Domain%z0
  Root%P1%x=Domain%x1
  Root%P1%y=Domain%y1
  Root%P1%z=Domain%z1
  Root%TypeCut='x'
  Root%Index=0
  CALL CreateTree(Root,Depth)
 
  DO i=1,NumberHaus
    IF (Haus(i)%Faces(Haus(i)%NumberOfFaces)%Normal%z>0.d0) THEN
      CALL InsertHaus(Root,Haus(i),i)
    END IF
  END DO

  WRITE(*,*) '    Objects:   ',NumberHaus
  WRITE(*,*) '      faulty:  ',nr_h ,' (not analyzed)'
  WRITE(*,*) '      outside: ',nr_h1,' (not analyzed)'
  WRITE(*,9) '    Expansion x:',Xmin,'...',Xmax
  WRITE(*,9) '              y:',Ymin,'...',Ymax
  WRITE(*,9) '              z:',Zmin,'...',Zmax
    IF (nr_wahlfkt==nr_fkt_gml_oro) THEN
  WRITE(*,9) '      incl. ORO:',ZminOro,'...',ZmaxOro
    END IF
9 FORMAT(1x,a16,f11.2,8x,a3,7x,f11.2)

  IF (     Xmin   >=Domain%x0.AND.Xmax   <=Domain%x1.AND. &
           Ymin   >=Domain%y0.AND.Ymax   <=Domain%y1.AND. &
           ZminOro>=Domain%z0.AND.ZmaxOro<=Domain%z1) THEN ! totally inside
    WRITE(*,*) '      inside of domain'
  ELSE IF (Xmin   >=Domain%x1.OR. Xmax   <=Domain%x0.OR. &
           Ymin   >=Domain%y1.OR. Ymax   <=Domain%y0) THEN ! beyond xy-domain
    STOP '***** Stop: beyond xy-domain *****'
  ELSE IF (ZminOro>=Domain%z1.OR. ZmaxOro<=Domain%z0) THEN ! beyond z-domain
    STOP '***** Stop: incl. ORO beyond z-domain *****'
  ELSE IF (ZminOro>=Domain%z0.AND.ZmaxOro<=Domain%z1) THEN ! partially inside z-domain
    WRITE(*,*) '      partially outside of xy-domain'
  ELSE
    WRITE(*,*) '      partially outside of z-domain (incl. ORO)'
  END IF

END SUBROUTINE Read_Val_Haus_GML


SUBROUTINE Read_Oro

  INTEGER :: i,j
  REAL(8) :: Dummy
  REAL(8) :: Hmin,Hmax
  INTEGER:: InUnitOro=8

  Hmin= 1.d20
  Hmax=-1.d20
  nxH=0

  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InUnitOro,FILE=TRIM(file_namefkt),STATUS='OLD')
  READ(InUnitOro,*) nxH
  READ(InUnitOro,*) nyH
  READ(InUnitOro,*) x0H
  READ(InUnitOro,*) y0H
  READ(InUnitOro,*) dxH
  READ(InUnitOro,*) dyH
  ALLOCATE(xPH(nxH))
  ALLOCATE(yPH(nyH))
  ALLOCATE(Height(nxH,nyH))
  DO i=1,nxH
    DO j=1,nyH
      READ(InUnitOro,*) Dummy,Dummy,Height(i,j)
      Hmin=MIN(Height(i,j),Hmin)
      Hmax=MAX(Height(i,j),Hmax)
    END DO
  END DO
  xPH(1)=x0H
  DO i=2,nxH
    xPH(i)=xPH(i-1)+dxH
  END DO
  yPH(1)=y0H
  DO j=2,nyH
    yPH(j)=yPH(j-1)+dyH
  END DO

  WRITE(*,9) '    Expansion x:',xPH(1),'...',xPH(nxH)
  WRITE(*,9) '              y:',yPH(1),'...',yPH(nyH)
  WRITE(*,9) '              z:',Hmin,  '...',Hmax
9 FORMAT(1x,a16,f11.2,8x,a3,7x,f11.2)

  IF (     xPH(1)>=Domain%x1.OR.xPH(nxH)<=Domain%x0.OR. &
           yPH(1)>=Domain%y1.OR.yPH(nyH)<=Domain%y0) THEN
    STOP '***** Stop: OROGRAPHY FIELD beyond xy-domain *****'
  ELSE IF (xPH(1)> Domain%x0.OR.xPH(nxH)< Domain%x1.OR. &
           yPH(1)> Domain%y0.OR.yPH(nyH)< Domain%y1) THEN
    STOP '***** Stop: OROGRAPHY FIELD not capturing whole xy-domain *****'
  ELSE IF (Hmin  >=Domain%z1.OR.Hmax    < Domain%z0) THEN
    STOP '***** Stop: OROGRAPHY FIELD beyond z-domain *****'
  ELSE IF (Hmin  < Domain%z0.OR.Hmax    > Domain%z1) THEN
    STOP '***** Stop: OROGRAPHY FIELD partially outside of z-domain'
  ELSE
    WRITE(*,*) '      capturing whole xy-domain'
  END IF
  CLOSE(UNIT=InUnitOro)
  RETURN

1 WRITE(*,*) '  Orography=Zero by default'

  CLOSE(UNIT=InUnitOro)
END SUBROUTINE Read_Oro


SUBROUTINE Read_Val_FktKugel    ! for mountain
  REAL(8) :: Pi

  Pi=ATAN(1.0d0)*4.0d0

  READ(10,*) H
  READ(10,*) d
  READ(10,*) lam0
  READ(10,*) phi0
  READ(10,*) NameFun
  phi0=phi0/180.0d0*Pi
  lam0=lam0/180.0d0*Pi

END SUBROUTINE Read_Val_FktKugel


SUBROUTINE  Read_Val_Sierra      ! Sierra
    READ(10,*) iopt
    READ(10,*) k
    READ(10,*) s 
   
    CALL Read_Val_FktSierra 

END SUBROUTINE Read_Val_Sierra

SUBROUTINE  Read_Val_FktSierra   ! Sierra
  INTEGER :: ix
  REAL(8)    :: ai,fp,xb,xe
  INTEGER :: i,ier,is,j,l,lwrk,l1,l2,m,nest,nk1
  INTEGER :: InUnitV2=7

  READ(10,*) iopt   ! specify whether a weighted spline variation
  READ(10,*) k      ! degree of the spline
  READ(10,*) s      ! smoothing factor
 
  IF(ALLOCATED(b_coef)) THEN
    RETURN
  END IF

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')

! m denotes the number of data points
      read(InUnitV2,*) m
      ALLOCATE(x(1:m))
      ALLOCATE(y(1:m))
      ALLOCATE(w(1:m))
!      ALLOCATE(sp(1:m))

! we fetch the co-ordinate and function values of each data point.
      DO  i=1,m
        read(InUnitV2,*) ix,y(i)
        x(i)=REAL(ix)
      END DO 

! we set up the abscissae and weights of the data points
      DO  i=1,m
        w(i) = 1.0
      END DO 

! we set up the boundaries of the approximation interval
      xb = x(1)
      xe = x(m)

! we set up the dimension information.
      !k=3    ! degree of the spline, als read im grid-file definieren
      nest = m+k+1
      lwrk = m*(k+1)+nest*(7+3*k)
      !Berechnung fuer k=5
      !nest = 450  !407 berechnet ->aufgerundet
      !lwrk = 11400 !11360 berechnet ->aufgerundet
      ALLOCATE(t(1:nest))
      ALLOCATE(b_coef(1:nest))
      ALLOCATE(iwrk(1:nest))
      ALLOCATE(wrk(1:lwrk))
! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1) or a smoothing spline (iopt=
! 0 or 1) must be determined.
      !s=30.   !read
      !iopt=0  !read
      IF(iopt==-1) THEN
        j=k+2
        DO l=1,7
          ai=3*l
          t(j)=ai
          j=j+1
        END DO
        n=9+2*k
      END IF

! spline approximations of degree k
      call curfit(iopt,m,x,y,w,xb,xe,k,s,nest,n,t,b_coef,fp,wrk,lwrk,iwrk,ier)

! evaluation of the spline approximation
!      call splev(t,n,b_coef,k,x,sp,m,ier)

   CLOSE(UNIT=InUnitV2)
END SUBROUTINE Read_Val_FktSierra


SUBROUTINE  Read_Val_FktSurfit   !  dasurf.htm

  REAL(8) :: xx(11),yy(11),zz(121) !Bsp:
  REAL(8) :: x_datmin,y_datmin,min_px,min_py
  INTEGER :: i,ier,iopt,is,j,lwrk1,m,mx,my,nc
  INTEGER :: nmax,nxest,nyest
  INTEGER :: u,v,km,ne,bx,by,b1,b2
  REAL(8) :: ai,delta,eps,fp,s,ww,xb,xe,yb,ye
  REAL(8) :: m_half,xmin,xmax,ymin,ymax,breite,laenge
  REAL(8) :: diffmx,diffmy
  INTEGER :: InUnitV2=7

  IF(ALLOCATED(c_spl)) THEN
    RETURN
  END IF

  READ(10,*) iopt   ! specify whether a weighted spline variation
  READ(10,*) kx     ! degree of the spline x-direction
  READ(10,*) ky     ! degree of the spline y-direction
  READ(10,*) s      ! smoothing factor

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')
  !open(6,file='surf.out',status='UNKNOWN')

! m denotes the number of data points
! m>=(kx+1)*(ky+1) !!!Wichtig!!!
      read(InUnitV2,*) m
      ALLOCATE(x(1:m))
      ALLOCATE(y(1:m))
      ALLOCATE(z(1:m))
      ALLOCATE(w(1:m))

! fetch the co-ordinate and function values of each data point.
      DO  i=1,m
        read(InUnitV2,*) y(i),x(i),z(i)
      END DO 
      x_datmin=x(1)
      y_datmin=y(m)
 
      IF(TRIM(file_namefkt).NE. 'dasurf.htm') THEN
        IF(conv_gk=='s'.OR.out_surf=='G') THEN
          DO  i=1,m
            Call OGKGEO(y(i)/1000, x(i)/1000, breite, laenge)
            y(i)=(breite*4.0d0*ATAN(1.0d0))/180
            x(i)=(laenge*4.0d0*ATAN(1.0d0))/180
          END DO
        ELSE 
          !!Betreff Nullpunkt ermitteln 
          !!Wenn kleinste Koordinate aus *.dat beachtet werden soll
          !!z.Zt nur berechnet aber nicht verwendet
          !!gauss[x,y]_min aus grid-Domain nur Bezug, 
          !!Nullpunkt wird erst beim OutputGMV verschoben
          IF(gaussx_min<x_datmin) THEN
            min_px=gaussx_min
            diffgx=0.0
          ELSE
            min_px=x_datmin
            diffgx=gaussx_min-x_datmin
          END IF
          IF(gaussy_min<y_datmin) THEN
            min_py=gaussy_min
            diffgy=0.0
          ELSE
            min_py=y_datmin
            diffgy=gaussy_min-y_datmin
          END IF
          !DO  i=1,m
            !!Nullpunkt verschieben
            !  x(i)=x(i)-min_px
            !  y(i)=y(i)-min_py  
          !END DO
        END IF
      END IF ! nicht Input dasurf.htm

 
! fetch an estimate of the standard deviation of the data values.
      read(InUnitV2,*) delta

! the weights are set equal to delta**(-1)
      ww = 1./delta
      DO  i=1,m
        w(i) = ww 
      END DO 

! set up the boundaries of the approximation domain.
      !xb = -2.   !bedingt Bsp.: dasurf.htm
      !xe = 2.    ! ''
      !yb = -2.    ! ''
      !ye = 2.     ! ''
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      DO i=2,m
        xmin=MIN(xmin,x(i))
        xmax=MAX(xmax,x(i))
        ymin=MIN(ymin,y(i))
        ymax=MAX(ymax,y(i))
      END DO 
      diffmx=xmax-xmin
      diffmy=ymax-ymin
      IF(conv_gk=='s'.OR.out_surf=='G') THEN
        xb=xmin-0.01*diffmx
        xe=xmax+0.01*diffmx
        yb=ymin-0.01*diffmy
        ye=ymax+0.01*diffmy
      ELSE
        xb=REAL(FLOOR(xmin))   ! nicht fuer Bogenmass
        xe=REAL(CEILING(xmax)) ! nicht fuer Bogenmass
        yb=REAL(FLOOR(ymin))   ! nicht fuer Bogenmass
        ye=REAL(CEILING(ymax)) ! nicht fuer Bogenmass
      END IF
 
! generate a rectangular grid for evaluating the splines.
!      mx = 11
!      my = 11
!      DO i=1,11
!        ai = i-6
!        xx(i) = ai*0.4
!        yy(i) = xx(i)
!      END DO 


! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1) or a smoothing spline (iopt=
! 0 or 1) must be determined.
     ! iopt = 0    ! read
     ! kx = 3      ! read
     ! ky = 3      ! read
     ! s = 900000.   !s= 30 !read


! set up the dimension information.
!    -must specify an upper bound for the number of knots required
!     in the x- and y-directions respect.
!    -nxest >= 2*(kx+1), nyest >= 2*(ky+1)
!    -in most practical situation 
!     nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient.
      m_half=m/2
      nxest = kx+1+CEILING(SQRT(m_half)) 
      nyest = ky+1+CEILING(SQRT(m_half))
      IF (nxest<(2*kx+2)) THEN
        nxest=2*kx+2
      END IF
      IF (nyest<(2*ky+2)) THEN
        nyest=2*ky+2
      END IF
      nmax = MAX(nxest,nyest)
      ALLOCATE(tx(1:nmax))
      ALLOCATE(ty(1:nmax))
      ALLOCATE(c_spl((nxest-kx-1)*(nyest-ky-1)))
      kwrk = m+(nxest-2*kx-1)*(nyest-2*ky-1)
      ALLOCATE(iwrk(1:kwrk))
      !computation for 'lwrk1'-> must specify the actual dimension of wrk1
      u = nxest-kx-1
      v = nyest-ky-1
      km = max(kx,ky)+1
      ne = nmax
      bx = kx*v+ky+1
      by = ky*u+kx+1
      if(bx.le.by) THEN
           b1 = bx
           b2 = b1+v-ky
      END IF
      if(bx.gt.by) THEN
           b1 = by
           b2 = b1+u-kx
      END IF
      lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      ALLOCATE(wrk1(1:lwrk1))
      lwrk2=u*v*(b2+1)+b2 
      ALLOCATE(wrk2(1:lwrk2))
      !im Bsp surfit def. : nxest = 15,nyest = 15,nmax = 15
      !                     !kwrk = 300, lwrk1 = 12000, lwrk2 = 6000
 
! choose a value for eps
      eps=0.1e-05

! integer flag. on entry iopt must specify whether a weighted
! least-squares spline (iopt=-1)
 
      IF(iopt==-1) THEN
        !kx = 3
        !ky = 3
        tx_nx=11
        ty_ny=11
        j=kx+2
        DO i=1,3
          ai = i-2
          tx(j)=ai
          ty(j)=ai
          j=j+1
        END DO
      ELSE
        tx_nx=nmax
        ty_ny=nmax
      END IF

! spline approximations of degree k
      call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
        nmax,eps,tx_nx,tx,ty_ny,ty,c_spl,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

! evaluation of the spline approximation
!      call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xx,mx,yy,my,zz,   &
!         wrk2,lwrk2,iwrk,kwrk,ier)

    CLOSE(UNIT=InUnitV2)
END SUBROUTINE Read_Val_FktSurfit


SUBROUTINE  Read_Hill_Froeh

   Write(*,*) "Unnecessary Input, Simulation 'Flow over periodic hill',"
   Write(*,*) "                          (Jochen Froelich)"

END SUBROUTINE  Read_Hill_Froeh

SUBROUTINE Read_WGS84_Koord !old, direkte geo-Angabe 
  INTEGER       :: InUnitOro=8
  CHARACTER(20) :: def_str

  OPEN(UNIT=InUnitOro,FILE=TRIM(file_namefkt),STATUS='OLD')
  READ(InUnitOro,*) def_str,grdy1,miny1,seky1   
  READ(InUnitOro,*) def_str,grdy0,miny0,seky0
  READ(InUnitOro,*) def_str,grdx1,minx1,sekx1
  READ(InUnitOro,*) def_str,grdx0,minx0,sekx0 
  READ(InUnitOro,*) def_str,nyWgs
  READ(InUnitOro,*) def_str,nxWgs

  Domain%x0=grdx0+minx0/60+sekx0/3600
  Domain%x1=grdx1+minx1/60+sekx1/3600
  Domain%y0=grdy0+miny0/60+seky0/3600
  Domain%y1=grdy1+miny1/60+seky1/3600
  Domain%nx=nxWgs-1  !nxWGS = number value direction
  Domain%ny=nyWgs-1  !nyWGS = number value direction
  CLOSE(UNIT=InUnitOro)
END SUBROUTINE Read_WGS84_Koord

SUBROUTINE Read_WGS84_Oro

  INTEGER   :: iy,ix,i,j
  REAL(8)   :: Dummy
  REAL(8)   :: grdy1,miny1,seky1,gy1H,gk_y1
  REAL(8)   :: grdy0,miny0,seky0,gy0H,gk_y0
  REAL(8)   :: grdx1,minx1,sekx1,gx1H,gk_x1
  REAL(8)   :: grdx0,minx0,sekx0,gx0H,gk_x0
  REAL(8)   :: breite,laenge
  INTEGER   :: rows,cols
  CHARACTER (LEN=20)  :: def_str,l_str
  CHARACTER (LEN=512) :: read_line
  INTEGER   :: InUnitOro=8

  IF (ALLOCATED(Height)) THEN
    RETURN
  END IF
  OPEN(UNIT=InUnitOro,FILE=TRIM(file_namefkt),STATUS='OLD')
  IF(TRIM(findata_type)=="gk") THEN
     READ(InUnitOro,*) def_str,gk_y1
     READ(InUnitOro,*) def_str,gk_y0
     READ(InUnitOro,*) def_str,gk_x1
     READ(InUnitOro,*) def_str,gk_x0
     READ(InUnitOro,*) def_str,nyH
     READ(InUnitOro,*) def_str,nxH
     IF(fconv_ingrid=='spherical'.OR.fout_wahl=='G') THEN
        Call OGKGEO(gk_y0/1000, gk_x0/1000, breite, laenge)
        y0H=(breite*4.0d0*ATAN(1.0d0))/180
        x0H=(laenge*4.0d0*ATAN(1.0d0))/180
        Call OGKGEO(gk_y1/1000, gk_x1/1000, breite, laenge)
        y1H=(breite*4.0d0*ATAN(1.0d0))/180
        x1H=(laenge*4.0d0*ATAN(1.0d0))/180
     ELSE
        y1H=gk_y1; y0H=gk_y0
        x1H=gk_x1; x0H=gk_x0
     END IF
  ELSE 
     !IF(TRIM(indata_type)=="wgs84") THEN
     IF (TRIM(findata_type)=="geo") THEN
       READ(InUnitOro,*) def_str,grdy1,miny1,seky1   
       READ(InUnitOro,*) def_str,grdy0,miny0,seky0
       READ(InUnitOro,*) def_str,grdx1,minx1,sekx1
       READ(InUnitOro,*) def_str,grdx0,minx0,sekx0 
       READ(InUnitOro,*) def_str,nyH
       READ(InUnitOro,*) def_str,nxH
       !wgs84-->geo 
       gy1H=grdy1+miny1/60+seky1/3600
       gy0H=grdy0+miny0/60+seky0/3600
       gx1H=grdx1+minx1/60+sekx1/3600
       gx0H=grdx0+minx0/60+sekx0/3600
       !geo-->rad
       y1H=(gy1H*4.0d0*ATAN(1.0d0))/180.0d0
       y0H=(gy0H*4.0d0*ATAN(1.0d0))/180.0d0
       x1H=(gx1H*4.0d0*ATAN(1.0d0))/180.0d0
       x0H=(gx0H*4.0d0*ATAN(1.0d0))/180.0d0
     ELSE
       Write(*,*) "Funktion Oro: Spezifiziert zu 'gk', 'geo', 'wgs84'"
       STOP "Input Koordinate nicht spezifiziert zur Funktion Oro!"
     END IF
  END IF       

  IF(Domain%x0<x0H.OR.Domain%x1>x1H.OR. &
     Domain%y0<y0H.OR.Domain%y1>y1H )THEN
     Write(*,*) "Wertebereich Def.Domain und Koord-Input-WGS"
     Write(*,*) "Domain%x0=",Domain%x0,"  Domain%x1=",Domain%x1
     Write(*,*) "x0H      =",x0H,"  x1H      =",x1H
     Write(*,*) "Domain%y0=",Domain%y0,"  Domain%y1=",Domain%y1
     Write(*,*) "y0H      =",y0H,"  y1H      =",y1H
     STOP "STOP Read_WGS84_Oro"
  ELSE
    dyH=(y1H-y0H)/(nyH-1)
    dxH=(x1H-x0H)/(nxH-1)
  
    ALLOCATE(xPH(nxH))
    ALLOCATE(yPH(nyH))
    ALLOCATE(Height(nxH,nyH))
    DO iy=nyH,1,-1
        READ(InUnitOro,*) (Height(ix,iy),ix=1,nxH)
    END DO
    maxHOro=MAXVAL(Height)
    Write(*,*) "z-Max-Orographie=",maxHOro, " ----> aus Read_WGS84_Oro"
  
    xPH(1)=x0H
    DO i=2,nxH
      xPH(i)=xPH(i-1)+dxH
    END DO
    yPH(1)=y0H
    DO j=2,nyH
      yPH(j)=yPH(j-1)+dyH
    END DO
  END IF

  CLOSE(UNIT=InUnitOro)
END SUBROUTINE Read_WGS84_Oro


SUBROUTINE Read_Corine_List
  INTEGER             :: InUnitCorine=8
  CHARACTER(100)      :: read_line
  CHARACTER (LEN=20)  :: nroutine     ! name of this subroutine
  CHARACTER (LEN=40)  :: nerrmsg      ! error message
  INTEGER             :: nstat,cstat  ! error status variables
  INTEGER             :: i

  ! Datei-Orig: CLC.csv geändert!!!, Label 1-3 teilweise Hochkomma Nachtrag !!!
  ! Datei     : CLC_IfT.csv , erweitert um Spalte mat_IfT (1-9mat)
  OPEN(UNIT=InUnitCorine,FILE=TRIM(file_clc_list),STATUS='OLD',IOSTAT=nstat)
  IF(nstat/=0) THEN
     nroutine='READ_Corine_List'
     nerrmsg = "OPENING OF FILE \'file_clc_list\' FAILED"
     Write(*,*) nerrmsg, nroutine
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
  END IF

  ALLOCATE(v_clc(1:Nr_Land_Cl),STAT=cstat)
  IF(cstat/=0) THEN
     nroutine='READ_Corine_List'
     nerrmsg = 'Allocate-Error:  "v_clc" '
     Write(*,*) nerrmsg, nroutine
     !CALL tool_break (my_cart_id, 1003, nerrmsg, nroutine)
  END IF
  READ(InUnitCorine,*) read_line
  DO i=1,Nr_Land_Cl
      READ(InUnitCorine,*) v_clc(i)
  END DO
  CLOSE(UNIT=InUnitCorine)
END SUBROUTINE Read_Corine_List

SUBROUTINE Read_Corine_Data

  REAL(8)    :: grdy1,miny1,seky1  ! read 
  REAL(8)    :: grdy0,miny0,seky0  ! read
  REAL(8)    :: grdx1,minx1,sekx1  ! read
  REAL(8)    :: grdx0,minx0,sekx0  ! read
  INTEGER    :: nxLC,nyLC          ! read
  INTEGER    :: InUnitCorine=8
  CHARACTER (LEN=20)  :: def_str,l_str
  CHARACTER (LEN=512) :: read_line
  CHARACTER (LEN=20)  :: nroutine     ! name of this subroutine
  CHARACTER (LEN=40)  :: nerrmsg      ! error message
  INTEGER             :: nstat,cstat  ! error status variables
  INTEGER             :: ix,iy

  IF (ALLOCATED(LandCover)) THEN
    RETURN
  END IF
  OPEN(UNIT=InUnitCorine,FILE=TRIM(file_ncorine),STATUS='OLD',IOSTAT=nstat)
  IF(nstat/=0) THEN
     nroutine='READ_Corine_Data'
     nerrmsg = "OPENING OF FILE \'file_ncorine\' FAILED"
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
  ELSE 
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdy1,miny1,seky1,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdy0,miny0,seky0,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdx1,minx1,sekx1,l_str
    READ(InUnitCorine,*) read_line !to skip sth.: def_str,grdx0,minx0,sekx0,l_str
    READ(InUnitCorine,*) def_str,nyLC
    READ(InUnitCorine,*) def_str,nxLC

    IF(nxH/=nxLC.OR.nyH/=nyLC)THEN
      Write(*,*) "Keine Übereinstimmung Anzahl: 'nxH'  'nyH'   aus Read_WGS84_Oro   "
      Write(*,*) "                  mit Anzahl: 'nxLC' 'nyLC'  aus Read_Corine_Data ! " 
      !CALL tool_message(my_cart_id,2001,nwarning,nroutine))
    END IF

    ALLOCATE(LandCover(nxLC,nyLC),STAT=cstat)
    DO iy=nyLC,1,-1
      READ(InUnitCorine,*) (LandCover(ix,iy),ix=1,nxLC)
    END DO
    CLOSE(UNIT=InUnitCorine)
  END IF
END SUBROUTINE Read_Corine_Data



SUBROUTINE  Read_WGS84_Surf
  INTEGER :: i,j,ix,iy
  !optional surfit/bispev
  !gleiche Nutzung global def.: kx,ky,tx_nx,tx,ty_ny,ty,c_spl,wrk2,lwrk2,iwrk,kwrk
  INTEGER :: iopt,m,nxest,nyest,nmax,is,lwrk1,ier
  INTEGER :: mx,my 
  REAL(8) :: xb,xe,yb,ye,s,eps,fp
  REAL(8) :: xx(11),yy(11),zz(121) !Bsp:
  !
  INTEGER :: u,v,km,ne,bx,by,b1,b2
  REAL(8) :: ai,delta,ww,yp
  REAL(8) :: m_half,xmin,xmax,ymin,ymax
  REAL(8) :: diffmx,diffmy
  INTEGER :: InUnitV2=7
  REAL(8) :: grdy1,miny1,seky1,gy1H,gk_y1
  REAL(8) :: grdy0,miny0,seky0,gy0H,gk_y0
  REAL(8) :: grdx1,minx1,sekx1,gx1H,gk_x1
  REAL(8) :: grdx0,minx0,sekx0,gx0H,gk_x0
  REAL(8) :: breite,laenge
  REAL(8) :: x_datmin,y_datmin,min_px,min_py
  CHARACTER (LEN=20)  :: def_str

  IF(ALLOCATED(c_spl)) THEN
    RETURN
  END IF

  READ(10,*) iopt   ! specify whether a weighted spline variation
  READ(10,*) kx     ! degree of the spline x-direction
  READ(10,*) ky     ! degree of the spline y-direction
  READ(10,*) s      ! smoothing factor

  OPEN(UNIT=InUnitV2,FILE=TRIM(file_namefkt),STATUS='OLD')
  IF(TRIM(findata_type)=="gk") THEN
     READ(InUnitV2,*) def_str,gk_y1
     READ(InUnitV2,*) def_str,gk_y0
     READ(InUnitV2,*) def_str,gk_x1
     READ(InUnitV2,*) def_str,gk_x0
     READ(InUnitV2,*) def_str,nyH
     READ(InUnitV2,*) def_str,nxH
     IF(fconv_ingrid=='spherical'.OR.fout_wahl=='G') THEN
        Call OGKGEO(gk_y0/1000, gk_x0/1000, breite, laenge)
        y0H=(breite*4.0d0*ATAN(1.0d0))/180
        x0H=(laenge*4.0d0*ATAN(1.0d0))/180
        Call OGKGEO(gk_y1/1000, gk_x1/1000, breite, laenge)
        y1H=(breite*4.0d0*ATAN(1.0d0))/180
        x1H=(laenge*4.0d0*ATAN(1.0d0))/180
     ELSE
        y1H=gk_y1; y0H=gk_y0
        x1H=gk_x1; x0H=gk_x0
     END IF
  ELSE 
     !IF(TRIM(indata_type)=="wgs84") THEN
     IF (TRIM(findata_type)=="geo") THEN
       READ(InUnitV2,*) def_str,grdy1,miny1,seky1   
       READ(InUnitV2,*) def_str,grdy0,miny0,seky0
       READ(InUnitV2,*) def_str,grdx1,minx1,sekx1
       READ(InUnitV2,*) def_str,grdx0,minx0,sekx0 
       READ(InUnitV2,*) def_str,nyH
       READ(InUnitV2,*) def_str,nxH
       !wgs84-->geo 
       gy1H=grdy1+miny1/60+seky1/3600
       gy0H=grdy0+miny0/60+seky0/3600
       gx1H=grdx1+minx1/60+sekx1/3600
       gx0H=grdx0+minx0/60+sekx0/3600
       !geo-->rad
       y1H=(gy1H*4.0d0*ATAN(1.0d0))/180.0d0
       y0H=(gy0H*4.0d0*ATAN(1.0d0))/180.0d0
       x1H=(gx1H*4.0d0*ATAN(1.0d0))/180.0d0
       x0H=(gx0H*4.0d0*ATAN(1.0d0))/180.0d0
     ELSE
       Write(*,*) "Funktion Oro: Spezifiziert zu 'gk', 'geo', 'wgs84'"
       STOP "Input Koordinate nicht spezifiziert zur Funktion Oro!"
     END IF
  END IF       

  IF(Domain%x0<x0H.OR.Domain%x1>x1H.OR. &
     Domain%y0<y0H.OR.Domain%y1>y1H )THEN
     Write(*,*) "Wertebereich Def.Domain und Koord-Input-WGS"
     Write(*,*) "Domain%x0=",Domain%x0,"  Domain%x1=",Domain%x1
     Write(*,*) "x0H      =",x0H,"  x1H      =",x1H
     Write(*,*) "Domain%y0=",Domain%y0,"  Domain%y1=",Domain%y1
     Write(*,*) "y0H      =",y0H,"  y1H      =",y1H
     STOP "STOP Read_WGS84_Surf"
  ELSE

  ! m denotes the number of data points
  ! m>=(kx+1)*(ky+1) !!!Wichtig!!!
      !read(InUnitV2,*) m
      m=nyH*nxH
      ALLOCATE(x(1:m))
      ALLOCATE(y(1:m))
      ALLOCATE(z(1:m))
      ALLOCATE(w(1:m))

  ! fetch the co-ordinate and function values of each data point.
      !DO  i=1,m
      !  read(InUnitV2,*) y(i),x(i),z(i)
      !END DO 
      !In Datei ist 1.Zeile/1.Zahl linke obere Ecke (north/west)des 
      !Koordinatensystems
      DO  iy=nyH,1,-1
         i=(iy-1)*nxH
         read(InUnitV2,*) (z(i+ix),ix=1,nxH)
      END DO 
      maxHOro=MAXVAL(z)
      Write(*,*) "z-Max-Orographie=",maxHOro, " ----> aus Read_WGS84_Surf"
      dyH=(y1H-y0H)/(nyH-1)
      dxH=(x1H-x0H)/(nxH-1)
      ! nur Zwischenergebnis 
      ! ......................
    ALLOCATE(xPH(nxH))
    ALLOCATE(yPH(nyH))
      xPH(1)=x0H
      DO i=2,nxH
        xPH(i)=xPH(i-1)+dxH
      END DO
      yPH(1)=y0H
      DO j=2,nyH
       yPH(j)=yPH(j-1)+dyH
      END DO
      ! ......................

      DO iy=1,nyH
        DO  ix=1,nxH
         i=(iy-1)*nxH+ix
         x(i)=x0H+dxH*(ix-1)
        END DO 
      END DO 
      DO iy=1,nyH
        yp=y0H+(dyH*(iy-1))
        DO  ix=1,nxH
         i=(iy-1)*nxH+ix
         y(i)=yp
        END DO 
      END DO 
      x_datmin=x(1)
      y_datmin=y(1)
 
!!      IF(TRIM(file_namefkt).NE. 'dasurf.htm') THEN
!!        IF(conv_gk=='s'.OR.out_surf=='G') THEN
!!          DO  i=1,m
!!            Call OGKGEO(y(i)/1000, x(i)/1000, breite, laenge)
!!            y(i)=(breite*4.0d0*ATAN(1.0d0))/180
!!            x(i)=(laenge*4.0d0*ATAN(1.0d0))/180
!!          END DO
!!        ELSE 
!!          !!Betreff Nullpunkt ermitteln 
!!          !!Wenn kleinste Koordinate aus *.dat beachtet werden soll
!!          !!z.Zt nur berechnet aber nicht verwendet
!!          !!gauss[x,y]_min aus grid-Domain nur Bezug, 
!!          !!Nullpunkt wird erst beim OutputGMV verschoben
!!          IF(gaussx_min<x_datmin) THEN
!!            min_px=gaussx_min
!!            diffgx=0.0
!!          ELSE
!!            min_px=x_datmin
!!            diffgx=gaussx_min-x_datmin
!!          END IF
!!          IF(gaussy_min<y_datmin) THEN
!!            min_py=gaussy_min
!!            diffgy=0.0
!!          ELSE
!!            min_py=y_datmin
!!            diffgy=gaussy_min-y_datmin
!!          END IF
!!          !DO  i=1,m
!!            !!Nullpunkt verschieben
!!            !  x(i)=x(i)-min_px
!!            !  y(i)=y(i)-min_py  
!!          !END DO
!!        END IF
!!      END IF ! nicht Input dasurf.htm

 
  ! fetch an estimate of the standard deviation of the data values.
      !read(InUnitV2,*) delta
      delta=1.0d-12
  ! the weights are set equal to delta**(-1)
      ww= 97.49 !ww = 1./delta
      DO  i=1,m
        w(i) = ww 
      END DO 

  ! set up the boundaries of the approximation domain.
      !xb = -2.   !bedingt Bsp.: dasurf.htm
      !xe = 2.    ! ''
      !yb = -2.    ! ''
      !ye = 2.     ! ''
      xmin=x(1)
      xmax=x(1)
      ymin=y(1)
      ymax=y(1)
      DO i=2,m
        xmin=MIN(xmin,x(i))
        xmax=MAX(xmax,x(i))
        ymin=MIN(ymin,y(i))
        ymax=MAX(ymax,y(i))
      END DO 
      diffmx=xmax-xmin
      diffmy=ymax-ymin
      IF(conv_gk=='s'.OR.out_surf=='G') THEN
        xb=xmin-0.01*diffmx
        xe=xmax+0.01*diffmx
        yb=ymin-0.01*diffmy
        ye=ymax+0.01*diffmy
      ELSE IF (TRIM(findata_type)=="geo") THEN
        xb=xmin-1.0d-12*diffmx
        xe=xmax+1.0d-12*diffmx
        yb=ymin-1.0d-12*diffmy
        ye=ymax+1.0d-12*diffmy
      ELSE
        xb=REAL(FLOOR(xmin))   ! nicht fuer Bogenmass
        xe=REAL(CEILING(xmax)) ! nicht fuer Bogenmass
        yb=REAL(FLOOR(ymin))   ! nicht fuer Bogenmass
        ye=REAL(CEILING(ymax)) ! nicht fuer Bogenmass
      END IF
  ! generate a rectangular grid for evaluating the splines.
  !      mx = 11
  !      my = 11
  !      DO i=1,11
  !        ai = i-6
  !        xx(i) = ai*0.4
  !        yy(i) = xx(i)
  !      END DO 


  ! integer flag. on entry iopt must specify whether a weighted
  ! least-squares spline (iopt=-1) or a smoothing spline (iopt=
  ! 0 or 1) must be determined.
     ! iopt = 0    ! read
     ! kx = 3      ! read
     ! ky = 3      ! read
     ! s = 900000.   !s= 30 !read


  ! set up the dimension information.
  !    -must specify an upper bound for the number of knots required
  !     in the x- and y-directions respect.
  !    -nxest >= 2*(kx+1), nyest >= 2*(ky+1)
  !    -in most practical situation 
  !     nxest = kx+1+sqrt(m/2), nyest = ky+1+sqrt(m/2) will be sufficient.
      m_half=m/2
      nxest = kx+1+CEILING(SQRT(m_half)) 
      nyest = ky+1+CEILING(SQRT(m_half))
      IF (nxest<(2*kx+2)) THEN
        nxest=2*kx+2
      END IF
      IF (nyest<(2*ky+2)) THEN
        nyest=2*ky+2
      END IF
      nmax = MAX(nxest,nyest)
      ALLOCATE(tx(1:nmax))
      ALLOCATE(ty(1:nmax))
      ALLOCATE(c_spl((nxest-kx-1)*(nyest-ky-1)))
      kwrk = m+(nxest-2*kx-1)*(nyest-2*ky-1)
      ALLOCATE(iwrk(1:kwrk))
      !computation for 'lwrk1'-> must specify the actual dimension of wrk1
      u = nxest-kx-1
      v = nyest-ky-1
      km = max(kx,ky)+1
      ne = nmax
      bx = kx*v+ky+1
      by = ky*u+kx+1
      if(bx.le.by) THEN
           b1 = bx
           b2 = b1+v-ky
      END IF
      if(bx.gt.by) THEN
           b1 = by
           b2 = b1+u-kx
      END IF
      lwrk1=u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
      ALLOCATE(wrk1(1:lwrk1))
      lwrk2=u*v*(b2+1)+b2 
      ALLOCATE(wrk2(1:lwrk2))
      !im Bsp surfit def. : nxest = 15,nyest = 15,nmax = 15
      !                     !kwrk = 300, lwrk1 = 12000, lwrk2 = 6000
 
  ! choose a value for eps
      eps=0.1e-05

  ! integer flag. on entry iopt must specify whether a weighted
  ! least-squares spline (iopt=-1)
 
      IF(iopt==-1) THEN
        !kx = 3
        !ky = 3
        tx_nx=11
        ty_ny=11
        j=kx+2
        DO i=1,3
          ai = i-2
          tx(j)=ai
          ty(j)=ai
          j=j+1
        END DO
      ELSE
        tx_nx=nmax
        ty_ny=nmax
      END IF

  ! spline approximations of degree k
      Write(*,*) "spline approximations of degree k=",kx
      call surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest, &
        nmax,eps,tx_nx,tx,ty_ny,ty,c_spl,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)

  ! evaluation of the spline approximation
  !      call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xx,mx,yy,my,zz,   &
  !         wrk2,lwrk2,iwrk,kwrk,ier)
  END IF

  CLOSE(UNIT=InUnitV2)
END SUBROUTINE Read_WGS84_Surf


SUBROUTINE SelectReadFkt(nr_wahlfkt)
INTEGER nr_wahlfkt
  SELECT CASE (nr_wahlfkt)
    CASE (1)
      CALL Read_Val_Fkt
    CASE (2)
      CALL Read_Val_Fkt2
    CASE (3)
      CALL Read_Val_Fkt3
    CASE (4)
      CALL Read_Val_Fkt4
    CASE (5)
      CALL Read_Val_Fkt5
    CASE (6)
      CALL Read_Val_Fkt6
    CASE (7)
      CALL Read_Val_Fkt7
    CASE (8)
      CALL Read_Val_Fkt8
    CASE (9)
      CALL Read_Val_Fkt9
    CASE (10)
      CALL Read_Val_Haus
    CASE (11)
      CALL Read_Val_Haus_GML
    CASE (12)
      CALL Read_Oro
    CASE (13)
      CALL Read_Val_FktKugel
    CASE (14)
      CALL Read_Val_FktSierra
    CASE (15)
      CALL Read_Val_FktSurfit
    CASE (16)
      CALL Read_Hill_Froeh
    CASE (17)
      CALL Read_WGS84_Oro
    CASE (18)
      CALL Read_WGS84_Surf
    CASE (nr_fkt_gml_oro)
      CALL Read_Val_Haus_GML
    CASE DEFAULT
      H=200.0
      aUP=1000.0
      aDOWN=1000.0
      WRITE(*,*) 'Funktion Agnesi als Standard definiert'
      WRITE(*,*) 'H = ',H
      WRITE(*,*) 'aUP = ',aUP
      WRITE(*,*) 'aDOWN = ',aDOWN
      !CALL Read_Val_Fkt      
  END SELECT
END SUBROUTINE SelectReadFkt


FUNCTION f1(x,y,z)  ! Agnesi
    REAL(8) :: f1
    REAL(8) :: x,y,z
    REAL(8) :: x0=3000.0 

    IF (x<0) THEN
      f1 = z-H/(1.0+ (x/aUP)**2)
    ELSE 
      f1 = z-H/(1.0+ (x/aDOWN)**2)
    END IF

    ! f1 = z-H/(1.0+ (x/a)**2)
    ! f1 = z-(H/(1.0+ ((x-x0)/a)**2))
END FUNCTION f1

FUNCTION f2(x,y,z)  ! Agnesi 3D 
    REAL(8) :: f2
    REAL(8) :: x,y,z

    f2 = z-(H/(1.0+ (x/a)**2+ (y/b)**2))
    !f2 = z-(H/(1.0+ ((x-y)/a)**2))
END FUNCTION f2

FUNCTION f3(x,y,z)  ! Agnesi 2sided 2D
    REAL(8) :: f3
    REAL(8) :: x,y,z

    IF (x<0) THEN
      f3 = z-H/(1.0+ (x/aUP)**2)
    ELSE
      f3 = z-H/(1.0+ (x/aDOWN)**2)
    END IF
END FUNCTION f3

FUNCTION f4(x,y,z)  ! Agnesi Double 
    REAL(8) :: f4
    REAL(8) :: x,y,z
     f4=MAX(z-1000.0/(1.0+(x/10000.0)**2.0+(y/1.0d15)**2.0) &
             -1000.0/(1+(y/10000.0)**2.0),0.0d0)

!    f3 = z-(H/(1.0+ ((x-d/2.0)/a)**2.0+ (y/b)**2.0) &
!           +H/(1.0+ ((x+d/2.0)/a)**2.0+ (y/b)**2.0))
END FUNCTION f4

FUNCTION f5(x,y,z)  ! Schaer
    REAL(8) :: f5
    REAL(8) :: x,y,z
    REAL(8) :: Pi
    Pi=ATAN(1.0d0)*4.0d0
    f5=z-H*EXP(-x**2/a**2)*(COS(Pi*x/Lambda))**2
END FUNCTION f5


FUNCTION f6(x,y,z)  ! Concave 
    REAL(8) :: f6
    REAL(8) :: x,y,z
    REAL(8) :: x1,y1,xc
    REAL(8) :: Pi

    IF (y>=d) THEN
      x1=x*COS(alpha)+(y-b)*SIN(alpha)
      y1=(y-b)*COS(alpha)-x*SIN(alpha)
    ELSE IF (y<=-d) THEN
      x1=x*COS(alpha)-(y+b)*SIN(alpha)
      y1=(y+b)*COS(alpha)+x*SIN(alpha)
    ELSE
      xc=(b-d/2-y**2/(2*d))*TAN(alpha)
      x1=(x-xc)*COS(alpha)
      y1=0
    END IF

    IF (ABS(y)>c) THEN
      IF (y*y1>=0) THEN
        f6=H/(1+x1**2/a**2+y1**2/a**2)
      ELSE
        f6=H/(1+x1**2/a**2)
      END IF
    ELSE
      Pi=ATAN(1.0d0)*4.0d0
      f6=H*(1+beta*COS(Pi*y/c)**2)/(1+x1**2/a**2)
    END IF
END FUNCTION f6

FUNCTION f7(x,y,z)    ! für Bannon 2D
    REAL(8) :: f7
    REAL(8) :: x,y,z
    REAL(8) :: Pi

    Pi=ATAN(1.0d0)*4.0d0
!    f7=z-H*(1.0d0+mu*cos(n*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0)  !!!ORIGINAL
!    f7=z-H*(1.0d0+mu*cos(n*Pi*x/a+(phi*Pi)))/(1.0d0+(x/a)**Pot)

    f7=z-H*(1.0d0+mu*cos(nn*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0)
END FUNCTION f7

FUNCTION f8(x,y,z)    ! für Bannon 3D
    REAL(8) :: f8
    REAL(8) :: x,y,z
    REAL(8) :: Pi

    Pi=ATAN(1.0d0)*4.0d0
!    f8=z-H*(1.0d0+mu*cos(n*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0+(y/b)**2.0d0) !!!ORIGINAL
!    f8=z-H*(1.0d0+mu*cos(n*Pi*x/a+(phi*Pi)))/(1.0d0+(x/a)**2.0d0+(y/b)**Pot)

    f8=z-H*(1.0d0+mu*cos(nn*Pi*x/a+phi))/(1.0d0+(x/a)**2.0d0+(y/b)**2.0d0)
END FUNCTION f8

FUNCTION f9(x,y,z)    ! Kegel 
    REAL(8) :: f9
    REAL(8) :: x,y,z

    REAL(8) :: r
    REAL(8) :: hx,hy
    REAL(8) :: fOut
    CHARACTER*20 :: Casex,Casey
    REAL(8) :: Grav
    REAL(8) :: eta,eta0,U0,eta_vertical,SurfGeoPot,SurfGeoPot0
    REAL(8) :: lam,phi,rot_lam,rot_phi


    !Pi=ATAN(1.0d0)*4.0d0

!    f9=z-(x-xStart)/(xEnd-xStart)*h  !Rampe
  SELECT CASE (NameFun)
    CASE ('Annulus')
      r=SQRT(x*x+y*y)
      f9=MIN(r-RIn,Rout-r)
    CASE ('Leer')
      f9=MIN(y-r1,r2-y)
    CASE ('ValleyTwo')
      IF (ABS(x)<=Vx) THEN
        hx=0.0d0
        Casex='1x'
      ELSE IF (ABS(x)<=Vx+Sx) THEN
        hx=0.5d0-0.5d0*COS(Pi*(ABS(x)-Vx)/Sx)
        Casex='2x'
      ELSE IF (ABS(x)<=Vx+Sx+Px) THEN
        hx=1.0d0
        Casex='3x'
      ELSE IF (ABS(x)<=Vx+2.0d0*Sx+Px) THEN
        hx=0.5d0+0.5d0*COS(Pi*(ABS(x)-(Vx+Sx+Px))/Sx)
        Casex='4x'
      ELSE
        hx=0.0d0
        Casex='5x'
      END IF
      IF (ABS(y)<=Py) THEN
        hy=1.0d0
        Casey='1y'
      ELSE IF (ABS(y)<=Py+Sy) THEN
        hy=0.5d0+0.5d0*COS(Pi*(ABS(y)-Py)/Sy)
        Casey='2y'
      ELSE
        hy=0.0d0
        Casey='3y'
      END IF
      f9=z-hP*hx*hy
!     IF (f9>0.0d0) THEN
!       WRITE(*,*) casex,Casey
!       WRITE(*,*) f9,hp,hx,hy
!     END IF
    CASE ('Valley3D')
      IF (ABS(x)<=Vx) THEN
        hx=0.0d0
      ELSE IF (ABS(x)<=Vx+Sx) THEN
        hx=0.5d0-0.5d0*COS(Pi*(ABS(x)-Vx)/Sx)
      ELSE
        hx=1.0d0
      END IF
      hy=0.5d0+0.5d0*TANH(y/Sy)
      f9=z-hP*hx*hy
    CASE ('Valley2D')
      hy=0.5d0+0.5d0*TANH(y/Sy)
      f9=z-hP*hy
    CASE ('Kegel')
      IF (x<=xL) THEN
        f9=z
      ELSE IF (x<=xC.AND.xC>xL) THEN
        f9=z-MIN((x-xL)/(xC-xL)*h,hMin)
      ELSE IF (x<=xR.AND.xR>xC) THEN
        f9=z-MIN((xR-x)/(xR-xC)*h,hMin)
      ELSE
        f9=z
      END IF
    CASE ('Inlet')
      IF (x<=x1Inlet) THEN
        f9=0.0d0
        f9=z-f9
      ELSE IF (x<=x2Inlet) THEN
        f9=h2Inlet/(x2Inlet-x1Inlet)**2.0d0*(x-x1Inlet)**2.0d0
        f9=z-f9
      ELSE IF (x<=x4Inlet) THEN
        f9=h3Inlet-(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x-x3Inlet)**2.0d0
        f9=z-f9
      ELSE   
        f9=h3Inlet-(h3Inlet-h2Inlet)/(x2Inlet-x3Inlet)**2.0d0*(x4Inlet-x3Inlet)**2.0d0 &
           -h4Inlet*(x-x4Inlet)/(x5Inlet-x4Inlet)
        f9=z-f9
      END IF
      IF (x>=xOut1.AND.x<=xOut3) THEN
        fOut=hOut1
        f9=MIN(f9,fOut-z)
      END IF  
      IF (x>=xOut2.AND.x<=xOut3) THEN
        fOut=hOut2
        f9=MAX(f9,z-fOut)
      END IF  
      f9=MIN(f9,h5Inlet-z)
    CASE ('BaroIn')
      Grav=9.81d0
      eta=1.0d0
      eta0=0.252d0
      U0=35.0d0
      phi=0.5d0*Pi
      SurfGeoPot0   = U0*1.5d0*(COS(eta_vertical))**1.5d0*                                          &
                     ((-2.d0*(SIN(phi))**6.0d0*((COS(phi))**2.0d0+1.d0/3.d0)+10.d0/63.d0)* &
                      U0*(COS(eta_vertical))**1.5d0  +                                             &
                     (8.d0/5.d0*(COS(phi))**3.0d0*((SIN(phi))**2.0d0+2.d0/3.d0)            &
                     - pi/4.d0)*Omega*RadEarth)
      lam=x
      phi=y
      CALL Rotate(lam,phi,rot_lam,rot_phi,RotAngle)
      eta_vertical = (eta - eta0) * 0.5d0*pi
      SurfGeoPot   = U0*1.5d0*(COS(eta_vertical))**1.5d0*                                          &
                     ((-2.d0*(SIN(rot_phi))**6.0d0*((COS(rot_phi))**2.0d0+1.d0/3.d0)+10.d0/63.d0)* &
                      U0*(COS(eta_vertical))**1.5d0  +                                             &
                     (8.d0/5.d0*(COS(rot_phi))**3.0d0*((SIN(rot_phi))**2.0d0+2.d0/3.d0)            &
                     - pi/4.d0)*Omega*RadEarth)-SurfGeoPot0
      f9=z-SurfGeoPot/Grav               
    CASE Default
       !f9=1.0d0
       !f9=0.0d0
       !f9=-10.0d0
       f9=z+1.d20 ! 
  END SELECT
END FUNCTION f9


FUNCTION DistAll(x,y,z)
  REAL(8) :: DistAll
  REAL(8) :: x,y,z
  INTEGER :: i,h
  TYPE(Point_T) ::Pdis
  TYPE(Box_T), POINTER :: Box
  REAL(8):: disp,distz0

  Pdis%x=x
  Pdis%y=y
  Pdis%z=z
  DistAll=+1.0d99
  CALL FindBoxPoint(Root,Box,PDis)
  !Distance 
  IF (z==0.0d0) THEN
    !DistAll=-1.0d-8   ! Test ContainerT_Verschiebe.grid V8.9.3.1
                       ! wenn inout=0 sein soll wird -1 gesetzt, 
                       ! Im Rücksprung in  CheckVertex auf > -1.0d-12 getestet 
                       ! kommt dann immer in den Else-Zweig --> in_out = -1
    DistAll=-1.0d-13  ! Test ContainerT_Verschiebe.grid V8.9.3.1 ok.
    DistAll=-0.48d-12 ! Anpassung um näher -1.0d-12 (CheckVertex), wenn nahe 0 ist 
                      ! Test ContainerT_Verschiebe.grid V8.9.3.1  ok.
    DistAll=-(dist_fscv*0.48) ! DistAll=-(48%dist_fscv) ,negativ value is stringent!
    distz0=DistAll
  END IF
  DO
    IF (ASSOCIATED(Box)) THEN
      IF (ASSOCIATED(Box%List)) THEN
        DO i=1,SIZE(Box%List)
          disp=Dist(Pdis,Haus(Box%List(i)))
          DistAll=MIN(DistAll,disp)
        END DO
      END IF
      Box=>Box%Parent
    ELSE
      EXIT
    END IF
  END DO
  !DistAll=DistAll-1.0d-8
  !DistAll=DistAll-1.0d-13
  if(z==0.0d0) THEN
    DistAll=DistAll-(dist_fscv*0.48)
    !Write(*,*) "DistAll=", DistAll
    
  end if
END FUNCTION DistAll


FUNCTION fOro(x,y,z)
  REAL(8) :: fOro
  REAL(8) :: x,y,z

  INTEGER :: iPos,jPos
  REAL(8) :: x1,x2,y1,y2
  REAL(8) :: h11,h12,h21,h22
  REAL(8) :: Dist

  IF      (nxH<=0) THEN
    fOro=0.d0
    GOTO 1
  ELSE IF (nxH==1) THEN
    fOro=Height(1,1)
    GOTO 1
  END IF

  ! 2-dimensionale Interpolation
  iPos=(x-x0H)/dxH+1
  jPos=(y-y0H)/dyH+1

  x1=xPH(iPos)
  x2=xPH(iPos+1)
  y1=yPH(jPos)
  y2=yPH(jPos+1)
  h11=Height(iPos,jPos)
  h21=Height(iPos+1,jPos)
  h12=Height(iPos,jPos+1)
  h22=Height(iPos+1,jPos+1)
  !wenn iPos,jPos Grenze: +dxH/2 od. +dyH/2 
  ! ---> für Skalierung damit Dist != 0 werden kann
  !wenn im grid-File nx,ny Grenzen '0 bis n(x,y)H-1' angegeben wird 
  !-->  'IF' nicht benutzt
  IF(iPos==nxH) THEN
     x2=xPH(iPos)+dxH/2
     h21=Height(iPos,jPos)
     IF(jPos==nyH) THEN 
        h22=Height(iPos,jPos)
     ELSE 
        h22=Height(iPos,jPos+1)
     END IF
  END IF
  IF(jPos==nyH) THEN
     y2=yPH(jPos)+dyH/2
     h12=Height(iPos,jPos)   
     IF(iPos/=nxH)  h22=Height(iPos+1,jPos)
  END IF

  Dist=(x2-x1)*(y2-y1)
  fOro = &
        (x2-x)*(y2-y)*h11 &
       -(x1-x)*(y2-y)*h21 &
       -(x2-x)*(y1-y)*h12 &
       +(x1-x)*(y1-y)*h22
  fOro=fOro/Dist

1 CONTINUE

  fOro=z-fOro

END FUNCTION fOro


FUNCTION fOroWgs(x,y,z)
  REAL(8) :: fOroWgs
  REAL(8) :: x,y,z

  INTEGER :: iPos,jPos
  REAL(8) :: x1,x2,y1,y2
  REAL(8) :: h11,h12,h21,h22
  REAL(8) :: Dist
  !search pos into Height
  !IF(INT(x-x0H)==0) THEN
  !  iPos=1
  !ELSE
  !  iPos=INT((x-x0H)/dxH)
  !END IF
  !IF(INT(y-y0H)==0)THEN
  !  jPos=1
  !ELSE
  !  jPos=INT((y-y0H)/dyH)
  !END IF
  iPos=(x-x0H)/dxH+1
  jPos=(y-y0H)/dyH+1

  x1=xPH(iPos)
  x2=xPH(iPos+1)
  y1=yPH(jPos)
  y2=yPH(jPos+1)
  h11=Height(iPos,jPos)
  h21=Height(iPos+1,jPos)
  h12=Height(iPos,jPos+1)
  h22=Height(iPos+1,jPos+1)
  !wenn iPos,jPos Grenze: +dxH/2 od. +dyH/2 
  ! ---> für Skalierung damit Dist != 0 werden kann
  !wenn im grid-File nx,ny Grenzen '0 bis n(x,y)WGS-1' angegeben wird 
  !-->  'IF' nicht benutzt
  IF(iPos==nxH) THEN
     x2=xPH(iPos)+dxH/2
     h21=Height(iPos,jPos)
     IF(jPos==nyH) THEN
        h22=Height(iPos,jPos)
     ELSE
        h22=Height(iPos,jPos+1)
     END IF
  END IF
  IF(jPos==nyH) THEN
     y2=yPH(jPos)+dyH/2
     h12=Height(iPos,jPos)
     IF(iPos/=nxH)  h22=Height(iPos+1,jPos)
  END IF

  ! 2-dimensionale Interpolation
  Dist=(x2-x1)*(y2-y1)
  fOroWgs = &
        (x2-x)*(y2-y)*h11 &
       -(x1-x)*(y2-y)*h21 &
       -(x2-x)*(y1-y)*h12 &
       +(x1-x)*(y1-y)*h22
  fOroWgs = z-fOroWgs/Dist

END FUNCTION fOroWgs


FUNCTION f_r(lambda0,lambda1,phi0,phi1)
    REAL(8) :: f_r
    REAL(8) :: lambda0,lambda1,phi0,phi1

    f_r = RadEarth*ACOS(SIN(phi0)*SIN(phi1)+COS(phi0)*COS(phi1)*COS(lambda1-lambda0))
END FUNCTION f_r


FUNCTION f_h(lam,phi,z)  !shape of mountain
     REAL(8) :: f_h
     REAL(8) :: lam,phi,z
     REAL(8) :: r,PiNinth
     !h0, height at the center of the mountain
     !d, half-width of the mountain
     !r, distance from the center

     IF (TRIM(NameFun)=='Agnesi') THEN
       r=f_r(lam0,lam,phi0,phi)
       f_h = z-H/(1.0d0+(r/d)**2) !h-h0
     ELSE IF (TRIM(NameFun)=='Cone') THEN
       !PiNinth=ATAN(1.0d0)*4.0d0/9.0d0
       PiNinth=Pi/9.0d0
       r=MIN(SQRT((lam-lam0)**2+(phi-phi0)**2),PiNinth)
       f_h=z-H*(PiNinth-r)/PiNinth
     ELSE
       STOP 'False NameFun'
     END IF
END FUNCTION f_h


FUNCTION f_splev(x,y,z)
   REAL(8)  :: f_splev
   REAL(8)  :: x,y,z
   INTEGER, PARAMETER :: m=1
   INTEGER  :: ier
   REAL(8)  :: xs(m),sp(m)
   xs(1)=x
   call splev(t,n,b_coef,k,xs,sp,m,ier)
   f_splev=z-sp(m)
END FUNCTION f_splev


FUNCTION f_bispev(x,y,z)
   REAL(8)  :: f_bispev
   REAL(8)  :: x,y,z
   INTEGER, PARAMETER :: m=1
   INTEGER  :: ier
   REAL(8)  :: xg(1),yg(1),b_sp(1)
   xg(1)=x
   yg(1)=y
   call bispev(tx,tx_nx,ty,ty_ny,c_spl,kx,ky,xg,m,yg,m,b_sp, &
        wrk2,lwrk2,iwrk,kwrk,ier)
   f_bispev=z-b_sp(m)
END FUNCTION f_bispev


FUNCTION f_hill_froeh(x,y,z)
   REAL(8)  :: f_hill_froeh
   REAL(8)  :: x,y,z
   REAL(8)  :: zb,xLoc
   REAL(8), PARAMETER  :: HillHeight=28.0d0

   xLoc=x
   IF (xLoc>=4.5d0*HillHeight) THEN
     xLoc=MAX(9.0d0*HillHeight-xLoc,0.0d0)
   END IF

   IF (xLoc>=0.0d0.AND.xLoc<9.0d0) THEN
     zb=2.800000000000d+01 &
        +0.000000000000d+00*xLoc &
        +6.775070969851d-03*(xLoc**2.0) &
        -2.124527775800d-03*(xLoc**3.0)
     zb=MIN(zb,28.0d0)
   ELSE IF (xLoc>=9.0d0.AND.xLoc<14.0d0) THEN
     zb=2.507355893131d+01 &
      +9.754803562315d-01*xLoc &
      -1.016116352781d-01*(xLoc**2.0d0) &
      +1.889794677828d-03*(xLoc**3.0d0)
   ELSE IF (xLoc>=14.0d0.AND.xLoc<20.0d0) THEN
     zb= 2.579601052357d+01 &
      +8.206693007457d-01*xLoc &
      -9.055370274339d-02*(xLoc**2.0d0) &
      +1.626510569859d-03*(xLoc**3.0d0)
   ELSE IF (xLoc>=20.0d0.AND.xLoc<30.0d0) THEN
     zb= 4.046435022819d+01 &
      -1.379581654948d+00*xLoc &
      +1.945884504128d-02*(xLoc**2.0d0) &
      -2.070318932190d-04*(xLoc**3.0d0)
   ELSE IF (xLoc>=30.0d0.AND.xLoc<40.0d0) THEN
     zb= 1.792461334664d+01  &
      +8.743920332081d-01*xLoc &
      -5.567361123058d-02*(xLoc**2.0d0) &
      +6.277731764683d-04*(xLoc**3.0d0)
   ELSE IF (xLoc>=40.0d0.AND.xLoc<54.0d0) THEN
     zb=5.639011190988d+01 &
        -2.010520359035d+00*xLoc &
        +1.644919857549d-02*(xLoc**2.0d0) &
        +2.674976141766d-05*(xLoc**3.0d0)
     zb=MAX(0.0d0,zb)
   ELSE
     zb=0.0
   END IF
   f_hill_froeh=z-zb
END FUNCTION f_hill_froeh


FUNCTION sf(wahlfkt, x,y,z)
    REAL(8) :: sf 
    REAL(8) :: x,y,z
    INTEGER :: wahlfkt

    SELECT CASE (wahlfkt)
      CASE (1)
        sf=f1(x,y,z)
      CASE (2)
        sf=f2(x,y,z)
      CASE (3)
        sf=f3(x,y,z)
      CASE (4)
        sf=f4(x,y,z)
      CASE (5)
        sf=f5(x,y,z)
      CASE (6)
        sf=f6(x,y,z)
      CASE (7)
        sf=f7(x,y,z)
      CASE (8)
        sf=f8(x,y,z)
      CASE (9)
        sf=f9(x,y,z)
      CASE (10)
        sf=DistAll(x,y,z)
      CASE (11)
        sf=DistAll(x,y,z)
      CASE (12)
        sf=fOro(x,y,z)
      CASE (13)
        sf=f_h(x,y,z)
      CASE (14)
        sf=f_splev(x,y,z)
      CASE (15)
        sf=f_bispev(x,y,z)
      CASE (16)
        sf=f_hill_froeh(x,y,z)
      CASE (17)
        sf=foroWgs(x,y,z)
      CASE (18)
        sf=f_bispev(x,y,z)
      CASE (nr_fkt_gml_oro)
        sf=fOro(x,y,z)
        IF (sf>=-1.d-2) THEN
          sf=MIN(sf,DistAll(x,y,z))
        END IF
      CASE DEFAULT
        sf=f1(x,y,z)
    END SELECT
END FUNCTION sf
   

!#####################################################################
! Deallocate Routine
!#####################################################################

SUBROUTINE DEALLOC_Var_Sierra(nr_wahlfkt)
 INTEGER :: nr_wahlfkt 
  IF(nr_wahlfkt==nr_fkt_sierra) THEN
    DEALLOCATE(x)
    DEALLOCATE(y)
    DEALLOCATE(w)
    !DEALLOCATE(sp)
    DEALLOCATE(iwrk)
    DEALLOCATE(wrk)
    DEALLOCATE(t)
    DEALLOCATE(b_coef)
  END IF
END SUBROUTINE DEALLOC_Var_Sierra

SUBROUTINE DEALLOC_Var_SurfitProx(nr_wahlfkt)
 INTEGER :: nr_wahlfkt 
  IF(nr_wahlfkt==nr_fkt_surfit) THEN
    DEALLOCATE(x)
    DEALLOCATE(y)
    DEALLOCATE(z)
    DEALLOCATE(w)
    DEALLOCATE(tx)
    DEALLOCATE(ty)
    DEALLOCATE(c_spl)
    DEALLOCATE(iwrk)
    DEALLOCATE(wrk1)
    DEALLOCATE(wrk2)
  END IF
END SUBROUTINE DEALLOC_Var_SurfitProx

END MODULE f_mod

