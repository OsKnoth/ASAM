MODULE Pitzer_Mod

  IMPLICIT NONE

!---  dimensions
  INTEGER, PARAMETER :: RealKind=8
  INTEGER :: nc = 4, na = 5
  INTEGER :: nac, nacc, naca, nc2, na2

!---  INTEGER variable arrays
  INTEGER, ALLOCATABLE :: nzc(:),nza(:)          ! charges
  INTEGER, ALLOCATABLE :: nuec(:,:),nuea(:,:)    ! ion number per salt molecule

!--------------------------------------------------------
!---  REAL(RealKind) constants

!---  REAL(RealKind) parameter arrays
  REAL(RealKind), PRIVATE, ALLOCATABLE :: b0(:,:),b1(:,:),c(:,:),        &
                                          psim(:,:,:),tetax(:,:),        &
                                          psix(:,:,:),tetam(:,:)         

!---  REAL(RealKind) variables 
  REAL(RealKind) :: stri,aw
  REAL(RealKind), ALLOCATABLE :: gmx(:,:), gm(:), gx(:)
  REAL(RealKind), ALLOCATABLE :: cma(:), cmc(:)

!--------------------------------------------------------
!---  CHARACTER  arrays with names
  CHARACTER(20), ALLOCATABLE :: mname(:), xname(:)
  CHARACTER(20), ALLOCATABLE :: mxname(:,:)

!--------------------------------------------------------
   INTEGER, ALLOCATABLE :: PitzInd(:)

!  HILFSGROESSEN

   REAL(RealKind), ALLOCATABLE :: Charge(:)
   CHARACTER*20, ALLOCATABLE :: SpeciesNameAqua(:)
   INTEGER :: nAqua

CONTAINS

! -----------------------------------------------------------------------
REAL(RealKind) FUNCTION  aphi(t,p,eps,rdhl,el)
! -----------------------------------------------------------------------
!     aphi calculates the debye-hueckel parameter aphi, the reciprocal
!     debey-hueckel length rdhl, and the electrostatic length el
!
!     t    - temperature  [K]
!     p    - pressure    [hpa]      
!     eps  - dielectric constant at t,p
!     rdhl - debey-hueckel reciprocal length 
!     el   - electrostatic length el
!
!     eps-validity: 0-355°C, <1kbar
!
!     ref: aphi
!     -Pitzer, K.S., 1991: ion interaction approach: theory and data 
!      correlation; inactivity coefficients in electrolyte solutions,
!      ed. Pitzer, pp. 87, 130
!     -clegg & whitfield, ibid p. 297
!
!     ref. rdhl, el
!     -Pitzer, K.S., 1991: ion interaction approach: theory and data 
!      correlation; inactivity coefficients in electrolyte solutions,
!      ed. Pitzer, p. 22, 122
!
! -----------------------------------------------------------------------

!---  external variables
      REAL(RealKind) ::  t,p,eps,rdhl,el                                     

!---  external variables
      REAL(RealKind) ::  pi,pp,b,c1

      REAL(RealKind), PARAMETER ::     u1=3.4279d2, u2=-5.0866d-3, u3=9.46900d-7  &
&                ,u4=-2.0525d0, u5=3.1159d3, u6=-1.8289d+2, u7=-8.0325d+3  &
&                ,u8=4.21142d6, u9=2.1417d0

      REAL(RealKind), PARAMETER :: t0 = 273.15d0
      REAL(RealKind), PARAMETER :: densw = 1.,            & ! water density [g/cm**3]
&                           avo   = 6.02257d23,    & ! avogadro's number [molec./mol]
&                           bol   = 1.3804d-16,    & ! bolzmann's constant [erg/K]
&                           e     = 4.80298d-10      ! electronic charge [esu]
! -----------------------------------------------------------------------
!
      pi  = 4.d0 * ATAN(1.d0)
      pp  = p * 1.d-3    				! [mbar] -> [bar]  
      b   = u7 + u8/t + u9*t
      c1  = u4 + u5/(u6+t)
      eps = u1 * exp (u2*t+u3*t*t)
      eps = eps + c1*log((b+pp)/(b+1.d3))
!      
      aphi = 1.d0/3.d0 * sqrt (2*pi*avo*densw*1.d-3)*(e*e/eps/bol/t)**(1.5d0)

!---------------------------------------------------------------------------
!fm   aphi1=1.83e6*(eps*t)**(-1.5)
!fm   if (t.lt.t0) then
!fm     aphi=0.13422*(0.0368329*t-14.62718*log(t)-1530.1474/t+80.40631)
!fm   else
!fm     aphi=0.13422*(4.1725332-0.1481291*sqrt(t)+1.5188505e-5*t*t-
!fm  &       1.8016317e-8*t**3+9.3816144e-10*t**3.5)
!fm   endif 
!fm   print*,aphi,aphi1,aphi2
!---------------------------------------------------------------------------
!
!     debey-hueckel reciprocal length 
      rdhl = 4.d0 * pi * e * e /eps/bol/t
!
!     electrostatic length el
      el = e * e /eps/bol/t
!
! -----------------------------------------------------------------------
END FUNCTION aphi

! ----------------------------------------------------------------------
SUBROUTINE eteta (a,stri1,nzi,nzj,etetamn,etetaprmn)
! ----------------------------------------------------------------------
!
!     calculation of additional terms for unsymmetrical mixing:

!     considers the electrostatic forces between like signed ions but
!     unlike charges (mx & nx or mx & my)
!
!     ref: pitzer, k.s.: ion interaction approach: theory and data
!     correlation, in electrolyte solutions..., pp. 76-147, chap. 3
!
! ----------------------------------------------------------------------
      
!---  external variables
      INTEGER :: nzi,nzj
      REAL(RealKind) :: a,stri1,etetamn,etetaprmn

!---  internal variables
      INTEGER :: k, n
      REAL(RealKind) :: xj, zj, dzjdxj

      REAL(RealKind) :: ak(0:20,2),bk(0:22),dk(0:22),xx(3),xjxj(3),xjprxj(3)
      DATA ak /    &  ! ak(1)
&               1.925154014814667d0,-0.060076477753119d0,-0.029779077456514d0, &
&              -0.007299499690937d0, 0.000388260636404d0, 0.000636874599598d0, &
&               0.000036583601823d0,-0.000045036975204d0,-0.000004537895710d0, &
&               0.000002937706971d0, 0.000000396566462d0,-0.000000202099617d0, &
&              -0.000000025267769d0, 0.000000013522610d0, 0.000000001229405d0, &
&              -0.000000000821969d0,-0.000000000050847d0, 0.000000000046333d0, &
&               0.000000000001943d0,-0.000000000002563d0,-0.000000000010991d0, &
!                ! ak(2)
&               0.628023320520852d0, 0.462762985338493d0, 0.150044637187895d0, &
&              -0.028796057604906d0,-0.036552745910311d0,-0.001668087945272d0, &
&               0.006519840398744d0, 0.001130378079086d0,-0.000887171310131d0, &
&              -0.000242107641309d0, 0.000087294451594d0, 0.000034682122751d0, &
&              -0.000004583768938d0,-0.000003548684306d0,-0.000000250453880d0, &
&               0.000000216991779d0, 0.000000080779570d0, 0.000000004558555d0, &
&              -0.000000006944757d0,-0.000000002849257d0, 0.000000000237816d0/
      SAVE ak
! ----------------------------------------------------------------------
!
      bk = 0.d0
      dk = 0.d0
!
      xx(1) = 6.d0*nzi*nzi*a*sqrt(stri1)
      xx(2) = 6.d0*nzi*nzj*a*sqrt(stri1)
      xx(3) = 6.d0*nzj*nzj*a*sqrt(stri1)
!
!---  evaluate j-integrals
      DO n=1,3
        xj=xx(n)
        IF (xj > 1.d-30) THEN  ! xj>0.
          IF (xj <= 1.) THEN
            zj=4.d0*xj**(0.2d0)-2.d0
            dzjdxj=0.8d0*xj**(-0.8d0)
            DO k=20,0,-1
              bk(k)=zj*bk(k+1)-bk(k+2)+ak(k,1)
              dk(k)=bk(k+1)+zj*dk(k+1)-dk(k+2)
            END DO
          ELSE
            zj=40.d0/9.d0*xj**(-0.1d0)-22.d0/9.d0
            dzjdxj=-4.d0/9.d0*xj**(-1.1d0)
            DO k=20,0,-1
              bk(k)=zj*bk(k+1)-bk(k+2)+ak(k,2)
              dk(k)=bk(k+1)+zj*dk(k+1)-dk(k+2)             
            END DO
          END IF
          xjxj(n)=0.25d0*xj-1.d0+0.5d0*(bk(0)-bk(2))
          xjprxj(n)=0.25d0+0.5d0*dzjdxj*(dk(0)-dk(2))
        ELSE             
          xjxj(n)=0.d0
          xjprxj(n)=0.d0
        END IF     
      END DO                    ! n-loop
!     
      etetamn   = FLOAT(nzi*nzj)/4.d0/stri1*(xjxj(2)-0.5d0*(xjxj(1) + xjxj(3)))
      etetaprmn = -etetamn/stri1+float(nzi*nzj)/8.d0/stri1/stri1              &
&               * (xx(2)*xjprxj(2)-0.5d0*(xx(1)*xjprxj(1)+xx(3)*xjprxj(3)))
!     
! ----------------------------------------------------------------------
END SUBROUTINE eteta

! ----------------------------------------------------------------------
REAL(RealKind) FUNCTION strength( )
! ----------------------------------------------------------------------
!     calculation of the total ionic strenght
!
! ----------------------------------------------------------------------
!
      INTEGER :: ia, ic
! ----------------------------------------------------------------------
!
      strength = 0.d0
      DO ia=1,na
        strength = strength+cma(ia)*nza(ia)*nza(ia)
      END DO
      DO ic=1,nc
        strength = strength+cmc(ic)*nzc(ic)*nzc(ic)
      END DO      
      strength = 0.5d0*strength
!      
! ----------------------------------------------------------------------
END FUNCTION strength

SUBROUTINE act_coeff(t,p)      

!---  external variables
      REAL(RealKind) :: t,p

!---  external variables
      INTEGER :: ia,ic,icc, iaa,ica,isp,iccp,iaap
      REAL(RealKind) :: rdhl,el, etetpr, etet, xw, phi, xlngpr,    &
&                eps
      REAL(RealKind) :: a, a1a, a2a, str, s1, s2, s3, s4, s5, s6, s7,   &
&                s8, s9, s10, s11, s11a, s12, s12a, s13, s13a , s22

      REAL(RealKind), PARAMETER :: a1=2.d0,a2=0.d0,b=1.2d0
      REAL(RealKind), PARAMETER :: cmw=1000.d0/18.015d0

!     OSSI
      REAL(RealKind) :: t1
      REAL(RealKind) :: ga1(na),gc1(nc)
      REAL(RealKind) :: ga2(na),gc2(nc)
      REAL(RealKind) :: ga3(na),gc3(nc)
      REAL(RealKind) :: ga4(na),gc4(nc)
!
! ----------------------------------------------------------------------
! --- definition of internal functions
      REAL(RealKind) :: f1, f2, bmx, bmxpr, bmxaw, nue
      f1(a,str)=1.-(1.+a*SQRT(str))*EXP(-a*SQRT(str))
      f2(a,str)=(1+a*SQRT(str)+0.5*a*a*str)*EXP(-a*SQRT(str))-1.
      bmx(ic,ia,a1a,a2a,str)=b0(ic,ia)+2.0*b1(ic,ia)/a1a/a1a/str    &
&                           *f1(a1a,str)                            
!    &                    +2.0*b2(ic,ia)/a2a/a2a*stri*f1(a2a,str)
      bmxpr(ic,ia,a1a,a2a,str)=2.0*b1(ic,ia)/a1a/a1a/str/str        &
&                             *f2(a1a,str)                          
!    &                    +2.0*b2(ic,ia)/a2a/a2a*str*str*f2(a2a,str)
      bmxaw(ic,ia,a1a,str)=b0(ic,ia)+b1(ic,ia)*EXP(-a1a*SQRT(str))
      nue(ic,ia)=nuec(ic,ia)+nuea(ic,ia)
! ----------------------------------------------------------------------
!
! --- ionic strength, etc.

      a=aphi(t,p,eps,rdhl,el)
!     stri=strength(nzc,nza,cmc,cma) !OSSI
      stri=strength()
      rdhl=sqrt(2.*rdhl*stri)
!
! ----------------------------------------------------------------------
! ---- loop over the different salts
      DO ic=1,nc
         DO ia=1,na
!             
           xlngpr=-a*( SQRT(stri)/(1.+b*SQRT(stri))+2./b * LOG(1.+b*SQRT(stri)) )
!             
! ----- 1st term
           s1=0.
           s11=0.
           DO icc=1,nc
             s11=s11+cmc(icc)*nzc(icc)
           END DO
           DO iaa=1,na
             s11=s11-cma(iaa)*nza(iaa)
           END DO
           DO iaa=1,na
             t1=bmx(ic,iaa,a1,a2,stri)
             s1=s1+cma(iaa)*(2.*bmx(ic,iaa,a1,a2,stri)+s11*c(ic,iaa))
           END DO
!             
! ----- 2nd term
           s2=0.
           DO icc=1,nc
             s11=0.
             DO ica=1,nc
               s11=s11+cmc(ica)*nzc(ica)
             END DO
             DO ica=1,na
               s11=s11-cma(ica)*nza(ica)
             END DO
             s2=s2+cmc(icc)*(2.*bmx(icc,ia,a1,a2,stri)+s11*c(icc,ia))
           END DO
!             
! ----- 3rd, & 4th term
           s3=0.
           s4=0.
           DO icc=1,nc
             DO iaa=1,na
                s3 = s3+cmc(icc)*cma(iaa) *     &
&                    (bmxpr(icc,iaa,a1,a2,stri)*nzc(ic)*nzc(ic)+ABS(nzc(ic))*c(icc,iaa))
                s4 = s4+cmc(icc)*cma(iaa) *     &
&                  (nza(ia)*nza(ia)*bmxpr(icc,iaa,a1,a2,stri)+ABS(nza(ia))*c(icc,iaa))
             END DO
           END DO
!             
! ----- 5th, & 6th term
           s5=0.
           s6=0.
           DO icc=1,nc
             s11=0.
             DO iaa=1,na
               s11=s11+cma(iaa)*psim(ic,icc,iaa) ! (ic,ia,iaa)
             END DO
             s5=s5+cmc(icc)*(2.*tetam(ic,icc)+s11)
           END DO
!                      
           DO iaa=1,na
             s22=0.
             DO icc=1,nc
               s22=s22+cmc(icc)*psix(ia,iaa,icc) ! (icc,ia,iaa)
             END DO
             s6=s6+cma(iaa)*(2.*tetax(ia,iaa)+s22)
           END DO
!                      
! ----- 7th, & 8th term
           s7=0.
           s8=0.
           DO iaa=1,na-1
             DO iaap=iaa,na
               s7=s7+cma(iaa)*cma(iaap)*psix(iaa,iaap,ic)
             END DO
           END DO
           DO icc=1,nc-1
             DO iccp=icc,nc
               s8=s8+cmc(icc)*cmc(iccp)*psim(icc,iccp,ia)
             END DO
           END DO
!                      
           IF (nzc(ic).ne.0) THEN
             gm(ic)=EXP(nzc(ic)*nzc(ic)*xlngpr+s1+s3+s5+s7)
             gc1(ic)=nzc(ic)*nzc(ic)*xlngpr
             gc2(ic)=nzc(ic)*nzc(ic)*xlngpr+s1+s3
             gc3(ic)=s5/2.0d0
             gc4(ic)=nzc(ic)*nzc(ic)*xlngpr+s1+s3+s5+s7
           END IF 
           IF (nza(ia).ne.0) THEN
             gx(ia)=EXP(nza(ia)*nza(ia)*xlngpr+s2+s4+s6+s8)
             ga1(ia)=nza(ia)*nza(ia)*xlngpr
             ga2(ia)=nza(ia)*nza(ia)*xlngpr+s2+s4
             ga3(ia)=s6/2.0d0
             ga4(ia)=nza(ia)*nza(ia)*xlngpr+s2+s4+s6+s8
           END IF 
           IF (nue(ic,ia).ne.0)                                           &
&             gmx(ic,ia)=exp ( 1./nue(ic,ia)*( nuec(ic,ia)*LOG(gm(ic))+   &
&             nuea(ic,ia)*LOG(gx(ia)) ))
         END DO
      END DO
!             
! --- water activity
! --- 9th, 10th & part of the 11th term
      s9=0.
      s11a=0.
      DO ia=1,na
        s9=s9+cma(ia)
        s11a=s11a+cma(ia)*ABS(nza(ia))
      END DO 
      DO ic=1,nc
        s9=s9+cmc(ic)
        s11a=s11a+cmc(ic)*nzc(ic)
      END DO
      s10=-a*stri**1.5/(1.+b*sqrt(stri))
      t1=-2.0d0*s10 
!     
! --- 11th term
      s11=0.
      DO ic=1,nc
        DO ia=1,na
          s11=s11+cmc(ic)*cma(ia)*(bmxaw(ic,ia,a1,stri)+s11a*c(ic,ia))
!FM          s11=s11+cmc(ic)*cma(ia)*(bmx(ic,ia,a1,a2,stri)+s11a*c(ic,ia))
        END DO 
      END DO
!     
! --- 12th term
      s12=0.
      DO ic=1,nc-1
        DO icc=ic+1,nc
          s12a=0.
          DO ia=1,na
            s12a=s12a+cma(ia)*psim(ic,icc,ia)
          END DO
!         correct for unsymmetrical mixing
          IF (nzc(icc)*nzc(ic).gt.0.and.nzc(icc).ne.nzc(ic)) THEN
            call eteta (a,stri,nzc(icc),nzc(ic),etet,etetpr)
          ELSE
            etet=0.
            etetpr=0.
          END IF 
!         write (6,'(2i3,1p,5e13.5)') ic,icc,stri,tetam(ic,icc)
!    &         ,tetam(ic,icc)+etet+stri*etetpr,etet+stri*etetpr
!
          s12=s12+cmc(ic)*cmc(icc)*(tetam(ic,icc)+etet+stri*etetpr+s12a)
        END DO 
      END DO     
!     
! --- 13th term
      s13=0.
      DO ia=1,na-1
        DO iaa=ia+1,na
          s13a=0.
          DO ic=1,nc
            s13a=s13a+cmc(ic)*psix(ia,iaa,ic)
          END DO
!         correct for unsymmetrical mixing
          IF (nza(iaa)*nza(ia).gt.0.and.nza(iaa).ne.nza(ia)) THEN
            call eteta (a,stri,nza(iaa),nza(ia),etet,etetpr)
          ELSE
            etet=0.
            etetpr=0.
          END IF 
!         write (6,'(2i3,1p,5e13.5)') ia,iaa,stri,tetax(ia,iaa)
!    &         ,tetax(ia,iaa)+etet+stri*etetpr,etet+stri*etetpr
!
          s13=s13+cma(ia)*cma(iaa)*(tetax(ia,iaa)+etet+stri*etetpr+s13a)
        END DO 
      END DO     
!     
! ----------------------------------------------------------------------
! --- water activity
      phi=1.+2./s9*(s10+s11+s12+s13)
      aw=exp (-s9/cmw*phi)
!
! ----------------------------------------------------------------------
! --- water mixing ratio
!     xw=cmw/(cmw+s9)
!     fw=aw/xw
!     
! ----------------------------------------------------------------------
END SUBROUTINE act_coeff 

! ------------------------------------------------------------------------
SUBROUTINE InitPitzer(nPitz)
! ------------------------------------------------------------------------
!
!    nPitz - flag for the input data set;
!            1:  Judit
!            2:  harvie (1984)
!            3:  pitzer (1973)
!
! -------------------------------------------------------------------------

  INTEGER :: nPitz

!---  internal variables
  INTEGER :: ia, ic, ios, jt, indspc, na1, nc1, pos
  CHARACTER(20):: string

!---  external function
      INTEGER :: ifind
      EXTERNAL   ifind

! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
!---  Define index array PitzInd, array for activity coefficients
      ALLOCATE (PitzInd(nAqua))

      PitzInd(:)      = 0
!
!---  Read Pitzer tables
      WRITE(*,801)
      call init_act(nPitz)

! -------------------------------------------------------------------------
!---  Set index transformation 
! -------------------------------------------------------------------------
!
      na1 = 0
      nc1 = 0
!
      PitzInd(:) = 0 
      DO jt=1,nc
         indspc=Position(mname(jt))
         IF (indspc > 0) THEN
            PitzInd(indspc)=jt
         END IF
      END DO
      DO jt=1,na
         indspc=Position(xname(jt))
         IF (indspc>0) THEN
            PitzInd(indspc)=-jt
         END IF
      END DO

! -------------------------------------------------------------------------
!---  Print Pitzer system
! -------------------------------------------------------------------------
!
      WRITE(*,802) nc1,na1, '  No.','  PitzInd','   Charge','     Species Name'
      DO jt=1,nAqua
         IF (PitzInd(jt) <= 0)  CYCLE
         WRITE(*,803)  jt, PitzInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      DO jt=1,nAqua
         IF (PitzInd(jt) >= 0)  CYCLE
         WRITE(*,803)  jt, PitzInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      WRITE(*,804)

! -------------------------------------------------------------------------
801   FORMAT(1x/1x,75('=')/ ' Pitzer Initialization:' )
802   FORMAT(1x/1x,75('-')/           &
&            ' Considered Pitzer System: Cations =',i3,'     Anions =',i3 // &
&            a5,a9,a9,a18 / 1x,40('-'))
803   FORMAT(i4,i9,f9.1,7x,a18)
804   FORMAT(1x/1x,75('=')/1x)
!
! -------------------------------------------------------------------------
END SUBROUTINE InitPitzer

! ------------------------------------------------------------------------
SUBROUTINE init_act (nPitz)
! ------------------------------------------------------------------------
!
!    nPitz - flag for the input data set;
!            1:  Judit
!            2:  harvie (1984)
!            3:  pitzer (1973)
!
! -------------------------------------------------------------------------
!---  exinternal variables
      INTEGER :: nPitz

!---  internal variables
      INTEGER :: ia, ic, ios
      INTEGER :: ir0=9, ir1=10, ir2=11, ir3=12

      CHARACTER(7)  ::  set
      CHARACTER(12) ::  Path
!
! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
      IF (nPitz <= 0) RETURN
      IF (nPitz == 1) THEN
         write (6,'(/a/)') 'use standard pitzer input data...'
         path = 'Data_IfT/'
         set  = '.judit'
      ELSE IF (nPitz == 2) THEN
         write (6,'(/a/)') 'use harvie & standard pitzer input data...'
         path = 'Data_Harvie/'
         set  = '.harvie'
      ELSE IF (nPitz == 3) THEN
         write (6,'(/a/)') 'use standard pitzer input data...'
         path = 'Data_Pitzer/'
         set  = '.pitzer'
      END IF
!
!------------------------------------------------------------------
! --- Read Pitzer System: Dimensions, Names, Ions, Charges
!------------------------------------------------------------------
!
      OPEN (ir0,FILE=TRIM(path)//'name'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'name'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Name= File !!'
      END IF

!---  read and set dimensions
      READ(ir0,*)
      READ(ir0,*)  nc, na

      nac  = na*nc
      nacc = na * nac
      naca = nc * nac
      nc2  = nc * nc
      na2  = na * na

!---  allocate arrays
      ALLOCATE ( nzc(nc),nza(na),nuec(nc,na),nuea(nc,na) )
      nza(:) = 0
      nzc(:) = 0
      nuea(:,:) = 0
      nuec(:,:) = 0

      ALLOCATE ( b0(nc,na),b1(nc,na),c(nc,na) )
      ALLOCATE ( psim(nc,nc,na),tetax(na,na) )
      ALLOCATE ( psix(na,na,nc),tetam(nc,nc) )
      b0(:,:) = 0.d0
      b1(:,:) = 0.d0
      c(:,:)  = 0.d0
      psim(:,:,:) = 0.d0
      psix(:,:,:) = 0.d0
      tetam(:,:)  = 0.d0
      tetax(:,:)  = 0.d0

      ALLOCATE ( gmx(nc,na),gm(nc),gx(na) )
      ALLOCATE ( cmc(nc),cma(na) )
      gm(:)  = 0.d0
      gx(:)  = 0.d0
      cma(:) = 0.d0
      cmc(:) = 0.d0
      gmx(:,:) = 0.d0

      ALLOCATE ( mname(nc),xname(na) )
      ALLOCATE ( mxname(nc,na) )

!---  read names
      REWIND (ir0)
      CALL read_name(ir0)
      CLOSE (ir0)
!
!------------------------------------------------------------------
! --- determine the stoichiometric coefficients
      DO ic=1,nc
         DO ia=1,na
            nuec(ic,ia) = ABS(nza(ia))
            nuea(ic,ia) = nzc(ic)
         END DO
      END DO
!
!------------------------------------------------------------------
! --- Read ion interaction characteristics: Pitzer coefficients
!------------------------------------------------------------------
!
! ---  B0, B1, C
      OPEN (ir1,FILE=TRIM(path)//'b01c'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'b01c'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =B01C= File !!'
      END IF

!---  read values
      CALL read_b01c(ir1)
      CLOSE (ir1)
!
!------------------------------------------------------------------
! --- convert c's according to their charge
      DO ic=1,nc
         DO ia=1,na
            c(ic,ia) = 0.5d0*c(ic,ia)/sqrt(float(abs(nza(ia)*nzc(ic))))
         END DO
      END DO
!
!------------------------------------------------------------------
! ---  TetaM, TetaX
      OPEN (ir2,FILE=TRIM(path)//'teta'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'teta'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Teta= File !!'
      END IF

      CALL read_teta(ir2)
      CLOSE (ir2)
!
!------------------------------------------------------------------
! ---  PsiM, PsiX
      OPEN (ir3,FILE=TRIM(path)//'psinew'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'psinew'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Psi= File !!'
      END IF

      CALL read_psi1(ir3)
      CLOSE (ir3)
!
! -------------------------------------------------------------------------
END SUBROUTINE init_act

! ----------------------------------------------------------------------
SUBROUTINE read_psi (iread)
! ----------------------------------------------------------------------
!--- external variables
      INTEGER :: iread
!
!--- internal variables
      INTEGER :: ia, ia1, ia2, ic, ic1, ic2
! ----------------------------------------------------------------------
!
      read (iread,*)
      DO ic1=1,nc
      DO ic2=1,nc
        read (iread,*) (psim(ic1,ic2,ia),ia=1,na)
      END DO
      END DO
      read (iread,*)
      read (iread,*)
      DO ia1=1,na
      DO ia2=1,na
        read (iread,*) (psix(ia1,ia2,ic),ic=1,nc)
      END DO
      END DO
!
      close (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_psi 
!
! ======================================================================
SUBROUTINE read_psi1 (iread)
! ----------------------------------------------------------------------
!--- external variables
      INTEGER :: iread
!
!--- internal variables
      INTEGER :: ia, ia1, ia2, ic, ic1, ic2
! ----------------------------------------------------------------------
!
      WRITE(*,801)  nc,mname
      WRITE(*,802)  na,xname
801   FORMAT(' Cations:',I3,/(11x,5a14))
802   FORMAT(' Anions :',I3,/(11x,5a14))
     
      READ (iread,*)
      DO ic1=1,nc-1
        DO ic2=ic1+1,nc
          READ (iread,*) (psim(ic1,ic2,ia),ia=1,na)
          DO ia=1,na
            psim(ic2,ic1,ia) = psim(ic1,ic2,ia)
          END DO 
        END DO
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ia1=1,na-1
        DO ia2=ia1+1,na
          READ (iread,*) (psix(ia1,ia2,ic),ic=1,nc)
          DO ic=1,nc
            psix(ia2,ia1,ic) = psix(ia1,ia2,ic)
          END DO 
        END DO
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_psi1
! ----------------------------------------------------------------------
SUBROUTINE read_name (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread
!
!---  internal variables
      INTEGER :: ia, ic, nc1, na1
! ----------------------------------------------------------------------

      REWIND (iread)
      READ (iread,*)
      READ (iread,*) nc1,na1
      IF (nc1 /= nc .OR. na1 /= na) THEN
        print*,'dimensions of the parameter space does not coincident'  &
&              //' with the model dimenssions!'
        print*,'nc1: ',nc1,'  nc: ',nc
        print*,'na1: ',na1,'  na: ',na
        stop 'READ_nmae...'
      END IF

      READ (iread,*)
      READ (iread,*) (mname(ic),ic=1,nc)
      READ (iread,*) (nzc(ic),ic=1,nc)
      READ (iread,*)
      READ (iread,*) (xname(ia),ia=1,na)
      READ (iread,*) (nza(ia),ia=1,na)
      READ (iread,*)
      DO ic=1,nc
        READ (iread,*) (mxname(ic,ia),ia=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_name
! ----------------------------------------------------------------------
SUBROUTINE read_b01c (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread
!
!---  internal variables
      INTEGER :: ia, ic

! ----------------------------------------------------------------------
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (b0(ic,ia),ia=1,na)
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (b1(ic,ia),ia=1,na)
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (c(ic,ia),ia=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_b01c 
! ----------------------------------------------------------------------
SUBROUTINE read_teta (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread

!---  internal variables
      INTEGER :: ia1, ic1, ia2, ic2
! ----------------------------------------------------------------------
!
      READ (iread,*)
      DO ic1=1,nc
        READ (iread,*) (tetam(ic1,ic2),ic2=1,nc)
      END DO
      READ (iread,*)
!
      READ (iread,*)
      DO ia1=1,na
        READ (iread,*) (tetax(ia1,ia2),ia2=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_teta

FUNCTION Position(Species)

  INTEGER :: Position
  CHARACTER(*) :: Species

  INTEGER :: i

  Position=0
  DO i=1,nAqua
    IF (TRIM(Species)==TRIM(SpeciesNameAqua(i))) THEN
      Position=i
      EXIT
    END IF
  END DO

END FUNCTION

END MODULE Pitzer_Mod

PROGRAM TPitzer

  USE Pitzer_Mod
  
  IMPLICIT NONE

  INTEGER :: i
  REAL(RealKind) :: SumA,SumC
  REAL(RealKind) :: TempLoc,Pre
  REAL(RealKind), ALLOCATABLE :: xC(:)

  nAqua=17
  ALLOCATE(Charge(nAqua))
  ALLOCATE(SpeciesNameAqua(nAqua))
  SpeciesNameAqua(1)='NAp'
  Charge(1)=1
  SpeciesNameAqua(2)='NH4p'
  Charge(2)=1
  SpeciesNameAqua(3)='Hp'
  Charge(3)=1
  SpeciesNameAqua(4)='Kp'
  Charge(4)=1
  SpeciesNameAqua(5)='CApp'
  Charge(5)=2
  SpeciesNameAqua(6)='MGpp'
  Charge(6)=2
  SpeciesNameAqua(7)='FEpp'
  Charge(7)=2
  SpeciesNameAqua(8)='MNpp'
  Charge(8)=2
  SpeciesNameAqua(9)='CUpp'
  Charge(9)=2
  SpeciesNameAqua(10)='CLm'
  Charge(10)=-1
  SpeciesNameAqua(11)='SO4mm'
  Charge(11)=-2
  SpeciesNameAqua(12)='OHm'
  Charge(12)=-1
  SpeciesNameAqua(13)='NO3m'
  Charge(13)=-1
  SpeciesNameAqua(14)='HSO4m'
  Charge(14)=-1
  SpeciesNameAqua(15)='HCO3m'
  Charge(15)=-1
  SpeciesNameAqua(16)='CO3mm'
  Charge(16)=-2
  SpeciesNameAqua(17)='BRm'
  Charge(17)=-1
 
  CALL InitPitzer(1)
  ALLOCATE(xC(nAqua))
  SumA=0.0d0
  SumC=0.0d0
  DO i=1,nA
    xC(i)=i
    SumA=SumA+Charge(i)*xC(i)
  END DO
  DO i=1,nA
    xC(i)=xC(i)/SumA
  END DO
  DO i=nA+1,nAqua
    xC(i)=i
    SumC=SumC-Charge(i)*xC(i)
  END DO
  DO i=nA+1,nAqua
    xC(i)=xC(i)/SumC
  END DO
  SumA=0.0d0
  DO i=1,nAqua
    SumA=SumA+Charge(i)*xC(i)
  END DO

  DO i=1,nAqua
    IF (PitzInd(i)>0) THEN
      cmc(PitzInd(i))=xC(i)
    ELSE IF (PitzInd(i)<0) THEN
      cma(-PitzInd(i))=xC(i)
    END IF
  END DO
  WRITE(*,*) cma
  WRITE(*,*) cmc
  gm(:)=0.d0
  gx(:)=0.d0
  gmx(:,:)=0.d0
  TempLoc=289.15D0
  Pre=1000.0d0
! Null setzen
  
  CALL act_coeff (TempLoc,Pre)
  WRITE(*,*) 'gm',gm
  WRITE(*,*) 'gx',gx
  WRITE(*,*) 'aw',aw
  
END PROGRAM TPitzer
