SUBROUTINE ReadWeights(InputFile)

  USE Kind_Mod
  USE Floor_Mod

  IMPLICIT NONE

  CHARACTER*80 :: InputFile
  CHARACTER(300) :: Line
  INTEGER :: ib,i,j,k
  INTEGER :: ix,iy,iz
  INTEGER :: ixref, iyref, izref
  INTEGER :: levX,levY,levZ,shX,shY,shZ
  INTEGER :: levNX,levNY,levNZ,shNX,shNY,shNZ
  INTEGER :: InputUnit=10
  INTEGER :: nxx,nyy,nzz,nlines
  REAL(8) :: xw0,xw1,yw0,yw1,zw0,zw1
 
  OPEN(UNIT=InputUnit,FILE=TRIM(InputFile),STATUS='old')

  ALLOCATE(domain%dx(domain%ix0+1:domain%ix1))
  ALLOCATE(domain%dy(domain%iy0+1:domain%iy1))
  ALLOCATE(domain%dz(domain%iz0+1:domain%iz1))
  ALLOCATE(domain%xP(domain%ix0:domain%ix1))
  ALLOCATE(domain%yP(domain%iy0:domain%iy1))
  ALLOCATE(domain%zP(domain%iz0:domain%iz1))

  DO 
    READ(InputUnit,*,END=1) Line
    !.......................................................................
    IF (INDEX(Line,'#xGrid')>0) THEN   ! set different distance x-direction
      READ(InputUnit,*) nlines
      nxx=0
      domain%xP(domain%ix0)=domain%x0
      DO i=1,nlines
        READ(InputUnit,*) xw0,xw1,domain%nx
        DO ix=nxx+1,nxx+domain%nx
          domain%dx(ix)=(xw1-xw0)/domain%nx
          domain%xP(ix)=domain%xP(ix-1)+domain%dx(ix)
        END DO
        nxx=nxx+domain%nx
      END DO
      domain%nx=nxx
      EXIT
    ELSE
      domain%xP(domain%ix0)=domain%x0
      DO ix=domain%ix0+1,domain%ix1
        domain%dx(ix)=(domain%x1-domain%x0)/domain%nx
        domain%xP(ix)=domain%xP(ix-1)+domain%dx(ix)
      END DO
    END IF
  END DO
1 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=2) Line
    !.......................................................................
    IF (INDEX(Line,'#yGrid')>0) THEN   ! set different distance y-direction
      READ(InputUnit,*) nlines
      nyy=0
      domain%yP(domain%iy0)=domain%y0
      DO j=1,nlines
        READ(InputUnit,*) yw0,yw1,domain%ny
        DO iy=nyy+1,nyy+domain%ny
          domain%dy(iy)=(yw1-yw0)/domain%ny
          domain%yP(iy)=domain%yP(iy-1)+domain%dy(iy)
        END DO
        nyy=nyy+domain%ny
      END DO
      domain%ny=nyy
      EXIT
    ELSE
      domain%yP(domain%iy0)=domain%y0
      DO iy=domain%iy0+1,domain%iy1
        domain%dy(iy)=(domain%y1-domain%y0)/domain%ny
        domain%yP(iy)=domain%yP(iy-1)+domain%dy(iy)
      END DO
    END IF
  END DO
2 CONTINUE
  REWIND(InputUnit)
  DO 
    READ(InputUnit,*,END=3) Line
    !........................................................................
    IF (INDEX(Line,'#zGrid')>0) THEN   ! set different distance z-direction
      READ(InputUnit,*) nlines
      nzz=0
      domain%zP(domain%iz0)=domain%z0
      DO k=1,nlines
        READ(InputUnit,*) zw0,zw1,domain%nz
        DO iz=nzz+1,nzz+domain%nz
          domain%dz(iz)=(zw1-zw0)/domain%nz
          domain%zP(iz)=domain%zP(iz-1)+domain%dz(iz)
        END DO
        nzz=nzz+domain%nz
      END DO
      domain%nz=nzz
      EXIT
    ELSE
      domain%zP(domain%iz0)=domain%z0
      DO iz=domain%iz0+1,domain%iz1
        domain%dz(iz)=(domain%z1-domain%z0)/domain%nz
        domain%zP(iz)=domain%zP(iz-1)+domain%dz(iz)
      END DO
    END IF
    !........................................................................
  END DO
3 CONTINUE
  CLOSE(UNIT=InputUnit)
  !..........................................................................
  ! Verschoben-Kopie  InitAnalyzeArea
  domain%x0View=domain%xP(domain%ix0)
  domain%x1View=domain%xP(domain%ix1)
  domain%y0View=domain%yP(domain%iy0)
  domain%y1View=domain%yP(domain%iy1)
  domain%z0View=domain%zP(domain%iz0)
  domain%z1View=domain%zP(domain%iz1)
  !..........................................................................

  DO ib=1,nb
    CALL Set(Floor(ib))
    levX=-RefineX
    shX=MAX(1,2**levX)-1
    levY=-RefineY
    shY=MAX(1,2**levY)-1
    levZ=-RefineZ
    shZ=MAX(1,2**levZ)-1
    !......................................
    ix=ix0+1
    ixref = ix * 2.e0**levX
    xP(ix0)=domain%xP(ixref-shX-1)
    DO ix=ix0+1,ix1
      ixref = ix * 2.e0**levX
      dx(ix)=SUM(domain%dx(ixref-shX:ixref))
      xP(ix)=xP(ix-1)+dx(ix)
    END DO
    dx(ix0)=dx(ix0+1)  !Rand
    dx(ix1+1)=dx(ix1)
    xP(ix0-1)=xP(ix0)-dx(ix0)
    xP(ix1+1)=xP(ix1)+dx(ix1)
    !......................................
    iy=iy0+1
    iyref = iy * 2.e0**levY
    yP(iy0)=domain%yP(iyref-shY-1)
    DO iy=iy0+1,iy1
      iyref = iy * 2.e0**levY
      dy(iy)=SUM(domain%dy(iyref-shY:iyref))
      yP(iy)=yP(iy-1)+dy(iy)
    END DO
    dy(iy0)=dy(iy0+1)  !Rand
    dy(iy1+1)=dy(iy1)
    yP(iy0-1)=yP(iy0)-dy(iy0)
    yP(iy1+1)=yP(iy1)+dy(iy1)
    !......................................
    iz=iz0+1
    izref = iz * 2.e0**levZ
    zP(iz0)=domain%zP(izref-shZ-1)
    DO iz=iz0+1,iz1
      izref = iz * 2.e0**levZ
      dz(iz)=SUM(domain%dz(izref-shZ:izref))
      zP(iz)=zP(iz-1)+dz(iz)
    END DO
    dz(iz0)=dz(iz0+1)  !Rand
    dz(iz1+1)=dz(iz1)
    zP(iz0-1)=zP(iz0)-dz(iz0)
    zP(iz1+1)=zP(iz1)+dz(iz1)
  END DO

END SUBROUTINE ReadWeights

