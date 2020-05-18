PROGRAM Grid_Built
!  Programm, um die Bloecke im Grid-File zu erzeugen.

  IMPLICIT NONE
  INTEGER :: NX,NY,NZ,BX,BY,BZ,xi,yi,zi
  INTEGER :: xRes,yRes,zRes
  INTEGER :: xIncr,yIncr,zIncr
  INTEGER :: ix,ix0,ix1
  INTEGER :: iy,iy0,iy1
  INTEGER :: iz,iz0,iz1
  INTEGER :: InputUnit
  INTEGER :: OutputUnit
  CHARACTER :: lc
  CHARACTER(len = 100) :: mls = ''

  TYPE Block_T
    INTEGER :: ix0,ix1,xC=1
    INTEGER :: iy0,iy1,yC=1
    INTEGER :: iz0,iz1,zC=1
  END TYPE Block_T

  TYPE(Block_T), ALLOCATABLE :: Blocks(:,:,:)
  INTEGER :: xCMax,yCMax,zCMax
  INTEGER :: xBMax,yBMax,zBMax

  CHARACTER(80) :: InputFileName
  CHARACTER(300) :: Line
  INTEGER :: i,j,k,iF
  INTEGER :: NumFineBlocks
  REAL(8) :: xC,yC,zC

  InputUnit=98
  OutputUnit=99
  CALL GET_COMMAND_ARGUMENT(1,InputFileName)
  OPEN(UNIT=InputUnit,FILE=InputFileName,STATUS='UNKNOWN')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'#SelfMultiblock')>0) THEN
      READ(InputUnit,*) ix0,ix1,BX
      READ(InputUnit,*) iy0,iy1,BY
      READ(InputUnit,*) iz0,iz1,BZ
      NX=ix1-ix0
      NY=iy1-iy0
      NZ=iz1-iz0
      ALLOCATE(Blocks(BX,BY,BZ))
      READ(InputUnit,*) NumFineBlocks
      READ(InputUnit,*) xCMax,yCMax,zCMax
      READ(InputUnit,*) xBMax,yBMax,zBMax
      DO iF=1,NumFineBlocks
        READ(InputUnit,*) i,j,k,xC,yC,zC
        Blocks(i,j,k)%xC=xC
        Blocks(i,j,k)%yC=yC
        Blocks(i,j,k)%zC=zC
      END DO
      EXIT
    END IF  
  END DO
1 CONTINUE  

  OPEN(UNIT=OutputUnit,FILE='Block.grid',STATUS='REPLACE')

! WRITE(*,*)'X - Domain'
! READ(*,*) NX
! WRITE(*,*)'X - blocks'
! READ(*,*) BX
! WRITE(*,*)'Y - Domain'
! READ(*,*) NY
! WRITE(*,*)'Y - blocks'
! READ(*,*) BY
! WRITE(*,*)'Z - Domain'
! READ(*,*) NZ
! WRITE(*,*)'Z - blocks'
! READ(*,*) BZ
! ALLOCATE(Blocks(BX,BY,BZ))
  WRITE(*,*)'MultiLayerSoil,example:"1 2 3 4 5 6 7 8"'
  WRITE(*,*)'1=ice; 2=rock; 3=sand; 4=sandy loam'
  WRITE(*,*)'5=loam; 6=clay loam; 7=clay; 8=peat'
  READ(*,'(A)') mls
  WRITE(*,*)'LandClass'
  WRITE(*,*)'0=bare soil; 1=urban area; 2=savannah; 3=deciduous forest;'
  WRITE(*,*)'4=coniferous forest; 5=mixed forest; 6=shrubland;'
  WRITE(*,*)'7=annual land; 8=grass land; 9=sea'
  READ(*,'(A)') lc

  xRes=MOD(NX,BX)
  xIncr=NX/BX
  yRes=MOD(NY,BY)
  yIncr=NY/BY
  zRes=MOD(NZ,BZ)
  zIncr=NZ/BZ

  WRITE(OutputUnit,'(a11)') '#Multiblock'
  WRITE(OutputUnit,*) 0,NX,'     ! X ... nxa, nxe'
  WRITE(OutputUnit,*) 0,NY,'     ! Y ... nya, nye'
  WRITE(OutputUnit,*) 0,NZ,'     ! Z ... nza, nze'
  WRITE(OutputUnit,*) ''
  WRITE(OutputUnit,*) BX*BY*BZ, '       !Number of blocks ... nb'

  iz0=0
  DO iz=1,BZ
    IF (iz<=zRes) THEN
      iz1=iz0+zIncr+1
    ELSE  
      iz1=iz0+zIncr
    END IF  
    iy0=0
    DO iy=1,BY
      IF (iy<=yRes) THEN
        iy1=iy0+yIncr+1
      ELSE  
        iy1=iy0+yIncr
      END IF  
      ix0=0 
      DO ix=1,BX
        IF (ix<=xRes) THEN
          ix1=ix0+xIncr+1
        ELSE  
          ix1=ix0+xIncr
        END IF  
        WRITE(OutputUnit,*)'                             !Block',ix+(iy-1)*BX+(iz-1)*BX*BY
        WRITE(OutputUnit,*) ix0,ix1,'     ! ixa, ixe'
        WRITE(OutputUnit,*) iy0,iy1,'     ! iya, iye'
        WRITE(OutputUnit,*) iz0,iz1,'     ! iza, ize'
        WRITE(OutputUnit,*)'0 0 0 0 4 T T'
        Blocks(ix,iy,iz)%ix0=ix0
        Blocks(ix,iy,iz)%ix1=ix1
        Blocks(ix,iy,iz)%iy0=iy0
        Blocks(ix,iy,iz)%iy1=iy1
        Blocks(ix,iy,iz)%iz0=iz0
        Blocks(ix,iy,iz)%iz1=iz1
        ix0=ix1
      END DO  
      iy0=iy1
    END DO  
    iz0=iz1
  END DO  
  WRITE(OutputUnit,'(a15)') '#MultiLayerSoil'
  WRITE(OutputUnit,*) BX*BY*BZ, '       !Number of blocks ... nb'

  iz0=0
  DO iz=1,BZ
    IF (iz<=zRes) THEN
      iz1=iz0+zIncr+1
    ELSE  
      iz1=iz0+zIncr
    END IF  
    iy0=0
    DO iy=1,BY
      IF (iy<=yRes) THEN
        iy1=iy0+yIncr+1
      ELSE  
        iy1=iy0+yIncr
      END IF  
      ix0=0 
      DO ix=1,BX
        IF (ix<=xRes) THEN
          ix1=ix0+xIncr+1
        ELSE  
          ix1=ix0+xIncr
        END IF  
        WRITE(OutputUnit,*) '!.............................................'
        WRITE(OutputUnit,*) ix+(iy-1)*BX+(iz-1)*BX*BY,'                    !Block',ix+(iy-1)*BX+(iz-1)*BX*BY
        WRITE(OutputUnit,*) '1                ! nr_soildef'
        WRITE(OutputUnit,*) ix0,ix1,'     ! ixa, ixe'
        WRITE(OutputUnit,*) iy0,iy1,'     ! iya, iye'
        WRITE(OutputUnit,*) iz0,iz1,'     ! iza, ize'
        WRITE(OutputUnit,*) TRIM(mls),'   ! soil layers'
        ix0=ix1
      END DO  
      iy0=iy1
    END DO  
    iz0=iz1
  END DO  
  WRITE(OutputUnit,'(a13)') '#LandClassDef'
  WRITE(OutputUnit,*) BX*BY*BZ, '       !Number of blocks ... nb'

  iz0=0
  DO iz=1,BZ
    IF (iz<=zRes) THEN
      iz1=iz0+zIncr+1
    ELSE  
      iz1=iz0+zIncr
    END IF  
    iy0=0
    DO iy=1,BY
      IF (iy<=yRes) THEN
        iy1=iy0+yIncr+1
      ELSE  
        iy1=iy0+yIncr
      END IF  
      ix0=0 
      DO ix=1,BX
        IF (ix<=xRes) THEN
          ix1=ix0+xIncr+1
        ELSE  
          ix1=ix0+xIncr
        END IF  
        WRITE(OutputUnit,*) '!.............................................'
!       WRITE(OutputUnit,*) '                             !Block',ix+(iy-1)*BX+(iz-1)*BX*BY
        WRITE(OutputUnit,*) ix0,ix1,'     ! ixa, ixe'
        WRITE(OutputUnit,*) iy0,iy1,'     ! iya, iye'
        WRITE(OutputUnit,*) lc,'              ! LandClass'
        ix0=ix1
      END DO  
      iy0=iy1
    END DO  
    iz0=iz1
  END DO  
  CLOSE(OutputUnit)
  CLOSE(InputUnit)
! OPEN(UNIT=OutputUnit,FILE='BlockNeu.grid',STATUS='REPLACE')
! WRITE(OutputUnit,'(a11)') '#Multiblock'
! WRITE(OutputUnit,*) 0,NX,'     ! X ... nxa, nxe'
! WRITE(OutputUnit,*) 0,NY,'     ! Y ... nya, nye'
! WRITE(OutputUnit,*) 0,NZ,'     ! Z ... nza, nze'
! WRITE(OutputUnit,*) ''
! WRITE(OutputUnit,*) BX*BY*BZ, '       !Number of blocks ... nb'
  CALL Coarse(Blocks,xCMax,yCMax,zCMax)
! CLOSE(OutputUnit)
CONTAINS
SUBROUTINE Coarse(Blocks,xCMax,yCMax,zCMax)
  TYPE(Block_T) :: Blocks(:,:,:)
  INTEGER :: xCMax,yCMax,zCMax

  INTEGER :: i,j,k,iC
  INTEGER :: ii,jj,kk
  INTEGER :: iMaxC
  INTEGER :: NumberCells
  INTEGER :: NumberBlocks
  INTEGER :: xCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                   LBOUND(Blocks,2):UBOUND(Blocks,2), &
                   LBOUND(Blocks,3):UBOUND(Blocks,3))  
  INTEGER :: yCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                   LBOUND(Blocks,2):UBOUND(Blocks,2), &
                   LBOUND(Blocks,3):UBOUND(Blocks,3))  
  INTEGER :: zCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                   LBOUND(Blocks,2):UBOUND(Blocks,2), &
                   LBOUND(Blocks,3):UBOUND(Blocks,3))  

  DO iC=1,xCMax
    xCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%xC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i-1,j,k)%xC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j-1,k)%xC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j,k-1)%xC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i+1,j,k)%xC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j+1,k)%xC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j,k+1)%xC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (xCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%xC=xCLoc(i,j,k)-1
          END IF  
        END DO
      END DO
    END DO
  END DO
  DO iC=1,yCMax
    yCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%yC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i-1,j,k)%yC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j-1,k)%yC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j,k-1)%yC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i+1,j,k)%yC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j+1,k)%yC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j,k+1)%yC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (yCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%yC=yCLoc(i,j,k)-1
          END IF  
        END DO
      END DO
    END DO
  END DO
  DO iC=1,zCMax
    zCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%zC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i-1,j,k)%zC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j-1,k)%zC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j,k-1)%zC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i+1,j,k)%zC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j+1,k)%zC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j,k+1)%zC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (zCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%zC=zCLoc(i,j,k)-1
          END IF  
        END DO
      END DO
    END DO
  END DO
  NumberCells=0
  NumberBlocks=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        IF (Blocks(i,j,k)%xC==1) Blocks(i,j,k)%xC=-xCMax
        IF (Blocks(i,j,k)%yC==1) Blocks(i,j,k)%yC=-yCMax
        IF (Blocks(i,j,k)%zC==1) Blocks(i,j,k)%zC=-zCMax
        NumberBlocks=NumberBlocks+2**MAX(xBMax+Blocks(i,j,k)%xC,0) &
                                 *2**MAX(yBMax+Blocks(i,j,k)%yC,0) &
                                 *2**MAX(zBMax+Blocks(i,j,k)%zC,0) 
!       WRITE(*,*) 'Block',i,j,k                         
!       WRITE(*,*) Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC
!       WRITE(*,*) 2**MAX(xBMax+Blocks(i,j,k)%xC,0) &
!                 *2**MAX(yBMax+Blocks(i,j,k)%yC,0) &
!                 *2**MAX(zBMax+Blocks(i,j,k)%zC,0) 
!       WRITE(OutputUnit,*)'                             !Block',i,j,k
!       WRITE(OutputUnit,*) Blocks(i,j,k)%ix0,Blocks(i,j,k)%ix1,'     ! ixa, ixe'
!       WRITE(OutputUnit,*) Blocks(i,j,k)%iy0,Blocks(i,j,k)%iy1,'     ! iya, iye'
!       WRITE(OutputUnit,*) Blocks(i,j,k)%iz0,Blocks(i,j,k)%iz1,'     ! iza, ize'
        NumberCells=NumberCells+(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/(2**(-Blocks(i,j,k)%xC))* &
                                (Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/(2**(-Blocks(i,j,k)%yC))* &
                                (Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/(2**(-Blocks(i,j,k)%zC))
        iMaxC=MIN(Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC)
!       WRITE(OutputUnit,*) iMaxC,Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC,'4 T T'
      END DO
    END DO
  END DO
  WRITE(*,*) 'NumberCells ', NumberCells
  WRITE(*,*) 'NumberBlocks', NumberBlocks
  OPEN(UNIT=OutputUnit,FILE='BlockNeu.grid',STATUS='REPLACE')
  WRITE(OutputUnit,'(a11)') '#Multiblock'
  WRITE(OutputUnit,*) 0,NX,'     ! X ... nxa, nxe'
  WRITE(OutputUnit,*) 0,NY,'     ! Y ... nya, nye'
  WRITE(OutputUnit,*) 0,NZ,'     ! Z ... nza, nze'
  WRITE(OutputUnit,*) ''
  WRITE(OutputUnit,*) NumberBlocks, '       !Number of blocks ... nb'
  NumberBlocks=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        DO ii=1,2**MAX(xBMax+Blocks(i,j,k)%xC,0)
          DO jj=1,2**MAX(yBMax+Blocks(i,j,k)%yC,0)
            DO kk=1,2**MAX(xBMax+Blocks(i,j,k)%zC,0)
              NumberBlocks=NumberBlocks+1
              WRITE(OutputUnit,*)'                             !Block',NumberBlocks
              WRITE(OutputUnit,*) Blocks(i,j,k)%ix0+(ii-1)*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0) &
                                 ,Blocks(i,j,k)%ix0+(ii  )*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0)  
              WRITE(OutputUnit,*) Blocks(i,j,k)%iy0+(jj-1)*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0) &
                                 ,Blocks(i,j,k)%iy0+(jj  )*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0)  
              WRITE(OutputUnit,*) Blocks(i,j,k)%iz0+(kk-1)*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0) &
                                 ,Blocks(i,j,k)%iz0+(kk  )*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0)  
              iMaxC=MIN(Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC)
              WRITE(OutputUnit,*) iMaxC,Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC,'4 T T'
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  CLOSE(OutputUnit)
END SUBROUTINE Coarse
END PROGRAM Grid_Built

