SUBROUTINE SelfMultiblockCompute

  INTEGER :: NX,NY,NZ,BX,BY,BZ,xi,yi,zi
  INTEGER :: xRes,yRes,zRes
  INTEGER :: xIncr,yIncr,zIncr
  INTEGER :: ix,ix0,ix1
  INTEGER :: iy,iy0,iy1
  INTEGER :: iz,iz0,iz1
  CHARACTER :: lc

  TYPE Block_T
    INTEGER :: ix0,ix1,xC=1
    INTEGER :: iy0,iy1,yC=1
    INTEGER :: iz0,iz1,zC=1
    INTEGER :: Coarse=1
    LOGICAL :: Active=.TRUE.
  END TYPE Block_T

  TYPE(Block_T), ALLOCATABLE :: Blocks(:,:,:)
  INTEGER :: xCMax,yCMax,zCMax
  INTEGER :: xBMax,yBMax,zBMax

  INTEGER :: i,j,k
  INTEGER :: NumFineBlocks
  REAL(8) :: xC,yC,zC
  INTEGER :: ib,iF,iC
  INTEGER :: ii,jj,kk
  INTEGER :: iMaxC
  INTEGER :: NumberCells
  INTEGER :: iCoarse,jCoarse,kCoarse
  INTEGER :: iCoarseAct,jCoarseAct,kCoarseAct
  INTEGER, ALLOCATABLE :: xCLoc(:,:,:)
  INTEGER, ALLOCATABLE :: yCLoc(:,:,:)
  INTEGER, ALLOCATABLE :: zCLoc(:,:,:)

  READ(InputUnit,*) Domain%ix0,Domain%ix1,BX
  READ(InputUnit,*) Domain%iy0,Domain%iy1,BY
  READ(InputUnit,*) Domain%iz0,Domain%iz1,BZ
  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  Domain%igx0=Domain%ix0
  Domain%igx1=Domain%ix1
  Domain%igy0=Domain%iy0
  Domain%igy1=Domain%iy1
  Domain%igz0=Domain%iz0
  Domain%igz1=Domain%iz1
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

  xRes=MOD(NX,BX)
  xIncr=NX/BX
  yRes=MOD(NY,BY)
  yIncr=NY/BY
  zRes=MOD(NZ,BZ)
  zIncr=NZ/BZ

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


  ALLOCATE(xCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  
  ALLOCATE(yCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  
  ALLOCATE(zCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  

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
            Blocks(i,j,k)%Coarse=MAX(-xCLoc(i,j,k),Blocks(i,j,k)%Coarse)
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
            Blocks(i,j,k)%Coarse=MAX(-yCLoc(i,j,k),Blocks(i,j,k)%Coarse)
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
  nb=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        IF (Blocks(i,j,k)%xC==1) Blocks(i,j,k)%xC=-xCMax
        IF (Blocks(i,j,k)%yC==1) Blocks(i,j,k)%yC=-yCMax
        IF (Blocks(i,j,k)%zC==1) Blocks(i,j,k)%zC=-zCMax
        Blocks(i,j,k)%Coarse=MAX(-Blocks(i,j,k)%xC,-Blocks(i,j,k)%yC,-Blocks(i,j,k)%zC)
        nb=nb+2**MAX(xBMax+Blocks(i,j,k)%xC,0) &
             *2**MAX(yBMax+Blocks(i,j,k)%yC,0) &
             *2**MAX(zBMax+Blocks(i,j,k)%zC,0) 
!       NumberCells=NumberCells+(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/(2**(-Blocks(i,j,k)%xC))* &
!                               (Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/(2**(-Blocks(i,j,k)%yC))* &
!                               (Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/(2**(-Blocks(i,j,k)%zC))
      END DO
    END DO
  END DO
! Coarsening 
  
  iCoarse=2 
  jCoarse=1 
  kCoarse=2 
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1),iCoarse
    IF (i+iCoarse-1>UBOUND(Blocks,1)) THEN
      iCoarseAct=UBOUND(Blocks,1)-i+1
    ELSE
      iCoarseAct=iCoarse
    END IF  
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2),jCoarse
      IF (j+jCoarse-1>UBOUND(Blocks,2)) THEN
        jCoarseAct=UBOUND(Blocks,2)-j+1
      ELSE
        jCoarseAct=jCoarse
      END IF  
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3),kCoarse
        IF (k+kCoarse-1>UBOUND(Blocks,3)) THEN
          kCoarseAct=UBOUND(Blocks,3)-k+1
        ELSE
          kCoarseAct=kCoarse
        END IF  
        IF (MAXVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)>0.AND. &
          MAXVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)==          &
          MINVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)) THEN
          DO ii=1,iCoarseAct
            DO jj=1,jCoarseAct
              DO kk=1,kCoarseAct
                Blocks(i+ii-1,j+jj-1,k+kk-1)%Active=.FALSE.
              END DO  
            END DO  
          END DO  
          Blocks(i,j,k)%Active=.TRUE.
          nb=nb+1-iCoarseAct*jCoarseAct*kCoarseAct
          Blocks(i,j,k)%ix1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%ix1
          Blocks(i,j,k)%iy1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%iy1
          Blocks(i,j,k)%iz1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%iz1
        END IF
      END DO
    END DO
  END DO

  ALLOCATE(Floor(nb))
  ib=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        IF (Blocks(i,j,k)%Active) THEN
          DO ii=1,2**MAX(xBMax+Blocks(i,j,k)%xC,0)
            DO jj=1,2**MAX(yBMax+Blocks(i,j,k)%yC,0)
              DO kk=1,2**MAX(xBMax+Blocks(i,j,k)%zC,0)
                ib=ib+1
                Floor(ib)%igx0=Blocks(i,j,k)%ix0+(ii-1)*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0) 
                Floor(ib)%igx1=Blocks(i,j,k)%ix0+(ii  )*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0)  
                Floor(ib)%igy0=Blocks(i,j,k)%iy0+(jj-1)*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0)
                Floor(ib)%igy1=Blocks(i,j,k)%iy0+(jj  )*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0)  
                Floor(ib)%igz0=Blocks(i,j,k)%iz0+(kk-1)*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0)
                Floor(ib)%igz1=Blocks(i,j,k)%iz0+(kk  )*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0)  
                iMaxC=MIN(Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC)
                Floor(ib)%Refine=iMaxC
                Floor(ib)%RefineX=Blocks(i,j,k)%xC
                Floor(ib)%RefineY=Blocks(i,j,k)%yC
                Floor(ib)%RefineZ=Blocks(i,j,k)%zC
                Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
                Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
                Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
                Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
                Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
                Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
                Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
                Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineY)) 
                Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineZ)) 
                Floor(ib)%ib=ib
                Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
                Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
                Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
                Floor(ib)%RefLevel=2

                IF (Floor(ib)%Refine<0) THEN
                  Floor(ib)%JacAdvection=.FALSE.
                  Floor(ib)%JacDiffusion=.FALSE.
                ELSE
                  Floor(ib)%JacAdvection=.TRUE.
                  Floor(ib)%JacDiffusion=.TRUE.
                END IF

                Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz

                Floor(ib)%TypeW='iw'
                IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
                  Floor(ib)%TypeW='ow'
                ELSE IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West=='Period') THEN
                   Floor(ib)%TypeW='pw'
                END IF
                Floor(ib)%TypeE='ie'
                IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
                  Floor(ib)%TypeE='oe'
                ELSE IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East=='Period') THEN
                 Floor(ib)%TypeE='pe'
                END IF
                Floor(ib)%TypeS='is'
                IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
                  Floor(ib)%TypeS='os'
                ELSE IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South=='Period') THEN
                   Floor(ib)%TypeS='ps'
                END IF
                Floor(ib)%TypeN='in'
                IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
                  Floor(ib)%TypeN='on'
                ELSE IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North=='Period') THEN
                 Floor(ib)%TypeN='pn'
                END IF
                Floor(ib)%TypeB='ib'
                IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
                  Floor(ib)%TypeB='ob'
                ELSE IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom=='Period') THEN
                  Floor(ib)%TypeB='pb'
                END IF
                Floor(ib)%TypeT='it'
                IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
                  Floor(ib)%TypeT='ot'
                ELSE IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top=='Period') THEN
                  Floor(ib)%TypeT='pt'
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE SelfMultiblockCompute
END MODULE Grid_Mod

