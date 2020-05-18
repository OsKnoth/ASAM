SUBROUTINE sconstxd(J12,R0W)

! Definition of the sectional grid, ...
!
! called by wmain.f
!
! calling none
!
USE data_parallel,      ONLY: my_cart_id

INCLUDE 'HEADER2'

! Definition der Variablen 
IMPLICIT NONE


INTEGER :: II,JJ,J1,J2,J12,KK
DOUBLE PRECISION :: msum1,msum2,R0W,M0W,J0W,R0S,M0S,J0S


!       write(*,*) "Hier sconst"
! grid resolution and lower boundary
J0W   = (NMAX-1.D0)/LOG(2.D0)
J0S   = (TMAX-1.D0)/LOG(2.D0)
R0W  = 1.D-9
M0W  = FACT*R0W**3
R0S  = 1.D-9
M0S  = fact_AS*R0S**3
! Calculation of the bin boundaries and bin centers (aerosol)
DO JJ=1,SIMAX+1
  SGRENZ(JJ)=M0S*EXP((JJ-1.D0)/J0S)
  RSGRENZ(JJ)=(SGRENZ(JJ)/fact_AS)**QU1D3
ENDDO
DO JJ=1,SIMAX
  SMITTE(JJ)=(SGRENZ(JJ)+SGRENZ(JJ+1))/2.D0
  SDIFF21(JJ)=SGRENZ(JJ+1)-SGRENZ(JJ)
  RSMITTE(JJ)=(SMITTE(JJ)/fact_AS)**QU1D3
ENDDO
! Calculation of the bin boundaries and bin centers (water)
J12=0
DO JJ=1,JMAX+1
  MGRENZ(JJ)=M0W*EXP((JJ-1.D0)/J0W)
  RGRENZ(JJ)=(MGRENZ(JJ)/FACT)**QU1D3
ENDDO
DO JJ=1,JMAX
  MMITTE(JJ)=(MGRENZ(JJ)+MGRENZ(JJ+1))/2.D0
  DIFF21(JJ)=MGRENZ(JJ+1)-MGRENZ(JJ)
  RMITTE(JJ)=(MMITTE(JJ)/FACT)**QU1D3
  IF(J12.EQ.0.AND.RMITTE(JJ).GT.12.D-6) J12=JJ
ENDDO
! allocation of collision products to target bins (water)
DO J1=1,JMAX
  DO J2=1,JMAX
    IF(J1.EQ.J2) THEN
      KZIEL(J1,J2)=MIN(JMAX,J1+NMAX-1)
    ELSE 
      msum1=MGRENZ(J1)+MGRENZ(J2)
      msum2=MGRENZ(J1+1)+MGRENZ(J2+1)
      DO II=MAX(J1,J2),JMAX
        IF(MGRENZ(II).LE.msum1.AND.MGRENZ(II+1).GT.msum1) THEN
          KZIEL(J1,J2)=II
          GOTO 50
        ENDIF
      ENDDO
      KZIEL(J1,J2)=JMAX
50    CONTINUE
    ENDIF
  ENDDO
ENDDO
! allocation of collision products to target bins (aerosol)
do JJ=1,SIMAX
  do II=1,SIMAX
    IF(JJ.EQ.II) THEN
      SZIEL(JJ,II)=MIN(JJ+TMAX-1,SIMAX)
    ELSE
      msum1=SGRENZ(JJ)+SGRENZ(II)
      msum2=SGRENZ(JJ+1)+SGRENZ(II+1)
      do KK=MAX(JJ,II),SIMAX
        if(SGRENZ(KK).le.msum1.and.SGRENZ(KK+1).gt.msum1) then
          SZIEL(JJ,II)=KK
          goto 51
        endif
      enddo
      SZIEL(JJ,II)=SIMAX
51    continue
    ENDIF
  enddo
enddo
! Estimation of the activated share
! 'special' indices
do JJ=JMAX-1,1,-1
  if(rgrenz(JJ+1).gt.0.5D-6) J_D01=JJ
  if(rgrenz(JJ+1).gt.2.5D-6) J_D05=JJ
  if(rgrenz(JJ+1).gt.5.D-6)  J_D10=JJ
enddo

IF (my_cart_id == 0) THEN 
! print out some constants
OPEN(20,FILE="RESULTATE/kon")
REWIND 20
WRITE(20,*) "PI      =   ", PI    
WRITE(20,*) "R0W     =   ", R0W   
WRITE(20,*) "M0W     =   ", M0W   
WRITE(20,*) "R0S     =   ", R0S   
WRITE(20,*) "M0S     =   ", M0S   
WRITE(20,*) "J12     =   ", J12
WRITE(20,*) "R(J12)  =   ", RMITTE(J12)   
CLOSE(20)

OPEN(20,FILE="RESULTATE/mrel")
REWIND 20
WRITE(20,*) 1,MGRENZ(1)/MMITTE(1)
WRITE(20,*) JMAX,MGRENZ(JMAX)/MMITTE(JMAX)
WRITE(20,*) 
WRITE(20,*) 1,MGRENZ(2)/MMITTE(1)
WRITE(20,*) JMAX,MGRENZ(JMAX+1)/MMITTE(JMAX)
CLOSE(20)

!!       RETURN
!!       OPEN(20,FILE="RESULTATE/kziel")
!!       REWIND 20
!!       DO J1=1,JMAX
!!         DO J2=1,JMAX
!!           WRITE(20,*) J1,J2,KZIEL(J1,J2)
!!         ENDDO
!!       ENDDO
!!       CLOSE(20)
ENDIF
  
RETURN
END 
