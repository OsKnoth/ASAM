! Calculation of the collisions between ice particles and ice particles
! only if a considerable water shell is present for both collision partners 
! using the two-dimensional Linear Discrete Method LDM
!
! calling koll_s.f

SUBROUTINE koll_eis_eisxd(knf,kqf,kqwf,kqsf,kqaf,nf,qf,qwf,qsf,qaf,      &
                          rho,TT,miquer,siquer,ip)

INCLUDE 'HEADER3'

IMPLICIT NONE

INTEGER :: ip
DOUBLE PRECISION :: rho,TT
DOUBLE PRECISION :: nf(jmax,smax,ipmax),qf(jmax,smax,ipmax)
DOUBLE PRECISION :: qsf(jmax,smax,ipmax),qaf(jmax,smax,ipmax)
DOUBLE PRECISION :: qwf(jmax,smax,ipmax),kqwf(jmax,smax)
DOUBLE PRECISION :: knf(jmax,smax),kqf(jmax,smax),kqsf(jmax,smax)
DOUBLE PRECISION :: kqaf(jmax,smax),miquer(jmax,smax),siquer(jmax,smax)

INTEGER :: i,i1,i2,j,j1,j2,s,k,ikolk      
DOUBLE PRECISION :: ddqf12,ddqf21,ddnf12,ddnf21,ddqsf12,ddqsf21
DOUBLE PRECISION :: ddqaf12,ddqaf21,nhilf1,nhilf2,mneu,sneu,aneu
DOUBLE PRECISION :: xsum1,xsum2,nneu,xgre,xmit,xdif,f1n,f2n,f1m,f2m
DOUBLE PRECISION :: f1s,f2s,small,kolk_eff
DOUBLE PRECISION :: ddqwf12,ddqwf21,mwfneu


! small constant  
! VG 
!small = 1.0D-5
small = 1.0D0

! checking the collisions
DO i1=1,smax
  DO j1=1,jmax
    IF(nf(j1,i1,ip).LE.small) GOTO 200
    DO i2=1,i1
      DO j2=1,j1
        IF(nf(j2,i2,ip).LE.small) GOTO 100
! For a sufficient water shell, collision occurs
        IF(DABS(qf(j1,i1,ip)/qwf(j1,i1,ip)).LE.0.7D0.OR.                &
           DABS(qf(j2,i2,ip)/qwf(j2,i2,ip)).LE.0.7D0) THEN
          ikolk = 1
          CALL kolk_avexd(kolk_eff,j1,i1,j2,i2,miquer(j1,i1),           &
                          siquer(j1,i1),miquer(j2,i2),siquer(j2,i2),    &
                          ikolk,1)
          IF(kolk_eff.EQ.0.0D0) GOTO 100

! target bins
          k = KZIEL(j1,j2)
          s = SZIEL(i1,i2)
          IF(smax.EQ.1) s = 1
          IF(s.GT.smax) WRITE(*,*) 'koll_contactxd S = ', s, smax

! interaction terms
          nhilf1 = nf(j1,i1,ip)*kolk_eff
          nhilf2 = nf(j2,i2,ip)*kolk_eff

! Unit volume instead of unit mass => change of air density, kg^-1 => m^-3
! NW, QW, QS in m^-3 and DD?? also in m^-3, therefore, only one RHOTOT
          nhilf1 = nhilf1*(rho)
          nhilf2 = nhilf2*(rho)

! correction due to increasing fall velocity with decreasing density
          nhilf1 = nhilf1*SQRT(RHOTOT0/rho)
          nhilf2 = nhilf2*SQRT(RHOTOT0/rho)

          ddnf12 = nf(j1,i1,ip)*nhilf2
          ddnf21 = ddnf12
          ddqf12 = qf(j1,i1,ip)*nhilf2
          ddqf21 = qf(j2,i2,ip)*nhilf1
          ddqwf12 = qwf(j1,i1,ip)*nhilf2
          ddqwf21 = qwf(j2,i2,ip)*nhilf1
          ddqsf12 = qsf(j1,i1,ip)*nhilf2
          ddqsf21 = qsf(j2,i2,ip)*nhilf1
          ddqaf12 = qaf(j1,i1,ip)*nhilf2
          ddqaf21 = qaf(j2,i2,ip)*nhilf1



! collision with drops out of the same size and aerosol bin J1=J2 and I1=I2
          IF((j1.EQ.j2.OR.k.EQ.jmax).AND.(i1.EQ.i2.OR.s.EQ.smax)) THEN

! source for (k,s)
            knf(k,s) = knf(k,s) + ddnf12/2.D0
            kqf(k,s) = kqf(k,s) + ddqf12
            kqwf(k,s) = kqwf(k,s) + ddqwf12
            kqsf(k,s) = kqsf(k,s) + ddqsf12
            kqaf(k,s) = kqaf(k,s) + ddqaf12
  
! sink for (j1,i1) = (j2,i2)
            knf(j1,i1) = knf(j1,i1) - ddnf12
            kqf(j1,i1) = kqf(j1,i1) - ddqf12
            kqwf(j1,i1) = kqwf(j1,i1) - ddqwf12
            kqsf(j1,i1) = kqsf(j1,i1) - ddqsf12
            kqaf(j1,i1) = kqaf(j1,i1) - ddqaf12

! two drops form one new drop.
! Therefore, qf, qwf, qsf, qaf are conserved, nf is divided by 2
          END IF

! collision of two drops with the same size j1 = j2, but different aerosol
! content i1 > i2
          IF((j1.EQ.j2.OR.k.EQ.jmax).AND.(i1.GT.i2.AND.s.LT.smax)) THEN

            mneu = ddqf12 + ddqf21
            mwfneu = ddqwf12 + ddqwf21
! parameters for linear approximation (aerosol)
            xsum1 = sgrenz(i1) + sgrenz(i2)
            xsum2 = sgrenz(i1+1) + sgrenz(i2+1)
            sneu = ddqsf12 + ddqsf21
            aneu = ddqaf12 + ddqaf21
            nneu = ddnf12
            xdif = sdiff21(i1) + sdiff21(i2)
            xmit = smitte(i1) + smitte(i2)
            xgre = sgrenz(s+1)

            CALL koll_sxd(xsum1,xsum2,xgre,aneu,nneu,xdif,xmit,f2s,f2n,10)

            f2m = f2n

! source for (k,s) and (k,s+1)
            knf(k,s) = knf(k,s) + f2n*nneu
            kqf(k,s) = kqf(k,s) + f2m*mneu
            kqwf(k,s) = kqwf(k,s) + f2m*mwfneu
            kqsf(k,s) = kqsf(k,s) + f2s*sneu
            kqaf(k,s) = kqaf(k,s) + f2s*aneu
           
            knf(k,s+1) = knf(k,s+1) + (1.0D0 - f2n)*nneu
            kqf(k,s+1) = kqf(k,s+1) + (1.0D0 - f2m)*mneu
            kqwf(k,s+1) = kqwf(k,s+1) + (1.0D0 - f2m)*mwfneu
            kqsf(k,s+1) = kqsf(k,s+1) + (1.0D0 - f2s)*sneu
            kqaf(k,s+1) = kqaf(k,s+1) + (1.0D0 - f2s)*aneu

! sinks for (j1,i1) and (j2,i2)
            knf(j1,i1) = knf(j1,i1) - ddnf12
            kqf(j1,i1) = kqf(j1,i1) - ddqf12
            kqwf(j1,i1) = kqwf(j1,i1) - ddqwf12
            kqsf(j1,i1) = kqsf(j1,i1) - ddqsf12
            kqaf(j1,i1) = kqaf(j1,i1) - ddqaf12
           
            knf(j2,i2) = knf(j2,i2) - ddnf21
            kqf(j2,i2) = kqf(j2,i2) - ddqf21
            kqwf(j2,i2) = kqwf(j2,i2) - ddqwf21
            kqsf(j2,i2) = kqsf(j2,i2) - ddqsf21
            kqaf(j2,i2) = kqaf(j2,i2) - ddqaf21

          END IF

! collision of two drops with same aerosol content i1 = i2, but of
! different size j1 > j2
          IF((j1.GT.j2.AND.k.LT.jmax).AND.(i1.EQ.i2.OR.s.EQ.smax)) THEN

            mneu = ddqf12 + ddqf21
            mwfneu = ddqwf12 + ddqwf21
! parameters for linear approximation (size)
            xsum1 = mgrenz(j1) + mgrenz(j2)
            xsum2 = mgrenz(j1+1) + mgrenz(j2+1)
            sneu = ddqsf12 + ddqsf21
            aneu = ddqaf12 + ddqaf21
            nneu = ddnf12
            xdif = diff21(j1) + diff21(j2)
            xmit = mmitte(j1) + mmitte(j2)
            xgre = mgrenz(k+1)

! CMS Muss es in der ‹bergabe mwfneu anstatt mneu heiﬂen???
!            CALL koll_sxd(xsum1,xsum2,xgre,mneu,nneu,xdif,xmit,f1m,f1n,11)
            CALL koll_sxd(xsum1,xsum2,xgre,mwfneu,nneu,xdif,xmit,f1m,f1n,11)
 
            f1s = f1n

! source for (k,s) and (k+1,s)
            knf(k,s) = knf(k,s) + f1n*nneu
            kqf(k,s) = kqf(k,s) + f1m*mneu
            kqwf(k,s) = kqwf(k,s) + f1m*mwfneu
            kqsf(k,s) = kqsf(k,s) + f1s*sneu
            kqaf(k,s) = kqaf(k,s) + f1s*aneu

            knf(k+1,s) = knf(k+1,s) + (1.0D0 - f1n)*nneu
            kqf(k+1,s) = kqf(k+1,s) + (1.0D0 - f1m)*mneu
            kqwf(k+1,s) = kqwf(k+1,s) + (1.0D0 - f1m)*mwfneu
            kqsf(k+1,s) = kqsf(k+1,s) + (1.0D0 - f1s)*sneu
            kqaf(k+1,s) = kqaf(k+1,s) + (1.0D0 - f1s)*aneu

! sinks for(j1,i1) and (j2,i2)
            knf(j1,i1) = knf(j1,i1) - ddnf12
            kqf(j1,i1) = kqf(j1,i1) - ddqf12
            kqwf(j1,i1) = kqwf(j1,i1) - ddqwf12
            kqsf(j1,i1) = kqsf(j1,i1) - ddqsf12
            kqaf(j1,i1) = kqaf(j1,i1) - ddqaf12

            knf(j2,i2) = knf(j2,i2) - ddnf21
            kqf(j2,i2) = kqf(j2,i2) - ddqf21
            kqwf(j2,i2) = kqwf(j2,i2) - ddqwf21
            kqsf(j2,i2) = kqsf(j2,i2) - ddqsf21
            kqaf(j2,i2) = kqaf(j2,i2) - ddqaf21

        END IF




! collision of two drops with different size I1>I2 and different aerosol
! content J1>J2
          IF((j1.GT.j2.AND.k.LT.jmax).AND.(i1.GT.i2.AND.s.LT.smax)) THEN

            mneu = ddqf12 + ddqf21
            mwfneu = ddqwf12 + ddqwf21
! parameters for linear approximation (size)
            xsum1 = mgrenz(j1) + mgrenz(j2)
            xsum2 = mgrenz(j1+1) + mgrenz(j2+1)
            nneu = ddnf12
            xdif = diff21(j1) + diff21(j2)
            xmit = mmitte(j1) + mmitte(j2)
            xgre = mgrenz(k+1)

! CMS Muss es in der ‹bergabe mwfneu anstatt mneu heiﬂen???
!            CALL koll_sxd(xsum1,xsum2,xgre,mneu,nneu,xdif,xmit,f1m,f1n,12)
            CALL koll_sxd(xsum1,xsum2,xgre,mwfneu,nneu,xdif,xmit,f1m,f1n,12)

            f1s = f1n

! parameters for linear approximation (aerosol)
            xsum1 = sgrenz(i1) + sgrenz(i2)
            xsum2 = sgrenz(i1+1) + sgrenz(i2+1)
            sneu = ddqsf12 + ddqsf21
            aneu = ddqaf12 + ddqaf21
            nneu = ddnf12
            xdif = sdiff21(i1) + sdiff21(i2)
            xmit = smitte(i1) + smitte(i2)
            xgre = sgrenz(s+1)

            CALL koll_sxd(xsum1,xsum2,xgre,aneu,nneu,xdif,xmit,f2s,f2n,13)

            f2m = f2n

! sources for (k,s), (k+1,s), (k,s+1) and (k+1,s+1)
            knf(k,s) = knf(k,s) + f1n*f2n*nneu
            kqf(k,s) = kqf(k,s) + f1m*f2m*mneu
            kqwf(k,s) = kqwf(k,s) + f1m*f2m*mwfneu
            kqsf(k,s) = kqsf(k,s) + f1s*f2s*sneu
            kqaf(k,s) = kqaf(k,s) + f1s*f2s*aneu

            knf(k+1,s) = knf(k+1,s) +(1.0D0 - f1n)*f2n*nneu
            kqf(k+1,s) = kqf(k+1,s) +(1.0D0 - f1m)*f2m*mneu
            kqwf(k+1,s) = kqwf(k+1,s) + (1.0D0 - f1m)*f2m*mwfneu
            kqsf(k+1,s) = kqsf(k+1,s) + (1.0D0 - f1s)*f2s*sneu
            kqaf(k+1,s) = kqaf(k+1,s) + (1.0D0 - f1s)*f2s*aneu

            knf(k,s+1) = knf(k,s+1) + f1n*(1.0D0 - f2n)*nneu
            kqf(k,s+1) = kqf(k,s+1) + f1m*(1.0D0 - f2m)*mneu
            kqwf(k,s+1) = kqwf(k,s+1) + f1m*(1.0D0 - f2m)*mwfneu
            kqsf(k,s+1) = kqsf(k,s+1) + f1s*(1.0D0 - f2s)*sneu
            kqaf(k,s+1) = kqaf(k,s+1) + f1s*(1.0D0 - f2s)*aneu

            knf(k+1,s+1) = knf(k+1,s+1) + (1.D0 - f1n)*(1.D0 - f2n)*nneu
            kqf(k+1,s+1) = kqf(k+1,s+1) + (1.D0 - f1m)*(1.D0 - f2m)*mneu
            kqwf(k+1,s+1) = kqwf(k+1,s+1) + (1.D0 - f1m)*(1.D0 - f2m)*mwfneu
            kqsf(k+1,s+1) = kqsf(k+1,s+1) + (1.D0 - f1s)*(1.D0 - f2s)*sneu
            kqaf(k+1,s+1) = kqaf(k+1,s+1) + (1.D0 - f1s)*(1.D0 - f2s)*aneu

! sinks for (j1,i1) and (j2,i2)
            knf(j1,i1) = knf(j1,i1) - ddnf12
            kqf(j1,i1) = kqf(j1,i1) - ddqf12
            kqwf(j1,i1) = kqwf(j1,i1) - ddqwf12
            kqsf(j1,i1) = kqsf(j1,i1) - ddqsf12
            kqaf(j1,i1) = kqaf(j1,i1) - ddqaf12
 
            knf(j2,i2) = knf(j2,i2) - ddnf21
            kqf(j2,i2) = kqf(j2,i2) - ddqf21
            kqwf(j2,i2) = kqwf(j2,i2) - ddqwf21
            kqsf(j2,i2) = kqsf(j2,i2) - ddqsf21
            kqaf(j2,i2) = kqaf(j2,i2) - ddqaf21

        END IF



        END IF
100     CONTINUE
      END DO
    END DO

200 CONTINUE
  END DO
END DO


RETURN
END 

