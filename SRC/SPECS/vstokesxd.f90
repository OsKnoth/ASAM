! Fallgeschwindigkeiten nach Pruppacher und Klett (1997) S.415ff

SUBROUTINE vstokesxd(vtw,vtf,kz,dz,pp,rho,TT,T0_local)

INCLUDE 'HEADER2'

USE data_parallel,      ONLY: my_cart_id

IMPLICIT NONE

INTEGER :: kz
DOUBLE PRECISION :: vtw(kz,jmax),vtf(kz,jmax),pp(kz),rho(kz)
DOUBLE PRECISION :: TT(kz),T0_local,dz

INTEGER :: j,k
DOUBLE PRECISION :: cdnre2,xi,yi,NBo,NP,NRe,z,lambda

DOUBLE PRECISION :: bi0,bi1,bi2,bi3,bi4,bi5,bi6
PARAMETER(bi0 = -0.318657D1, bi1 = 0.992696, bi2 = -0.153193D-2)
PARAMETER(bi3 = -0.987059D-3, bi4 = -0.578878D-3)
PARAMETER(bi5 = 0.855176D-4, bi6 = -0.327815D-5)
DOUBLE PRECISION :: bw0,bw1,bw2,bw3,bw4,bw5
PARAMETER(bw0 = -0.500015D1, bw1 = 0.523778D1, bw2 = -0.204914D1)
PARAMETER(bw3 = 0.475294, bw4 = -0.542819D-1, bw5 = 0.238449D-2)

! Function
DOUBLE PRECISION :: etaxd,sigmaxd

IF (my_cart_id == 0) THEN 
  OPEN(11,FILE='RESULTATE/vtw')
  OPEN(12,FILE='RESULTATE/vtf')
ENDIF

! VG reverted arrays for LM
! pp(1) --> pp(40)

z = -dz
DO k=1,kz
  z = z + dz      
  DO j=1,jmax

! Fallgeschwindigkeiten fuer sphaerische Eisteilchen (Hagel, Graupel)
    IF(rmitte(j).LE.10.0D-6) THEN
      lambda = 6.6D-8*(pp(40)/pp(k))*(TT(k)/T0_local)
      vtf(k,j) = -(1.0D0 + 1.26D0*lambda/rmitte(j))*2.0D0             &
                 *rmitte(j)**2.0D0*grav*(rhoi - rho(k))/(9.0D0*etaxd(TT(k)))
    ELSE IF(rmitte(j).GT.10.0D-6.AND.rmitte(j).LE.535.0D-6) THEN
      cdnre2 = 32.0D0*rmitte(j)**3.0D0*(rhoi - rho(k))*rho(k)         &
               *grav/(3.0D0*etaxd(TT(k))**2.0D0)
      xi = DLOG(cdnre2)
      yi = bi0 + bi1*xi + bi2*xi**2.0D0 + bi3*xi**3.0D0               &
         + bi4*xi**4.0D0 + bi5*xi**5.0D0 + bi6*xi**6.0D0
      vtf(k,j) = -etaxd(TT(k))*DEXP(yi)/(2.0D0*rho(k)*rmitte(j))
    ELSE IF(rmitte(j).GT.535.0D-6) THEN
! Fuer Hagel und Graupel nach Pruppacher und Klett (1997) S.441
! VG diese Formel "knickt" die Geschwindigkeitsverteilung
!    vor allem ist sie Hoehenunabhaengig. Neue Berechnung analog zu Tropfen
!      vtf(k,j) = -9.0D0*(2.0D0*rmitte(j)*100.0D0)**0.8D0
      NBo = grav*(rhoi - rho(k))*rmitte(j)**2.0D0/sigmaxd(TT(k))
      NP = sigmaxd(TT(k))**3.0D0*rho(k)**2.0D0/(etaxd(TT(k))**4.0D0       &
           *grav*(rhoi - rho(k)))
      xi = DLOG(16.0D0*NBo*NP**(1.0D0/6.0D0)/3.0D0)
      yi = bw0 + bw1*xi + bw2*xi**2.0D0 + bw3*xi**3.0D0               &
         + bw4*xi**4.0D0 + bw5*xi**5.0D0
      NRe = NP**(1.0D0/6.0D0)*DEXP(yi)
      vtf(k,j) = -etaxd(TT(k))*NRe/(2.0D0*rho(k)*rmitte(j)) 
   END IF

! Fallgeschwindigkeiten fuer Wassertropfen
    IF(rmitte(j).LE.10.0D-6) THEN
      lambda = 6.6D-8*(pp(40)/pp(k))*(TT(k)/T0_local)
      vtw(k,j) = -(1.0D0 + 1.26D0*lambda/rmitte(j))*2.0D0             &
                 *rmitte(j)**2.0D0*grav*(rhow - rho(k))/(9.0D0*etaxd(TT(k)))
    ELSE IF(rmitte(j).GT.10.0D-6.AND.rmitte(j).LE.535.0D-6) THEN
      cdnre2 = 32.0D0*rmitte(j)**3.0D0*(rhow - rho(k))*rho(k)         &
               *grav/(3.0D0*etaxd(TT(k))**2.0D0)
      xi = DLOG(cdnre2)
      yi = bi0 + bi1*xi + bi2*xi**2.0D0 + bi3*xi**3.0D0               &
         + bi4*xi**4.0D0 + bi5*xi**5.0D0 + bi6*xi**6.0D0
      vtw(k,j) = -etaxd(TT(k))*DEXP(yi)/(2.0D0*rho(k)*rmitte(j))
    ELSE IF(rmitte(j).GT.535.0D-6) THEN
      NBo = grav*(rhow - rho(k))*rmitte(j)**2.0D0/sigmaxd(TT(k))
      NP = sigmaxd(TT(k))**3.0D0*rho(k)**2.0D0/(etaxd(TT(k))**4.0D0       &
           *grav*(rhow - rho(k)))
      xi = DLOG(16.0D0*NBo*NP**(1.0D0/6.0D0)/3.0D0)
      yi = bw0 + bw1*xi + bw2*xi**2.0D0 + bw3*xi**3.0D0               &
         + bw4*xi**4.0D0 + bw5*xi**5.0D0
      NRe = NP**(1.0D0/6.0D0)*DEXP(yi)
      vtw(k,j) = -etaxd(TT(k))*NRe/(2.0D0*rho(k)*rmitte(j))
    END IF

IF (my_cart_id == 0) THEN 
    WRITE(11,100) rmitte(j),z,vtw(k,j)
    WRITE(12,100) rmitte(j),z,vtf(k,j)
ENDIF

  END DO

IF (my_cart_id == 0) THEN 
  WRITE(11,*)
  WRITE(12,*)
ENDIF

END DO

IF (my_cart_id == 0) THEN 
CLOSE(11)
CLOSE(12)
ENDIF

100  FORMAT(D13.6,4X,D13.6,4X,D13.6)

           
RETURN
END
