MODULE DynCore_Mod
  USE Kind_Mod
  USE Parameter_Mod

  IMPLICIT NONE

!-----------------------------------------------------------------------
! steady-state and baroclinic wave tuning parameter
!----------------------------------------------------------------------- 

  REAL(RealKind), PARAMETER :: radius                 = 10._RealKind  ! reciprocal radius of the perturbation without 'a'
  REAL(RealKind), PARAMETER :: perturbation_amplitude =  1._RealKind  ! amplitude of u perturbation 1 m/s
  REAL(RealKind), PARAMETER :: perturbation_longitude = 20._RealKind  ! longitudinal position, 20E
  REAL(RealKind), PARAMETER :: perturbation_latitude  = 40._RealKind  ! latitudinal position, 40N
  REAL(RealKind), PARAMETER :: eta_sfc                = 1._RealKind  ! hybrid value at surface
  REAL(RealKind), PARAMETER :: delta_T                = 480000._RealKind  ! in K, for T mean calculation
    !
  REAL(RealKind), PARAMETER :: perturbation_latitude_tracer = 55.d0


CONTAINS

FUNCTION F_lon_lat_eta(lon,lat,eta,height,u0)
  REAL(RealKind) :: F_lon_lat_eta
  REAL(RealKind) :: lon,lat,eta,u0,height
  F_lon_lat_eta=-Grav*height+surface_geopotential(lon,lat,eta,u0)
END FUNCTION F_lon_lat_eta

FUNCTION dFdeta(lon,lat,eta,u0)
  REAL(RealKind) :: dFdeta
  REAL(RealKind) :: lon,lat,eta,u0
  dFdeta=-Rd/eta*temperature(lon,lat,eta,u0)
END FUNCTION dFdeta

FUNCTION temperature(lon,lat,eta,u0) !,rotation_angle)
  REAL(RealKind) :: temperature
  REAL(RealKind) :: eta,u0
  REAL(RealKind) :: lon, lat !, rotation_angle
  REAL(RealKind) :: rot_lon, rot_lat

!  IF (ABS(rotation_angle)<1.0E-8) THEN
     rot_lon = lon
     rot_lat = lat
!  ELSE
!    CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
!  ENDIF

  temperature  = t_mean(eta) + t_deviation(rot_lon,rot_lat,eta,u0)
END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
FUNCTION t_mean(eta)
  REAL(RealKind) :: t_mean
  REAL(RealKind) :: eta

  IF ((eta).gt.(eta_tropo)) THEN
     t_mean = Temp0*eta**exponent1                                ! mean temperature at each level (troposphere)
  ELSE
     t_mean = Temp0*eta**exponent1 + delta_T*(eta_tropo-eta)**5  ! mean temperature at each level (stratosphere)
  ENDIF
END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean 
  ! (equation (6) minus horizontally averaged temperature)
  !
FUNCTION t_deviation(lon,lat,eta,u0)
  REAL(RealKind) :: t_deviation
  REAL(RealKind) :: eta,u0
  REAL(RealKind) :: lon, lat
  REAL(RealKind) :: factor, eta_vertical, rot_lon, rot_lat

  factor       = eta*pi*u0/Rd             
  eta_vertical = (eta - eta0) * 0.5d0*pi

  rot_lon = lon
  rot_lat = lat

  t_deviation = factor * 1.5d0 * SIN(eta_vertical) * (cos(eta_vertical))**0.5d0 *                        &
                ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                u0 * (COS(eta_vertical))**1.5d0  +                                                       &
                (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
END FUNCTION t_deviation

!**************************************************************************
!
! Surface geopotential (equation (7) in Jablonowski and Williamson, 2006)
!
!**************************************************************************  
FUNCTION geo_deviation(lon,lat,eta,u0) !,rotation_angle)
  REAL(RealKind) :: geo_deviation
  REAL(RealKind) :: eta,u0
  REAL(RealKind) :: lon, lat!, rotation_angle
  REAL(RealKind) :: cos_tmp, rot_lon, rot_lat

!  IF (ABS(rotation_angle)<1.0E-8) THEN
     rot_lon = lon
     rot_lat = lat
!  ELSE
!     CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
!  ENDIF
  cos_tmp    = u0 * (cos((eta-eta0)*pi*0.5d0))**1.5d0

  geo_deviation = ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp
END FUNCTION geo_deviation

FUNCTION geo_mean(eta)
  REAL(RealKind) :: geo_mean
  REAL(RealKind) :: factor,eta

  factor=Temp0*Grav/gamma1
  IF (eta.GE.eta_tropo) THEN
    geo_mean = factor*(1.d0-eta**exponent1)
  ELSE
    geo_mean = factor*(1.d0-eta**exponent1)    &
              - Rd*delta_T*( (LOG(eta/eta_tropo)+137.d0/60.d0)*eta_tropo**5.d0   &
                              - 5.d0*eta_tropo**4.d0*eta                &
                              + 5.d0*eta_tropo**3.d0*eta**2.d0          &
                              - 10.d0/3.d0*eta_tropo**2.d0*eta**3.d0    &
                              + 5.d0/4.d0*eta_tropo*eta**4.d0           &
                              - 1.d0/5.d0*eta**5.d0      )
  END IF
END FUNCTION geo_mean

FUNCTION surface_geopotential(lon,lat,eta,u0)
  REAL(RealKind) :: surface_geopotential
  REAL(RealKind) :: lon,lat,eta,u0
  surface_geopotential=geo_mean(eta) + geo_deviation(lon,lat,eta,u0)

END FUNCTION surface_geopotential


END MODULE DynCore_Mod

