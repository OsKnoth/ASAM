MODULE TurbulentInflow_Mod

  USE Kind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: Amp
  REAL(RealKind) :: uMax
  REAL(RealKind) :: vMax
  REAL(RealKind) :: Inflow_lenx=1.0d0  
  REAL(RealKind) :: Inflow_leny=1.0d0  
  REAL(RealKind) :: Inflow_lenz=1.0d0  
  REAL(RealKind) :: offset_x=1.0d0  
  REAL(RealKind) :: offset_x1=1.0d0  
  REAL(RealKind) :: offset_y=1.0d0  
  REAL(RealKind) :: offset_y1=1.0d0  
  REAL(RealKind) :: offset_z=1.0d0  
  REAL(RealKind) :: offset_z1=1.0d0  
  REAL(RealKind) :: intenz=0.1d0


CONTAINS

FUNCTION vTurb(x,y,z)
  REAL(RealKind) :: vTurb
  REAL(RealKind) :: x,y,z,Time

  INTEGER :: s,ss
  REAL(RealKind) :: k,kNull,k_End
  REAL(RealKind) :: wx,wy,wz
  REAL(RealKind) :: VNorm
  REAL(RealKind) :: Pi

  REAL(RealKind) :: uTurb

  Pi=ATAN(1.0d0)*4.0d0

  knull=18.0d0
  k_end=INT(Inflow_leny/4.0)
  DO s=1,10 !k_end
    k=((10.0)*(0.9d0**s))
    DO ss=1,2
      VNorm=SQRT(uMax*uMax+vMax*vMax)
      wx=(-One)**ss*k*VNorm/(offset_x1-offset_x)
      wy=(-One)**ss*k*VNorm/(offset_y1-offset_y)
      wz=(-One)**ss*k*VNorm/(offset_z1-offset_z)

      vTurb=AMP*COS(k*((x-offset_x)/(offset_x1-offset_x))*Two*Pi+ &
            Time*(wx)*Two*Pi+k**(2.0))  &
           +AMP*COS(k*((y-offset_y)/(offset_y1-offset_y))*Two*Pi+ &
            Time*(wy)*Two*Pi+k**(2.0))  &
           +AMP*COS(k*((z-offset_z)/(offset_z1-offset_z))*Two*Pi+ &
            Time*(wz)*Two*Pi+k**(2.0))
      uTurb=AMP*SIN(k*((x-offset_x)/(offset_x1-offset_x))*Two*Pi + &
            Time*(wx)*Two*Pi+k**(2.0))  &
           +AMP*SIN(k*((y-offset_y)/(offset_y1-offset_y))*Two*Pi + &
            Time*(wy)*Two*Pi+k**(2.0))  &
           +AMP*SIN(k*((z-offset_z)/(offset_z1-offset_z))*Two*Pi + &
            Time*(wz)*Two*Pi+k**(2.0))
    END DO
  END DO
END FUNCTION vTurb
END MODULE TurbulentInflow_Mod
