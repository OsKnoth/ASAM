SUBROUTINE sammelnxd(nw,qw,qws,qwa,nf,qf,qfs,qfa,qfw,rrnw,rrqw,rrw,rrnf,   &
                     rrqf,rrf,kz,jmax,smax,ipmax,r1,dt,dz)

IMPLICIT NONE
INTEGER :: kz,jmax,smax,ipmax
DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qws(kz,jmax,smax,ipmax),qwa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: nf(kz,jmax,smax,ipmax),qf(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfs(kz,jmax,smax,ipmax),qfa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: rrnw(jmax),rrqw(jmax),rrw
DOUBLE PRECISION :: rrnf(jmax),rrqf(jmax),rrf
DOUBLE PRECISION :: r1,dt,dz

INTEGER :: i,j,k,ip
DOUBLE PRECISION :: vgz,vquer,pi,vgesw,vgesf,rlinks,rrechts,rmitte,rrtw,rrtf


pi = 2.0D0*ASIN(1.0D0)

! Volumen der unteren inneren Gitterzelle berechnen.
vgz = pi*dz*r1**2.0D0

! Regentropfen haben Radien >= 0.1 mm.
! Der Klasse 49 entspricht ein Radius von 0.065 mm.
! vquer ist das mittlere Tropfenvolumen und vges ist das gesamte
! Tropfenvolumen in der unteren Gitterzelle.
vgesw = 0.0D0
vgesf = 0.0D0
DO j=49,jmax
  rlinks = 1.0D-9*2.0D0**((DBLE(FLOAT(j)) - 1.0D0)/3.0D0)
  rrechts = 1.0D-9*2.0D0**(DBLE(FLOAT(j))/3.0D0)
  rmitte = (rlinks + rrechts)/2.0D0
  vquer = 4.0D0*pi*rmitte**3.0D0/3.0D0
  DO i=1,smax
    vgesw = vgesw + vquer*nw(1,j,i,1)*vgz
    vgesf = vgesf + vquer*nf(1,j,i,1)*vgz
    rrnw(j) = rrnw(j) + nw(1,j,i,1)
    rrqw(j) = rrqw(j) + qw(1,j,i,1)
    rrnf(j) = rrnf(j) + nf(1,j,i,1)
    rrqf(j) = rrqf(j) + qf(1,j,i,1)
    DO ip=1,ipmax
      nw(1,j,i,ip) = 0.0D0
      qw(1,j,i,ip) = 0.0D0
      qws(1,j,i,ip) = 0.0D0
      qwa(1,j,i,ip) = 0.0D0
      nf(1,j,i,ip) = 0.0D0
      qf(1,j,i,ip) = 0.0D0 
      qfs(1,j,i,ip) = 0.0D0
      qfa(1,j,i,ip) = 0.0D0
      qfw(1,j,i,ip) = 0.0D0
    END DO
  END DO
END DO
! Niederschlag in [mm] berechnen.
rrtw = 1000.0D0*vgesw/(pi*r1**2.0D0)
rrtf = 1000.0D0*vgesf/(pi*r1**2.0D0)
rrw = rrw + rrtw
rrf = rrf + rrtf
! Niederschlag in [mm/h]
! rr = rr*3600.0D0*1000.0D0/dt


RETURN
END
