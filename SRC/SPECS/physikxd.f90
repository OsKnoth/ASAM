!*****************************************************************
!*                                                               *
!*    subroutine  physikxd                                       *
!*                                                               *
!*                                                               *
!*****************************************************************
SUBROUTINE physikxd(teta,pp,T0_local,qq,usatt,nw,qw,qws,qwa,nf,qf,qfs,qfa,qfw,   &
                    ni,qi,rho,kz,mquer,rquer,squer,mfquer,rfquer,sfquer,   &
                    Talt,anz,dteta,dtetadyn,dqq,dqqdyn,dnw,dqw,dqws,dqwa,  &
                    dnf,dqf,dqfs,dqfa,dqfw,dni,dqi,ifreeze,miv,vtw,vtf)

INCLUDE 'HEADER2'


IMPLICIT NONE

INTEGER :: kz,anz,ifreeze(itmax)
DOUBLE PRECISION :: T0_local,Talt(kz,ipmax),rho(kz)
DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qws(kz,jmax,smax,ipmax),qwa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: nf(kz,jmax,smax,ipmax),qf(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfs(kz,jmax,smax,ipmax),qfa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: ni(kz,simax,itmax,ipmax),qi(kz,simax,itmax,ipmax)
DOUBLE PRECISION :: teta(kz,ipmax),pp(kz),qq(kz,ipmax),usatt(kz,ipmax)
DOUBLE PRECISION :: mquer(kz,jmax,smax,ipmax),rquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: squer(kz,jmax,smax,ipmax),mfquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: rfquer(kz,jmax,smax,ipmax),sfquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: miv(itmax),vtw(kz,jmax),vtf(kz,jmax)
! dx Groessen      
DOUBLE PRECISION :: dteta(kz,ipmax),dtetadyn(kz,ipmax),dqq(kz,ipmax)
DOUBLE PRECISION :: dqqdyn(kz,ipmax)
DOUBLE PRECISION :: dnw(kz,jmax,smax,ipmax),dqw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dqws(kz,jmax,smax,ipmax),dqwa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dnf(kz,jmax,smax,ipmax),dqf(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dqfs(kz,jmax,smax,ipmax),dqfa(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dqfw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dni(kz,simax,itmax,ipmax),dqi(kz,simax,itmax,ipmax)      

INTEGER :: i,j,k,si,it,ip
DOUBLE PRECISION :: nnw(jmax,smax,ipmax),qqw(jmax,smax,ipmax)
DOUBLE PRECISION :: qqws(jmax,smax,ipmax),qqwa(jmax,smax,ipmax)
DOUBLE PRECISION :: nnf(jmax,smax,ipmax),qqf(jmax,smax,ipmax)
DOUBLE PRECISION :: qqfs(jmax,smax,ipmax),qqfa(jmax,smax,ipmax)
DOUBLE PRECISION :: qqfw(jmax,smax,ipmax)
DOUBLE PRECISION :: nni(simax,itmax,ipmax),qqi(simax,itmax,ipmax)
DOUBLE PRECISION :: ddnw(jmax,smax),ddqw(jmax,smax),ddqws(jmax,smax)
DOUBLE PRECISION :: ddqwa(jmax,smax),ddnf(jmax,smax),ddqf(jmax,smax)
DOUBLE PRECISION :: ddqfs(jmax,smax),ddqfa(jmax,smax),ddqfw(jmax,smax)
DOUBLE PRECISION :: ddni(simax,itmax),ddqi(simax,itmax)
DOUBLE PRECISION :: mmquer(jmax,smax),rrquer(jmax,smax),ssquer(jmax,smax)
DOUBLE PRECISION :: mmfquer(jmax,smax),rrfquer(jmax,smax),ssfquer(jmax,smax)
DOUBLE PRECISION :: vvtw(jmax),vvtf(jmax)

DOUBLE PRECISION :: TT(kz,ipmax)
DOUBLE PRECISION :: ddqq,ddqqdyn,dTT,dTTdyn,dpp,drho

DOUBLE PRECISION :: r_dry(jmax,smax)


! Temperaturen aus teta berechnen.
CALL teta2txd(teta,kz,ipmax,TT,pp,cappa,1)

DO ip=1,ipmax
  DO k=1,kz
!    WRITE(*,*) 'Schicht = ', k
    ddqq = 0.0D0
    ddqqdyn = dqqdyn(k,ip)
    dTT = 0.0D0
    dTTdyn = dtetadyn(k,ip)/(pp(1)/pp(k))**cappa  ! da p zeitl. konst
    drho = 0.0D0
    dpp = 0.0D0
    DO j=1,jmax
      vvtw(j) = vtw(k,j)
      vvtf(j) = vtf(k,j)
      DO i=1,smax
        nnw(j,i,ip) = nw(k,j,i,ip)
        qqw(j,i,ip) = qw(k,j,i,ip)
        qqws(j,i,ip) = qws(k,j,i,ip)
        qqwa(j,i,ip) = qwa(k,j,i,ip)
        nnf(j,i,ip) = nf(k,j,i,ip)
        qqf(j,i,ip) = qf(k,j,i,ip)
        qqfs(j,i,ip) = qfs(k,j,i,ip)
        qqfa(j,i,ip) = qfa(k,j,i,ip)
        qqfw(j,i,ip) = qfw(k,j,i,ip)
        ddnw(j,i) = 0.0D0
        ddqw(j,i) = 0.0D0
        ddqws(j,i) = 0.0D0
        ddqwa(j,i) = 0.0D0
        ddnf(j,i) = 0.0D0
        ddqf(j,i) = 0.0D0
        ddqfs(j,i) = 0.0D0
        ddqfa(j,i) = 0.0D0
        ddqfw(j,i) = 0.0D0
        mmquer(j,i) = mquer(k,j,i,ip)
        rrquer(j,i) = rquer(k,j,i,ip)
        ssquer(j,i) = squer(k,j,i,ip)
        mmfquer(j,i) = mfquer(k,j,i,ip)
        rrfquer(j,i) = rfquer(k,j,i,ip)
        ssfquer(j,i) = sfquer(k,j,i,ip)
      END DO
    END DO
    DO si=1,simax
      DO it=1,itmax
        nni(si,it,ip) = ni(k,si,it,ip)
        qqi(si,it,ip) = qi(k,si,it,ip)
        ddni(si,it) = 0.0D0
        ddqi(si,it) = 0.0D0
      END DO
    END DO
    CALL cloudxd(dTT,     &
                 dTTdyn,  &
                 dpp,     &
                 drho,    &
                 ddnw,    &
                 ddqw,    &
                 ddqws,   &
                 ddqwa,   &
                 ddqq,    &
                 ddqqdyn, &
                 usatt(k,ip), &
                 nnw,     &
                 qqw,     &
                 qqws,    &
                 qqwa,    &
                 qq(k,ip),&
                 pp(k),   &
                 rho(k),  &
                 TT(k,ip),&
                 mmquer,  &
                 rrquer,  &
                 ssquer,  &
                 mmfquer, &
                 rrfquer, &
                 ssfquer, &
                 ddnf,    &
                 ddqf,    &
                 ddqfs,   &
                 ddqfa,   &
                 ddqfw,   &
                 nnf,     &
                 qqf,     &
                 qqfs,    &
                 qqfa,    &
                 qqfw,    &
                 ddni,    &
                 ddqi,    &
                 nni,     &
                 qqi,     &
                 Talt(k,ip),  &
                 ifreeze, &
                 miv,     &
                 ip,      &
                 k,       &
                 vvtw,    &
                 vvtf)
    dqq(k,ip) = ddqq
    dteta(k,ip) = dTT*(pp(1)/pp(k))**cappa      ! da p zeitl. konst
    Talt(k,ip) = TT(k,ip)
    DO j=1,jmax
      DO i=1,smax
        dnw(k,j,i,ip) = ddnw(j,i)/deltat 
        dqw(k,j,i,ip) = ddqw(j,i)/deltat
        dqws(k,j,i,ip) = ddqws(j,i)/deltat
        dqwa(k,j,i,ip) = ddqwa(j,i)/deltat
        dnf(k,j,i,ip) = ddnf(j,i)/deltat
        dqf(k,j,i,ip) = ddqf(j,i)/deltat
        dqfs(k,j,i,ip) = ddqfs(j,i)/deltat
        dqfa(k,j,i,ip) = ddqfa(j,i)/deltat
        dqfw(k,j,i,ip) = ddqfw(j,i)/deltat
      END DO
    END DO
    DO si=1,simax
      DO it=1,itmax
        dni(k,si,it,ip) = ddni(si,it)/deltat
        dqi(k,si,it,ip) = ddqi(si,it)/deltat
      END DO
    END DO
  END DO
END DO

      
RETURN
END
