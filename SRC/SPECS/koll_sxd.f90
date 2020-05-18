SUBROUTINE koll_sxd(XSUM1,XSUM2,XGRE,XNEU,NNEU,XDIF,XMIT,F1X,F1N,icall)

! calculation of the linear approximation and the factors for the 
! partition of the interaction terms to the target bins
!
! called by koll.f
!
! calling none

implicit none
integer , INTENT (IN) :: icall
DOUBLE PRECISION :: ug,og,n0,n1,n2,n0kx0,steig,steigstar,mstar,FN,FX
DOUBLE PRECISION , INTENT(OUT) :: F1N,F1X
DOUBLE PRECISION , INTENT(IN)  :: NNEU,XNEU,XDIF,XMIT,XSUM1,XSUM2,XGRE


! parameters of the linear distribution
n0=NNEU/XDIF
steig=12.D0*(XNEU-XMIT*NNEU)/XDIF**3
n0kx0=n0-steig*XMIT
n1=n0-steig*XDIF/2.D0
n2=n0+steig*XDIF/2.D0
! Looking for the target intervals
! n1 > 0 und n2 > 0
if(n1.ge.0.D0.and.n2.ge.0.D0) then
  ug=XSUM1
  og=XGRE
  FX=(0.5D0*n0kx0*(og**2-ug**2)+steig*(og**3-ug**3)/3.D0)
  FN=n0kx0*(og-ug)+0.5D0*steig*(og**2-ug**2)
! share of the source term for target bin K and S (between 0 und 1)
  F1X=FX/XNEU
  F1N=FN/NNEU
endif
! n1 < 0
if(n1.lt.0.D0) then
  mstar=3.D0*(XNEU/NNEU)-2.D0*XSUM2
  if(mstar.gt.XGRE) then
    F1X=0.D0
    F1N=0.D0
  else
    steigstar=2.D0*NNEU/(XSUM2-mstar)**2
    ug=mstar
    og=XGRE
    FX=steigstar*((og**3-ug**3)/3.D0-0.5D0*mstar*(og**2-ug**2))
    FN=steigstar*((og**2-ug**2)/2.D0-mstar*(og-ug))
! share of the source term for target bin K and S (between 0 und 1)
    F1X=FX/XNEU
    F1N=FN/NNEU
  endif
endif
! n2 < 0
if(n2.lt.0.D0) then
  mstar=3.D0*(XNEU/NNEU)-2.D0*XSUM1
  if(mstar.lt.XGRE) then
    F1X=1.D0
    F1N=1.D0
  else
    steigstar=-2.D0*NNEU/(XSUM1-mstar)**2
    ug=XSUM1
    og=XGRE
    FX=steigstar*((og**3-ug**3)/3.D0-0.5D0*mstar*(og**2-ug**2))
    FN=steigstar*((og**2-ug**2)/2.D0-mstar*(og-ug))
    F1X=FX/XNEU
    F1N=FN/NNEU
  endif
endif
! Test of F1X and F1N
if(F1X.lt.0.D0.or.F1X.gt.1.D0.or.F1N.lt.0.D0.or.F1N.gt.1.D0) then
!!!        write(*,*) "Fehler in koll_s",F1X,F1N,n1/abs(n1),n2/abs(n2)
!!!     &             ,icall
  F1X=MIN(F1X,1.D0)
  F1X=MAX(F1X,0.D0)
  F1N=MIN(F1N,1.D0)
  F1N=MAX(F1N,0.D0)
endif
      
!      if(F1X.gt.0.9D0.or.F1X.lt.0.1D0) write(81,*) F1X,F1N
!      if(F1X.le.0.9D0.and.F1X.ge.0.1D0) write(82,*) F1X,F1N

 
return
end
