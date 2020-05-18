subroutine kolk_avexd(KOLK_eff,J1,I1,J2,I2,MQ1,SQ1,MQ2,SQ2,ikolk,icall)
  
  INCLUDE 'HEADER3'
  
  implicit none
  
  integer :: KK,II,JJ,IL,IS
  INTEGER , INTENT (IN) :: ikolk,icall,J1,I1,J2,I2
  DOUBLE PRECISION , INTENT(OUT) :: KOLK_eff
  DOUBLE PRECISION , INTENT( IN) :: MQ1,SQ1,MQ2,SQ2
  DOUBLE PRECISION :: MS1,MS2,MSL,MSS,PBIN,PBIN_INV
  DOUBLE PRECISION :: KOLK11,KOLK12,KOLK21,KOLK22
  DOUBLE PRECISION :: f11,f12,f21,f22,f1L,f2L,f1S,f2S
  
  KOLK_eff = 0.D0

  if(ikoll.lt.3) then
     ! mass bin of 1st particle
     MS1 = MQ1+SQ1
     !        do KK=J1+1,JMAX+1
     do KK = 2,JMAX+1
        if(MS1.lt.MGRENZ(KK)) then
           JJ = KK-1
           !            write(55,*) JJ,J1
           goto 100
        endif
     enddo
     JJ = JMAX

100  continue
!!!        if(MS1.lt.MGRENZ(JJ).or.MS1.gt.MGRENZ(JJ+1)) then
!!!          write(*,*) "kolk_ave.f: MS1 ",J1,I1,JJ,MS1,MQ1,SQ1
     !     &        ,MGRENZ(JJ),MGRENZ(JJ+1),MQ1/MGRENZ(J1),SQ1/SGRENZ(I1)
!!!     &        ,MS1/MGRENZ(JJ),icall
!!!        endif
     ! mass bin of 2nd particle
     MS2 = MQ2+SQ2
     !        do KK=J2+1,JMAX+1
     do KK=2,JMAX+1
        if(MS2.lt.MGRENZ(KK)) then
           II = KK-1
           !            write(56,*) II,J2
           goto 200
        endif
     enddo
     II=JMAX
200  continue
!!!        if(MS2.lt.MGRENZ(II).or.MS2.gt.MGRENZ(II+1)) then
!!!          write(*,*) "kolk_ave.f: MS2 ",J2,I2,II,MS2,MQ2,SQ2
!!!     &       ,MS2/MGRENZ(II),icall
!!!        endif
     ! Calculation of weighting factors
     if(ikoll.eq.1) then
        ! bin resolution
        PBIN     = 2.D0**(1.D0/(NMAX-1.D0))
        PBIN_INV = 1.D0/(PBIN-1.D0)
        ! linear weighting
        f1L = PBIN_INV*(PBIN-MS1/MGRENZ(JJ))
        f2L = 1.D0-f1L
        f1S = PBIN_INV*(PBIN-MS2/MGRENZ(II))
        f2S = 1.D0-f1S
        ! area weighting
        f11 = f1L*f1S
        f21 = f2L*f1S
        f12 = f1L*f2S
        f22 = f2L*f2S
        IF(f11.lt.0.D0.or.f21.lt.0.D0.or.f12.lt.0.D0.or.f22.lt.0.D0) THEN 
           write(*,*) JJ,II,"kolk_ave.f: f?? ",f11,f21,f12,f22           &
                ,f1L,f2L,f1S,f2S,icall
           f11 = 0.25D0
           f21 = 0.25D0
           f12 = 0.25D0
           f22 = 0.25D0
        ENDIF


        IF((f11+f21+f12+f22).gt.1.00001D0) THEN 
           write(*,*) "f sum !!Fehler!!!!!!!",f11+f21+f12+f22            &
                ,f11,f21,f12,f22,icall
           f11 = 0.25D0
           f21 = 0.25D0
           f12 = 0.25D0
           f22 = 0.25D0
        ENDIF

     else
        f11 = 0.25D0
        f21 = 0.25D0
        f12 = 0.25D0
        f22 = 0.25D0
        if(ikoll.ne.2) write(*,*) "kolk_ave.f: Fehler in ikoll"
     endif
     ! kernel at the corners of the averaging area
     if(ikolk.eq.0) then
        KOLK11 = KOLK2D(JJ,II) 
        KOLK12 = KOLK2D(JJ,II+1) 
        KOLK21 = KOLK2D(JJ+1,II) 
        KOLK22 = KOLK2D(JJ+1,II+1) 
     endif
     if(ikolk.eq.1) then
        KOLK11 = KOLKI2D(JJ,II) 
        KOLK12 = KOLKI2D(JJ,II+1) 
        KOLK21 = KOLKI2D(JJ+1,II) 
        KOLK22 = KOLKI2D(JJ+1,II+1) 
     endif
     ! result of area averaging
     KOLK_eff = f11*KOLK11+f12*KOLK12+f21*KOLK21+f22*KOLK22


  endif
  
  if(ikoll.eq.3) then
     if(ikolk.eq.0) KOLK_eff = KOLK(J1,J2)
     if(ikolk.eq.1) KOLK_eff = KOLKI(J1,J2)
  endif


!if (icall==2) print *,'kolk_avexd.f90, kolk_eff',kolk_eff
  
  return
end subroutine kolk_avexd
