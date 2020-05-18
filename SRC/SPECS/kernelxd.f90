SUBROUTINE kernelxd
  
  ! Initialisation of the collision/coalescence kernel
  ! choice between Golovin and Kerkweg et al. 
  ! called by wmain.f
  ! calling none
  ! reading files "kernel*"
  ! to do: implement the different kernel averaging versions?
  ! => done
  
  ! Definition der Variablen 
  
  INCLUDE 'HEADER3'
  
  IMPLICIT NONE
  
  INTEGER :: II,JJ,igolov
  double precision :: sabi(JMAX,JMAX)
  
  igolov=2  ! Golovin kernel
  if(igolov.eq.1) then
     DO JJ=1,JMAX
        DO II=1,JMAX
           KOLK(JJ,II)  = 1.5D0*(MMITTE(JJ)+MMITTE(II))
           KOLKI(JJ,II) = 1.5D0*(MMITTE(JJ)+MMITTE(II))
        ENDDO
     ENDDO
     DO JJ=1,JMAX+1
        DO II=1,JMAX+1
           KOLK2D(JJ,II)  = 1.5D0*(MGRENZ(JJ)+MGRENZ(II))
           KOLKI2D(JJ,II) = 1.5D0*(MGRENZ(JJ)+MGRENZ(II))
        ENDDO
     ENDDO
  endif

  ! read in the kernel given as file (here: Kerkweg et al.) for drop-drop,ap
  if(igolov.eq.2) then

!    print *,'reading kernel drops'
     ! 2-dimensional
     if(JMAX.eq.66)  OPEN(21,FILE="DATEN/kernel_66")
     if(JMAX.eq.132) OPEN(21,FILE="DATEN/kernel_132")
     if(JMAX.eq.264) OPEN(21,FILE="DATEN/kernel_264")

     REWIND 21
     READ(21,*) KOLK
     CLOSE(21)

     DO JJ=1,JMAX
        DO II=1,JMAX
           if(KOLK(JJ,II).lt.0.D0) then
              write(*,*) "KOLK wrong ",JJ,II,KOLK(JJ,II)
              KOLK(JJ,II)=0.D0
           endif
        ENDDO
     ENDDO

     ! 2-dimensional at bin corners
     if(JMAX.eq.66)  OPEN(21,FILE="DATEN/kernel2d67")
     if(JMAX.eq.132) OPEN(21,FILE="DATEN/kernel2d133")
     if(JMAX.eq.264) OPEN(21,FILE="DATEN/kernel2d265")

     REWIND 21
     READ(21,*) KOLK2D
     CLOSE(21)
     
     DO JJ=1,JMAX+1
        DO II=1,JMAX+1
           if(KOLK2D(JJ,II).lt.0.D0) then
              write(*,*) "KOLK2D wrong ",JJ,II,KOLK2D(JJ,II)
              KOLK2D(JJ,II)=0.D0
           endif
        ENDDO
     ENDDO

!     print *,'reading kernel ice'
     ! 2-dimensional ice
     if(JMAX.eq.66)  OPEN(21,FILE="DATEN/kernel_66ice")
     if(JMAX.eq.132) OPEN(21,FILE="DATEN/kernel_132ice")
     if(JMAX.eq.264) OPEN(21,FILE="DATEN/kernel_264ice")

     REWIND 21
     READ(21,*) KOLKI
     CLOSE(21)

     DO JJ=1,JMAX
        DO II=1,JMAX
           if(KOLKI(JJ,II).lt.0.D0) then
              write(*,*) "KOLKI wrong ",JJ,II,KOLKI(JJ,II)
              KOLKI(JJ,II)=0.D0
           endif
           IF(jj.GE.61) kolki(jj,ii) = 0.0D0
           IF(ii.GE.61) kolki(jj,ii) = 0.0D0
           !!            IF(jj.GE.61) kolki(jj,ii) = kolki(jj,ii)/50.0D0
           !!            IF(ii.GE.61) kolki(jj,ii) = kolki(jj,ii)/50.0D0
        ENDDO
     ENDDO

!     print *,'reading kernel 2D'
     ! 2-dimensional at bin corners
     if(JMAX.eq.66) OPEN(21,FILE="DATEN/kernel2d67ice")
     if(JMAX.eq.132) OPEN(21,FILE="DATEN/kernel2d133ice")
     if(JMAX.eq.264) OPEN(21,FILE="DATEN/kernel2d265ice")

     REWIND 21
     READ(21,*) KOLKI2D
     CLOSE(21)

     DO JJ=1,JMAX+1
        DO II=1,JMAX+1
           if(KOLKI2D(JJ,II).lt.0.D0) then
              write(*,*) "KOLKI2D wrong ",JJ,II,KOLKI2D(JJ,II)
              KOLKI2D(JJ,II)=0.D0
           endif
           !!            IF(jj.GE.61) kolki2d(jj,ii) = 0.0D0
           !!            IF(ii.GE.61) kolki2d(jj,ii) = 0.0D0
           !!            IF(jj.GE.61) kolki2d(jj,ii) = kolki2d(jj,ii)/50.0D0
           !!            IF(ii.GE.61) kolki2d(jj,ii) = kolki2d(jj,ii)/50.0D0
        ENDDO
     ENDDO
  endif

  RETURN
END SUBROUTINE kernelxd
