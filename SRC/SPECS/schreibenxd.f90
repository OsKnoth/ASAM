SUBROUTINE schreibenxd(ww,teta,qq,u2,usatt,pp,nw,qw,nf,qf,qfw,ni,qi,ini,   &
                       rrnw,rrqw,rrnf,rrqf,rrw,rrf,kz,anz,jmax,smax,       &
                       ipmax,simax,itmax,dz,tk,zeit,cappa)

IMPLICIT NONE
INTEGER :: kz,anz,jmax,smax,ipmax,simax,itmax,tk,zeit
DOUBLE PRECISION :: dz,ini(kz,anz),cappa
DOUBLE PRECISION :: ww(kz,ipmax),teta(kz,ipmax),qq(kz,ipmax)
DOUBLE PRECISION :: u2(kz),usatt(kz,ipmax),pp(kz),TT(kz,ipmax)
DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: nf(kz,jmax,smax,ipmax),qf(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: ni(kz,simax,itmax,ipmax),qi(kz,simax,itmax,ipmax)
DOUBLE PRECISION :: rrnw(jmax),rrqw(jmax),rrnf(jmax),rrqf(jmax),rrw,rrf

INTEGER :: k,j,i,ip,si,it
DOUBLE PRECISION :: nwsum(kz,ipmax),qwsum(kz,ipmax),nwsum28(kz,ipmax)
DOUBLE PRECISION :: qwsum28(kz,ipmax),nfsum(kz,ipmax),qfsum(kz,ipmax)
DOUBLE PRECISION :: qfwsum(kz,ipmax),qfwsum28(kz,ipmax)
DOUBLE PRECISION :: nfsum28(kz,ipmax),qfsum28(kz,ipmax),nisum(kz,ipmax)
DOUBLE PRECISION :: qisum(kz,ipmax)
DOUBLE PRECISION :: rrnwsum,rrqwsum,rrnfsum,rrqfsum


OPEN(9,FILE='RESULTATE/Ta',POSITION='APPEND')
OPEN(10,FILE='RESULTATE/Tb',POSITION='APPEND')
OPEN(11,FILE='RESULTATE/wa',POSITION='APPEND')
OPEN(12,FILE='RESULTATE/tetaa',POSITION='APPEND')
OPEN(13,FILE='RESULTATE/qa',POSITION='APPEND')
OPEN(14,FILE='RESULTATE/u2',POSITION='APPEND')
OPEN(15,FILE='RESULTATE/dtetaa',POSITION='APPEND')
OPEN(16,FILE='RESULTATE/dqa',POSITION='APPEND')
OPEN(17,FILE='RESULTATE/satta',POSITION='APPEND')
OPEN(18,FILE='RESULTATE/wb',POSITION='APPEND')
OPEN(19,FILE='RESULTATE/nwas',POSITION='APPEND')
OPEN(20,FILE='RESULTATE/qwas',POSITION='APPEND')
OPEN(21,FILE='RESULTATE/nwas28',POSITION='APPEND')
OPEN(22,FILE='RESULTATE/qwas28',POSITION='APPEND')
OPEN(23,FILE='RESULTATE/nfas',POSITION='APPEND')
OPEN(24,FILE='RESULTATE/qfas',POSITION='APPEND')
OPEN(25,FILE='RESULTATE/nfas28',POSITION='APPEND')
OPEN(26,FILE='RESULTATE/qfas28',POSITION='APPEND')
OPEN(27,FILE='RESULTATE/rr',POSITION='APPEND')
OPEN(28,FILE='RESULTATE/nwa',POSITION='APPEND')
OPEN(29,FILE='RESULTATE/qwa',POSITION='APPEND')
OPEN(30,FILE='RESULTATE/nfa',POSITION='APPEND')
OPEN(31,FILE='RESULTATE/qfa',POSITION='APPEND')
OPEN(32,FILE='RESULTATE/nwb',POSITION='APPEND')
OPEN(33,FILE='RESULTATE/qwb',POSITION='APPEND')
OPEN(34,FILE='RESULTATE/nfb',POSITION='APPEND')
OPEN(35,FILE='RESULTATE/qfb',POSITION='APPEND')
OPEN(36,FILE='RESULTATE/nwbs',POSITION='APPEND')
OPEN(37,FILE='RESULTATE/qwbs',POSITION='APPEND')
OPEN(38,FILE='RESULTATE/nwbs28',POSITION='APPEND')
OPEN(39,FILE='RESULTATE/qwbs28',POSITION='APPEND')
OPEN(40,FILE='RESULTATE/nfbs',POSITION='APPEND')
OPEN(41,FILE='RESULTATE/qfbs',POSITION='APPEND')
OPEN(42,FILE='RESULTATE/nfbs28',POSITION='APPEND')
OPEN(43,FILE='RESULTATE/qfbs28',POSITION='APPEND')
OPEN(44,FILE='RESULTATE/sattb',POSITION='APPEND')
OPEN(45,FILE='RESULTATE/nias',POSITION='APPEND')
OPEN(46,FILE='RESULTATE/qias',POSITION='APPEND')
OPEN(47,FILE='RESULTATE/nibs',POSITION='APPEND')
OPEN(48,FILE='RESULTATE/qibs',POSITION='APPEND')
OPEN(49,FILE='RESULTATE/nia',POSITION='APPEND')
OPEN(50,FILE='RESULTATE/qia',POSITION='APPEND')
OPEN(51,FILE='RESULTATE/nib',POSITION='APPEND')
OPEN(52,FILE='RESULTATE/qib',POSITION='APPEND')
OPEN(53,FILE='RESULTATE/qfwas',POSITION='APPEND')
OPEN(54,FILE='RESULTATE/qfwas28',POSITION='APPEND')
OPEN(55,FILE='RESULTATE/qfwa',POSITION='APPEND')
OPEN(56,FILE='RESULTATE/qfwb',POSITION='APPEND')


! Anfangswerte schreiben.
IF(tk.EQ.1) THEN
  DO k=1,kz
    WRITE(9,100) 0, FLOAT(k)*dz, ini(k,10) - 273.16D0
    WRITE(10,100) 0, FLOAT(k)*dz, ini(k,10) - 273.16D0
    WRITE(11,100) 0, FLOAT(k)*dz, ini(k,1)
    WRITE(12,100) 0, FLOAT(k)*dz, ini(k,2)
    WRITE(13,100) 0, FLOAT(k)*dz, ini(k,3)*1000.0D0
    WRITE(14,100) 0, FLOAT(k)*dz, ini(k,5)
  END DO
  WRITE(9,*)
  WRITE(10,*)
  WRITE(11,*)
  WRITE(12,*)
  WRITE(13,*)
  WRITE(14,*)
END IF

CALL teta2txd(teta,kz,ipmax,TT,pp,cappa,1)

! Aktuelle Daten schreiben.
DO k=1,kz
  DO ip=1,ipmax
    nwsum(k,ip) = 0.0D0
    qwsum(k,ip) = 0.0D0
    nwsum28(k,ip) = 0.0D0
    qwsum28(k,ip) = 0.0D0
    nfsum(k,ip) = 0.0D0
    qfsum(k,ip) = 0.0D0
    qfwsum(k,ip) = 0.0D0
    nfsum28(k,ip) = 0.0D0
    qfsum28(k,ip) = 0.0D0
    qfwsum28(k,ip) = 0.0D0
    DO j=1,jmax    
      DO i=1,smax
        nwsum(k,ip) = nwsum(k,ip) + nw(k,j,i,ip)
        qwsum(k,ip) = qwsum(k,ip) + qw(k,j,i,ip)
        nfsum(k,ip) = nfsum(k,ip) + nf(k,j,i,ip)
        qfsum(k,ip) = qfsum(k,ip) + qf(k,j,i,ip)
        qfwsum(k,ip) = qfwsum(k,ip) + qfw(k,j,i,ip)
        IF(j.GE.28) THEN
          nwsum28(k,ip) = nwsum28(k,ip) + nw(k,j,i,ip)
          qwsum28(k,ip) = qwsum28(k,ip) + qw(k,j,i,ip)
          nfsum28(k,ip) = nfsum28(k,ip) + nf(k,j,i,ip)
          qfsum28(k,ip) = qfsum28(k,ip) + qf(k,j,i,ip)
          qfwsum28(k,ip) = qfwsum28(k,ip) + qfw(k,j,i,ip)
        END IF
      END DO
    END DO
    nisum(k,ip) = 0.0D0
    qisum(k,ip) = 0.0D0
    DO si=1,simax
      nisum(k,ip) = nisum(k,ip) + ni(k,si,1,ip)
      qisum(k,ip) = qisum(k,ip) + qi(k,si,1,ip)
    END DO
  END DO
  WRITE(9,100) tk, FLOAT(k)*dz, TT(k,1) - 273.16D0
  WRITE(10,100) tk, FLOAT(k)*dz, TT(k,2) - 273.16D0
  WRITE(11,100) tk, FLOAT(k)*dz, ww(k,1)
  WRITE(12,100) tk, FLOAT(k)*dz, teta(k,1)
  WRITE(13,100) tk, FLOAT(k)*dz, qq(k,1)*1000.0D0
  WRITE(14,100) tk-1, FLOAT(k)*dz, u2(k)
  WRITE(15,100) tk, FLOAT(k)*dz, teta(k,1) - ini(k,2)
  WRITE(16,100) tk, FLOAT(k)*dz, (qq(k,1) - ini(k,3))*1000.0D0
  WRITE(17,100) tk, FLOAT(k)*dz, usatt(k,1)*100.0D0
  WRITE(44,100) tk, FLOAT(k)*dz, usatt(k,2)*100.0D0
  WRITE(18,100) tk, FLOAT(k)*dz, ww(k,2)
  WRITE(19,100) tk, FLOAT(k)*dz, nwsum(k,1)*1.0D-6
  WRITE(20,100) tk, FLOAT(k)*dz, qwsum(k,1)*1000.0D0
  WRITE(21,100) tk, FLOAT(k)*dz, nwsum28(k,1)*1.0D-6
  WRITE(22,100) tk, FLOAT(k)*dz, qwsum28(k,1)*1000.0D0
  WRITE(23,100) tk, FLOAT(k)*dz, nfsum(k,1)*1.0D-6
  WRITE(24,100) tk, FLOAT(k)*dz, qfsum(k,1)*1000.0D0
  WRITE(25,100) tk, FLOAT(k)*dz, nfsum28(k,1)*1.0D-6
  WRITE(26,100) tk, FLOAT(k)*dz, qfsum28(k,1)*1000.0D0
  WRITE(36,100) tk, FLOAT(k)*dz, nwsum(k,2)*1.0D-6
  WRITE(37,100) tk, FLOAT(k)*dz, qwsum(k,2)*1000.0D0
  WRITE(38,100) tk, FLOAT(k)*dz, nwsum28(k,2)*1.0D-6
  WRITE(39,100) tk, FLOAT(k)*dz, qwsum28(k,2)*1000.0D0
  WRITE(40,100) tk, FLOAT(k)*dz, nfsum(k,2)*1.0D-6
  WRITE(41,100) tk, FLOAT(k)*dz, qfsum(k,2)*1000.0D0
  WRITE(42,100) tk, FLOAT(k)*dz, nfsum28(k,2)*1.0D-6
  WRITE(43,100) tk, FLOAT(k)*dz, qfsum28(k,2)*1000.0D0
  WRITE(45,100) tk, FLOAT(k)*dz, nisum(k,1)*1.0D-6
  WRITE(46,100) tk, FLOAT(k)*dz, qisum(k,1)*1000.0D0
  WRITE(47,100) tk, FLOAT(k)*dz, nisum(k,2)*1.0D-6
  WRITE(48,100) tk, FLOAT(k)*dz, qisum(k,2)*1000.0D0
  WRITE(53,100) tk, FLOAT(k)*dz, qfwsum(k,1)*1000.0D0
  WRITE(54,100) tk, FLOAT(k)*dz, qfwsum28(k,1)*1000.0D0
END DO
DO k=9,26
  WRITE(k,*)
  CLOSE(k)
END DO
DO k=36,48
  WRITE(k,*)
  CLOSE(k)
END DO
DO k=53,54
  WRITE(k,*)
  CLOSE(k)
END DO

! Niederschlag schreiben
rrnwsum = 0.0D0
rrqwsum = 0.0D0
rrnfsum = 0.0D0
rrqfsum = 0.0D0
DO j=1,jmax
  rrnwsum = rrnwsum + rrnw(j)
  rrqwsum = rrqwsum + rrqw(j)
  rrnfsum = rrnfsum + rrnf(j)
  rrqfsum = rrqfsum + rrqf(j)
END DO
WRITE(27,200) tk,rrnwsum*1.0D-6,rrqwsum*1000.0D0,rrw,rrnfsum*1.0D-6,   &
              rrqfsum*1000.0D0,rrf
CLOSE(27)

! Spektren schreiben
IF(tk.EQ.1) THEN
  DO k=28,35
    WRITE(k,300) zeit,kz,jmax,smax,dz
  END DO
  DO k=49,52
    WRITE(k,300) zeit,kz,simax,itmax,dz
  END DO
END IF
WRITE(28,400) (((nw(k,j,i,1), i=1,smax), j=1,jmax), k=1,kz)
WRITE(29,400) (((qw(k,j,i,1), i=1,smax), j=1,jmax), k=1,kz)
WRITE(30,400) (((nf(k,j,i,1), i=1,smax), j=1,jmax), k=1,kz)
WRITE(31,400) (((qf(k,j,i,1), i=1,smax), j=1,jmax), k=1,kz)
WRITE(55,400) (((qfw(k,j,i,1), i=1,smax), j=1,jmax), k=1,kz)
WRITE(32,400) (((nw(k,j,i,2), i=1,smax), j=1,jmax), k=1,kz)
WRITE(33,400) (((qw(k,j,i,2), i=1,smax), j=1,jmax), k=1,kz)
WRITE(34,400) (((nf(k,j,i,2), i=1,smax), j=1,jmax), k=1,kz)
WRITE(35,400) (((qf(k,j,i,2), i=1,smax), j=1,jmax), k=1,kz)
WRITE(56,400) (((qfw(k,j,i,2), i=1,smax), j=1,jmax), k=1,kz)
WRITE(49,400) (((ni(k,si,it,1), it=1,itmax), si=1,simax), k=1,kz)
WRITE(50,400) (((qi(k,si,it,1), it=1,itmax), si=1,simax), k=1,kz)
WRITE(51,400) (((ni(k,si,it,2), it=1,itmax), si=1,simax), k=1,kz)
WRITE(52,400) (((qi(k,si,it,2), it=1,itmax), si=1,simax), k=1,kz)
DO k=28,35
  CLOSE(k)
END DO
DO k=49,52
  CLOSE(k)
END DO
CLOSE(55)
CLOSE(56)

100 FORMAT(I10,4X,E13.6,4X,E14.6)
200 FORMAT(I10,6(4X,E13.6))
300 FORMAT(4(I4,2X),E13.6)
400 FORMAT(10(E13.6,2X))


RETURN
END
