      SUBROUTINE DGEFA (A, LDA, N,INFO)
      IMPLICIT NONE
      INTEGER LDA,N,INFO
      REAL(8) :: A(LDA,N)
!
      REAL(8) :: T
      INTEGER :: J,K,KP1,NM1
      NM1 = N - 1
      DO K = 1, NM1
         KP1 = K + 1
         T = 1.0D0/A(K,K)
         CALL DSCAL(N-K,T,A(K+1,K),1)
         DO J = KP1, N
           T = -A(K,J)
           CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
         END DO
      END DO
      END SUBROUTINE DGEFA


      SUBROUTINE DGEFAP (A, LDA, N,IP, INFO)
      IMPLICIT NONE
      INTEGER LDA,N,IP(N),INFO
      REAL(8) :: A(LDA,N)
!
      REAL(8) ::t,Temp 
      INTEGER :: ii,jj,kk,kp,iTemp
      DO kk=1,n
        ip(kk)=kk
      END DO
      DO kk=1,n
!       Bestimme Pivot Element
        temp=ABS(A(ip(kk),ip(kk)))
        kp=kk
        DO jj=kk+1,n
          IF (temp<ABS(A(ip(jj),ip(jj)))) THEN
            temp=ABS(A(ip(jj),ip(jj)))
            kp=jj
          END IF
        END DO
        IF (temp>0.0d0) THEN
          itemp=ip(kk)
          ip(kk)=ip(kp)
          ip(kp)=itemp
          t=A(ip(kk),ip(kk))
          DO jj=kk+1,n
            A(ip(jj),ip(kk))=A(ip(jj),ip(kk))/t
          END DO
          DO jj=kk+1,n
            t=-A(ip(kk),ip(jj)) 
            DO ii=kk+1,n
              a(ip(ii),ip(jj))=a(ip(ii),ip(jj))+t*a(ip(ii),ip(kk))
            END DO
          END DO
        ELSE
          A(ip(kk),ip(kk))=1.0d0 
        END IF
      END DO
      END SUBROUTINE DGEFAP


