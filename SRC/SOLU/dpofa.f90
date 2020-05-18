      SUBROUTINE DPOFA(A,LDA,N,INFO)                                    

      USE Kind_Mod

      IMPLICIT NONE

      INTEGER LDA,N,INFO                                                
      REAL(RealKind) A(LDA,N)                                         
!                                                                       
!     DPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      
!     MATRIX.                                                           
!                                                                       
!     DPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED            
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          
!     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR DPOFA) .                  
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE         
!                DIAGONAL AND UPPER TRIANGLE ARE USED.                  
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R 
!                WHERE  TRANS(R)  IS THE TRANSPOSE.                     
!                THE STRICT LOWER TRIANGLE IS UNALTERED.                
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   
!                                                                       
!        INFO    INTEGER                                                
!                = 0  FOR NORMAL RETURN.                                
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DDOT                                                         
!     FORTRAN DSQRT                                                     
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      REAL(RealKind) DDOT,T                                           
      REAL(RealKind) S                                                
      INTEGER J,K                                                   
!     BEGIN BLOCK WITH ...EXITS TO 40                                   
!                                                                       
!                                                                       
         DO J = 1, N                                                 
            INFO = J                                                    
            S = 0.0d0                                                   
            DO K = 1, J-1                                            
               T = A(K,J) - DDOT(K-1,A(1,K),1,A(1,J),1)                 
               T = T/A(K,K)                                             
               A(K,J) = T                                               
               S = S + T*T                                              
            END DO
            S = A(J,J) - S                                              
!     ......EXIT                                                        
            IF (S <= 0.0d0) EXIT
            INFO = 0
            A(J,J) = SQRT(S)                                           
         END DO
      END SUBROUTINE dpofa                                                               
