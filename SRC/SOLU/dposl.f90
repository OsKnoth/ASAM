      SUBROUTINE DPOSL(A,LDA,N,B)                                       

      USE Kind_Mod

      IMPLICIT NONE

      INTEGER LDA,N                                                     
      REAL(RealKind) A(LDA,N),B(N)                                    
!                                                                       
!     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     
!     SYSTEM A * X = B                                                  
!     USING THE FACTORS COMPUTED BY DPOCO OR DPOFA.                     
!                                                                       
!     ON ENTRY                                                          
!                                                                       
!        A       DOUBLE PRECISION(LDA, N)                               
!                THE OUTPUT FROM DPOCO OR DPOFA.                        
!                                                                       
!        LDA     INTEGER                                                
!                THE LEADING DIMENSION OF THE ARRAY  A .                
!                                                                       
!        N       INTEGER                                                
!                THE ORDER OF THE MATRIX  A .                           
!                                                                       
!        B       DOUBLE PRECISION(N)                                    
!                THE RIGHT HAND SIDE VECTOR.                            
!                                                                       
!     ON RETURN                                                         
!                                                                       
!        B       THE SOLUTION VECTOR  X .                               
!                                                                       
!     ERROR CONDITION                                                   
!                                                                       
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES            
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED    
!        CORRECTLY AND  INFO .EQ. 0 .                                   
!                                                                       
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 
!     WITH  P  COLUMNS                                                  
!           CALL DPOCO(A,LDA,N,RCOND,Z,INFO)                            
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          
!           DO 10 J = 1, P                                              
!              CALL DPOSL(A,LDA,N,C(1,J))                               
!        10 CONTINUE                                                    
!                                                                       
!     LINPACK.  THIS VERSION DATED 08/14/78 .                           
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      
!                                                                       
!     SUBROUTINES AND FUNCTIONS                                         
!                                                                       
!     BLAS DAXPY,DDOT                                                   
!                                                                       
!     INTERNAL VARIABLES                                                
!                                                                       
      REAL(RealKind) DDOT,T                                           
      INTEGER K,KB                                                      
!                                                                       
!     SOLVE TRANS(R)*Y = B                                              
!                                                                       
      DO K = 1, N                                                    
         T = DDOT(K-1,A(1,K),1,B(1),1)                                  
         B(K) = (B(K) - T)/(A(K,K)+1.0d-40)                                       
      END DO
!                                                                       
!     SOLVE R*X = Y                                                     
!                                                                       
      DO KB = 1, N                                                   
         K = N + 1 - KB                                                 
         B(K) = B(K)/(A(K,K)+1.0d-40)                                             
         T = -B(K)                                                      
         CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                              
      END DO
      END SUBROUTINE dposl                                                               
