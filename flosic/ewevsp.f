C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE EWEVSP 
C     ================= 
C 
C *************************************************************************
C 
C  Ewevsp calculates eigenvalues and eigenvectors of the real special 
C  symmetric eigenvalue problem. 
C  In the explanation of the method, primes denote an exchange of rows 
C  and columns.
C 
C  Method:  *  A*C = E*C 
C           *  Solve INV(R')*A*INV(R)*Y = E*Y  (Householder + IQL) 
C           *  Sorting of eigenvalues and eigenvectors 
C 
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix  
C       EW      (O) :  Eigenvalues 
C       H       (-) :  Auxiliary vector 
C       IEV     (I) :  0: No eigenvectors  
C       IORD    (I) :  1: Descending order of eigenvalues 
C                     -1: Ascending order of eigenvalues
C                      otherwise: no sorting 
C       IER     (O) :  Error indication  
C                      0: No error 
C                     -1: N is larger than NA 
C                      K: Convergence failure for eigenvalue
C                         K, (K-1) eigenvalues are correct
C 
C *************************************************************************
C
      SUBROUTINE EWEVSP(NA,N,A,EW,H,IEV,IORD,IER) 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  A(NA,N),EW(N),H(N)
       SAVE
C
       IER = 0
       IF (N .LE. 0) RETURN
       IER = -1 
       IF (N .GT. NA) RETURN
       IER = 0 
       JEV = 0 
       IF (IEV .NE. 0) JEV = 1 
       CALL TRIDIA(NA,N,EW,H,A,JEV) 
       CALL IQLDIA(NA,N,EW,H,A,JEV,IER) 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC(NA,N,N,EW,H,A,IORD,JEV)
       RETURN 
      END
