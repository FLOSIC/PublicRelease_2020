C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE EWEVGE 
C     ================= 
C 
C *************************************************************************
C 
C  Ewevge calculates eigenvalues and eigenvectors of the real general 
C  symmetric eigenvalue problem. The Matrix B must be positiv definite.
C  In the explanation of the method, primes denote an exchange of rows 
C  and columns.
C 
C  Method:  *  A*C = E*B*C 
C           *  Choleski decomposition  B = R'*R
C           *  A*C = E*R'*R*C  ->  INV(R')*A*C = E*R*C
C           *  Transformation Y = R*C  ->  C = INV(R)*Y
C           *  Solve INV(R')*A*INV(R)*Y = E*Y  (Householder + IQL) 
C           *  Sorting of eigenvalues and eigenvectors 
C           *  Back transformation C = INV(R)*Y
C 
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       NB      (I) :  Dimension of B 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix  
C       B       (I) :  Matrix B (lower triangle)
C               (O) :  R where B = R'*R (upper triangle) 
C       EW      (O) :  Eigenvalues 
C       H       (-) :  Auxiliary vector 
C       IEV     (I) :  0: No eigenvectors  
C       IORD    (I) :  1: Descending order of eigenvalues 
C                     -1: Ascending order of eigenvalues
C                      otherwise: no sorting 
C       IER     (O) :  Error indication  
C                      0: No error 
C                     -1: N is larger than NA or NB
C                      K: (K <= N)  B is not positive definite
C                      K: (K > N) Convergence failure for eigenvalue
C                                 (K-N), (K-N-1) eigenvalues are correct
C 
C *************************************************************************
C
      SUBROUTINE EWEVGE(NA,NB,N,A,B,EW,H,IEV,IORD,IER) 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  A(NA,N),B(NB,N),EW(N),H(N)
       SAVE
C
       IER = 0
       IF (N .LE. 0) RETURN
       IER = -1 
       IF ((N .GT. NA) .OR. (N .GT. NB)) RETURN
       IER = 0 
       JEV = 0 
       IF (IEV .NE. 0) JEV = 1 
       CALL CHOLES(NB,N,B,IER) 
       IF (IER .NE. 0) RETURN 
       CALL MATRAF(NA,NB,N,H,A,B) 
       CALL TRIDIA(NA,N,EW,H,A,JEV) 
       CALL IQLDIA(NA,N,EW,H,A,JEV,IER) 
       IF (IER .NE. 0) IER = IER+N 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC(NA,N,N,EW,H,A,IORD,JEV)
       IF (JEV .EQ. 0) RETURN
       CALL BACKTR(NB,NA,NA,N,N,H,B,A,A) 
       RETURN 
      END
