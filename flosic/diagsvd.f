C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C 
C     SUBROUTINE DIAGSVD
C     ================== 
C     
C *************************************************************************
C
C  Diagsvd creates an interface to the diagonalization routine ewevsvd.
C
C     Parameters: 
C 
C       NA      (I) :  Dimension of A and B 
C       N       (I) :  Dimension of Problem  
C       NE      (I) :  Returned number of eigenvalues
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix (eigenvectors in columns) 
C       B       (I) :  Matrix B (lower triangle)
C               (O) :  Rubbish  
C       EVAL    (O) :  Eigenvalues 
C       SC1     (-) :  Auxiliary vector 1
C       SC2     (-) :  Auxiliary vector 2
C       DELTA   (I) :  Smallest allowed eigenvalue for matrix B
C       IEV     (I) :  should be 0 if eigenvectors are not required
C 
C *************************************************************************
C
      SUBROUTINE DIAGSVD(NA,N,NE,A,B,EVAL,SC1,SC2,DELTA,IEV)
       IMPLICIT  REAL*8 (A-H,O-Z)
       DIMENSION A(NA,N),B(NA,N),EVAL(N),SC1(N),SC2(N)
       SAVE
C
       IF (N .GT. NA) THEN
        PRINT *,'DIAGSVD: PROBLEM SIZE EXCEEDS MATRIX DIMENSIONS'
        CALL STOPIT
       END IF
C
C  CALL EWEVSVD WITH EIGENVECTORS AND ASCENDING SORT
C
       IER = 0
       IORD = -1
       CALL EWEVSVD(NA,NA,N,NE,A,B,EVAL,SC1,SC2,
     &              DELTA,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGSVD: FAILED IN EWEVSVD, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
