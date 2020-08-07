C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C 
C     SUBROUTINE DIAGSP
C     =================
C     
C *************************************************************************
C
C  Diagsp creates an interface to the diagonalization routine ewevsp.
C
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix (eigenvectors in columns) 
C       EVAL    (O) :  Eigenvalues 
C       SC      (-) :  Auxiliary vector 
C       IEV     (I) :  should be 0 if eigenvectors are not required
C 
C *************************************************************************
C
      SUBROUTINE DIAGSP_D(NA,N,A,EVAL,SC,IEV)
       IMPLICIT  REAL*8 (A-H,O-Z)
       DIMENSION A(NA,N),EVAL(N),SC(N)
       SAVE
C
       IF (N .GT. NA) THEN
        PRINT *,'DIAGSP: PROBLEM SIZE EXCEEDS MATRIX DIMENSIONS'
        CALL STOPIT
       END IF
C
C  CALL EWEVGE WITH EIGENVECTORS AND ASCENDING SORT
C
       IER = 0
       IORD = -1
       CALL EWEVSP(NA,N,A,EVAL,SC,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGSP: FAILED IN EWEVGE, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
