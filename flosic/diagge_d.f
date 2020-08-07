C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C  This is a collection of subroutines designated to solve the real
C  general symmetric eigenvalue problem with or without eigenvectors.
C  Some of the routines have been taken from different freeware FORTRAN
C  libraries and optimized by hand (or eye ?! ;-)). The routine ewevsvd
C  has been written by Dirk Porezag. Most of the optimizations have been
C  done with respect to stride minimization for the innermost loops of 
C  the subroutines. Problems with bugs, roaches and other livestock 
C  please report to
C
C  Dirk Porezag   porezag@physik.tu-chemnitz.de
C
C  or to your nearest pest control agency (I doubt they will help).
C  Have fun !!
C
C  Copyright for this file by Dirk Porezag
C  Washington, DC, May 12th, 1995 
C
C *************************************************************************
C 
C     SUBROUTINE DIAGGE
C     =================
C     
C *************************************************************************
C
C  Diagge creates an interface to the diagonalization routine ewevge.
C
C     Parameters: 
C 
C       NA      (I) :  Dimension of A and B 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix (eigenvectors in columns) 
C       B       (I) :  Matrix B (lower triangle)
C               (O) :  R where B = R'*R (upper triangle) 
C       EVAL    (O) :  Eigenvalues 
C       SC      (-) :  Auxiliary vector 
C       IEV     (I) :  should be 0 if eigenvectors are not required
C 
C *************************************************************************
C
      SUBROUTINE DIAGGE_D(NA,N,A,B,EVAL,SC,IEV)
!     SUBROUTINE DIAGGE(NA,N,A,B,EVAL,SC,IEV)
       IMPLICIT  REAL*8 (A-H,O-Z)
       DIMENSION A(NA,N),B(NA,N),EVAL(N),SC(N)
       SAVE
C
       IF (N .GT. NA) THEN
        PRINT *,'DIAGGE: PROBLEM SIZE EXCEEDS MATRIX DIMENSIONS'
        CALL STOPIT
       END IF
C
C  CALL EWEVGE WITH EIGENVECTORS AND ASCENDING SORT
C
       IER = 0
       IORD = -1
       CALL EWEVGE(NA,NA,N,A,B,EVAL,SC,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGGE: FAILED IN EWEVGE, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
