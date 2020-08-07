C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE MATRAF
C     =================
C
C *************************************************************************
C
C  Matraf calculates out of the symmetric matrix A and the 
C  upper triangular matrix R the product INV(R') * A * INV(R), 
C  where the main diagonal of R is given inverted.
C
C     Parameters:
C
C       NA      (I) :  Dimension of A
C       NB      (I) :  Dimension of B
C       N       (I) :  Dimension of problem
C       H       (-) :  Auxiliary vector
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Transformed matrix (lower triangle)
C       B       (I) :  Matrix R (upper triangle), inverted main diagonal
C
C *************************************************************************
C
      SUBROUTINE MATRAF(NA,NB,N,H,A,B)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  A(NA,N),B(NB,N),H(N)
       SAVE
C
C  FILL MATRIX 
C
       DO 20 I = 1,N
        DO 10 J = I+1,N
         A(I,J) = A(J,I)
   10   CONTINUE
   20  CONTINUE
C
C  CALCULATION OF A = INV(R') * A
C
       DO 60 I = 1,N
        I1 = I-1
        D = B(I,I)
        DO 50 J = 1,N
         S = A(I,J)
         DO 30 K = 1,I1
          S = S - B(K,I) * A(K,J)
   30    CONTINUE
         A(I,J) = S * D
   50   CONTINUE  
   60  CONTINUE
C
C  CALCULATION OF A = A * INV(R) (USE BUFFER FOR STRIDE OPTIMIZATION)
C
       DO 160 I = 1,N
        I1 = I-1
        D = B(I,I)
        DO 110 J = I,N
         H(J) = A(J,I)
  110   CONTINUE              
        DO 130 K = 1,I1
         S = B(K,I)
         DO 120 J = I,N
          H(J) = H(J) - S * A(J,K)
  120    CONTINUE
  130   CONTINUE
        DO 140 J = I,N
         A(J,I) = H(J) * D
  140   CONTINUE              
  160  CONTINUE
       RETURN
      END
