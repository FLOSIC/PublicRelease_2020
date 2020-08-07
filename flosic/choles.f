C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE CHOLES
C     =================
C
C *************************************************************************
C
C  Choles calculates the Choleski decomposition B = R' * R of B
C  into an upper triangle matrix R for the symmetric positive
C  definite Matrix B. The elements of the main diagonal are
C  stored inverted.
C
C     Parameters:
C
C       NB      (I) :  Dimension of B 
C       N       (I) :  Dimension of problem 
C       B       (I) :  Matrix B (lower triangle)
C               (O) :  Matrix R (upper triangle), inverted main diagonal
C       ICHO    (I) :  ICHO - 1 is the dimension of the submatrix that
C                      is available as Choleski decomposition ( < 1 = 1)
C               (O) :  Row number where decomposition failed (0 if success)
C
C *************************************************************************
C
      SUBROUTINE CHOLES(NB,N,B,ICHO)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  B(NB,N)
       SAVE
C
       IF (ICHO .GT. N)  GOTO 200
       IF (ICHO .LT. 1)  ICHO = 1
       D=0.0D0
       DO 80 I = ICHO,N
        I1 = I - 1
        DO 70 J = I,N
         S = B(J,I)
         DO 20 K = 1,I1
          S = S - B(K,I) * B(K,J)
   20    CONTINUE
         IF (I .NE .J) GOTO 40
         IF (S .LE. 0.0D0) GOTO 100
         S = 1.0D0 / SQRT(S)
         D = S
         GOTO 60
   40    S = S * D
   60    B(I,J) = S
   70   CONTINUE
   80  CONTINUE
       ICHO = 0
       GOTO 200
  100  ICHO = I
  200  RETURN
      END
