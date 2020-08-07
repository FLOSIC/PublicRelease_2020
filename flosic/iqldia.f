C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE IQLDIA
C     =================
C
C *************************************************************************
C
C  Iqldia calculates eigenvalues and eigenvectors of a tridiagonal
C  matrix using the QL algorithm with implicit shifting.
C
C     Parameters:
C
C       NM      (I) :  Dimension of Z
C       N       (I) :  Dimension of the problem 
C       D       (I) :  Diagonal of tridiagonal matrix
C               (O) :  Eigenvalues
C       E       (I) :  Subdiagonal of tridiagonal matrix
C       Z       (I) :  Transformation matrix 
C               (O) :  Eigenvectors according to Z
C       IEV     (I) :  0: No eigenvectors
C       IER     (O) :  Error indication
C                      0: no error
C                      K: Convergence failure for the eigenvalue 
C                         number k, k-1 eigenvalues are correct
C
C *************************************************************************
C
      SUBROUTINE IQLDIA(NM,N,D,E,Z,IEV,IER)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  D(N),E(N),Z(NM,N)
       SAVE
       DATA ACCU /1.0D-15/
C
       IER = 0
       IF (N .EQ. 1) RETURN
C
C  GET MACHINE EPSILON AND BIG
C
       CALL MACHEPS(EPS)
       EPS = MAX(EPS,ACCU)
       EPSS = SQRT(EPS)
       EPS4 = EPS * 1.0D-4
       BIG = 1.0D0/EPS4
C
       ANORM = 0.0D0  
       R = 0.0D0
       DO 30 I = 2, N
        S = E(I)
        E(I-1) = S
        S = ABS(S)
        P = ABS(D(I-1)) + R + S
        IF (P .GT. ANORM) ANORM = P
        R = S
   30  CONTINUE
       P = ABS(D(N)) + R
       IF (P .GT. ANORM) ANORM = P
       E(N) = 0.0D0
       DO 250 L = 1, N
        J = 0
C
C  LOOK FOR SMALL SUBDIAGONAL ELEMENT 
C 
   50   DO 60 M = L, N-1
         DD = ABS(D(M)) + ABS(D(M+1))
         IF (ABS(E(M)) .LE. (EPS * DD)) GOTO 70
         IF (ABS(E(M)) .LE. (EPS4 * ANORM)) GOTO 70
   60   CONTINUE
        M = N
   70   P = D(L)
        MM1 = M - 1
        IF (M .EQ. L) GOTO 250
        IF (J .EQ. 30) GOTO 900
        J = J + 1
C
C  FORM SHIFT. THIS IS A SLIGHTLY ADVANCED FORM OF SHIFTING MAKING
C  THE ROUTINE ABOUT 20 PERCENT FASTER THAN THE USUAL STUFF.
C
        G = (D(L+1) - P) / (2.0D0 * E(L))
        R = SQRT (G * G + 1.0D0)
        S = P - E(L) / (G + SIGN (R, G))
        IF (M .EQ. L+1) GOTO 120
        T = S
        R = MAX(ABS(S),(ANORM / N))
        DO 100 I = 1, 6
         PSI = D(M) - T
         PSJ = -1.0D0
         DO 90 KK = L, MM1
          K = L + MM1 - KK
          IF (ABS(PSI) .GE. (EPS * ABS(E(K)))) GOTO 80
          PSI = BIG
          PSJ = BIG * BIG
          GOTO 90
   80     P = E(K) / PSI
          PSI = D(K) - T - P * E(K)
          PSJ = P * P * PSJ - 1.0D0
   90    CONTINUE
         IF (ABS(PSJ) .LE. EPS4) GOTO 120
         P = PSI / PSJ
         C = P
         IF (ABS(P) .GT. (0.5D0 * R)) C = SIGN(R,P)
         T = T - C
         IF (ABS(P) .LE. (EPSS * R)) GOTO 110
  100   CONTINUE
        GOTO 120
  110   S = T
  120   G = D(M) - S
        S = 1.0D0
        C = 1.0D0
        P = 0.0D0
        MML = M - L
C
C  FOR I = M - 1 STEP -1 UNTIL L DO 
C
        DO 200 II = 1, MML
         I = M - II
         F = S * E(I)
         B = C * E(I)
C
C  SAFE CALCULATION OF SQRT(G * G + F * F) AND SIMILAR STUFF
C
         IF (ABS(F) .LT. ABS(G)) GOTO 150
         C = G / F
         R = SQRT(1.0D0 + C * C)
         E(I+1) = F * R
         S = 1.0D0 / R
         C = C * S
         GOTO 160
  150    S = F / G
         R = SQRT (1.0D0 + S * S)
         E(I+1) = G * R
         C = 1.0D0 / R
         S = S * C
  160    G = D(I+1) - P
         R = (D(I) - G) * S + 2.0D0 * C * B
         P = S * R
         D(I+1) = G + P
         G = C * R - B
         IF (IEV .EQ. 0) GOTO 200
C
C  FORM VECTOR
C
         DO 180 K = 1,N
          F = Z(K,I+1)
          B = Z(K,I)
          Z(K,I+1) = S * B + C * F
          Z(K,I)   = C * B - S * F
  180    CONTINUE
  200   CONTINUE
        D(L) = D(L) - P
        E(L) = G
        E(M) = 0.0D0
        GOTO 50
  250  CONTINUE
       RETURN
  900  IER = L
       RETURN
      END
