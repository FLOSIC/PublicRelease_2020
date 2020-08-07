C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C  This is another version of Iqldia using a less sophisticated 
C  shifting algorithm. It is much simpler but 20 percent slower.
C
C *************************************************************************
C
C     SUBROUTINE IQLDIA(NM,N,D,E,Z,IEV,IER)
C      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION  D(N),E(N),Z(NM,N)
C      SAVE
C
C      IER = 0
C      IF (N .EQ. 1) RETURN
C      DO 10 I = 2, N
C       E(I-1) = E(I)
C  10  CONTINUE
C      E(N) = 0.0D0
C      DO 250 L = 1, N
C       ITER = 0
C
C  LOOK FOR SMALL SUBDIAGONAL ELEMENT 
C 
C 100   DO 110 M = L, N-1
C        DD = ABS(D(M)) + ABS(D(M+1))
C        IF ((ABS(E(M)) + DD) .EQ. DD) GOTO 120
C 110   CONTINUE
C       M = N
C 120   IF (M .EQ. L) GOTO 250
C       IF (ITER .EQ. 30) GOTO 900
C       ITER = ITER + 1
C
C  FORM SHIFT 
C
C       G = (D(L+1) - D(L)) / (2.0D0 * E(L))
C       R = SQRT (G * G + 1.0D0)
C       G = D(M) - D(L) + E(L) / (G + SIGN(R,G))
C       S = 1.0D0
C       C = 1.0D0
C       P = 0.0D0
C
C  FOR I = M - 1 STEP -1 UNTIL L DO 
C
C       DO 200 II = 1, M-L
C        I = M - II
C        F = S * E(I)
C        B = C * E(I)
C
C  SAFE CALCULATION OF SQRT(G * G + F * F) AND SIMILAR STUFF
C
C        IF (ABS(F) .LT. ABS(G)) GOTO 150
C        C = G / F
C        R = SQRT(1.0D0 + C * C)
C        E(I+1) = F * R
C        S = 1.0D0 / R
C        C = C * S
C        GOTO 160
C 150    S = F / G
C        R = SQRT (1.0D0 + S * S)
C        E(I+1) = G * R
C        C = 1.0D0 / R
C        S = S * C
C 160    G = D(I+1) - P
C        R = (D(I) - G) * S + 2.0D0 * C * B
C        P = S * R
C        D(I+1) = G + P
C        G = C * R - B
C        IF (IEV .EQ. 0) GOTO 200
C
C  FORM VECTOR
C
C        DO 180 K = 1, N
C         F = Z(K,I+1)
C         Z(K,I+1) = S * Z(K,I) + C * F
C         Z(K,I) =   C * Z(K,I) - S * F
C 180    CONTINUE
C 200   CONTINUE
C       D(L) = D(L) - P
C       E(L) = G
C       E(M) = 0.0D0
C       GOTO 100
C 250  CONTINUE
C      RETURN
C 900  IER = L
C      RETURN
C     END
C
C *************************************************************************
C
C     SUBROUTINE BACKTR
C     =================
C
C *************************************************************************
C
C  Backtr solves the system R * X = Y (R upper triangular matrix),
C  where the main diagonal of R is given inverted.
C
C     Parameters:
C
C       NR      (I) :  Dimension of R
C       NX      (I) :  Dimension of X 
C       NY      (I) :  Dimension of Y
C       N       (I) :  Dimension of problem 
C       M       (I) :  Number of columns in X and Y
C       H       (I) :  Auxiliary vector
C       R       (I) :  Matrix R (upper triangle)
C       X       (O) :  Matrix X (solution of system)
C       Y       (I) :  Matrix Y (right side)
C
C *************************************************************************
C
      SUBROUTINE BACKTR(NR,NX,NY,N,M,H,R,X,Y)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  R(NR,N),X(NX,M),Y(NY,M),H(N)
       SAVE
C
C  CALCULATION OF X = INV(R) * Y 
C
       DO 40 II = 1,N
        I = N + 1 - II
        I1 = I + 1
        D = R(I,I)
        DO 10 J= I,N
         H(J)= R(I,J)
   10   CONTINUE
        DO 30 J = 1,M
         S = Y(I,J)
         DO 20 K = I1,N
          S = S - H(K) * X(K,J)
   20    CONTINUE
         X(I,J) = S * D
   30   CONTINUE
   40  CONTINUE
       RETURN
      END
