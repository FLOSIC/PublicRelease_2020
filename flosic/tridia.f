C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE TRIDIA
C     =================
C
C *************************************************************************
C
C  Tridiagonalization of a given symmetric matrix A using Householder
C
C     Parameters:
C
C       NM      (I) :  Dimension of A 
C       N       (I) :  Dimension of problem
C       D       (O) :  Diagonal of tridiagonal matrix
C       E       (O) :  Subdiagonal of tridiagonal matrix (E(1) = 0.0)
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Transformation Matrix
C       IEV     (I) :  0: No eigenvectors
C
C *************************************************************************
C
      SUBROUTINE TRIDIA(NM,N,D,E,A,IEV)
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  A(NM,N),D(N),E(N)
       SAVE
C
       DO 100 I = 1,N
        D(I) = A(N,I)
  100  CONTINUE
       IF (N .EQ. 1) GOTO 510
C
C  FOR I = N STEP -1 UNTIL 2 DO
C
       DO 300 II = 2,N
        I = N + 2 - II
        L = I - 1
        H = 0.0D0
        SCALE = 0.0D0
        IF (L .LT. 2) GOTO 130
C
C  SCALE ROW
C
        DO 120 K = 1,L
         SCALE = SCALE + ABS(D(K))
  120   CONTINUE
C
        IF (SCALE .NE. 0.0D0) GOTO 140
  130   E(I) = D(L)
        DO 135 J = 1,L
         D(J) = A(L,J)
         A(I,J) = 0.0D0
         A(J,I) = 0.0D0
  135   CONTINUE
        GOTO 290
C
  140   DO 150 K = 1,L
         D(K) = D(K) / SCALE
         H = H + D(K) * D(K)
  150   CONTINUE 
        F = D(L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE * G
        H = H - F * G
        D(L) = F - G
C
C  FORM A * U
C
        DO 170 J = 1,L
         E(J) = 0.0D0
  170   CONTINUE
        DO 240 J = 1,L
         F = D(J)
         A(J,I) = F
         G = E(J) + A(J,J) * F
         JP1 = J + 1
         DO 200 K = JP1,L
          G = G + A(K,J) * D(K)
          E(K) = E(K) + A(K,J) * F
  200    CONTINUE 
         E(J) = G
  240   CONTINUE
C
C  FORM P
C
        F = 0.0D0
        DO 245 J = 1,L
         E(J) = E(J) / H
         F = F + E(J) * D(J)
  245   CONTINUE
        HH = F / (H + H)
C
C  FORM Q
C
        DO 250 J = 1,L
         E(J) = E(J) - HH * D(J)
  250   CONTINUE
C
C  FORM REDUCED A
C
        DO 280 J = 1,L
         F = D(J)
         G = E(J)
         DO 260 K = J,L
          A(K,J) = A(K,J) - F * E(K) - G * D(K)
  260    CONTINUE
         D(J) = A(L,J)
         A(I,J) = 0.0D0 
  280   CONTINUE 
C
C  DONE WITH THIS TRANSFORMATION
C 
  290   D(I) = H
  300  CONTINUE 
C
C  ACCUMULATION OF TRANSFORMATION MATRICES
C
       IF (IEV .EQ. 0) GOTO 600
       DO 500 I = 2,N
        L = I - 1
        A(N,L) = A(L,L)
        A(L,L) = 1.0D0
        H = D(I)
        IF (H .EQ. 0.0D0) GOTO 380
        DO 330 K = 1,L
         D(K) = A(K,I) / H
  330   CONTINUE
        DO 360 J = 1,L
         G = 0.0D0
         DO 340 K = 1,L
          G = G + A(K,I) * A(K,J)
  340    CONTINUE
         DO 350 K = 1,L
          A(K,J) = A(K,J) - G * D(K)
  350    CONTINUE
  360   CONTINUE
C        
  380   DO 400 K = 1,L
         A(K,I) = 0.0D0
  400   CONTINUE
  500  CONTINUE
  510  DO 520 I = 1,N
        D(I) = A(N,I)
        A(N,I) = 0.0D0
  520  CONTINUE
       GOTO 700
C
C  ONLY EIGENVALUES REQUIRED
C
  600  DO 610 I = 1,N
        D(I) = A(I,I)
  610  CONTINUE
C
  700  A(N,N) = 1.0D0
       E(1) = 0.0D0
       RETURN
      END
