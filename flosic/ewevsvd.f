C UTEP Electronic Structure Lab (2020)
C
C *************************************************************************
C
C     SUBROUTINE EWEVSVD 
C     ================== 
C 
C *************************************************************************
C 
C  Ewevsvd calculates eigenvalues and eigenvectors of the real general 
C  symmetric eigenvalue problem. The diagonalization is performed in
C  the subspace where B is positiv definite. This scheme is called
C  Singular Value Decomposition (SVD). B can have negative eigenvalues, 
C  but the matrix must have positive diagonal elements.
C  For the solution of the real symmetric eigenvalue problem using 
C  Householder and implicit QL algorithms see the documentation for 
C  routine ewevge. In the explanation of the method, simple vectors
C  (such as the vector of eigenvalues) are considered to be diagonal
C  matrices. Primes represent an exchange of rows and columns.
C  
C
C  Method:  *  A*C = B*C*E 
C           *  Set P = diag(1/SQRT(S(I,I))) and C = P*D
C              -> (P*A*P)*D = (P*B*P)*D*E
C           *  Set H = P*A*P and S = P*B*P -> H*D = S*D*E
C              -> S is now a normalized matrix
C           *  Find eigenvalues Z of S: S*Q = Q*Z
C           *  Use all eigenvectors Q with Z > Delta as basis to expand
C              D = Q*F ->  F = D'*Q  (Q is an orthogonal matrix)
C              -> H*Q*F = S*Q*F*E = Q*Z*F*E 
C              -> Q'*H*Q*F = Z*F*E
C           *  Set Z = Y*Y and X = Y**(-1) and G = Y*F -> F = X*G  
C              -> X*Q'*H*Q*X*G = G*E
C           *  Set U = X*Q'*H*Q*X*G -> U*G = G*E
C           *  Find eigenvalues of U: U*G = G*E
C           *  Back transformation: C = P*D = P*Q*F = P*Q*X*G
C 
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       NB      (I) :  Dimension of B 
C       N       (I) :  Dimension of Problem  
C       NE      (O) :  Number of calculated eigenvalues
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix  
C       B       (I) :  Matrix B (lower triangle)
C               (O) :  Rubbish
C       EW      (O) :  Eigenvalues 
C       H1      (-) :  Auxiliary vector 1
C       H2      (-) :  Auxiliary vector 2
C       DELTA   (I) :  Smallest allowed eigenvalue for matrix B
C       IEV     (I) :  0: No eigenvectors  
C                      1: Calculate eigenvectors
C       IORD    (I) :  1: Descending order of eigenvalues 
C                     -1: Ascending order of eigenvalues
C                      otherwise: no sorting 
C       IER     (O) :  Error indication  
C                      0: No error 
C                     -1: N larger than NA or NB
C                     -2: Delta is too small
C                     -3: Diagonal elements of B are not positive
C                      K <= N: Convergence failure for eigenvalue K of B
C                      K > N:  Convergence failure for eigenvalue K-N
C 
C *************************************************************************
C
      SUBROUTINE EWEVSVD(NA,NB,N,NE,A,B,EW,H1,H2,DELTA,IEV,IORD,IER) 
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER(ZERO=0.0D0)
       DIMENSION A(NA,N),B(NB,N),EW(N),H1(N),H2(N)
       SAVE
C
C  CHECK INPUT AND SETUP H1
C
       NE = 0
       IER = 0
       IF (N .LE. 0) RETURN
       IER = -1 
       IF ((N .GT. NA) .OR. (N .GT. NB)) RETURN
       IER = -2 
       IF (DELTA .LE. ZERO) RETURN
       IER = -3
       DO 10 I = 1,N
        IF (B(I,I) .LE. ZERO) RETURN
        H1(I) = 1.0D0 / SQRT(B(I,I))
   10  CONTINUE
       IER = 0
C
C  TRANSFORM A = (H1 * A * H1) AND B = (H1 * A * H1) 
C  AFTER COMPLETION, ALL DIAGONAL ELEMENTS OF S ARE NORMALIZED
C
       DO 30 I = 1,N
        FAC = H1(I)
        DO 20 J = I,N
         A(J,I) = A(J,I) * FAC * H1(J)
         B(J,I) = B(J,I) * FAC * H1(J)
   20   CONTINUE
   30  CONTINUE
C
C  GET EIGENVALUES OF B AND SORT EIGENSTUFF
C  AFTER COMPLETION, EIGENVECTORS ARE STORED IN B
C
       JEV = 1
       CALL TRIDIA(NB,N,EW,H2,B,JEV) 
       CALL IQLDIA(NB,N,EW,H2,B,JEV,IER) 
       IF (IER .NE. 0) RETURN 
       JORD = 1 
       CALL SORTVC(NB,N,N,EW,H2,B,JORD,JEV)
C
C  DEFINE NE
C
       NE = 0
       DO 40 I = 1,N
        IF (EW(I) .LT. DELTA) GOTO 50
        NE = NE+1
   40  CONTINUE
   50  IF (NE .EQ. 0) RETURN
C
C  TRANSFORM A. FIRST FILL MATRIX
C
       DO 70 I = 1,N
        DO 60 J = I+1,N
         A(I,J) = A(J,I)
   60   CONTINUE
   70  CONTINUE
C
C  TRANSFORM A = (B' * A) 
C
       DO 110 I = 1,N
        DO 80 K = 1,N
         H2(K) = A(K,I)
   80   CONTINUE
        DO 100 J = 1,NE
         SUM = 0.0D0
         DO 90 K = 1,N
          SUM = SUM + B(K,J) * H2(K)
   90    CONTINUE
         A(J,I) = SUM
  100   CONTINUE
  110  CONTINUE
C
C  TRANSFORM A = (A * B) (DEFINE ONLY THE LOWER TRIANGLE)
C
       DO 150 I = 1,NE
        FAC = 1.0D0
        DO 120 K = 1,N
         H2(K) = A(I,K)
  120   CONTINUE
        DO 140 J = 1,NE
         SUM = 0.0D0
         DO 130 K = 1,N
          SUM = SUM + H2(K) * B(K,J)
  130    CONTINUE
         A(I,J) = SUM * FAC
  140   CONTINUE
  150  CONTINUE
C
C  TRANSFORM EW = EW**(-1/2)
C
       DO 160 I = 1,NE
        EW(I) = 1.0D0 / SQRT(EW(I))
  160  CONTINUE
C
C  TRANSFORM A = (EW * A * EW)
C
       DO 180 I = 1,NE
        FAC = EW(I)
        DO 170 J = I,NE
         A(J,I) = A(J,I) * FAC * EW(J)
  170   CONTINUE
  180  CONTINUE
C
C  TRANSFORM B = (H1 * B * EW)
C
       DO 200 I = 1,NE
        FAC = EW(I)
        DO 190 J = 1,N
         B(J,I) = B(J,I) * FAC * H1(J)
  190   CONTINUE
  200  CONTINUE
C
C  DIAGONALIZE A AND SORT EIGENSTUFF
C  AFTER COMPLETION, EIGENVECTORS ARE STORED IN A
C  IF ONLY EIGENVALUES ARE NEEDED, RETURN
C
       JEV = 0
       IF (IEV .NE. 0) JEV = 1
       CALL TRIDIA(NA,NE,EW,H2,A,JEV) 
       CALL IQLDIA(NA,NE,EW,H2,A,JEV,IER) 
       IF (IER .NE. 0) IER = IER+N 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC(NA,NE,NE,EW,H2,A,IORD,JEV)
       IF (JEV .EQ. 0) RETURN
C
C  IN ORDER TO MINIMIZE STRIDE IN THE INNER LOOP FOR THE FOLLOWING
C  THE BACK TRANSFORMATION, SWITCH ROWS AND COLUMNS IN B
C
       DO 220 I = 1,NE
        DO 210 J = I+1,N
         FAC = B(J,I)
         B(J,I) = B(I,J)
         B(I,J) = FAC
  210   CONTINUE
  220  CONTINUE
C
C  BACK TRANSFORMATION. FORM (B * A) IN A
C
       DO 260 I = 1,NE
        DO 230 K = 1,NE
         H2(K) = A(K,I)
  230   CONTINUE
        DO 250 J = 1,N
         SUM = 0.0D0
         DO 240 K = 1,NE
          SUM = SUM + B(K,J) * H2(K)
  240    CONTINUE
         A(J,I) = SUM
  250   CONTINUE
  260  CONTINUE
       RETURN 
      END
