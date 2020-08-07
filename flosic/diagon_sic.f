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
      SUBROUTINE DIAGGE_FO(NA,N,A,B,EVAL,SC,IEV)
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
       CALL EWEVGE_FO(NA,NA,N,A,B,EVAL,SC,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGGE: FAILED IN EWEVGE, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
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
      SUBROUTINE DIAGSP_FO(NA,N,A,EVAL,SC,IEV)
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
       CALL EWEVSP_FO(NA,N,A,EVAL,SC,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGSP: FAILED IN EWEVGE, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
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
      SUBROUTINE DIAGSVD_FO(NA,N,NE,A,B,EVAL,SC1,SC2,DELTA,IEV)
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
       CALL EWEVSVD_FO(NA,NA,N,NE,A,B,EVAL,SC1,SC2,
     &              DELTA,IEV,IORD,IER) 
       IF (IER .NE. 0) THEN
        PRINT *,'DIAGSVD: FAILED IN EWEVSVD, ERROR CODE: ',IER
        PRINT *,'MATRIX SIZE: ',NA,' NDIM: ',N
        CALL STOPIT
       END IF
       RETURN
      END
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
      SUBROUTINE EWEVGE_FO(NA,NB,N,A,B,EW,H,IEV,IORD,IER) 
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
       CALL CHOLES_FO(NB,N,B,IER) 
       IF (IER .NE. 0) RETURN 
       CALL MATRAF_FO(NA,NB,N,H,A,B) 
       CALL TRIDIA_FO(NA,N,EW,H,A,JEV) 
       CALL IQLDIA_FO(NA,N,EW,H,A,JEV,IER) 
       IF (IER .NE. 0) IER = IER+N 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC_FO(NA,N,N,EW,H,A,IORD,JEV)
       IF (JEV .EQ. 0) RETURN
       CALL BACKTR_FO(NB,NA,NA,N,N,H,B,A,A) 
       RETURN 
      END
C
C *************************************************************************
C
C     SUBROUTINE EWEVSP 
C     ================= 
C 
C *************************************************************************
C 
C  Ewevsp calculates eigenvalues and eigenvectors of the real special 
C  symmetric eigenvalue problem. 
C  In the explanation of the method, primes denote an exchange of rows 
C  and columns.
C 
C  Method:  *  A*C = E*C 
C           *  Solve INV(R')*A*INV(R)*Y = E*Y  (Householder + IQL) 
C           *  Sorting of eigenvalues and eigenvectors 
C 
C     Parameters: 
C 
C       NA      (I) :  Dimension of A 
C       N       (I) :  Dimension of Problem  
C       A       (I) :  Matrix A (lower triangle)
C               (O) :  Eigenvector matrix  
C       EW      (O) :  Eigenvalues 
C       H       (-) :  Auxiliary vector 
C       IEV     (I) :  0: No eigenvectors  
C       IORD    (I) :  1: Descending order of eigenvalues 
C                     -1: Ascending order of eigenvalues
C                      otherwise: no sorting 
C       IER     (O) :  Error indication  
C                      0: No error 
C                     -1: N is larger than NA 
C                      K: Convergence failure for eigenvalue
C                         K, (K-1) eigenvalues are correct
C 
C *************************************************************************
C
      SUBROUTINE EWEVSP_FO(NA,N,A,EW,H,IEV,IORD,IER) 
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION  A(NA,N),EW(N),H(N)
       SAVE
C
       IER = 0
       IF (N .LE. 0) RETURN
       IER = -1 
       IF (N .GT. NA) RETURN
       IER = 0 
       JEV = 0 
       IF (IEV .NE. 0) JEV = 1 
       CALL TRIDIA_FO(NA,N,EW,H,A,JEV) 
       CALL IQLDIA_FO(NA,N,EW,H,A,JEV,IER) 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC_FO(NA,N,N,EW,H,A,IORD,JEV)
       RETURN 
      END
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
      SUBROUTINE EWEVSVD_FO(NA,NB,N,NE,A,B,EW,H1,H2,DELTA,IEV,IORD,IER) 
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
       CALL TRIDIA_FO(NB,N,EW,H2,B,JEV) 
       CALL IQLDIA_FO(NB,N,EW,H2,B,JEV,IER) 
       IF (IER .NE. 0) RETURN 
       JORD = 1 
       CALL SORTVC_FO(NB,N,N,EW,H2,B,JORD,JEV)
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
       CALL TRIDIA_FO(NA,NE,EW,H2,A,JEV) 
       CALL IQLDIA_FO(NA,NE,EW,H2,A,JEV,IER) 
       IF (IER .NE. 0) IER = IER+N 
       IF (IER .NE. 0) RETURN 
       CALL SORTVC_FO(NA,NE,NE,EW,H2,A,IORD,JEV)
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
      SUBROUTINE CHOLES_FO(NB,N,B,ICHO)
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
      SUBROUTINE MATRAF_FO(NA,NB,N,H,A,B)
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
      SUBROUTINE TRIDIA_FO(NM,N,D,E,A,IEV)
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
      SUBROUTINE IQLDIA_FO(NM,N,D,E,Z,IEV,IER)
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
       CALL MACHEPS_FO(EPS)
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
      SUBROUTINE BACKTR_FO(NR,NX,NY,N,M,H,R,X,Y)
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
C
C *************************************************************************
C
C     SUBROUTINE SORTVC
C     =================
C
C *************************************************************************
C
C  Sortvc sorts D and (if required) E and the columns of Q.
C
C     Parameters:
C
C       NM      (I) :  Dimension of Q
C       N       (I) :  Dimension of problem (size of one vector in Q)
C       NQ      (I) :  Number of elements in D (or columns in Q)
C       D       (I) :  Vector to sort
C               (O) :  Sorted vector 
C       E       (I) :  Additional Vector to sort
C               (O) :  Sorted additional vector
C       Q       (I) :  Matrix to sort (vectors in columns)
C               (O) :  Sorted matrix (vectors in columns)
C       M       (I) :  1: Descending order in D
C                     -1: Ascending order in D
C                      otherwise: no sorting
C       IEV     (I) :  0: No sorting of Q and E
C                      1: Sorting of Q, no sorting of E
C                      2: Sorting of Q and E
C
C *************************************************************************
C
      SUBROUTINE SORTVC_FO(NM,N,NQ,D,E,Q,M,IEV)
       IMPLICIT REAL*8 (A-H,O-Z)
       LOGICAL    LMIN,LMAX
       DIMENSION  D(NQ),E(NQ),Q(NM,NQ)
       SAVE
C
       IF (NQ .LT. 2) RETURN
       LMAX = (M .EQ.  1)
       LMIN = (M .EQ. -1)
       IF (.NOT. (LMAX .OR. LMIN)) RETURN
       DO 40 KK = 2,NQ
        K = KK - 1
        J = K
        H = D(K)
C
C  FIND EXTREMUM
C
        DO 10 I = KK,NQ
         S = D(I)
         IF (LMIN .AND. (S .GE. H)) GOTO 10
         IF (LMAX .AND. (S .LE. H)) GOTO 10
         J = I
         H = S
   10   CONTINUE
        IF (J .EQ. K) GOTO 40
C
C  SORT D
C
        D(J) = D(K)
        D(K) = H
        IF (IEV .EQ. 0) GOTO 40
C
C  SORT Q
C
        DO 20 I = 1,N
         H = Q(I,K)
         Q(I,K) = Q(I,J)
         Q(I,J) = H
   20   CONTINUE
        IF (IEV .LT. 2) GOTO 40
C
C  SORT E
C
        H    = E(K)
        E(K) = E(J)
        E(J) = H
   40  CONTINUE
       RETURN
      END
C
C ***********************************************************************
C
      SUBROUTINE MACHEPS_FO(EPS)
C
C DETERMINES THE MACHINE ACCURACY, SHOULD WORK AS LONG AS EPS >= 10^(-48)
C
       IMPLICIT NONE
       CHARACTER*50 STR
       REAL*8 EPS,XHF
       SAVE
C
       EPS= 1.0D0
   10  CONTINUE
        EPS=  EPS*0.5D0
        WRITE(STR,1000) EPS+1.0D0
        READ (STR,1000) XHF
        IF (XHF .GT. 1.0D0) GOTO 10
       CONTINUE
       EPS= EPS*2.0D0
       RETURN
 1000  FORMAT(F50.48)
      END
