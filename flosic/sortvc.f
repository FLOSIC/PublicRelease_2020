C UTEP Electronic Structure Lab (2020)
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
      SUBROUTINE SORTVC(NM,N,NQ,D,E,Q,M,IEV)
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
