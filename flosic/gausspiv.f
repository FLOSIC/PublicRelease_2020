C UTEP Electronic Structure Lab (2020)
c       
c ***********************************************************
c
        SUBROUTINE GAUSSPIV(NA,NB,NDIM,NCOL,A,B,SCR,ISCR,IER)
c
c ***********************************************************
c  Copyright by Dirk V. Porezag (porezag@dave.nrl.navy.mil) *
c ***********************************************************/
c 
c Gausspiv performs a Gaussian elimination to solve the matrix 
c equation a*x = b. Complete line and row pivoting is applied.
c The result x is stored in b. Gausspiv modifies a and b.
c
c Dimension of the matrices: a(na rows, na columns) 
c                            b(na rows, nb columns)
c
c Input:   na,nb   see above matrix dimensions
c          ndim    actual dimension of the problem (must be <= na)
c          ncol    actual number of vectors in b (must be <= nb)
c          a       matrix a
c          b       matrix b
c          scr     scratch vector of size ndim
c          iscr    scratch vector of size ndim
c Output:  b       matrix x with a*x = b
c          ier     error code
c
c Error codes:      0: no error
c                  -1: na, nb, ncol, and ndim are incompatible
c                  -2: should not happen, can only be caused by
c                      a bug or broken compiler.
c                   n: (n > 0) means that the matrix a has a zero
c                      determinant and that the gauss algorithm
c                      failed in the nth line
c
c *************************************************************** */
c
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         PARAMETER(ZERO=1.0D-30)
         DIMENSION A(NA,NA),B(NA,NB),SCR(NDIM),ISCR(NDIM)
         SAVE
c
c special cases and setup
c
         IER= 0
         IF ((NDIM .LT. 1) .OR. (NCOL .LT. 1)) RETURN
         IER= -1
         IF ((NDIM .GT. NA) .OR. (NCOL .GT. NB)) RETURN
         DO I= 1,NDIM
          ISCR(I)= I
         END DO
c
c Gauss elimination: construction of triangular form
c
         DO I= 1,NDIM
c
c look for largest element in submatrix
c
          DMAX= ABS(A(I,I))
          INDJ= I
          INDK= I
          DO J= I,NDIM
           DO K= I,NDIM
            SAV= ABS(A(K,J))
            IF (SAV .GT. DMAX) THEN
             DMAX= SAV
             INDJ= J
             INDK= K
            END IF
           END DO
          END DO
          XREF= ZERO*ABS(A(1,1))
          IF (I .EQ. 1) XREF= 0.0D0
          IF (DMAX .LE. XREF) THEN
           IER= I
           RETURN
          END IF
c
c exchange rows
c
          IF (I .NE. INDK) THEN
           DO J= I,NDIM
            SAV= A(I,J)
            A(I,J)= A(INDK,J)
            A(INDK,J)= SAV
           END DO
           DO J= 1,NCOL
            SAV= B(I,J)
            B(I,J)= B(INDK,J)
            B(INDK,J)= SAV
           END DO
          END IF
c
c exchange columns
c
          IF (I .NE. INDJ) THEN
           DO J= 1,NDIM
            SAV= A(J,I)
            A(J,I)= A(J,INDJ)
            A(J,INDJ)= SAV
           END DO
           ISAV= ISCR(I)
           ISCR(I)= ISCR(INDJ)
           ISCR(INDJ)= ISAV 
          END IF
c
c one step towards triangular form, invert diagonal element
c
          AREC= 1.0D0/A(I,I)
          A(I,I)= AREC
          DO J= I+1,NDIM
           SCR(J)= AREC*A(J,I)
          END DO
          DO K= I+1,NDIM
           FAC= A(I,K)
           DO J= I+1,NDIM
            A(J,K)= A(J,K)-FAC*SCR(J)
           END DO
          END DO
          DO K= 1,NCOL
           FAC= B(I,K)
           DO J= I+1,NDIM
            B(J,K)= B(J,K)-FAC*SCR(J)
           END DO
          END DO
         END DO
c
c Gauss elimination: backward substitution 
c
         DO I= NDIM,1,-1
          AREC= A(I,I)
          DO K= I+1,NDIM
           SCR(K)= A(I,K)
          END DO
          DO J= 1,NCOL
           DO K= I+1,NDIM
            B(I,J)= B(I,J)-SCR(K)*B(K,J)
           END DO
           B(I,J)= B(I,J)*AREC
          END DO
         END DO
c
c restore old ordering
c 
         IER= -2
         DO I= 1,NDIM
          DO INDI= I,NDIM
           IF (ISCR(INDI) .EQ. I) GOTO 10
          END DO
          RETURN
   10     ISCR(INDI)= ISCR(I)
          DO J= 1,NCOL
           SAV= B(I,J)
           B(I,J)= B(INDI,J)
           B(INDI,J)= SAV
          END DO
         END DO
         IER= 0
         RETURN
        END
