C UTEP Electronic Structure Lab (2020)
      SUBROUTINE SYMMAT(MATRIX,DMAT,N,OPTION)
C
C     SYMMATRISATION AND ANTISYMMETRIZATION OF A LOWER OR UPPER 
C     TRIANGLE MATRIX.
C
C     BY: ULISES REVELES, JUNE 2013.
C
C     ------------------------------------------------------------------
C
C     LOCAL DIMENSIONS:
C
C     DMAT: Dimension of MATRIX.
C
C     LOCAL VARIABLES:
C
C     MATRIX: Matrix to be symmetrised.
C     N     : Number of columns or rows of matrix.
C     OPTION: a) LOWUP = Symmetrisation of lower triangular matrix.
C             b) UPLOW = Symmetrisation of upper triangular matrix.
C             c) -LOWUP = Antisymmetrisation of lower triangular matrix.
C             d) -UPLOW = Antisymmetrisation of upper triangular matrix.
C
C     ------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER DMAT
      REAL*8 MATRIX(DMAT,DMAT)
C
      CHARACTER OPTION*(*)
      INTEGER I,N
C
C     ------------------------------------------------------------------
C
      IF (OPTION.EQ.'LOWUP') THEN
        DO I=1,N-1
          CALL DCOPY(N-I,MATRIX(I+1,I),1,MATRIX(I,I+1),DMAT)
        END DO
      ELSE IF (OPTION.EQ.'UPLOW') THEN
        DO I=1,N-1
          CALL DCOPY(N-I,MATRIX(I,I+1),DMAT,MATRIX(I+1,I),1)
        END DO
      ELSE IF (OPTION.EQ.'-LOWUP') THEN
        DO I=1,N-1
          CALL DCOPY(N-I,MATRIX(I+1,I),1,MATRIX(I,I+1),DMAT)
        END DO
        DO I=2,N
          CALL DSCAL(I-1,-1.0,MATRIX(1,I),1)
        END DO
      ELSE IF (OPTION.EQ.'-UPLOW') THEN
        DO I=1,N-1
          CALL DCOPY(N-I,MATRIX(I,I+1),DMAT,MATRIX(I+1,I),1)
        END DO
        DO I=1,N-1
          CALL DSCAL(N-I,-1.0,MATRIX(I+1,I),1)
        END DO
      ELSE
        WRITE(*,*)'ERROR IN SYMMAT: UNKNOWN OPTION FOR SYMMAT'
      END IF
C
C     ------------------------------------------------------------------
C
      END
