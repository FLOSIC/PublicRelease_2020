C UTEP Electronic Structure Lab (2020)
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE INVERSE(A,B,M)
C     =============================================================
      IMPLICIT REAL*8  (A-H,O-Z)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PARAMETER (IMATSZ=40)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DIMENSION A(IMATSZ,IMATSZ),B(IMATSZ,IMATSZ)
      DIMENSION TD(IMATSZ),AD(IMATSZ),BD(IMATSZ)
      SAVE
C
C SUBROUTINE TO PREFORM GAUSSIAN ELIMINATION
C            NO ZEROS ALONG THE DIAGONAL
C
      N=M
      IF(N.GT.IMATSZ)THEN
       PRINT *,'INVERT: MATRIX A TOO LARGE'
       CALL STOPIT
      END IF
C
      DO 14 I=1,N
      ATMP=A(I,I)
      IF(ABS(ATMP) .LT. 1.0D-08)THEN
        WRITE(66,'(2A,I4)') 'INVERT: MATRIX HAS ZERO DIAGONAL ',
     &                      'ELEMENT IN ROW: ',I
        CALL STOPIT
      ENDIF
  14  CONTINUE
C
      IF(N.EQ.1) GO TO 605
C
      DO 23 I=1,N
C
      DO 35 J=1,N
 35      TD(J)=A(J,I)/A(I,I)
C
C     TD(I)=(0.0E+00,0.0E+00)
      TD(I)=0.0D0
C
      DO 71 K=1,N
         BD(K)=B(I,K)
 71      AD(K)=A(I,K)
C
      DO 601 K=1,N
      DO 601 J=1,N
         B(J,K)=B(J,K)-(TD(J)*BD(K))
 601     A(J,K)=A(J,K)-(TD(J)*AD(K))
C
 23   CONTINUE
C
      DO 603 I=1,N
      DO 603 J=1,N
 603     B(J,I)=B(J,I)/A(J,J)
C
      RETURN
C
 605  B(1,1)=1.0D0/A(1,1)
      RETURN
      END
