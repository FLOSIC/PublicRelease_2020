C UTEP Electronic Structure Lab (2020)
      INTEGER FUNCTION STREXT(STRING)
C
C     GET STRING EXTENSION
C
C     BY ULISES REVELES, JUNE 2013
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      CHARACTER STRING*(*)
C
      INTEGER I
C
C     ------------------------------------------------------------------
C
      STREXT = 1
C
      DO I=LEN(STRING),1,-1
        IF (STRING(I:I).NE.' ') THEN
          STREXT = I
          GO TO 999
        END IF
      END DO
C
  999 CONTINUE
C
C     ------------------------------------------------------------------
      END
