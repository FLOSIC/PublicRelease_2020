C UTEP Electronic Structure Lab (2020)
      CHARACTER*(*) FUNCTION STRCOMP(STRING)
C
C     ELIMINATE MULTIPLE SPACES IN A STRING
C
C     BY ULISES REVELES, JUNE 2013
C
C     ******************************************************************
C
      IMPLICIT NONE
C
      CHARACTER STRING*(*)
      INTEGER I,Z
C
C     ------------------------------------------------------------------
C
      STRCOMP = ' '
C
      Z = 1
C
      DO I=1,LEN(STRING)
        IF (Z.EQ.1) THEN
          IF (STRING(I:I).NE.' ') THEN
            STRCOMP(Z:Z) = STRING(I:I)
            Z = Z + 1
          END IF
        ELSE
          IF ((STRING(I:I).NE.' ').OR.(STRCOMP(Z-1:Z-1).NE.' ')) THEN
            STRCOMP(Z:Z) = STRING(I:I)
            Z = Z + 1
          END IF
        END IF
      END DO
C
C     ------------------------------------------------------------------
C
      END
