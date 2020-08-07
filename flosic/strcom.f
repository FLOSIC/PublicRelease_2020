C UTEP Electronic Structure Lab (2020)
      CHARACTER*(*) FUNCTION STRCOM(STRING)
C
C     ELIMINATES SPACE CHARACTERS AROUND ASSIGNMENT AND COMBINATION
C     CHARACTERS IN AN INPUT STRING.
C
C     BY ULISES REVELES, JULY 2013.
C
C     ------------------------------------------------------------------
C
      CHARACTER STRING*(*)
C
      CHARACTER*1 ASSIGNS,COMBINE,LOLIMIT,UPLIMIT
C
      LOGICAL SYMBOL
      INTEGER I,IM,IP,Z
C
C     ------------------------------------------------------------------
C
C     --- INITIALIZATION ---
C
      ASSIGNS = '='
      COMBINE = '-'
      LOLIMIT = '<'
      UPLIMIT = '>'
C
      IF ((INDEX(STRING,ASSIGNS).EQ.0).AND.
     $    (INDEX(STRING,COMBINE).EQ.0).AND.
     $    (INDEX(STRING,LOLIMIT).EQ.0).AND.
     $    (INDEX(STRING,UPLIMIT).EQ.0)) THEN
        STRCOM = STRING
        RETURN
      ELSE
        Z = 1
        STRCOM = ' '
      END IF
C
      DO 10 I=1,LEN(STRING)
        SYMBOL = .FALSE.
        IM = MAX(1,I-1)
        IP = MIN(LEN(STRING),I+1)
        IF (STRING(IM:IM).EQ.ASSIGNS) SYMBOL = .TRUE.
        IF (STRING(IP:IP).EQ.ASSIGNS) SYMBOL = .TRUE.
        IF (STRING(IM:IM).EQ.COMBINE) SYMBOL = .TRUE.
        IF (STRING(IP:IP).EQ.COMBINE) SYMBOL = .TRUE.
        IF (STRING(IM:IM).EQ.LOLIMIT) SYMBOL = .TRUE.
        IF (STRING(IP:IP).EQ.LOLIMIT) SYMBOL = .TRUE.
        IF (STRING(IM:IM).EQ.UPLIMIT) SYMBOL = .TRUE.
        IF (STRING(IP:IP).EQ.UPLIMIT) SYMBOL = .TRUE.
        IF (SYMBOL.AND.STRING(I:I).EQ.' ') THEN
          GO TO 10
        ELSE
          STRCOM(Z:Z) = STRING(I:I)
          Z = Z + 1
        END IF
   10 CONTINUE
C
C     ------------------------------------------------------------------
C
      END
