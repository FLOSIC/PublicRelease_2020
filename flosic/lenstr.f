C UTEP Electronic Structure Lab (2020)
C
C MINIMUM FUNCTIONALITY PREPROCESSOR FOR CONDITIONAL COMPILATION
C READS FROM STDIN AND WRITES TO STDOUT
C DIRK POREZAG, MAY 1998
C 
C CALL: condcomp -Ddef1 -Ddef2 ... 
C
C *****************************************************************
C
C SUPPLEMENTARY FUNCTION TO FIND NONBLANK LENGTH OF STRING
C
      INTEGER FUNCTION LENSTR(S)
       CHARACTER*300 S
       LENSTR= 0
       DO I= 300, 1, -1
        IF (S(I:I) .NE. ' ') THEN
         LENSTR= I
         RETURN
        END IF
       END DO
       RETURN
      END
