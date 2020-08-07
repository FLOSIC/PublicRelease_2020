C UTEP Electronic Structure Lab (2020)
C
C ***********************************************************************
C
      SUBROUTINE MACHEPS(EPS)
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
