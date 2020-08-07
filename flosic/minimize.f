C UTEP Electronic Structure Lab (2020)
C
C *********************************************************************
C
       SUBROUTINE MINIMIZE(IMN,XINIT,X,F)
C
C FINDS A MINIMUM OF THE FUNCTION F(X) WITH X > 0.
C XINIT IS THE INITIAL VALUE (MUST BE > 0)
C
        IMPLICIT REAL*8 (A-H,O-Z)
        LOGICAL BRAK,GOUP,GODN
        DIMENSION XTB(3),FTB(3)
        SAVE
        DATA FACUP,FACDN /10.0D0,0.1D0/
C
C INITIALIZATION
C
        IF (XINIT .LE. 0.0D0) THEN
         PRINT *,'MINIMIZE: XINIT MUST BE > 0'
         CALL STOPIT
        END IF
        IF (IMN .LE. 0) RETURN
        IF (IMN .EQ. 1) THEN
         BRAK= .FALSE.
         GOUP= .FALSE.
         GODN= .FALSE.
         X= XINIT
         RETURN
        ELSE IF (IMN .EQ. 2) THEN
         XTB(1)= X 
         FTB(1)= F
         X= X*FACUP
        ELSE IF (IMN .EQ. 3) THEN
         IF (F .LT. FTB(1)) THEN
          GOUP= .TRUE.
          XTB(2)= X
          FTB(2)= F
          X= X*FACUP
         ELSE
          GODN= .TRUE.
          XTB(3)= X
          FTB(3)= F
          XTB(2)= XTB(1)
          FTB(2)= FTB(1)
          X= XTB(2)*FACDN
         END IF
        ELSE
C
C WE ARE BEYOND THE THIRD OPTIMIZATION STEP
C DEFINE INDICES. THE TRIPLE XTB DEFINES TWO INTERVALS. TO FIND THE
C MINIMUM, WE HAVE TO CHECK BOTH BY TURNS.
C EVEN IMN -> NEW X WILL BE BETWEEN XTB(1) AND XTB(2) -> IND3=1
C ODD  IMN -> NEW X WILL BE BETWEEN XTB(2) AND XTB(3) -> IND3=3
C
         IF (MOD(IMN,2) .EQ. 1) THEN
          IND1=1
          IND3=3
         ELSE
          IND1=3
          IND3=1
         END IF
C
C IF BRACKETS ARE AVAILABLE, MODIFIED BISECTION
C
         IF (BRAK) THEN
          IF (F .LT. FTB(2)) THEN
           XTB(IND3)= XTB(2)
           FTB(IND3)= FTB(2)
           XTB(2)= X
           FTB(2)= F
          ELSE
           XTB(IND1)= X
           FTB(IND1)= F
          END IF
          X= SQRT(XTB(2)*XTB(IND3))
         END IF
C
C WE ARE LOOKING FOR BRACKETS BY GOING UP
C
         IF (GOUP) THEN
          IF (F .LT. FTB(2))THEN
           XTB(1)= XTB(2)
           FTB(1)= FTB(2)
           XTB(2)= X
           FTB(2)= F
           X= XTB(2)*FACUP
          ELSE
           BRAK= .TRUE.
           GOUP= .FALSE.
           XTB(3)= X
           FTB(3)= F
           X= SQRT(XTB(2)*XTB(IND3))
          END IF
         END IF
C
C WE ARE LOOKING FOR BRACKETS BY GOING DOWN
C
         IF (GODN) THEN
          IF (F .LT. FTB(2)) THEN
           XTB(3)= XTB(2)
           FTB(3)= FTB(2)
           XTB(2)= X
           FTB(2)= F
           X= XTB(2)*FACDN
          ELSE
           BRAK= .TRUE.
           GODN= .FALSE.
           XTB(1)= X
           FTB(1)= F
           X= SQRT(XTB(2)*XTB(IND3))
          END IF
         END IF
        END IF
       END
