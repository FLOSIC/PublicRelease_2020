C UTEP Electronic Structure Lab (2020)
        
        FUNCTION FCTRL(NN)
        IMPLICIT REAL*8 (A-H,O-Z)
           IF(NN.GE.0)THEN
           FCTRL=1.0D0
           ELSE
           FCTRL=0.0D0
           END IF
           DO I=1,NN
           FCTRL=FCTRL*I
           END DO
        RETURN
        END
