C UTEP Electronic Structure Lab (2020)
c
c ********************************************************************
c
c determines where to cut a box between two atoms of type ifnct 
c and jfnct and distance d
c
      FUNCTION SPLITHERE(IFNCT,JFNCT,D)
       IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER(EPS=1.0D-15)
       SAVE
c      
       IF (D .LE. 0.0D0) THEN
        SPLITHERE= 0.0D0
        RETURN
       END IF
c 
       R1= 0.0D0
       R2= D
       DO I=1,100
        R= 0.5D0*(R1+R2)
        IF(R.LT.EPS) THEN
         R=0.5D0*D
         GOTO 10
        ENDIF
        CALL RPFIT(IFNCT,R,1.0D0/R,RHOF,P1)
        IF(D-R.LT.EPS) THEN
         R=0.5D0*D
         GOTO 10
        ENDIF
        CALL RPFIT(JFNCT,D-R,1.0D0/(D-R),RHOF,P2)
        IF (ABS(P2-P1) .LT. EPS*ABS(P2+P1)) GOTO 10
        IF (P1 .GT. P2) THEN
         R2= R
        ELSE
         R1= R
        ENDIF
       END DO
   10  SPLITHERE= R
c
c set reasonable limits
c
       SPLITHERE= MIN(SPLITHERE,D*0.8D0)
       SPLITHERE= MAX(SPLITHERE,D*0.2D0)
       RETURN
      END
