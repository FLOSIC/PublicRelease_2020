C UTEP Electronic Structure Lab (2020)
C
       FUNCTION DET(NDMREP,REP)
       IMPLICIT REAL*8 (A-H,O-Z)
       SAVE
       DIMENSION REP(3,3)
       IF(NDMREP.EQ.1)THEN
       DET=REP(1,1)
       RETURN
       END IF
       IF(NDMREP.EQ.2)THEN
       DET=REP(1,1)*REP(2,2)-REP(1,2)*REP(2,1)
       RETURN
       END IF
       DET=REP(1,1)*(REP(2,2)*REP(3,3)-REP(2,3)*REP(3,2))-
     &     REP(1,2)*(REP(2,1)*REP(3,3)-REP(2,3)*REP(3,1))+
     &     REP(1,3)*(REP(2,1)*REP(3,2)-REP(2,2)*REP(3,1))
       RETURN
       END
