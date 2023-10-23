C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE GETTID(TID) 
c
c get id for free node, stop if all nodes are busy
c 02/13/97 David Clay Patton
c changed 10/98 by DVP
c
#ifdef MPI
       use mpidat1,only : INUSE, NPROC
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      SAVE
      INTEGER :: TID,I
      TID=0
      DO I=1,NPROC
       IF (INUSE(I).EQ.0) THEN
        TID=I
        GOTO 10
       END IF
      END DO
c
   10 IF (TID .EQ. 0) THEN
       write(6,*)'gettid: called while all nodes are busy'
       CALL STOPIT
      ELSE
       INUSE(TID)=1
      ENDIF
#endif
      RETURN
      END
