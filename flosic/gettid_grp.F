C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE GETTID_GRP(TID) 
c
c get id for free node, stop if all nodes are busy
c 02/13/97 David Clay Patton
c changed 10/98 by DVP
! Group version by Luis Basurto 10/19/2017
c
#ifdef GROUP
       use debug1
       use mpidat1,only : INUSE_GRP, SHM_SIZE
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      SAVE
      INTEGER :: TID,I
      TID=0
!      CALL TRACER('GETTID_GRP:SHM_SIZE',SHM_SIZE)
!      CALL TRACER('GETTID_GRP:SIZE INUSE_GRP',SHM_SIZE)
      DO I=1,SHM_SIZE
!       CALL TRACER('I',I)
!       CALL TRACER('INUSE_GRP',INUSE_GRP(I))
       IF (INUSE_GRP(I).EQ.0) THEN
        TID=I
        GOTO 10
       END IF
      END DO
c
   10 IF (TID .EQ. 0) THEN
!       CALL TRACER('gettid_grp: called while all nodes are busy')
       CALL STOPIT
      ELSE
       INUSE_GRP(TID)=1
      ENDIF
#endif
      RETURN
      END
