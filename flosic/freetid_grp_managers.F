C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE FREETID_GRP_MANAGER(TID)
c
c 02/13/97 David Clay Patton
c changed 10/98 by DVP
c
#ifdef MPI
       use mpidat1,only : INUSE_MANAGER, MANAGER_SIZE
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      SAVE
      INTEGER :: TID
      IF ((TID .GE. 1) .AND. (TID.LE.MANAGER_SIZE)) THEN
       INUSE_MANAGER(TID)=0
      ELSE
       write(6,*)'freetid: called with invalid tid number: ',TID
       CALL STOPIT
      ENDIF
#endif
      RETURN
      END
