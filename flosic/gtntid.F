C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE GTNTID(TID) 
c
c get tid for free node or next available node     
c 12/02/99 Mark R. Pederson 
c
#ifdef MPI
       use mpidat1,only : INUSE, NPROC, MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:52 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IERR, IFLAG
       REAL*8 :: TAG 
      INCLUDE 'mpif.h'
      SAVE
      INTEGER :: TID,I,IRECVSTAT(MPI_STATUS_SIZE)
      TID=0
      DO I=1,NPROC
       IF (INUSE(I).EQ.0) THEN
        TID=I
        INUSE(TID)=1
        GOTO 10
       END IF
      END DO
   10 IF (TID .EQ. 0) THEN
C WAIT FOR NEXT PROCESSOR TO FINISH
       TAG=1
       CALL MPI_IPROBE(MPI_ANY_SOURCE,TAG,MPI_COMM_WORLD,
     &                 IFLAG,IRECVSTAT,IERR)
        CALL MPI_RECV(TID,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
      ENDIF
#endif
      RETURN
      END
