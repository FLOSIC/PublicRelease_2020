C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE CKWORKER(MODE,TID)
c
c 02/13/97 David Clay Patton
c 04/17/97 converted from PVM to MPI (DCP)
c 04/26/97 fixed call to mpi_iprobe so that mode 1 now works (DCP)
c 07/23/97 revised (DCP)
c
       use mpidat1,only : NPROC, NCALLED, MYCOMM
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:38 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      INCLUDE 'mpif.h'
      SAVE
      INTEGER TID,NPROC,MODE,NCALLED,TAG,I
      INTEGER IRECVSTAT(MPI_STATUS_SIZE),IFLAG,IERR
c
      
      integer profiler(2) / 0, 0 /
      save profiler

      call TAU_PROFILE_TIMER(profiler, '                                &
     &CKWORKER [{ckworker.f} {4,7}-{72,9}]')
      
      call TAU_PROFILE_START(profiler)
      TID=0
      IF (NCALLED.EQ.0) then
      call TAU_PROFILE_STOP(profiler)
                        RETURN
         endif
c
c check to see if any workers are finished playing
c
      TAG=1
      IF (MODE.EQ.1) THEN
       CALL MPI_IPROBE(MPI_ANY_SOURCE,TAG,MYCOMM,
     &                 IFLAG,IRECVSTAT,IERR)
c
c if so, get their tid
c
       IF (ABS(IFLAG).NE.0) THEN
        CALL MPI_RECV(TID,1,MPI_INTEGER,IRECVSTAT(MPI_SOURCE),TAG,
     &                MYCOMM,IRECVSTAT,IERR)
        NCALLED=NCALLED-1
        CALL FREETID(TID)
       END IF
      END IF
c
c if all workers are out playing then get result back from first to finish
c otherwise return immediately
c
      IF (MODE.EQ.2) THEN
       IF (NCALLED.EQ.NPROC) THEN
        CALL MPI_RECV(TID,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,
     &                MYCOMM,IRECVSTAT,IERR)
        NCALLED=NCALLED-1
        CALL FREETID(TID)
       END IF
      END IF
c
c let all workers finish
c
      IF (MODE.EQ.3) THEN
       IF (NCALLED.NE.0) THEN
        DO I=1,NCALLED
         CALL MPI_RECV(TID,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,
     &                 MYCOMM,IRECVSTAT,IERR)
         CALL FREETID(TID)
        END DO
       END IF
       NCALLED=0
      END IF
c
c wait for the first worker to finish
c
      IF (MODE.EQ.4) THEN
       IF (NCALLED.EQ.0) then
      call TAU_PROFILE_STOP(profiler)
                         RETURN
         endif
       CALL MPI_RECV(TID,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,
     &               MYCOMM,IRECVSTAT,IERR)
       NCALLED=NCALLED-1
       CALL FREETID(TID)
      END IF
      call TAU_PROFILE_STOP(profiler)
      RETURN
      call TAU_PROFILE_STOP(profiler)
      END
