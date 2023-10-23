C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PAMMESH(MODE,BOX,ERRBOX,NMULT)
c
c 12/98 Dirk V. Porezag 
c
c mode 1: send box data to workers 
c mode 2: get mesh points from workers
c
#ifdef MPI
      use mesh1,only : wmsh,rmsh,nmsh
      use mpidat1,only : NPROC, NCALLED, IRANK, MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IPROC, NTID
       REAL*8 :: BOX , ERRBOX
      INCLUDE 'mpif.h'
      SAVE
      DIMENSION BOX(2,3)
      INTEGER MODE,NMULT
      INTEGER TID,JOB,TAG,IERR,MPTS,IRECVSTAT(MPI_STATUS_SIZE)
c
c produce fatal error if I am a worker
c
      IF (IRANK.NE.0) THEN
       write(6,*)'FATAL: PAMMESH CALLED BY WORKER'
       CALL STOPIT
      END IF
c
c mode= 1: send data for interstitial mesh creation to worker
c
      IF (MODE.EQ.1) THEN
       NCALLED=NCALLED+1
       CALL GETTID(TID)
       JOB= 14
       TAG= 0
       CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 1401
       CALL MPI_SSEND(BOX(1,1),6,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 1402
       CALL MPI_SSEND(ERRBOX,1,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 1403
       CALL MPI_SSEND(NMULT,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
      END IF
c
c mode= 2: receive mesh points from workers
c
      IF (MODE.EQ.2) THEN
       CALL CKWORKER(3,NTID)
       DO 100 IPROC=1,NPROC
        JOB= 15
        TAG= 0
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
        TAG= 1501
        CALL MPI_RECV(MPTS,1,MPI_INTEGER,IPROC,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
        IF (MPTS .LE. 0) GOTO 100
        IF (NMSH+MPTS .GT. MAX_PTS) THEN
         write(6,*)'PAMMESH: MAX_PTS MUST BE AT LEAST: ',NMSH+MPTS
         CALL STOPIT
        END IF
        TAG= 1502
        CALL MPI_RECV(RMSH(1,NMSH+1),3*MPTS,MPI_DOUBLE_PRECISION,
     &                IPROC,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        TAG= 1503
        CALL MPI_RECV(WMSH(NMSH+1),MPTS,MPI_DOUBLE_PRECISION,
     &                IPROC,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        NMSH=NMSH+MPTS
  100  CONTINUE
      END IF
#endif
      RETURN
      END
