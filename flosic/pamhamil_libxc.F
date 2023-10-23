C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PAMHAMIL_LIBXC(MODE,IFNCT,JFNCT,AI,AJ,HSUB,MYIP)
c
c Subroutine for parallel hamiltonian called by manager
c 20/03/98 Jens Kortus   kortus@theo.physik.tu-freiberg.de
c modifications by DV Porezag 08/99
c
c Mode 1: looking for a free node and sending all data
c Mode 2: receiving data and doing a call do gethold
c
#ifdef MPI
      use mpidat1,only : NCALLED, IRANK, MYCOMM
      INCLUDE 'PARAMA2'
      INCLUDE 'mpif.h'
      SAVE
      INTEGER MODE,IFNCT,JFNCT,MYIP,MXXX
      REAL*8 :: AI,AJ,HSUB
      DIMENSION AI(3),AJ(3)
      DIMENSION HSUB(MAXUNSYM,MAXUNSYM)
      INTEGER TID,JOB,TAG,IERR,ITRANS(3),IRECVSTAT(MPI_STATUS_SIZE)
c
c produce fatal error if I am a worker
c
      IF (IRANK.NE.0) THEN
        write(6,*)'FATAL: PAMHAMIL CALLED BY WORKER'
        CALL STOPIT
      END IF         
c
c send data for gethold execution
c
      IF (MODE.EQ.1) THEN
       NCALLED=NCALLED+1
       CALL GETTID(TID)
!       JOB= 5
       JOB= 28
       TAG= 0
       CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
       ITRANS(1)= IFNCT
       ITRANS(2)= JFNCT
       ITRANS(3)= MYIP 
       TAG= 501
       CALL MPI_SSEND(ITRANS(1),3,MPI_INTEGER,
     &                TID,TAG,MPI_COMM_WORLD,IERR)   
       TAG= 502
       CALL MPI_SSEND(AI(1),3,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)   
       TAG= 503
       CALL MPI_SSEND(AJ(1),3,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)   
      END IF  
c
c receive gethold data from worker
C
      IF (MODE.EQ.2) THEN
       TAG=1
       CALL MPI_RECV(TID,1,MPI_INTEGER,MPI_ANY_SOURCE,TAG,
     &               MPI_COMM_WORLD,IRECVSTAT,IERR)
!       JOB= 6
       JOB= 29
       TAG= 0
       CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
       MXXX=MAXUNSYM*MAXUNSYM
       TAG= 601
       CALL MPI_RECV(HSUB(1,1),MXXX,MPI_DOUBLE_PRECISION,
     &               TID,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       TAG= 602
       CALL MPI_RECV(MYIP,1,MPI_INTEGER,
     &               TID,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
 
       NCALLED=NCALLED-1
       CALL FREETID(TID)
      END IF
#endif  
      RETURN
      END
