C UTEP Electronic Structure Lab (2020)
c 
c *******************************************************************
c
      SUBROUTINE PASWAVE(MODE)
C
C  RECEIVING BLOCKHAMILTONIANS FOR DIAGGE
C  JK 09/99
C  MODE 1: RECEIVE DATA AND CALL DIAGGE
C  MODE 2: SEND RESULTS BACK
C
#ifdef MPI
      use common8,only : REP
      use mpidat1,only : IRANK, MPI_COMM_WORLD
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MPI_COMM_WORLD, NDH2, NEIG
       REAL*8 :: DSINGV , SC1, SC2, HAM
      INCLUDE 'mpif.h'
      SAVE       
      INTEGER MODE
      INTEGER TID,JOB,TAG,IERR,ITRANS(3),IRECVSTAT(MPI_STATUS_SIZE)         
      INTEGER IREP,ISVD,NBAS
c
C  ITRANS(1)=IREP
C  ITRANS(2)=ISVD
C  ITRANS(3)=NBAS
c
c produce fatal error if I am manager
c
      IF (IRANK.EQ.0) THEN
       write(6,*)'FATAL: PASWAVE CALLED BY MANAGER'
       CALL STOPIT
      END IF
C
      IF (MODE.EQ.1) THEN
       TAG= 2000
       CALL MPI_RECV(ITRANS(1),3,MPI_INTEGER,0,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
       IREP=ITRANS(1)
       ISVD=ITRANS(2)
       NBAS=ITRANS(3)
       TAG= 2001
       NDH2=NDH*NBAS
       CALL MPI_RECV(HAM(1,1),NDH2,MPI_DOUBLE_PRECISION,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       TAG= 2002
       CALL MPI_RECV(OVER(1,1),NDH2,MPI_DOUBLE_PRECISION,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
C
C call DIAGGE or SVD depending on ISVD, tell mom I am done
        IF(ISVD.EQ.0) THEN
          CALL DIAGSVD(NDH,NBAS,NEIG,HAM,OVER,EVAL,
     &                  SC1,SC2,DSINGV,1)
          NBAS=NEIG 
        ELSE
          CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
        END IF                
       TAG=1
       CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,
     &                 MPI_COMM_WORLD,IERR)                   
C
      END IF
c
c MODE = 2: SEND RESULTS BACK
c
      IF (MODE.EQ.2) THEN
       TAG=2100
       ITRANS(1)=IREP
       ITRANS(2)=NBAS
       CALL MPI_SSEND(ITRANS(1),2,MPI_INTEGER,
     &               0,TAG,MPI_COMM_WORLD,IERR)
       NDH2=NDH*NBAS
       TAG= 2101
       CALL MPI_SSEND(HAM(1,1),NDH2,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IERR)
       TAG= 2102
       CALL MPI_SSEND(OVER(1,1),NDH2,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IERR)
       TAG= 2103
       CALL MPI_SSEND(EVAL(1),NBAS,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IERR)
 
      END IF
#endif
      RETURN
      END     
