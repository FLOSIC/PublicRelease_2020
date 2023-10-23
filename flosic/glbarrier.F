C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE GLBARRIER
c     
c this routine establishes a global synchronization barrier
c
#ifdef MPI
       use mpidat1,only : NPROC, IRANK, MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INCLUDE 'mpif.h'
       SAVE
       INTEGER IPROC,JOB,TAG,IERR
C
       JOB=100
       TAG=0
       IF (IRANK .EQ. 0) THEN
        DO IPROC= 1,NPROC
         CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
        END DO
       END IF
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
       RETURN
      END
