C UTEP Electronic Structure Lab (2020)
      SUBROUTINE STOPIT
C
C     STOP THE PROGRAM IN A CLEAN WAY
C
C     UPDATE BY ULISES REVELES, JULY 2013.
C
C     ------------------------------------------------------------------
C
C
      use mpidat1,only : MYCOMM
      LOGICAL EXIST
#ifdef MPI
       INCLUDE 'mpif.h'
#endif
C     ------------------------------------------------------------------
C
      OPEN (UNIT=11, FILE='ERROR_NRLMOL')
      WRITE(11,'(A)') 'ERROR ...'
      CLOSE(11)
C LB: REMOVE LOCK FILE
C      INQUIRE(FILE='RUNNING',EXIST=EXIST)
C      IF(EXIST)THEN
C        WRITE(6,*)'REMOVING RUNNING FILE'
C        CALL SYSTEM('rm RUNNING')
C      ELSE
C        WRITE(6,*)'SOMETHING WRONG HERE, RUNNING FILE NOT FOUND!'
C      ENDIF
C LB
C
#ifdef MPI
      PRINT '(A)','Aborting MPI ...'
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR)
#endif
C
      STOP
C
C     ------------------------------------------------------------------
C
      END
