! UTEP Electronic Structure Lab (2020)
SUBROUTINE SCALA_CALL(JOB)
  use global_inputs,only : inbas,iiev,iimesh
  use for_diag1
  use hstor1,only : MAXCLUSTERSIZE
  use mpidat1
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  INTEGER,INTENT(IN) :: JOB
  INTEGER :: IPROC,IERR
  INTEGER,PARAMETER :: ROOT=0,TAG=0

  IF(IRANK==ROOT) THEN
    DO IPROC = 1, NPROC
      CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
    END DO
  END IF
  IF(JOB==1) THEN
!
! Perform parallel diagonalization (send parameters first)
!
    CALL MPI_BCAST(MAXCLUSTERSIZE,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(INBAS,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(IIEV,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(IIMESH,1,MPI_LOGICAL,ROOT,MPI_COMM_WORLD,IERR)
! worker nodes allocate local array for eigenvalues
    IF(IRANK.NE.0)THEN
      ALLOCATE(AEVAL(INBAS),STAT=IERR)
      IF(IERR.NE.0)then
        WRITE(6,*)'scala_call:Error allocating eval',irank
        CALL FLUSH(6)
      ENDIF
    ELSE
!      IF(INBAS/=SIZE(AEVAL)) THEN
!        DEALLOCATE(AEVAL)
!        ALLOCATE(AEVAL(INBAS))
!      ENDIF
!      IF(INBAS/=SIZE(AHAM,1))THEN
!        DEALLOCATE(AHAM)
!        ALLOCATE(AHAM(INBAS,INBAS))
!      ENDIF
      write(6,*)'Calling Scalapack'
    ENDIF
! call scalapack routine
    CALL DIAGGES(INBAS,AHAM,AEVAL,IIEV)
    IF(IIEV.eq.2) THEN
      IIEV=1
      CALL DIAGGES(INBAS,AHAM,AEVAL,IIEV)
    ENDIF
! worker nodes deallocate local array for eigenvalues
    IF(IRANK.NE.0)THEN
      DEALLOCATE(AEVAL,STAT=IERR)
      IF(IERR.NE.0)then
        WRITE(6,*)'scala_call:Error deallocating eval',irank
        CALL FLUSH(6)
      ENDIF
    ELSE
      CALL SYSTEM('rm HAMTOT OVLTOT')
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE SCALA_CALL
