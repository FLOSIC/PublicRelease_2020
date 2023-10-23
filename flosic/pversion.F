C UTEP Electronic Structure Lab (2020)
c
c **********************************************************************
c
      SUBROUTINE PVERSION()
c
c print statement if serial or parallel version is used
c
#ifdef MPI
c
       use mpidat1,only : NPROC
       PRINT '(A)','MPI VERSION'
C
#ifdef MPI_3
C
       PRINT '(A)','USING MPI-3 SHARED MEMORY'
C
#endif
C
#ifdef GROUP
C
       PRINT '(A)','CALCULATION OF SIC MATRIX IN GROUPS'
       IF(NPROC < 2) THEN
        PRINT '(A)','GROUP CALCULATION REQUIRES 3 OR MORE PROCS'
        !CALL STOPIT
       END IF
C
#endif
c
#ifdef SCALAPACK
c
       PRINT '(A)','USING SCALAPACK FOR DIAGONALIZATION'
c
#endif
c
#else
c
       PRINT '(A)','SERIAL VERSION'
c
#endif
c
       RETURN
      END
