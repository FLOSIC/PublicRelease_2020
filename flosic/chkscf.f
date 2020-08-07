C UTEP Electronic Structure Lab (2020)
      SUBROUTINE CHKSCF(ERGCNV,MODE_RUN,OLDERG1,OLDERG2,ITSCF)
C
C     CHECK FOR SCF CONVERGENCE
C
C     ORIGINAL VERSION BY MARK R PEDERSON (1990)
C     UPDATE BY ULISES REVELES (2013)
C
C     ------------------------------------------------------------------
C
C     --- GLOBAL VARIABLES ---
C  
      use global_inputs,only : SCFTOL
      use common2,only : ETOTAL, EKINONL 
      use common5,only : CONVERGENCE
C
C     --- LOCAL VARIABLES ---
C
      INTEGER MODE_RUN,ITSCF
      REAL*8 ERGCNV,OLDERG1,OLDERG2
C 
C     ------------------------------------------------------------------
C
C     --- INITIALIZATION ---
C
      CONVERGENCE=.FALSE.
C
      IF (SCFTOL .LT. 0.0D0) THEN
        ERGCNV=EKINONL
      ELSE
        ERGCNV=ETOTAL
      END IF
C
C     --- CHECK SCF CONVERGENCE ---
C
      IF(MODE_RUN.EQ.3) THEN
        IF ((ABS(ERGCNV-OLDERG1) .LT. ABS(SCFTOL)) .AND.
     &     (ABS(ERGCNV-OLDERG2) .LT. ABS(SCFTOL))) CONVERGENCE=.TRUE.
          OLDERG2=OLDERG1
C
      ELSE
        IF (ABS(ERGCNV-OLDERG1) .LT. ABS(SCFTOL))  CONVERGENCE=.TRUE.
C
      ENDIF
C
C     --- SPECIAL CASES ---
C
      IF (ABS(SCFTOL) .GE. 1.0D0) CONVERGENCE=.TRUE.
      IF (ITSCF.EQ.1) CONVERGENCE=.FALSE.
      OLDERG1=ERGCNV
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END 
