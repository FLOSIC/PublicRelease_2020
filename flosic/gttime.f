C UTEP Electronic Structure Lab (2020)
C
C ************************************************************
C
      SUBROUTINE GTTIME(TIME)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL FIRST
!DIR$ NAME (MYCLOCK="myclock")
      SAVE
      DATA FIRST/.TRUE./
      CALL MYCLOCK(XTIME)
      TIMEX=0.01D0*XTIME
      IF(FIRST)THEN
      FIRST=.FALSE.
      TIME0=TIMEX
      ENDIF
      TIME=TIMEX-TIME0
      RETURN
      END
