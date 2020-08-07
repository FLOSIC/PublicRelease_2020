C UTEP Electronic Structure Lab (2020)
      REAL*8  FUNCTION FFERMI(E,EFRM,TEMP)
C
C     ------------------------------------------------------------------
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:44 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: E , EFRM, TEMP, EDIFF
      SAVE
C
C     ------------------------------------------------------------------
C
      EDIFF=E-EFRM
C
      IF (ABS(EDIFF) .GT. 50*ABS(TEMP)) THEN
        IF (EDIFF*TEMP .LT. 0.0D0) THEN
          FFERMI=1.0D0
        ELSE
          FFERMI=0.0D0
        END IF
      ELSE
        FFERMI=1.0D0/(EXP(EDIFF/TEMP)+1.0D0)
      END IF
C
      RETURN
C
C     ------------------------------------------------------------------
C
       END
