C UTEP Electronic Structure Lab (2020)
c
c ************************************************************
c
c logcgr version dirk porezag august 1994
c
       SUBROUTINE LOGCGR(IWARN,IMODE,ICHRA,ICHRB,ILINE,ISTEP)
        CHARACTER ICHRA,ICHRB
        CHARACTER*23 LINE
        SAVE
        OPEN(45,FILE='CGRLOG',FORM='formatted',STATUS='unknown')
        REWIND(45)
   10    READ(45,'(a23)',END=20) LINE
         GOTO 10
   20   CONTINUE
        BACKSPACE(45)
        IF (IWARN .EQ. 1) THEN
         WRITE(45,*) '1d-minimization resulted in converged gradients'
         WRITE(45,*) 'but did not lead to the best function value'
         WRITE(45,*) 'check accuracy of function value and gradients'
         GOTO 30
        END IF
        IF (IWARN .EQ. 2) THEN
         WRITE(45,*) 'maximum between ',ICHRA,' and ',ICHRB,
     &               ' in line search detected'
         WRITE(45,*) 'there is either an error in the ',
     &               'function value / gradients'
         WRITE(45,*) 'or the stepwidth is too big'
        END IF
        IF (IMODE .EQ. 0) THEN
         LINE= 'error: imode=0         '
        ELSE IF (IMODE.EQ.1) THEN
         LINE= 'expanding interval     '
        ELSE IF (IMODE.EQ.2) THEN
         LINE= 'quadratic interpolation'
        ELSE IF (IMODE.EQ.3) THEN
         LINE= 'linear interpolation   '
        ELSE IF (IMODE.EQ.4) THEN
         LINE= 'bisection              '
        END IF
        WRITE(45,100) ISTEP,ILINE,LINE
  100   FORMAT('istep= ',I2,' iline= ',I2,'  ',A23)
   30   CLOSE(45)
        RETURN
       END
