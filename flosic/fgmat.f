C UTEP Electronic Structure Lab (2020)
C
C *********************************************************************
C
       SUBROUTINE FGMAT
C WRITTEN BY MARK R PEDERSON (1985)
C
C READS THE GRPMAT FILE AND CHECKS FOR A LARGER GROUP
C
C       use debug1
       use common3,only : RMAT, NGRP
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:44 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IGRP, IIIG, ITIMES, J, JGRP, K, KGRP, MGRP, NGRERR
       REAL*8 :: ERROR , PMAT, TOLER
       SAVE
       LOGICAL BEEN_CALLED,EXIST
       DIMENSION PMAT(3,3)
       DATA BEEN_CALLED,TOLER/.FALSE.,1.0D-8/
C
       IIIG=0
 15    CONTINUE
       IIIG=IIIG+1
       IF(IIIG.GT.2) CALL STOPIT
       IF (BEEN_CALLED) RETURN
       PRINT '(A)','READING AND ANALYZING GROUP MATRICES'
       BEEN_CALLED=.TRUE.
       NGRERR=0
C
C DEFAULT
C
       NGRP=1
       IF (NGRP .GT. MX_GRP) GOTO 500
       DO I=1,3
        DO J=1,3
         RMAT(J,I,1)= 0.0D0
        END DO
        RMAT(I,I,1)= 1.0D0
       END DO
C
C READ GRPMAT
C IF UNAVAILABLE, WRITE DEFAULT GRPMAT
C
       INQUIRE(FILE='GRPMAT',EXIST=EXIST)
       IF (.NOT. EXIST) THEN
        OPEN(77,FILE='GRPMAT',FORM='FORMATTED',STATUS='NEW')
        REWIND(77)
        WRITE(77,*) NGRP
        DO IGRP=1,NGRP
         WRITE(77,1000)((RMAT(J,I,IGRP),J=1,3),I=1,3)
         WRITE(77,*)' '
        END DO
       ELSE
        OPEN(77,FILE='GRPMAT',FORM='FORMATTED',STATUS='OLD')
        REWIND(77)
        READ(77,*,END=600) NGRP
        IF (NGRP .GT. MX_GRP) THEN
         CLOSE(77)
         GOTO 500
        END IF
        IF (NGRP .LT. 1) GOTO 600
        DO IGRP=1,NGRP
         READ(77,*,END=600)((RMAT(J,I,IGRP),J=1,3),I=1,3)
        END DO
C
C CHECK IF IDENTITY MATRIX IS FIRST IN TABLE
C
        ERROR=0.0D0
        DO I=1,3
         DO J=1,3
          IF (I .EQ. J) THEN
           ERROR=ERROR+ABS(RMAT(I,I,1)-1.0D0)
          ELSE
           ERROR=ERROR+ABS(RMAT(J,I,1))
          END IF
         END DO
        END DO
        IF (ERROR .GT. TOLER) THEN
         PRINT *,'FGMAT: IDENTITY MATRIX MUST BE FIRST IN GRPMAT'
         CLOSE(77)
         CALL STOPIT
        END IF
C
C CHECK TO MAKE SURE THAT MATRICES FORM A GROUP
C
  10    CONTINUE
        MGRP=NGRP
        DO 100 IGRP=1,NGRP
         DO 90 JGRP=1,NGRP
          DO I=1,3
           DO J=1,3
            PMAT(J,I)=0
            DO K=1,3
             PMAT(J,I)=RMAT(J,K,IGRP)*RMAT(K,I,JGRP)+PMAT(J,I)
            END DO
           END DO
          END DO
          ITIMES=0
          DO 20 KGRP=1,MGRP
           ERROR=0.0D0
           DO I=1,3
            DO J=1,3
             ERROR=ERROR+ABS(PMAT(J,I)-RMAT(J,I,KGRP))
            END DO
           END DO
           IF (ERROR .LT. TOLER) ITIMES=ITIMES+1
   20     CONTINUE
C
          IF (ITIMES.EQ.0) THEN
           PRINT *,'FGMAT: PRODUCT NOT IN GROUP: ',IGRP,JGRP
           MGRP=MGRP+1
           IF (MGRP .LE. MX_GRP) THEN
            DO I=1,3
             DO J=1,3
              RMAT(J,I,MGRP)=PMAT(J,I)
             END DO
            END DO
           ELSE 
            NGRP=MGRP
            CLOSE(77)
            GOTO 500
           END IF
          ELSE IF (ITIMES .GT. 1) THEN
           PRINT *,'BIZARRE ERROR IN REPRESENTATION MATRICES'
           PRINT *,'ARE THERE TWO IDENTICAL MATRICES IN INPUT FILE ?'
           PRINT *,'PROGRAM CRASHED WITH: '
           PRINT *,'IGRP=',IGRP
           DO I=1,3
            PRINT *,(RMAT(J,I,IGRP),J=1,3)
           END DO
           PRINT *,'JGRP=',JGRP
           DO I=1,3
            PRINT *,(RMAT(J,I,JGRP),J=1,3)
           END DO
           CALL STOPIT
          END IF
   90    CONTINUE
  100   CONTINUE
C
        IF (MGRP.NE.NGRP) THEN
         NGRERR=NGRERR+1
         PRINT *,'FGMAT: ERROR IN GRPMAT: NERR,MGRP=',NGRERR,MGRP
         NGRP=MGRP
         IF (NGRERR .GT. 100) CALL STOPIT
         GOTO 10
        END IF
C
        IF (NGRERR.NE.0) THEN
        REWIND(77)
        BEEN_CALLED=.FALSE.
C        WRITE(77,*) 'A LARGER GROUP OF ORDER',NGRP,' HAS BEEN FOUND'
C        WRITE(77,*)'GROUP REPRESENTATION:'
         WRITE(77,*) NGRP
         DO IGRP=1,NGRP
          WRITE(77,1000)((RMAT(J,I,IGRP),J=1,3),I=1,3)
          WRITE(77,*)' '
         END DO
         PRINT *,'FGMAT: FOUND LARGER GROUP: CHECK GRPMAT FILE'
C        CALL STOPIT
         REWIND(77)
         GO TO 15
 1000    FORMAT(' ',3G25.16)
        END IF
       END IF
       CLOSE(77)
       RETURN
C
  500  PRINT *,'FGMAT: MX_GRP MUST BE AT LEAST:',NGRP
       CALL STOPIT
  600  PRINT *,'FGMAT: GRPMAT IS BROKEN'
       CLOSE(77)
       CALL STOPIT
       END
