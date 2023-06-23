C UTEP Electronic Structure Lab (2020)
C> @file mldecomp.f
C> NEWWAVES
C> Path diverges depending on ATOMSIC and FRMORB files
       SUBROUTINE NEWWAVES(ITTOT,TRACE)
       use global_inputs,only : symmetrymodule1
       IMPLICIT REAL*8 (A-H,O-Z)
       CHARACTER*50 LINE
       LOGICAL EXIST
       INQUIRE(FILE='ATOMSIC', EXIST = EXIST)
       IF(.NOT.EXIST)THEN
         INQUIRE(FILE='FRMORB',EXIST=EXIST)
         IF(EXIST) THEN
           IF(symmetrymodule1) THEN
             CALL NEWWAVE_2020(ITTOT,TRACE)
           ELSE
             CALL SCFSICW(ITTOT,TRACE)
           END IF
         ELSE
           CALL NEWWAVE_serial(ITTOT,TRACE)
           !Temporary modifiction
           !CALL NEWWAVE(ITTOT,TRACE)
         END IF
       ELSE
         INQUIRE(FILE='PURGRSQ',EXIST=EXIST)
         OPEN(54,FILE='PURGRSQ')
         IF(EXIST)READ(54,*,END=10)EXIST
 10      CONTINUE
         IF(.NOT.EXIST)THEN
           REWIND(54)
           WRITE(54,*)'.TRUE.'
           PRINT*,'RERUN WITH NEW PURGRSQ'
           CALL STOPIT
         END IF
         CLOSE(54)
         OPEN(54,FILE='XMOL.DAT')
         READ(54,*)NATM
         CLOSE(54)
         OPEN(54,FILE='GRPMAT')
         READ(54,*)MXXG
         CLOSE(54)
         IF(NATM.NE.1.OR.MXXG.NE.1)THEN
           PRINT*,' ATOMS:',NATM
           PRINT*,'GRPMAT:',MXXG
           PRINT*,'SICWAVE SHOULD NOT BE CALLED'
           CALL STOPIT
         END IF
         !CALL SICWAVE(ITTOT,TRACE)
         print *,"SICWAVE is not merged yet in this version"
         call stopit  
       END IF
       RETURN
       END
