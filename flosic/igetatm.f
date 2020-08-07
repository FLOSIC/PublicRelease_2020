C UTEP Electronic Structure Lab (2020)
C
C ******************************************************************
C
       SUBROUTINE IGETATM(IUNIT,R,IZ,NAME,IREAD)
         IMPLICIT REAL*8 (A-H,O-Z)
         DIMENSION    R(3)
         SAVE
C
         CHARACTER*6  NAME
         CHARACTER*1  CCH
         CHARACTER*80 LINE 
         IREAD=1
         READ(IUNIT,'(A80)',END=100) LINE
         K=81
         DO I=80,1,-1
           IF(LINE(I:I).EQ.'!')K=I
         END DO
         DO I=K,80
           LINE(I:I)=' '
         END DO
         READ(LINE,*,END=100) R,IZ
         IPOS=0
         IBLANK=1
         DO IJUMP=1,5
   10      CONTINUE
           IPOS=IPOS+1
           IF (IPOS .GT. 80) GOTO 100
           CCH=LINE(IPOS:IPOS)
           IF (IBLANK .EQ. 1) THEN
             IF ((CCH .NE. ' ') .AND. (CCH .NE. '\t') 
     &                          .AND. (CCH .NE. ',')) THEN
               IF (IJUMP .EQ. 5) GOTO 20
                 IBLANK=0
               END IF
             GOTO 10
           ELSE
             IF ((CCH .NE. ' ') .AND. (CCH .NE. '\t') 
     &                         .AND. (CCH .NE. ',')) THEN
               GOTO 10
             ELSE
               IBLANK=1
             END IF
           END IF
           CONTINUE
         END DO
   20    CONTINUE 
         K=0
         DO I=IPOS,80
           IF(LINE(I:I).NE.' ')THEN
             K=K+1
             LINE(K:K)=LINE(I:I)
           END IF
         END DO  
         NAME=LINE(1:6)
         IREAD=0
  100    RETURN
       END
