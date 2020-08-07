C UTEP Electronic Structure Lab (2020)
C
C
C
         SUBROUTINE TESTNSB
         IMPLICIT REAL*8 (A-H,O-Z)
         PARAMETER (NDIM=20)
         DIMENSION SS(NDIM,NDIM),WW(10,10),A(3),B(3) 
 1       CONTINUE
         write(6,*)'ALP,BET=?'
         READ*,ALP,BET,A,B
         CALL KINM3D(NDIM,ALP,BET,A,B,SS)
         CALL KNMXSF(ALP,BET,A,B,WW)
c        call ovlp3d(ndim,alp,bet,a,b,ss)
c        call ovmxsf(alp,bet,a,b,ww)
            PRINT 20,(SS(I,I),I=1,10)
            PRINT 20,(WW(I,I),I=1,10)
 20      FORMAT(' ',5G15.6)
         ERROR=0.0D0
                 DO I=1,10
                    IF(I.LE.4)THEN
                    PRINT 10,(WW(I,J),J=1,4)
                    PRINT 11,(SS(I,J),J=1,4)
 10      FORMAT(' OLD:',4G15.6)
 11      FORMAT(' NEW:',4G15.6)
                    END IF
                 DO J=1,10
                 DIFF=       ABS(WW(I,J)-SS(I,J))
                 IF(DIFF.GT.0.0001D0)THEN
                 write(6,*)I,J,WW(I,J),SS(I,J)
                 END IF
                 ERROR=ERROR+ABS(WW(I,J)-SS(I,J))
                 END DO
                 END DO
                 write(6,*)'ERROR:',ERROR
        GO TO 1
        END
