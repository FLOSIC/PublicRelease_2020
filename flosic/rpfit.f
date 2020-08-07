C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE RPFIT(IFNCT,R,RRC,RHO,POT)
       use common4,only : RPFALP, RPFCMX, RPFCOF, NRPFIT, LDIVR
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:01 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFNCT, I
       REAL*8 :: R , RRC, RHO, POT, ACCU, FAC
        SAVE
        DATA ACCU /1.0D-20/
C
        FAC= EXP(-RPFALP(IFNCT)*R*R)
        RHO= 0.0D0
        POT= 0.0D0
        DO I= 1,NRPFIT(IFNCT)
         RHO= RHO+RPFCOF(1,I,IFNCT)*FAC
         POT= POT+RPFCOF(2,I,IFNCT)*FAC
         FAC= FAC*FAC
         IF (FAC*RPFCMX(IFNCT) .LT. ACCU*MIN(RHO,POT)) GOTO 100
        END DO
  100   CONTINUE
        IF (LDIVR(IFNCT) .EQ. 1) THEN
         POT= -POT*RRC
        ELSE
         POT= -POT
        END IF
        RETURN
       END
