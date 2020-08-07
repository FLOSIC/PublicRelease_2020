C UTEP Electronic Structure Lab (2020)
C
C *******************************************************
C 
      SUBROUTINE FRCSYM(RVEC,FVEC,FDIFF)
C ORIGINAL VERSION BY MARK R PEDERSON (1990) 
C 
C DIRK POREZAG, AUGUST 1997
C PROJECT OUT FORCE COMPONENTS THAT VIOLATE POINT GROUP SYMMETRY 
C INPUT:  COORDINATES OF THE ATOMIC LOCATION IN RVEC
C         UNSYMMETRIZED FORCE IN FVEC
C OUTPUT: SYMMETRIZED FORCE IN FVEC 
C         |F(INPUT)-F(OUTPUT)| IN FDIFF
C
       use common3,only : RMAT, NGRP
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IGRP, J, MGRP
       REAL*8 :: RVEC , FVEC, FDIFF, COORD, DIFF, FAC, FSYMM, ZER2,
     & ZERO
       SAVE
       DIMENSION RVEC(3),FVEC(3)
       DIMENSION FSYMM(3)
       DATA ZERO/5.0D-4/
C
C SUM OVER ALL GROUP OPERATIONS THAT LEAVE RVEC UNCHANGED
C
       ZER2=ZERO**2
       DO I=1,3
        FSYMM(I)= 0.0D0
       END DO
       MGRP=0
       DO IGRP=1,NGRP
        DIFF= 0.0D0
        DO I=1,3
         COORD= 0.0D0
         DO J=1,3
          COORD= COORD+RMAT(J,I,IGRP)*RVEC(J)
         END DO
         DIFF= DIFF+(COORD-RVEC(I))**2
        END DO
        IF (DIFF .LT. ZER2) THEN
         MGRP= MGRP+1
         DO I=1,3
          DO J=1,3
           FSYMM(I)= FSYMM(I)+RMAT(J,I,IGRP)*FVEC(J)
          END DO
         END DO
        END IF
       END DO
C
C FINAL CORRECTION, DEFINE FDIFF 
C
       IF (MGRP .LT. 1) THEN
        write(6,*)'FRCSYM: WHERE THE HECK IS THE IDENTITY MATRIX ?'
        CALL STOPIT
       END IF
       FAC= 1.0D0/MGRP
       FDIFF= 0.0D0
       DO I=1,3
        FSYMM(I)= FSYMM(I)*FAC
        FDIFF= FDIFF+(FVEC(I)-FSYMM(I))**2
        FVEC(I)= FSYMM(I)
       END DO
       FDIFF= SQRT(FDIFF)
       RETURN
      END  
