C UTEP Electronic Structure Lab (2020)
C
C ***********************************************************
C
       SUBROUTINE ADDPTS(NSITES,R_NUC,DELR,WTR,ANGLE,DOMEGA,
     &                   NR,N_ANGLE)
C ORIGINAL VERSION BY MARK R PEDERSON (1989)
       use mesh1,only : wmsh,rmsh,nmsh
       use common3,only : RMAT, NGRP
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:46 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NSITES, NR, N_ANGLE, IANG, IGRP, INDEX, IRAD, J,
     & JANG, JGRP, K, M_NUC, MANG, MGRP, MSITES, MXANG
       REAL*8 :: SYMBOL , R_NUC, DELR, WTR, ANGLE, DOMEGA, DIST, DVJ,
     & DVOL, RR, RVEC, TOLER, VEC
       SAVE
       PARAMETER (MXANG=400) !DFT:200 SIC: 400
       DIMENSION R_NUC(3),ANGLE(3,N_ANGLE),DOMEGA(N_ANGLE)
       DIMENSION DELR(NR),WTR(NR)
       DIMENSION VEC(3),RVEC(3,MX_GRP),MSITES(1)
       DIMENSION INDEX(MX_GRP)
       DATA TOLER/1.0D-8/
C
C CALL GASITES AND CHECK PARAMETERS
C
       CALL GASITES(1,R_NUC,M_NUC,RVEC,MSITES)
       IF(MSITES(1).NE.NSITES)THEN
        write(6,*)'ADDPTS: MSITES AND NSITES ARE NOT EQUAL'
        CALL STOPIT
       END IF
       IF(N_ANGLE.GT.MXANG)THEN
        write(6,*)'ADDPTS: MXANG MUST BE AT LEAST: ',N_ANGLE
        CALL STOPIT
       END IF
C
C FIND GROUP ELEMENTS WHICH LEAVE R_NUC INVARIANT:
C
       MGRP=0
       DO 30 IGRP=1,NGRP
        DO 20 J=1,3
         VEC(J)=0.0D0
         DO 10 K=1,3
          VEC(J)=VEC(J)+RMAT(J,K,IGRP)*R_NUC(K)
   10    CONTINUE
   20   CONTINUE
        DVJ=SQRT((VEC(1)-R_NUC(1))**2+
     &           (VEC(2)-R_NUC(2))**2+
     &           (VEC(3)-R_NUC(3))**2)
        IF(DVJ .LT. TOLER)THEN
         MGRP=MGRP+1
         INDEX(MGRP)=IGRP
        END IF
   30  CONTINUE
C
C REDUCE ANGULAR MESH TO INEQUIVALENT POINTS
C
       DO 60 IANG=1,N_ANGLE
        IF (DOMEGA(IANG) .NE. 0.0D0) THEN
         DO 50 JANG=IANG+1,N_ANGLE
          IF (DOMEGA(JANG) .NE. 0.0D0) THEN
           DO 40 JGRP=1,MGRP
            IGRP=INDEX(JGRP)
            DO J=1,3
             VEC(J)=0.0D0
             DO K=1,3
              VEC(J)=VEC(J)+RMAT(J,K,IGRP)*ANGLE(K,JANG)
             END DO
            END DO
            DIST=(VEC(1)-ANGLE(1,IANG))**2
     &          +(VEC(2)-ANGLE(2,IANG))**2
     &          +(VEC(3)-ANGLE(3,IANG))**2
            IF (DIST .LT. TOLER*TOLER) THEN
             DOMEGA(IANG)=DOMEGA(IANG)+DOMEGA(JANG)
             DOMEGA(JANG)=0.0D0
             GOTO 50
            END IF
   40      CONTINUE
          END IF
   50    CONTINUE
        END IF
   60  CONTINUE
       MANG=0
       DO IANG=1,N_ANGLE
        IF (DOMEGA(IANG) .NE. 0.0D0) THEN
         MANG=MANG+1
         DOMEGA(MANG)=DOMEGA(IANG)
         ANGLE(1,MANG)=ANGLE(1,IANG)
         ANGLE(2,MANG)=ANGLE(2,IANG)
         ANGLE(3,MANG)=ANGLE(3,IANG)
        END IF
       END DO
C
C FILL RMSH AND WMSH
C
       DO 100 IRAD=1,NR
        RR=DELR(IRAD)
        DVOL=WTR(IRAD)
        DO J=1,MANG
         NMSH=NMSH+1
         RMSH(1,NMSH)= RR*ANGLE(1,J)+R_NUC(1)
         RMSH(2,NMSH)= RR*ANGLE(2,J)+R_NUC(2)
         RMSH(3,NMSH)= RR*ANGLE(3,J)+R_NUC(3)
         WMSH(NMSH) = DOMEGA(J)*DVOL*NSITES
        END DO  
  100  CONTINUE
       RETURN 
       END
