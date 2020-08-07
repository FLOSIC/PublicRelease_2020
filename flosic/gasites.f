C UTEP Electronic Structure Lab (2020)
C
C *******************************************************
C 
       SUBROUTINE GASITES(N_NUC,R_NUC,M_NUC,R_NUCA,MSITES)
C WRITTEN BY MARK R PEDERSON (1985)
       use common3,only : RMAT, NGRP
       use common8,only : IGEN
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:46 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: N_NUC, M_NUC, MSITES, I_BEG, I_NUC, IC, IGRP, ITOT,
     & J_NUC, JC
       REAL*8 :: R_NUC , R_NUCA, DISTANCE, ERROR, TOL
       SAVE
       DIMENSION R_NUC(3,N_NUC),R_NUCA(3,N_NUC*MX_GRP),MSITES(N_NUC)
       DATA TOL/1.0D-8/
       M_NUC=0
       DO 50 I_NUC=1,N_NUC
        DISTANCE=R_NUC(1,I_NUC)**2+R_NUC(2,I_NUC)**2+R_NUC(3,I_NUC)**2
        DISTANCE=SQRT(DISTANCE)
        IF (DISTANCE .GT. TOL) THEN
         I_BEG=M_NUC+1
         MSITES(I_NUC)=0
         DO 25 IGRP=1,NGRP
          M_NUC=M_NUC+1
          MSITES(I_NUC)=MSITES(I_NUC)+1
          DO IC=1,3
           R_NUCA(IC,M_NUC)=0.0D0
           DO JC=1,3
            R_NUCA(IC,M_NUC)=R_NUCA(IC,M_NUC)+RMAT(IC,JC,IGRP)
     &                      *R_NUC(JC,I_NUC)
           END DO
          END DO
C
C IS THIS A NEW SITE?
C
          ITOT=0
          DO J_NUC=I_BEG,M_NUC-1
           ERROR=ABS(R_NUCA(1,M_NUC)-R_NUCA(1,J_NUC))
     &          +ABS(R_NUCA(2,M_NUC)-R_NUCA(2,J_NUC))
     &          +ABS(R_NUCA(3,M_NUC)-R_NUCA(3,J_NUC))
           IF (ERROR .LE. TOL*DISTANCE) ITOT=ITOT+1
          END DO
          IF (ITOT .NE. 0) THEN
           M_NUC=M_NUC-1
           MSITES(I_NUC)=MSITES(I_NUC)-1
          ELSE
           IF (M_NUC .LE. MX_GRP) IGEN(M_NUC)=IGRP
          END IF
   25    CONTINUE
        ELSE
         M_NUC=M_NUC+1
         R_NUCA(1,M_NUC)=R_NUC(1,I_NUC)
         R_NUCA(2,M_NUC)=R_NUC(2,I_NUC)
         R_NUCA(3,M_NUC)=R_NUC(3,I_NUC)
         IF (M_NUC .LE. MX_GRP) IGEN(M_NUC)=1
         MSITES(I_NUC)=1
        END IF
  50   CONTINUE
       RETURN
       END
