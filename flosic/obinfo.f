C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE OBINFO(N_NUC,R_NUC,R_NUCA,M_NUC,KSHELL)
C ORIGINALLY WRITTEN BY M.R. PEDERSON (1985)
       use debug1
       use common3,only : RMAT, NGRP
       use common8,only : IGEN, RDENT, NUMSITES, IGGEN, N_IDNT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:55 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: N_NUC, M_NUC, KSHELL, I, I_NUC, IC, IGRP, ISHELL,
     & ITOT, J, J_NUC, JC, NEW
       REAL*8 :: R_NUC , R_NUCA, DISTANCE, ERROR, TOL
       SAVE
       LOGICAL FIRST
       DIMENSION R_NUC(3,N_NUC),R_NUCA(3,N_NUC*MX_GRP)
C
C THIS PROGRAM RECEIVES:
C
C N_NUC = TOTAL NUMBER OF LATTICE POINTS
C R_NUC(J,I_NUC) J=1,3 I_NUC=1,N_NUC  : LATTICE POINTS
C R_NUCA                              : STORAGE
C
C IT FINDS ALL THE DIFFERENT KINDS OF IDENTITY MEMBERS AND STORES THEM
C IN RDENT
C
C N_IDNT = TOTAL NUMBER OF IDENTITY MEMBERS
C NUMSITES(ISHELL) = NUMBER OF EQUIVALENT LATTICE POINTS
C
C PROGRAM RETURNS:
C
C     M_NUC NUMBER OF SITES ASSOCIATED WITH SHELL N_NUC
C     KSHELL = SHELL TYPE FOR SHELL N_NUC
C
       DATA TOL/1.0D-8/
       DATA FIRST/.TRUE./
C
       IF(FIRST)THEN
        FIRST=.FALSE.
        N_IDNT=0
       END IF
       DO 55 I_NUC=1,N_NUC
        DISTANCE=R_NUC(1,I_NUC)**2+R_NUC(2,I_NUC)**2+R_NUC(3,I_NUC)**2
        DISTANCE=SQRT(DISTANCE)
        IF(DISTANCE.GT.TOL)THEN
         M_NUC=0
         DO 25 IGRP=1,NGRP
          M_NUC=M_NUC+1
          DO 10 IC=1,3
           R_NUCA(IC,M_NUC)=0.0D0
           DO JC=1,3
            R_NUCA(IC,M_NUC)=R_NUCA(IC,M_NUC)+RMAT(IC,JC,IGRP)
     &                      *R_NUC(JC,I_NUC)
           END DO
   10     CONTINUE
C
C IS THIS A NEW SITE?
C
          ITOT=0
          DO 15 J_NUC=1,M_NUC-1
           ERROR=ABS(R_NUCA(1,M_NUC)-R_NUCA(1,J_NUC))
     &          +ABS(R_NUCA(2,M_NUC)-R_NUCA(2,J_NUC))
     &          +ABS(R_NUCA(3,M_NUC)-R_NUCA(3,J_NUC))
           IF(ERROR.LE.TOL*DISTANCE)THEN
            ITOT=ITOT+1
           END IF
 15       CONTINUE
          IF(ITOT.NE.0)THEN
           M_NUC=M_NUC-1
          ELSE
           IGEN(M_NUC)=IGRP
          END IF
 25      CONTINUE
        ELSE
         M_NUC=1
         R_NUCA(1,M_NUC)=R_NUC(1,I_NUC)
         R_NUCA(2,M_NUC)=R_NUC(2,I_NUC)
         R_NUCA(3,M_NUC)=R_NUC(3,I_NUC)
         IGEN(M_NUC)=1
        END IF
C
C CHECK TO SEE IF THIS IS A NEW SHELL TYPE
C
        NEW=0
        DO 45 ISHELL=1,N_IDNT
         ERROR=0.0D0
         DO 40 I=1,3
          ERROR=ERROR+ABS(R_NUC(I,I_NUC)-RDENT(I,ISHELL))
 40      CONTINUE
         IF(ERROR .LT. 1.0D-6)THEN
          KSHELL=ISHELL
          NEW=NEW+1
         END IF
 45     CONTINUE
        IF(NEW.GT.1)THEN
         PRINT *,'OBINFO: NEW > 1'
         CALL STOPIT
        ELSE IF(NEW.EQ.0)THEN
         N_IDNT=N_IDNT+1
         KSHELL=N_IDNT
         RDENT(1,N_IDNT)=R_NUC(1,I_NUC)
         RDENT(2,N_IDNT)=R_NUC(2,I_NUC)
         RDENT(3,N_IDNT)=R_NUC(3,I_NUC)
         NUMSITES(N_IDNT)=M_NUC
         DO 50 IGRP=1,M_NUC
          IGGEN(IGRP,N_IDNT)=IGEN(IGRP)
 50      CONTINUE
         IF (DEBUG) THEN
          write(6,*)'NEW SHELL, N_IDNT=',N_IDNT
          write(6,*)(RDENT(J,N_IDNT),J=1,3)
          write(6,*)'GENERATORS:'
          PRINT 105,(IGGEN(IGRP,N_IDNT),IGRP=1,M_NUC)
         END IF
        END IF
 55    CONTINUE
 100   FORMAT(2(1X,I5),3(1X,G15.6))
 105   FORMAT(18(1X,I3))
       RETURN
       END
