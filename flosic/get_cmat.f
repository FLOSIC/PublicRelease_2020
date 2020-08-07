C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
       SUBROUTINE GET_CMAT(ISYM,MSYM,ISHELL,NBTOT,CMAT)
C WRITTEN BY MARK R PEDERSON (1985)
       use common3,only : RMAT, NGRP
       use common8,only : S_REP, P_REP, D_REP, REP, RDENT,
     &   NUMSITES, IGGEN
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:47 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: ISYM, MSYM, ISHELL, NBTOT, I, IBEG, IGP, ISITE, J,
     & NBASE
       REAL*8 :: CMAT , ERROR, SITE, VEC
       SAVE
       DIMENSION VEC(3,MX_GRP),NBASE(3),CMAT(MX_GRP,6*MX_GRP),SITE(3)
       DATA NBASE/1,3,6/
C
C      THIS PROGRAM ACCEPTS:
C       1)A SHELL NUMBER:   ISHELL
C       2)A SYMMETRY TYPE: 0=S,1=P,2=D
C       3)A SUB SYMMETRY TYPE: 1   FOR S
C                          1 2 3   FOR PX,PY,PZ
C                      1 2 3 4 5 6 FOR XX,YY,ZZ,XY,XZ,YZ
C
C      IT RETURNS:
C
C      CMAT(I,J)
C         J=1,# OF FUNCTIONS  (1*NSITE FOR S)
C              (3*NSITE FOR P)
C              (6*NSITE FOR D)
C         I=1,NGRP            (DIMENSION OF GROUP)
       NBTOT=NBASE(ISYM+1)*NUMSITES(ISHELL)
       DO 1 J=1,NBTOT
       DO 1 I=1,NGRP
 1     CMAT(I,J)=0.0D0
       DO 10 ISITE=1,NUMSITES(ISHELL)
       DO 5 I=1,3
       VEC(I,ISITE)=0.0D0
       DO 4 J=1,3
       VEC(I,ISITE)=VEC(I,ISITE)+
     &   RMAT(I,J,IGGEN(ISITE,ISHELL))*RDENT(J,ISHELL)
 4     CONTINUE
 5     CONTINUE
 10    CONTINUE
       DO 50 IGP=1,NGRP
       DO 15 I=1,3
       SITE(I)=0.0D0
       DO 14 J=1,3
       SITE(I)=SITE(I)+RMAT(J,I,IGP)*RDENT(J,ISHELL)
 14    CONTINUE
 15    CONTINUE
       DO 25 ISITE=1,NUMSITES(ISHELL)
       ERROR=0.0D0
       DO 20 I=1,3
 20    ERROR=ERROR+ABS(SITE(I)-VEC(I,ISITE))
       IF(ERROR .LE. 1.0D-4)GO TO 30
 25    CONTINUE
       write(6,*)'MISTAKE IN GET_CMAT'
       write(6,*)'SITE GENERATED IS NOT IN SHELL:',ISHELL
       DO 35 ISITE=1,NUMSITES(ISHELL)
        PRINT 26,(VEC(I,ISITE),I=1,3)
   35  CONTINUE
       write(6,*)'GET_CMAT: BAD VECTOR:'
       PRINT 26,SITE
 26    FORMAT(' ',3G15.6)
       write(6,*)'CALLING CKRMAT AND QUITTING ...'
       CALL CKRMAT
       CALL STOPIT
 30    CONTINUE
       IBEG=(ISITE-1)*NBASE(ISYM+1)
       DO 40 J=1,NBASE(ISYM+1)
       IF(ISYM.EQ.0)THEN
       CMAT(IGP,J+IBEG)=S_REP(IGP)
       ELSE IF(ISYM.EQ.1)THEN
       CMAT(IGP,J+IBEG)=P_REP(MSYM,J,IGP)
       ELSE
       CMAT(IGP,J+IBEG)=D_REP(MSYM,J,IGP)
       END IF
 40    CONTINUE
 50    CONTINUE
       RETURN
       END
