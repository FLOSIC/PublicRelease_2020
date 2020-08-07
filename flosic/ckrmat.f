C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
       SUBROUTINE CKRMAT
C WRITTEN BY MARK R PEDERSON (1985)
       use common3,only : RMAT, NGRP, MULTAB
       use common8,only : S_REP, P_REP, D_REP, REP, N_REP, NDMREP
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:38 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, ICOL, IFAIL, IGP, II, IJ, INDEX, IREP, IROW, J,
     & JI, JJ, JREP, K, KP, L, M, N
       REAL*8 :: ADD , ANS, D, ERROR, FLTDIM, FLTNGP
       SAVE
       LOGICAL IFCHK,OK
       DIMENSION ANS(6,6),INDEX(3,3)
       real *8 :: det
       external det
! Small changes - Raja Oct. 2015.

C
C CONSTRUCT S,P AND D REPRESENTATION MATRICES:
C
       DO 10 IGP=1,NGRP
       S_REP(IGP)=1.0D0
   10  CONTINUE
C
       DO 25 IGP=1,NGRP
       DO 20 I=1,3
       DO 15 J=1,3
       P_REP(J,I,IGP)=RMAT(J,I,IGP)
   15  CONTINUE
   20  CONTINUE
   25  CONTINUE
C
       INDEX(1,1)=1
       INDEX(2,2)=2
       INDEX(3,3)=3
       INDEX(2,1)=4
       INDEX(1,2)=4
       INDEX(3,1)=5
       INDEX(1,3)=5
       INDEX(3,2)=6
       INDEX(2,3)=6
C CREATE "D"-REPRESENTATION FROM P REPRESENTATION:
       DO 110 IGP=1,NGRP
       DO 35 I=1,6
       DO 30 J=1,6
        D_REP(J,I,IGP)=0.0D0
   30  CONTINUE
   35  CONTINUE
       DO 80 I=1,3
       DO 75 J=I,3
       IROW=INDEX(I,J)
       DO 55 L=1,3
       DO 45 M=1,3
        ICOL=INDEX(L,M)
        D_REP(IROW,ICOL,IGP)=D_REP(IROW,ICOL,IGP)+
     &     P_REP(I,L,IGP)*P_REP(J,M,IGP)
   45  CONTINUE
   55  CONTINUE
   75  CONTINUE
   80  CONTINUE
  110  CONTINUE
C
C      CHECK DETERMINANTS:
C
       IFAIL=0
       IFCHK=.FALSE.
       IF(IFCHK)THEN
        PRINT *,'DETERMINANT CHECK:'
        DO 5 IREP=1,N_REP
         OK=.TRUE.
         DO 3 K=1,NGRP
          D=DET(NDMREP(IREP),REP(1,1,K,IREP))
          IF(ABS(ABS(D)-1.0D0).GT.1.0D-5)THEN
           IFAIL=IFAIL+1
           write(6,*)'IREP,K,DET=',IREP,K,D
           OK=.FALSE.
          END IF
 3       CONTINUE
         IF(OK)THEN
          write(6,*)'DETERMINANTS ARE OK FOR REP:',IREP
         END IF
 5      CONTINUE
       END IF
C      CHECK GREAT ORTHOGONALITY THEOREM:
       IFCHK=.TRUE.
       IF(IFCHK)THEN
       PRINT '(A)','CHECKING GREAT ORTHOGONALITY THEOREM'
       DO 50 IREP=1,N_REP
       DO 50 JREP=IREP,N_REP
       DO 50 II=1,NDMREP(IREP)
       DO 50 JI=1,NDMREP(IREP)
       DO 50 IJ=1,NDMREP(JREP)
       DO 50 JJ=1,NDMREP(JREP)
       ADD=0.0D0
       DO 40 K=1,NGRP
   40  ADD=ADD+REP(II,JI,K,IREP)*REP(IJ,JJ,K,JREP)
       IF(IREP.EQ.JREP.AND.II.EQ.IJ.AND.JI.EQ.JJ)THEN
       FLTNGP=NGRP
       FLTDIM=NDMREP(IREP)
       ERROR=ABS(ADD-FLTNGP/FLTDIM)
       ELSE
       ERROR=ABS(ADD)
       END IF
       IF(ERROR.GT.1.0D-4)THEN
        write(6,*)'CKRMAT: FAILURE OF GREAT ORTHONGONALITY THEOREM:'
        PRINT 60,IREP,JREP,II,IJ,JI,JJ,ADD
        IFAIL=IFAIL+1
       END IF
 50    CONTINUE
 60    FORMAT(' ',2I3,'   ',2I3,'   ',2I3,'   ',G15.6)
       END IF
C      CHECK OF MULTIPLICATION TABLE
       IFCHK=.TRUE.
       IF(IFCHK)THEN
        PRINT '(A)','CHECKING MULTIPLICATION TABLE'
        DO 90 IREP=1,N_REP
         OK=.TRUE.
         DO 85 I=1,NGRP
         DO 85 J=1,NGRP
          DO 65 K=1,NDMREP(IREP)
          DO 65 L=1,NDMREP(IREP)
          ANS(K,L)=0.0D0
          DO 65 N=1,NDMREP(IREP)
           ANS(K,L)=REP(K,N,I,IREP)*REP(N,L,J,IREP)
     &             +ANS(K,L)
 65       CONTINUE
          KP=MULTAB(I,J)
C      COMPARE TO ANSWER IN MULTIPLICATION TABLE
          ERROR=0.0D0
          DO 70 K=1,NDMREP(IREP)
          DO 70 L=1,NDMREP(IREP)
           ERROR=ERROR+ABS(ANS(K,L)-REP(K,L,KP,IREP))
 70       CONTINUE
          IF(ERROR.GT.1.0D-5)THEN
           IFAIL=IFAIL+1
           PRINT *,'CKRMAT: ANS AND REP DIFFER'
           write(6,*)I,J,ERROR
           OK=.FALSE.
          END IF
  85     CONTINUE
         IF (OK) PRINT '(A,I2,A)','REPRESENTATION ',IREP,' IS OK'
  90    CONTINUE
C
        PRINT '(A)','CHECKING GROUP REPRESENTATION'
        OK=.TRUE.
        DO 105 I=1,NGRP
        DO 105 J=1,NGRP
         DO 95 K=1,3
         DO 95 L=1,3
         ANS(K,L)=0.0D0
         DO 95 N=1,3
          ANS(K,L)=RMAT(K,N,I)*RMAT(N,L,J)+ANS(K,L)
 95      CONTINUE
         KP=MULTAB(I,J)
C      COMPARE TO ANSWER IN MULTIPLICATION TABLE
         ERROR=0.0D0
         DO 100 K=1,3
         DO 100 L=1,3
          ERROR=ERROR+ABS(ANS(K,L)-RMAT(K,L,KP))
 100     CONTINUE
         IF(ERROR.GT.1.0D-5)THEN
          IFAIL=IFAIL+1
          PRINT *,'CKRMAT: ANS AND RMAT DIFFER'
          write(6,*)I,J,ERROR
          OK=.FALSE.
         END IF
 105    CONTINUE
        IF (OK) PRINT '(A)','GROUP REPRESENTATION IS OK'
        OK=.TRUE.
        PRINT '(A)','CHECKING S REPRESENTATION'
        DO 205 I=1,NGRP
        DO 205 J=1,NGRP
         DO 195 K=1,1
         DO 195 L=1,1
         ANS(K,L)=0.0D0
         DO 195 N=1,1
          ANS(K,L)=S_REP(I)*S_REP(J)
     &            +ANS(K,L)
 195     CONTINUE
         KP=MULTAB(I,J)
C      COMPARE TO ANSWER IN MULTIPLICATION TABLE
         ERROR=0.0D0
         DO 200 K=1,1
         DO 200 L=1,1
          ERROR=ERROR+ABS(ANS(K,L)-S_REP(KP))
 200     CONTINUE
         IF(ERROR.GT.1.0D-5)THEN
          IFAIL=IFAIL+1
          PRINT *,'CKRMAT: ANS AND S_REP DIFFER'
          write(6,*)I,J,ERROR
          OK=.FALSE.
         END IF
 205    CONTINUE
        IF (OK) PRINT '(A)','S REPRESENTATION IS OK'
        OK=.TRUE.
        PRINT '(A)','CHECKING P REPRESENTATION'
        DO 305 I=1,NGRP
        DO 305 J=1,NGRP
         DO 295 K=1,3
         DO 295 L=1,3
         ANS(K,L)=0.0D0
         DO 295 N=1,3
          ANS(K,L)=P_REP(K,N,I)*P_REP(N,L,J)
     &            +ANS(K,L)
 295     CONTINUE
         KP=MULTAB(I,J)
C      COMPARE TO ANSWER IN MULTIPLICATION TABLE
         ERROR=0.0D0
         DO 300 K=1,3
         DO 300 L=1,3
          ERROR=ERROR+ABS(ANS(K,L)-P_REP(K,L,KP))
 300     CONTINUE
         IF(ERROR.GT.1.0D-5)THEN
          IFAIL=IFAIL+1
          PRINT *,'CKRMAT: ANS AND P_REP DIFFER'
          write(6,*)I,J,ERROR
          OK=.FALSE.
         END IF
 305    CONTINUE
        IF (OK) PRINT '(A)','P REPRESENTATION IS OK'
        OK=.TRUE.
        PRINT '(A)','CHECKING D REPRESENTATION:'
        DO 405 I=1,NGRP
        DO 405 J=1,NGRP
         DO 395 K=1,6
         DO 395 L=1,6
         ANS(K,L)=0.0D0
         DO 395 N=1,6
          ANS(K,L)=D_REP(K,N,I)*D_REP(N,L,J)
     &            +ANS(K,L)
 395     CONTINUE
         KP=MULTAB(I,J)
C      COMPARE TO ANSWER IN MULTIPLICATION TABLE
         ERROR=0.0D0
         DO 400 K=1,6
         DO 400 L=1,6
          ERROR=ERROR+ABS(ANS(K,L)-D_REP(K,L,KP))
 400     CONTINUE
         IF(ERROR.GT.1.0D-5)THEN
          IFAIL=IFAIL+1
          PRINT *,'CKRMAT: ANS AND D_REP DIFFER'
          write(6,*)I,J,ERROR
          OK=.FALSE.
         END IF
 405    CONTINUE
        IF (OK) PRINT '(A)','D REPRESENTATION IS OK'
       END IF
       PRINT '(A,I3)','TOTAL NUMBER OF FAILED TESTS: ',IFAIL
       IF (IFAIL.NE.0)THEN
C      CALL STOPIT
       write(6,*)' ATTEMPING TO RUN IN REDUCIBLE REPRESENTATION MODE!'
       END IF
       RETURN
       END
