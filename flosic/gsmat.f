C UTEP Electronic Structure Lab (2020)
C
C *******************************************************
C
C DVP 10/98: GSMAT EXPECTS NOW THAT CREPMAT HAS BEEN CALLED BEFORE
C
       SUBROUTINE GSMAT(ISHELL,IPT)
C WRITTEN BY MARK R PEDERSON (1985)
        use debug1
       use common3,only : NGRP !RMAT
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, NUMSITES
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: ISHELL, IPT, I, I_REP, IB, IBEG, ICOL, IEND, !ICO,
     & IGP, IND, IROW, ISYM, ISYMMAX, J, JB, K, KB, KSALC, KSALCBEG,
     & MCOUNT, MSYM, NBTOT, NBTOTS, NDEG, NDSV, NIND, NINP, NTOT
       REAL(8) :: ADD , CMAT, DOT, T_MAT, ONEOSQRT3
        SAVE
        DATA ONEOSQRT3  /0.57735026918962576D0/
        LOGICAL LININD,EXIST,PURGRSQ
        DIMENSION NTOT(MAX_REP),NDEG(3)
        DIMENSION T_MAT(ISMAX,6*MX_GRP,5)
        DIMENSION CMAT(MX_GRP,6*MX_GRP)
               INQUIRE(FILE='PURGRSQ',EXIST=EXIST)
               IF(EXIST)THEN
               OPEN(90,FILE='PURGRSQ',FORM='FORMATTED')
               READ(90,*,END=5)PURGRSQ
 5             CONTINUE
               REWIND(90)
               WRITE(90,*)PURGRSQ,' T/F TO/NOT PURGE R^2 GAUSSIANS'
               CLOSE(90)
               END IF
C             
C
        NDEG(1)=1
        NDEG(2)=3
        NDEG(3)=6
        ISYMMAX=2
        DO 150 ISYM=0,ISYMMAX
         IF(ISYM.EQ.0)THEN
          IBEG=1
          IEND=1
         ELSE IF(ISYM.EQ.1)THEN
          IBEG=1
          IEND=3
         ELSE IF(ISYM.EQ.2)THEN
          IBEG=1
          IEND=6
         END IF
         KSALC=0
         DO I_REP=1,N_REP
          KSALC=KSALC+NDMREP(I_REP)
         END DO
         IF(KSALC.GT.MAXSYMSALC)THEN
          write(6,*)'MAXSYMSALC MUST BE AT LEAST:',KSALC
          CALL STOPIT
         END IF
         KSALCBEG=0
         DO 110 I_REP=1,N_REP
          NIND=0
          IF(ISYM.EQ.2.AND.PURGRSQ)THEN
C DEFINE FIRST SALC'S OF EACH REP TO BE S-LIKE...
             KSALC=KSALCBEG
             DO  IROW=1,NDMREP(I_REP)
              KSALC=KSALC+1
             DO NIND=1,N_SALC(KSALC,1,ISHELL)
              DO IB=1,NBTOTS
              DO KB=1,3
                JB=(IB-1)*6+KB
              T_MAT(NIND,JB  ,IROW)=U_MAT(IB,NIND,KSALC,1,IPT)
!     &                               /SQRT(3.0D0)
     &                               *ONEOSQRT3
              T_MAT(NIND,JB+3,IROW)=0.0D0
              END DO
              END DO
              DO IB=1,NBTOTS*6
               U_MAT(IB,NIND,KSALC,ISYM+1,IPT)=T_MAT(NIND,IB,IROW)
              END DO      
             END DO
              NIND=N_SALC(KSALC,1,ISHELL)
              N_SALC(KSALC,ISYM+1,ISHELL)=NIND
              NDSV=NIND
!101   FORMAT(' ',6g12.3)
             END DO      
          END IF
          DO 90 MSYM=IBEG,IEND
           CALL GET_CMAT(ISYM,MSYM,ISHELL,NBTOT,CMAT)
                      IF(ISYM.EQ.0)NBTOTS=NBTOT 
           IF(NBTOT.GT.LOCMAX)THEN
            write(6,*)'LOCMAX MUST BE AT LEAST:',NBTOT
            CALL STOPIT
           END IF
           DO 80 ICOL=1,NDMREP(I_REP)
            DO 15 IROW=1,NDMREP(I_REP)
C
C FIND A (ICOL,MSYM) SALC OF REP(I_REP,IROW)
C
             DO 16 K=1,NBTOT
              ADD=0.0D0
              DO 10 IGP=1,NGRP
               ADD=ADD+REP(ICOL,IROW,IGP,I_REP)*CMAT(IGP,K)
 10           CONTINUE
              T_MAT(NIND+1,K,IROW)=ADD
 16          CONTINUE
 15         CONTINUE
C CHECK FOR LINEAR DEPENDENCIES:
            DO IND=1,NIND
             DOT=0.0D0
             DO K=1,NBTOT
              DOT=DOT+T_MAT(NIND+1,K,1)*T_MAT(IND,K,1)
             END DO
             DO IROW=1,NDMREP(I_REP)
              DO K=1,NBTOT
               T_MAT(NIND+1,K,IROW)=T_MAT(NIND+1,K,IROW)
     &                             -DOT*T_MAT( IND  ,K,IROW)
              END DO
             END DO
            END DO
C FIND NORM OF NEW SALC:
            DOT=0.0D0
            DO K=1,NBTOT
             DOT=DOT+T_MAT(NIND+1,K,1)*T_MAT(NIND+1,K,1)
            END DO
            IF(DOT.GE.1.0D-4)THEN
             DOT=1.0D0/SQRT(DOT)
             DO IROW=1,NDMREP(I_REP)
              DO K=1,NBTOT
               T_MAT(NIND+1,K,IROW)=T_MAT(NIND+1,K,IROW)*DOT
              END DO
             END DO
C CHECK ORTHOGONALITY:
C            DO IND=1,NIND+1
C             DOT=0.0D0
C             DO K=1,NBTOT
C              DOT=DOT+T_MAT(NIND+1,K,1)*T_MAT(IND,K,1)
C             END DO
C             write(6,*)IND,NIND+1,DOT
C            END DO
             LININD=.TRUE.
             NIND=NIND+1
            ELSE
             LININD=.FALSE.
            END IF
            IF(LININD)THEN
             IF(DEBUG) write(6,*)'I_REP=',I_REP,' MSYM, ICOL=',MSYM,ICOL
C      MOVE SALC TO APPROPRIATE PLACE IN U_MAT:
             IF(NIND.GT.ISMAX)THEN
              write(6,*)'ISMAX MUST BE AT LEAST: ', NIND
              CALL STOPIT
             END IF
             KSALC=KSALCBEG
                IF(ISYM.EQ.2.AND.PURGRSQ)THEN
                  NINP=NIND-NDSV
                ELSE
                  NINP=NIND
                END IF
             DO 25 IROW=1,NDMREP(I_REP)
              KSALC=KSALC+1
              DO 20 IB=1,NBTOT
               U_MAT(IB,NINP,KSALC,ISYM+1,IPT)=T_MAT(NIND,IB,IROW)
 20           CONTINUE
 25          CONTINUE
            END IF
  80       CONTINUE
  90      CONTINUE
          KSALC=KSALCBEG 
                IF(ISYM.EQ.2.AND.PURGRSQ)THEN
                  NINP=NIND-NDSV
                ELSE
                  NINP=NIND
                END IF
          DO IROW=1,NDMREP(I_REP)
           KSALC=KSALC+1
           N_SALC(KSALC,ISYM+1,ISHELL)=NINP
           NTOT(I_REP)=NINP
          END DO
          KSALCBEG=KSALC
 110     CONTINUE
! 111     FORMAT(' SHELL:',I3,' R=',3G15.6)
         MCOUNT=0
         DO 105 I=1,N_REP
          DO 102 J=1,NDMREP(I)
           MCOUNT=MCOUNT+NTOT(I)
 102      CONTINUE
          IF (DEBUG) PRINT *,I,ISYM+1,(NTOT(I),J=1,NDMREP(I))
! 103      FORMAT(' REPRESENTATION, SYM:',2I3,' PROJECTIONS:',10I3)
 105     CONTINUE
         IF (DEBUG) PRINT *,'TOTAL # OF PROJECTIONS:',MCOUNT,
     &                      ' EXPECTED:',NUMSITES(ISHELL)*NDEG(ISYM+1)
 150    CONTINUE
        RETURN
        END
