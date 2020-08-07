C UTEP Electronic Structure Lab (2020)
       SUBROUTINE MACYS_SCISSOR(NBAS)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use common2, only : BFCON,N_BARE,N_CON
       use for_diag1
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:52 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NBAS, MAX_TOT, I, IBARE, IBAS, INDREP, ISPN, J, JBAS,
     & JSPN, K, KOUNT, L, ML, N22, N23, N42, N8, N_VIRT, NDEG, NN1,
     & NSAV, NSPN, NTEMP
       REAL*8 :: SYMBOL , CUTOCC, DELTA, DIAG, DLT, DOT, DSINGV, ERROR,
     & EVALSAV, HSSA, OCCTMP, POTIN, POTOUT, SCIS, TEMP, UMAT
!YY moved from common.inc
!       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!     &  EVAL(NDH),SC1(NDH),SC2(NDH)
!       INCLUDE 'commons.inc'
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       CHARACTER*9 QUANTUM
       DIMENSION NAMES(3)
       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
     &  NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
       DIMENSION DIAG(NDH,MAX_REP)
       DIMENSION NTEMP(MAX_TOT*MXSPN)
       COMMON/MIXPOT1/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN)
       DIMENSION UMAT(NDH,NDH),DELTA(NDH,2)
       COMMON/DELSCI/SCIS(200,200)

C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS IS DER-YOU's CONSTRUCTION ZONE...
       INQUIRE(FILE='DELTA_TEMPLATE',EXIST=EXIST)
       OPEN(40,FILE='DELTA_TEMPLATE')
       IF(.NOT.EXIST) THEN
            DO JSPN=1,NSPN
            IBAS=0
             DO L=0,2
                  IF(L.EQ.0)QUANTUM(1:1)='S'
                  IF(L.EQ.1)QUANTUM(1:1)='P'
                  IF(L.EQ.2)QUANTUM(1:1)='D'
                  IF(L.EQ.3)QUANTUM(1:1)='F'
               DO I=1,N_CON(L+1,1)
C SINGLE OR CONTRACTED GAUSSIAN?
                  NN1=0
                  DO IBARE=1,N_BARE(1)
                   IF(ABS(BFCON(IBARE,I,L+1,1)).GT.1.0D-5)NN1=NN1+1
                  END DO
                  IF(I.EQ.N_CON(L+1,1))NN1=2
               DO ML=1,2*L+1
                                DO K=2,8
                                QUANTUM(K:K)=' '
                                END DO
                   IF(L.EQ.1.AND.ML.EQ.1)QUANTUM(2:2)='X'
                   IF(L.EQ.1.AND.ML.EQ.2)QUANTUM(2:2)='Y'
                   IF(L.EQ.1.AND.ML.EQ.3)QUANTUM(2:2)='Z'
                   IF(L.EQ.2.AND.ML.EQ.1)QUANTUM(2:3)='XY'
                   IF(L.EQ.2.AND.ML.EQ.2)QUANTUM(2:3)='XZ'
                   IF(L.EQ.2.AND.ML.EQ.3)QUANTUM(2:3)='YZ'
                   IF(L.EQ.2.AND.ML.EQ.4)QUANTUM(2:8)='ZZ-RR/3'
                   IF(L.EQ.2.AND.ML.EQ.5)QUANTUM(2:6)='XX-YY'
                IBAS=IBAS+1
                IF(NN1.NE.1)THEN
                WRITE(40,40)JSPN,IBAS,I+L,QUANTUM(1:8)
                END IF
 40             FORMAT('0.000',2I4,I2,A8,I3)
               END DO
              END DO
             END DO
            END DO
       ENDIF
       DO JSPN=1,2
       DO KOUNT=1,NDH
       DELTA(KOUNT,JSPN)=0.0D0
       END DO
       END DO
       REWIND(40)
       DO KOUNT=1,NDH*2
       READ(40,*,END=45)DLT,JSPN,JBAS
              DELTA(JBAS,JSPN)=DLT
       END DO
 45    CONTINUE

        PRINT*,"CONSTRUCTION DIAGNOSTICS"
        DO I=1,N_CON(1,1)
        PRINT 1011,(OVER(I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=1,N_CON(1,1))
        END DO
        PRINT*,' '
        DO I=1,N_CON(1,1)
        PRINT 1011,(HAM (I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=1,N_CON(1,1))
        END DO
        PRINT*,' '
        N8=N_CON(1,1)+1
        N22=N_CON(1,1)+3*N_CON(2,1)
        DO I=N8,N22
        PRINT 1011,(OVER(I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=N8,N22)
        END DO
        PRINT*,' '
        DO I=N8,N22
        PRINT 1011,(HAM(I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=N8,N22)
        END DO
        PRINT*,' '
        N23=N22+1
        N42=N22+5*N_CON(3,1)
        DO I=N23,N42
        PRINT 1011,(OVER(I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=N23,N42)
        END DO
        PRINT*,' '
        DO I=N23,N42
        PRINT 1011,(HAM (I,J)/SQRT(OVER(I,I)*OVER(J,J)),J=N23,N42)
        END DO
        PRINT*,' '
        PRINT 1012,(HAM(I,I)/OVER(I,I),I=1,NBAS)
        DO I=1,NBAS
           EVAL(I)=SQRT(OVER(I,I))
        END DO
C       DO I=1,NBAS
C          DO J=1,NBAS
C             OVER(J,I)=OVER(J,I)/SQRT(EVAL(I)*EVAL(J))
C          END DO
C       END DO
        DO I=1,NBAS
           DO J=1,NBAS
              FILO(J,I)=OVER(J,I)/EVAL(I)/EVAL(J)
           END DO
        END DO
                   CALL LOWDEN(NDH,NBAS,FILO,UMAT,EVAL,SC1)
                            DO I=1,NBAS
                              DO J=1,NBAS
                              UMAT(J,I)=UMAT(J,I)/SQRT(OVER(J,J))
                              END DO
                            END DO
C       PRINT*, 'UNITARY?'
C       PRINT 1012, (UMAT(I,I),I=1,NBAS)
 1011   FORMAT(20F7.2)
 1012   FORMAT(10F10.4)
C      PRINT*,'LOWDENOVERLAPS'
       ERROR=0.0D0
       DO I=1,NBAS
                 DOT=0.0D0
                 DO K=1,NBAS
                 DO L=1,NBAS
                 DOT=DOT+UMAT(K,I)*UMAT(L,I)*OVER(K,L)
                 END DO
                 END DO
c           PRINT*,'I,DOT:',I,DOT
       ERROR=ERROR+ABS(DOT-1.0D0)
          DO J=1,NBAS
              SC1(J)=0.0D0
             DO K=1,NBAS
                DO L=1,NBAS
                    SC1(J)=SC1(J)+UMAT(K,I)*UMAT(L,J)*OVER(K,L)
                END DO
             END DO
          IF(J.NE.I)ERROR=ERROR+ABS(SC1(J))
          END DO
C      PRINT 1012,(SC1(J),J=I,NBAS)
c      PRINT*, ' '
       END DO
       PRINT*,'LOWDIN ERROR IN MACYS:',ERROR
       IF(ERROR.GT.1.0D-8)CALL STOPIT
            DO I=1,NBAS !OLD ORBITALS
            DO K=1,NBAS !NEW ORBITAL
            FILO(I,K)=0.0D0
              DO L=1,NBAS
              FILO(I,K)=FILO(I,K)+UMAT(L,K)*OVER(L,I)
              END DO
            END DO
            END DO

        DO I=1,NBAS !OLD
           HSSA=HAM(I,I)
        DO J=1,NBAS !OLD
           SCIS(I,J)=0.0D0
           DO K=1,NBAS !NEW
c           HAM(I,J)=HAM(I,J)+DELTA(K,ISPN)*FILO(I,K)*FILO(J,K)
           SCIS(I,J)=SCIS(I,J)+DELTA(K,ISPN)*FILO(I,K)*FILO(J,K)
           END DO
           HAM(I,J)=HAM(I,J)+SCIS(I,J)
        END DO
           PRINT 8383,ISPN,I,DELTA(I,ISPN),HAM(I,I)-HSSA
 8383      format('scissor1:',2i3,3f12.4)
        END DO

       RETURN
       END
