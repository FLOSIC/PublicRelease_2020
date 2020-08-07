C UTEP Electronic Structure Lab (2020)
C/
C ******************************************************************\
C
       SUBROUTINE GETBAS(LSETUP,MPTS,IGP,R,PSIBR,JREP,NBAS)
C
C CALCULATES THE VALUE OF ALL SALC BASIS FUNCTIONS FOR A
C BUNCH OF MESH POINTS AND STORES THEM IN PSIBR
C WRITTEN BY MARK R PEDERSON
C DVP: THIS ROUTINE MAY NOT BE CALLED WITH ARBITRARY MPTS/R/IGP/JREP
C ARGUMENTS. FOR A PARTICULAR SET MPTS/R OF POINTS, IT MUST BE CALLED
C WITH JREP=IGP=1 FIRST AND CAN ONLY BE CALLED WITH OTHER JREP/IGP
C ARGUMENTS IF MPTS/R DOES NOT CHANGE.
C
       use debug1
       use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT
       use common5,only : PSI
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MPTS, IGP, JREP, NBAS, I_POS, IBEG, IC, IEND, IFNCT,
     & ILOC, INS, IPTS, IQ, IR, IREP, IROW, IS, ISALC, ISHELLA, ISS,
     & JR, K_POS, KSHELL, LI, M_NUC, MU, NBI, NBI_BEG, NDEG, NSC
       REAL*8 :: R , PSIBR, FAC, T1
       SAVE
       LOGICAL LSETUP
       DIMENSION PSIBR(MPBLOCK,5*NDH),R(3,MPBLOCK)
C       COMMON/TMP1/ORBT(MPBLOCK,10,MAX_CON,3)
C     &  ,PSIBS(MPBLOCK,MAXSLC*MAX_IDENT)
C     &  ,CONORB(MPBLOCK,MAX_CON,3),ANGLE(MPBLOCK,10),DTS(3,MPBLOCK)
C     &  ,RVEC(3),RVECA(3,MX_GRP),NSALC(MAX_REP)
C     &  ,INDEX(ISMAX,MAX_CON,3,MAX_IDENT,MAXSYMSALC)
       REAL*8 :: ORBT(MPBLOCK,10,MAX_CON,3),
     &    PSIBS(MPBLOCK,MAXSLC*MAX_IDENT),
     &    CONORB(MPBLOCK,MAX_CON,3),ANGLE(MPBLOCK,10),
     *    DTS(3,MPBLOCK),RVECA(3,MX_GRP),RVEC(3)
       INTEGER :: NSALC(MAX_REP),
     &    INDEX(ISMAX,MAX_CON,3,MAX_IDENT,MAXSYMSALC)
       DIMENSION IBEG(3),IEND(3),NDEG(3)
       DATA IBEG,IEND,NDEG/1,2,5,1,4,10,1,3,6/
C
C IF (LSETUP), SET UP INDEX ARRAY FOR SALCS
C
       IF (LSETUP) THEN
        INS=0
        ISALC=0
        DO IREP=1,N_REP
         DO IROW=1,NDMREP(IREP)
          ISALC=ISALC+1
          ISHELLA=0
          NSC=0
          DO IFNCT=1,NFNCT
           DO I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            RVEC(1)=RIDT(1,ISHELLA)
            RVEC(2)=RIDT(2,ISHELLA)
            RVEC(3)=RIDT(3,ISHELLA)
            CALL OBINFO(1,RVEC,RVECA,M_NUC,KSHELL)
            CALL GSMAT(KSHELL,1)
            DO LI=1,LSYMMAX(IFNCT)+1
             DO IC=1,N_CON(LI,IFNCT)
              DO IQ=1,N_SALC(ISALC,LI,KSHELL)
               INS=INS+1
               NSC=NSC+1
               INDEX(IQ,IC,LI,ISHELLA,ISALC)=INS
              END DO
             END DO
            END DO
           END DO
          END DO
         END DO
         IF (DEBUG) write(6,*)'IREP,IROW-1,NSC',IREP,IROW-1,NSC
         NSALC(IREP)=NSC
        END DO
        IF (DEBUG) THEN
         write(6,*)'GETBAS: TOTAL NUMBER OF BASIS FUNCTIONS:',INS
        END IF

        RETURN
       END IF
C
       IF (MPTS.GT.MPBLOCK) THEN
        write(6,*)'MPTS MUST BE:',MPBLOCK
        CALL STOPIT
       END IF
C
C INIT PSIBS
C
       IF ((IGP.EQ.1).AND.(JREP.EQ.1)) THEN
        DO IS=1,INS
         DO IPTS=1,MPTS
          PSIBS(IPTS,IS)=0.0D0
         END DO
        END DO
        ISHELLA=0
        CALL GTTIME(T1)
        DO 250 IFNCT=1,NFNCT
         DO 240 I_POS=1,N_POS(IFNCT)
          ISHELLA=ISHELLA+1
          RVEC(1)=RIDT(1,ISHELLA)
          RVEC(2)=RIDT(2,ISHELLA)
          RVEC(3)=RIDT(3,ISHELLA)
          CALL OBINFO(1,RVEC,RVECA,M_NUC,KSHELL)
          CALL GSMAT(KSHELL,1)
          DO 230 K_POS=1,M_NUC
           DO IPTS=1,MPTS
            DTS(1,IPTS)=R(1,IPTS)-RVECA(1,K_POS)
            DTS(2,IPTS)=R(2,IPTS)-RVECA(2,K_POS)
            DTS(3,IPTS)=R(3,IPTS)-RVECA(3,K_POS)
           END DO
           CALL GTORBNH(MPTS,DTS,IFNCT,ANGLE,CONORB)
           DO  LI=1,LSYMMAX(IFNCT)+1
            DO  IC=1,N_CON(LI,IFNCT)
             DO  MU=IBEG(LI),IEND(LI)
              DO  IPTS=1,MPTS
               ORBT(IPTS,MU,IC,LI)=CONORB(IPTS,IC,LI)*ANGLE(IPTS,MU)
              END DO
             END DO
            END DO
           END DO
C
           ISALC=0
           NBI_BEG=0
           DO IREP=1,N_REP
            DO IROW=1,NDMREP(IREP)
             ISALC=ISALC+1
             DO LI=1,LSYMMAX(IFNCT)+1
              ILOC=NDEG(LI)*(K_POS-1)-IBEG(LI)+1
              DO MU=IBEG(LI),IEND(LI)
               NBI=NBI_BEG
               DO IC=1,N_CON(LI,IFNCT)
                DO IQ=1,N_SALC(ISALC,LI,KSHELL)
                 INS=INDEX(IQ,IC,LI,ISHELLA,ISALC)
                 FAC= U_MAT(ILOC+MU,IQ,ISALC,LI,1)
                 NBI=NBI+1
                 DO IPTS=1,MPTS
                  PSIBS(IPTS,INS)=PSIBS(IPTS,INS)
     &            +FAC*ORBT(IPTS,MU,IC,LI)
                 END DO
                END DO
               END DO
              END DO
              NBI_BEG=NBI
             END DO
            END DO
           END DO
           IF (NBI.GT.MAXSLC) THEN
            write(6,*)'NUMHAM: MAXSLC MUST BE AT LEAST: ',NBI
            CALL STOPIT
           END IF
 230      CONTINUE
 240     CONTINUE
 250    CONTINUE
       END IF
C
       INS=0
       DO IREP=1,N_REP
        IF (IREP.EQ.JREP) THEN
         NBAS=NSALC(IREP)
         DO ISS=1,NSALC(IREP)
          INS=INS+1
          DO IR=1,NDMREP(IREP)
           DO IPTS=1,MPTS
            PSIBR(IPTS,ISS+(IR-1)*NSALC(IREP))=0.0D0
           END DO
           DO JR=1,NDMREP(IREP)
            DO IPTS=1,MPTS
             PSIBR(IPTS,ISS+(IR-1)*NSALC(IREP))=
     &       PSIBR(IPTS,ISS+(IR-1)*NSALC(IREP))+
     &       PSIBS(IPTS,INS+(JR-1)*NSALC(IREP))*REP(JR,IR,IGP,IREP)
            END DO
           END DO
          END DO
         END DO
         INS=INS+(NDMREP(IREP)-1)*NSALC(IREP)
        ELSE
         INS=INS+ NDMREP(IREP)   *NSALC(IREP)
        END IF
       END DO

       RETURN
       END
