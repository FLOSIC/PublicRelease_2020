C UTEP Electronic Structure Lab (2020)
       SUBROUTINE CHARGES         
C MARK R. PEDERSON 6 FEBRUARY 2001 
C LOWDIN CHARGES .... EXAMPLE OUTPUT FOR CARBON:
C  1  0  0.00  0.00 1S ATOMIC ORBITAL
C  2  0  0.00  0.00 2S ATOMIC ORBITAL
C  2  1  0.00  0.00 2P ATOMIC ORBITAL
C
       use common2,only : RIDT, IFUIDT, NIDENT, ZNUC, BFCON, N_BARE,
     &   N_CON, LSYMMAX, N_POS, NFNCT
       use common3,only : RMAT
       use common8,only : REP, NDMREP, N_SALC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:36 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IADD, IBAS, IBEG, ICON, IDUM, IEND, IFNCT, IID,
     & IND_SALC, IPOINT, IPOS, IREP, ISALC, ISHELL, JBAS, KBAS, KREP,
     & KSALC, L, LNDX, MNUC, MSITES, NADD, NBAS, NBSMAX, NDEG, NORBS,
     & NSITE
       REAL*8 :: AI , AJ, DEBUG, FACT, RNUC, RNUCI, RNUCJ, SHFT, SS
       SAVE
       CHARACTER*7 FNAME
       CHARACTER*1 NAM(4)
       LOGICAL EXIST,SCISS 
       DIMENSION NDEG(3),AI(3),AJ(3),IPOINT(MAX_REP)
       DIMENSION SS(10,10)
       DIMENSION RNUCI(3,MX_GRP),RNUCJ(3,MX_GRP)
       DIMENSION IBEG(3),IEND(3)
       DIMENSION IND_SALC(ISMAX,MAX_CON,3,2)
       DIMENSION SHFT(MAX_CON,3,2,MAX_FUSET)
       DIMENSION LNDX(6,MAX_CON,3,2)
       DIMENSION RNUC(3,MX_GRP)
       DATA NDEG/1,3,6/
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA NAM/'S','P','D','F'/
       DATA SCISS/.FALSE./
        DO 10 IFNCT=1,NFNCT
         NSITE=0
         IF(.NOT.EXIST)THEN
         WRITE(70,*)IFNCT,ZNUC(IFNCT),' FUNCTION SET'
         ELSE
         READ(70,*)IDUM
         END IF
         DO IPOS=1,N_POS(IFNCT)
          IID=IID+1
          CALL GASITES(1,RIDT(1,IID),MNUC,RNUC,MSITES)
          NSITE=NSITE+MNUC
          IF(.NOT.EXIST)THEN
          WRITE(70,1010) IFNCT,IPOS,(RIDT(I,IID),I=1,3)
          ELSE
          READ(70,*)IDUM
          END IF
 1010     FORMAT('    ',I3,' IP:',I3,' POSITION:',3F15.6)
         END DO
         NADD=N_CON(1,IFNCT)+3*N_CON(2,IFNCT)+6*N_CON(3,IFNCT)
         NORBS=NORBS+NADD*NSITE
         NBSMAX=MAX(NBSMAX,NADD)
         DO L=0,LSYMMAX(IFNCT)
          DO ICON=1,N_CON(L+1,IFNCT)
          IADD=0
           DO I=1,N_BARE(IFNCT)
           IF(ABS(BFCON(I,ICON,L+1,IFNCT)).GT.1D-5) IADD=IADD+1
           END DO
          IF(IADD.GT.1)THEN
                 KEEP(L+1,ICON,IFNCT)=1  
          ELSE
                 KEEP(L+1,ICON,INFCT)=0
          END IF
          END DO
         END DO
   10   CONTINUE
        NBAS=0.0D0
        DO 20 IID=1,NIDENT
         IF(DEBUG) write(6,*)'IID:',IID
         IFNCT=IFUIDT(IID)
         CALL OBINFO(1,RIDT(1,IID),RNUC,MNUC,ISHELL)
         CALL GSMAT(ISHELL,2)
         KSALC=0
         DO KREP=1,IREP  
          KSALC=KSALC+NDMREP(KREP)
         END DO
          DO L=0,LSYMMAX(IFNCT)
           DO ICON=1,N_CON(L+1,IFNCT)
            DO ISALC=1,N_SALC(KSALC,L+1,ISHELL)
            NBAS=NBAS+1                   
            KPR(NBAS)=KEEP(L+1,ICON,IFNCT)                     
            END DO
           END DO
          END DO
   20   CONTINUE
                   DO IBAS=1   ,NBAS
                   DO JBAS=IBAS,NBAS
                   OVER(JBAS,IBAS)=OVER(JBAS,IBAS)
                   END DO                                 
                   END DO
       DO IBAS=1,NBAS
        IF(ABS(EVAL(IBAS)).GE.0.000001)THEN
        DO JBAS=1   ,NBAS
         FACT=EVAL(IBAS)*OVER(IBAS,JBAS)
         DO KBAS=JBAS,NBAS
          HAM(KBAS,JBAS)=HAM(KBAS,JBAS)+FACT*OVER(KBAS,IBAS)
         END DO
        END DO
        END IF
       END DO
       CLOSE(70)
       RETURN
       END
