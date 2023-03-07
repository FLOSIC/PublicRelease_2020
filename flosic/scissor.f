C UTEP Electronic Structure Lab (2020)
       SUBROUTINE SCISSOR(IREP)          
C MARK R. PEDERSON 20 NOVEMBER 2000
C
C NOTE: FORCES NOT NECESSARILY CORRECT....
C OK FOR SINGLE GEOMETRIES....
C *********************************************************************
C
C THIS PROGRAM CREATES "SCISSOR" WHICH LISTS EACH BASIS FUNCTION 
C ON EACH ATOM AS SOMETHING LIKE...
C
C
C1  6. FUNCTION SET
C     1 IP:  1 POSITION:       0.000000       0.000000       0.000000
C  1  0  0.00  0.00 1S ATOMIC ORBITAL
C  2  0  0.00  0.00 2S ATOMIC ORBITAL
C  3  0  0.00  0.00  S(BARE GAUSSIAN)
C  4  0  0.00  0.00  S(BARE GAUSSIAN)
C  5  0  0.00  0.00  S(BARE GAUSSIAN)
C  2  1  0.00  0.00 2P ATOMIC ORBITAL
C  3  1  0.00  0.00  P(BARE GAUSSIAN)
C  4  1  0.00  0.00  P(BARE GAUSSIAN)
C  5  1  0.00  0.00  P(BARE GAUSSIAN)
C  3  2  0.00  0.00  D(BARE GAUSSIAN)
C  4  2  0.00  0.00  D(BARE GAUSSIAN)
C  5  2  0.00  0.00  D(BARE GAUSSIAN)                                           
C
C
C TO ADD A SCISSOR SHIFT OF 27.2 eV TO, FOR EXAMPLE THE CARBON 1S UP
C  ORBITAL WE CHANGE:
C
C  1  0  0.00  0.00 1S ATOMIC ORBITAL
C TO:
C  1  0 -1.00  0.00 1S ATOMIC ORBITAL
C
C
       use debug1
       use for_diag1 !added
       use common2,only : RIDT, IFUIDT, NIDENT, ZNUC, BFCON, N_BARE,
     &   N_CON, LSYMMAX, N_POS, NFNCT, ISPN
       use common3,only : RMAT
       use common8,only : REP, NDMREP, N_SALC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:02 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IREP, I, IADD, IBAS, IBEG, ICON, IDUM, IEND, IFNCT,
     & IID, IND_SALC, IPOINT, IPOS, ISALC, ISHELL, JBAS, JDUM, KBAS,
     & KREP, KSALC, L, LNDX, MNUC, MSITES, NADD, NBAS, NBSMAX, NDEG,
     & NORBS, NSITE
       REAL*8 :: AI , AJ, FACT, RNUC, RNUCI, RNUCJ, SHFT, SS
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
       INQUIRE(FILE='SCISSOR',EXIST=EXIST)
       IF(.NOT.EXIST)RETURN
          OPEN(70,FILE='SCISSOR')
          READ(70,*,END=420)SCISS
          IF(.NOT.SCISS) THEN
            CLOSE(70)
            RETURN
           ENDIF
          GOTO 430
 420   CONTINUE
       EXIST=.FALSE.
       IF(.NOT.EXIST) WRITE(70,*)SCISS,' IF SCISSOR=T PLEASE RUN WITH
     & PURGRSQ= T !!!'
 430   CONTINUE
        NORBS=0
        NBSMAX=0
        IID=0
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
          IF(.NOT.EXIST)THEN
          WRITE(70,1020)ICON+L,L,ICON+L,NAM(L+1)
          ELSE
          READ(70,*   )IDUM,JDUM,SHFT(ICON,L+1,1,IFNCT)
     &                          ,SHFT(ICON,L+1,2,IFNCT)
          END IF
 1020     FORMAT(' ',2I3,'  0.00  0.00 ',I1,A1,' ATOMIC ORBITAL')
          ELSE
          IF(.NOT.EXIST)THEN
          WRITE(70,1021)ICON+L,L,       NAM(L+1)
          ELSE
          READ(70,*   )IDUM,JDUM,SHFT(ICON,L+1,1,IFNCT)
     &                          ,SHFT(ICON,L+1,2,IFNCT)
          END IF
 1021     FORMAT(' ',2I3,'  0.00  0.00  ',   A1,'(BARE GAUSSIAN)')
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
C NOTE THEN NEXT TWO DO LOOPS MAY NEED TO BE PERMUTED.
          DO ICON=1,N_CON(L+1,IFNCT)
          DO ISALC=1,N_SALC(KSALC,L+1,ISHELL)
          NBAS=NBAS+1                   
          AEVAL(NBAS)=SHFT(ICON,L+1,ISPN,IFNCT)/AOVER(NBAS,NBAS)
          IF(DEBUG) THEN
          PRINT 1030,NBAS,AHAM(NBAS,NBAS),AOVER(NBAS,NBAS),
     &             SHFT(ICON,L+1,ISPN,IFNCT),AEVAL(NBAS)
          ENDIF
          END DO
          END DO
          END DO
 1030   FORMAT(' ',I5,4F12.4)
   20   CONTINUE
                         DO IBAS=1,NBAS
                         DO JBAS=IBAS,NBAS
                         AOVER(IBAS,JBAS)=AOVER(JBAS,IBAS)
                         END DO
                         END DO
       DO IBAS=1,NBAS
        IF(ABS(AEVAL(IBAS)).GE.0.000001)THEN
        DO JBAS=1   ,NBAS
         FACT=AEVAL(IBAS)*AOVER(IBAS,JBAS)
         DO KBAS=JBAS,NBAS
          AHAM(KBAS,JBAS)=AHAM(KBAS,JBAS)+FACT*AOVER(KBAS,IBAS)
         END DO
        END DO
        END IF
       END DO
       CLOSE(70)
       RETURN
       END
