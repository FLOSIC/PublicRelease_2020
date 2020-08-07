C UTEP Electronic Structure Lab (2020)
       SUBROUTINE WANNIER
C MARK R. PEDERSON 20 NOVEMBER 2000
C
C *********************************************************************
C CALCULATES WANNIER FUNCTIONS - AO BASED
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
C
       use debug1
       use common2,only : RIDT, RCNT, IFUIDT, NIDENT, ZNUC, BFCON,
     &   N_BARE, N_CON, LSYMMAX, N_POS, NFNCT, IGGA, ISPN, NSPN, SYMATM
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI
       use common8,only : REP, N_REP, NDMREP, N_SALC, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IADD, IBAS, IBEG, IBLW, ICON, IDUM, IEND, IF,
     & IFNCT, IID, ILOW, IND, IND_SALC, IPOINT, IPOS, IREP, ISALC,
     & ISHELL, ITER, IWF, J, JBAS, JDUM, JLOW, JSPN, K, KBAS, KLOW,
     & KREP, KSALC, KVRT, L, L1, LBAS, LNDX, LTRAN, MNUC, MODE, MSITES,
     & NADD, NBAS, NBSMAX, NDEG, NLOW, NORBS, NSITE
       REAL*8 :: AI , AJ, CGT, CHRGL, DERIV, DT, EPS, ERROR, ESIC,
     & ESICATM, EVAL, FN, RNUC, RNUCI, RNUCJ, SC1, SHFT, SS, ST, TMP1,
     & TMP2, TOTC, X, XX, Y
       SAVE
       CHARACTER*7 FNAME
       CHARACTER*1 NAM(4)
       CHARACTER*1 LORB
       CHARACTER*2 UD(2)
       CHARACTER*12 QFILE
       CHARACTER*3 ATNAME
       LOGICAL EXIST,SICON,SLOW,FIRST
       DIMENSION NDEG(3),AI(3),AJ(3),IPOINT(MAX_REP)
       DIMENSION SS(10,10),X(100),Y(100),TMP1(100),TMP2(100)
       DIMENSION RNUCI(3,MX_GRP),RNUCJ(3,MX_GRP)
       DIMENSION IBEG(3),IEND(3)
       DIMENSION IND_SALC(ISMAX,MAX_CON,3,2)
       DIMENSION SHFT(MAX_CON,3,2,MAX_FUSET)
       DIMENSION LORB(MAX_CON,3,  MAX_FUSET)
       DIMENSION LNDX(6,MAX_CON,3,2)
       DIMENSION RNUC(3,MX_GRP)
C********
       DIMENSION CHRGL(MAX_IDENT,MAX_CON,3,2) !FIX THIS DIMENSION STATEMENT...
       DIMENSION IND(100,3),ST(100,100),IBLW(100)
       COMMON/SICENRGY/ESIC,ESICATM(MAX_FUSET)
       COMMON/TMP1/OVLTAB(MAX_CON,MAX_CON,LDIM),OVLDIA(MAX_CON,LDIM)
     &  ,HKNTAB(MAX_CON,MAX_CON,LDIM),HNLTAB(MAX_CON,MAX_CON,LDIM)
     &  ,BASFU(3,MAX_CON),TABEXP(MAX_BARE)
     &  ,RRAD(MAXRAD),WRAD(MAXRAD),VLOC(MAXRAD),VPOT(2,MAXRAD)
     &  ,VXCSUB(2),VXC(MAXRAD),RHO(5,MAXRAD),RHCL(NCOUL,3),RHOGRD(10)
     &  ,HMAT(MAX_CON,MAX_CON,LDIM),OMAT(MAX_CON,MAX_CON,LDIM)
     &  ,EVALA(MAX_CON,LDIM),SCRV(MAX_CON),RPW(5),ORCNT(3),EXCVEC(4)
     &  ,EVTB(MAX_CON*LDIM),OCCU(MAX_CON*LDIM)
     &  ,AMAT(MATSZ,MATSZ),AVEC(MATSZ),FCOF(MATSZ)
     &  ,NDEG(MAX_CON*LDIM),IVEC(MATSZ),NIGGA(2)

       DATA NDEG/1,3,6/
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA NAM/'S','P','D','F'/
       DATA UD/'UP','DN'/
       DATA SICON/.FALSE./
       DATA SLOW/.FALSE./
       DATA FIRST/.TRUE./
C      FIRST=.TRUE.
       IF(MODE.EQ.2)ESIC=0.0D0
       OPEN(60,FILE='LOWDEN',FORM='UNFORMATTED')
       REWIND(60)
              DO L=1,2                        !spin index
              DO K=1,3
              DO J=1,MAX_CON
              DO I=1,100                      !FIX THIS
              CHRGL(I,J,K,L)=0.0D0
              END DO
              END DO
              END DO
              END DO
                  DO IF=1,MAX_FUSET
                  DO L1=1,3
                  DO ICON=1,MAX_CON
                  LORB(ICON,L1,IF)='N'        
                  END DO
                  END DO
                  END DO
       INQUIRE(FILE='WNF',EXIST=EXIST)
       IF(.NOT.EXIST)RETURN
C MODE=1 GENERATE WANNIER FUNCTIONS FROM ATOMIC ORBITALS.
C MODE=2 CALCULATE LOCAL CHARGES.
C MODE=3 CONSTRUCT NONLOCAL SIC POTENTIAL.
          OPEN(70,FILE='WNF')
          READ(70,*,END=430)LTRAN
          IF(.NOT.LTRAN) THEN
            CLOSE(70)
            RETURN
          ENDIF
 420   EXIST=.FALSE.
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
                  LORB(ICON,L+1,IFNCT)='Y'        
          IF(.NOT.EXIST)THEN
          WRITE(70,1020)ICON+L,L,OCCU(ICON+L),ICON+L,NAM(L+1)
          ELSE
          READ(70,*   )IDUM,JDUM,SHFT(ICON,L+1,1,IFNCT)
     &                          ,SHFT(ICON,L+1,2,IFNCT)
          END IF
 1020     FORMAT(' ',2I3,'  0.00  0.00 ',I1,A1,' ATOMIC ORBITAL')
          ELSE
          IF(.NOT.EXIST)THEN
          WRITE(70,1021)ICON+L,L, OCCU(ICON+L),       NAM(L+1)
          ELSE
          READ(70,*   )IDUM,JDUM,SHFT(ICON,L+1,1,IFNCT)
     &                          ,SHFT(ICON,L+1,2,IFNCT)
          END IF
 1021     FORMAT(' ',2I3,'  0.00  0.00  ', F10.5,2X,  A1,'(BARE GAUSSIAN)')
          END IF
          END DO
         END DO
   10   CONTINUE


       TOTC=0.0D0
       CALL OVERLAP(1)
       KVRT=0
       DO 2000 JSPN=1,NSPN
       KIND=0
       DO 2000 IREP=1,N_REP
       NBAS=NS_TOT(IREP)
       DO I=1,NBAS            
         DO J=I,NBAS            
          KIND=KIND+1
          OVER(J,I)=HSTOR(KIND,1)
         END DO
        END DO     
        
                         DO IBAS=1,NBAS            
                         DO JBAS=IBAS,NBAS              
                         OVER(IBAS,JBAS)=OVER(JBAS,IBAS)
                         END DO
                         END DO
        DO I=1,3
        PRINT 53,(OVER(I,J),J=4,9)
        END DO
 53     FORMAT('SP OVL:',6G15.6)
                         DO IBAS=1,NBAS
                         DO JBAS=1,NBAS
                         HAM(JBAS,IBAS)=0.0D0
                         END DO
                         END DO
        NLOW=0
        NBAS=0
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
                 IF(LORB(ICON,L+1,IFNCT).EQ.'Y')THEN
                 NLOW=NLOW+1
                 IBLW(NLOW)=NBAS
                 IND(NLOW,1)=IID
                 IND(NLOW,2)=ICON
                 IND(NLOW,3)=L+1
                 HAM(NBAS,NLOW)=1.0D0
                 END IF        
          END DO
          END DO
          END DO
 1030   FORMAT(' ',I5,4F12.4)
   20   CONTINUE
C       write(6,*)'REP:',IREP,' NBAS:',NBAS,' NLOW:',NLOW
        IF(FIRST)THEN
        IF(SLOW)THEN
C BEGIN LOWDIN ITERATIVE LOWDIN ORTHNORMALIZATION:
        ITER=0
        DO ITER=1,8
        ERROR=0.0D0
        DO ILOW=1   ,NLOW
          DT=0.0D0
          DO KBAS=1,NBAS
          DO JBAS=1,NBAS
          DT=DT+HAM(JBAS,ILOW)*HAM(KBAS,ILOW)*OVER(JBAS,KBAS)
          END DO
          END DO
        ERROR=ERROR+ABS(1.0D0-DT)
          DO KBAS=1,NBAS
          HAM(KBAS,ILOW)=HAM(KBAS,ILOW)/SQRT(DT)
          END DO
        END DO
C THIS NEEDS TO BE REWRITTEN AS AN N**3 ALGORITHM...
        DO ILOW=1   ,NLOW
        DO JLOW=ILOW,NLOW
             DT=0.0D0
             DO IBAS=1,NBAS
             DO JBAS=1,NBAS
             DT=DT+HAM(IBAS,ILOW)*HAM(JBAS,JLOW)*
     &               OVER(JBAS,IBAS)
             END DO
             END DO
        ST(JLOW,ILOW)=DT                        
        ST(ILOW,JLOW)=DT
        END DO
        PRINT 64,(ST(JLOW,ILOW),JLOW=1,NLOW)
        END DO
 64     FORMAT(8G15.6,' OVLP')
                   DO ILOW=1,NLOW
                   DO JBAS=1,NBAS
                   HAM(JBAS,ILOW+NLOW)=HAM(JBAS,ILOW)*1.5
                   END DO
                   END DO
        DO ILOW=1,NLOW
        DO JLOW=1,NLOW
C       IF(JLOW.NE.ILOW)THEN
        DT=ST(JLOW,ILOW)*0.5D0
        IF(ILOW.NE.JLOW)THEN
           ERROR=ERROR+ABS(DT)
        ELSE
           ERROR=ERROR+ABS(DT-0.5D0)
        END IF
           DO JBAS=1,NBAS
           HAM(JBAS,ILOW+NLOW)=HAM(JBAS,ILOW+NLOW) 
     &                        -HAM(JBAS,JLOW     )*DT
           END DO
C       END IF
        END DO
        END DO
                   DO ILOW=1,NLOW
                   DO JBAS=1,NBAS
                   HAM(JBAS,ILOW)=HAM(JBAS,ILOW+NLOW)
                   END DO
                   END DO
        write(6,*)'REP,ITER,ERROR:',IREP,ITER,ERROR
C              DO ILOW=1,NLOW
C              write(6,*)'ILOW:',ILOW
C              PRINT 23,(HAM(JBAS,ILOW),JBAS=1,NBAS)
C              END DO
        END DO
       ELSE
C EXACT LOWDIN...
        OPEN(68,FILE='SAVOV',FORM='UNFORMATTED')
        WRITE(68)((OVER(J,I),J=1,NBAS),I=1,NBAS)
        REWIND(68)
            DO I=1,NLOW
            DO J=I,NLOW
            OVER(J,I)=OVER(IBLW(J),IBLW(I))
            END DO
            END DO
            DO I=1,NLOW
            DO J=I,NLOW
            OVER(I,J)=OVER(J,I)
            END DO
            END DO
C USE N**3 version of LOWDEN:
        CALL LOWDEN(NDH,NLOW,OVER,HAM,EVAL,SC1)
          DO I=1,NLOW
            DO J=1,NLOW
            EVAL(J)=HAM(J,I)
            END DO
            DO J=1,NBAS
            HAM(J,I)=0.0D0
            END DO
            DO J=1,NLOW
            HAM(IBLW(J),I)=EVAL(J)
            END DO
          END DO
        READ(68)((OVER(J,I),J=1,NBAS),I=1,NBAS)
        CLOSE(68,STATUS='DELETE')
c           IF(MODE.EQ.3) THEN
c              DO ILOW=1,NLOW
c              WRITE(43,*)'ILOW:',ILOW
c              WRITE(43, 23)(HAM(JBAS,ILOW),JBAS=1,NBAS)
c              END DO
c           END IF
       END IF
C CALCULATE (JBAS|ILOW):
       DO ILOW=1,NLOW
       DO JBAS=1,NBAS
       EVAL(JBAS)=0.0D0
       DO IBAS=1,NBAS
       EVAL(JBAS)=EVAL(JBAS)+HAM(IBAS,ILOW)*OVER(IBAS,JBAS)
       END DO
       END DO
       DO JBAS=1,NBAS
       HAM(JBAS,ILOW)=EVAL(JBAS)
       END DO 
       END DO
       write(6,*)'WRITING REP:',IREP
       WRITE(60)NBAS,NLOW
       WRITE(60)((HAM(JBAS,ILOW),JBAS=1,NBAS),ILOW=1,NLOW)
       WRITE(60)((IND(KLOW,J),KLOW=1,NLOW),J=1,3)
       ELSE
       write(6,*)'READING REP:',IREP
       READ(60)NBAS,NLOW
       READ(60)((HAM(JBAS,ILOW),JBAS=1,NBAS),ILOW=1,NLOW)
       READ(60)((IND(KLOW,J),KLOW=1,NLOW),J=1,3)
       END IF
     
C UPDATE TOTAL CHARGES PROJECTED ONTO LOWDEN STATES:
              DO IWF =1,N_OCC(IREP,JSPN)
              KVRT=KVRT+1
              DO ILOW=1,NLOW
                 DT=0.0D0
                 DO IBAS=1,NBAS
                 DT=DT+PSI_COEF(IBAS,IWF,IREP,JSPN)*HAM(IBAS,ILOW)
                 END DO
              CHRGL(IND(ILOW,1),IND(ILOW,2),IND(ILOW,3),JSPN)=
     &        CHRGL(IND(ILOW,1),IND(ILOW,2),IND(ILOW,3),JSPN)+
     &                  NDMREP(IREP)*DT*DT*OCCUPANCY(KVRT)
              TOTC=TOTC+NDMREP(IREP)*DT*DT*OCCUPANCY(KVRT)
              END DO
              END DO
       write(6,*)'IREP,TOTC:',IREP,TOTC
 2000  CONTINUE
       write(6,*)'FINAL CHARGE:',TOTC
       write(6,*)'LOCAL CHARGES:'
        DO IID=1,NIDENT
        IFNCT=IFUIDT(IID)
        CALL GASITES(1,RIDT(1,IID),MNUC,RNUC,MSITES)
        ESIC=ESICATM(IFNCT)*MNUC+ESIC
        DO L=0,LSYMMAX(IFNCT)
        DO ICON=1,N_CON(L+1,IFNCT)
             DO JSPN=1,NSPN 
             CHRGL(IID,ICON,L+1,JSPN)=CHRGL(IID,ICON,L+1,JSPN)/MNUC
             END DO
           CGT=CHRGL(IID,ICON,L+1,1)+CHRGL(IID,ICON,L+1,NSPN)
        IF(CGT.GT.0.0001)THEN
            L1=L+1
            ATNAME=SYMATM(5,IID)//SYMATM(6,IID)//SYMATM(7,IID)
            WRITE(QFILE,'(A,I1.1,A,A3)')'EVSQ',ICON,NAM(L1),ATNAME
            PRINT *,QFILE
            XX=CGT
            CALL SICINTERPOL(XX,QFILE,FN,DERIV)
            ESIC=ESIC+FN*MNUC
            write(6,*) 'XX', XX, FN, ESIC
        PRINT 43,IID,ICON+L,NAM(L+1),(UD(JSPN),
     &     CHRGL(IID,ICON,L+1,JSPN),JSPN=1,NSPN)
        END IF
          
        END DO
        END DO
        END DO
       IF(.NOT.SICON)ESIC=0.0D0
       IF(MODE.EQ.3.AND.SICON) THEN
       REWIND(60)
       LBAS=0
       DO 3000 IREP=1,N_REP
       write(6,*)'READING REP:',IREP
       READ(60)NBAS,NLOW
       READ(60)((HAM(JBAS,ILOW),JBAS=1,NBAS),ILOW=1,NLOW)
       READ(60)((IND(KLOW,J),KLOW=1,NLOW),J=1,3)
        DO IBAS=1,NBAS
         DO JBAS=IBAS,NBAS
          LBAS=LBAS+1
           DT=0.0D0
           DO ILOW=1,NLOW
            IID =IND(ILOW,1)
            ICON=IND(ILOW,2)
            L1  =IND(ILOW,3)
            ATNAME=SYMATM(5,IID)//SYMATM(6,IID)//SYMATM(7,IID)
            WRITE(QFILE,'(A,I1.1,A,A3)')'EVSQ',ICON,NAM(L1),ATNAME
            PRINT *,QFILE
            XX=CHRGL(IID,ICON,L1,ISPN)
            CALL SICINTERPOL(XX,QFILE,FN,DERIV)
C           EPS=F(XX)
            EPS=DERIV
            XX=CHRGL(IID,ICON,L1,ISPN)
            WRITE(6,*)'EPSILON=', EPS,ICON,L1,XX
c                     OPEN(90,FILE=QFILE)
c                     NNN=1
c 101                 READ(90,*,END=201)X(NNN),DUM,DUM,Y(NNN)
c                     wRITE(6,*)X(NNN),Y(NNN),NNN
c                     NNN=NNN+1
c                     GO TO 101
c 201                 CLOSE(90)
c                     NNN=NNN-1
c                     jjj=0
c                     DO III=NNN,1,-1
c                      jjj=jjj+1
c                      tmp1(jjj)= X(III)
c                      tmp2(jjj)= Y(III)
c                     END DO
c                     DO III=1,NNN
c                       X(III)=TMP1(III)
c                       Y(III)=TMP2(III)
c                       write(6,*)X(III),Y(III)
c                     END DO
c
c            IFD1=0
c            IFDN=0
c            CALL PCIN(X,Y,FD,B,C,11,1,IFD1,IFDN)
c            XX=CHRGL(IID,ICON,L1,ISPN)
C            IF(XX.GT.1.0D0) XX=1.0D0
c            CALL EVLPC(X,Y,FD,B,C,11,XX,FF,FF1,IL)
c            EPSILON=FF
c            WRITE(6,*)'EPSILON=', EPS,ICON,L1,XX
             DT=DT+ HAM(IBAS,ILOW)*HAM(JBAS,ILOW)*EPS    
           END DO
            HSTOR(LBAS,2)=HSTOR(LBAS,2)+DT
         END DO
        END DO
 3000  CONTINUE
       END IF
 23    FORMAT(' ',6G15.6)
 43    FORMAT('QCALC, IID:',I2,' ',I1,A1,' ',2(A2,' ',G15.6))
       CLOSE(70)
       CLOSE(60)
       write(6,*)'LOWDEN CLOSED'
       FIRST=.FALSE.
       RETURN    
       END
