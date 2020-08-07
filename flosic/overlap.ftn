C UTEP Electronic Structure Lab (2020)
C *********************************************************************
C
       SUBROUTINE OVERLAP(MODE)
C
C  ORIGINALLY WRITTEN BY MARK R PEDERSON (1985), 
C  modified by Arlin Briley,
C  Dirk Porezag and God knows who else
C  CALLS TO REWRITE/INITZE REMOVED BY MRP 4-Feb 1998
C
C  MODE=1   CALCULATE OVERLAPS                -STORE IN HSTOR(I,1)
C  MODE=2   CALCULATE HAMILTONIAN             -STORE IN HSTOR(I,2)
C  MODE=3   CALCULATE KINETIC+NONLOCAL EOP    -STORE IN HSTOR(I,1)
C  MODE=4   CALCULATE KINETIC EOP             -STORE IN HSTOR(I,1)
C  EOP IS ABBREVIATION FOR ENERGY OPERATOR
C
       use hstor1,only : hstor
       use debug1
       use common1,only : ISITPSP
       use common2,only : RIDT, IFUIDT, NIDENT, BFCON, BFALP, N_BARE,
     &   N_CON, LSYMMAX, N_POS, NFNCT
       use common3,only : RMAT
       use common5,only : ISTSCF, HOLD, HTEMP
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC,
     &   INDBEG, NS_TOT
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MODE, I, IALP, IBASE, IBEG, ICON, IEND, IFNCT, II,
     & IID, IJ, IND_SALC, INDEX, IOFS, IPOINT, IPOS, IQ, IREC, IS, ISA,
     & ISHELL, ISIT, ISITE, ISTRT, ITIMES, J, JALP, JB, JBASE, JFNCT,
     & JID, JNDEX, JNDSV, JQ, JS, JSA, JSHELL, JSIT, JSITE, JSTRT,
     & JTEMP, K_REP, K_ROW, KNDEX, KOUNT, KREP, KSALC, L, LI, LJ, LNDX,
     & MAXI, MAXJ, MI, MJ, MNUC, MNUCI, MNUCJ, MODSTOR, MSITES, MTMAX,
     & MUI, MUJ, NADD, NBLOCK, NBSMAX, NDEG, NHTOT, NORBS, NREC, NSITE,
     & NTIMES
       REAL*8 :: AI , AJ, ALPHAI, ALPHAJ, ARG, FF, PROD, RNUC, RNUCI,
     & RNUCJ, SS
       SAVE
       CHARACTER*7 FNAME
       LOGICAL FIRST,OVLBABY,HAMBABY,HNLBABY,KINBABY,STARTUP
       DIMENSION NDEG(3),AI(3),AJ(3),IPOINT(MAX_REP)
       DIMENSION SS(10,10)
       DIMENSION RNUCI(3,MX_GRP),RNUCJ(3,MX_GRP)
       DIMENSION IBEG(3),IEND(3)
       DIMENSION IND_SALC(ISMAX,MAX_CON,3,2)
       DIMENSION LNDX(6,MAX_CON,3,2)
       DIMENSION RNUC(3,MX_GRP)
       DATA NDEG/1,3,6/
       DATA IBEG,IEND/1,2,5,1,4,10/
       DATA FIRST/.TRUE./
       DATA OVLBABY/.FALSE./
       DATA HAMBABY/.FALSE./
       DATA HNLBABY/.FALSE./
       DATA KINBABY/.FALSE./
C
C CHECK AND PRINT BASIS SET PROPERTIES
C
       IF (FIRST) THEN
        FIRST=.FALSE.
        NORBS=0
        NBSMAX=0
        WRITE(7,*) ' '
        WRITE(7,*) ' SUMMARY OF BASIS SET:'
        WRITE(7,*) ' NUMBER OF CONTRACTED FUNCTION SETS:',NFNCT
        IID=0
        DO 10 IFNCT=1,NFNCT
         NSITE=0
         WRITE(7,*) ' '
         WRITE(7,*) ' FUNCTION SET:',IFNCT
         DO IPOS=1,N_POS(IFNCT)
          IID=IID+1
          CALL GASITES(1,RIDT(1,IID),MNUC,RNUC,MSITES)
          NSITE=NSITE+MNUC
          WRITE(7,1010) IFNCT,IPOS,(RIDT(I,IID),I=1,3)
 1010     FORMAT(' FS:',I3,' IP:',I3,' POSITION:',3F15.6)
         END DO
         WRITE(7,*) 'CONTRACTED ORBITALS FOR ANGULAR TYPE S:',
     &               N_CON(1,IFNCT)
         WRITE(7,*) 'CONTRACTED ORBITALS FOR ANGULAR TYPE P:',
     &               N_CON(2,IFNCT)
         WRITE(7,*) 'CONTRACTED ORBITALS FOR ANGULAR TYPE D:',
     &               N_CON(3,IFNCT)
         WRITE(7,*) 'GAUSSIANS FOR FUNCTION SET:'
         WRITE(7,1020)(BFALP(I,IFNCT),I=1,N_BARE(IFNCT))
 1020    FORMAT(' ',7G11.3)
         NADD=N_CON(1,IFNCT)+3*N_CON(2,IFNCT)+6*N_CON(3,IFNCT)
         NORBS=NORBS+NADD*NSITE
         NBSMAX=MAX(NBSMAX,NADD)
         DO L=0,LSYMMAX(IFNCT)
          WRITE(7,*) 'CONTRACTED ORBITALS FOR ANGULAR MOMENTUM:',L
          DO ICON=1,N_CON(L+1,IFNCT)
           WRITE(7,1020)(BFCON(I,ICON,L+1,IFNCT),I=1,N_BARE(IFNCT))
          END DO
         END DO
   10   CONTINUE
        WRITE(7,*) 'TOTAL NUMBER OF CONTRACTED ORBITALS:',NORBS
        WRITE(7,*) ' '
        IF (IID.NE.NIDENT) THEN
         write(6,*)'OVERLAP: IID AND NIDENT ARE DIFFERENT: ',IID,NIDENT
         CALL STOPIT
        END IF
        IF (NBSMAX.GT.MAXUNSYM) THEN
         write(6,*)'OVERLAP: MAXUNSYM MUST BE AT LEAST: ',NBSMAX
         CALL STOPIT
        END IF
C
C INITITIALIZE SALC INDICES
C
        DO KREP=1,N_REP
         INDBEG(1,KREP)=0
         NS_TOT(KREP)  =0
        END DO
        DO 20 IID=1,NIDENT
         IFNCT=IFUIDT(IID)
         CALL OBINFO(1,RIDT(1,IID),RNUC,MNUC,ISHELL)
         CALL GSMAT(ISHELL,2)
         KSALC=0
         DO KREP=1,N_REP
          KSALC=KSALC+NDMREP(KREP)
          DO L=0,LSYMMAX(IFNCT)
           NS_TOT(KREP)=NS_TOT(KREP)
     &                 +N_CON(L+1,IFNCT)*N_SALC(KSALC,L+1,ISHELL)
          END DO
          IF (IID .NE. NIDENT) INDBEG(IID+1,KREP)=NS_TOT(KREP)
         END DO
   20   CONTINUE
C
C SQUISHING TO SINGLE INDICES
C
        NHTOT=0
        DO KREP=1,N_REP
         IPOINT(KREP)=NHTOT
         NHTOT=NHTOT+(NS_TOT(KREP)*(NS_TOT(KREP)+1))/2
         IF (NS_TOT(KREP).GT.NDH) THEN
          write(6,*)'OVERLAP: NDH MUST BE AT LEAST: ',NS_TOT(KREP)
          CALL STOPIT
         END IF
        END DO
        IF (NHTOT.GT.NDH_TOT) THEN
         write(6,*)'OVERLAP: NDH_TOT MUST BE AT LEAST: ',NHTOT
         CALL STOPIT
        END IF
C
C CHECK IF MTEMP IS LARGE ENOUGH
C
        MTMAX=0
        DO 40 IID=1,NIDENT
         DO JID=IID,NIDENT
          KOUNT=0
          DO KREP=1,N_REP
           IF (IID.EQ.NIDENT) THEN
            MI=NS_TOT(KREP)-INDBEG(IID,KREP)
           ELSE
            MI=INDBEG(IID+1,KREP)-INDBEG(IID,KREP)
           END IF
           IF (JID.EQ.NIDENT) THEN
            MJ=NS_TOT(KREP)-INDBEG(JID,KREP)
           ELSE
            MJ=INDBEG(JID+1,KREP)-INDBEG(JID,KREP)
           END IF
           KOUNT=KOUNT+MI*MJ
          END DO
          MTMAX=MAX(KOUNT,MTMAX)
         END DO
   40   CONTINUE
        IF (MTMAX.GT.MTEMP) THEN
         write(6,*)'OVERLAP: MTEMP MUST BE AT LEAST: ',MTMAX
         CALL STOPIT
        END IF
C
C PRINT DATA TO FILE OUTPUT
C
        DO KREP=1,N_REP
         WRITE(7,*) 'REPRESENTATION:',KREP,' HAS:',NS_TOT(KREP),
     &              ' BASES'
        END DO
        IF (DEBUG) write(6,*)'DONE WITH OVERLAP PRELIMINARIES'
       END IF
C
C END OF FIRST-TIME-ONLY STUFF
C NOW, DETERMINE WHETHER VSTART SHOULD BE CALLED AND DEFINE MODSTOR
C
       STARTUP=.FALSE.
       IF (((ISTSCF.EQ.0).OR.(ISTSCF.EQ.3)).AND.(MODE.EQ.2)) THEN
        STARTUP=.TRUE.
       END IF
       MODSTOR=1
       IF (MODE.EQ.2) MODSTOR=2
C
C DETERMINE NHTOT (TOTAL NUMBER OF NONZERO HAMILTONIAN MATRIX ELEMENTS)
C CHECK WHETHER DATA CAN BE READ FROM FILE
C
       NHTOT=0
       DO KREP=1,N_REP
        NHTOT=NHTOT+(NS_TOT(KREP)*(NS_TOT(KREP)+1))/2
       END DO
       IF ((MODE.EQ.1).AND.OVLBABY) THEN
        FNAME='OVLBABY'
        GOTO 100
       END IF
       IF (((MODE.EQ.2).OR.(MODE.EQ.3)).AND.HAMBABY) THEN
        FNAME='HAMBABY'
        GOTO 100
       END IF
       IF ((MODE.EQ.4).AND.KINBABY) THEN
        FNAME='KINBABY'
        GOTO 100
       ENDIF
       GOTO 150
C
C DATA IS AVAILABLE, READ IT
C
  100  CONTINUE
       OPEN(66,FILE=FNAME,FORM='UNFORMATTED',STATUS='OLD')
       REWIND(66)
       READ(66,END=110) NREC
       IF (NREC.NE.NHTOT) THEN
        write(6,*)'OVERLAP: READ: WRONG NUMBER OF DATA: ',NREC,NHTOT
        CALL STOPIT
       ENDIF
       IF (size(HSTOR,1).LT.NREC) THEN
        write(6,*)'OVERLAP: wrong size of HSTOR'
        write(6,*)'overlap:Number of records is',NREC
        write(6,*)'overlap:size of HSTOR is',size(HSTOR,1)
        CALL STOPIT
       ENDIF
       READ(66,END=110,ERR=110)(HSTOR(IREC,MODSTOR),IREC=1,NREC)
       CLOSE(66)
       GOTO 800
C
C READ ERROR
C
  110  write(6,*)'OVERLAP: COULD NOT READ STORED DATA'
       write(6,*)'TRIED TO READ FILE: ',FNAME
       CALL STOPIT
C
C END READING FROM FILES
C CALCULATE NONLOCAL PART OF HAMILTONIAN AND WRITE IT TO HAMBABY
C
  150  IF ((MODE.EQ.2).OR.(MODE.EQ.3)) THEN
        IF ((ISITPSP.EQ.1).AND.(.NOT.HNLBABY)) THEN
         DO IREC=1,NHTOT
          HSTOR(IREC,MODSTOR)= 0.0D0
         END DO
         CALL VNLHAM(MODSTOR)
         HNLBABY=.TRUE.
         HAMBABY=.FALSE.
         OPEN(66,FILE='HAMBABY',FORM='UNFORMATTED',STATUS='UNKNOWN')
         REWIND(66)
         WRITE(66) NHTOT,MTEMP
         DO IOFS=0,NHTOT-1,MTEMP
          NREC=MIN(MTEMP,NHTOT-IOFS)
          WRITE(66)(HSTOR(IREC+IOFS,MODSTOR), IREC=1,NREC)
         END DO
         CLOSE(66)
        END IF
       END IF
C
C CALCULATE KINETIC ENERGY OR OVERLAP MATRICES
C
       DO IREC=1,NHTOT
        HSTOR(IREC,MODSTOR)= 0.0D0
       END DO
       DO 400 IID=1,NIDENT
        CALL OBINFO(1,RIDT(1,IID),RNUCI,MNUCI,ISHELL)
        CALL GSMAT(ISHELL,1)
        DO 320 JID=IID,NIDENT
         CALL OBINFO(1,RIDT(1,JID),RNUCJ,MNUCJ,JSHELL)
         CALL GSMAT(JSHELL,2)
C
C LOOP OVER ALL SITES OF EACH SHELL
C
         ISA=IID
         JSA=JID
         IFNCT=IFUIDT(IID)
         JFNCT=IFUIDT(JID)
         IS=ISHELL
         JS=JSHELL
         MAXI=N_CON(1,IFNCT)+3*N_CON(2,IFNCT)+6*N_CON(3,IFNCT)
         MAXJ=N_CON(1,JFNCT)+3*N_CON(2,JFNCT)+6*N_CON(3,JFNCT)
         DO JTEMP=1,MTMAX
          HTEMP(JTEMP)=0.0D0
         END DO
         DO 260 ISITE=1 ,MNUCI
          AI(1)=RNUCI(1,ISITE)
          AI(2)=RNUCI(2,ISITE)
          AI(3)=RNUCI(3,ISITE)
          IF (ISA.EQ.JSA) THEN
           JB=ISITE
          ELSE
           JB=1
          END IF
          DO 258 JSITE=JB,MNUCJ
           AJ(1)=RNUCJ(1,JSITE)
           AJ(2)=RNUCJ(2,JSITE)
           AJ(3)=RNUCJ(3,JSITE)
           DO I=1,MAXI
            DO J=1,MAXJ
             HOLD(J,I)=0.0D0
            END DO
           END DO
c
C calculate local indices:
c
           INDEX=0
           DO LI=0,LSYMMAX(IFNCT)
            DO IBASE=1,N_CON(LI+1,IFNCT)
             DO MUI=1,NDEG(LI+1)
              INDEX=INDEX+1
              LNDX(MUI,IBASE,LI+1,1)=INDEX
             END DO
            END DO
           END DO
           INDEX=0
           DO LJ=0,LSYMMAX(JFNCT)
            DO JBASE=1,N_CON(LJ+1,JFNCT)
             DO MUJ=1,NDEG(LJ+1)
              INDEX=INDEX+1
              LNDX(MUJ,JBASE,LJ+1,2)=INDEX
             END DO
            END DO
           END DO
           DO 230 IALP=1,N_BARE(IFNCT)
            ALPHAI=BFALP(IALP,IFNCT)
            DO 220 JALP=1,N_BARE(JFNCT)
             ALPHAJ=BFALP(JALP,JFNCT)
             ARG=(ALPHAI*ALPHAJ/(ALPHAI+ALPHAJ))
     &          *((AI(1)-AJ(1))**2+(AI(2)-AJ(2))**2+(AI(3)-AJ(3))**2)
             IF (ARG .GT. CUTEXP) GOTO 220
             IF (MODE.EQ.1) THEN  ! Overlap
              CALL OVMXSF(ALPHAI,ALPHAJ,AI,AJ,SS)
             ELSE                                ! Kinetic energy
              CALL KNMXSF(ALPHAI,ALPHAJ,AI,AJ,SS)
             END IF
C
C VSTART IS NEEDED FOR A GOOD FIRST GUESS
C
             IF (STARTUP) CALL VSTART(ALPHAI,ALPHAJ,AI,AJ,SS)
             INDEX=0
             DO LI =0,LSYMMAX(IFNCT)
              DO IBASE=1,N_CON(LI+1,IFNCT)
               DO MUI=IBEG(LI+1),IEND(LI+1)
                INDEX=INDEX+1
                JNDEX=0
                DO LJ =0,LSYMMAX(JFNCT)
                 DO JBASE=1,N_CON(LJ+1,JFNCT)
                  PROD=BFCON(IALP,IBASE,LI+1,IFNCT)
     &                *BFCON(JALP,JBASE,LJ+1,JFNCT)
                  DO MUJ=IBEG(LJ+1),IEND(LJ+1)
                   JNDEX=JNDEX+1
                   HOLD(JNDEX,INDEX)=HOLD(JNDEX,INDEX)
     &                              +PROD*SS(MUI,MUJ)
                  END DO
                 END DO
                END DO
               END DO
              END DO
             END DO
  220       CONTINUE
  230      CONTINUE
C
C NOW UPDATE SALC MATRICES FOR EACH REPRESENTATION:
C
           NTIMES=2
           IF (ISA.NE.JSA)     NTIMES=1
           IF (ISITE.EQ.JSITE) NTIMES=1
           DO 256 ITIMES=1,NTIMES
            IF (ITIMES.EQ.1) THEN
             ISIT=ISITE
             JSIT=JSITE
             KSALC=0
             JNDEX=0
             DO 240 K_REP=1,N_REP
              FF=1.0D0/NDMREP(K_REP)
              DO K_ROW=1,NDMREP(K_REP)
               IF (K_ROW.EQ.1) THEN
                JNDSV=JNDEX
               ELSE
                JNDEX=JNDSV
               END IF
               KSALC=KSALC+1
               DO LI =0,LSYMMAX(IFNCT)
                ISTRT=(ISIT-1)*NDEG(LI+1)
                DO LJ =0,LSYMMAX(JFNCT)
                 JSTRT=(JSIT-1)*NDEG(LJ+1)
                 DO IBASE=1,N_CON(LI+1,IFNCT)
                  DO JBASE=1,N_CON(LJ+1,JFNCT)
                   DO MUJ=1,NDEG(LJ+1)
                    DO MUI=1,NDEG(LI+1)
                     KNDEX=JNDEX
                     DO IQ=1,N_SALC(KSALC,LI+1,IS)
                      DO JQ=1,N_SALC(KSALC,LJ+1,JS)
                       KNDEX=KNDEX+1
                       HTEMP(KNDEX)=HTEMP(KNDEX)+
     &                 U_MAT(MUI+ISTRT,IQ,KSALC,LI+1,1)*
     &                 U_MAT(MUJ+JSTRT,JQ,KSALC,LJ+1,2)*FF*
     &                 HOLD(LNDX(MUJ,JBASE,LJ+1,2),
     &                      LNDX(MUI,IBASE,LI+1,1))
                      END DO
                     END DO
                    END DO
                   END DO
                   JNDEX=KNDEX
                  END DO
                 END DO
                END DO
               END DO
              END DO
  240        CONTINUE
            ELSE
             ISIT=ISITE
             JSIT=JSITE
             KSALC=0
             JNDEX=0
             DO 250 K_REP=1,N_REP
              FF=1.0D0/NDMREP(K_REP)
              DO K_ROW=1,NDMREP(K_REP)
               IF (K_ROW.EQ.1) THEN
                JNDSV=JNDEX
               ELSE
                JNDEX=JNDSV
               END IF
               KSALC=KSALC+1
               DO LJ =0,LSYMMAX(JFNCT)
                JSTRT=(JSIT-1)*NDEG(LJ+1)
                DO LI =0,LSYMMAX(IFNCT)
                 ISTRT=(ISIT-1)*NDEG(LI+1)
                 DO JBASE=1,N_CON(LJ+1,JFNCT)
                  DO IBASE=1,N_CON(LI+1,IFNCT)
                   DO MUI=1,NDEG(LI+1)
                    DO MUJ=1,NDEG(LJ+1)
                     KNDEX=JNDEX
                     DO JQ=1,N_SALC(KSALC,LJ+1,JS)
                      DO IQ=1,N_SALC(KSALC,LI+1,IS)
                       KNDEX=KNDEX+1
                       HTEMP(KNDEX)=HTEMP(KNDEX)+
     &                 U_MAT(MUJ+JSTRT,JQ,KSALC,LJ+1,2)*
     &                 U_MAT(MUI+ISTRT,IQ,KSALC,LI+1,1)*FF*
     &                 HOLD(LNDX(MUJ,JBASE,LJ+1,2),
     &                      LNDX(MUI,IBASE,LI+1,1))
                      END DO
                     END DO
                    END DO
                   END DO
                   JNDEX=KNDEX
                  END DO
                 END DO
                END DO
               END DO
              END DO
  250        CONTINUE
            END IF
  256      CONTINUE
  258     CONTINUE
  260    CONTINUE
C
C MOVE THINGS TO THE CORRECT ARRAY LOCATION
C FIRST, CALCULATE SALC INDICES
C
         JNDEX=0
         KSALC=0
         DO 300 KREP=1,N_REP
          KSALC=KSALC+NDMREP(KREP)
          INDEX=INDBEG(ISA,KREP)
          DO LI =0,LSYMMAX(IFNCT)
           DO IBASE=1,N_CON(LI+1,IFNCT)
            DO IQ=1,N_SALC(KSALC,LI+1,IS)
             INDEX=INDEX+1
             IND_SALC(IQ,IBASE,LI+1,1)=INDEX
            END DO
           END DO
          END DO
          INDEX=INDBEG(JSA,KREP)
          DO LJ =0,LSYMMAX(JFNCT)
           DO JBASE=1,N_CON(LJ+1,JFNCT)
            DO JQ=1,N_SALC(KSALC,LJ+1,JS)
             INDEX=INDEX+1
             IND_SALC(JQ,JBASE,LJ+1,2)=INDEX
            END DO
           END DO
          END DO
C
C END CALCULATION OF SALC INDICES FOR KREP
C
!YY right here
          DO 280 LI=0,LSYMMAX(IFNCT)
           DO LJ=0,LSYMMAX(JFNCT)
            DO IBASE=1,N_CON(LI+1,IFNCT)
             DO JBASE=1,N_CON(LJ+1,JFNCT)
              DO IQ=1,N_SALC(KSALC,LI+1,IS)
               II=IND_SALC(IQ,IBASE,LI+1,1)
               IJ=IND_SALC(1 ,JBASE,LJ+1,2)-1
               DO JQ=1,N_SALC(KSALC,LJ+1,JS)
                IJ=IJ+1
                JNDEX=JNDEX+1
                IF (JNDEX.GT.MTEMP) THEN
                 write(6,*)'OVERLAP: MTEMP IS TOO SMALL'
                 CALL STOPIT
                END IF
                IF (IJ.GE.II) THEN
                 KNDEX=IPOINT(KREP)+1+(IJ-II)
     &                +(NS_TOT(KREP)      *(NS_TOT(KREP)+1)
     &                -(NS_TOT(KREP)-II+1)*(NS_TOT(KREP)-II+2))/2
                 HSTOR(KNDEX,MODSTOR)=HSTOR(KNDEX,MODSTOR)+HTEMP(JNDEX)
                ELSE
                 KNDEX=IPOINT(KREP)+1+(II-IJ)
     &                +(NS_TOT(KREP)      *(NS_TOT(KREP)+1)
     &                -(NS_TOT(KREP)-IJ+1)*(NS_TOT(KREP)-IJ+2))/2
                END IF
               END DO
              END DO
             END DO
            END DO
           END DO
  280     CONTINUE
  300    CONTINUE
  320   CONTINUE
  400  CONTINUE
C
C ADD STORED NONLOCAL HAMILTONIAN TO HSTOR
C
  500  IF (((MODE.EQ.2).OR.(MODE.EQ.3)).AND.(ISITPSP.EQ.1)) THEN
        IF (.NOT.HNLBABY) GOTO 550
        OPEN(66,FILE='HAMBABY',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(66)
        READ(66) NREC,NBLOCK
        IF ((NREC.NE.NHTOT).OR.(NBLOCK.NE.MTEMP)) GOTO 550
        DO IOFS=0,NHTOT-1,MTEMP
         NREC=MIN(MTEMP,NHTOT-IOFS)
         READ(66)(HTEMP(IREC), IREC=1,NREC)
         DO IREC=1,NREC
          HSTOR(IREC+IOFS,MODSTOR)=HSTOR(IREC+IOFS,MODSTOR)+HTEMP(IREC)
         END DO
        END DO
        CLOSE(66)
        GOTO 600
C
C READ ERROR
C
  550   write(6,*)'OVERLAP: NONLOCAL HAMBABY IS INCOMPATIBLE' 
        CALL STOPIT
       END IF
C
C CHECK IF RESULTS SHOULD BE WRITTEN TO FILE
C
  600  IF ((MODE.EQ.1).AND.(.NOT.OVLBABY)) THEN
        FNAME='OVLBABY'
        OVLBABY=.TRUE.
        GOTO 700
       END IF
       IF (((MODE.EQ.2).OR.(MODE.EQ.3)).AND.
     &     (.NOT.HAMBABY).AND.(.NOT.STARTUP)) THEN
        FNAME='HAMBABY'
        HAMBABY=.TRUE.
        HNLBABY=.FALSE.
        GOTO 700
       END IF
       IF ((MODE.EQ.4).AND.(.NOT.KINBABY)) THEN
        FNAME='KINBABY'
        KINBABY=.TRUE.
        GOTO 700
       ENDIF
       GOTO 800
C
C WRITE RESULTS
C
  700  OPEN(66,FILE=FNAME,FORM='UNFORMATTED',STATUS='UNKNOWN')
       REWIND(66)
       WRITE(66) NHTOT
       WRITE(66)(HSTOR(IREC,MODSTOR),IREC=1,NHTOT)
       CLOSE(66)
C
C GET NUMERICAL PART OF HAMILTONIAN INTEGRALS
C
  800  IF (MODE.EQ.2) CALL NUMHAM
       RETURN
      END
