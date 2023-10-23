C UTEP Electronic Structure Lab (2020)
C
C ************************************************************
C
      SUBROUTINE OVERNUM(IZERO) 
C
C WRITTEN BY MARK R PEDERSON (1998)
C CALCULATE OVERLAP-LIKE MATRIX AS AN  O(N) PROBLEM
C USE AS MUCH SYMMETRY AS POSSIBLE.
C   STORE RESULT IN HSTOR(I,2)
C PARALLEL VERSION OF JENS KORTUS, MARK PEDERSON (14-APRIL 1998)
C
       use hstor1,only : hstor
       use debug1
       use common2,only : RIDT, IFUIDT, NIDENT, N_CON, LSYMMAX,
     &   N_POS, NFNCT
       use common3,only : RMAT
       use common5,only : HOLD, HTEMP
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC,
     &   INDBEG, NS_TOT
       use mpidat1,only : NPROC
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IZERO, IBASE, IFNCT, IFNP, IH, II, IID, IIDK, IIDL,
     & IIDP, IIP, IJ, IND_SALC, INDEX, IPOINT, IPOS, IPW, IQ, IREC, IS,
     & ISA, ISHELL, ISIT, ISITE, ISITEK, ISITEP, ISTRT, ITIMES, JB,
     & JBASE, JFNCT, JFNP, JH, JID, JIDK, JIDL, JIDP, JNDEX, JNDSV, JQ,
     & JS, JSA, JSHELL, JSIT, JSITE, JSITEK, JSITEP, JSTRT, JTEMP,
     & K_REP, K_ROW, KNDEX, KOUNT, KREP, KSALC, L, LI, LJ, LNDX, MI,
     & MJ, MNUC, MNUCI, MNUCJ, MODSTOR, MPROC, MSITES, MTMAX, MUI, MUJ,
     & MYIP, MYRES, NADD, NBSMAX, NDEG, NEWIT, NHTOT, NMPROC, NORBS,
     & NSITE, NTIMES
       REAL*8 :: ADD , AI, AIP, AJ, AJP, FF, HTIM1, HTIM2, HTIM3,
     & HTIM4, PHOLD, RNUC, RNUCI, RNUCJ, TIME1, TIME2
       SAVE
       COMMON/HAMTIM/HTIM1,HTIM2,HTIM3,HTIM4
       LOGICAL DOIT
       DIMENSION AIP(3,NPROC+1), AJP(3,NPROC+1)
       DIMENSION IFNP(NPROC+1),JFNP(NPROC+1),ISITEP(NPROC+1)
     &          ,JSITEP(NPROC+1),IIDP(NPROC+1),JIDP(NPROC+1)
       DIMENSION PHOLD(MAXUNSYM,MAXUNSYM,NPROC+1)
C JK98
       LOGICAL FIRST
       DIMENSION NDEG(3),IPOINT(MAX_REP)
       DIMENSION RNUC(3,MX_GRP)
       DIMENSION RNUCI(3,MX_GRP),RNUCJ(3,MX_GRP)
       DIMENSION IND_SALC(ISMAX,MAX_CON,3,2)
       DIMENSION LNDX(6,MAX_CON,3,2)
       DATA NDEG/1,3,6/
       DATA FIRST/.TRUE./
       HTIM1=0.0D0
       HTIM2=0.0D0
       HTIM3=0.0D0
       HTIM4=0.0D0
C
C DEFINE NUMBER OF AVAILABLE WORKER PROCESSES
C
       MPROC=NPROC
c
#ifdef MPI
c
       IF (DEBUG) THEN
        write(6,*)'total number of processors in overnum:',
     &          NPROC+1
       END IF
c
#endif
c
       CALL GTTIME(TIME1)
C
C CHECK AND PRINT BASIS SET PROPERTIES
C
       IF (FIRST) THEN
        FIRST=.FALSE.
        NORBS=0
        NBSMAX=0
        IID=0
        DO 10 IFNCT=1,NFNCT
         NSITE=0
         DO IPOS=1,N_POS(IFNCT)
          IID=IID+1
          CALL GASITES(1,RIDT(1,IID),MNUC,RNUC,MSITES)
          NSITE=NSITE+MNUC
         END DO
         NADD=N_CON(1,IFNCT)+3*N_CON(2,IFNCT)+6*N_CON(3,IFNCT)
         NORBS=NORBS+NADD*NSITE
         NBSMAX=MAX(NBSMAX,NADD)
   10   CONTINUE
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
C
       MODSTOR=2
C
C JK98/DVP98: send general data to workers
C
c
#ifdef MPI
c
       CALL SENDDATA(104)
c
#endif
c
C DETERMINE NHTOT (TOTAL NUMBER OF NONZERO HAMILTONIAN MATRIX ELEMENTS)
C CHECK WHETHER DATA CAN BE READ FROM FILE
C
       NHTOT=0
       DO KREP=1,N_REP
        NHTOT=NHTOT+(NS_TOT(KREP)*(NS_TOT(KREP)+1))/2
       END DO
       IF (IZERO.EQ.0) THEN
        DO IREC=1,NHTOT
         HSTOR(IREC,MODSTOR)= 0.0D0
        END DO
       END IF
C
C CALCULATE OVERLAP/POTENTIAL MATRICES
C
       IIDL= -1
       JIDL= -1
       IPW=0
       DO 400 IID=1,NIDENT
        DO 380  JID=IID,NIDENT
         IF (IID.NE.IIDL) THEN
          CALL OBINFO(1,RIDT(1,IID),RNUCI,MNUCI,ISHELL)
          CALL GSMAT(ISHELL,1)
          IIDL=IID
         END IF
         IF (JID.NE.JIDL) THEN
          CALL OBINFO(1,RIDT(1,JID),RNUCJ,MNUCJ,JSHELL)
          CALL GSMAT(JSHELL,2)
          JIDL=JID
         END IF
         ISA=IID
         JSA=JID
         IFNCT=IFUIDT(IID)
         JFNCT=IFUIDT(JID)
 
         DO 370 ISITE=1,MNUCI
          IF (ISA.EQ.JSA) THEN
           JB=ISITE
          ELSE
           JB=1
          END IF
          DO 360 JSITE=JB,MNUCJ
           IPW=IPW+1
           AIP(1,IPW)=RNUCI(1,ISITE)
           AIP(2,IPW)=RNUCI(2,ISITE)
           AIP(3,IPW)=RNUCI(3,ISITE)
           IFNP(IPW)=IFNCT
           IIDP(IPW)=IID
           ISITEP(IPW)=ISITE
c
           AJP(1,IPW)=RNUCJ(1,JSITE)
           AJP(2,IPW)=RNUCJ(2,JSITE)
           AJP(3,IPW)=RNUCJ(3,JSITE)
           JFNP(IPW)=JFNCT
           JIDP(IPW)=JID
           JSITEP(IPW)=JSITE
           DOIT=.TRUE.
           IF (IID.NE.NIDENT ) DOIT=.FALSE.
           IF (JID.NE.NIDENT ) DOIT=.FALSE.
           IF (ISITE.NE.MNUCI) DOIT=.FALSE.
           IF (JSITE.NE.MNUCJ) DOIT=.FALSE.
           IF (IPW.EQ.MPROC+1) DOIT=.TRUE.
           IF (DOIT) THEN
            NMPROC=IPW
c
c send tasks to workers 
c
#ifdef MPI
c
            DO MYIP=1,NMPROC-1
             CALL PAMHAMIL(1,IFNP(MYIP),JFNP(MYIP),AIP(1,MYIP),
     &                     AJP(1,MYIP),HOLD,MYIP)
            ENDDO
c
#endif
c
c have the manager do some work too
c
            MYIP=NMPROC
            NEWIT=1
            CALL GETHOLD(NEWIT,IFNP(MYIP),JFNP(MYIP),AIP(1,MYIP),
     &                   AJP(1,MYIP),HOLD)
            DO JH=1,MAXUNSYM
             DO IH=1,MAXUNSYM
              PHOLD(IH,JH,MYIP)=HOLD(IH,JH)
             END DO
            END DO
c
c get data back from workers
c
#ifdef MPI
c
            DO IIP=1,NMPROC-1
             CALL PAMHAMIL(2,IFNCT,JFNCT,AI,AJ,HOLD,MYRES)
             ADD=0.0D0
             DO JH=1,MAXUNSYM
              DO IH=1,MAXUNSYM
               PHOLD(IH,JH,MYRES)=HOLD(IH,JH)
               ADD=ADD+ABS(PHOLD(IH,JH,MYRES))
              END DO
             END DO
             IF (DEBUG) write(6,*)'processor:',IIP,ADD
            END DO
c
#endif
c
            DO 350 IPW=1,NMPROC
             IIDK=IIDP(IPW)
             JIDK=JIDP(IPW)
             ISITEK=ISITEP(IPW)
             JSITEK=JSITEP(IPW)
             IF (IIDK.NE.IIDL) THEN
              CALL OBINFO(1,RIDT(1,IIDK),RNUCI,MNUCI,ISHELL)
              CALL GSMAT(ISHELL,1)
              IIDL=IIDK
             END IF
             IF (JIDK.NE.JIDL) THEN
              CALL OBINFO(1,RIDT(1,JIDK),RNUCJ,MNUCJ,JSHELL)
              CALL GSMAT(JSHELL,2)
              JIDL=JIDK
             END IF
             ISA=IIDK
             JSA=JIDK
             IFNCT=IFUIDT(IIDK)
             JFNCT=IFUIDT(JIDK)
             IS=ISHELL
             JS=JSHELL
             DO JTEMP=1,MTMAX
              HTEMP(JTEMP)=0.0D0
             END DO
             DO JH=1,MAXUNSYM
              DO IH=1,MAXUNSYM
               HOLD(IH,JH)=PHOLD(IH,JH,IPW)
              END DO
             END DO
C
C JK98
C
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
C
C NOW UPDATE SALC MATRICES FOR EACH REPRESENTATION
C
             NTIMES=2
             IF (ISA.NE.JSA) NTIMES=1
             IF (ISITEK.EQ.JSITEK) NTIMES=1
             DO 260 ITIMES=1,NTIMES
              IF (ITIMES.EQ.1) THEN
               ISIT=ISITEK
               JSIT=JSITEK
               KSALC=0
               JNDEX=0
               DO 150 K_REP=1,N_REP
                FF=1.0D0/NDMREP(K_REP)
                DO 140 K_ROW=1,NDMREP(K_REP)
                 IF (K_ROW.EQ.1) THEN
                  JNDSV=JNDEX
                 ELSE
                  JNDEX=JNDSV
                 END IF
                 KSALC=KSALC+1
                 DO 120 LI=0,LSYMMAX(IFNCT)
                  ISTRT=(ISIT-1)*NDEG(LI+1)
                  DO LJ=0,LSYMMAX(JFNCT)
                   JSTRT=(JSIT-1)*NDEG(LJ+1)
                   DO IBASE=1,N_CON(LI+1,IFNCT)
                    DO JBASE=1,N_CON(LJ+1,JFNCT)
                     DO MUJ=1,NDEG(LJ+1)
                      DO MUI=1,NDEG(LI+1)
                       KNDEX=JNDEX
                       DO IQ=1,N_SALC(KSALC,LI+1,IS)
                        DO JQ=1,N_SALC(KSALC,LJ+1,JS)
                         KNDEX=KNDEX+1
                         HTEMP(KNDEX)=HTEMP(KNDEX)
     &                   +U_MAT(MUI+ISTRT,IQ,KSALC,LI+1,1)
     &                   *U_MAT(MUJ+JSTRT,JQ,KSALC,LJ+1,2)*FF
     &                   *HOLD(LNDX(MUJ,JBASE,LJ+1,2),
     &                         LNDX(MUI,IBASE,LI+1,1))
                        END DO
                       END DO
                      END DO
                     END DO
                     JNDEX=KNDEX
                    END DO
                   END DO
                  END DO
  120            CONTINUE
  140           CONTINUE
  150          CONTINUE
              ELSE
               ISIT=ISITEK
               JSIT=JSITEK
               KSALC=0
               JNDEX=0
               DO 250 K_REP=1,N_REP
                FF=1.0D0/NDMREP(K_REP)
                DO 240 K_ROW=1,NDMREP(K_REP)
                 IF (K_ROW.EQ.1) THEN
                  JNDSV=JNDEX
                 ELSE
                  JNDEX=JNDSV
                 END IF
                 KSALC=KSALC+1
                 DO 220 LJ=0,LSYMMAX(JFNCT)
                  JSTRT=(JSIT-1)*NDEG(LJ+1)
                  DO LI=0,LSYMMAX(IFNCT)
                   ISTRT=(ISIT-1)*NDEG(LI+1)
                   DO JBASE=1,N_CON(LJ+1,JFNCT)
                    DO IBASE=1,N_CON(LI+1,IFNCT)
                     DO MUI=1,NDEG(LI+1)
                      DO MUJ=1,NDEG(LJ+1)
                       KNDEX=JNDEX
                       DO JQ=1,N_SALC(KSALC,LJ+1,JS)
                        DO IQ=1,N_SALC(KSALC,LI+1,IS)
                         KNDEX=KNDEX+1
                         HTEMP(KNDEX)=HTEMP(KNDEX)
     &                   +U_MAT(MUJ+JSTRT,JQ,KSALC,LJ+1,2)
     &                   *U_MAT(MUI+ISTRT,IQ,KSALC,LI+1,1)*FF
     &                   *HOLD(LNDX(MUJ,JBASE,LJ+1,2),
     &                         LNDX(MUI,IBASE,LI+1,1))
                        END DO
                       END DO
                      END DO
                     END DO
                     JNDEX=KNDEX
                    END DO
                   END DO
                  END DO
  220            CONTINUE
  240           CONTINUE
  250          CONTINUE
              END IF
  260        CONTINUE
C
C MOVE THINGS TO THE CORRECT ARRAY LOCATION
C FIRST, CALCULATE SALC INDICES
C
             JNDEX=0
             KSALC=0
             DO 300 KREP=1,N_REP
              KSALC=KSALC+NDMREP(KREP)
              INDEX=INDBEG(ISA,KREP)
              DO LI=0,LSYMMAX(IFNCT)
               DO IBASE=1,N_CON(LI+1,IFNCT)
                DO IQ=1,N_SALC(KSALC,LI+1,IS)
                 INDEX=INDEX+1
                 IND_SALC(IQ,IBASE,LI+1,1)=INDEX
                END DO
               END DO
              END DO
              INDEX=INDBEG(JSA,KREP)
              DO LJ=0,LSYMMAX(JFNCT)
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
     &               +(NS_TOT(KREP)      *(NS_TOT(KREP)+1)
     &               -(NS_TOT(KREP)-II+1)*(NS_TOT(KREP)-II+2))/2
                     HSTOR(KNDEX,MODSTOR)=HSTOR(KNDEX,MODSTOR)
     &                                   +HTEMP(JNDEX)
                    ELSE
                     KNDEX=IPOINT(KREP)+1+(II-IJ)
     &               +(NS_TOT(KREP)      *(NS_TOT(KREP)+1)
     &               -(NS_TOT(KREP)-IJ+1)*(NS_TOT(KREP)-IJ+2))/2
                    END IF
                   END DO
                  END DO
                 END DO
                END DO
               END DO
  280         CONTINUE
  300        CONTINUE
  350       CONTINUE
            IPW=0
           END IF
  360     CONTINUE
  370    CONTINUE
  380   CONTINUE
  400  CONTINUE
       CALL GTTIME(TIME2)
       PRINT 500,TIME2-TIME1
  500  FORMAT('TIME IN OVERNUM:',F12.4)
       IF (DEBUG) write(6,*)'HTIME : ',HTIM1,HTIM2,HTIM3,HTIM4
       RETURN
      END
