C UTEP Electronic Structure Lab (2020)
      SUBROUTINE WFWIND2(EMIN,EMAX,SPN1,SPN2,IFAIL,NBASTOT)
C
C     ------------------------------------------------------------------
C
C      VERSION DIRK POREZAG NOVEMBER 1996 
C      UPDATED BY ULISES REVELES, JUNE 2014
C
C     ------------------------------------------------------------------
C
      use global_inputs,only : idiag2
      use debug1
      use for_diag1
      use hstor1,only : hstor
      use debug1
      use common2,only : ISPN, NSPN
      use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &                   NWF, NWFS, EVLOCC
      use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFAIL, I, IB, IBAS, IDEG, IEIG, IERR, IND1, IND2,
     & IOFS, IREC, IREP, ISPBEG, ISPEND, JBAS, JSPN, KBAS, MSPN, NBAS,
     & NEIG, NREC, NSYMM, NUSYM
       REAL*8 :: EMIN , EMAX, TIME1, TIME2
      SAVE
      LOGICAL EXIST,SPN1,SPN2,DOIT
C
      INTEGER NBASTOT
C
C     ------------------------------------------------------------------
C
C     --- CHECKING AND SETTING UP SOME STUFF ---
C
      IFAIL = 0
      CALL GTTIME(TIME1)
C
      IF (DEBUG) PRINT '(A,2(1X,F15.5))',
     &            ' IN WFWIND: EMIN, EMAX:',EMIN,EMAX
      IF(N_REP.GT.MAX_REP)THEN
        PRINT *,'WFWIND: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
      END IF
C
C     --- INITIALIZATION ---
C
      NSYMM = 0
      NUSYM = 0
      ISPBEG = NSPN
      ISPEND = 1
      IF (SPN1) ISPBEG = 1
      IF (SPN2) ISPEND = NSPN
      NWF = 0
      NWFS(1) = 0
      NWFS(NSPN) = 0
C
      DO JSPN = 1,NSPN
        DO IREP = 1,N_REP
          N_OCC(IREP,JSPN) = 0
        END DO
      END DO
C
C     --- LOOP OVER SPIN ---
C
      DO 240 ISPN = ISPBEG,ISPEND
        DOIT = .FALSE.
        IF (ISPN.EQ.1.AND.SPN1) DOIT=.TRUE.
        IF (ISPN.EQ.2.AND.SPN2) DOIT=.TRUE.
        IF (DOIT) THEN
C
C           IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 1'
C           CALL OVERLAP(1)
C           IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 2'
C           CALL OVERLAP(2)
C
C     --- READ HAMILTONIAN ---
C
        INQUIRE(FILE='HAMOLD',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
          WRITE(6,*)'WFIWND2: HAMOLD DOES NOT EXIST'
          CALL STOPIT
        END IF
C
        OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='OLD')
        REWIND(99)
        READ(99,END=50) NREC,MSPN
C
        IF (MSPN.EQ.1) THEN
          READ(99,END=50)(HSTOR(IREC,2),IREC=1,NREC)
C
        ELSE IF (MSPN.EQ.2) THEN
          DO 40 I=1,ISPN
            READ(99,END=50)(HSTOR(IREC,2),IREC=1,NREC)
  40      CONTINUE
C
        ELSE
          GOTO 50
        END IF
C
        GOTO 60
  50    WRITE(6,*)'WFWIDN2: HAMOLD UNREADABLE'
        CALL STOPIT
  60    CLOSE(99)
C
C     --- READ OVERLAP ---
C
       OPEN(66,FILE='OVLBABY',FORM='UNFORMATTED',STATUS='OLD')
       REWIND(66)
       READ(66,END=65) NREC
C       IF (NREC.NE.NHTOT) THEN
C        write(6,*)'WFWIDN2: READ: WRONG NUMBER OF DATA: ',NREC,NHTOT
C        CALL STOPIT
C       ENDIF
       IF (size(HSTOR,1).LT.NREC) THEN
        write(6,*)'WFWIND2: wrong size of HSTOR'
        write(6,*)'WFWIND2:Number of records is',NREC
        write(6,*)'WFWIND2:size of HSTOR is',size(HSTOR,1)
        CALL STOPIT
       ENDIF
       READ(66,END=65,ERR=65)(HSTOR(IREC,1),IREC=1,NREC)
       CLOSE(66)
       GOTO 66
  65  write(6,*)'WFWIND2: COULD NOT READ STORED DATA'
       write(6,*)'TRIED TO READ FILE: OVLBABY'
       CALL STOPIT
C
C      --- LOOP OVER REPRESENTATIONS ---
C      --- GET MATRIX ELEMENTS       ---
C
  66   KBAS = 0
       NBASTOT = 0
       DO 130 IREP=1,N_REP
         NBAS=NS_TOT(IREP)
         NBASTOT = NBASTOT + NBAS
C
         IF(NBAS.GT.NDH)THEN
           PRINT *,'WFWIND: NDH MUST BE AT LEAST: ',NBAS
           CALL STOPIT
         END IF
C
C     --- ALLOCATE LOCAL FILEDS ---
C
         ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'wfwind:Error allocating Ham'
         ENDIF 
         ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'wfwind:Error allocating Overlap'
         ENDIF 
         ALLOCATE(AEVAL(NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'wfwind:Error allocating Eval'
         ENDIF 
         ALLOCATE(ASC1(NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'wfwind:Error allocating Sc1'
         ENDIF 
         DO 80 IBAS=1,NBAS
           DO 70 JBAS=IBAS,NBAS
             KBAS=KBAS+1
             AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
             AHAM (JBAS,IBAS)=HSTOR(KBAS,2)
   70      CONTINUE
   80    CONTINUE
CJK01/2001
c         CALL SCISSOR(IREP)
CJK01/2001       
C
C    --- GET EIGENVECTORS AND EIGENVALUES   ---
C    --- (WILL BE RETURNED IN HAM AND EVAL) ---
C
         IF(NBAS.GT.0)THEN
C           CALL TRACER('WFWIND DIAGONALIZATION')
            CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,ASC1,1)
C           idiag2=1
C           CALL DIAGGE3(NBAS,1,AHAM,AEVAL,1)
         ENDIF
C
C    --- DETERMINE WHICH EIGENSTATES ARE NEEDED ---
C
         IND1=0
         IND2=0
         DO I=1,NBAS
          IF (IND1.EQ.0) THEN
           IF (AEVAL(I) .GE. EMIN) IND1=I
          END IF
          IF (IND2.EQ.0) THEN
           IF (AEVAL(I) .GT. EMAX) IND2=I
          END IF
         END DO
         IF (IND1.EQ.0) THEN
          IOFS=0
          NEIG=0
         ELSE IF (IND2.EQ.0) THEN
          IOFS=IND1-1
          NEIG=NBAS-IOFS
         ELSE
          IOFS=IND1-1
          NEIG=IND2-IND1
         END IF
         IF (NEIG.GT.MAX_VIRT_PER_SYM) THEN
          PRINT *,'WFWIND: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',NEIG
          IFAIL=1
          RETURN
         END IF
C
C      --- STORE WAVEFUNCTIONS AND EIGENVALUES ---
C
         N_OCC(IREP,ISPN)=NEIG
         NWFS(ISPN)=NWFS(ISPN)+NEIG*NDMREP(IREP)
         NWF=NWF+NEIG*NDMREP(IREP)
         DO IEIG=1,NEIG
          NSYMM=NSYMM+1
          OCCUPANCY(NSYMM)=1.0D0
          DO IB=1,NBAS
           PSI_COEF(IB,IEIG,IREP,ISPN)=AHAM(IB,IEIG+IOFS)
          END DO
          DO IDEG=1,NDMREP(IREP)
            NUSYM=NUSYM+1
            IF (NUSYM.GT.MAX_OCC) THEN
              PRINT *,'WFWIND: MAX_OCC MUST BE AT LEAST ',NUSYM
              IFAIL=1
              RETURN
            END IF
            EVLOCC(NUSYM)=AEVAL(IEIG+IOFS)
           END DO
          END DO
C
C     --- DEALLOCATE LOCAL FIELDS ---      
C
          DEALLOCATE(AHAM,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Ham'
          ENDIF 
          DEALLOCATE(AOVER,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Overlap'
          ENDIF 
          DEALLOCATE(AEVAL,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Eval'
          ENDIF 
          DEALLOCATE(ASC1,STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'wfwind:Error deallocating Sc1'
          ENDIF 
  130   CONTINUE
        END IF
  240 CONTINUE
C
      CALL GTTIME(TIME2)
      CALL TIMOUT('WAVEFUNCTIONS IN ENERGY WINDOW:    ',TIME2-TIME1)
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
