C UTEP Electronic Structure Lab (2020)
      SUBROUTINE TESTBAS
C
C     ORIGINAL VERSION BY MARK R PEDERSON (1990)
C
C     ------------------------------------------------------------------
C
C     GLOBAL DIMENSIONS:
C
      use debug1
      use global_inputs,only : NBO1,inbas
      use for_diag1
      use hstor1,only : hstor
      use common3,only : RMAT
      use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:59 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I_REP, IERR, IOFFSET, IREPSIZE, ISTOP, KBAS, KREP,
     & MAXB, MINBAS, NBAS
       REAL*8 :: DUMMY 
      SAVE
C
C     LOCAL DIMENSIONS:
C
      INTEGER I,J,NTOT
      REAL*8 timea,timeb
C
C     ------------------------------------------------------------------
C
      PRINT '(A)','TESTING BASIS SET'
      IF (N_REP.GT.MAX_REP) THEN
        write(6,*)'TESTBAS: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
      END IF
C
      MINBAS=500     

      IF (DEBUG) write(6,*)'TESTBAS CALLS OVERLAP MODE: 1'
      CALL GTTIME(timea)
      CALL OVERLAP(1)
      CALL GTTIME(timeb)
      WRITE(6,*) 'TESTBAS: back from OVERLAP. TIME:',timeb-timea
      OPEN(97,FILE='OVLPTST',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND(97)
      WRITE(97,*) '******** TEST FOR LINEAR DEPENDENCIES *********'
      WRITE(97,*) '***** EIGENVALUES OF THE OVERLAP MATRICES *****'
C
C     --- DIAGONALIZE OVERLAP MATRIX FOR EACH REPRESENTATION ---
C     --- MOVE OVERLAP MATRIX TO HAM                         ---
C
      ISTOP=0
      DO 70 I_REP=1,N_REP
        WRITE(97,*) 'REPRESENTATION:',I_REP,' DIMEN:',NDMREP(I_REP)
        IF (NS_TOT(I_REP).GT.NDH) THEN
          write(6,*)'TESTBAS: NDH MUST BE AT LEAST: ',NS_TOT(I_REP)
          CALL STOPIT
        END IF
C
C     --- ALLOCATE LOCAL FIELDS ---
C
      NBAS=NS_TOT(I_REP)
C      ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
C      IF(IERR.NE.0) WRITE(6,*)'testbas:Error allocating Ham'
C
      ALLOCATE(AEVAL(NBAS),STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'testbas:Error allocating Eval'
C
C      ALLOCATE(ASC1(NBAS),STAT=IERR)
C      IF(IERR.NE.0) WRITE(6,*)'testbas:Error allocating Sc1'
C
C     --- NOW MOVE OVERLAP MATRIX TO HAM                         ---
C
C      DO I=1,NS_TOT(I_REP)
C        DO J=I,NS_TOT(I_REP)
C          KIND=KIND+1
C            write(6,*)'HSTOR=',HSTOR(KIND,1),kind,I,J
C          AHAM (J,I)=HSTOR(KIND,1)
C        END DO
C       END DO
C
CJUR
C     --- WRITE UNNORMALIZED OVERLAP MATRIX (LOWER TRIANGLE ONLY) ---
C
      IOFFSET =0
      DO KREP=1,I_REP-1
            IREPSIZE=NS_TOT(KREP)
            NTOT = (IREPSIZE*(IREPSIZE+1))/2
            IOFFSET = IOFFSET+NTOT 
      END DO

      KBAS=0
      CALL CHECK_INPUTS
      IF(NBO1) THEN
        OPEN(66,FILE='OVLNRM',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(66)
        WRITE(66) NBAS
        DO I=1,NBAS
          DO J=1,I
            KBAS=KBAS+1
            WRITE(66) HSTOR(KBAS+IOFFSET,1)
          END DO
        END DO
        CLOSE(66)
      ENDIF
CJUR
C     --- NORMALIZE BASIS FUNCTIONS ---
C
        
  
      KBAS=0
      DO 30 I=1,NBAS
        DO 35 J=I,NBAS
C
C     --- The value of first corresponds to the first element of the ---
C     --- given representation in HSTOR this value will be given to  ---
C     --- DIAGSP as the initial position to do the diagonalization   ---
C
C          AEVAL(I)=1.0D0/SQRT(HSTOR(iindex,1))
C
          KBAS=KBAS+1
          IF(I.EQ.J) THEN
           AEVAL(I)=1.0D0/SQRT(HSTOR(KBAS+IOFFSET,1))
          ENDIF
C
   35 CONTINUE
   30 CONTINUE
C
      KBAS=0
      DO 50 I=1,NBAS
        DO 40 J=I,NBAS
          KBAS=KBAS+1
          HSTOR(KBAS+IOFFSET,1)=HSTOR(KBAS+IOFFSET,1)*AEVAL(I)*AEVAL(J)
C
C          AHAM(J,I)=AHAM(J,I)*ASC1(I)*ASC1(J)
C
   40   CONTINUE
   50 CONTINUE

       CALL GTTIME(timea)
       WRITE(6,*)'TESTBAS: pre-diag time', timea-timeb
C
C     --- DIAGONALIZE ---
C
      MAXB=NS_TOT(I_REP)
      IF (MAXB.NE.0) THEN
C          CALL DIAGSP2(MAXB,ifirst,AEVAL,0)
C          CALL DIAGSP(MAXB,MAXB,AHAM,AEVAL,ASC1,0)

!          CALL DIAG_P_SYM(NBAS,IOFFSET+1,1,DUMMY,AEVAL,0)

#ifdef SCALAPACK
          IF(NBAS.GE.MINBAS)THEN
!          CALL DIAG_S_SYM(NBAS,IOFFSET,AEVAL,0)
           INBAS=NBAS
           CALL SENDDATA(113)
          ELSE
           CALL DIAG_P_SYM(NBAS,IOFFSET+1,1,DUMMY,AEVAL,0)
          END IF
#else
          CALL DIAG_P_SYM(NBAS,IOFFSET+1,1,DUMMY,AEVAL,0)
#endif
          CALL GTTIME(timeb)
          WRITE(6,*)'TESTBAS:diag time',timeb-timea
       
C
C
        WRITE(97,*) NDMREP(I_REP),MAXB
        WRITE(97,1000)(AEVAL(I),I=1,MAXB)
      END IF
C
      DO 60 I=1,MAXB
        IF (AEVAL(I) .LE. 0.0D0) ISTOP=ISTOP+1
   60 CONTINUE
C
C     --- DEALLOCATE LOCAL FIELDS ---
C
C      DEALLOCATE(AHAM,STAT=IERR)
C      IF(IERR.NE.0) WRITE(6,*)'testbas:Error deallocating Ham'
C
C        DEALLOCATE(AOVER,STAT=IERR)
C        IF(IERR.NE.0)THEN
C          WRITE(6,*)'testbas:Error deallocating Overlap'
C        ENDIF
C
      DEALLOCATE(AEVAL,STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'testbas:Error deallocating Eval'
C      DEALLOCATE(ASC1,STAT=IERR)
C      IF(IERR.NE.0) WRITE(6,*)'testbas:Error deallocating Sc1'
C
   70 CONTINUE
      IF (ISTOP.GT.0) THEN
        write(6,*)'TESTBAS: ERROR, ISTOP: ',ISTOP
        CALL STOPIT
      END IF
      CLOSE(97)
C
C
      RETURN
C
 1000 FORMAT(' ',5G15.7)
C
C     ------------------------------------------------------------------
C
      END
