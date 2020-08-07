C UTEP Electronic Structure Lab (2020)
C ****************************************************************
C
       SUBROUTINE READWFSIC(FAILED,LSIC)
C        implicit real*8(a-h,o-z)
C       INCLUDE 'PARAMS'
       use debug1
       use common2,only : N_CON, LSYMMAX, ISPN, NSPN, WFFILE
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF, NWFS,
     & EVLOCC
       use common7,only : T1UNRV, T2UNRV
       use for_diag1,only : HAM, OVER
       use hstor1,only : HSTOR
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, INDBEG
     &                   ,NS_TOT
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
C       INCLUDE 'commons.inc'
       INCLUDE  'PARAMA2'  
       INTEGER :: MD, KSPN, KORB, IFNCT, ISHELLA, I_SITE, N_NUC, ILOC,
     & I_LOCAL, I_SALC, IBASE, ICALL1, ICALL2, IL, IMS, IND_SALC,
     & INDEX, IOCC, IORB, IORBB, IORBE, IQ, IQ_BEG, IROW, ISHELL,
     & ISHELLV, ITOT, IWF, J_LOCAL, JORB, JWF, K_REP, K_ROW, KSALC, LI,
     & LSPN, LSPNB, LSPNE, MLOCAL, MO, MSTOT, MU, NDEG, NTEST
       REAL*8 :: RVEC , RVECI, ADD1, ADD2, ADD3, ADDTMAT, PHILOC,
     & TIMER1, TIMER2, TMATMAX
       INTEGER :: MSPN, I, MREP, I_REP,I_OCC,NBASF, KB, IB, JB, IREP
       INTEGER :: ITER, NBAS
       REAL*8 :: ERR, DOT,ER1
       
C
C READWF READS IN WAVEFUNCTIONS TO ALLOW RESTARTING A CALCULATION
C AT A SELF-CONSISTENT END POINT. THIS IS CONVENIENT FOR FINE-TUNING
C A SET OF FIT GAUSSIANS AND FOR GETTING GOOD NEW STARTING POINTS.
C IN ORDER TO MAKE SURE THE FILE IS COMPATIBLE, SOME CRITICAL
C VARIABLES ARE COMPARED TO THE CURRENT CALCULATION.
C IMPORTANT: IN ORDER TO HAVE THE CORRECT NS_TOT AVAILABLE, OVERLAP
C MUST HAVE BEEN CALLED BEFORE.
C
       LOGICAL FAILED
       LOGICAL EXIST, LSIC
       CHARACTER*12 EVALSTR
C
C READ IN ATOMIC LOCATIONS AND BASIS SET INFORMATION
C
       print *, 'from readwfsic'
       FAILED= .FALSE.
       IF (LSIC) THEN
         EXIST=.FALSE.
         INQUIRE(FILE='WFSIC',EXIST=EXIST) 
         IF(EXIST) WFFILE='WFSIC'
       END IF
       PRINT '(2A)','READING OLD WAVEFUNCTIONS FROM FILE ',WFFILE
       OPEN(99,FILE=WFFILE,FORM='UNFORMATTED',STATUS='UNKNOWN')
       READ(99,END=900) MSPN
       print *, 'MSPN',MSPN
       IF (MSPN.NE.NSPN) GOTO 900
       IF (NSPN.GT.MXSPN) THEN
        PRINT *,'READWF: NSPN > MXSPN'
        CALL STOPIT
       END IF
       READ(99,END=900) NWF,(NWFS(ISPN), ISPN=1,NSPN)
       write(6,*) NWF,(NWFS(ISPN), ISPN=1,NSPN)
         
c
c  KAJ 3-5-2020  If this is a restart from a DFT calculation
c                NWF and NWFS(i) will the the total number of occupied
c                wave functions and the total number of occ wave
c                 functions of spin i, respectively.  An SIC calculation
c                needs the occupied states and the virtual states in
c                order to build the SIC Hamiltonian matrix elements. 
c                Check below to see if NWF = the sum of n_bas for each
c                irreducible representation.  If it is, then the WFOUT
c                file is compatible with an SIC calculation.  If not
c                ORBSIC is called to get virtual state eigenvectors.
c
       IF (NWF.GT.MAX_OCC) THEN
        PRINT *,'READWF: MAX_OCC MUST BE AT LEAST: ',NWF
        CALL STOPIT
       END IF
       READ(99,END=900)(EVLOCC(IWF), IWF=1,NWF)
!      print *, 'readwfsic evlocc'
!      DO I=1,NWF
!        WRITE(6,*)EVLOCC(I)
!      END DO
       ITOT=0
       READ(99,END=900) MREP
         WRITE(6,*)'MREP =', MREP
       IF (MREP.NE.N_REP) GOTO 900
        ntest = 0
       DO 150 ISPN=1,NSPN
        DO I_REP=1,N_REP
         READ(99,END=900) N_OCC(I_REP,ISPN),NBASF
         WRITE(6,*)'I_REP N_OCC, NBASF NS_TOT=', I_REP, 
     &       N_OCC(I_REP,ISPN), NBASF,
     &       NS_TOT(I_REP)
         ntest = ntest + NS_TOT(i_rep)
         IF (NBASF.NE.NS_TOT(I_REP)) GOTO 900
         READ(99,END=900)(OCCUPANCY(I_OCC+ITOT),
     &                    I_OCC=1,N_OCC(I_REP,ISPN))
!        WRITE(6,*) 'OCCUPANCY  ITOT',itot
!        WRITE(6,*) (OCCUPANCY(I_OCC+ITOT),
!    &                    I_OCC=1,N_OCC(I_REP,ISPN))
         ITOT=ITOT+N_OCC(I_REP,ISPN)
         DO IWF=1,N_OCC(I_REP,ISPN)
          READ(99,END=900)(PSI_COEF(I,IWF,I_REP,ISPN),
     &                     I=1,NS_TOT(I_REP))
         END DO
        END DO
  150  CONTINUE
       CLOSE(99)
       print *, 'readwfsic  ntest, nwf, nwfs',ntest, nwf,
     &            nwfs(1),nwfs(2)
c
c  KAJ 4-1-2020  MOVE CALL TO ORBSIC until after
c     renormalization of evecs
c

C
C REORTHONORMALIZE OCCUPIED WAVEFUNCTIONS: 1st order
c  orthonormalization...
C
       IF (DEBUG) PRINT *,'READWF CALLS OVERLAP MODE: 1'
       CALL OVERLAP(1)
       DO 240 ISPN=1,NSPN
       KB=0
       DO 240 IREP=1,N_REP
        NBAS=NS_TOT(IREP)
        DO 160 IB=1 ,NBAS
        DO 160 JB=IB,NBAS
         KB=KB+1
         OVER(JB,IB)=HSTOR(KB,1)
         OVER(IB,JB)=HSTOR(KB,1)
  160   CONTINUE
        ERR=0.0D0
        DO 175 IWF=1,N_OCC(IREP,ISPN)
        DO 175 JWF=1,IWF
         DOT=0.0D0
         IF (JWF.EQ.IWF) THEN
          DO IB=1,NBAS
           DO JB=1,NBAS
            DOT=DOT+OVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
     &                         *PSI_COEF(IB,IWF,IREP,ISPN)
           END DO
          END DO
          ERR=ERR+ABS(DOT-1.0D0)
          DO IB=1,NBAS
           PSI_COEF(IB,IWF,IREP,ISPN)=
     &     PSI_COEF(IB,IWF,IREP,ISPN)/SQRT(DOT)
          END DO
          DO IB=1,NBAS
           HAM(IB,IWF)=0.0D0
           DO JB=1,NBAS
            HAM(IB,IWF)=HAM(IB,IWF)
     &                 +PSI_COEF(JB,IWF,IREP,ISPN)*OVER(JB,IB)
           END DO
          END DO
         ELSE
          DOT=0.0D0
          DO IB=1,NBAS
           DOT=DOT+PSI_COEF(IB,IWF,IREP,ISPN)*HAM(IB,JWF)
          END DO
          ERR=ERR+ABS(DOT      )
          DO IB=1,NBAS
           PSI_COEF(IB,IWF,IREP,ISPN)=PSI_COEF(IB,IWF,IREP,ISPN)
     &    -PSI_COEF(IB,JWF,IREP,ISPN)*DOT
          END DO
         END IF
  175   CONTINUE
        ER1=ERR
C
C TEST
C
c       IF(DEBUG) THEN
        ERR=0.0D0
        DO 185 IWF=1,N_OCC(IREP,ISPN)
        DO 185 JWF=1,IWF
         DOT=0.0D0
         DO IB=1,NBAS
          DO JB=1,NBAS
           DOT=DOT+OVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
     &                        *PSI_COEF(IB,IWF,IREP,ISPN)
          END DO
         END DO
         IF (IWF.EQ.JWF) THEN
          ERR=ERR+ABS(DOT-1.0D0)
         ELSE
          ERR=ERR+ABS(DOT)
         END IF
  185   CONTINUE
        PRINT *,'REORTHOGONALIZATION ERROR:',ERR,ER1
c       ENDIF
  240  CONTINUE
c
c  LSIC = .TRUE. if this is a FLOSIC calculation
c  CALL ORBSIC if WFOUT file only has eigenvecs for occupied states 
c
       IF (LSIC) THEN
        IF (.NOT.EXIST) THEN
           if(ntest.ne.nwf) then
             print *, 'calling orbsic'
             CALL ORBSIC(LSIC)
             RETURN
           end if
        END IF
       END IF

C
C REMOVE OLD EVALUE FILES AND CREATE DUMMY FILE FOR FIRST ITERATION
C
       ITER=0
  300   ITER=ITER+1
        WRITE(EVALSTR,'(A,I3.3)')'EVAL',ITER
        INQUIRE(FILE=EVALSTR,EXIST=EXIST)
        IF (EXIST) THEN
         OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='OLD')
         CLOSE(98,STATUS='DELETE')
         GOTO 300
        END IF
       CONTINUE
       OPEN(98,FILE='EVAL001',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(98)
       WRITE(98,*) 'NO EIGENVALUES FOR READ WAVEFUNCTIONS'
       CLOSE(98)
       RETURN
C
C FAILURE
C
  900  PRINT '(A)','UNABLE TO READ OR PROCESS OLD WAVEFUNCTIONS'
       PRINT '(A)','PROCEEDING WITH DEFAULT STARTING POINT'
       FAILED= .TRUE.
       RETURN
       END
C
