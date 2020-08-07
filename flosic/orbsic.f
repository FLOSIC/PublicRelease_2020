C UTEP Electronic Structure Lab (2020)
C ****************************************************************
C
       SUBROUTINE ORBSIC(LSIC)
C       INCLUDE 'PARAMS'
C       INCLUDE 'commons.inc'
C
C
       use debug1
       use common2,only : N_CON, LSYMMAX, ISPN, NSPN, WFFILE
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF, NWFS,
     & EVLOCC
       use common7,only : T1UNRV, T2UNRV
       use for_diag1,only : HAM, OVER, filo, SC1, eval
       use hstor1,only : HSTOR
       use common8,only : REP, N_REP, NDMREP, U_MAT, N_SALC, INDBEG
     &                   ,NS_TOT
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
       use ORBENG,only : EVALOCC
       INCLUDE  'PARAMA2'
       INTEGER,PARAMETER :: MAX_TOT=NDH*MAX_REP
       LOGICAL FAILED
       LOGICAL EXIST, LSIC
       CHARACTER*12 EVALSTR
!      REAL*8, DIMENSION (NDH):: EVALSAV(NDH), OCCTMP(MAX_TOT*MXSPN)
!      REAL*8, DIMENSION (NDH):: EVL(NDH)
       REAL*8, allocatable :: EVALSAV(:),OCCTMP(:),EVL(:)
C      COMMON/ORBENG/EVALOCC(MAX_OCC)
C      COMMON/EVLTMP/EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
C      INCLUDE 'commons.inc'
       INTEGER :: MD, KSPN, KORB, IFNCT, ISHELLA, I_SITE, N_NUC, ILOC,
     & I_LOCAL, I_SALC, IBASE, ICALL1, ICALL2, IL, IMS, IND_SALC,
     & INDEX, IOCC, IORB, IORBB, IORBE, IQ, IQ_BEG, IROW, ISHELL,
     & ISHELLV, ITOT, IWF, J_LOCAL, JORB, JWF, K_REP, K_ROW, KSALC, LI,
     & LSPN, LSPNB, LSPNE, MLOCAL, MO, MSTOT, MU, NDEG
       REAL*8 :: RVEC , RVECI, ADD1, ADD2, ADD3, ADDTMAT, PHILOC,
     & TIMER1, TIMER2, TMATMAX
       INTEGER :: MSPN, I, MREP, I_REP,I_OCC,NBASF, KB, IB, JB, IREP
       INTEGER :: ITER, NBAS, J, LBAS, KBAS, IBAS, JBAS,K
       REAL*8 :: ERR, DOT,ER1, ERROR, EVV

       ALLOCATE(EVALSAV(NDH))
       ALLOCATE(OCCTMP(MAX_TOT*MXSPN))
       ALLOCATE(EVL(NDH))
C
       write(6,*)'from orbsic'
       print *, 'NWFS', NWFS(1), NWFS(2)
       print *, 'NS_TOT NBAS', NS_TOT(1), NS_TOT(2), NBAS
c
c  KAJ  Here NWFS holds the number of occupied states of spin 1 and 2
c       NS_TOT should be the number of basis functions of irep 1 and 2
c       As noted below, this subroutine is currently using 
c         NS_TOT(i) to be NBAS for spin i.  This is a bug that will
c         probably make this routine fail if symmetry is used.
c
c  KAJ 2-27-2020  SHIFT occupied states down to insure
c  that anionic systems can start from wave functions
c
       DO 250 ISPN=1,NSPN
         DO I=1,NWFS(ISPN)
            J=I+(ISPN-1)*NWFS(1)
          EVL(I)=EVLOCC(J) - 10.0d0
         END DO

       CALL OVERLAP(1)
       KB=0
       DO IREP=1,N_REP
        NBAS=NS_TOT(IREP)
       DO IB=1,NBAS
        DO JB=1,NBAS
             HAM(IB,JB)=0.0D0
             OVER(IB,JB)=0.0D0
             FILO(IB,JB)=0.0D0
        END DO
       END DO
c
c  KAJ:  build "Hamiltonian" in spectral representation
c  based on number of states read in from WFOUT
c  If this number < NBAS, diagonalizing this Hamiltonian
c  will return good evecs for the number of states read in.
c  The rest will correspond to zero eigenvalues and will
c  be arbitrary linear combinations of the virtual states.
c

        DO 160 IB=1 ,NBAS
        DO 160 JB=IB,NBAS
         KB=KB+1
         OVER(JB,IB)=HSTOR(KB,1)
         OVER(IB,JB)=HSTOR(KB,1)
         HAM(JB,IB)=0.0D0
         HAM(IB,JB)=0.0D0
         FILO(IB,JB)=0.0D0
         FILO(JB,IB)=0.0D0
  160   CONTINUE
         DO KBAS=1,NBAS
          DO LBAS=1,NBAS
           DO IWF=1,NWFS(ISPN)
           HAM(KBAS,LBAS)=HAM(KBAS,LBAS)+
     &    EVL(IWF)*PSI_COEF(KBAS,IWF,IREP,ISPN)
     &             *PSI_COEF(LBAS,IWF,IREP,ISPN)
           END DO
          END DO
         END DO
         ERROR=0.0D0
         DO KBAS=1,NBAS
          DO LBAS=KBAS,NBAS
           ERROR=ERROR+ABS(HAM(KBAS,LBAS)-HAM(LBAS,KBAS))
          END DO
         END DO
         WRITE(6,*)'ERROR=',ERROR
        
         DO KBAS=1,NBAS
           DO JBAS=1,NBAS
             DO LBAS=1,NBAS
               FILO(KBAS,JBAS)=FILO(KBAS,JBAS)
     &                        +HAM(KBAS,LBAS)*OVER(LBAS,JBAS) 
             END DO
           END DO
         END DO
             
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             HAM(IBAS,JBAS)=0.0D0
           END DO
         END DO
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               HAM(IBAS,JBAS)=HAM(IBAS,JBAS)
     &                        +OVER(IBAS,KBAS)*FILO(KBAS,JBAS) 
             END DO
           END DO
         END DO
c
c  KAJ: diagonalize spectral Hamiltonian
c         
        CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
         WRITE(6,*)'AFTER DIAG'
         DO I=1,NBAS
             WRITE(6,*) EVAL(I)
         END DO
c           EVV=EVAL(NWFS(ISPN))
        DO I=1,NWFS(ISPN)
           J=I+(ISPN-1)*NS_TOT(1)
c
c  KAJ Note that NS_TOT is indexed by IREP, NOT SPIN
c       The logic above assumes spin! 
c       The current structure will likely only work for a system
c         with no symmetry.
c
c  KAJ unshift evalues -- restore them to their original values.
c
           EVALSAV(J)=EVAL(I) + 10.0d0
        END DO
c       DO I=NWFS(ISPN)+1,NS_TOT(1)
c          J=I+(ISPN-1)*NS_TOT(1)
c  KAJ Change NS_TOT to NBAS
        DO I=NWFS(ISPN)+1,NBAS
           J=I+(ISPN-1)*NBAS
c         EVALSAV(J)=EVV+0.1*I
          EVALSAV(J)= EVAL(I)
        END DO
    
        DO I=1,NS_TOT(IREP)
c           WRITE(6,*) EVALSAV(I)
          DO IBAS=1,NS_TOT(IREP)
            PSI_COEF(IBAS,I,IREP,ISPN)=HAM(IBAS,I)
          END DO
        END DO
       END DO
  250  CONTINUE



        DO I=1,MAX_TOT*MXSPN
           OCCTMP(I)=0.0D0
        END DO
        J=0
        DO ISPN=1,NSPN
           K=(ISPN-1)*NWFS(1)
           DO IREP=1,N_REP
             DO IBAS=1,N_OCC(IREP,ISPN)
              J=J+1
              K=K+1
              OCCTMP(J)=OCCUPANCY(K)
             END DO
             DO IBAS=N_OCC(IREP,ISPN)+1,NS_TOT(IREP)
              J=J+1
              OCCTMP(J)=0.0D0
             END DO
           END DO
        END DO
            WRITE(6,*) 'ISPN, IREP, J, OCCUPANCY(J)'
        J=0
        NWF=0
        DO ISPN=1,NSPN
         NWFS(ISPN)=NS_TOT(1)
         NWF=NWF+NWFS(ISPN)
         DO IREP=1,N_REP
           N_OCC(IREP,ISPN)=NWFS(ISPN)
           DO IBAS=1,NS_TOT(IREP)
            J=J+1
            OCCUPANCY(J)=OCCTMP(J)
            EVALOCC(J)=EVALSAV(J)
            WRITE(6,*) ISPN, IREP, J, OCCUPANCY(J),EVALOCC(J)
 600        FORMAT(I3,I3,I5,F8.4,F12.8)
           END DO
         END DO
        END DO
c           WRITE(6,*) 'ORBSIC:'
c           print *, 'nwf nwfs', NWF, NWFS(1),NWFS(2)

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

       DEALLOCATE(EVL)
       DEALLOCATE(OCCTMP)
       DEALLOCATE(EVALSAV)

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
