C UTEP Electronic Structure Lab (2020)
C
C ****************************************************************
C
       SUBROUTINE READWF2(FAILED)
       use for_diag1
       use hstor1,only : hstor
       use debug1
       use common2,only : WFFILE, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF,
     &   NWFS, EVLOCC, ISTSCF
       use common8,only : REP, N_REP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_OCC, I_REP, IB, IERR, IREP, ITER, ITOT, ITOT2,
     & IWF, J, JB, JWF, KB, MREP, MSPN, NBAS, NBASF
       REAL*8 :: DOT , ER1, ERR, SQRT2, TMP
       SAVE
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
       LOGICAL EXIST
       CHARACTER*12 EVALSTR
C
C READ IN ATOMIC LOCATIONS AND BASIS SET INFORMATION
C
       FAILED= .FALSE.
       PRINT '(2A)','READING OLD WAVEFUNCTIONS FROM FILE ',WFFILE
       OPEN(99,FILE=WFFILE,FORM='UNFORMATTED',STATUS='UNKNOWN')
       READ(99,END=900) MSPN
       IF (NSPN.GT.MXSPN) THEN
        write(6,*)'READWF: NSPN > MXSPN'
        CALL STOPIT
       END IF
       READ(99,END=900) NWF,NWFS(1)
       NWFS(NSPN)=NWFS(1)
       IF (NWF.GT.MAX_OCC) THEN
        write(6,*)'READWF: MAX_OCC MUST BE AT LEAST: ',NWF
        CALL STOPIT
       END IF
C       WRITE(6,*) 'Read from WFOUT'
C       WRITE(6,*) 'NSPN=',NSPN
C       WRITE(6,*) 'MSPN=',MSPN
C       WRITE(6,*) 'NWF=',NWF
C       WRITE(6,*) 'NWFS=',NWFS(1),NWFS(2)
       READ(99,END=900)(EVLOCC(IWF), IWF=1,NWF)
       J=0
       DO I=NWF+1,2*NWF
          J=J+1
          EVLOCC(I)=EVLOCC(J)
       END DO
       NWF=2*NWF
       DO I=1,NWF
         DO J=I+1,NWF
          IF(EVLOCC(J).LE.EVLOCC(I)) THEN
           TMP=EVLOCC(I)
           EVLOCC(I)=EVLOCC(J)
           EVLOCC(J)=TMP
          END IF
         END DO 
       END DO 
       ITOT=0
       READ(99,END=900) MREP
C       WRITE(6,*) 'MREP=',MREP
C       WRITE(6,*) 'N_PREP=',N_REP
       IF (MREP.NE.N_REP) GOTO 900
C       DO 150 ISPN=1,NSPN
C        DO I_REP=1,N_REP
C         READ(99,END=900) N_OCC(I_REP,ISPN),NBASF
C         IF (NBASF.NE.NS_TOT(I_REP)) GOTO 900
C         READ(99,END=900)(OCCUPANCY(I_OCC+ITOT),
C     &                    I_OCC=1,N_OCC(I_REP,ISPN))
C         ITOT=ITOT+N_OCC(I_REP,ISPN)
C         DO IWF=1,N_OCC(I_REP,ISPN)
C          READ(99,END=900)(PSI_COEF(I,IWF,I_REP,ISPN),
C     &                     I=1,NS_TOT(I_REP))
C         END DO
C        END DO
        WRITE(6,*)'NSPN : ',NSPN
        ISPN=1
        SQRT2=1.0D0/DSQRT(2.0D0)
        DO I_REP=1,N_REP
         READ(99,END=900) N_OCC(I_REP,ISPN),NBASF
         IF (NBASF.NE.NS_TOT(I_REP)) GOTO 900
         READ(99,END=900)(OCCUPANCY(I_OCC+ITOT),
     &                    I_OCC=1,N_OCC(I_REP,ISPN))
         DO I_OCC =1,N_OCC(I_REP,ISPN)
            OCCUPANCY(I_OCC+ITOT)=OCCUPANCY(I_OCC+ITOT)
         ENDDO
         ITOT=ITOT+N_OCC(I_REP,ISPN)
         DO IWF=1,N_OCC(I_REP,ISPN)
          READ(99,END=900)(PSI_COEF(I,IWF,I_REP,ISPN),
     &                     I=1,NS_TOT(I_REP))
          DO I=1,NS_TOT(I_REP)
             PSI_COEF(I,IWF,I_REP,ISPN)=PSI_COEF(I,IWF,I_REP,ISPN)*SQRT2
          END DO 
         END DO
       ENDDO
       WRITE(6,*) 'FINISHED PROCESSING SPIN 1'
       ISPN=2
       ITOT2=0
        DO I_REP=1,N_REP
         N_OCC(I_REP,NSPN)=N_OCC(I_REP,1)
c         IF (NBASF.NE.NS_TOT(I_REP)) GOTO 900
         DO I_OCC=1,N_OCC(I_REP,1)
           OCCUPANCY(I_OCC+ITOT)=OCCUPANCY(I_OCC+ITOT2)
         END DO
         ITOT=ITOT+N_OCC(I_REP,NSPN)
         ITOT2=ITOT2+N_OCC(I_REP,1)
         DO IWF=1,N_OCC(I_REP,1)
          DO I=1,NS_TOT(I_REP)
            PSI_COEF(I,IWF,I_REP,ISPN)=PSI_COEF(I,IWF,I_REP,1)
          END DO
         END DO
       ENDDO
C THE LINES BELOW ARE TO BE USED ONLY FOR IP CALCULATIONS
       IF(ISTSCF==12)THEN
         OCCUPANCY(I_OCC+ITOT2)= OCCUPANCY(I_OCC+ITOT2-1)-1
       ENDIF
C THE LINES BELOW ARE TO BE USED ONLY FOR EA CALCULATIONS
       IF(ISTSCF==13)THEN
         OCCUPANCY(I_OCC+ITOT2)= OCCUPANCY(I_OCC+ITOT2-1)+1
       ENDIF

C  150  CONTINUE
       CLOSE(99)
C
C Print out WF for verification
C
C      CALL WFOUT2
C      CALL STOPIT
C
C REORTHONORMALIZE OCCUPIED WAVEFUNCTIONS:
C
       IF (DEBUG) write(6,*)'READWF CALLS OVERLAP MODE: 1'
       CALL OVERLAP(1)
       DO 240 ISPN=1,NSPN
       KB=0
       DO 240 IREP=1,N_REP
        NBAS=NS_TOT(IREP)
        ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'readwf:Error allocating Ham'
        ENDIF
        ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'readwf:Error allocating Overlap'
        ENDIF
        DO 160 IB=1 ,NBAS
        DO 160 JB=IB,NBAS
         KB=KB+1
         AOVER(JB,IB)=HSTOR(KB,1)
         AOVER(IB,JB)=HSTOR(KB,1)
  160   CONTINUE
        ERR=0.0D0
        DO 175 IWF=1,N_OCC(IREP,ISPN)
        DO 175 JWF=1,IWF
         DOT=0.0D0
         IF (JWF.EQ.IWF) THEN
          DO IB=1,NBAS
           DO JB=1,NBAS
            DOT=DOT+AOVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
     &                         *PSI_COEF(IB,IWF,IREP,ISPN)
           END DO
          END DO
          ERR=ERR+ABS(DOT-1.0D0)
          DO IB=1,NBAS
           PSI_COEF(IB,IWF,IREP,ISPN)=
     &     PSI_COEF(IB,IWF,IREP,ISPN)/SQRT(DOT)
          END DO
           
          DO IB=1,NBAS
           AHAM(IB,IWF)=0.0D0
           DO JB=1,NBAS
            AHAM(IB,IWF)=AHAM(IB,IWF)
     &                 +PSI_COEF(JB,IWF,IREP,ISPN)*AOVER(JB,IB)
           END DO
          END DO
         ELSE
          DOT=0.0D0
          DO IB=1,NBAS
           DOT=DOT+PSI_COEF(IB,IWF,IREP,ISPN)*AHAM(IB,JWF)
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
        IF(DEBUG) THEN
        ERR=0.0D0
        DO 185 IWF=1,N_OCC(IREP,ISPN)
        DO 185 JWF=1,IWF
         DOT=0.0D0
         DO IB=1,NBAS
          DO JB=1,NBAS
           DOT=DOT+AOVER(JB,IB)*PSI_COEF(JB,JWF,IREP,ISPN)
     &                        *PSI_COEF(IB,IWF,IREP,ISPN)
          END DO
         END DO
         IF (IWF.EQ.JWF) THEN
          ERR=ERR+ABS(DOT-1.0D0)
         ELSE
          ERR=ERR+ABS(DOT)
         END IF
  185   CONTINUE
        write(6,*)'REORTHOGONALIZATION ERROR:',ERR,ER1
        ENDIF
        DEALLOCATE(AHAM,STAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'readwf:Error deallocating Ham'
        ENDIF
        DEALLOCATE(AOVER,STAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'readwf:Error deallocating Overlap'
        ENDIF
  240  CONTINUE
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
