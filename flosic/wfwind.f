C UTEP Electronic Structure Lab (2020)
      SUBROUTINE WFWIND(EMIN,EMAX,SPN1,SPN2,IFAIL)
C
C     WFWIND VERSION KAJ April 2020
c     In this revised version, 
c     READ eigenvector coefficients from WFOUT  (canonical wf's)
c     READ evalues from EVALUES file
C
c     Earlier version of WFWIND got EVALUES, EVECs from diagonalization of
c     Hamiltonian stored in HSTOR.  This is the DFT Hamiltonian.
C     ------------------------------------------------------------------
C
       use for_diag1
       use hstor1,only : hstor
       use debug1
       use common2,only : ISPN, NSPN
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:07 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IFAIL, I, IB, IBAS, IDEG, idum, IEIG, 
     & IERR, IND1, IND2,
     & I_OCC,IOFS, IREP, ISPBEG, ISPEND, itot, iwf,j,
     & JBAS, JSPN, KBAS, mrep,mspn,
     & NBAS, nbasf, NEIG, nstest, ntest, NSYMM,
     & NUSYM
       REAL*8 :: EMIN , EMAX, TIME1, TIME2
       real*8 :: evalread(ndh,max_rep,mxspn)
       SAVE
       LOGICAL SPN1,SPN2,DOIT
C
C     ------------------------------------------------------------------
C
C     --- CHECKING AND SETTING UP SOME STUFF ---
C
       IFAIL=0
       CALL GTTIME(TIME1)
       IF (DEBUG) PRINT '(A,2(1X,F15.5))',
     &            ' IN WFWIND: EMIN, EMAX:',EMIN,EMAX
       IF(N_REP.GT.MAX_REP)THEN
         PRINT *,'WFWIND: MAX_REP MUST BE AT LEAST: ',N_REP
         CALL STOPIT
       END IF
       print *, 'WFWIND EMIN EMAX', EMIN, EMAX
c
c     read wave functions
c
       open(99,file='WFOUT',form='unformatted')
       read(99) mspn
       if(mspn.ne.nspn) then
              print *, 'spin mixup in wfwind'
              call stopit
       end if
       READ(99) NWF,(NWFS(ISPN), ISPN=1,NSPN)
       write(6,*) NWF,(NWFS(ISPN), ISPN=1,NSPN)
       IF (NWF.GT.MAX_OCC) THEN
        PRINT *,'READWF: MAX_OCC MUST BE AT LEAST: ',NWF
        CALL STOPIT
       END IF
       READ(99)(EVLOCC(IWF), IWF=1,NWF)
       ITOT=0
       READ(99) MREP
       IF (MREP.NE.N_REP) then
              print *, 'mrep wrong'
              call STOPIT
       end if
        ntest = 0
       DO 150 ISPN=1,NSPN
        DO IREP=1,N_REP
         READ(99) N_OCC(IREP,ISPN),NBASF
         ntest = ntest + NS_TOT(irep)
         IF (NBASF.NE.NS_TOT(IREP)) then
                 print *, 'nbasf wrong'
                 call STOPIT
         end if
         READ(99)(OCCUPANCY(I_OCC+ITOT),
     &                    I_OCC=1,N_OCC(IREP,ISPN))
         ITOT=ITOT+N_OCC(IREP,ISPN)
         DO IWF=1,N_OCC(IREP,ISPN)
          READ(99)(PSI_COEF(I,IWF,IREP,ISPN),
     &                     I=1,NS_TOT(IREP))
         END DO
        END DO
  150  CONTINUE
       CLOSE(99)
c
c  Read EVALUES
c
       open(99,file='EVALUES',form='FORMATTED')
       do ISPN = 1,nspn
         read(99,*)
         do irep = 1,n_rep
           read(99,*)
           read(99,*) idum, nstest
           if(nstest.ne.ns_tot(irep)) then
              print *, 'wfind: nstest not equal to ns_tot'
              call stopit
           end if       
           read(99,*) (evalread(i,irep,ispn),i=1,nstest)
         end do
       end do
       close(99)
C
C     --- LOOP OVER SPIN ---
C
      NSYMM=0
      NUSYM=0
      ISPBEG=NSPN
      ISPEND=1
      IF(SPN1)ISPBEG=1
      IF(SPN2)ISPEND=NSPN
C
      NWF=0
      NWFS   (1)=0
      NWFS(NSPN)=0
      DO JSPN=1,NSPN
        DO IREP=1,N_REP
          N_OCC(IREP,JSPN)=0
        END DO
      END DO
C
      DO 240 ISPN=ISPBEG,ISPEND
        DOIT=.FALSE.
        IF (ISPN.EQ.1.AND.SPN1) DOIT=.TRUE.
        IF (ISPN.EQ.2.AND.SPN2) DOIT=.TRUE.
C
        IF(DOIT)THEN
          IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 1'
          CALL OVERLAP(1)
          IF (DEBUG) PRINT *,'WFWIND CALLS OVERLAP MODE: 2'
          CALL OVERLAP(2)
C
C     --- LOOP OVER REPRESENTATIONS ---
C     --- GET MATRIX ELEMENTS       ---
C
          KBAS=0
          DO 130 IREP=1,N_REP
            NBAS=NS_TOT(IREP)
            IF(NBAS.GT.NDH)THEN
              PRINT *,'WFWIND: NDH MUST BE AT LEAST: ',NBAS
              CALL STOPIT
            END IF
            ALLOCATE(AEVAL(NBAS),STAT=IERR)
            IF(IERR.NE.0)THEN
              WRITE(6,*)'wfwind:Error allocating Eval'
            ENDIF 
            ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
            IF(IERR.NE.0)THEN
              WRITE(6,*)'wfwind:Error allocating AHAM'
            ENDIF 
C
C     --- fill in eigenvalues and eigenvectors read from EVALUES and
C     WFOUT
C
            DO i = 1,nbas
              AEVAL(i) = evalread(i,irep,ispn)
              do j = 1,nbas
                AHAM(j,i) = psi_coef(j,i,irep,ispn)
              end do
            end do
C
C     --- DETERMINE WHICH EIGENSTATES ARE NEEDED ---
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
C
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
C
            IF (NEIG.GT.MAX_VIRT_PER_SYM) THEN
              PRINT *,'WFWIND: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',NEIG
              IFAIL=1
              RETURN
            END IF
C
C     --- STORE WAVEFUNCTIONS AND EIGENVALUES ---
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
C
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
            END IF 
            DEALLOCATE(AEVAL,STAT=IERR)
            IF(IERR.NE.0)THEN
              WRITE(6,*)'wfwind:Error deallocating Eval'
            ENDIF 
C
  130     CONTINUE
C
        END IF
  240 CONTINUE
C
       CALL GTTIME(TIME2)
       CALL TIMOUT('WAVEFUNCTIONS IN ENERGY WINDOW:    ',TIME2-TIME1)
C
       RETURN
C
       END
