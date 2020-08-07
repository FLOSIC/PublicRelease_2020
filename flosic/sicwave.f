C UTEP Electronic Structure Lab (2020)
C> @file sicwave.f
       SUBROUTINE SICWAVE(NITER,TRACE)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use global_inputs,only : inbas,iiev,iimesh,iinitial,mpi_io1,
     &        fixm1
       use for_diag1
       use hstor1, only : hstor
       use mixpot1,only : POTIN,POTOUT
       use debug1
       use common2,only : E_UP, E_DN, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : ISTSCF, IHIPOL, PSI_COEF, OCCUPANCY,
     &   N_OCC, PSI, NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT, LDMREP
       use mpidat1,only : NPROC, NCALLED
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:02 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NITER, MAX_TOT, I, IB, IBAS, ICOUNT, IDEG, IEIG, IL,
     & INDREP, IOCC, IOFS, IREP, ISAV, ISORT, ISPFAC, IVIRT, J, JBAS,
     & JNDXSAV, JREP, JSORT, JVIRT, K, KBAS, N_VIRT, NBAS, NDEG, NEIG,
     & NOCCU, NSAV, NTEMP, NTMP, NTMP2, NVIRT_DN, NVIRT_UP, NVIRTTOT
       REAL*8 :: SYMBOL , TRACE, CUTOCC, DELTA, DIAG, DSINGV, ECUT, EF,
     & ELEC, EVALSAV, OCCTMP, SC3, SCENG, SCIS, SCVLSAV, SDL, TEMP,
     & TIME1, TIME2, UMAT
!       INCLUDE 'commons.inc'
! YY moved from commons.inc
       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
     &  EVAL(NDH),SC1(NDH),SC2(NDH)
       save
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       DIMENSION NAMES(3)
       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION JNDXSAV(MAX_TOT*MXSPN)
       DIMENSION SCVLSAV(MAX_TOT*MXSPN)
       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
     &  NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
       DIMENSION DIAG(NDH,MAX_REP)
       DIMENSION NTEMP(MAX_TOT*MXSPN)
!       COMMON/MIXPOT1/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN)
       COMMON/DELSCI/SCIS(200,200)
       COMMON/SCNENG/SCENG
       DIMENSION UMAT(NDH,NDH),DELTA(NDH,2),SC3(NDH,2)

C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/
cccccccccccccccccccccfirst executabel
C
C CHECKING AND SETTING UP SOME STUFF
C
       IF (N_REP.GT.MAX_REP) THEN
        PRINT *,'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       TRACE=0.0D0
       SCENG=0.0D0
**************************************************************
*	Set default mode for calculation of E_F
        EF_MODE=.FALSE.
        FERMISTAT=.TRUE.
C
       OPEN(39,FILE='TMPTRE',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(39)
       READ(39,*,END=10)TEMP
   10  REWIND(39)
       WRITE(39,*)TEMP,' KT IN HARTREES'
       CLOSE(39)
       IF (DEBUG) PRINT*,'NEWWAVE CALLS OVERLAP MODE: 1'
       CALL OVERLAP(1)



       FERMISTAT=.TRUE.
       IF (NSPN.GT.MXSPN) THEN
        PRINT *,'NEWWAVE: MXSPN MUST BE AT LEAST: ',NSPN
        CALL STOPIT
       END IF
C
C REMOVE OLDER EVALXXX FILES IF NITER=1
C
       IF (NITER.EQ.1) call system('rm EVAL0*')
       WRITE(EVALSTR,'(A,I3.3)')'EVAL',NITER
       OPEN(97,FILE='EVALUES',FORM='FORMATTED',STATUS='UNKNOWN')
       OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='UNKNOWN')
 1000  FORMAT(A4)
************************************************************************
C      IF ((FLINE.EQ.'FIXM' ).OR.(FLINE.EQ.'fixm'))
C    &   EF_MODE=.FALSE.
C
C
C
C DIAGONALIZE AND GET OCCUPANCIES
C LOOP OVER SPIN
C
       NOCCU=0
       NWF=0
       PRINT '(A)','CONSTRUCTING NEW WAVEFUNCTIONS'
       CALL GTTIME(TIME1)
**************************TB***********************************
       ELEC=E_UP+E_DN
**************************TB***********************************
       ISPFAC=2/NSPN
       DO 240 ISPN=1,NSPN
        NWFS(ISPN)=0
        WRITE(97,*)'********* NEW TRY ************, SPIN: ',ISPN
        WRITE(98,*)'********* NEW TRY ************, SPIN: ',ISPN

        PRINT '(A,I1,A)','SPIN ',ISPN,':'
        CALL OVERLAP(2)
        CALL OVERLAP(1)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS
C
        KBAS=0
c        NVIRTTOT=0
        DO 130 IREP=1,N_REP
         N_VIRT(IREP,ISPN)=0
         NBAS=NS_TOT(IREP)
         IF (NBAS.GT.NDH) THEN
          PRINT *,'NEWWAVE: NDH MUST BE AT LEAST: ',NBAS
          CALL STOPIT
         END IF

         DO 80 IBAS=1,NBAS
          DO 70 JBAS=IBAS,NBAS
           KBAS=KBAS+1
           OVER(JBAS,IBAS)=HSTOR(KBAS,1)
           OVER(IBAS,JBAS)=HSTOR(KBAS,1)
           HAM (JBAS,IBAS)=HSTOR(KBAS,2)
           HAM (IBAS,JBAS)=HSTOR(KBAS,2)
   70     CONTINUE
   80    CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         CALL MACYS_SCISSOR(NBAS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
           DO I=1,NBAS
           SC3(I,ISPN)=0.0D0
           DO J=1,NBAS
           DO K=1,NBAS
           SC3(I,ISPN)=SC3(I,ISPN)+HAM(J,I)*HAM(K,I)*SCIS(J,K)
           END DO
           END DO
           PRINT*,I,EVAL(I),SC3(I,ISPN),' SCISSOR EXPECTATION'
           END DO
           NEIG=NBAS
         WRITE(97,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
         WRITE(98,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
          N_VIRT(IREP,ISPN)=NEIG
*********************************TB*******************************
          IF(ISPN.EQ.1.AND.IREP.EQ.1) NVIRTTOT=0
*********************************TB*******************************
          DO IEIG=1,NEIG
           EVALSAV(NVIRTTOT+IEIG)=EVAL(IEIG)
           JNDXSAV(NVIRTTOT+IEIG)=IEIG
           SCVLSAV(NVIRTTOT+IEIG)=SC3(IEIG,ISPN)
           NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)/LDMREP(IREP)!MRP 19OCT
           INDREP(NVIRTTOT+IEIG)=IREP
          END DO
           NVIRTTOT=NVIRTTOT+NEIG
          WRITE(97,*)NDMREP(IREP),NEIG
          WRITE(97,1300)(EVAL(IEIG),IEIG=1,NEIG)
          WRITE(98,*)NDMREP(IREP),NEIG
          WRITE(98,1300)(EVAL(IEIG),IEIG=1,NEIG)
 1300     FORMAT(' ',5G15.7)
C
C STORE ALL EIGENVECTORS THAT FIT INTO PSI_COEF
C
         NSAV(IREP,ISPN)=MIN(NEIG,MAX_VIRT_PER_SYM)
         DO 120 ISAV=1,NSAV(IREP,ISPN)
          DO 110 IB=1,NBAS
           PSI_COEF(IB,ISAV,IREP,ISPN)=HAM(IB,ISAV)
  110     CONTINUE
  120    CONTINUE
  130   CONTINUE

C
C FERMI STATISTICS: DEFINE OCCUPATION NUMBERS
C
CJK99 FIXING CASES WITH NEIG.NE.NBAS
C     NEEDED TO BE CHANGED FOR PARALLEL VERSION
C looks messy and probably is never really needed :-(
C only NVIRTTOT has to be right
*********************************TB*******************************
        IF(ISPN.EQ.1) THEN
           NTMP=0
           NTMP2=0
        ELSE
           NTMP=NVIRT_UP
           NTMP2=NVIRT_UP
        ENDIF
*********************************TB*******************************
        DO JREP=1,N_REP-1
         NTMP=NTMP+N_VIRT(JREP,ISPN)
         NTMP2=NTMP2+NS_TOT(JREP)
         IF(NS_TOT(JREP).NE.N_VIRT(JREP,ISPN)) THEN
          DO IL=1,N_VIRT(JREP,ISPN)
           NDEG(NTMP+IL)=NDEG(NTMP2+IL)
           EVALSAV(NTMP+IL)=EVALSAV(NTMP2+IL)
          ENDDO
         ENDIF
        ENDDO
        NVIRTTOT=NTMP+N_VIRT(N_REP,ISPN)
C
********************************TB**************************************

        IF(ISPN.EQ.1) NVIRT_UP=NVIRTTOT
  240  CONTINUE

        IF(EF_MODE) PRINT *, 'OPTIMIZED MOMENT  MODE :'
       IF(.NOT.EF_MODE) THEN
        DO ISPN=1,NSPN
           DO IL=1,MAX_PTS*MXSPN
             POTIN(IL)=0.0D0
             POTOUT(IL)=0.0D0
           ENDDO
           IF (ISPN.EQ.1) THEN
             ELEC=E_UP
             DO IL=1,NVIRT_UP
               POTIN(IL)=EVALSAV(IL)
               NTEMP(IL)=NDEG(IL)
             ENDDO
             CALL FERMILV(NVIRT_UP,ELEC,EF,TEMP,POTIN,POTOUT,NTEMP)
             EFERMI(ISPN)=EF
             DO IL=1,NVIRT_UP
               OCCTMP(IL)=POTOUT(IL)
             ENDDO
           ELSE
             ELEC=E_DN
             IB=0
             DO IL=NVIRT_UP+1,NVIRTTOT
               IB=IB+1
               POTIN(IB)=EVALSAV(IL)
               NTEMP(IB)=NDEG(IL)
c            PRINT *, 'CHECK POINT6', IB, IL, NDEG(IL)
             ENDDO
             NVIRT_DN=NVIRTTOT-NVIRT_UP
             CALL FERMILV(NVIRT_DN,ELEC,EF,TEMP,POTIN,POTOUT,NTEMP)
             EFERMI(ISPN)=EF
             IB=0
             DO IL=NVIRT_UP+1,NVIRTTOT
               IB=IB+1
               OCCTMP(IL)=POTOUT(IB)
             ENDDO
           ENDIF
        END DO !ISPN
       END IF !EF_MODE
*********************************TB*******************************

        DO IL = 1,MAX_PTS*MXSPN
          POTIN(IL)=0.0D0
          POTOUT(IL)=0.0D0
        ENDDO


        ISAV=0
        IOFS=0
        IF (FERMISTAT) THEN
         IF(NSPN.EQ.1)ELEC=E_UP
C        IF(EF_MODE)
C    &     CALL FERMILV(NVIRTTOT,ELEC,EF,TEMP,EVALSAV,OCCTMP,NDEG)
*********************************TB*******************************
        ICOUNT=0

        DO  ISPN=1,NSPN
        IF (EF_MODE) EFERMI(ISPN)=EF
         DO 150 IREP=1,N_REP
          N_OCC(IREP,ISPN)=0
          DO IVIRT=1,N_VIRT(IREP,ISPN)
           ICOUNT=ICOUNT+1
C           IF (OCCTMP(IOFS+IVIRT) .LT. CUTOCC) GOTO 140
           IF (OCCTMP(ICOUNT) .LT. CUTOCC) GOTO 140

           IF (IVIRT .GT. NSAV(IREP,ISPN)) THEN
            ISAV= MAX(ISAV,IVIRT)
           ELSE
            NOCCU=NOCCU+1
            OCCUPANCY(NOCCU)=OCCTMP(ICOUNT)/LDMREP(IREP) ! MRP 19Oct98
            N_OCC(IREP,ISPN)=N_OCC(IREP,ISPN)+1
            NTEMP(ICOUNT)=ISPN
           END IF
  140       CONTINUE
          END DO
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  150    CONTINUE
         IF (ISAV .NE. 0) THEN
          PRINT *,'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',ISAV
          CALL STOPIT
         END IF
         ENDDO
        END IF
C
C GET TRACE AND EVLOCC. FILL N_VIRT WITH THE NUMBER OF
C WAVEFUNCTIONS THAT IS ACTUALLY STORED IN PSI_COEF (NSAV -> N_VIRT)
C
        ISAV=0
        IOFS=0
        JVIRT=0

*************************************TB*******************************
        SDL=0.0D0
        DO ISPN=1,NSPN

        ELEC=0.0D0
*************************************TB*******************************
        DO 200 IREP=1,N_REP
         DO 190 IVIRT=1,N_OCC(IREP,ISPN)
          JVIRT=IOFS+IVIRT
c          JVIRT=JVIRT+1
          TRACE=TRACE+EVALSAV(JVIRT)*OCCTMP(JVIRT)*NDEG(JVIRT)
          SCENG=SCENG+SCVLSAV(JVIRT)*OCCTMP(JVIRT)*NDEG(JVIRT)
*************************************TB*******************************
          ELEC=ELEC+OCCTMP(JVIRT)*NDEG(JVIRT)
          IF(OCCTMP(JVIRT).GT.0.0001)THEN
          PRINT 1401,ISPN,IREP,IVIRT,
     &OCCTMP(JVIRT),EVALSAV(JVIRT),SCVLSAV(JVIRT),SCENG
 1401     FORMAT(' SCENG FLAG:',3I4,4F12.4)
          END IF
          DO 180 IDEG=1,NDMREP(IREP)
           NWF=NWF+1
           NWFS(ISPN)=NWFS(ISPN)+1
           IF (NWF .GT. MAX_OCC) THEN
            ISAV=NWF
           ELSE
            EVLOCC(NWF)=EVALSAV(JVIRT)
           END IF
  180     CONTINUE
          IF (ISAV .EQ. 0) THEN
           PRINT 1400,NWF,EVLOCC(NWF),OCCTMP(JVIRT)*NDEG(JVIRT)*ISPFAC
 1400      FORMAT('STATE ',I4,', EV= ',F15.6,', OCCUP= ',F12.6)
          END IF
  190    CONTINUE
         IOFS=IOFS+N_VIRT(IREP,ISPN)
         N_VIRT(IREP,ISPN)=NSAV(IREP,ISPN)
  200   CONTINUE
*************************************TB*******************************
        IF(ISPN.EQ.1) THEN
           E_UP=ELEC
        ELSE
           E_DN=ELEC
        ENDIF
*************************************TB*******************************
        IF (ISAV .NE. 0) THEN
         PRINT *,'NEWWAVE: MAX_OCC MUST BE AT LEAST ',ISAV
         CALL STOPIT
        END IF
       ENDDO
        IF(NSPN.EQ.1.AND.EF_MODE) THEN
           E_DN=E_UP
        ENDIF
       Write(6,*)'ELECTRON :',E_UP,E_DN
C
C SORT EIGENVALUES
C
        DO ISORT=1,NVIRTTOT
         PRINT 1499,ISORT,EVALSAV(ISORT),NTEMP(ISORT),
     &    INDREP(ISORT),NDEG(ISORT),JNDXSAV(ISORT)
        END DO
        DO 220 ISORT=1,NVIRTTOT
 1499   FORMAT(' ',I3,F12.4,4I4,' <LOOK AT THIS')
         DO 210 JSORT=ISORT+1,NVIRTTOT
          IF (EVALSAV(JSORT).LT.EVALSAV(ISORT)) THEN
           CALL ISWAP(JNDXSAV(ISORT),JNDXSAV(JSORT))
           CALL SWAP(EVALSAV(ISORT),EVALSAV(JSORT))
           CALL SWAP(OCCTMP(ISORT),OCCTMP(JSORT))
           CALL ISWAP(NDEG(ISORT),NDEG(JSORT))
           CALL ISWAP(INDREP(ISORT),INDREP(JSORT))
           CALL ISWAP(NTEMP(ISORT),NTEMP(JSORT))
          END IF
  210    CONTINUE
  220   CONTINUE
C
C OUTPUT
C

        WRITE(6,*)'ELECTRONS OF SPIN UP : ', E_UP
        WRITE(6,*)'ELECTRONS OF SPIN DN : ', E_DN
        WRITE(97,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(97,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        WRITE(98,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(98,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        ECUT=MAX(0.0D0,EF+2.0D1*TEMP)
        DO 230 IOCC=1,NVIRTTOT
         WRITE(97,1500) IOCC,NTEMP(IOCC),JNDXSAV(IOCC),INDREP(IOCC),
     &                  NDEG(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         WRITE(98,1500) IOCC,NTEMP(IOCC),JNDXSAV(IOCC),INDREP(IOCC),
     &                  NDEG(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         IF (EVALSAV(IOCC).GT.ECUT) GOTO 290
  230   CONTINUE
  290  CONTINUE
************************************************************************
 1500  FORMAT(I5,2X,'SPIN/ORB:',I1,'-', I3.3, 2X,
     &'REP: ',I2,2X,'DEG: ',I2,2X,
     &'ENERGY: ',  G14.6,2X,'OCC: ',G14.6)
       ISPFAC=2/NSPN
       TRACE=TRACE*ISPFAC
       CLOSE(97)
       CLOSE(98)
       CALL GTTIME(TIME2)
       CALL TIMOUT('CONSTRUCTION OF NEW WAVEFUNCTIONS: ',TIME2-TIME1)
       RETURN
       END
