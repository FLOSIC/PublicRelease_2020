C
C
C **************************************************************
C
C NEWWAVE VERSION DIRK POREZAG AUGUST 1994 
C PARALLEL VERSION JENS KORTUS SEPTEMBER 1999 
C NEWWAVE DIAGONALIZES THE HAMILTONIAN AND CALCULATES THE 
C OCCUPATION NUMBERS
C
       SUBROUTINE NEWWAVE(NITER,TRACE)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
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

!      INCLUDE 'PARAMS'
       INCLUDE 'PARAMA2'
       !INCLUDE 'commons.inc'
       INTEGER :: NITER, MAX_TOT, I_REP, IB, IBAS, ICOUNT, IDEG, IEIG,
     & IERR, IL, INITIAL, IOCC, IOFS, IR, IREP, ISAV, ISORT, ISPFAC,
     & ITER, IVIRT, JBAS, JREP, JSORT, JVIRT, KBAS, L, N_VIRT, NBAS,
     & NEIG, NHTOT, NOCCU, NSAV, NTMP, NTMP2, NVIRT_DN, NVIRT_UP,
     & NVIRTTOT
       INTEGER :: IPTS,MPTS,MSPN
       REAL*8 :: SYMBOL , TRACE, AUX, CUTOCC, DSINGV, ECUT, EF, ELEC,
     & EVALSAV, TEMP, TIME1, TIME2, TIME_A, TIME_B
       SAVE
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       DIMENSION NAMES(3)
C      DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION EVALSAV(MAX_TOT*MXSPN)
C       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
C     &  NSAV(MAX_REP,MXSPN)
       DIMENSION NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
C       DIMENSION DIAG(NDH,MAX_REP)
C       DIMENSION NTEMP(MAX_TOT*MXSPN)
       INTEGER :: NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
     &  NTEMP(MAX_TOT*MXSPN)                                   
       REAL*8 :: OCCTMP(MAX_TOT*MXSPN) 
C       COMMON/MIXPOT1/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN) 
C       COMMON/SICSHIFT/SICSCI(NDH_TOT,2)
! additional 
       INTEGER :: J,KSPN,imove
       REAL*8 :: SICSCI(NDH_TOT,2)
C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/ 
       SICSCI=0.0D0
       INQUIRE(FILE='SICSCI',EXIST=EXIST)
       IF(EXIST)THEN
       OPEN(50,FILE='SICSCI',FORM='UNFORMATTED')
       DO KSPN=1,NSPN
       READ(50)(SICSCI(J,KSPN),J=1,NDH_TOT)
       END DO
       CLOSE(50)
       END IF
C
C CHECKING AND SETTING UP SOME STUFF
C
       IF (N_REP.GT.MAX_REP) THEN
        PRINT *,'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       TRACE=0.0D0
**************************************************************
*	Set default mode for calculation of E_F
        EF_MODE=.TRUE.
        IF(ISTSCF.EQ.3.AND.NITER.EQ.1) EF_MODE=.FALSE.
**************************************************************
C
C READ IN TEMPERATURE AND CALL OVERLAP
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
       IF (NITER.EQ.1) THEN
        ITER=0
   15    ITER=ITER+1
         WRITE(EVALSTR,'(A,I3.3)')'EVAL',ITER 
         INQUIRE(FILE=EVALSTR,EXIST=EXIST)
         IF (EXIST) THEN
          OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='OLD')
          CLOSE(98,STATUS='DELETE')
          GOTO 15
         END IF  
        CONTINUE
       END IF
C
C CHECK IF FROZEN OCCUPATION MODE  
C
       WRITE(EVALSTR,'(A,I3.3)')'EVAL',NITER 
       OPEN(97,FILE='EVALUES',FORM='FORMATTED',STATUS='UNKNOWN')
       OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='UNKNOWN')
       READ(97,1000,END=60,ERR=60)FLINE
 1000  FORMAT(A4)
************************************************************************
       IF ((FLINE.EQ.'FIXM' ).OR.(FLINE.EQ.'fixm')) 
     &   EF_MODE=.FALSE.
       IF ((FLINE.EQ.'OCCU').OR.(FLINE.EQ.'occu')) THEN
             CLOSE(97)
             CLOSE(98)
        CALL NEWOCC(NITER,TRACE)
        RETURN
C       FERMISTAT=.FALSE.
C       IOCC=0
C       DO 30 ISPN=1,NSPN
C        DO 20 IREP=1,N_REP
C         READ(97,*)N_OCC(IREP,ISPN)
C         IF (N_OCC(IREP,ISPN).GT.MAX_VIRT_PER_SYM) THEN
C          PRINT *,'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',
C    &              N_OCC(IREP,ISPN)
C          CALL STOPIT
C         END IF
C         IF (N_OCC(IREP,ISPN).GT.0) THEN
C          READ(97,*)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
C          IOCC=IOCC+N_OCC(IREP,ISPN)
C         END IF
C  20    CONTINUE
C  30   CONTINUE
C       REWIND(97)
C       WRITE(97,1100)FLINE
C       WRITE(98,1100)FLINE
 1100   FORMAT(A4)
C       IOCC=0
C       DO 50 ISPN=1,NSPN
C        DO 40 IREP=1,N_REP
C         WRITE(97,*)N_OCC(IREP,ISPN)
C         WRITE(97,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
C         WRITE(98,*)N_OCC(IREP,ISPN)
C         WRITE(98,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
C1200     FORMAT(' ',5G15.7)
C         IOCC=IOCC+N_OCC(IREP,ISPN)
C  40    CONTINUE
C  50   CONTINUE
C       GOTO 65
       END IF
   60  REWIND(97)
   65  CONTINUE

       
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
       ELEC=E_UP+E_DN
       ISPFAC=2/NSPN
        NVIRTTOT=0
       DO 240 ISPN=1,NSPN
        NWFS(ISPN)=0
        IF(FLINE.EQ.'FIXM' .OR. FLINE.EQ.'fixm') THEN
          WRITE(97,1100)FLINE
          WRITE(6,1100)FLINE
        ENDIF
        WRITE(97,*)'********* NEW TRY ************, SPIN: ',ISPN
        WRITE(98,*)'********* NEW TRY ************, SPIN: ',ISPN

        PRINT '(A,I1,A)','SPIN ',ISPN,':'
            IF (DEBUG) PRINT *,'NEWWAVE CALLS OVERLAP MODE: 2'
            CALL OVERLAP(2)
        IF (DEBUG) PRINT *,'NEWWAVE CALLS OVERLAP MODE: 1'
        CALL OVERLAP(1)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS
C
        KBAS=0
#ifdef MPI
        NLOOP=N_REP/NPROC+1
        MRUN=NPROC
        DO ILOOP=1,NLOOP
         IF(ILOOP.EQ.NLOOP) MRUN=N_REP-(ILOOP-1)*MRUN 
C
        DO IREPP=1,MRUN
         IREP=IREPP+(ILOOP-1)*NPROC
#else
        DO 130 IREP=1,N_REP
#endif
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
           HAM (JBAS,IBAS)=HSTOR(KBAS,2)-SICSCI(KBAS,ISPN)
   70     CONTINUE
   80    CONTINUE

CJK11/2000
         CALL SCISSOR(IREP)
CJK11/2000

C
C GET EIGENVECTORS AND EIGENVALUES
C IF WE START A NEW GEOMETRY FROM AN OLD HAMILTONIAN, USE SINGULAR
C VALUE DECOMPOSITION TO AVOID PROBLEMS WITH SCREWED UP EIGENVALUES
C USE N_VIRT TO TEMPORARILY STORE THE NUMBER OF AVAILABLE EIGENSTATES 
C
         IF (NBAS.NE.0) THEN
#ifdef MPI
         ISVD=1
         IF ((NITER.EQ.1).AND.(ISTSCF.EQ.1).AND.(IHIPOL.EQ.1))ISVD=0 
          CALL PAMWAVE(1,ISVD,IREP,NBAS) 
#else
         IF (NBAS.NE.0) THEN
          IF ((NITER.EQ.1).AND.(ISTSCF.EQ.1).AND.(IHIPOL.EQ.1)) THEN
           PRINT '(A)','NEWWAVE: USING SINGULAR VALUE DECOMPOSITION'
           CALL DIAGSVD(NDH,NBAS,NEIG,HAM,OVER,EVAL,
     &                  SC1,SC2,DSINGV,1)
          ELSE
           CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
           NEIG=NBAS
          END IF
#endif
         ENDIF
C
#ifdef MPI
         ENDDO          ! IREPP
C
         DO IREPP=1,MRUN
         IF (NCALLED.NE.0) THEN
          CALL PAMWAVE(2,ISVD,IREP,NBAS)
          NEIG=NBAS
#endif         
         WRITE(97,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
         WRITE(98,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
          N_VIRT(IREP,ISPN)=NEIG
*********************************TB*******************************
C          IF(ISPN.EQ.1.AND.IREP.EQ.1) NVIRTTOT=0
*********************************TB*******************************
c         DO JREP=1,IREP-1
c          NVIRTTOT=NVIRTTOT+NS_TOT(JREP)
c         ENDDO
#ifdef MPI
         print*,'nvirt : ',ispn,irep,nvirttot,neig,N_VIRT(IREP,ISPN)
          imove=0
          do jrep=1,irep-1
           imove=imove+ns_tot(jrep)
          enddo
          if(ispn.eq.2) then
          do jrep=1,N_rep
           imove=imove+ns_tot(jrep)
          enddo
          endif
          print*,'imove : ',imove,nvirttot,ispn,irep
#else
          imove=nvirttot
#endif
          DO IEIG=1,NEIG
           EVALSAV(imove+IEIG)=EVAL(IEIG) 
C           EVALSAV(NVIRTTOT+IEIG)=EVAL(IEIG) 
C          NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)
C           NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)/LDMREP(IREP)!MRP 19OCT
           NDEG(imove+IEIG)=NDMREP(IREP)/LDMREP(IREP)!MRP 19OCT
C           INDREP(NVIRTTOT+IEIG)=IREP
           INDREP(imove+IEIG)=IREP
          END DO
           NVIRTTOT=NVIRTTOT+NEIG
          WRITE(97,*)NDMREP(IREP),NEIG
          WRITE(97,1300)(EVAL(IEIG),IEIG=1,NEIG)
          WRITE(98,*)NDMREP(IREP),NEIG
          WRITE(98,1300)(EVAL(IEIG),IEIG=1,NEIG)
 1300     FORMAT(' ',5G15.7)
         END IF
C
C STORE ALL EIGENVECTORS THAT FIT INTO PSI_COEF
C
         NSAV(IREP,ISPN)=MIN(NEIG,MAX_VIRT_PER_SYM)
         DO 120 ISAV=1,NSAV(IREP,ISPN)
          DO 110 IB=1,NBAS
           PSI_COEF(IB,ISAV,IREP,ISPN)=HAM(IB,ISAV)
  110     CONTINUE
  120    CONTINUE
#ifdef MPI
         END DO   !IREPP
         END DO   !ILOOP
#else
  130   CONTINUE
#endif

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
       print*,'NTMP : ',ispn,ntmp,N_VIRT(N_REP,ISPN),NVIRTTOT
C
********************************TB**************************************

        IF(ISPN.EQ.1) NVIRT_UP=NVIRTTOT
  240  CONTINUE    

        IF(EF_MODE) PRINT *, 'OPTIMIZED MOMENT  MODE :'
       IF(.NOT.EF_MODE) THEN
        print*,'entering .NOT.EF_MODE :',NVIRT_UP,ispn
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
             print*,'fermi up : ',nvirt_up,elec,ef
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
             ENDDO
             NVIRT_DN=NVIRTTOT-NVIRT_UP
             CALL FERMILV(NVIRT_DN,ELEC,EF,TEMP,POTIN,POTOUT,NTEMP)
             EFERMI(ISPN)=EF
             print*,'fermi dn : ',nvirt_dn,elec,ef
             IB=0
             DO IL=NVIRT_UP+1,NVIRTTOT
               IB=IB+1
               OCCTMP(IL)=POTOUT(IB)
             ENDDO

            ENDIF
          ENDDO 
        ENDIF 
*********************************TB*******************************

        DO IL = 1,MAX_PTS*MXSPN
          POTIN(IL)=0.0D0
          POTOUT(IL)=0.0D0
        ENDDO


        ISAV=0
        IOFS=0
        IF (FERMISTAT) THEN
         IF(NSPN.EQ.1)ELEC=E_UP
         IF(EF_MODE) then
           print*,'should not be called for fixm'
           CALL FERMILV(NVIRTTOT,ELEC,EF,TEMP,EVALSAV,OCCTMP,NDEG)
         endif
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
C           OCCUPANCY(NOCCU)=OCCTMP(IOFS+IVIRT)                
C            OCCUPANCY(NOCCU)=OCCTMP(IOFS+IVIRT)/LDMREP(IREP) ! MRP 19Oct98
            OCCUPANCY(NOCCU)=OCCTMP(ICOUNT)/LDMREP(IREP) ! MRP 19Oct98
            N_OCC(IREP,ISPN)=N_OCC(IREP,ISPN)+1
            NTEMP(ICOUNT)=ISPN
           END IF
  140       CONTINUE
          END DO
c  140     CONTINUE
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  150    CONTINUE
          print*,'BUG suspected:', icount,noccu
           print*,'N_VIRT : ',ispn, (N_VIRT(IREP,ISPN),IREP=1,N_REP)
           print*,'N_SAV  : ',ispn, (NSAV (IREP,ISPN),IREP=1,N_REP)
           print*,'N_OCC  : ',ispn, (N_OCC(IREP,ISPN),IREP=1,N_REP)

         IF (ISAV .NE. 0) THEN
          PRINT *,'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',ISAV
          CALL STOPIT
         END IF
         ENDDO

        ELSE
C
C FROZEN OCCUPATION MODE. N_OCC IS ALREADY DEFINED
C
        DO  ISPN=1,NSPN
         EFERMI(ISPN)= -1.0D20
         DO IVIRT=1,NVIRTTOT
          OCCTMP(IVIRT)=0.0D0
         END DO
         DO 170 IREP=1,N_REP
          IF (N_OCC(IREP,ISPN) .GT. NSAV(IREP,ISPN)) THEN
           ISAV= IREP
          ELSE
           DO 160 IVIRT=1,N_OCC(IREP,ISPN)
            NOCCU=NOCCU+1
c            OCCTMP(IOFS+IVIRT)=OCCUPANCY(NOCCU)
            OCCTMP(NOCCU)=OCCUPANCY(NOCCU)
            IF (OCCUPANCY(NOCCU) .GT. CUTOCC) THEN
C             EFERMI(ISPN)=MAX(EFERMI(ISPN),EVALSAV(IOFS+IVIRT))
             EFERMI(ISPN)=MAX(EFERMI(ISPN),EVALSAV(NOCCU))
            END IF
  160      CONTINUE
          END IF
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  170    CONTINUE
         IF (ISAV .NE. 0) THEN
          PRINT *,'NEWWAVE: NOT ENOUGH STATES FOR GIVEN OCCUPATION'
          PRINT *,'         IN REPRESENTATION: ',ISAV
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
        DO ISPN=1,NSPN

        ELEC=0.0D0
*************************************TB*******************************
        DO 200 IREP=1,N_REP
         DO 190 IVIRT=1,N_OCC(IREP,ISPN)
          JVIRT=IOFS+IVIRT
c          JVIRT=JVIRT+1
          TRACE=TRACE+EVALSAV(JVIRT)*OCCTMP(JVIRT)*NDEG(JVIRT)
*************************************TB*******************************
          ELEC=ELEC+OCCTMP(JVIRT)*NDEG(JVIRT)
*************************************TB*******************************
C         DO 180 IDEG=1,NDEG(JVIRT)
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
        DO 220 ISORT=1,NVIRTTOT
         DO 210 JSORT=ISORT+1,NVIRTTOT
          IF (EVALSAV(JSORT).LT.EVALSAV(ISORT)) THEN
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
        WRITE(6,*)'ELECTRONS OF SPIN UP : ', E_UP
        WRITE(6,*)'ELECTRONS OF SPIN DN : ', E_DN
        WRITE(97,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(97,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        WRITE(98,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(98,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        ECUT=MAX(0.0D0,EF+2.0D1*TEMP)
        DO 230 IOCC=1,NVIRTTOT
         WRITE(97,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         WRITE(98,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         IF (EVALSAV(IOCC).GT.ECUT) GOTO 290
  230   CONTINUE
  290  CONTINUE
************************************************************************
 1500  FORMAT(I5,2X,'REP: ',I2,2X,'DEG: ',I2,2X,'SPIN: ',I2,2X,
     &'ENERGY: ',  G14.6,2X,'OCC: ',G14.6)
       ISPFAC=2/NSPN
       TRACE=TRACE*ISPFAC
       CLOSE(97)
       CLOSE(98)
       CALL GTTIME(TIME2)
       CALL TIMOUT('CONSTRUCTION OF NEW WAVEFUNCTIONS: ',TIME2-TIME1)
       RETURN
       END

!***********************************************************************

       SUBROUTINE NEWOCC(NITER,TRACE)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use for_diag1
       use hstor1, only : hstor
       use mixpot1,only : POTIN,POTOUT
       use debug1
       use common2,only : E_UP, E_DN, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : ISTSCF, IHIPOL, PSI_COEF, OCCUPANCY,
     &   N_OCC, PSI, NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT, LDMREP
       INCLUDE 'PARAMA2'
       !INCLUDE 'PARAMS'
       !INCLUDE 'commons.inc'
       !PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM
       INTEGER :: IB,IBAS,IDEG,IEIG,INDREP,IOCC,IOFS,IREP,ISAV,ISORT,
     &            ISPFAC,IVIRT,JBAS,JSORT,JVIRT,KBAS,KVIRT,
     &            L,MAX_TOT,NBAS,NDEG,NEIG,NITER,NOCCU,NSAV,NTEMP,
     &            N_VIRT,NVIRT_UP,NVIRTTOT,OCCTMP
       REAL*8 :: CUTOCC,DIAG,DSINGV,ECUT,EF,ELEC,ESPN,
     &           EVALSAV,TEMP,TIME1,TIME2,TRACE
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       DIMENSION NAMES(3)
       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
     &  NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
       DIMENSION DIAG(NDH,MAX_REP)
       DIMENSION NTEMP(MAX_TOT*MXSPN)
!       COMMON/MIXPOT1/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN)
C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/
       IF (N_REP.GT.MAX_REP) THEN
        PRINT *,'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       TRACE=0.0D0
       CALL OVERLAP(1)
       IF (NSPN.GT.MXSPN) THEN
        PRINT *,'NEWWAVE: MXSPN MUST BE AT LEAST: ',NSPN
        CALL STOPIT
       END IF
C
C CHECK IF FROZEN OCCUPATION MODE
C
       WRITE(EVALSTR,'(A,I3.3)')'EVAL',NITER
       OPEN(97,FILE='EVALUES',FORM='FORMATTED',STATUS='UNKNOWN')
       OPEN(98,FILE=EVALSTR,FORM='FORMATTED',STATUS='UNKNOWN')
       READ(97,1000,END=60,ERR=60)FLINE
 1000  FORMAT(A4)
        FERMISTAT=.FALSE.
        IOCC=0
        DO 30 ISPN=1,NSPN
         DO 20 IREP=1,N_REP
          READ(97,*)N_OCC(IREP,ISPN)
          IF (N_OCC(IREP,ISPN).GT.MAX_VIRT_PER_SYM) THEN
           PRINT *,'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',
     &              N_OCC(IREP,ISPN)
           CALL STOPIT
          END IF
          IF (N_OCC(IREP,ISPN).GT.0) THEN
           READ(97,*)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
           IOCC=IOCC+N_OCC(IREP,ISPN)
          END IF
   20    CONTINUE
   30   CONTINUE
        REWIND(97)
        WRITE(97,1100)FLINE
        WRITE(98,1100)FLINE
 1100   FORMAT(A4)
        IOCC=0
        DO 50 ISPN=1,NSPN
        ESPN=0.0D0
         DO 40 IREP=1,N_REP
          WRITE(97,*)N_OCC(IREP,ISPN)
          WRITE(97,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
          WRITE(98,*)N_OCC(IREP,ISPN)
          WRITE(98,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
                    DO L=IOCC+1,IOCC+N_OCC(IREP,ISPN)
                    ESPN=ESPN+OCCUPANCY(L)*NDMREP(IREP)
                    END DO
 1200     FORMAT(' ',5G15.7)
          IOCC=IOCC+N_OCC(IREP,ISPN)
   40    CONTINUE
              IF(ISPN.EQ.1)E_UP=ESPN
              IF(ISPN.EQ.2)E_DN=ESPN
   50   CONTINUE
        GOTO 65
   60  REWIND(97)
   65  CONTINUE
       NOCCU=0
       NWF=0
       PRINT '(A)','CONSTRUCTING NEW WAVEFUNCTIONS'
       CALL GTTIME(TIME1)
       ELEC=E_UP+E_DN
       ISPFAC=2/NSPN
        NVIRTTOT=0
       L=0
       DO 240 ISPN=1,NSPN
        NWFS(ISPN)=0
        IF(FLINE.EQ.'FIXM' .OR. FLINE.EQ.'fixm') THEN
          WRITE(97,1100)FLINE
          WRITE(6,1100)FLINE
        ENDIF
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
           HAM (JBAS,IBAS)=HSTOR(KBAS,2)
   70     CONTINUE
   80    CONTINUE
         IF (NBAS.NE.0) THEN
           CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
           NEIG=NBAS
C
         WRITE(97,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
         WRITE(98,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
          N_VIRT(IREP,ISPN)=NEIG
          NVIRTTOT=NVIRTTOT+NEIG
          WRITE(97,*)NDMREP(IREP),NEIG
          WRITE(97,1300)(EVAL(IEIG),IEIG=1,NEIG)
          WRITE(98,*)NDMREP(IREP),NEIG
          WRITE(98,1300)(EVAL(IEIG),IEIG=1,NEIG)
 1300     FORMAT(' ',5G15.7)
         END IF
C
C STORE ALL EIGENVECTORS THAT FIT INTO PSI_COEF
C
         NSAV(IREP,ISPN)=MIN(NEIG,MAX_VIRT_PER_SYM)
         DO 120 ISAV=1,NSAV(IREP,ISPN)
           IF(ISAV.LE.N_OCC(IREP,ISPN))THEN
              L=L+1
              TRACE=TRACE+OCCUPANCY(L)*EVAL(ISAV)*NDMREP(IREP)
           EVALSAV(L)=EVAL(ISAV)
           OCCTMP (L)=OCCUPANCY(L)       
           NDEG   (L)=NDMREP(IREP)
           INDREP (L)=IREP           
           NTEMP  (L)=ISPN           
           END IF
          DO 110 IB=1,NBAS
           PSI_COEF(IB,ISAV,IREP,ISPN)=HAM(IB,ISAV)
  110     CONTINUE
  120    CONTINUE
  130   CONTINUE
        IF(ISPN.EQ.1) NVIRT_UP=NVIRTTOT
  240  CONTINUE   
        ISAV=0
        IOFS=0
        JVIRT=0
        KVIRT=0
        DO ISPN=1,NSPN
        ELEC=0.0D0
        DO 200 IREP=1,N_REP
         DO 190 IVIRT=1,N_OCC(IREP,ISPN)
          KVIRT=KVIRT+1
          ELEC=ELEC+OCCTMP(KVIRT)*NDEG(KVIRT)
          DO 180 IDEG=1,NDMREP(IREP)
           NWF=NWF+1
           NWFS(ISPN)=NWFS(ISPN)+1
           IF (NWF .GT. MAX_OCC) THEN
            CALL STOPIT
           ELSE
            EVLOCC(NWF)=EVALSAV(KVIRT)
           END IF
  180     CONTINUE
  190    CONTINUE
         N_VIRT(IREP,ISPN)=NSAV(IREP,ISPN)
  200   CONTINUE
        IF (ISAV .NE. 0) THEN
         PRINT *,'NEWWAVE: MAX_OCC MUST BE AT LEAST ',ISAV
         CALL STOPIT
        END IF
       ENDDO
       Write(6,*)'ELECTRON :',E_UP,E_DN
C
C SORT EIGENVALUES
C
        DO 220 ISORT=1,L           
         DO 210 JSORT=ISORT+1,L          
          IF (EVALSAV(JSORT).LT.EVALSAV(ISORT)) THEN
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
        DO 230 IOCC=1,L           
         WRITE(97,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         WRITE(98,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                  EVALSAV(IOCC),OCCTMP(IOCC)
         IF (EVALSAV(IOCC).GT.ECUT) GOTO 290
  230   CONTINUE
  290  CONTINUE
 1500  FORMAT(I5,2X,'REP: ',I2,2X,'DEG: ',I2,2X,'SPIN: ',I2,2X,
     &'ENERGY: ',  G14.6,2X,'OCC: ',G14.6)
       ISPFAC=2/NSPN
       TRACE=TRACE*ISPFAC
       CLOSE(97)
       CLOSE(98)
       CALL GTTIME(TIME2)
       CALL TIMOUT('CONSTRUCTION OF NEW WAVEFUNCTIONS: ',TIME2-TIME1)
       RETURN
       END
