C UTEP Electronic Structure Lab (2020)
C
C
C **************************************************************
C NEWWAVE VERSION DIRK POREZAG AUGUST 1994
C PARALLEL VERSION JENS KORTUS SEPTEMBER 1999 
C NEWWAVE DIAGONALIZES THE HAMILTONIAN AND CALCULATES THE 
C OCCUPATION NUMBERS
C
       SUBROUTINE NEWWAVE_serial(NITER,TRACE)
C
C WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use global_inputs,only : inbas,iiev,iimesh,iinitial,mpi_io1,
     &        fixm1,mixing1
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
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
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
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM,sichammix
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
c       CHARACTER*7 NAMES
c       DIMENSION NAMES(3)
C       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION EVALSAV(MAX_TOT*MXSPN)
C       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
C     &  NSAV(MAX_REP,MXSPN)
       DIMENSION NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
C       DIMENSION DIAG(NDH,MAX_REP)
C       DIMENSION NTEMP(MAX_TOT*MXSPN)
       INTEGER,ALLOCATABLE :: NDEG(:),INDREP(:),NTEMP(:)
       REAL*8, ALLOCATABLE :: OCCTMP(:)
C       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN) 
! additional
       INTEGER :: J,KSPN,ISUCCEED
       REAL*8 :: SICSCI(NDH_TOT,2)
       character*10 :: control_wave
C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
c       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/

!      hammixing option H_DFA only or H_DFA+SIC
       sichammix=.true. 
!reading SICSCI  !if sichammix=.true., add it in the sichammixdrv
       CALL CHECK_INPUTS
       if(sichammix .and. (MIXING1.EQ.2)) then
        continue !do nothing
       else
        SICSCI=0.0D0
        INQUIRE(FILE='SICSCI',EXIST=EXIST)
        IF(EXIST)THEN
         OPEN(50,FILE='SICSCI',FORM='UNFORMATTED')
         DO KSPN=1,NSPN
          READ(50)(SICSCI(J,KSPN),J=1,NDH_TOT)
         END DO
         CLOSE(50)
        END IF
       control_wave='NO-FLIP'
        ISUCCEED=1
        IF(MAX_VIRT_PER_SYM.NE.NDH)THEN
        PRINT*,'IF USING NEWWAVE, MAX_VIRT_PER_SYM MUST BE:',NDH
        ISUCCEED=0
        RETURN
        END IF 
       end if
!end reading 
C
C CHECKING AND SETTING UP SOME STUFF
C
       IF (N_REP.GT.MAX_REP) THEN
        write(6,*)'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       TRACE=0.0D0

C       Set default mode for calculation of E_F

        EF_MODE=.TRUE.
        IF(ISTSCF.EQ.3.AND.NITER.EQ.1) EF_MODE=.FALSE.
C
C READ IN TEMPERATURE AND CALL OVERLAP
C
       OPEN(39,FILE='TMPTRE',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(39)
       READ(39,*,END=10)TEMP
   10  REWIND(39)
       WRITE(39,*)TEMP,' KT IN HARTREES'
       CLOSE(39)
       IF (DEBUG) write(6,*)'NEWWAVE CALLS OVERLAP MODE: 1'
C       write(6,*)'NEWWAVE CALLS OVERLAP MODE: 1'
       CALL OVERLAP(1)



       FERMISTAT=.TRUE.
       IF (NSPN.GT.MXSPN) THEN
        write(6,*)'NEWWAVE: MXSPN MUST BE AT LEAST: ',NSPN
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
       CALL CHECK_INPUTS
       IF(FIXM1) EF_MODE=.FALSE.
       IF ((FLINE.EQ.'FIXM' ).OR.(FLINE.EQ.'fixm')) 
     &   EF_MODE=.FALSE.
       IF ((FLINE.EQ.'OCCU').OR.(FLINE.EQ.'occu')) THEN
        FERMISTAT=.FALSE.
        EF_MODE=.FALSE.
        IOCC=0
        DO 30 ISPN=1,NSPN
         DO 20 IREP=1,N_REP
          READ(97,*)N_OCC(IREP,ISPN)
          IF (N_OCC(IREP,ISPN).GT.MAX_VIRT_PER_SYM) THEN
           write(6,*)'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',
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
         DO 40 IREP=1,N_REP
          WRITE(97,*)N_OCC(IREP,ISPN)
          WRITE(97,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
          WRITE(98,*)N_OCC(IREP,ISPN)
          WRITE(98,1200)(OCCUPANCY(L),L=IOCC+1,IOCC+N_OCC(IREP,ISPN))
 1200     FORMAT(' ',5G15.7)
          IOCC=IOCC+N_OCC(IREP,ISPN)
   40    CONTINUE
   50   CONTINUE
        GOTO 65
       END IF
   60  REWIND(97)
   65  CONTINUE

       NHTOT=0
       DO ISPN=1,NSPN
          DO I_REP=1,N_REP
           NHTOT=NHTOT+NS_TOT(I_REP)
          END DO
       END DO
C       write(6,*)'newwave:NHTOT=',NHTOT,NSPN,N_REP,NS_TOT(1) 

       NHTOT=NHTOT*MXSPN
       ALLOCATE(NDEG(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating NDEG'
       ENDIF
       ALLOCATE(INDREP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating INDREP'
       ENDIF
       ALLOCATE(NTEMP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating NTEMP'
       ENDIF
       ALLOCATE(OCCTMP(NHTOT),STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Allocating OCCTMP'
       ENDIF

C       write(6,*)'newwave:NHTOT=',NHTOT
C       write(6,*)'newwave:MAX_TOT*MXSPN=',MAX_TOT*MXSPN

C       DO IR=1,MAX_TOT*MXSPN
       DO IR=1,NHTOT
           OCCTMP(IR)=0.0D0
       END DO
C
C
C CHECK TO SEE IF WE ARE DOING HAMILTONIAN MIXING
C
       CALL CHECK_INPUTS
       IF(MIXING1.EQ.2)THEN
         DO ISPN=1,NSPN
           !YY. Calculate hamiltonian. This will create HAMOLD. 
           CALL OVERLAP(2)
         END DO
         OPEN(99,FILE='HAMOLD', FORM='UNFORMATTED',STATUS='UNKNOWN')
         READ(99)MPTS,MSPN
         CLOSE(99)
         IF(sichammix) THEN
          CALL SICHAMMIXDRV(MPTS)
         ELSE
          CALL HAMMIXDRV(MPTS)
         END IF
       ENDIF
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
       write(6,*)'NSPN=',NSPN
       ISPFAC=2/NSPN
       DO 240 ISPN=1,NSPN
        NWFS(ISPN)=0
c        IF (ISPN.EQ.1) ELEC=E_UP
c        IF (ISPN.EQ.2) ELEC=E_DN
!YY uncommented 4 lines below for backward compatibility
        IF(FLINE.EQ.'FIXM' .OR. FLINE.EQ.'fixm') THEN
          WRITE(97,1100)FLINE
          WRITE(6,1100)FLINE
        ENDIF
        WRITE(97,*)'********* NEW TRY ************, SPIN: ',ISPN
        WRITE(98,*)'********* NEW TRY ************, SPIN: ',ISPN

        PRINT '(A,I1,A)','SPIN ',ISPN,':'
        IF (DEBUG) write(6,*)'NEWWAVE CALLS OVERLAP MODE: 2'
        IF(MIXING1.NE.2) THEN
          CALL OVERLAP(2)
        ELSE
          OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
          READ(99)MPTS,MSPN
          READ(99)(HSTOR(IPTS,2),IPTS=1,MPTS)
          IF(ISPN.EQ.2)READ(99)(HSTOR(IPTS,2),IPTS=1,MPTS)
          CLOSE(99)
        ENDIF
        IF (DEBUG) write(6,*)'NEWWAVE CALLS OVERLAP MODE: 1'
        CALL OVERLAP(1)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS
C
C        write(6,*)'NEWWAVE:Back from OVERLAP 1 & 2'
        KBAS=0
        IF(ISPN.EQ.1)NVIRTTOT=0
c#ifdef MPI
c        NLOOP=N_REP/NPROC+1
c        MRUN=NPROC
c        DO ILOOP=1,NLOOP
c         IF(ILOOP.EQ.NLOOP) MRUN=N_REP-(ILOOP-1)*MRUN 
C
c        DO IREPP=1,MRUN
c         IREP=IREPP+(ILOOP-1)*NPROC
c#else
        DO 130 IREP=1,N_REP
c#endif
         initial=kbas+1
         write(6,*)'IREP=',IREP,'ISPN=',ISPN
         N_VIRT(IREP,ISPN)=0
         NBAS=NS_TOT(IREP)
         IF (NBAS.GT.NDH) THEN
          write(6,*)'NEWWAVE: NDH MUST BE AT LEAST: ',NBAS
          CALL STOPIT
         END IF
#ifndef SCALAPACK
         ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Ham'
         ENDIF
#else
         IF(NBAS.LT.200)THEN
          ALLOCATE(AHAM(NBAS,NBAS),STAT=IERR)
          IF(IERR.NE.0)THEN
            WRITE(6,*)'newwave_serial:Error allocating Ham'
          ENDIF
         END IF
#endif
!Uncommented aover alloc. for scissor operation
         ALLOCATE(AOVER(NBAS,NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Overlap'
         ENDIF
         ALLOCATE(AEVAL(NBAS),STAT=IERR)
         IF(IERR.NE.0)THEN
           WRITE(6,*)'newwave_serial:Error allocating Eval'
         ENDIF
C         ALLOCATE(ASC1(NBAS),STAT=IERR)
C         IF(IERR.NE.0)THEN
C           WRITE(6,*)'newwave_serial:Error allocating Sc1'
C         ENDIF
C         IF (NITER.GT.1)CALL SICSET(3,IREP)
         if(.true.)then !symm code
           DO 80 IBAS=1,NBAS
            AEVAL(IBAS)=0.0
            DO 70 JBAS=IBAS,NBAS
             KBAS=KBAS+1
!C           AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
!C           AHAM (JBAS,IBAS)=HSTOR(KBAS,2)
!YY quick fix
             AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
             if(sichammix .and. (MIXING1.EQ.2)) then
               AHAM (JBAS,IBAS)=HSTOR(KBAS,2)
             else
               AHAM (JBAS,IBAS)=HSTOR(KBAS,2)-SICSCI(KBAS,ISPN)
             endif
   70       CONTINUE
   80      CONTINUE
CJK11/2000
c          CALL SCISSOR(IREP,2,1)
           CALL SCISSOR(IREP)
         else
           AEVAL(:)=0.0
         end if
CJK11/2000

C
C GET EIGENVECTORS AND EIGENVALUES
C IF WE START A NEW GEOMETRY FROM AN OLD HAMILTONIAN, USE SINGULAR
C VALUE DECOMPOSITION TO AVOID PROBLEMS WITH SCREWED UP EIGENVALUES
C USE N_VIRT TO TEMPORARILY STORE THE NUMBER OF AVAILABLE EIGENSTATES 
C
         IF (NBAS.NE.0) THEN
c#ifdef MPI
c         ISVD=1
c         IF ((NITER.EQ.1).AND.(ISTSCF.EQ.1).AND.(IHIPOL.EQ.1))ISVD=0 
c          CALL PAMWAVE(1,ISVD,IREP,NBAS) 
c#else
         IF (NBAS.NE.0) THEN
          IF ((NITER.EQ.1).AND.(ISTSCF.EQ.1).AND.(IHIPOL.EQ.1)) THEN
           PRINT '(A)','NEWWAVE: USING SINGULAR VALUE DECOMPOSITION'
       
            call GTTIME(time_A)  ! Raja
           write(6,*)'newwave_serial:Function DIAGSVD not available' 
           CALL STOPIT

C           CALL DIAGSVD(NBAS,NBAS,NEIG,AHAM,AOVER,AEVAL,
C     &                  ASC1,SC2,DSINGV,1)
            call GTTIME(time_B)
             aux = time_B-time_A
            write(6,*) 'Time: Singular value decomposition diagonlz',aux
          ELSE
              call GTTIME(time_A)  ! Raja
C Perform diagonalization in parallel with ScaLAPACK with a call to senddata
           
!replacing #ifdef statements with allocation check of AHAM -cmd
!C           if(NBAS.GE.2000) then
!           if(NBAS.GE.200) then !NBAS limit set in AHAM allocation above -cmd
           if(.not.allocated(AHAM))then
             iirep=irep
             iispn=ispn
             inbas=nbas
             mvps=MAX_VIRT_PER_SYM
             iiev=1
             iinitial=initial
             iimesh=.FALSE.
C the call to scalapack is in senddata(111)
C             write(6,*)'newwave_serial calls Scalapack'
             IF(MPI_IO1)THEN
               CALL PREP_SCALA
             ENDIF
             CALL SENDDATA(111)
C             WRITE(6,*)'newwave_serial comes back from Scalapack'
           else
            if(.true.) then !for symmetry code
             CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,AUX,1)
            else
             WRITE(6,*)'CALLING DIAGGE3'
             CALL DIAGGE3(NBAS,initial,AHAM,AEVAL,1)
            end if
           endif
           NEIG=NBAS
              call GTTIME(time_B)
              aux = time_B-time_A
            write(6,*) 'Time required for DIAGGE',aux
          END IF
c#endif
         ENDIF
C
c#ifdef MPI
c         ENDDO          ! IREPP
C
c         DO IREPP=1,MRUN
c         IF (NCALLED.NE.0) THEN
c          CALL PAMWAVE(2,ISVD,IREP,NBAS)
c          NEIG=NBAS
c#endif         
         WRITE(97,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
         WRITE(98,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &              ' NUMBER OF BASES: ',NBAS
          N_VIRT(IREP,ISPN)=NEIG
*********************************TB*******************************
c          IF(ISPN.EQ.1.AND.IREP.EQ.1) NVIRTTOT=0
*********************************TB*******************************
c         DO JREP=1,IREP-1
c          NVIRTTOT=NVIRTTOT+NS_TOT(JREP)
c         ENDDO
          DO IEIG=1,NEIG
           EVALSAV(NVIRTTOT+IEIG)=AEVAL(IEIG) 
C          NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)
           NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)/LDMREP(IREP)!MRP 19OCT
           INDREP(NVIRTTOT+IEIG)=IREP
          END DO
           NVIRTTOT=NVIRTTOT+NEIG
          WRITE(97,*)NDMREP(IREP),NEIG
          WRITE(97,1300)(AEVAL(IEIG),IEIG=1,NEIG)
          WRITE(98,*)NDMREP(IREP),NEIG
          WRITE(98,1300)(AEVAL(IEIG),IEIG=1,NEIG)
 1300     FORMAT(' ',5G15.7)
         END IF
C
C WRITE SOME EIGENVECTORS
C
C        OPEN(69,FILE='EIGVECT')
C        DO I=1,100
C          DO J=1,100
C            WRITE(69,*)AHAM(I,J)
C          ENDDO
C        ENDDO
C        CLOSE(69)
C        STOP
C
C STORE ALL EIGENVECTORS THAT FIT INTO PSI_COEF
C
         NSAV(IREP,ISPN)=MIN(NEIG,MAX_VIRT_PER_SYM)
         IF(allocated(AHAM))THEN
          DO 120 ISAV=1,NSAV(IREP,ISPN)
           DO 110 IB=1,NBAS
            PSI_COEF(IB,ISAV,IREP,ISPN)=AHAM(IB,ISAV)
  110      CONTINUE
  120     CONTINUE
         END IF
c#ifdef MPI
c         END DO   !IREPP
c         END DO   !ILOOP
c#else
       IF(allocated(AHAM))THEN
        DEALLOCATE(AHAM,STAT=IERR)
        IF(IERR.NE.0)THEN
          WRITE(6,*)'newwave_serial:Error deallocating Ham'
        ENDIF
       END IF
       DEALLOCATE(AOVER,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'newwave_serial:Error deallocating Overlap'
       ENDIF
c       DEALLOCATE(ASC1,STAT=IERR)
c       IF(IERR.NE.0)THEN
c         WRITE(6,*)'newwave_serial:Error deallocating Sc1'
c       ENDIF
       DEALLOCATE(AEVAL,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'newwave_serial:Error deallocating Eval'
       ENDIF
  130   CONTINUE
c#endif

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

        IF(EF_MODE) write(6,*) 'OPTIMIZED MOMENT  MODE :'
        WRITE(6,*)'NVIRTTOT= ',NVIRTTOT
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
         IF(EF_MODE)
     &     CALL FERMILV(NVIRTTOT,ELEC,EF,TEMP,EVALSAV,OCCTMP,NDEG)
*********************************TB*******************************
C        DO ICOUNT=1,MAX_TOT*MXSPN
        DO ICOUNT=1,NHTOT
         NTEMP(ICOUNT)=0
        ENDDO

        ICOUNT=0
        
        DO  ISPN=1,NSPN

        IF (EF_MODE) EFERMI(ISPN)=EF
         DO 150 IREP=1,N_REP
          N_OCC(IREP,ISPN)=0
          DO IVIRT=1,N_VIRT(IREP,ISPN)
           ICOUNT=ICOUNT+1
            NTEMP(ICOUNT)=ISPN
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
           END IF
  140       CONTINUE
          END DO
c  140     CONTINUE
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  150    CONTINUE
         IF (ISAV .NE. 0) THEN
          write(6,*)'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',ISAV
          CALL STOPIT
         END IF
         ENDDO

        ELSE
C
C FROZEN OCCUPATION MODE. N_OCC IS ALREADY DEFINED
C
        IOFS=0
        NOCCU=0
         DO IVIRT=1,NVIRTTOT
          OCCTMP(IVIRT)=0.0D0
         END DO
        DO  ISPN=1,NSPN
         EFERMI(ISPN)= -1.0D20
         DO 170 IREP=1,N_REP
          IF (N_OCC(IREP,ISPN) .GT. NSAV(IREP,ISPN)) THEN
           ISAV= IREP
          ELSE
           DO 160 IVIRT=1,N_OCC(IREP,ISPN)
            NOCCU=NOCCU+1
            OCCTMP(IOFS+IVIRT)=OCCUPANCY(NOCCU)
            IF (OCCUPANCY(NOCCU) .GT. CUTOCC) THEN
C             EFERMI(ISPN)=MAX(EFERMI(ISPN),EVALSAV(IOFS+IVIRT))
             EFERMI(ISPN)=MAX(EFERMI(ISPN),EVALSAV(NOCCU))
            END IF
  160      CONTINUE
          END IF
          DO  JVIRT=1,N_VIRT(IREP,ISPN)
           NTEMP(IOFS+JVIRT)=ISPN
          END DO
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  170    CONTINUE
         IF (ISAV .NE. 0) THEN
          write(6,*)'NEWWAVE: NOT ENOUGH STATES FOR GIVEN OCCUPATION'
          write(6,*)'         IN REPRESENTATION: ',ISAV
          CALL STOPIT
         END IF
         WRITE(97,*)'FERMI LEVEL: ',EFERMI(ISPN),' TEMP: ',TEMP
         ENDDO
        END IF
C
C GET TRACE AND EVLOCC. FILL N_VIRT WITH THE NUMBER OF
C WAVEFUNCTIONS THAT IS ACTUALLY STORED IN PSI_COEF (NSAV -> N_VIRT)
C
        WRITE(6,*)'NVIRTTOT before trace = ',NVIRTTOT
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
         write(6,*)'NEWWAVE: MAX_OCC MUST BE AT LEAST ',ISAV
         CALL STOPIT
        END IF
       ENDDO
        IF(NSPN.EQ.1.AND.EF_MODE) THEN 
           E_DN=E_UP
        ENDIF
       Write(6,*)'ELECTRON :',E_UP,E_DN
        WRITE(6,*)'NVIRTTOT after trace = ',NVIRTTOT
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
C
        
        WRITE(6,*)'ELECTRONS OF SPIN UP : ', E_UP
        WRITE(6,*)'ELECTRONS OF SPIN DN : ', E_DN
        WRITE(97,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(97,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        WRITE(98,*)'FERMI LEVEL: ',EF,' TEMP: ',TEMP
        WRITE(98,*)'SUMMARY OF EVALUES AND THEIR OCCUPANCIES:'
        ECUT=MAX(0.0D0,EF+2.0D1*TEMP)
        WRITE(6,*)'NVIRTTOT before EVALXXX = ',NVIRTTOT
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
C deallocate local arrays that were created for diagonalization

       DEALLOCATE(NDEG,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating NDEG'
       ENDIF
       DEALLOCATE(INDREP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating INDREP'
       ENDIF
       DEALLOCATE(NTEMP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating NTEMP'
       ENDIF
       DEALLOCATE(OCCTMP,STAT=IERR)
       IF(IERR.NE.0)THEN
         WRITE(6,*)'NEWWAVE:Error Deallocating OCCTMP'
       ENDIF
       RETURN
       END
