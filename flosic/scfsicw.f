C UTEP Electronic Structure Lab (2020)
C> @detail Diagonalization routine for SIC
C> @param[in] SIC SIC 2D potential array
C
C **************************************************************
C
C> @author NEWWAVE VERSION DIRK POREZAG AUGUST 1994
C> @author PARALLEL VERSION JENS KORTUS SEPTEMBER 1999
C> @brief NEWWAVE DIAGONALIZES THE HAMILTONIAN AND CALCULATES THE
C>        OCCUPATION NUMBERS
C
       SUBROUTINE SCFSICW(NITER,TRACE)
C
C> @author WRITTEN BY MARK R PEDERSON (1986-1989)
c
       use global_inputs,only : inbas,iiev,iimesh,iinitial,mpi_io1,
     &        fixm1,mixing1,uniham1
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
!SIC module
       use LOCORB,only : TMAT,MORB,ZSIC,IRBSIC
       use MOCORB,only : SLAT,NFRM,ZTZL,JJJJJJ
       use SICMAT,only : SIC
       use diagv1,only : NORB,PHIRES
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:02 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NITER, MAX_TOT, I, IB, IBAS, ICOUNT, IDEG, IEIG, IL,
     & INDREP, INEW, IOCC, IOFS, IOLD, IRB, IREP, ISAV, ISORT, ISPFAC,
     & ITER, IVIRT, IWF, J, JBAS, JNEW, JOLD, JREP, JSORT, JVIRT, JWF,
     & KBAS, L, N_VIRT, NBAS, NDEG, NEIG, NOCCU, NSAV, NTEMP, NTMP,
     & NTMP2, NVIRT_DN, NVIRT_UP, NVIRTTOT
       INTEGER :: MPTS,MSPN
       REAL*8 :: SYMBOL , TRACE, CUTOCC, DIAG, DSINGV, ECUT, EF, ELEC,
     & ERR1, ERR2, ERROR, EVALSAV, OCCTMP, TEMP, TIME1,
     & TIME2
       REAL(8) :: HAM1,PSII
!       INCLUDE 'commons.inc'
! YY moved from commons.inc
!       COMMON/FOR_DIAG/OVER(NDH,NDH),HAM(NDH,NDH),FILO(NDH,NDH),
!     &  EVAL(NDH),SC1(NDH),SC2(NDH)
!       COMMON/LOCORB/TMAT(MAX_OCC,MAX_OCC,MXSPN),MORB(2),ZSIC,IRBSIC
!       COMMON/MOCORB/SLAT(MAX_OCC,MAX_OCC,MXSPN),NFRM(2),ZTZL,JJJJJJ
!       SAVE
       PARAMETER (MAX_TOT=NDH*MAX_REP)
       LOGICAL EXIST,FERMISTAT
       LOGICAL AVERAGE,EF_MODE,HAMAVG,RENORM,printmore
       CHARACTER*4 FLINE
       CHARACTER*12 EVALSTR
       CHARACTER*7 NAMES
       DIMENSION NAMES(3)
       DIMENSION EVALSAV(MAX_TOT*MXSPN),OCCTMP(MAX_TOT*MXSPN)
       DIMENSION NDEG(MAX_TOT*MXSPN),INDREP(MAX_TOT*MXSPN),
     &  NSAV(MAX_REP,MXSPN)
       DIMENSION N_VIRT(MAX_REP,MXSPN)
       !DIMENSION DIAG(NDH,MAX_REP) !<<<< YY Not used
       DIMENSION NTEMP(MAX_TOT*MXSPN)

!       DIMENSION PHIRES(NDH,NDH)
!       DIMENSION FOOVR2(NDH,NDH)
!       real*8, allocatable :: PHIRES(:,:)
       real*8, allocatable :: FOOVR2(:,:)
       real*8, allocatable :: M(:,:)
       real*8, allocatable :: UNIHAM(:,:)

!       COMMON/MIXPOT1/POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN)
!       COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)

C
C DEFINE TEMPERATURE, MINIMUM OCCUPANCY AND SMALLEST ALLOWED
C EIGENVALUE OF OVERLAP MATRIX FOR SINGULAR VALUE DECOMPOSITION
C
       DATA TEMP  /1.0D-4/
       DATA CUTOCC/1.0D-10/
       DATA DSINGV/2.0D-4/
       DATA NAMES/'BROYDEN','KBROY1','KBROY2'/
C
C CHECKING AND SETTING UP SOME STUFF
C
       printmore=.false. !diagnostic prints and error checks

       IF (N_REP.GT.MAX_REP) THEN
        PRINT *,'NEWWAVE: MAX_REP MUST BE AT LEAST: ',N_REP
        CALL STOPIT
       END IF
       TRACE=0.0D0
**************************************************************
*      Set default mode for calculation of E_F
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
       CALL CHECK_INPUTS
       IF(FIXM1) EF_MODE=.FALSE.
       IF ((FLINE.EQ.'FIXM' ).OR.(FLINE.EQ.'fixm'))
     &   EF_MODE=.FALSE.
       IF ((FLINE.EQ.'OCCU').OR.(FLINE.EQ.'occu')) THEN
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
!YY Hack fix for the frozen occupation mode 1/24/18
           if (N_OCC(IREP,ISPN) .LT. NS_TOT(IREP))  then
             do L=IOCC+1,IOCC+NS_TOT(IREP)-N_OCC(IREP,ISPN)
              OCCUPANCY(L)=0.0d0
             end do
             IOCC=IOCC+NS_TOT(IREP)-N_OCC(IREP,ISPN)
             N_OCC(IREP,ISPN) = NS_TOT(IREP)
           end if
!end Hack fix
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


!#############################
!#    HAM MIXING SECTION     #
!#############################
       CALL CHECK_INPUTS
       print *,"MIXING",MIXING1," 1:POTENTIAL,2:HAMILT,3:DENMAT"
       IF(MIXING1.EQ.2)THEN !copied from DFT newwave
         DO ISPN=1,NSPN
           !YY. Calculate hamiltonian. This will create HAMOLD.
           CALL OVERLAP(2)
         END DO
         OPEN(99,FILE='HAMOLD', FORM='UNFORMATTED',STATUS='UNKNOWN')
         READ(99)MPTS,MSPN
         CLOSE(99)
         CALL HAMMIXDRV(MPTS)
       ENDIF

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
!<<<<<<<<<<<<<
c  Force 0.0
c  open file for storing lagrange multiplier matrix
c
       open(unit=37,file='LAGMULT.DAT',form='formatted')
!>>>>>>>>>>>>>
       ISPFAC=2/NSPN
       DO 240 ISPN=1,NSPN
        NWFS(ISPN)=0
c       IF (ISPN.EQ.1) ELEC=E_UP
c       IF (ISPN.EQ.2) ELEC=E_DN
        IF(FLINE.EQ.'FIXM' .OR. FLINE.EQ.'fixm') THEN
         WRITE(97,1100)FLINE
         WRITE(6,1100)FLINE
        ENDIF
        WRITE(97,*)'********* NEW TRY ************, SPIN: ',ISPN
        WRITE(98,*)'********* NEW TRY ************, SPIN: ',ISPN

        PRINT '(A,I1,A)','SPIN ',ISPN,':'
        IF (DEBUG) PRINT *,'NEWWAVE CALLS OVERLAP MODE: 2'
        CALL CHECK_INPUTS
!        if(MIXING1.EQ.0) CALL OVERLAP(2)
        IF(MIXING1.EQ.0) THEN
          CALL OVERLAP(2)
        ELSE !copied from DFT newwave
          OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
          READ(99)MPTS,MSPN
          READ(99)(HSTOR(KBAS,2),KBAS=1,MPTS)
          IF(ISPN.EQ.2)READ(99)(HSTOR(KBAS,2),KBAS=1,MPTS)
          CLOSE(99)
        ENDIF
        IF (DEBUG) PRINT *,'NEWWAVE CALLS OVERLAP MODE: 1'
        CALL OVERLAP(1)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS
C
        KBAS=0
c       NVIRTTOT=0
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

           OVER(IBAS,JBAS)=OVER(JBAS,IBAS)
           HAM(IBAS,JBAS)=HAM(JBAS,IBAS)

   70     CONTINUE
   80    CONTINUE

C
C GET EIGENVECTORS AND EIGENVALUES
C IF WE START A NEW GEOMETRY FROM AN OLD HAMILTONIAN, USE SINGULAR
C VALUE DECOMPOSITION TO AVOID PROBLEMS WITH SCREWED UP EIGENVALUES
C USE N_VIRT TO TEMPORARILY STORE THE NUMBER OF AVAILABLE EIGENSTATES
C
         IF (NBAS.NE.0) THEN
          IF (NBAS.NE.0) THEN
           IF ((NITER.EQ.1).AND.(ISTSCF.EQ.1).AND.(IHIPOL.EQ.1)) THEN
            PRINT '(A)','NEWWAVE: USING SINGULAR VALUE DECOMPOSITION'
            CALL DIAGSVD(NDH,NBAS,NEIG,HAM,OVER,EVAL,
     &                   SC1,SC2,DSINGV,1)
           ELSE IF(NITER.EQ.1) THEN
            OPEN(49,FILE='STORE',FORM='UNFORMATTED')
            WRITE(49)OVER
!           WRITE(6,*)(HAM(I,I),I=1,NBAS)
!           WRITE(6,*)(OVER(I,I),I=1,NBAS)
            PRINT*,'CALLING DIAGGE SIC'

            CALL FLUSH(6)

!           print *,"BEFORE DIAG",NDH,NBAS !YY
!           do I=1,NBAS
!            print *,'H',HAM(1:NBAS,I) !YY
!           end do
!           do I=1,NBAS
!            print *,'O',OVER(1:NBAS,I)
!           end do

            !CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
            CALL DIAGGE_FO(NDH,NBAS,HAM,OVER,EVAL,SC1,1)

!            print *,"AFTER DIAG"
!            do I=1,NBAS
!            print *,'H',HAM(1:NBAS,I) !YY
!!            print *,'O',OVER(1:NBAS,I)
!            end do

            NEIG=NBAS
            REWIND(49)
            READ(49)OVER
            CLOSE(49,STATUS='DELETE')
!            do I=1,NBAS !YY
!             print *,'o',OVER(1:NBAS,I)
!            end do
           ELSE  !NITER.NE.1

            OPEN(49,FILE='STORE',FORM='UNFORMATTED')
            WRITE(49)OVER
C           WRITE(6,*)(HAM(I,I),I=1,NBAS)
C           WRITE(6,*)(OVER(I,I),I=1,NBAS)
            PRINT*,'CALLING DIAGV'
C           CALL FLUSH(6)

C NEW WORKING AREA...
            DO J=1,NBAS  !
             DO IRB=1,NFRM(ISPN)   ! THIS IS A FERMI ORBITAL
              SC1(IRB)=0.0D0
              DO IWF=1,NFRM(ISPN)    !THIS IS A SUM OVER OCCUPIED WAVEFUNCTION
               SC1(IRB)=SC1(IRB)
     &                 +TMAT(IRB,IWF,ISPN)*PSI_COEF(J,IWF,IREP,ISPN)
              END DO
             END DO
             DO IRB=1,NFRM(ISPN)
              PSI_COEF(J,IRB,IREP,ISPN)=SC1(IRB)
             END DO
            END DO
            !IF(NITER.NE.1)THEN  NITER is.NE.1 already
             DO IWF=1,NBAS
              DO JWF=1,NBAS
               OVER(JWF,IWF)=HAM(JWF,IWF)
              END DO
             END DO
!
! CMD translating O(N^4) loops to 2 O(N^3) mat mults
!
! HAM = PSI(transpose) * OVER * PSI 
!  M  = PSI(transpose) * OVER
             allocate(M(NBAS,NBAS))
             call dgemm('T','N',NBAS,NBAS,NBAS,1.0d0,
     &         PSI_COEF(1,1,IREP,ISPN),NDH,
     &         OVER,NDH, 0.0d0,
     &         M, NBAS)
!  HAM= M * PSI
             call dgemm('N','N',NBAS,NBAS,NBAS,1.0d0,
     &         M,  NBAS,
     &         PSI_COEF(1,1,IREP,ISPN),NDH, 0.0d0,
     &         HAM,NDH)
             deallocate(M)

! ORIGINAL CODE
!             DO IWF=1,NBAS!NWFS(ISPN)
!              DO JWF=1,NBAS!NWFS(ISPN)
!               HAM(JWF,IWF)=0.0D0
!               DO I=1,NBAS
!                PSII=PSI_COEF(I,IWF,IREP,ISPN)
!                DO J=1,NBAS
!                 HAM(JWF,IWF)=HAM(JWF,IWF)+OVER(J,I)*
!     &            PSI_COEF(J,JWF,IREP,ISPN)*PSI_COEF(I,IWF,IREP,ISPN)
!                END DO
!               END DO
!              END DO
!!              PRINT 2016,(HAM(JWF,IWF),JWF=1,NBAS)
!             END DO

 2016         FORMAT(20G11.3)
!YY
!            print*,'HAM before'
!            DO IWF=1,NBAS
!             PRINT 2016,(HAM(JWF,IWF),JWF=1,NBAS)
!            END DO
!            PRINT*,'SIC MATRIX:'
!            DO IWF=1,NBAS
!             PRINT 2016,(SIC(JWF,IWF,ISPN),JWF=1,NBAS)
!            END DO
             PRINT*,'DFA-SIC MATRIX:'
             DO IWF=1,NBAS
              DO JWF=1,NBAS
               HAM(JWF,IWF)=HAM(JWF,IWF)+SIC(IWF,JWF,ISPN)
              END DO
!             PRINT 2016,(HAM(JWF,IWF),JWF=1,NBAS)
             END DO
!<<<<<<<<<
c  Force 0.0
c  store lagrange multiplier matrix to disc
c
             do iwf = 1,nfrm(ispn)
               write(37,*) (HAM(JWF,IWF),JWF=1,NFRM(ISPN))
             end do
!>>>>>>>>
             err1=0.0d0
             err2=0.0d0
             do iwf=1,nfrm(ispn)
              do jwf=nfrm(ispn)+1,nbas
               err1=err1+abs(ham(jwf,iwf))
               err2=err2+abs(ham(iwf,jwf))
              end do
             end do
             print*,'diagv hamiltonian error:',err1,err2
!YY At this point Hamiltonian mixing should be done.

!FOR Unified Hamiltonian SCF-SIC    Kamal Sharkas CMU Jun 2019 BEG
            IF (UNIHAM1) THEN
             print *,"Unified Hamiltonian is used in SCFSICW"
             allocate(UNIHAM(NBAS,NBAS))
             do iwf = 1,nbas
              do jwf = 1,nbas
               UNIHAM(jwf,iwf)=0.0d0
              end do
             end do
c  Occupied-occupied block
c
             do IWF = 1,nfrm(ispn)
              UNIHAM(iwf,iwf) = HAM(iwf,iwf)
             end do
c
c occupied-virtual blocks
             do iwf = nfrm(ispn)+1,nbas
              do jwf = 1,nfrm(ispn)
               UNIHAM(jwf,iwf) = HAM(jwf,iwf)
               UNIHAM(iwf,jwf) = UNIHAM(jwf,iwf)
              end do
             end do
c  virtual-virtual block
             do iwf = nfrm(ispn)+1,nbas
              do jwf = nfrm(ispn)+1,nbas
               UNIHAM(iwf,jwf) = HAM(iwf,jwf)
              end do
             end do
c
c  print unified hamiltonian
c
C            do iwf = 1,nbas
C             PRINT 2016,(UNIHAM(JWF,IWF),JWF=1,NBAS)
C            end do
            ENDIF
!FOR Unified Hamiltonian SCF-SIC    Kamal Sharkas CMU Jun 2019 END
!YY allocate PHIRES
!           allocate(PHIRES(NDH,NDH))
            allocate(PHIRES(NBAS,NBAS))
            PRINT*,'NFRM,NDH',NFRM(ISPN),NBAS
! LB HAM,PHIRES shared over modules now
!           call diagv_fo(NDH,NFRM(ISPN),NBAS,HAM,PHIRES,0)
!FOR Unified Hamiltonian SCF-SIC    Kamal Sharkas CMU Jun 2019 BEG
            IF (UNIHAM1) THEN
             call diagv_fo_uni(NBAS,NFRM(ISPN),NBAS,UNIHAM,PHIRES,0)

!            PRINT*,'NEW ORBITALS:'
!            DO INEW=1,NBAS
!             PRINT 2016,(PHIRES(INEW,JOLD),JOLD=1,NBAS)
!            END DO
!            PRINT*,'"Diagonalized Matrix:"'

!            DO INEW=1,NBAS
!             DO JNEW=1,NBAS
!              OVER(JNEW,INEW)=0.0D0
!              DO IOLD=1,NBAS
!               DO JOLD=1,NBAS
!                OVER(JNEW,INEW)=OVER(JNEW,INEW)+
!     &           PHIRES(INEW,IOLD)*PHIRES(JNEW,JOLD)*UNIHAM(JOLD,IOLD)
!!               OVER(JNEW,INEW)=OVER(JNEW,INEW)+
!!    &           PHIRES(INEW,IOLD)*PHIRES(JNEW,JOLD)*HAM(JOLD,IOLD)
!               END DO
!              END DO
!             END DO
!             EVAL(INEW)=OVER(INEW,INEW)

!             IF(INEW.LE.NFRM(ISPN))THEN
!              PRINT 2016,(OVER(JNEW,INEW),JNEW=1,NFRM(ISPN))
!             END IF

!            END DO
!            NEIG=NBAS

             DO INEW=1,NBAS
              OVER(INEW,INEW)=0.0D0
              DO IOLD=1,NBAS
               DO JOLD=1,NBAS
                OVER(INEW,INEW)=OVER(INEW,INEW)+
     &           PHIRES(INEW,IOLD)*PHIRES(INEW,JOLD)*UNIHAM(JOLD,IOLD)
!               OVER(JNEW,INEW)=OVER(JNEW,INEW)+
!    &           PHIRES(INEW,IOLD)*PHIRES(JNEW,JOLD)*HAM(JOLD,IOLD)
               END DO
              END DO
              EVAL(INEW)=OVER(INEW,INEW)
!             IF(INEW.LE.NFRM(ISPN))THEN
!              PRINT 2016,(OVER(JNEW,INEW),JNEW=1,NFRM(ISPN))
!             END IF
             END DO
!            NEIG=NBAS

C ROTATE BACK INTO BASIS SET REPRESENTATION





!            print*,'Eigenvalues: ',EVAL(1:NBAS)

             NEIG=NBAS
!            PRINT*, 'BACKROTATION: O(N^3)', NBAS
C ROTATE BACK INTO BASIS SET REPRESENTATION
! Not quite yet merging the TWOSTEP !YY. 8.13.19
!            IF(TWOSTEP1) THEN
!             DO INEW=1,NBAS
!              DO IBAS=1,NBAS
!               HAM(IBAS,INEW)=0.0D0
!               DO IOLD=1,NBAS
!                HAM(IBAS,INEW)=HAM(IBAS,INEW)
!    &           +PHIRES(INEW,IOLD)*PSI_COEF_FLO(IBAS,IOLD,IREP,ISPN)
!               END DO
!              END DO
!             END DO
!            ELSE
              DO INEW=1,NBAS
               DO IBAS=1,NBAS
                HAM(IBAS,INEW)=0.0D0
                DO IOLD=1,NBAS
                 HAM(IBAS,INEW)=HAM(IBAS,INEW)
     &           +PHIRES(INEW,IOLD)*PSI_COEF(IBAS,IOLD,IREP,ISPN)
                END DO
               END DO
              END DO
!            ENDIF
! UNIHAM no longer needed
             deallocate(UNIHAM)
            ELSE  
!FOR Unified Hamiltonian and TWO-STEP SCF-SIC    Kamal Sharkas CMU Jun 2019 ELSE             

! Jacobi sweep
             INBAS=NBAS
             NORB=NFRM(ISPN)
             call diagv_fo(0)
!            CALL TRACER('STOPPING FOR NOW')
!            CALL STOPIT
! LB

             PRINT*,'NEW ORBITALS:'
!            DO INEW=1,NBAS
!             PRINT 2016,(PHIRES(INEW,JOLD),JOLD=1,NBAS)
!            END DO

             PRINT*,'"Diagonalized Matrix:"'
!            if (printmore) then
!             DO INEW=1,NBAS
!              DO JNEW=1,NBAS
!               OVER(JNEW,INEW)=0.0D0
!               DO IOLD=1,NBAS
!                DO JOLD=1,NBAS
!                 OVER(JNEW,INEW)=OVER(JNEW,INEW)+
!    &           PHIRES(INEW,IOLD)*PHIRES(JNEW,JOLD)*HAM(JOLD,IOLD)
!                END DO
!               END DO
!              END DO
!              EVAL(INEW)=OVER(INEW,INEW)
!              IF(INEW.LE.NFRM(ISPN))THEN
!               PRINT 2016,(OVER(JNEW,INEW),JNEW=1,NFRM(ISPN))
!              END IF
!             END DO
!             NEIG=NBAS
C ROTATE BACK INTO BASIS SET REPRESENTATION
!             print*,'Eigenvalues: ',EVAL(1:NBAS)
!            else
              !only calculate diagonal elements.
              DO INEW=1,NBAS
!              DO JNEW=1,NBAS
               OVER(INEW,INEW)=0.0D0
               DO IOLD=1,NBAS
                DO JOLD=1,NBAS
                 OVER(INEW,INEW)=OVER(INEW,INEW)+
     &           PHIRES(INEW,IOLD)*PHIRES(INEW,JOLD)*HAM(JOLD,IOLD)
                END DO
               END DO
!              END DO
               EVAL(INEW)=OVER(INEW,INEW)
              END DO
             
!            end if

             NEIG=NBAS
             PRINT*, 'BACKROTATION: O(N^3)', NBAS
C ROTATE BACK INTO BASIS SET REPRESENTATION
             DO INEW=1,NBAS
              DO IBAS=1,NBAS
               HAM(IBAS,INEW)=0.0D0
               DO IOLD=1,NBAS
                HAM(IBAS,INEW)=HAM(IBAS,INEW)
     &          +PHIRES(INEW,IOLD)*PSI_COEF(IBAS,IOLD,IREP,ISPN)
               END DO
              END DO
             END DO
!!FOR Unified Hamiltonian Kamal Sharkas CMU Jun 2019 END
            ENDIF
!YY PHIRES no longer needed
             deallocate(PHIRES)
             DO INEW=1,NBAS
              DO IBAS=1,NBAS
               PSI_COEF(IBAS,INEW,IREP,ISPN)=HAM(IBAS,INEW)
              END DO
             END DO

!            print*,'New orbitals:'
!            do inew=1,nbas
!             print '(100F10.5)',(HAM(ibas,inew),ibas=1,nbas)
!            end do

             REWIND(49)
             READ(49)OVER

C            PRINT *, 'CALCULATE FOOVR2 O(N^4)', NBAS
!YY allocatable FOOVR2
             if(printmore)then
              allocate(FOOVR2(NDH,NDH))
              FOOVR2=0.0d0
              do inew=1,nbas
               do jnew=1,nbas
                do ibas=1,nbas
                 do jbas=1,nbas
                  FOOVR2(inew,jnew)=FOOVR2(inew,jnew)+
     &            PSI_COEF(ibas,inew,IREP,ISPN)
     &           *PSI_COEF(jbas,jnew,IREP,ISPN)*OVER(ibas,jbas)
                 end do
                end do
               end do
              end do
              print*,'Overlap of new orbitals:'
              error=0.0d0
              do inew=1,nbas
               do jnew=inew+1,nbas
                error=error+abs(foovr2(inew,jnew))
               end do
c              print '(100F10.5)',(FOOVR2(inew,jnew),jnew=1,nbas)
              end do
              deallocate(FOOVR2)
              print*,' diav overlap error:',error
             end if !printmore
           
            !END IF
            PRINT*,'IS IT OK?',NWFS(ISPN)

C        IF(DUMMY.LT.0)RETURN
C END OF NEW WORKING AREA...
            REWIND(49)
            READ(49)OVER
            CLOSE(49,STATUS='DELETE')
           END IF
          ENDIF
C
          WRITE(97,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &               ' NUMBER OF BASES: ',NBAS
          WRITE(98,*)'REP: ',IREP,' DIM: ',NDMREP(IREP),
     &               ' NUMBER OF BASES: ',NBAS
          N_VIRT(IREP,ISPN)=NEIG
*********************************TB*******************************
          IF(ISPN.EQ.1.AND.IREP.EQ.1) NVIRTTOT=0
*********************************TB*******************************
c         DO JREP=1,IREP-1
c          NVIRTTOT=NVIRTTOT+NS_TOT(JREP)
c         ENDDO
          DO IEIG=1,NEIG
           EVALSAV(NVIRTTOT+IEIG)=EVAL(IEIG)
C          NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)
           NDEG(NVIRTTOT+IEIG)=NDMREP(IREP)/LDMREP(IREP)!MRP 19OCT
           INDREP(NVIRTTOT+IEIG)=IREP
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
!<<<<<<<<<<
c  Force 0.0
c  close file for lagrange multiplier matrix
c
       close(37)
!>>>>>>>>>>

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
          CALL FERMILVOCC(NVIRT_UP,ELEC,EF,TEMP,POTIN,POTOUT,NTEMP)
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
          CALL FERMILVOCC(NVIRT_DN,ELEC,EF,TEMP,POTIN,POTOUT,NTEMP)
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
     &     CALL FERMILVOCC(NVIRTTOT,ELEC,EF,TEMP,EVALSAV,OCCTMP,NDEG)
*********************************TB*******************************
        ICOUNT=0

        DO  ISPN=1,NSPN

         IF (EF_MODE) EFERMI(ISPN)=EF
         DO 150 IREP=1,N_REP
          N_OCC(IREP,ISPN)=0
          DO IVIRT=1,N_VIRT(IREP,ISPN)
           ICOUNT=ICOUNT+1
C          IF (OCCTMP(IOFS+IVIRT) .LT. CUTOCC) GOTO 140
C          IF (OCCTMP(ICOUNT) .LT. CUTOCC) GOTO 140

           IF (IVIRT .GT. NSAV(IREP,ISPN)) THEN
            ISAV= MAX(ISAV,IVIRT)
           ELSE
            NOCCU=NOCCU+1
C           OCCUPANCY(NOCCU)=OCCTMP(IOFS+IVIRT)
C           OCCUPANCY(NOCCU)=OCCTMP(IOFS+IVIRT)/LDMREP(IREP) ! MRP 19Oct98
            OCCUPANCY(NOCCU)=OCCTMP(ICOUNT)/LDMREP(IREP) ! MRP 19Oct98
            N_OCC(IREP,ISPN)=N_OCC(IREP,ISPN)+1
            NTEMP(ICOUNT)=ISPN
           END IF
  140      CONTINUE
          END DO
c  140    CONTINUE
          IOFS=IOFS+N_VIRT(IREP,ISPN)
  150    CONTINUE
         IF (ISAV .NE. 0) THEN
          PRINT *,'NEWWAVE: MAX_VIRT_PER_SYM MUST BE AT LEAST: ',ISAV
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
C
        DO  ISPN=1,NSPN
         EFERMI(ISPN)= -1.0D20
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
!The next three lines will keep track of spin up (1)
!and spin down (2) occupations.
          DO  JVIRT=1,N_VIRT(IREP,ISPN)
           NTEMP(IOFS+JVIRT)=ISPN
          END DO
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
c         JVIRT=JVIRT+1
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
  210   CONTINUE
  220  CONTINUE
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
        WRITE(97,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                 EVALSAV(IOCC),OCCTMP(IOCC)
        WRITE(98,1500) IOCC,INDREP(IOCC),NDEG(IOCC),NTEMP(IOCC),
     &                 EVALSAV(IOCC),OCCTMP(IOCC)
        IF (EVALSAV(IOCC).GT.ECUT) GOTO 290
  230  CONTINUE
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

