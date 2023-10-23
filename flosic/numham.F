C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
       SUBROUTINE NUMHAM
C
C WRITTEN BY MARK R PEDERSON (1988-1995)
C

       use hstor1,only : hstor
       use mixpot1,only : POTIN,POTDV=>POTOUT
       use mesh1,only : wmsh,nmsh
       use common2,only : ISPN, NSPN, IGGA !Need IGGA to check if XC is MGGA
       use common3,only : RMAT
       use common5,only : HAVEHAM, ISTSCF
       use common8,only : REP, N_REP, NS_TOT
       use mpidat1,only : NPROC
       use global_inputs, only : libxc1
       use xtmp2a,only : HSTORSIZE
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IMODE, IPTR, IPTS, IREC, KPTS, KREP, MPTS, MSPN,
     & NDUM, NREC
       REAL*8 :: DUMMY , SUM, TIME1, TIME2
       SAVE
       LOGICAL EXIST
! YY. Additional local variables
       real(8),allocatable :: HXC(:)
       integer :: HSTORSIZEF,NHSTORSIZE,NSPNF
       logical isMGGA
C       COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POTDV(MAX_PTS*MXSPN)
C
C IMODE=1 USE O(N)+PARALLEL APPROACH
C IMODE=2 USE STANDARD APPROACH - FASTER FOR SMALL SYSTEMS / FEW CPUS
C      
       IMODE=2
C
#ifdef MPI
C
       IF (NPROC.GE.3) IMODE=1
C
#endif
C
       INQUIRE(FILE='HAMMODE',EXIST=EXIST)
       IF (EXIST) THEN
        OPEN(99,FILE='HAMMODE',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        READ(99,*) IMODE
        REWIND(99)
        WRITE(99,*) IMODE,' 1=O(N)/PARALLEL OR 2=STANDARD'
        CLOSE(99)
       END IF
       IF (IMODE .EQ. 1) THEN
        PRINT '(A)','USING O(N)/PARALLEL HAMILTONIAN ROUTINE'
       ELSE
        PRINT '(A)','USING STANDARD HAMILTONIAN ROUTINE'
       END IF
C
C ZERO POTOLD. READ POTOLD IF POTENTIAL IS USED AS STARTING POINT
C
       CALL GTTIME(TIME1)
       DO 10 IPTS=1,NMSH*NSPN
        POTDV(IPTS)=0.0D0
   10  CONTINUE
       IF (ISTSCF.EQ.2) THEN
        INQUIRE(FILE='POTOLD',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
         write(6,*)'NUMHAM: POTOLD DOES NOT EXIST'
         CALL STOPIT
        END IF
        OPEN(99,FILE='POTOLD',FORM='UNFORMATTED',STATUS='OLD')
        REWIND(99)
        READ(99,END=20)MPTS,MSPN,NDUM
        IF (MPTS.NE.NMSH) THEN
         write(6,*)'NUMHAM: NMSH IN POTOLD DIFFERS FROM CURRENT JOB'
         write(6,*)'MPTS,NMSH=',MPTS,NMSH,MSPN,NDUM
         CALL STOPIT
        END IF
        KPTS=NMSH*MIN(MSPN,NSPN)
        READ(99,END=20)(POTDV(IPTS), IPTS=1,KPTS)
        GOTO 30
   20   write(6,*)'NUMHAM: POTOLD UNREADABLE'
        CALL STOPIT
   30   CLOSE(99)
C
C DEFINE SPIN DOWN POTENTIAL IF POTOLD POTENTIAL IS SPIN UNPOLARIZED
C
        IF ((MSPN.EQ.1).AND.(NSPN.EQ.2)) THEN
         DO 40 IPTS=1,NMSH
          POTDV(IPTS+NMSH)=POTDV(IPTS)
   40    CONTINUE
        END IF
C
C MULTIPLY POTENTIAL BY WMSH AND CALL CORRECT SUBROUTINE
C
        IPTR=(ISPN-1)*NMSH
        DO IPTS=1,NMSH
         IF(POTDV(IPTS).NE.POTDV(IPTS))WRITE(6,*)'POTDV NAN'
         IF(WMSH (IPTS).NE.WMSH (IPTS))WRITE(6,*)'WMSH  NAN'
         POTDV(IPTS)=WMSH(IPTS)*POTDV(IPTS+IPTR)
        END DO
        WRITE(6,*)'NMSH:',NMSH,' IN NUMHAM'
        CALL FLUSH(6) !YY copied from SIC code
        SUM=0.0D0
        DO I=1,NMSH
        SUM=SUM+POTDV(I)
        END DO
        IF (IMODE.EQ.1) THEN
C         write(6,*)'numham:calling overnum'
         CALL OVERNUM(1)
C         write(6,*)'numham:done overnum'
        ELSE
         CALL PATCH(POTDV)
        END IF
       END IF

! YY. Adding Hamiltonian elements from Libxc
! (Turn this on when you are using libxc or MGGA)
       if(IGGA(1)==3 .or. IGGA(2)==3) isMGGA = .true.

       call CHECK_INPUTS
       !write(*,*) "LIBXC1", libxc1, isMGGA,IGGA(1),IGGA(2)

       if(libxc1 .or. isMGGA) then
         inquire(file='HAMXC',EXIST=EXIST)
         if (EXIST) then
         !if(.TRUE.) then
! YY. read and add Hamiltonian mixing elements from patchmod to hstor
           open(99999,file='HAMXC',FORM='UNFORMATTED')
           rewind(99999)
           read(99999) NHSTORSIZE,NSPNF
!          CALL TRACER('HSTORSIZEF',NHSTORSIZE)
!          CALL TRACER('NSPNF',NSPNF)
!          read(99999,*) HSTORSIZEF,SPNF
           IF(HSTORSIZE/=NHSTORSIZE) THEN
!            CALL TRACER('HSTOR SIZES DIFFER!')
             !Skip adding extra Hamiltonian elements if HSTOR sizes are
             !different. Using goto statement rather than stopit since this
             !happens in restarting calculation using RUNS file.
             write(*,*) "NOT adding extra XC Ham. elements to HSTOR"
             close(99999) !Close the file before proceeding
             goto 111
             CALL STOPIT
           ENDIF
           IF(NSPNF/=NSPN) THEN
!            CALL TRACER('SPINS DIFFER!')
             print *,"SPINS DIFFER!"
             CALL STOPIT
           ENDIF
           allocate(HXC(HSTORSIZE))
           read(99999) (HXC(IPTS),IPTS=1,HSTORSIZE)
!          CALL TRACER('READ FIRST PART')
           if(ISPN == 1) then
           !if(ISPN == 1 .AND. N_OCC(1,1) > 0) then
           ! check N_OCC is non zero. -> RHOV array need to be zeroed in getvlxc
           ! instead.
             do IPTS=1,HSTORSIZE
!              READ(99999,*)HXC(IPTS)
               HSTOR(IPTS,2)=HSTOR(IPTS,2) + HXC(IPTS)
             enddo
           end if
           if(ISPN == 2) then
           !if(ISPN == 2 .AND. N_OCC(1,2) > 0) then
             read(99999) (HXC(IPTS),IPTS=1,HSTORSIZE)
!            CALL TRACER('READ SECOND PART')
             do IPTS=1,HSTORSIZE
!              READ(99999,*)HXC(IPTS)
               HSTOR(IPTS,2)=HSTOR(IPTS,2) + HXC(IPTS)
             end do
           end if
           close(99999)
           deallocate(HXC)
         end if
  111    continue
       end if
C
C USE OLD HAMILTONIAN OR WRITE NEW ONE
C
       NREC=0
       DO 510 KREP=1,N_REP
        NREC=NREC+(NS_TOT(KREP)*(NS_TOT(KREP)+1))/2
  510  CONTINUE
       IF (ISTSCF.EQ.1) THEN
        PRINT '(A)','USING OLD HAMILTONIAN'
        INQUIRE(FILE='HAMOLD',EXIST=EXIST)
        IF (.NOT.EXIST) THEN
         write(6,*)'NUMHAM: HAMOLD DOES NOT EXIST'
         CALL STOPIT
        END IF
        OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='OLD')
        REWIND(99)
        READ(99,END=530)NREC,MSPN
        IF (MSPN.EQ.1) THEN
         READ(99,END=530)(HSTOR(IREC,2),IREC=1,NREC)
        ELSE IF (MSPN.EQ.2) THEN
         DO 520 I=1,ISPN
          READ(99,END=530)(HSTOR(IREC,2),IREC=1,NREC)
  520    CONTINUE
        ELSE
         GOTO 530
        END IF
        GOTO 540
  530   write(6,*)'NUMHAM: HAMOLD UNREADABLE'
        CALL STOPIT
  540   CLOSE(99)
       ELSE
        OPEN(99,FILE='HAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        IF (ISPN.EQ.1) THEN
         WRITE(99)NREC,NSPN
         WRITE(99)(HSTOR(IREC,2),IREC=1,NREC)
        ELSE
         READ(99)NDUM
         READ(99)DUMMY
         WRITE(99)(HSTOR(IREC,2),IREC=1,NREC)
        END IF
        CLOSE(99)
       END IF
       HAVEHAM=.TRUE.
       CALL GTTIME(TIME2)
       CALL TIMOUT('POTENTIAL MATRIX ELEMENTS:         ',TIME2-TIME1)
       CLOSE(45,STATUS='DELETE')
       RETURN
       END
