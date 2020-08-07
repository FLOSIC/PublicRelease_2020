! UTEP Electronic Structure Lab (2020)
! Subroutine to perform density matrix mixing
!
! Luis Basurto (2014)
!
SUBROUTINE DMATMIX
 use debug1
 use mixpot1,only  : DMATMIXIN,DMATMIXOUT
 use common2,only : NSPN,NIDENT 
 use common5,only : NDH_TOT
 use common7,only : GAUSS_CUT
 use den_mat,only : DMAT,DMAT_1D

 IMPLICIT NONE
 SAVE
 CHARACTER(7) :: NAMES(3)
 LOGICAL      :: FIRST,AVERAGE,EXIST1
 REAL(8)      :: AVG
 INTEGER      :: INDX,XSPN,MPTS,ITER,IPTS,NSPTS,IERR,NHTOT,SIZE1,SIZE2
 DATA FIRST/.TRUE./
 DATA NAMES/'BROYDEN','KBROY1','KBROY2'/

    NHTOT=size(dmat,1)*size(dmat,2)
    call tracer('in dmat mix',NHTOT)
!    call tracer('from module',int(NDH_TOT,4))
    ITER=0
    CALL TRACER('DMATMIX:STARTING')
    WRITE(6,*)'DMATMIX:FIRST',FIRST
    IF(FIRST)THEN
      FIRST=.FALSE.
!
! READ AVERAGE
!
      CALL TRACER('DMATMIX:IN FIRST SECTION')
      AVG=0.70D0
      AVERAGE=.TRUE.
      OPEN(99,FILE='AVRGDAT',FORM='FORMATTED',STATUS='UNKNOWN')
      REWIND 99
      READ(99,*,END=30)AVG,AVERAGE
30    CONTINUE
      REWIND(99)
      WRITE(99,*)AVG,AVERAGE,' AVG, AVERAGE'
      CLOSE(99)
!
! SET GAUSS_CUT (IN THIS VERSION ONLY NEEDED BY COULOMB1)
!   IS IT NEEDED FOR HAMILTONIAN MIXING? (LB)
!   OR EVEN WORSE: IS IT NEEDED AT ALL?
!   YES. It is neeed (in Coupot I think) (RRZ).
!   not needed in combined_mpi. coupot_section.F assumes 1.0D30 (CMD)
      DO IPTS=1,NIDENT
        GAUSS_CUT(IPTS)=1.0D30
      END DO
!YY. Commenting out return on the ITER=0.  
!    Make an initial DMATMIXOLD 
!      return
      open(99,file='DMATMIXOLD',status='UNKNOWN',form='UNFORMATTED')
      close(99,status='DELETE')

      return
      !NSPTS=NHTOT*NSPN
!      go to 501
    END IF
!
! ALLOCATE MIXING ARRAYS
!
    allocate(DMATMIXIN(NHTOT*NSPN),STAT=IERR)
     if(IERR.NE.0)write(6,*)'DMATMIX:Error allocating DMATMIXIN'
    allocate(DMATMIXOUT(NHTOT*NSPN),STAT=IERR)
     if(IERR.NE.0)write(6,*)'DMATMIX:Error allocating DMATMIXOUT'
! not set up for sparse yet. fills HSTOR using sparse arrays -CMD
!    if(SPARSE1) allocate(HSTOR(NHTOT))

!
! PERFORM DENSITY MATRIX MIXING
!
    NSPTS=NHTOT*NSPN
    IF(AVERAGE) THEN
! The DMATMIXOLD may not exist for the first iteration.
      INQUIRE(FILE='DMATMIXOLD',EXIST=EXIST1)
      IF (.NOT. EXIST1) GO TO 501 
      OPEN(99,FILE='DMATMIXOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(99)
      ITER=0
      READ(99,END=500)MPTS,XSPN,ITER
      IF((MPTS.NE.NSPTS).OR.(NSPN.NE.XSPN))THEN
        WRITE(6,*)'DMATMIX: FILE DMATMIXOLD IS UNUSABLE'
        WRITE(6,*)MPTS, NSPTS, NSPN, XSPN
        CALL STOPIT
      END IF
!     OLD DMAT IS IN DMATMIXIN ARRAY STORED IN DMATMIXOLD FILE
      READ(99)(DMATMIXIN(IPTS), IPTS=1,NSPTS)
!
!     NEW DMAT IS IN DMATMIXOUT ARRAY
      OPEN(96,FILE='DMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
      READ(96) MPTS,XSPN
      IF((MPTS.NE.NHTOT).OR.(NSPN.NE.XSPN))THEN
        WRITE(6,*)'DMATMIX: FILE DMAT IS INCOMPATIBLE'
        WRITE(6,*)MPTS,NHTOT,XSPN,NSPN
        CALL STOPIT
      END IF
      ALLOCATE(DMAT_1D(MPTS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'DMATMIX:ERROR ALLOCATING DMAT_1D'
      DO INDX=1,NSPN
!        WRITE(6,*)'DMATMIX:READING SPIN',INDX,MPTS
        READ(96)(DMAT_1D(IPTS),IPTS=1,MPTS)
        DO IPTS=1,MPTS
           DMATMIXOUT(IPTS+(INDX-1)*MPTS)=DMAT_1D(IPTS)
        END DO
      END DO
      CLOSE(96)
!YY. Trying something here: if(avg .gt. 0.0d0 .AND. ITER .gt. 2)
      IF(AVG.GT.0.0D0)THEN
        PRINT '(A)','BROYDEN MIXING OF DENSITY MATRIX'
        CALL MIXING(3,ITER,AVG,NSPTS,NAMES)
        !write(*,*) 'dmatmix mixing', ITER
!YY. do simple mixing
      ELSE
        PRINT '(A)','SIMPLE LINEAR MIXING OF DENSITY MATRIX'
        DO IPTS=1,NSPTS
          DMATMIXIN(IPTS)=(1.0D0+AVG)*DMATMIXIN(IPTS)- &
                               AVG*DMATMIXOUT(IPTS)
        END DO
      END IF

!       MIXED DENSITY MATRIX IS IN DMAT
      OPEN(96,FILE='DMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(96)MPTS,XSPN
      DO INDX=1,NSPN
        DO IPTS=1,MPTS
          DMAT_1D(IPTS)=DMATMIXIN(IPTS+(INDX-1)*MPTS)
        ENDDO
        WRITE(96)(DMAT_1D(IPTS),IPTS=1,MPTS)
      END DO
      CLOSE(96)
      DEALLOCATE(DMAT_1D)
      ALLOCATE(DMAT_1D(NSPTS))
      DMAT_1D(:)=DMATMIXIN(:)
!
!
500   ITER=ITER+1
!       MIXED DENSITY MATRIX IS IN DMATMIXIN FOR THE NEXT ITERATION
      REWIND(99)
      WRITE(99) NSPTS,NSPN,ITER
      write(*,*) 'DMATMIX: writing', ITER
      WRITE(99)(DMATMIXIN(IPTS), IPTS=1,NSPTS)
      CLOSE(99)

      CALL CONVERT_DMAT(1)
     GOTO 600 !deallocate and write DMAT
!      return
    ENDIF
!     
!  The following should be exectuted when iter is zero (actually 1) since we
!   do not have anything in the DMATMIXOLD yet. So copy the DMAT of the
!   first iteration to DMATMIXOLD.
 501   ITER = 1
      write(6,*) 'after normh 1 in dmatmix(501)'
      OPEN(99,FILE='DMATMIXOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(103,FILE='DMAT',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(99)
      READ(103) MPTS,XSPN
      write(99) NSPTS,NSPN,ITER
      IF(ALLOCATED(DMAT_1D)) DEALLOCATE(DMAT_1D)
      ALLOCATE(DMAT_1D(NSPTS),STAT=IERR)
      IF(IERR/=0) WRITE(6,*)'ERROR ALLOCATING DMAT_1D (501)'
      write(*,*) 'DMATMIX writing1', ITER
      DO INDX=1,NSPN
        READ(103)(DMAT_1D(IPTS),IPTS=1,MPTS)
          DO IPTS=1,MPTS
           DMATMIXIN(IPTS+(INDX-1)*MPTS)=DMAT_1D(IPTS)
!           DMATMIXIN(IPTS)=DMAT_1D(IPTS)
          END DO
      ENDDO
      WRITE(99)(DMATMIXIN(IPTS),IPTS=1,NSPTS)
      CLOSE(99)
      CLOSE(103)

! deallocate arrays
 600  deallocate(DMATMIXOUT,DMATMIXIN)

      call tracer('DMATMIX:finished dmatmix')
    RETURN
END SUBROUTINE DMATMIX
