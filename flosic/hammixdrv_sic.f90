! UTEP Electronic Structure Lab (2020)
!> @brief Subroutine to perform hamiltonian mixing
!
!> @author Luis Basurto (2014)
!
!YY modifying this driver to mix DFA-SIC matrix.
!1. We will mix HAM(NBAS,NBAS) instead of HSTOR(NDHTOT). Symmetry of the matrix
!   is not guaranteed.
!2. We save DFA-SIC in a file SICHAMOLD in the subroutine above.
!3. Save the old iteration data in  HAMMIXOLD.
!4. At the subroutine above, read SICHAMOLD back to HAM.
!
SUBROUTINE HAMMIXDRV_SIC(NHTOT)
 use mixpot1,only  : HAMMIXIN,HAMMIXOUT
! use hstor1,only  : HSTOR
 use common2,only : NSPN,NIDENT
 use common7,only : GAUSS_CUT

 IMPLICIT NONE
 SAVE
 INTEGER,INTENT(IN) :: NHTOT !Should be same as NBAS*NBAS
 CHARACTER*7        :: NAMES(3)
 LOGICAL            :: FIRST,AVERAGE,EXIST
 REAL*8             :: AVG
 INTEGER            :: INDX,XSPN,MPTS,ITER,IPTS,NSPTS
 REAL*8,ALLOCATABLE :: AHAM(:) !Temporary HAM array but in 1D
 DATA FIRST/.TRUE./
 DATA NAMES/'BROYDEN','KBROY1','KBROY2'/

    ITER=0

    IF(FIRST)THEN
      FIRST=.FALSE.
!
! READ AVERAGE
!
      AVG=0.15D0
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
!
! IS IT NEEDED FOR HAMILTONIAN MIXING? (LB)
! OR EVEN WORSE: IS IT NEEDED AT ALL?
! YES. It is neeed (in Coupot I think) (RRZ).
      DO IPTS=1,NIDENT
        GAUSS_CUT(IPTS)=1.0D30
      END DO
          
!YY. Commenting out return on the ITER=0.  
!    Make an initial HAMMIXOLD 
!      return
      open(99,file='HAMMIXOLD',status='UNKNOWN',form='UNFORMATTED')
      close(99,status='DELETE')
      return
      !NSPTS=NHTOT*NSPN
      !go to 501
    END IF
!
! PERFORM HAMILTONIAN MIXING
!
    NSPTS=NHTOT*NSPN
!
! ALLOCATE TEMP HAM ARRAY AND MIXING ARRAYS
!
    ALLOCATE(AHAM(NSPTS))
    ALLOCATE(HAMMIXIN(NSPTS))
    ALLOCATE(HAMMIXOUT(NSPTS))
    AHAM(:)=0.0d0
!
    IF(AVERAGE) THEN
! The HAMMIXOLD may not exist for the first iteration.
      INQUIRE(FILE='HAMMIXOLD',EXIST=EXIST) 
      IF (.NOT. EXIST) GO TO 501 
      OPEN(99,FILE='HAMMIXOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
      REWIND(99)
      ITER=0
      READ(99,END=500)MPTS,XSPN,ITER
      write(*,*) 'hammix reading', MPTS,XSPN,ITER
      write(*,*) 'debug', NSPTS,NSPN
      IF((MPTS.NE.NSPTS).OR.(NSPN.NE.XSPN))THEN
        WRITE(6,*)'HAMMIXDRV: FILE HAMIXOLD IS UNUSABLE'
        !WRITE(6,*)MPTS, NSPTS, NSPN, XSPN
        CALL STOPIT
      END IF
!     OLD HAM IS IN HAMMIXIN ARRAY STORED IN HAMMIXOLD FILE
      READ(99)(HAMMIXIN(IPTS), IPTS=1,NSPTS)
!
! NORMALIZE NEW HAMILTONIAN
!
!      CALL NORMH(1)
!     NEW HAM (NORMALIZED) IS IN HAMMIXOUT ARRAY STORED IN HAMOLD FILE
      OPEN(96,FILE='SICHAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
      READ(96) MPTS,XSPN
      IF((MPTS.NE.NHTOT).OR.(NSPN.NE.XSPN))THEN
        WRITE(6,*)'HAMMIXDRV: FILE SICHAMOLD IS INCOMPATIBLE'
        CALL STOPIT
      END IF
      DO INDX=1,NSPN
        READ(96)(AHAM(IPTS),IPTS=1,MPTS)
        DO IPTS=1,MPTS
          HAMMIXOUT(IPTS+(INDX-1)*MPTS)=AHAM(IPTS)
        END DO
      END DO
      CLOSE(96)
!YY. Trying something here: if(avg .gt. 0.0d0 .AND. ITER .gt. 2)
      IF(AVG.GT.0.0D0)THEN
        PRINT '(A)','BROYDEN MIXING OF HAMILTONIAN'
        CALL MIXING(2,ITER,AVG,NSPTS,NAMES)
        !write(*,*) 'hammix mixing', ITER
!YY. do simple mixing
!     elseif( AVG .GT. 0.0d0 ) then
!       do IPTS=1,NSPTS
!         HAMMIXIN(IPTS)=(1.0d0-AVG)*HAMMIXIN(IPTS)+ &
!                        AVG*HAMMIXOUT(IPTS)
!       end do
      ELSE
        PRINT '(A)','SIMPLE LINEAR MIXING OF HAMILTONIAN'
        DO IPTS=1,NSPTS
          HAMMIXIN(IPTS)=(1.0D0+AVG)*HAMMIXIN(IPTS)- &
                               AVG*HAMMIXOUT(IPTS)
        END DO
      END IF

!       NORMLIZED MIXED HAM IS IN HAMOLD
      OPEN(96,FILE='SICHAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(96)MPTS,XSPN
      DO INDX=1,NSPN
        DO IPTS=1,MPTS
           AHAM(IPTS)=HAMMIXIN(IPTS+(INDX-1)*MPTS)
        ENDDO
        WRITE(96)(AHAM(IPTS),IPTS=1,MPTS)
      END DO
      CLOSE(96)
!
!
500   ITER=ITER+1
!       NORMLIZED MIXED HAM IS IN HAMMIXIN FOR THE NEXT ITERATION
      REWIND(99)
      WRITE(99) NSPTS,NSPN,ITER
      write(*,*) 'hammix writing', ITER
      WRITE(99)(HAMMIXIN(IPTS), IPTS=1,NSPTS)
      CLOSE(99)
!     UNNORMLIZE THE  MIXED HAM IN SICHAMOLD FOR DIAGONALIZATION
!      CALL NORMH(2)
!
! DEALLOCATE
!
      DEALLOCATE(HAMMIXOUT)
      DEALLOCATE(HAMMIXIN)
      DEALLOCATE(AHAM)
      return
    endif
!     
!  THe following should be exectuted when iter is zero (actually 1) since we
!   do not have anything in the HAMIXOLD yet. So copy the HAMOLD of the
!   first iteration to HAMMIXOLD.
 501 ITER = 1
!   NORMALIZE BEFORE STORING IN THE FIRST ITERATION      
!    CALL NORMH(1)
    write(6,*) 'after normh 1'
    OPEN(99,FILE='HAMMIXOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
    OPEN(100,FILE='SICHAMOLD',FORM='UNFORMATTED',STATUS='UNKNOWN')
    REWIND(99)
    READ(100) MPTS,XSPN
    write(99) MPTS*XSPN,XSPN,ITER
    write(*,*) 'hammix writing1', MPTS,XSPN,ITER
    DO INDX=1,NSPN
      READ(100)(AHAM(IPTS),IPTS=1,MPTS)
      DO IPTS=1,MPTS
        HAMMIXIN(IPTS+(INDX-1)*MPTS)=AHAM(IPTS)
      END DO
    ENDDO
    WRITE(99)(HAMMIXIN(IPTS),IPTS=1,NSPTS)
    CLOSE(99)
    CLOSE(100)
!    CALL NORMH(2)
!
! DEALLOCATE
!
    DEALLOCATE(HAMMIXOUT)
    DEALLOCATE(HAMMIXIN)
    DEALLOCATE(AHAM)
    RETURN
END SUBROUTINE HAMMIXDRV_SIC
