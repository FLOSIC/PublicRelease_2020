C UTEP Electronic Structure Lab (2020)
      SUBROUTINE PAMLMSIC_DER(MODE,NSPN,LFM,MPTS,WMSA,SICP,RPTS)
! LFM: change (orbital index -> total number of orbitals)
! SICP: change ( 1d NMAX -> 2d NMAXxMAX_OCC )
#ifdef MPI
      use mpidat1,only : NPROC,NCALLED,IRANK
      use SICMAT,only : SIC,DERSIC,ZMGGAS
      use FRM,only    : BFRM,RESULTS,LFRM,DEBDAX
      USE XTMP2A,ONLY : ISMGGA
      use global_inputs,only : libxc1
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:56 MDT 2017

!     INCLUDE  'PARAMAS'  
      INCLUDE  'PARAMA2'  
      INCLUDE 'mpif.h'
!     INCLUDE 'commons.inc' 
      INTEGER :: IERR,IPROC,IWF,JOB,JWF,LFM,MPTS,NMAX,NSPN
      REAL*8 :: BUF,ADD0,ADD1,ADD2,ADD3,RPTS,SICP,WMSA
      PARAMETER (NMAX=MPBLOCK) 
      INTEGER TID,ITRANS(3),MODE,TAG,IID,LF,IX
      INTEGER IRECVSTAT(MPI_STATUS_SIZE)
      DIMENSION WMSA(NMAX),SICP(NMAX,MAX_OCC),RPTS(3,NMAX)
      DIMENSION ADD1(NMAX),ADD2(NMAX),ADD3(3,NMAX),ADD0(3)

      IF (IRANK.NE.0) THEN
        PRINT *,'FATAL: PAMLMSIC_DER CALLED BY WORKER'
        CALL STOPIT
      END IF
      IF(MODE.EQ.1)THEN
        NCALLED=NCALLED+1
        CALL GTNTID(TID)
        JOB=74 ! CHECK THIS
        TAG=0  ! I AM THE MANAGER (0)
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)

        ITRANS(1)=NSPN
        ITRANS(2)=LFM
        ITRANS(3)=MPTS
        TAG=7401
!       print *, 'calling mpi_ssend ', TAG, ITRANS
        CALL MPI_SSEND(ITRANS(1),3,MPI_INTEGER,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        TAG=7402
!       print *, 'calling mpi_ssend NMAX ', TAG, NMAX
        CALL MPI_SSEND(WMSA(1),NMAX,MPI_DOUBLE_PRECISION,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        TAG=7403
!       print *, 'calling mpi_ssend ', TAG
        CALL MPI_SSEND(SICP(1,1),NMAX*MAX_OCC,MPI_DOUBLE_PRECISION,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        TAG=7404
!       print *, 'calling mpi_ssend ', TAG
        CALL MPI_SSEND(RPTS(1,1),3*NMAX,MPI_DOUBLE_PRECISION,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        TAG=7405
!       print *, 'calling mpi_ssend ', TAG
        CALL MPI_SSEND(ISMGGA,1,MPI_LOGICAL,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        TAG=7406
!       print *, 'calling mpi_ssend ', TAG
        CALL MPI_SSEND(LIBXC1,1,MPI_LOGICAL,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        IF(ISMGGA.OR.LIBXC1)THEN
         TAG=7407
         CALL MPI_SSEND(ZMGGAS,4*NMAX*MAX_OCC,MPI_DOUBLE_PRECISION,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
        ENDIF
!       PRINT*,'I AM RETURNING TO SICLAG_DER',IRANK
!       print *, 'pamlmsic_der ncalled NPROC', ncalled, NPROC
        IF(NCALLED.EQ.NPROC)THEN
C         PRINT*,'WAITING FOR SOMEONE TO FINISH'
          CALL CKWORKER(4,TID)
C         PRINT*,'DONE WAITING'
        END IF  
      END IF
c=============================================
!     PRINT*,'MODE:',MODE
      IF (MODE .EQ. 2 ) THEN
C     PRINT*,'NCALLED:',NCALLED,'SHOULD BE ZERO'
        IF(NCALLED.NE.0)CALL STOPIT
c       IF(LFM.EQ.0)THEN
c         PRINT*,'LFM IS ZERO:',LFM
c         CALL STOPIT
c       END IF
        JOB=75 ! CHECK THIS
        TAG=0  ! I AM THE MANAGER (0)
        DO IPROC=1,NPROC
C        PRINT*,'IPROC',IPROC
         CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
        END DO
        CALL MPI_REDUCE(MPI_IN_PLACE,DERSIC(1,1,1),3*MAX_OCC*MX_CNT,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR)
      END IF
#endif
      RETURN
      END
