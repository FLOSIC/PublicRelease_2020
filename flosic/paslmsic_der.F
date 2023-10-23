C UTEP Electronic Structure Lab (2020)
      SUBROUTINE PASLMSIC_DER(MODE)
#ifdef MPI
!     use common2,only : NSPN
      use mpidat1, only : IRANK
      use SICMAT, only : SIC,DERSIC,ZMGGAS
      use FRM, only : BFRM,RESULTS,DEBDAX
      USE XTMP2A,ONLY : ISMGGA
      use global_inputs,only : libxc1
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:56 MDT 2017

!     INCLUDE  'PARAMAS'  
      INCLUDE  'PARAMA2'  
      INCLUDE 'mpif.h'
!     INCLUDE 'commons.inc'
      save
      INTEGER :: IERR,IWF,JWF,LFM,LFN,LSPN,MISTAKES,MPTS,NMAX
      REAL*8 :: APT1,APT2,BUF,RPTS,SICP,WMSA
      PARAMETER (NMAX=MPBLOCK)
      INTEGER TID,ITRANS(3),MODE,TAG,IID,IX
      INTEGER IRECVSTAT(MPI_STATUS_SIZE)
      DIMENSION WMSA(NMAX),SICP(NMAX,MAX_OCC),RPTS(3,NMAX)
!     DIMENSION BUF(MAX_OCC,MAX_OCC)
      DIMENSION MISTAKES(2)
C     COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C     COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
C     COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
C    &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
C    &  ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!      COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!      COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!     &   LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)
       DATA MISTAKES/0,0/
       IF (IRANK.EQ.0) THEN
        PRINT *,'FATAL: PASLMSIC CALLED BY MANAGER'
        CALL STOPIT
       END IF
!      print *, 'in PASLMSIC_DER MODE', MODE

      IF (MODE .EQ. 1 ) THEN
       TAG=7401
       CALL MPI_RECV(ITRANS(1),3,MPI_INTEGER,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
c
c  In this version LSPN = NSPN (total spin) and LFN = LFRM(1)+LFRM(2) = total
c    occupied orbitals
       LSPN = ITRANS(1)
       LFN  = ITRANS(2)
       LFM  =LFN
       MPTS = ITRANS(3)
c      PRINT*,IRANK,'PAS RECV LSPN LFN MPTS:',LSPN,LFN, MPTS
       TAG=7402
       CALL MPI_RECV(WMSA(1),NMAX,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
C      PRINT*,IRANK,'PAS RECV   LFM IS ZERO:',LFM
       TAG=7403
       CALL MPI_RECV(SICP(1,1),NMAX*MAX_OCC,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
C      PRINT*,IRANK,'PAS RECV   LFM IS ZERO:',LFM
       TAG=7404
       CALL MPI_RECV(RPTS(1,1),3*NMAX,MPI_DOUBLE_PRECISION,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
C      PRINT*,IRANK,'PAS RECV   LFM IS ZERO:',LFM
       TAG=7405
!MGGA flag
       CALL MPI_RECV(ISMGGA,1,MPI_LOGICAL,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
!LIBXC flag
       TAG=7406
       CALL MPI_RECV(LIBXC1,1,MPI_LOGICAL,
     &               0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       IF(ISMGGA.OR.LIBXC1)THEN
        TAG=7407
        !allocate array for mgga
        if(.not.allocated(ZMGGAS)) then
         ALLOCATE(ZMGGAS(4,NMAX,MAX_OCC),STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING ZMGGAS'
        end if
        ZMGGAS(:,:,:)=0.0d0
        CALL MPI_RECV(ZMGGAS(1,1,1),4*NMAX*MAX_OCC,MPI_DOUBLE_PRECISION,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       END IF
       CALL GTTIME(APT1)
       IF(LFM.EQ.0)THEN
C        PRINT*,IRANK,'PAS BEFORE LFM IS ZERO:',LFM,LFN
         LFM=MAX(LFN,LFM)
         MISTAKES(1)=MISTAKES(1)+1
         IF(MISTAKES(1).GT.MISTAKES(2))THEN
           MISTAKES(2)=MISTAKES(1)*2
           PRINT*,'LFM MISTAKES:',IRANK,MISTAKES
         END IF
       END IF
!      print *,'paslmsic_der calling siclm_derslv LSPN MPTS', LSPN,MPTS
       CALL SICLM_DERSLV(MPTS,WMSA,SICP,RPTS)
C                 PRINT*,IRANK,'PAS AFTER  LFM IS ZERO:',LFM
       CALL GTTIME(APT2)
C      PRINT*,'I AM RANK:',IRANK,' AND I AM DONE',APT2-APT1,MPTS
       TAG=1
!      print *, 'paslmsic_der calling mpi_ssend irank', IRANK
       CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,IERR)
      END IF

      IF (MODE .EQ. 2 ) THEN
c      PRINT*,IRANK,' REDUCE',LFM,LSPN
       TAG=7501
!      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       CALL MPI_REDUCE(DERSIC(1,1,1),DERSIC(1,1,1),
     &                  3*MAX_OCC*MX_CNT,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                  MPI_COMM_WORLD,IERR)

! Reset the DERSIC array of worker processes for the next iteration or next set
! of points
c        DERSIC(:,:,:)=0.0d0
         DERSIC=0.0d0
       !Deallocate an array for mgga
       IF(allocated(ZMGGAS))THEN
        DEALLOCATE(ZMGGAS,STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING ZMGGAS'
       ENDIF
      END IF
#endif
      RETURN
      END
