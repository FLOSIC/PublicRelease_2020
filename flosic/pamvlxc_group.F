C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PAMVLXC_GROUP(NGRAD,NDIM,RHOG,VXC,VLO,
     &                   LDFTF,MXXD,KXXS,KXXO)
c
c subroutine for parallel evaluation of local and xc potential
c called by manager 
c DV Porezag, 08/99
c
#ifdef MPI
       use global_inputs,only : LIBXC1
       use mesh1,only : nmsh
       use common2,only : ISPN, NSPN, ERGXL, ERGXN, ERGCL, ERGCN
       use mpidat1,only : IRANK,SHM_SIZE,NCALLED_GRP,SHMRANK,SHMCOMM
       use XTMP2A,only : MIXINS,MIXIN,ISMGGA
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:55 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IOFS1, IOFS2, JOB
       REAL*8 :: RHOG , VXC, VLO, EXCVEC
       INCLUDE 'mpif.h'
       SAVE
       INTEGER   NGRAD,NDIM
       DIMENSION RHOG(NDIM,10,MXSPN),VXC(NDIM*MXSPN),VLO(NDIM)
       REAL*8,ALLOCATABLE :: RHOV(:),VXCS(:),VLOS(:)
       DIMENSION EXCVEC(4)
       INTEGER TAG,TID,IERR,IRECVSTAT(MPI_STATUS_SIZE)
       INTEGER NMAST,NSMALL,LPTS,MPTS,MODE,LREF,MDIM,NSIZE
       INTEGER IFAC,IGRAD,IOFS,IPTS,ITRANS(6)
       LOGICAL,INTENT(IN) :: LDFTF
       INTEGER,INTENT(IN) :: MXXD,KXXS,KXXO
c
c produce fatal error if I am a worker
c
       IF (SHMRANK.NE.0) THEN
        write(6,*)'FATAL: PAMVLXC_GRP CALLED BY WORKER'
        CALL STOPIT
       END IF 
!
! NSPN=2 always for SIC orbital loop.
!
!       IF(.NOT. LDFTF) NSPN=2
C
C ALLOCATE LOCAL ARRAYS
C
       ALLOCATE(RHOV(10*MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING RHOV'
       ALLOCATE(VXCS(MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING VXCS'
       ALLOCATE(VLOS(MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING VLOS'
c
c determine maximum number of points for manager, smallest allowed packet size
c     
       NMAST=  MPBLOCK/(3*SHM_SIZE+1)
       !IF(LIBXC1.OR.ISMGGA) NMAST = MPBLOCKXC !YY Possibly needed
       NSMALL= MPBLOCKXC/10 
C
C start looping over mesh points
C mode 1: all points sent out, wait for next worker to complete
C         and continue with mode 5
C      2: workers are available for work
C      3: workers are busy, do some work myself
C      4: workers are busy, wait for next one to finish
C         and continue with mode 5
C      5: get data from workers
C                        
       LPTS= 0
   10  CONTINUE
        MPTS= MIN(MPBLOCKXC,NMSH-LPTS)
        MODE= 0
        IF (MPTS .EQ. 0) THEN
         MODE= 1
         IF (NCALLED_GRP .EQ. 0) GOTO 20
         CALL CKWORKER_GRP(4,TID)
         MODE= 5
        ELSE
         IF (NCALLED_GRP.NE.SHM_SIZE) THEN
          MODE= 2
          NCALLED_GRP=NCALLED_GRP+1
          CALL GETTID_GRP(TID)
         ELSE
          MODE= 3
          MPTS= MIN(NMAST,MPTS)
          !YY test. It may be faster if the manager proc. doesn't do any work.
          IF((LIBXC1.OR.ISMGGA).AND.(SHM_SIZE .NE. 0) ) MPTS= 0

          IF ((MPTS .LT. NSMALL) .AND. (MPTS .LT. NMSH-LPTS)) THEN
           MODE= 4
           CALL CKWORKER_GRP(4,TID)
           MODE= 5
          END IF
         END IF
        END IF
c
c mode 2,3: setup points
c
        IF ((MODE .EQ. 2) .OR. (MODE .EQ. 3)) THEN
         IFAC=NSPN*NGRAD
         DO ISPN=1,NSPN
          DO IGRAD=1,NGRAD
           IOFS= IGRAD+(ISPN-1)*NGRAD
           DO IPTS=1,MPTS
            RHOV(IOFS+IFAC*(IPTS-1))=RHOG(LPTS+IPTS,IGRAD,ISPN)
           END DO
          END DO
         END DO
        END IF
c
c mode 2: send job id and points, update lpts and return to beginning
c                                        
        IF (MODE .EQ. 2) THEN
         JOB= 40
         TAG= 0
         CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,SHMCOMM,IERR)
         NSIZE= MPTS*NSPN*NGRAD
         ITRANS(1)= MPTS
         ITRANS(2)= LPTS
         ITRANS(3)= NSIZE
         ITRANS(4)= MXXD
         ITRANS(5)= KXXS
         ITRANS(6)= KXXO
         TAG= 701
         CALL MPI_SSEND(ITRANS(1),6,MPI_INTEGER,TID,
     &                  TAG,SHMCOMM,IERR)
         TAG= 702
         CALL MPI_SSEND(RHOV(1),NSIZE,MPI_DOUBLE_PRECISION,TID,
     &                  TAG,SHMCOMM,IERR)
         LPTS= LPTS+MPTS
         TAG= 703
         CALL MPI_SSEND(LDFTF,1,MPI_LOGICAL,TID,
     &                  TAG,SHMCOMM,IERR)
         GOTO 10
        END IF
c
c mode 3: call SUBVLXC
c
        LREF= LPTS
        MDIM= MPTS
        IF (MODE .EQ. 3) THEN
         CALL SUBVLXC(2,LPTS,MPTS,RHOV,VXCS,VLOS,EXCVEC,
     &                LDFTF,MXXD,KXXS,KXXO)
         LPTS=LPTS+MPTS
        END IF
c
c mode 5: receive data from worker
c
        IF (MODE .EQ. 5) THEN
         TAG= 705
         CALL MPI_RECV(ITRANS(1),2,MPI_INTEGER,TID,TAG,
     &                 SHMCOMM,IRECVSTAT,IERR)
         MDIM= ITRANS(1)
         LREF= ITRANS(2)
         TAG= 706
         CALL MPI_RECV(VXCS(1),MDIM*NSPN,MPI_DOUBLE_PRECISION,TID,
     &                 TAG,SHMCOMM,IRECVSTAT,IERR)
         TAG= 707
         CALL MPI_RECV(VLOS(1),MDIM,MPI_DOUBLE_PRECISION,TID,
     &                 TAG,SHMCOMM,IRECVSTAT,IERR)
         TAG= 708
         CALL MPI_RECV(EXCVEC(1),4,MPI_DOUBLE_PRECISION,TID,
     &                 TAG,SHMCOMM,IRECVSTAT,IERR)
         IF(LIBXC1.OR.ISMGGA)THEN
           TAG= 709
           CALL MPI_RECV(MIXINS,MDIM*MXSPN*4,MPI_DOUBLE_PRECISION,TID,
     &                 TAG,SHMCOMM,IRECVSTAT,IERR)
         ENDIF
         IF ((LREF .LT. 0) .OR. (LREF+MDIM .GT. NMSH)) THEN
          write(6+IRANK,*)'PAMVLXC: RETURNED INDICES ARE OUT OF BOUNDS'
          write(6+IRANK,*)'LREF, MDIM, IRANK: ',LREF,MDIM,TID
          CALL STOPIT
         END IF
        END IF                               
c
c update arrays
c
        ERGXL=ERGXL+EXCVEC(1)
        ERGXN=ERGXN+EXCVEC(2)
        ERGCL=ERGCL+EXCVEC(3)
        ERGCN=ERGCN+EXCVEC(4)
        IOFS1=(NSPN-1)*NMSH
        IOFS2=(NSPN-1)*MDIM
        DO IPTS=1,MDIM
         VXC(LREF+IPTS)=VXCS(IPTS)
         VXC(LREF+IPTS+IOFS1)=VXCS(IPTS+IOFS2)
         VLO(LREF+IPTS)=VLOS(IPTS)
         IF(LIBXC1.OR.ISMGGA)THEN
! GGA
           MIXIN(1,LREF+IPTS)=MIXINS(1,IPTS)          !Spin up mixin elems
           MIXIN(2,LREF+IPTS)=MIXINS(2,IPTS)
           MIXIN(3,LREF+IPTS)=MIXINS(3,IPTS)
           MIXIN(1,LREF+IPTS+IOFS1)=MIXINS(1,IPTS+IOFS2) !Spin dn mixinelems
           MIXIN(2,LREF+IPTS+IOFS1)=MIXINS(2,IPTS+IOFS2)
           MIXIN(3,LREF+IPTS+IOFS1)=MIXINS(3,IPTS+IOFS2)
! meta-GGA
           MIXIN(4,LREF+IPTS)=MIXINS(4,IPTS)          !Spin up mixin elems
           MIXIN(4,LREF+IPTS+IOFS1)=MIXINS(4,IPTS+IOFS2) !Spin dn mixin elems
         ENDIF

        END DO
        GOTO 10
   20  CONTINUE
C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(RHOV,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING RHOV'
       DEALLOCATE(VXCS,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING VXCS'
       DEALLOCATE(VLOS,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING VLOS'

#endif                      
       RETURN
      END
