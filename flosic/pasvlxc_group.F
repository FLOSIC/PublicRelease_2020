C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASVLXC_GROUP
c
c 08/99 Dirk V. Porezag
c 
c Pasvlxc is called by a worker process to deal with the calculation
c of the exchange-correlation and coulomb potential
c 
#ifdef MPI
       use global_inputs,only : LIBXC1
       use common2,only : NSPN
       use common8,only : REP
       use mpidat1,only : IRANK,SHMRANK,SHMCOMM
       use XTMP2A,only : MIXINS,ISMGGA
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: DTMP , EXCVEC, RHOC, RTMP, VLOC, XTMP
!      REAL*8 :: DTMP , EXCVEC, RHOC, RHOV, RTMP, VLOC, VLOS, VXCS,
!    & XTMP
       INCLUDE 'mpif.h'
       SAVE
       !Replaced common
       REAL*8,ALLOCATABLE :: RHOV(:),VXCS(:),VLOS(:)
!      COMMON/TMP2/RHOV(10*MXSPN*MPBLOCK),VXCS(MXSPN*MPBLOCK)
!    &  ,VLOS(MPBLOCK),RHOC(10,MPBLOCK),XTMP(3,MPBLOCK),DTMP(3,MPBLOCK)
!    &  ,RTMP(MPBLOCK),VLOC(NSPEED)
       INTEGER TAG,MPTS,LPTS,NSIZE,IERR
       INTEGER ITRANS(6),IRECVSTAT(MPI_STATUS_SIZE)
       DIMENSION EXCVEC(4)
       LOGICAL :: LDFTF
       INTEGER :: MXXD,KXXS,KXXO
!
! Allocate local arrays
!
       ALLOCATE(RHOV(10*MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING RHOV'
       ALLOCATE(VXCS(MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING VXCS'
       ALLOCATE(VLOS(MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR ALLOCATING VLOS'
c
c produce fatal error if I am the manager
c
       IF (SHMRANK.EQ.0) THEN
        write(6,*)'FATAL: PASVLXC CALLED BY MANAGER'
        CALL STOPIT
       END IF               
c
c get data and execute subvlxc
c
       TAG= 701
       CALL MPI_RECV(ITRANS(1),6,MPI_INTEGER,0,TAG,
     &               SHMCOMM,IRECVSTAT,IERR)
       MPTS=  ITRANS(1)
       LPTS=  ITRANS(2)
       NSIZE= ITRANS(3)
       MXXD=  ITRANS(4)
       KXXS=  ITRANS(5)
       KXXO=  ITRANS(6)
       ALLOCATE(MIXINS(4,MPTS*MXSPN),STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'PASVLXC:ERROR ALLOCATING MIXINS'
       TAG= 702
       CALL MPI_RECV(RHOV(1),NSIZE,MPI_DOUBLE_PRECISION,0,TAG,
     &               SHMCOMM,IRECVSTAT,IERR)
       TAG= 703
       CALL MPI_RECV(LDFTF,1,MPI_LOGICAL,0,TAG,
     &               SHMCOMM,IRECVSTAT,IERR)
       !SIC orbital calculaton uses NSPN=2 always.
       IF(.NOT. LDFTF) NSPN=2
       CALL SUBVLXC(2,LPTS,MPTS,RHOV,VXCS,VLOS,EXCVEC,
     &              LDFTF,MXXD,KXXS,KXXO)
c
c tell mom that I'm done and report results
c
       TAG=1
       CALL MPI_SSEND(SHMRANK,1,MPI_INTEGER,0,TAG,SHMCOMM,IERR)
       ITRANS(1)= MPTS
       ITRANS(2)= LPTS
       TAG= 705
       CALL MPI_SSEND(ITRANS(1),2,MPI_INTEGER,
     &                0,TAG,SHMCOMM,IERR)
       TAG= 706
       CALL MPI_SSEND(VXCS(1),MPTS*NSPN,MPI_DOUBLE_PRECISION,
     &                0,TAG,SHMCOMM,IERR)
       TAG= 707
       CALL MPI_SSEND(VLOS(1),MPTS,MPI_DOUBLE_PRECISION,
     &                0,TAG,SHMCOMM,IERR)
       TAG= 708
       CALL MPI_SSEND(EXCVEC(1),4,MPI_DOUBLE_PRECISION,
     &                0,TAG,SHMCOMM,IERR)
       IF(LIBXC1.OR.ISMGGA)THEN
         TAG= 709
         CALL MPI_SSEND(MIXINS,MPTS*MXSPN*4,MPI_DOUBLE_PRECISION,
     &                0,TAG,SHMCOMM,IERR)
       ENDIF
       DEALLOCATE(MIXINS,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'PASVLXC:ERROR DEALLOCATING MIXINS'
!
! Deallocate local arrays
!
       DEALLOCATE(VLOS,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING VLOS'
       DEALLOCATE(VXCS,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING VXCS'
       DEALLOCATE(RHOV,STAT=IERR)
       IF(IERR/=0)WRITE(6+IRANK,*)'GETVLXC:ERROR DEALLOCATING RHOV'

#endif
       RETURN
      END
