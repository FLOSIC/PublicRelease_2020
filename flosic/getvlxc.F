C UTEP Electronic Structure Lab (2020)
C
C *************************************************************
C
C GETVLXC DIRK POREZAG AUGUST 1999
C CALCULATES THE LOCAL AND EXCHANGE-CORRELATION POTENTIALS ON THE MESH
C
       SUBROUTINE GETVLXC(NDIM,RHOG,VXC,VLO,LDFTF,MXXD,KXXS,KXXO)
       use mesh1,only : nmsh
       use common2,only : IGGA, ISPN, NSPN, ERGXL, ERGXN, ERGCL, ERGCN
! YY. compute ham. elements if convergence is false. Otherwise, write
! mixin to a file for force calculations later in frcslv.
       use common5,only : CONVERGENCE
       use common8,only : REP, N_REP, NDMREP, NS_TOT
       use global_inputs, only : libxc1
       USE XTMP2A,ONLY : ISMGGA,HSTORSIZE,HSTOR_LIBXC,MIXIN,MIXINS
       USE scaledpzsic,only : SICEXCS,SICEXC,scaledsic
!      use SICMAT,only : ZMGGA
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:51 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'
       INTEGER :: KREP
       REAL*8 :: TIMEA, TIMEB, timing1, timing2
       LOGICAL :: EXIST
 
       INTEGER :: IFAC, IGRAD, IOFS, IOFS1, IOFS2, IPTS, LPTS, MPTS
       INTEGER :: NDIM, IERR, NGRAD
       REAL*8 :: RHOG , VXC, VLO, EXCVEC
       CHARACTER*12 :: MGGASTR
       SAVE
       DIMENSION RHOG(NDIM,KRHOG,MXSPN)
       DIMENSION VXC(NDIM*MXSPN),VLO(NDIM)
C      COMMON/TMP2/RHOV(10*MXSPN*MPBLOCK),VXCS(MXSPN*MPBLOCK)
C    &  ,VLOS(MPBLOCK),RHOC(10,MPBLOCK),XTMP(3,MPBLOCK),DTMP(3,MPBLOCK)
C    &  ,RTMP(MPBLOCK),VLOC(NSPEED)
       REAL*8,ALLOCATABLE :: RHOV(:),VXCS(:),VLOS(:)
       DIMENSION EXCVEC(4)
       LOGICAL ISGGA    ! .true. for gga and mgga
       LOGICAL,INTENT(IN) :: LDFTF
       INTEGER,INTENT(IN) :: MXXD,KXXS,KXXO
!
! SIC orbital calculation uses NSPN=2 always.
!
!       IF(.NOT. LDFTF) NSPN=2
C
C ALLOCATE LOCAL ARRAYS
C
       ALLOCATE(RHOV(10*MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING RHOV'
       RHOV(:)=0.0d0
       ALLOCATE(VXCS(MXSPN*MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING VXCS'
       VXCS(:)=0.0d0
       ALLOCATE(VLOS(MPBLOCKXC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING VLOS'
       VLOS(:)=0.0d0
C
C FIGURE OUT SIZE OF NGRAD
C INITIALIZE ENERGIES
C
! YY. ISGGA is set to .true. for GGA and above the Jacob's ladder
! FYI: IGGA(IXC)= 0:LDA  1:GGA  3:MGGA
       !ISGGA= ((IGGA(1) .EQ. 1) .OR. (IGGA(2) .EQ. 1))
       ISGGA= ((IGGA(1) >= 1) .OR. (IGGA(2) >= 1))
       ISMGGA= ((IGGA(1) >= 3) .OR. (IGGA(2) >= 3))
       NGRAD=1
       IF (ISGGA) NGRAD=KRHOG
       ERGXL= 0.0D0
       ERGXN= 0.0D0
       ERGCL= 0.0D0
       ERGCN= 0.0D0

C LB GET TIMINGS
       CALL GTTIME(TIMEA)
C LB
! Allocate HSTOR_LIBXC with size of upper triangle of a Hamiltonian matrix
       IF(LIBXC1.OR.ISMGGA)THEN
        ALLOCATE(MIXIN(4,NMSH*MXSPN),STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING MIXIN',LDFTF
        MIXIN(:,:)=0.0d0
        ALLOCATE(MIXINS(4,MPBLOCKXC*MXSPN),STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING MIXINS',LDFTF
        MIXINS(:,:)=0.0d0
        HSTORSIZE=0
        DO KREP=1,N_REP
          HSTORSIZE= HSTORSIZE+NS_TOT(KREP)*(NS_TOT(KREP)+1)/2
        ENDDO
        !ALLOCATE(HSTOR_LIBXC(HSTORSIZE*NSPN,2))
        !IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR ALLOCATING HSTOR_LIBXC'
        !HSTOR_LIBXC(:,:)=0.0D0
       ENDIF

       IF(scaledsic)THEN
        ALLOCATE(SICEXC(NMSH),STAT=IERR)
        SICEXC(:)=0.0d0
        ALLOCATE(SICEXCS(MPBLOCKXC),STAT=IERR)
        SICEXCS(:)=0.0d0
       END IF
C
C MPI: CALL PARALLEL VERSION 
C
!YY Since there is mesh screening in apotnl_sic, mpi in getvlxc
!   is disabled.
#ifdef MPI
#ifdef GROUP
       IF(LDFTF) THEN
        IF(LIBXC1.OR.ISMGGA) CALL SENDDATA(202)
        CALL PAMVLXC(NGRAD,NDIM,RHOG,VXC,VLO,
     &               LDFTF,MXXD,KXXS,KXXO)
       ELSE
        IF(LIBXC1.OR.ISMGGA) CALL SENDDATA_GRP(202) 
        CALL PAMVLXC_GROUP(NGRAD,NDIM,RHOG,VXC,VLO,
     &                     LDFTF,MXXD,KXXS,KXXO)
       END IF
#else
       !print *,"PAMVLXC"
       !Broadcast LIBXC1 and ISMGGA if nproc >= 1
       IF(LIBXC1.OR.ISMGGA) CALL SENDDATA(202)  
       CALL PAMVLXC(NGRAD,NDIM,RHOG,VXC,VLO,
     &              LDFTF,MXXD,KXXS,KXXO)
#endif
#else
C
C MOVE RHOG DATA INTO TEMPORARY BLOCK
C CALL SUBVLXC FOR CURRENT SET OF POINTS AND UPDATE ARRAYS
C
       DO LPTS=0,NMSH-1,MPBLOCKXC
        MPTS=MIN(NMSH-LPTS,MPBLOCKXC)
        IFAC= NSPN*NGRAD
        DO ISPN=1,NSPN
         DO IGRAD=1,NGRAD
          IOFS= IGRAD+(ISPN-1)*NGRAD
          DO IPTS=1,MPTS
           RHOV(IOFS+IFAC*(IPTS-1))=RHOG(LPTS+IPTS,IGRAD,ISPN)
          END DO
         END DO
        END DO
        CALL SUBVLXC(2,LPTS,MPTS,RHOV,VXCS,VLOS,EXCVEC,LDFTF,
     &               MXXD,KXXS,KXXO)
        ERGXL=ERGXL+EXCVEC(1)
        ERGXN=ERGXN+EXCVEC(2)
        ERGCL=ERGCL+EXCVEC(3)
        ERGCN=ERGCN+EXCVEC(4)
        IOFS1=(NSPN-1)*NMSH   ! Spin down offset of the entire mesh
        IOFS2=(NSPN-1)*MPTS   ! Spin down offset of chopped mesh
        DO IPTS=1,MPTS
         VXC(LPTS+IPTS)=VXCS(IPTS)  ! Spin up
         VXC(LPTS+IPTS+IOFS1)=VXCS(IPTS+IOFS2)  ! Spin down
         VLO(LPTS+IPTS)=VLOS(IPTS)
! EXC for LSIC 
         IF(scaledsic) SICEXC(LPTS+IPTS)=SICEXCS(IPTS)
         IF(LIBXC1.OR.ISMGGA)THEN
! GGA
           MIXIN(1,LPTS+IPTS)=MIXINS(1,IPTS)  !Spin up mixin elems
           MIXIN(2,LPTS+IPTS)=MIXINS(2,IPTS)
           MIXIN(3,LPTS+IPTS)=MIXINS(3,IPTS)
           MIXIN(1,LPTS+IPTS+IOFS1)=MIXINS(1,IPTS+IOFS2)  !Spin dn mixin elems
           MIXIN(2,LPTS+IPTS+IOFS1)=MIXINS(2,IPTS+IOFS2)
           MIXIN(3,LPTS+IPTS+IOFS1)=MIXINS(3,IPTS+IOFS2)
         
! meta-GGA
           MIXIN(4,LPTS+IPTS)=MIXINS(4,IPTS)  !Spin up mixin elems
           MIXIN(4,LPTS+IPTS+IOFS1)=MIXINS(4,IPTS+IOFS2)  !Spin dn mixin elems
         ENDIF
!         write(*,*) LPTS+IPTS, VXC(LPTS+IPTS), LPTS+IPTS+IOFS1, 
!     &                  VXC(LPTS+IPTS+IOFS1), VLO(LPTS+IPTS)
        END DO
       END DO
C
#endif
!      write(6,*)'getvlxc:libxc1',libxc1,'ismgga',ismgga,
!    &      'CONVERGENCE',CONVERGENCE,'NSPN',nspn

! Compute HAMXC file if you are running a DFT calculation.
! We need a different approach for SIC MGGA.
       IF(LDFTF) THEN
        IF(.NOT.CONVERGENCE) THEN
         IF(LIBXC1 .OR. ISMGGA) THEN
          call gttime(timing1)
! compute hamiltonian elements and save them to HAMXC
          do ispn=1,nspn
           call OVERNUM_LIBXC
          end do
          call gttime(timing2)
          call timout('OVERNUM_LIBXC EXECUTION:           ',
     &         timing2-timing1)
! file HAMXC written in OVERNUM_LIBXC
         ENDIF
! If CONVERGENCE is true save mixin to a file for force calculation.
        ELSE 
         !If the LIBXCFRC already exists, delete it.
         INQUIRE(FILE='LIBXCFRC',EXIST=EXIST)
         IF(EXIST) THEN
          OPEN(12001,FILE='LIBXCFRC',FORM='UNFORMATTED',STATUS='OLD')
          CLOSE(12001,STATUS='DELETE')
         END IF

         !Write MIXIN arrays to file. You use these in Pulay force
         !calculation.
         if(LIBXC1.OR.ISMGGA) then
           open(12001,FILE='LIBXCFRC',FORM='UNFORMATTED',STATUS='NEW')
           do IPTS=1,NMSH*NSPN
             write(12001) MIXIN(1,IPTS),MIXIN(2,IPTS),
     &                    MIXIN(3,IPTS),MIXIN(4,IPTS)
           end do
           close(12001)
         end if
        
        END IF
       ELSE !SIC KXXO
        IF(LIBXC1 .OR. ISMGGA) THEN
         call gttime(timing1)
!        ZMGGA(1,:,KXXO) = MIXIN(1,:)
!        ZMGGA(2,:,KXXO) = MIXIN(2,:)
!        ZMGGA(3,:,KXXO) = MIXIN(3,:)
!        ZMGGA(4,:,KXXO) = MIXIN(4,:)
!File name string
         WRITE(MGGASTR,'(A,I4.4)')'MGGA',KXXO
!Clean up
         INQUIRE(FILE=MGGASTR,EXIST=EXIST)
         IF(EXIST) THEN
          OPEN(106+KXXO,FILE=MGGASTR,FORM='UNFORMATTED',STATUS='OLD')
          CLOSE(106+KXXO,STATUS='DELETE')
         END IF
!Write file
         OPEN(106+KXXO,FILE=MGGASTR,FORM='UNFORMATTED',STATUS='NEW')
         print *,"ZMGGA created",KXXO
         DO IPTS=1,NMSH*NSPN
          WRITE(106+KXXO) MIXIN(1:4,IPTS)
         END DO
         CLOSE(106+KXXO)
         call gttime(timing2)
         call timout('MGGA FILE: ', timing2-timing1)
        END IF
       ENDIF

C LB PRINT TIME
       CALL GTTIME(TIMEB)
       CALL TIMOUT('EXCHANGE-CORRELATION:              ',TIMEB-TIMEA)
C LB
C
C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(RHOV,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING RHOV'
       DEALLOCATE(VXCS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING VXCS'
       DEALLOCATE(VLOS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING VLOS'
       IF(LIBXC1.OR.ISMGGA)THEN
!      IF(LDFTF) THEN
         DEALLOCATE(MIXIN,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING MIXIN'
         DEALLOCATE(MIXINS,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING MIXINS'
!      ENDIF
       !DEALLOCATE(HSTOR_LIBXC,STAT=IERR)
       !IF(IERR/=0)WRITE(6,*)'GETVLXC:ERROR DEALLOCATING HSTOR_LIBXC'
       IF(scaledsic) DEALLOCATE(SICEXCS,STAT=IERR)
       ENDIF
C
       RETURN
       END
