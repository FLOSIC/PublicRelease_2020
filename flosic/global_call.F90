! UTEP Electronic Structure Lab (2020)
      subroutine global_call(idx,ARG1)
#ifdef MPI
      use debug1
      use mpidat1,only : IRANK,NPROC,&
                         SHMRANK,SHM_SIZE,NCALLED_GRP,SHMCOMM,&
                         MANAGER_SIZE,NCALLED_MANAGERS,SHM_MANAGER_COMM,&
                         MYGROUP
      use mesh1,only   : NMSH,RMSH,WMSH
      use common2,only : NSPN,IGGA
      use common5,only : NWFS,PSI_COEF,CONVERGENCE
      use common8,only : U_MAT
      use locorb,only  : TMAT,IRBSIC
      use mocorb,only  : NFRM
      use SICMAT,only  : SIC !,ZPOT
      use FRM,only     : RESULTS,DEBDAX
      use coupotmod,only : USE_DDOT
      use mpi
      include 'PARAMA2'
      integer,intent(in) :: idx
      integer,intent(in),optional :: arg1
      integer,parameter :: ROOT=0
      integer :: I,J,K,TAG,IERR,ISPX,IORBX,NTID,EFP_RESULT
      real(8) :: TOTQNUM,TIME1,TIME2

!
! Bring in other processors to this subroutine
! in global communicator mode or group communicator mode
!
      IF(IDX==70) THEN
        IF(SHMRANK==ROOT) THEN
          TAG=0
          DO I=1,SHM_SIZE
!            CALL TRACER('SENDING NOTICE TO',I)
            CALL MPI_SSEND(IDX,1,MPI_INTEGER,I,TAG,SHMCOMM,IERR)
!            CALL TRACER('SENT NOTICE TO',I)
          ENDDO
!          CALL TRACER('FINISHED CALLING PROCESSORS')
        ENDIF
      ELSE
        IF(IRANK==ROOT) THEN
          TAG=0
          DO I=1,NPROC
!            CALL TRACER('SENDING NOTICE TO',I)
            CALL MPI_SSEND(IDX,1,MPI_INTEGER,I,TAG,MPI_COMM_WORLD,IERR)
!            CALL TRACER('SENT NOTICE TO',I)
          ENDDO
!          CALL TRACER('FINISHED CALLING PROCESSORS')
        ENDIF
      ENDIF
!
! Execute accordingly
!
      select case(idx)
!
! Allocate FOMESH and read cutoff value
!
      case (60)
!         CALL TRACER('ALLOCATING FO_MESH')
!         call allocate_fomesh
!         CALL TRACER('FINISHED ALLOCATING FO_MESH')
!
! Read mesh for given orbital or original mesh
!
      case (61)
!        CALL MPI_BCAST(ARG1,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
!        call read_fomesh(arg1)
!
! Initialize fo_mesh array to false at start of mesh calculation
      case(62)
!        CALL TRACER('GLOBAL CALL 62')
!        IF(ALLOCATED(FO_MESH)) THEN
!          CALL TRACER('INITIALIZING FO_MESH TO ZERO')
!          fo_mesh(:)=.false.
!          CALL TRACER('FOMESH_CUTOFF',IRANK,FOMESH_CUTOFF)
!        ELSE
!          CALL TRACER('ERROR:FO_MESH IS NOT ALLOCATED')
!          CALL STOPIT
!        ENDIF
!
! Broadcast data needed for siclm_msh
!
      case(63)
        CALL MPI_BCAST(NWFS,MXSPN,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(PSI_COEF,NDH*MAX_VIRT_PER_SYM*MAX_REP*MXSPN, &
             MPI_DOUBLE_PRECISION,ROOT,MPI_COMM_WORLD,IERR)
        CALL MPI_BCAST(U_MAT,LOCMAX*ISMAX*MAXSYMSALC*3*2,&
             MPI_DOUBLE_PRECISION,ROOT,MPI_COMM_WORLD,IERR)
#ifdef GROUP
!
! Split communicator into groups
!
      case(64)
!        call SHM_COMM_SPLIT
        CALL READ_GROUPS
        CALL GROUP_SPLIT
!        CALL TRACER('IRANK',IRANK)
!        CALL TRACER('SHMRANK',SHMRANK)
!        CALL TRACER('MYGROUP',MYGROUP)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        CALL STOPIT
#endif
!
! Allocate memory for SIC matrix
!
      case(65)
!        call tracer('IN GLOBAL CALL 65')
        if(irank/=ROOT) call allocate_sic_shm
#ifdef GROUP
!
! Split off to work in groups
!
      case(66)
!            CALL TRACER('SHMRANK',SHMRANK)
! only global manager continues here other managers go to waiting room
            IF(SHMRANK==ROOT) THEN
              CALL MPI_BCAST(CONVERGENCE,1,MPI_LOGICAL,ROOT,SHM_MANAGER_COMM,IERR)
              IF(IRANK==ROOT) THEN
!                CALL TRACER('I AM GLOBAL MANAGER')
!                CALL TRACER('STARTING SENNDATA')
                CALL SENDDATA_MANAGERS(101)
                CALL SENDDATA_MANAGERS(102)
                CALL SENDDATA_MANAGERS(103)
                CALL SENDDATA_MANAGERS(104)
                CALL SENDDATA_MANAGERS(105)
                CALL SENDDATA_MANAGERS(208)
                CALL SENDDATA_MANAGERS(209)
                CALL SENDDATA_MANAGERS(202) !Broadcast of LIBXC1 and ISMGGA
                IF(CONVERGENCE) THEN
                  CALL SENDDATA_MANAGERS(211)
                ENDIF
!                CALL TRACER('DONE SENNDATA')
! Loop over orbitals to hand out work
!                CALL TRACER('STARTING ORBITAL LOOP')
                DO IORBX=1,NFRM(ARG1)
!                  CALL GTTIME(TIME1)
!                  CALL TRACER('HANDING OUT ORBITAL',IORBX)
                  CALL CKWORKER_GRP_MANAGER(1,NTID)
!                  CALL TRACER('AFTER CKWORKER_GRP_MANAGER',NTID)
                  IF (NCALLED_MANAGERS .NE. MANAGER_SIZE) THEN
!                    CALL TRACER('CALLING PAMAPOTNL(1)')
!                    CALL GTTIME(TIME2)
!                    CALL TIMOUT('HANDING OUT',TIME2-TIME1)
                    CALL PAMAPOTNL_GRP(1,ARG1,IORBX)
                  ELSE
! No free workers found, manager does this iteration
!                    CALL TRACER('NONE FOUND I WILL DO',IORBX)
!                    CALL APOTNL_SIC(TOTQNUM,ARG1,IORBX)
!                    CALL TRACER('WAITING FOR SOMEONE TO FINISH')
                    CALL CKWORKER_GRP_MANAGER(4,NTID)
!                    CALL GTTIME(TIME2)
!                    CALL TIMOUT('HANDING OUT 2',TIME2-TIME1)
                    CALL PAMAPOTNL_GRP(1,ARG1,IORBX)
                  END IF
                ENDDO
                CALL PAMAPOTNL_GRP(2,ARG1,1)
!                CALL WRITE_SIC
! LOOP finished now tell other managers and own workers to go back to global communicatior
                TAG=0
!                CALL TRACER('TELLING OTHER MANAGERS TO FINISH')
                DO I=1,MANAGER_SIZE
                  CALL MPI_SEND(-1,1,MPI_INTEGER,I,TAG,SHM_MANAGER_COMM,IERR)
                ENDDO
!                CALL TRACER('FINISHED TELLING OTHER MANAGERS')
              ELSE
!                CALL TRACER('I AM A MANAGER,GOING TO CPUHOG_GROUP_MANAGER')
                CALL CPUHOG_GROUP_MANAGER
!                CALL TRACER('FINISHED CPUHOG_GROUP_MANAGER')
              ENDIF
              TAG=0
!              CALL TRACER('TELLING MY WORKERS TO FINISH')
              DO I=1,SHM_SIZE
                CALL MPI_SEND(-1,1,MPI_INTEGER,I,TAG,SHMCOMM,IERR)
              ENDDO
!              CALL TRACER('FINISHED TELLING MY WORKERS')
            ELSE
! All workers go to group waiting room to get orders from their respective manager
!              CALL TRACER('I AM A WORKER,GOING TO CPUHOG_GROUP_WORKER')
              CALL CPUHOG_GROUP_WORKER
!              CALL TRACER('FINISHED CPUHOG_WORKER')
            ENDIF
!            CALL TRACER('FINISHED GLOBAL CALL 66')
#endif
!
! Initialize RESULTS array
!
      CASE(67)
!           WRITE(6+IRANK,*)'ALLOCATED SIC MATRIX',ALLOCATED(SIC)
!           CALL FLUSH(6+IRANK)
!           CALL TRACER('IN GLOBAL_CALL(67)')
           DO I=1,NSPN
            !IF(SHMRANK==0)THEN
!             CALL TRACER('INITIALIZING SIC MATRIX')
              DO J=1,MAX_OCC
               DO K=1,MAX_OCC
                 SIC(K,J,I)=0.0D0
               END DO
              END DO
              !DO J=1,MAX_OCC
              ! ZPOT(:,J,I)=0.0D0
              !END DO
            !ENDIF
!           CALL TRACER('INITIALIZING DEBDAX')
            DO IORBX=1,MAX_OCC
             DO J=1,3
               DEBDAX(J,IORBX,I)=0.0d0
             END DO
            END DO
!           CALL TRACER('INITIALIZING RESULTS')
            DO IORBX=1,MAX_OCC
              DO J=1,13
                RESULTS(J,IORBX,I)=0.0d0
              END DO
            END DO
           END DO
!
! Deallocate memory for SIC matrix
!
      case(68)
        if(IRANK/=ROOT) call deallocate_sic_shm
!
! Allocate shared memory array ACOUL_SHARED
!
      case(69)
           CALL TRACER('ACOUL NOT AVAILABLE WITH THIS COMPILATION')
           CALL STOPIT
      case(70)
!
! Produce error if job number does not exist
!
      case default
           CALL TRACER('ERROR: JOB DOES NOT EXIST')
           call stopit
      end select
!      CALL TRACER('FINISHED GLOBAL_CALL',IDX)
#endif
      return

      end subroutine global_call
