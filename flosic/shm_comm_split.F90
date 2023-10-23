! UTEP Electronic Structure Lab (2020)
subroutine shm_comm_split
!5/2017 split from allocate_hstor
! used for shared memory:
!   allocate_hstor:    hams,iah,jah - sparse hamiltonian arrays
!   allocate_psi_coef: psi_coef
#ifdef GROUP
  use debug1
  use mpidat1,only : shmcomm,shmrank,shm_manager_comm,shm_manager_rank,irank,SHM_SIZE,MANAGER_SIZE,&
                     inuse_grp,inuse_manager !,nproc
  use mpi
  implicit none
!  include 'mpif.h'

  integer :: ierr
!  integer :: shm_size,manager_size

  shm_manager_rank=-1
  shmrank=-1
  shm_size=-1
  manager_size=-1 
! split communicator for shared memory access

  call mpi_comm_split_type(mpi_comm_world,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,shmcomm,ierr)
  call mpi_comm_rank(shmcomm,shmrank,ierr)
  call mpi_comm_size(shmcomm,shm_size,ierr)
  shm_size=shm_size-1
!  CALL TRACER('SHMRANK',SHMRANK)
!  CALL TRACER('SHM_SIZE',SHM_SIZE)
!  if(irank==0)call tracer('size of shared memory group',shm_size)
  if(shmrank==0) then
    allocate(inuse_grp(shm_size))
    inuse_grp(:)=0
  endif
 
  call mpi_Barrier(mpi_comm_world,ierr)
! 
! create communicator for shm managers
!
  call mpi_comm_split(mpi_comm_world,shmrank,irank,shm_manager_comm,ierr)
  call mpi_comm_rank(shm_manager_comm,shm_manager_rank,ierr)
  call mpi_comm_size(shm_manager_comm,manager_size,ierr)
  manager_size=manager_size-1
!  CALL TRACER('SHM_MANAGER_RANK',SHM_MANAGER_RANK)
!  CALL TRACER('MANAGER_SIZE',MANAGER_SIZE)
!  if(irank==0)call tracer('# of shared memory groups',manager_size)
  if(irank==0) then
    allocate(inuse_manager(manager_size))
    inuse_manager(:)=0
  endif
  
! 
! Sync 
!
  call MPI_Barrier(MPI_COMM_WORLD,ierr) !just to sync prints
  call MPI_Barrier(shmcomm,ierr) !time barrier to make sure all ranks have updated their info

  return
#endif
end subroutine

!--Closing calls--

! call mpi_comm_free(shm_manager_comm,ierr)
! call mpi_comm_free(shmcomm,ierr)

