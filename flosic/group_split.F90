! UTEP Electronic Structure Lab (2020)
subroutine group_split
#ifdef MPI
!#ifdef GROUP
use debug1
use mpidat1,only : IRANK,NGROUPS,MEMBERS,MYGROUP, &
                   SHMRANK,SHM_SIZE,SHMCOMM,INUSE_GRP, &
                   SHM_MANAGER_RANK,MANAGER_SIZE,SHM_MANAGER_COMM,INUSE_MANAGER
use mpi
implicit none
integer, parameter :: dp = kind(1.d0)
integer   :: key, ierr
!logical   :: exist_dir
real(dp)    :: part

NGROUPS=NGROUPS+1
part=IRANK/MEMBERS
MYGROUP=floor(part)
key=mod(IRANK,MEMBERS)

if(irank/=0) then
  MYGROUP=MYGROUP+1
endif

if(MYGROUP==1)THEN
  key=key-1
endif
!write(6,*) 'after barrier This is node',IRANK,'of',NPROC,'I will be',key,'in group',MYGROUP
!call flush(6)
call MPI_Comm_split(MPI_COMM_WORLD,MYGROUP,key,SHMCOMM,ierr)

call MPI_Comm_size(SHMCOMM,SHM_SIZE,ierr)
call MPI_Comm_rank(SHMCOMM,SHMRANK,ierr)
SHM_SIZE=SHM_SIZE-1
call MPI_Errhandler_set(SHMCOMM,MPI_ERRORS_RETURN,ierr)

! Create Communicator for Managers
IF(SHMRANK==0)THEN
  allocate(inuse_grp(shm_size))
  inuse_grp(:)=0
ENDIF
call mpi_Barrier(mpi_comm_world,ierr)
! 
! create communicator for shm managers
!
call mpi_comm_split(mpi_comm_world,shmrank,irank,shm_manager_comm,ierr)
call mpi_comm_rank(shm_manager_comm,shm_manager_rank,ierr)
call mpi_comm_size(shm_manager_comm,manager_size,ierr)
manager_size=manager_size-1

if(irank==0) then
  allocate(inuse_manager(manager_size))
  inuse_manager(:)=0
endif

return
#endif
end subroutine group_split
