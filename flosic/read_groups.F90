! UTEP Electronic Structure Lab (2020)
subroutine read_groups
#ifdef MPI
use mpidat1,only : IRANK,NPROC,NGROUPS,MEMBERS
use mpi
implicit none
integer :: ierr
logical :: exist1
!LB: These are just default values
NGROUPS=1
MEMBERS=NPROC+1
!LB
INQUIRE(file='igroup',exist=exist1)
IF(EXIST1) THEN
  open(2,file='igroup')
  read(2,*) NGROUPS
  read(2,*) MEMBERS
  close(2)
  if(NGROUPS*MEMBERS/=NPROC+1) write(6,*)'PROBLEM IN GROUP_SPLIT',NPROC+1,NGROUPS*MEMBERS
ELSE
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if(IRANK.eq.0) then
    open(3,file='igroup')
    rewind(3)
    write(3,*) NGROUPS ! ,'Number of groups'
    write(3,*) MEMBERS !,'Members per group'
    close(3)
  end if
ENDIF

return

#endif
end subroutine read_groups
