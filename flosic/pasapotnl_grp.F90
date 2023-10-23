! UTEP Electronic Structure Lab (2020)
subroutine pasapotnl_grp(mode)
#ifdef MPI
use debug1
use common2,only : NSPN
use common5,only : CONVERGENCE
use mpidat1,only : SHM_MANAGER_COMM,SHM_MANAGER_RANK
use SICMAT,only  : SIC,DERSIC!,ZPOT
use FRM,only     : RESULTS,DEBDAX
use mesh1,only    : NMSH
use mpi
include 'PARAMA2'
integer,intent(in) :: mode
integer :: tag,ierr,IRECVSTAT(MPI_STATUS_SIZE)
integer :: mpidata(2),spnx,iorbx
real(8) :: TOTQNUM
if(mode==1)then
  tag=301
! Receive data from global manager
  call mpi_recv(mpidata,2,mpi_integer,0,TAG,shm_manager_comm,irecvstat,ierr)
  spnx=mpidata(1)
  iorbx=mpidata(2)
! Calculate for this orbital/spin
  CALL APOTNL_SIC(TOTQNUM,spnx,IORBX)
! Tell global manager I am done and get ready for next task
  TAG=1
  CALL MPI_SSEND(SHM_MANAGER_RANK,1,MPI_INTEGER,0,TAG,SHM_MANAGER_COMM,IERR)
endif

if(mode==2)then
  tag=0
  call mpi_recv(spnx,1,mpi_integer,0,TAG,shm_manager_comm,irecvstat,ierr)
! CALL MPI_REDUCE(SIC(1,1,SPNX),SIC(1,1,SPNX),MAX_OCC*MAX_OCC,&
!      MPI_DOUBLE_PRECISION,MPI_SUM,0,SHM_MANAGER_COMM,IERR)
! CALL MPI_REDUCE(ZPOT(1,1,SPNX),ZPOT(1,1,SPNX),NMSH*MAX_OCC,&
!      MPI_DOUBLE_PRECISION,MPI_SUM,0,SHM_MANAGER_COMM,IERR)
  CALL MPI_REDUCE(RESULTS(1,1,SPNX),RESULTS(1,1,SPNX),13*MAX_OCC,&
       MPI_DOUBLE_PRECISION,MPI_SUM,0,SHM_MANAGER_COMM,IERR)
  IF(CONVERGENCE) THEN
     CALL TRACER('SIC CONVERGENCE')
     CALL MPI_REDUCE(DEBDAX(1,1,SPNX),DEBDAX(1,1,SPNX),3*MAX_OCC,&
         MPI_DOUBLE_PRECISION,MPI_SUM,0,SHM_MANAGER_COMM,IERR)
     IF(SPNX.EQ.NSPN) THEN
     CALL MPI_REDUCE(DERSIC(1,1,1),DERSIC(1,1,1),3*MAX_OCC*MX_CNT,&
         MPI_DOUBLE_PRECISION,MPI_SUM,0,SHM_MANAGER_COMM,IERR)
     ENDIF
  ENDIF
endif
#endif
return
end subroutine pasapotnl_grp
