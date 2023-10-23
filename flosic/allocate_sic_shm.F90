! UTEP Electronic Structure Lab (2020)
!5/2017
!routine allocates shared memory window to store SIC matrix
!--shared memory communicators split in shm_comm_split

subroutine allocate_sic_shm

  use common2,only : NSPN
!  use common5,only : PSI_COEF
  use sicmat,only : SIC,SIC_COL!,ZPOT,ZMGGA
  use common8,only : NS_TOT, N_REP
  use mesh1,only : nmsh
#ifdef MPI
  use debug1
  use global_inputs,only : mb_size
  use mpidat1,only : shmcomm,shmrank,shm_manager_comm,shm_manager_rank, &
                     sic_win,irank !,nproc

  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
  use mpi
#endif
  
  INCLUDE 'PARAMA2'
!  implicit none

  integer :: NH,VIRT,I_REP,ierr,MY_COMM,MPI3_SIZE
#ifdef MPI_3
!  integer :: shm_size,manager_size
  real(dp):: mem
 
  integer(KIND=MPI_ADDRESS_KIND) :: winsize
  integer :: disp_unit
  type(C_PTR) :: baseptr
  integer :: arrayshape(1)
#endif

#ifndef MPI_3
! allocate SIC matrix on root
! other processors will allocate private copies in later call
!  NH=0
!  do I_REP=1,N_REP !allocate to largest rep size
!   NH=max(NH,NS_TOT(I_REP)) 
!  end do
!  VIRT=min(NH,MAX_VIRT_PER_SYM)
! VIRT=min(VIRT, max(int(E_UP),int(E_DN))+500 )
#ifdef MPI
  if(irank==0)then
    call global_call(65)
  end if
#endif
#ifdef GROUP
!Previous
! if(shmrank==0)then
!   write(6+IRANK,*)'allocating sic matrix',MAX_OCC,MAX_OCC,NSPN
!   allocate(SIC(MAX_OCC,MAX_OCC,NSPN),stat=ierr)
!   SIC(:,:,:)=0.0d0
!   allocate(ZPOT(NMSH,MAX_OCC,NSPN),stat=ierr)
!   ZPOT(:,:,:)=0.0d0
! endif
!After - Manager proc need ZPOT and everybody needs SIC
  if(shmrank==0) then
   !Let manager proc print out
   write(6+IRANK,*)'allocating sic matrix',MAX_OCC,MAX_OCC,NSPN
  endif
  allocate(SIC(MAX_OCC,MAX_OCC,NSPN),stat=ierr)
  SIC(:,:,:)=0.0d0
! if(shmrank==0)then
!   allocate(ZPOT(NMSH,MAX_OCC,NSPN),stat=ierr)
!   ZPOT(:,:,:)=0.0d0
!   allocate(ZMGGA(4,NMSH*MXSPN,MAX_OCC),stat=ierr)
!   ZMGGA(:,:,:)=0.0d0
! endif
#else
    allocate(SIC(MAX_OCC,MAX_OCC,NSPN),stat=ierr)
    SIC(:,:,:)=0.0d0
!   allocate(ZPOT(NMSH,MAX_OCC,NSPN),stat=ierr)
!   ZPOT(:,:,:)=0.0d0
!   allocate(ZMGGA(4,NMSH*MXSPN,MAX_OCC),stat=ierr)
!   ZMGGA(:,:,:)=0.0d0
#endif
! allocate(SIC_COL(MAX_OCC),stat=ierr)
! SIC_COL(:)=0.0d0
  return
#else 

! allocate shared memory SIC matrix for all shared memory groups

  if(irank==0)then
   !brings everyone else in
   call global_call(65)
  end if
!  call tracer('--- allocating shared memory SIC matrix ---')

!  call MPI_Bcast(NS_TOT,N_REP,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  NH=0
!  do I_REP=1,N_REP !allocate to largest rep size
!   NH=max(NH,NS_TOT(I_REP)) 
!  end do
!  VIRT=min(NH,MAX_VIRT_PER_SYM)

! 
! allocate shared memory window
!
  MPI3_SIZE=MAX_OCC
  mem=(8.0d0*MPI3_SIZE) /mb_size 
!  if(irank==0) write(6+irank,1010) mem
 1010  FORMAT('allocating memory: ',F10.2,' MB per node')

  if(shmrank==0)then !SIC matrix(NH,VIRT,N_REP,NSPN), allocated on rank 0
    winsize=8_MPI_ADDRESS_KIND*MPI3_SIZE
    allocate(SIC(MAX_OCC,MAX_OCC,NSPN),stat=ierr)
    SIC(:,:,:)=0.0d0
!   allocate(ZPOT(NMSH,MAX_OCC,NSPN),stat=ierr)
!   ZPOT(:,:,:)=0.0d0
  else
    !remaining ranks do not allocate memory
    winsize=0_MPI_ADDRESS_KIND
  end if

#ifdef GROUP
  MY_COMM=shmcomm
#else
  MY_COMM=MPI_COMM_WORLD
#endif

!  write(6+irank,*)'allocating shared',winsize,MAX_OCC,NSPN
  call mpi_win_allocate_shared(winsize,1,MPI_INFO_NULL,MY_COMM,baseptr,sic_win,ierr)
!  write(6+irank,*)'after allocating shared',winsize,MAX_OCC,NSPN
! 
! query for location of SIC matrix array, on rank 0
!                                   0
  call MPI_Win_shared_query(sic_win,0,winsize,disp_unit,baseptr,ierr)
!  write(6+irank,*)'post allocate_shared',shmrank,winsize
! 
! associate Fortran pointer with shared memory
!
  arrayshape=(/ MAX_OCC/)
  call C_F_POINTER(baseptr, SIC_COL, arrayshape) 
!  write(6+irank,*)'SIC matrix allocated and set up'

  call mpi_win_lock_all(MPI_MODE_NOCHECK,sic_win,ierr)
! 
! fill memory
!
!  if(shmrank==0) then
!    write(6+irank,*)'zero matrix',irank
!    SIC(:,:,:)=0.0d0
!    write(6+irank,*)'end zero matrix',irank
!  end if
!
! Sync 
!
  call MPI_Barrier(MPI_COMM_WORLD,ierr) !just to sync prints
!  write(6+irank,*)'pre win sync'
  call MPI_Win_sync(sic_win,ierr) !memory fence to sync node exchanges
  call MPI_Barrier (MY_COMM,ierr) !time barrier to make sure all ranks have updated their info
  call MPI_Win_sync(sic_win,ierr) !additional memory fence maybe needed on some platforms
!  write(6+irank,*)'post sync'
!  if(irank==0)call tracer('--- end allocating shared memory SIC matrix ')

  return
#endif
end subroutine

!--Closing calls--
! call mpi_win_unlock_all(psi_win,ierr)
! call MPI_win_free (psi_win,ierr)

! call mpi_comm_free(shm_manager_comm,ierr)
! call mpi_comm_free(shmcomm,ierr)

subroutine deallocate_sic_shm
  use sicmat,only : SIC,SIC_COL,ZPOT,ZMGGA
#ifdef MPI
  if(irank==0)then
    call global_call(68)
  end if
#endif
#ifdef GROUP
  deallocate(SIC,stat=ierr)
! if(shmrank==0)then
!   deallocate(ZPOT,stat=ierr)
!   deallocate(ZMGGA,stat=ierr)
! endif
#else
    deallocate(SIC,stat=ierr)
!   deallocate(ZPOT,stat=ierr)
!   deallocate(ZMGGA,stat=ierr)
#endif
! deallocate(SIC_COL,stat=ierr)
  return

end subroutine
