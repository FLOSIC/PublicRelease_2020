C UTEP Electronic Structure Lab (2020)
#ifdef GROUP
      program cluster
      use mpidat1,only : IRANK,MPI_COMM_WORLD,MYGROUP,NPROC
      implicit none
      include 'mpif.h'
      integer   :: key, ierr
      integer   :: i, total_groups, total_members
      integer   :: istatus(MPI_STATUS_SIZE)
      character :: group_txt*10
      character :: dir_name*40
      logical   :: exist_dir
      real*8    :: part

      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,NPROC,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,IRANK,ierr)

      call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)

!LB: These are just default values
      total_groups=1
      total_members=12
!LB
      open(2,file='igroup')
      read(2,*) total_groups
      read(2,*) total_members
      close(2)
  
      if(IRANK.eq.0) then
        do i=1,total_groups
          MYGROUP=i-1
          write(group_txt,'(1I0)') MYGROUP
          group_txt=adjustl(group_txt)
          dir_name='g'//trim(group_txt)//'/'
          inquire(FILE=dir_name,EXIST=exist_dir)
          if(exist_dir)then
            write(6,*) 'Directory ',dir_name,' exists'
          else
            dir_name='mkdir g'//trim(group_txt)
!      write(6,*) 'Creating directory ',dir_name
!      call flush(6)
            call system(dir_name)

! move files to that directory
            dir_name='cp * g'//trim(group_txt)//'/' 
            call system(dir_name)
          endif
        enddo
      endif

      part=IRANK/total_members
      MYGROUP=floor(part)
      key=mod(IRANK,total_members)
!call MPI_Barrier(MPI_COMM_WORLD,ierr)

!write(6,*) 'after barrier This is node',IRANK,'of',NPROC,'I will be',key,'in group',MYGROUP
!call flush(6)
      call MPI_Comm_split(MPI_COMM_WORLD,MYGROUP,key,MPI_COMM_WORLD,ierr)

      call MPI_Comm_size(MPI_COMM_WORLD,NPROC,ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,IRANK,ierr)

      call MPI_Errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN,ierr)

!write(6,*) 'Now split node',IRANK,'of',NPROC,'in group',MYGROUP
!call flush(6)

      write(group_txt,'(1I0)') MYGROUP
      group_txt=adjustl(group_txt)
      dir_name='g'//trim(group_txt)
      call chdir(dir_name)
      if(IRANK.eq.0)then
!  write(6,*) 'Group ',MYGROUP,'will do calculation',MYGROUP+1
!  call flush(6)
        call create_run(MYGROUP+1)
      endif

!if(group.eq.0)then
      call mpiclus
!endif

!write(6,*) 'Now split node',IRANK,'of',NPROC,'in group',MYGROUP
!call flush(6)
!call MPI_Barrier(MPI_COMM_WORLD,ierr)

!if(rank.eq.0)then
      close(6)
!endif

      call MPI_Comm_free(MPI_COMM_WORLD,ierr)

      call MPI_Finalize(ierr)
#else
      subroutine group
#endif
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                   !
! This subroutine is run by manager of each group    !
! to create the RUNS file for that group            !
!                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine create_run(jgroup)
       integer :: jgroup
       open(2,file='RUNS')
       write(2,*) 0,jgroup
       write(2,*) 3,4
       write(2,*) 0
       close(2)
       return
      end
