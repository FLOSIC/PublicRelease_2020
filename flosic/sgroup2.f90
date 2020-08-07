! UTEP Electronic Structure Lab (2020)
!	implicit none
	program cluster
	include 'mpif.h'
	integer :: rank, size, key, group, ierr,igroup
	integer :: i, j, k, total_groups, total_members
	integer istatus(MPI_STATUS_SIZE)
	integer :: MPI_COMM_WORLD, test
	character :: group_txt*10
	character :: dir_name*40
	logical :: exist_dir

	call MPI_Init(ierr)
	call MPI_Comm_size(MY_NEW_COM,size,ierr)
	call MPI_Comm_rank(MY_NEW_COM,rank,ierr)

	write(6,*) 'This is initial node',rank,'of',size
	call flush(6)
	group=0
	i=0
	igroup=1
	if(rank.eq.0) then
!	  open(5,file='assign_groups')
!	  read(5) total_groups,total_members
!	  close(5)
	  total_groups=3
	  total_members=4
	  do j=1,total_groups
	    group=j-1
	    write(group_txt,'1I0') group
	    group_txt=adjustl(group_txt)
	    dir_name='g'//trim(group_txt)//'/'
	    inquire(FILE=dir_name,EXIST=exist_dir)
	    if(exist_dir)then
!	      write(6,*) 'Directory ',dir_name,' exists'
	    else
	      dir_name='mkdir g'//trim(group_txt)
	      write(6,*) 'Creating directory ',dir_name
	      call flush(6)
	      call system(dir_name)

! move file to that directory
              dir_name='cp * g'//trim(group_txt)//'/' 
	      call system(dir_name)
	    endif
	    do k=1,total_members
	      key=k-1
	      if(key.ne.0.and.group.ne.0)then
	       write(6,*)'curently at member',key,' of group',group
	       call flush(6)
 	       call MPI_Send(key,1,MPI_INTEGER,i,1,MY_NEW_COM,ierr)
	       call MPI_Send(group,1,MPI_INTEGER,i,1,MY_NEW_COM,ierr)
	      endif
	      i=i+1
	    enddo
	  enddo
	  group=0
	  key=0
	else
	  call MPI_Recv(key,1,MPI_INTEGER,0,1,MY_NEW_COM,istatus,ierr)
	  call MPI_Recv(group,1,MPI_INTEGER,0,1,MY_NEW_COM,istatus,ierr)
	  write(6,*)'This is node',rank,'i got key=',key,'group=',group
	  call flush(6)
	endif

	
	write(6,*) 'before barrier This is node',rank,'of',size
	call flush(6)
	call MPI_Barrier(MY_NEW_COM,ierr)

	write(6,*) 'after barrier This is node',rank,'of',size
	call flush(6)
	call MPI_Comm_split(MY_NEW_COM,group,key,MPI_COMM_WORLD,ierr)

	call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

	write(6,*) 'Now split node',rank,'of',size,'in group',group,'ierr=',ierr,'comm=',MPI_COMM_WORLD
	call flush(6)

	write(group_txt,'1I0') group
	group_txt=adjustl(group_txt)
	dir_name='g'//trim(group_txt)
	call chdir(dir_name)
	if(rank.eq.0)then
	  igroup=group+igroup
	  write(6,*) 'Group ',group,'will do calculation',igroup
	  call create_run(igroup)
	endif

	write(6,*) 'node',rank,'of',size,'group',group,'test',test

!	call mpiclus(group,rank,size,MPI_COMM_WORLD)

	call MPI_Barrier(MPI_COMM_WORLD,ierr)

	call MPI_Comm_free(MPI_COMM_WORLD,ierr)

	call MPI_Finalize(ierr)
	end

	subroutine create_run(jgroup)
	integer :: jgroup

	
	open(2,file='RUNS')
	write(2,*) 0,jgroup
	write(2,*) 3,4
	write(2,*) 0
	close(2)
	return

	end
