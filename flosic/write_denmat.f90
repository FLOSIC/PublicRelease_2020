! UTEP Electronic Structure Lab (2020)
subroutine write_denmat
use debug1
use den_mat,only : DMAT,DMAT_1D
use common2,only : NSPN
implicit none
integer :: i,j,size1d,ierr,offset

!call tracer('calling convert_mat in write_denmat')
call convert_dmat(0)
size1d=size(dmat,1)*size(dmat,2)
if(allocated(dmat_1d)) then
  write(6,*)'IN WRITE_DENMAT:DMAT_1D ALLOCATED'
else
  write(6,*)'IN WRITE_DENMAT:DMAT_1D NOT ALLOCATED'
endif
call flush(6)
open(96,file='DMAT',form='unformatted',status='unknown')
write(96)size1d,nspn
do i=1,NSPN
  offset=(i-1)*size1d
  write(6,*) 'WRITE_DENMAT:WRITING SPIN',i,size1d,offset 
  call flush(6)
  write(96)(dmat_1d(j),j=1+offset,size1d+offset)
enddo
close(96)

deallocate(dmat_1d,stat=ierr)
if(ierr/=0) write(6,*)'ERROR DEALLOCATING DMAT_1D'
return
end subroutine write_denmat
