! UTEP Electronic Structure Lab (2020)
subroutine convert_dmat(mode)
use den_mat,only : DMAT,DMAT_1D
use common2,only : NSPN
implicit none
integer,intent(in) :: mode
integer :: i,j,k,l,dmat_size1,dmat_size2,ierr
dmat_size1=size(dmat,1)
dmat_size2=size(dmat,2)

if(mode==0)then
  allocate(dmat_1d(dmat_size1*dmat_size2*nspn),stat=ierr)
  if(ierr/=0) write(6,*)'ERROR ALLOCATING DMAT_1D'
  write(6,*)'CONVERT_DMAT:ALLOCATED SIZE',dmat_size1*dmat_size2*nspn
endif

l=0
do k=1,nspn
  do j=1,dmat_size2
    do i=1,dmat_size1
      l=l+1
      if(mode==0) then
        dmat_1d(l)=dmat(i,j,k)
      else
        dmat(i,j,k)=dmat_1d(l)
      endif
    enddo
  enddo
enddo

!if(mode==0)then
!  deallocate(dmat,stat=ierr)
!  if(ierr/=0) write(6,*)'ERROR DEALLOCATING DMAT'
!endif

if(mode==1)then
  deallocate(dmat_1d,stat=ierr)
  if(ierr/=0) write(6,*)'ERROR DEALLOCATING DMAT_1D'
endif

return
end subroutine convert_dmat
