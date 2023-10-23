! UTEP Electronic Structure Lab (2020)
subroutine write_phires
use diagv1,only : PHIRES
use for_diag1,only : HAM
use global_inputs,only : INBAS
use mpidat1,only : IRANK
integer :: i,j
character(3) :: rnktxt
character(14) :: filename
#ifdef MPI
write(rnktxt,'(i3.3)')irank
filename='PHIRES_TXT'//trim(rnktxt)
#else
irank=0
filename='PHIRES_TXT'
#endif
open(98+irank,file=filename,form='formatted')
do i=1,INBAS
  do j=1,INBAS
    write(98+irank,*)j,i,phires(j,i)
  end do
enddo
close(98+irank)
!filename='HAM_TXT'//trim(rnktxt)
!open(98+irank,file=filename,form='formatted')
!do i=1,INBAS
!  do j=1,INBAS
!    write(98+irank,*)i,j,ham(i,j)
!  enddo
!enddo
!close(98+irank)
end subroutine write_phires
