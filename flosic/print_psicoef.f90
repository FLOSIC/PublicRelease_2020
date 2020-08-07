! UTEP Electronic Structure Lab (2020)
subroutine print_psicoef

use common2,only : NSPN
use common5,only : psi_coef
use common8,only : ns_tot

implicit none
integer :: n,i,j,ispx
n=ns_tot(1)
do ispx=1,NSPN
  do j=1,n
    do i=1,n
      write(20,*)i,j,ispx,psi_coef(i,j,1,ispx)
    enddo
  enddo
enddo

return
end subroutine print_psicoef
