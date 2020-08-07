! UTEP Electronic Structure Lab (2020)
subroutine get_denmat_atoms(atom1,atom2)
use debug1
use den_mat,only : DMAT,DMAT_LOCAL
use common8,only : INDBEG,NS_TOT
use common2,only : NIDENT,NSPN
implicit none
integer,intent(in) :: atom1,atom2
integer :: i,j,k,init_row,init_col,mat_size1,mat_size2

init_row=indbeg(atom1,1)
init_col=indbeg(atom2,1)
if(atom1==NIDENT)THEN
  mat_size1=NS_TOT(1)-INDBEG(atom1,1)
else
  mat_size1=INDBEG(atom1+1,1)-INDBEG(atom1,1)
endif
if(atom2==NIDENT)THEN
  mat_size2=NS_TOT(1)-INDBEG(atom2,1)
else
  mat_size2=INDBEG(atom2+1,1)-INDBEG(atom2,1)
endif

!call tracer('init_row',init_row)
!call tracer('init_col',init_col)
!call tracer('mat_size1',mat_size1)
!call tracer('mat_size2',mat_size2)
!call tracer('NSPN',nspn)

do k=1,nspn
  do j=1,mat_size2
    do i=1,mat_size1
      dmat_local(i,j,k)=DMAT(I+INIT_ROW,J+INIT_COL,k)
!      write(17,*)i,j,k,dmat(i,j,k)
    enddo
  enddo
enddo

return

end subroutine get_denmat_atoms
