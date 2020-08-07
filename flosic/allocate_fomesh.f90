! UTEP Electronic Structure Lab (2020)
subroutine allocate_fomesh
use mesh1,only : NMSH,FO_MESH,FOMESH_CUTOFF
use frm,only : LFRM
implicit none
integer :: i,max_orb

OPEN(90,file='FOMSH_CUTOFF',form='FORMATTED')
read(90,*)FOMESH_CUTOFF
close(90)

allocate(fo_mesh(nmsh),stat=i)
if(i/=0) write(6,*)'INIT_FOMSH: ERROR ALLOCATING FO_MESH'
return
end subroutine allocate_fomesh
