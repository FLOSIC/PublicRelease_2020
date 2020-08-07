! UTEP Electronic Structure Lab (2020)
! convert global index to local index in block-cyclic distribution

subroutine global_2_local(i,np,nb,p,il)

implicit none
integer, intent(in) :: i    ! global array index, input
integer, intent(in) :: np   ! processor array dimension, input
integer, intent(in) :: nb   ! block size, input
integer, intent(out):: p    ! processor array index, output
integer, intent(out):: il   ! local array index, output
integer :: im1   

im1 = i-1
p   = mod((im1/nb),np)
il  = (im1/(np*nb))*nb + mod(im1,nb) + 1

return
end subroutine global_2_local
