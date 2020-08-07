! UTEP Electronic Structure Lab (2020)

! packed to full matrices
subroutine packed_2_full(n,P,A)
  implicit none
  integer,parameter :: dp=kind(1.0d0)
  integer,intent(in)  :: n
  real(dp),intent(in) :: P((n+1)*n/2)
  real(dp),intent(out):: A(n,n)

  integer :: i,j,ii

  ii=1
  do i=1,n
    do j=i,n
      A(j,i)=P(ii)
      A(i,j)=A(j,i)
      ii=ii+1
    enddo
  enddo

end subroutine

subroutine full_2_packed(n,A,P)
  implicit none
  integer,parameter :: dp=kind(1.0d0)
  integer,intent(in)  :: n
  real(dp),intent(in) :: A(n,n)
  real(dp),intent(out):: P((n+1)*n/2)

  integer :: i,j,ii
  
  ii=1
  do i=1,n
    do j=i,n
      P(ii)=A(j,i)
      ii=ii+1
    end do
  end do

end subroutine
