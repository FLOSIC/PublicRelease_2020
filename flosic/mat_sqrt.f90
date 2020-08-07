! UTEP Electronic Structure Lab (2020)
!> @author Carlos M Diaz, UTEP
!> @brief
!> calculates the square root of a symmetric matrix: A^(1/2)
!> @detail
!> 1: decompose A=Z*S*Z^T
!>  where: S is a diagonal matrix of the eigenvalues of A
!>     and Z contains the eigenvectors of A
!> 2: sqrt(A) = Z*sqrt(S)*Z^T
!> param[in] n matrix dimension
!> param[inout] A square root of A on output

  subroutine mat_sqrt(n,A)
    use blas_module,only : mat_mult
    implicit none
    integer, PARAMETER :: dp = KIND(1.0d0)
    integer ,intent(IN)    :: n
    real(dp),intent(INOUT) :: A(n,n)

    integer  :: i,ierr
    !work arrays
    real(dp) :: Z(n,n),ZT(n,n),W(n)
    !lapack variables
    character :: jobz,uplo
    integer :: lwork,info
    real(dp),allocatable :: work(:)
    real(dp) :: alpha !,beta

    Z(:,:)=A(:,:)

    !A=Z*S*Z^T
    jobz='V'
    uplo='U'
    !allocate work array
    allocate(work(1),stat=ierr)
    
    !query for optimal work size
    lwork=-1

    call dsyev(jobz,uplo,n,Z,n,W,work,lwork,info) 
    if(info/=0)write(6,*)'dsyev query: error',info
    lwork=int(work(1))
    deallocate(work)
    allocate(work(lwork))

    call dsyev(jobz,uplo,n,Z,n,W,work,lwork,info)
    if(info/=0)write(6,*)'dsyev: error',info
    
    deallocate(work)

    !eigenvectors stored in Z
    !copy to Z^T
    ZT(:,:)=transpose(Z(:,:))

    !Z*sqrt(S): scale column(i) by sqrt(W(i))
    do i=1,n
      alpha=sqrt(W(i))
      call dscal(n,alpha,Z(:,i),1)
    end do

    ![Z*sqrt(S)]*Z^T = sqrt(A)
    call mat_mult(n,Z,ZT,A)

    !sqrt(A) returned in A

  end


