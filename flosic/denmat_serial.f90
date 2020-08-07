! UTEP Electronic Structure Lab (2020)
!> @file denmat_serial.f90
!> Creates the density matrix in serial mode
!> Electronic Structure Lab @ UTEP
!> 07/2016
subroutine denmat_serial
!for testing
use debug1,only : tracer
use global_inputs,only : INBAS
use den_mat,only : dmat
use common2,only : NSPN
use common5,only : n_occ,occupancy,psi_coef
use common8,only : ns_tot

implicit none 

integer, parameter :: dp=kind(1.d0)
integer :: n,nuwf

integer :: k, ierr, k_virt
real(dp) :: FF
real(dp),allocatable :: coef(:,:),t_dmat(:,:) !,S(:,:)

character(1) ::  transa,transb
real(dp) :: alpha,beta
!real(dp) :: sum1 
integer  :: ispx,i_rep
!real(dp) :: timea,timeb,time1

k_virt=0
i_rep=1
n=ns_tot(i_rep)
inbas=n
if(.not.allocated(dmat))then
  allocate(dmat(n,n,NSPN),stat=ierr)
  if(ierr/=0) write(6,*)'ERROR ALLOCATING DMAT'
  dmat=0.0d0
end if

do ispx=1,NSPN

  nuwf=N_OCC(i_rep,ispx)
  allocate(coef(n,nuwf),stat=ierr)
  if(ierr/=0) write(6,*)'ERROR ALLOCATING COEF'
  allocate(t_dmat(n,n),stat=ierr)
  if(ierr/=0) write(6,*)'ERROR ALLOCATING T_DMAT'
  COEF=0.0d0
  t_dmat=0.0d0

  COEF(1:n,1:nuwf)=psi_coef(1:n,1:nuwf,i_rep,ispx)

  do k=1,nuwf
    k_virt=k_virt+1
    ff =sqrt(occupancy(k_virt))
    call dscal(n, ff, coef(:,k), 1)
  enddo

  alpha=1.d0
  beta=0.d0
  transa='N'
  transb='C'
  call DGEMM(transa,transb,n,n,nuwf,alpha, &
               coef,n,       &
               coef,n,beta, &
               t_dmat,n)

! hold real part of density.
  dmat(:,:,ispx)=t_dmat(:,:)

! ------------------------------------  
! testing dmat
!
!  call overlap(1)
!  allocate(S(n,n))
!  call packed_2_full(n,hstor,S)
!  S=matmul(dmat,S)
!  do i=1,n
!   sum1=sum1+S(i,i)
!  end do
!  call tracer('trace of D*S',0,sum1)
!  deallocate(S)
! ------------------------------------  

  deallocate(t_dmat,stat=ierr)
  if(ierr/=0)write(6,*)'ERROR DEALLOCATING T_DMAT'
  deallocate(COEF,stat=ierr)
  if(ierr/=0)write(6,*)'ERROR DEALLOCATING COEF'
! enddo
enddo
! write to file
call write_denmat
return

end subroutine 
! end Raja test

