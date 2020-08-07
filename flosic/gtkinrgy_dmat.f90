! UTEP Electronic Structure Lab (2020)
subroutine gtkinrgy_dmat
use debug1
use for_diag1,only : AHAM
!use global_inputs,only : SPARSE1
use common2,only : NSPN, EKINONL2
use common8,only : NDMREP, NS_TOT 
use hstor1,only  : hstor
use den_mat,only : dmat 
implicit none

integer, parameter :: dp=kind(1.d0)
integer :: ISPN,N, i_rep, ierr 
real(dp) :: EREP, kintest 
real(dp),allocatable :: COEF(:,:),t_dmat(:,:)

character(1) ::  transa,transb
real(dp) :: alpha,beta
real(dp) :: timeb,time1 
integer  :: i,j,k
!real(8)  :: mat_trace
call gttime(time1)

N=NS_TOT(1)
!allocate(AHAM(n,n),COEF(n,n),t_dmat(n,n),stat=ierr)
allocate(AHAM(n,n),COEF(n,n),stat=ierr)
COEF(:,:)=0.0d0
AHAM(:,:)=0.0d0

!call tracer('distributing hstor')
k=1
do j=1,n
  do i=j,n
    aham(i,j)=hstor(k,1)
    aham(j,i)=aham(i,j)
    k=k+1
  enddo
enddo

EKINONL2=0.0d0

!call tracer('gtkinrgy_dmat:done distributing HAM')

alpha=1.d0
beta=0.d0
transa='N'
transb='T'
I_REP=1
N=NS_TOT(1)
do ispn=1,NSPN
! do I_REP=1,N_REP
!  t_dmat(:,:)=dmat(:,:,ispn)
!  call tracer('moved dmat',ispn)
!  call DGEMM(transa,transb,n,n,n,alpha,t_dmat,n,AHAM,n,beta,COEF,n) ! This works
  call DGEMM(transa,transb,n,n,n,alpha,dmat(1,1,ispn),n,AHAM,n,beta,COEF,n) ! This works
!  call tracer('after dgemm')
!  EREP=mat_trace(COEF) !calculate trace
  EREP=0.0
  DO I=1,N
    EREP=EREP+COEF(I,I)
  ENDDO
!  CALL TRACER('AFTER TRACE',ISPN,EREP)
  KINTEST=EREP*NDMREP(I_REP)
  EKINONL2=EKINONL2+KINTEST
!  CALL TRACER('KINRGY',ISPN,KINTEST)
  
! end do
end do

!deallocate(coef,aham,t_dmat)
deallocate(coef,aham)
call gttime(timeb)
!call timout('GTKINRGY_DMAT',timeb-time1)

return
end subroutine  gtkinrgy_dmat


