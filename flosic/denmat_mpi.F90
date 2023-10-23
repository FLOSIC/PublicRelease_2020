! UTEP Electronic Structure Lab (2020)
subroutine denmat_mpi
#ifdef MPI
#ifdef SCALAPACK
! Raja- ELP, Feb 4, 2017. This is an improvement over the Novemeber 2016 implementation of use of 
! density matrix in Coupot. It is still restricted to C1 symmetry (i.e. NREP =1) 
! It now uses LAPACK routines to build the density matrix, resulting in speed
! up. Also this should allow straightforward implementation of parallel
! computation of Density Matrix and Kinetic Energy using SCALAPACK routines
! instead of LAPACK routines.
! 
!  ToDo - parallelize using SCALAPACK
!  Rajendra Zope - Sat Feb  4 08:34:28 MST 2017
!  
!  Parallelized with SCALPACK 2/15/2017 - CMD

use debug1,only  : tracer

use mpidat1,only : IRANK,LDR,LDC, DESCA, ICONTXT
use common2,only : NSPN
use common5,only : N_OCC, OCCUPANCY, psi_coef
use common8,only : N_REP, NS_TOT
use hstor1,only  : distr_psi_coef  !,local_over
use den_mat,only : distr_dmat !,T_DMAT,tgem_dmat
use mpi

implicit none 

integer, parameter :: dp=kind(1.d0)
integer :: ISPN,N,NUWF,k, i_rep, ierr, k_virt
real(dp) :: FF
real(dp),allocatable :: coef(:,:),t_dmat(:,:)

character(1) ::  transa,transb
real(dp) :: alpha,beta
real(dp) :: timea,timeb,time1

integer :: descfull(10)

IF(N_REP>1) THEN
  WRITE(6,*) 'USE of DENSITY MATRIX in Coupot works only for single representation'
  RETURN
ENDIF

call gttime(time1)

if(IRANK==0)then
  call global_call(38) !bring other processors in
  write(6,*)'IN DENMAT_MPI'
  call tracer('psi_coef sum',0,sum(psi_coef))
end if

N= NS_TOT(1)

!call overlap(3)
if(.not.allocated(distr_dmat))then
  if(IRANK==0)call tracer('allocating distributed dmat')
  allocate(distr_dmat(LDR,LDC,NSPN), stat=ierr)
  if(ierr/=0) call tracer('denmat:error allocating distr_dmat')
end if
          
allocate(coef(LDR,LDC), t_dmat(LDR,LDC), stat=ierr)
if(ierr/=0)call tracer('denmat: error allocating')
coef(:,:)=0.0d0
t_dmat(:,:)=0.0d0

call tracer('done allocating')

!EKINONL2=0.0d0

k_virt=0
do ispn=1,NSPN

  i_rep=1
! do I_REP=1,N_REP

  N=NS_TOT(I_REP)
  nuwf=N_OCC(i_rep,ISPN)

  if(allocated(distr_psi_coef))then
  !  call tracer('distr_psi_coef allocated')
  else
    call tracer('distr_psi_coef not allocated')
    allocate(distr_psi_coef(LDR,LDC,N_REP,NSPN),stat=ierr)
    if(ierr/=0) call tracer('denmat:error allocating distr_psi_coef')
    distr_psi_coef=0.0d0
    call tracer('denmat:calling pdgemr2d',n,real(nuwf,dp))
    call descinit(descfull,n,n,n,n,0,0,icontxt,n,ierr)
    call pdgemr2d(n,nuwf,psi_coef(:,:,ispn,i_rep),1,1,descfull, &
                   distr_psi_coef(:,:,ispn,i_rep),1,1,desca,icontxt)
    call tracer('denmat:post pdgemr2d',ierr,sum(distr_psi_coef))
  endif

  coef(:,:)=distr_psi_coef(:,:,i_rep,ISPN)

!  call tracer('pre p?scal',nuwf)
  do k=1,nuwf
    k_virt=k_virt+1
    ff = dsqrt(occupancy(k_virt))
    call pdscal(N, ff, coef, 1,k, DESCA, 1)
    !call dscal(nbas, ff, PSI_coef(:,k,1,1), 1)
  enddo

  call gttime(timea)

  alpha=1.d0
  beta=0.d0
  transa='N'
  transb='T'
!  call tracer('pre-pdgemm',n)
  
  call PDGEMM(transa,transb,N,N,nuwf,alpha, &
               coef,1,1,DESCA,       &
               coef,1,1,DESCA, BETA, &
             t_dmat,1,1,DESCA)
!  call tracer('post-pdgemm',nuwf)
   !call DGEMM('N','C',Nbas,Nbas,Nwf,ALPHA,PSI_coef(1:nbas,1:nwf,1,1),Nbas, &
                     ! PSI_coef(1:nbas,1:nwf,1,1),Nbas,BETA,tgem_dmat,Nbas)
  distr_dmat(:,:,ISPN)=t_dmat(:,:)

  call gttime(timeb)
  call ptimout('denmat pdgemm 1:',timeb-timea)

! Calculate sparsity of density matrix
  CALL CALC_DMAT_NNZ(ISPN)
enddo
deallocate(coef,T_DMAT,distr_psi_coef)

call gttime(timeb)
call ptimout('TOTAL IN DENMAT_MPI',timeb-time1)

return

#endif !SCALAPACK
#endif !MPI
end subroutine DENMAT_MPI
! end Raja test

