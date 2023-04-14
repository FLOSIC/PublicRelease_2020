! UTEP Electronic Structure Lab (2020)
!> @file population.f90
!> @ESLAB at UTEP 01/2018
!> Carlos M Diaz
!
!  <S**2> code by Rajendra Joshi (CMU)
! 
! DOI 10.1002/qua
! 1. Lowdin Population
! Charge on atom A (NA)
!  TOTAL N = tr(S^(1/2) * P * S^(1/2))
!  NA=sum over k (S^(1/2) * P * S^(1/2))kk


! 2. Mulliken (later)
!  TOTAL N = tr(P*S)
!  NA=sum over k (P*S)kk

!> @brief calculates Lowdin population matrix
!> @param[in] n matrix dimension
!> @param[in] D Density matrix
!> @param[inout] S Overlap matrix
!> @param[out] P Population matrix
  SUBROUTINE lowdin_population(n,ispx,D,S,P)
    use blas_module,only : mat_mult
    use common2,only : NSPN
    implicit none
    integer, PARAMETER :: dp = KIND(1.0d0)
    integer,intent(IN)    :: n
    integer,intent(IN)    :: ispx
    real(dp),intent(IN)   :: D(n,n)
    real(dp),intent(INOUT):: S(n,n)
    real(dp),intent(OUT)  :: P(n,n)
    !work matrix
    real(dp) :: T(n,n)

    integer  :: i, j
    real(dp) :: SSq, S1, S2
 
    P(:,:)=0.0_dp
    
!Lowdin
  ! calculate S^(1/2)
  ! LB: Need to do it once, only for spin up, since S is going back with the new value
    if(ispx==1) CALL mat_sqrt(n,S)

  ! S^(1/2) *D
    CALL mat_mult(n,S,D,T)
    !saved in temp matrix T

  ! [S^(1/2) *D] *S^(1/2)
    CALL mat_mult(n,T,S,P)

  END SUBROUTINE lowdin_population
!
!> @brief
!> interface to lowdin_population to write partial charges to file
!
!> param[in] popcase 1 for Mulliken, 2 for Lowdin
  SUBROUTINE write_population(popcase)
    use blas_module,only : mat_mult
    use debug1, only : mat_trace,tracer
    use common2,only : nident,IFUCNT,ZELC,RCNT,nspn
    use common8,only : ns_tot,indbeg
    use hstor1, only : hstor
    use den_mat,only : dmat
    implicit none
    integer, PARAMETER :: dp = KIND(1.0d0)
    integer,intent(IN) :: popcase
    real(dp),allocatable :: S(:,:),P(:,:,:)

    integer :: n,ia,k,istart,iend,ifnct,ierr,ispx
    real(dp):: sumA
    real(dp):: NA(nident,nspn+1)
    real(dp):: Nst(nident,nspn)
    character(LEN=30) :: format1
!
!   Rajendra
!
    integer  :: i, j
    real(dp) :: T(ns_tot(1),ns_tot(1))
    real(dp) :: SSq, S1, S2
!

! Calculate dimension of problem
    n=ns_tot(1)

    IF (.not.allocated(dmat)) THEN
      call denmat_serial
      if (.not.allocated(dmat)) return !something went wrong in denmat
    END IF

    WRITE(6,*)''

! ALLOCATE Overlap and Population matrices
    ALLOCATE(S(n,n),P(n,n,nspn),stat=ierr)
    IF (ierr/=0) THEN
      WRITE(6,*)'write_population: error allocating'
      goto 10   !deallocate and exit
    END IF
    P(:,:,:)=0.0_dp

! fill overlap matrix
    CALL overlap(1)
    CALL packed_2_full(n,hstor,S)

! fill population matrix 
    do ispx=1,nspn
      select case(popcase)
      case(1)
        if(ispx==1) then
          write(6,*)'MULLIKEN POPULATION ANALYSIS'
          call flush(6)
        endif
        CALL mat_mult(n,dmat(1,1,ispx),S,P(1,1,ispx))
    
      case(2)
        if(ispx==1) then
          write(6,*)'LOWDIN POPULATION ANALYSIS'
          call flush(6)
        endif
        CALL lowdin_population(n,ispx,dmat(1,1,ispx),S,P(1,1,ispx))
      end select
      write(6,'(A,I2,F10.5)')' trace of Population matrix spin ',ispx,mat_trace(P(:,:,ispx))
    enddo

!
!   Rajendra Joshi
!   <S**2> values from 
!   peralta et. al, J. Chem. Theory Comput., 2017, 13(12), 6101–6107
!
      S1  = 0.0d0
      S2  = 0.0d0
       Do i = 1,n
         S1 = S1 + 0.75d0*(P(i,i,1) + P(i,i,2))
         Do j = 1,n
           S2 = S2 + 0.25d0*(P(i,i,1)*P(j,j,1) + P(i,i,2)*P(j,j,2)) &
                   - 0.25d0*(P(i,j,1)*P(j,i,1) + P(i,j,2)*P(j,i,2)) &
                   - 0.25d0*(P(i,i,1)*P(j,j,2) + P(i,i,2)*P(j,j,1)) &
                   - 0.50d0*(P(i,j,1)*P(j,i,2) + P(i,j,2)*P(j,i,1))
         End do
       End do
       SSq= S1+S2
       Write(*,*)' <S**2> = ', SSq

!
 
! calculate partial charge
    ! Loop over spin
    DO ISPX=1,NSPN
    ! loop over atoms
      DO ia=1,nident
        istart=indbeg(ia,1)+1 
        IF (ia==nident) THEN
          iend=NS_TOT(1)
        ELSE
          iend=indbeg(ia+1,1)
        END IF
        !sum diagonal elements for atom ia
        sumA=0.0_dp
        DO k=istart,iend
          sumA=sumA+P(k,k,ISPX)
        END DO
        !number of electrons on atom A
        IF (nspn==1) sumA=2*sumA
        NA(ia,ISPX)=sumA

        !partial charge on atom A
        ifnct=IFUCNT(ia)
        NA(ia,ISPX)=abs(ZELC(ifnct))-NA(ia,ISPX)
      END DO
    END DO
    IF(NSPN==2)THEN
      DO IA=1,NIDENT
        IFNCT=IFUCNT(IA)
        Nst(IA,1) = NA(IA,1)
        Nst(IA,2) = NA(IA,2)
        NA(IA,3)=NA(IA,1)+NA(IA,2)
        NA(IA,2)=-NA(IA,1)+ NA(IA,2)
        NA(IA,1)= NA(IA,3) - abs(ZELC(IFNCT))
      ENDDO
    ELSE
      DO IA=1,NIDENT
        IFNCT=IFUCNT(IA)
        Nst(IA,1) = NA(IA,1)
      END DO
    ENDIF

!output partial charges to file
    ! atom x y z charge
    OPEN(12,FILE='POPULATION_ANALYSIS',FORM='FORMATTED')
    WRITE(12,*) 'Natoms = ',nident
    select case(popcase)
    case(1)
      IF(NSPN==1)THEN
        WRITE(12,*)'MULLIKEN POPULATION: ATOM, Z     <n>'
      ELSE
        WRITE(12,*)'MULLIKEN POPULATION: ATOM  Z    <n_up>    <n_down>'
      ENDIF
    case(2) 
      IF(NSPN==1)THEN
        WRITE(12,*)'LOWDIN POPULATION  : ATOM  Z     <n>'
      ELSE
        WRITE(12,*)'LOWDIN POPULATION  : ATOM  Z    <n_up>    <n_down>' 
      ENDIF
    end select
    IF(NSPN==1)THEN
      format1="(20X,I4,1x,I4,1(1X,F10.5))"
    ELSE
      format1="(20X,I4,1x,I4,2(1X,F10.5))"
    ENDIF
    If (NSPN==1) Then 
       Write(6,*)'*******************************'
       Write(6,*)'   Atom         Charge'
       Write(6,*)'*******************************'
    Else 
       Write(6,*)'***********************************'
       Write(6,*)'   Atom         Charge   Spin'
       Write(6,*)'***********************************'
    EndIf
    DO ia=1,nident
      IF(NSPN==1) THEN
        WRITE(12,format1)ia,int(ZELC(IFUCNT(ia))),Nst(ia,1)
        WRITE(6,249)ia, int(ZELC(IFUCNT(ia))),NA(ia,1)
      ELSE
!
        WRITE(12,format1)ia,int(ZELC(IFUCNT(ia))),Nst(ia,2), Nst(ia,1)
        WRITE(6,250)ia, int(ZELC(IFUCNT(ia))),NA(ia,1:2)
      ENDIF
    END DO
249 Format(1x,I4,' (Z=',I4,') ',(1X,F8.5))
250 Format(1x,I4,' (Z=',I4,') ',2(1X,F8.5))
! write total charge
    If (NSPN==1) Then 
       Write(6,*)'*******************************'
    Else 
       Write(6,*)'***********************************'
    EndIf
    WRITE(6,'(A,F10.5)')' Total Charge',sum(NA(:,1))
    WRITE(6,'(A,F10.5)')' Total Spin  ',sum(NA(:,2))
    WRITE(6,*)''
    call flush(6)

!DEALLOCATE work arrays
    DEALLOCATE(S,P,stat=ierr)
10  continue
    deallocate(dmat)

  END SUBROUTINE write_population

