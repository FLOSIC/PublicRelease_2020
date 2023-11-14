! UTEP Electronic Structure Lab (2020)
module xmol

  type Ion
    integer :: anum
    real*8  :: Rx,Ry,Rz
  end type Ion

  type Ion2
    integer     :: atnum
    character*2 :: atlab
    character*5 :: orblab
  end type Ion2

  type(Ion),allocatable  :: xmol_list(:)
  type(Ion2),allocatable :: atomorb_lab(:)
  integer                :: num_atms
  real*8,parameter       :: au2ang=0.52917721067d0 !CODATA 2014 !0.529177
  real*8,parameter       :: ang2au=1.0d0/au2ang    !1.889726878

 contains
      SUBROUTINE GET_LETTER(NUMBER,LETTER)
!
!     RETURNS THE ATOMIC LABEL OF THE SUPPLIED ATOMIC NUMBER.
!
!     BY RAJENDRA ZOPE AND LUIS BASURTO DEC. 2012.
!
!     -----------------------------------------------------------------
!
      INTEGER,INTENT(IN)      :: NUMBER
      CHARACTER*2,INTENT(OUT) :: LETTER
      CHARACTER*2             :: ELEMENTSLETTER(0:112)
!
      DATA ELEMENTSLETTER  /"X", "H", "He", "Li", "Be", "B", "C",      &
           "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",&
           "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co",   &
           "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", &
           "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",  &
           "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",  &
           "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", &
           "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",  &
           "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", &
           "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",  &
           "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", &
           "Hs", "Mt", "Ds", "Rg", "Cn" /
!
      LETTER = ELEMENTSLETTER(NUMBER)
!
!     -----------------------------------------------------------------
!
      END SUBROUTINE GET_LETTER

end module xmol

module mpidat1

! MPIDAT: used by the MPI version to keep track of CPU usage
! and store communicator in group calculations (MYCOMM)
  integer :: NPROC,NCALLED,IRANK,MYCOMM,MYGROUP,IOUT
  integer,allocatable :: INUSE(:)

! SHARED MEMORY COMMUNICATION
  INTEGER :: NGROUPS,MEMBERS
  INTEGER :: SHMRANK,SHM_SIZE,SHMCOMM,NCALLED_GRP
  INTEGER,ALLOCATABLE :: INUSE_GRP(:)
! SHARED WINDOW
  INTEGER :: SIC_WIN
! SHARE MEMORY MANAGERS COMMUNICATION
  INTEGER :: NCALLED_MANAGERS
  INTEGER :: SHM_MANAGER_RANK,MANAGER_SIZE,SHM_MANAGER_COMM
  INTEGER,ALLOCATABLE :: INUSE_MANAGER(:)
! SHARED ACOUL ARRAY WINDOW
  INTEGER :: ACOUL_WIN
! SHARED ACOUL_DFT ARRAY WINDOW AMONGST MANAGERS
  INTEGER :: ACOUL_DFT_WIN

end module mpidat1

module scaladat1
! Adding for excrewrite. only needs setup once. -cmd
  integer :: nprow,npcol,ICONTXT,myrow,mycol
  integer :: NBU,NB,LDR,LDC,LDC1D
  integer,allocatable :: DESCA(:),DESCB(:),DESCZ(:),descfull(:)
  real*8, allocatable :: Z(:,:),AA(:,:),BB(:,:)
end module scaladat1

module global_inputs

  real*8  :: SCFTOL
  integer :: inbas,iinitial,iiev,idiag1,idiag2,idiag3,calctype1,calc_basis,MAXSCF,mixing1,ITTOT,nwfout,&
             scanmesh1,fod_opt3
  logical :: iimesh,excited1,dosjnt1,wfgrid1,dosoccu1,formfak1,atomsph1,matdipole1,nonscf1,nonscfforces1,&
             dftd31,fragment1,polarizability1,rhogrid1,solvent1,molden1,nbo1,symmetry1,mpi_io1,spnorb1,&
             fixm1,pcm1,wffrm1,dmat1,libxc1,spnpol1,efp1,population1,scaledlbfgs1, &
             uniham1,frozendensity1,fod_loop1,fod_opt1,fod_opt2,symmetrymodule1
  logical :: cexcited,cdosjnt,cwfgrid,cdosoccu,cformfak,catomsph,cmatdipole,cnonscf,cnonscfforces,&
             cdftd3,ccalctype,cfragment,cbasis,cpolarizability,crhogrid,csolvent,cmaxscf,cscftol,cmolden,&
             cnbo,csymmetry,cmpi_io,cspnorb,cfixm,cpcm,cdmat,cmixing,clibxc,cefp
  character*30 :: basis_filename

  real(8),parameter :: mb_size=1024.0*1024.0
  real(8),parameter :: gb_size=1024.0*1024.0*1024.0

end module global_inputs

module pot_dens

  real*8, allocatable :: coulomb(:),rhog(:,:,:)

end module pot_dens

module xtmp1a

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! Nedd PARAMS for MAX_PTS
  real*8  :: phig(max_pts,2)

end module xtmp1a

module xtmp2

  real*8,allocatable :: psig(:,:)

end module xtmp2

module xtmp2a

  REAL*8,ALLOCATABLE :: PSIBR(:,:),PSIB(:),VOL(:),QR(:,:)
! The next variables are used when running libxc
  REAL*8,ALLOCATABLE :: PSIBRA(:,:,:),VOL1(:,:),POTDV(:,:),HSTOR_LIBXC(:,:),PROD(:)
  INTEGER :: HSTORSIZE
  REAL*8,ALLOCATABLE :: TAUTOT(:,:),TAUCHOP(:,:),MIXINS(:,:),MIXIN(:,:),TAUW_TAU(:),BETA(:)
  LOGICAL :: ISMGGA

end module xtmp2a

module xtmp3

  integer,parameter :: MAXTST=2000
  integer,parameter :: MAXNP=6
  real*8,allocatable :: ALTAB(:),PREFAC(:,:),RTST(:,:),&
      XYZWGT(:,:,:),CORRECT(:,:,:),YZDECAY(:,:),SUMMESH(:,:),&
      TABINT(:,:),TABSUM(:,:),TABNRM(:)
  integer,allocatable :: NPTAB(:)
  integer :: NTST

end module xtmp3

module for_diag1

  real*8, allocatable :: aham(:,:),aover(:,:),aeval(:),asc1(:)
  real*8, allocatable :: ham(:,:),over(:,:),EVAL(:),SC1(:),filo(:,:),sc2(:)
  integer             :: iirep,iispn,mvps

end module for_diag1

module hstor1

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! NEED PARAMS FOR NDH_TOT
  real*8  :: HSTOR(NDH_TOT,2)
  integer :: MAXCLUSTERSIZE
  integer,allocatable :: JAH(:),IAH(:),JAO(:),IAO(:)
  integer :: NNZH,NNZO

end module hstor1

module coupdata1

  REAL*8,ALLOCATABLE :: AIV(:,:,:,:),AJV(:,:,:,:),DMTV(:,:,:,:,:),&
                      ALPIV(:,:,:),ALPJV(:,:,:),CENTER(:,:),&
                      ADD(:,:,:),RVECI(:,:),RVECJ(:,:)
  INTEGER,ALLOCATABLE :: NPAIRS(:),IP(:,:)
  INTEGER :: XTRA

end module coupdata1

module formfakm

  real*8,allocatable :: rkpt(:,:),rhokpt(:,:,:)
  integer :: nkpt

end module formfakm

module mixpot1

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! NEED PARAMS FOR MAX_PTS AND MXSPN
  real*8 :: POTIN(MAX_PTS*MXSPN),POTOUT(MAX_PTS*MXSPN)
  real(8),allocatable :: dmatmixin(:),dmatmixout(:)
!YY changed HAMMIXIN/OUT to allocatable to experiment SIC-HAM mixing
!  real*8 :: HAMMIXIN(NDH_TOT*MXSPN),HAMMIXOUT(NDH_TOT*MXSPN)
   real*8,allocatable :: HAMMIXIN(:),HAMMIXOUT(:)

end module mixpot1

module mesh1

! Conversion to implicit none.  Raja Zope Sun Aug 20 10:02:04 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
! NEED PARAMS FOR MAX_PTS
!  real*8,allocatable :: wmsh(:),rmsh(:,:)
  real*8 :: wmsh(max_pts),rmsh(3,max_pts)
  integer :: nmsh,nxmsh
  integer :: nmsh_save

end module mesh1

module zero1

!This module is used to zero out arrays of different types and dimensions
  CONTAINS
  SUBROUTINE ZEROINTARRAY(I_A)

!$  use omp_lib
    INTEGER :: I_A(:)
    INTEGER :: I
!$OMP PARALLEL
!$OMP DO
    DO I=1,SIZE(I_A)
      I_A(I)=0
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE ZEROINTARRAY

  SUBROUTINE ZEROREALARRAY(R_A)

!$  use omp_lib
    REAL*8  :: R_A(:)
    INTEGER :: I
!$OMP PARALLEL
!$OMP DO
    DO I=1,SIZE(R_A)
      R_A(I)=0.0
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE ZEROREALARRAY

  SUBROUTINE ZEROINTMATRIX(I_M,IR,IC)

!$ use omp_lib
   INTEGER :: I_M(:,:)
   INTEGER :: IR,IC,I,J
!$OMP PARALLEL
!$OMP DO
    DO I=1,IR
      DO J=1,IC
        I_M(I,J)=0
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE ZEROINTMATRIX

  SUBROUTINE ZEROREALMATRIX(R_M,IR,IC)

!$ use omp_lib
   REAL*8  :: R_M(:,:)
   INTEGER :: IR,IC,I,J
!$OMP PARALLEL
!$OMP DO
    DO I=1,IR
      DO J=1,IC
        R_M(I,J)=0.0
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL
  END SUBROUTINE ZEROREALMATRIX

  SUBROUTINE ZEROREALCUBE(R_C,I1,I2,I3)
!$ use omp_lib
  REAL*8 :: R_C(:,:,:)
  INTEGER :: I1,I2,I3,I,J,K
!$OMP PARALLEL
!$OMP DO
  DO I=1,I1
    DO J=1,I2
      DO K=1,I3
        R_C(I,J,K)=0.0
      ENDDO
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
  END SUBROUTINE ZEROREALCUBE

end module zero1

module blas_module

  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine DOT_PROD
! implments a vector dot prdouct using BLAS
! dot=X.Y
!
! INPUTS
! N size of problem
! X vector of size (N)
! Y vector of size (N)
! OUTPUT
! dot scalar result for dot product
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dot_prod(N,X,Y,dot)

    implicit none
    integer,intent(in)     :: N
    real*8,intent(in)      :: X(N),Y(N)
    real*8,intent(out)     :: dot
    real*8,external        :: ddot
    integer                :: incx,incy

    incx=1
    incy=1
    dot=ddot(N,X,incx,Y,incy)

  end subroutine dot_prod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine MAT_VET
! implments a matrix vector multiplication using BLAS
! Y=A*X
!
! INPUTS
! N size of problem
! A matrix of size (N,N)
! X vector of size (N)
! OUTPUT
! Y vector of size(N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mat_vet(N,A,X,Y)

    implicit none
    integer,intent(in) :: N
    real*8,intent(in)  :: A(N,N)
    real*8,intent(in)  :: X(N)
    real*8,intent(out) :: Y(N)

    integer   :: lda,incx,incy
    character :: trans
    real*8    :: alpha,beta

    trans='N'
    lda=N
    alpha=1.0
    beta=0.0
    incx=1
    incy=1
    call dgemv(trans,N,N,alpha,A,lda,X,incx,beta,Y,incy)

  end subroutine mat_vet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine MAT_MULT
! implments a matrix-matrix multiplication using BLAS
! vector B is not transposed
! C=alpha*A*B+beta*C
!
! INPUTS
! N size of problem
! A matrix of size (N,N)
! B matrix of size (N,N)
! C matrix of size (N,N)
! OUTPUT
! C matrix of size (N,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mat_mult(N,A,B,C)

    implicit none
    integer,intent(in)   :: N
    real*8,intent(in)    :: A(N,N),B(N,N)
    real*8,intent(inout) :: C(N,N)

    character :: transa,transb
    integer :: lda,ldb,ldc
    real*8 :: alpha,beta

    alpha=1.0
!    beta=1.0
    beta=0.0
    transa='N'
!    transb='T'
    transb='N'
    lda=N
    ldb=N
    ldc=N
!    call PMAT_MUL(N,A,B,C)
    call DGEMM(TRANSA,TRANSB,N,N,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

  end subroutine mat_mult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine TRASPOSE
! implments transpose of a matrix
! A=AT
!
! INPUTS
! N size of problem
! A matrix of size (N,N)
! OUTPUT
! A transposed matrix of size (N,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine traspose(N,A)

    implicit none
    integer,intent(in)   :: N
    real*8,intent(inout) :: A(N,N)
    real*8               :: temp
    integer              :: i,j

    do i=1,N
      do j=i,N
        temp=A(j,i)
        A(j,i)=A(i,j)
        A(i,j)=temp
      enddo
    enddo

  end subroutine traspose
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine TRIPLE_PROD
! implments dot product between two matrices using BLAS
!  index1=something
!  index2=something
!  do i=1,N
!    do i=1,N
!      dor_r=dot_r+A(J,index1)*A(I,index2)*B(J,I)
!    enddo
!  enddo
!
! INPUTS
! N size of problem
! index1
! index2
! A matrix of size (N,N)
! B matrix of size (N,N)
! OUTPUT
! dot_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine triple_prod(N,index1,index2,A,B,dot_r)

!$    use omp_lib
    implicit none
    integer,intent(in) :: N,index1,index2
    real*8,intent(in)  :: A(N,N),B(N,N)
    real*8,intent(out) :: dot_r
    integer            :: i,j,err
    real*8,allocatable :: C(:),D(:)
    real*8             :: dot2

    dot_r=0.0
!$OMP PARALLEL SHARED(dot_r) PRIVATE(dot2)
!$OMP DO REDUCTION(+:dot_r)
    do i=1,N
      call dot_prod(N,A(1:N,index1),B(1:N,i),dot2)
      dot_r=dot_r+dot2*A(i,index2)
    enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine triple_prod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine TRIPLE_PROD
! implments dot product between two matrices using BLAS
!  index1=something
!  index2=something
!  do i=1,N
!    do i=1,N
!      dor_r=dot_r+A(J,index1)*A(I,index2)*B(J,I)
!    enddo
!  enddo
!
! INPUTS
! N size of problem
! index1
! index2
! A matrix of size (N,N)
! B matrix of size (N,N)
! OUTPUT
! dot_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine triple_prod2(N,M,A,B,dot_r)

    implicit none
    integer,intent(in) :: N,M
    real*8,intent(in)  :: A(N,N),B(N,N)
    real*8,intent(out) :: dot_r
    real*8,allocatable :: X(:),Y(:),Z(:)
    integer            :: i,j
    real*8             :: dot2

    allocate(X(N),Y(N),Z(N))
    dot_r=0.0
    do i=1,M
      do j=1,N
        X(j)=A(j,i)
        Y(j)=X(j)
      enddo
      call mat_vet(N,B,X,Z)
      call dot_prod(N,Y,Z,dot2)
      dot_r=dot_r+dot2
    enddo
    deallocate(X,Y,Z)
  end subroutine triple_prod2

end module blas_module

module den_mat
 integer, parameter :: dp = kind(1.d0)
 real(8)              :: TIME_DMT
 logical              :: USE_DMAT
 real(dp),allocatable :: distr_dmat(:,:,:)
 real(dp),allocatable :: dmat(:,:,:), dmat_1d(:)
 real(dp)             :: dmat_local(50,50,2)
! complex(dp),allocatable :: zmat(:,:),d2(:,:)
! real(dp),allocatable :: DMAT_SPARSE(:)
! real(dp),allocatable :: DMAT_SPARSE_UP(:),DMAT_SPARSE_DN(:)
! integer,allocatable  :: DMAT_IA(:),DMAT_JA(:)
! integer,allocatable  :: DMAT_IA_UP(:),DMAT_JA_UP(:)
! integer,allocatable  :: DMAT_IA_DN(:),DMAT_JA_DN(:)

end module den_mat

module normham

  real*8, allocatable :: diag(:)
  integer             :: ndiag

end module normham

module istitl
  real*8,allocatable :: RNUC(:,:)
  real*8,allocatable :: ZALP(:,:)
  real*8 :: AFUDIS,ALONG
  integer, allocatable :: IFNU(:)
  integer, allocatable :: NPOW(:)
  integer :: NNUC
  integer :: MX1D

end module istitl

module acoulmod

include 'PARAMA2'

integer :: MUI,IBASE,LI,IALP,IFNCT,I_LOC_MAX
integer :: MUJ,JBASE,LJ,JALP,JFNCT,J_LOC_MAX
REAL(dp),PARAMETER :: ZED=1.0D-30
logical :: STORAGE_FIRST
end module acoulmod

module coupotmod

include 'PARAMA2'

real(dp),pointer     :: acoul_shared(:,:,:)
real(dp)             :: dm_atm_pair(MAXUNSYM,MAXUNSYM,MXSPN)
real(dp),allocatable :: ACOUL_PAIR(:,:,:)
real(dp),allocatable :: COUL_ALL(:,:,:)
real(dp),pointer     :: ACOUL_DFT(:)
real(dp),allocatable :: ACOUL_ORBITALS(:,:,:,:)
logical              :: use_ddot

end module coupotmod

module densold1

  real(8),allocatable :: PHIG(:,:),PSIG(:,:,:),PTS(:,:),GRAD(:,:,:,:,:),RVECA(:,:)
  LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
  LOGICAL             :: LGGA

end module densold1
