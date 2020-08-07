! UTEP Electronic Structure Lab (2020)
!     ------------------------------------------------------------------
!
!     The subroutine interface for LAPACK routines. The interface is
!     created following old rouitne style for easy comptability and 
!     understanding of routines.
!
!     Raja, El Paso, TX, Wed Jul 15 13:36:29 CDT 2009
!
!     ------------------------------------------------------------------
! 
!     SUBROUTINE DIAGSP
!     =================
!     
!     ------------------------------------------------------------------
!
!     Diagge creates an interface to the diagonalization routine ewevge.
!
!     Parameters: 
! 
!       N       (I) :  Dimension of Problem  
!       A       (I) :  Matrix A (lower triangle)
!               (O) :  Eigenvector matrix (eigenvectors in columns) 
!       Eig     (O) :  Eigenvalues 
!       SC      (-) :  Auxiliary vector 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
!     ------------------------------------------------------------------
!
      SUBROUTINE DIAGSP2(N,initial,Eig,IEV)
!
!     ------------------------------------------------------------------
!
      use hstor1,only : hstor
      implicit none
      integer, intent(IN)     ::  N,initial,IEV
!
!     real*8, intent(INOUT)   ::  A,AUX
!
      real*8, intent(OUT)     ::  Eig(N)
!
!     --- LOCAL VARIABLES ---
!
      CHARACTER*1          :: JOBZ, UPLO
      integer              :: info,lda,ldz,liwork,lwork,idiag,IERR
      integer, allocatable :: IWORK(:)
      real*8, allocatable  :: WORK(:),Z(:,:)
      logical              :: IMEM1,imem2
!
!     ------------------------------------------------------------------
!
!     --- SOLVING:  A*C = E*C ---
!
      idiag=2
      IMEM1=.true.
      IMEM2=.true.
!
!     --- Z MATRIX IS ALWAYS NEEDED, ALLOCATE HERE ---
!
      ALLOCATE (Z(N,N))
!
      IF (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
      IF (iev .eq. 1) THEN
         jobz='V'  ! Both eigenvectors and eigenvalues are calculated.
!JUR        ALLOCATE(Z(N,N))
      END IF
!
      UPLO='L'      ! Lower triangular matrix
      LDZ=N
!
!     --- LWORK, and LIWORK compute size of temperary arrays ---
!
      LWORK = 1 + 6*N + 2*N*N + 64  
      ALLOCATE(WORK(LWORK),STAT=IERR);
      IF (IERR.ne.0) IMEM1=.false.
!
      LIWORK = 3 + 5*N  + 10
      ALLOCATE(IWORK(LIWORK),STAT=IERR);
      IF(IERR.ne.0) IMEM2=.false.
      INFO =999
!
!     call DSYEVD( JOBZ, UPLO, N, A, LDA, Eig, WORK, &
!       &                   LWORK, IWORK, LIWORK, INFO )
!
      IF (IMEM1.and.imem2) THEN
        call DSPEVD( JOBZ, UPLO, N, HSTOR(initial,1), Eig, Z, LDZ, WORK, LWORK,&
       &              IWORK, LIWORK, INFO )
         DEALLOCATE(WORK); deALLOCATE(IWORK);
!
      ELSE
        ALLOCATE(WORK(3*N));
!
        call  DSPEV( JOBZ, UPLO, N, HSTOR(initial,1), Eig, Z, LDZ, WORK, INFO )
        DEALLOCATE(WORK) 

      END IF
!
!     WRITE(6,*)'Eigenvalues'
!     do i=1,N
!       WRITE(6,*)eig(i)
!     end do
!
      IF(ALLOCATEd(z)) DEALLOCATE(z)
      IF (INFO .EQ. 0) then
        WRITE(6,*) 'Diagonalization is succesful'  ! Remove later
!
      ELSE
        WRITE(6,*) 'DIAGSP2: Diagonalization is not succesful'
        WRITE(6,*) 'Err: INFO=',INFO
!
      END IF
!
      RETURN
!
!     ------------------------------------------------------------------
!
      END SUBROUTINE DIAGSP2
