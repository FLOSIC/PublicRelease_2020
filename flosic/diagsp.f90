! UTEP Electronic Structure Lab (2020)
!   The subroutine interface for LAPACK routines. The interface is created 
! following old rouitne style for easy comptability and understanding of 
! routines.
!    Raja, El Paso, TX, Wed Jul 15 13:36:29 CDT 2009
!
! *************************************************************************
! 
!     SUBROUTINE DIAGSP
!     =================
!     
! *************************************************************************
!
!  Diagge creates an interface to the diagonalization routine ewevge.
!
!     Parameters: 
! 
!       NA      (I) :  Dimension of A and B 
!       N       (I) :  Dimension of Problem  
!       A       (I) :  Matrix A (lower triangle)
!               (O) :  Eigenvector matrix (eigenvectors in columns) 
!       Eig    (O) :  Eigenvalues 
!       SC      (-) :  Auxiliary vector 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
! *************************************************************************

     SUBROUTINE DIAGSP(NA,N,A,Eig,AUX,IEV)
       implicit none
       integer, intent(IN)     ::  NA,N,IEV
       real*8, intent(INOUT)   ::  A,AUX
       real*8, intent(OUT)     ::  Eig

       ! local variables
       CHARACTER          JOBZ, UPLO
       integer              :: info,lda,liwork,lwork
       integer, allocatable :: IWORK(:)
       real*8, allocatable  :: WORK(:)
! Solving:  A*C = E*C 


        if (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
        if (iev .eq. 1) jobz='V'  ! Both eigenvectors and eigenvalues are to be calculated.
        ! Remove folllowing block later
!       if (NA   .ne. N) then
!          write(6,*) 'DIAGSP: Note Dimension of matrix and that of problem are different!' 
!       endif

        UPLO='L'      ! Lower triangular matrix
        LDA=NA        

        ! Following LWORK, and LIWORK is to compute size of temperary arrays

        LWORK = 1 + 6*Na + 2*Na*Na + 64  
        allocate(WORK(LWORK));
        LIWORK = 3 + 5*Na  + 10
        allocate(IWORK(LIWORK));

        INFO =999

        call DSYEVD( JOBZ, UPLO, N, A, LDA, Eig, WORK, &
       &                   LWORK, IWORK, LIWORK, INFO )


      if (INFO .EQ. 0) then
!          write(6,*) 'Diagonalization is succesful'  ! Remove later
      else
           write(6,*) 'DIAGSP: Diagonalization is not succesful'
           write(6,*) 'Err: INFO=',INFO
      endif

       deallocate(WORK); deallocate(IWORK);
        return

      END SUBROUTINE DIAGSP
!-----------------------------------------------------------------------------------------
