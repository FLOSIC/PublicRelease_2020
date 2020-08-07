! UTEP Electronic Structure Lab (2020)
! The subroutine interfaces LAPACK Packed storage routine.
! It solves the standard eigenvalue problem 
! following old rouitne style for easy comptability and understanding of 
! routines.
!    Luis Basurto, El Paso, TX, Wed Jun 26 2015
!
! *************************************************************************
! 
!     SUBROUTINE DIAG_P_SYM
!     =======================
!     
! *************************************************************************
!
!  Diagge creates an interface to the diagonalization routine ewevge.
!
!     Parameters: 
! 
!     N       (I) :  Dimension of Problem 
!     INITIAL (I) :  Starting index of Hamiltonian (used in representations)
!     IDIM    (I) :  Dimension of HSTOR to diagonalize (1 for Overlap, 2 for Hamiltonian)
!                    (Overlap needs to be diagonalize in testbas.ftn)
!     A       (O) :  Eigenvector matrix (eigenvectors in columns) 
!     Eig     (O) :  Eigenvalues 
!     IEV     (I) :  Should be 0 if eigenvectors are not required
! 
! *************************************************************************
      SUBROUTINE DIAG_P_SYM(N,initial,IDIM,A,Eig,IEV)
       use hstor1,only : hstor
       implicit none
       integer, intent(IN)     ::  N,IEV
       integer, intent(INOUT)  ::  initial
       integer, intent(in)     ::  idim
       real*8, intent(OUT)     ::  A(N,N),Eig(N)

       ! local variables
       CHARACTER           :: JOBZ, UPLO,RANGE
       integer             :: info,lda,liwork,lwork,itype,ierr,il,iu,m,&
                              ngrp,nspn
       real*8              :: vl,vu,abstol,ele_up,ele_dn
       integer,allocatable :: IWORK(:)
       real*8,allocatable  :: WORK(:)
       logical             :: imem1,imem2,exist
! Remove this variable once its global
       integer             :: idiag4

       idiag4=0
       if(initial==0)initial=1

! Solving:  A*C = E*C
       itype=1
       imem1=.TRUE.
       imem2=.TRUE.
       if (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
       if (iev .eq. 1) jobz='V'  ! Both eigenvectors and eigenvalues are to be calculated.

       UPLO='L'      ! Lower triangular matrix
       LDA=N        
       INFO =999

        SELECT CASE(idiag4)
        CASE(0)
            ! Perform work size query
          allocate(WORK(1));
          allocate(IWORK(1));
          LWORk=-1
          LIWORK=-1
          CALL DSPEVD( JOBZ, UPLO, N, HSTOR(initial,idim), Eig, A, LDA, &
                        WORK, LWORK, IWORK, LIWORK, INFO )
          LWORK= INT(WORK(1))
          LIWORK = IWORK(1)
          deallocate(work)
          deallocate(iwork)
          allocate(WORK(LWORK),STAT=ierr)
          allocate(IWORK(LIWORK),STAT=ierr)
          CALL DSPEVD( JOBZ, UPLO, N, HSTOR(initial,idim), Eig, A, LDA, &
                        WORK, LWORK, IWORK, LIWORK, INFO )
        END SELECT
      if (INFO .EQ. 0) then
!          write(6,*) 'Diagonalization is succesful'  ! Remove later
      else
           write(6,*) 'Diagonalization is not succesful'
           write(6,*) 'Err: INFO=',INFO
      endif

      deallocate(WORK); deallocate(IWORK);
      return

      END SUBROUTINE DIAG_P_SYM
!-----------------------------------------------------------------------------------------
