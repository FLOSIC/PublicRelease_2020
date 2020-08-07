! UTEP Electronic Structure Lab (2020)
!   The subroutine interface for LAPACK routines. The interface is created 
! following old rouitne style for easy comptability and understanding of 
! routines.
!    Raja, El Paso, TX, Wed Jul 15 13:36:29 CDT 2009
!
! *************************************************************************
! 
!     SUBROUTINE DIAGGE
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
!       B       (I) :  Matrix B (lower triangle)
!               (O) :  R where B = R'*R (upper triangle) 
!       Eig    (O) :  Eigenvalues 
!       SC      (-) :  Auxiliary vector 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
! *************************************************************************
      SUBROUTINE DIAGGE(NA,N,A,B,Eig,AUX,IEV)
       use global_inputs,only : idiag1
       implicit none
       integer, intent(IN)     ::  NA,N,IEV
       real(8), intent(INOUT)   ::  A(Na,Na),B(Na,Na),AUX
       real(8), intent(OUT)     ::  Eig(NA)

       ! local variables
       CHARACTER          JOBZ, UPLO
       CHARACTER(1)        RANGE
       integer        :: info,lda,ldb,liwork,lwork,itype,ngrp,nspn
!  if idiag == 1 use divide and conquer
       integer, allocatable :: IWORK(:)
       real(8), allocatable  :: WORK(:)
! Following variable & two arrarys needed for selected calculation of eigenvalues using dsygvx 
       integer   :: IL,IU,LDZ,M,i,j
       real(8)    :: ABSTOL,VL,VU, ele_up, ele_dn,timeA,timeB
       real(8), allocatable ::    PartialEgv(:,:), PEig(:)
       integer, allocatable ::    IFAIL(:)
       logical :: exist, first=.true.
! Solving:  A*C = E*B*C  for itype =1
        itype=1; 

         CALL GTTIME(TIMEA)
         CALL CHECK_INPUTS
!         INQUIRE(FILE='input',EXIST=EXIST)
!         IF (EXIST) THEN
!           open(9,file='input')
!           read(9,*) idiag
!           close(9)
!         ENDIF
           if (NA < 80) idiag1 = 0
          if (idiag1 < 0 .OR. idiag1 > 3) idiag1 = 1

!            write(6,*) 'Current value of idiag is',idiag

        if (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
        if (iev .eq. 1) jobz='V'  ! Both eigenvectors and eigenvalues are to be calculated.
        ! Remove folllowing block later
        if (NA   .ne. N) then
!          write(6,*) 'Note Matrix and problem dimensions are different!' 
        endif

        UPLO='L'      ! Lower triangular matrix
        LDA=NA        
        LDB=NA        
        ! Following LWORK, and LIWORK is to compute size of temperary arrays
        INFO =99

            call flush(6)
        select case(idiag1)
        case(0)
!          write(6,*) 'Calling DSYGV, QR factorization'
          LWORK=3*NA-1
          allocate(WORK(LWORK));
          call  DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, Eig, WORK,  & 
      &                 LWORK, INFO )
           deallocate(WORK)

!          allocate(work(NA))
!           write(6,*) 'Using old dirk_porezag diagonalization routine'
!           write(6,*) 'LB: Not currently available'
!          call DIAGGE_D(NA,N,A,B,Eig,WORK,IEV)
!          deallocate(work)
        case(1)
           write(6,*) 'Calling DSYGVD, divide and conquer'
           LWORK= 1 + 6*Na + 2*Na*Na + 64  
           allocate(WORK(LWORK));
           LIWORK = 3 + 5*Na  + 10
           allocate(IWORK(LIWORK));
         call DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, Eig, WORK, &
       &                   LWORK, IWORK, LIWORK, INFO )
           deallocate(IWORK);
           deallocate(WORK);
        case(2)
!          write(6,*) 'Calling DSYGV, QR factorization'
          LWORK=3*NA-1
          allocate(WORK(LWORK));
          call  DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, Eig, WORK,  & 
      &                 LWORK, INFO )
           deallocate(WORK)
        case(3)
       
           RANGE='I'
           IL=1
           call ELE_INFO(ngrp,nspn,ele_up,ele_dn)
           IU= max(ele_up,ele_dn)+60
           if (IU > N) IU=N
            ABSTOL=1.0e-05
           LWORK=8*NA
           allocate(WORK(LWORK));
           allocate(IWORK(5*NA));
           allocate(IFAIL(NA));
           LDZ=LDA;
           allocate(PartialEgv(LDZ,IU));
           write(6,*) 'No. of eigenvalues requested',IU
          call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,             &
     &                   VL, VU, IL, IU, ABSTOL, M,  Eig, PartialEgv, LDZ, WORK,             &
     &                   LWORK, IWORK, IFAIL, INFO )
 
           write(6,*) 'No. of eigenvalues found',M
!          do i=1, M;write(6,100) i, PEig(i); enddo
               A=0; 
                do i=1,LDZ
                   do j=1,IU
                       A(i,j) = PartialEgv(i,j)
                   enddo
                 enddo
           deallocate(PartialEgv); 
!          deallocate(PEig);
           deallocate(IFAIL);
           deallocate(IWORK);
           deallocate(WORK);
       end select



      if (INFO .EQ. 0) then
!           write(6,*) 'Diagonalization is succesful'  ! Remove later
      else
           write(6,*) 'Diagonalization is not succesful'
           write(6,*) 'Err: INFO=',INFO
      endif


!          write(6,*)'Returning from diagge_z'
         CALL GTTIME(TIMEB)
          call flush(6)
!          write (6,*) 'Time required for diagonalization', timeB-timeA                                                                          
        return
 100       format(1x,i3,f13.6)
      END SUBROUTINE DIAGGE
!-----------------------------------------------------------------------------------------
