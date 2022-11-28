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
!     N       (I) :  Dimension of Problem 
!     INITIAL (I) :  Starting index of Hamiltonian (used in representations)
!     A       (O) :  Eigenvector matrix (eigenvectors in columns) 
!     Eig     (O) :  Eigenvalues 
!     IEV     (I) :  Should be 0 if eigenvectors are not required
! 
! *************************************************************************
      SUBROUTINE DIAGGE3(N,initial,A,Eig,IEV)
       use global_inputs,only : idiag2
       use hstor1,only : hstor
       use common9,only : old_mode
       implicit none
       integer, intent(IN)     ::  N,IEV
       integer, intent(INOUT)  ::  initial
       real*8, intent(OUT)     ::  A(N,N),Eig(N)

       ! local variables
       CHARACTER           :: JOBZ, UPLO,RANGE
       integer             :: info,lda,liwork,lwork,itype,ierr,il,iu,m,&
                               ngrp,nspn
       real*8              :: vl,vu,abstol,ele_up,ele_dn
       integer,allocatable :: IWORK(:),ifail(:)
       real*8,allocatable  :: WORK(:)
       logical             :: imem1,imem2,exist

       if (old_mode) call check_inputs
!       INQUIRE(FILE='input2',EXIST=EXIST)
!       IF (EXIST) THEN
!         open(9,file='input2')
!         read(9,*) idiag
!         close(9)
!       ENDIF
       if (N < 80) idiag2 = 1
       if (idiag2 < 0 .OR. idiag2 > 3) idiag2 = 1
!       write(6,*) 'Current value of idiag2 is',idiag2
       if(initial==0)initial=1

! Solving:  A*C = E*B*C  for itype =1
        itype=1
        imem1=.TRUE.
        imem2=.TRUE.
        if (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
        if (iev .eq. 1) jobz='V'  ! Both eigenvectors and eigenvalues are to be calculated.

        UPLO='L'      ! Lower triangular matrix
        LDA=N        
        INFO =999

        SELECT CASE(idiag2)
        CASE(0)
          allocate(work(8*N),stat=ierr)
          if(ierr.ne.0)then
            write(6,*)'diagge3:no space for WORK for DSPGVX'
            write(6,*)'DIAGGE3:will not execute DSPGVX'
            write(6,*)'SORRY,BYE'
            CALL STOPIT
          endif
            allocate(iwork(5*N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge2:no space for IWORK DSPGVX'
            allocate(ifail(N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge2:no space for IFAIL DSPGVX'
            range='I'
            vl=0.0
            vu=0.0
            il=1
            call ELE_INFO(ngrp,nspn,ele_up,ele_dn)
            IU=INT( max(ele_up,ele_dn))+400
!             IU=N
            if(iu.gt.n)iu=n
!            write(6,*)'N=',N,'IU=',IU
!            write(6,*)'initial=',initial
!            write(6,*)'size of HSTOR',size(HSTOR,1)
!            write(6,*)'size of A',size(A,1)
!            write(6,*)'size of Eig',size(Eig,1)
            abstol=1.0e-10
!            write(6,*)'DIAGGE2:Calling DSPGVX'
            CALL DSPGVX(ITYPE, JOBZ, RANGE, UPLO, N, HSTOR(initial,2), &
               HSTOR(initial,1), VL, VU, IL, IU, ABSTOL, M, Eig, A, &
               LDA, WORK, IWORK, IFAIL, INFO)
        CASE(1)
        ! Following LWORK, and LIWORK is to compute size of temperary arrays  
          allocate(WORK(1));
          allocate(IWORK(1));
          LWORk=-1
          LIWORK=-1
          CALL DSPGVD( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
     &                HSTOR(initial,1), Eig, A, LDA, &
     &                   WORK, LWORK, IWORK, LIWORK, INFO )
        
          LWORK= INT(WORK(1))
          LIWORK = IWORK(1)
          deallocate(work)
          deallocate(iwork)
          allocate(WORK(LWORK),STAT=ierr)
          if(ierr.ne.0) then
            write(6,*)'DIAGGE3:could not allocate memory for WORK'
            imem1=.FALSE.
          endif
          allocate(IWORK(LIWORK),STAT=ierr)
          if(ierr.ne.0) then
            imem2=.FALSE.
            write(6,*)'DIAGGE3:could not allocate memory for IWORK'
          endif
          if(imem1.and.imem2)then
            write(6,*)'DIAGGE3:Calling DSPGVD'
            CALL DSPGVD( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
     &                HSTOR(initial,1), Eig, A, LDA, &
     &                   WORK, LWORK, IWORK, LIWORK, INFO )
          else
            write(6,*)'DIAGGE3:not enough memory for WORK'
            write(6,*)'DIAGGE3:will not execute DSPGVD'
            write(6,*)'SORRY,BYE'
            CALL STOPIT
          endif
        CASE(2)
            write(6,*)'DIAGGE3:Calling DSPGV'
            LWORK=3*N
            allocate(WORK(LWORK),STAT=ierr)
            if(ierr.ne.0)then
              write(6,*)'diagge3:no space for WORK for DSPGV'
              write(6,*)'DIAGGE3:will not execute DSPGV'
              write(6,*)'SORRY,BYE'
              CALL STOPIT
            endif
            CALL DSPGV( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
                       HSTOR(initial,1), Eig, A, LDA, WORK, INFO )
!    call DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, eig, WORK, &
!    call DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, eig, WORK, &
!    &                   LWORK,  INFO )
        CASE DEFAULT
          write(6,*)'DIAGGE3:Running default case:DSPGVX'
          allocate(work(8*N),stat=ierr)
          if(ierr.ne.0)then
            write(6,*)'diagge3:no space for WORK for DSPGVX'
            write(6,*)'DIAGGE3:will not execute DSPGVX'
            write(6,*)'SORRY,BYE'
            CALL STOPIT
          endif
            LWORK=3*N
            allocate(WORK(LWORK),STAT=ierr)
            allocate(iwork(5*N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge3:no space for IWORK DSPGVX'
            allocate(ifail(N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge3:no space for IFAIL DSPGVX'
            range='I'
            vl=0.0
            vu=0.0
            il=1
            call ELE_INFO(ngrp,nspn,ele_up,ele_dn)
            IU= max(ele_up,ele_dn)+50
            abstol=1.0e-10
            write(6,*)'DIAGGE3:Calling DSPGVX'
            CALL DSPGVX(ITYPE, JOBZ, RANGE, UPLO, N, HSTOR(initial,2), &
               HSTOR(initial,1), VL, VU, IL, IU, ABSTOL, M, Eig, A, &
               LDA, WORK, IWORK, IFAIL, INFO)

        END SELECT
      if (INFO .EQ. 0) then
!          write(6,*) 'Diagonalization is succesful'  ! Remove later
      else
           write(6,*) 'Diagonalization is not succesful'
           write(6,*) 'Err: INFO=',INFO
      endif

      deallocate(WORK); deallocate(IWORK);
      return

      END SUBROUTINE DIAGGE3
!-----------------------------------------------------------------------------------------
