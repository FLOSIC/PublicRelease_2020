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
!       N       (I) :  Dimension of Problem  
!       A       (O) :  Eigenvector matrix (eigenvectors in columns) 
!       Eig    (O) :  Eigenvalues 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
! *************************************************************************
      SUBROUTINE DIAGGE2(N,initial,A,Eig,IEV)
       use hstor1,only : hstor
       implicit none
       integer, intent(IN)     ::  N,IEV
       integer, intent(INOUT)  ::  initial
!       real*8, intent(INOUT)   ::  A(NA,NA),B(NA,NA),AUX(NA)
       real*8, intent(OUT)     ::  A(N,N),Eig(N)

       ! local variables
       CHARACTER            :: JOBZ, UPLO,RANGE
       integer              :: info,lda,liwork,lwork,itype,ierr,il,iu,m,&
                               ngrp,nspn,ele_up,ele_dn
       real*8               :: vl,vu,abstol
       integer, allocatable :: IWORK(:),ifail(:)
       real*8, allocatable  :: WORK(:)
       logical              :: imem1,imem2
! Solving:  A*C = E*B*C  for itype =1


        itype=1; 
        imem1=.TRUE.
        imem2=.TRUE.
        if (iev .eq. 0) jobz='N'  ! Only eigenvalues are needed.
        if (iev .eq. 1) jobz='V'  ! Both eigenvectors and eigenvalues are to be calculated.
        ! Remove folllowing block later
!       if (NA   .ne. N) then
!           write(6,*) 'Note Dimension of matrix and that of problem are different!' 
!       endif

        UPLO='L'      ! Lower triangular matrix
        LDA=N        
        ! Following LWORK, and LIWORK is to compute size of temperary arrays
        allocate(WORK(1));
        allocate(IWORK(1));
        LWORk=-1
        LIWORK=-1
        if(initial.eq.0) initial=1
        CALL DSPGVD( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
     &                HSTOR(initial,1), Eig, A, LDA, &
     &                   WORK, LWORK, IWORK, LIWORK, INFO )
        

        LWORK= INT(WORK(1))
        LIWORK = IWORK(1)
!        write(6,*)'diagge2: ', Lwork, LIWORK, N
        deallocate(work)
        deallocate(iwork)
        allocate(WORK(LWORK),STAT=ierr)
        if(ierr.ne.0) then
           write(6,*)'DIAGGE2:could not allocate memory for WORK'
           imem1=.FALSE.
        endif
        allocate(IWORK(LIWORK),STAT=ierr)
        if(ierr.ne.0) then
           imem2=.FALSE.
           write(6,*)'DIAGGE2:could not allocate memory for IWORK'
           itype=2
        endif
        INFO =999
        if(imem1.and.imem2)then
          write(6,*)'DIAGGE2:Calling DSPGVD'
          CALL DSPGVD( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
     &                HSTOR(initial,1), Eig, A, LDA, &
     &                   WORK, LWORK, IWORK, LIWORK, INFO )
!    call DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, eig, WORK, &
!    call DSYGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, eig, WORK, &
!    &                   LWORK,  INFO )
        else
          write(6,*)'DIAGGE2:not enough memory for WORK'
          write(6,*)'DIAGGE2:will not execute DSPGV'
          allocate(work(8*N),stat=ierr)
          if(ierr.ne.0)then
            write(6,*)'diagge2:no space for WORK for DSPGVX'
            write(6,*)'DIAGGE2:will not execute DSPGVX'
            LWORK=3*N
            allocate(WORK(LWORK),STAT=ierr)
            if(ierr.ne.0)then
              write(6,*)'DIAGGE2:not space for WORK for DSPGV'
              write(6,*)'Not enough memory for anything, forget it!'
              call stopit
            else
              write(6,*)'DIAGGE2:Calling DSPGV'
              CALL DSPGV( ITYPE, JOBZ, UPLO, N, HSTOR(initial,2), &
                       HSTOR(initial,1), Eig, A, LDA, WORK, INFO )
            endif
          else 
            allocate(iwork(5*N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge2:no space for IWORK DSPGVX'
            allocate(ifail(N),stat=ierr)
            if(ierr.ne.0)write(6,*)'diagge2:no space for IFAIL DSPGVX'
            range='I'
            vl=0.0
            vu=0.0
            il=1
            call ELE_INFO(ngrp,nspn,ele_up,ele_dn)
            IU= max(ele_up,ele_dn)+20
            abstol=0.0
            write(6,*)'DIAGGE2:Calling DSPGVX'
            CALL DSPGVX(ITYPE, JOBZ, RANGE, UPLO, N, HSTOR(initial,2), &
               HSTOR(initial,1), VL, VU, IL, IU, ABSTOL, M, Eig, A, &
               LDA, WORK, IWORK, IFAIL, INFO)
          endif
        endif
      if (INFO .EQ. 0) then
!          write(6,*) 'Diagonalization is succesful'  ! Remove later
      else
           write(6,*) 'Diagonalization is not succesful'
           write(6,*) 'DSPGVD Err: INFO=',INFO
      endif

       deallocate(WORK); deallocate(IWORK);
        return

      END SUBROUTINE DIAGGE2
!-----------------------------------------------------------------------------------------
