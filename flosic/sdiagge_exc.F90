! UTEP Electronic Structure Lab (2020)
!   The subroutine interface for LAPACK routines. The interface is created 
! following old rouitne style for easy comptability and understanding of 
! routines.
!    Raja, El Paso, TX, Wed Jul 15 13:36:29 CDT 2009
!
! *************************************************************************
! 
!     SUBROUTINE DIAGGES_EXC
!     =================
!     
! *************************************************************************
!
!  Diagge creates an interface to the diagonalization routine ewevge.
!
!  modified for diagonalization in excited state code
! diagonalizes full AHAM and AOVER matrices
!     - cmd, may 2016
!
!     Parameters: 
! 
!       N       (I) :  Dimension of Problem  
!       A       (O) :  Eigenvector matrix (eigenvectors in columns) 
!       Eig     (O) :  Eigenvalues 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
! *************************************************************************
!      SUBROUTINE DIAGGES(N,A,EIG,IEV)
    SUBROUTINE DIAGGES_EXC
#ifdef SCALAPACK
       use mpidat1,only : NPROC,IRANK
       use global_inputs,only : inbas,iiev,iimesh,mpi_io1,idiag3
       use hstor1,only : hstor,maxclustersize
       use for_diag1
       use debug1
       use common5,only : PSI_COEF
!      use general
!      use mpi
      implicit none
      include 'mpif.h'
!      integer, intent(IN)    :: N
!      integer, intent(inout) :: IEV
!      real*8,  intent(out)   :: A(:,:),EIG(:)
       
       ! local variables
      integer              :: N
      integer              :: NBU
      CHARACTER            :: JOBZ, UPLO
      CHARACTER*1          :: RANGE,ORDER
      integer              :: info,lda,ldb,ngrp,nspn,ITYPE
      integer              :: LWORK,LIWORK,CLUSTERSIZE
      real*8               :: ele_up,ele_dn
      integer, allocatable :: IWORK(:)
      real*8,  allocatable :: WORK(:)
! Following variable & two arrarys needed for selected calculation of eigenvalues using dsygvx 
      integer :: IL,IU,LDZ,M,i1,i,j,ig,jg,kindex,auxi,auxj
      integer :: IA,JA,IB,JB,IZ,JZ,NZ,NB,LDR,LDC
      integer :: NPROW,NPCOL,ICONTXT,MYROW,MYCOL
      integer :: LDR2,LDC2
      integer :: IAM,NODE,NPROC2
      integer :: NUMROC
      real*8  :: ABSTOL,VL,VU
      real*8  :: ORFAC
      integer, allocatable :: IFAIL(:)
      integer, allocatable :: ICLUSTR(:)
      integer, allocatable :: DESCA(:), DESCB(:), DESCZ(:), descfull(:)
      real*8,  allocatable :: W(:), Z(:,:), GAP(:),AA(:,:),BB(:,:),PAA(:),PBB(:)
      integer :: dest_prow,dest_pcol,local_row,local_col
      integer :: block_l,block_m,ik,jk
      integer :: ierr,istatus(MPI_STATUS_SIZE)
      integer :: nbrows,nbcolumns,col_max,row_max,ibrows,ibcolumns
      integer :: rows,columns
      integer :: oldclustersize
      integer,dimension(2) :: pdims,dims,distribs,dargs
      integer :: darray
      integer :: locsize, nelements
      integer(kind=MPI_ADDRESS_KIND) :: lb, locextent
      integer(kind=MPI_OFFSET_KIND) :: disp

      integer :: infile,up_limit
      real*8  :: PDLAMCH
      integer :: DSCALE

!      integer :: clock
!      real*8  :: comptime
! Solving:  A*X = E*B*X  for itype =1

      N=inbas

!broadcast matrix size for MPI setup
      CALL MPI_BCAST(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  
      ITYPE=1
! BLACS_PINFO is the same as MPI_COMM_SIZE
      CALL BLACS_PINFO(IAM,NPROC2)
      if(iam.eq.0)then
        if(NPROC+1/=NPROC2)write(6,*)'Number of processors problem'
        call flush(6)
      endif
! Generate  processor grid from MPI_Dims_create
      pdims = 0
      call MPI_Dims_create(nproc2, 2, pdims, ierr)
      nprow = pdims(1)
      npcol = pdims(2)
!      if(iam==0) write(6,*)'prow=',nprow,'pcol=',npcol

      JOBZ='V'                   ! Both eigenvectors and eigenvalues are to be calculated.
      if (IDIAG3 .eq. 0) JOBZ='V'! Both eigenvectors and eigenvalues are to be calculated.

      UPLO='U'      ! Upper/Lower triangular matrix
      INFO =99
      ORDER='R'
! Generate context
      CALL BLACS_GET(-1, 0, ICONTXT)
! Initialize context
      CALL BLACS_GRIDINIT(ICONTXT, ORDER, NPROW, NPCOL)
! Get context info
      CALL BLACS_GRIDINFO(ICONTXT, NPROW, NPCOL, MYROW, MYCOL)
!         write(iam+7,*)"myrow=",myrow,"mycol=",mycol
!         write(iam+7,*)"N=",N,"myrow=",myrow,"mycol=",mycol
      if(MYROW.eq.-1) go to 20

      RANGE='A'

      IA=1
      JA=1
      IB=IA
      JB=JA
      IZ=IA
      JZ=JA
      VL=0.0
      VU=0.0
      IL=1

      IU=N

      ABSTOL=PDLAMCH(ICONTXT,'U')

      ORFAC=-1
!      NBU=INT(N/(5*NPROC2))
      NBU=64
! Determine Global matrix blocking factor
      call BLOCKSET(NB, NBU, N, NPROW, NPCOL)
!        write(6,*)'IU=',IU
!        write(6,*)'IIEV=',IIEV
!        write(6,*)'NB=',NB

! Calculate size of local matrix
      LDR=NUMROC(N,NB,MYROW,0,NPROW)
      LDC=NUMROC(N,NB,MYCOL,0,NPCOL)

!      write(iam+7,*)"diagges_exc:N,NBU,LDR,LDC",N,NB,LDR,LDC
!      call flush(iam+7)

! Initialize array descriptors
      ALLOCATE(DESCA(10),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating DESCA'
      ALLOCATE(DESCB(10),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating DESCB'
      ALLOCATE(DESCZ(10),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating DESCZ'
      allocate(descfull(10),stat=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating descfull'
! Fill array descriptors
      CALL DESCINIT(DESCA, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL DESCINIT(DESCB, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL DESCINIT(DESCZ, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL descinit(descfull, N, N, N, N, 0, 0, ICONTXT, N , INFO)
! Allocate local arrays
      allocate(Z(LDR, LDC),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating Z'
      allocate(AA(LDR, LDC),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating Z'
      allocate(BB(LDR, LDC),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating Z'

!      write(iam+7,*)"diagges_exc:post local allocate"
!      call flush(iam+7)

! Get data into local arrays by either MPI read or Previous Broadcast of HSTOR
! or by scattering AHAM, AOVER with pdgemr2d

!      if(N<7000)then
      if(.false.)then
       if(iam==0)write(6,*)'Beginning distribution via pdgemr2d'
       CALL pdgemr2d(N,N,aham, 1,1,descfull,AA,1,1,DESCA,ICONTXT)
       CALL pdgemr2d(N,N,aover,1,1,descfull,BB,1,1,DESCB,ICONTXT)
      else
!        WRITE(IAM+7,*)'Beggining distribution'
!        CALL FLUSH(IAM+7)
!        if(iam==0) call tick(clock)
        auxi=0
        do ig=1,N
          auxi=auxi+ig
          call global_2_local(ig,nprow,NB,dest_prow,ik)
!          auxj=0
          do jg=ig,N
!            auxj=auxj+jg
            call global_2_local(jg,npcol,NB,dest_pcol,jk)
!            dest_prow=mod(floor(real((ig-1)/NB)),NPROW)
!            dest_pcol=mod(floor(real((jg-1)/NB)),NPCOL)
            if(myrow.eq.dest_prow.and.mycol.eq.dest_pcol)then
!              block_l=floor(real((ig-1)/(NPROW*NB)))
!              block_m=floor(real((jg-1)/(NPCOL*NB)))
!              local_row=mod(ig-1,NB)+1
!              local_col=mod(jg-1,NB)+1
!              ik=block_l*NB+local_row
!              jk=block_m*NB+local_col
              if(jg>=ig)then
                kindex=(ig-1)*N+ig+jg-auxi
              else
!                kindex=(jg-1)*N+ig+jg-auxj
                stop
              endif
! First dimension of HSTOR is Overlap, second is Hamiltonian
              AA(ik,jk)=hstor(kindex,2)
              BB(ik,jk)=hstor(kindex,1)
!            write(iam+7,*)'---------------------------'
            endif
          enddo
        enddo
! End traverse distribution
      end if

!*****************************************************************************************


! Allocate driver arrays
      allocate(WORK(3),stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error allocating WORK(3)',iam
      allocate(IWORK(3),stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error allocating IWORK(3)',iam
      LIWORK=-1
      LWORK=-1
! These arrays are not used by divide and conquer
      if(idiag3.ne.0)then
        allocate(IFAIL(N),stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error allocating IFAIL',iam
        allocate(ICLUSTR(2*NPROW*NPCOL),stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error allocating ICLUSTR',iam
        allocate(GAP(NPROC2),stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error allocating GAP',iam
      endif

!      write(iam+7,*)"diagges_exc:post work/ifail allocate"
!      call flush(iam+7)

! Do transformation for other solvers (standard solver and divide and conquer)
      if(idiag3.ne.1)then
        INFO=99
        call pdpotrf(UPLO,N,BB,1,1,DESCB,INFO)
        if(iam==0.and.info/=0)write(6,*)'From Cholesky factorization INFO=',INFO
        call pdsyngst(1,UPLO,N,AA,1,1,DESCA,BB,1,1,DESCB,dscale,work,&
                      lwork,info)
        if(iam==0.and.INFO/=0)write(6,*)'From job query INFO=',INFO
        lwork=INT(work(1))
        deallocate(work,stat=ierr)
        allocate(work(lwork),stat=ierr)
        call pdsyngst(1,UPLO,N,AA,1,1,DESCA,BB,1,1,DESCB,dscale,work,&
                      lwork,info)
        if(iam==0.and.INFO/=0)write(6,*)'From transformation INFO=',INFO
        deallocate(work,stat=ierr)
        lwork=-1
        liwork=-1
        allocate(work(3),stat=ierr)
      endif

! Select diagonalization to use
      select case(idiag3)
        case (0)
          if(iam==0) write(6,*)'Solving with Divide and Conquer'
          call pdsyevd(JOBZ,UPLO,N,AA,1,1,DESCA,AEVAL,Z,1,1,DESCZ,WORK,&
                     LWORK,IWORK,LIWORK,INFO)
          lwork=INT(WORK(1))+(maxclustersize-1)*N
          liwork=7*N+8*NPCOL+2         
!          if(iam==0)write(6,*)'Second job query INFO=',INFO
          deallocate(work,stat=ierr)
          deallocate(iwork,stat=ierr)
          allocate(work(lwork),stat=ierr)
          allocate(iwork(liwork),stat=ierr)
          call pdsyevd(JOBZ,UPLO,N,AA,1,1,DESCA,AEVAL,Z,1,1,DESCZ,WORK,&
                       LWORK,IWORK,LIWORK,INFO)
          if(iam==0.and.INFO/=0)write(6,*)'Divide and Conquer Solver INFO=',INFO

        case (1)
          if(iam==0)write(6,*)'Solving using PDSYGVX'
          if(myrow.eq.0.and.mycol.eq.0)then
            write(6,*)'enter PDSYGVX 1'
            WRITE(6,*)'UPLO=',UPLO
            WRITE(6,*)'ITYPE=',ITYPE
            WRITE(6,*)'JOBZ=',JOBZ
            WRITE(6,*)'IL=',IL
          endif
          call PDSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, AA, IA, JA, DESCA, &
      &                  BB, IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, &
      &                  NZ,AEVAL,ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, &
      &                  IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )

          if(INFO.ne.0)then
            write(6,*)'sdiagge_n:Error in work request INFO=',INFO
            return 
          endif
          LWORK=INT(WORK(1))+(maxclustersize-1)*N
          LIWORK=IWORK(1)
          if(myrow.eq.0.and.mycol.eq.0)then
           write(6,*)'LWORK2=',real(LWORK*8)/1048576.0,'Mb',' LIWORK=',real(LIWORK*8/1024.0),'Kb'
          endif
          deallocate(WORK,stat=ierr)
          if(ierr/=0)write(6,*)'sdiagge_n:error deallocating WORK',iam
          deallocate(IWORK,stat=ierr)
          if(ierr/=0)write(6,*)'sdiagge_n:error deallocating IWORK',iam
          allocate(WORK(LWORK),stat=ierr)
          if(ierr/=0)write(6,*)'sdiagge_n:error allocating WORK',iam
          allocate(IWORK(LIWORK),stat=ierr)
          if(ierr/=0)write(6,*)'sdiagge_n:error allocating IWORK',iam
          if(myrow.eq.0.and.mycol.eq.0)then
            write(6,*)'enter PDSYGVX 2'
          endif
          call PDSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, AA, IA, JA, DESCA, &
      &                  BB, IB, JB, DESCB, VL, VU, IL, IU, ABSTOL, M, &
      &                  NZ,AEVAL,ORFAC, Z, IZ, JZ, DESCZ, WORK, LWORK, &
      &                  IWORK, LIWORK, IFAIL, ICLUSTR, GAP, INFO )


          if(iam==0)then
            write(6,*)'exit PDSYGVX, INFO=',INFO
          endif
      
          if(info==16.and.iam==0) write(6,*)'IFAIL(1)=',IFAIL(1)
          if(info.eq.2.)then
            oldclustersize=maxclustersize
            do i=1,nproc2-1
              clustersize=iclustr(2*i)-iclustr(2*i-1)
              maxclustersize=max(maxclustersize,clustersize)
            enddo
            if(maxclustersize>oldclustersize.and.iam==0) &
              write(6,*)'maxclustersize',maxclustersize
          endif

        case (2)
          if(iam==0)write(6,*)'Solving standard egivenvalue problem'
          call pdsyevx(JOBZ,RANGE,UPLO,N,AA,1,1,DESCA,VL,VU,IL,IU,ABSTOL,&
                       M,NZ,AEVAL,ORFAC,Z,1,1,DESCZ,WORK,LWORK,IWORK,&
                       LIWORK,IFAIL,ICLUSTR,GAP,INFO)
          if(iam==0)write(6,*)'Second job query INFO=',INFO
          LWORK=INT(WORK(1))+(maxclustersize-1)*N
          liwork=iwork(1)
          deallocate(work)
          deallocate(iwork)
          allocate(work(lwork))
          allocate(iwork(liwork))
          call pdsyevx(JOBZ,RANGE,UPLO,n,AA,1,1,DESCA,VL,VU,IL,IU,ABSTOL,&
                       M,NZ,AEVAL,ORFAC,Z,1,1,DESCZ,work,lwork,iwork,&
                       liwork,IFAIL,ICLUSTR,GAP,info)
          if(iam==0)write(6,*)'Standard Solver INFO=',INFO
          if(info.eq.2.)then
            oldclustersize=maxclustersize
            do i=1,nproc2-1
              clustersize=iclustr(2*i)-iclustr(2*i-1)
              maxclustersize=max(maxclustersize,clustersize)
            enddo
            if(maxclustersize>oldclustersize.and.iam==0) &
              write(6,*)'maxclustersize',maxclustersize
          endif

      end select

! back transformation
      if(idiag3.ne.1)then
        if(JOBZ=='V') then
          call pdtrsm('L',UPLO,'N','N',n,n,1.0D0,BB,1,1,DESCB,Z,1,1,DESCZ)
        end if
      endif
!_________________________________________________________________________________________
!*****************************************************************************************

! collect all the orthonormal eigenvectors of the matrix
! store in AHAM
!       CALL pdgemr2d(N,N,Z,1,1,DESCZ,C,1,1,descfull,ICONTXT)

!      if(N<7000)then
      if(.false.)then

       CALL pdgemr2d(N,N,Z,1,1,DESCZ,aham,1,1,descfull,ICONTXT)

      else

!       if(IIEV.eq.1.and.info.eq.0)then
       if(info.eq.0)then
!Begin manager part
        IF (iam.EQ.0) THEN
!         up_limit=MIN(N,MVPS)
!         write(6,*)'Scalapack:Performing collection of eigenvectors'
!obtain dimensions of the local array on manager node
         LDR2=LDR
         LDC2=LDC
!k is the rank of the processor that sent the data
!nprocs is the number of processors in the communicator
         DO NODE = 0, NPROC2 -1
           IF (node.EQ.0) THEN
             i = myrow
             j = mycol
           ELSE

            !node -- is the rank of the processor that sent the data
            !status -- message status array(of type integer)
             CALL MPI_Recv(i, 1, MPI_Integer, NODE, 10, MPI_COMM_WORLD,ISTATUS, IERR)
             CALL MPI_Recv(j, 1, MPI_Integer, NODE, 20, MPI_COMM_WORLD,ISTATUS, IERR)
         
             CALL MPI_Recv(LDR2, 1, MPI_Integer, NODE, 40, MPI_COMM_WORLD,ISTATUS, IERR)
             CALL MPI_Recv(LDC2, 1, MPI_Integer, NODE, 50, MPI_COMM_WORLD,ISTATUS,IERR )
         
            !create a local matrix with the size of the matrix passed
             DEALLOCATE(Z,stat=ierr)
             if(ierr/=0)write(6,*)'sdiagge_n:error deallocating Z',iam
             ALLOCATE( Z(LDR2,LDC2),stat=ierr)
             if(ierr/=0)write(6,*)'sdiagge_n:error allocating Z',iam
         
             !recieve the local matrix sent from node
             CALL MPI_Recv(Z, LDR2*LDC2, MPI_DOUBLE_PRECISION, NODE,30,MPI_COMM_WORLD, ISTATUS, IERR)
           ENDIF

!compute the number of blocks in each local array
           nbrows = CEILING(REAL(LDR2)/REAL(NB))    !number of blocks in the row
           nbcolumns = CEILING(REAL(LDC2)/REAL(NB)) !number of blocks in the columns

     !loop over each block in the column
!         write(6,*)'working with submatrix from',i,j
            DO ibcolumns = 1, nbcolumns
            !special case - number of columns is less than NBCOL     
              IF(ibcolumns.EQ.nbcolumns) THEN
                col_max = LDC2-(nbcolumns-1)*NB
              ELSE
                col_max = NB !number of columns in a block        
              ENDIF
          
            !loop over each block in the row
              DO ibrows = 1, nbrows
               !special case - number of columns is less than NBCO
                IF (ibrows.EQ.nbrows) THEN
                  row_max = LDR2-(nbrows-1)*NB
                ELSE
                  row_max = NB !number of rows in a block
                ENDIF
               !for each column in the block, loop over each row in the
               !block
                DO columns = 1, col_max
                  DO rows = 1, row_max
                     ig=(i+(ibrows-1)*nprow)*NB + rows
                     jg=(j+(ibcolumns-1)*npcol)*NB + columns
!                     if(jg.LE.up_limit)then
!                       PSI_COEF(ig,jg,iirep,iispn)=Z((ibrows-1)*NB + rows, &
!                               (ibcolumns-1)*NB +columns)
!                     endif
                     AHAM(ig,jg)= Z((ibrows-1)*NB + rows, (ibcolumns-1)*NB +columns)
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
         ENDDO

!End manager part.
! Begin worker part.
        ELSE
!send the processor coordinates in grid
         CALL MPI_Send(myrow, 1, MPI_Integer, 0, 10, MPI_COMM_WORLD,IERR)
         CALL MPI_Send(mycol, 1, MPI_Integer, 0, 20, MPI_COMM_WORLD,IERR)
!send number rows and columns
         CALL MPI_Send(LDR, 1, MPI_Integer, 0, 40, MPI_COMM_WORLD, IERR)
         CALL MPI_Send(LDC, 1, MPI_Integer, 0, 50, MPI_COMM_WORLD, IERR)
!send the local matrix
         CALL MPI_Send(Z,LDR*LDC , MPI_DOUBLE_PRECISION, 0, 30,MPI_COMM_WORLD, IERR)

        ENDIF
! End worker part.
!------------------------------------------------------------------------------------------
       end if
      end if

!------------------------------------------------------------------------------------------

! Deallocate local arrays
      deallocate(DESCZ,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating DESCZ',iam
      deallocate(DESCB,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating DESCB',iam
      deallocate(DESCA,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating DESCA',iam
      deallocate(Z,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating Z',iam
      deallocate(BB,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating BB',iam
      deallocate(AA,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating AA',iam
      deallocate(IWORK,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating IWORK',iam
      deallocate(WORK,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating WORK',iam
      if(idiag3.ne.0)then
        deallocate(GAP,stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error deallocating GAP',iam
        deallocate(ICLUSTR,stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error deallocating ICLUSTR',iam
        deallocate(IFAIL,stat=ierr)
        if(ierr/=0)write(6,*)'sdiagge_n:error deallocating IFAIL',iam
      endif
      IIEV=INFO
      if(INFO.ne.0)then
        if (IAM.EQ.0) then
           write(6,*) 'Diagonalization is not succesful'
        endif 
      endif
      call BLACS_GRIDEXIT(ICONTXT)
      return
20    write(6,*)'PROBLEM: A PROCESSOR IS NOT IN SCALAPACK GRID'
      stop

 100       format(1x,i5,f13.6)
 200       format(1x,i5,i5,f13.6)
#endif
      END SUBROUTINE DIAGGES_EXC


!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
