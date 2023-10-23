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
!       Eig     (O) :  Eigenvalues 
!       IEV     (I) :  should be 0 if eigenvectors are not required
! 
! *************************************************************************
!      SUBROUTINE DIAG_S_SYM
      SUBROUTINE DIAG_S_SYM(initial,IAEVAL)
#ifdef SCALAPACK
!      SUBROUTINE DIAGGES
       use mpidat1,only : NPROC,IRANK
       use global_inputs,only : inbas,iiev,iimesh,mpi_io1,idiag3
!       use hstor1,only : maxclustersize,NNZO,OVERS,JAO,IAO
       use hstor1,only : hstor, maxclustersize
       use for_diag1
       use debug1
!       use common2,only : ISPN
!      use general
!      use mpi
      implicit none
      include 'mpif.h'
!      integer, intent(IN)    :: N
!      integer, intent(inout) :: IEV
!       real*8,  intent(in)  :: HSTOR(:,:)
!      real*8,  intent(out)   :: A(:,:),EIG(:)
      integer, intent(INOUT) :: initial
      real*8, intent(OUT)    :: IAEVAL(inbas)
       
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
      integer :: IL,IU,LDZ,M,i,j,ig,jg
      integer :: kindex,auxi,auxj
      integer :: IA,JA,IB,JB,IZ,JZ,NZ,NB,LDR,LDC
      integer :: NPROW,NPCOL,ICONTXT,MYROW,MYCOL
      integer :: LDR2,LDC2
      integer :: IAM,NODE,NPROC2
      integer :: NUMROC
      real*8  :: ABSTOL,VL,VU
      real*8  :: ORFAC
      integer, allocatable :: DESCA(:), DESCZ(:)
      real*8,  allocatable :: W(:), Z(:,:),AA(:,:)
      integer :: dest_prow,dest_pcol,local_row,local_col
      integer :: block_l,block_m,ik,jk
      integer :: ierr,istatus(MPI_STATUS_SIZE)
      integer :: nbrows,nbcolumns,col_max,row_max,ibrows,ibcolumns
      integer :: rows,columns
      integer,dimension(2) :: pdims,dims,distribs,dargs
      integer :: darray
      integer :: locsize, nelements
      integer(kind=MPI_ADDRESS_KIND) :: lb, locextent
      integer(kind=MPI_OFFSET_KIND) :: disp

      integer :: infile,up_limit
      real*8  :: PDLAMCH
      integer :: DSCALE

      integer :: row
!      integer :: clock
!      real*8  :: comptime
! Solving:  A*X = E*B*X  for itype =1
      if(initial==0)initial=1
      N=INBAS
      ITYPE=1
!      IF(IRANK.eq.0)THEN
!        WRITE(6,*)'BEFORE BLACS_PINFO'
!        CALL FLUSH(6)
!      ENDIF
! BLACS_PINFO is the same as MPI_COMM_SIZE
      CALL BLACS_PINFO(IAM,NPROC2)
      if(iam.eq.0)then
        up_limit=MIN(N,MVPS)
!         write(6,*)'MPI: Number of processes=',NPROC2
        if(NPROC+1/=NPROC2)write(6,*)'Number of processors problem'
        call flush(6)
      endif
! Generate  processor grid from MPI_Dims_create
      pdims = 0
      call MPI_Dims_create(nproc2, 2, pdims, ierr)
      nprow = pdims(1)
      npcol = pdims(2)
!      if(iam==0) write(6,*)'prow=',nprow,'pcol=',npcol
      IIEV=0
      if (IIEV   .eq. 0) JOBZ='N' ! Only eigenvalues are needed.
      if (IIEV   .eq. 1) JOBZ='V' ! Both eigenvectors and eigenvalues are to be calculated.
      JOBZ='V'
!      if (IDIAG3 .eq. 0) JOBZ='V'! Both eigenvectors and eigenvalues are to be calculated.
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
      if(MYROW.eq.-1) go to 20
!         if(iam==0)then
!           write(6,*)'icontxt3',ICONTXT,NPROW,NPCOL,MYROW,MYCOL
!           call flush(6)
!         endif
!         if(MYROW.ne.-1)then
!         write(6,*)'Myrow',MYROW,'Mycol',MYCOL
!         end if
      if(IIEV.eq.0)then
        RANGE='A'
      else
        RANGE='I'
      endif
      IA=1
      JA=1
      IB=IA
      JB=JA
      IZ=IA
      JZ=JA
      VL=0.0
      VU=0.0
      IL=1
      IF(IIMESH)THEN
        IU=N
      ELSE
        IF(IAM.EQ.0)THEN
          call ELE_INFO(ngrp,nspn,ele_up,ele_dn)
          IU= INT(max(ele_up,ele_dn))+100
        ENDIF
        CALL MPI_BCAST(IU,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
      ENDIF
!      IU=N
!         ABSTOL=1.0e-10
!         ABSTOL=2*PDLAMCH('S')
      ABSTOL=PDLAMCH(ICONTXT,'U')
!         ABSTOL=-1.0
      ORFAC=-1
!      NBU=INT(N/(5*NPROC2))
      NBU=64
! Determine Global matrix blocking factor
      call BLOCKSET(NB, NBU, N, NPROW, NPCOL)
!         NB=4
!      if(iam==0) then
!        write(6,*)'IU=',IU
!        write(6,*)'IIEV=',IIEV
!        write(6,*)'NB=',NB
!      endif
! Calculate size of local matrix
      LDR=NUMROC(N,NB,MYROW,0,NPROW)
      LDC=NUMROC(N,NB,MYCOL,0,NPCOL)
!      WRITE(IAM+7,*)LDR,LDC
!      WRITE(IAM+7,*)N,NB,MYROW,MYCOL,NPROW,NPCOL
!      CALL FLUSH(IAM+7)

! Initialize array descriptors
      ALLOCATE(DESCA(10),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating DESCA'
      ALLOCATE(DESCZ(10),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating DESCZ'

! Fill array descriptors
      CALL DESCINIT(DESCA, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL DESCINIT(DESCZ, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)

! Allocate local arrays
      allocate(Z(LDR, LDC),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating Z'
      allocate(AA(LDR, LDC),STAT=IERR)
      IF(IERR/=0)WRITE(6,*)'sdiagge_n:error allocating AA'
      AA(:,:)=0.0D0

! Get data into local arrays by either MPI read or Previous Broadcast of HSTOR
!      WRITE(IAM+7,*)'Made it this far'
!      CALL FLUSH(IAM+7)
! Start MPI read
      if(mpi_io1)then
! Define MPI data type to read from file
      else
!        WRITE(IAM+7,*)'Beggining distribution'
!        CALL FLUSH(IAM+7)
!        if(iam==0) call tick(clock)

! Fill local HAM
!        IF(iam==0) WRITE(6,*) 'sdiagge: fill local HAM, ISPN:',ISPN
    
! Fill local OVER
        auxi=0
        do ig=1,N
          auxi=auxi+ig
          do jg=ig,N

            dest_prow=mod(floor(real((ig-1)/NB)),NPROW)
            dest_pcol=mod(floor(real((jg-1)/NB)),NPCOL)
            if(myrow.eq.dest_prow.and.mycol.eq.dest_pcol)then
              block_l=floor(real((ig-1)/(NPROW*NB)))
              block_m=floor(real((jg-1)/(NPCOL*NB)))
              local_row=mod(ig-1,NB)+1
              local_col=mod(jg-1,NB)+1
              ik=block_l*NB+local_row
              jk=block_m*NB+local_col
              if(jg>=ig)then
                kindex=(ig-1)*N+ig+jg-auxi
              else
                write(6,*)'ERROR: Problem in distribution',jg,ig
!                kindex=(jg-1)*N+ig+jg-auxj
                stop
              endif
! First dimension of HSTOR is Overlap, second is Hamiltonian
               AA(ik,jk)=hstor(kindex,1)
!              AA(ik,jk)=OVERS(i)
            endif
          end do
        enddo
! End traverse distribution
      endif
! Allocate driver arrays
      allocate(WORK(3),stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error allocating WORK(3)',iam
      allocate(IWORK(3),stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error allocating IWORK(3)',iam
      LIWORK=-1
      LWORK=-1
! These arrays are not used by divide and conquer
!      if(idiag3.ne.0)then
!      endif
! Do transformation for other solvers (standard solver and divide and conquer)
!      if(idiag3.ne.1)then
!      endif
! Select diagonalization to use
          if(iam==0) write(6,*)'Solving with Divide and Conquer'
          call pdsyevd(JOBZ,UPLO,N,AA,1,1,DESCA,IAEVAL,Z,1,1,DESCZ,WORK,&
                     LWORK,IWORK,LIWORK,INFO)
          lwork=INT(WORK(1))+(maxclustersize-1)*N
          liwork=7*N+8*NPCOL+2         
          if(iam==0)write(6,*)'Second job query INFO=',INFO
!          if(iam==0)then
!            write(6,*)'LWORK2=',LWORK,'LIWORK=',LIWORK
!           write(6,*)'LWORK2=',real(LWORK*8)/1048576.0,'Mb',' LIWORK=',real(LIWORK*8/1024.0),'Kb'
!          endif
          deallocate(work,stat=ierr)
          deallocate(iwork,stat=ierr)
          allocate(work(lwork),stat=ierr)
          allocate(iwork(liwork),stat=ierr)
          call pdsyevd(JOBZ,UPLO,N,AA,1,1,DESCA,IAEVAL,Z,1,1,DESCZ,WORK,&
                       LWORK,IWORK,LIWORK,INFO)
          if(iam==0)write(6,*)'Divide and Conquer Solver INFO=',INFO
! back transformation
! write eigenvalues

! collect all the orthonormal eigenvectors of the matrix
! collection is done only when eigenvectors were requested (IIEV=1)
       if(IIEV.eq.1.and.info.eq.0)then
!Begin manager part
!Begin worker part.
       else
!send the processor coordinates in grid
!send number rows and columns
!send the local matrix
!        CALL MPI_Send(Z,LDR*LDC , MPI_DOUBLE_PRECISION, 0, 30,MPI_COMM_WORLD, IERR)
! End worker part.
      endif!if(IIEV.eq.1 ...)

! Deallocate local arrays
      deallocate(DESCZ,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating DESCZ',iam
      deallocate(DESCA,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating DESCA',iam
      deallocate(Z,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating Z',iam
      deallocate(AA,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating AA',iam
      deallocate(IWORK,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating IWORK',iam
      deallocate(WORK,stat=ierr)
      if(ierr/=0)write(6,*)'sdiagge_n:error deallocating WORK',iam
      if(idiag3.ne.0)then
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
      END SUBROUTINE DIAG_S_SYM


!-----------------------------------------------------------------------------------------
!      subroutine blockset( nb, nbuser, n, nprow, npcol)
!
!      end subroutine blockset
!---------------------------------------------------------------------------------------
