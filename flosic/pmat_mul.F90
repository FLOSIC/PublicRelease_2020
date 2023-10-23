! UTEP Electronic Structure Lab (2020)
!   The subroutine interface for LAPACK routines. The interface is created 
! following old routine style for easy comptability and understanding of 
! routines.
!    Raja, El Paso, TX, Wed Jul 15 13:36:29 CDT 2009
!
!   Same MPI setup from diagges. Scatters and gathers matrices to block-cyclic 
!  using pdgemr2d. performs matrix-multiplication using pdgemm.  
!    cmd, May 2016
!
! *************************************************************************
! 
!     SUBROUTINE PMAT_MUL
!     =================
!     
! *************************************************************************
!
!  PMAT_MUL creates an interface to pdgemm matrix multiplication
!  C=alpha*A*B+beta*C
!
!     Parameters: 
! 
!       N       (I) :  Dimension of Problem  
!       A,B     (O) :  Input Matrices 
!       C       (O) :  Output matrix
! 
! *************************************************************************
    SUBROUTINE PMAT_MUL(N,A,B,C)
#ifdef SCLAPACK
!      SUBROUTINE PMAT_MUL
       use mpidat1,only : NPROC,IRANK
       use debug1
!      use mpi
      use global_inputs, only : inbas
      implicit none
      include 'mpif.h'
      integer, intent(IN)    :: N
      real*8,  intent(in)  :: A(N,N),B(N,N)
      real*8,  intent(out)   :: C(N,N)
       
      ! local variables
!      integer :: N
      integer :: info,ierr
      integer :: i,j,ig,jg
      integer :: NBU,NB,LDR,LDC,NUMROC
      integer :: nprow,npcol,ICONTXT,myrow,mycol
      integer :: IAM,NPROC2
      integer, allocatable :: DESCA(:), DESCB(:), DESCZ(:), descfull(:)
      real*8,  allocatable :: Z(:,:),AA(:,:),BB(:,:)
      integer,dimension(2) :: pdims
   
      character*1  :: ORDER
      character    :: transa,transb
      real*8 :: alpha, beta

!for scattering/gathering local matrices
!      integer :: IL,IU,LDZ,M
!      integer :: IA,JA,IB,JB,IZ,JZ,NZ
!      integer :: dest_prow,dest_pcol,local_row,local_col
!      integer :: block_l,block_m,ik,jk
      integer :: nbrows,nbcolumns,col_max,row_max,ibrows,ibcolumns
      integer :: NODE
      integer :: LDR2,LDC2
      integer :: rows,columns
      integer :: istatus(MPI_STATUS_SIZE)

      integer :: sendr,sendc,recvr,recvc,row,col,nr,nc
!      integer, external :: indxg2l, indxg2p

      INFO =99

! BLACS_PINFO is the same as MPI_COMM_SIZE
      CALL BLACS_PINFO(IAM,NPROC2)
!      if(iam.eq.0) write(6,*)'in pmat_mul',N
    
!Bring everyone else in
      if(iam==0)then
       call senddata(115)
      end if

!broadcast matrix size for MPI setup
!      write(iam+7,*)'in pmat_mul. matrix size',N
!      call flush(iam+7)
      CALL MPI_Bcast(N,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!      write(iam+7,*)'in pmat_mul. matrix size',N
!      call flush(iam+7)

! Generate  processor grid from MPI_Dims_create
      pdims = 0
      CALL MPI_Dims_create(nproc2, 2, pdims, ierr)
      nprow = pdims(1)
      npcol = pdims(2)
!      if(iam==0) write(6,*)'prow=',nprow,'pcol=',npcol

! Generate context
      CALL BLACS_GET(-1, 0, ICONTXT)

! Initialize context
      ORDER='R'
      CALL BLACS_GRIDINIT(ICONTXT, ORDER, nprow, npcol)

! Get context info
      CALL BLACS_GRIDINFO(ICONTXT, nprow, npcol, myrow, mycol)
!         write(iam+7,*)"myrow=",myrow,"mycol=",mycol

      if(myrow.eq.-1) go to 20

!       if(iam==0)then
!         write(6,*)'icontxt3',ICONTXT,nprow,npcol,myrow,mycol
!         CALL flush(6)
!       endif

! Determine Global matrix blocking factor
!      NBU=INT(N/(5*NPROC2))
      NBU=64
      CALL BLOCKSET(NB, NBU, N, nprow, npcol)
!      if(iam==0) write(6,*)'NB=',NB

! Calculate size of local matrix
      LDR=NUMROC(N,NB,myrow,0,nprow)
      LDC=NUMROC(N,NB,mycol,0,npcol)
!      write(IAM+7,*)LDR,LDC
!      CALL FLUSH(IAM+7)

! Initialize array descriptors
      allocate(DESCA(10),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating DESCA'
      allocate(DESCB(10),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating DESCB'
      allocate(DESCZ(10),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating DESCZ'

! Fill array descriptors
      CALL DESCINIT(DESCA, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL DESCINIT(DESCB, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)
      CALL DESCINIT(DESCZ, N, N, NB, NB, 0, 0, ICONTXT, LDR, INFO)

! Allocate local arrays
      allocate(AA(LDR, LDC),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating Z'
      allocate(BB(LDR, LDC),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating Z'
      allocate(Z(LDR, LDC),STAT=ierr)
       if(ierr/=0)write(6,*)'pmat_mul:error allocating Z'


       AA(:,:)=0.0D0
       BB(:,:)=0.0D0
       Z (:,:)=0.0D0
!       C (:,:)=0.0D0

! Send A and B matrices to block distributed matrices
      sendr=0
      sendc=0
      recvr=1
      recvc=1
      do row=1,N ,NB
       sendc=0
       ! Number of rows to be sent
       ! Is this the last row block?
       nr=NB
       if(N-row+1<NB) nr=N-row+1
       do col=1,N,NB
        ! Number of cols to be sent
        ! Is this the last col block?
        nc=NB
        if(N-col+1<NB) nc=N-col+1
        if(iam==0)then
         ! Send  a nr-by-nc submatrix to process (sendr,sendc)
         call dgesd2d(ICONTXT,nr,nc,A(row,col),N,sendr,sendc)
         call dgesd2d(ICONTXT,nr,nc,B(row,col),N,sendr,sendc)
        end if
        if(myrow==sendr .and. mycol==sendc)then
         ! Receive the same data
         ! The leading dimension of the local matrix is LDR(nrows)!
         call dgerv2d(ICONTXT,nr,nc,AA(recvr,recvc),LDR,0,0)
         call dgerv2d(ICONTXT,nr,nc,BB(recvr,recvc),LDR,0,0)
         recvc=mod(recvc+nc,LDC)
        end if
        sendc=mod(sendc+1,npcol)
       end do
       if(myrow==sendr) recvr=mod(recvr+nr,LDR)
       sendr=mod(sendr+1,nprow)
      end do
       
      if(.false.)then
! pdgemr2d does not work if N*N > max 32bit int (2^31-1)
       allocate(descfull(10),stat=ierr)
        if(ierr/=0)write(6,*)'pmat_mul:error allocating descfull'
       CALL descinit(descfull, N, N, N, N, 0, 0, ICONTXT, N , INFO)
       CALL pdgemr2d(N,N,A,1,1,descfull,AA,1,1,DESCA,ICONTXT)
       CALL pdgemr2d(N,N,B,1,1,descfull,BB,1,1,DESCB,ICONTXT)
      end if

!*****************************************************************************************

! setup matrix multiplication Z=alpha*A*B+beta*Z
       alpha=1.0
       beta=0.0
       transa='N'
       transb='N'

!       if(iam==0)write(6,*)'pre pdgemm',iam,N
!        CALL flush(6)

       CALL PDGEMM(transa,transb,N,N,N,alpha, &
                   AA,1,1,DESCA, &
                   BB,1,1,DESCB, beta, &
                    Z,1,1,DESCZ)
!       if(iam==0)write(6,*)'post pdgemm',iam
!        CALL flush(6)

!gather matrix solution back to C 
      if(.false.)then
       CALL pdgemr2d(N,N,Z,1,1,DESCZ,C,1,1,descfull,ICONTXT)
       deallocate(descfull)
      else
!        WRITE(IAM+7,*)'Beggining distribution to C'
!        CALL FLUSH(IAM+7)

!       if(IIEV.eq.1.and.info.eq.0)then
!       if(info.eq.0)then
!Begin manager part
        IF (iam.EQ.0) THEN
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
!                     AHAM(ig,jg)= Z((ibrows-1)*NB + rows, (ibcolumns-1)*NB +columns)
                     C(ig,jg)= Z((ibrows-1)*NB + rows, (ibcolumns-1)*NB +columns)
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
!       end if   !if(info==0)
      end if   !ifpdgemr2d


!       if(iam==0)write(6,*)'end of pmat_mul',iam
!        CALL flush(6)
!------------------------------------------------------------------------------------------

! Deallocate local arrays
!      if(allocated(A_S))deallocate(A_S,B_S,JAA,IAA,JAB,IAB,stat=ierr)
!      if(ierr/=0)write(6,*)'pmat_mul:error deallocating sparse arrays',iam
      deallocate(DESCZ,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating DESCZ',iam
      deallocate(DESCB,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating DESCB',iam
      deallocate(DESCA,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating DESCA',iam
      deallocate(Z,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating Z',iam
      deallocate(BB,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating BB',iam
      deallocate(AA,stat=ierr)
      if(ierr/=0)write(6,*)'pmat_mul:error deallocating AA',iam

!      if(INFO.ne.0)then
!        if (IAM.EQ.0) then
!           write(6,*) 'Multiplication is not succesful'
!        endif 
!      endif
      CALL BLACS_GRIDEXIT(ICONTXT)
      return
20    write(6,*)'PROBLEM: A PROCESSOR IS NOT IN SCALAPACK GRID'
      stop
#endif
      END SUBROUTINE PMAT_MUL


!-----------------------------------------------------------------------------------------
!      subroutine blockset( nb, nbuser, n, nprow, npcol)
!
!     This subroutine try to choose an optimal block size
!     for the distributd matrix.
!     Written by Carlo Cavazzoni, CINECA
!
!      implicit none
!      integer :: n,nb, nprow, npcol, nbuser
 
!      nb = min ( n/nprow, n/npcol )
!      if(nbuser.gt.0) then
!        nb = min ( nb, nbuser )
!      endif
!      nb = max(nb,1)
 
!      return
!      end subroutine blockset
!!---------------------------------------------------------------------------------------
