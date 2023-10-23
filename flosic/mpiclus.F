C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
c
c main routine for MPI
c
c 
c 09/11/97 David C. Patton
c 11/06/97 put in statistics (DCP)
c numerous fixes and extensions by DVP and MRP
c 09/08/11 Luis Basurto
c implemented second level parallelism
c for group calculations
c
c MPI Implementation of cluster code
c
c (note: specific MPI names are in caps)
c
#ifdef MPI
!#ifdef GROUP
!      SUBROUTINE MPICLUS
!#else
      PROGRAM MPICLUS
!#endif
      use debug1
      use common3,only : RMAT
      use mpidat1,only : INUSE,NPROC,NCALLED,IRANK,MYCOMM,MYGROUP,IOUT
!     &                   NPROC_GLOBAL,IRANK_GLOBAL
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      INCLUDE 'mpif.h'
      SAVE
      INTEGER IERR,ERRCODE,LENGTH,JOB
      INTEGER I,TAG
      INTEGER IRECVSTAT(MPI_STATUS_SIZE)
      REAL*8  TIME1,TIME2
      CHARACTER*(MPI_MAX_PROCESSOR_NAME) NAME
      REAL*8,ALLOCATABLE :: WCTIME(:)
      CHARACTER*20 group_txt
      CHARACTER*20 dir_name
      iout=6
      group_txt='0'

c
c setup MPI environment
c note: nproc is the number of processes-1 
c 
!#ifdef GROUP
!      write(group_txt,'(1I0)') irank
!      group_txt=adjustl(group_txt)
!      dir_name='print-g'//trim(group_txt)//'.log'
!      MYCOMM=MPI_COMM_WORLD 
!#else
      CALL MPI_INIT(IERR)
      IF (IERR.NE.0) THEN
        write(iout,*)'Error starting MPI. Terminating'
        ERRCODE=0
        CALL MPI_ABORT(MPI_COMM_WORLD,ERRCODE,IERR)
      END IF
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
!      IRANK_GLOBAL=IRANK
!#endif
      if(irank.eq.0) then
        iout=6
C        open(unit=iout,file=dir_name)
C        write(iout,*)'MPI: Number of processes=',NPROC
      else
C        iout=6+irank
C        open(unit=iout,file=dir_name)
C        iout=6
      endif
!#ifdef GROUP
!      write(iout,*)'MPI: Number of processes=',NPROC,' My rank=',
!     &                            IRANK,' My group=',MYGROUP
!#endif
      NPROC=NPROC-1
!      NPROC_GLOBAL=NPROC
      ALLOCATE(INUSE(NPROC),STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'MPICLUS: ERROR ALLOCATING INUSE'
      ALLOCATE(WCTIME(NPROC),STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'MPICLUS: ERROR ALLOCATING WCTIME'
      call flush(iout)
      CALL MPI_GET_PROCESSOR_NAME(NAME,LENGTH,IERR)
C      write(iout,*)'MPI: I am processor ',NAME(1:LENGTH)
      NCALLED=0
      DO I=1,NPROC
       INUSE(I)=0
      END DO
      call flush(iout)
c
c startup cluster code on one node and cpuhog on the rest
c
      IF (IRANK.EQ.0) THEN
!#ifdef GROUP
!       write(iout,*)'THIS IS MANAGER NODE FROM GROUP',MYGROUP
!#else
       write(iout,*)'MPI: Number of worker processes=',NPROC
!#endif
       TIME1=MPI_WTIME()
       CALL CLUSTER2()
       JOB = -1
       TAG = 0
       DO I = 1,NPROC
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,I,TAG,MPI_COMM_WORLD,IERR)
       END DO
       TIME2 = MPI_WTIME()
       write(iout,*)' '
       write(iout,*)'MPI: Total Wall clock TIME = ',TIME2-TIME1
       write(iout,*)' '
c
c get some statistics and print them 
c
       DO I=1,NPROC
        TAG= 12
        CALL MPI_RECV(WCTIME(I),1,MPI_DOUBLE_PRECISION,I,
     &                TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       END DO 
       write(iout,*)'***********************'
       write(iout,*)'MPI Node       CPU time'
       write(iout,*)'-----------------------'
       call flush(iout)
       DO I=1,NPROC
        PRINT 1000,I,WCTIME(I)
1000    FORMAT(1X,I5,5X,F12.2)
       END DO
       write(iout,*)'************************************************'
      ELSE
       CALL CPUHOG()
      END IF
c
        call flush(iout)
c shutdown MPI
c
      DEALLOCATE(INUSE,STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'MPICLUS: ERROR DEALLOCATING INUSE'
      DEALLOCATE(WCTIME,STAT=IERR)
      IF(IERR.NE.0) WRITE(6,*)'MPICLUS: ERROR DEALLOCATING WCTIME'
!#ifdef GROUP
!      RETURN
!#else
      CALL MPI_FINALIZE(IERR)
!#endif

#else
      SUBROUTINE MPICLUS
#endif
      END
