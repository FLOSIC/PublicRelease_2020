C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE CPUHOG()
c
c 07/23/97 David C. Patton, modified by DVP
c
#ifdef MPI
!      use debug1
      use mpidat1, only : MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:48 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: PRTIME , PRTIME1, PRTIME2
      INCLUDE 'mpif.h'
      SAVE
      INTEGER TAG,JOB,IERR,IRECVSTAT(MPI_STATUS_SIZE) 

c
      PRTIME=0.0D0
c
c check to see what type of task to perform
c
      JOB=0
      DO WHILE (JOB.GE.0)
       TAG=0
       CALL MPI_RECV(JOB,1,MPI_INTEGER,0,
     &               TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       PRTIME1=MPI_WTIME()
!      CALL TRACER('CPUHOG regular:GOT JOB',JOB)
c
c global synchronization barrier
c
       IF (JOB.EQ.100) THEN
        CALL GLBARRIER
       END IF
c
c any kind of general data transfer
c
       IF (JOB.GT.100) THEN
        CALL SENDDATA(JOB)
       END IF
c
c get data and do cpu intensive work in poisson calculation
c     
       IF (JOB.EQ.3) CALL PASPOISS(1)
c 
c send in coulomb and rhog arrays 
c
       IF (JOB.EQ.4) CALL PASPOISS(2)
c
c receive data for gethold execution
c
       IF (JOB.EQ.5) CALL PASHAMIL(1)
c
c send gethold data to manager
c
       IF (JOB.EQ.6) CALL PASHAMIL(2)
c
c deal with parallel calculation of XC and local potentials
c
       IF (JOB.EQ.7) CALL PASVLXC
c
c receive data for parallel forces and call frcslv
c
       IF (JOB.EQ.11) CALL PASFPUL(1)
c
c send force data created by frcslv
c
       IF (JOB.EQ.12) CALL PASFPUL(2) 
c
c receive input data for interstitial mesh and call fillist
c
       IF (JOB.EQ.14) CALL PASMESH(1)
c
c send mesh points to manager
c
       IF (JOB.EQ.15) CALL PASMESH(2)
c
c recieive data  from manager
c      IF (JOB.EQ.16) CALL PASLS(1)
c
c send data  to  manager
c      IF (JOB.EQ.17) CALL PASLS(2)
c
c receive input data for diagge from NEWWAVE
c
C       IF(JOB.EQ.20) CALL PASWAVE(1)
c
c send eigenvalues and vectors back
c
C       IF(JOB.EQ.21) CALL PASWAVE(2)
C
c send data to workers for rhofft  (fourier trans)
c
c      IF(JOB.EQ.22) CALL PASFORM(1)
c
c get the results back
c
c      IF(JOB.EQ.23) CALL PASFORM(2)

C Receive data to do joint density of states section 3

       IF(JOB.EQ.24) CALL PASDOSJNT(3,1)

C Send data back to manager

       IF(JOB.EQ.25) CALL PASDOSJNT(3,2)

C Receive data to do joint density of states section 1

       IF(JOB.EQ.26) CALL PASDOSJNT(1,1)

C Send data back to manager

       IF(JOB.EQ.27) CALL PASDOSJNT(1,2)

c
c receive data for gethold2 execution for libxc/mgga
c
       IF (JOB.EQ.28) CALL PASHAMIL_LIBXC(1)
c
c send gethold2 data to manager
c
       IF (JOB.EQ.29) CALL PASHAMIL_LIBXC(2)
c
c send data to workers for rhofft  (fourier trans)
c
       IF(JOB.EQ.34)THEN
C      PRINT*,'IRANK IS CALLING  PASLMSIC'
            CALL PASLMSIC(1)
       END IF
                           
c
c get the results back
c
       IF(JOB.EQ.35) CALL PASLMSIC(2)
c
c send data to workers for siclag_der calculation
c
c  KAJ 1-23-2023
c
       IF(JOB.EQ.74)THEN
C      PRINT*,'IRANK IS CALLING  PASLMSIC'
            CALL PASLMSIC_DER(1)
       END IF
                           
c
c get the results back
c
       IF(JOB.EQ.75) CALL PASLMSIC_DER(2)
C    
c send data to workers for mesh determination
c
       IF(JOB.EQ.36)THEN
!            CALL PASLM_MSH(1)
       END IF
c
c get the results back
c
!       IF(JOB.EQ.37) CALL PASLM_MSH(2)
C    
#ifdef GROUP
C
C Allocate MESH for fo_mesh and read cutoff value
C
       IF(JOB.EQ.60) CALL GLOBAL_CALL(60)
C
C Read mesh for a given orbital
C
       IF(JOB.EQ.61) CALL GLOBAL_CALL(61)
C
C Initialize FO_MESH to zero
C
       IF(JOB.EQ.62) CALL GLOBAL_CALL(62)
C
C Initialize FO_MESH to zero
C
       IF(JOB.EQ.63) CALL GLOBAL_CALL(63)
C
C Split communicator into groups
C
       IF(JOB.EQ.64) CALL GLOBAL_CALL(64)
#endif
C
C Allocate SIC matrix
C
       IF(JOB.EQ.65) CALL GLOBAL_CALL(65)
#ifdef GROUP
C
C Run APOTNL in group mode
C
       IF(JOB.EQ.66) CALL GLOBAL_CALL(66)
C
C Initialize RESULTS array to zero
C
       IF(JOB.EQ.67) CALL GLOBAL_CALL(67)
C
C Allocate ACOUL_SHARED array
C
       IF(JOB.EQ.69) CALL GLOBAL_CALL(69)
C
C RUN DENSOLD_SIC_GROUP
C
       IF(JOB.EQ.70) CALL GLOBAL_CALL(70)
#endif

       PRTIME2=MPI_WTIME()
       PRTIME=PRTIME+PRTIME2-PRTIME1
      END DO
c
c send total CPU time 
c
      TAG= 12
      CALL MPI_SSEND(PRTIME,1,MPI_DOUBLE_PRECISION,0,
     &               TAG,MPI_COMM_WORLD,IERR)   
#endif
      RETURN
      END
