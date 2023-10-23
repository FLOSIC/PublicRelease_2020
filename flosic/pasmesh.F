C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASMESH(MODE)
c
c 12/98 Dirk V. Porezag 
c
c mode 1: get box data from manager and call fillist
c mode 2: send mesh points to manager
c
#ifdef MPI
       use mesh1,only : wmsh,rmsh,nmsh
       use mpidat1,only : IRANK, MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: BOX , ERRBOX
       INCLUDE 'mpif.h'
       SAVE
       INTEGER MODE
       INTEGER TAG,IERR,IRECVSTAT(MPI_STATUS_SIZE)
       DIMENSION BOX(2,3)
       INTEGER NMULT
c
c produce fatal error if I am the manager
c
       IF (IRANK.EQ.0) THEN
        write(6,*)'FATAL: PASMESH CALLED BY MANAGER'
        CALL STOPIT
       END IF    
      
c
c receive input data for interstitial mesh and call fillist 
c
       IF (MODE.EQ.1) THEN
        TAG= 1401
        CALL MPI_RECV(BOX(1,1),6,MPI_DOUBLE_PRECISION,0,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
        TAG= 1402
        CALL MPI_RECV(ERRBOX,1,MPI_DOUBLE_PRECISION,0,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
        TAG= 1403
        CALL MPI_RECV(NMULT,1,MPI_INTEGER,0,TAG,
     &                MPI_COMM_WORLD,IRECVSTAT,IERR)
        CALL FILLIST(BOX,ERRBOX,NMULT)
        TAG=1
        CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,IERR)
       END IF                      
c
c send mesh points to manager
c  
       IF (MODE.EQ.2) THEN
        TAG= 1501
        CALL MPI_SSEND(NMSH,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,IERR)
        IF (NMSH .GT. 0) THEN
         TAG= 1502
         CALL MPI_SSEND(RMSH(1,1),3*NMSH,MPI_DOUBLE_PRECISION,0,TAG,
     &                  MPI_COMM_WORLD,IERR)
         TAG= 1503
         CALL MPI_SSEND(WMSH(1),NMSH,MPI_DOUBLE_PRECISION,0,TAG,
     &                  MPI_COMM_WORLD,IERR)
        END IF
       END IF
#endif                            
       RETURN
      END
