C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASFPUL(MODE)
c
c calculate forces in parallel
c called by workers DVP 10/98  
c
#ifdef MPI
       use debug1
       use global_inputs,only : LIBXC1
       use common2,only : NIDENT, FRC1, FRC2,
     &                    over1, ek, dftV, allH !<<<<
       use mpidat1,only : IRANK, MYCOMM
       use xtmp2a,only : MIXINS,ISMGGA
       use SICFLAG,only : LSICF !<<<<
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:55 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IID, IX
       REAL*8 :: BUF1 , BUF2, 
     &           BUF3, BUF4, BUF5, BUF6 !<<<<
       INCLUDE 'mpif.h'
       SAVE
       INTEGER MODE
       INTEGER TAG,IERR,MPTS,LPTS,ITRANS(2),IRECVSTAT(MPI_STATUS_SIZE)
       COMMON/TMP3/BUF1(3,MAX_IDENT),BUF2(3,MAX_IDENT),
     &  BUF3(MAX_OCC,MAX_OCC),BUF4(MAX_OCC,MAX_OCC),  !<<<<
     &  BUF5(MAX_OCC,MAX_OCC),BUF6(MAX_OCC,MAX_OCC)   !<<<<
c
c produce fatal error if I am the manager
c
       IF (IRANK.EQ.0) THEN
        write(6,*)'FATAL: PASFPUL CALLED BY MANAGER'
        CALL STOPIT
       END IF         
c
c mode= 1: receive data for parallel forces and call frcslv 
c
       IF (MODE .EQ. 1) THEN
        TAG= 1101
        CALL MPI_RECV(ITRANS(1),2,MPI_INTEGER,0,
     &                TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        MPTS= ITRANS(1)
        LPTS= ITRANS(2)
        IF(LIBXC1.OR.ISMGGA)THEN
          ALLOCATE(MIXINS(4,MPTS*MXSPN),STAT=IERR)
          IF(IERR.NE.0)WRITE(6,*)'PASFPUL:ERROR ALLOCATING MIXINS'
          TAG=1102
          CALL MPI_RECV(MIXINS,4*MPTS*MXSPN,MPI_DOUBLE_PRECISION,0,
     &                TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        ENDIF
        IF(LSICF) THEN
         CALL FRCSLV0(MPTS,LPTS)
        ELSE
         CALL FRCSLV(MPTS,LPTS)
        END IF
        IF(LIBXC1.OR.ISMGGA)THEN
          DEALLOCATE(MIXINS,STAT=IERR)
          IF(IERR.NE.0)WRITE(6,*)'PASFPUL:ERROR DEALLOCATING MIXINS'
        ENDIF
        TAG=1
        CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,IERR)
       END IF
c
c mode= 2: send force data created by frcslv 
c
       IF (MODE.EQ.2) THEN
C       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        CALL MPI_REDUCE(FRC1(1,1),BUF1(1,1),3*NIDENT,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                  MPI_COMM_WORLD,IERR) 
        CALL MPI_REDUCE(FRC2(1,1),BUF2(1,1),3*NIDENT,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                  MPI_COMM_WORLD,IERR) 
!<<<<<<<<<<<<<<<<<<<
        IF(LSICF) THEN
         CALL MPI_REDUCE(over1(1,1),BUF3(1,1),max_occ**2,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(ek(1,1),BUF4(1,1),max_occ**2,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(dftV(1,1),BUF5(1,1),max_occ**2,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(allH(1,1),BUF6(1,1),max_occ**2,
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
        END IF 
!>>>>>>>>>>>>>>>>>>>
        DO IX=1,3
         DO IID=1,NIDENT
          FRC1(IX,IID)=BUF1(IX,IID)
          FRC2(IX,IID)=BUF2(IX,IID)
         END DO
        END DO
       END IF
#endif
      END
