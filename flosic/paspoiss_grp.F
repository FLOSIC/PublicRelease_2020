C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASPOISS_GRP(MODE)
c
c subroutine for parallel poisson evaluation called by worker
c 09/16/97 David C. Patton
c moved out of cpuhog into PSPOISS (D. Porezag, 08/1999)
c
c Mode 1: The data is received from the manager and poisson is
c         called  (job 2)
c Mode 2: The resultant data (pot & rhog) are sent to the manager (job 3)
c
#ifdef MPI
      use debug1
      use mesh1,only : nmsh
      use pot_dens,only : COULOMB, RHOG
      use common2,only : IGGA, NSPN
      use common7,only : MODDEN
      use mpidat1,only : SHMRANK, SHMCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IB, J, KSPN,DATA1(12)
       REAL*8 :: A , ALPHA, B, BETA, RHO, TMPRECV
      INCLUDE 'mpif.h'
      SAVE
C      COMMON/TMP1/ACOULOMB(MAX_PTS),ARHOG(MAX_PTS,KRHOG,MXSPN)
      LOGICAL   NWRD
      INTEGER   MODE,TAG,NGRAD,NPAIR,ND,MD,IERR
      INTEGER   ITRANS(3),IRECVSTAT(MPI_STATUS_SIZE)
      DIMENSION ALPHA(MXPOISS),BETA(MXPOISS)
      DIMENSION A(3*MXPOISS),B(3*MXPOISS)
      DIMENSION RHO(10*10*MXPOISS)       
      COMMON/TMP3/TMPRECV(MAX_PTS)
c
c produce fatal error if I am the manager
c
      IF (SHMRANK.EQ.0) THEN
       write(6,*)'FATAL: PASPOISS CALLED BY MANAGER'
       CALL STOPIT
      END IF
c
c worker mode= 1: receive data from manager, call poisson
c
      IF (MODE.EQ.1) THEN
       TAG= 301
       CALL MPI_RECV(NWRD,1,MPI_LOGICAL,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       TAG= 302
       CALL MPI_RECV(ITRANS(1),3,MPI_INTEGER,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       NPAIR= ITRANS(1)
       ND=    ITRANS(2)
       MD=    ITRANS(3)
       TAG= 303
       CALL MPI_RECV(ALPHA(1),NPAIR,MPI_DOUBLE_PRECISION,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       TAG= 304
       CALL MPI_RECV(A(1),3*NPAIR,MPI_DOUBLE_PRECISION,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       TAG= 305
       CALL MPI_RECV(BETA(1),NPAIR,MPI_DOUBLE_PRECISION,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)                 
       TAG= 306
       CALL MPI_RECV(B(1),3*NPAIR,MPI_DOUBLE_PRECISION,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       TAG= 307
       CALL MPI_RECV(RHO(1),10*10*NPAIR,MPI_DOUBLE_PRECISION,
     &               0,TAG,SHMCOMM,IRECVSTAT,IERR)
       CALL POISSON2(NWRD,NPAIR,ND,MD,ALPHA,A,BETA,B,RHO)
c
c tell mom I am done
c
!       CALL TRACER('I AM DONE',SHMRANK)
       TAG=1
       CALL MPI_SSEND(SHMRANK,1,MPI_INTEGER,0,TAG,SHMCOMM,IERR)
      END IF                    
c
c worker mode= 2: send potential data to manager via reduce
c
      IF (MODE.EQ.2) THEN
       NGRAD=1
       IF ((IGGA(1).GT.0) .OR. (IGGA(2).GT.0)) NGRAD=10
       TAG=401
C       CALL MPI_BARRIER(SHMCOMM,IERR)
       CALL MPI_REDUCE(COULOMB(1),TMPRECV(1),NMSH,MPI_DOUBLE_PRECISION,
     &                 MPI_SUM,0,SHMCOMM,IERR)
           DO IB=1,NMSH
            COULOMB(IB)=TMPRECV(IB)
           END DO
       IF (MODDEN .EQ. 1) THEN
        DO KSPN=1,NSPN
         DO J=1,NGRAD
C          CALL MPI_BARRIER(SHMCOMM,IERR)
          CALL MPI_REDUCE(RHOG(1,J,KSPN),TMPRECV(1),NMSH,
     &                    MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                    SHMCOMM,IERR)
          DO IB=1,NMSH
           RHOG(IB,J,KSPN)=TMPRECV(IB)
          END DO
         END DO
        END DO
       END IF
      END IF
#endif                 
      RETURN
      END
