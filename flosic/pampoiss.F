C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************

c
      SUBROUTINE PAMPOISS(MODE,ND,MD,MINC,LNG,NWRD,POT,RHOG2)
c
c subroutine for parallel poisson evaluation called by manager 
c 09/16/97 David C. Patton
c updates and name change by D. Porezag
c
c Mode 1: A free node is found via a call to gettid and the 
c         data from COUPOT1 is sent to the node for processing
c         (a call of POISSON2) (job2)
c Mode 2: The resultant data (pot & rhog) is collected from the nodes
c         and accumulated in the existing array pot (job 3)
c
#ifdef MPI
      use coupdata1
      use mesh1,only : nmsh
      use common2,only : IGGA, NSPN
      use common7,only : MODDEN
      use mpidat1,only : NPROC, NCALLED, IRANK, MYCOMM
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       REAL*8 :: SYMBOL , POT, RHOG2
       INTEGER :: MXPR,MXLG
      INCLUDE 'mpif.h'
      SAVE
      PARAMETER (MXPR=MXPOISS)
      PARAMETER (MXLG=3)
      LOGICAL   NWRD
      INTEGER   MODE,ND,MD,MINC,LNG
      DIMENSION POT(NMSH),RHOG2(NMSH,KRHOG,MXSPN)
C      COMMON/COUPDATA/
C     &  AIV(3,MXPR,MX_CNT+2,MXLG),AJV(3,MXPR,MX_CNT+2,MXLG)
C     &  ,DMTV(10,10,MXPR,MX_CNT+2,MXLG),ALPIV(MXPR,MX_CNT+2,MXLG)
C     &  ,ALPJV(MXPR,MX_CNT+2,MXLG),CENTER(3,MX_CNT+2)
C     &  ,ADD(MAXUNSYM,MAXUNSYM,2),RVECI(3,MX_GRP),RVECJ(3,MX_GRP)
C     &  ,NPAIRS(MX_CNT+2),IP(MX_CNT+2,MXLG)
C      COMMON/TMP2/ADDBUF(MAX_PTS)
      REAL*8,ALLOCATABLE :: ADDBUF(:)
      INTEGER TID,JOB,TAG,NPAIR,IERR,NGRAD,IPROC,I,J,KSPN,ITRANS(3)
c
c produce fatal error if I am a worker
c
      IF (IRANK.NE.0) THEN
       write(6,*)'FATAL: PAMPOISS CALLED BY WORKER'
       CALL STOPIT
      END IF 
c
c mode= 1: send data to workers
c
      IF (MODE.EQ.1) THEN
       NCALLED=NCALLED+1
       CALL GETTID(TID)
       JOB= 3
       TAG= 0
       CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 301
       CALL MPI_SSEND(NWRD,1,MPI_LOGICAL,TID,TAG,MPI_COMM_WORLD,IERR)
       NPAIR=IP(MINC,LNG)
       ITRANS(1)= NPAIR
       ITRANS(2)= ND
       ITRANS(3)= MD
       TAG= 302
       CALL MPI_SSEND(ITRANS(1),3,MPI_INTEGER,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 303
       CALL MPI_SSEND(ALPIV(1,MINC,LNG),NPAIR,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 304
       CALL MPI_SSEND(AIV(1,1,MINC,LNG),3*NPAIR,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 305
       CALL MPI_SSEND(ALPJV(1,MINC,LNG),NPAIR,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 306
       CALL MPI_SSEND(AJV(1,1,MINC,LNG),3*NPAIR,MPI_DOUBLE_PRECISION,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
       TAG= 307
       CALL MPI_SSEND(DMTV(1,1,1,MINC,LNG),10*10*NPAIR,
     &                MPI_DOUBLE_PRECISION,TID,TAG,MPI_COMM_WORLD,IERR)
      END IF
c
c manager mode= 2: receive results from workers
c
      IF (MODE.EQ.2) THEN
       ALLOCATE(ADDBUF(MAX_PTS),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'PAMPOISS:ERROR ALLOCATING ADDBUF'
       CALL CKWORKER(3,TID)
       NGRAD=1
       IF ((IGGA(1).GT.0).OR.(IGGA(2).GT.0)) NGRAD=10
       JOB= 4
       TAG= 0
       DO IPROC=1,NPROC
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
       END DO
C       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR) 
       CALL MPI_REDUCE(POT(1),ADDBUF(1),NMSH,MPI_DOUBLE_PRECISION,
     &                 MPI_SUM,0,MPI_COMM_WORLD,IERR)
       DO I=1,NMSH
        POT(I)=ADDBUF(I)
       END DO        
       IF (MODDEN .EQ. 1) THEN
        DO KSPN=1,NSPN
         DO J=1,NGRAD
C          CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          CALL MPI_REDUCE(RHOG2(1,J,KSPN),ADDBUF(1),NMSH,
     &                    MPI_DOUBLE_PRECISION,MPI_SUM,
     &                    0,MPI_COMM_WORLD,IERR)
          DO I=1,NMSH
           RHOG2(I,J,KSPN)=ADDBUF(I)
          END DO
         END DO
        END DO
       END IF
       DEALLOCATE(ADDBUF,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'PAMPOISS:ERROR DEALLOCATING ADDBUF'
      END IF
#endif
      RETURN
      END
