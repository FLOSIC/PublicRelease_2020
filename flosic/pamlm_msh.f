C UTEP Electronic Structure Lab (2020)
      SUBROUTINE PAMLM_MSH(MODE,LSPN,LFM,MPTS,LPTS)
      use debug1
      use mpidat1,only : NPROC,NCALLED,IRANK
      use SICMAT,only : SIC
      use FRM,only    : BFRM,RESULTS,LFRM
      use mesh1,only : NMSH,FO_MESH
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:56 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
      INCLUDE 'mpif.h'
      INCLUDE 'commons.inc' 
      INTEGER :: IERR,IPROC,JOB,JWF,LFM,MPTS,NMAX,LSPN,LPTS
      REAL*8 :: ADD0,ADD1,ADD2,ADD3,RPTS,SICP,WMSA
      PARAMETER (NMAX=MPBLOCK) 
      INTEGER TID,ITRANS(4),MODE,TAG,I,SUM1
      INTEGER IRECVSTAT(MPI_STATUS_SIZE)
!     LOGICAL BUF
!     DIMENSION WMSA(NMAX),SICP(NMAX),RPTS(3,NMAX)
      DIMENSION ADD1(NMAX),ADD2(NMAX),ADD3(3,NMAX),ADD0(3)
!     DIMENSION BUF(NMSH)
CCC PUT ALL OF SICLGMATs specific common blocks here.
C     COMMON/MIXPOT/POTIN(MAX_PTS*MXSPN),POT(MAX_PTS*MXSPN)
C     COMMON/TMP1/COULOMB(MAX_PTS),RHOG(MAX_PTS,10,MXSPN)
C     COMMON/TMP2/PSIG(NMAX,10,MAX_OCC)
C    & ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
C    & ,RVECA(3,MX_GRP),ICOUNT(MAX_CON,3)
!      COMMON/SICMAT/SIC(MAX_OCC,MAX_OCC,MXSPN)
!      COMMON/FRM/BFRM(3,MAX_OCC,MXSPN),RESULTS(13, MAX_OCC,MXSPN),
!     &  LFRM(MXSPN),DEBDAX(3,MAX_OCC,MXSPN)

      IF (IRANK.NE.0) THEN
        PRINT *,'FATAL: PAMLMSIC CALLED BY WORKER'
        CALL STOPIT
      END IF
!      PRINT*,'IN PAMLM_MSH MODE=',MODE
      IF(MODE.EQ.1)THEN
        NCALLED=NCALLED+1
        CALL GTNTID(TID)
!       PRINT*,'WORKER:',TID,' IS ABOUT TO BE CALLED'
        JOB=36 ! CHECK THIS
        TAG=0  ! I AM THE MANAGER (0)
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)

        ITRANS(1)=LSPN
        ITRANS(2)=LFM
        ITRANS(3)=MPTS
        ITRANS(4)=LPTS
        TAG=2401
C       PRINT*,IRANK,'PAM LFM IS ZERO:',LFM
        CALL MPI_SSEND(ITRANS(1),4,MPI_INTEGER,
     &                 TID,TAG,MPI_COMM_WORLD,IERR)
!       PRINT*,'I AM RETURNING TO FRMORB_MSH',IRANK
        IF(NCALLED.EQ.NPROC)THEN
C         PRINT*,'WAITING FOR SOMEONE TO FINISH'
          CALL CKWORKER(4,TID)
C         PRINT*,'DONE WAITING'
        END IF  
      END IF
c=============================================
C     PRINT*,'MODE:',MODE
      IF (MODE .EQ. 2 ) THEN
C     PRINT*,'NCALLED:',NCALLED,'SHOULD BE ZERO'
        IF(NCALLED.NE.0)CALL STOPIT
        IF(LFM.EQ.0)THEN
          PRINT*,'LFM IS ZERO:',LFM
          CALL STOPIT
        END IF
        DO IPROC=1,NPROC
C        PRINT*,'IPROC',IPROC
         JOB=37 ! CHECK THIS
         TAG=0  ! I AM THE MANAGER (0)
         CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
        END DO
!        DO JWF=1,MAX_OCC
!          SIC(JWF,LFM,LSPN)= 0.0D0          
!        END DO
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!        BUF(:)=0
        CALL MPI_REDUCE(MPI_IN_PLACE,FO_MESH,NMSH,
     &                   MPI_LOGICAL,MPI_LOR,0,
     &                   MPI_COMM_WORLD,IERR)
!        SUM1=0
!        DO I=1,NMSH
!          FO_MESH(I)=BUF(I)
!        END DO
!        DO I=1,NMSH
!          IF(FO_MESH(I)==1) SUM1=SUM1+1
!        ENDDO
!        CALL TRACER('SUM',SUM1)
      END IF

      RETURN
      END
