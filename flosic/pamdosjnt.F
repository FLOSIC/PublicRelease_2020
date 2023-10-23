C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PAMDOSJNT(SECTION,MODE,ITER)
c
c Subroutine for parallel joint density of states called by manager
C Luis Basurto 06/2013
c taken from
c 20/03/98 Jens Kortus   kortus@theo.physik.tu-freiberg.de
c modifications by DV Porezag 08/99
c
c Mode 1: looking for a free node and sending all data
c Mode 2: receiving data and doing a call do gethold
c
#ifdef MPI
      use debug1
      use dosjnt_mod,only: H,TEMP_RCV,TEMP_RCV2,SOS_POL,
     &     SOS_FREQ,SPDIP,SPTOT,CHARGE,NSPEC
      use common5,only : NWF
      use mpidat1, only : NCALLED,NPROC,MYCOMM
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER,INTENT(IN) :: SECTION,MODE,ITER
      INTEGER TID,JOB,TAG,IERR,IRECVSTAT(MPI_STATUS_SIZE)
      INTEGER I,J,K,IIP,IPTS,NTID
      REAL*8 POL_TEMP,TIME1,TIME2
c
c produce fatal error if I am a worker
c
      IF (IRANK.NE.0) THEN
       write(6,*)'FATAL: PAMiDOSJNT CALLED BY WORKER'
       CALL STOPIT
      END IF         

      SELECT CASE(SECTION)
      CASE (1)
c
c send data for joint dos execution section 1
c
      IF (MODE.EQ.1) THEN
        NCALLED=NCALLED+1
        CALL GETTID(TID)
C        CALL TRACER('SENDING',ITER)
C        CALL TRACER('TO',TID)
        JOB= 26
        TAG= 0
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
        TAG= 2601
        CALL MPI_SSEND(ITER,1,MPI_INTEGER,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
      END IF  
c
c receive data from workers
C
      IF (MODE.EQ.2) THEN
        ALLOCATE(TEMP_RCV2(NWF,NWF,3),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*) 'ERROR ALLOCATING TEMP_RCV2'
        CALL CKWORKER(3,NTID)
        JOB= 27
        TAG= 0
C        CALL TRACER('TELLING ALL TO SEND RESULTS')
        DO IIP=1,NPROC
          CALL MPI_SSEND(JOB,1,MPI_INTEGER,IIP,TAG,MPI_COMM_WORLD,IERR)
        ENDDO
C        CALL TRACER('REDUCING CHARGE')
        CALL MPI_REDUCE(CHARGE,POL_TEMP,1,MPI_DOUBLE_PRECISION,
     &               MPI_SUM,0,MPI_COMM_WORLD,IERR)
        CHARGE=POL_TEMP
C        CALL TRACER('REDUCING H')
C        CALL TRACER('NWF=',NWF)
        CALL MPI_ALLREDUCE(H,TEMP_RCV2,NWF*NWF*3,MPI_DOUBLE_PRECISION,
     &               MPI_SUM,MPI_COMM_WORLD,IERR)
        DO K=1,3
          DO J=1,NWF
            DO I=1,NWF
              H(I,J,K)=TEMP_RCV2(I,J,K)
            ENDDO
          ENDDO
        ENDDO
C        CALL TRACER('DONE REDUCING SECTION 1')
        DEALLOCATE(TEMP_RCV2,STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*) 'ERROR DEALLOCATING TEMP_RCV2'
      END IF
      RETURN
      CASE (3)
c
c send data for joint dos execution section 3
c
      IF (MODE.EQ.1) THEN
        NCALLED=NCALLED+1
        CALL GETTID(TID)
C        CALL TRACER('SENDING',ITER)
C        CALL TRACER('TO',TID)
        JOB= 24
        TAG= 0
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
        TAG= 2401
        CALL MPI_SSEND(ITER,1,MPI_INTEGER,
     &                TID,TAG,MPI_COMM_WORLD,IERR)
      END IF  
c
c receive data from workers
C
      IF (MODE.EQ.2) THEN
        ALLOCATE(TEMP_RCV(NSPEC),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'ERROR ALLOCATING TEMP_RCV'
        CALL CKWORKER(3,NTID)
        JOB= 25
        TAG= 0
        DO IIP=1,NPROC
          CALL MPI_SSEND(JOB,1,MPI_INTEGER,IIP,TAG,MPI_COMM_WORLD,IERR)
        ENDDO
        CALL MPI_REDUCE(SOS_POL,POL_TEMP,1,MPI_DOUBLE_PRECISION,
     &               MPI_SUM,0,MPI_COMM_WORLD,IERR)
        SOS_POL=POL_TEMP
        CALL MPI_REDUCE(SPDIP(1),TEMP_RCV(1),NSPEC,
     &       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
        DO IPTS=1,NSPEC
          SPDIP(IPTS)=TEMP_RCV(IPTS)
        ENDDO
        CALL MPI_REDUCE(SOS_FREQ(1),TEMP_RCV(1),NSPEC,
     &       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
        DO IPTS=1,NSPEC
          SOS_FREQ(IPTS)=TEMP_RCV(IPTS)
        ENDDO
        CALL MPI_REDUCE(SPTOT(1),TEMP_RCV(1),NSPEC,
     &       MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
        DO IPTS=1,NSPEC
          SPTOT(IPTS)=TEMP_RCV(IPTS)
        ENDDO
        DEALLOCATE(TEMP_RCV,STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*)'PASERROR DEALLOCATING TEMP_RCV'
C        DEALLOCATE(H,STAT=IERR)
C        IF(IERR.NE.0) WRITE(6,*)'ERROR DEALLOCATING TEMP_RCV'
      END IF
      RETURN
      END SELECT
#endif  
      END
