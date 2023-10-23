C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASDOSJNT(SECTION,MODE)
c
C Luis Basurto 06/2013
c Taken from
c J Kortus 1998
c
c Mode 1: receive data 
c Mode 2: send data 
c
#ifdef MPI
      USE DEBUG1
      use dosjnt_mod,only: H,TEMP_RCV,TEMP_RCV2,SOS_POL,
     &     SOS_FREQ,SPDIP,SPTOT,CHARGE,NSPEC
      use common5,only : NWF
      use mpidat1,only: MYCOMM
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER SECTION,MODE,ITER,I,J,K
      INTEGER TAG,IERR,IRECVSTAT(MPI_STATUS_SIZE)
      REAL*8 POL_TEMP
c
c produce fatal error if I am the manager
c
       IF (IRANK.EQ.0) THEN
        write(6,*)'FATAL: PASDOSJNT CALLED BY MANAGER'
        CALL STOPIT
       END IF         
       SELECT CASE (SECTION)
       CASE(1)
c
c receive data for parallel joint dos section 1
c
       IF (MODE.EQ.1) THEN
         TAG= 2601
         CALL MPI_RECV(ITER,1,MPI_INTEGER,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       
C
C Perform share of work
C
         CALL DOSJNT_S1(ITER)
C
C Tell mom I am done
C
         TAG=1
         CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,
     &                 MPI_COMM_WORLD,IERR)   
       ENDIF
c
c send joint dos data to manager
c
       IF (MODE.EQ.2) THEN
        ALLOCATE(TEMP_RCV2(NWF,NWF,3),STAT=IERR)
        IF(IERR.NE.0) WRITE(6,*) 'ERROR ALLOCATING TMP_RCV2'
        CALL MPI_REDUCE(CHARGE,POL_TEMP,1,MPI_DOUBLE_PRECISION,
     &               MPI_SUM,0,MPI_COMM_WORLD,IERR)
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
        IF(IERR.NE.0) WRITE(6,*) 'ERROR DEALLOCATING TMP_RCV2'
       ENDIF
       RETURN
       CASE(3)
c
c receive data for parallel joint dos section 3
c
       IF (MODE.EQ.1) THEN
         TAG= 2401
         CALL MPI_RECV(ITER,1,MPI_INTEGER,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
       
C
C Perform share of work
C
         CALL DOSJNT3(ITER)
C
C Tell mom I am done
C
         TAG=1
         CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,
     &                 MPI_COMM_WORLD,IERR)   
       ENDIF
c
c send joint dos data to manager
c
       IF (MODE.EQ.2) THEN
         ALLOCATE(TEMP_RCV(NSPEC),STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*) 'ERROR ALLOCATING TEMP_RCV'
         CALL MPI_REDUCE(SOS_POL,POL_TEMP,1,MPI_DOUBLE_PRECISION,
     &               MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SPDIP(1),TEMP_RCV(1),NSPEC,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SOS_FREQ(1),TEMP_RCV(1),NSPEC,
     &        MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         CALL MPI_REDUCE(SPTOT(1),TEMP_RCV(1),NSPEC,
     &         MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IERR)
         DEALLOCATE(TEMP_RCV,STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*)'PASDOSJNT:ERROR DEALLOCATING TEMP_RCV'
         DEALLOCATE(H,STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*)'PASDOSJNT:ERROR DEALLOCATING H'
         DEALLOCATE(SPTOT,STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*)'PASDOSJNT:ERROR DEALLOCATING SPTOT'
         DEALLOCATE(SPDIP,STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*)'PASDOSJNT:ERROR DEALLOCATING SPDIP'
         DEALLOCATE(SOS_FREQ,STAT=IERR)
         IF(IERR.NE.0) WRITE(6,*)'PASDOSJNT:ERROR DEALLOCATING SOS_FREQ'
       ENDIF
       RETURN
       END SELECT
#endif
       END
