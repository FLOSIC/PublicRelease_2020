C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PASHAMIL_LIBXC(MODE)
c
c J Kortus 1998
c code taken from cpuhog and moved into PASHAMIL by DV Porezag 08/1999
c
c Mode 1: receive data for call to gethold
c Mode 2: send data from call to gethold
c
#ifdef MPI
       use mpidat1,only : IRANK, MYCOMM
       use XTMP2A,only  : ismgga
       INCLUDE 'PARAMA2'
       INCLUDE 'mpif.h'
       SAVE
       INTEGER MODE
       INTEGER TAG,IERR,IRECVSTAT(MPI_STATUS_SIZE)
       INTEGER IFNCT,JFNCT,MYIP,NEWIT,MXXX,ITRANS(3)
       REAL*8 :: AI,AJ,HSUB
       DIMENSION AI(3),AJ(3),HSUB(MAXUNSYM,MAXUNSYM)
c
c produce fatal error if I am the manager
c
       IF (IRANK.EQ.0) THEN
        write(6,*)'FATAL: PASHAMIL CALLED BY MANAGER'
        CALL STOPIT
       END IF         
c
c receive data for parallel hamiltonian
c
       IF (MODE.EQ.1) THEN
        TAG= 501
        CALL MPI_RECV(ITRANS(1),3,MPI_INTEGER,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        IFNCT= ITRANS(1)
        JFNCT= ITRANS(2)
        MYIP=  ITRANS(3)
        TAG= 502
        CALL MPI_RECV(AI(1),3,MPI_DOUBLE_PRECISION,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
        TAG= 503
        CALL MPI_RECV(AJ(1),3,MPI_DOUBLE_PRECISION,
     &                0,TAG,MPI_COMM_WORLD,IRECVSTAT,IERR)
c
c call gethold and tell mom I am done
c
        NEWIT=1
!        CALL GETHOLD2(NEWIT,IFNCT,JFNCT,AI,AJ,HSUB)
        if(ismgga)then
         CALL GETHOLD_MGGA(NEWIT,IFNCT,JFNCT,AI,AJ,HSUB)
        else
         CALL GETHOLD_LIBXC(NEWIT,IFNCT,JFNCT,AI,AJ,HSUB)
        end if
        TAG=1
        CALL MPI_SSEND(IRANK,1,MPI_INTEGER,0,TAG,
     &                 MPI_COMM_WORLD,IERR)   
       ENDIF
c
c send gethold data to manager
c
       IF (MODE.EQ.2) THEN
        MXXX=MAXUNSYM*MAXUNSYM
        TAG= 601
        CALL MPI_SSEND(HSUB(1,1),MXXX,MPI_DOUBLE_PRECISION,
     &                 0,TAG,MPI_COMM_WORLD,IERR)
        TAG= 602
        CALL MPI_SSEND(MYIP,1,MPI_INTEGER,0,TAG,MPI_COMM_WORLD,IERR)
       ENDIF
#endif
       RETURN
       END
