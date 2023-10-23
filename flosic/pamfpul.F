C UTEP Electronic Structure Lab (2020)
c
c *******************************************************************
c
      SUBROUTINE PAMFPUL(MODE,MPTS,LPTS)
c
c calculate forces in parallel
c called by the manager DVP 10/98
c
#ifdef MPI
       use global_inputs,only : LIBXC1
       use common2,only : NIDENT, FRC1, FRC2,
     &                    over1, ek, dftV, allH !<<<<
       use common5,only : NWF !<<<<
       use mpidat1,only : NPROC, NCALLED, IRANK, MYCOMM
       use xtmp2a,only  : MIXINS,ISMGGA
       use SICFLAG,only : LSICF
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:54 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IID, IPROC, IX, NTID
       INCLUDE 'mpif.h'
       SAVE
       INTEGER MODE,MPTS,LPTS
C       COMMON/TMP2/FADD1(3,MAX_IDENT),FADD2(3,MAX_IDENT)
       REAL*8,ALLOCATABLE :: FADD1(:,:),FADD2(:,:),
     &        over1ADD(:,:),ekADD(:,:),dftVADD(:,:),allHADD(:,:) !<<<<
       INTEGER :: TID,JOB,TAG,IERR,ITRANS(2),
     &            IWF,JWF !<<<<
c
c produce fatal error if I am a worker
c
       IF (IRANK.NE.0) THEN
        write(6,*)'FATAL: PAMFPUL CALLED BY WORKER'
        CALL STOPIT
       END IF         
c
c mode= 1: send force data to workers
c
       IF (MODE .EQ. 1) THEN
        NCALLED= NCALLED+1
        CALL GETTID(TID)
        JOB= 11
        TAG= 0
        CALL MPI_SSEND(JOB,1,MPI_INTEGER,TID,TAG,MPI_COMM_WORLD,IERR)
        ITRANS(1)= MPTS
        ITRANS(2)= LPTS
        TAG= 1101
        CALL MPI_SSEND(ITRANS(1),2,MPI_INTEGER,TID,
     &                 TAG,MPI_COMM_WORLD,IERR)
        IF(LIBXC1.OR.ISMGGA)THEN
          TAG=1102
          CALL MPI_SSEND(MIXINS,4*MPTS*MXSPN,MPI_DOUBLE_PRECISION,TID,
     &        TAG,MPI_COMM_WORLD,IERR)
        ENDIF
       END IF
c
c mode= 2: receive results from workers
c
       IF (MODE.EQ.2) THEN
        ALLOCATE(FADD1(3,MAX_IDENT),STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING FADD1'
        ALLOCATE(FADD2(3,MAX_IDENT),STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING FADD2'
!<<<<<<<<<<<
        IF(LSICF) THEN
         ALLOCATE(over1ADD(max_occ,max_occ),STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING over1'
         ALLOCATE(ekADD(max_occ,max_occ),STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING ek'
         ALLOCATE(dftVADD(max_occ,max_occ),STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING dftV'
         ALLOCATE(allHADD(max_occ,max_occ),STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR ALLOCATING allH'
        END IF
!>>>>>>>>>>>>
        CALL CKWORKER(3,NTID)
        DO IPROC=1,NPROC
         JOB= 12
         TAG= 0
         CALL MPI_SSEND(JOB,1,MPI_INTEGER,IPROC,TAG,MPI_COMM_WORLD,IERR)
        END DO
C       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        CALL MPI_REDUCE(FRC1(1,1),FADD1(1,1),3*NIDENT,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                  MPI_COMM_WORLD,IERR) 
        CALL MPI_REDUCE(FRC2(1,1),FADD2(1,1),3*NIDENT,
     &                  MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                  MPI_COMM_WORLD,IERR) 
!<<<<<<<<<<<
        IF(LSICF) THEN
         CALL MPI_REDUCE(over1(1,1),over1ADD(1,1),max_occ**2, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(ek(1,1),ekADD(1,1),max_occ**2, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(dftV(1,1),dftVADD(1,1),max_occ**2, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
         CALL MPI_REDUCE(allH(1,1),allHADD(1,1),max_occ**2, 
     &                   MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                   MPI_COMM_WORLD,IERR) 
        END IF
!>>>>>>>>>>>
        DO IID= 1,NIDENT
         DO IX= 1,3
          FRC1(IX,IID)= FADD1(IX,IID)
          FRC2(IX,IID)= FADD2(IX,IID)
         END DO
        END DO
!<<<<<<<<<<<<<<<<<
        IF(LSICF) THEN
c
c  transfer overadd to over1
c
         do iwf = 1,NWF
          do jwf = 1, NWF
           over1(jwf,iwf)=over1ADD(jwf,iwf)
           ek(jwf,iwf)=ekADD(jwf,iwf)
           dftV(jwf,iwf)=dftVADD(jwf,iwf)
           allH(jwf,iwf)=allHADD(jwf,iwf)
          end do
         end do
        END IF
!>>>>>>>>>>>>>>>>>
        DEALLOCATE(FADD1,STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING FADD1'
        DEALLOCATE(FADD2,STAT=IERR)
        IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING FADD2'
!<<<<<<<<<<<<<<
        IF(LSICF) THEN
         DEALLOCATE(over1ADD,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING over1ADD'
         DEALLOCATE(ekADD,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING ekADD'
         DEALLOCATE(dftVADD,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING dftVADD'
         DEALLOCATE(allHADD,STAT=IERR)
         IF(IERR/=0)WRITE(6,*)'PAMFPUL:ERROR DEALLOCATING allHADD'
        END IF
!>>>>>>>>>>>>>>
       END IF
#endif
      END
