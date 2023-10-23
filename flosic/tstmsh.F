C UTEP Electronic Structure Lab (2020)
C
      SUBROUTINE TSTMSH
C
C     DETERMINES THE LARGEST DEVIATION OF THE NUMERICAL OVERLAP MATRIX
C     FROM THE CORRECT OVERLAP MATRIX
C
C     WRITTEN BY MARK R PEDERSON (1997)
C
C     ------------------------------------------------------------------
C
      use debug1
      use global_inputs,only : inbas,iiev,iimesh,iinitial,mpi_io1
      use for_diag1
      use hstor1,only : hstor
      use mixpot1,only : POTIN,POTDV=>POTOUT
      use mesh1,only : WMSH,NMSH
      use common3,only : RMAT
      use common8,only : REP, N_REP, IGEN, NS_TOT
      use mpidat1,only : NPROC
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:59 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IBAS, IERR, IMODE, IOFFEST, IOFFSET, IPTS, IREP,
     & ISIZE, JBAS, KBAS1, KBAS2, NBAS, NDHTOT
       REAL*8 :: ERR1 , ERR2, ERR3
      SAVE
      LOGICAL EXIST
C
C     ------------------------------------------------------------------
C
C     --- CALCULATION OF NUMERICAL OVERLAP MATRIX ---
C
C     IMODE=1 USE O(N)+PARALLEL APPROACH
C     IMODE=2 USE STANDARD APPROACH - FASTER FOR SMALL
C             SYSTEMS / FEW CPUS
C      
      IMODE=2
C
#ifdef MPI
      IF (NPROC .GE. 3) IMODE=1
#endif
C
      INQUIRE(FILE='HAMMODE',EXIST=EXIST)
      IF (EXIST) THEN
        OPEN(99,FILE='HAMMODE',FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(99)
        READ(99,*) IMODE
        REWIND(99)
        WRITE(99,*) IMODE,' 1=O(N)/PARALLEL OR 2=STANDARD'
        CLOSE(99)
      END IF
C
      IF (DEBUG) write(6,*)'TSTMSH CALLS OVERLAP IN MODE: 1'
      CALL OVERLAP(1)
      DO IPTS=1,NMSH
        POTDV(IPTS)=WMSH(IPTS)
      END DO   
C
      DO I=1,NDHTOT
        HSTOR(I,2)=0.0D0
      END DO
C
      IF (IMODE.EQ.1) THEN
        CALL OVERNUM(0)
      ELSE
        CALL PATCH(POTDV)
      END IF
C
      ERR1=0.0D0
      ERR2=0.0D0
      ERR3=0.0D0
      KBAS1=0
      KBAS2=0
C
      DO 100 IREP=1,N_REP
        NBAS=NS_TOT(IREP)
        IF(NBAS.GT.NDH)THEN
          write(6,*)'TSTMSH: NDH MUST BE AT LEAST: ',NBAS
        CALL STOPIT
      END IF
C
      ISIZE=NBAS
      if(NBAS.eq.0)then
        ISIZE=1
      endif
C
CJUR  --- HERE WE ALLOCATE AHAM ---
C
C     ALLOCATE(AHAM(ISIZE,ISIZE),STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error allocating Ham'
C     ENDIF
CJUR
C     ALLOCATE(AOVER(ISIZE,ISIZE),STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error allocating Overlap'
C     ENDIF
      ALLOCATE(AEVAL(ISIZE),STAT=IERR)
      IF(IERR.NE.0)THEN
        WRITE(6,*)'tstmsh:Error allocating Eval'
      ENDIF
C     ALLOCATE(ASC1(ISIZE),STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error allocating Sc1'
C     ENDIF

      IOFFSET=KBAS1+1
C      KBAS1=0
      DO IBAS=1,NBAS
        DO JBAS=IBAS,NBAS
          KBAS1=KBAS1+1
C         AOVER(JBAS,IBAS)=HSTOR(KBAS1,1)
C         AHAM (JBAS,IBAS)=HSTOR(KBAS1,2)
          IF(IBAS.EQ.JBAS) THEN
            AEVAL(IBAS)=HSTOR(KBAS1,1)
            AEVAL(IBAS)= 1.0D0/SQRT(AEVAL(IBAS))
          ENDIF
        END DO  
C       AEVAL(IBAS)=AOVER(IBAS,IBAS)
C       ASC1(IBAS)=AOVER(IBAS,IBAS)
C       write(6,*)'AEVAL',AEVAL(IBAS),'ASC1',ASC1(IBAS)
      END DO 

C     DO IBAS=1,NBAS
C      AEVAL(IBAS)= 1.0D0/SQRT(AEVAL(IBAS))
C      ASC1(IBAS)= 1.0D0/SQRT(ASC1(IBAS))
C     END DO

C     KBAS=0
      DO IBAS=1,NBAS
C       ERR11=MAX(ERR11,ABS(AHAM(IBAS,IBAS)-AOVER(IBAS,IBAS))
C    &       *ASC1(IBAS)*ASC1(IBAS))
        DO JBAS=IBAS,NBAS
          KBAS2=KBAS2+1
          IF(IBAS.EQ.JBAS)THEN
            ERR1=MAX(ERR1,ABS(HSTOR(KBAS2,2)-
     &           HSTOR(KBAS2,1))*AEVAL(IBAS)*AEVAL(IBAS))
C           write(6,*)'ERR1=',ERR1,'ERR11',ERR11
C           write(6,*) 'HAM=',AHAM(IBAS,IBAS),'HSTOR=',
C    &                  HSTOR(KBAS2,2)
C           write(6,*) 'OVER=',AOVER(IBAS,IBAS),
C    &                 'HSTOR=',HSTOR(KBAS2,1)
          ENDIF
C         ERR22=MAX(ERR22,ABS(AHAM(JBAS,IBAS)-AOVER(JBAS,IBAS))
C    &         *ASC1(IBAS)*ASC1(JBAS))
          ERR2=MAX(ERR2,ABS(HSTOR(KBAS2,2)-
     &         HSTOR(KBAS2,1))*AEVAL(IBAS)*AEVAL(JBAS))
C         write(6,*) 'HAM=',AHAM(JBAS,IBAS),'HSTOR=',
C    &                HSTOR(KBAS2,2)
C         write(6,*) 'OVER=',AOVER(JBAS,IBAS),
C    &               'HSTOR=',HSTOR(KBAS2,1)
C         write(6,*)'ERR2=',ERR2,'ERR22',ERR22
        END DO  
      END DO  
      if(NBAS.NE.0)then
#ifdef SCALAPACK
        IF(NBAS>2000)THEN
          IINITIAL=IOFFEST
          iiev=0
          inbas=NBAS
          iimesh=.TRUE.
C         write(6,*)'tstmsh:before scalapack,NBAS=',NBAS
          IF(MPI_IO1)THEN
            CALL PREP_SCALA
          ENDIF
          CALL SENDDATA(111)
        ELSE
C         write(6,*)'tstmsh:before diagge3,NBAS=',NBAS,'JKIND=',JKIND
          CALL DIAGGE3(NBAS,IOFFSET,AHAM,AEVAL,0)
        ENDIF

#else
C       write(6,*)'tstmsh:before diagge3,NBAS=',NBAS,'JKIND=',JKIND
        CALL DIAGGE3(NBAS,IOFFSET,AHAM,AEVAL,0)
#endif
      END IF
C
C The decision block is needed for certain symmetric cases
C to avoid calling DIAGGE with NBAS=0
      IF (NBAS .GT. 0) THEN
        ERR3=MAX(ERR3,ABS(AEVAL(1)   -1.0D0))
        ERR3=MAX(ERR3,ABS(AEVAL(NBAS)-1.0D0))
      END IF
      IF (DEBUG) THEN
        IF (NBAS.GT.0) THEN
          PRINT 1010,IREP,NBAS,AEVAL(1),AEVAL(NBAS)
          PRINT 1020,(AEVAL(I),I=1,NBAS)
        ELSE
          PRINT 1030,IREP,NBAS
        END IF
      END IF
 1010   FORMAT(1X,2(1X,I4),2(1X,F12.6),' OVERLAP EIGENVALUES - MESH')
 1020   FORMAT(1X,6(1X,F10.7))
 1030   FORMAT(1X,2(1X,I4),' NO OVERLAP EIGENVALUES')
CJUR
C     --- AHAM IS NEEDED AND WE DEALLOCATE IT HERE ---
C
C     DEALLOCATE(AHAM,STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error deallocating Ham'
C     ENDIF
CJUR
C     DEALLOCATE(AOVER,STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error deallocating Overlap'
C     ENDIF
      DEALLOCATE(AEVAL,STAT=IERR)
      IF(IERR.NE.0)THEN
        WRITE(6,*)'tstmsh:Error deallocating Eval'
      ENDIF
C     DEALLOCATE(ASC1,STAT=IERR)
C     IF(IERR.NE.0)THEN
C       WRITE(6,*)'tstmsh:Error deallocating Sc1'
C     ENDIF
  100  CONTINUE
C
C     --- PRINT AND FORMATS ---
C
      PRINT 1040,ERR1
      PRINT 1050,ERR2
      PRINT 1060,ERR3
 1040 FORMAT('LARGEST S-DIAGONAL ERROR: ',D12.4)
 1050 FORMAT('LARGEST S-GENERAL  ERROR: ',D12.4)
 1060 FORMAT('LARGEST S-EIGEN    ERROR: ',D12.4)
C
C     call stopit
C
C     ------------------------------------------------------------------
C
      END
