C UTEP Electronic Structure Lab (2020)
! This program must be compiled first wiht condcomp
! ./condcomp -DMPI -DSCALAPACK < fragment.ftn > fragment.f
! mpif90 -c fragment.f
! In order to use Scalapack
      program clusterfrag
        use global_inputs,only : idiag1,inbas,iiev,iinitial,
     &         iimesh,mpi_io1
        use for_diag1
        use zero1,only : ZEROREALMATRIX
        use mpidat1
        use hstor1, only : MAXCLUSTERSIZE
        implicit none
        include 'mpif.h'
        INTEGER :: NBAS,NREC,NTOT,ISPN,NSPN,NFRAG
        INTEGER :: I,J,ISPFAC,KBAS,IBAS,JBAS,IWFS,KTOT
        INTEGER :: MBAS,IFRAG,IB,JB,MSPN,MXNMO
        INTEGER :: IWF,JFRAG,IA,JA,MREC
        INTEGER :: IOC,JOC,JWFS
        INTEGER :: IERR
        REAL*8 :: A,ETOT,ERR,SUMC
        LOGICAL :: RESTART
        CHARACTER*4 NAME
        CHARACTER*5 FILENAME
        INTEGER, allocatable,dimension(:) :: BAS, NMO
        INTEGER, allocatable,dimension(:,:) :: INDMO
        INTEGER, allocatable,dimension(:,:,:) :: NFRMO
!        real*8, allocatable,dimension(:):: AEVAL
        real*8, allocatable,dimension(:):: EVAL,SC1
!        real*8, allocatable,dimension(:,:):: AHAM,AOVER
        real*8, allocatable,dimension(:,:):: HAM,OVER,EVALFRAG
        INTEGER, allocatable,dimension(:,:):: NWFS
        INTEGER, allocatable,dimension(:):: NWF
        real*8, allocatable,dimension(:,:,:,:):: CHARGE
        CHARACTER*7, allocatable,dimension(:,:):: MOLABEL
        CHARACTER*7, allocatable,dimension(:):: COMLABEL
        
        real*8, allocatable,dimension(:,:,:):: PSI_COEFTOT,PSI_COEF
        real*8, allocatable,dimension(:,:):: HSTOR
! Use divide and conquer for diagonalization
#ifdef SCALAPACK
        CALL CHECK_INPUTS
        CALL MPI_INIT(IERR)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,IERR)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
        NPROC=NPROC-1
        IOUT=6
        MAXCLUSTERSIZE=200
        IF(IRANK.EQ.0) THEN
#endif
        idiag1=1
        OPEN(5,FILE='FRAGMENT.IN')
        READ(5,*) NSPN
        READ(5,*) NFRAG
        ALLOCATE(BAS(NFRAG))
        ALLOCATE(NWFS(NFRAG,NSPN))
        ALLOCATE(NWF(NSPN))
        ALLOCATE(CHARGE(NFRAG,2,NFRAG,2))

        NTOT=0
        ETOT=0.0D0
        READ(5,*)(NWF(ISPN),ISPN=1,NSPN)
        DO  I=1,NFRAG
          READ(5,*) BAS(I),(NWFS(I,ISPN),ISPN=1,NSPN)
          NTOT=NTOT+BAS(I)
          DO ISPN=1,NSPN
            ETOT=ETOT+NWFS(I,ISPN)
          END DO
        END DO
         
        READ(5,*)RESTART
        ALLOCATE(NMO(NSPN))
!        READ IN THE MOs FOR ANALYZING FOR EACH SPIN
    
        READ(5,*,END=20) (NMO(ISPN),ISPN=1,NSPN)
        MXNMO=MAX(NMO(1),NMO(NSPN))
        ALLOCATE(INDMO(MXNMO,NSPN))
        ALLOCATE(NFRMO(MXNMO,NFRAG,NSPN))

        DO ISPN=1,NSPN
          DO I=1,NMO(ISPN)
            READ(5,*,END=20)INDMO(I,ISPN),(NFRMO(I,J,ISPN),J=1,NFRAG)
          END DO
        END DO
        GO TO 30

 20     DO ISPN=1,NSPN
          NMO(ISPN)=11
          IF(NMO(ISPN).GT.NWF(ISPN)) NMO(ISPN)=NWF(ISPN)
        END DO
        MXNMO=MAX(NMO(1),NMO(NSPN))
        ALLOCATE(INDMO(MXNMO,NSPN))
        ALLOCATE(NFRMO(MXNMO,NFRAG,NSPN))

                  
        DO ISPN=1,NSPN
          I=0 
          JB=NMO(ISPN)/2
          DO IBAS=-JB, JB
            I=I+1
            INDMO(I,ISPN)=NWF(ISPN)+IBAS
            IF (INDMO(I,ISPN).LE.0) THEN
              I=I-1
            ELSE
              DO J=1,NFRAG
                NFRMO(I,J,ISPN)=10
                IF (NFRMO(I,J,ISPN).GT.NWFS(J,ISPN)) THEN
                   NFRMO(I,J,ISPN)=NWFS(J,ISPN)
                ENDIF
              END DO
            END IF
          END DO
               
        END DO
 30     CLOSE(5)

!        WRITE(6,*) 'DONE WITH INPUTS'
         
        ISPFAC=2/NSPN
        ETOT=ETOT*ISPFAC
!        WRITE(6,*) 'NUMBER OF ELECTRONS:', ETOT

        OPEN(20,FILE='OVLBABY',FORM='UNFORMATTED')

        READ(20)NREC
        ALLOCATE(HSTOR(NREC,2))
        READ(20) (HSTOR(I,1),I=1,NREC)
        A=DSQRT(1.0D0+8.0D0*NREC)
        NBAS=(INT(A)-1)/2
!        WRITE(6,*) 'TOTAL BASIS SET SIZE :',NBAS,NREC
        IF(NBAS.NE.NTOT) THEN
         WRITE(6,*) 'NTOT.NE.NBAS'
         STOP
        END IF
        CLOSE(20)

        ALLOCATE(AOVER(NBAS,NBAS))
        ALLOCATE(AHAM(NBAS,NBAS))
        ALLOCATE(PSI_COEFTOT(NBAS,NBAS,NSPN))
        ALLOCATE(AEVAL(NBAS))
        ALLOCATE(SC1(NBAS))
        ALLOCATE(PSI_COEF(NBAS,NBAS,NSPN))
        ALLOCATE(EVALFRAG(NBAS,NSPN))

        DO IBAS=1,NBAS
          DO JBAS=1,NBAS
            AOVER(JBAS,IBAS)=0.0D0
          END DO
        END DO

        KBAS=0
        ETOT=0
        DO IBAS=1,NBAS
          DO JBAS=IBAS,NBAS
            KBAS=KBAS+1
            AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
          END DO
          ETOT=ETOT+AOVER(IBAS,IBAS)
        END DO
!        WRITE(6,*)'OVER TRACE =', ETOT

        IF(KBAS.NE.NREC) THEN
          WRITE(6,*)'ATTENTION : KBAS.NE.NREC'
          STOP
        END IF
        OPEN(25,FILE='HAMOLD',FORM='UNFORMATTED')
        OPEN(60,FILE='OVERLAP',FORM='UNFORMATTED')
        OPEN(65,FILE='WAVEMO',FORM='UNFORMATTED')
        REWIND(60)
        REWIND(65)

        DO 240 ISPN=1,NSPN

        IF (RESTART) GO TO 300
        READ(25)NREC,MSPN
!        WRITE(6,*) 'NREC :',NREC
        IF(MSPN.EQ.ISPN) THEN
          READ(25) (HSTOR(I,2),I=1,NREC)
        END IF 

        DO IBAS=1,NBAS
          EVALFRAG(IBAS,ISPN)=0.0D0
          AEVAL(IBAS)=0.0D0
          DO JBAS=1,NBAS
            PSI_COEF(JBAS,IBAS,ISPN)=0.0D0
            AHAM(JBAS,IBAS)=0.0D0
          END DO
        END DO

        KBAS=0
        ETOT=0.0D0
        DO IBAS=1,NBAS
          DO JBAS=IBAS,NBAS
            KBAS=KBAS+1
            AHAM(JBAS,IBAS)=HSTOR(KBAS,2)
          END DO
        END DO

#ifdef SCALAPACK
        IF(NBAS>2000)THEN
          KTOT=KBAS
          inbas=nbas
          iiev=1
          iinitial=1
          iimesh=.TRUE.
          if(mpi_io1)CALL PREP_SCALA(KTOT,HSTOR)
          CALL SCALA_CALL(1)
        ELSE
          CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,SC1,1)
        ENDIF
#else
        CALL DIAGGE(NBAS,NBAS,AHAM,AOVER,AEVAL,SC1,1)
#endif
        DEALLOCATE(SC1)

        DO JBAS=1,NBAS
          DO IBAS=1,NBAS
            PSI_COEFTOT(IBAS,JBAS,ISPN)=AHAM(IBAS,JBAS)
          END DO
        END DO
           
!        WRITE(6,*) 'OCCUPIED EIGENVALUES OF THE COMPLEX FOR SPIN =',ISPN
!        WRITE(6,100)(AEVAL(I),I=1,NBAS)

100     format(5F14.6)


          KBAS=0
          DO IBAS=1,NBAS
            DO JBAS=IBAS,NBAS
               KBAS=KBAS+1
               AHAM(JBAS,IBAS)=HSTOR(KBAS,2)
               AOVER(JBAS,IBAS)=HSTOR(KBAS,1)
               AHAM(IBAS,JBAS)=AHAM(JBAS,IBAS)
               AOVER(IBAS,JBAS)=AOVER(JBAS,IBAS)
            END DO
          END DO

!       NOW DO THE FRAGMENTS
             
         IWFS=0
         ETOT=0.0D0

          ALLOCATE(HAM(NBAS,NBAS))
          ALLOCATE(OVER(NBAS,NBAS))
          ALLOCATE(EVAL(NBAS))
          ALLOCATE(SC1(NBAS))

! LB: Begin fragment loop         
         DO IFRAG=1,NFRAG
!          WRITE(6,*)'FRAGMENT', IFRAG
!           MBAS=BAS(IFRAG)
          DO IBAS=1,NREC
            HSTOR(IBAS,2)=0.0D0
          END DO
          NAME='HAM.'
          WRITE(FILENAME,'(A,I1.1)')NAME,IFRAG
!          WRITE(6,*)filename
          OPEN(35,FILE=FILENAME,FORM='UNFORMATTED')
           READ(35)MREC
           READ(35) (HSTOR(I,2),I=1,MREC)
           A=DSQRT(1.0D0+8.0D0*MREC)
           MBAS=(INT(A)-1)/2
!           WRITE(6,*) 'FRAGMENT BASIS SET SIZE :',MBAS,MREC,IFRAG
           CLOSE(35)

 

          DO IBAS=1,NBAS
            EVAL(IBAS)=0.0D0
            SC1(IBAS)=0.0D0
            DO JBAS=1,NBAS
              HAM(IBAS,JBAS)=0.0D0
              OVER(IBAS,JBAS)=0.0D0
            END DO
          END DO

          KBAS=0
          DO IBAS=1,MBAS
            DO JBAS=IBAS,MBAS
              KBAS=KBAS+1
              HAM(JBAS,IBAS)=HSTOR(KBAS,2)
            END DO
          END DO
        
          IB=0
          DO IBAS=IWFS+1,IWFS+MBAS
            IB=IB+1
            JB=0
            DO JBAS=IWFS+1,IWFS+MBAS
              JB=JB+1
!              HAM(JB,IB)=AHAM(JBAS,IBAS)
              OVER(JB,IB)=AOVER(JBAS,IBAS)
            END DO
          END DO
!          WRITE(6,*) IWFS+1, IWFS+MBAS, ':RANGE'

#ifdef SCALAPACK
        IF(NBAS>2000)THEN
          KBAS=0
          DO IBAS=1,MBAS
            DO JBAS=IBAS,MBAS
              KBAS=KBAS+1
              HSTOR(KBAS,1)=OVER(IBAS,JBAS)
            ENDDO
          ENDDO 
          INBAS=MBAS
          IIEV=1
          IINITIAL=1
          IIMESH=.TRUE.
          if(mpi_io1)CALL PREP_SCALA(KTOT,HSTOR)
          CALL SCALA_CALL(1)
          DO IBAS=1,MBAS
            EVAL(IBAS)=AEVAL(IBAS)
            DO JBAS=1,MBAS
               HAM(IBAS,JBAS)=AHAM(IBAS,JBAS)
            END DO
          END DO
        ELSE
          CALL DIAGGE(NBAS,MBAS,HAM,OVER,EVAL,SC1,1)
        ENDIF
#else
          CALL DIAGGE(NBAS,MBAS,HAM,OVER,EVAL,SC1,1)
#endif
            
          IB=0
          DO IBAS=IWFS+1,IWFS+MBAS
            IB=IB+1
            EVALFRAG(IBAS,ISPN)=EVAL(IB)
            JB=0
            DO JBAS=IWFS+1,IWFS+MBAS
              JB=JB+1
              PSI_COEF(JBAS,IBAS,ISPN)=HAM(JB,IB)
            END DO
          END DO
!          WRITE(6,100)(EVAL(I),I=1,MBAS)

          IWFS=IWFS+MBAS 
! LB: End fragment loop
         END DO
        
         DEALLOCATE(SC1)
         DEALLOCATE(HSTOR)

!     NOW FIND  MIXING CHARACTER FOR MOs OF INTEREST

!         ALLOCATE (HAM(NBAS,NBAS))  
!         ALLOCATE (OVER(NBAS,NBAS))  
!         CALL OMP_MAX_THREADS(I)
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             HAM(JBAS,IBAS)=0.0D0
             OVER(JBAS,IBAS)=0.0D0
             AHAM(JBAS,IBAS)=0.0D0
           END DO
         END DO
       
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               AHAM(IBAS,JBAS)=AHAM(IBAS,JBAS)+
     &         AOVER(IBAS,KBAS)*PSI_COEFTOT(KBAS,JBAS,ISPN)
             END DO
           END DO
         END DO

         
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               HAM(IBAS,JBAS)=HAM(IBAS,JBAS)+
     &         AHAM(KBAS,JBAS)*PSI_COEF(KBAS,IBAS,ISPN)
             END DO
           END DO
         END DO

!         WRITE(6,*) 'EXPANSION COEFFICIENTS : DONE'
!         WRITE(6,100) (HAM(JBAS,1),JBAS=1,NBAS)
!        HAM CONTAINS THE COEFFICIENTS FOR EXPANSION IN THE FRAGMENT MO
!        BASIS

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             AHAM(IBAS,JBAS)=0.0D0
             OVER(IBAS,JBAS)=0.0D0
           END DO
         END DO

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               AHAM(JBAS,IBAS)=AHAM(JBAS,IBAS)+
     &         AOVER(KBAS,IBAS)*PSI_COEF(KBAS,JBAS,ISPN)
             END DO
           END DO
         END DO

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               OVER(JBAS,IBAS)=OVER(JBAS,IBAS)+
     &         AHAM(JBAS,KBAS)*PSI_COEF(KBAS,IBAS,ISPN)
             END DO
           END DO
         END DO

         ERR=0.0D0
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
               PSI_COEFTOT(JBAS,IBAS,ISPN)=HAM(JBAS,IBAS)
               ERR=ERR+OVER(IBAS,JBAS)-OVER(JBAS,IBAS)
           END DO
         END DO

!         WRITE(6,*) 'OVERLAP OF FRAGMENT MOs: DONE'
!         WRITE(6,*)'ERR IN OVERLAP = ', ERR
!        OVER NOW CONTAINS THE OVERLAP OF THE FRAGMENT MOs
         DEALLOCATE(PSI_COEF) 

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             AHAM(JBAS,IBAS)=0.0D0
             AOVER(JBAS,IBAS)=0.0D0
           END DO
         END DO
 
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               AHAM(JBAS,IBAS)=AHAM(JBAS,IBAS)+
     &           OVER(KBAS,IBAS)*HAM(KBAS,JBAS)
             END DO
           END DO
         END DO

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               AOVER(IBAS,JBAS)=AOVER(IBAS,JBAS)+
     &         AHAM(JBAS,KBAS)*HAM(KBAS,IBAS)
             END DO
           END DO
         END DO
!
         ERR=0.0D0
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             ERR=ERR+ABS(AOVER(IBAS,JBAS)-AOVER(JBAS,IBAS))
             AHAM(JBAS,IBAS)=0.0D0
             PSI_COEFTOT(JBAS,IBAS,ISPN)=0.0D0
           END DO
         END DO
!         WRITE(6,*)'ERR IN PROJECTION OVERLAP : ', ERR
!         WRITE(6,*)(AOVER(IBAS,IBAS),IBAS=1,NBAS)
    
!         OVER=AOVER
         
!         ALLOCATE(EVAL(NBAS))
         ALLOCATE(SC1(2*NBAS))
         CALL LOWDEN(NBAS,NBAS)
          
!         DO IBAS=1,NBAS
!           DO JBAS=1,NBAS
!             AOVER(JBAS,IBAS)=0.0D0
!             DO KBAS=1,NBAS
!               AHAM(JBAS,IBAS)=AHAM(JBAS,IBAS)+ OVER(JBAS,KBAS)*HAM(KBAS,IBAS)
!             END DO
!           END DO
!         END DO
          
!         DO IBAS=1,NBAS
!           DO JBAS=1,NBAS
!             DO KBAS=1,NBAS
!               AOVER(JBAS,IBAS)=AOVER(JBAS,IBAS)+ AHAM(KBAS,IBAS)*HAM(KBAS,JBAS)
!             END DO
!           END DO
!         END DO

! ROTATE BACK 

         AOVER=0.0D0
         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               PSI_COEFTOT(JBAS,IBAS,ISPN)=PSI_COEFTOT(JBAS,IBAS,ISPN)+
     &         HAM(JBAS,KBAS)*AHAM(KBAS,IBAS)
             END DO
           END DO
         END DO

         WRITE(60) NBAS,ISPN
         DO IBAS=1,NBAS
            WRITE(60) (OVER(JBAS,IBAS),JBAS=1,NBAS)
         END DO
         WRITE(65) NBAS,ISPN
         DO IBAS=1,NBAS
            WRITE(65) (PSI_COEFTOT(JBAS,IBAS,ISPN),JBAS=1,NBAS)
         END DO
        GO TO 310

 300     CONTINUE
         IF (RESTART) THEN
            ALLOCATE(OVER(NBAS,NBAS))
            ALLOCATE(HAM(NBAS,NBAS))
            ALLOCATE(EVAL(NBAS))
            REWIND(60)
            READ(60) NBAS,MSPN
!            WRITE(6,*) NBAS,MSPN, '60'
            IF(MSPN.LT.ISPN) THEN
              DO IBAS=1,NBAS
                 READ(60) (OVER(JBAS,IBAS),JBAS=1,NBAS)
              END DO
            END IF
!            write(6,*) 'HERE DONE'
            DO IBAS=1,NBAS
              READ(60) (OVER(JBAS,IBAS),JBAS=1,NBAS)
            END DO
!            write(6,*) '1 HERE DONE'

            REWIND(65)
            READ(65) NBAS,MSPN
            IF(MSPN.LT.ISPN) THEN
              DO IBAS=1,NBAS
                READ(65) (PSI_COEFTOT(JBAS,IBAS,MSPN),JBAS=1,NBAS)
              END DO
            END IF
!            write(6,*) '2 HERE DONE'
            DO IBAS=1,NBAS
               READ(65) (PSI_COEFTOT(JBAS,IBAS,ISPN),JBAS=1,NBAS)
            END DO
         END IF
!                write(6,*) '3 HERE DONE'


 310     CONTINUE
         AOVER=0.0D0
         HAM=0.0D0

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               AOVER(JBAS,IBAS)=AOVER(JBAS,IBAS)+
     &         OVER(JBAS,KBAS)*PSI_COEFTOT(KBAS,IBAS,ISPN)
             END DO
           END DO
         END DO

         DO IBAS=1,NBAS
           DO JBAS=1,NBAS
             DO KBAS=1,NBAS
               HAM(JBAS,IBAS)=HAM(JBAS,IBAS)+
     &         AOVER(KBAS,JBAS)*PSI_COEFTOT(KBAS,IBAS,ISPN)
             END DO
           END DO
         END DO

!         WRITE(6,*) 'AOVER AFTER  LOWDEN'
         ERR=0.0D0
         DO IBAS=1,NBAS
!         WRITE(6,*) IBAS
!         WRITE(6,100) (PSI_COEFTOT(JBAS,IBAS,ISPN),JBAS=1,NBAS)
!         WRITE(6,100) (HAM(JBAS,IBAS),JBAS=1,NBAS)
           ERR=ERR+HAM(IBAS,IBAS)
         END DO
!          WRITE(6,*)'ERR = ', ERR
!          WRITE(6,*)(AOVER(IBAS,IBAS),IBAS=1,NBAS)

         ISPFAC=2/NSPN
         DO IBAS=1,NMO(ISPN)
           I=INDMO(IBAS,ISPN)
           IB=0
           EVAL=0.0D0
           AEVAL=0.0D0
           AHAM=0.0D0
           AOVER=OVER
           DO JBAS=1,NBAS
             EVAL(JBAS)=PSI_COEFTOT( JBAS,I,ISPN) 
           END DO
           CALL MULLIKEN(NBAS,EVAL,AOVER,AHAM,AEVAL)

           
!           WRITE(6,*)'WF =', I
!           WRITE(6,100) SUM(AEVAL)
!           WRITE(6,100) (AEVAL(JB),JB=1,NBAS)
           DO IFRAG=1,NFRAG
             KBAS=NFRMO(IBAS,IFRAG,ISPN)/2
             IF(KBAS.EQ.0) KBAS=1
             JB=IB+NWFS(IFRAG,ISPN)-KBAS
            DO JBAS=1,NFRMO(IBAS,IFRAG,ISPN)
               JB=JB+1
!               WRITE(6,*) JB, IFRAG, I, AEVAL(JB)
            END DO
            IB=IB+BAS(IFRAG)
           END DO
         END DO
        
         
!    ATTEMPT AT ELECTRON DONATION BETWEEN FRAGMENTS
!
        CHARGE=0.0D0
        DO IWF=1,NWF(ISPN)
            IWFS=0
            DO IFRAG=1, NFRAG
             DO IOC=1,2
                  IF(IOC.EQ.1) THEN
                    IA=IWFS+1
                    IB=IWFS+NWFS(IFRAG,ISPN)
                   ELSE
                    IA=IWFS+NWFS(IFRAG,ISPN)+1
                    IB=IWFS+BAS(IFRAG)
                   END IF
               JWFS=0
               DO  JFRAG=1, NFRAG
                 DO JOC=1,2
                  SUMC=0.0D0
                  IF(JOC.EQ.1) THEN
                    JA=JWFS+1
                    JB=JWFS+NWFS(JFRAG,ISPN)
                   ELSE
                    JA=JWFS+NWFS(JFRAG,ISPN)+1
                    JB=JWFS+BAS(JFRAG)
                   END IF
                   EVAL=0.0D0
                   AEVAL=0.0D0
                   AHAM=0.0D0
                   AOVER=OVER
                   DO JBAS=1,NBAS
                      EVAL(JBAS)=PSI_COEFTOT(JBAS,IWF,ISPN)
                   END DO
                   CALL MULLIKEN2(NBAS,IA,IB,JA,JB,EVAL,
     &             AOVER,AHAM,AEVAL)
                  DO JBAS=IA,IB
                    SUMC=SUMC+AEVAL(JBAS)
                  END DO
                  CHARGE(IFRAG,IOC,JFRAG,JOC)=
     &            CHARGE(IFRAG,IOC,JFRAG,JOC)+SUMC
                 END DO
                  JWFS=JWFS+BAS(JFRAG)
               END DO
              END DO
             IWFS=IWFS+BAS(IFRAG)
            END DO
         END DO
        
       WRITE(6,*)'FRAGEMENT,OCC/VIR,FRAGMENT,OCC/VIR,CHARGE 1:OCC,2:VIR'

            ISPFAC=2/NSPN
            SUMC=0.0D0
            DO IFRAG=1,NFRAG
              DO IOC=1,2
               DO JFRAG=1,NFRAG
                DO JOC=1,2
                 WRITE(6,150)IFRAG,IOC,JFRAG,JOC,ISPFAC*
     &                 CHARGE(IFRAG,IOC,JFRAG,JOC)
                 SUMC=SUMC+CHARGE(IFRAG,IOC,JFRAG,JOC)
                END DO
               END DO
              END DO
             END DO
!             WRITE(6,*) 'TOTAL  = ', SUMC*ISPFAC
               
 150    FORMAT(4I10,F12.6)
 240     CONTINUE
         CLOSE(25)
         CLOSE(60)
         CLOSE(65)
#ifdef SCALAPACK
         CALL SCALA_CALL(-1)
       ELSE
! CODE FOR REST OF PROCESSORS
! HERE THEY WAIT TO BE SUMMONED
         CALL TRIGGER
       ENDIF
       CALL MPI_FINALIZE(IERR)
#endif
        STOP
        END

        SUBROUTINE MULLIKEN(NBAS,EVAL,AOVER,AHAM,AEVAL)
        IMPLICIT NONE
        INTEGER :: NBAS,IBAS,JBAS
        REAL*8, DIMENSION(NBAS):: EVAL,AEVAL
        REAL*8, DIMENSION(NBAS,NBAS):: AOVER,AHAM
        AHAM=0.0D0
        AEVAL=0.0D0

         DO IBAS=1,NBAS
          DO JBAS=1,NBAS
             AHAM(JBAS,IBAS)=AHAM(JBAS,IBAS)+EVAL(JBAS)*EVAL(IBAS) 
          END DO
         END DO

         DO IBAS=1,NBAS
          DO JBAS=1,NBAS
             AEVAL(IBAS)=AEVAL(IBAS)+
     &              AHAM(IBAS,JBAS)*AOVER(JBAS,IBAS)
          END DO
         END DO

        RETURN
        END

        SUBROUTINE MULLIKEN2(NBAS,IA,IB,JA,JB,EVAL,AOVER,AHAM,AEVAL)
        IMPLICIT NONE
        INTEGER :: NBAS,IBAS,JBAS,IA,IB,JA,JB
        REAL*8, DIMENSION(NBAS):: EVAL,AEVAL
        REAL*8, DIMENSION(NBAS,NBAS):: AOVER,AHAM
        AHAM=0.0D0
        AEVAL=0.0D0
         DO IBAS=1,NBAS
          DO JBAS=1,NBAS
             AHAM(JBAS,IBAS)=EVAL(JBAS)*EVAL(IBAS) 
          END DO
         END DO

         DO IBAS=IA,IB
          DO JBAS=JA,JB
             AEVAL(IBAS)=AEVAL(IBAS)+AHAM(IBAS,JBAS)*
     &                     AOVER(JBAS,IBAS)
          END DO
         END DO

        RETURN
        END
