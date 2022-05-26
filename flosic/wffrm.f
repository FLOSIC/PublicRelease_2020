C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
C Fermi_orbitals are written out on a cubic grid for visualization
C purpose
C
       SUBROUTINE WFFRM
       use global_inputs,only : wffrm1
       use MOCORB,only : NFRM
       use common2,only : LSYMMAX,NSPN,NATOM,N_CON,N_POS,
     &                    NCNT,NFNCT,RCNT,RIDT
       use common5,only : NWF,NWFS,EVLOCC,PSI  ,N_OCC
       use common8,only : N_REP,NDMREP,NS_TOT
       use mesh1,only : NMSH
       INCLUDE 'PARAMA2'
!      INCLUDE 'commons.inc'
!      INTEGER,PARAMETER :: NTOTWF=200
       LOGICAL EXIST,IUPDAT,WFF
       CHARACTER*20 FORMSTR
       CHARACTER*20 BAND, TRASH
       CHARACTER*6 FILEN
!      CHARACTER*9 FILENAME(NTOTWF)
       CHARACTER*4 caltype
       INTEGER :: I,II,IW,IX,IY,IZ,IPV,IWF,
     &            IBAS,ICNT,ICON,IERR,ILOC,IOFS,IPTS,IREP,ISPN,
     &            IFNCT,IFORM,ISIZE,ITYPE,IUNIT,I_POS,ISHDUM,ISHELLA,
     &            ITOTWF,
     &            J,JJ,JOFS,JRP,JSP,JST,J_POS,
     &            K,KSPN,K_REP,
     &            LI,LMAX1,L_NUC,
     &            MPTS,MU,M_NUC,
     &            NDIM,NDM,NPV,NSTORE,NUNSYM,NWAVF
     &            ,IWOFS,IDIM,IORBX,LFN
       REAL*8 :: CHG(MAX_OCC),CHR,FACT,FACTOR,X,Y,Z
!       COMMON/TMP2/PSIG(MPBLOCK,MAX_OCC)
!     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
!     &  ,RVECA(3,MX_GRP),RBAS(3,4),RGRID(3,MPBLOCK)
!     &  ,NGRID(3),INFOWF(4,MAX_OCC),ICOUNT(MAX_CON,3)
       REAL*8,ALLOCATABLE :: PSIG(:,:),PTS(:,:)
     &                      ,GRAD(:,:,:,:,:),RVECA(:,:),RGRID(:,:)
       REAL*8 :: RBAS(3,4)
       INTEGER :: NGRID(3)
       LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
       INTEGER,ALLOCATABLE :: INFOWF(:,:)
       DIMENSION ISIZE(3)
!      DIMENSION INDWRTHOMO(NTOTWF),ISP(NTOTWF),IRP(NTOTWF),IBS(NTOTWF)
       INTEGER,ALLOCATABLE :: INDWRTHOMO(:),ISP(:),IRP(:),IBS(:)
       CHARACTER*9,ALLOCATABLE :: FILENAME(:)
       DATA ISIZE/1,3,6/
       DATA FACT/0.148203857088/
C
C READ FILE WFGRID WHICH CONTAINS: 
C * MODE (1=UNFORMATTED, 2=FORMATTED)
C * NUMBER OF POINTS FOR THE THREE BASIS VECTORS
C * ORIGIN
C * BASIS VECTORS
C * NUMBER OF STATES
c * COR/VAL STATES; VAL REQUIRES LEVELS ORDERED WRT HOMO: E.G. +1, -1
C * SPIN, REPRESENTATION, INDEX FOR EACH STATE (ORDERED) FOR OPTION COR
C
       CALL CHECK_INPUTS
!       INQUIRE(FILE='WFRMGRID',EXIST=EXIST)
       FORMSTR= ' '
       IF (.NOT.WFFRM1) FORMSTR= ' --> NOTHING TO DO'
       PRINT '(2A)','WAVEFUNCTION GRID',FORMSTR
       IF (.NOT. WFFRM1) RETURN

       NWF=0
       DO ISPN=1,NSPN
         NWF=NWF+NWFS(ISPN)
       END DO

       ALLOCATE(INDWRTHOMO(NWF))
       ALLOCATE(ISP(NWF))
       ALLOCATE(IRP(NWF))
       ALLOCATE(IBS(NWF))
       ALLOCATE(FILENAME(NWF))

       NWAVF=21
       NSTORE=NMSH
!       OPEN(80,FILE='WFRMGRID',FORM='FORMATTED',STATUS='OLD')
       OPEN(80,FILE='WFRMGRID',FORM='FORMATTED')
       REWIND(80)
       READ(80,*,END=110) ITYPE,IFORM
       GO TO 120
 110   CONTINUE
       REWIND(80) 
       

       ITYPE=1
       IFORM=2
       DO I=1,3
        RBAS(I,1)=1.0D30
        RBAS(I,2)=-1.0D30
       END DO
       DO ICNT=1,NCNT
        DO I=1,3
         RBAS(I,1)=MIN(RBAS(I,1),RCNT(I,ICNT))
         RBAS(I,2)=MAX(RBAS(I,2),RCNT(I,ICNT))
        END DO
       END DO
       DO I=1,3
        RBAS(I,1)=RBAS(I,1)-5.0D0
        RBAS(I,2)=RBAS(I,2)+5.0D0
        NGRID(I)=(RBAS(I,2)-RBAS(I,1))/0.5D0+2
        RBAS(I,2)=(RBAS(I,2)-RBAS(I,1))/(NGRID(I)-1)
       END DO
       WRITE(80,1010) 1,2, 'GRID MODE, FORMATTED FILE IN CUBE FORMAT'
       WRITE(80,1020) (NGRID(I),I=1,3),   'NUMBER OF GRID POINTS'
       WRITE(80,1030) (RBAS(I,1),I=1,3),   'ORIGIN'
       WRITE(80,1030) RBAS(1,2), 0.0D0, 0.0D0,  'BASIS VECTOR 1'
       WRITE(80,1030) 0.0D0, RBAS(2,2), 0.0D0,  'BASIS VECTOR 2'
       WRITE(80,1030) 0.0D0, 0.0D0, RBAS(1,2),  'BASIS VECTOR 1'
       NWAVF=NWF    !sum(NFRM(1:NSPN))
!       WRITE(80,'(I3,A)') NWAVF, '    NUMBER OF ORBITALS '
       WRITE(80,'(I3,A)') sum(NFRM(1:NSPN)), '    NUMBER OF ORBITALS '
       BAND='ALL' 
!       DO ii=1,NWF
!        INDWRTHOMO(ii)= ii
!       END DO
       DO ii=1,NFRM(1)
        INDWRTHOMO(ii)= ii
       END DO
       DO ii=1,NFRM(2)
        INDWRTHOMO(NFRM(1)+ii)=NWFS(1)+ii
       END DO
       WRITE(80,'(A)') 'VAL    IF VALENCE GIVE 
     &   INDEX WITH RESPECT TO HOMO'
!       WRITE(80,'(21I4,2X)') (INDWRTHOMO(jj),jj=1,NWAVF) 
       WRITE(80,'(21I4,2X)') (INDWRTHOMO(jj),jj=1,sum(NFRM(1:NSPN)))
! YY Can't skip to 230 since you need to assign a format string
!    FORMSTR='FORMATTED' and reread basis vectors.
!       GO TO 230
 1010  FORMAT( 2(I6,1X), 8X,A)
 1020  FORMAT( 3(I6,1X), 1X,A)
 1030  FORMAT( 3(F7.3,1X), 1X,A)
       REWIND(80)
       READ(80,*) ITYPE,IFORM
 120   CONTINUE

       IF (ITYPE .GT. 2) ITYPE= 2
       IF (ITYPE .LT. 1) ITYPE= 1
       IF (IFORM .GT. 2) IFORM= 2
       IF (IFORM .LT. 1) IFORM= 1
       IF (ITYPE .EQ. 2) IFORM= 2
       FORMSTR='UNFORMATTED'
       IF (IFORM .EQ. 2) FORMSTR='FORMATTED'
C
C GET NUMBER OF GRID POINTS, ORIGIN, AND BASIS VECTORS
C
       READ(80,*,END=110)(NGRID(J), J=1,3)
       IF (ITYPE .EQ. 1) THEN
        DO I=1,4
         READ(80,*,END=110)(RBAS(J,I), J=1,3)
        END DO
       END IF
       DO I=1,3
        IF (NGRID(I) .LT. 1) THEN
         PRINT *,'WFGRID: NUMBER OF GRID POINTS MUST BE >= 1'
         GOTO 920
        END IF
       END DO
C
C
       READ(80,*,END=110) NWAVF
!      IF (NWAVF.GT.NTOTWF) THEN
!       WRITE(6,*)'WFGRID: CHANGE NTOTWF TO ',NWAVF
!       RETURN
!      END IF
       if (NWAVF.GT.NWF) NWAVF=NWF
 220   READ(80,225,END=110) BAND
 225   FORMAT(A20)
C
       NUNSYM=0
       READ(80,*,END=910) (INDWRTHOMO(I),I=1,NWAVF)
       NWAVF=NWF !YY Swapping NWAVF to NWF
 230   CONTINUE
C Allocate INFOWF
       ALLOCATE(INFOWF(4,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING INFOWF'
       DO I=1,NWAVF
!YY alternative loop
!       DO II=1,NWAVF
!        I=INDWRTHOMO(II)
!exclude unoccupied states/non Fermi orbitals
        !if(I .GT. NWFS(1)+NFRM(2)) cycle
        !if(I .GT. NFRM(1) .AND. I .LE. NWFS(1)) cycle

        IRP(I)=1
        IF (I.LE.NWFS(1)) THEN !YY Changed .LT. to .LE.
         ISP(I)=1
        ELSE
         ISP(I)=NSPN
        END IF
        IF(ISP(I).EQ.1) THEN
         IBS(I)=I
        ELSE
         IBS(I)=I-NWFS(1)
        ENDIF
!        WRITE(6,*)'STATE :', ISP(I),IRP(I),IBS(I)
        NDIM=NDMREP(IRP(I))
        IF ((ISP(I) .LT. 1) .OR. (ISP(I) .GT. NSPN)) THEN
         PRINT *,'WFGRID: SPIN FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IRP(I).LT. 1) .OR. (IRP(I).GT. N_REP)) THEN
         PRINT *,'WFGRID: REPRESENTATION FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IBS(I) .LT. 1) .OR. (IBS(I) .GT. NS_TOT(IRP(I)))) THEN
         PRINT *,'WFGRID: INDEX FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        NDIM=NDMREP(IRP(I))
        NUNSYM=NUNSYM+NDIM
        IF (NUNSYM .GT. MAX_OCC) THEN
         PRINT *,'WFGRID: NUMBER OF STATES IS TOO LARGE'
         PRINT *,'        INCREASE MAX_OCC TO AT LEAST: ',NUNSYM
         GOTO 920
        END IF
        INFOWF(1,I)=ISP(I)
        INFOWF(2,I)=IRP(I)
        INFOWF(3,I)=IBS(I)
        INFOWF(4,I)=NDIM
       END DO
         
************************************************************************
C
C
   30  PRINT '(A)','SETUP OF WAVEFUNCTIONS ... '
c
C
       FILEN='WFFRMI'
!        print *, "NWAVF PUSKAR",NWAVF
       DO I=1,NWAVF
!YY exclude non Fermi orbitals
        if(I .GT. NWFS(1)+NFRM(2)) cycle
        if(I .GT. NFRM(1) .AND. I .LE. NWFS(1)) cycle

        II=I
!        II=INDWRTHOMO(I)
        IF (II.LT.100) THEN
         WRITE(FILENAME(I),'(A,I2.2)') FILEN, II
        ELSE
         WRITE(FILENAME(I),'(A,I3.3)') FILEN, II
        END IF

C
C WRITE FILE HEADER
C Grid data RBAS converted to angstrom PUSKAR 01/24/2022
        IUNIT=90+I
        OPEN(IUNIT,FILE=FILENAME(I),FORM=FORMSTR,STATUS='UNKNOWN')
        REWIND(IUNIT)
        IF (IFORM .EQ. 1) THEN
         WRITE(IUNIT) ITYPE,NSPN
         WRITE(IUNIT)(NGRID(J), J=1,3),MPBLOCK
         IF (ITYPE .EQ. 1) THEN
          WRITE(IUNIT)((RBAS(J,K)*0.529177d0, J=1,3), K=1,4)
         END IF
         WRITE(IUNIT) NWAVF
         WRITE(IUNIT)(INFOWF(J,I), J=1,4),EVLOCC(I)
        ELSE
         JSP=INFOWF(1,I)
         JRP=INFOWF(2,I)
         JST=INFOWF(3,I)
         WRITE(IUNIT,*) 'ORBITAL DENSITY FOR WF'
         WRITE(IUNIT,'(A,I1,A,I2,A,I5)')'SPIN ',JSP,
     &        ' REPRESENTATION ',JRP, ' STATE ',JST
         OPEN(77,FILE='XMOL.DAT')
         REWIND(77)
         READ(77,*) NATOM
         READ(77,*)
         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,
     &   (0.529177d0*RBAS(J,1),J=1,3)
         DO K=1,3
          WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(K),
     &   (0.529177d0*RBAS(J,K+1),J=1,3)
         ENDDO
         DO K=1,NATOM
          READ(77,*)IZ, X, Y, Z
          CHR=REAL(IZ) 
!PB converted Angstrom to Bohr
          WRITE(IUNIT,2002)IZ, CHR, X, Y, Z

         END DO
         CLOSE(77)
 2002    FORMAT(I6,4F16.10)
        END IF
       END DO
C
C ALLOCATE ARRAYS
C
       ALLOCATE(PSIG(MPBLOCK,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING PSIG'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING NSPEED'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING RVECA'
       ALLOCATE(RGRID(3,MPBLOCK),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING RGRID'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR ALLOCATING ICOUNT'
C
C LOOP OVER ALL POINTS IN BLOCKS
C
       ITOTWF=0
       DO IW=1,NWAVF
        K_REP=INFOWF(2,IW)
        NDM=NDMREP(K_REP)
        ITOTWF=ITOTWF+NDM
       END DO
!        print *,"ITOTWF PUSKAR", ITOTWF
!PB changed
!        CHG=0.0d0
       DO 800 IX=1,NGRID(1)
        DO 790 IY=1,NGRID(2)
         NMSH=NGRID(3)
         DO 780 IOFS=0,NMSH-1,MPBLOCK
          MPTS=MIN(MPBLOCK,NMSH-IOFS)
C
C SETUP GRID POINTS AND INITIALIZE PSIG
C
          IF (ITYPE .EQ. 1) THEN
           DO IPTS=1,MPTS
            IZ=IOFS+IPTS
            DO I=1,3
             RGRID(I,IPTS)=RBAS(I,1)+(IX-1)*RBAS(I,2)
     &                    +(IY-1)*RBAS(I,3)+(IZ-1)*RBAS(I,4)
            END DO
           END DO
          ELSE
           DO IPTS=1,MPTS
            READ(80,*,END=900)(RGRID(I,IPTS), I=1,3)
           END DO
          END IF
          DO IWF=1,MAX_OCC
           DO IPTS=1,MPTS
            PSIG(IPTS,IWF)=0.0D0
           END DO
          END DO
C
C LOOP OVER ALL FUNCTION SETS, THEIR POSITIONS, EQUIVALENT SITES
C
          ISHELLA=0
          DO IFNCT=1,NFNCT
           LMAX1=LSYMMAX(IFNCT)+1
           DO I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUM)
            DO J_POS=1,M_NUC
!            CALL LORAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
!     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
! YY. According to Mark, it is safer to use unravel2 than loravel.
! In MD=2 mode, you don't really need spin index since all occupied
! and unoccupied orbitals are computed. KSPN=NSPN is used as a 
! dummy index here.
             KSPN=NSPN
             CALL UNRAVEL2(2,KSPN,0,IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &            RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
             IF(L_NUC .NE. M_NUC)THEN
              PRINT *,'WFGRID: PROBLEM IN UNRAVEL'
              CALL STOPIT
             END IF
C
C FOR ALL MESHPOINTS IN BLOCK DO A SMALLER BLOCK
C
             DO JOFS=0,MPTS-1,NSPEED
              NPV=MIN(NSPEED,MPTS-JOFS)
              DO IPV=1,NPV
               PTS(IPV,1)=RGRID(1,JOFS+IPV)-RVECA(1,J_POS)
               PTS(IPV,2)=RGRID(2,JOFS+IPV)-RVECA(2,J_POS)
               PTS(IPV,3)=RGRID(3,JOFS+IPV)-RVECA(3,J_POS)
              END DO
C
C GET VALUE OF BASIS FUNCTIONS
C
              CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
C
C UPDATE PSIG
C

              IF (IUPDAT) THEN
               DO IWF=1,ITOTWF
                ILOC=0
                DO LI=1,LMAX1
                 DO MU=1,ISIZE(LI)
                  DO ICON=1,N_CON(LI,IFNCT)
                   ILOC=ILOC+1
                   FACTOR=PSI(ILOC,IWF,1)
                   if(abs(FACTOR).GT. 1.0D-10) then
                    DO IPV=1,NPV
                     PSIG(JOFS+IPV,IWF)=PSIG(JOFS+IPV,IWF)
     &                                 +FACTOR*GRAD(IPV,1,MU,ICON,LI)
                    END DO
                   end if
                  END DO
                 END DO
                END DO
               END DO
              END IF
             END DO
            END DO
           END DO
          END DO
C
C GET ORBITAL DENSITY FROM WAVEFUNCTION 
C

          IWF=0
          DO IW=1,NWAVF
           ISPN=INFOWF(1,IW)
           IREP=INFOWF(2,IW)
           IBAS=INFOWF(3,IW)
           NDIM=NDMREP(IREP)
           IWF=IWF+1
           DO IPTS=1,MPTS
             PSIG(IPTS,IW)=PSIG(IPTS,IWF)**2
             CHG(IW)=CHG(IW)+PSIG(IPTS,IW)*RBAS(3,4)**3
           END DO
C          DO IDIM=2,NDIM
C           IWF=IWF+1
C           DO IPTS=1,MPTS
C            PSIG(IPTS,IW)=PSIG(IPTS,IW)
C     &                   +PSIG(IPTS,IWF)**2
C           END DO
C          END DO
          END DO
!          print *,"NWF",NWF
          IF(IWF.GT.ITOTWF) THEN
           WRITE(6,*) 'WFGRID: THE TOTAL NUMBER OF STATES DO NOT MATCH'
           WRITE(6,*)  IWF, ' VS ', NWF
           CALL STOPIT
          END IF
C
C WRITE OUTPUT
C
          DO IWF=1,NWAVF
!YY alternative loop
!          DO II=1,NWAVF
!           IWF=INDWRTHOMO(II)
!exclude unoccupied states/non Fermi orbitals
           if(IWF .GT. NWFS(1)+NFRM(2)) cycle
           if(IWF .GT. NFRM(1) .AND. IWF .LE. NWFS(1)) cycle

           IUNIT=90+IWF
!PSIG is converted from Bohr^-3 to Angstrom^-3 
           IF (ITYPE .EQ. 1) THEN
            IF (IFORM .EQ. 1) THEN
             WRITE(IUNIT)(PSIG(IPTS,IWF)/FACT, IPTS=1,MPTS)
            ELSE
             WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF)/FACT,
     &        IPTS=1,MPTS)
            END IF
           ELSE
            DO IPTS=1,MPTS
             WRITE(IUNIT,'(3(1X,F20.12))')(RGRID(I,IPTS), I=1,3)
             WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF))/FACT
            END DO
           END IF
          END DO
  780    CONTINUE
  790   CONTINUE
  800  CONTINUE
!       print *,"CHARGE PUSKAR"
!       print 737,(CHG(IW),IW=1,NFRM(1))
!  737  FORMAT(5(F10.6))
C
C Deallocate arrays
C
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING PSIG'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING NSPEED'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING GRAD'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING RVECA'
       DEALLOCATE(RGRID,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING RGRID'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING ICOUNT'
       DEALLOCATE(INFOWF,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'ERROR DEALLOCATING INFOWF'

       DEALLOCATE(INDWRTHOMO)
       DEALLOCATE(ISP)
       DEALLOCATE(IRP)
       DEALLOCATE(IBS)
       DEALLOCATE(FILENAME)

       CLOSE(80)
       DO I=1,NWAVF
        IUNIT=90+I
        CLOSE(IUNIT)
       END DO
       NMSH=NSTORE
       RETURN
C
C ERROR HANDLING
C
  900  CLOSE(80)
  910  PRINT *,'WFGRID: ERROR READING FILE WFRMGRID' 
  920  CLOSE(90)
       NMSH=NSTORE
       RETURN
      END

