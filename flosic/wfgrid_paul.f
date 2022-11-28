C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
C WFGRID DIRK POREZAG JUNE 1998 
C CALCULATION OF WAVEFUNCTION DENSITY ON A GRID OF POINTS
C
C MODIFIED TO WRITE THE OUTPUT IN THE GAUSSIAN CUBE FORMAT
C THE HOMO DENSITY IS PRINTED OUT BY DEFAULT  -- TB 04/03
C
       SUBROUTINE WFGRID
       use global_inputs,only : wfgrid1
       use mesh1,only : nmsh
       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
       use common9, only: old_mode
       use xmol, only: au2ang, num_atms, xmol_list
!       use xmol,only : AU2ANG
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NTOTWF, I, I_POS, IBAS, IBS, IBSO, ICNT, ICON, IDIM,
     & IERR, IFNCT, IFORM, II, ILOC, INDWRTHOMO, IOFS, IPTS, IPV, IREP,
     & IRP, IRPO, ISHDUM, ISHELLA, ISIZE, ISP, ISPO, ITOTWF, ITYPE,
     & IUNIT, IW, IWF, IX, IY, IZ, J, J_POS, JJ, JOFS, JRP, JSP, JST,
     & K, K_REP, L_NUC, LI, LMAX1, M_NUC, MPTS, MU, NATOM, NDIM, NDM,
     & NPV, NSTORE, NUNSYM, NWAVF, INDEX, N, IS, IR, IB, IDEG, IREPR,
     & ISPIN, IVIRT, IWFHOMO, JVIRT, KWF, IV, NBASE, NS
       REAL*8 :: SYMBOL , CHR, FACT, FACTOR, X, Y, Z, EVL, DIFF,
     & DIFFNEW, EF, EV
       SAVE
       PARAMETER (NTOTWF=200)
       LOGICAL EXIST,IUPDAT,WFF
       CHARACTER*20 FORMSTR
       CHARACTER*20 BAND, TRASH
       CHARACTER*6  FILE
       CHARACTER*15 FILENAME(NTOTWF)
       CHARACTER*5  EXT
       CHARACTER*4  caltype
       CHARACTER*1  spn
C       COMMON/TMP2/PSIG(MPBLOCK,MAX_OCC)
C     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
C     &  ,RVECA(3,MX_GRP),RBAS(3,4),RGRID(3,MPBLOCK)
C     &  ,NGRID(3),INFOWF(4,MAX_OCC),ICOUNT(MAX_CON,3)
       LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
       REAL*8,ALLOCATABLE :: PSIG(:,:),PTS(:,:),GRAD(:,:,:,:,:)
     &  ,RVECA(:,:),RGRID(:,:)
       INTEGER,ALLOCATABLE :: INFOWF(:,:)
       REAL*8 :: RBAS(3,4)
       INTEGER :: NGRID(3)
       DIMENSION ISIZE(3)
       DIMENSION INDWRTHOMO(NTOTWF),ISP(NTOTWF),IRP(NTOTWF),IBS(NTOTWF)
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
       if (old_mode) call check_inputs
C       INQUIRE(FILE='WFGRID',EXIST=EXIST)
       FORMSTR= ' '
       IF (.NOT.WFGRID1) FORMSTR= ' --> NOTHING TO DO'
       PRINT '(2A)','WAVEFUNCTION GRID',FORMSTR
       IF (.NOT.WFGRID1) RETURN
C
C ALLOCATE LOCAL ARRAYS
C
       ALLOCATE(PSIG(MPBLOCK,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING PSIG'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING PTS'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING GRAD'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING RVECA'
       ALLOCATE(RGRID(3,MPBLOCK),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING RGRID'
       ALLOCATE(INFOWF(4,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING INFOWF'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR ALLOCATING ICOUNT'

        NWF=0
        DO ISPN=1,NSPN
         NWF=NWF+NWFS(ISPN)
        END DO

       NWAVF=21
       NSTORE=NMSH
C       OPEN(80,FILE='WFGRID',FORM='FORMATTED',STATUS='OLD')
       OPEN(80,FILE='WFGRID',FORM='FORMATTED')
       REWIND(80)
       READ(80,*,END=110) ITYPE,IFORM
       GO TO 120
***********************TB******************************
 110   CONTINUE
       REWIND(80) 
       
C      CREATE A DEFAULT INPUT FILE   - TB
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
 1010  FORMAT( 2(I6,1X), 8X,A)
 1020  FORMAT( 3(I6,1X), 1X,A)
 1030  FORMAT( 3(F7.3,1X), 1X,A)
       REWIND(80)
       READ(80,*) ITYPE,IFORM
 120   CONTINUE
***********************TB******************************

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
C READ DATA FOR WAVEFUNCTIONS
C PRINT OUT THE HOMO DENSITY BY DEFAULT
C
       READ(80,*,END=210) NWAVF
       IF (NWAVF.GT.NTOTWF) THEN
        WRITE(6,*)'WFGRID: CHANGE NTOTWF TO ',NWAVF
        WRITE(6,*)'AND RECOMPILE AND RERUN'
        WRITE(6,*) 'TILL THEN, BYE BYE FROM WGRID'
        RETURN
       END IF
       if (NWAVF.GT.NWF) NWAVF=NWF
       GO TO 220
 210   CONTINUE
       BACKSPACE(80)
       WRITE(80,'(I3,A)') NWAVF, '    NUMBER OF ORBITALS '
C LB:
C 210   BACKSPACE(80)
C       WRITE(80,'(I3,A)') NWAVF, '    NUMBER OF ORBITALS '
C LB:
 220   READ(80,240,END=310) BAND
 240   FORMAT(A20)
       GO TO 320
 310   BAND='VAL' 
C LB:
C 310   BACKSPACE(80)
C       BAND='VAL'
C LB:
       jj=0
       DO ii=-10,10
       jj=jj+1
       INDWRTHOMO(jj)= ii
       END DO
       BACKSPACE(80)
       WRITE(80,'(2A)') 'VAL    IF VALENCE GIVE ',
     &   'INDEX WITH RESPECT TO HOMO'
       WRITE(80,'(21I4,2X,A)') (INDWRTHOMO(jj),jj=1,NWAVF), 
     &   '   FOR CORE STATES GIVE SPIN, REPRESENATION, BASIS INDEX'
       GOTO 230
 320   CONTINUE

       IF(BAND(1:3).NE.'VAL') THEN
       NUNSYM=0
       ISPO=0
       IRPO=0
       IBSO=0
       DO I=1,NWAVF
        READ(80,*,END=910) ISPN,IREP,IBAS
        IF ((ISPN .LT. 1) .OR. (ISPN .GT. NSPN)) THEN
         PRINT *,'WFGRID: SPIN FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IREP .LT. 1) .OR. (IREP .GT. N_REP)) THEN
         PRINT *,'WFGRID: REPRESENTATION FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        IF ((IBAS .LT. 1) .OR. (IBAS .GT. NS_TOT(IREP))) THEN
         PRINT *,'WFGRID: INDEX FOR STATE ',I,' IS INVALID'
         GOTO 920
        END IF
        NDIM=NDMREP(IREP)
        NUNSYM=NUNSYM+NDIM
        IF (NUNSYM .GT. MAX_OCC) THEN
         PRINT *,'WFGRID: NUMBER OF STATES IS TOO LARGE'
         PRINT *,'        INCREASE MAX_OCC TO AT LEAST: ',NUNSYM
         GOTO 920
        END IF


C
C CHECK FOR ORDERING
C
        IF (ISPN .LT. ISPO) GOTO 10
        IF (ISPN .GT. ISPO) THEN
         ISPO=ISPN
         IRPO=0
        END IF
        IF (IREP .LT. IRPO) GOTO 10
        IF (IREP .GT. IRPO) THEN
         IRPO=IREP
         IBSO=0
        END IF
        IF (IBAS .LT. IBSO) GOTO 10
        IBSO=IBAS
        INFOWF(1,I)=ISPN
        INFOWF(2,I)=IREP
        INFOWF(3,I)=IBAS
        INFOWF(4,I)=NDIM
       END DO
       GOTO 30
C
C WRONG ORDERING
C
   10  PRINT *,'WFGRID: STATES MUST BE ORDERED ACCORDING TO: '
       PRINT *,'      * SPIN           (HIGHEST PRIORITY)'
       PRINT *,'      * REPRESENTATION (SECOND HIGHEST PRIORITY)'
       PRINT *,'      * INDEX          (LOWEST PRIORITY)'
       GOTO 920

       ELSE
         NUNSYM=0
         READ(80,*,END=910) (INDWRTHOMO(I),I=1,NWAVF)
 230      CALL FINDSTATE(INDWRTHOMO,NWAVF,ISP,IRP,IBS,IERR)
         DO I=1,NWAVF
           WRITE(6,*)'STATE :', ISP(I),IRP(I),IBS(I)
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
       END IF
         
************************************************************************
C
C DIAGONALIZE FOR EACH SPIN AND REPRESENTATION 
C STORE RELEVANT STATES IN PSI_COEF 
C
C  THIS PART IS SKIPPED, THE STORED WF ARE USED -- TB
   30  PRINT '(A)','SETUP OF WAVEFUNCTIONS ... '
c       GO TO 200
c
c      IWOFS=0
c      NWF=0
c      DO 200 ISPN=1,NSPN
c        NWFS(ISPN)=0
c        IF (IWOFS .GE. NWAVF) GOTO 190
c        IF (INFOWF(1,IWOFS+1) .GT. ISPN) GOTO 190
c        CALL OVERLAP(1)
c        CALL OVERLAP(2)
C
C LOOP OVER REPRESENTATIONS
C GET MATRIX ELEMENTS AND DIAGONALIZE
C
c        KBAS=0
c        DO 100 IREP=1,N_REP
c         NBAS=NS_TOT(IREP)
c         IF ((N_OCC(IREP,ISPN) .LT. 1) .OR. (NBAS .LT. 1)) THEN
c          KBAS=KBAS+(NBAS*(NBAS+1))/2
c          GOTO 90
c         END IF
c         IF (NBAS .GT. NDH) THEN
c          PRINT *,'WFGRID: NDH MUST BE AT LEAST: ',NBAS
c          CALL STOPIT
c         END IF
c         DO IBAS=1,NBAS
c          DO JBAS=IBAS,NBAS
c           KBAS=KBAS+1
c           OVER(JBAS,IBAS)=HSTOR(KBAS,1)
c           HAM (JBAS,IBAS)=HSTOR(KBAS,2)
c          END DO
c         END DO
c         CALL DIAGGE(NDH,NBAS,HAM,OVER,EVAL,SC1,1)
c
C
C SAVE CORRECT EIGENSTATES
C
c         DO I=1,N_OCC(IREP,ISPN)
c          INDX=INFOWF(3,IWOFS+I)
c          EVLOCC(IWOFS+I)=EVAL(INDX)
c          OCCUPANCY(IWOFS+I)= 1.0D0
c          DO IB=1,NBAS
c           PSI_COEF(IB,I,IREP,ISPN)=HAM(IB,INDX)
c          END DO
c          PRINT '(A,I1,A,I2,A,I5,A,F15.6)',
c     &          'SPIN ',ISPN,', REPRESENTATION ',IREP,', STATE ',
c     &           INDX,', EIGENVALUE: ',EVAL(INDX)
c         END DO
c         IWOFS=IWOFS+N_OCC(IREP,ISPN)
c         NWFS(ISPN)=NWFS(ISPN)+N_OCC(IREP,ISPN)*NDMREP(IREP)
c   90    CONTINUE
cc  100   CONTINUE
c  190   CONTINUE
c        NWF=NWF+NWFS(ISPN)
c
c  200  CONTINUE

C
C OPEN OUTPUT FILE
C
       WRITE(6,*)'ISP(I)',ISP(I)
       IF(ISP(I).EQ.1)THEN
          SPN='a'
       ELSE
          SPN='b'
       ENDIF
       EXT='.cube'
       IF  (BAND(1:3).EQ.'VAL') THEN
          FILE='WFHOMO'
       ELSE
          FILE='WFCORE'
       END IF

       DO I=1,NWAVF
        IF(ISP(I).EQ.1)THEN
          SPN='a'
        ELSE
          SPN='b'
        ENDIF
        IF(BAND(1:3).EQ.'VAL') THEN
           II=INDWRTHOMO(I)
        ELSE
           II=I
        END IF
        IF (II.LT.0) THEN
           WRITE(FILENAME(I),'(A,I3.2,A,A)') FILE, II, SPN, EXT
        ELSE
           WRITE(FILENAME(I),'(A,I2.2,A,A)') FILE, II, SPN, EXT
        END IF

C
C WRITE FILE HEADER
C
       IUNIT=90+I
       OPEN(IUNIT,FILE=FILENAME(I),FORM=FORMSTR,STATUS='UNKNOWN')
       REWIND(IUNIT)
       IF (IFORM .EQ. 1) THEN
        WRITE(IUNIT) ITYPE,NSPN
        WRITE(IUNIT)(NGRID(J), J=1,3),MPBLOCK
        IF (ITYPE .EQ. 1) THEN
         WRITE(IUNIT)((RBAS(J,K), J=1,3), K=1,4)
        END IF
        WRITE(IUNIT) NWAVF
        WRITE(IUNIT)(INFOWF(J,I), J=1,4),EVLOCC(I)
       ELSE
        JSP=INFOWF(1,I)
        JRP=INFOWF(2,I)
        JST=INFOWF(3,I)
        WRITE(IUNIT,*) 'ORBITAL DENSITY FOR WF'
        WRITE(IUNIT,'(A,I1,A,I2,A,I5)')'SPIN ',JSP,
     &     ' REPRESENTATION ',JRP, ' STATE ',JST
        !<LA: Updating here to use xmol mod 
!         OPEN(77,FILE='XMOL.DAT')
!         REWIND(77)
!         READ(77,*) NATOM
!         READ(77,*)
!         WRITE(IUNIT,'(1X,I10,3F20.12)') NATOM,(RBAS(J,1),J=1,3)
         WRITE(IUNIT,'(1X,I10,3F20.12)') num_atms,(RBAS(J,1),J=1,3)
         DO K=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(K),(RBAS(J,K+1),J=1,3)
         ENDDO
!         DO K=1,NATOM
         DO K=1,num_atms
!           CHR=REAL(IZ)
           iz = xmol_list(k)%anum
           x = xmol_list(k)%rx
           y = xmol_list(k)%ry
           z = xmol_list(k)%rz
!           READ(77,*)IZ, X, Y, Z
           CHR=REAL(IZ)
C LB: WRITE COORDINATES IN ANGSTROMS
!           WRITE(IUNIT,2002)IZ, CHR, X/AU2ANG, Y/AU2ANG, Z/AU2ANG
           WRITE(IUNIT,2002)IZ, CHR, X, Y, Z
         END DO
!         CLOSE(77)
 2002    FORMAT(I6,4F16.10)
       END IF
       END DO
C
C LOOP OVER ALL POINTS IN BLOCKS
C
       ITOTWF=0
       DO IW=1,NWAVF
        K_REP=INFOWF(2,IW)
        NDM=NDMREP(K_REP)
        ITOTWF=ITOTWF+NDM
       END DO

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
             CALL WFRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &           RVECA,L_NUC,1,INFOWF,NWAVF,ITOTWF)
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
c                  IF (ICOUNT(ICON,LI)) THEN
                    FACTOR=PSI(ILOC,IWF,1)
                    DO IPV=1,NPV
                     PSIG(JOFS+IPV,IWF)=PSIG(JOFS+IPV,IWF)
     &                                 +FACTOR*GRAD(IPV,1,MU,ICON,LI)
                    END DO
c                  END IF
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
             END DO
             DO IDIM=2,NDIM
              IWF=IWF+1
              DO IPTS=1,MPTS
               PSIG(IPTS,IW)=PSIG(IPTS,IW)
     &                              +PSIG(IPTS,IWF)**2
              END DO
             END DO
           END DO
           IF(IWF.GT.ITOTWF) THEN
            WRITE(6,*) 'WFGRID: THE TOTAL NUMBER OF STATES DO NOT MATCH'
            WRITE(6,*)  IWF, ' VS ', NWF
            CALL STOPIT
           END IF
C
C WRITE OUTPUT
C
                DO IWF=1,NWAVF
                 IUNIT=90+IWF

          IF (ITYPE .EQ. 1) THEN
           IF (IFORM .EQ. 1) THEN
            WRITE(IUNIT)(PSIG(IPTS,IWF)*FACT, IPTS=1,MPTS)
           ELSE
             WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF)*FACT,
     &        IPTS=1,MPTS)
           END IF
          ELSE
           DO IPTS=1,MPTS
            WRITE(IUNIT,'(3(1X,F20.12))')(RGRID(I,IPTS), I=1,3)
            WRITE(IUNIT,'(3(1X,E20.12))')(PSIG(IPTS,IWF))
           END DO
          END IF

c                   GO TO 770
c                  END IF
c                END DO
c             END DO
c  770        CONTINUE

           END DO
  780    CONTINUE
  790   CONTINUE
  800  CONTINUE
       CLOSE(80)
       DO I=1,NWAVF
         IUNIT=90+I
         CLOSE(IUNIT)
       END DO
       NMSH=NSTORE

C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING PSIG'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING PTS'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING GRAD'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING RVECA'
       DEALLOCATE(RGRID,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING RGRID'
       DEALLOCATE(INFOWF,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING INFOWF'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'WFGRID:ERROR DEALLOCATING ICOUNT'

       RETURN
C
C ERROR HANDLING
C
  900  CLOSE(80)
  910  PRINT *,'WFGRID: ERROR READING FILE WFGRID' 
  920  CLOSE(90)
       PRINT *,'WFGRID:ERROR NOT DEALLOCATING LOCAL ARRAYS'
       NMSH=NSTORE
       RETURN
      END


********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE_old(INDEX,N,IS,IR,IB)
       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NTOTWF, I, I_POS, IBAS, IBS, IBSO, ICNT, ICON, IDIM,
     & IERR, IFNCT, IFORM, II, ILOC, INDWRTHOMO, IOFS, IPTS, IPV, IREP,
     & IRP, IRPO, ISHDUM, ISHELLA, ISIZE, ISP, ISPO, ITOTWF, ITYPE,
     & IUNIT, IW, IWF, IX, IY, IZ, J, J_POS, JJ, JOFS, JRP, JSP, JST,
     & K, K_REP, L_NUC, LI, LMAX1, M_NUC, MPTS, MU, NATOM, NDIM, NDM,
     & NPV, NSTORE, NUNSYM, NWAVF, INDEX, N, IS, IR, IB, IDEG, IREPR,
     & ISPIN, IVIRT, IWFHOMO, JVIRT, KWF, IV, NBASE, NS
       REAL*8 :: SYMBOL , CHR, FACT, FACTOR, X, Y, Z, EVL, DIFF,
     & DIFFNEW, EF, EV
      SAVE 
      DIMENSION INDEX(N),IS(N),IR(N),IB(N)
      DIMENSION ISP(MAX_OCC),IRP(MAX_OCC),IBS(MAX_OCC),EVL(MAX_OCC)
      
        WRITE(6,*)'FINDSTATE:' ,NSPN,N_REP
        IWF=0
        DO ISPIN=1,NSPN
          DO IREPR=1,N_REP
            JVIRT=0
            WRITE(6,*)'N_OCC', IREPR,ISPIN,N_OCC(IREPR,ISPIN)
            DO IVIRT=1,N_OCC(IREPR,ISPIN)
               DO IDEG=1,NDMREP(IREPR)
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EVLOCC(IWF)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
               END DO
            END DO
          END DO
        END DO
        WRITE(6,*) IWF
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J).GT.EVL(I)) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
        DO I=1,N
           KWF=IWFHOMO+INDEX(I)
           IS(I)=ISP(KWF)
           IR(I)=IRP(KWF)
           IB(I)=IBS(KWF)
        END DO
      RETURN
      END
     
********************************************************************

      SUBROUTINE FINDHOMO(IWF,EVL,IWFHOMO)
       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NTOTWF, I, I_POS, IBAS, IBS, IBSO, ICNT, ICON, IDIM,
     & IERR, IFNCT, IFORM, II, ILOC, INDWRTHOMO, IOFS, IPTS, IPV, IREP,
     & IRP, IRPO, ISHDUM, ISHELLA, ISIZE, ISP, ISPO, ITOTWF, ITYPE,
     & IUNIT, IW, IWF, IX, IY, IZ, J, J_POS, JJ, JOFS, JRP, JSP, JST,
     & K, K_REP, L_NUC, LI, LMAX1, M_NUC, MPTS, MU, NATOM, NDIM, NDM,
     & NPV, NSTORE, NUNSYM, NWAVF, INDEX, N, IS, IR, IB, IDEG, IREPR,
     & ISPIN, IVIRT, IWFHOMO, JVIRT, KWF, IV, NBASE, NS
       REAL*8 :: SYMBOL , CHR, FACT, FACTOR, X, Y, Z, EVL, DIFF,
     & DIFFNEW, EF, EV
      SAVE
      DIMENSION EVL(IWF) 
      EF=EFERMI(1)
      IF(NSPN.EQ.2) EF=MAX(EF,EFERMI(NSPN))
      DIFF=1.0D30
       DO I=1,IWF
         DIFFNEW=EF-EVL(I)
         IF(DIFFNEW.GE.-1.0D-4) THEN
           IF(DIFFNEW.LT.DIFF) THEN
            DIFF=DIFFNEW
            IWFHOMO=I
           END IF 
         END IF 
       END DO 
        WRITE(6,*)'IWFHOMO  = ', IWFHOMO,EF
      RETURN
      END

********************************************************************
C
C      THIS ROUTINE FINDS SPIN, REPRESENTATION AND THE INDEX OF THE BASIS
C      HOMO OF THE MOLECULE AND THE HOMO+INDEX STATES.
C                                          TB 04/03
C
      SUBROUTINE FINDSTATE(INDEX,N,IS,IR,IB,IERR)
       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI,
     &   NWF, NWFS, EFERMI, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: NTOTWF, I, I_POS, IBAS, IBS, IBSO, ICNT, ICON, IDIM,
     & IERR, IFNCT, IFORM, II, ILOC, INDWRTHOMO, IOFS, IPTS, IPV, IREP,
     & IRP, IRPO, ISHDUM, ISHELLA, ISIZE, ISP, ISPO, ITOTWF, ITYPE,
     & IUNIT, IW, IWF, IX, IY, IZ, J, J_POS, JJ, JOFS, JRP, JSP, JST,
     & K, K_REP, L_NUC, LI, LMAX1, M_NUC, MPTS, MU, NATOM, NDIM, NDM,
     & NPV, NSTORE, NUNSYM, NWAVF, INDEX, N, IS, IR, IB, IDEG, IREPR,
     & ISPIN, IVIRT, IWFHOMO, JVIRT, KWF, IV, NBASE, NS
       REAL*8 :: SYMBOL , CHR, FACT, FACTOR, X, Y, Z, EVL, DIFF,
     & DIFFNEW, EF, EV
      SAVE 
      LOGICAL EXIST
      CHARACTER*20 TRASH
      CHARACTER*4 caltype
      DIMENSION INDEX(N),IS(N),IR(N),IB(N),EV(MXSPN*NDH)
C      DIMENSION IOCC(N)
      INTEGER,ALLOCATABLE :: IOCC(:)
      DIMENSION ISP(MXSPN*NDH),IRP(MXSPN*NDH),
     &           IBS(MXSPN*NDH),EVL(MXSPN*NDH)
 
        IERR=0
        INQUIRE(FILE='EVALUES',EXIST=EXIST) 
        IF(EXIST) THEN
          OPEN(90,FILE='EVALUES',STATUS='OLD')
          REWIND(90)
        ELSE
          WRITE(6,*)'EVALUES DOES NOT EXIST, CALLING'
          WRITE(6,*)'OLD ROUTINE TO FIND HOMO AND OCCUPIED STATES'
          WRITE(6,*)'NEED EVALUES TO FIND LUMO'
          IERR=1
          RETURN
        END IF
        IWF=0
        READ(90,'(A4)') caltype
        write(6,'(A4)') caltype
        call flush(6)
        ALLOCATE(IOCC(N))        
        IF((caltype.eq.'OCCU').or.(caltype.eq.'occu')) then
        do ispin=1,nspn
         read(90,*)Ns
         read(90,*)(IOCC(I),I=1,Ns)
        end do
        elseif((caltype=='FIXM').or.(caltype=='fixm')) then
        else
         rewind(90)
        end if
        DEALLOCATE(IOCC)        
        DO 300 ISPIN=1,NSPN
          READ(90,240) TRASH
          DO 200 IREPR=1,N_REP
            IF(NS_TOT(IREPR).LE.0) GO TO 200
            READ(90,240) TRASH
            READ(90,*) NDIM, NBASE
             JVIRT=0
             READ(90,*) (EV(IVIRT),IVIRT=1,NBASE)
              DO IV=1,NBASE
                IWF=IWF+1
                JVIRT=JVIRT+1
                EVL(IWF)=EV(IV)
                ISP(IWF)=ISPIN
                IRP(IWF)=IREPR
                IBS(IWF)=JVIRT
              END DO
 200      CONTINUE
 300    CONTINUE
 240    FORMAT(A20)
        DO I=1,IWF
          DO J=1,IWF
           IF(EVL(J)-EVL(I).GT.1.0D-6) THEN  
             CALL SWAP(EVL(I),EVL(J))
             CALL ISWAP(ISP(I),ISP(J))
             CALL ISWAP(IRP(I),IRP(J))
             CALL ISWAP(IBS(I),IBS(J))
           END IF
          END DO
        END DO
        CLOSE(90)
        CALL FINDHOMO(IWF,EVL,IWFHOMO) 
         WRITE(6,*)'FINDSTATE : IWFHOMO', IWFHOMO
        DO I=1,N
            KWF=IWFHOMO+INDEX(I)
            IS(I)=ISP(KWF)
            IR(I)=IRP(KWF)
            IB(I)=IBS(KWF)
        END DO
      RETURN
      END

C
C******************************************************************
