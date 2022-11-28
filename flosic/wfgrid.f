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
       use mesh1,only : nmsh
       PARAMETER (NTOTWF=200)
       use common2,only : RIDT, RCNT, NCNT, N_CON, LSYMMAX, N_POS, NFNCT, ISPN, NSPN
       use common3,only : RMAT
       use common5,only : PSI_COEF, OCCUPANCY, N_OCC, PSI, NWF, NWFS, EVLOCC
       use common8,only : REP, N_REP, NDMREP, IGEN, NS_TOT
       use xmol, only: au2ang, num_atms, xmol_list
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:06 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_POS, IBAS, IBS, IBSO, ICNT, ICON, IDIM, IERR,
     & IFNCT, IFORM, II, ILOC, INDWRTHOMO, INFOWF, IOFS, IPTS, IPV,
     & IREP, IRP, IRPO, ISHDUM, ISHELLA, ISIZE, ISP, ISPN, ISPO,
     & ITOTWF, ITYPE, IUNIT, IW, IWF, IX, IY, IZ, J, J_POS, JOFS, JRP,
     & JSP, JST, K, K_REP, L_NUC, LI, LMAX1, M_NUC, MPTS, MU, NATOM,
     & NDIM, NDM, NGRID, NPV, NSPN, NSTORE, NUNSYM, NWAVF
       REAL*8 :: EV , CHR, FACT, FACTOR, GRAD, PSIG, PTS, RBAS, RGRID,
     & RVECA, X, Y, Z
       SAVE
       LOGICAL EXIST,IUPDAT,ICOUNT,WFF
       CHARACTER*20 FORMSTR
       CHARACTER*20 BAND, TRASH
       CHARACTER*6 FILE
       CHARACTER*9 FILENAME(20)
       COMMON/TMP2/PSIG(MPBLOCK,MAX_OCC)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,RVECA(3,MX_GRP),RBAS(3,4),RGRID(3,MPBLOCK)
     &  ,NGRID(3),INFOWF(4,MAX_OCC),ICOUNT(MAX_CON,3)
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
       INQUIRE(FILE='WFGRID',EXIST=EXIST)
       FORMSTR= ' '
       IF (.NOT.EXIST) FORMSTR= ' --> NOTHING TO DO'
       PRINT '(2A)','WAVEFUNCTION GRID',FORMSTR
       IF (.NOT. EXIST) RETURN

       NWAVF=21
       NSTORE=NMSH
       OPEN(80,FILE='WFGRID',FORM='FORMATTED',STATUS='OLD')
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
       GO TO 220
 210   WRITE(80,'(I3,A)') NWAVF, '    NUMBER OF ORBITALS '
 220   READ(80,240,END=310) BAND
 240   FORMAT(A20)
       GO TO 320
 310   BAND='VAL' 
       INDWRTHOMO(NWAVF)=0
       WRITE(80,'(A)') 'VAL    IF VALENCE GIVE 
     &   INDEX WITH RESPECT TO HOMO'
       WRITE(80,'(I3,2X,A)') INDWRTHOMO(NWAVF), 
     &   '   FOR CORE STATES GIVE SPIN, REPRESENATION, BASIS INDEX'
       GO TO 230
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
 230     CALL FINDSTATE(INDWRTHOMO,NWAVF,ISP,IRP,IBS,IERR)
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

        NWF=0
        DO ISPN=1,NSPN
         NWF=NWF+NWFS(ISPN)
        END DO
C
C OPEN OUTPUT FILE
C
       IF  (BAND(1:3).EQ.'VAL') THEN
          FILE='WFHOMO'
       ELSE
          FILE='WFCORE'
       END IF

       DO I=1,NWAVF
        IF(BAND(1:3).EQ.'VAL') THEN
           II=INDWRTHOMO(I)
        ELSE
           II=I
        END IF
        IF (II.LT.0) THEN
           WRITE(FILENAME(I),'(A,I3.2)') FILE, II
        ELSE
           WRITE(FILENAME(I),'(A,I2.2)') FILE, II
        END IF

c       OPEN(90,FILE='WFGROUT',FORM=FORMSTR,STATUS='UNKNOWN')
c       REWIND(90)
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
!         OPEN(77,FILE='XMOL.DAT')
!         REWIND(77)
!         READ(77,*) NATOM
!         READ(77,*)
         WRITE(IUNIT,'(1X,I10,3F20.12)') num_atms,(RBAS(J,1),J=1,3)
         DO K=1,3
         WRITE(IUNIT,'(1X,I10,3F20.12)') NGRID(K),(RBAS(J,K+1),J=1,3)
         ENDDO
         DO k=1,num_atms
!           READ(77,*)IZ, X, Y, Z
           iz = xmol_list(k)%anum
           x = xmol_list(k)%rx
           y = xmol_list(k)%ry
           z = xmol_list(k)%rz
           CHR=REAL(IZ)
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
              PSIG(IPTS,IW)=0.0d0
             END DO
             WRITE(6,*)'IW =', IW, 'NDIM=',NDIM
             DO IDIM=2,NDIM
              IWF=IWF+1
              if (IDIM.EQ.3) THEN
              DO IPTS=1,MPTS
               PSIG(IPTS,IW)=PSIG(IPTS,IW)
     &                              +PSIG(IPTS,IWF) 
              END DO
              end if
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
       RETURN
C
C ERROR HANDLING
C
  900  CLOSE(80)
  910  PRINT *,'WFGRID: ERROR READING FILE WFGRID' 
  920  CLOSE(90)
       NMSH=NSTORE
       RETURN
      END
