C UTEP Electronic Structure Lab (2020)
C
C *********************************************************************
C
      SUBROUTINE FRCNONL
C
C DIRK POREZAG, AUGUST 1997
C CALCULATE NONLOCAL PSP CONTRIBUTION TO THE FORCES AND NONLOCAL ENERGY
C ATTENTION: THIS ROUTINE WILL ONLY WORK CORRECTLY IF THE IDENTITY
C MATRIX IS THE FIRST ELEMENT OF RMAT
C
C THEORY: The total energy contribution due to the nonlocal PSP is:
C
C        occ rad  2 lmx  A    / ang      _   ___         \ 2
C E= SUM SUM INT r  SUM V (r)(  INT PSI (R + rom) Y  (om) )
C     A   i  dr     lm   l    \ dom    i  A        lm    / 
C       ___
C where rom is the vector based on radius r and angular coordinate om.
C The gradient of this expression can be written as:
C
C dE            A     B
C -- = 2 SUM ( g   - g   ) (k= x,y,z), where
C dR      B     Bk    Ak
C   Ak
C
C  A   rad  2 lmx  A       B    occ                        B
C g  = INT r  SUM V (r)   SUM   SUM c     PSILM     GRADKLM
C  Bk  dr     lm   l    mu of B  i   mu_i      i_lm        k_mu_lm
C
C                   B          ang    B  _   _    ___
C PSILM    = SUM   SUM   c     INT PHI  (R - R  + rom) Y  (om)
C      i_lm   B  mu of B  mu_i dom    mu  A   B         lm
C
C        B        ang       B  _   _    ___
C GRADKLM       = INT (d PHI  (R - R  + rom))/(dk) Y  (om)
C        k mu lm  dom       mu  A   B               lm
C
C                              B
C where all basis functions PHI  are centered on atom B.
C                              mu
C
       use common1,only : PSRZONE, LMXPSRZ, PSPSYM, ISITPSP, RPSNLO,
     &    WPSNLO, VPSNLO, LMAXNLO, NRPSP
       use common2,only : RIDT, IFUIDT, NIDENT, N_CON, LSYMMAX, N_POS,
     &   NFNCT, NSPN, ENONLO, FNONL
       use common3,only : RMAT, NGRP
       use common5,only : PSI, NWF
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:45 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IANG, IATOM, IATOMS, IBAS, IBG, IBSSAV, ICON,
     & IDUM, IERR, IFACSP, IFNCT, IFUPSP, IGRP, IID, IIDAT, IIDPSP,
     & ILOC, INDX, INUC, IOFS, IPOS, IPTS, IRAD, IWF, IYLM, IZONE, J,
     & JID, JIDAT, JNUC, K, KNUC, L, L1, LMAX1, LMDIM, LMDMAX, LMULT,
     & LPSPMX, M, MAXANG, MAXBAS, MGRP, MNUC, MPTS, MU, MXYLM, NATOMS,
     & NBAS
       REAL*8 :: SYMBOL , DIF, FACT, FDIFF, RRAD, SMALL, SUM, TIME1,
     & TIME2, ZERO
       SAVE
       PARAMETER(MAXBAS=MAXUNSYM*MX_CNT)
       PARAMETER(MAXANG=300)
       PARAMETER(MXYLM=(MXLPSP+1)**2)
C       COMMON/TMP2/ANGLE(3,MAXANG,5),DOMEGA(MAXANG)
C     &  ,YLMP(MAXANG,MXYLM,5),RNUC(3,MX_GRP),VNLSAV(MXYLM)
C     &  ,RTAB(3,MX_CNT),PTS(NSPEED,3),GRDKBA(3,MX_CNT,MAX_IDENT)
C     &  ,GRAD(NSPEED,10,6,MAX_CON,3),RVEC(3),FVEC(3)
C     &  ,PSILM(MAX_OCC,MXYLM),PSISAV(NSPEED,MAX_OCC)
C     &  ,GRADKLM(3,MAXBAS,MXYLM)
C     &  ,NMULTI(MAX_IDENT),NANG(5),NYLM(5)

       REAL*8,ALLOCATABLE :: ANGLE(:,:,:),DOMEGA(:),YLMP(:,:,:)
     &  ,RNUC(:,:),VNLSAV(:),RTAB(:,:),PTS(:,:),GRDKBA(:,:,:)
     &  ,GRAD(:,:,:,:,:),PSILM(:,:),PSISAV(:,:),GRADKLM(:,:,:)
       INTEGER,ALLOCATABLE :: NMULTI(:)
       REAL*8 :: RVEC(3),FVEC(3)
       INTEGER :: NANG(5),NYLM(5)
       LOGICAL ICOUNT,IUPDAT,SKIPIT
       DIMENSION LMULT(3)
       DIMENSION ICOUNT(MAX_CON,3)
       DIMENSION SKIPIT(MX_CNT)
       DATA LMULT/1,3,6/
       DATA SMALL/1.0D-3/
       DATA ZERO /1.0D-4/
C
       CALL GTTIME(TIME1)
C
C ALLOCATE LOCAL ARRAYS
C
       ALLOCATE(ANGLE(3,MAXANG,5),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING ANGLE'
       ALLOCATE(DOMEGA(MAXANG),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING DOMEGA'
       ALLOCATE(YLMP(MAXANG,MXYLM,5),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING YLMP'
       ALLOCATE(RNUC(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING RNUC'
       ALLOCATE(VNLSAV(MXYLM),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING VNLSAV'
       ALLOCATE(RTAB(3,MX_CNT),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING RTAB'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING PTS'
       ALLOCATE(GRDKBA(3,MX_CNT,MAX_IDENT),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING GRDKBA'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING GRAD'
       ALLOCATE(PSILM(MAX_OCC,MXYLM),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING PSILM'
       ALLOCATE(PSISAV(NSPEED,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING PSISAV'
       ALLOCATE(GRADKLM(3,MAXBAS,MXYLM),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING GRADKLM'
       ALLOCATE(NMULTI(MAX_IDENT),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR ALLOCATING NMULTI'
C
C INITIALIZE ENONLO AND THE ARRAY WITH NONLOCAL FORCES
C DEFINE SPIN PREFACTOR
C LEAVE IF THIS IS NOT A PSP CALCULATION 
C
       ENONLO= 0.0D0
       IFACSP= 2/NSPN
       DO IID=1,NIDENT
        FNONL(1,IID)= 0.0D0
        FNONL(2,IID)= 0.0D0
        FNONL(3,IID)= 0.0D0
       END DO
       IF (ISITPSP .NE. 1) GOTO 900
C
C DETERMINE TOTAL NUMBER OF ATOMS AND BASIS FUNCTIONS
C SAVE MULTIPLICITY OF IDENTITY MEMBERS IN NMULTI
c SAVE COORDINATES OF ALL ATOMS IN RTAB
C
       NATOMS=0
       NBAS=0
       DO IID=1,NIDENT
        CALL OBINFO(1,RIDT(1,IID),RNUC,MNUC,IDUM)
        NMULTI(IID)=MNUC
        IF (NATOMS+MNUC .GT. MX_CNT) THEN
         PRINT *,'FRCNONL: MX_CNT MUST BE AT LEAST: ',NATOMS+MNUC
         CALL STOPIT
        END IF
        DO INUC=1,MNUC
         RTAB(1,NATOMS+INUC)= RNUC(1,INUC) 
         RTAB(2,NATOMS+INUC)= RNUC(2,INUC) 
         RTAB(3,NATOMS+INUC)= RNUC(3,INUC) 
        END DO
        NATOMS=NATOMS+MNUC
        IFNCT=IFUIDT(IID)
        DO L1=1,LSYMMAX(IFNCT)+1
         NBAS=NBAS+MNUC*N_CON(L1,IFNCT)*LMULT(L1)
        END DO
       END DO
       IF (NBAS .GT. MAXBAS) THEN
        PRINT *,'FRCNONL: MAXBAS MUST BE AT LEAST: ',NBAS
        CALL STOPIT
       END IF
C
C DEFINE LPSPMX (MAXIMUM L FOR PSEUDO-POTENTIAL)
C
       LPSPMX=0
       DO IIDPSP=1,NIDENT
        IFUPSP=IFUIDT(IIDPSP)
        IF (PSPSYM(IFUPSP)(1:3) .NE. 'ALL') THEN
         LPSPMX=MAX(LPSPMX,LMAXNLO(IFUPSP))
        END IF
       END DO
C
C DEFINE YLMP: CONTAINS (SPHERICAL HARMONICS)*(ANGULAR WEIGHT)
C
       MPTS=0
       CALL HARMONICS(MAXANG,MPTS,LPSPMX,ANGLE(1,1,1),
     &                YLMP(1,1,1),LMDMAX)
       IF (LMDMAX .GT. MXYLM) THEN
        PRINT *,'VNLHAM: MXYLM MUST BE AT LEAST: ',LMDMAX
        CALL STOPIT
       END IF
       DO IZONE=1,5
        CALL ANGMSH(MAXANG,LMXPSRZ(IZONE),NANG(IZONE),
     &              ANGLE(1,1,IZONE),DOMEGA)
        CALL HARMONICS(MAXANG,NANG(IZONE),LPSPMX,ANGLE(1,1,IZONE),
     &                 YLMP(1,1,IZONE),NYLM(IZONE))
        DO IYLM=1,NYLM(IZONE)
         DO IANG=1,NANG(IZONE)
          YLMP(IANG,IYLM,IZONE)= YLMP(IANG,IYLM,IZONE)*DOMEGA(IANG)
         END DO
        END DO
       END DO
C
C HERE STARTS THE BIG LOOP
C INITIALIZE GRDKBA
C SET UP RADIUS, TYPE OF ANGULAR MESH, AND NONLOCAL POTENTIAL
C
       DO 500 IIDPSP=1,NIDENT 
        DO IATOMS=1,NATOMS
         GRDKBA(1,IATOMS,IIDPSP)= 0.0D0
         GRDKBA(2,IATOMS,IIDPSP)= 0.0D0
         GRDKBA(3,IATOMS,IIDPSP)= 0.0D0
        END DO
        IFUPSP=IFUIDT(IIDPSP)
        IF (PSPSYM(IFUPSP)(1:3) .EQ. 'ALL') GOTO 500
        DO 490 IRAD=1,NRPSP(IFUPSP)
         RRAD= RPSNLO(IRAD,IFUPSP)
         IZONE=5
         DO I= 4, 1, -1
          IF (RRAD .LT. PSRZONE(I)) THEN
           IZONE=I
          ELSE
           GOTO 20
          END IF
         END DO
   20    CONTINUE
C
C SAVE NONLOCAL POTENTIAL FOR LATER
C
         LMDIM=(LMAXNLO(IFUPSP)+1)**2
         DO L=0,LMAXNLO(IFUPSP)
          IBG=L*L
          DO M=1,2*L+1
           VNLSAV(IBG+M)=VPSNLO(L+1,IRAD,IFUPSP)*WPSNLO(IRAD,IFUPSP)
          END DO
         END DO
C
C INITIALIZE PSILM AND GRADKLM 
C
         DO IYLM=1,LMDIM
          DO IWF=1,NWF
           PSILM(IWF,IYLM)= 0.0D0
          END DO
          DO IBAS=1,NBAS
           GRADKLM(1,IBAS,IYLM)= 0.0D0
           GRADKLM(2,IBAS,IYLM)= 0.0D0
           GRADKLM(3,IBAS,IYLM)= 0.0D0
          END DO
         END DO
C
C LOOP ONE: CALCULATE ANGULAR PROJECTION OF PSI AND THE 
C BASIS FUNCTION DERIVATIVES
C
         IID=0
         IATOM=0
         IBAS=0
         DO 200 IFNCT=1,NFNCT
          LMAX1=LSYMMAX(IFNCT)+1
          DO IPOS=1,N_POS(IFNCT)
           IID=IID+1
           CALL OBINFO(1,RIDT(1,IID),RNUC,MNUC,IDUM)
           DO 150 INUC=1,MNUC
            IATOM=IATOM+1
            SKIPIT(IATOM)= .TRUE.
            CALL UNRAVEL(IFNCT,IID,INUC,RIDT(1,IID),
     &                   RNUC,KNUC,1)
            IF (KNUC .NE. MNUC) THEN
             PRINT *,'APTSLV: PROBLEM IN UNRAVEL'
             CALL STOPIT
            END IF
C
C SUM OVER ALL ANGULAR POINTS 
C
            RVEC(1)= RIDT(1,IIDPSP)-RNUC(1,INUC)
            RVEC(2)= RIDT(2,IIDPSP)-RNUC(2,INUC)
            RVEC(3)= RIDT(3,IIDPSP)-RNUC(3,INUC)
            IBSSAV=IBAS
            DO IOFS=0,NANG(IZONE)-1,NSPEED
             IBAS=IBSSAV
             MPTS=MIN(NSPEED,NANG(IZONE)-IOFS)
             DO I=1,3
              DO IPTS=1,MPTS
               PTS(IPTS,I)= RVEC(I)+RRAD*ANGLE(I,IOFS+IPTS,IZONE)
              END DO
             END DO
             CALL GORBDRV(1,IUPDAT,ICOUNT,MPTS,PTS,IFNCT,GRAD)
             IF (IUPDAT) THEN
              SKIPIT(IATOM)= .FALSE.
              DO IWF=1,NWF
               DO IPTS=1,MPTS
                PSISAV(IPTS,IWF)= 0.0D0
               END DO
              END DO
              ILOC=0
              DO L1=1,LMAX1
               DO MU=1,LMULT(L1)
                DO ICON=1,N_CON(L1,IFNCT)
                 ILOC=ILOC+1
                 IBAS=IBAS+1
                 IF (ICOUNT(ICON,L1)) THEN
                  DO IWF=1,NWF
                   FACT=PSI(ILOC,IWF,1) 
                   DO IPTS=1,MPTS
                    PSISAV(IPTS,IWF)= PSISAV(IPTS,IWF)
     &              +FACT*GRAD(IPTS,1,MU,ICON,L1)
                   END DO
                  END DO
                  DO IYLM=1,LMDIM
                   DO K=1,3
                    SUM= 0.0D0
                    DO IPTS=1,MPTS
                     SUM= SUM+GRAD(IPTS,K+1,MU,ICON,L1)
     &                    *YLMP(IOFS+IPTS,IYLM,IZONE)
                    END DO
                    GRADKLM(K,IBAS,IYLM)= GRADKLM(K,IBAS,IYLM)+SUM
                   END DO
                  END DO    ! IYLM
                 END IF
                END DO      ! ICON
               END DO       ! MU
              END DO        ! L1
              DO IYLM=1,LMDIM
               DO IWF=1,NWF
                SUM= 0.0D0
                DO IPTS=1,MPTS
                 SUM= SUM+PSISAV(IPTS,IWF)*YLMP(IOFS+IPTS,IYLM,IZONE)
                END DO
                PSILM(IWF,IYLM)= PSILM(IWF,IYLM)+SUM
               END DO
              END DO
             ELSE
              DO L1=1,LMAX1
               IBAS=IBAS+N_CON(L1,IFNCT)*LMULT(L1)
              END DO
             END IF
            END DO          ! IOFS 
  150      CONTINUE         ! INUC
          END DO            ! IPOS
  200    CONTINUE           ! IFNCT
C
C UPDATE NONLOCAL ENERGY
C
         DO IYLM=1,LMDIM
          SUM= 0.0D0
          DO IWF=1,NWF
           SUM= SUM+PSILM(IWF,IYLM)*PSILM(IWF,IYLM)
          END DO
          ENONLO= ENONLO+NMULTI(IIDPSP)*VNLSAV(IYLM)*SUM
         END DO
C
C LOOP TWO: UPDATE GRDKBA
C
         IID=0
         IATOM=0
         IBAS=0
         DO 300 IFNCT=1,NFNCT
          LMAX1=LSYMMAX(IFNCT)+1
          DO IPOS=1,N_POS(IFNCT)
           IID=IID+1
           CALL OBINFO(1,RIDT(1,IID),RNUC,MNUC,IDUM)
           DO 250 INUC=1,MNUC
            IATOM=IATOM+1
            IF (SKIPIT(IATOM)) THEN
             DO L1=1,LMAX1
              IBAS=IBAS+N_CON(L1,IFNCT)*LMULT(L1)
             END DO
             GOTO 250
            END IF
            CALL UNRAVEL(IFNCT,IID,INUC,RIDT(1,IID),
     &                 RNUC,KNUC,1)
            ILOC=0
            DO L1=1,LMAX1
             DO MU=1,LMULT(L1)
              DO ICON=1,N_CON(L1,IFNCT)
               ILOC=ILOC+1
               IBAS=IBAS+1
               DO IYLM=1,LMDIM
                SUM= 0.0D0
                DO IWF=1,NWF
                 SUM= SUM+PSI(ILOC,IWF,1)*PSILM(IWF,IYLM)
                END DO
                DO K=1,3
                 GRDKBA(K,IATOM,IIDPSP)= GRDKBA(K,IATOM,IIDPSP)
     &           +SUM*GRADKLM(K,IBAS,IYLM)*VNLSAV(IYLM) 
                END DO
               END DO       ! IYLM
              END DO        ! ICON
             END DO         ! MU
            END DO          ! L1
  250      CONTINUE         ! INUC
          END DO            ! IPOS
  300    CONTINUE           ! IFNCT
  490   CONTINUE            ! IRAD
  500  CONTINUE             ! IIDPSP
C
C NOW THAT WE HAVE GRDKBA, THE FORCES ARE NOT A BIG DEAL ANY MORE
C THE ONLY THING WE STILL NEED TO DEAL WITH IS THE SYMMETRY STUFF
C FIRST, THE SIMPLE PART
C
       IIDAT=1
       DO IID=1,NIDENT
        DO IATOMS= 1,NATOMS 
         FNONL(1,IID)= FNONL(1,IID)+GRDKBA(1,IATOMS,IID)
         FNONL(2,IID)= FNONL(2,IID)+GRDKBA(2,IATOMS,IID)
         FNONL(3,IID)= FNONL(3,IID)+GRDKBA(3,IATOMS,IID)
        END DO
C
C PART TWO IS A LITTLE MORE TRICKY: WE NEED TO ROTATE OUR DATA 
C
        JIDAT=1
        DO JID=1,NIDENT
         DO JNUC=1,NMULTI(JID)
C
C SUM OVER ALL SYMMETRY OPERATIONS THAT CREATE THIS SITE
C
          FVEC(1)= 0.0D0   
          FVEC(2)= 0.0D0   
          FVEC(3)= 0.0D0   
          MGRP=0
          DO IGRP= 1,NGRP
           DO I=1,3
            RVEC(I)= 0.0D0
            DO J=1,3
             RVEC(I)= RVEC(I)+RMAT(J,I,IGRP)*RIDT(J,JID)
            END DO
           END DO
           DIF= ABS(RVEC(1)-RTAB(1,JIDAT+JNUC-1))
     &         +ABS(RVEC(2)-RTAB(2,JIDAT+JNUC-1))
     &         +ABS(RVEC(3)-RTAB(3,JIDAT+JNUC-1))
           IF (DIF .LE. ZERO) THEN
            MGRP=MGRP+1
C
C LOOK FOR THE ATOM THAT CREATES IDENTITY MEMBER IID WHEN USING
C SYMMETRY OPERATION IGRP
C
            INDX=0
            DO INUC=1,NMULTI(IID)
             DO I=1,3
              RVEC(I)= 0.0D0
              DO J=1,3
               RVEC(I)= RVEC(I)+RMAT(J,I,IGRP)*RTAB(J,IIDAT+INUC-1)
              END DO
             END DO
             DIF= ABS(RVEC(1)-RIDT(1,IID))
     &           +ABS(RVEC(2)-RIDT(2,IID))
     &           +ABS(RVEC(3)-RIDT(3,IID))
             IF (DIF .LE. ZERO) INDX=IIDAT+INUC-1
            END DO
            IF (INDX .EQ. 0) THEN
             PRINT *,'FRCNONL: ERROR 1 WHILE ROTATING FORCES'
             CALL STOPIT
            END IF
C
C ROTATE GRDKBA AND ADD TO FVEC
C
            DO I=1,3
             SUM= 0.0D0
             DO J=1,3
              SUM= SUM+RMAT(J,I,IGRP)*GRDKBA(J,INDX,JID)
             END DO
             FVEC(I)= FVEC(I)+SUM
            END DO
           END IF
          END DO  
          IF (MGRP .EQ. 0) THEN
           PRINT *,'FRCNONL: ERROR 2 WHILE ROTATING FORCES'
           CALL STOPIT
          END IF
          FNONL(1,IID)= FNONL(1,IID)-FVEC(1)/MGRP
          FNONL(2,IID)= FNONL(2,IID)-FVEC(2)/MGRP
          FNONL(3,IID)= FNONL(3,IID)-FVEC(3)/MGRP
         END DO
         JIDAT=JIDAT+NMULTI(JID)
        END DO
        IIDAT=IIDAT+NMULTI(IID)
       END DO
       DO IID=1,NIDENT
        FNONL(1,IID)= -2*IFACSP*FNONL(1,IID)
        FNONL(2,IID)= -2*IFACSP*FNONL(2,IID)
        FNONL(3,IID)= -2*IFACSP*FNONL(3,IID)
        CALL FRCSYM(RIDT(1,IID),FNONL(1,IID),FDIFF)
        IF (FDIFF .GT. SMALL) THEN
         PRINT 1010,IID,FDIFF
 1010    FORMAT(' WARNING: NONLOCAL FORCE OF ATOM ',I3,
     &          ' VIOLATES SYMMETRY BY ',D12.4)
        END IF
       END DO
       ENONLO= ENONLO*IFACSP
C
C THE END - DEPENDS ON YOU WHETHER IT IS A GOOD ONE :-)
C

C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(ANGLE,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING ANGLE'
       DEALLOCATE(DOMEGA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING DOMEGA'
       DEALLOCATE(YLMP,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING YLMP'
       DEALLOCATE(RNUC,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING RNUC'
       DEALLOCATE(VNLSAV,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING VNLSAV'
       DEALLOCATE(RTAB,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING RTAB'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING PTS'
       DEALLOCATE(GRDKBA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING GRDKBA'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING GRAD'
       DEALLOCATE(PSILM,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING PSILM'
       DEALLOCATE(PSISAV,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING PSISAV'
       DEALLOCATE(GRADKLM,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING GRADKLM'
       DEALLOCATE(NMULTI,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'FRCNONL:ERROR DEALLOCATING NMULTI'
C
  900  CALL GTTIME(TIME2)
       CALL TIMOUT('NONLOCAL CONTRIBUTION TO FORCES:   ',TIME2-TIME1)
 1020  FORMAT(1X,A,F12.3)
       RETURN
      END
