C UTEP Electronic Structure Lab (2020)
       SUBROUTINE DOSJNT_ORIG 
C
C DOSJNT BY MRP/DVP NOV 1996
C CALCULATES THE JOINT DENSITY OF STATES WITH AND WITHOUT DIPOLE WEIGHTS
C
       use debug1
       use dosjnt_mod,only : H,PSIG,SPTOT,SPDIP,SOS_FREQ,
     &   RVECA,PTS,GRAD,ICOUNT,ESTEP,EALP,SOS_POL,VFAC,
     &   ENJD,EXJD,TEMP,HA2EV,NSPEC,ISIZE,P,Q,V,FCGRP,CHARGE
       use global_inputs,only : dosjnt1
       use mesh1,only : wmsh,rmsh,nmsh
       use common2,only : RIDT, N_CON, LSYMMAX, N_POS, NFNCT, E_UP,
     &   E_DN, ISPN, NSPN, DIPOLE
       use common3,only : RMAT, NGRP
       use common5,only : PSI_COEF, PSI, NWF, NWFS, EFERMI, EVLOCC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:42 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I_POS, ICON, IERR, IFAIL, IFNCT, IGP, ILOC, IPILE,
     & IPTS, ISHDUMMY, ISHELLA, ISWITCH, IWF, IX, J, J_POS, JPTS, JSPN,
     & JWF, JWMAX, K, L_NUC, LI, LMAX1, LMOM, LPV, M_NUC, MPTS, MU,
     & MXSPEC, NMAX, NOFS, NPILE, NPV
       REAL*8 :: SYMBOL , BACKW, DIP, EDIFF, EDIFF27, EFMAX, EFMIN,
     & EFRMI, EMAX, EMIN, ERG, EV1, EV2, EWIND, FACT, FACTOR, FCT,
     & FUNC, FWHM, OCCI, OCCJ, PI, SDMAX, STMAX, TIME1, TIME2, TIMEA,
     & TIMEB, TIMEC, TIMEOLD
       SAVE
       PARAMETER (MXSPEC=10000)
       PARAMETER (NMAX=MPBLOCK)
C
C FOOL THE COMPILER FOR MXSPN=1 TO SUPRESS WARNING MESSAGES
C THAT ARE REALLY IRRELEVANT
C
       LOGICAL IUPDAT,EXIST,LMKFIL,USEMPI,DMOM
C       COMMON/TMP2/PSIG(NMAX,MAX_OCC)
C
C SCRATCH COMMON BLOCK FOR LOCAL ARRAYS
C
C       COMMON/TMP1/H(MAX_OCC,MAX_OCC,3)
C     &  ,SPDIP(MXSPEC),SPTOT(MXSPEC),RVECA(3,MX_GRP)
C     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
C     &  ,ICOUNT(MAX_CON,3)
C       REAL*8,ALLOCATABLE :: H(:,:,:),SPDIP(:),SPTOT(:),RVECA(:,:),
C     &                       PTS(:,:),GRAD(:,:,:,:,:)
C       REAL*8,ALLOCATABLE :: RVECA(:,:),PTS(:,:),GRAD(:,:,:,:,:)
C       LOGICAL,ALLOCATABLE :: ICOUNT(:,:)
C       DIMENSION ISIZE(3)
C       DIMENSION P(NMAX,3),Q(NMAX,3),V(NMAX)
       DATA USEMPI/.FALSE./
C       DATA ISIZE/1,3,6/
C       DATA HA2EV/27.2116D0/
C       DATA TEMP/1.0D-4/
       DATA EMIN,EMAX/-10000.,1.0D0/
C       DATA ENJD,EXJD/ 0.0,1.0D0/
       DATA FWHM/0.05D0/
C
C RETURN IF INPUT FILE DOES NOT EXIST
C
       PRINT '(A)','CALCULATING JOINT DENSITY OF STATES'
       CALL CHECK_INPUTS
       IF (.NOT.DOSJNT1) THEN
        PRINT '(2A)','DOSJNT: CHECK NRLMOL_INPUT.DAT ',
     &               '--> NOTHING TO DO'
        RETURN
       END IF
C
C ALLOCATE LOCAL ARRAYS
C
       ALLOCATE(H(MAX_OCC,MAX_OCC,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING H'
       ALLOCATE(SPDIP(MXSPEC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING SPDIP'
       ALLOCATE(SPTOT(MXSPEC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING SPTOT'
       ALLOCATE(SOS_FREQ(MXSPEC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING SOS_FREQ'
       ALLOCATE(RVECA(3,MX_GRP),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING RVECA'
       ALLOCATE(PTS(NSPEED,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING PTS'
       ALLOCATE(GRAD(NSPEED,10,6,MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING GRAD'
       ALLOCATE(ICOUNT(MAX_CON,3),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING ICOUNT'
       ALLOCATE(PSIG(NMAX,MAX_OCC),STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR ALLOCATING PSIG'
       ENJD=0.0
       EXJD=1.0D0
       TEMP=1.0D-4
C
C CREATE A STANDARD INPUT FILE IF THE CURRENT INPUT FILE IS EMPTY
C
C       CALL GTTIME(TIME1)
       LMKFIL=.TRUE.
C       OPEN(74,FILE='DOSJNT',FORM='FORMATTED',STATUS='OLD')
       OPEN(74,FILE='DOSJNT',FORM='FORMATTED')
       REWIND(74)
       READ(74,*,END=5,ERR=5) ISWITCH
       IF (ISWITCH .NE. 0) LMKFIL=.FALSE.
    5  CLOSE(74)
C
       IF (LMKFIL) THEN
C        OPEN(74,FILE='DOSJNT',FORM='FORMATTED',STATUS='OLD')
        OPEN(74,FILE='DOSJNT',FORM='FORMATTED')
        REWIND(74)
        WRITE(74,*) '0  auto=0, otherwise user-defined'
        WRITE(74,99) EMIN,EMAX,FWHM,'EMIN, EMAX, FWHM IN HARTREES'
        WRITE(74,98) ENJD,EXJD,'PLOTTED ENERGY WINDOW IN HARTREES'
        CLOSE(74)
   98  FORMAT(2F14.5,16X,A)
   99  FORMAT(3F14.5,2X,A)
       END IF
C
C READ IN TEMPERATURE
C
       OPEN(39,FILE='TMPTRE',FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(39)
       READ(39,*,END=10) TEMP
   10  CLOSE(39)
C
C READ IN NECESSARY INPUT DATA
C
C       OPEN(74,FILE='DOSJNT',FORM='FORMATTED',STATUS='OLD')
       OPEN(74,FILE='DOSJNT',FORM='FORMATTED')
       REWIND(74)
       READ(74,*,END=20) ISWITCH
       READ(74,*,END=20) EMIN,EMAX,FWHM
       READ(74,*,END=20) ENJD,EXJD
       GOTO 25
   20  CLOSE(74)
       PRINT *,'DOSJNT: FILE DOSJNT IS BROKEN'
       GOTO 900
C
   25  CLOSE(74)
       EWIND=(EXJD-ENJD)
       IF ((EWIND .LE. 0.0D0) .OR. (FWHM .LE. 0.0D0)) THEN
        PRINT '(A)','CALCULATION OF JDOS HAS BEEN SKIPPED'
        GOTO 900
       END IF
       NSPEC=INT(10*EWIND/FWHM)+2
       EALP=4*LOG(2.0D0)/FWHM**2
       PI=4*ATAN(1.0D0)
       VFAC=SQRT(EALP/PI)/HA2EV
       IF (NSPEC.GT.MXSPEC) THEN
        PRINT *,'DOSJNT: MXSPEC MUST BE AT LEAST: ',NSPEC
        CALL STOPIT
       END IF
      IF(.NOT.USEMPI)THEN
C
C CALL WFWIND TO GET THE CORRECT PSI_COEF
C WFWIND DEFINES NWF AND NWFS
C
       CALL GTTIME(TIME1)
       TIMEOLD=TIME1
       EFMIN=MIN(EFERMI(1),EFERMI(NSPN))
       EFMAX=MAX(EFERMI(1),EFERMI(NSPN))
C      EMIN=EFMIN-EWIND-(4*FWHM+20*TEMP)
C      EMAX=EFMAX+EWIND+(4*FWHM+20*TEMP)
       CALL WFWIND(EMIN,EMAX,.TRUE.,.TRUE.,IFAIL)
       CALL TRACER('NWF=',NWF)
       IF (IFAIL .EQ. 1) THEN
        PRINT *,'DOSJNT: WFWIND FAILED, ABORTING JDOS CALCULATION'
        RETURN
       END IF
C
C ZERO H (CONTAINS THE DIPOLE MATRIX ELEMENTS)
C
       DO IX=1,3
        DO IWF=1,NWF
         DO JWF=1,NWF
          H(JWF,IWF,IX)=0.0D0
         END DO
        END DO
       END DO
C
C CALCULATE DIPOLE MATRIX ELEMENTS BY MESH INTEGRATION
C
       CHARGE=0.0D0
       NPILE=NMSH/NMAX
       FCGRP=1.0D0/NGRP
       DO 850 IPILE=0,NPILE
        NOFS=IPILE*NMAX
        MPTS=MIN(NMAX,NMSH-NOFS)
        CALL TRACER('MPTS',MPTS)
        DO IPTS=1,MPTS
         Q(IPTS,1)=RMSH(1,IPTS+NOFS)
         Q(IPTS,2)=RMSH(2,IPTS+NOFS)
         Q(IPTS,3)=RMSH(3,IPTS+NOFS)
         V(  IPTS)=WMSH(IPTS+NOFS)*FCGRP
        END DO
        CALL TRACER('NGRP=',NGRP)
        DO 800 IGP=1,NGRP
         DO J=1,3
          DO IPTS=1,MPTS
           P(IPTS,J)=0.0D0
          END DO
          DO K=1,3
           DO IPTS=1,MPTS
            P(IPTS,J)=P(IPTS,J)+RMAT(K,J,IGP)*Q(IPTS,K)
           END DO
          END DO
         END DO
C
C INITIALIZE PSIG 
C
         DO IWF=1,NWF
          DO IPTS=1,MPTS
           PSIG(IPTS,IWF)=0.0D0
          END DO
         END DO  
         ISHELLA=0
         DO 86 IFNCT=1,NFNCT
          LMAX1=LSYMMAX(IFNCT)+1
          DO 84 I_POS=1,N_POS(IFNCT)
           ISHELLA=ISHELLA+1
           CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
           DO 82 J_POS=1,M_NUC
            CALL UNRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                   RVECA,L_NUC,1)
            IF(L_NUC.NE.M_NUC)THEN
             PRINT *,'DOSJNT: PROBLEM IN UNRAVEL'
             CALL STOPIT
            END IF
            DO 80 JPTS=1,MPTS,NSPEED
             NPV=MIN(NSPEED,MPTS-JPTS+1)
             IPTS=JPTS-1
             DO LPV=1,NPV
              PTS(LPV,1)=P(IPTS+LPV,1)-RVECA(1,J_POS)
              PTS(LPV,2)=P(IPTS+LPV,2)-RVECA(2,J_POS)
              PTS(LPV,3)=P(IPTS+LPV,3)-RVECA(3,J_POS)
             END DO
             CALL GORBDRV(0,IUPDAT,ICOUNT,NPV,PTS,IFNCT,GRAD)
             IF (IUPDAT) THEN
              ILOC=0
              DO 78 LI=1,LMAX1
               DO MU=1,ISIZE(LI)
                DO ICON=1,N_CON(LI,IFNCT)
                 ILOC=ILOC+1
                 IF (ICOUNT(ICON,LI)) THEN
                  DO IWF=1,NWF
                   FACTOR=PSI(ILOC,IWF,1)
                   DO LPV=1,NPV
                    PSIG(IPTS+LPV,IWF)=PSIG(IPTS+LPV,IWF)
     &              +FACTOR*GRAD(LPV,1,MU,ICON,LI)
                   END DO  
                  END DO  
                 END IF
                END DO  
               END DO  
   78         CONTINUE
             END IF
   80       CONTINUE
   82      CONTINUE
   84     CONTINUE
   86    CONTINUE
        CALL GTTIME(TIMEA)
        CALL TIMOUT('DOSJNT Marker 1                   :',TIMEA-TIMEOLD)
         TIMEOLD=TIMEA
C
C UPDATE CHARGE AND DIPOLE MATRICES
C
         DO IWF=1,NWF
          DO IPTS=1,MPTS
           CHARGE=CHARGE+V(IPTS)*PSIG(IPTS,IWF)**2
          END DO
         END DO
         DO IX=1,3
          DO IWF=1,NWF
           JWMAX=NWFS(1)
           IF (IWF.GT.NWFS(1)) JWMAX=NWF
           DO JWF=IWF,JWMAX
            DO IPTS=1,MPTS
             H(JWF,IWF,IX)=H(JWF,IWF,IX)
     &       +PSIG(IPTS,IWF)*PSIG(IPTS,JWF)*V(IPTS)*P(IPTS,IX)
            END DO
           END DO
          END DO
         END DO
        CALL GTTIME(TIMEC)
        CALL TIMOUT('DOSJNT Marker 2                   :',TIMEC-TIMEOLD)
  800   CONTINUE
  850  CONTINUE
         CALL TRACER('CHARGE=',1,CHARGE)
         CALL GTTIME(TIMEB)
         CALL TIMOUT('DOSJNT Marker 3                   :',TIMEB-TIME1)
       ELSE
       write(6,*)'READING DIPOLE MATRIX ELEMENTS FROM SPNDAT'
           INQUIRE(FILE='SPNDAT',EXIST=EXIST)
           IF(EXIST)THEN
           OPEN(74,FILE='SPNDAT',FORM='UNFORMATTED')
           READ (74)NWF,NWFS,NSPN,DMOM,LMOM
           READ (74)EVLOCC
           READ (74)TEMP
           READ (74)EFRMI 
                        EFERMI(1)   =EFRMI
                        EFERMI(NSPN)=EFRMI
           READ (74)(((H(JWF,IWF,IX),JWF=1,NWF),IWF=1,NWF),IX=1,3)
           READ (74)CHARGE
           READ (74)E_UP,E_DN
           CLOSE(74)
           ELSE
           IF(.NOT.DMOM.OR..NOT.EXIST)THEN
           write(6,*)'DIPOLE MATRIX ELEMENTS UNAVAILABLE'
           write(6,*)'CHANGE DMOM=.TRUE. IN SPNORB'
           write(6,*)'ABANDONING JOINT DOS CALCULATION'
           RETURN
           END IF
       END IF
                    DO IWF=1,NWF
                    DO JWF=IWF+1,NWF
                    H(JWF,IWF,IX)=H(IWF,JWF,IX)
                    END DO
                    END DO
         CALL GTTIME(TIMEC)
         CALL TIMOUT('DOSJNT Marker 4                   :',TIMEC-TIMEB)
       END IF
C
C JOINT DOS CALCULATION
C CALCULATE <PSI-I | DIPOL | PSI-J>**2 AND SPECTRUM
C
       ESTEP=EWIND/(NSPEC-1)
       sos_pol = 0.0d0
       DO IPTS=1,NSPEC
        sos_freq(ipts) = 0.0d0
        SPTOT(IPTS)=0.0D0
        SPDIP(IPTS)=0.0D0
       END DO
       CALL GTTIME(TIMEA)
       DO IWF=1,NWF
        ISPN=1
        IF (IWF.GT.NWFS(1)) ISPN=2
        OCCI=FFERMI(EVLOCC(IWF),EFERMI(ISPN),TEMP)
        DO JWF=IWF,NWF
         JSPN=1
         IF (JWF.GT.NWFS(1)) JSPN=2
         OCCJ=FFERMI(EVLOCC(JWF),EFERMI(JSPN),TEMP)
         DIP=H(JWF,IWF,1)**2
     &      +H(JWF,IWF,2)**2
     &      +H(JWF,IWF,3)**2
          write(11,121)evlocc(iwf),evlocc(jwf),h(jwf,iwf,1),
     &             h(jwf,iwf,2),h(jwf,iwf,3),dip,occi,occj 
 121      format(1x,2(f7.4,1x), 4(f13.6), 2(f6.3,1x))
             IF(ISPN.NE.JSPN)DIP=0.0D0
         EDIFF=EVLOCC(JWF)-EVLOCC(IWF)
         EDIFF27=EDIFF*27.2118
         FCT=MAX(OCCI*(1.-OCCJ),OCCJ*(1.-OCCI))
         IF(FCT*DIP.GT.0.00001)THEN
             EV1=MIN(EVLOCC(IWF),EVLOCC(JWF))
             EV2=MAX(EVLOCC(IWF),EVLOCC(JWF))
             EDIFF=EV1-EV2
             sos_pol = sos_pol + dip/ediff
             
         PRINT 203,IWF,JWF,EV1,EV2,EDIFF,FCT,DIP
         END IF
 203     FORMAT(' ',2i5,6g15.6)
C
C LOOP OVER ENERGY GRID TO GET INTENSITIES
C
         BACKW= 1.0D0
         IF (IWF.EQ.JWF) BACKW= 0.0D0
          fact = occi*(1.0d0-occj)
         DO IPTS=1,NSPEC
          ERG=(IPTS-1)*ESTEP+ENJD
          FUNC=VFAC*EXP(-EALP*(ERG+EDIFF)**2)
          SPTOT(IPTS)=SPTOT(IPTS)+fact*FUNC
          SPDIP(IPTS)=SPDIP(IPTS)+fact*FUNC*DIP
          sos_freq(ipts) = sos_freq(ipts)+fact*DIP/erg
          FUNC=BACKW*VFAC*EXP(-EALP*(ERG+EDIFF)**2)
          SPTOT(IPTS)=SPTOT(IPTS)+OCCJ*(1.0D0-OCCI)*FUNC
          SPDIP(IPTS)=SPDIP(IPTS)+OCCJ*(1.0D0-OCCI)*FUNC*DIP
         END DO
        END DO
       END DO
       CALL GTTIME(TIMEB)
       CALL TIMOUT('SECTION 3                         :',TIMEB-TIMEA)
C
C CORRECT FOR SPIN EFFECTS AND PRINT SPECTRA
C
              DO IPTS=1,NSPEC
              SDMAX=MAX(SDMAX,SPDIP(IPTS))
              STMAX=MAX(STMAX,SPTOT(IPTS))
              END DO 
                ! Normalize the spectrum.
              DO IPTS=1,NSPEC
                SPDIP(IPTS)=SPDIP(IPTS)/SDMAX
                SPTOT(IPTS)=SPTOT(IPTS)/STMAX
              END DO
       write(6,*)'Unscreened polarizability (Bohr^3)', sos_pol
       write(6,*)'Unscreened polarizability ang^3',sos_pol*(0.5292**3)
       OPEN(73,FILE='JNTOUT',STATUS='UNKNOWN')
       REWIND(73)
       DO IPTS=1,NSPEC
        ERG=((IPTS-1)*ESTEP+ENJD)
        WRITE(73,1010)-ERG*HA2EV,SPDIP(IPTS),SPTOT(IPTS),sos_freq(ipts)
 1010   FORMAT(4(1X,F15.6))
       END DO
       CLOSE(73)
  900  CONTINUE
       CALL GTTIME(TIME2)
       CALL TIMOUT('JOINT DENSITY OF STATES:           ',TIME2-TIME1)
C
C DEALLOCATE LOCAL ARRAYS
C
       DEALLOCATE(H,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING H'
       DEALLOCATE(SPDIP,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING SPDIP'
       DEALLOCATE(SPTOT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING SPTOT'
       DEALLOCATE(SOS_FREQ,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING SOS_FREQ'
       DEALLOCATE(RVECA,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING RVECA'
       DEALLOCATE(PTS,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING PTS'
       DEALLOCATE(GRAD,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING GRAD'
       DEALLOCATE(ICOUNT,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING ICOUNT'
       DEALLOCATE(PSIG,STAT=IERR)
       IF(IERR/=0)WRITE(6,*)'DOSJNT:ERROR DEALLOCATING PSIG'

       RETURN
       END
