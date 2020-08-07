C UTEP Electronic Structure Lab (2020)
       SUBROUTINE JDOS 
C WRITTEN BY MARK R. PEDERSON
C
C MAGNETOANISOTROPY ENERGY BY MRP 18-MARCH 1999
C MAJOR CHANGE ON 2-AUGUST 1999
C CALCULATES SPIN-ORBIT COUPLING TERMS
C
       use for_diag1
       use common2,only : E_UP, E_DN, NSPN, DIPOLE
       use common3,only : RMAT
       use common5,only : PSI, NWF, NWFS, EFERMI, EVLOCC
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:34:50 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MXSPEC, NMAX, MDH, I, IBET, IDEBUG, IE, IENG, INFO,
     & ITHE, IWF, IWSPN, IX, J, JBAS, JWF, K, KWF, L, LL, MU, MWF,
     & NBET, NBETA, NDEG, NMAXJNT, NTHE, NTHET, NTMP, NU
       REAL*8 :: SYMBOL , AJDOS, AMX, ARG, AX, AXIS, AY, AZ, BETA,
     & BETMAX, BRD, BT, BXIS, CHARGE, CHGE, CLIGHT, CP, CSPD, DE, DIP,
     & DIP1, DIP2, DIP3, DIRLIST, DOS, DOSMAX, DOWNCHAR, ECORE, EDIF,
     & EDIFF, EF, EFRMI, ELEC, EMAX, EMAXJNT, EMIN, EMINJNT, EN, ENG,
     & EVLS, FCTLS, GRAD, H, HA2EV, HA2KEL, HT, OCC1, OCC2, ORBM, PI,
     & POT, PSIG, PTS, Q, REALTMP, RVECA, SC1R, SPDIP, SPEV, SPHM,
     & SPMT, SPOV, SPSC, SPTOT, TEMP, TH, THETA, THETMAX, TIME1, TIME2,
     & TRACE, UPCHAR, V, ZEXX, ZPVCHG
       SAVE
       COMPLEX*16 SPN,SPT,SPB,HAMEL,HAMLS,EVECC,SC1C,SC2C
       COMPLEX*16 FCT,FIJ
       COMPLEX*16 DIPMAT,OVLP,DIPC
       DIMENSION  ORBM(3),HT(3)
       PARAMETER (MXSPEC=10000)
       PARAMETER (NMAX=MPBLOCK)
       PARAMETER (MDH=MAX_OCC)
C
       LOGICAL ICOUNT,EXIST,FAST,DOIT,RDIT
       LOGICAL LMOM,LDIR,DMOM,LSQMOM
       DIMENSION SPN(2,2,3),SPT(2,2,3),SPB(2,2,3)
       CHARACTER*20 FNAME
       INTEGER IERR
       COMMON/TMP2/PSIG(4,NMAX,MAX_OCC),H(MAX_OCC,MAX_OCC,3)
       COMMON/TMP1/POT(MAX_PTS*MXSPN)
     &  ,SPDIP(MXSPEC),SPTOT(MXSPEC),RVECA(3,MX_GRP)
     &  ,PTS(NSPEED,3),GRAD(NSPEED,10,6,MAX_CON,3)
     &  ,ICOUNT(MAX_CON,3)
     &  ,NDEG(MAX_OCC)
     &  ,HAMLS(MDH*MDH),EVLS(MDH),SC1R(MDH)
     &  ,SC1C(2*MDH),SC2C(2*MDH)
     &  ,EVECC(MDH,MDH),DIRLIST(3,100),DIPMAT(MDH,MDH,3)
       DIMENSION Q(NMAX,3),V(NMAX),CP(NMAX,2),AXIS(3,2),BXIS(3,2)
       DIMENSION SPMT(2,2,3,3),SPHM(3,3),SPOV(3,3),SPEV(3),SPSC(6)
       DIMENSION ENG(MXSPEC),AJDOS(MXSPEC),DOS(MXSPEC)
       DIMENSION  IWSPN(MDH)
       DIMENSION OVLP(MDH,MDH)
       EQUIVALENCE (OVLP(1,1),HAMLS(1))
       CHARACTER*1 CHDUM  
       CHARACTER*80 CHLIN  
       SAVE
       DATA HA2EV/27.2116D0/
       DATA CLIGHT/137.088D0/
       DATA HA2KEL/315891.1826D0/        !27.2116*1.602/1.38E-04   
       DATA PI /3.141592654D0/
       DATA TEMP/1.0D-4/
       DATA IDEBUG/0/
C
       DOIT=.FALSE.
C      DOIT=.TRUE.
       RDIT=.FALSE.
       FNAME='SPNORB'
       OPEN(74,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
       REWIND(74)
       READ(74,*,END=33)DOIT,RDIT
  33   CONTINUE
       IF(.NOT.DOIT)THEN
        CLOSE(74)
        PRINT '(A)','MIX: CALCULATION OF SPNORB MATRIX HAS BEEN SKIPPED'
        RETURN
       END IF         
C
C
C READ IN NECESSARY INPUT DATA
C
       REWIND(74)
       LDIR=.FALSE.
       EMIN=-1.5D0 
       EMAX= 0.5D0
       NTHET=3
       NBETA=1
       THETMAX=0.5D0
       BETMAX =2.0D0
       REALTMP=1.0D-5 
       CSPD=1.0D0
       ZEXX=0.0D0
       CALL CORESPLIT(1,ECORE,EMIN,EMAX,ZPVCHG)
       DMOM=.TRUE.  !Dipole matrix elements
       EMINJNT=-1.
       EMAXJNT= 1.
       NMAXJNT=1000.
       READ(74,*,END=20)DOIT,RDIT
              INQUIRE(FILE='SPNDAT',EXIST=EXIST)
              IF(RDIT)THEN
                 IF(.NOT.EXIST)RDIT=.FALSE.
              END IF
       READ(74,*,END=20)EMIN,EMAX 
       READ(74,*,END=20)NTHET,NBETA,THETMAX,BETMAX
        IF(NTHET.LT.0) LDIR=.TRUE.
        IF(LDIR) THEN 
         NTHET=-NTHET
         NBETA=1
         IF(NTHET.GT.100) THEN
          write(6,*)'ERROR : MAXIMUM NUMBER OF DIRECTION < 100'
          CALL STOPIT
         ENDIF
         DO I=1,NTHET
          READ(74,*,END=22,ERR=22) (DIRLIST(J,I),J=1,3) 
         ENDDO
         GOTO 25
   22    CONTINUE
         write(6,*)'PROBLEM READING DIRECTIONS :-(  ',I
         CALL STOPIT
   25    CONTINUE
c	 ELSE
c         NTHET=5
c         NBETA=5     
         ENDIF

       READ(74,*,END=20)REALTMP,CSPD
       READ(74,*,END=20)ZEXX
       READ(74,*,END=20)LMOM,LSQMOM,DMOM
       READ(74,*,END=20)EMINJNT,EMAXJNT,NMAXJNT
   20  REWIND(74)
       WRITE(74,740)DOIT,RDIT
       WRITE(74,741) EMIN,EMAX
       IF(LDIR) THEN
        WRITE(74,742)-NTHET,NBETA,THETMAX,BETMAX  
        DO I=1,NTHET
         WRITE(74,*) (DIRLIST(J,I),J=1,3)
        ENDDO
       ELSE
        WRITE(74,742)NTHET,NBETA,THETMAX,BETMAX  
       ENDIF
       WRITE(74,743)REALTMP,CSPD
       WRITE(74,744)ZEXX
       WRITE(74,745)LMOM,LSQMOM,DMOM
       WRITE(74,755)EMINJNT,EMAXJNT,NMAXJNT
          CSPD=CLIGHT*CSPD
 740   FORMAT(2L8,14X,' T/F TO DO/(NOT DO) SPIN-ORBIT')
 741   FORMAT(2F15.6,' ABSOLUTE ENERGY WINDOW (HARTREES)') 
 742   FORMAT(2I5,2F9.4,3X,'NTHETA,NBETA, THETMAX,BETMAX (PI UNITS)')
 743   FORMAT(2F15.6,    ' REAL TEMPERATURE, SPEED OF LIGHT C=1')
 744   FORMAT(F15.6,16X,'EXTRA CHARGE ADDED')
 745   FORMAT(3L8,15X,'T/F CALCULATE <L> <LxLx>, <I|GRAD|J> (JDOS)')
 750   FORMAT(' ',2I5,2F14.4,G15.6)
 755   FORMAT(' ',2F14.4,I5,' EMINJNT,EMAXJNT,NMAXJNT (AU)')
       CLOSE(74)

       IF(EMIN.LT.-999.9) ECORE=-9999.9
       CALL CORESPLIT(2,ECORE,EMIN,EMAX,ZPVCHG)
       ZPVCHG=ZPVCHG+ZEXX
C
       IF(.NOT.RDIT) THEN
	 WRITE(6,*)'MAKE RDIT 1 IN SPNORB'
	 WRITE(6,*)'BYE BYE'
         RETURN
       END IF
C
 75    CONTINUE
           INQUIRE(FILE='SPNDAT',EXIST=EXIST)
           OPEN(74,FILE='SPNDAT',FORM='UNFORMATTED')
           READ (74)NWF,NWFS,NSPN
           READ (74)EVLOCC
           READ (74)TEMP
           READ (74)EFERMI 
           READ (74)(((H(JWF,IWF,IX),JWF=1,NWF),IWF=1,NWF),IX=1,3)
           READ (74)CHARGE
           READ (74)E_UP,E_DN
          CLOSE(74)
C
C  NOW CALCULATE 
C  SUM(I,J)<PHI_I S| VX | PHI_J S'><PHI_J S'| VY | PHI_I S>]/[E(IS)-E(JS))]
C
       OPEN(78,FILE='SO-DOS',STATUS='UNKNOWN')
       REWIND(78)
  110  READ(78,'(A)',END=120) CHDUM
       GOTO 110
  120  CONTINUE               
C
C
  900  CONTINUE
       DO I=1,MXSPEC
         AJDOS(I)=0.0D0
       END DO
       ZPVCHG=CHARGE

C NOW FOR EXACT DIAGONALIZATION:
c       DO ICALL=0,1    
c         IF(ICALL.EQ.0)THEN
c         FCTLS=0.0D0  
c         THETA=0.0D0
c         BETA=0.0D0
c         NTHE=1
c         NBET=1
c         ELSE
         FCTLS=1.0D0  
         NTHE=NTHET
         NBET=NBETA
         THETMAX=THETMAX*PI
         BETMAX=BETMAX*PI
c         END IF
       DO ITHE=1,NTHE
       DO IBET=1,NBET
c         IF(ICALL.EQ.1)THEN
          IF(LDIR)THEN
           CALL AX2THET(DIRLIST(1,ITHE),THETA,BETA)
          ELSE
           IF(NTHE.EQ.1)THEN
             THETA=0.0D0
           ELSE
             THETA=(ITHE-1)*THETMAX/(NTHE-1)
           ENDIF
           IF(NBET.EQ.1)THEN
             BETA =0.0D0
           ELSE
             BETA =(IBET-1)*BETMAX/(NBET-1)
           ENDIF
          ENDIF 
c         END IF
C
         CALL SPNMAT(SPN, SPT, SPB, THETA, BETA )
C
C CONSTRUCT HAMILTONIAN MATRICES:
         FCT=1.0D0/2.0D0/CSPD**2
         FCT=FCT*DCMPLX(0.0D0,-1.0D0)*FCTLS
         LL=0
         DO IWF=1,NWF
           IF(IWF.LE.NWFS(1))THEN
             NU=1
           ELSE
             NU=2
           END IF
           DO JWF=IWF,NWF
             HAMEL=DCMPLX(0.0D0,0.0D0)
             LL=LL+1
             IF(IWF.EQ.JWF)THEN
               HAMLS(LL)=EVLOCC(IWF)
             ELSE
               HAMLS(LL)=DCMPLX(0.0D0,0.0D0)
             END IF
             IF(JWF.LE.NWFS(1))THEN
               MU=1
             ELSE
               MU=2
             END IF
             IF(IWF.NE.JWF) THEN
              DO IX=1,3
               HAMEL=HAMEL+H(JWF,IWF,IX)*SPN(MU,NU,IX)
              END DO
             ENDIF
             HAMEL=HAMEL*FCT
             HAMLS(LL)=HAMLS(LL)+HAMEL
           END DO
         END DO
C DIAGONALIZE COMPLEX HAMILTONIAN
C OVERLAP=IDENTITY
      write(6,*)'NWF,NWFS:',NWF,NWFS
      OPEN(93,FILE='WHEREAMI')
      REWIND(93)
      WRITE(93,*)'CALLING ZHPEV'
      WRITE(93,*)'NWF,MDH:',NWF,MDH
      MWF=NWF
      CLOSE(93)
C     IF(LMOM.OR.DMOM.OR.LSQMOM) THEN
      CALL ZHPEV('V','L',MWF,HAMLS,EVLS,EVECC,MDH,SC1C,SC2C,INFO)
C for the Cray T3E 
C      CALL CHPEV('V','L',MWF,HAMLS,EVLS,EVECC,MDH,SC1C,SC2C,INFO)
C     ELSE
C     CALL ZHPEV('N','L',MWF,HAMLS,EVLS,EVECC,MDH,SC1C,SC2C,INFO)
C for the Cray T3E 
C      CALL CHPEV('N','L',MWF,HAMLS,EVLS,EVECC,MDH,SC1C,SC2C,INFO)
C     ENDIF
      IF(INFO.NE.0) THEN
       write(6,*)'PROBLEM IN SPIN-ORB: ZHPEV  ',INFO
       CALL STOPIT
      ENDIF
C
      OPEN(93,FILE='WHEREAMI')
      REWIND(93)
      WRITE(93,*)'DONE WITH ZHPEV'
      CLOSE(93)
         DO IWF=1,NWF 
           NDEG(IWF)=1
         END DO
         ELEC=DBLE(NINT(CHARGE))-2.0D0
         IF(ZPVCHG.GT.0.001D0)ELEC=ZPVCHG
         write(6,*)'ELEC:',ELEC
         CALL FERMILV(NWF,ELEC,EF,REALTMP,EVLS,SC1R,NDEG)

c          IF(DMOM.AND.ICALL.NE.0)THEN
c          DO IX=1,3
c                  DO LBAS=1,NWF
c                  DO IWF =1,NWF
c                  PSILI(IWF,LBAS)=CMPLX(0.0D0,0.0D0)
c                  END DO
c                  END DO
c             DO KBAS=1,NWF
c                IF(KBAS.LE.NWFS(1))THEN
c                 LBEG=1
c                 LEND=NWFS(1)
cc                ELSE
c                 LBEG=NWFS(1)+1
c                 LEND=NWF    
c                END IF
c             DO LBAS=LBEG,LEND
c                  HHH=H(MIN(LBAS,KBAS),MAX(LBAS,KBAS),IX)
c                  DO IWF=1,NWF
c                  PSILI(IWF,LBAS)=PSILI(IWF,LBAS)+EVECC(KBAS,IWF)*HHH
cc                  END DO
c             END DO
c             END DO
c             DO IWF =1  ,NWF
c             DO JWF =IWF,NWF
c                  WTO=SC1R(IWF)*(1.0D0-SC1R(JWF)) 
cc            IF(WTO.GE.1.0D-4)THEN
c              FIJ=CMPLX(0.0D0,0.0D0)
c             DO LBAS=1,NWF
c               FIJ=FIJ+CONJG(EVECC(LBAS,JWF))*PSILI(IWF,LBAS)
cc             END DO
c               WTO=WTO*FIJ*CONJG(FIJ)
c                DE=(EMAXJNT-EMINJNT)/NMAXJNT
c                EN=EMINJNT-DE
c                RD=0.05  
c               DO IENG=1,NMAXJNT
c                EN=EN+DE
c                ARG=(EVLS(IWF)-EVLS(JWF)-EN)**2
c                ARG=ARG/BRD**2
c                ARG=EXP(-ARG)*WTO
c                ENG(IENG)=EN!*27.2118
c                DOS(IENG)=DOS(IENG)+ARG
c               END DO
c             END IF
c             END DO
c             END DO
c          END DO
c          END IF
c         IF(ICALL.EQ.0)THEN
c           TRACE=0.0D0
c         ELSE
c           TRACE=-TNLS
c         END IF

         MWF=NINT(CHARGE)
         TRACE=0.0D0
         CHGE=0.0D0
         NTMP=0 
         DO IWF=1,NWF
           IF(SC1R(IWF).GT.1.0D-12)NTMP=NTMP+1
         ENDDO
         NTMP=NTMP+10
         DO KWF=1,NTMP
             CHGE=CHGE+SC1R(KWF)
             TRACE=TRACE+EVLS(KWF)*SC1R(KWF)
             UPCHAR=0.0D0
             DOWNCHAR=0.0D0
             DO JBAS=1,NWF 
               IF(JBAS.LE.NWFS(1)) THEN
               UPCHAR=UPCHAR+DCONJG(EVECC(JBAS,KWF))*EVECC(JBAS,KWF)
               ELSE
               DOWNCHAR=DOWNCHAR+DCONJG(EVECC(JBAS,KWF))*EVECC(JBAS,KWF)
               END IF
              END DO
              PRINT 8384,KWF,EVLS(KWF),SC1R(KWF),UPCHAR,DOWNCHAR
              WRITE(78,8384)KWF,EVLS(KWF),SC1R(KWF),UPCHAR,DOWNCHAR
 8384        FORMAT(' ',I5,G20.12,1F12.5,2F16.7)
 8385        FORMAT(' ',I5,G20.12,4F12.5)
         END DO

           IF(DMOM)   THEN  ! TUNNA
C          READ IN THE DIPOLE MATRIX FROM SPNDAT
           INQUIRE(FILE='SPNDAT',EXIST=EXIST)
           IF(EXIST) THEN
              OPEN(74,FILE='SPNDAT',FORM='UNFORMATTED')
              READ (74)NWF,NWFS,NSPN
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
              write(6,*)'DIPOLE MATRIX ELEMENTS UNAVAILABLE'
              write(6,*)'CHANGE DMOM=.TRUE. IN SPNORB'
              write(6,*)'ABANDONING JOINT DOS CALCULATION FROM SPNORB'
              RETURN
           END IF
              DO IWF=1,NWF
                IF (IWF.LE.NWFS(1)) THEN
                   IWSPN(IWF)=1
                ELSE
                   IWSPN(IWF)=2
                END IF
              END DO
              DO IX=1,3
               DO IWF=1,NWF
                DO JWF=1,NWF
                 DIPMAT(JWF,IWF,IX)=DCMPLX(0.0D0,0.0D0)
                END DO
               END DO
              END DO
C         STORE THE UPPER TRIANGLE OF H IN HAM
              ALLOCATE(AHAM(NWF,NWF),STAT=IERR)
              IF(IERR.NE.0)THEN
                WRITE(6,*)'jdos:Error allocating Ham'
              ENDIF 
              DO IX=1,3
                  DO IWF=1,NWF
                   DO JWF=1,NWF
                    AHAM(JWF,IWF)=0.0D0
                    OVLP(JWF,IWF)=DCMPLX(0.0D0,0.0D0)
                   END DO
                  END DO
                    DO IWF=1,NWF
                     DO JWF=1,IWF
                      AHAM(JWF,IWF)=H(JWF,IWF,IX)
                      AHAM(IWF,JWF)=AHAM(JWF,IWF)
                     END DO
                    END DO
                  DO IWF=1,NWF
                   DO K=1,NWF
                     DO L=1,NWF
                       IF(IWSPN(L).EQ.IWSPN(K)) THEN
                        OVLP(K,IWF)=OVLP(K,IWF)
     &                  +AHAM(K,L)*DCONJG(EVECC(L,IWF))
                       END IF
                     END DO    
                   END DO    
                  END DO

                 DO IWF=1,NWF
                   DO JWF=1,NWF
                     DO K=1,NWF
                        DIPMAT(JWF,IWF,IX)=DIPMAT(JWF,IWF,IX)
     &                           +OVLP(K,IWF)*EVECC(K,JWF)
                     END DO    
                   END DO    
                  END DO    
              END DO  ! IX
           DEALLOCATE(AHAM,STAT=IERR) 
           IF(IERR.NE.0)THEN
             WRITE(6,*)'jdos:Error deallocating Ham'
           ENDIF 
C          MAKE THE ENERGY GRID
           BRD=0.05
           DE=(EMAXJNT - EMINJNT)/NMAXJNT
           EN=EMINJNT-DE
           DO IENG=1,NMAXJNT
            EN=EN+DE
            ENG(IENG)=EN*27.2118
                 DO IWF=1,NWF
                   DO JWF=IWF,NWF
                    OCC1=SC1R(IWF)*(1.0D0-SC1R(JWF))
c                    OCC2=SC1R(IWF)*(1.0D0-SC1R(JWF))
                    EDIFF=EVLS(IWF)-EVLS(JWF)
                    ARG=(EDIFF-EN)**2
                    ARG=EXP(-ARG/BRD**2)
c                    OCM=MAX(OCC1,OCC2)
                    DIP=0.0D0
	       IF(OCC1.GT.0.001D0) THEN
               DIP=DCONJG(DIPMAT(JWF,IWF,1))*DIPMAT(JWF,IWF,1)+
     &             DCONJG(DIPMAT(JWF,IWF,2))*DIPMAT(JWF,IWF,2)+
     &             DCONJG(DIPMAT(JWF,IWF,3))*DIPMAT(JWF,IWF,3)
               END IF
               IF (IWF.EQ.JWF) DIP=0.0D0
               AJDOS(IENG)=AJDOS(IENG)+DIP*OCC1*ARG
              END DO
             END DO
            END DO
C      NOW WRITE OUT THE DIPOLE ELEMENTS
        WRITE(6,*) 'DIPOLE MATRIX'
          DO IWF=1,NWF
           DO JWF=IWF,NWF
            EDIF=DABS(EVLS(JWF)-EVLS(IWF))
            IF (EDIF.LE.EMAX) THEN
            OCC1=SC1R(IWF)*(1.0D0-SC1R(JWF))
            OCC2=SC1R(JWF)*(1.0D0-SC1R(IWF))
            DIP1=DCONJG(DIPMAT(JWF,IWF,1))*DIPMAT(JWF,IWF,1)
            DIP2=DCONJG(DIPMAT(JWF,IWF,2))*DIPMAT(JWF,IWF,2)
            DIP3=DCONJG(DIPMAT(JWF,IWF,3))*DIPMAT(JWF,IWF,3)
            DIP=DIP1+DIP2+DIP3
            IF (IWF.EQ.JWF)DIP=0.0D0
            WRITE(6,200) IWF,JWF,EDIF,DIP,MAX(OCC1,OCC2)
            END IF
           END DO
          END DO
          ENDIF  ! DMOM IF

 100     FORMAT(F13.7,4X,F16.8)
 200     FORMAT(2I5,3F16.8)
         TRACE=TRACE!*HA2KEL
         TH=THETA/PI
         BT=BETA /PI
         AX=SIN(THETA)*COS(BETA)
         AY=SIN(THETA)*SIN(BETA)
         AZ=COS(THETA)
         PRINT 92,TH,BT,TRACE,AX,AY,AZ
         WRITE(78,92)TH,BT,TRACE,AX,AY,AZ
         PRINT*
         WRITE(78,*)
       END DO  ! IBET
       END DO  ! ITHE
   
c       END DO  ! ICALL
         WRITE(6,*) 'JDOS'
         AMX= -1.0D+30
         DO IE=1,NMAXJNT
           AMX=MAX(AMX,AJDOS(IE))
         END DO
         WRITE(6,*) 'AMX = ', AMX
         DO IE=1,NMAXJNT
          AJDOS(IE)=AJDOS(IE)/AMX
          WRITE(6,100)ENG(IE),AJDOS(IE)
          WRITE(78,100)ENG(IE),AJDOS(IE)
         END DO
         CLOSE(75)
       
C WRITE OUT JNT DOS:
       IF(DMOM)THEN
             DOSMAX=0.0D0
              DO I=1,NMAXJNT
              IF(DOS(I).GT.DOSMAX)DOSMAX=DOS(I)
              END DO
              DO I=1,NMAXJNT
c              DOS(I)=DOS(I)/DOSMAX
              END DO
c             OPEN(78,FILE='JNTSPN')
             DO I=1,NMAXJNT
             IF(DOS(I).GT.0.0001)THEN
c             WRITE(78,96)ENG(I),DOS(I)
             END IF
             END DO
c       CLOSE(78)
       END IF
 90       FORMAT(' ',2G20.10,'   ',2G20.10)
 91       FORMAT(' ',10G15.6)
 92       FORMAT(2F10.5,1G12.4,3F8.5,3X,'TH,BT,SO,AXIS')
 940      FORMAT(3F10.5,28X,' <Lx>,<Ly>,<Lz>')
 941      FORMAT(3F10.5,28X,' <LxLx>,<LyLy>,<LzLz>')
 942      FORMAT(3F14.5,28X,' <L.Easy>,<L.Hard>,<L(L+1)>')
 95       FORMAT(3G12.4,' IMAG PART <LTOT>')
 96       FORMAT(2G15.6)
C
       CALL GTTIME(TIME2)
       CALL TIMOUT('SPIN ORBIT COUPLING:               ',TIME2-TIME1)
       RETURN        
       END
