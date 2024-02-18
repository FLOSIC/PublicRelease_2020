C UTEP Electronic Structure Lab (2020)
       SUBROUTINE WRITEUS(ISTEP)
C
C     ------------------------------------------------------------------
C
C      This subroutine does the writing of the MO coefficients 
C      into the molden file and allows molden to graphically display
C      the results.
C
C      BY MARK PEDERSON AND ULISES REVELES (JUNE, 2013)
C
C      INTRODUCED A NEW PARAMETER 'IHOMOP'
C      TOTAL NUMBER OF MO IS 2*IHOMOP+1, HOMO STATE IS IN THE MIDDLE
C
C     ------------------------------------------------------------------
C 
C      COMMON VARIABLES
C
      use xmol,only    : NUM_ATMS,XMOL_LIST,GET_LETTER,ATOMORB_LAB
      use common2,only : EGMIN, EGMAX, RIDT, LSYMMAX, N_CON, NFNCT,
     &                   N_POS, NSPN, NATOMS, INDICES
      use common5,only : PSI, EVLOCC, NWF, NWFS, OCCUPANCY, EFERMI
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:07 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: ISTEP, I, I_POS, ICON, IDIM, IFAIL, IFNCT, ILOC,
     & INUC, ISHDUMMY, ISHELLA, ISIZE, IWF, J, J_POS, JCON, JLOC, JSPN,
     & KLOC, KNT, L_NUC, LI, LMX1, M_NUC, MCON, MU, NBAS, NU
       REAL*8 :: AU2ANG , EFMAX, EFMIN, ERGMAX, ERGMIN, OCCP, PSIKH,
     & SPN
C
      REAL*8 EPS
      PARAMETER (EPS = 1.E-6)
      real*8 :: ffermi
      external ffermi
C
C     --- LOCAL WORK FIELDS ---
C
       DIMENSION ISIZE(3),SPN(2)
       DIMENSION KLOC(10,6,3),PSIKH(6)
C
       CHARACTER*1 PRN
       CHARACTER*2 AZI,LETTER
       CHARACTER*3 IATOM
       CHARACTER*15 ORBNME(6)
       CHARACTER*80 LINE,STRCOMP
       INTEGER K,STREXT,NBASTOT
       REAL*8 TEMP
C
       LOGICAL BETAMOS
C
       DATA ISIZE/1,3,6/
       DATA AU2ANG/0.529177D0/
C
C      LOCAL DYNAMICAL FIELDS 
C
       INTEGER ALLOCATION
       REAL*8,ALLOCATABLE  :: RVECA(:,:)
       LOGICAL EXIST
C
C     ------------------------------------------------------------------
C
C     --- CALCULATE ENERGY TRESHOLDS ---
C
      EFMIN=MIN(EFERMI(1),EFERMI(NSPN))
      EFMAX=MAX(EFERMI(1),EFERMI(NSPN))
C
C     --- SET THRESHOLDS HERE FOR NOW ---
C 
      ERGMIN=EFMIN-(EGMIN)
      ERGMAX=EFMAX+(EGMAX)
C
      CALL WFWIND2(ERGMIN,ERGMAX,.TRUE.,.TRUE.,IFAIL,NBASTOT)
C
      BETAMOS = .FALSE.
C
C     --- ALLOCATE LOCAL FIELDS ---
C
      IF (NWF.GT.NWFS(1)) THEN
        IDIM = 2*NWFS(1)
        BETAMOS = .TRUE.
      ELSE
        IDIM = NWFS(1)
      END IF
C
      ALLOCATE(RVECA(3,MX_GRP),ATOMORB_LAB(NDH),STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: ALLOCATION FAILED IN ','WRITE_UNSYM'
      END IF
C
C     --- OPEN AND APPEND MOLDEN FILE ---
C
      OPEN(43,FILE='CLUSTER.MOLDEN',FORM='FORMATTED',STATUS='OLD',
     &     ACCESS='APPEND')
C
C     --- WRITE ENERGY, SPIN, AND OCCUPATIONS ---
C
      WRITE(43,37)
 37   FORMAT('[MO]')
C
C     --- RUN OVER REQUESTED MOLECULAR ORBITALS ---
C
      NBAS = NWFS(1)
C
      TEMP=0.0001D0
      DO 87 IWF=1,NWF
C
C     --- FIND OUT IF BETA ORBITALS ARE NEEDED ---
C 
        JSPN=1
        IF(IWF.GT.NWFS(1))JSPN=2
C 
        KNT = 0
C
C     --- PRINT LINE WITH ORBITAL INFORMATION:             ---
C     --- INDICES(1-4,IWF)= NUMBER, K_REP, IROW, AND ISPN  ---
C
        INDICES(1,IWF) =  IWF
        LINE = ''
        LINE(1:4) = 'Sym='
        DO J=1,4
          K = STREXT(STRCOMP(LINE))
          WRITE(LINE(K+1:K+6),'(I5)')INDICES(J,IWF)
          IF (J.LT.4) WRITE(LINE(K+7:K+8),'(A1)')'-'
          LINE = STRCOMP(LINE) 
        END DO
        WRITE(43,'(A80)') LINE
C        
        WRITE(43,39)EVLOCC(IWF),EFERMI(JSPN)
C
        IF(IWF.LE.NWFS(1))THEN
          OCCP=FFERMI(EVLOCC(IWF),EFERMI(1),TEMP)
          WRITE(43,40)'Alpha'
C
C     --- OCCUPATION FOR CLOSED SHELL SYSTEMS ---
C
          IF (.NOT.BETAMOS) OCCP = 2.0*OCCP
C
        ELSE
          OCCP=FFERMI(EVLOCC(IWF),EFERMI(JSPN),TEMP)
          WRITE(43,40)'Beta '
        END IF
C
C     --- NUMERICAL STABILITY TEST ---
C
        IF (ABS(OCCP).LT.EPS) OCCP = 0.0
C
        WRITE(43,41)OCCP!OCCUPANCY(IWF)
 39     FORMAT('Ene=',2G14.6)
 40     FORMAT('Spin= ',A5) 
 41     FORMAT('Occup= ',G14.6)
C
        ISHELLA=0
        DO 86 IFNCT=1,NFNCT
          LMX1=LSYMMAX(IFNCT)+1
          DO 84 I_POS=1,N_POS(IFNCT)
            ISHELLA=ISHELLA+1
            CALL OBINFO(1,RIDT(1,ISHELLA),RVECA,M_NUC,ISHDUMMY)
            DO 82 J_POS=1,M_NUC
              CALL UTRAVEL(IFNCT,ISHELLA,J_POS,RIDT(1,ISHELLA),
     &                     RVECA,L_NUC,1)
              ILOC=0
              DO LI=1,LMX1
                DO MU=1,ISIZE(LI)
                  DO ICON=1,N_CON(LI,IFNCT)
                    ILOC=ILOC+1
                    KLOC(ICON,MU,LI)=ILOC
                  END DO
                END DO
              END DO
C
              DO LI=1,LMX1
                DO MCON=1,N_CON(LI,IFNCT)
                  DO MU=1,ISIZE(LI)
                    NU=0
                    DO ICON=1,N_CON(LI,IFNCT)
                      JCON=ICON
C 
C     --- LABELS FOR 'S' ATOMIC ORBITALS ---
C 
                      IF(LI.EQ.1)THEN
                        PRN='S'
                        AZI=' '
C
C     --- LABELS FOR 'P' ATOMIC ORBITALS ---
C
                      ELSE  IF(LI.EQ.2)THEN
                        PRN='P'
                        IF(MU.EQ.1)AZI='X '
                        IF(MU.EQ.2)AZI='Y '
                        IF(MU.EQ.3)AZI='Z '
                        JCON=JCON+1
C
C     --- LABELS FOR 'D' ATOMIC ORBITALS ---
C
                      ELSE IF(LI.EQ.3)THEN
                        JCON=JCON+2
                        PRN='D'
                        IF(MU.EQ.1)AZI='XX'
                        IF(MU.EQ.2)AZI='YY'
                        IF(MU.EQ.3)AZI='ZZ'
                        IF(MU.EQ.4)AZI='XY'
                        IF(MU.EQ.5)AZI='XZ'
                        IF(MU.EQ.6)AZI='YZ'
                      END IF
C
                      IF(ICON.EQ.MCON)THEN
                        JLOC=KLOC(ICON,MU,LI)
                        NU=MU   
                        PSIKH(NU)=PSI(JLOC,IWF,1)
                        WRITE(ORBNME(NU),42)JCON,PRN,AZI,JLOC
C                       WRITE(43,43)JCON,PRN,AZI,JLOC,PSI(JLOC,IWF,1)
                      END IF
 42                   FORMAT(I2,A1,A2,I5)        
 43                   FORMAT(I2,A1,A2,I5,F12.6)
                    END DO !ICON
                  END DO !MU
C
C     --- NOW WRITE 'S', 'P', AND 'D' ORBIALS ---
C
                  DO MU=1,ISIZE(LI)
                    KNT=KNT+1
                    ATOMORB_LAB(KNT)%ORBLAB = ORBNME(MU)(1:5)
                    WRITE(43,44)KNT,PSIKH(MU)
 44                 FORMAT(I10,F12.6)
                  END DO
C
                END DO
              END DO
   82       CONTINUE
   84     CONTINUE
   86   CONTINUE
   87 CONTINUE
C
C     --- CLOSE MOLDEN FILE ---
C
      CLOSE(43)
C
C     --- ADD ATOMIC LABEL TO ATOMIC ORBITAL LABEL ---
C
      INUC = 0
      DO I=1,NBASTOT
        IF (ATOMORB_LAB(I)%ORBLAB(1:3).EQ.' 1S') INUC = INUC + 1
        CALL GET_LETTER(XMOL_LIST(INUC)%ANUM,LETTER)
        ATOMORB_LAB(I)%ATNUM = INUC
        ATOMORB_LAB(I)%ATLAB = LETTER
      END DO
CJUR
C     DO I=1,NBASTOT
C       WRITE(*,*) I,'ATOMORB_LAB(I)%ORBLAB(3:5):',
C    %              ATOMORB_LAB(I)%ORBLAB(3:5)
C     END DO
CJUR
C     --- DEALLOCATE LOCAL FIELDS ---
C
      DEALLOCATE(RVECA,STAT=ALLOCATION)
      IF (ALLOCATION.GT.0) THEN
        WRITE(*,*) 'ERROR: DEALLOCATION FAILED IN ','WRITE_UNSYM'
      END IF
C
      RETURN
C
C     ------------------------------------------------------------------
C
      END
