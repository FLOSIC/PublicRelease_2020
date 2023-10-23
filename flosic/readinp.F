C UTEP Electronic Structure Lab (2020)
C
C **************************************************************
C
      SUBROUTINE READINP
C 
C DIRK POREZAG, JULY 1997
C READ ATOM DATA FROM INPUT AND ELECTRIC FIELD FROM EFIELD
C
C       use debug1
       use common1,only : ALPCOR, PSPSYM, ISITPSP, ISNLCC, NLCC
       use common2,only : RIDT, RCNT, IFUIDT, IFUCNT, NIDENT, NCNT,
     &   ZELC, ZNUC, BFCON, BFALP, N_BARE, N_CON, LSYMMAX, N_POS, NFNCT,
     &   E_UP, E_DN, WFFILE, NSPN, EFIELD
       use common3,only : RMAT
       use common4,only : ISPIDT
       use inputvars, only: simplified_input 
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:57 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, I_CON, IEXTRA, IFADD, IFNCT, IFOUND, IID, INUC,
     & IPOS, ISAV, J, L, MNUC
       REAL*8 :: ATMCL , DR, DUMMY, SAV, TOLSP, ZERO
       SAVE
C       COMMON/TMP1/RVEC(3,MX_GRP),NFADD(3),MSITES(1),IORD(2,MAX_BARE)
       CHARACTER*80 LINE
       LOGICAL EXIST,ANTIFERRO
       REAL*8 :: RVEC(3,MX_GRP)
       INTEGER :: NFADD(3),MSITES(1),IORD(2,MAX_BARE)
       DATA ZERO /1.0D-5/
       DATA TOLSP/1.0D-8/
       DATA ATMCL/1.0D-1/
C
C READ/WRITE EFIELD 
C
       !<LA: old_mode, use input module
       if(.not. simplified_input) then
         OPEN(21,FILE='EFIELD',FORM='FORMATTED',STATUS='UNKNOWN')
         REWIND(21)
         READ(21,*,END=20) EFIELD
   20    CONTINUE
         REWIND(21)
         WRITE(21,'(3(1X,D24.16))') EFIELD
         CLOSE(21)
       endif
C
C OPEN INPUT, READ NUMBER OF FUNCTION SETS 
C
       PRINT '(A)','READING ATOMIC AND BASIS SET DATA'
       INQUIRE(FILE='INPUT',EXIST=EXIST)
       IF (.NOT. EXIST) GOTO 500
       OPEN(21,FILE='INPUT',FORM='FORMATTED',STATUS='OLD')
       REWIND(21)
       READ(21,*,END=500) NFNCT
       IF (NFNCT .GT. MAX_FUSET) THEN
        write(6,*)'READINP: MAX_FUSET MUST BE AT LEAST: ',NFNCT
        GOTO 600
       END IF
       
C
C CHARGE, PSEUDOPOTENTIAL TYPE, AND LOCATIONS
C
       NIDENT= 0
       NCNT= 0
       ISITPSP= 0 
       ISNLCC= 0
       ANTIFERRO= .FALSE.
       DO 100 IFNCT=1,NFNCT
        READ(21,*,END=500) ZELC(IFNCT),ZNUC(IFNCT)
        READ(21,'(A7)',END=500) PSPSYM(IFNCT)
        IF (PSPSYM(IFNCT)(1:3) .NE. 'ALL') ISITPSP=1
        READ(21,*,END=500) N_POS(IFNCT)
        DO 50 IPOS=1,N_POS(IFNCT)
         NIDENT=NIDENT+1
         IF (NIDENT .GT. MAX_IDENT) THEN
          write(6,*)'READINP: MAX_IDENT MUST BE AT LEAST: ',NIDENT
          GOTO 600
         END IF
         IFUIDT(NIDENT)=IFNCT
         ISPIDT(NIDENT)=0
         READ(21,'(A80)') LINE
         DO I=1,78
          IF (LINE(I:I+2) .EQ. 'SUP') THEN
           ISPIDT(NIDENT)= 1
           ANTIFERRO= .TRUE.
           GOTO 30
          ENDIF
          IF (LINE(I:I+2) .EQ. 'SDN') THEN
           ISPIDT(NIDENT)= -1
           ANTIFERRO= .TRUE.
           GOTO 30
          ENDIF
         END DO
   30    READ(LINE,*,END=500)(RIDT(J,NIDENT),J=1,3)
C
C SET UP COORDINATES AND TYPE FOR EQUIVALENT ATOMS
C
         CALL GASITES(1,RIDT(1,NIDENT),MNUC,RVEC,MSITES)
         IF (NCNT+MNUC .GT. MX_CNT) THEN
          write(6,*)'READINP: MX_CNT MUST BE AT LEAST: ',NCNT+MNUC
          GOTO 600
         END IF
         DO INUC=1,MNUC
          IFUCNT(NCNT+INUC)= IFNCT
          RCNT(1,NCNT+INUC)= RVEC(1,INUC)
          RCNT(2,NCNT+INUC)= RVEC(2,INUC)
          RCNT(3,NCNT+INUC)= RVEC(3,INUC)
         END DO
         NCNT=NCNT+MNUC
C
C CHECK IF ATOMS ARE TOO CLOSE
C
         DO IID=1,NIDENT
          IFOUND=0
          DO INUC=1,MNUC
           DR= ABS(RVEC(1,INUC)-RIDT(1,IID))
     &        +ABS(RVEC(2,INUC)-RIDT(2,IID))
     &        +ABS(RVEC(3,INUC)-RIDT(3,IID))
           IF (DR .LT. ATMCL) IFOUND=IFOUND+1
          END DO
          IF (IID .EQ. NIDENT) IFOUND= IFOUND-1
          IF (IFOUND .NE. 0) THEN
           write(6,*)'READINP: IDENTITY MEMBERS ARE TOO CLOSE TO ',
     &             'EACH OTHER:'
           write(6,*)'ID MEMBER 1: ',IID,', ID MEMBER 2: ',NIDENT
           PRINT 1000,'LOCATION 1: ',(RIDT(J,NIDENT), J=1,3)
           PRINT 1000,'LOCATION 2: ',(RIDT(J,IID),    J=1,3)
           WRITE(6,*)'IFOUND = ', IFOUND
 1000      FORMAT(A,3(1X,F15.5))
           GOTO 600
          END IF
         END DO
   50   CONTINUE
C
C BASIS SET INFO
C
        READ(21,*,END=500) IEXTRA
        READ(21,*,END=500) N_BARE(IFNCT)
        IF (N_BARE(IFNCT) .GT. MAX_BARE) THEN
        write(6,*)'READINP: MAX_BARE SHOULD BE AT LEAST: ',N_BARE(IFNCT)
         GOTO 600
        END IF
        READ(21,*,END=500)(N_CON(L,IFNCT),L=1,3)
        READ(21,*,END=500)(NFADD(L),L=1,3)
        READ(21,*,END=500)(BFALP(I,IFNCT),I=1,N_BARE(IFNCT))
C
C SORT ALPHA FROM BIG TO SMALL
C
        DO I=1,N_BARE(IFNCT)
         IORD(1,I)=I
        END DO
        DO I=1,N_BARE(IFNCT)
         DO J=I+1,N_BARE(IFNCT)
          IF (BFALP(I,IFNCT) .LT. BFALP(J,IFNCT)) THEN
           SAV= BFALP(I,IFNCT)
           BFALP(I,IFNCT)= BFALP(J,IFNCT)
           BFALP(J,IFNCT)= SAV
           ISAV= IORD(1,I)
           IORD(1,I)= IORD(1,J)
           IORD(1,J)= ISAV
          END IF
         END DO
        END DO
        DO I=1,N_BARE(IFNCT)
         DO J=1,N_BARE(IFNCT)
          IF (IORD(1,J) .EQ. I) IORD(2,I)=J
         END DO
        END DO
        LSYMMAX(IFNCT)=0
        DO L=1,3
         IF (IEXTRA .GT. 0) N_CON(L,IFNCT)=N_CON(L,IFNCT)+NFADD(L)
         IF (N_CON(L,IFNCT) .GT. 0) LSYMMAX(IFNCT)=L-1
         IF (N_CON(L,IFNCT) .GT. MAX_CON) THEN
         write(6,*)'READINP: MAX_CON MUST BE AT LEAST: ',N_CON(L,IFNCT)
          GOTO 600
         END IF
         DO I_CON=1,N_CON(L,IFNCT)
          READ(21,*,END=500)
     &        (BFCON(IORD(2,I),I_CON,L,IFNCT),I=1,N_BARE(IFNCT))
         END DO
         IF (IEXTRA .LE. 0) THEN
          DO IFADD=1,NFADD(L)
           READ(21,*,END=500)(DUMMY, I=1,N_BARE(IFNCT))
          END DO
         END IF
        END DO 
        ALPCOR(IFNCT)=0.0D0
  100  CONTINUE
C
C ELECTRONS AND NAME FOR FILE FOR READING/WRITING WAVEFUNCTIONS
C DEFINE NSPN
C
       READ(21,*,END=500) E_UP,E_DN
       NSPN=2
       IF ((ABS(E_UP-E_DN) .LT. TOLSP) .AND. (.NOT. ANTIFERRO)) NSPN=1
       IF (NSPN .GT. MXSPN) THEN
        write(6,*)'READINP: MXSPN MUST BE AT LEAST: ',NSPN
        CLOSE(21)
        CALL STOPIT
       END IF
       E_UP=ABS(E_UP)
       E_DN=ABS(E_DN)
       IF ((E_UP+E_DN) .LT. ZERO) THEN
        write(6,*)'READINP: THERE IS NO POINT IN WASTING CPU TIME FOR'
        write(6,*)'CALCULATIONS WITHOUT ELECTRONS. GET A LIFE !'
        GOTO 600
       END IF
       WFFILE='WFOUT'
       READ(21,'(A)') WFFILE
       CLOSE(21)
C
C INITIALIZE PSEUDOPOTENTIAL PARAMETERS
C
       CALL READPSP
C
C GET ATOMIC POTENTIALS AND DENSITIES
C
       DO IFNCT= 1,NFNCT
        CALL ATOMSCF(IFNCT)
       END DO
       RETURN
C
  500  write(6,*)'READINP: FILE INPUT IS MISSING OR BROKEN'
  600  CLOSE(21)
       CALL STOPIT
C
      END
