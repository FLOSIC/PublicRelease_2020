C UTEP Electronic Structure Lab (2020)
C
C *****************************************************************
C
C SUBROUTINE VSTART based on VHAM/KANDV by Briley/Pederson/Porezag
C Gets the potential part of the Hamiltonian: a fit contribution
C based on the expansion of the free atom potentials 
C or similar data read from STPOTUP/DN file. If (ISTSCF.EQ.0), a simple 
C superposition of the atomic start potentials will be used. 
C If (ISTSCF.EQ.3), the same superposition is used, but three-center 
C and crystal field integrals are neglected. In all other cases, 
C the routine does nothing.
C The Potential expansion is defined as follows:
C
C        NALPFIT           -AFIT(I)*R*R    -LDRFIT
C POT= -  SUM  CFIT(I) * E              * R     
C          I
C
C Default potentials are provided in COMMON/RHOPFIT/ (see routine
C ATOMSCF for details).
C For calculations employing nonlocal PSPs, the above expansion defines
C the local part of the atomic potential.
C If the default starting potentials are used, the starting potential 
C for each atom can be made spin-dependent by corresponding non-zero
C entries in COMMON/STPOT/ (1 for spin-up and -1 for spin-down).
C 
       SUBROUTINE VSTART(ALPHAI,ALPHAJ,AI,AJ,SS)
       use common2,only : RIDT, IFUIDT, NIDENT, ISPN
       use common3,only : RMAT
       use common4,only : ISPIDT, RPFALP, RPFCOF, NRPFIT, LDIVR
       use common5,only : ISTSCF
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:02:00 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: I, IFNCT, IID, IP, ISPOLD, IVEC, J, LDRFIT, MSITES,
     & NEQUIV, NFIT, NIDFIT
       REAL*8 :: ALPHAI , ALPHAJ, AI, AJ, SS, AIMAJSQ, AIMD1SQ,
     & AJMD1SQ, ARG, ATMP, B, BSQR, D, D1, GAUSSP, POTCOF, POTEXP,
     & RFIT, RVECA, TOLER, VDEL, WW
       SAVE
       DIMENSION AI(3),AJ(3),SS(10,10)
C
       LOGICAL NEGLECT,EXIST
       CHARACTER*7 NAMST
       DIMENSION POTEXP(MAXLSQF,MAX_IDENT),POTCOF(MAXLSQF,MAX_IDENT)
       DIMENSION RFIT(3,MAX_IDENT),LDRFIT(MAX_IDENT)
       DIMENSION RVECA(3,MX_GRP),MSITES(1)
       DIMENSION NFIT(MAX_IDENT)
       DIMENSION B(3),D(3),D1(3)
       DIMENSION WW(10,10)
C
C TOLER DEFINES WHEN TWO ATOMS ARE EQUIVALENT
C ISPOLD IS USED TO SAVE THE LAST SPIN
C
       DATA TOLER    /1.0D-5/
       DATA ISPOLD   /100000/
C
C DO NOT ZERO SS, WW ARRAYS: SS ALREADY CONTAINS KINETIC ENERGY INFO,
C WW IS ZEROED WITHIN THCNOV AND THCNDR
C
C RETURN IF NOTHING TO DO
C
       IF ((ISTSCF.NE.0).AND.(ISTSCF.NE.3)) RETURN
C
C START WITH SETUP
C
       IF (ISTSCF.EQ.0) THEN
        NEGLECT=.FALSE.
       ELSE
        NEGLECT=.TRUE.
       END IF
C
C SETUP DEFAULT POTENTIAL EXPANSION IF CALLED WITH NEW ISPN
C
       IF (ISPN.NE.ISPOLD) THEN
        ISPOLD=ISPN
        IF (ISPN.EQ.1) THEN
         NAMST='STPOTUP'
        ELSE
         NAMST='STPOTDN'
        END IF
        INQUIRE(FILE=NAMST,EXIST=EXIST)
        IF (EXIST) THEN
         OPEN(21,FILE=NAMST,FORM='FORMATTED',STATUS='OLD')
         REWIND(21)
        END IF
C
C READ START POTENTIAL INFO FROM FILE OR TAKE DEFAULT VALUES
C
        IF (EXIST) THEN
         PRINT '(4A,I1)','USING FILE ',NAMST,' TO DEFINE ',
     &                   'STARTING POTENTIAL FOR SPIN: ',ISPN
         READ(21,*,END=90) NIDFIT
         IF (NIDFIT .GT. MAX_IDENT) THEN
          write(6,*)'VSTART: MAX_IDENT MUST BE AT LEAST ',NIDFIT
          CALL STOPIT
         END IF
        ELSE
         NIDFIT=NIDENT
        END IF
C
C LOOP OVER ALL ATOMS
C FIRST, READ LOCATION OF FIT CENTER
C
        DO 20 IID=1,NIDFIT
         IF (EXIST) THEN
          READ(21,*,END=90)(RFIT(J,IID),J=1,3)
         ELSE
          DO J=1,3
           RFIT(J,IID)=RIDT(J,IID)
          END DO
         END IF
C
C EITHER USE DEFAULT POTENTIAL OR READ ZFIT, DECAY, GAUSSIAN DATA 
C
         IF (.NOT.EXIST) THEN
          IFNCT=IFUIDT(IID)
          NFIT(IID)=NRPFIT(IFNCT)
          LDRFIT(IID)=LDIVR(IFNCT)
          ATMP=RPFALP(IFNCT)
          DO I=1,NFIT(IID)
           POTEXP(I,IID)= ATMP
           POTCOF(I,IID)= -RPFCOF(2,I,IFNCT)
           ATMP= ATMP*2
          END DO
C
C CORRECT FOR SPIN-DEPENDENCE
C
          IF (NFIT(IID) .GT. 0) THEN
           I= MIN(2,NFIT(IID))
           VDEL= -0.1D0
           POTCOF(I,IID)= POTCOF(I,IID)+VDEL*(3-2*ISPN)*ISPIDT(IID)
          END IF
         ELSE 
          READ(21,*,END=90) NFIT(IID),LDRFIT(IID)
         END IF
         IF (NFIT(IID) .GT. MAXLSQF) THEN
          write(6,*)'VSTART: MAXLSQF MUST BE AT LEAST: ',NFIT(IID)
          CALL STOPIT
         END IF
         IF ((LDRFIT(IID) .LT. 0) .OR. (LDRFIT(IID) .GT. 1)) THEN
          write(6,*)'VSTART: LDRFIT MUST BE 0 OR 1'
          CALL STOPIT
         END IF
         IF (EXIST) THEN
          DO I=1,NFIT(IID)
           READ(21,*,END=90) POTEXP(I,IID),POTCOF(I,IID)
          END DO 
         END IF
   20   CONTINUE
        IF (EXIST) CLOSE(21)
        GOTO 100
C
C ERROR
C
   90   write(6,*)'VSTART: FILE ',NAMST,' IS BROKEN'
        CALL STOPIT
C
C WRITE STRTS(P)
C
  100   IF (ISPN.EQ.1) THEN
         NAMST='STRTS'
        ELSE
         NAMST='STRTSP'
        END IF
        OPEN(21,FILE=NAMST,FORM='FORMATTED',STATUS='UNKNOWN')
        REWIND(21)
        WRITE(21,*) NIDFIT
        DO 120 IID=1,NIDFIT
         WRITE(21,*) (RFIT(I,IID), I=1,3)
         WRITE(21,*) NFIT(IID),LDRFIT(IID)
         DO I=1,NFIT(IID)
          WRITE(21,*) POTEXP(I,IID),POTCOF(I,IID)
         END DO   
         WRITE(21,*) ' '
  120   CONTINUE
        WRITE(21,*) 'TO USE THE STARTING POTENTIAL FEATURE, RENAME'
        WRITE(21,*) 'STRTS/STRTSP TO STPOTUP/STPOTDN AND RERUN'
        CLOSE(21)
       END IF
C
C NOW THE REAL ACTION
C
       B(1)=AJ(1)-AI(1)
       B(2)=AJ(2)-AI(2)
       B(3)=AJ(3)-AI(3)
       BSQR=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
       AIMAJSQ=BSQR
       DO 280 IID=1,NIDFIT
        CALL GASITES(1,RFIT(1,IID),NEQUIV,RVECA,MSITES)
        DO 260 IP=1,NFIT(IID)
         GAUSSP=POTEXP(IP,IID)
         DO 250 IVEC=1,NEQUIV
          AIMD1SQ=0.0D0
          AJMD1SQ=0.0D0
          DO I=1,3
           D1(I)=RVECA(I,IVEC)
           D(I) =RVECA(I,IVEC)-AI(I)
           AJMD1SQ=AJMD1SQ+(D1(I)-AJ(I))**2
           AIMD1SQ=AIMD1SQ+D(I)*D(I)
          END DO  
          ARG=(ALPHAI*ALPHAJ*AIMAJSQ+ALPHAI*GAUSSP*AIMD1SQ
     &        +ALPHAJ*GAUSSP*AJMD1SQ)/(ALPHAI+ALPHAJ+GAUSSP)
C
C IF (NEGLECT), SKIP EVALUATION OF THREE-CENTER AND
C CRYSTAL FIELD INTEGRALS
C
          IF ((.NOT.NEGLECT).OR.(AIMD1SQ.LT.TOLER)
     &                      .OR.(AJMD1SQ.LT.TOLER)) THEN
           IF (ARG .LT. CUTEXP) THEN
            IF (LDRFIT(IID) .EQ. 1) THEN
             CALL THCNDR(ALPHAI,ALPHAJ,GAUSSP,B,D,WW)
            ELSE
             CALL THCNOV(ALPHAI,ALPHAJ,GAUSSP,B,D,WW)
            END IF
            DO I=1,10
             DO J=1,10
              SS(J,I)=SS(J,I)+WW(J,I)*POTCOF(IP,IID)
             END DO
            END DO
           END IF
          END IF
  250    CONTINUE
  260   CONTINUE
  280  CONTINUE
       RETURN
       END
