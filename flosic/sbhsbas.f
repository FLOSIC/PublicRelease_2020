C UTEP Electronic Structure Lab (2020)
C
C ******************************************************************
C
       SUBROUTINE SBHSBAS(IZNUC,ALP,CON,NALP,NBASF)
C
C DIRK POREZAG, AUGUST 1998
C BASIS SET GENERATOR FOR BHS PSEUDOPOTENTIAL CALCULATIONS
C DATA STORED IN COMMON BLOCK BDBHSBAS
C
C ATTENTION: PARAMETERS MXCHG,MXBAR,MXCON MUST BE IDENTICAL IN
C BDBHSBAS AND SBHSBAS
C
! Conversion to implicit none.  Raja Zope Thu Aug 17 14:35:01 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IZNUC, NALP, NBASF, MXBAR, MXCHG, MXCON, I, IADD,
     & IALP, IBOFS, IOFS, L, NADD, NFULL, NG, NMAX, NSPD
       REAL*8 :: SYMBOL , ALP, CON, ALPH, COEF
       PARAMETER (MXCHG=36)
       PARAMETER (MXBAR=10)
       PARAMETER (MXCON=1)
       COMMON/BHSBAS/ALPH(MXBAR,MXCHG),COEF(MXBAR,MXCON,3,MXCHG)
     &  ,NSPD(3,MXCHG),NG(MXCHG)
       DIMENSION ALP(MAX_BARE),CON(MAX_BARE,MAX_CON,3),NBASF(2,3)
       DIMENSION NFULL(3),NADD(2,3),IBOFS(2,3)
       SAVE
C
        IF (IZNUC.GT.MXCHG) THEN
         PRINT *,'SBHSBAS: AUTOMATED BASIS SETS NOT AVAILABLE '
     &          ,'FOR Z > ',MXCHG
         CALL STOPIT
        END IF
        IF (IZNUC.LT.0) THEN 
         PRINT *,'SBHSBAS: Z < 0 IS NOT SUPPORTED'
         CALL STOPIT
        END IF
C
C IF NUCLEAR CHARGE IS ZERO, WE JUST ADD A STANDARD BASIS
C
        IF (IZNUC.EQ.0) THEN
         NALP=4
         IF (NALP.GT.MAX_BARE) THEN
          PRINT *,'SBHSBAS: MAX_BARE MUST BE AT LEAST: ',NALP
          CALL STOPIT
         END IF
         ALP(1)=1.50D0
         ALP(2)=0.60D0
         ALP(3)=0.25D0
         ALP(4)=0.10D0
         DO L=1,3
          NADD(1,L)=NALP
          NADD(2,L)=0
          IBOFS(1,L)=0
          IBOFS(2,L)=0
         END DO
         NADD (1,3)=1
         IBOFS(1,3)=2
         NADD (2,3)=2
         IBOFS(2,3)=0
         NMAX=0
         DO L=1,3
          NFULL(L)=0
          NBASF(1,L)=NADD(1,L)
          NBASF(2,L)=NADD(2,L)
          NMAX=MAX(NMAX,NADD(1,L)+NADD(2,L))
         END DO
         IF (NMAX.GT.MAX_CON) THEN
          PRINT *,'SBHSBAS: MAX_CON MUST BE AT LEAST: ',NMAX
          CALL STOPIT
         END IF
         GOTO 100
        END IF
C
C HERE STARTS THE DEFINITION FOR REAL ATOMS
C
        NALP=NG(IZNUC)
        IF (NALP.GT.MAX_BARE) THEN
         PRINT *,'SBHSBAS: MAX_BARE MUST BE AT LEAST: ',NALP
         CALL STOPIT
        END IF
C
C DEFINE IBOFS AND NADD
C
        DO L=1,3
         NADD (1,L)=3
         IBOFS(1,L)=0
         NADD (2,L)=0
         IBOFS(2,L)=0
        END DO
        IBOFS(1,3)=1
        NADD (2,3)=1
        IBOFS(2,3)=0
        IF (NSPD(2,IZNUC).EQ.0) THEN
         IBOFS(1,2)=1
         NADD (2,2)=1
         IBOFS(2,2)=0
         NADD (1,3)=1
         IBOFS(1,3)=2
         NADD (2,3)=1
         IBOFS(2,3)=1
        END IF
C
C DEFINE NBASF, RANGE CHECKING
C
        NMAX=0
        DO L=1,3
         NFULL(L)=NSPD(L,IZNUC)
         NBASF(1,L)=NFULL(L)+NADD(1,L)
         NBASF(2,L)=NADD(2,L)
         NMAX=MAX(NMAX,NBASF(1,L)+NBASF(2,L))
        END DO
        IF (NMAX.GT.MAX_CON) THEN
         PRINT *,'SBHSBAS: MAX_CON MUST BE AT LEAST: ',NMAX
         CALL STOPIT
        END IF
        IF (NMAX.GT.NALP) THEN
         PRINT *,'SBHSBAS: INTERNAL ERROR: LINEAR DEPENDECY'
         CALL STOPIT
        END IF
        DO IALP=1,NALP
         ALP(IALP)=ALPH(IALP,IZNUC)
        END DO
        DO L=1,3
         DO I=1,NSPD(L,IZNUC)
          DO IALP=1,NALP
           CON(IALP,I,L)=COEF(IALP,I,L,IZNUC)
          END DO
         END DO
        END DO
C
C ADD ADDITIONAL FUNCTIONS TO OBTAIN BETTER BASIS
C
  100   DO L=1,3
         IOFS=NFULL(L)
         DO IADD=1,2
          IF ((NADD(IADD,L)+IBOFS(IADD,L)).GT.NALP) THEN
           PRINT *,'SBHSBAS: INTERNAL ERROR: NADD+IBOFS TOO LARGE'
           CALL STOPIT
          END IF
          DO I=1,NADD(IADD,L)
           DO IALP=1,NALP
            CON(IALP,IOFS+I,L)=0.0D0
           END DO
           IALP=NALP-IBOFS(IADD,L)-NADD(IADD,L)+I
           CON(IALP,IOFS+I,L)=1.0D0
          END DO
          IOFS=IOFS+NADD(IADD,L)
         END DO
        END DO
        RETURN
        END
