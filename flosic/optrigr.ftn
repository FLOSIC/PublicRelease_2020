C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE OPTRIGR(MODE,IFNCT,NRAD,RRAD,WRAD,SUMTAB)
C
C SUBROUTINE THAT DETERMINES THE NUMERICAL INTEGRALS NEEDED BY OPTRMSH,
C OPTVMSH, AND OPTMMSH
C MODE: 1 IF WEIGHTS NEED TO BE MULTIPLIED BY 4*PI*RRAD**2
C       2 IF WEIGHTS NEED TO BE MULTIPLIED BY 4*PI
C
       use common2,only : BFCON, BFALP, N_BARE, N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: MODE, IFNCT, NRAD, IBARE, ICON, IR, IRAD, L, L1, NPV
       REAL*8 :: RRAD , WRAD, SUMTAB, FAC, PI, R2, SUM, TABEXP, VLSAV,
     & W
        SAVE
        DIMENSION RRAD(NRAD),WRAD(NRAD),SUMTAB(2,MAX_CON,LDIM)
        DIMENSION VLSAV(NSPEED),TABEXP(MAX_BARE)
C
        PI= 4*ATAN(1.0D0)
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          SUMTAB(1,ICON,L1)= 0.0D0
          SUMTAB(2,ICON,L1)= 0.0D0
         END DO
        END DO
C
        DO IRAD= 0,NRAD-1,NSPEED
         NPV= MIN(NSPEED,NRAD-IRAD)
         CALL VLOCAL(1,NPV,IFNCT,RRAD(IRAD+1),VLSAV)
         DO IR= 1,NPV
          R2= RRAD(IRAD+IR)**2
          DO IBARE= 1,N_BARE(IFNCT)
           FAC= BFALP(IBARE,IFNCT)*R2
           IF (FAC .GT. CUTEXP) THEN
            TABEXP(IBARE)= 0.0D0
           ELSE
            TABEXP(IBARE)= EXP(-FAC)
           END IF
          END DO
          W= 4*PI*WRAD(IRAD+IR)
          IF (MODE .EQ. 1) W= W*R2
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           DO ICON= 1,N_CON(L1,IFNCT)
            SUM= 0.0D0
            DO IBARE= 1,N_BARE(IFNCT)
             SUM= SUM+BFCON(IBARE,ICON,L1,IFNCT)*TABEXP(IBARE)
            END DO
            SUM= W*SUM*SUM
            SUMTAB(1,ICON,L1)= SUMTAB(1,ICON,L1)+SUM
            SUMTAB(2,ICON,L1)= SUMTAB(2,ICON,L1)+SUM*VLSAV(IR)
           END DO
           W= W*R2
          END DO
         END DO
        END DO
        RETURN
       END
