C UTEP Electronic Structure Lab (2020)
C
C ********************************************************************
C
       SUBROUTINE OPTVMSH(IOPT,IMSH,NRSIZ,IFNCT,RST,OVLTAB,ERRMAX,
     &                    NRAD,RRAD,WRAD,RTRAFO)
C
C CREATE OPTIMIZED RADIAL MESH FOR INTEGRATION IN THE INTERVAL (RST,INF)
C USE RGCMSH TO GET DATA 
C IOPT: 1: OPTIMIZE NRAD AND RTRAFO
C       2: LOOK FOR SMALLEST RTRAFO FOR WHICH THE ERROR IS LESS THAN
C          ERRMAX IF NRAD POINTS ARE USED (START AT INPUT RTRAFO)
C
       use common2,only : N_CON, LSYMMAX
! Conversion to implicit none.  Raja Zope Sun Aug 20 09:01:53 MDT 2017

!      INCLUDE  'PARAMAS'  
       INCLUDE  'PARAMA2'  
       INTEGER :: IOPT, IMSH, NRSIZ, IFNCT, NRAD, I, ICON, IMN, IRAD,
     & L, L1, MUL, NSEARCH, NTRUE
       REAL*8 :: SYMBOL , RST, OVLTAB, ERRMAX, RRAD, WRAD, RTRAFO, DIF,
     & ERRO, ERT1, ERT2, RMID, SUMTAB, TINY, WOVL
        SAVE
        PARAMETER (NSEARCH=50)
        DIMENSION OVLTAB(MAX_CON,LDIM)
        DIMENSION RRAD(NRSIZ),WRAD(NRSIZ)
        DIMENSION SUMTAB(2,MAX_CON,LDIM,2),WOVL(MAX_CON,LDIM)
        DATA TINY /1.0D-10/
C
C CHECK INPUT PARAMETERS
C
        IF (RST .LT. 0.0D0) THEN
         write(6,*)'OPTRMSH: RST MUST BE >= 0'
         CALL STOPIT
        END IF
        IF (ERRMAX .LT. TINY) THEN
         write(6,*)'OPTRMSH: DESIRED ACCURACY IS TOO SMALL: ',ERRMAX
         write(6,*)'         MUST NOT BE SMALLER THAN: ',TINY
         CALL STOPIT
        END IF
C
C SETUP OF WOVL 
C
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          WOVL(ICON,L1)= 1.0D0/OVLTAB(ICON,L1)
         END DO
        END DO
C
C GET TRUE VALUE OF INTEGRALS
C
        DO L= 0,LSYMMAX(IFNCT)
         L1= L+1
         DO ICON= 1,N_CON(L1,IFNCT)
          DO I= 1,2
           SUMTAB(I,ICON,L1,2)= 1.0D30
          END DO
         END DO
        END DO
C
        RMID= MAX(1.0D0,5*RST)
        NTRUE= 0
  100   CONTINUE
         NTRUE= NTRUE+10
         IF (NTRUE .GT. NRSIZ) THEN
          write(6,*)'OPTVMSH: NRSIZ MUST BE AT LEAST: ',NTRUE
          CALL STOPIT
         END IF
         DO L= 0,LSYMMAX(IFNCT)
          L1= L+1
          DO ICON= 1,N_CON(L1,IFNCT)
           DO I= 1,2
            SUMTAB(I,ICON,L1,1)= SUMTAB(I,ICON,L1,2)
           END DO
          END DO
         END DO
         CALL RGCMSH(IMSH,NTRUE,RST,RMID,RRAD,WRAD)
         CALL OPTRIGR(2,IFNCT,NTRUE,RRAD,WRAD,SUMTAB(1,1,1,2))
         DO L= 0,LSYMMAX(IFNCT)
          L1= L+1
          MUL= 4*L+2
          DO ICON= 1,N_CON(L1,IFNCT)
           DO I= 1,2
            ERRO= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &            *WOVL(ICON,L1)*MUL
            IF (ABS(ERRO) .GT. ERRMAX) GOTO 100
           END DO
          END DO
         END DO
        CONTINUE
        IF (IOPT .EQ. 2) GOTO 600
C
C MESH CREATION (IOPT .EQ. 1)
C
        ERRO= 0.0D0
        ERT1= 1.0D30
        ERT2= 1.0D30
        NTRUE= MAX(NTRUE,NRAD)
        DO 400 IRAD= NRAD,NTRUE
         DO 200 IMN= 1,NSEARCH
          CALL MINIMIZE(IMN,1.0D0,RMID,ERRO)
          CALL RGCMSH(IMSH,IRAD,RST,RST+RMID,RRAD,WRAD)
          CALL OPTRIGR(2,IFNCT,IRAD,RRAD,WRAD,SUMTAB(1,1,1,1))
          ERRO= 0.0D0
          DO L= 0,LSYMMAX(IFNCT)
           L1= L+1
           MUL= 4*L+2
           DO ICON= 1,N_CON(L1,IFNCT)
            DO I= 1,2
             DIF= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &            *WOVL(ICON,L1)*MUL
             ERRO= MAX(ERRO,ABS(DIF))
            END DO
           END DO
          END DO
C 
C IF THE ERROR FOR THE BEST MESH HAS CHANGED LESS THAN 1 PERCENT
C DURING THE LAST 5 ITERATIONS, STOP OPTIMIZING ALOPT
C
          ERT2= MIN(ERT2,ERRO) 
          IF (MOD(IMN,5) .EQ. 0) THEN 
           IF (ERT2 .GT. ERT1*0.99D0) GOTO 300
           ERT1= ERT2
          END IF
  200    CONTINUE
  300    IF (ERRO .LT. ERRMAX) GOTO 500
  400   CONTINUE
  500   CONTINUE
        NRAD= IRAD
        RTRAFO= RST+RMID
        RETURN
C
C MESH CREATION (IOPT .EQ. 2)
C
  600   IF (RTRAFO .LE. RST) THEN
         write(6,*)'OPTVMSH: RTRAFO MUST BE > RST'
         CALL STOPIT
        END IF
        RMID= RTRAFO-RST
  610   CONTINUE
         
         RMID= RMID*0.95D0
         CALL RGCMSH(IMSH,NRAD,RST,RST+RMID,RRAD,WRAD)
         CALL OPTRIGR(2,IFNCT,NRAD,RRAD,WRAD,SUMTAB(1,1,1,1))
         DO L= 0,LSYMMAX(IFNCT)
          L1= L+1
          MUL= 4*L+2
          DO ICON= 1,N_CON(L1,IFNCT)
           DO I= 1,2
            DIF= (SUMTAB(I,ICON,L1,1)-SUMTAB(I,ICON,L1,2))
     &           *WOVL(ICON,L1)*MUL
            IF (ABS(DIF) .GT. ERRMAX) GOTO 700
           END DO
          END DO
         END DO
         GOTO 610
  700   CONTINUE
        RTRAFO= RST+RMID/0.95D0
        CALL RGCMSH(IMSH,NRAD,RST,RTRAFO,RRAD,WRAD)
        RETURN
       END
